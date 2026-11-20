
/* ikchain: analytic chain IK with twist distribution and joint-limit clamping.
 *
 * ik_solve pipeline:
 *   1. validate definition (reject invalid length/limit params)
 *   2. project target onto the reachable shell
 *   3. analytic positions: recursive binary split of the chain, each split
 *      solved as a two-bone triangle by the law of cosines, pole vector
 *      picks the bend side
 *   4. swing limits: project each bone direction onto the cone around the
 *      previous bone (sequential root->tip, positions rebuilt)
 *   5. twist: distribute a total twist angle across joints; joint k rotates
 *      its sub-chain about the hinge->end axis (P[k-1] and the end both lie
 *      on the axis, so bone lengths and the end position cannot drift),
 *      per-joint twist clamped to its range
 *   6. final swing clamp so no joint is out of range after twisting
 */
#include <math.h>
#include <stdio.h>
#include <string.h>

#define IK_PI 3.14159265358979323846f
#define IK_MAX_BONES 16

typedef struct { float x, y, z; } V3;

static V3 v3(float x, float y, float z) { V3 v = { x, y, z }; return v; }
static V3 vadd(V3 a, V3 b) { return v3(a.x + b.x, a.y + b.y, a.z + b.z); }
static V3 vsub(V3 a, V3 b) { return v3(a.x - b.x, a.y - b.y, a.z - b.z); }
static V3 vscale(V3 a, float s) { return v3(a.x * s, a.y * s, a.z * s); }
static float vdot(V3 a, V3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static V3 vcross(V3 a, V3 b) { return v3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x); }
static float vlen(V3 a) { return sqrtf(vdot(a, a)); }
static V3 vnorm(V3 a) { float l = vlen(a); return l > 1e-12f ? vscale(a, 1.0f / l) : v3(0, 0, 0); }
static float clampf(float x, float lo, float hi) { return x < lo ? lo : (x > hi ? hi : x); }

/* Rodrigues rotation of v about a unit axis */
static V3 vrot(V3 v, V3 axis, float ang) {
    float c = cosf(ang), s = sinf(ang);
    return vadd(vadd(vscale(v, c), vscale(vcross(axis, v), s)),
                vscale(axis, vdot(axis, v) * (1.0f - c)));
}
static V3 rot_about_line(V3 p, V3 lp, V3 axis, float ang) {
    return vadd(lp, vrot(vsub(p, lp), axis, ang));
}
static V3 any_perp(V3 d) {
    V3 up = fabsf(d.y) < 0.9f ? v3(0, 1, 0) : v3(1, 0, 0);
    return vnorm(vcross(d, up));
}
static float dir_angle(V3 a, V3 b) {
    return acosf(clampf(vdot(vnorm(a), vnorm(b)), -1.0f, 1.0f));
}

typedef struct {
    float len;                  /* bone length, > 0 */
    float swing_max;            /* max bend vs previous bone, [0, pi] */
    float twist_min, twist_max; /* twist range of the joint at this bone's start */
    float twist_w;              /* twist distribution weight, >= 0 */
} IkBone;

typedef struct {
    int n;
    IkBone bone[IK_MAX_BONES];
} IkChainDef;

typedef struct {
    V3 p[IK_MAX_BONES + 1];     /* joint positions, p[0] = root, p[n] = end */
    float twist[IK_MAX_BONES];  /* twist applied per joint (rad) */
} IkPose;

enum { IK_OK = 0, IK_ERR_N, IK_ERR_LEN, IK_ERR_SWING, IK_ERR_TWIST, IK_ERR_WEIGHT };

int ik_validate(const IkChainDef* d) {
    if (!d || d->n < 1 || d->n > IK_MAX_BONES) return IK_ERR_N;
    float wsum = 0.0f;
    for (int i = 0; i < d->n; i++) {
        const IkBone* b = &d->bone[i];
        if (!isfinite(b->len) || !(b->len > 0.0f)) return IK_ERR_LEN;
        if (!isfinite(b->swing_max) || !(b->swing_max >= 0.0f && b->swing_max <= IK_PI)) return IK_ERR_SWING;
        if (!isfinite(b->twist_min) || !isfinite(b->twist_max) || b->twist_min > b->twist_max) return IK_ERR_TWIST;
        if (!isfinite(b->twist_w) || !(b->twist_w >= 0.0f)) return IK_ERR_WEIGHT;
        if (i >= 1) wsum += b->twist_w;
    }
    if (d->n > 1 && !(wsum > 0.0f)) return IK_ERR_WEIGHT;
    return IK_OK;
}

/* two-bone triangle by the law of cosines; pole picks the bend side */
static V3 two_bone_mid(V3 root, V3 target, float l1, float l2, V3 pole) {
    V3 d = vsub(target, root);
    float dist = vlen(d);
    V3 dir = dist > 1e-9f ? vscale(d, 1.0f / dist) : v3(1, 0, 0);
    float eps = 1e-6f * (l1 + l2);
    dist = clampf(dist, fabsf(l1 - l2) + eps, l1 + l2 - eps);
    float cos_a = clampf((l1 * l1 + dist * dist - l2 * l2) / (2.0f * l1 * dist), -1.0f, 1.0f);
    V3 n = vcross(dir, pole);
    if (vlen(n) < 1e-6f) n = any_perp(dir);
    V3 jdir = vrot(dir, vnorm(n), acosf(cos_a));
    return vadd(root, vscale(jdir, l1));
}

static float span_len(const IkChainDef* d, int i, int j) {
    float s = 0.0f;
    for (int k = i; k <= j; k++) s += d->bone[k].len;
    return s;
}

/* solve joints i..j given P[i] = root; fills P[i+1 .. j+1] */
static void chain_positions(const IkChainDef* d, int i, int j, V3 root, V3 target, V3 pole, V3* P) {
    if (i == j) {
        V3 dlt = vsub(target, root);
        V3 dir = vlen(dlt) > 1e-9f ? vnorm(dlt) : v3(1, 0, 0);
        P[j + 1] = vadd(root, vscale(dir, d->bone[i].len));
        return;
    }
    int m = i + (j - i + 1) / 2;
    V3 mid = two_bone_mid(root, target, span_len(d, i, m - 1), span_len(d, m, j), pole);
    chain_positions(d, i, m - 1, root, mid, pole, P);
    chain_positions(d, m, j, mid, target, pole, P);
}

/* project each bone direction onto the swing cone of the previous bone */
static void clamp_swing(const IkChainDef* d, V3* P) {
    int n = d->n;
    V3 dir[IK_MAX_BONES];
    dir[0] = vnorm(vsub(P[1], P[0]));
    if (vlen(dir[0]) < 0.5f) dir[0] = v3(1, 0, 0);
    for (int k = 1; k < n; k++) {
        V3 dk = vnorm(vsub(P[k + 1], P[k]));
        float smax = d->bone[k].swing_max;
        float c = clampf(vdot(dir[k - 1], dk), -1.0f, 1.0f);
        if (c < cosf(smax)) {
            V3 perp = vsub(dk, vscale(dir[k - 1], c));
            if (vlen(perp) > 1e-9f)
                dk = vadd(vscale(dir[k - 1], cosf(smax)), vscale(vnorm(perp), sinf(smax)));
            else
                dk = dir[k - 1]; /* anti-parallel: fold straight */
        }
        dir[k] = dk;
    }
    for (int k = 0; k < n; k++)
        P[k + 1] = vadd(P[k], vscale(dir[k], d->bone[k].len));
}

/* distribute total twist across joints 1..n-1; joint k rotates its sub-chain
 * about the axis through hinge P[k-1] and end P[n]: both stay on the axis,
 * so bone lengths are preserved and the end effector cannot drift */
static void apply_twist(const IkChainDef* d, V3* P, float total, float* applied) {
    int n = d->n;
    applied[0] = 0.0f;
    if (n < 2) return;
    float wsum = 0.0f;
    for (int k = 1; k < n; k++) wsum += d->bone[k].twist_w;
    if (!(wsum > 0.0f)) {
        for (int k = 1; k < n; k++) applied[k] = 0.0f;
        return;
    }
    for (int k = 1; k < n; k++) {
        float tk = clampf(total * d->bone[k].twist_w / wsum,
                          d->bone[k].twist_min, d->bone[k].twist_max);
        V3 axis = vsub(P[n], P[k - 1]);
        if (vlen(axis) < 1e-9f) { applied[k] = 0.0f; continue; }
        applied[k] = tk;
        if (fabsf(tk) < 1e-12f) continue;
        axis = vnorm(axis);
        for (int j = k; j <= n; j++)
            P[j] = rot_about_line(P[j], P[k - 1], axis, tk);
    }
}

int ik_solve(const IkChainDef* d, V3 root, V3 target, V3 pole, float twist_total, IkPose* out) {
    int v = ik_validate(d);
    if (v != IK_OK) return v;
    float maxR = span_len(d, 0, d->n - 1);
    V3 toT = vsub(target, root);
    float dist = vlen(toT);
    V3 dirT = dist > 1e-9f ? vscale(toT, 1.0f / dist) : v3(1, 0, 0);
    float r = clampf(dist, maxR * 1e-6f, maxR * (1.0f - 1e-6f));
    V3 T = vadd(root, vscale(dirT, r));
    out->p[0] = root;
    chain_positions(d, 0, d->n - 1, root, T, pole, out->p);
    clamp_swing(d, out->p);
    apply_twist(d, out->p, twist_total, out->twist);
    clamp_swing(d, out->p);
    return IK_OK;
}

/* numerical reference: cyclic coordinate descent, per-frame joint rotation.
 * records the largest single joint rotation in *max_step */
static float ccd_solve(int n, V3* P, V3 target, int iters, float tol, float* max_step) {
    float err = vlen(vsub(P[n], target));
    for (int it = 0; it < iters && err > tol; it++) {
        for (int j = n - 1; j >= 0; j--) {
            V3 toEnd = vsub(P[n], P[j]);
            V3 toTgt = vsub(target, P[j]);
            if (vlen(toEnd) < 1e-9f || vlen(toTgt) < 1e-9f) continue;
            V3 a = vnorm(toEnd), b = vnorm(toTgt);
            V3 ax = vcross(a, b);
            if (vlen(ax) < 1e-9f) continue;
            float ang = acosf(clampf(vdot(a, b), -1.0f, 1.0f));
            if (ang > 0.5f) ang = 0.5f; /* damped step */
            ax = vnorm(ax);
            for (int k = j + 1; k <= n; k++)
                P[k] = rot_about_line(P[k], P[j], ax, ang);
            if (ang > *max_step) *max_step = ang;
        }
        err = vlen(vsub(P[n], target));
    }
    return err;
}

/* ---------------- self test ---------------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static IkChainDef make_def(int n, const float* lens, float swing, float tmin, float tmax) {
    IkChainDef d;
    d.n = n;
    for (int i = 0; i < n; i++) {
        d.bone[i].len = lens[i];
        d.bone[i].swing_max = swing;
        d.bone[i].twist_min = tmin;
        d.bone[i].twist_max = tmax;
        d.bone[i].twist_w = 1.0f;
    }
    return d;
}

static int swing_violations(const IkChainDef* d, const V3* P) {
    int c = 0;
    for (int k = 1; k < d->n; k++) {
        float ang = dir_angle(vsub(P[k], P[k - 1]), vsub(P[k + 1], P[k]));
        if (ang > d->bone[k].swing_max + 1e-4f) c++;
    }
    return c;
}
static int twist_violations(const IkChainDef* d, const float* tw) {
    int c = 0;
    for (int k = 1; k < d->n; k++)
        if (tw[k] < d->bone[k].twist_min - 1e-5f || tw[k] > d->bone[k].twist_max + 1e-5f) c++;
    return c;
}

int main(void) {
    float L[5] = { 0.9f, 0.7f, 0.5f, 0.4f, 0.3f };
    IkPose pose;

    /* invalid limit params must be rejected */
    IkChainDef def5 = make_def(5, L, 0.6f, -0.8f, 0.8f);
    check("validate_ok", ik_validate(&def5) == IK_OK);
    IkChainDef bad;
    bad = def5; bad.n = 0;                    check("reject_n_zero", ik_validate(&bad) == IK_ERR_N);
    bad = def5; bad.n = IK_MAX_BONES + 1;     check("reject_n_over", ik_validate(&bad) == IK_ERR_N);
    bad = def5; bad.bone[2].len = 0.0f;       check("reject_len_zero", ik_validate(&bad) == IK_ERR_LEN);
    bad = def5; bad.bone[0].len = -1.0f;      check("reject_len_neg", ik_validate(&bad) == IK_ERR_LEN);
    bad = def5; bad.bone[1].swing_max = -0.1f; check("reject_swing_neg", ik_validate(&bad) == IK_ERR_SWING);
    bad = def5; bad.bone[1].swing_max = 3.2f; check("reject_swing_over_pi", ik_validate(&bad) == IK_ERR_SWING);
    bad = def5; bad.bone[3].twist_min = 0.5f; bad.bone[3].twist_max = -0.5f;
    check("reject_twist_inverted", ik_validate(&bad) == IK_ERR_TWIST);
    bad = def5; bad.bone[1].twist_w = -1.0f;  check("reject_weight_neg", ik_validate(&bad) == IK_ERR_WEIGHT);
    bad = def5; for (int i = 0; i < bad.n; i++) bad.bone[i].twist_w = 0.0f;
    check("reject_weight_all_zero", ik_validate(&bad) == IK_ERR_WEIGHT);
    bad = def5; bad.bone[0].len = -1.0f;
    check("solve_rejects_bad_len", ik_solve(&bad, v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0), 0.0f, &pose) == IK_ERR_LEN);
    bad = def5; bad.bone[2].twist_min = 1.0f; bad.bone[2].twist_max = -1.0f;
    check("solve_rejects_bad_twist", ik_solve(&bad, v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0), 0.0f, &pose) == IK_ERR_TWIST);

    /* two-bone triangle by the law of cosines: l1=3, l2=4, dist=4 */
    V3 mid = two_bone_mid(v3(0, 0, 0), v3(4, 0, 0), 3.0f, 4.0f, v3(0, 0, 1));
    check("tb_on_sphere_a", fabsf(vlen(mid) - 3.0f) < 1e-4f);
    check("tb_on_sphere_b", fabsf(vlen(vsub(mid, v3(4, 0, 0))) - 4.0f) < 1e-4f);
    check("tb_cosine_x", fabsf(mid.x - 1.125f) < 1e-3f); /* 3*cos(acos(0.375)) */
    check("tb_pole_side", mid.z > 0.0f);

    /* target beyond reach is projected onto the reachable shell */
    float L2[2] = { 1.0f, 1.0f };
    IkChainDef def2 = make_def(2, L2, IK_PI, -IK_PI, IK_PI);
    check("solve_reach_ok", ik_solve(&def2, v3(0, 0, 0), v3(5, 0, 0), v3(0, 0, 1), 0.0f, &pose) == IK_OK);
    check("reach_clamped", fabsf(vlen(pose.p[2]) - 2.0f) < 1e-3f);
    check("reach_dir", fabsf(pose.p[2].y) < 1e-4f && fabsf(pose.p[2].z) < 1e-4f);

    /* analytic vs numerical iteration, two-bone: CCD started from the
     * analytic pose must barely move; per-joint angle error has an upper bound */
    V3 T2 = v3(1.2f, 0.3f, 0.0f);
    ik_solve(&def2, v3(0, 0, 0), T2, v3(0, 0, 1), 0.0f, &pose);
    V3 q[IK_MAX_BONES + 1];
    memcpy(q, pose.p, sizeof(q));
    float step = 0.0f;
    float err = ccd_solve(2, q, T2, 100, 1e-6f, &step);
    printf("tb_angle_err=%.2e bound=2.0e-02\n", (double)step);
    check("tb_ccd_end", err < 1e-3f);
    check("tb_angle_bound", step < 0.02f);

    /* the numerical reference itself converges from a bent start */
    V3 s[3] = { v3(0, 0, 0), v3(0.7f, 0, 0.5f), v3(1.3f, 0, 0.8f) };
    float step_s = 0.0f;
    check("ccd_converges", ccd_solve(2, s, T2, 200, 1e-6f, &step_s) < 1e-3f);

    /* chain of 5, no active limits: end accuracy + angle error bound vs CCD */
    IkChainDef free5 = make_def(5, L, IK_PI, -IK_PI, IK_PI);
    free5.bone[1].twist_w = 2.0f;
    free5.bone[3].twist_w = 0.5f;
    V3 T5 = v3(1.4f, 0.8f, 0.5f);
    check("solve_chain5_ok", ik_solve(&free5, v3(0, 0, 0), T5, v3(0, 1, 0), 0.0f, &pose) == IK_OK);
    check("chain5_end", vlen(vsub(pose.p[5], T5)) < 1e-3f);
    memcpy(q, pose.p, sizeof(q));
    step = 0.0f;
    err = ccd_solve(5, q, T5, 100, 1e-6f, &step);
    float maxd = 0.0f;
    for (int k = 0; k < 5; k++) {
        float a = dir_angle(vsub(pose.p[k + 1], pose.p[k]), vsub(q[k + 1], q[k]));
        if (a > maxd) maxd = a;
    }
    printf("chain5_angle_err=%.2e bound=2.0e-02\n", (double)(step > maxd ? step : maxd));
    check("chain5_ccd_end", err < 1e-3f);
    check("chain5_angle_bound", step < 0.02f && maxd < 0.02f);

    /* twist distribution must not move the end effector */
    V3 P[IK_MAX_BONES + 1];
    P[0] = v3(0, 0, 0);
    chain_positions(&free5, 0, 4, v3(0, 0, 0), T5, v3(0, 1, 0), P);
    V3 E0 = P[5];
    float tw[IK_MAX_BONES];
    apply_twist(&free5, P, 0.7f, tw);
    check("twist_no_drift", vlen(vsub(P[5], E0)) < 1e-4f);
    float tsum = 0.0f;
    for (int k = 1; k < 5; k++) tsum += tw[k];
    check("twist_sum", fabsf(tsum - 0.7f) < 1e-4f);
    int lok = 1;
    for (int k = 0; k < 5; k++)
        if (fabsf(vlen(vsub(P[k + 1], P[k])) - L[k]) > 1e-4f) lok = 0;
    check("twist_keeps_lengths", lok);

    /* same through the full solver: twist on/off, end stays */
    IkPose p0, p1;
    ik_solve(&free5, v3(0, 0, 0), T5, v3(0, 1, 0), 0.0f, &p0);
    ik_solve(&free5, v3(0, 0, 0), T5, v3(0, 1, 0), 0.7f, &p1);
    check("solve_twist_no_drift", vlen(vsub(p0.p[5], p1.p[5])) < 1e-4f);

    /* tight limits: out-of-range joint count must be zero after solve */
    float L3[3] = { 1.0f, 1.0f, 1.0f };
    IkChainDef lim3 = make_def(3, L3, 0.25f, -0.2f, 0.2f);
    check("solve_lim3_ok", ik_solve(&lim3, v3(0, 0, 0), v3(0.5f, 1.8f, 0.6f), v3(0, 0, 1), 0.9f, &pose) == IK_OK);
    int sv = swing_violations(&lim3, pose.p);
    int tv = twist_violations(&lim3, pose.twist);
    printf("swing_violations=%d twist_violations=%d\n", sv, tv);
    check("swing_violations_zero", sv == 0);
    check("twist_violations_zero", tv == 0);

    /* twist clamp engages: requested 3.0 rad, range +-0.2 */
    IkChainDef tc = make_def(5, L, IK_PI, -0.2f, 0.2f);
    ik_solve(&tc, v3(0, 0, 0), T5, v3(0, 1, 0), 3.0f, &pose);
    check("twist_clamped_zero_violations", twist_violations(&tc, pose.twist) == 0);
    check("twist_clamp_active", fabsf(pose.twist[1] - 0.2f) < 1e-5f);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
