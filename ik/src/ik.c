
/* ik: analytic two-bone IK with pole vector and joint travel limits.
 *
 * Analytic path: the law of cosines gives the opening angle gamma between
 * the first bone and the root->target line; gamma is projected onto the
 * joint travel [min_bend, max_bend] and the pole vector picks the bend
 * side. A golden-section numerical solver is kept as an independent
 * frame-by-frame reference for cross-checking. */
#include <stdio.h>
#include <string.h>
#include <math.h>

typedef struct { float x, y; } V2;

enum {
    IK_OK = 0,
    IK_ERR_ZERO_BONE = 1,     /* bone length is zero or negative */
    IK_ERR_TARGET_AT_ROOT = 2,
    IK_ERR_UNREACHABLE = 3,   /* only raised when limits are disabled */
    IK_ERR_BAD_LIMIT = 4
};

typedef struct {
    int ok;
    int err;
    V2 joint;
    V2 end;
    float a1;   /* absolute angle of first bone */
    float a2;   /* absolute angle of second bone */
    float bend; /* opening angle of first bone vs root->target line, [0, PI] */
} IkResult;

typedef struct {
    int enable_limit;
    float min_bend;   /* joint travel lower bound, radians */
    float max_bend;   /* joint travel upper bound, radians */
} IkLimit;

static const float EPS_DIST = 1e-6f;
static const float PI_F = 3.14159265358979323846f;

static float vlen(V2 v) { return sqrtf(v.x * v.x + v.y * v.y); }
static V2 vsub(V2 a, V2 b) { V2 r = { a.x - b.x, a.y - b.y }; return r; }
static V2 vadd(V2 a, V2 b) { V2 r = { a.x + b.x, a.y + b.y }; return r; }
static V2 vscale(V2 v, float s) { V2 r = { v.x * s, v.y * s }; return r; }
static float vdot(V2 a, V2 b) { return a.x * b.x + a.y * b.y; }
static float vcross(V2 a, V2 b) { return a.x * b.y - a.y * b.x; }
static V2 vnorm(V2 v) {
    float l = vlen(v);
    V2 r = { 0, 0 };
    if (l > 0) { r.x = v.x / l; r.y = v.y / l; }
    return r;
}
static V2 vdir(float ang) { V2 r = { cosf(ang), sinf(ang) }; return r; }
static float clampf(float v, float lo, float hi) { return v < lo ? lo : (v > hi ? hi : v); }

/* pole decides which side of the root->target line the joint sits on */
static float bend_side(V2 dir, V2 pole) {
    return vcross(dir, pole) >= 0.0f ? 1.0f : -1.0f;
}

/* assemble the pose from an opening angle gamma: the angle between the
 * first bone and the root->target line (0 = fully straight, PI = folded
 * back onto itself). The joint lies on the root circle and the end bone
 * always points at (root + dist * dir), i.e. the target ray at the
 * target's distance; for an unreachable target this is the closest point
 * on that ray the chain can aim at. */
static IkResult build_pose(V2 root, V2 dir, float dist, float l1, float l2,
                           float gamma, float side) {
    IkResult r;
    memset(&r, 0, sizeof(r));
    float base = atan2f(dir.y, dir.x);
    r.a1 = base + side * gamma;
    r.joint = vadd(root, vscale(vdir(r.a1), l1));
    V2 aim = vadd(root, vscale(dir, dist));
    V2 d2 = vsub(aim, r.joint);
    r.a2 = atan2f(d2.y, d2.x);
    /* end follows the actual second bone length; when the target is out of
     * reach the gap between end and `aim` is the residual miss distance */
    r.end = vadd(r.joint, vscale(vnorm(d2), l2));
    r.bend = gamma;
    r.ok = 1;
    return r;
}

static int limit_valid(const IkLimit* lim) {
    if (!lim->enable_limit) return 1;
    if (isnan(lim->min_bend) || isinf(lim->min_bend)) return 0;
    if (isnan(lim->max_bend) || isinf(lim->max_bend)) return 0;
    if (lim->min_bend < 0.0f || lim->min_bend > PI_F) return 0;
    if (lim->max_bend < 0.0f || lim->max_bend > PI_F) return 0;
    if (lim->min_bend > lim->max_bend) return 0;
    return 1;
}

/* analytic two-bone solve: cosine law, project gamma into joint travel,
 * pole vector selects the bend side */
static int solve_two_bone(V2 root, V2 target, float l1, float l2,
                          V2 pole, const IkLimit* lim, IkResult* out) {
    memset(out, 0, sizeof(*out));
    if (!(l1 > 0.0f) || !(l2 > 0.0f)) { out->err = IK_ERR_ZERO_BONE; return out->err; }
    if (!limit_valid(lim)) { out->err = IK_ERR_BAD_LIMIT; return out->err; }

    V2 d = vsub(target, root);
    float dist = vlen(d);
    if (dist < EPS_DIST) { out->err = IK_ERR_TARGET_AT_ROOT; return out->err; }

    /* law of cosines on the root triangle:
     * l2^2 = l1^2 + dist^2 - 2*l1*dist*cos(gamma), so gamma is the
     * opening angle between the first bone and the root->target line */
    float cos_g = (l1 * l1 + dist * dist - l2 * l2) / (2.0f * l1 * dist);
    cos_g = clampf(cos_g, -1.0f, 1.0f);
    float gamma = acosf(cos_g);

    if (!lim->enable_limit) {
        /* identical reachability contract to the original solver */
        if (dist > l1 + l2 + 1e-5f) { out->err = IK_ERR_UNREACHABLE; return out->err; }
        if (dist < fabsf(l1 - l2) - 1e-5f) { out->err = IK_ERR_UNREACHABLE; return out->err; }
    } else {
        /* project the desired opening into the joint travel; an out-of-
         * reach target then yields the closest in-travel pose instead of
         * the chain folding straight through itself */
        gamma = clampf(gamma, lim->min_bend, lim->max_bend);
    }

    V2 dir = vnorm(d);
    *out = build_pose(root, dir, dist, l1, l2, gamma, bend_side(dir, pole));
    out->err = IK_OK;
    return IK_OK;
}

/* double-precision geometry used only inside the numeric reference, so
 * its golden-section interval can shrink well past float resolution */
typedef struct { double x, y; } V2d;

static double vlen_d(V2d v) { return sqrt(v.x * v.x + v.y * v.y); }
static V2d vsub_d(V2d a, V2d b) { V2d r = { a.x - b.x, a.y - b.y }; return r; }
static V2d vadd_d(V2d a, V2d b) { V2d r = { a.x + b.x, a.y + b.y }; return r; }
static V2d vscale_d(V2d v, double s) { V2d r = { v.x * s, v.y * s }; return r; }
static V2d vnorm_d(V2d v) {
    double l = vlen_d(v);
    V2d r = { 0.0, 0.0 };
    if (l > 0.0) { r.x = v.x / l; r.y = v.y / l; }
    return r;
}

/* double-precision mirror of build_pose for the reference solver */
static void build_pose_d(V2d root, V2d dir, double dist, double l1, double l2,
                         double gamma, double side, V2d* joint, V2d* end) {
    double base = atan2(dir.y, dir.x);
    double a1 = base + side * gamma;
    V2d u1 = { cos(a1), sin(a1) };
    *joint = vadd_d(root, vscale_d(u1, l1));
    V2d aim = vadd_d(root, vscale_d(dir, dist));
    V2d d2 = vnorm_d(vsub_d(aim, *joint));
    *end = vadd_d(*joint, vscale_d(d2, l2));
}

typedef struct {
    V2d root, dir, target;
    double dist, l1, l2, side;
} NumericCtx;

static double numeric_err(const NumericCtx* cx, double g) {
    V2d j, e;
    build_pose_d(cx->root, cx->dir, cx->dist, cx->l1, cx->l2,
                 g, cx->side, &j, &e);
    return vlen_d(vsub_d(e, cx->target));
}

/* independent reference: golden-section search over the opening angle on
 * a fixed bend side, minimizing end-effector distance frame by frame.
 * The search interval is the joint travel when limits are enabled. */
static IkResult solve_numeric(V2 root, V2 target, float l1, float l2,
                              V2 pole, float gamma_lo, float gamma_hi) {
    const double GR = 0.6180339887498948482;
    NumericCtx cx;
    cx.root = (V2d){ root.x, root.y };
    cx.target = (V2d){ target.x, target.y };
    V2d d = vsub_d(cx.target, cx.root);
    cx.dist = vlen_d(d);
    cx.dir = vnorm_d(d);
    cx.l1 = l1; cx.l2 = l2;
    cx.side = (cx.dir.x * pole.y - cx.dir.y * pole.x) >= 0.0 ? 1.0 : -1.0;

    double a = gamma_lo, b = gamma_hi;
    double c = b - GR * (b - a);
    double e = a + GR * (b - a);
    for (int i = 0; i < 100; i++) {
        if (numeric_err(&cx, c) < numeric_err(&cx, e)) b = e;
        else a = c;
        c = b - GR * (b - a);
        e = a + GR * (b - a);
    }
    /* check the interval ends as well as the interior estimate */
    double best = (a + b) * 0.5;
    double best_err = numeric_err(&cx, best);
    double lo_err = numeric_err(&cx, gamma_lo);
    double hi_err = numeric_err(&cx, gamma_hi);
    if (lo_err < best_err) best = gamma_lo;
    if (hi_err < best_err) best = gamma_hi;

    V2d jd, ed;
    build_pose_d(cx.root, cx.dir, cx.dist, cx.l1, cx.l2, best, cx.side, &jd, &ed);
    IkResult r;
    memset(&r, 0, sizeof(r));
    r.joint.x = (float)jd.x; r.joint.y = (float)jd.y;
    r.end.x = (float)ed.x; r.end.y = (float)ed.y;
    r.bend = (float)best;
    r.a1 = (float)(atan2(cx.dir.y, cx.dir.x) + cx.side * best);
    r.a2 = (float)atan2(ed.y - jd.y, ed.x - jd.x);
    r.ok = 1;
    r.err = IK_OK;
    return r;
}

/* opening angle actually produced by a pose, measured from geometry:
 * angle between the root->target line and the first bone */
static float measured_bend(const IkResult* r, V2 root, V2 target) {
    V2 u = vnorm(vsub(r->joint, root));
    V2 d = vnorm(vsub(target, root));
    return acosf(clampf(vdot(u, d), -1.0f, 1.0f));
}

/* ---- original solver, kept only to prove limit-off poses are unchanged ---- */
static IkResult legacy_solve(V2 root, V2 target, float l1, float l2, V2 pole) {
    IkResult r;
    memset(&r, 0, sizeof(r));
    V2 d = vsub(target, root);
    float dist = vlen(d);
    if (l1 <= 0 || l2 <= 0) return r;
    if (dist > l1 + l2 + 1e-5f) return r;
    if (dist < fabsf(l1 - l2) - 1e-5f) return r;
    if (dist < 1e-6f) return r;
    float cos_a = (l1 * l1 + dist * dist - l2 * l2) / (2 * l1 * dist);
    cos_a = clampf(cos_a, -1.0f, 1.0f);
    float a = acosf(cos_a);
    V2 dir = vnorm(d);
    float side = bend_side(dir, pole);
    float ang = atan2f(dir.y, dir.x) + side * a;
    r.joint.x = root.x + l1 * cosf(ang);
    r.joint.y = root.y + l1 * sinf(ang);
    r.a1 = ang;
    V2 d2 = vsub(target, r.joint);
    r.a2 = atan2f(d2.y, d2.x);
    r.ok = 1;
    return r;
}

/* deterministic LCG so the cross-check is reproducible */
static unsigned long rng_state = 1;
static float randf(float lo, float hi) {
    rng_state = rng_state * 1103515245UL + 12345UL;
    float t = (float)((rng_state >> 8) & 0xffffff) / (float)0xffffff;
    return lo + t * (hi - lo);
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    IkLimit off = { 0, 0.0f, 0.0f };
    V2 root = { 0, 0 }, target = { 1, 0 }, pole = { 0, 1 };

    /* ---- original contract, unchanged ---- */
    IkResult r;
    solve_two_bone(root, target, 1.0f, 1.0f, pole, &off, &r);
    check("reachable", r.err == IK_OK);
    check("joint_on_circle", fabsf(vlen(vsub(r.joint, root)) - 1.0f) < 1e-4f);
    check("end_matches", fabsf(vlen(vsub(r.end, target))) < 1e-4f);
    check("pole_up", r.joint.y > 0.0f);

    solve_two_bone(root, target, 1.0f, 1.0f, (V2){ 0, -1 }, &off, &r);
    check("pole_down", r.joint.y < 0.0f);

    solve_two_bone(root, (V2){ 5, 0 }, 1.0f, 1.0f, pole, &off, &r);
    check("unreachable", r.err == IK_ERR_UNREACHABLE);

    solve_two_bone(root, (V2){ 0.1f, 0 }, 1.0f, 0.2f, pole, &off, &r);
    check("too_close", r.err == IK_ERR_UNREACHABLE);

    solve_two_bone(root, target, 0.0f, 1.0f, pole, &off, &r);
    check("zero_len", r.err == IK_ERR_ZERO_BONE);
    solve_two_bone(root, target, 1.0f, -0.5f, pole, &off, &r);
    check("neg_len", r.err == IK_ERR_ZERO_BONE);

    solve_two_bone(root, root, 1.0f, 1.0f, pole, &off, &r);
    check("target_at_root", r.err == IK_ERR_TARGET_AT_ROOT);

    IkLimit bad1 = { 1, 2.0f, 1.0f };
    IkLimit bad2 = { 1, -0.1f, 1.0f };
    IkLimit bad3 = { 1, 0.0f, 4.0f };
    check("bad_limit_order", solve_two_bone(root, target, 1, 1, pole, &bad1, &r) == IK_ERR_BAD_LIMIT);
    check("bad_limit_min", solve_two_bone(root, target, 1, 1, pole, &bad2, &r) == IK_ERR_BAD_LIMIT);
    check("bad_limit_max", solve_two_bone(root, target, 1, 1, pole, &bad3, &r) == IK_ERR_BAD_LIMIT);

    /* ---- limit-off poses match the original solver exactly ---- */
    int legacy_mismatch = 0;
    rng_state = 1;
    for (int i = 0; i < 200; i++) {
        float l1 = randf(0.3f, 2.0f), l2 = randf(0.3f, 2.0f);
        float lo = fabsf(l1 - l2) + 0.05f, hi = (l1 + l2) * 0.95f;
        float dist = randf(lo, hi);
        float ang = randf(-PI_F, PI_F);
        V2 t = { root.x + dist * cosf(ang), root.y + dist * sinf(ang) };
        V2 p = { randf(-1, 1), randf(-1, 1) };
        IkResult a, b;
        solve_two_bone(root, t, l1, l2, p, &off, &a);
        b = legacy_solve(root, t, l1, l2, p);
        if (vlen(vsub(a.joint, b.joint)) > 1e-5f || fabsf(a.a1 - b.a1) > 1e-5f)
            legacy_mismatch = 1;
    }
    check("limit_off_matches_legacy", !legacy_mismatch);

    /* ---- analytic vs frame-by-frame numeric, with and without limits ---- */
    const struct { float lo, hi; } travels[] = {
        { 0.0f, PI_F },          /* limit enabled, full travel */
        { 0.3f, 2.4f },          /* shoulder-like clamp */
        { 0.8f, 1.8f },          /* narrow travel */
        { 1.2f, 1.2f }           /* fully locked joint */
    };
    float max_joint_err = 0.0f;
    int violations = 0;
    int frames = 0;

    for (int mode = 0; mode < 2; mode++) {
        for (unsigned k = 0; k < sizeof(travels) / sizeof(travels[0]); k++) {
            IkLimit lim = { 1, travels[k].lo, travels[k].hi };
            if (mode == 0) lim.enable_limit = 0;
            rng_state = 7 + k * 101 + mode * 1009;
            for (int i = 0; i < 300; i++) {
                float l1 = randf(0.3f, 2.0f), l2 = randf(0.3f, 2.0f);
                /* mode 0 stays inside the feasible annulus (unreachable
                 * targets are contractually rejected without limits);
                 * mode 1 reaches beyond it so the clamp path is exercised */
                float dist;
                if (mode == 0)
                    dist = randf(fabsf(l1 - l2) + 0.02f, (l1 + l2) * 0.98f);
                else
                    dist = randf(0.05f, (l1 + l2) * 1.3f);
                float ang = randf(-PI_F, PI_F);
                V2 t = { root.x + dist * cosf(ang), root.y + dist * sinf(ang) };
                int side_up = (i & 1);
                /* pole as a point offset to one side of the root->target
                 * line, so the bend side never degenerates when the
                 * target happens to align with a world axis */
                V2 dir_t = vnorm(vsub(t, root));
                V2 perp = { -dir_t.y, dir_t.x };
                V2 p = vadd(t, vscale(perp, side_up ? 1.0f : -1.0f));

                IkResult an;
                int rc = solve_two_bone(root, t, l1, l2, p, &lim, &an);
                if (rc != IK_OK) { check("random_frame_ok", 0); continue; }

                float glo = lim.enable_limit ? lim.min_bend : 0.0f;
                float ghi = lim.enable_limit ? lim.max_bend : PI_F;
                IkResult num = solve_numeric(root, t, l1, l2, p, glo, ghi);

                float je = vlen(vsub(an.joint, num.joint));
                if (je > max_joint_err) max_joint_err = je;
                frames++;

                if (lim.enable_limit) {
                    float beta = measured_bend(&an, root, t);
                    if (beta < lim.min_bend - 1e-4f || beta > lim.max_bend + 1e-4f)
                        violations++;
                    /* pole still decides the side even while clamped:
                     * side=+1 puts the joint left of the root->target
                     * line (cross > 0). Fully straight poses have a
                     * zero cross product by construction and are exempt. */
                    if (an.bend > 1e-3f) {
                        float sgn = vcross(vnorm(vsub(t, root)), vsub(an.joint, root));
                        if ((side_up && sgn < -1e-5f) || (!side_up && sgn > 1e-5f))
                            violations++;
                    }
                }
            }
        }
    }
    printf("frames=%d max_joint_err=%.2e violations=%d\n", frames, max_joint_err, violations);
    /* bound on joint position disagreement with the double-precision
     * numeric reference: float32 geometry is the only remaining error
     * source, so the bound stays an order of magnitude above it */
    check("joint_err_bound", max_joint_err < 1e-4f);
    check("zero_limit_violations", violations == 0);

    /* ---- the reported failure mode: near-root target no longer flips ---- */
    IkResult free_pose, clamped_pose;
    /* arm cannot bend back past 120 degrees at the joint */
    IkLimit elbow = { 1, 0.0f, 2.0944f };
    V2 near_t = { 0.05f, 0.0f };
    solve_two_bone(root, near_t, 1.0f, 1.0f, pole, &off, &free_pose);
    solve_two_bone(root, near_t, 1.0f, 1.0f, pole, &elbow, &clamped_pose);
    printf("free_bend=%.3f clamped_bend=%.3f max_travel=%.3f\n",
           free_pose.bend, clamped_pose.bend, elbow.max_bend);
    /* near-root target wants a near-flat fold-back (opening close to PI/2
     * for equal bones); the free pose is already near the edge and the
     * 120-degree travel does not bind, so both poses must agree */
    check("free_pose_near_foldback", free_pose.bend > 1.5f);
    check("unbound_pose_unchanged", fabsf(clamped_pose.bend - free_pose.bend) < 1e-4f);
    check("clamped_in_travel",
          measured_bend(&clamped_pose, root, near_t) <= elbow.max_bend + 1e-4f);

    /* a reachable target whose opening exceeds the travel: without a
     * limit the arm bends back past it; with the limit it stops at the
     * travel end instead of flipping through. Unequal bones are used on
     * purpose: for equal bones the opening is always PI/2. */
    IkLimit tight = { 1, 0.0f, 1.2f };
    IkResult free_tight, tight_pose;
    V2 close_t = { 0.55f, 0.0f };
    solve_two_bone(root, close_t, 0.5f, 1.0f, pole, &off, &free_tight);
    solve_two_bone(root, close_t, 0.5f, 1.0f, pole, &tight, &tight_pose);
    printf("free_open=%.3f tight_bend=%.3f travel_max=%.3f\n",
           free_tight.bend, tight_pose.bend, tight.max_bend);
    check("free_opening_exceeds_travel", free_tight.bend > tight.max_bend);
    check("tight_clamp_stops_flip", fabsf(tight_pose.bend - tight.max_bend) < 1e-4f);
    check("tight_in_travel", measured_bend(&tight_pose, root, close_t) <= tight.max_bend + 1e-4f);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
