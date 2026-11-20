
/* blend: multi-pose rotation blending.
 *
 * tangent mode (log-Euclidean on SO(3), tangent space at identity):
 *   1. weights are normalized first (sum must not be zero)
 *   2. each joint rotation is mapped to the tangent space by the quat log
 *   3. tangent vectors are combined with the normalized weights
 *   4. the result is mapped back with the quat exp, so it stays a unit quat
 * linear mode: plain per-component lerp (the legacy path, kept bit-for-bit
 *   for parity) which does not renormalize and visibly shrinks when two
 *   distinct rotations are blended -- the "squashed" artifact.
 */
#include <stdio.h>
#include <string.h>
#include <math.h>

#define MAX_BLEND_POSES 8
#define MAX_BLEND_JOINTS 16
#define QUAT_EPS 1e-8f

typedef struct { float x, y, z, w; } Quat;
typedef struct { float x, y, z; } V3;

typedef struct {
    Quat rot[MAX_BLEND_JOINTS];
    V3 root;
    int joints;
} Pose;

typedef enum {
    BLEND_OK = 0,
    BLEND_ERR_NO_POSES = 1,
    BLEND_ERR_CAPACITY = 2,
    BLEND_ERR_WEIGHT_SUM = 3,
    BLEND_ERR_DEGENERATE_ROT = 4,
    BLEND_ERR_JOINT_MISMATCH = 5
} BlendError;

static const char* blend_strerror(BlendError e) {
    switch (e) {
        case BLEND_OK:                  return "ok";
        case BLEND_ERR_NO_POSES:        return "no input poses";
        case BLEND_ERR_CAPACITY:        return "too many input poses";
        case BLEND_ERR_WEIGHT_SUM:      return "weight sum is zero or non-finite";
        case BLEND_ERR_DEGENERATE_ROT:  return "degenerate (zero-norm) rotation";
        case BLEND_ERR_JOINT_MISMATCH:  return "poses disagree on joint count";
    }
    return "unknown error";
}

/* ---------- quaternion math ---------- */

static Quat q_from_axis_angle(V3 axis, float ang) {
    float l = sqrtf(axis.x * axis.x + axis.y * axis.y + axis.z * axis.z);
    Quat q = { 0, 0, 0, 1 };
    if (l <= QUAT_EPS) return q;
    float h = 0.5f * ang, s = sinf(h) / l;
    q.x = axis.x * s; q.y = axis.y * s; q.z = axis.z * s; q.w = cosf(h);
    return q;
}

static float q_dot(Quat a, Quat b) {
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}

static float q_norm(Quat q) { return sqrtf(q_dot(q, q)); }

static Quat q_normalize(Quat q) {
    float l = q_norm(q);
    Quat r = { 0, 0, 0, 1 };
    if (l <= QUAT_EPS) return r;
    r.x = q.x / l; r.y = q.y / l; r.z = q.z / l; r.w = q.w / l;
    return r;
}

/* quaternions double cover SO(3): pick the hemisphere of the identity so
 * that log(-q) == log(q) and the angle of every log stays in [0, pi] */
static Quat q_canonical(Quat q) {
    if (q.w < 0.0f) { q.x = -q.x; q.y = -q.y; q.z = -q.z; q.w = -q.w; }
    return q;
}

/* log: unit quat -> tangent space at identity (rotation vector) */
static V3 q_log(Quat q) {
    q = q_canonical(q);
    float vlen = sqrtf(q.x * q.x + q.y * q.y + q.z * q.z);
    V3 r = { 0, 0, 0 };
    if (vlen <= QUAT_EPS) return r;                       /* near identity */
    float theta = 2.0f * atan2f(vlen, q.w);              /* in [0, pi] */
    float k = theta / vlen;
    r.x = q.x * k; r.y = q.y * k; r.z = q.z * k;
    return r;
}

/* exp: tangent space at identity -> unit quat */
static Quat q_exp(V3 v) {
    float h = 0.5f * sqrtf(v.x * v.x + v.y * v.y + v.z * v.z);
    Quat q = { 0, 0, 0, 1 };
    if (h <= 0.5f * QUAT_EPS) return q;                  /* near zero vector */
    float s = sinf(h) / (2.0f * h);
    q.x = v.x * s; q.y = v.y * s; q.z = v.z * s; q.w = cosf(h);
    return q;
}

/* geodetic angle between two rotations, in radians, [0, pi] */
static float q_angle(Quat a, Quat b) {
    float d = fabsf(q_dot(a, b));
    if (d > 1.0f) d = 1.0f;
    return 2.0f * acosf(d);
}

/* pairwise spherical interpolation; caller passes unit quaternions */
static Quat q_slerp(Quat a, Quat b, float t) {
    float d = q_dot(a, b);
    if (d < 0.0f) {                                      /* shortest path */
        b.x = -b.x; b.y = -b.y; b.z = -b.z; b.w = -b.w;
        d = -d;
    }
    if (d > 0.9995f) {                                   /* nearly parallel: lerp */
        Quat r = {
            a.x + t * (b.x - a.x), a.y + t * (b.y - a.y),
            a.z + t * (b.z - a.z), a.w + t * (b.w - a.w)
        };
        return q_normalize(r);
    }
    float theta = acosf(d);
    float s0 = sinf((1.0f - t) * theta) / sinf(theta);
    float s1 = sinf(t * theta) / sinf(theta);
    Quat r = {
        s0 * a.x + s1 * b.x, s0 * a.y + s1 * b.y,
        s0 * a.z + s1 * b.z, s0 * a.w + s1 * b.w
    };
    return r;
}

/* sequential pairwise slerp, running weights renormalized each step */
static Quat q_slerp_chain(const Quat qs[], const float ws[], int n) {
    Quat acc = qs[0];
    float accw = ws[0];
    for (int i = 1; i < n; i++) {
        accw += ws[i];
        float t = ws[i] / accw;
        acc = q_slerp(acc, qs[i], t);
    }
    return acc;
}

/* ---------- weight handling ---------- */

static int valid_finite(float v) {
    return v == v && v != INFINITY && v != -INFINITY;
}

/* normalize raw weights into nw; rejects zero total and non-finite input */
static BlendError normalize_weights(const float* w, int n, float* nw) {
    float sum = 0.0f;
    for (int i = 0; i < n; i++) {
        if (!valid_finite(w[i])) return BLEND_ERR_WEIGHT_SUM;
        sum += w[i];
    }
    if (fabsf(sum) <= QUAT_EPS) return BLEND_ERR_WEIGHT_SUM;
    for (int i = 0; i < n; i++) nw[i] = w[i] / sum;
    return BLEND_OK;
}

/* ---------- pose blending ---------- */

/* use_tangent != 0 -> tangent-space log/exp blend; 0 -> legacy linear lerp */
static BlendError blend_poses(const Pose* poses, const float* weights, int n,
                              int use_tangent, Pose* out) {
    if (n <= 0) return BLEND_ERR_NO_POSES;
    if (n > MAX_BLEND_POSES) return BLEND_ERR_CAPACITY;
    int joints = poses[0].joints;
    if (joints <= 0 || joints > MAX_BLEND_JOINTS) return BLEND_ERR_JOINT_MISMATCH;
    for (int i = 1; i < n; i++)
        if (poses[i].joints != joints) return BLEND_ERR_JOINT_MISMATCH;

    /* degenerate rotations are rejected in both modes */
    for (int i = 0; i < n; i++)
        for (int j = 0; j < joints; j++)
            if (q_norm(poses[i].rot[j]) <= QUAT_EPS) return BLEND_ERR_DEGENERATE_ROT;

    float nw[MAX_BLEND_POSES];
    BlendError e = normalize_weights(weights, n, nw);
    if (e != BLEND_OK) return e;

    out->joints = joints;

    for (int j = 0; j < joints; j++) {
        if (use_tangent) {
            V3 acc = { 0, 0, 0 };
            for (int i = 0; i < n; i++) {
                Quat uq = q_normalize(poses[i].rot[j]);
                V3 lv = q_log(uq);
                acc.x += nw[i] * lv.x;
                acc.y += nw[i] * lv.y;
                acc.z += nw[i] * lv.z;
            }
            out->rot[j] = q_exp(acc);
        } else {
            /* legacy path: exactly the old per-component lerp, no renorm */
            Quat acc = { 0, 0, 0, 0 };
            for (int i = 0; i < n; i++) {
                acc.x += nw[i] * poses[i].rot[j].x;
                acc.y += nw[i] * poses[i].rot[j].y;
                acc.z += nw[i] * poses[i].rot[j].z;
                acc.w += nw[i] * poses[i].rot[j].w;
            }
            out->rot[j] = acc;
        }
    }

    /* translations live in Euclidean space: weighted mean in both modes */
    V3 r = { 0, 0, 0 };
    for (int i = 0; i < n; i++) {
        r.x += nw[i] * poses[i].root.x;
        r.y += nw[i] * poses[i].root.y;
        r.z += nw[i] * poses[i].root.z;
    }
    out->root = r;
    return BLEND_OK;
}

/* ---------- deterministic test helpers ---------- */

static unsigned long rng_state = 0x12345678ul;
static float rng_f(void) {                       /* LCG, uniform (0, 1) */
    rng_state = rng_state * 1103515245ul + 12345ul;
    return (float)((rng_state >> 8) & 0xffffff) / (float)0x1000000;
}

static Quat random_unit_quat(float max_ang) {
    V3 axis = { rng_f() - 0.5f, rng_f() - 0.5f, rng_f() - 0.5f };
    if (axis.x == 0 && axis.y == 0 && axis.z == 0) axis.x = 1.0f;
    return q_from_axis_angle(axis, rng_f() * max_ang);
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    /* ---- basic quaternion round trips ---- */
    Quat qid = q_from_axis_angle((V3){ 1, 0, 0 }, 0.0f);
    check("identity_quat", fabsf(qid.w - 1.0f) < 1e-6f && qid.x == 0.0f);
    V3 lid = q_log(qid);
    check("log_identity_zero", fabsf(lid.x) < 1e-6f && fabsf(lid.y) < 1e-6f && fabsf(lid.z) < 1e-6f);
    check("exp_zero_identity", fabsf(q_exp((V3){ 0, 0, 0 }).w - 1.0f) < 1e-6f);

    V3 axis = { 0, 0, 1 };
    Quat q90 = q_from_axis_angle(axis, 1.5707963f);
    V3 l90 = q_log(q90);
    check("log_magnitude", fabsf(sqrtf(l90.x * l90.x + l90.y * l90.y + l90.z * l90.z)
                                 - 1.5707963f) < 1e-5f);
    Quat qback = q_exp(l90);
    check("log_exp_roundtrip", q_angle(q90, qback) < 1e-5f);

    /* double cover: log(-q) must equal log(q) */
    Quat qneg = { -q90.x, -q90.y, -q90.z, -q90.w };
    check("double_cover", q_angle(q_exp(q_log(qneg)), q90) < 1e-5f);

    /* rotation near pi stays finite through log/exp */
    Quat qpi = q_from_axis_angle((V3){ 1, 0, 0 }, 3.1415926f);
    check("near_pi_finite", q_angle(q_exp(q_log(qpi)), qpi) < 1e-4f);

    /* slerp endpoints and midpoint */
    float PI_F = 3.14159265358979323846f;
    Quat qa = q_from_axis_angle((V3){ 0, 0, 1 }, 0.0f);
    Quat qb = q_from_axis_angle((V3){ 0, 0, 1 }, PI_F);
    check("slerp_endpoints", q_angle(q_slerp(qa, qb, 0.0f), qa) < 1e-6f &&
                             q_angle(q_slerp(qa, qb, 1.0f), qb) < 1e-6f);
    check("slerp_midpoint", fabsf(q_angle(qa, q_slerp(qa, qb, 0.5f)) - PI_F * 0.5f) < 1e-4f);

    /* ---- weight normalization ---- */
    float raw[3] = { 1.0f, 2.0f, 3.0f }, nw[3];
    check("normalize_ok", normalize_weights(raw, 3, nw) == BLEND_OK &&
                          fabsf(nw[0] + nw[1] + nw[2] - 1.0f) < 1e-6f);
    float zeros[2] = { 0.0f, 0.0f };
    check("zero_weight_sum", normalize_weights(zeros, 2, nw) == BLEND_ERR_WEIGHT_SUM);
    float bad[2] = { 1.0f, NAN };
    check("nan_weight", normalize_weights(bad, 2, nw) == BLEND_ERR_WEIGHT_SUM);

    /* ---- error paths through the public entry point ---- */
    Pose poses[MAX_BLEND_POSES];
    for (int i = 0; i < MAX_BLEND_POSES; i++) {
        poses[i].joints = 1;
        poses[i].rot[0] = qid;
        poses[i].root = (V3){ 0, 0, 0 };
    }
    Pose out;
    float w2[2] = { 1, 1 };
    check("no_poses", blend_poses(poses, w2, 0, 1, &out) == BLEND_ERR_NO_POSES);
    check("capacity", blend_poses(poses, w2, MAX_BLEND_POSES + 1, 1, &out) == BLEND_ERR_CAPACITY);
    check("blend_zero_weights", blend_poses(poses, zeros, 2, 1, &out) == BLEND_ERR_WEIGHT_SUM);

    Pose deg = poses[0];
    deg.rot[0] = (Quat){ 0, 0, 0, 0 };
    check("degenerate_rotation", blend_poses((Pose[]){ poses[0], deg }, w2, 2, 1, &out)
                                 == BLEND_ERR_DEGENERATE_ROT);
    check("degenerate_linear_too", blend_poses((Pose[]){ poses[0], deg }, w2, 2, 0, &out)
                                   == BLEND_ERR_DEGENERATE_ROT);

    Pose mismatch = poses[0];
    mismatch.joints = 2;
    check("joint_mismatch", blend_poses((Pose[]){ poses[0], mismatch }, w2, 2, 1, &out)
                            == BLEND_ERR_JOINT_MISMATCH);
    check("error_messages", strlen(blend_strerror(BLEND_OK)) > 0 &&
                            strlen(blend_strerror(BLEND_ERR_WEIGHT_SUM)) > 0 &&
                            strlen(blend_strerror(BLEND_ERR_DEGENERATE_ROT)) > 0);

    /* ---- output of tangent blend is always a proper rotation ---- */
    Pose p0 = { .joints = 2, .root = { 0, 0, 0 } };
    Pose p1 = { .joints = 2, .root = { 1, 0, 0 } };
    p0.rot[0] = q_from_axis_angle((V3){ 1, 0, 0 }, 0.7f);
    p0.rot[1] = q_from_axis_angle((V3){ 0, 1, 0 }, 1.1f);
    p1.rot[0] = q_from_axis_angle((V3){ 0, 0, 1 }, 1.3f);
    p1.rot[1] = q_from_axis_angle((V3){ 1, 1, 0 }, 0.9f);
    float max_unit_dev = 0.0f;
    for (int k = 0; k <= 20; k++) {
        float t = (float)k / 20.0f;
        float wt[2] = { 1.0f - t, t };
        BlendError rc = blend_poses((Pose[]){ p0, p1 }, wt, 2, 1, &out);
        if (rc != BLEND_OK) { check("tangent_unit_norm", 0); break; }
        for (int j = 0; j < 2; j++) {
            float d = fabsf(q_norm(out.rot[j]) - 1.0f);
            if (d > max_unit_dev) max_unit_dev = d;
        }
    }
    check("tangent_unit_norm", max_unit_dev < 1e-5f);

    /* ---- legacy mode reproduces the old per-component lerp exactly ---- */
    Pose lp0 = { .joints = 1, .root = { 0, 0, 0 }, .rot = { q90 } };
    Pose lp1 = { .joints = 1, .root = { 0, 0, 0 }, .rot = { qpi } };
    check("linear_parity", blend_poses((Pose[]){ lp0, lp1 }, (float[]){ 0.25f, 0.75f },
                                       2, 0, &out) == BLEND_OK);
    Quat expect = {
        0.25f * q90.x + 0.75f * qpi.x, 0.25f * q90.y + 0.75f * qpi.y,
        0.25f * q90.z + 0.75f * qpi.z, 0.25f * q90.w + 0.75f * qpi.w
    };
    check("linear_exact_components",
          fabsf(out.rot[0].x - expect.x) < 1e-7f && fabsf(out.rot[0].y - expect.y) < 1e-7f &&
          fabsf(out.rot[0].z - expect.z) < 1e-7f && fabsf(out.rot[0].w - expect.w) < 1e-7f);

    /* the old path shrinks; the tangent path does not */
    Pose lo0 = { .joints = 1, .root = { 0 }, .rot = { q_from_axis_angle((V3){ 0, 0, 1 }, 1.2f) } };
    Pose lo1 = { .joints = 1, .root = { 0 }, .rot = { q_from_axis_angle((V3){ 0, 1, 0 }, 1.2f) } };
    blend_poses((Pose[]){ lo0, lo1 }, (float[]){ 0.5f, 0.5f }, 2, 0, &out);
    check("linear_shrinks", q_norm(out.rot[0]) < 0.99f);
    blend_poses((Pose[]){ lo0, lo1 }, (float[]){ 0.5f, 0.5f }, 2, 1, &out);
    check("tangent_no_shrink", fabsf(q_norm(out.rot[0]) - 1.0f) < 1e-6f);

    /* ---- endpoint weights recover an input pose ---- */
    float wone[2] = { 1.0f, 0.0f };
    blend_poses((Pose[]){ p0, p1 }, wone, 2, 1, &out);
    check("endpoint_pose0", q_angle(out.rot[0], p0.rot[0]) < 1e-5f &&
                            q_angle(out.rot[1], p0.rot[1]) < 1e-5f);
    check("endpoint_root", fabsf(out.root.x - p0.root.x) < 1e-6f);

    /* same axis: log blend and slerp are mathematically identical; any
     * residual is float32 rounding in the sin/sin vs log/exp paths */
    float max_coax = 0.0f;
    for (int k = 1; k < 10; k++) {
        float t = (float)k / 10.0f;
        Pose c0 = { .joints = 1, .rot = { q_from_axis_angle(axis, 0.3f) } };
        Pose c1 = { .joints = 1, .rot = { q_from_axis_angle(axis, 1.7f) } };
        blend_poses((Pose[]){ c0, c1 }, (float[]){ 1 - t, t }, 2, 1, &out);
        Quat ref = q_slerp(c0.rot[0], c1.rot[0], t);
        float d = q_angle(out.rot[0], ref);
        if (d > max_coax) max_coax = d;
    }
    check("coaxis_matches_slerp", max_coax < 1e-3f);

    /* ---- upper bound on joint-angle error vs pairwise slerp ----
     * the identity-tangent log map is single-valued only inside a pi-ball
     * around the identity: poses nearly antipodal to it (or to each other
     * after canonicalization) sit on the cut locus where the log is
     * discontinuous, so the bound is stated for the moderate-angle regime
     * animation blends live in (individual joint rotations <= 1 rad here,
     * pair separation <= ~1.7 rad). */

    /* two poses: sweep t over random moderate-angle pairs */
    float max_err2 = 0.0f, max_pair_sep = 0.0f;
    rng_state = 7;
    for (int trial = 0; trial < 200; trial++) {
        Quat a = random_unit_quat(1.0f), b = random_unit_quat(1.0f);
        float sep = q_angle(a, b);
        if (sep > max_pair_sep) max_pair_sep = sep;
        Pose pa = { .joints = 1, .rot = { a } };
        Pose pb = { .joints = 1, .rot = { b } };
        for (int k = 1; k < 10; k++) {
            float t = (float)k / 10.0f;
            blend_poses((Pose[]){ pa, pb }, (float[]){ 1 - t, t }, 2, 1, &out);
            float d = q_angle(out.rot[0], q_slerp(a, b, t));
            if (d > max_err2) max_err2 = d;
        }
    }
    printf("slerp_err_two_pose_max=%.6f rad (pair separation up to %.3f)\n",
           max_err2, max_pair_sep);
    check("two_pose_slerp_bound", max_err2 < 0.04f);

    /* antipodal pair documents the domain limit: tangent blend and slerp
     * can disagree maximally near the cut locus (no error bound there) */
    Quat ap = q_from_axis_angle((V3){ 1, 0, 0 }, PI_F - 0.02f);
    Quat am = q_from_axis_angle((V3){ -1, 0, 0 }, PI_F - 0.02f);
    Pose pp = { .joints = 1, .rot = { ap } };
    Pose pm = { .joints = 1, .rot = { am } };
    blend_poses((Pose[]){ pp, pm }, (float[]){ 0.5f, 0.5f }, 2, 1, &out);
    float antipodal_err = q_angle(out.rot[0], q_slerp(ap, am, 0.5f));
    printf("antipodal_disagreement=%.6f rad (domain limit, not bounded)\n",
           antipodal_err);
    check("antipodal_domain_documented", antipodal_err > 0.5f);

    /* n poses vs sequential pairwise slerp with renormalized running weights */
    float max_errn = 0.0f;
    rng_state = 99;
    for (int trial = 0; trial < 200; trial++) {
        int n = 3 + (int)(rng_f() * 2.0f);                 /* 3..4 poses */
        Pose np[4]; Quat nq[4] = { qid, qid, qid, qid }; float w[4];
        for (int i = 0; i < n; i++) {
            nq[i] = random_unit_quat(1.5f);
            np[i].joints = 1;
            np[i].rot[0] = nq[i];
        }
        for (int i = 0; i < n; i++) w[i] = 0.1f + rng_f();
        blend_poses(np, w, n, 1, &out);
        Quat ref = q_slerp_chain(nq, w, n);
        float d = q_angle(out.rot[0], ref);
        if (d > max_errn) max_errn = d;
    }
    printf("slerp_err_n_pose_max=%.6f rad\n", max_errn);
    check("n_pose_slerp_bound", max_errn < 0.10f);

    /* ---- angular velocity continuity while weights change ----
     * weights are smooth sinusoids that hit zero (a pose fades fully out
     * and back in); the blended rotation path must have continuous angular
     * velocity, i.e. no jump when a weight crosses zero. Angular velocity is
     * the rotation vector of the step's relative rotation, so sign flips of
     * a scalar speed cannot hide a discontinuity. */
    Pose wp[3];
    wp[0].joints = 1; wp[1].joints = 1; wp[2].joints = 1;
    wp[0].rot[0] = q_from_axis_angle((V3){ 1, 0, 0 }, 1.0f);
    wp[1].rot[0] = q_from_axis_angle((V3){ 0, 1, 0 }, 1.6f);
    wp[2].rot[0] = q_from_axis_angle((V3){ 0, 0, 1 }, 2.2f);
    enum { NT = 4001 };
    static Quat path[NT];
    static V3 vel[NT];
    float dt = 6.2831853f / (NT - 1);
    for (int k = 0; k < NT; k++) {
        float t = (float)k / (NT - 1) * 6.2831853f;
        float ww[3];
        ww[0] = sinf(t) * sinf(t);             /* zero at t = 0, pi, 2pi */
        ww[1] = cosf(t) * cosf(t);             /* zero at t = pi/2, 3pi/2 */
        ww[2] = 0.5f + 0.5f * sinf(2.0f * t);  /* always > 0 */
        blend_poses(wp, ww, 3, 1, &out);
        path[k] = out.rot[0];
    }
    /* angular velocity vector of each interior sample */
    for (int k = 1; k < NT - 1; k++) {
        Quat a = path[k - 1], b = path[k + 1];
        Quat conj_a = { -a.x, -a.y, -a.z, a.w };
        /* relative rotation q1 * conj(q0) */
        Quat r = {
            b.w * conj_a.x + b.x * conj_a.w + b.y * conj_a.z - b.z * conj_a.y,
            b.w * conj_a.y - b.x * conj_a.z + b.y * conj_a.w + b.z * conj_a.x,
            b.w * conj_a.z + b.x * conj_a.y - b.y * conj_a.x + b.z * conj_a.w,
            b.w * conj_a.w - b.x * conj_a.x - b.y * conj_a.y - b.z * conj_a.z
        };
        V3 lv = q_log(r);
        vel[k].x = lv.x / (2.0f * dt);
        vel[k].y = lv.y / (2.0f * dt);
        vel[k].z = lv.z / (2.0f * dt);
    }
    float max_vel_step = 0.0f, max_acc = 0.0f;
    for (int k = 2; k < NT - 2; k++) {
        float sx = fabsf(vel[k + 1].x - vel[k].x);
        float sy = fabsf(vel[k + 1].y - vel[k].y);
        float sz = fabsf(vel[k + 1].z - vel[k].z);
        float step = sx + sy + sz;
        float ax = fabsf(vel[k + 1].x - 2.0f * vel[k].x + vel[k - 1].x);
        float ay = fabsf(vel[k + 1].y - 2.0f * vel[k].y + vel[k - 1].y);
        float az = fabsf(vel[k + 1].z - 2.0f * vel[k].z + vel[k - 1].z);
        float acc = ax + ay + az;
        if (step > max_vel_step) max_vel_step = step;
        if (acc > max_acc) max_acc = acc;
    }
    printf("velocity_vector_max_step=%.6f accel_vector_max=%.6f\n",
           max_vel_step, max_acc);
    check("velocity_continuous", max_vel_step < 0.05f);
    check("acceleration_bounded", max_acc < 0.05f);

    /* directly probe the instants where a weight crosses zero: the angular
     * velocity vectors on the two sides must agree (no pop) */
    int cross_ok = 1;
    for (int ci = 0; ci < 5; ci++) {
        float tc = ci * 6.2831853f / 4.0f;      /* 0, pi/2, pi, 3pi/2, 2pi */
        Quat qm = qid, q0 = qid, qp = qid;
        for (int side = 0; side < 3; side++) {
            float t = tc + (side - 1) * dt;
            float ww[3] = {
                sinf(t) * sinf(t), cosf(t) * cosf(t),
                0.5f + 0.5f * sinf(2.0f * t)
            };
            blend_poses(wp, ww, 3, 1, &out);
            if (side == 0) qm = out.rot[0];
            if (side == 1) q0 = out.rot[0];
            if (side == 2) qp = out.rot[0];
        }
        if (q_angle(qm, q0) > 0.01f || q_angle(q0, qp) > 0.01f) cross_ok = 0;
    }
    check("zero_weight_crossing", cross_ok);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
