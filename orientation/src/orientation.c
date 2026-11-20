
/* orientation: ground-normal alignment with lagged normal filtering,
   two-stage shortest-arc rotation (ground fit + in-plane yaw), angular
   speed cap and exponential damping. */
#include <stdio.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---------------- vectors ---------------- */

typedef struct { float x, y, z; } V3;

static V3 v3(float x, float y, float z) { V3 r = { x, y, z }; return r; }
static V3 vadd(V3 a, V3 b) { return v3(a.x + b.x, a.y + b.y, a.z + b.z); }
static V3 vsub(V3 a, V3 b) { return v3(a.x - b.x, a.y - b.y, a.z - b.z); }
static V3 vmul(V3 a, float s) { return v3(a.x * s, a.y * s, a.z * s); }
static float vdot(V3 a, V3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static V3 vcross(V3 a, V3 b) {
    return v3(a.y * b.z - a.z * b.y,
              a.z * b.x - a.x * b.z,
              a.x * b.y - a.y * b.x);
}
static float vlen(V3 a) { return sqrtf(vdot(a, a)); }
static V3 vnorm(V3 a) {
    float l = vlen(a);
    return l > 0.0f ? vmul(a, 1.0f / l) : v3(0, 0, 0);
}
static V3 vlerp(V3 a, V3 b, float t) { return vadd(vmul(a, 1.0f - t), vmul(b, t)); }

static float clampf(float x, float lo, float hi) {
    return x < lo ? lo : (x > hi ? hi : x);
}

/* angle between two unit vectors, always in [0, pi] */
static float vangle(V3 a, V3 b) {
    return acosf(clampf(vdot(a, b), -1.0f, 1.0f));
}

/* ---------------- quaternions (Hamilton; local +y up, -z forward) ---------------- */

typedef struct { float x, y, z, w; } Quat;

static Quat qid(void) { Quat q = { 0, 0, 0, 1 }; return q; }

static Quat qnorm(Quat q) {
    float l = sqrtf(q.x * q.x + q.y * q.y + q.z * q.z + q.w * q.w);
    if (l <= 0.0f) return qid();
    float s = 1.0f / l;
    Quat r = { q.x * s, q.y * s, q.z * s, q.w * s };
    return r;
}

static Quat qmul(Quat a, Quat b) {
    Quat r;
    r.w = a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z;
    r.x = a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y;
    r.y = a.w * b.y - a.x * b.z + a.y * b.w + a.z * b.x;
    r.z = a.w * b.z + a.x * b.y - a.y * b.x + a.z * b.w;
    return r;
}

static Quat qconj(Quat q) { Quat r = { -q.x, -q.y, -q.z, q.w }; return r; }

static V3 qrot(Quat q, V3 v) {
    Quat vq = { v.x, v.y, v.z, 0 };
    Quat r = qmul(qmul(q, vq), qconj(q));
    return v3(r.x, r.y, r.z);
}

/* rotation of angle a about unit axis n */
static Quat qaxis(V3 n, float a) {
    float h = a * 0.5f;
    float s = sinf(h), c = cosf(h);
    Quat r = { n.x * s, n.y * s, n.z * s, c };
    return r;
}

/* rotation matrix, columns are the rotated basis vectors */
typedef struct { V3 c0, c1, c2; } Mat3;

static Mat3 qtomat(Quat q) {
    float xx = q.x * q.x, yy = q.y * q.y, zz = q.z * q.z;
    float xy = q.x * q.y, xz = q.x * q.z, yz = q.y * q.z;
    float wx = q.w * q.x, wy = q.w * q.y, wz = q.w * q.z;
    Mat3 m;
    m.c0 = v3(1 - 2 * (yy + zz), 2 * (xy + wz),     2 * (xz - wy));
    m.c1 = v3(2 * (xy - wz),     1 - 2 * (xx + zz), 2 * (yz + wx));
    m.c2 = v3(2 * (xz + wy),     2 * (yz - wx),     1 - 2 * (xx + yy));
    return m;
}

static float mdet(Mat3 m) { return vdot(m.c0, vcross(m.c1, m.c2)); }

static Quat mattoquat(Mat3 m) {
    float trace = m.c0.x + m.c1.y + m.c2.z;
    Quat q;
    if (trace > 0.0f) {
        float s = sqrtf(trace + 1.0f) * 2.0f;
        q.w = 0.25f * s;
        q.x = (m.c2.y - m.c1.z) / s;
        q.y = (m.c0.z - m.c2.x) / s;
        q.z = (m.c1.x - m.c0.y) / s;
    } else if (m.c0.x > m.c1.y && m.c0.x > m.c2.z) {
        float s = sqrtf(1.0f + m.c0.x - m.c1.y - m.c2.z) * 2.0f;
        q.w = (m.c2.y - m.c1.z) / s;
        q.x = 0.25f * s;
        q.y = (m.c0.y + m.c1.x) / s;
        q.z = (m.c0.z + m.c2.x) / s;
    } else if (m.c1.y > m.c2.z) {
        float s = sqrtf(1.0f + m.c1.y - m.c0.x - m.c2.z) * 2.0f;
        q.w = (m.c0.z - m.c2.x) / s;
        q.x = (m.c0.y + m.c1.x) / s;
        q.y = 0.25f * s;
        q.z = (m.c1.z + m.c2.y) / s;
    } else {
        float s = sqrtf(1.0f + m.c2.z - m.c0.x - m.c1.y) * 2.0f;
        q.w = (m.c1.x - m.c0.y) / s;
        q.x = (m.c0.z + m.c2.x) / s;
        q.y = (m.c1.z + m.c2.y) / s;
        q.z = 0.25f * s;
    }
    return qnorm(q);
}

/* shortest-arc rotation from unit a to unit b (angle never exceeds pi) */
static Quat qbetween(V3 a, V3 b) {
    float d = vdot(a, b);
    if (d <= -1.0f + 1e-6f) {
        /* antipodal: any perpendicular axis, 180 degrees */
        V3 helper = fabsf(a.x) < 0.9f ? v3(1, 0, 0) : v3(0, 1, 0);
        V3 axis = vnorm(vcross(a, helper));
        return qaxis(axis, (float)M_PI);
    }
    V3 axis = vcross(a, b);
    float s = sqrtf((1.0f + d) * 2.0f);
    Quat q = { axis.x / s, axis.y / s, axis.z / s, s * 0.5f };
    return qnorm(q);
}

/* ---------------- normal filter: EMA with lag, rejects zero normals ---------------- */

typedef struct {
    V3 n;        /* filtered unit normal */
    int valid;   /* has been initialized with a valid sample */
    float tau;   /* time constant (s); larger = smoother but more lag */
} NormalFilter;

static void nf_init(NormalFilter* f, float tau) {
    memset(f, 0, sizeof(*f));
    f->tau = tau;
}

/* returns 0 and keeps the previous estimate when the raw normal is ~zero */
static int nf_push(NormalFilter* f, V3 raw, float dt) {
    float l = vlen(raw);
    if (l < 1e-6f) return 0;                       /* reject: no ground */
    V3 n = vmul(raw, 1.0f / l);
    if (!f->valid) {
        f->n = n;
        f->valid = 1;
        return 1;
    }
    float a = dt / (f->tau + dt);                  /* EMA weight */
    f->n = vnorm(vlerp(f->n, n, a));
    return 1;
}

/* ---------------- target orientation (analytic reference) ---------------- */

/*
 * Desired orientation built analytically from a unit ground normal and a
 * world-space heading hint. Local basis: +y up, -z forward, right-handed.
 * heading is projected onto the ground plane; this matrix is the ground
 * truth used to validate the incremental two-stage path.
 */
static Quat orient_target(V3 normal, V3 heading) {
    V3 up = vnorm(normal);
    V3 f = vsub(heading, vmul(up, vdot(heading, up)));
    if (vlen(f) < 1e-6f) f = v3(0, 0, 1);          /* heading parallel to up */
    f = vnorm(f);
    V3 right = vnorm(vcross(f, up));               /* right = f x up */
    Mat3 m = { right, up, vmul(f, -1.0f) };        /* columns: right, up, -f */
    return mattoquat(m);
}

/* ---------------- alignment step: shortest arc + speed cap + damping ---------------- */

typedef struct {
    float wmax;     /* angular speed cap (rad/s) */
    float tau;      /* damping time constant (s); 0 disables damping */
    float eps;      /* convergence threshold (rad) */
} AlignCfg;

/*
 * Rotate cur toward tgt by one fixed step. The relative rotation is always
 * the shortest arc (|angle| <= pi). The step angle is the smaller of the
 * damped first-order step and the hard cap wmax*dt, and never overshoots.
 */
static Quat align_step(Quat cur, Quat tgt, float dt, const AlignCfg* cfg) {
    Quat rel = qmul(qconj(cur), tgt);
    rel = qnorm(rel);
    if (rel.w < 0.0f) {                            /* enforce shortest arc */
        rel.x = -rel.x; rel.y = -rel.y; rel.z = -rel.z; rel.w = -rel.w;
    }
    float ang = 2.0f * acosf(clampf(rel.w, -1.0f, 1.0f));
    if (ang <= cfg->eps) return tgt;

    float step = ang;
    if (cfg->tau > 0.0f) {
        float damped = ang * (1.0f - expf(-dt / cfg->tau));
        step = damped;
    }
    float cap = cfg->wmax * dt;
    if (cap < step) step = cap;
    if (step > ang) step = ang;                   /* never overshoot */

    V3 axis = v3(rel.x, rel.y, rel.z);
    float al = vlen(axis);
    Quat dq = al > 1e-9f ? qaxis(vmul(axis, 1.0f / al), step) : qid();
    return qnorm(qmul(dq, cur));
}

/*
 * Two-stage per-frame orientation update, robust to a moving normal:
 *   stage 1 ground fit: shortest arc taking local +y onto the filtered normal
 *                       (axis lies in the ground plane, never flips at pi);
 *   stage 2 yaw:        shortest-arc rotation about the new up axis, driving
 *                       local forward (-z) toward the projected heading.
 * Building one full target matrix per frame is avoided: a fixed world heading
 * becomes parallel to a rotating normal twice per revolution, which makes the
 * projected basis degenerate and flip. The two stages stay well defined.
 */
static Quat orient_update(Quat cur, V3 normal, V3 heading, float dt,
                          const AlignCfg* cfg) {
    V3 up = vnorm(normal);

    /* stage 1: fit up to the ground normal along the shortest arc */
    V3 cur_up = qrot(cur, v3(0, 1, 0));
    Quat ground = qbetween(cur_up, up);
    Quat pose = qnorm(qmul(ground, cur));

    /* stage 2: yaw about the new up toward the projected heading */
    V3 f = vsub(heading, vmul(up, vdot(heading, up)));
    if (vlen(f) >= 1e-6f) {
        f = vnorm(f);
        V3 cur_fwd = qrot(pose, v3(0, 0, -1.0f));  /* local forward = local -z */
        V3 want_fwd = f;                           /* drive it onto the heading */
        /* shortest in-plane angle: both vectors are (nearly) perpendicular to up */
        float yaw = atan2f(vdot(up, vcross(cur_fwd, want_fwd)),
                           vdot(cur_fwd, want_fwd));
        if (fabsf(yaw) > cfg->eps) {
            float step = cfg->tau > 0.0f
                       ? yaw * (1.0f - expf(-dt / cfg->tau)) : yaw;
            float cap = cfg->wmax * dt;
            if (fabsf(cap) < fabsf(step)) step = copysignf(cap, step);
            if (fabsf(step) > fabsf(yaw)) step = yaw;
            pose = qnorm(qmul(qaxis(up, step), pose));
        }
    }
    return pose;
}

/* ---------------- self tests ---------------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

/* analytic ground truth rotation built directly from normal + heading */
static Mat3 analytic_matrix(V3 normal, V3 heading) {
    V3 up = vnorm(normal);
    V3 f = vsub(heading, vmul(up, vdot(heading, up)));
    f = vnorm(f);
    V3 right = vnorm(vcross(f, up));
    Mat3 m = { right, up, vmul(f, -1.0f) };
    return m;
}

/* signed phase of a unit vector rotating about the x axis: atan2(z, y) */
static float phase_x(V3 n) { return atan2f(n.z, n.y); }

static float wrap_pi(float a) {
    while (a > (float)M_PI) a -= 2.0f * (float)M_PI;
    while (a < -(float)M_PI) a += 2.0f * (float)M_PI;
    return a;
}

int main(void) {
    const float dt = 1.0f / 60.0f;
    const float EPS = 1e-4f;
    V3 local_up = v3(0, 1, 0);

    /* ---- 1. analytic construction vs. quaternion path; up error bound ---- */
    V3 slope_n = vnorm(v3(0, 1, 0.5f));
    V3 heading = v3(0, 0, 1);
    Mat3 ref = analytic_matrix(slope_n, heading);
    Quat qt = orient_target(slope_n, heading);
    Mat3 mt = qtomat(qt);

    check("matrix_col0", vlen(vsub(mt.c0, ref.c0)) < EPS);
    check("matrix_col1", vlen(vsub(mt.c1, ref.c1)) < EPS);
    check("matrix_col2", vlen(vsub(mt.c2, ref.c2)) < EPS);
    check("orthonormal", fabsf(vdot(mt.c0, mt.c1)) < EPS &&
                         fabsf(vdot(mt.c1, mt.c2)) < EPS &&
                         fabsf(vlen(mt.c0) - 1.0f) < EPS);
    check("right_handed", fabsf(mdet(mt) - 1.0f) < EPS);
    check("up_matches_normal", vangle(qrot(qt, local_up), slope_n) < EPS);

    /* ---- 2. slope walk: convergence, cap, jitter frames = 0 ---- */
    AlignCfg cfg = { 3.0f, 0.15f, 1e-5f };
    Quat cur = qid();                              /* flat ground start: feet level */
    Quat tgt = orient_target(slope_n, heading);
    float init_up_err = vangle(qrot(cur, local_up), slope_n);
    float max_up_err = 0.0f, max_speed = 0.0f;
    int jitter_frames = 0;
    int sign = 0;
    V3 last_axis = v3(0, 0, 0);
    for (int i = 0; i < 600; i++) {
        Quat prev = cur;
        cur = align_step(cur, tgt, dt, &cfg);
        float e = vangle(qrot(cur, local_up), slope_n);
        if (e > max_up_err) max_up_err = e;

        Quat dq = qnorm(qmul(qconj(prev), cur));
        float da = fabsf(2.0f * acosf(clampf(dq.w, -1.0f, 1.0f)));
        float speed = da / dt;
        if (speed > max_speed) max_speed = speed;

        /* jitter: an appreciable step whose axis flips against the prior one */
        V3 ax = v3(dq.x, dq.y, dq.z);
        if (vlen(ax) > 1e-7f && da > 1e-5f) {
            ax = vnorm(ax);
            int s = vdot(ax, last_axis) < 0.0f ? -1 : 1;
            if (sign != 0 && s != sign) jitter_frames++;
            sign = s;
            last_axis = ax;
        }
    }
    float final_up_err = vangle(qrot(cur, local_up), slope_n);
    check("up_err_bounded", max_up_err <= init_up_err + EPS);
    check("up_err_converges", final_up_err < 1e-3f);
    check("speed_capped", max_speed <= cfg.wmax + 1e-3f);
    check("jitter_frames_zero", jitter_frames == 0);
    printf("up_err_init_rad=%.6f up_err_max_rad=%.6f up_err_final_rad=%.6f\n",
           init_up_err, max_up_err, final_up_err);
    printf("max_angular_speed_radps=%.6f cap_radps=%.6f\n", max_speed, cfg.wmax);
    printf("jitter_frames=%d\n", jitter_frames);

    /* ---- 3. filter lag: numeric bias on a continuously rotating normal ---- */
    const float tau_f = 0.25f;
    const float omega = 1.0f;                      /* rad/s about x */
    NormalFilter nf;
    nf_init(&nf, tau_f);
    float t = 0.0f;
    for (int i = 0; i < 4000; i++) {               /* reach steady state */
        V3 raw = v3(0, cosf(omega * t), sinf(omega * t));
        nf_push(&nf, raw, dt);
        t += dt;
    }
    /* collect steady-state phase lag over two periods */
    float lag_sum = 0.0f;
    int lag_n = 0;
    for (int i = 0; i < 720; i++) {
        V3 raw = v3(0, cosf(omega * t), sinf(omega * t));
        nf_push(&nf, raw, dt);
        lag_sum += wrap_pi(phase_x(raw) - phase_x(nf.n));
        lag_n++;
        t += dt;
    }
    float lag_meas = lag_sum / lag_n;

    /* closed-form references for the *linear* (unnormalized) EMA: the exact
       discrete transfer phase and the continuous first-order approximation.
       Renormalization adds a small extra lag, quantified by measured-vs-linear. */
    float a = dt / (tau_f + dt);
    float r = 1.0f - a;
    float theta = omega * dt;
    float lag_linear = atanf(r * sinf(theta) / (1.0f - r * cosf(theta)));
    float lag_cont = atanf(omega * tau_f);        /* continuous first-order approx */

    /* independent 2D re-implementation of the same filter (component EMA +
       renormalization) as the numeric reference for the normalized lag. */
    float ry = 1.0f, rz = 0.0f, rt = 0.0f;
    for (int i = 0; i < 4000; i++) {
        float uy = cosf(omega * rt), uz = sinf(omega * rt);
        float ny = r * ry + a * uy, nz = r * rz + a * uz;
        float nl = sqrtf(ny * ny + nz * nz);
        ry = ny / nl; rz = nz / nl;
        rt += dt;
    }
    float ref_sum = 0.0f;
    for (int i = 0; i < 720; i++) {
        float uy = cosf(omega * rt), uz = sinf(omega * rt);
        float ny = r * ry + a * uy, nz = r * rz + a * uz;
        float nl = sqrtf(ny * ny + nz * nz);
        ry = ny / nl; rz = nz / nl;
        ref_sum += wrap_pi(atan2f(uz, uy) - atan2f(rz, ry));
        rt += dt;
    }
    float lag_ref2d = ref_sum / 720.0f;

    printf("filter_lag_rad=%.6f reference_2d_rad=%.6f linear_ema_rad=%.6f "
           "continuous_approx_rad=%.6f\n",
           lag_meas, lag_ref2d, lag_linear, lag_cont);
    printf("renorm_extra_lag_rad=%.6f\n", lag_meas - lag_linear);
    check("filter_lag_matches_reference", fabsf(lag_meas - lag_ref2d) < 2e-6f);
    check("lag_positive", lag_meas > 0.0f);

    /* ---- 4. end-to-end: two-stage update tracks a moving normal; the only
       residual up bias is the filter lag ---- */
    NormalFilter nf2;
    nf_init(&nf2, tau_f);
    AlignCfg track = { 30.0f, 0.02f, 1e-6f };
    /* heading co-rotates 90 deg ahead of the normal inside the rotation plane,
       so it is never parallel to the normal; ground fit stays a pure x turn. */
    Quat pose = qid();
    t = 0.0f;
    float bias_sum = 0.0f, fwd_err_sum = 0.0f, e2e_max_speed = 0.0f;
    int bias_n = 0, e2e_jitter = 0, e2e_sign = 0;
    V3 e2e_last_axis = v3(0, 0, 0);
    for (int i = 0; i < 4000; i++) {
        V3 raw = v3(0, cosf(omega * t), sinf(omega * t));
        nf_push(&nf2, raw, dt);
        V3 head = v3(0, -raw.z, raw.y);           /* perpendicular, co-rotating */
        Quat prev = pose;
        pose = orient_update(pose, nf2.n, head, dt, &track);
        if (i == 0) prev = pose;                  /* first frame snaps to ground */

        Quat dq = qnorm(qmul(qconj(prev), pose));
        float da = fabsf(2.0f * acosf(clampf(dq.w, -1.0f, 1.0f)));
        float spd = da / dt;
        if (spd > e2e_max_speed) e2e_max_speed = spd;
        V3 ax = v3(dq.x, dq.y, dq.z);
        if (vlen(ax) > 1e-7f && da > 1e-5f) {
            ax = vnorm(ax);
            int s = vdot(ax, e2e_last_axis) < 0.0f ? -1 : 1;
            if (e2e_sign != 0 && s != e2e_sign) e2e_jitter++;
            e2e_sign = s;
            e2e_last_axis = ax;
        }

        if (i >= 3000) {
            bias_sum += vangle(qrot(pose, local_up), raw);
            /* forward (-z) should lie along the projected heading on the ground */
            V3 fwd = qrot(pose, v3(0, 0, -1.0f));
            V3 fp = vsub(head, vmul(nf2.n, vdot(head, nf2.n)));
            fwd_err_sum += vangle(fwd, vnorm(fp));
            bias_n++;
        }
        t += dt;
    }
    float e2e_bias = bias_sum / bias_n;
    float e2e_fwd_err = fwd_err_sum / bias_n;
    printf("end_to_end_up_bias_rad=%.6f filter_lag_rad=%.6f\n", e2e_bias, lag_meas);
    printf("end_to_end_fwd_err_rad=%.6f max_speed_radps=%.6f jitter_frames=%d\n",
           e2e_fwd_err, e2e_max_speed, e2e_jitter);
    check("lag_is_the_bias", fabsf(e2e_bias - lag_meas) < 5e-3f);
    check("forward_tracks_heading", e2e_fwd_err < 1e-2f);
    check("e2e_speed_capped", e2e_max_speed <= track.wmax + 1e-3f);
    check("e2e_jitter_zero", e2e_jitter == 0);

    /* ---- 5. zero normal rejection ---- */
    NormalFilter nfz;
    nf_init(&nfz, 0.2f);
    check("reject_before_init", nf_push(&nfz, v3(0, 0, 0), dt) == 0);
    check("accept_first_valid", nf_push(&nfz, v3(0, 1, 0), dt) == 1);
    V3 kept = nfz.n;
    check("reject_zero_keeps_state", nf_push(&nfz, v3(1e-9f, 0, 0), dt) == 0 &&
                                     nfz.valid == 1 &&
                                     vlen(vsub(nfz.n, kept)) < EPS);
    check("recover_after_bad_sample", nf_push(&nfz, v3(0, 0, 1), dt) == 1);

    /* ---- 6. shortest-arc primitives and a planar 180-degree turn ---- */
    V3 from = v3(0, 0, 1), to = vnorm(v3(1, 0, 0.3f));
    Quat qb0 = qbetween(from, to);
    check("qbetween_rotates", vlen(vsub(qrot(qb0, from), to)) < EPS);
    check("qbetween_short", qb0.w >= 0.0f);
    Quat qanti = qbetween(v3(1, 0, 0), v3(-1, 0, 0));
    check("qbetween_antipodal", fabsf(2.0f * acosf(clampf(qanti.w, -1, 1)) -
                                      (float)M_PI) < 1e-4f);

    Quat qa = orient_target(v3(0, 1, 0), v3(0, 0, 1));
    Quat qb = orient_target(v3(0, 1, 0), v3(0, 0, -1));
    Quat rel = qmul(qconj(qa), qb);
    if (rel.w < 0.0f) { rel.x = -rel.x; rel.y = -rel.y; rel.z = -rel.z; rel.w = -rel.w; }
    float turn = 2.0f * acosf(clampf(rel.w, -1.0f, 1.0f));
    printf("turn_180_angle_rad=%.6f\n", turn);
    check("turn_is_pi", fabsf(turn - (float)M_PI) < 1e-4f);
    check("turn_about_up", fabsf(rel.x) < EPS && fabsf(rel.z) < EPS &&
                          fabsf(fabsf(rel.y) - 1.0f) < EPS);

    /* ---- 7. frame-rate independence of the damped, capped response ---- */
    Quat p1 = qid();
    for (int i = 0; i < 300; i++) p1 = align_step(p1, tgt, dt, &cfg);       /* 60 Hz, 5 s */
    Quat p2 = qid();
    for (int i = 0; i < 150; i++) p2 = align_step(p2, tgt, 2.0f * dt, &cfg); /* 30 Hz, 5 s */
    Quat rr = qmul(qconj(p1), p2);
    float drift = 2.0f * acosf(clampf(fabsf(rr.w), -1.0f, 1.0f));
    printf("framerate_drift_rad=%.6f\n", drift);
    check("framerate_independent", drift < 1e-2f);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
