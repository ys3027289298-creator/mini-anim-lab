
/* blend: distance-decayed pose blending in rotation space with time-varying weights.
 *
 * Rotations are angles (unit complex numbers).  A weighted rotation average is
 * taken on the circle (geodesic / Karcher mean), not by component-wise linear
 * interpolation, so the result is always a rotation (no length shrink) and the
 * +/-pi seam is handled (no jump).  Translation is interpolated linearly and
 * independently.  Weights are decay(distance) * user curve(time); a total that
 * is zero (or below a noise floor) is rejected.  The places where the weight
 * functions' first derivatives are discontinuous are reported, as is any time
 * at which the total weight vanishes.
 */
#include <stdio.h>
#include <string.h>
#include <math.h>

#define MAX_POSES 8
#define MAX_WKEYS 8
#define MAX_KINKS 24

#define PI_F 3.14159265358979323846f
#define TWO_PI_F 6.28318530717958647693f
#define WEIGHT_EPS 1e-6f      /* total weight at or below this is "zero" */
#define ROT_EPS 1e-6f         /* resultant below this: rotation mean ambiguous */
#define SLOPE_EPS 1e-6f

typedef struct { float x, y, ang; } Pose2;

typedef enum { DECAY_SMOOTH = 0, DECAY_LINEAR = 1, DECAY_GAUSS = 2 } DecayKernel;

typedef struct {
    DecayKernel kernel;
    float radius;             /* support radius (SMOOTH/LINEAR), falloff scale (GAUSS) */
} DecayCfg;

/* piecewise-linear time weight; count == 0 means constant 1 */
typedef struct { float t, v; } WKey;
typedef struct { WKey k[MAX_WKEYS]; int count; } WeightCurve;

typedef struct {
    int count;
    Pose2 pose[MAX_POSES];
    float dist[MAX_POSES];          /* distance used by the decay kernel */
    WeightCurve time[MAX_POSES];    /* per-pose time factor */
} PoseSet;

typedef enum {
    KINK_SPATIAL_RADIUS = 0,   /* decay kernel's slope changes at d = radius */
    KINK_TIME_SLOPE = 1,       /* time weight curve's slope changes */
    KINK_TIME_SUM_ZERO = 2     /* total weight is zero at this time */
} KinkKind;

typedef struct { KinkKind kind; float at; } Kink;
typedef struct { Kink items[MAX_KINKS]; int count; } KinkReport;

enum {
    BLEND_OK = 0,
    BLEND_ERR_BAD_COUNT = -1,
    BLEND_ERR_ZERO_WEIGHT = -2,
    BLEND_ERR_AMBIGUOUS_ROT = -3,
    BLEND_ERR_NEG_WEIGHT = -4,
    BLEND_ERR_ZERO_RADIUS = -5
};

typedef struct {
    int ok;
    int err;
    Pose2 pose;
    float w[MAX_POSES];    /* effective decay * time weights at the query */
    float sum;             /* sum of w[]; > WEIGHT_EPS when ok */
    float rho;             /* |sum w_i * z_i| / sum: chordal resultant, in [0,1] */
    float shrink;          /* 1 - rho: magnitude a raw linear average loses */
    int iters;             /* Karcher iteration count */
    float residual;        /* |last Karcher step| */
} BlendResult;

static float wrap_pi(float a) {
    a = fmodf(a + PI_F, TWO_PI_F);
    if (a < 0.0f) a += TWO_PI_F;
    return a - PI_F;       /* [-pi, pi) */
}

/* ---- distance decay ---------------------------------------------------- */

static float decay_eval(const DecayCfg* dc, float d) {
    if (!dc || dc->radius <= 0.0f) return 0.0f;
    float x = fabsf(d) / dc->radius;
    switch (dc->kernel) {
        case DECAY_SMOOTH:                   /* C1 compact support (1-x^2)^2 */
            if (x >= 1.0f) return 0.0f;
            return (1.0f - x * x) * (1.0f - x * x);
        case DECAY_LINEAR:                   /* C0: derivative jumps at x = 1 */
            if (x >= 1.0f) return 0.0f;
            return 1.0f - x;
        case DECAY_GAUSS:                    /* smooth, no compact support */
        default:
            return expf(-x * x);
    }
}

static int kink_push(KinkReport* r, KinkKind kind, float at) {
    for (int i = 0; i < r->count; i++) {
        if (r->items[i].kind == kind && fabsf(r->items[i].at - at) < 1e-5f) return i;
    }
    if (r->count >= MAX_KINKS) return -1;
    r->items[r->count].kind = kind;
    r->items[r->count].at = at;
    return r->count++;
}

static void decay_kinks(const DecayCfg* dc, KinkReport* rep) {
    /* LINEAR falls off with slope -1/r inside and 0 outside, so its first
     * derivative is discontinuous at d = radius.  SMOOTH and GAUSS are C1. */
    if (dc && dc->kernel == DECAY_LINEAR && dc->radius > 0.0f) {
        kink_push(rep, KINK_SPATIAL_RADIUS, dc->radius);
    }
}

/* ---- time weight curves ------------------------------------------------ */

static int wcurve_add(WeightCurve* c, float t, float v) {
    if (!c || c->count >= MAX_WKEYS) return -1;
    if (c->count > 0 && !(t > c->k[c->count - 1].t)) return -1;  /* strictly increasing */
    c->k[c->count].t = t;
    c->k[c->count].v = v;
    return c->count++;
}

static float wcurve_eval(const WeightCurve* c, float t) {
    if (!c || c->count == 0) return 1.0f;
    if (t <= c->k[0].t) return c->k[0].v;
    if (t >= c->k[c->count - 1].t) return c->k[c->count - 1].v;
    for (int i = 1; i < c->count; i++) {
        if (t <= c->k[i].t) {
            float dt = c->k[i].t - c->k[i - 1].t;
            float u = (t - c->k[i - 1].t) / dt;
            return c->k[i - 1].v + u * (c->k[i].v - c->k[i - 1].v);
        }
    }
    return c->k[c->count - 1].v;
}

static float wcurve_slope(const WeightCurve* c, int seg) {
    float dt = c->k[seg + 1].t - c->k[seg].t;
    if (dt <= 0.0f) return 0.0f;
    return (c->k[seg + 1].v - c->k[seg].v) / dt;
}

static void curve_kinks(const WeightCurve* c, float lo, float hi, KinkReport* rep) {
    if (!c || c->count < 2) return;   /* constant: no slope change */
    int n = c->count;
    for (int j = 1; j < n - 1; j++) { /* interior corners */
        float sl = wcurve_slope(c, j - 1);
        float sr = wcurve_slope(c, j);
        if (fabsf(sr - sl) > SLOPE_EPS && c->k[j].t >= lo && c->k[j].t <= hi) {
            kink_push(rep, KINK_TIME_SLOPE, c->k[j].t);
        }
    }
    /* outside the key range the curve is clamped, so slope 0 on that side */
    if (fabsf(wcurve_slope(c, 0)) > SLOPE_EPS && c->k[0].t >= lo && c->k[0].t <= hi) {
        kink_push(rep, KINK_TIME_SLOPE, c->k[0].t);
    }
    if (fabsf(wcurve_slope(c, n - 2)) > SLOPE_EPS && c->k[n - 1].t >= lo && c->k[n - 1].t <= hi) {
        kink_push(rep, KINK_TIME_SLOPE, c->k[n - 1].t);
    }
}

/* ---- weights ----------------------------------------------------------- */

static float pose_weight(const PoseSet* ps, const DecayCfg* dc, int i, float t) {
    return decay_eval(dc, ps->dist[i]) * wcurve_eval(&ps->time[i], t);
}

static float weight_sum_at(const PoseSet* ps, const DecayCfg* dc, float t) {
    float s = 0.0f;
    for (int i = 0; i < ps->count; i++) s += pose_weight(ps, dc, i, t);
    return s;
}

static int compute_weights(const PoseSet* ps, const DecayCfg* dc, float t,
                           float w[], float* sum) {
    if (!ps || !dc || ps->count < 1 || ps->count > MAX_POSES) return BLEND_ERR_BAD_COUNT;
    if (!(dc->radius > 0.0f)) return BLEND_ERR_ZERO_RADIUS;
    float s = 0.0f;
    for (int i = 0; i < ps->count; i++) {
        float wi = pose_weight(ps, dc, i, t);
        if (!(wi >= 0.0f) || !isfinite(wi)) return BLEND_ERR_NEG_WEIGHT;  /* NaN too */
        w[i] = wi;
        s += wi;
    }
    if (!(s > WEIGHT_EPS)) return BLEND_ERR_ZERO_WEIGHT;
    *sum = s;
    return BLEND_OK;
}

/* ---- rotation averaging ------------------------------------------------ */

static int pick_heaviest(const float w[], int n) {
    int best = -1;
    float bw = -1.0f;
    for (int i = 0; i < n; i++) {
        if (w[i] > 0.0f && w[i] > bw) { bw = w[i]; best = i; }
    }
    return best;
}

/* Geodesic (Karcher) weighted mean on the circle.  The step is the weighted
 * mean of the tangent-space logs of the inputs relative to the current
 * estimate; when all inputs lie in an arc shorter than pi this is a linear
 * problem and the iteration reaches the lifted weighted mean in one step.
 * Otherwise the mean may not be unique; iterations and residual are reported
 * so callers can tell. */
static int rot_mean(const PoseSet* ps, const float w[], int n, float sum,
                    float* out, float* rho, int* iters, float* residual) {
    int best = pick_heaviest(w, n);
    if (best < 0) return BLEND_ERR_ZERO_WEIGHT;

    float cx = 0.0f, cy = 0.0f;
    for (int i = 0; i < n; i++) {
        if (w[i] <= 0.0f) continue;
        cx += w[i] * cosf(ps->pose[i].ang);
        cy += w[i] * sinf(ps->pose[i].ang);
    }
    float r = hypotf(cx, cy) / sum;
    if (r < ROT_EPS) return BLEND_ERR_AMBIGUOUS_ROT;   /* directions cancel */

    float theta = ps->pose[best].ang;   /* start on the data, inside any arc */
    float step = 0.0f;
    int it = 0;
    for (; it < 64; it++) {
        float g = 0.0f;
        for (int i = 0; i < n; i++) {
            if (w[i] <= 0.0f) continue;
            g += w[i] * wrap_pi(ps->pose[i].ang - theta);
        }
        step = g / sum;
        theta = wrap_pi(theta + step);
        if (fabsf(step) < 1e-7f) { it++; break; }
    }
    *out = theta;
    *rho = r;
    *iters = it;
    *residual = fabsf(step);
    return BLEND_OK;
}

/* ---- blending ---------------------------------------------------------- */

static BlendResult blend_poses(const PoseSet* ps, const DecayCfg* dc, float t) {
    BlendResult r;
    memset(&r, 0, sizeof(r));
    r.err = compute_weights(ps, dc, t, r.w, &r.sum);
    if (r.err != BLEND_OK) { r.ok = 0; return r; }

    float sx = 0.0f, sy = 0.0f;
    for (int i = 0; i < ps->count; i++) {
        if (r.w[i] <= 0.0f) continue;
        sx += r.w[i] * ps->pose[i].x;
        sy += r.w[i] * ps->pose[i].y;
    }
    r.pose.x = sx / r.sum;
    r.pose.y = sy / r.sum;              /* translation: linear, independent */

    r.err = rot_mean(ps, r.w, ps->count, r.sum, &r.pose.ang, &r.rho, &r.iters, &r.residual);
    if (r.err != BLEND_OK) { r.ok = 0; return r; }
    r.shrink = 1.0f - r.rho;
    r.ok = 1;
    return r;
}

/* Sequential pairwise reference: merge pose i into the running result with
 * u = w_i / (sum so far).  Translation telescopes to the weighted mean; the
 * rotation walks the shortest arc from the running estimate to each pose,
 * so the result depends on the merge order. */
static BlendResult blend_pairwise(const PoseSet* ps, const DecayCfg* dc, float t) {
    BlendResult r;
    memset(&r, 0, sizeof(r));
    r.err = compute_weights(ps, dc, t, r.w, &r.sum);
    if (r.err != BLEND_OK) { r.ok = 0; return r; }

    float acc = 0.0f;
    int have = 0;
    int merges = 0;
    for (int i = 0; i < ps->count; i++) {
        if (r.w[i] <= 0.0f) continue;
        if (!have) {
            r.pose = ps->pose[i];
            acc = r.w[i];
            have = 1;
            continue;
        }
        float nw = acc + r.w[i];
        float u = r.w[i] / nw;
        r.pose.x += u * (ps->pose[i].x - r.pose.x);
        r.pose.y += u * (ps->pose[i].y - r.pose.y);
        r.pose.ang = wrap_pi(r.pose.ang + u * wrap_pi(ps->pose[i].ang - r.pose.ang));
        acc = nw;
        merges++;
    }
    if (!have) { r.ok = 0; r.err = BLEND_ERR_ZERO_WEIGHT; return r; }

    float cx = 0.0f, cy = 0.0f;
    for (int i = 0; i < ps->count; i++) {
        if (r.w[i] <= 0.0f) continue;
        cx += r.w[i] * cosf(ps->pose[i].ang);
        cy += r.w[i] * sinf(ps->pose[i].ang);
    }
    r.rho = hypotf(cx, cy) / r.sum;
    r.shrink = 1.0f - r.rho;
    r.iters = merges;
    r.residual = 0.0f;
    r.ok = 1;
    return r;
}

static float angle_error(float a, float b) {
    return fabsf(wrap_pi(a - b));
}

/* ---- angle error bound ------------------------------------------------- */

/* Minimal arc containing all positively-weighted input angles. */
static float pose_arc_span(const PoseSet* ps, const float w[]) {
    float a[MAX_POSES];
    int m = 0;
    for (int i = 0; i < ps->count; i++) {
        if (w[i] > 0.0f) a[m++] = wrap_pi(ps->pose[i].ang);
    }
    if (m < 2) return 0.0f;
    for (int i = 1; i < m; i++) {           /* insertion sort, n <= 8 */
        float v = a[i];
        int j = i - 1;
        while (j >= 0 && a[j] > v) { a[j + 1] = a[j]; j--; }
        a[j + 1] = v;
    }
    float max_gap = a[0] + TWO_PI_F - a[m - 1];
    for (int i = 1; i < m; i++) {
        float g = a[i] - a[i - 1];
        if (g > max_gap) max_gap = g;
    }
    float span = TWO_PI_F - max_gap;
    return span < 0.0f ? 0.0f : span;
}

static float max_angular_dev(const PoseSet* ps, const float w[], float theta) {
    float dmax = 0.0f;
    for (int i = 0; i < ps->count; i++) {
        if (w[i] <= 0.0f) continue;
        float d = angle_error(ps->pose[i].ang, theta);
        if (d > dmax) dmax = d;
    }
    return dmax;
}

/* Upper bound on |theta_a - theta_b| for two blends of the same weighted
 * poses.  If the poses fit in an arc of span <= pi, positive weights keep
 * both the geodesic mean and the sequential pairwise walk inside that arc,
 * so their separation is at most the span.  Otherwise fall back to the
 * triangle inequality through the inputs (capped at pi, the diameter). */
static float blend_angle_error_bound(const PoseSet* ps, const float w[],
                                     float theta_a, float theta_b) {
    float span = pose_arc_span(ps, w);
    if (span <= PI_F) return span;
    float b = max_angular_dev(ps, w, theta_a) + max_angular_dev(ps, w, theta_b);
    return b > PI_F ? PI_F : b;
}

typedef struct {
    BlendResult geodesic;
    BlendResult pairwise;
    float angle_err;      /* |geodesic.ang - pairwise.ang| */
    float angle_bound;    /* certified upper bound on angle_err */
} BlendCompare;

/* Blend both ways over the same weights and certify the rotation gap. */
static BlendCompare blend_compare(const PoseSet* ps, const DecayCfg* dc, float t) {
    BlendCompare c;
    memset(&c, 0, sizeof(c));
    c.geodesic = blend_poses(ps, dc, t);
    c.pairwise = blend_pairwise(ps, dc, t);
    if (!c.geodesic.ok || !c.pairwise.ok) return c;
    c.angle_err = angle_error(c.geodesic.pose.ang, c.pairwise.pose.ang);
    c.angle_bound = blend_angle_error_bound(ps, c.geodesic.w,
                                            c.geodesic.pose.ang, c.pairwise.pose.ang);
    return c;
}

/* ---- kink reporting ---------------------------------------------------- */

static int collect_breaks(const PoseSet* ps, float lo, float hi, float out[], int cap) {
    int n = 0;
    out[n++] = lo;
    for (int i = 0; i < ps->count; i++) {
        const WeightCurve* c = &ps->time[i];
        for (int j = 0; j < c->count; j++) {
            if (c->k[j].t > lo && c->k[j].t < hi && n < cap) out[n++] = c->k[j].t;
        }
    }
    if (n < cap && hi > lo) out[n++] = hi;
    for (int i = 1; i < n; i++) {           /* insertion sort */
        float v = out[i];
        int j = i - 1;
        while (j >= 0 && out[j] > v) { out[j + 1] = out[j]; j--; }
        out[j + 1] = v;
    }
    int m = 0;
    for (int i = 0; i < n; i++) {           /* drop duplicates */
        if (m == 0 || fabsf(out[i] - out[m - 1]) > 1e-6f) out[m++] = out[i];
    }
    return m;
}

/* The total weight is piecewise linear between the curve keys (decay is
 * constant in time), so a zero total over time is located exactly: check the
 * breakpoints, then solve any sign change inside a piece. */
static void sum_zero_kinks(const PoseSet* ps, const DecayCfg* dc, float lo, float hi,
                           KinkReport* rep) {
    float bps[MAX_POSES * MAX_WKEYS + 2];
    int nb = collect_breaks(ps, lo, hi, bps, (int)(sizeof(bps) / sizeof(bps[0])));
    if (nb < 1) return;
    for (int j = 0; j < nb; j++) {
        float s = weight_sum_at(ps, dc, bps[j]);
        if (fabsf(s) <= WEIGHT_EPS) kink_push(rep, KINK_TIME_SUM_ZERO, bps[j]);
    }
    for (int j = 0; j + 1 < nb; j++) {
        float sa = weight_sum_at(ps, dc, bps[j]);
        float sb = weight_sum_at(ps, dc, bps[j + 1]);
        if (sa * sb < 0.0f) {
            float tz = bps[j] + (bps[j + 1] - bps[j]) * (sa / (sa - sb));
            kink_push(rep, KINK_TIME_SUM_ZERO, tz);
        }
    }
}

/* Reports all first-derivative discontinuities of the weight functions over
 * [lo, hi]: the decay kernel's break at the support radius, every corner of
 * the time weight curves (and their clamped ends), and every time at which
 * the total weight is zero. */
static void blend_kinks(const PoseSet* ps, const DecayCfg* dc, float lo, float hi,
                        KinkReport* rep) {
    memset(rep, 0, sizeof(*rep));
    decay_kinks(dc, rep);
    for (int i = 0; i < ps->count; i++) {
        /* a pose with zero spatial decay has a constant zero weight: no kinks */
        if (decay_eval(dc, ps->dist[i]) <= 0.0f) continue;
        curve_kinks(&ps->time[i], lo, hi, rep);
    }
    sum_zero_kinks(ps, dc, lo, hi, rep);
}

/* ---- tests ------------------------------------------------------------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static int report_has(const KinkReport* r, KinkKind kind, float at) {
    for (int i = 0; i < r->count; i++) {
        if (r->items[i].kind == kind && fabsf(r->items[i].at - at) < 1e-4f) return 1;
    }
    return 0;
}

int main(void) {
    /* --- decay kernels --- */
    DecayCfg smooth = { DECAY_SMOOTH, 1.0f };
    DecayCfg linear = { DECAY_LINEAR, 1.0f };
    DecayCfg gauss = { DECAY_GAUSS, 1.0f };
    check("decay_center", fabsf(decay_eval(&smooth, 0.0f) - 1.0f) < 1e-6f &&
                          fabsf(decay_eval(&linear, 0.0f) - 1.0f) < 1e-6f &&
                          fabsf(decay_eval(&gauss, 0.0f) - 1.0f) < 1e-6f);
    check("decay_at_radius", fabsf(decay_eval(&smooth, 1.0f)) < 1e-6f &&
                             fabsf(decay_eval(&linear, 1.0f)) < 1e-6f);
    check("decay_beyond_radius", decay_eval(&smooth, 2.0f) == 0.0f &&
                                 decay_eval(&linear, 2.0f) == 0.0f);
    check("decay_monotone", decay_eval(&smooth, 0.3f) > decay_eval(&smooth, 0.7f) &&
                            decay_eval(&gauss, 0.3f) > decay_eval(&gauss, 0.7f));

    /* --- derivative discontinuity report --- */
    PoseSet kp;
    memset(&kp, 0, sizeof(kp));
    kp.count = 2;
    kp.pose[0] = (Pose2){ 0.0f, 0.0f, 0.0f };  kp.dist[0] = 0.2f;
    kp.pose[1] = (Pose2){ 1.0f, 0.0f, 1.0f };  kp.dist[1] = 0.4f;
    wcurve_add(&kp.time[1], 0.0f, 1.0f);
    wcurve_add(&kp.time[1], 0.5f, 0.25f);
    wcurve_add(&kp.time[1], 1.0f, 1.0f);

    KinkReport kr;
    blend_kinks(&kp, &linear, 0.0f, 1.0f, &kr);
    check("kink_spatial_radius", report_has(&kr, KINK_SPATIAL_RADIUS, 1.0f));
    check("kink_time_corner", report_has(&kr, KINK_TIME_SLOPE, 0.5f));
    check("kink_time_ends", report_has(&kr, KINK_TIME_SLOPE, 0.0f) &&
                            report_has(&kr, KINK_TIME_SLOPE, 1.0f));

    KinkReport kg;
    blend_kinks(&kp, &smooth, 0.0f, 1.0f, &kg);
    check("smooth_no_spatial_kink", !report_has(&kg, KINK_SPATIAL_RADIUS, 1.0f));

    PoseSet constant = kp;
    memset(&constant.time[1], 0, sizeof(constant.time[1]));
    KinkReport kc;
    blend_kinks(&constant, &smooth, 0.0f, 1.0f, &kc);
    check("constant_curve_no_time_kink", kc.count == 0);

    /* numeric confirmation of both reported jumps */
    float h = 1e-3f;
    float f0 = pose_weight(&kp, &linear, 1, 0.5f);
    float fl = pose_weight(&kp, &linear, 1, 0.5f - h);
    float fr = pose_weight(&kp, &linear, 1, 0.5f + h);
    check("time_kink_numeric_jump", fabsf((fr - f0) / h - (f0 - fl) / h) > 0.05f);

    float e = 1e-3f;
    float lj = (decay_eval(&linear, 1.0f) - decay_eval(&linear, 1.0f - e)) / e
             - (decay_eval(&linear, 1.0f + e) - decay_eval(&linear, 1.0f)) / e;
    float sj = (decay_eval(&smooth, 1.0f) - decay_eval(&smooth, 1.0f - e)) / e
             - (decay_eval(&smooth, 1.0f + e) - decay_eval(&smooth, 1.0f)) / e;
    check("spatial_kink_numeric", fabsf(lj) > 0.5f && fabsf(sj) < 0.05f);

    /* --- zero total weight is rejected --- */
    PoseSet zr;
    memset(&zr, 0, sizeof(zr));
    zr.count = 2;
    zr.pose[0] = (Pose2){ 0.0f, 0.0f, 0.0f };  zr.dist[0] = 2.0f;
    zr.pose[1] = (Pose2){ 1.0f, 0.0f, 1.0f };  zr.dist[1] = 3.0f;
    BlendResult rz = blend_poses(&zr, &linear, 0.0f);
    check("reject_zero_sum", rz.ok == 0 && rz.err == BLEND_ERR_ZERO_WEIGHT && rz.sum == 0.0f);
    BlendResult rzp = blend_pairwise(&zr, &linear, 0.0f);
    check("reject_zero_sum_pairwise", rzp.ok == 0 && rzp.err == BLEND_ERR_ZERO_WEIGHT);

    /* --- a single pose is reproduced exactly --- */
    PoseSet one;
    memset(&one, 0, sizeof(one));
    one.count = 1;
    one.pose[0] = (Pose2){ 1.3f, -0.7f, 0.4f };
    one.dist[0] = 0.1f;
    BlendResult r1 = blend_poses(&one, &smooth, 0.5f);
    check("single_pose_exact", r1.ok &&
          fabsf(r1.pose.x - 1.3f) < 1e-6f && fabsf(r1.pose.y + 0.7f) < 1e-6f &&
          angle_error(r1.pose.ang, 0.4f) < 1e-6f && fabsf(r1.rho - 1.0f) < 1e-6f);

    /* --- rotations across the +/-pi seam --- */
    PoseSet seam;
    memset(&seam, 0, sizeof(seam));
    seam.count = 2;
    seam.pose[0] = (Pose2){ 0.0f, 0.0f, 3.0f };   seam.dist[0] = 0.0f;
    seam.pose[1] = (Pose2){ 1.0f, 0.0f, -3.0f };  seam.dist[1] = 0.0f;
    BlendResult rs = blend_poses(&seam, &smooth, 0.0f);
    BlendResult rsp = blend_pairwise(&seam, &smooth, 0.0f);
    float naive = 0.5f * 3.0f + 0.5f * (-3.0f);  /* linear angle lerp */
    check("seam_geodesic_mean", rs.ok && fabsf(fabsf(rs.pose.ang) - PI_F) < 1e-3f);
    check("seam_pairwise_agrees", angle_error(rs.pose.ang, rsp.pose.ang) < 1e-3f);
    check("seam_naive_far_off", angle_error(naive, rs.pose.ang) > 3.0f);
    check("seam_span_is_small", pose_arc_span(&seam, rs.w) < 0.4f);
    check("seam_bound_tight", blend_angle_error_bound(&seam, rs.w, rs.pose.ang, rsp.pose.ang)
                              < 0.4f);

    /* --- 0 and 120 degrees, equal weights --- */
    PoseSet spread;
    memset(&spread, 0, sizeof(spread));
    spread.count = 2;
    spread.pose[0] = (Pose2){ 0.0f, 0.0f, 0.0f };                 spread.dist[0] = 0.0f;
    spread.pose[1] = (Pose2){ 0.0f, 0.0f, 2.0f * PI_F / 3.0f };    spread.dist[1] = 0.0f;
    BlendResult rsp2 = blend_poses(&spread, &smooth, 0.0f);
    check("chordal_shrinks_to_half", fabsf(rsp2.shrink - 0.5f) < 0.01f);
    check("rotation_mean_unit", angle_error(rsp2.pose.ang, PI_F / 3.0f) < 1e-3f);
    check("rotation_mean_converges", rsp2.iters <= 3 && rsp2.residual < 1e-5f);

    /* --- three poses: translation and rotation comparison --- */
    PoseSet tri;
    memset(&tri, 0, sizeof(tri));
    tri.count = 3;
    tri.pose[0] = (Pose2){ 0.0f, 0.0f, 0.0f };   tri.dist[0] = 0.1f;
    tri.pose[1] = (Pose2){ 1.0f, 2.0f, 0.6f };   tri.dist[1] = 0.6f;
    tri.pose[2] = (Pose2){ -1.0f, 1.0f, -0.3f }; tri.dist[2] = 0.8f;
    BlendResult rt = blend_poses(&tri, &smooth, 0.0f);
    BlendResult rtp = blend_pairwise(&tri, &smooth, 0.0f);
    float wx = 0.0f, wy = 0.0f, ws = 0.0f;
    for (int i = 0; i < tri.count; i++) {
        wx += rt.w[i] * tri.pose[i].x;
        wy += rt.w[i] * tri.pose[i].y;
        ws += rt.w[i];
    }
    check("translation_weighted_mean",
          fabsf(rt.pose.x - wx / ws) < 1e-5f && fabsf(rt.pose.y - wy / ws) < 1e-5f);
    check("translation_pairwise_agrees",
          fabsf(rt.pose.x - rtp.pose.x) < 1e-5f && fabsf(rt.pose.y - rtp.pose.y) < 1e-5f);

    float err_t = angle_error(rt.pose.ang, rtp.pose.ang);
    float bound_t = blend_angle_error_bound(&tri, rt.w, rt.pose.ang, rtp.pose.ang);
    BlendCompare ctri = blend_compare(&tri, &smooth, 0.0f);
    check("angle_error_within_bound", err_t <= bound_t + 1e-5f);
    check("angle_error_bound_is_span", fabsf(bound_t - pose_arc_span(&tri, rt.w)) < 1e-5f);
    check("compare_exposes_same_err", ctri.geodesic.ok && ctri.pairwise.ok &&
          fabsf(ctri.angle_err - err_t) < 1e-6f &&
          fabsf(ctri.angle_bound - bound_t) < 1e-6f);

    /* rot mean is independent of input order when the span is below pi */
    PoseSet tri_rev = tri;
    tri_rev.pose[0] = tri.pose[2]; tri_rev.dist[0] = tri.dist[2];
    tri_rev.pose[2] = tri.pose[0]; tri_rev.dist[2] = tri.dist[0];
    BlendResult rtr = blend_poses(&tri_rev, &smooth, 0.0f);
    check("narrow_mean_order_invariant", angle_error(rt.pose.ang, rtr.pose.ang) < 1e-4f);

    /* --- wide spread (> pi): pairwise merges are order dependent --- */
    PoseSet wide;
    memset(&wide, 0, sizeof(wide));
    wide.count = 3;
    wide.pose[0] = (Pose2){ 0.0f, 0.0f, 0.0f };  wide.dist[0] = 0.0f;
    wide.pose[1] = (Pose2){ 1.0f, 0.0f, 2.5f };  wide.dist[1] = 0.0f;
    wide.pose[2] = (Pose2){ 2.0f, 0.0f, -2.0f }; wide.dist[2] = 0.0f;
    BlendResult rw = blend_poses(&wide, &smooth, 0.0f);
    BlendResult rwp = blend_pairwise(&wide, &smooth, 0.0f);

    PoseSet wide_rev = wide;
    wide_rev.pose[0] = wide.pose[2];
    wide_rev.pose[2] = wide.pose[0];
    BlendResult rwpr = blend_pairwise(&wide_rev, &smooth, 0.0f);

    BlendCompare cwide = blend_compare(&wide, &smooth, 0.0f);
    check("wide_arc_span", pose_arc_span(&wide, rw.w) > PI_F);
    check("pairwise_order_dependent", angle_error(rwp.pose.ang, rwpr.pose.ang) > 0.5f);
    float bound_w = blend_angle_error_bound(&wide, rw.w, rw.pose.ang, rwp.pose.ang);
    check("wide_bound_is_diameter", fabsf(bound_w - PI_F) < 1e-5f);
    check("wide_angle_error_within_bound",
          angle_error(rw.pose.ang, rwp.pose.ang) <= bound_w + 1e-5f);
    check("wide_compare_within_bound", cwide.angle_err <= cwide.angle_bound + 1e-5f);

    /* --- time-varying weights --- */
    PoseSet tv;
    memset(&tv, 0, sizeof(tv));
    tv.count = 2;
    tv.pose[0] = (Pose2){ 0.0f, 0.0f, 0.0f };  tv.dist[0] = 0.0f;
    tv.pose[1] = (Pose2){ 2.0f, 0.0f, 1.5f };  tv.dist[1] = 0.0f;
    wcurve_add(&tv.time[0], 0.0f, 1.0f);
    wcurve_add(&tv.time[0], 1.0f, 0.0f);
    wcurve_add(&tv.time[1], 0.0f, 0.0f);
    wcurve_add(&tv.time[1], 1.0f, 1.0f);
    BlendResult rt0 = blend_poses(&tv, &smooth, 0.0f);
    BlendResult rt1 = blend_poses(&tv, &smooth, 1.0f);
    BlendResult rth = blend_poses(&tv, &smooth, 0.5f);
    check("time_gate_a", rt0.ok && fabsf(rt0.pose.x) < 1e-6f &&
                         angle_error(rt0.pose.ang, 0.0f) < 1e-6f && rt0.w[0] == 1.0f);
    check("time_gate_b", rt1.ok && fabsf(rt1.pose.x - 2.0f) < 1e-6f &&
                         angle_error(rt1.pose.ang, 1.5f) < 1e-6f && rt1.w[0] == 0.0f);
    check("time_half", rth.ok && fabsf(rth.pose.x - 1.0f) < 1e-6f &&
                       angle_error(rth.pose.ang, 0.75f) < 1e-3f);

    /* a zero total weight can also appear as a function of time: report it */
    PoseSet zt = tv;
    memset(&zt.time[0], 0, sizeof(zt.time[0]));
    memset(&zt.time[1], 0, sizeof(zt.time[1]));
    wcurve_add(&zt.time[0], 0.0f, 1.0f);
    wcurve_add(&zt.time[0], 0.5f, 0.0f);
    wcurve_add(&zt.time[0], 1.0f, 1.0f);
    wcurve_add(&zt.time[1], 0.0f, 1.0f);
    wcurve_add(&zt.time[1], 0.5f, 0.0f);
    wcurve_add(&zt.time[1], 1.0f, 1.0f);
    KinkReport kz;
    blend_kinks(&zt, &smooth, 0.0f, 1.0f, &kz);
    check("kink_sum_zero_reported", report_has(&kz, KINK_TIME_SUM_ZERO, 0.5f));
    check("reject_zero_sum_at_time", blend_poses(&zt, &smooth, 0.5f).err == BLEND_ERR_ZERO_WEIGHT);

    /* --- bad inputs --- */
    PoseSet bad;
    memset(&bad, 0, sizeof(bad));
    bad.count = MAX_POSES + 1;
    check("reject_capacity", blend_poses(&bad, &smooth, 0.0f).err == BLEND_ERR_BAD_COUNT &&
                             blend_pairwise(&bad, &smooth, 0.0f).err == BLEND_ERR_BAD_COUNT);
    PoseSet empty;
    memset(&empty, 0, sizeof(empty));
    check("reject_empty", blend_poses(&empty, &smooth, 0.0f).err == BLEND_ERR_BAD_COUNT);
    DecayCfg badradius = { DECAY_SMOOTH, 0.0f };
    check("reject_zero_radius", blend_poses(&one, &badradius, 0.0f).err == BLEND_ERR_ZERO_RADIUS);

    PoseSet neg;
    memset(&neg, 0, sizeof(neg));
    neg.count = 2;
    neg.pose[0] = (Pose2){ 0.0f, 0.0f, 0.0f };  neg.dist[0] = 0.0f;
    neg.pose[1] = (Pose2){ 1.0f, 0.0f, 0.0f };  neg.dist[1] = 0.0f;
    wcurve_add(&neg.time[1], 0.0f, -1.0f);
    check("reject_negative_weight", blend_poses(&neg, &smooth, 0.0f).err == BLEND_ERR_NEG_WEIGHT);

    /* two opposite directions with equal weight: no well-defined mean */
    PoseSet anti;
    memset(&anti, 0, sizeof(anti));
    anti.count = 2;
    anti.pose[0] = (Pose2){ 0.0f, 0.0f, 0.0f };  anti.dist[0] = 0.0f;
    anti.pose[1] = (Pose2){ 1.0f, 0.0f, PI_F };  anti.dist[1] = 0.0f;
    BlendResult ra = blend_poses(&anti, &smooth, 0.0f);
    check("reject_ambiguous_rotation", !ra.ok && ra.err == BLEND_ERR_AMBIGUOUS_ROT);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
