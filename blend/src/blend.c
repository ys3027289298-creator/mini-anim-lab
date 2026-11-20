
/* blend: animation/physics blend continuity analysis.

   The blended pose is

       p(t) = (1 - w(t)) * anim(t) + w(t) * phys(t)

   with anim/phys as sample tables (linear between samples) and w an explicit
   piecewise weight curve. Differentiating:

       p'(t) = (1 - w) * anim' + w * phys' + w' * (phys - anim)

   So a velocity discontinuity can only come from a break in w' , a break in
   anim' / phys' (input velocity jumps), or from a jump of w itself. A jump of
   w (a hold segment ending) makes p' impulsive: the one-sided derivatives do
   not exist, and this is reported as undefined rather than blended silently.

   Each knot is a linear, smoothstep or hold segment. With u = (t-t0)/(t1-t0):

       linear    w = w0 + dw*u                w' = dw/dt
       smooth    w = w0 + dw*(3u^2 - 2u^3)    w' = dw*6u(1-u)/dt
       hold      w = w0 on [t0,t1)            w' = 0, then w jumps to w1

   Continuity conditions at an interior knot, with wl/wr the one-sided slopes:

   - smoothstep ends have w' = 0, so smooth|smooth is always C1;
   - linear|linear is C1 iff their secant slopes are equal;
   - linear|smoothstep is C1 iff that one side's slope is 0;
   - hold|anything jumps in w' (bounded velocity jump) unless the hold is
     flat, and anything|hold that ends the hold jumps in w (unbounded).
     These are reported, never averaged away.

   Boundary knots (t == t0 or t == T) are not breaks: the blend is clamped
   outside the knot domain. */
#include <stdio.h>
#include <string.h>
#include <math.h>

#define BLEND_MAX_KEYS 32
#define BLEND_MAX_BREAKS (2 * BLEND_MAX_KEYS)
#define BLEND_EPS 1e-6f
#define BLEND_ND_H 1e-3f
#define BLEND_R2D 57.29577951308232f

#define BB_WEIGHT_NOT_DIFF 1u   /* w' jumps: the weight has a corner        */
#define BB_WEIGHT_DISCONT  2u   /* w itself jumps (hold segment exit)        */
#define BB_INPUT_ANIM      4u   /* anim' jumps at a sample                   */
#define BB_INPUT_PHYS      8u   /* phys' jumps at a sample                   */
#define BB_DIRECTION       16u  /* motion direction angle > threshold        */
#define BB_WEIGHT_RANGE    32u  /* w outside [0, 1]                          */
#define BB_OUTSIDE         64u  /* t outside the blend domain                */

typedef enum {
    WSEG_LINEAR = 0,
    WSEG_SMOOTH,
    WSEG_HOLD
} WSegKind;

typedef struct { float t; float w; WSegKind kind; } WKey;
typedef struct { WKey k[BLEND_MAX_KEYS]; int n; } WeightCurve;

typedef struct {
    float dt;
    int n;
    float pos[BLEND_MAX_KEYS];
    float dx[BLEND_MAX_KEYS];   /* motion direction, need not be unit */
    float dy[BLEND_MAX_KEYS];
} Track;

typedef struct {
    float t;
    unsigned flags;
    float w, dw_lo, dw_hi;   /* w and its one-sided derivatives at t */
    float d;                 /* signed phys - anim at t             */
    float jump;              /* analytic |p'(t+) - p'(t-)|, INF if undefined */
    float bound;             /* analytic upper bound for |jump|     */
    float pose_jump;         /* hold only: |dw| * |d|               */
} Break;

typedef struct {
    int n_breaks;
    int weight_corners;      /* w' jumps, velocity bounded           */
    int weight_undefined;    /* w jumps, velocity undefined          */
    int input_breaks;        /* anim/phys velocity breaks            */
    int direction_warnings;
    float max_jump;
    float max_bound;
} Report;

typedef struct {
    float t0, t1;
    float w_min, w_max;
    int w_in_range;
} Interval;

typedef struct {
    int compared;
    int skipped;
    int n;
    float t[BLEND_MAX_KEYS];
    float exact[BLEND_MAX_KEYS];
    float numeric[BLEND_MAX_KEYS];
    int used[BLEND_MAX_KEYS];
    float max_abs_err;
    float t_worst;
} Check;

typedef struct {
    float measured;
    float jump;
    float bound;
    int undefined;
    int ok;                  /* measured within [jump-slack, bound+slack] */
} StepCheck;

typedef struct {
    int ok;
    float pose;
    float w;
    float angle_deg;
    unsigned flags;
} Gate;

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

/* --------------------------------------------------------------- tracks */

static float track_duration(const Track* tr) {
    return tr->n > 1 ? tr->dt * (float)(tr->n - 1) : 0.0f;
}

static int track_at_sample(const Track* tr, float t) {
    float x = t / tr->dt;
    int i = (int)floorf(x + 0.5f);
    if (fabsf(x - (float)i) <= 1e-4f) return i;
    return -1;
}

static float track_eval(const Track* tr, float t) {
    if (tr->n <= 0) return 0.0f;
    if (t <= 0.0f) return tr->pos[0];
    float x = t / tr->dt;
    int i = (int)floorf(x);
    if (i >= tr->n - 1) return tr->pos[tr->n - 1];
    float u = x - (float)i;
    return tr->pos[i] + (tr->pos[i + 1] - tr->pos[i]) * u;
}

static float track_slope(const Track* tr, int seg) {
    return (tr->pos[seg + 1] - tr->pos[seg]) / tr->dt;
}

static int track_seg_right(const Track* tr, float t) {
    float x = t / tr->dt;
    int i = (int)floorf(x);
    if (i < 0) i = 0;
    if (i > tr->n - 2) i = tr->n - 2;
    return i;
}

static int track_seg_left(const Track* tr, float t) {
    int s = track_at_sample(tr, t);
    if (s >= 1 && s <= tr->n - 2) return s - 1;
    return track_seg_right(tr, t);
}

/* Index of an interior sample whose slope jumps, else -1. Phantom bounds at
   t = 0 and t = duration are not breaks: the tables are clamped there. */
static int track_break_at(const Track* tr, float t) {
    int i = track_at_sample(tr, t);
    if (i <= 0 || i >= tr->n - 1) return -1;
    if (fabsf(track_slope(tr, i) - track_slope(tr, i - 1)) > 1e-5f) return i;
    return -1;
}

static float track_angle_deg(const Track* tr, float t) {
    int i = (int)floorf(t / tr->dt + 0.5f);
    if (i < 0) i = 0;
    if (i >= tr->n) i = tr->n - 1;
    if (tr->dx[i] * tr->dx[i] + tr->dy[i] * tr->dy[i] < 1e-12f) return NAN;
    return atan2f(tr->dy[i], tr->dx[i]) * BLEND_R2D;
}

/* Smallest angle between the two motion directions at the nearest sample.
   Returns 0 when either track is not moving (no direction to disagree on). */
static float blend_direction_deg(const Track* anim, const Track* phys, float t) {
    float a = track_angle_deg(anim, t);
    float b = track_angle_deg(phys, t);
    if (isnan(a) || isnan(b)) return 0.0f;
    float d = fmodf(a - b + 180.0f, 360.0f);
    if (d < 0.0f) d += 360.0f;
    return fabsf(d - 180.0f);
}

/* --------------------------------------------------------- weight curve */

static void weight_init(WeightCurve* c) { memset(c, 0, sizeof(*c)); }

static int weight_add(WeightCurve* c, float t, float w, WSegKind kind) {
    if (c->n >= BLEND_MAX_KEYS) return -1;
    if (c->n > 0 && !(t > c->k[c->n - 1].t)) return -1;   /* strictly increasing */
    c->k[c->n].t = t;
    c->k[c->n].w = w;
    c->k[c->n].kind = kind;
    return ++c->n - 1;
}

static int weight_segs(const WeightCurve* c) { return c->n > 1 ? c->n - 1 : 0; }

static float weight_eval_seg(const WeightCurve* c, int seg, float t) {
    const WKey* a = &c->k[seg];
    const WKey* b = &c->k[seg + 1];
    float u = (t - a->t) / (b->t - a->t);
    if (u <= 0.0f) return a->w;
    if (u >= 1.0f) return b->w;
    if (a->kind == WSEG_SMOOTH) return a->w + (b->w - a->w) * (u * u * (3.0f - 2.0f * u));
    if (a->kind == WSEG_HOLD) return a->w;
    return a->w + (b->w - a->w) * u;
}

/* Analytic derivative at t; on the closed segment end this is the left
   derivative for a segment ending at t and the right one for a segment
   starting at t. */
static float weight_slope_seg(const WeightCurve* c, int seg, float t) {
    const WKey* a = &c->k[seg];
    const WKey* b = &c->k[seg + 1];
    float d = b->t - a->t;
    float u = (t - a->t) / d;
    if (u < 0.0f) u = 0.0f;
    if (u > 1.0f) u = 1.0f;
    if (a->kind == WSEG_SMOOTH) return (b->w - a->w) * 6.0f * u * (1.0f - u) / d;
    if (a->kind == WSEG_HOLD) return 0.0f;
    return (b->w - a->w) / d;
}

static float weight_eval(const WeightCurve* c, float t) {
    if (c->n <= 0) return 0.0f;
    if (t <= c->k[0].t) return c->k[0].w;
    if (t >= c->k[c->n - 1].t) return c->k[c->n - 1].w;
    for (int s = 0; s < weight_segs(c); s++)
        if (t < c->k[s + 1].t) return weight_eval_seg(c, s, t);
    return c->k[c->n - 1].w;
}

/* Weight-side break at t. Fills the one-sided slopes and returns the flags:
   BB_WEIGHT_NOT_DIFF for a corner, BB_WEIGHT_DISCONT for a hold exit whose
   value jumps (pose_jump = |dw|*|d| is the finite pose jump). Zero means w is
   differentiable at t (possibly with w' = 0 inside a hold). */
static unsigned weight_break_at(const WeightCurve* c, float t,
                                float* dw_lo, float* dw_hi, float* dw_pose) {
    *dw_lo = 0.0f;
    *dw_hi = 0.0f;
    *dw_pose = 0.0f;
    int k = -1;
    for (int i = 1; i <= c->n - 2; i++)
        if (fabsf(t - c->k[i].t) <= BLEND_EPS) { k = i; break; }
    if (k < 0) return 0;

    if (c->k[k - 1].kind == WSEG_HOLD) {          /* hold ends at this knot */
        *dw_lo = 0.0f;
        *dw_hi = weight_slope_seg(c, k, t);
        *dw_pose = c->k[k].w - c->k[k - 1].w;
        if (fabsf(*dw_pose) <= BLEND_EPS) return 0;   /* flat hold: no jump */
        return BB_WEIGHT_DISCONT;
    }
    *dw_lo = weight_slope_seg(c, k - 1, t);
    *dw_hi = weight_slope_seg(c, k, t);
    if (fabsf(*dw_hi - *dw_lo) > BLEND_EPS) return BB_WEIGHT_NOT_DIFF;
    return 0;
}

/* --------------------------------------------------- derivative analysis */

/* Build the Break record at an arbitrary time t. Returns 1 if p' is
   discontinuous or undefined there. Both the exact signed jump and a
   triangle-inequality upper bound are produced:

     jump  = |(1-w)*[anim'] + w*[phys'] + [w']*(phys-anim)|
     bound = |[w']|*|phys-anim| + |[anim']|*|1-w| + |[phys']|*|w|

   where [x] is the one-sided change. bound >= |jump| always, and the two
   coincide when no input break shares the time. */
static int break_at_time(const WeightCurve* c, const Track* anim, const Track* phys,
                         float t, Break* br) {
    memset(br, 0, sizeof(*br));
    br->t = t;
    float dw_lo, dw_hi, dw_pose;
    unsigned flags = weight_break_at(c, t, &dw_lo, &dw_hi, &dw_pose);

    int al = track_seg_left(anim, t), ar = track_seg_right(anim, t);
    int pl = track_seg_left(phys, t), pr = track_seg_right(phys, t);
    float apL = track_slope(anim, al), apR = track_slope(anim, ar);
    float bpL = track_slope(phys, pl), bpR = track_slope(phys, pr);
    float dap = fabsf(apR - apL), dbp = fabsf(bpR - bpL);
    if (dap > 1e-5f) flags |= BB_INPUT_ANIM;
    if (dbp > 1e-5f) flags |= BB_INPUT_PHYS;
    if (flags == 0) return 0;

    float w = weight_eval(c, t);
    float a = track_eval(anim, t);
    float b = track_eval(phys, t);
    float d = b - a;

    br->flags = flags;
    br->w = w;
    br->dw_lo = dw_lo;
    br->dw_hi = dw_hi;
    br->d = d;

    if (flags & BB_WEIGHT_DISCONT) {
        br->jump = INFINITY;
        br->bound = INFINITY;
        br->pose_jump = fabsf(dw_pose) * fabsf(d);
        return 1;
    }

    float vL = (1.0f - w) * apL + w * bpL + dw_lo * d;
    float vR = (1.0f - w) * apR + w * bpR + dw_hi * d;
    br->jump = fabsf(vR - vL);
    br->bound = fabsf(dw_hi - dw_lo) * fabsf(d) +
                dap * fabsf(1.0f - w) + dbp * fabsf(w);
    return 1;
}

static float blend_pose(const WeightCurve* c, const Track* anim,
                        const Track* phys, float t) {
    float w = weight_eval(c, t);
    return (1.0f - w) * track_eval(anim, t) + w * track_eval(phys, t);
}

/* Closed-form p'(t); 0 when t is outside the domain or p' is not two-sided
   (weight corner, hold exit, or input velocity break at exactly t). */
static int blend_deriv_exact(const WeightCurve* c, const Track* anim,
                             const Track* phys, float t, float* out) {
    if (c->n < 2 || anim->n < 2 || phys->n < 2) return 0;
    if (t <= c->k[0].t || t >= c->k[c->n - 1].t) return 0;
    if (t >= track_duration(anim) || t >= track_duration(phys)) return 0;

    Break br;
    if (break_at_time(c, anim, phys, t, &br)) return 0;

    int s = -1;
    for (int i = 0; i < weight_segs(c); i++)
        if (t >= c->k[i].t && t < c->k[i + 1].t) { s = i; break; }
    if (s < 0) return 0;
    float w = weight_eval_seg(c, s, t);
    float dw = weight_slope_seg(c, s, t);

    int ia = track_seg_right(anim, t);
    int ip = track_seg_right(phys, t);
    float ua = (t - (float)ia * anim->dt) / anim->dt;
    float up = (t - (float)ip * phys->dt) / phys->dt;
    if (ua < 0.0f) ua = 0.0f;
    if (ua > 1.0f) ua = 1.0f;
    if (up < 0.0f) up = 0.0f;
    if (up > 1.0f) up = 1.0f;

    float a = anim->pos[ia] + (anim->pos[ia + 1] - anim->pos[ia]) * ua;
    float b = phys->pos[ip] + (phys->pos[ip + 1] - phys->pos[ip]) * up;
    *out = (1.0f - w) * track_slope(anim, ia) + w * track_slope(phys, ip) + dw * (b - a);
    return 1;
}

/* Central difference of the blended pose: the velocity the runtime sees. */
static float blend_deriv_numeric(const WeightCurve* c, const Track* anim,
                                 const Track* phys, float t) {
    const float h = BLEND_ND_H;
    return (blend_pose(c, anim, phys, t + h) -
            blend_pose(c, anim, phys, t - h)) / (2.0f * h);
}

/* ------------------------------------------------------------- reports */

static int time_exists(const Break* out, int n, float t) {
    for (int i = 0; i < n; i++)
        if (fabsf(out[i].t - t) <= 1e-4f) return 1;
    return 0;
}

/* All interior times where p' jumps or does not exist, sorted. Sources are
   weight knots and interior sample times of both inputs. */
static int blend_breaks(const WeightCurve* c, const Track* anim, const Track* phys,
                        Break* out, int cap) {
    int n = 0;
    float t0, t1;
    if (c->n < 2) return 0;
    t0 = c->k[0].t;
    t1 = c->k[c->n - 1].t;

    for (int k = 1; k <= c->n - 2 && n < cap; k++) {
        Break br;
        if (break_at_time(c, anim, phys, c->k[k].t, &br)) out[n++] = br;
    }
    for (int side = 0; side < 2 && n < cap; side++) {
        const Track* tr = side == 0 ? anim : phys;
        for (int i = 1; i < tr->n - 1 && n < cap; i++) {
            float t = tr->dt * (float)i;
            if (t <= t0 + 1e-4f || t >= t1 - 1e-4f) continue;
            if (time_exists(out, n, t)) continue;
            Break br;
            if (break_at_time(c, anim, phys, t, &br)) out[n++] = br;
        }
    }
    for (int i = 1; i < n; i++) {          /* insertion sort by time */
        Break b = out[i];
        int j = i - 1;
        while (j >= 0 && out[j].t > b.t) { out[j + 1] = out[j]; j--; }
        out[j + 1] = b;
    }
    return n;
}

static int blend_direction_warnings(const WeightCurve* c, const Track* anim,
                                    const Track* phys, float thresh_deg) {
    int count = 0;
    if (c->n < 2) return 0;
    float t0 = c->k[0].t, t1 = c->k[c->n - 1].t;
    float dur = track_duration(anim) < track_duration(phys)
              ? track_duration(anim) : track_duration(phys);
    if (t1 > dur) t1 = dur;
    for (float t = t0; t <= t1 + 1e-4f; t += anim->dt)
        if (blend_direction_deg(anim, phys, t) > thresh_deg) count++;
    return count;
}

static Report blend_report(const WeightCurve* c, const Track* anim, const Track* phys,
                           float thresh_deg) {
    Break br[BLEND_MAX_BREAKS];
    Report r;
    memset(&r, 0, sizeof(r));
    int n = blend_breaks(c, anim, phys, br, BLEND_MAX_BREAKS);
    r.n_breaks = n;
    for (int i = 0; i < n; i++) {
        if (br[i].flags & BB_WEIGHT_NOT_DIFF) {
            r.weight_corners++;
            if (br[i].jump > r.max_jump) r.max_jump = br[i].jump;
            if (br[i].bound > r.max_bound) r.max_bound = br[i].bound;
        }
        if (br[i].flags & BB_WEIGHT_DISCONT) r.weight_undefined++;
        if (br[i].flags & (BB_INPUT_ANIM | BB_INPUT_PHYS)) r.input_breaks++;
    }
    r.direction_warnings = blend_direction_warnings(c, anim, phys, thresh_deg);
    return r;
}

/* Maximal spans on which p' exists and is continuous (the complement of the
   break set). w_min/w_max and w_in_range are the additional conditions for
   the weight to be a real blend rather than an extrapolation. */
static int blend_intervals(const WeightCurve* c, const Track* anim, const Track* phys,
                           Interval* out, int cap) {
    Break br[BLEND_MAX_BREAKS];
    if (c->n < 2) return 0;
    int nb = blend_breaks(c, anim, phys, br, BLEND_MAX_BREAKS);

    float cuts[BLEND_MAX_BREAKS + 2];
    int nc = 0;
    cuts[nc++] = c->k[0].t;
    for (int i = 0; i < nb; i++) cuts[nc++] = br[i].t;
    cuts[nc++] = c->k[c->n - 1].t;

    int n = 0;
    for (int i = 0; i + 1 < nc && n < cap; i++) {
        float a = cuts[i], b = cuts[i + 1];
        if (b - a <= 1e-4f) continue;
        Interval iv;
        iv.t0 = a;
        iv.t1 = b;
        iv.w_min = weight_eval(c, a);
        iv.w_max = iv.w_min;
        float wb = weight_eval(c, b);
        if (wb < iv.w_min) iv.w_min = wb;
        if (wb > iv.w_max) iv.w_max = wb;
        for (int k = 0; k < c->n; k++) {
            if (c->k[k].t < a - 1e-4f || c->k[k].t > b + 1e-4f) continue;
            if (c->k[k].w < iv.w_min) iv.w_min = c->k[k].w;
            if (c->k[k].w > iv.w_max) iv.w_max = c->k[k].w;
        }
        iv.w_in_range = iv.w_min >= -BLEND_EPS && iv.w_max <= 1.0f + BLEND_EPS;
        out[n++] = iv;
    }
    return n;
}

/* ------------------------------------------------- numeric cross-check */

static int near_break(const Break* br, int nb, float t, float margin) {
    for (int i = 0; i < nb; i++)
        if (fabsf(br[i].t - t) < margin) return 1;
    return 0;
}

/* A C1 weight join still kinks w'', and a central-difference stencil that
   straddles it picks up O(h * |dw''| * |d|) error. Those samples are not
   evidence against the closed form, so the comparison skips them too. */
static int near_weight_kink(const WeightCurve* c, float t, float margin) {
    for (int k = 1; k <= c->n - 2; k++)
        if (fabsf(c->k[k].t - t) < margin) return 1;
    return 0;
}

/* Exact vs central-difference velocity, sampled uniformly. Samples within
   4h of a break or a weight knot (C1 or not) are skipped because the numeric
   stencil smears across the jump or the second-derivative kink; everything
   else must agree. */
static Check blend_compare(const WeightCurve* c, const Track* anim, const Track* phys,
                           float t0, float t1, int nsamp) {
    Check ck;
    memset(&ck, 0, sizeof(ck));
    Break br[BLEND_MAX_BREAKS];
    int nb = blend_breaks(c, anim, phys, br, BLEND_MAX_BREAKS);
    float margin = 4.0f * BLEND_ND_H;
    float durA = track_duration(anim), durP = track_duration(phys);
    float k0 = c->k[0].t, k1 = c->k[c->n - 1].t;
    if (nsamp < 2) return ck;
    if (nsamp > BLEND_MAX_KEYS) nsamp = BLEND_MAX_KEYS;
    ck.n = nsamp;

    for (int i = 0; i < nsamp; i++) {
        float t = t0 + (t1 - t0) * (float)i / (float)(nsamp - 1);
        ck.t[i] = t;
        int skip = near_break(br, nb, t, margin) || near_weight_kink(c, t, margin);
        if (t < k0 + margin || t > k1 - margin) skip = 1;
        if (t + margin >= durA || t + margin >= durP) skip = 1;
        float ex;
        if (!skip && !blend_deriv_exact(c, anim, phys, t, &ex)) skip = 1;
        if (skip) {
            ck.skipped++;
            ck.used[i] = 0;
            continue;
        }
        ck.used[i] = 1;
        ck.compared++;
        ck.exact[i] = ex;
        ck.numeric[i] = blend_deriv_numeric(c, anim, phys, t);
        float err = fabsf(ck.exact[i] - ck.numeric[i]);
        if (err > ck.max_abs_err) { ck.max_abs_err = err; ck.t_worst = t; }
    }
    return ck;
}

/* Numeric one-frame jump across a break: sample the numeric velocity at
   probe distance on each side and compare with the analytic jump/bound. */
static StepCheck blend_step_check(const WeightCurve* c, const Track* anim,
                                  const Track* phys, const Break* br, float probe) {
    StepCheck sc;
    memset(&sc, 0, sizeof(sc));
    sc.jump = br->jump;
    sc.bound = br->bound;
    sc.undefined = (br->flags & BB_WEIGHT_DISCONT) != 0;
    float vl = blend_deriv_numeric(c, anim, phys, br->t - probe);
    float vr = blend_deriv_numeric(c, anim, phys, br->t + probe);
    sc.measured = fabsf(vr - vl);
    if (sc.undefined) {
        sc.ok = isinf(sc.measured) || sc.measured > 1.0f;  /* clearly impulsive */
    } else {
        float slack = 0.5f;
        sc.ok = sc.measured >= sc.jump - slack && sc.measured <= sc.bound + slack;
    }
    return sc;
}

/* Gate the blend at one time: a non-differentiable weight, a weight outside
   [0,1], or a direction disagreement is reported (ok = 0) instead of being
   blended silently. Input velocity breaks are reported by blend_report but
   do not block evaluation, since p is still well defined there. */
static Gate blend_gate(const WeightCurve* c, const Track* anim, const Track* phys,
                       float t, float thresh_deg) {
    Gate g;
    memset(&g, 0, sizeof(g));
    if (c->n < 2 || t < c->k[0].t || t > c->k[c->n - 1].t) {
        g.flags |= BB_OUTSIDE;
        g.ok = 0;
        return g;
    }
    g.pose = blend_pose(c, anim, phys, t);
    g.w = weight_eval(c, t);
    if (g.w < -BLEND_EPS || g.w > 1.0f + BLEND_EPS) g.flags |= BB_WEIGHT_RANGE;

    float dw_lo, dw_hi, dw_pose;
    unsigned wf = weight_break_at(c, t, &dw_lo, &dw_hi, &dw_pose);
    g.flags |= wf;

    g.angle_deg = blend_direction_deg(anim, phys, t);
    if (g.angle_deg > thresh_deg) g.flags |= BB_DIRECTION;

    g.ok = (g.flags == 0) ? 1 : 0;
    return g;
}

/* --------------------------------------------------------------- tests */

static Track make_track(float dt, int n, const float* pos, float dx, float dy) {
    Track tr;
    memset(&tr, 0, sizeof(tr));
    tr.dt = dt;
    tr.n = n;
    for (int i = 0; i < n && i < BLEND_MAX_KEYS; i++) {
        tr.pos[i] = pos[i];
        tr.dx[i] = dx;
        tr.dy[i] = dy;
    }
    return tr;
}

int main(void) {
    const float ramp[5] = { 0, 1, 2, 3, 4 };        /* anim: +4/s */
    const float fast[5] = { 10, 12, 14, 16, 18 };   /* phys: +8/s */
    Track anim = make_track(0.25f, 5, ramp, 1, 0);
    Track phys = make_track(0.25f, 5, fast, 1, 0);
    Break br[BLEND_MAX_BREAKS];

    /* --- closed-form derivative and numeric agreement --- */
    WeightCurve lin;
    weight_init(&lin);
    check("add_knot_0", weight_add(&lin, 0.0f, 0.0f, WSEG_LINEAR) == 0);
    check("add_knot_1", weight_add(&lin, 1.0f, 1.0f, WSEG_LINEAR) == 1);

    float v = 0.0f;
    int ok = blend_deriv_exact(&lin, &anim, &phys, 0.5f, &v);
    /* 0.5*4 + 0.5*8 + 1*(14-2) = 18 */
    check("deriv_closed_form", ok && fabsf(v - 18.0f) < 1e-3f);
    check("deriv_numeric_agrees",
          fabsf(blend_deriv_numeric(&lin, &anim, &phys, 0.5f) - 18.0f) < 0.05f);

    Check ck = blend_compare(&lin, &anim, &phys, 0.1f, 0.9f, 9);
    check("compare_samples", ck.compared >= 7 && ck.skipped <= 2);
    check("compare_max_err", ck.max_abs_err < 0.05f);

    /* --- linear|linear corner: velocity jump is reported and bounded --- */
    WeightCurve corner;
    weight_init(&corner);
    weight_add(&corner, 0.0f, 0.0f, WSEG_LINEAR);
    weight_add(&corner, 0.5f, 0.2f, WSEG_LINEAR);
    weight_add(&corner, 1.0f, 1.0f, WSEG_LINEAR);
    int nb = blend_breaks(&corner, &anim, &phys, br, BLEND_MAX_BREAKS);
    check("corner_count", nb == 1);
    check("corner_position", nb == 1 && fabsf(br[0].t - 0.5f) < 1e-5f);
    check("corner_flagged", nb == 1 && (br[0].flags & BB_WEIGHT_NOT_DIFF));
    check("corner_one_sided", nb == 1 &&
          fabsf(br[0].dw_lo - 0.4f) < 1e-4f && fabsf(br[0].dw_hi - 1.6f) < 1e-4f);
    /* 1.2 * |14-2| = 14.4 */
    check("corner_jump", nb == 1 && fabsf(br[0].jump - 14.4f) < 0.01f);
    check("corner_bound_at_jump", nb == 1 && br[0].bound >= br[0].jump - 1e-4f);
    StepCheck sc = nb == 1 ? blend_step_check(&corner, &anim, &phys, &br[0], 0.01f)
                           : (StepCheck){ 0 };
    check("corner_step_measurable", nb == 1 && fabsf(sc.measured - 14.4f) < 0.3f);
    check("corner_step_within_bound", nb == 1 && sc.ok);
    check("corner_gate_blocks", blend_gate(&corner, &anim, &phys, 0.5f, 45.0f).ok == 0);
    float cv;
    check("corner_deriv_undefined", blend_deriv_exact(&corner, &anim, &phys, 0.5f, &cv) == 0);
    Interval civ[BLEND_MAX_BREAKS];
    int cni = blend_intervals(&corner, &anim, &phys, civ, BLEND_MAX_BREAKS);
    check("corner_splits_intervals", cni == 2 &&
          fabsf(civ[0].t0) < 1e-5f && fabsf(civ[0].t1 - 0.5f) < 1e-5f &&
          fabsf(civ[1].t0 - 0.5f) < 1e-5f && fabsf(civ[1].t1 - 1.0f) < 1e-5f);
    check("corner_intervals_in_range", cni == 2 && civ[0].w_in_range && civ[1].w_in_range);

    /* --- smooth|smooth: w' = 0 at the join, one C1 interval --- */
    WeightCurve smooth;
    weight_init(&smooth);
    weight_add(&smooth, 0.0f, 0.0f, WSEG_SMOOTH);
    weight_add(&smooth, 0.5f, 0.5f, WSEG_SMOOTH);
    weight_add(&smooth, 1.0f, 1.0f, WSEG_SMOOTH);
    check("smooth_no_break", blend_breaks(&smooth, &anim, &phys, br, BLEND_MAX_BREAKS) == 0);
    Interval iv[BLEND_MAX_BREAKS];
    int ni = blend_intervals(&smooth, &anim, &phys, iv, BLEND_MAX_BREAKS);
    check("smooth_single_interval", ni == 1);
    check("smooth_interval_covers", ni == 1 &&
          fabsf(iv[0].t0) < 1e-5f && fabsf(iv[0].t1 - 1.0f) < 1e-5f);
    check("smooth_weight_range", ni == 1 &&
          iv[0].w_in_range && fabsf(iv[0].w_min) < 1e-5f && fabsf(iv[0].w_max - 1.0f) < 1e-5f);
    check("smooth_gate_allows", blend_gate(&smooth, &anim, &phys, 0.5f, 45.0f).ok == 1);
    check("smooth_slope_at_join",
          fabsf(weight_slope_seg(&smooth, 0, 0.5f)) < 1e-6f &&
          fabsf(weight_slope_seg(&smooth, 1, 0.5f)) < 1e-6f);
    ck = blend_compare(&smooth, &anim, &phys, 0.05f, 0.95f, 19);
    check("smooth_compare_err", ck.compared >= 15 && ck.max_abs_err < 0.05f);

    /* --- hold: unbounded velocity, finite pose jump, gated --- */
    WeightCurve hold;
    weight_init(&hold);
    weight_add(&hold, 0.0f, 0.0f, WSEG_LINEAR);
    weight_add(&hold, 0.25f, 0.2f, WSEG_HOLD);
    weight_add(&hold, 0.75f, 0.8f, WSEG_LINEAR);
    weight_add(&hold, 1.0f, 1.0f, WSEG_LINEAR);
    nb = blend_breaks(&hold, &anim, &phys, br, BLEND_MAX_BREAKS);
    check("hold_two_breaks", nb == 2);
    int hold_exit = -1;
    for (int i = 0; i < nb; i++)
        if (br[i].flags & BB_WEIGHT_DISCONT) hold_exit = i;
    check("hold_undefined", hold_exit >= 0);
    check("hold_undefined_at_exit", hold_exit >= 0 && fabsf(br[hold_exit].t - 0.75f) < 1e-5f);
    check("hold_jump_infinite", hold_exit >= 0 && isinf(br[hold_exit].jump));
    /* |0.8-0.2| * |16-3| = 7.8 */
    check("hold_pose_jump", hold_exit >= 0 && fabsf(br[hold_exit].pose_jump - 7.8f) < 0.01f);
    check("hold_enter_is_corner", nb == 2 && (br[0].t < br[1].t ? (br[0].flags & BB_WEIGHT_NOT_DIFF)
                                                                 : (br[1].flags & BB_WEIGHT_NOT_DIFF)));
    check("hold_gate_exit_blocks",
          blend_gate(&hold, &anim, &phys, 0.75f, 45.0f).ok == 0);
    check("hold_gate_inside_allows",
          blend_gate(&hold, &anim, &phys, 0.5f, 45.0f).ok == 1);
    Report hr = blend_report(&hold, &anim, &phys, 45.0f);
    check("hold_report_undefined", hr.weight_undefined == 1 && hr.weight_corners == 1);

    /* --- input-only break: no weight corner, still split and flagged --- */
    const float stepped[4] = { 0, 1, 3, 4 };        /* slopes 1, 2, 1 */
    const float zero4[4] = { 0, 0, 0, 0 };
    Track anim_s = make_track(1.0f, 4, stepped, 1, 0);
    Track phys_z = make_track(1.0f, 4, zero4, 1, 0);
    WeightCurve wsm;
    weight_init(&wsm);
    weight_add(&wsm, 0.0f, 0.0f, WSEG_SMOOTH);
    weight_add(&wsm, 2.0f, 1.0f, WSEG_SMOOTH);
    nb = blend_breaks(&wsm, &anim_s, &phys_z, br, BLEND_MAX_BREAKS);
    check("input_break_found", nb == 1 && fabsf(br[0].t - 1.0f) < 1e-5f);
    check("input_break_flag", nb == 1 && br[0].flags == BB_INPUT_ANIM);
    /* w(1) = 0.5, anim' 1 -> 2: (1-0.5)*1 = 0.5 */
    check("input_break_jump", nb == 1 && fabsf(br[0].jump - 0.5f) < 1e-4f);
    check("input_break_bound", nb == 1 && fabsf(br[0].bound - 0.5f) < 1e-4f);
    ni = blend_intervals(&wsm, &anim_s, &phys_z, iv, BLEND_MAX_BREAKS);
    check("input_break_splits", ni == 2 &&
          fabsf(iv[0].t1 - 1.0f) < 1e-5f && fabsf(iv[1].t0 - 1.0f) < 1e-5f);
    StepCheck isc = nb == 1 ? blend_step_check(&wsm, &anim_s, &phys_z, &br[0], 0.01f)
                            : (StepCheck){ 0 };
    check("input_step_within_bound", nb == 1 && isc.ok);
    /* the pose is continuous at an input-only break, so single-frame
       evaluation is still allowed; only the velocity is split */
    check("input_break_gate_allows",
          blend_gate(&wsm, &anim_s, &phys_z, 1.0f, 45.0f).ok == 1);
    check("input_break_deriv_undefined",
          blend_deriv_exact(&wsm, &anim_s, &phys_z, 1.0f, &v) == 0);

    /* --- bound is strict when weight and input jumps cancel --- */
    WeightCurve mix;
    weight_init(&mix);
    weight_add(&mix, 0.0f, 0.0f, WSEG_LINEAR);
    weight_add(&mix, 1.0f, 0.25f, WSEG_LINEAR);
    weight_add(&mix, 2.0f, 0.75f, WSEG_LINEAR);
    nb = blend_breaks(&mix, &anim_s, &phys_z, br, BLEND_MAX_BREAKS);
    check("mixed_break_count", nb == 1);
    check("mixed_flags", nb == 1 &&
          (br[0].flags & BB_WEIGHT_NOT_DIFF) && (br[0].flags & BB_INPUT_ANIM));
    /* jump = 0.5, bound = 1.0: the terms partly cancel */
    check("mixed_jump", nb == 1 && fabsf(br[0].jump - 0.5f) < 1e-4f);
    check("mixed_bound_strict", nb == 1 && fabsf(br[0].bound - 1.0f) < 1e-4f &&
          br[0].bound > br[0].jump + 0.4f);
    StepCheck msc = nb == 1 ? blend_step_check(&mix, &anim_s, &phys_z, &br[0], 0.01f)
                            : (StepCheck){ 0 };
    check("mixed_step_in_bound", nb == 1 && msc.ok && msc.measured <= msc.bound);

    /* --- direction disagreement: report, do not blend --- */
    Track anim_x = make_track(0.25f, 5, ramp, 1, 0);   /* 0 deg   */
    Track phys_y = make_track(0.25f, 5, fast, 0, 1);   /* 90 deg  */
    check("direction_angle_90",
          fabsf(blend_direction_deg(&anim_x, &phys_y, 0.5f) - 90.0f) < 1e-4f);
    Gate dg = blend_gate(&smooth, &anim_x, &phys_y, 0.5f, 45.0f);
    check("direction_gate_blocks", dg.ok == 0 && (dg.flags & BB_DIRECTION));
    check("direction_gate_reports_angle", fabsf(dg.angle_deg - 90.0f) < 1e-4f);
    check("direction_threshold_passes",
          blend_gate(&smooth, &anim_x, &phys_y, 0.5f, 135.0f).ok == 1);

    float c170 = cosf(170.0f / BLEND_R2D), s170 = sinf(170.0f / BLEND_R2D);
    float cm170 = cosf(-170.0f / BLEND_R2D), sm170 = sinf(-170.0f / BLEND_R2D);
    Track anim_w = make_track(0.25f, 5, ramp, c170, s170);
    Track phys_w = make_track(0.25f, 5, fast, cm170, sm170);
    check("direction_wraps", fabsf(blend_direction_deg(&anim_w, &phys_w, 0.5f) - 20.0f) < 1e-3f);
    Track phys_still = make_track(0.25f, 5, fast, 0, 0);
    check("direction_zero_is_neutral",
          fabsf(blend_direction_deg(&anim_x, &phys_still, 0.5f)) < 1e-6f);

    /* --- weight range and domain gate --- */
    WeightCurve over;
    weight_init(&over);
    weight_add(&over, 0.0f, 0.0f, WSEG_SMOOTH);
    weight_add(&over, 0.5f, 1.5f, WSEG_SMOOTH);
    weight_add(&over, 1.0f, 0.0f, WSEG_SMOOTH);
    Gate og = blend_gate(&over, &anim, &phys, 0.5f, 45.0f);
    check("range_gate_blocks", og.ok == 0 && (og.flags & BB_WEIGHT_RANGE));
    check("range_not_a_corner", (og.flags & BB_WEIGHT_NOT_DIFF) == 0);
    Interval oiv[BLEND_MAX_BREAKS];
    ni = blend_intervals(&over, &anim, &phys, oiv, BLEND_MAX_BREAKS);
    check("range_interval_flag", ni == 1 && oiv[0].w_in_range == 0 &&
          fabsf(oiv[0].w_max - 1.5f) < 1e-4f);
    check("outside_gate_blocks",
          blend_gate(&smooth, &anim, &phys, 1.5f, 45.0f).flags & BB_OUTSIDE);

    /* --- rejections and capacity --- */
    check("reject_non_increasing", weight_add(&lin, 1.0f, 0.5f, WSEG_LINEAR) == -1);
    WeightCurve fullw;
    weight_init(&fullw);
    int last = 0;
    for (int i = 0; i < BLEND_MAX_KEYS; i++) last = weight_add(&fullw, (float)i, 0.5f, WSEG_LINEAR);
    check("knot_capacity", last == BLEND_MAX_KEYS - 1 &&
          weight_add(&fullw, 99.0f, 0.0f, WSEG_LINEAR) == -1);

    /* --- table break boundaries are not interior breaks --- */
    const float two[3] = { 0, 1, 3 };
    Track t2 = make_track(1.0f, 3, two, 1, 0);
    check("break_at_left_edge", track_break_at(&t2, 0.0f) == -1);
    check("break_at_interior", track_break_at(&t2, 1.0f) == 1);
    check("break_at_right_edge", track_break_at(&t2, 2.0f) == -1);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
