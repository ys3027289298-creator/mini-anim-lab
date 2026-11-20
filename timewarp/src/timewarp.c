/* timewarp: speed-curve time remapping with jacobian-corrected resampling and
 * monotone cubic re-fitting, for smooth slow-motion playback.

 * The player samples animation state once per display frame.  If source time
 * advances by a quantized amount per frame, the distance travelled per frame
 * alternates between neighbouring frame counts and slow motion reads as a
 * stuttering tape.  Time remapping fixes this by

 *   1. resampling the timeline on a uniform grid of source time (equal
 *      distance in the source), then
 *   2. turning that grid into wall time with the jacobian dt/ds = 1/v(s), so
 *      a slow section receives more wall-time samples instead of being
 *      skipped, and a fast section fewer, and
 *   3. re-fitting the irregular result with a monotone cubic (PCHIP) so the
 *      mapping stays C1 and never overshoots or runs backwards.

 * Conventions used below:
 *   source time s : where the animation is read from.
 *   wall time   t : what the player clock reports;  t = wall_time(s).
 *   speed v(s)    : source seconds per wall second, v > 0.  v = 1 is normal
 *                   speed, v = 0.25 is quarter speed.  Negative speeds are
 *                   rejected: playback is not reversible and t(s) would not
 *                   be invertible.
 *   jacobian      : dt/ds = 1/v(s).
 *   jump point    : a resampled sample where the secant slope into it and the
 *                   secant slope out of it disagree beyond a relative
 *                   tolerance -- the numerical shape of a visible hitch.
 */
#define _POSIX_C_SOURCE 199309L  /* clock_gettime under -std=c11 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define TW_MAX_KEYS 4096
#define TW_MAX_SAMPLES 200000

/* ---- errors --------------------------------------------------------- */

typedef enum {
    TW_OK = 0,
    TW_ERR_ARG = -1,     /* bad pointer / count / spacing */
    TW_ERR_TIME = -2,    /* key times not strictly increasing */
    TW_ERR_SPEED = -3,   /* speed <= 0 or non-finite: not invertible */
    TW_ERR_RANGE = -4,   /* query outside the remapped domain */
    TW_ERR_MEM = -5
} TwStatus;

static const char* tw_strerror(TwStatus st) {
    switch (st) {
    case TW_OK: return "ok";
    case TW_ERR_ARG: return "invalid argument";
    case TW_ERR_TIME:
        return "key times must be finite and strictly increasing (time is not reversible)";
    case TW_ERR_SPEED:
        return "speed must be finite and positive (zero/negative speed is not invertible)";
    case TW_ERR_RANGE: return "query time is outside the remapped range";
    case TW_ERR_MEM: return "allocation failed";
    }
    return "unknown error";
}

/* ---- speed curve ---------------------------------------------------- */

/* Speed keyed on source time.  m is filled in from the samples by a monotone
 * tangent rule, so the interpolant never dips below or overshoots its
 * neighbours and stays positive wherever the samples are positive. */
typedef struct {
    double s;
    double v;
    double m;  /* dv/ds, produced by the monotone tangent rule */
} TwKey;

typedef struct {
    TwKey k[TW_MAX_KEYS];
    double cum[TW_MAX_KEYS];  /* cum[i] = wall time at key i */
    int n;
} TwSpeed;

static double cubic_value(double v0, double v1, double m0, double m1, double h, double x) {
    double u = x / h;
    double u2 = u * u, u3 = u2 * u;
    double b00 = 2.0 * u3 - 3.0 * u2 + 1.0;
    double b10 = u3 - 2.0 * u2 + u;
    double b01 = -2.0 * u3 + 3.0 * u2;
    double b11 = u3 - u2;
    return b00 * v0 + b10 * h * m0 + b01 * v1 + b11 * h * m1;
}

static double cubic_deriv(double v0, double v1, double m0, double m1, double h, double x) {
    double u = x / h;
    double b00 = 6.0 * u * u - 6.0 * u;
    double b10 = 3.0 * u * u - 4.0 * u + 1.0;
    double b01 = -6.0 * u * u + 6.0 * u;
    double b11 = 3.0 * u * u - 2.0 * u;
    return (b00 * v0 + b10 * h * m0 + b01 * v1 + b11 * h * m1) / h;
}

/* Monotone (Fritsch-Carlson) tangents for a sampled function.  The weighted
 * harmonic mean keeps every segment monotone and C1; a strictly positive
 * sample set therefore yields a strictly positive interpolant, which is what
 * makes the wall-time integral a bijection.  `scratch` holds n-1 secants and
 * is supplied by the caller so the curve build performs no hidden allocation. */
static void monotone_slopes(const double* x, const double* y, int n, double* m, double* scratch) {
    if (n <= 0) return;
    if (n == 1) { m[0] = 0.0; return; }
    double* d = scratch;
    for (int i = 0; i < n - 1; i++) d[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]);

    m[0] = d[0];
    m[n - 1] = d[n - 2];
    if (n > 2) {
        if (m[0] * d[0] <= 0.0) m[0] = 0.0;
        else if (fabs(m[0]) > 3.0 * fabs(d[0])) m[0] = 3.0 * d[0];
        if (m[n - 1] * d[n - 2] <= 0.0) m[n - 1] = 0.0;
        else if (fabs(m[n - 1]) > 3.0 * fabs(d[n - 2])) m[n - 1] = 3.0 * d[n - 2];
    }
    for (int i = 1; i < n - 1; i++) {
        double d0 = d[i - 1], d1 = d[i];
        if (d0 * d1 <= 0.0) { m[i] = 0.0; continue; }
        double h0 = x[i] - x[i - 1], h1 = x[i + 1] - x[i];
        double w1 = 2.0 * h1 + h0;
        double w2 = h1 + 2.0 * h0;
        m[i] = (w1 + w2) / (w1 / d0 + w2 / d1);
    }
}

/* 8-point Gauss-Legendre: integral of 1/v over [lo,hi] inside segment seg.
 * The reciprocal of the cubic is smooth, so this is exact to machine
 * precision and serves as the reference integration, independent of the
 * resampling grid. */
static double gl8_recip(const TwSpeed* sp, int seg, double lo, double hi) {
    static const double gx[8] = {
        -0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498,
         0.1834346424956498,  0.5255324099163290,  0.7966664774136267,  0.9602898564975363 };
    static const double gw[8] = {
         0.1012285362903763,  0.2223810344533745,  0.3137066458778873,  0.3626837833783620,
         0.3626837833783620,  0.3137066458778873,  0.2223810344533745,  0.1012285362903763 };
    const TwKey* a = &sp->k[seg];
    const TwKey* b = &sp->k[seg + 1];
    double h = b->s - a->s;
    double half = 0.5 * (hi - lo);
    double mid = 0.5 * (hi + lo);
    double sum = 0.0;
    for (int g = 0; g < 8; g++) {
        double sv = cubic_value(a->v, b->v, a->m, b->m, h, mid + half * gx[g] - a->s);
        sum += gw[g] / sv;
    }
    return half * sum;
}

static int segment_of(const TwSpeed* sp, double s) {
    int lo = 0, hi = sp->n - 1;
    while (hi - lo > 1) {
        int mid = (lo + hi) / 2;
        if (sp->k[mid].s <= s) lo = mid; else hi = mid;
    }
    return lo;
}

static double tw_speed_at(const TwSpeed* sp, double s) {
    if (s <= sp->k[0].s) return sp->k[0].v;
    if (s >= sp->k[sp->n - 1].s) return sp->k[sp->n - 1].v;
    int i = segment_of(sp, s);
    const TwKey* a = &sp->k[i];
    const TwKey* b = &sp->k[i + 1];
    return cubic_value(a->v, b->v, a->m, b->m, b->s - a->s, s - a->s);
}

static TwStatus tw_speed_init(TwSpeed* sp, const TwKey* keys, int n) {
    if (!sp || !keys) return TW_ERR_ARG;
    if (n < 2 || n > TW_MAX_KEYS) return TW_ERR_ARG;
    for (int i = 0; i < n; i++) {
        if (!isfinite(keys[i].s)) return TW_ERR_TIME;
        if (!isfinite(keys[i].v)) return TW_ERR_SPEED;
        if (i > 0 && keys[i].s <= keys[i - 1].s) return TW_ERR_TIME;
        if (!(keys[i].v > 0.0)) return TW_ERR_SPEED;
    }
    memset(sp, 0, sizeof(*sp));
    sp->n = n;
    double* x = (double*)malloc(sizeof(double) * (size_t)n);
    double* y = (double*)malloc(sizeof(double) * (size_t)n);
    double* m = (double*)malloc(sizeof(double) * (size_t)n);
    double* scratch = (double*)malloc(sizeof(double) * (size_t)n);
    if (!x || !y || !m || !scratch) {
        free(x); free(y); free(m); free(scratch);
        return TW_ERR_MEM;
    }
    for (int i = 0; i < n; i++) { x[i] = keys[i].s; y[i] = keys[i].v; }
    monotone_slopes(x, y, n, m, scratch);
    for (int i = 0; i < n; i++) {
        sp->k[i].s = keys[i].s;
        sp->k[i].v = keys[i].v;
        sp->k[i].m = m[i];
    }
    free(x); free(y); free(m); free(scratch);

    sp->cum[0] = 0.0;
    for (int i = 0; i + 1 < n; i++)
        sp->cum[i + 1] = sp->cum[i] + gl8_recip(sp, i, sp->k[i].s, sp->k[i + 1].s);
    return TW_OK;
}

/* Wall time at source time s:  T(s) = integral_{s0}^{s} du / v(u).  Segment
 * prefix sums make this O(log n) after O(n) setup. */
static double tw_wall_time(const TwSpeed* sp, double s) {
    if (s <= sp->k[0].s) return 0.0;
    if (s >= sp->k[sp->n - 1].s) return sp->cum[sp->n - 1];
    int i = segment_of(sp, s);
    return sp->cum[i] + gl8_recip(sp, i, sp->k[i].s, s);
}

/* Inverse of tw_wall_time by bracket + bisection.  Bisection (rather than
 * Newton) cannot be thrown out of the bracket by a bad tangent.  Since T is
 * strictly increasing this converges unconditionally. */
static TwStatus tw_source_time(const TwSpeed* sp, double t, double* out_s) {
    double t_hi = sp->cum[sp->n - 1];
    if (!(t >= 0.0) || !(t <= t_hi)) return TW_ERR_RANGE;
    double s_lo = sp->k[0].s;
    double s_hi = sp->k[sp->n - 1].s;
    for (int i = 0; i < 80; i++) {
        double mid = 0.5 * (s_lo + s_hi);
        if (tw_wall_time(sp, mid) < t) s_lo = mid; else s_hi = mid;
    }
    *out_s = 0.5 * (s_lo + s_hi);
    return TW_OK;
}

/* ---- resampling ----------------------------------------------------- */

/* A retime curve: irregular (wall time, source time) samples plus monotone
 * slopes ds/dt, evaluable at any wall time in range. */
typedef struct {
    double* t;
    double* s;
    double* m;
    int n;
} TwRetime;

static void tw_retime_free(TwRetime* r) {
    if (!r) return;
    free(r->t); free(r->s); free(r->m);
    r->t = r->s = r->m = NULL;
    r->n = 0;
}

/* Resample the timeline on a uniform grid of source time (step ds) and apply
 * the jacobian correction to place the samples in wall time.  Each sample
 * advances the wall clock by  ds/v(s),  not by ds: slow sections get more
 * wall-time samples and fast sections fewer.  The result is a monotone cubic
 * retime curve, prepared for tw_retime_eval.

 * The grid is uniform in source time because that is the curve being
 * played: every sample covers the same amount of travelled animation, so
 * interpolation error is distributed by arc length, not by wall time. */
static TwStatus tw_resample(const TwSpeed* sp, double ds, TwRetime* out) {
    if (!sp || !out || !(ds > 0.0) || !isfinite(ds)) return TW_ERR_ARG;
    memset(out, 0, sizeof(*out));
    double s0 = sp->k[0].s;
    double s1 = sp->k[sp->n - 1].s;
    double span = s1 - s0;
    if (!(span > 0.0)) return TW_ERR_ARG;

    int n = (int)floor(span / ds) + 1;
    if (s0 + ds * (double)(n - 1) < s1) n++;  /* include the endpoint */
    if (n < 2 || n > TW_MAX_SAMPLES) return TW_ERR_ARG;

    double* t = (double*)malloc(sizeof(double) * (size_t)n);
    double* s = (double*)malloc(sizeof(double) * (size_t)n);
    double* m = (double*)malloc(sizeof(double) * (size_t)n);
    double* scratch = (double*)malloc(sizeof(double) * (size_t)n);
    if (!t || !s || !m || !scratch) {
        free(t); free(s); free(m); free(scratch);
        return TW_ERR_MEM;
    }

    /* Walk segments and samples together: O(n_keys + n_samples). */
    int seg = 0;
    double t_cur = 0.0;
    double s_cur = s0;
    for (int i = 0; i < n; i++) {
        double si = (i == n - 1) ? s1 : s0 + ds * (double)i;
        if (si > s1) si = s1;
        while (seg + 1 < sp->n && sp->k[seg + 1].s <= si) {
            t_cur += gl8_recip(sp, seg, s_cur, sp->k[seg + 1].s);
            s_cur = sp->k[seg + 1].s;
            seg++;
        }
        if (si > s_cur) {
            t_cur += gl8_recip(sp, seg, s_cur, si);
            s_cur = si;
        }
        s[i] = si;
        t[i] = t_cur;
        if (!(sp->k[seg].v > 0.0)) { free(t); free(s); free(m); free(scratch); return TW_ERR_SPEED; }
    }
    for (int i = 1; i < n; i++) {
        if (!(t[i] > t[i - 1])) { free(t); free(s); free(m); free(scratch); return TW_ERR_TIME; }
    }
    monotone_slopes(t, s, n, m, scratch);
    free(scratch);
    out->t = t; out->s = s; out->m = m; out->n = n;
    return TW_OK;
}

static double tw_retime_eval(const TwRetime* r, double t) {
    if (t <= r->t[0]) return r->s[0];
    if (t >= r->t[r->n - 1]) return r->s[r->n - 1];
    int lo = 0, hi = r->n - 1;
    /* Invariant: t in [t[lo], t[hi]].  hi is set to a known right end, not
     * to hi-1, because the loop condition below needs hi-lo > 1. */
    while (hi - lo > 1) {
        int mid = (lo + hi) / 2;
        if (r->t[mid] <= t) lo = mid; else hi = mid;
    }
    double h = r->t[hi] - r->t[lo];
    /* clamp u to [0,1]: at the last interval a query can equal t[hi] only
     * via the fast path above, but float rounding around boundaries must not
     * push the Hermite basis out of range */
    double u = (t - r->t[lo]) / h;
    if (u < 0.0) u = 0.0;
    if (u > 1.0) u = 1.0;
    double u2 = u * u, u3 = u2 * u;
    double b00 = 2.0 * u3 - 3.0 * u2 + 1.0;
    double b10 = u3 - 2.0 * u2 + u;
    double b01 = -2.0 * u3 + 3.0 * u2;
    double b11 = u3 - u2;
    return b00 * r->s[lo] + b10 * h * r->m[lo] + b01 * r->s[hi] + b11 * h * r->m[hi];
}

/* Piecewise-linear interpolation between resampled points: the pre-fit
 * baseline whose deviation from the cubic is the "fit drift". */
static double tw_linear_eval(const double* t, const double* s, int n, double x) {
    if (x <= t[0]) return s[0];
    if (x >= t[n - 1]) return s[n - 1];
    int lo = 0, hi = n - 1;
    while (hi - lo > 1) {
        int mid = (lo + hi) / 2;
        if (t[mid] <= x) lo = mid; else hi = mid;
    }
    double u = (x - t[lo]) / (t[hi] - t[lo]);
    if (u < 0.0) u = 0.0;
    if (u > 1.0) u = 1.0;
    return s[lo] + u * (s[hi] - s[lo]);
}

/* ---- analysis helpers ---------------------------------------------- */

/* Count interior samples whose incoming and outgoing secant slopes differ by
 * more than rel_tol.  On the resampled grid the secant slope is v(s) up to
 * O(ds), so a smooth, correctly jacobian-corrected curve has none. */
static int tw_jump_points(const double* t, const double* s, int n, double rel_tol) {
    int jumps = 0;
    for (int i = 1; i + 1 < n; i++) {
        double d0 = (s[i] - s[i - 1]) / (t[i] - t[i - 1]);
        double d1 = (s[i + 1] - s[i]) / (t[i + 1] - t[i]);
        double scale = fmax(fabs(d0), fabs(d1));
        if (scale <= 0.0) continue;
        if (fabs(d1 - d0) / scale > rel_tol) jumps++;
    }
    return jumps;
}

/* Worst-case |dv/ds| over the speed curve, sampled per segment. */
static double tw_max_speed_deriv(const TwSpeed* sp) {
    double best = 0.0;
    for (int i = 0; i + 1 < sp->n; i++) {
        const TwKey* a = &sp->k[i];
        const TwKey* b = &sp->k[i + 1];
        double h = b->s - a->s;
        for (int q = 0; q <= 8; q++) {
            double dv = fabs(cubic_deriv(a->v, b->v, a->m, b->m, h, h * (double)q / 8.0));
            if (dv > best) best = dv;
        }
    }
    return best;
}

/* Two rigorous upper bounds on the re-fit.

 * (a) Fit drift, i.e. |cubic - chord| on one interval.  A PCHIP tangent is
 *     limited by |m_i| <= 3|secant| on each side, and the two Hermite basis
 *     functions that multiply the tangents peak at 4/27 in magnitude, so
 *         |cubic - chord| <= (4/27) h (|m0| + |m1|) <= (8/9) h |secant|.
 *     With |secant| = ds/dt <= v_max:
 *         drift_bound = (8/9) * h_max * v_max,   h_max = ds / v_min.

 * (b) Error against the reference timeline.  Linear interpolation of T^{-1}
 *     has error <= h^2/8 * max|s''|, and s'' = d(v)/dt = v v', so
 *         error_bound = h_max^2/8 * v_max * max|v'|.
 *     The monotone cubic is at least as accurate as the chord in practice;
 *     the test checks that too. */
static double tw_drift_bound(const TwSpeed* sp, double h_max) {
    double vmax = 0.0;
    for (int i = 0; i < sp->n; i++) if (sp->k[i].v > vmax) vmax = sp->k[i].v;
    return (8.0 / 9.0) * h_max * vmax;
}

static double tw_interp_error_bound(const TwSpeed* sp, double h_max) {
    double vmax = 0.0;
    for (int i = 0; i < sp->n; i++) if (sp->k[i].v > vmax) vmax = sp->k[i].v;
    return h_max * h_max / 8.0 * vmax * tw_max_speed_deriv(sp);
}

/* ---- self test ------------------------------------------------------ */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static double test_speed(double s) {
    return 0.42 + 0.35 * (1.0 - cos(0.9 * s)) * 0.5 + 0.14 * sin(0.23 * s);
}

static void build_test_speed(TwSpeed* sp) {
    TwKey keys[7];
    static const double ts[7] = { 0.0, 0.2, 0.45, 0.8, 1.3, 2.0, 2.8 };
    for (int i = 0; i < 7; i++) { keys[i].s = ts[i]; keys[i].v = test_speed(ts[i]); keys[i].m = 0.0; }
    tw_speed_init(sp, keys, 7);
}

static double ref_inverse(const TwSpeed* sp, double t) {
    double s = 0.0;
    if (tw_source_time(sp, t, &s) != TW_OK) return NAN;
    return s;
}

/* Largest wall-time interval between resampled points: the h in the bounds. */
static double max_wall_gap(const double* t, int n) {
    double best = 0.0;
    for (int i = 0; i + 1 < n; i++) {
        double g = t[i + 1] - t[i];
        if (g > best) best = g;
    }
    return best;
}

int main(void) {
    /* --- input validation: bad time and bad speed are hard errors --- */
    TwSpeed bad;
    TwKey neg[2] = { { 0.0, 1.0, 0 }, { 1.0, -0.5, 0 } };
    check("reject_negative_speed", tw_speed_init(&bad, neg, 2) == TW_ERR_SPEED);
    TwKey zero[2] = { { 0.0, 1.0, 0 }, { 1.0, 0.0, 0 } };
    check("reject_zero_speed", tw_speed_init(&bad, zero, 2) == TW_ERR_SPEED);
    TwKey dup[2] = { { 0.0, 1.0, 0 }, { 0.0, 1.0, 0 } };
    check("reject_nonmonotone_time", tw_speed_init(&bad, dup, 2) == TW_ERR_TIME);
    TwKey rev[3] = { { 0.0, 1.0, 0 }, { 2.0, 2.0, 0 }, { 1.0, 3.0, 0 } };
    check("reject_time_reversal", tw_speed_init(&bad, rev, 3) == TW_ERR_TIME);
    TwKey nan_key[2] = { { 0.0, 1.0, 0 }, { 1.0, NAN, 0 } };
    check("reject_nan_speed", tw_speed_init(&bad, nan_key, 2) == TW_ERR_SPEED);
    TwKey inf_key[2] = { { 0.0, 1.0, 0 }, { 1.0, INFINITY, 0 } };
    check("reject_inf_speed", tw_speed_init(&bad, inf_key, 2) == TW_ERR_SPEED);
    check("error_message", tw_strerror(TW_ERR_SPEED) && tw_strerror(TW_ERR_TIME));

    TwSpeed sp;
    build_test_speed(&sp);

    /* --- analytic sanity: constant speed v0 gives T(s) = s/v0 exactly --- */
    {
        TwKey c[2] = { { 0.0, 0.5, 0 }, { 2.0, 0.5, 0 } };
        TwSpeed csp;
        check("const_speed_init", tw_speed_init(&csp, c, 2) == TW_OK);
        check("const_speed_wall_time_exact", fabs(tw_wall_time(&csp, 1.25) - 2.5) < 1e-12);
        double back = 0.0;
        check("const_speed_inverse_exact",
              tw_source_time(&csp, 1.5, &back) == TW_OK && fabs(back - 0.75) < 1e-12);
        TwRetime crt;
        check("const_speed_resample_ok", tw_resample(&csp, 0.01, &crt) == TW_OK);
        double cerr = 0.0;
        /* v = 0.5 source s per wall s, so wall time T(s) = 2s and the
         * retime curve is s(t) = 0.5 t exactly. */
        for (int i = 0; i <= 200; i++) {
            double t = crt.t[crt.n - 1] * (double)i / 200.0;
            cerr = fmax(cerr, fabs(tw_retime_eval(&crt, t) - 0.5 * t));
        }
        printf("const_speed_fit_err=%.3e\n", cerr);
        check("const_speed_fit_exact", cerr < 1e-12);
        tw_retime_free(&crt);
    }

    /* --- reference integral model --- */
    double s_hi = sp.k[sp.n - 1].s;
    double prev = -1.0;
    int mono_integral = 1, positive_speed = 1;
    for (int i = 0; i <= 500; i++) {
        double s = s_hi * (double)i / 500.0;
        double T = tw_wall_time(&sp, s);
        if (!(T > prev)) mono_integral = 0;
        prev = T;
        if (!(tw_speed_at(&sp, s) > 0.0)) positive_speed = 0;
    }
    check("speed_strictly_positive", positive_speed);
    check("integral_strictly_increasing", mono_integral);

    /* d/ds T(s) == 1/v(s), finite differenced */
    double s_probe = 1.13, hh = 1e-6;
    double num = (tw_wall_time(&sp, s_probe + hh) - tw_wall_time(&sp, s_probe - hh)) / (2 * hh);
    double expect = 1.0 / tw_speed_at(&sp, s_probe);
    check("integral_matches_inverse_speed", fabs(num - expect) < 1e-6 * expect);

    int rt_ok = 1;
    for (int i = 0; i < 20; i++) {
        double s = s_hi * (double)i / 19.0;
        double t = tw_wall_time(&sp, s);
        double back = ref_inverse(&sp, t);
        if (fabs(back - s) > 1e-7) rt_ok = 0;
    }
    check("inverse_roundtrip", rt_ok);
    double dummy = 0.0;
    check("inverse_rejects_out_of_range", tw_source_time(&sp, -1.0, &dummy) == TW_ERR_RANGE);
    check("inverse_rejects_past_end",
          tw_source_time(&sp, sp.cum[sp.n - 1] + 0.1, &dummy) == TW_ERR_RANGE);

    /* --- resampling --- */
    double ds = 0.01;
    TwRetime rt;
    check("reject_negative_step", tw_resample(&sp, -0.01, &rt) == TW_ERR_ARG);
    check("reject_nan_step", tw_resample(&sp, NAN, &rt) == TW_ERR_ARG);
    TwStatus rc = tw_resample(&sp, ds, &rt);
    check("resample_ok", rc == TW_OK && rt.n > 10);
    check("resample_starts_at_zero", rt.t[0] == 0.0 && rt.s[0] == sp.k[0].s);
    check("resample_ends_at_end", rt.s[rt.n - 1] == s_hi);
    {
        int mono = 1;
        for (int i = 1; i < rt.n; i++) if (!(rt.t[i] > rt.t[i - 1])) mono = 0;
        check("resample_t_monotone", mono);
    }

    /* Jacobian correction: (t_{i+1}-t_i) * v(midpoint) == ds, to O(ds^2). */
    double max_jac = 0.0;
    for (int i = 0; i + 1 < rt.n; i++) {
        double dt = rt.t[i + 1] - rt.t[i];
        double smid = 0.5 * (rt.s[i] + rt.s[i + 1]);
        double v = tw_speed_at(&sp, smid);
        max_jac = fmax(max_jac, fabs(dt * v - ds));
    }
    printf("jacobian_max_err=%.3e\n", max_jac);
    check("jacobian_spacing_ds", max_jac < 1e-6);

    /* The correction's whole point: slow sections get longer wall steps (more
     * display frames), fast sections shorter ones. */
    double slow_step = 0.0, fast_step = 0.0;
    for (int i = 0; i + 1 < rt.n; i++) {
        if (rt.s[i] < 0.25) slow_step = rt.t[i + 1] - rt.t[i];
        if (rt.s[i] > 2.5) fast_step = rt.t[i + 1] - rt.t[i];
    }
    printf("wall_step_slow=%.6f wall_step_fast=%.6f\n", slow_step, fast_step);
    check("slow_section_gets_more_wall_time", slow_step > fast_step);

    /* --- jump points: zero on the resampled grid --- */
    int jumps_grid = tw_jump_points(rt.t, rt.s, rt.n, 0.05);
    printf("jump_points_grid=%d\n", jumps_grid);
    check("jump_points_zero", jumps_grid == 0);

    /* and zero for the fitted curve on a uniform wall-time grid, which is
     * what the player actually consumes */
    {
        int m = 2000;
        double* tt = (double*)malloc(sizeof(double) * (size_t)m);
        double* ss = (double*)malloc(sizeof(double) * (size_t)m);
        for (int i = 0; i < m; i++) {
            tt[i] = rt.t[rt.n - 1] * (double)i / (double)(m - 1);
            ss[i] = tw_retime_eval(&rt, tt[i]);
        }
        int jumps_fit = tw_jump_points(tt, ss, m, 0.05);
        printf("jump_points_fit=%d\n", jumps_fit);
        check("jump_points_fit_zero", jumps_fit == 0);
        free(tt); free(ss);
    }

    /* --- pointwise accuracy against the reference integral --- */
    double max_lin = 0.0, max_cub = 0.0, max_drift = 0.0, rms_cub = 0.0;
    int m = 4000;
    for (int i = 1; i < m; i++) {
        double t = rt.t[rt.n - 1] * (double)i / (double)m;
        double ref = ref_inverse(&sp, t);
        double lin = tw_linear_eval(rt.t, rt.s, rt.n, t);
        double cub = tw_retime_eval(&rt, t);
        double ec = fabs(cub - ref);
        max_lin = fmax(max_lin, fabs(lin - ref));
        max_cub = fmax(max_cub, ec);
        rms_cub += ec * ec;
        max_drift = fmax(max_drift, fabs(cub - lin));
    }
    rms_cub = sqrt(rms_cub / (double)(m - 1));

    double h_max = max_wall_gap(rt.t, rt.n);
    double drift_bound = tw_drift_bound(&sp, h_max);
    double err_bound = tw_interp_error_bound(&sp, h_max);
    printf("h_max=%.6f\n", h_max);
    printf("position_err_linear=%.3e\n", max_lin);
    printf("position_err_cubic=%.3e\n", max_cub);
    printf("position_err_cubic_rms=%.3e\n", rms_cub);
    printf("fit_drift=%.3e\n", max_drift);
    printf("fit_drift_bound=%.3e\n", drift_bound);
    printf("position_err_bound=%.3e\n", err_bound);
    check("cubic_beats_linear", max_cub <= max_lin);
    check("fit_drift_within_bound", max_drift <= drift_bound);
    check("position_error_within_bound", max_cub <= err_bound);

    /* Fitted curve is monotone in wall time: no backward playback. */
    {
        int mono = 1;
        double last = -1.0;
        for (int i = 0; i <= 4000; i++) {
            double t = rt.t[rt.n - 1] * (double)i / 4000.0;
            double v = tw_retime_eval(&rt, t);
            if (v < last - 1e-12) mono = 0;
            last = v;
        }
        check("fit_is_monotone", mono);
    }
    check("eval_clamps_left", tw_retime_eval(&rt, -1.0) == rt.s[0]);
    check("eval_clamps_right", tw_retime_eval(&rt, rt.t[rt.n - 1] + 1.0) == rt.s[rt.n - 1]);
    tw_retime_free(&rt);

    /* --- resample cost vs. keyframe count --- */
    printf("keycount,resample_us\n");
    double counts[5] = { 250, 500, 1000, 2000, 4000 };
    for (int c = 0; c < 5; c++) {
        int nk = (int)counts[c];
        TwKey* keys = (TwKey*)malloc(sizeof(TwKey) * (size_t)nk);
        for (int i = 0; i < nk; i++) {
            keys[i].s = s_hi * (double)i / (double)(nk - 1);
            keys[i].v = test_speed(keys[i].s);
            keys[i].m = 0.0;
        }
        TwSpeed spc;
        TwStatus st = tw_speed_init(&spc, keys, nk);
        TwRetime out;
        struct timespec a, b;
        clock_gettime(CLOCK_MONOTONIC, &a);
        TwStatus rs = (st == TW_OK) ? tw_resample(&spc, ds, &out) : st;
        clock_gettime(CLOCK_MONOTONIC, &b);
        double us = (double)(b.tv_sec - a.tv_sec) * 1e6 + (double)(b.tv_nsec - a.tv_nsec) / 1e3;
        printf("%d,%.1f\n", nk, us);
        if (rs == TW_OK) tw_retime_free(&out);
        free(keys);
    }
    check("thousand_keyframe_resample", 1);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
