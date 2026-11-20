
/* curve: keyframe sampling with cubic hermite tangents and clamping.
 *
 * Two tangent modes:
 *   - default: caller-supplied tan_in / tan_out per key (unchanged behaviour);
 *   - monotone: Fritsch-Carlson tangents derived from the key values.
 *     Secant slopes of adjacent intervals clip the endpoint tangents
 *     (alpha^2 + beta^2 <= 9 per interval), so every cubic stays monotone
 *     on its interval and cannot overshoot. Boundary keys use one-sided
 *     differences. The derived slope is shared as tan_in == tan_out, so the
 *     piecewise curve is C1: there are no velocity discontinuities at keys.
 */
#include <stdio.h>
#include <string.h>
#include <math.h>

#define MAX_KEYS 16
#define PI_F 3.14159265358979323846f  /* M_PI is not exposed under -std=c11 */

typedef struct { float t; float v; float tan_in; float tan_out; } Key;

typedef struct {
    Key keys[MAX_KEYS];
    int count;
    int monotone;              /* set while Fritsch-Carlson mode is on */
    float mono_tan[MAX_KEYS];  /* shared in/out slope per key, mono mode */
} Curve;

static void curve_init(Curve* c) { memset(c, 0, sizeof(*c)); }

static int curve_build_monotone(Curve* c);  /* defined below, used by curve_add */

/* Timestamps must be finite and strictly increasing; duplicate or reversed
 * timestamps are rejected, as are non-finite values. */
static int curve_validate_times(const Curve* c) {
    for (int i = 0; i < c->count; i++) {
        if (!isfinite(c->keys[i].t) || !isfinite(c->keys[i].v)) return -1;
    }
    for (int i = 1; i < c->count; i++) {
        if (!(c->keys[i].t > c->keys[i - 1].t)) return -1;
    }
    return 0;
}

static int curve_add(Curve* c, float t, float v, float tin, float tout) {
    if (c->count >= MAX_KEYS) return -1;
    if (!isfinite(t) || !isfinite(v) || !isfinite(tin) || !isfinite(tout))
        return -1;
    if (c->count > 0 && t <= c->keys[c->count - 1].t) return -1;  /* strictly increasing */
    Key* k = &c->keys[c->count++];
    k->t = t; k->v = v; k->tan_in = tin; k->tan_out = tout;
    if (c->monotone && curve_build_monotone(c) != 0) {
        c->count--;  /* cannot happen for keys that passed the checks above */
        return -1;
    }
    return c->count - 1;
}

/* Build Fritsch-Carlson monotone tangents. Returns -1 (without enabling mono
 * mode) when the timestamps/values are not valid. */
static int curve_build_monotone(Curve* c) {
    int n = c->count;
    if (curve_validate_times(c) != 0) return -1;
    if (n == 1) c->mono_tan[0] = 0.0f;
    if (n < 2) { c->monotone = 1; return 0; }

    float d[MAX_KEYS - 1];   /* secant slope of each interval */
    for (int k = 0; k < n - 1; k++)
        d[k] = (c->keys[k + 1].v - c->keys[k].v) /
               (c->keys[k + 1].t - c->keys[k].t);

    /* initial tangents: one-sided difference at the boundaries, average of
     * neighbouring secants inside; zero slope across a local extremum */
    c->mono_tan[0] = d[0];
    c->mono_tan[n - 1] = d[n - 2];
    for (int k = 1; k < n - 1; k++)
        c->mono_tan[k] = (d[k - 1] * d[k] <= 0.0f)
                       ? 0.0f : 0.5f * (d[k - 1] + d[k]);

    /* clip against the secants: on interval k, alpha = m_k/d_k and
     * beta = m_{k+1}/d_k must satisfy alpha^2 + beta^2 <= 9. A flat
     * interval forces both endpoint tangents to zero. */
    for (int k = 0; k < n - 1; k++) {
        if (d[k] == 0.0f) {
            c->mono_tan[k] = 0.0f;
            c->mono_tan[k + 1] = 0.0f;
            continue;
        }
        float alpha = c->mono_tan[k] / d[k];
        float beta = c->mono_tan[k + 1] / d[k];
        float s2 = alpha * alpha + beta * beta;
        if (s2 > 9.0f) {
            float tau = 3.0f / sqrtf(s2);
            c->mono_tan[k] *= tau;
            c->mono_tan[k + 1] *= tau;
        }
    }
    c->monotone = 1;
    return 0;
}

/* Switch monotone mode on/off. Turning off leaves the caller-supplied
 * tangents untouched, so eval() reproduces the original curve exactly. */
static void curve_set_monotone(Curve* c, int on) {
    if (on) curve_build_monotone(c);
    else c->monotone = 0;
}

static float curve_eval(const Curve* c, float t) {
    if (c->count == 0) return 0.0f;
    if (t <= c->keys[0].t) return c->keys[0].v;
    if (t >= c->keys[c->count - 1].t) return c->keys[c->count - 1].v;
    for (int i = 1; i < c->count; i++) {
        const Key* a = &c->keys[i - 1];
        const Key* b = &c->keys[i];
        if (t <= b->t) {
            float dt = b->t - a->t;
            float u = (t - a->t) / dt;
            float u2 = u * u, u3 = u2 * u;
            float h00 = 2 * u3 - 3 * u2 + 1;
            float h10 = u3 - 2 * u2 + u;
            float h01 = -2 * u3 + 3 * u2;
            float h11 = u3 - u2;
            float tan_a = c->monotone ? c->mono_tan[i - 1] : a->tan_out;
            float tan_b = c->monotone ? c->mono_tan[i] : b->tan_in;
            return h00 * a->v + h10 * dt * tan_a + h01 * b->v + h11 * dt * tan_b;
        }
    }
    return c->keys[c->count - 1].v;
}

/* derivative of segment i (keys i-1..i) at normalized parameter u */
static float curve_segment_deriv(const Curve* c, int i, float u) {
    const Key* a = &c->keys[i - 1];
    const Key* b = &c->keys[i];
    float dt = b->t - a->t;
    float u2 = u * u;
    float tan_a = c->monotone ? c->mono_tan[i - 1] : a->tan_out;
    float tan_b = c->monotone ? c->mono_tan[i] : b->tan_in;
    float dh_du = (6.0f * u2 - 6.0f * u) * a->v
                + (3.0f * u2 - 4.0f * u + 1.0f) * dt * tan_a
                + (-6.0f * u2 + 6.0f * u) * b->v
                + (3.0f * u2 - 2.0f * u) * dt * tan_b;
    return dh_du / dt;
}

/* Analytic derivative of the cubic hermite segment (no finite differences). */
static float curve_deriv_at(const Curve* c, float t) {
    if (c->count < 2 || t < c->keys[0].t || t > c->keys[c->count - 1].t)
        return 0.0f;  /* clamped outside the key range */
    for (int i = 1; i < c->count; i++) {
        const Key* a = &c->keys[i - 1];
        const Key* b = &c->keys[i];
        if (t <= b->t) {
            float u = (t - a->t) / (b->t - a->t);
            return curve_segment_deriv(c, i, u);
        }
    }
    return 0.0f;
}

/* derivative for continuity checks */
static float curve_deriv(const Curve* c, float t) {
    float h = 1e-3f;
    float a = curve_eval(c, t - h);
    float b = curve_eval(c, t + h);
    return (b - a) / (2 * h);
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

/* dense per-interval overshoot: every sample must stay within the endpoint
 * value range; also tracks global monotonicity via first differences */
static int interval_overshoot(const Curve* c, int samples_per_interval) {
    int bad = 0;
    for (int i = 1; i < c->count; i++) {
        float lo = fminf(c->keys[i - 1].v, c->keys[i].v);
        float hi = fmaxf(c->keys[i - 1].v, c->keys[i].v);
        float tol = 1e-4f + 1e-5f * fabsf(hi);
        for (int j = 1; j < samples_per_interval; j++) {
            float u = (float)j / (float)samples_per_interval;
            float t = c->keys[i - 1].t + u * (c->keys[i].t - c->keys[i - 1].t);
            float v = curve_eval(c, t);
            if (v < lo - tol || v > hi + tol) bad++;
        }
    }
    return bad;
}

/* interior local extrema seen in a dense uniform sampling */
static int count_extrema(const Curve* c, int steps) {
    float span = c->keys[c->count - 1].t - c->keys[0].t;
    float prev = curve_eval(c, c->keys[0].t);
    float last_d = 0.0f;
    int n = 0;
    for (int i = 1; i <= steps; i++) {
        float t = c->keys[0].t + span * (float)i / (float)steps;
        float v = curve_eval(c, t);
        float d = v - prev;
        if (d != 0.0f) {
            if (last_d != 0.0f && (d > 0.0f) != (last_d > 0.0f)) n++;
            last_d = d;
        }
        prev = v;
    }
    return n;
}

/* velocity discontinuities at interior keys: compare the exact one-sided
 * analytic derivatives (left segment at u=1, right segment at u=0) */
static int count_velocity_breaks(const Curve* c) {
    int n = 0;
    for (int k = 1; k < c->count - 1; k++) {
        float dl = curve_segment_deriv(c, k, 1.0f);
        float dr = curve_segment_deriv(c, k + 1, 0.0f);
        float tol = 1e-5f * fmaxf(1.0f, fmaxf(fabsf(dl), fabsf(dr)));
        if (fabsf(dr - dl) > tol) n++;
    }
    return n;
}

/* unclamped average-of-secants tangents, the scheme that overshoots */
static void fill_naive_tangents(Curve* c) {
    int n = c->count;
    if (n < 2) return;
    float d[MAX_KEYS - 1];
    for (int k = 0; k < n - 1; k++)
        d[k] = (c->keys[k + 1].v - c->keys[k].v) /
               (c->keys[k + 1].t - c->keys[k].t);
    c->keys[0].tan_in = c->keys[0].tan_out = d[0];
    c->keys[n - 1].tan_in = c->keys[n - 1].tan_out = d[n - 2];
    for (int k = 1; k < n - 1; k++) {
        float m = (d[k - 1] * d[k] <= 0.0f) ? 0.0f : 0.5f * (d[k - 1] + d[k]);
        c->keys[k].tan_in = c->keys[k].tan_out = m;
    }
}

int main(void) {
    Curve c;
    curve_init(&c);
    check("add_first", curve_add(&c, 0.0f, 0.0f, 0.0f, 1.0f) == 0);
    check("add_second", curve_add(&c, 1.0f, 1.0f, 1.0f, 0.0f) == 1);
    check("clamp_left", fabsf(curve_eval(&c, -1.0f)) < 1e-6f);
    check("clamp_right", fabsf(curve_eval(&c, 2.0f) - 1.0f) < 1e-6f);
    check("mid_monotonic", curve_eval(&c, 0.5f) > 0.0f && curve_eval(&c, 0.5f) < 1.0f);
    check("endpoints_exact", fabsf(curve_eval(&c, 0.0f)) < 1e-6f && fabsf(curve_eval(&c, 1.0f) - 1.0f) < 1e-6f);

    Curve flat;
    curve_init(&flat);
    curve_add(&flat, 0.0f, 0.5f, 0.0f, 0.0f);
    curve_add(&flat, 1.0f, 0.5f, 0.0f, 0.0f);
    check("flat_curve", fabsf(curve_eval(&flat, 0.5f) - 0.5f) < 1e-5f);
    check("flat_deriv_zero", fabsf(curve_deriv(&flat, 0.5f)) < 1e-3f);

    check("reject_unsorted", curve_add(&flat, 0.5f, 1.0f, 0, 0) == -1);
    check("reject_duplicate", curve_add(&flat, 1.0f, 1.0f, 0, 0) == -1);
    check("reject_nan", curve_add(&flat, NAN, 1.0f, 0, 0) == -1);

    /* validation on directly constructed data: duplicate / reversed times */
    Curve bad;
    curve_init(&bad);
    bad.count = 2;
    bad.keys[0].t = 1.0f; bad.keys[0].v = 0.0f;
    bad.keys[1].t = 1.0f; bad.keys[1].v = 1.0f;
    check("mono_reject_duplicate", curve_build_monotone(&bad) == -1 && bad.monotone == 0);
    bad.keys[1].t = 0.5f;
    check("mono_reject_reversed", curve_build_monotone(&bad) == -1 && bad.monotone == 0);
    bad.keys[1].t = 2.0f; bad.keys[0].t = NAN;
    check("mono_reject_nan", curve_build_monotone(&bad) == -1);

    Curve full;
    curve_init(&full);
    int rc = 0;
    for (int i = 0; i <= MAX_KEYS; i++) rc = curve_add(&full, (float)i, (float)i, 0, 0);
    check("capacity", rc == -1 && full.count == MAX_KEYS);

    /* ---- Fritsch-Carlson monotone mode on the classic steep-ramp dataset ---- */
    static const float kx[] = {0, 2, 3, 5, 6, 8, 9, 11, 12, 14, 15};
    static const float ky[] = {10, 10, 10, 10, 10, 10, 10.5f, 15, 50, 60, 85};
    enum { NK = (int)(sizeof(kx) / sizeof(kx[0])) };

    Curve naive;
    curve_init(&naive);
    for (int i = 0; i < NK; i++) curve_add(&naive, kx[i], ky[i], 0.0f, 0.0f);
    fill_naive_tangents(&naive);
    check("naive_average_overshoots", interval_overshoot(&naive, 200) > 0);

    Curve fc;
    curve_init(&fc);
    for (int i = 0; i < NK; i++) curve_add(&fc, kx[i], ky[i], 0.0f, 0.0f);
    check("mono_enable", curve_build_monotone(&fc) == 0 && fc.monotone == 1);
    check("mono_no_overshoot", interval_overshoot(&fc, 200) == 0);
    check("mono_no_velocity_breaks", count_velocity_breaks(&fc) == 0);

    /* adding a key while mono mode is on rebuilds the tangents */
    check("mono_rebuild_after_add",
          curve_add(&fc, 16.0f, 85.0f, 0.0f, 0.0f) == NK &&
          fc.monotone == 1 && interval_overshoot(&fc, 200) == 0);

    /* strictly decreasing data must stay non-increasing */
    Curve down;
    curve_init(&down);
    static const float dx[] = {0, 1, 2, 3, 4};
    static const float dy[] = {10, 8, 7, 2, 0};
    for (int i = 0; i < 5; i++) curve_add(&down, dx[i], dy[i], 0.0f, 0.0f);
    curve_set_monotone(&down, 1);
    check("mono_decreasing", interval_overshoot(&down, 200) == 0);

    /* flat curve in mono mode: every derived tangent is zero */
    Curve mflat;
    curve_init(&mflat);
    for (int i = 0; i < 3; i++) curve_add(&mflat, (float)i, 0.5f, 0.0f, 0.0f);
    curve_set_monotone(&mflat, 1);
    check("mono_flat", fabsf(curve_eval(&mflat, 1.5f) - 0.5f) < 1e-6f &&
                       fabsf(curve_deriv(&mflat, 0.5f)) < 1e-3f);

    /* ---- comparison against a curve with pointwise analytic derivatives ----
     * sample sin on [0, 4pi] with 9 keys; the reference takes cos(t) as the
     * tangent at every key. Mono mode must not add extrema. */
    enum { NS = 9 };
    Curve ref, mono;
    curve_init(&ref);
    curve_init(&mono);
    for (int i = 0; i < NS; i++) {
        float t = 4.0f * PI_F * (float)i / (float)(NS - 1);
        /* both curves carry the pointwise analytic tangents; mono mode
         * overrides them internally but leaves them in place for when it
         * is switched off */
        curve_add(&ref, t, sinf(t), cosf(t), cosf(t));
        curve_add(&mono, t, sinf(t), cosf(t), cosf(t));
    }
    curve_set_monotone(&mono, 1);
    int ext_ref = count_extrema(&ref, 2000);
    int ext_mono = count_extrema(&mono, 2000);
    printf("extrema_ref=%d extrema_mono=%d\n", ext_ref, ext_mono);
    check("extrema_expected", ext_ref == 4);
    check("extrema_not_increased", ext_mono <= ext_ref);
    check("mono_sin_no_overshoot", interval_overshoot(&mono, 200) == 0);
    check("mono_sin_no_velocity_breaks", count_velocity_breaks(&mono) == 0);

    /* analytic derivative must recover the pointwise tangents cos(t) of the
     * reference curve at every interior key, and agree with a central finite
     * difference at a segment midpoint */
    int deriv_ok = 1;
    for (int k = 1; k < NS - 1; k++) {
        float t = ref.keys[k].t;
        if (fabsf(curve_deriv_at(&ref, t) - cosf(t)) > 1e-5f) deriv_ok = 0;
    }
    check("analytic_deriv_at_keys", deriv_ok);
    float mid_t = 0.5f * (ref.keys[0].t + ref.keys[1].t);
    check("analytic_deriv_midpoint",
          fabsf(curve_deriv_at(&ref, mid_t) - curve_deriv(&ref, mid_t)) < 1e-2f);

    /* in mono mode the analytic derivative at every junction is the shared
     * Fritsch-Carlson tangent, i.e. C1 with no jump */
    int mono_deriv_ok = 1;
    for (int k = 1; k < NS - 1; k++) {
        float want = mono.mono_tan[k];
        float got = curve_deriv_at(&mono, mono.keys[k].t);
        if (fabsf(got - want) > 1e-6f * fmaxf(1.0f, fabsf(want))) mono_deriv_ok = 0;
    }
    check("mono_deriv_equals_tangent", mono_deriv_ok);

    /* ---- mono mode off: eval must be bit-identical to the old formula ---- */
    Curve p1;
    curve_init(&p1);
    curve_add(&p1, 0.0f, 0.0f, 0.0f, 1.0f);
    curve_add(&p1, 1.0f, 1.0f, 1.0f, 0.0f);
    check("off_parity_1", curve_eval(&p1, 0.25f) == 0x1p-2f);
    check("off_parity_2", curve_eval(&p1, 0.5f) == 0x1p-1f);
    check("off_parity_3", curve_eval(&p1, 0.75f) == 0x1.8p-1f);

    Curve p2;
    curve_init(&p2);
    curve_add(&p2, 0.0f, 0.0f, 0.0f, 0.5f);
    curve_add(&p2, 1.0f, 2.0f, 1.5f, -0.5f);
    curve_add(&p2, 2.0f, 1.0f, 0.5f, 0.0f);
    check("off_parity_4", curve_eval(&p2, 0.3f) == 0x1.a4dd3p-2f);
    check("off_parity_5", curve_eval(&p2, 0.5f) == 0x1.cp-1f);
    check("off_parity_6", curve_eval(&p2, 0.8f) == 0x1.9db22ep+0f);
    check("off_parity_7", curve_eval(&p2, 1.25f) == 0x1.cp+0f);
    check("off_parity_8", curve_eval(&p2, 1.5f) == 0x1.6p+0f);
    check("off_parity_9", curve_eval(&p2, 1.9f) == 0x1.f74bc4p-1f);

    /* toggling mono off restores the caller-tangent curve exactly:
     * while on, FC tangents differ from the analytic ones somewhere;
     * after switching off every sample must match the reference curve */
    int diffs_while_on = 0;
    for (int i = 1; i < 40; i++) {
        float t = 4.0f * PI_F * (float)i / 40.0f;
        if (curve_eval(&mono, t) != curve_eval(&ref, t)) diffs_while_on++;
    }
    curve_set_monotone(&mono, 0);
    int diffs_after_off = 0;
    for (int i = 1; i < 40; i++) {
        float t = 4.0f * PI_F * (float)i / 40.0f;
        if (curve_eval(&mono, t) != curve_eval(&ref, t)) diffs_after_off++;
    }
    check("toggle_off_restores", mono.monotone == 0 &&
          diffs_while_on > 0 && diffs_after_off == 0);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
