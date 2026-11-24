
/* curve: keyframe sampling with cubic hermite tangents and clamping.
   also: tangent smoothing at extrema with a computable curvature bound. */
#include <stdio.h>
#include <string.h>
#include <math.h>

#define MAX_KEYS 64

typedef struct { float t; float v; float tan_in; float tan_out; } Key;

typedef struct {
    Key keys[MAX_KEYS];
    int count;
} Curve;

static void curve_init(Curve* c) { memset(c, 0, sizeof(*c)); }

static int curve_add(Curve* c, float t, float v, float tin, float tout) {
    if (c->count >= MAX_KEYS) return -1;
    if (c->count > 0 && t <= c->keys[c->count - 1].t) return -1;  /* strictly increasing */
    Key* k = &c->keys[c->count++];
    k->t = t; k->v = v; k->tan_in = tin; k->tan_out = tout;
    return c->count - 1;
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
            return h00 * a->v + h10 * dt * a->tan_out + h01 * b->v + h11 * dt * b->tan_in;
        }
    }
    return c->keys[c->count - 1].v;
}

/* derivative for continuity checks */
#ifndef CURVE_NO_MAIN
static float curve_deriv(const Curve* c, float t) {
    float h = 1e-3f;
    float a = curve_eval(c, t - h);
    float b = curve_eval(c, t + h);
    return (b - a) / (2 * h);
}
#endif

/* extremum key: neighbor height differences change sign */
static int curve_is_extremum(const Curve* c, int i) {
    if (i <= 0 || i >= c->count - 1) return 0;
    float dl = c->keys[i].v - c->keys[i - 1].v;
    float dr = c->keys[i + 1].v - c->keys[i].v;
    return dl * dr < 0.0f;
}

/* exact per-segment curvature bound for the graph (t, v(t)):
   d2v/du2 is linear in u, so max|d2v| on the segment is at an endpoint;
   curvature kappa = |v''| / (1 + v'^2)^1.5 <= |v''|. */
static float curve_seg_curvature_bound(const Curve* c, int i) {
    const Key* a = &c->keys[i];
    const Key* b = &c->keys[i + 1];
    float dt = b->t - a->t;
    float dv = b->v - a->v;
    float d2_0 = 6.0f * dv - dt * (4.0f * a->tan_out + 2.0f * b->tan_in);
    float d2_1 = -6.0f * dv + dt * (2.0f * a->tan_out + 4.0f * b->tan_in);
    float m = fabsf(d2_0) > fabsf(d2_1) ? fabsf(d2_0) : fabsf(d2_1);
    return m / (dt * dt);
}

/* whole-curve curvature bound. a tangent kink at an interior key means
   unbounded curvature there, reported as INFINITY. */
static float curve_curvature_bound(const Curve* c) {
    float bound = 0.0f;
    for (int i = 0; i < c->count - 1; i++) {
        float b = curve_seg_curvature_bound(c, i);
        if (b > bound) bound = b;
    }
    for (int i = 1; i < c->count - 1; i++)
        if (fabsf(c->keys[i].tan_in - c->keys[i].tan_out) > 1e-4f) return INFINITY;
    return bound;
}

/* smooth tangents at extremum keys. each extremum gets a window of +-w
   seconds; if another extremum falls within 2w the windows overlap and
   the key is left untouched (window too wide for the local cadence --
   counted in out_skipped). otherwise tan_in/tan_out are set to a common
   tangent averaged from the secants crossing the window. if the new
   tangent would push either adjacent segment past kappa_max, it is
   clamped into the feasible interval (kappa_max enters linearly, so the
   interval is exact); infeasible extrema are flattened to 0. clamped or
   infeasible keys are the "handled separately" set, returned as count. */
static int curve_smooth_extrema(Curve* c, float w, float kappa_max, int* out_skipped) {
    int special = 0, skipped = 0;
    for (int i = 1; i < c->count - 1; i++) {
        if (!curve_is_extremum(c, i)) continue;
        float ti = c->keys[i].t;
        int overlap = 0;
        for (int j = 1; j < c->count - 1; j++) {
            if (j != i && curve_is_extremum(c, j) && fabsf(c->keys[j].t - ti) < 2.0f * w)
                overlap = 1;
        }
        if (overlap) { skipped++; continue; }

        float sum = 0.0f;
        int n = 0;
        for (int j = 0; j < c->count - 1; j++) {
            float t0 = c->keys[j].t, t1 = c->keys[j + 1].t;
            if (t1 < ti - w || t0 > ti + w) continue;
            sum += (c->keys[j + 1].v - c->keys[j].v) / (t1 - t0);
            n++;
        }
        float m = n > 0 ? sum / n : 0.0f;

        /* feasible interval for m from the two adjacent segment ends */
        float lo = -1e30f, hi = 1e30f;
        const Key* a = &c->keys[i - 1];
        const Key* b = &c->keys[i + 1];
        float dt_l = ti - a->t, dt_r = b->t - ti;
        float p = -6.0f * (c->keys[i].v - a->v) + dt_l * 2.0f * a->tan_out;
        float kl = kappa_max * dt_l * dt_l;
        float l1 = (-kl - p) / (4.0f * dt_l), h1 = (kl - p) / (4.0f * dt_l);
        float q = 6.0f * (b->v - c->keys[i].v) - dt_r * 2.0f * b->tan_in;
        float kr = kappa_max * dt_r * dt_r;
        float l2 = (q - kr) / (4.0f * dt_r), h2 = (q + kr) / (4.0f * dt_r);
        if (l1 > lo) lo = l1;
        if (l2 > lo) lo = l2;
        if (h1 < hi) hi = h1;
        if (h2 < hi) hi = h2;

        if (lo > hi) { m = 0.0f; special++; }        /* infeasible: flatten, report */
        else if (m < lo) { m = lo; special++; }
        else if (m > hi) { m = hi; special++; }
        c->keys[i].tan_in = m;
        c->keys[i].tan_out = m;
    }
    if (out_skipped) *out_skipped = skipped;
    return special;
}

#ifndef CURVE_NO_MAIN
static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
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

    Curve full;
    curve_init(&full);
    int rc = 0;
    for (int i = 0; i <= MAX_KEYS; i++) rc = curve_add(&full, (float)i, (float)i, 0, 0);
    check("capacity", rc == -1 && full.count == MAX_KEYS);

    /* extremum smoothing: kinked peak -> C1, finite curvature bound */
    Curve sw;
    curve_init(&sw);
    curve_add(&sw, 0.0f, 0.0f, 0.0f, 1.0f);
    curve_add(&sw, 1.0f, 1.0f, 0.0f, -1.0f);   /* kinked peak */
    curve_add(&sw, 2.0f, 0.0f, 0.0f, 0.0f);
    check("extremum_detected", curve_is_extremum(&sw, 1) == 1);
    check("extremum_not_endpoint", curve_is_extremum(&sw, 0) == 0);
    check("kink_bound_infinite", isinf(curve_curvature_bound(&sw)));
    int sp = curve_smooth_extrema(&sw, 0.25f, 8.0f, NULL);
    check("smooth_c1", fabsf(sw.keys[1].tan_in - sw.keys[1].tan_out) < 1e-6f);
    check("smooth_no_special", sp == 0);
    check("smooth_within_bound", curve_curvature_bound(&sw) <= 8.0f + 1e-4f);

    /* over-bound extremum handled separately: spike cannot fit kappa_max */
    Curve spike;
    curve_init(&spike);
    curve_add(&spike, 0.0f, 0.0f, 0.0f, 0.0f);
    curve_add(&spike, 1.0f, 2.0f, 0.0f, 0.0f);
    curve_add(&spike, 2.0f, 0.0f, 0.0f, 0.0f);
    sp = curve_smooth_extrema(&spike, 0.25f, 5.0f, NULL);
    check("spike_special", sp == 1);
    check("spike_still_c1", fabsf(spike.keys[1].tan_in - spike.keys[1].tan_out) < 1e-6f);

    /* feasible clamp: averaged tangent pulled back into the bound */
    Curve cl;
    curve_init(&cl);
    curve_add(&cl, 0.0f, 0.0f, 0.0f, -4.0f);
    curve_add(&cl, 1.0f, -2.0f, 0.0f, 1.0f);   /* kinked trough */
    curve_add(&cl, 2.0f, -1.0f, 0.2f, 0.0f);
    sp = curve_smooth_extrema(&cl, 0.25f, 5.0f, NULL);
    check("clamp_special", sp == 1);
    check("clamp_within_bound", curve_curvature_bound(&cl) <= 5.0f + 1e-4f);

    /* window too wide for the spacing: extrema skipped, kinks remain */
    Curve dense;
    curve_init(&dense);
    curve_add(&dense, 0.0f, 0.0f, 0.0f, 5.0f);
    curve_add(&dense, 0.2f, 1.0f, 0.0f, -5.0f);
    curve_add(&dense, 0.4f, 0.0f, 0.0f, 5.0f);
    curve_add(&dense, 0.6f, 1.0f, 0.0f, -5.0f);
    curve_add(&dense, 0.8f, 0.0f, 0.0f, 0.0f);
    int skipped = 0;
    sp = curve_smooth_extrema(&dense, 0.25f, 1e6f, &skipped);
    check("overlap_skipped", skipped == 3 && sp == 0);
    check("overlap_kink_remains", isinf(curve_curvature_bound(&dense)));

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
#endif
