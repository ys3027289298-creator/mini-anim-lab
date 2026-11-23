
/* curve: keyframe sampling with cubic hermite tangents and clamping. */
#include <stdio.h>
#include <string.h>
#include <math.h>

#ifndef MAL_MODULE
#define MAL_MODULE "curve"
#endif
#define MAL_TIME_ITERS 1200000L
#include "../../common/harness.h"

#define MAX_KEYS 16

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
static float curve_deriv(const Curve* c, float t) {
    float h = 1e-3f;
    float a = curve_eval(c, t - h);
    float b = curve_eval(c, t + h);
    return (b - a) / (2 * h);
}

/* the rise curve from the sampling checks, reused by the evidence paths */
MAL_UNUSED static void curve_build(Curve* c) {
    curve_init(c);
    curve_add(c, 0.0f, 0.0f, 0.0f, 1.0f);
    curve_add(c, 1.0f, 1.0f, 1.0f, 0.0f);
}

static void run_checks(void) {
    Curve c;
    curve_init(&c);
    mal_check("add_first", curve_add(&c, 0.0f, 0.0f, 0.0f, 1.0f) == 0);
    mal_check("add_second", curve_add(&c, 1.0f, 1.0f, 1.0f, 0.0f) == 1);
    mal_check("clamp_left", fabsf(curve_eval(&c, -1.0f)) < 1e-6f);
    mal_check("clamp_right", fabsf(curve_eval(&c, 2.0f) - 1.0f) < 1e-6f);
    mal_check("mid_monotonic", curve_eval(&c, 0.5f) > 0.0f && curve_eval(&c, 0.5f) < 1.0f);
    mal_check("endpoints_exact", fabsf(curve_eval(&c, 0.0f)) < 1e-6f && fabsf(curve_eval(&c, 1.0f) - 1.0f) < 1e-6f);

    Curve flat;
    curve_init(&flat);
    curve_add(&flat, 0.0f, 0.5f, 0.0f, 0.0f);
    curve_add(&flat, 1.0f, 0.5f, 0.0f, 0.0f);
    mal_check("flat_curve", fabsf(curve_eval(&flat, 0.5f) - 0.5f) < 1e-5f);
    mal_check("flat_deriv_zero", fabsf(curve_deriv(&flat, 0.5f)) < 1e-3f);

    mal_check("reject_unsorted", curve_add(&flat, 0.5f, 1.0f, 0, 0) == -1);

    Curve full;
    curve_init(&full);
    int rc = 0;
    for (int i = 0; i <= MAX_KEYS; i++) rc = curve_add(&full, (float)i, (float)i, 0, 0);
    mal_check("capacity", rc == -1 && full.count == MAX_KEYS);
}

int main(void) {
    run_checks();

#ifdef MAL_EVIDENCE_POSE
    if (mal_failures == 0) {
        Curve c;
        curve_build(&c);
        mal_pose("keys", c.count);
        for (int i = 0; i <= 8; i++) {
            char k[16];
            snprintf(k, sizeof k, "t%d", i);
            mal_pose(k, curve_eval(&c, (float)i / 8.0f));
        }
        mal_pose("deriv_050", curve_deriv(&c, 0.5f));
        mal_pose("clamp_lo", curve_eval(&c, -1.0f));
        mal_pose("clamp_hi", curve_eval(&c, 2.0f));
        mal_emit_pose();
    }
#endif

#ifdef MAL_EVIDENCE_IMAGE
    if (mal_failures == 0) {
        MalImage im;
        Curve c;
        curve_build(&c);
        mal_img_clear(&im);
        double px = 0.0, py = 0.0;
        for (int i = 0; i <= 256; i++) {
            double t = (double)i / 256.0;
            double x = 16.0 + t * 96.0;
            double y = 112.0 - curve_eval(&c, (float)t) * 96.0;
            if (i > 0) mal_img_line(&im, px, py, x, y, 255, 180, 90);
            px = x; py = y;
        }
        mal_emit_png(&im, MAL_MODULE);
    }
#endif

#ifdef MAL_EVIDENCE_TIMING
    if (mal_failures == 0) {
        Curve c;
        curve_build(&c);
        double acc = 0.0;
        double t0 = mal_now();
        for (long it = 0; it < MAL_TIME_ITERS; it++) {
            float t = (float)(it % 1000) / 1000.0f;
            acc += curve_eval(&c, t) + curve_deriv(&c, t);
        }
        double t1 = mal_now();
        if (!isfinite(acc)) mal_check("time_finite", 0);
        mal_emit_time(MAL_TIME_TIER, t1 - t0, MAL_TIME_ITERS);
    }
#endif

    printf("all_ok=%d\n", mal_failures == 0);
    return mal_failures == 0 ? 0 : 1;
}
