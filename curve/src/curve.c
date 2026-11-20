
/* curve: keyframe sampling with cubic hermite tangents and clamping. */
#include <stdio.h>
#include <string.h>
#include <math.h>

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

/* ---- error-band fitting ----
 * curve_fit: 把稠密采样 (ts, vs) 按 L∞ 误差带分段，段内用三次埃尔米特
 * 拟合（端点精确插值，切线由采样差分估计，段界共享切线保证 C1）。
 * 贪心分段：段在参数空间尽量向右延伸，内部采样点偏差超过 band 即切断。
 * 拒绝（返回 -1，out 不变）：band <= 0、n < 2、ts 非严格递增、
 * max_segments 越界，或拟合所需段数超过 max_segments。
 * 成功时逐点复核：拟合曲线在每个原始采样点处的偏差都不超过 band，
 * 最大偏差写入 info->max_dev。 */

typedef struct {
    int segments;   /* 拟合曲线的三次段数（关键帧数 - 1） */
    float max_dev;  /* 全部采样点上 |拟合值 - 原始值| 的最大值 */
} FitInfo;

/* tangent at sample i by central difference (one-sided at the ends) */
static float fit_tangent(const float* ts, const float* vs, int n, int i) {
    if (i <= 0)     return (vs[1] - vs[0]) / (ts[1] - ts[0]);
    if (i >= n - 1) return (vs[n - 1] - vs[n - 2]) / (ts[n - 1] - ts[n - 2]);
    return (vs[i + 1] - vs[i - 1]) / (ts[i + 1] - ts[i - 1]);
}

/* same cubic hermite form as curve_eval, on an explicit span */
static float fit_span_eval(float t0, float v0, float m0,
                           float t1, float v1, float m1, float t) {
    float dt = t1 - t0;
    float u = (t - t0) / dt;
    float u2 = u * u, u3 = u2 * u;
    float h00 = 2 * u3 - 3 * u2 + 1;
    float h10 = u3 - 2 * u2 + u;
    float h01 = -2 * u3 + 3 * u2;
    float h11 = u3 - u2;
    return h00 * v0 + h10 * dt * m0 + h01 * v1 + h11 * dt * m1;
}

/* max deviation of the span s..e at its interior samples */
static float fit_span_err(const float* ts, const float* vs, int n, int s, int e) {
    float m0 = fit_tangent(ts, vs, n, s);
    float m1 = fit_tangent(ts, vs, n, e);
    float err = 0.0f;
    for (int i = s + 1; i < e; i++) {
        float d = fabsf(fit_span_eval(ts[s], vs[s], m0, ts[e], vs[e], m1, ts[i]) - vs[i]);
        if (d > err) err = d;
    }
    return err;
}

static int curve_fit(const float* ts, const float* vs, int n,
                     float band, int max_segments, Curve* out, FitInfo* info) {
    if (!ts || !vs || !out || n < 2) return -1;
    if (!(band > 0.0f)) return -1;                      /* zero/negative band */
    if (max_segments < 1 || max_segments > MAX_KEYS - 1) return -1;
    for (int i = 1; i < n; i++)
        if (!(ts[i] > ts[i - 1])) return -1;            /* strictly increasing */

    Curve tmp;
    curve_init(&tmp);
    int segments = 0;
    int s = 0;
    curve_add(&tmp, ts[0], vs[0], fit_tangent(ts, vs, n, 0), fit_tangent(ts, vs, n, 0));
    while (s < n - 1) {
        int e = s + 1;                                   /* one interval always fits */
        while (e < n - 1 && fit_span_err(ts, vs, n, s, e + 1) <= band) e++;
        float m = fit_tangent(ts, vs, n, e);
        curve_add(&tmp, ts[e], vs[e], m, m);             /* shared tangent: C1 */
        segments++;
        if (segments > max_segments) return -1;          /* too many segments */
        s = e;
    }

    /* pointwise verification against the original samples */
    float max_dev = 0.0f;
    for (int i = 0; i < n; i++) {
        float d = fabsf(curve_eval(&tmp, ts[i]) - vs[i]);
        if (d > max_dev) max_dev = d;
    }
    if (max_dev > band) return -1;

    *out = tmp;
    if (info) { info->segments = segments; info->max_dev = max_dev; }
    return 0;
}

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

    /* ---- error-band fitting ---- */
#define FIT_N 241
    float ts[FIT_N], vs[FIT_N];
    for (int i = 0; i < FIT_N; i++) {
        ts[i] = (float)i * 0.01f;                       /* 0 .. 2.4 */
        vs[i] = sinf(3.0f * ts[i]) + 0.3f * sinf(17.0f * ts[i]);
    }

    Curve fc;
    curve_init(&fc);
    FitInfo fi = { 0, 0.0f };
    float band = 0.05f;
    check("fit_ok", curve_fit(ts, vs, FIT_N, band, MAX_KEYS - 1, &fc, &fi) == 0);
    check("fit_reports", fi.segments >= 1 && fi.segments <= MAX_KEYS - 1 &&
          fi.max_dev >= 0.0f && fi.max_dev <= band);
    printf("fit band=%.3f segments=%d max_dev=%.6f\n", band, fi.segments, fi.max_dev);

    /* pointwise: fitted curve stays within the band at every sample */
    int within = 1;
    for (int i = 0; i < FIT_N; i++)
        if (fabsf(curve_eval(&fc, ts[i]) - vs[i]) > band) within = 0;
    check("fit_pointwise", within);

    /* rejections */
    Curve guard;
    curve_init(&guard);
    curve_add(&guard, 0.0f, 7.0f, 0.0f, 0.0f);
    check("fit_reject_zero_band", curve_fit(ts, vs, FIT_N, 0.0f, MAX_KEYS - 1, &guard, 0) == -1);
    check("fit_reject_neg_band", curve_fit(ts, vs, FIT_N, -0.1f, MAX_KEYS - 1, &guard, 0) == -1);
    check("fit_reject_seg_limit", curve_fit(ts, vs, FIT_N, 0.001f, 2, &guard, 0) == -1);
    check("fit_reject_keeps_out", guard.count == 1 && guard.keys[0].v == 7.0f);

    /* segments vs band: non-increasing as the band widens;
     * bands too tight for the segment cap are rejected, not truncated */
    float sweep[] = { 0.005f, 0.01f, 0.02f, 0.05f, 0.1f, 0.2f, 0.5f };
    int prev = MAX_KEYS;
    int monotone = 1;
    int rejected = 0;
    for (int i = 0; i < 7; i++) {
        Curve sc;
        FitInfo si = { 0, 0.0f };
        if (curve_fit(ts, vs, FIT_N, sweep[i], MAX_KEYS - 1, &sc, &si) != 0) {
            printf("sweep band=%.3f rejected(segments>%d)\n", sweep[i], MAX_KEYS - 1);
            rejected++;
            continue;
        }
        printf("sweep band=%.3f segments=%d max_dev=%.6f\n", sweep[i], si.segments, si.max_dev);
        if (si.segments > prev || si.max_dev > sweep[i]) monotone = 0;
        prev = si.segments;
    }
    check("fit_sweep_monotone", monotone);
    check("fit_sweep_some_ok", rejected < 7);

    /* fitting is opt-in: hand-built curves evaluate exactly as before */
    check("fit_disabled_unchanged", fabsf(curve_eval(&c, 0.5f) - 0.5f) < 1e-6f);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
