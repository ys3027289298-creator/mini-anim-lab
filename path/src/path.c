
/* path: polyline arc-length table, constant-speed sampling, curvature and
   closest point queries. */
#include <stdio.h>
#include <math.h>
#include <string.h>

#ifndef MAL_MODULE
#define MAL_MODULE "path"
#endif
#define MAL_TIME_ITERS 1200000L
#include "../../common/harness.h"

#define MAX_PTS 32

typedef struct { float x, y; } P2;

typedef struct {
    P2 p[MAX_PTS];
    float cum[MAX_PTS];    /* cumulative arc length at each point */
    int n;
    float total;
} Path;

static void path_init(Path* pa) { memset(pa, 0, sizeof(*pa)); }

static int path_add(Path* pa, float x, float y) {
    if (pa->n >= MAX_PTS) return -1;
    if (pa->n > 0) {
        P2 prev = pa->p[pa->n - 1];
        float dx = x - prev.x, dy = y - prev.y;
        float d = sqrtf(dx * dx + dy * dy);
        if (d < 1e-6f) return -1;               /* duplicate point rejected */
        pa->cum[pa->n] = pa->cum[pa->n - 1] + d;
        pa->total = pa->cum[pa->n];
    }
    pa->p[pa->n].x = x;
    pa->p[pa->n].y = y;
    pa->n++;
    return pa->n - 1;
}

/* constant speed sample: s in [0, total] */
static int path_sample(const Path* pa, float s, P2* out, float* seg_t) {
    if (pa->n < 2) return -1;
    if (s < 0.0f || s > pa->total) return -1;
    for (int i = 1; i < pa->n; i++) {
        if (s <= pa->cum[i]) {
            float seg = pa->cum[i] - pa->cum[i - 1];
            float t = seg > 0.0f ? (s - pa->cum[i - 1]) / seg : 0.0f;
            out->x = pa->p[i - 1].x + (pa->p[i].x - pa->p[i - 1].x) * t;
            out->y = pa->p[i - 1].y + (pa->p[i].y - pa->p[i - 1].y) * t;
            if (seg_t) *seg_t = t;
            return i - 1;
        }
    }
    return -1;
}

/* closest point on path to q, returns distance and arc length */
static float path_closest(const Path* pa, P2 q, float* out_s) {
    float best = 1e30f, best_s = 0.0f;
    for (int i = 1; i < pa->n; i++) {
        P2 a = pa->p[i - 1], b = pa->p[i];
        float vx = b.x - a.x, vy = b.y - a.y;
        float wx = q.x - a.x, wy = q.y - a.y;
        float vv = vx * vx + vy * vy;
        float t = vv > 0.0f ? (wx * vx + wy * vy) / vv : 0.0f;
        if (t < 0.0f) t = 0.0f;
        if (t > 1.0f) t = 1.0f;
        float px = a.x + vx * t, py = a.y + vy * t;
        float dx = q.x - px, dy = q.y - py;
        float d = sqrtf(dx * dx + dy * dy);
        if (d < best) {
            best = d;
            best_s = pa->cum[i - 1] + (pa->cum[i] - pa->cum[i - 1]) * t;
        }
    }
    if (out_s) *out_s = best_s;
    return best;
}

/* signed curvature at segment i (turn angle over length) */
static float path_curvature(const Path* pa, int i) {
    if (i < 1 || i >= pa->n - 1) return 0.0f;
    P2 a = pa->p[i - 1], b = pa->p[i], c = pa->p[i + 1];
    float abx = b.x - a.x, aby = b.y - a.y;
    float bcx = c.x - b.x, bcy = c.y - b.y;
    float cross = abx * bcy - aby * bcx;
    float la = sqrtf(abx * abx + aby * aby);
    float lb = sqrtf(bcx * bcx + bcy * bcy);
    float lc = sqrtf((c.x - a.x) * (c.x - a.x) + (c.y - a.y) * (c.y - a.y));
    if (la * lb * lc == 0.0f) return 0.0f;
    return 2.0f * cross / (la * lb * lc);
}

/* the L-shaped path used by the sampling checks, reused by the evidence */
MAL_UNUSED static void path_build(Path* pa) {
    path_init(pa);
    path_add(pa, 0.0f, 0.0f);
    path_add(pa, 3.0f, 0.0f);
    path_add(pa, 3.0f, 4.0f);
}

static void run_checks(void) {
    Path pa;
    path_init(&pa);
    mal_check("add0", path_add(&pa, 0, 0) == 0);
    mal_check("add1", path_add(&pa, 3, 0) == 1);
    mal_check("add2", path_add(&pa, 3, 4) == 2);
    mal_check("total_len", fabsf(pa.total - 7.0f) < 1e-5f);

    P2 s;
    mal_check("sample_mid", path_sample(&pa, 1.5f, &s, NULL) == 0 && fabsf(s.x - 1.5f) < 1e-5f);
    mal_check("sample_corner", path_sample(&pa, 3.0f, &s, NULL) >= 0 && fabsf(s.x - 3.0f) < 1e-5f && fabsf(s.y) < 1e-5f);
    mal_check("sample_end", path_sample(&pa, 7.0f, &s, NULL) >= 0 && fabsf(s.y - 4.0f) < 1e-4f);
    mal_check("sample_out", path_sample(&pa, 8.0f, &s, NULL) == -1);

    mal_check("dup_rejected", path_add(&pa, 3, 4) == -1);

    float best_s = 0;
    float d = path_closest(&pa, (P2){ 1.5f, 2.0f }, &best_s);
    mal_check("closest_dist", fabsf(d - 1.5f) < 1e-4f);
    mal_check("closest_s", fabsf(best_s - 5.0f) < 1e-4f);

    Path straight;
    path_init(&straight);
    path_add(&straight, 0, 0); path_add(&straight, 1, 0); path_add(&straight, 2, 0);
    mal_check("curvature_straight", fabsf(path_curvature(&straight, 1)) < 1e-5f);

    Path bend;
    path_init(&bend);
    path_add(&bend, 0, 0); path_add(&bend, 1, 0); path_add(&bend, 1, 1);
    mal_check("curvature_turn", path_curvature(&bend, 1) > 0.0f);
}

int main(void) {
    run_checks();

#ifdef MAL_EVIDENCE_POSE
    if (mal_failures == 0) {
        Path pa;
        path_build(&pa);
        mal_pose("n", pa.n);
        mal_pose("total", pa.total);
        for (int i = 1; i <= 7; i++) {
            P2 s;
            char k[16];
            snprintf(k, sizeof k, "s%d_x", i);
            if (path_sample(&pa, (float)i, &s, NULL) >= 0) mal_pose(k, s.x); else mal_pose(k, NAN);
            snprintf(k, sizeof k, "s%d_y", i);
            if (path_sample(&pa, (float)i, &s, NULL) >= 0) mal_pose(k, s.y); else mal_pose(k, NAN);
        }
        float best_s = 0.0f;
        mal_pose("closest", path_closest(&pa, (P2){ 1.5f, 2.0f }, &best_s));
        mal_pose("closest_s", best_s);
        mal_emit_pose();
    }
#endif

#ifdef MAL_EVIDENCE_IMAGE
    if (mal_failures == 0) {
        MalImage im;
        Path pa;
        path_build(&pa);
        mal_img_clear(&im);
        const double scale = 18.0, ox = 20.0, oy = 120.0;
        for (int i = 0; i < pa.n - 1; i++)
            mal_img_line(&im, ox + pa.p[i].x * scale, oy - pa.p[i].y * scale,
                         ox + pa.p[i + 1].x * scale, oy - pa.p[i + 1].y * scale,
                         80, 200, 255);
        for (int i = 1; i <= 7; i++) {
            P2 s;
            if (path_sample(&pa, (float)i, &s, NULL) >= 0)
                mal_img_line(&im, ox + s.x * scale - 2, oy - s.y * scale,
                             ox + s.x * scale + 2, oy - s.y * scale, 255, 180, 90);
        }
        mal_emit_png(&im, MAL_MODULE);
    }
#endif

#ifdef MAL_EVIDENCE_TIMING
    if (mal_failures == 0) {
        Path pa;
        path_build(&pa);
        double acc = 0.0;
        double t0 = mal_now();
        for (long it = 0; it < MAL_TIME_ITERS; it++) {
            float s = (float)(it % 70) / 10.0f;
            P2 out;
            if (path_sample(&pa, s, &out, NULL) >= 0) acc += out.x + out.y;
            acc += path_closest(&pa, (P2){ 1.5f, 2.0f }, NULL);
            acc += path_curvature(&pa, 1) + path_curvature(&pa, 2);
        }
        double t1 = mal_now();
        if (!isfinite(acc)) mal_check("time_finite", 0);
        mal_emit_time(MAL_TIME_TIER, t1 - t0, MAL_TIME_ITERS);
    }
#endif

    printf("all_ok=%d\n", mal_failures == 0);
    return mal_failures == 0 ? 0 : 1;
}
