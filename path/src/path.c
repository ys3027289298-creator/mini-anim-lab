
/* path: polyline arc-length table, constant-speed sampling, curvature and
   closest point queries. */
#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "../../shared/pose_export.h"

#define MAX_PTS 32
#define POSE_SAMPLES 15

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

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

/* ---------- unified pose export (new format, see shared/pose_export.h) ----- */

/* canonical pose: L-shaped polyline (0,0)->(3,0)->(3,4), total length 7 */
static void path_pose_build(Path* pa) {
    path_init(pa);
    path_add(pa, 0, 0);
    path_add(pa, 3, 0);
    path_add(pa, 3, 4);
}

static int path_pose_export(const char* outdir, int argc, char** argv) {
    if (pe_mkdir_p(outdir) != 0) return 2;

    P2 pts[POSE_SAMPLES];
    float ss[POSE_SAMPLES];
    /* timing tier 1: pose numeric computation, 1000 repeats */
    Path pa;
    long t0 = pe_now_us();
    for (int rep = 0; rep < 1000; rep++) {
        path_pose_build(&pa);
        for (int i = 0; i < POSE_SAMPLES; i++) {
            ss[i] = pa.total * (float)i / (POSE_SAMPLES - 1);
            if (path_sample(&pa, ss[i], &pts[i], NULL) < 0) return 2;
        }
        if (pts[0].x + pts[POSE_SAMPLES - 1].y > 1e30f) return 2;
    }
    long t1 = pe_now_us();
    long compute_us = t1 - t0;

    /* timing tier 2: image rendering, single pass */
    PeCanvas cv;
    pe_canvas_init(&cv, -0.5f, -0.5f, 4.0f, 4.5f);
    long t2 = pe_now_us();
    for (int i = 1; i < pa.n; i++)
        pe_line(&cv, pa.p[i - 1].x, pa.p[i - 1].y, pa.p[i].x, pa.p[i].y,
                40, 90, 180);
    for (int i = 0; i < POSE_SAMPLES; i++)
        pe_circle(&cv, pts[i].x, pts[i].y, 0.06f, 200, 60, 60);
    long t3 = pe_now_us();
    long render_us = t3 - t2;

    char ppm_path[600], json_path[600];
    snprintf(ppm_path, sizeof(ppm_path), "%s/pose.ppm", outdir);
    snprintf(json_path, sizeof(json_path), "%s/pose.json", outdir);
    if (pe_canvas_save(&cv, ppm_path) != 0) return 2;

    char hash[65];
    if (pe_sha256_file(ppm_path, hash) != 0) return 2;

    FILE* f = fopen(json_path, "wb");
    if (!f) return 2;
    fprintf(f, "{\n");
    fprintf(f, "  \"module\": \"path\",\n");
    fprintf(f, "  \"format_version\": %d,\n", POSE_FORMAT_VERSION);
    fprintf(f, "  \"pose\": {\n");
    fprintf(f, "    \"point_count\": %d,\n", pa.n);
    fprintf(f, "    \"total_length\": %.6f,\n", pa.total);
    fprintf(f, "    \"sample_count\": %d,\n", POSE_SAMPLES);
    fprintf(f, "    \"samples\": [\n");
    for (int i = 0; i < POSE_SAMPLES; i++) {
        fprintf(f, "      {\"s\": %.6f, \"x\": %.6f, \"y\": %.6f}%s\n",
                ss[i], pts[i].x, pts[i].y,
                i == POSE_SAMPLES - 1 ? "" : ",");
    }
    fprintf(f, "    ]\n");
    fprintf(f, "  },\n");
    fprintf(f, "  \"image\": {\n");
    fprintf(f, "    \"file\": \"pose.ppm\",\n");
    fprintf(f, "    \"format\": \"ppm-p6\",\n");
    fprintf(f, "    \"width\": %d,\n", POSE_CANVAS_W);
    fprintf(f, "    \"height\": %d,\n", POSE_CANVAS_H);
    fprintf(f, "    \"sha256\": \"%s\"\n", hash);
    fprintf(f, "  },\n");
    fprintf(f, "  \"timing\": {\n");
    fprintf(f, "    \"compute_pose_us\": %ld,\n", compute_us);
    fprintf(f, "    \"compute_repeats\": 1000,\n");
    fprintf(f, "    \"render_image_us\": %ld\n", render_us);
    fprintf(f, "  },\n");
    fprintf(f, "  \"reproduce\": ");
    fputc('[', f);
    for (int i = 0; i < argc; i++) {
        pe_write_json_string(f, argv[i]);
        fputs(i == argc - 1 ? "" : ", ", f);
    }
    fputs("],\n", f);
    fprintf(f, "  \"reproduce_hint\": \"run the command above from the repository root\"\n");
    fprintf(f, "}\n");
    fclose(f);
    return 0;
}

int main(int argc, char** argv) {
    if (strcmp(argc > 1 ? argv[1] : "", "--pose-export") == 0) {
        if (argc != 3) {
            fprintf(stderr, "usage: %s --pose-export <outdir>\n", argv[0]);
            return 2;
        }
        return path_pose_export(argv[2], argc, argv);
    }

    Path pa;
    path_init(&pa);
    check("add0", path_add(&pa, 0, 0) == 0);
    check("add1", path_add(&pa, 3, 0) == 1);
    check("add2", path_add(&pa, 3, 4) == 2);
    check("total_len", fabsf(pa.total - 7.0f) < 1e-5f);

    P2 s;
    check("sample_mid", path_sample(&pa, 1.5f, &s, NULL) == 0 && fabsf(s.x - 1.5f) < 1e-5f);
    check("sample_corner", path_sample(&pa, 3.0f, &s, NULL) >= 0 && fabsf(s.x - 3.0f) < 1e-5f && fabsf(s.y) < 1e-5f);
    check("sample_end", path_sample(&pa, 7.0f, &s, NULL) >= 0 && fabsf(s.y - 4.0f) < 1e-4f);
    check("sample_out", path_sample(&pa, 8.0f, &s, NULL) == -1);

    check("dup_rejected", path_add(&pa, 3, 4) == -1);

    float best_s = 0;
    float d = path_closest(&pa, (P2){ 1.5f, 2.0f }, &best_s);
    check("closest_dist", fabsf(d - 1.5f) < 1e-4f);
    check("closest_s", fabsf(best_s - 5.0f) < 1e-4f);

    Path straight;
    path_init(&straight);
    path_add(&straight, 0, 0); path_add(&straight, 1, 0); path_add(&straight, 2, 0);
    check("curvature_straight", fabsf(path_curvature(&straight, 1)) < 1e-5f);

    Path bend;
    path_init(&bend);
    path_add(&bend, 0, 0); path_add(&bend, 1, 0); path_add(&bend, 1, 1);
    check("curvature_turn", path_curvature(&bend, 1) > 0.0f);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
