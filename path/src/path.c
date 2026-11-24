
/* path: polyline arc-length table, constant-speed sampling, curvature and
   closest point queries. */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "../../shared/pose_format.h"

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

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

/* ===================== shared pose format (new) ==========================
 *
 * A "meta" record carries the total length, then one "pt" record per
 * polyline vertex with its cumulative arc length:
 *
 *   meta total=7
 *   pt x=0 y=0 cum=0
 *   pt x=3 y=0 cum=3
 *   pt x=3 y=4 cum=7
 *
 * Reading it back restores points, the arc-length table and total; the
 * roundtrip compares all three internal arrays.
 */

static Path path_pose_build(void) {
    Path pa;
    path_init(&pa);
    path_add(&pa, 0.0f, 0.0f);
    path_add(&pa, 3.0f, 0.0f);
    path_add(&pa, 3.0f, 4.0f);
    return pa;
}

static int path_export(const Path* pa, const char* file) {
    FILE* f = fopen(file, "w");
    if (!f) return -1;
    pf_write_header(f, "path", pa->n + 1);
    pf_record_begin(f, "meta");
    pf_field_float(f, "total", pa->total);
    pf_record_end(f);
    for (int i = 0; i < pa->n; i++) {
        pf_record_begin(f, "pt");
        pf_field_float(f, "x", pa->p[i].x);
        pf_field_float(f, "y", pa->p[i].y);
        pf_field_float(f, "cum", pa->cum[i]);
        pf_record_end(f);
    }
    fclose(f);
    return 0;
}

static int path_import(const char* file, Path* pa) {
    FILE* f = fopen(file, "r");
    if (!f) return -1;

    int n = 0;
    if (pf_read_header(f, "path", &n) != 0 || n < 1 || n > MAX_PTS + 1) {
        fclose(f);
        return -1;
    }

    path_init(pa);
    PfRecord r;
    if (pf_read_record(f, &r) != 0 || strcmp(r.name, "meta") != 0 ||
        pf_get_float(&r, "total", &pa->total) != 0) {
        fclose(f);
        return -1;
    }
    for (int k = 1; k < n; k++) {
        if (pf_read_record(f, &r) != 0 || strcmp(r.name, "pt") != 0 ||
            pf_get_float(&r, "x", &pa->p[k - 1].x) != 0 ||
            pf_get_float(&r, "y", &pa->p[k - 1].y) != 0 ||
            pf_get_float(&r, "cum", &pa->cum[k - 1]) != 0) {
            fclose(f);
            return -1;
        }
    }
    pa->n = n - 1;
    fclose(f);
    return 0;
}

/* ---- legacy path format (pre-unification): two headerless sections -----
 *
 *   total <total>
 *   <x>,<y>,<cum>
 *   ...
 */

static int path_export_legacy(const Path* pa, const char* file) {
    FILE* f = fopen(file, "w");
    if (!f) return -1;
    fprintf(f, "total %.9g\n", (double)pa->total);
    for (int i = 0; i < pa->n; i++)
        fprintf(f, "%.9g,%.9g,%.9g\n",
                (double)pa->p[i].x, (double)pa->p[i].y, (double)pa->cum[i]);
    fclose(f);
    return 0;
}

static int path_import_legacy(const char* file, Path* pa) {
    FILE* f = fopen(file, "r");
    if (!f) return -1;

    path_init(pa);
    char line[256];
    int have_total = 0, ok = 1;
    while (ok && fgets(line, sizeof line, f)) {
        char* p = line;
        while (*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n') p++;
        if (*p == '\0' || *p == '#') continue;

        if (strncmp(p, "total", 5) == 0) {
            if (sscanf(p, "total %f", &pa->total) != 1) { ok = 0; break; }
            have_total = 1;
            continue;
        }
        float x, y, cum;
        if (sscanf(p, "%f,%f,%f", &x, &y, &cum) != 3 || pa->n >= MAX_PTS) {
            ok = 0;
            break;
        }
        pa->p[pa->n].x = x;
        pa->p[pa->n].y = y;
        pa->cum[pa->n] = cum;
        pa->n++;
    }
    fclose(f);
    return (ok && have_total && pa->n > 0) ? 0 : -1;
}

static int path_same(const Path* a, const Path* b) {
    if (a->n != b->n || a->total != b->total) return 0;
    for (int i = 0; i < a->n; i++) {
        if (a->p[i].x != b->p[i].x || a->p[i].y != b->p[i].y ||
            a->cum[i] != b->cum[i]) return 0;
    }
    return 1;
}

/* ------------------------------------------------------------ CLI modes */

static int path_cli(int argc, char** argv) {
    const char* mode = argv[1];
    Path src = path_pose_build(), back;

    if (strcmp(mode, "--pose-export") == 0 && argc == 3)
        return path_export(&src, argv[2]) == 0 ? 0 : 1;
    if (strcmp(mode, "--pose-import") == 0 && argc == 3) {
        if (path_import(argv[2], &back) != 0) { printf("imported=0\n"); return 1; }
        printf("imported=1\nn=%d\n", back.n);
        return 0;
    }
    if (strcmp(mode, "--pose-roundtrip") == 0 && argc == 3) {
        int ok = path_export(&src, argv[2]) == 0 &&
                 path_import(argv[2], &back) == 0 &&
                 path_same(&src, &back);
        printf("roundtrip=%d\n", ok ? 1 : 0);
        return ok ? 0 : 1;
    }
    if (strcmp(mode, "--legacy-export") == 0 && argc == 3)
        return path_export_legacy(&src, argv[2]) == 0 ? 0 : 1;
    if (strcmp(mode, "--legacy-check") == 0 && argc == 3) {
        int ok = path_import_legacy(argv[2], &back) == 0 && path_same(&src, &back);
        printf("legacy=%d\n", ok ? 1 : 0);
        return ok ? 0 : 1;
    }

    fprintf(stderr, "usage: %s [--pose-export|--pose-import|--pose-roundtrip|"
                    "--legacy-export|--legacy-check] <file>\n", argv[0]);
    return 2;
}

int main(int argc, char** argv) {
    if (argc > 1) return path_cli(argc, argv);

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
