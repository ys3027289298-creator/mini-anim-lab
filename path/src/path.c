
/* path: polyline arc-length table, constant-speed sampling, curvature and
   closest point queries.
 * Self-test rewrite: every original check is retained; an arc-length field
 * is rasterized to a PGM "path length map" and constant-speed samples are
 * dumped to CSV, both re-read from disk and asserted on. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>

#define MAX_PTS 32
#define ART_PATH_MAX 256

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

/* ---------------- visible artifact helpers ---------------- */

static const char* art_dir(void) {
    const char* d = getenv("ARTIFACT_DIR");
    return (d && *d) ? d : "artifacts";
}

static void art_mkdir(void) {
    struct stat st;
    const char* d = art_dir();
    if (stat(d, &st) != 0) mkdir(d, 0755);
}

static void art_fullpath(char* buf, size_t n, const char* name) {
    snprintf(buf, n, "%s/%s", art_dir(), name);
}

static void art_announce(const char* name) {
    char path[ART_PATH_MAX];
    struct stat st;
    art_fullpath(path, sizeof path, name);
    if (stat(path, &st) != 0) { printf("artifact %s MISSING\n", name); return; }
    printf("artifact %s %lld\n", name, (long long)st.st_size);
}

typedef struct {
    int w, h;
    unsigned char* px;
    float xmin, ymin, scale;
} Canvas;

static void cv_init(Canvas* cv, int w, int h, float xmin, float ymin, float scale) {
    cv->w = w; cv->h = h; cv->xmin = xmin; cv->ymin = ymin; cv->scale = scale;
    cv->px = (unsigned char*)malloc((size_t)w * h);
    memset(cv->px, 0, (size_t)w * h);
}

static void cv_map(const Canvas* cv, float wx, float wy, int* px, int* py) {
    *px = (int)lroundf((wx - cv->xmin) * cv->scale);
    *py = cv->h - 1 - (int)lroundf((wy - cv->ymin) * cv->scale);
}

static void cv_dot(Canvas* cv, float wx, float wy, unsigned char v) {
    int x, y; cv_map(cv, wx, wy, &x, &y);
    for (int dy = -1; dy <= 1; dy++)
        for (int dx = -1; dx <= 1; dx++) {
            int xx = x + dx, yy = y + dy;
            if (xx >= 0 && yy >= 0 && xx < cv->w && yy < cv->h)
                cv->px[yy * cv->w + xx] = v;
        }
}

static void cv_line(Canvas* cv, float x0, float y0, float x1, float y1, unsigned char v) {
    int X0, Y0, X1, Y1;
    cv_map(cv, x0, y0, &X0, &Y0);
    cv_map(cv, x1, y1, &X1, &Y1);
    int n = abs(X1 - X0) > abs(Y1 - Y0) ? abs(X1 - X0) : abs(Y1 - Y0);
    if (n == 0) n = 1;
    for (int i = 0; i <= n; i++) {
        int x = X0 + (X1 - X0) * i / n;
        int y = Y0 + (Y1 - Y0) * i / n;
        if (x >= 0 && y >= 0 && x < cv->w && y < cv->h)
            cv->px[y * cv->w + x] = v;
    }
}

static void cv_write(const Canvas* cv, const char* name) {
    char path[ART_PATH_MAX];
    art_fullpath(path, sizeof path, name);
    FILE* f = fopen(path, "wb");
    fprintf(f, "P5\n%d %d\n255\n", cv->w, cv->h);
    fwrite(cv->px, 1, (size_t)cv->w * cv->h, f);
    fclose(f);
}

static int cv_count_lit(const Canvas* cv, unsigned char min_v) {
    int n = 0;
    for (int i = 0; i < cv->w * cv->h; i++)
        if (cv->px[i] >= min_v) n++;
    return n;
}

static unsigned char cv_world_pixel(const Canvas* cv, float wx, float wy) {
    int x, y; cv_map(cv, wx, wy, &x, &y);
    if (x < 0 || y < 0 || x >= cv->w || y >= cv->h) return 0;
    return cv->px[y * cv->w + x];
}

static int pgm_disk_header_ok(const char* name, int expect_w, int expect_h) {
    char path[ART_PATH_MAX];
    art_fullpath(path, sizeof path, name);
    FILE* f = fopen(path, "rb");
    if (!f) return 0;
    char magic[8] = {0};
    int w = 0, h = 0, maxv = 0;
    int ok = fscanf(f, "%7s %d %d %d", magic, &w, &h, &maxv) == 4;
    fclose(f);
    return ok && strcmp(magic, "P5") == 0 && w == expect_w && h == expect_h && maxv == 255;
}

static int art_count_data_lines(const char* name) {
    char path[ART_PATH_MAX];
    art_fullpath(path, sizeof path, name);
    FILE* f = fopen(path, "r");
    if (!f) return -1;
    int lines = 0;
    char line[256];
    if (!fgets(line, sizeof line, f)) { fclose(f); return 0; }   /* header */
    while (fgets(line, sizeof line, f))
        if (line[0] != '\n' && line[0] != '\0') lines++;
    fclose(f);
    return lines;
}

/* ---------------- checks ---------------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("check %s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    art_mkdir();
    const char* csv = "path_samples.csv";
    const char* img = "path_length_map.pgm";

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
    /* exact magnitude for a unit 90deg bend: 2*cross/(la*lb*lc)
       = 2/(1*1*sqrt(2)) = sqrt(2); a doubled factor yields 2*sqrt(2) */
    check("curvature_magnitude", fabsf(path_curvature(&bend, 1) - 1.4142136f) < 1e-4f);

    /* visible artifact 1: constant-speed samples every 0.5 length units */
    {
        char path[ART_PATH_MAX];
        art_fullpath(path, sizeof path, csv);
        FILE* f = fopen(path, "w");
        fprintf(f, "s,x,y\n");
        for (int k = 0; k <= 14; k++) {
            float sd = 0.5f * (float)k;
            P2 q;
            if (path_sample(&pa, sd, &q, NULL) < 0) continue;
            fprintf(f, "%.4f,%.6f,%.6f\n", sd, q.x, q.y);
        }
        fclose(f);
    }
    art_announce(csv);

    /* visible artifact 2: arc-length field map (pixel value = arc length of
       the closest path point, scaled to 0..252); path itself burned at 255 */
    Canvas cv;
    cv_init(&cv, 64, 80, -0.5f, -0.5f, 16.0f);   /* x in [-.5,3.5), y in [-.5,4.5) */
    for (int gy = 0; gy < cv.h; gy++) {
        for (int gx = 0; gx < cv.w; gx++) {
            float wx = cv.xmin + (gx + 0.5f) / cv.scale;
            float wy = cv.ymin + (cv.h - 1 - gy + 0.5f) / cv.scale;
            float as = 0.0f;
            path_closest(&pa, (P2){ wx, wy }, &as);
            unsigned char v = (unsigned char)(as * 36.0f);   /* s in [0,7] -> [0,252] */
            cv.px[gy * cv.w + gx] = v;
        }
    }
    for (int i = 1; i < pa.n; i++)
        cv_line(&cv, pa.p[i - 1].x, pa.p[i - 1].y, pa.p[i].x, pa.p[i].y, 255);
    for (int i = 0; i < pa.n; i++)
        cv_dot(&cv, pa.p[i].x, pa.p[i].y, 255);
    cv_write(&cv, img);
    art_announce(img);

    /* ---- assertions against the visible artifacts on disk ---- */
    char path[ART_PATH_MAX];
    art_fullpath(path, sizeof path, csv);
    int rows = 0, mid_ok = 0, vert_ok = 0, end_ok = 0;
    FILE* f = fopen(path, "r");
    if (f) {
        char line[256];
        if (!fgets(line, sizeof line, f)) { fclose(f); f = NULL; }   /* header consumed */
        float sd, x, y;
        while (f && fgets(line, sizeof line, f)) {
            if (sscanf(line, "%f,%f,%f", &sd, &x, &y) == 3) {
                rows++;
                if (fabsf(sd - 1.5f) < 1e-6f && fabsf(x - 1.5f) < 1e-4f && fabsf(y) < 1e-4f) mid_ok = 1;
                if (fabsf(sd - 5.0f) < 1e-6f && fabsf(x - 3.0f) < 1e-4f && fabsf(y - 2.0f) < 1e-4f) vert_ok = 1;
                if (fabsf(sd - 7.0f) < 1e-6f && fabsf(y - 4.0f) < 1e-4f) end_ok = 1;
            }
        }
        if (f) fclose(f);
    }
    check("art_csv_rows", art_count_data_lines(csv) == 15);
    check("art_sample_mid_disk", mid_ok);
    check("art_sample_vert_disk", vert_ok);
    check("art_sample_end_disk", end_ok);
    check("art_map_header", pgm_disk_header_ok(img, 64, 80));
    check("art_map_path_lit_band", cv_count_lit(&cv, 255) >= 40 && cv_count_lit(&cv, 255) <= 220);
    /* arc-length shading sampled OFF the burned 255 path line:
       (1,1.5) projects to horizontal seg s=1 -> ~36;
       (3.25,2) projects to vertical seg s=5 -> 180 */
    unsigned char v_h = cv_world_pixel(&cv, 1.0f, 1.5f);
    unsigned char v_v = cv_world_pixel(&cv, 3.25f, 2.0f);
    check("art_map_horiz_value", v_h >= 20 && v_h <= 60);
    check("art_map_vert_value", v_v >= 150 && v_v <= 210);
    check("art_map_monotonic", v_v > v_h);

    free(cv.px);
    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
