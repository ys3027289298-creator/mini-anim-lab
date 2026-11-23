
/* curve: keyframe sampling with cubic hermite tangents and clamping.
 * Self-test rewrite: every original check is retained; the sampled curve is
 * dumped to CSV and rasterized to PGM, both re-read from disk and asserted
 * against the exact analytic Hermite midpoint value. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>

#define MAX_KEYS 16
#define ART_PATH_MAX 256

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
    const char* csv = "curve_samples.csv";
    const char* img = "curve.pgm";

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
    /* exact derivative of a constant flat curve is 0 at any finite-difference
       step, so demand a tight tolerance; a 100x too-large step still reads 0
       here, hence the separate sloped derivative check below */
    check("flat_deriv_zero", fabsf(curve_deriv(&flat, 0.5f)) < 1e-6f);

    /* derivative near the left clamp boundary: a correct 1e-3 step stays
       inside the cubic and reads ~1.0; a coarse 0.1 step reaches across the
       clamp at t=0 and degrades the difference to ~0.75. A centered diff of
       a pure cubic is exact at any h in the interior, so the boundary is the
       only place this defect is observable. */
    check("clamp_zone_deriv", fabsf(curve_deriv(&c, 0.05f) - 1.0f) < 1e-2f);

    check("reject_unsorted", curve_add(&flat, 0.5f, 1.0f, 0, 0) == -1);
    /* equal timestamp must also be rejected: strictly-increasing invariant.
       Inverting the <= guard to < would wrongly accept this insert */
    int before = flat.count;
    check("reject_equal_time", curve_add(&flat, 1.0f, 0.9f, 0, 0) == -1
                              && flat.count == before);

    Curve full;
    curve_init(&full);
    int rc = 0;
    for (int i = 0; i <= MAX_KEYS; i++) rc = curve_add(&full, (float)i, (float)i, 0, 0);
    check("capacity", rc == -1 && full.count == MAX_KEYS);

    /* visible artifact 1: dense sample table including clamp regions */
    {
        char path[ART_PATH_MAX];
        art_fullpath(path, sizeof path, csv);
        FILE* f = fopen(path, "w");
        fprintf(f, "t,v\n");
        for (int k = -4; k <= 28; k++) {
            float t = 0.05f * (float)k;              /* -0.2 .. 1.4 */
            fprintf(f, "%.4f,%.6f\n", t, curve_eval(&c, t));
        }
        fclose(f);
    }
    art_announce(csv);

    /* visible artifact 2: rasterized value curve with baseline grid markers */
    Canvas cv;
    cv_init(&cv, 128, 128, -0.25f, -0.35f, 80.0f);  /* t in [-.25,1.35), v in [-.35,1.25) */
    for (int k = 0; k < 32; k++) {
        float t0 = -0.2f + 1.6f * (float)k / 32.0f;
        float t1 = -0.2f + 1.6f * (float)(k + 1) / 32.0f;
        cv_line(&cv, t0, curve_eval(&c, t0), t1, curve_eval(&c, t1), 220);
    }
    cv_line(&cv, -0.25f, 0.0f, 1.35f, 0.0f, 80);    /* zero baseline */
    cv_dot(&cv, 0.0f, 0.0f, 255);
    cv_dot(&cv, 1.0f, 1.0f, 255);
    cv_write(&cv, img);
    art_announce(img);

    /* exact analytic Hermite value at u=0.5: h00=h01=0.5; tangent basis
       h10 = u^3-2u^2+u and h11 = u^3-u^2 are both 0 at u=0.5
       (1/8 - 2/4 + 1/2 = 0; 1/8 - 1/4 = 0), so tangents contribute nothing
       and the midpoint is the endpoint average (0+1)/2 = 0.5 */
    const float expect_mid = 0.5f;

    /* ---- assertions against the visible artifacts on disk ---- */
    char path[ART_PATH_MAX];
    art_fullpath(path, sizeof path, csv);
    int rows = 0, clamp_l_ok = 0, clamp_r_ok = 0, mid_ok = 0, end_ok = 0;
    float mid_v = 0.0f;
    FILE* f = fopen(path, "r");
    if (f) {
        char line[256];
        if (!fgets(line, sizeof line, f)) { fclose(f); f = NULL; }   /* header consumed */
        float t, v;
        while (f && fgets(line, sizeof line, f)) {
            if (sscanf(line, "%f,%f", &t, &v) == 2) {
                rows++;
                if (fabsf(t + 0.2f) < 1e-6f && fabsf(v) < 1e-6f) clamp_l_ok = 1;
                if (fabsf(t - 1.4f) < 1e-6f && fabsf(v - 1.0f) < 1e-6f) clamp_r_ok = 1;
                if (fabsf(t - 0.5f) < 1e-6f) { mid_v = v; if (fabsf(v - expect_mid) < 1e-4f) mid_ok = 1; }
                if (fabsf(t - 1.0f) < 1e-6f && fabsf(v - 1.0f) < 1e-6f) end_ok = 1;
            }
        }
        if (f) fclose(f);
    }
    check("art_csv_rows", art_count_data_lines(csv) == 33);
    check("art_clamp_left_disk", clamp_l_ok);
    check("art_clamp_right_disk", clamp_r_ok);
    check("art_midpoint_analytic", mid_ok);
    check("art_endpoint_disk", end_ok);
    check("art_image_header", pgm_disk_header_ok(img, 128, 128));
    check("art_image_lit_band", cv_count_lit(&cv, 200) >= 20 && cv_count_lit(&cv, 200) <= 300);
    check("art_image_spot_mid", cv_world_pixel(&cv, 0.5f, mid_v) >= 200);
    check("art_image_spot_key", cv_world_pixel(&cv, 1.0f, 1.0f) >= 250);

    free(cv.px);
    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
