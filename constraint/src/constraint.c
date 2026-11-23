
/* constraint: verlet particle chain with distance constraints, iteration
   count and damping; used as the toy model for cloth/rope style solvers.
 * Self-test rewrite: every original check is retained; after stepping, the
 * full particle chain is dumped to CSV and rasterized to PGM, both re-read
 * from disk and asserted on. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>

#define MAX_P 16
#define ART_PATH_MAX 256

typedef struct { float x, y, px, py; float inv_mass; int pinned; } Particle;

typedef struct {
    Particle p[MAX_P];
    int n;
    float rest[MAX_P];      /* rest length between i and i+1 */
    int iterations;
    float damping;
} Chain;

static void chain_init(Chain* c, int n, float seg, float x0, float y0) {
    memset(c, 0, sizeof(*c));
    c->n = n;
    c->iterations = 8;
    c->damping = 0.99f;
    for (int i = 0; i < n; i++) {
        c->p[i].x = x0 + seg * i; c->p[i].y = y0;
        c->p[i].px = c->p[i].x; c->p[i].py = c->p[i].y;
        c->p[i].inv_mass = 1.0f;
        c->p[i].pinned = 0;
        c->rest[i] = seg;
    }
    c->p[0].pinned = 1;
    c->p[0].inv_mass = 0.0f;
}

static void chain_step(Chain* c, float dt, float gx, float gy) {
    for (int i = 0; i < c->n; i++) {
        Particle* p = &c->p[i];
        if (p->pinned) continue;
        float vx = (p->x - p->px) * c->damping;
        float vy = (p->y - p->py) * c->damping;
        p->px = p->x; p->py = p->y;
        p->x += vx + gx * dt * dt;
        p->y += vy + gy * dt * dt;
    }
    for (int it = 0; it < c->iterations; it++) {
        for (int i = 0; i < c->n - 1; i++) {
            Particle* a = &c->p[i];
            Particle* b = &c->p[i + 1];
            float dx = b->x - a->x, dy = b->y - a->y;
            float d = sqrtf(dx * dx + dy * dy);
            if (d < 1e-6f) continue;
            float diff = (d - c->rest[i]) / d;
            float wa = a->inv_mass, wb = b->inv_mass;
            float wsum = wa + wb;
            if (wsum <= 0) continue;
            a->x += dx * diff * (wa / wsum);
            a->y += dy * diff * (wa / wsum);
            b->x -= dx * diff * (wb / wsum);
            b->y -= dy * diff * (wb / wsum);
        }
    }
}

static float chain_seg_error(const Chain* c) {
    float worst = 0.0f;
    for (int i = 0; i < c->n - 1; i++) {
        float dx = c->p[i + 1].x - c->p[i].x;
        float dy = c->p[i + 1].y - c->p[i].y;
        float d = sqrtf(dx * dx + dy * dy);
        float e = fabsf(d - c->rest[i]);
        if (e > worst) worst = e;
    }
    return worst;
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

static void write_chain_csv(const char* name, const Chain* c, float worst) {
    char path[ART_PATH_MAX];
    art_fullpath(path, sizeof path, name);
    FILE* f = fopen(path, "w");
    fprintf(f, "i,x,y,pinned,inv_mass,seg_error_worst\n");
    for (int i = 0; i < c->n; i++)
        fprintf(f, "%d,%.6f,%.6f,%d,%.3f,%.6f\n",
                i, c->p[i].x, c->p[i].y, c->p[i].pinned, c->p[i].inv_mass, worst);
    fclose(f);
}

int main(void) {
    art_mkdir();
    const char* csv = "constraint_chain.csv";
    const char* img = "constraint.pgm";

    Chain c;
    chain_init(&c, 8, 1.0f, 0.0f, 0.0f);
    check("init_len", c.n == 8 && fabsf(chain_seg_error(&c)) < 1e-5f);

    /* let it fall for a while */
    for (int i = 0; i < 200; i++) chain_step(&c, 0.016f, 0.0f, -9.8f);
    check("stays_in_stretch", chain_seg_error(&c) < 0.02f);
    check("pinned_anchor", fabsf(c.p[0].x) < 1e-6f && fabsf(c.p[0].y) < 1e-6f);
    check("falls_down", c.p[7].y < -1.0f);

    /* more iterations must not increase error */
    Chain a, b;
    chain_init(&a, 8, 1.0f, 0.0f, 0.0f);
    chain_init(&b, 8, 1.0f, 0.0f, 0.0f);
    a.iterations = 1; b.iterations = 16;
    for (int i = 0; i < 100; i++) { chain_step(&a, 0.016f, 0.0f, -9.8f); chain_step(&b, 0.016f, 0.0f, -9.8f); }
    check("more_iters_better", chain_seg_error(&b) <= chain_seg_error(&a) + 1e-6f);

    Chain cap;
    chain_init(&cap, 16, 0.5f, 0.0f, 0.0f);
    check("capacity_ok", cap.n == 16);

    /* visible artifacts from the 200-step fallen chain */
    float worst = chain_seg_error(&c);
    write_chain_csv(csv, &c, worst);
    art_announce(csv);

    Canvas cv;
    cv_init(&cv, 160, 240, -3.5f, -7.0f, 32.0f);  /* x in [-3.5,1.5), y in [-7,0.5) */
    for (int i = 0; i < c.n - 1; i++)
        cv_line(&cv, c.p[i].x, c.p[i].y, c.p[i + 1].x, c.p[i + 1].y, 200);
    for (int i = 0; i < c.n; i++)
        cv_dot(&cv, c.p[i].x, c.p[i].y, 255);
    cv_write(&cv, img);
    art_announce(img);

    /* ---- assertions against the visible artifacts on disk ---- */
    char path[ART_PATH_MAX];
    art_fullpath(path, sizeof path, csv);
    int rows = 0, disk_pinned = 0, disk_tail_ok = 0, disk_err_ok = 0;
    float tail_y = 0.0f, disk_worst = -1.0f;
    FILE* f = fopen(path, "r");
    if (f) {
        char line[256];
        if (!fgets(line, sizeof line, f)) { fclose(f); f = NULL; }   /* header consumed */
        int i, pinned; float x, y, inv_mass, e;
        while (f && fgets(line, sizeof line, f)) {
            if (sscanf(line, "%d,%f,%f,%d,%f,%f", &i, &x, &y, &pinned, &inv_mass, &e) == 6) {
                rows++;
                disk_worst = e;
                if (i == 0 && pinned == 1 && fabsf(x) < 1e-6f && fabsf(y) < 1e-6f) disk_pinned = 1;
                if (i == 7) { tail_y = y; if (y < -1.0f) disk_tail_ok = 1; }
            }
        }
        if (f) fclose(f);
    }
    disk_err_ok = disk_worst >= 0.0f && disk_worst < 0.02f;
    check("art_csv_rows", rows == 8 && art_count_data_lines(csv) == 8);
    check("art_anchor_disk", disk_pinned);
    check("art_tail_disk", disk_tail_ok);
    check("art_seg_error_disk", disk_err_ok);
    check("art_image_header", pgm_disk_header_ok(img, 160, 240));
    check("art_image_lit_band", cv_count_lit(&cv, 200) >= 40 && cv_count_lit(&cv, 200) <= 400);
    check("art_image_spot_anchor", cv_world_pixel(&cv, 0.0f, 0.0f) >= 200);
    check("art_image_spot_tail", cv_world_pixel(&cv, c.p[7].x, c.p[7].y) >= 200);
    check("art_image_spot_blank", cv_world_pixel(&cv, 1.5f, -0.5f) == 0);
    (void)tail_y;

    free(cv.px);
    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
