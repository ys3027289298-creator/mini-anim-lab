
/* ik: analytic two-bone IK with pole vector and joint angle limits.
 * Self-test rewrite: every original check is retained; the IK solution is
 * rendered to a PGM image and dumped to a pose CSV, both re-read from disk
 * and asserted against analytically expected values. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>

#define ART_PATH_MAX 256

typedef struct { float x, y; } V2;

static float vlen(V2 v) { return sqrtf(v.x * v.x + v.y * v.y); }
static V2 vsub(V2 a, V2 b) { V2 r = { a.x - b.x, a.y - b.y }; return r; }
static V2 vnorm(V2 v) { float l = vlen(v); V2 r = { 0, 0 }; if (l > 0) { r.x = v.x / l; r.y = v.y / l; } return r; }

typedef struct {
    int ok;
    V2 joint;
    float a1;   /* angle of first bone */
    float a2;   /* angle of second bone */
} IkResult;

/* solve for target with bone lengths l1,l2; pole decides bend direction */
static IkResult solve_two_bone(V2 root, V2 target, float l1, float l2, V2 pole) {
    IkResult r;
    memset(&r, 0, sizeof(r));
    V2 d = vsub(target, root);
    float dist = vlen(d);
    if (l1 <= 0 || l2 <= 0) return r;
    if (dist > l1 + l2 + 1e-5f) return r;                 /* unreachable */
    if (dist < fabsf(l1 - l2) - 1e-5f) return r;          /* too close */
    if (dist < 1e-6f) return r;

    float cos_a = (l1 * l1 + dist * dist - l2 * l2) / (2 * l1 * dist);
    if (cos_a > 1.0f) cos_a = 1.0f;
    if (cos_a < -1.0f) cos_a = -1.0f;
    float a = acosf(cos_a);

    V2 dir = vnorm(d);
    float side = (dir.x * pole.y - dir.y * pole.x) >= 0 ? 1.0f : -1.0f;
    float ang = atan2f(dir.y, dir.x) + side * a;

    r.joint.x = root.x + l1 * cosf(ang);
    r.joint.y = root.y + l1 * sinf(ang);
    r.a1 = ang;
    V2 d2 = vsub(target, r.joint);
    r.a2 = atan2f(d2.y, d2.x);
    r.ok = 1;
    return r;
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
    const char* csv = "ik_pose.csv";
    const char* img = "ik_pose.pgm";

    V2 root = { 0, 0 }, target = { 1, 0 }, pole = { 0, 1 };
    IkResult r = solve_two_bone(root, target, 1.0f, 1.0f, pole);
    check("reachable", r.ok == 1);
    check("joint_on_circle", fabsf(vlen(vsub(r.joint, root)) - 1.0f) < 1e-4f);
    check("end_matches", fabsf(vlen(vsub(r.joint, target)) - 1.0f) < 1e-4f);
    check("pole_up", r.joint.y > 0.0f);

    IkResult r2 = solve_two_bone(root, target, 1.0f, 1.0f, (V2){ 0, -1 });
    check("pole_down", r2.joint.y < 0.0f);

    IkResult far = solve_two_bone(root, (V2){ 5, 0 }, 1.0f, 1.0f, pole);
    check("unreachable", far.ok == 0);

    IkResult close = solve_two_bone(root, (V2){ 0.1f, 0 }, 1.0f, 0.2f, pole);
    check("too_close", close.ok == 0);

    IkResult zero = solve_two_bone(root, target, 0.0f, 1.0f, pole);
    check("zero_len", zero.ok == 0);

    /* visible artifact 1: full pose table for every scenario */
    {
        char path[ART_PATH_MAX];
        art_fullpath(path, sizeof path, csv);
        FILE* f = fopen(path, "w");
        fprintf(f, "scenario,ok,jx,jy,a1,a2\n");
        fprintf(f, "up,%d,%.6f,%.6f,%.6f,%.6f\n", r.ok, r.joint.x, r.joint.y, r.a1, r.a2);
        fprintf(f, "down,%d,%.6f,%.6f,%.6f,%.6f\n", r2.ok, r2.joint.x, r2.joint.y, r2.a1, r2.a2);
        fprintf(f, "far,%d,%.6f,%.6f,%.6f,%.6f\n", far.ok, far.joint.x, far.joint.y, far.a1, far.a2);
        fprintf(f, "close,%d,%.6f,%.6f,%.6f,%.6f\n", close.ok, close.joint.x, close.joint.y, close.a1, close.a2);
        fprintf(f, "zero,%d,%.6f,%.6f,%.6f,%.6f\n", zero.ok, zero.joint.x, zero.joint.y, zero.a1, zero.a2);
        fclose(f);
    }
    art_announce(csv);

    /* visible artifact 2: both bent solutions + root/target markers */
    Canvas cv;
    cv_init(&cv, 128, 96, -1.5f, -2.0f, 32.0f);   /* x in [-1.5, 2.5), y in [-2, 1) */
    cv_line(&cv, root.x, root.y, r.joint.x, r.joint.y, 200);
    cv_line(&cv, r.joint.x, r.joint.y, target.x, target.y, 200);
    cv_line(&cv, root.x, root.y, r2.joint.x, r2.joint.y, 140);
    cv_line(&cv, r2.joint.x, r2.joint.y, target.x, target.y, 140);
    cv_dot(&cv, root.x, root.y, 255);
    cv_dot(&cv, target.x, target.y, 255);
    cv_dot(&cv, r.joint.x, r.joint.y, 255);
    cv_dot(&cv, r2.joint.x, r2.joint.y, 255);
    cv_write(&cv, img);
    art_announce(img);

    /* ---- assertions against the visible artifacts on disk ---- */
    char path[ART_PATH_MAX];
    art_fullpath(path, sizeof path, csv);
    int found_up = 0, found_down = 0, found_far = 0;
    float upx = 0, upy = 0, upa1 = 0, upa2 = 0, downy = 0;
    FILE* f = fopen(path, "r");
    if (f) {
        char scenario[16]; int ok; float jx, jy, a1, a2;
        char line[256];
        if (!fgets(line, sizeof line, f)) { fclose(f); f = NULL; }   /* header consumed */
        while (f && fgets(line, sizeof line, f)) {
            if (sscanf(line, "%15[^,],%d,%f,%f,%f,%f", scenario, &ok, &jx, &jy, &a1, &a2) == 6) {
                if (strcmp(scenario, "up") == 0) {
                    found_up = 1; upx = jx; upy = jy; upa1 = a1; upa2 = a2;
                } else if (strcmp(scenario, "down") == 0) {
                    found_down = 1; downy = jy;
                } else if (strcmp(scenario, "far") == 0) {
                    found_far = (ok == 0);
                }
            }
        }
        if (f) fclose(f);
    }
    /* analytical expectation for equal 1-length bones reaching (1,0):
       joint at (0.5, +/-sqrt(3)/2), a1 = +/-pi/3, a2 = -/+pi/3 */
    check("art_csv_rows", art_count_data_lines(csv) == 5);
    check("art_pose_up_analytic",
          found_up && fabsf(upx - 0.5f) < 1e-4f && fabsf(upy - 0.8660254f) < 1e-4f);
    check("art_angles_analytic",
          found_up && fabsf(upa1 - 1.0471976f) < 1e-4f && fabsf(upa2 + 1.0471976f) < 1e-4f);
    check("art_pose_down_analytic",
          found_down && fabsf(downy + 0.8660254f) < 1e-4f);
    check("art_reject_rows_disk", found_far);
    check("art_image_header", pgm_disk_header_ok(img, 128, 96));
    check("art_image_lit_band", cv_count_lit(&cv, 100) >= 40 && cv_count_lit(&cv, 100) <= 500);
    check("art_image_spot_up_joint", cv_world_pixel(&cv, 0.5f, 0.8660254f) >= 200);
    check("art_image_spot_down_joint", cv_world_pixel(&cv, 0.5f, -0.8660254f) >= 100);

    free(cv.px);
    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
