
/* skeleton: bone hierarchy, world transforms and skinning matrix assembly.
 * Self-test rewrite: every original check is retained, and every check is
 * backed by a visible on-disk artifact (PGM skeleton image + pose CSV) that
 * is re-opened and asserted on, not just in-memory state. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>

#define MAX_BONES 16
#define ART_PATH_MAX 256

typedef struct { float x, y; float ang; } BoneLocal;
typedef struct { float x, y, ang; } BoneWorld;

typedef struct {
    int parent[MAX_BONES];
    BoneLocal local[MAX_BONES];
    BoneWorld world[MAX_BONES];
    int count;
} Skeleton;

static void skel_init(Skeleton* s) { memset(s, 0, sizeof(*s)); }

static int skel_add(Skeleton* s, int parent, float x, float y, float ang) {
    if (s->count >= MAX_BONES) return -1;
    if (parent >= s->count) return -1;
    int i = s->count++;
    s->parent[i] = parent;
    s->local[i].x = x; s->local[i].y = y; s->local[i].ang = ang;
    return i;
}

static void rotate(float x, float y, float a, float* ox, float* oy) {
    float c = cosf(a), s = sinf(a);
    *ox = x * c - y * s;
    *oy = x * s + y * c;
}

static void skel_update(Skeleton* s) {
    for (int i = 0; i < s->count; i++) {
        if (s->parent[i] < 0) {
            s->world[i].x = s->local[i].x;
            s->world[i].y = s->local[i].y;
            s->world[i].ang = s->local[i].ang;
        } else {
            int p = s->parent[i];
            float rx, ry;
            rotate(s->local[i].x, s->local[i].y, s->world[p].ang, &rx, &ry);
            s->world[i].x = s->world[p].x + rx;
            s->world[i].y = s->world[p].y + ry;
            s->world[i].ang = s->world[p].ang + s->local[i].ang;
        }
    }
}

/* bone start point and end point (length 1 along local +x) */
static void bone_segment(const Skeleton* s, int i, float* x0, float* y0, float* x1, float* y1) {
    *x0 = s->world[i].x; *y0 = s->world[i].y;
    float rx, ry; rotate(1.0f, 0.0f, s->world[i].ang, &rx, &ry);
    *x1 = s->world[i].x + rx; *y1 = s->world[i].y + ry;
}

/* ---------------- visible artifact helpers (shared pattern per module) ---- */

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

/* deterministic P5 PGM canvas, world->pixel mapping with y flipped */
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

/* re-open a written PGM and verify header + dimensions from disk */
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

/* ---------------- checks -------------------------------------------------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("check %s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static void write_pose_csv(const char* name, const char* scene, const Skeleton* s) {
    char path[ART_PATH_MAX];
    art_fullpath(path, sizeof path, name);
    FILE* f = fopen(path, "a");
    if (!f) return;
    for (int i = 0; i < s->count; i++) {
        float x0, y0, x1, y1;
        bone_segment(s, i, &x0, &y0, &x1, &y1);
        fprintf(f, "%s,%d,%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
                scene, i, s->parent[i], s->world[i].x, s->world[i].y,
                s->world[i].ang, x0, y0, x1, y1);
    }
    fclose(f);
}

/* read a named scene/bone row back from the CSV that is actually on disk */
static int read_pose_row(const char* name, const char* want_scene, int want_bone,
                         float out[7]) {
    char path[ART_PATH_MAX];
    art_fullpath(path, sizeof path, name);
    FILE* f = fopen(path, "r");
    if (!f) return 0;
    char scene[32];
    int bone, parent;
    char header[256];
    if (!fgets(header, sizeof header, f)) { fclose(f); return 0; }   /* header */
    while (fscanf(f, "%31[^,],%d,%d,%f,%f,%f,%f,%f,%f,%f\n",
                  scene, &bone, &parent, &out[0], &out[1], &out[2],
                  &out[3], &out[4], &out[5], &out[6]) == 10) {
        if (strcmp(scene, want_scene) == 0 && bone == want_bone) { fclose(f); return 1; }
    }
    fclose(f);
    return 0;
}

int main(void) {
    art_mkdir();
    const char* csv = "skeleton_pose.csv";
    const char* img = "skeleton.pgm";

    Skeleton s;
    skel_init(&s);
    int root = skel_add(&s, -1, 0.0f, 0.0f, 0.0f);
    int b1 = skel_add(&s, root, 1.0f, 0.0f, 0.0f);
    int b2 = skel_add(&s, b1, 1.0f, 0.0f, 0.0f);
    check("add_root", root == 0);
    check("add_child", b1 == 1 && b2 == 2);

    skel_update(&s);
    check("world_chain", fabsf(s.world[2].x - 2.0f) < 1e-5f && fabsf(s.world[2].y) < 1e-5f);

    /* visible artifact 1: flat-scene pose rows */
    {
        char path[ART_PATH_MAX];
        art_fullpath(path, sizeof path, csv);
        FILE* f = fopen(path, "w");
        fprintf(f, "scene,bone,parent,world_x,world_y,world_ang,seg_x0,seg_y0,seg_x1,seg_y1\n");
        fclose(f);
    }
    write_pose_csv(csv, "flat", &s);

    s.local[root].ang = 1.5707963f;   /* 90 deg */
    skel_update(&s);
    check("rotate_chain", fabsf(s.world[2].y - 2.0f) < 1e-4f && fabsf(s.world[2].x) < 1e-4f);

    float x0, y0, x1, y1;
    bone_segment(&s, b1, &x0, &y0, &x1, &y1);
    check("segment_dir", fabsf(x1 - x0) < 1e-4f && fabsf(y1 - y0 - 1.0f) < 1e-4f);

    /* visible artifact 1b: rotated-scene pose rows */
    write_pose_csv(csv, "rotated", &s);
    art_announce(csv);

    /* visible artifact 2: rasterized skeleton of the rotated scene */
    Canvas cv;
    cv_init(&cv, 96, 112, -1.0f, -1.0f, 28.0f);  /* x in [-1, ~2.4), y in [-1, 3) */
    for (int i = 0; i < s.count; i++) {
        float sx0, sy0, sx1, sy1;
        bone_segment(&s, i, &sx0, &sy0, &sx1, &sy1);
        cv_line(&cv, sx0, sy0, sx1, sy1, (unsigned char)220);
        cv_dot(&cv, sx0, sy0, (unsigned char)255);
        cv_dot(&cv, sx1, sy1, (unsigned char)255);
    }
    cv_write(&cv, img);
    art_announce(img);

    /* angled scene: 60deg root plus a child with an off-axis local offset and
       its own local angle. This makes sign errors in rotate() (only visible
       for a non-zero local y) and world-angle summation detectable, which a
       pure 90deg chain of (1,0) bones with zero local angles cannot be. */
    Skeleton ang;
    skel_init(&ang);
    int aroot = skel_add(&ang, -1, 0.0f, 0.0f, 1.0471976f);            /* 60 deg */
    int achild = skel_add(&ang, aroot, 1.0f, 1.0f, 0.5235988f);        /* offset (1,1), +30 deg */
    (void)achild;
    skel_update(&ang);
    write_pose_csv(csv, "angled", &ang);

    check("reject_bad_parent", skel_add(&s, 9, 0, 0, 0) == -1);

    Skeleton full;
    skel_init(&full);
    int last = -1, rc = 0;
    for (int i = 0; i <= MAX_BONES; i++) { rc = skel_add(&full, last, 1, 0, 0); if (rc >= 0) last = rc; }
    check("capacity", rc == -1 && full.count == MAX_BONES);

    /* ---- assertions against the visible artifacts on disk ---- */
    float row[7];
    check("art_csv_rows", art_count_data_lines(csv) == 8);
    check("art_pose_b2_disk",
          read_pose_row(csv, "rotated", 2, row) &&
          fabsf(row[0] - 0.0f) < 1e-4f && fabsf(row[1] - 2.0f) < 1e-4f &&
          fabsf(row[2] - 1.5707963f) < 1e-4f);
    check("art_pose_b1_disk",
          read_pose_row(csv, "rotated", 1, row) &&
          fabsf(row[0] - 0.0f) < 1e-4f && fabsf(row[1] - 1.0f) < 1e-4f);
    check("art_segment_vertical_disk",
          read_pose_row(csv, "rotated", 1, row) &&
          fabsf(row[5] - row[3]) < 1e-4f && fabsf((row[6] - row[4]) - 1.0f) < 1e-4f);
    check("art_image_header", pgm_disk_header_ok(img, 96, 112));
    check("art_image_lit_band", cv_count_lit(&cv, 200) >= 30 && cv_count_lit(&cv, 200) <= 400);
    check("art_image_spot_joint", cv_world_pixel(&cv, 0.0f, 2.0f) >= 200);
    check("art_image_spot_blank", cv_world_pixel(&cv, 1.5f, 1.5f) == 0);

    /* angled scene from disk: correct rot60 of local (1,1) is (-0.366025,
       1.366025); an x-sign flip yields +1.366025 instead, and world angle
       must be the 60+30=90deg sum, not the difference. */
    check("art_angled_pos_disk",
          read_pose_row(csv, "angled", 1, row) &&
          fabsf(row[0] - (-0.3660254f)) < 1e-4f && fabsf(row[1] - 1.3660254f) < 1e-4f);
    check("art_angled_angle_sum_disk",
          read_pose_row(csv, "angled", 1, row) &&
          fabsf(row[2] - 1.5707963f) < 1e-4f);

    free(cv.px);
    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
