
/* skeleton: bone hierarchy, world transforms and skinning matrix assembly. */
#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "../../shared/pose_export.h"

#define MAX_BONES 16

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

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

/* ---------- unified pose export (new format, see shared/pose_export.h) ----- */

/* canonical pose: 3-bone chain with root rotated 90 degrees */
static void skeleton_pose_build(Skeleton* s) {
    skel_init(s);
    int root = skel_add(s, -1, 0.0f, 0.0f, 1.5707963f);
    int b1 = skel_add(s, root, 1.0f, 0.0f, 0.0f);
    int b2 = skel_add(s, b1, 1.0f, 0.0f, 0.0f);
    (void)b2;
    skel_update(s);
}

static int skeleton_pose_export(const char* outdir, int argc, char** argv) {
    if (pe_mkdir_p(outdir) != 0) return 2;

    /* timing tier 1: pose numeric computation, 1000 repeats */
    Skeleton s;
    long t0 = pe_now_us();
    for (int rep = 0; rep < 1000; rep++) {
        skeleton_pose_build(&s);
        /* touch results so the loop cannot be optimized away */
        if (s.world[s.count - 1].x + s.world[s.count - 1].y > 1e30f) return 2;
    }
    long t1 = pe_now_us();
    long compute_us = t1 - t0;

    /* timing tier 2: image rendering, single pass */
    PeCanvas cv;
    pe_canvas_init(&cv, -0.5f, -0.5f, 2.5f, 2.5f);
    long t2 = pe_now_us();
    for (int i = 0; i < s.count; i++) {
        float x0, y0, x1, y1;
        bone_segment(&s, i, &x0, &y0, &x1, &y1);
        pe_line(&cv, x0, y0, x1, y1, 40, 90, 180);
    }
    for (int i = 0; i < s.count; i++)
        pe_circle(&cv, s.world[i].x, s.world[i].y, 0.06f, 200, 60, 60);
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
    fprintf(f, "  \"module\": \"skeleton\",\n");
    fprintf(f, "  \"format_version\": %d,\n", POSE_FORMAT_VERSION);
    fprintf(f, "  \"pose\": {\n");
    fprintf(f, "    \"bone_count\": %d,\n", s.count);
    fprintf(f, "    \"bones\": [\n");
    for (int i = 0; i < s.count; i++) {
        float x0, y0, x1, y1;
        bone_segment(&s, i, &x0, &y0, &x1, &y1);
        fprintf(f,
                "      {\"index\": %d, \"parent\": %d, \"x\": %.6f, \"y\": %.6f, "
                "\"ang\": %.6f, \"end_x\": %.6f, \"end_y\": %.6f}%s\n",
                i, s.parent[i], s.world[i].x, s.world[i].y, s.world[i].ang,
                x1, y1, i == s.count - 1 ? "" : ",");
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
    /* exact command line this binary was invoked with, plus the output dir */
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
        return skeleton_pose_export(argv[2], argc, argv);
    }

    Skeleton s;
    skel_init(&s);
    int root = skel_add(&s, -1, 0.0f, 0.0f, 0.0f);
    int b1 = skel_add(&s, root, 1.0f, 0.0f, 0.0f);
    int b2 = skel_add(&s, b1, 1.0f, 0.0f, 0.0f);
    check("add_root", root == 0);
    check("add_child", b1 == 1 && b2 == 2);

    skel_update(&s);
    check("world_chain", fabsf(s.world[2].x - 2.0f) < 1e-5f && fabsf(s.world[2].y) < 1e-5f);

    s.local[root].ang = 1.5707963f;   /* 90 deg */
    skel_update(&s);
    check("rotate_chain", fabsf(s.world[2].y - 2.0f) < 1e-4f && fabsf(s.world[2].x) < 1e-4f);

    float x0, y0, x1, y1;
    bone_segment(&s, b1, &x0, &y0, &x1, &y1);
    check("segment_dir", fabsf(x1 - x0) < 1e-4f && fabsf(y1 - y0 - 1.0f) < 1e-4f);

    check("reject_bad_parent", skel_add(&s, 9, 0, 0, 0) == -1);

    Skeleton full;
    skel_init(&full);
    int last = -1, rc = 0;
    for (int i = 0; i <= MAX_BONES; i++) { rc = skel_add(&full, last, 1, 0, 0); if (rc >= 0) last = rc; }
    check("capacity", rc == -1 && full.count == MAX_BONES);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
