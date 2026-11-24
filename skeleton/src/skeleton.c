
/* skeleton: bone hierarchy, world transforms and skinning matrix assembly. */
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "../../shared/pose_format.h"

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

/* ===================== shared pose format (new) ==========================
 *
 * One "bone" record per bone, carrying the hierarchy (parent + local pose)
 * and the resulting world transform:
 *
 *   bone index=0 parent=-1 lx=0 ly=0 lang=1.57079637 wx=0 wy=0 ang=1.57079637
 *
 * Reading the file back restores both arrays, so a roundtrip compares the
 * whole internal Skeleton, not just the text.
 */

static void skel_pose_build(Skeleton* s) {
    skel_init(s);
    int root = skel_add(s, -1, 0.0f, 0.0f, 1.5707963f);  /* 90 deg root */
    int b1 = skel_add(s, root, 1.0f, 0.0f, 0.0f);
    skel_add(s, b1, 1.0f, 0.0f, 0.0f);
    skel_update(s);
}

static int skel_export(const Skeleton* s, const char* file) {
    FILE* f = fopen(file, "w");
    if (!f) return -1;
    pf_write_header(f, "skeleton", s->count);
    for (int i = 0; i < s->count; i++) {
        pf_record_begin(f, "bone");
        pf_field_int(f, "index", i);
        pf_field_int(f, "parent", s->parent[i]);
        pf_field_float(f, "lx", s->local[i].x);
        pf_field_float(f, "ly", s->local[i].y);
        pf_field_float(f, "lang", s->local[i].ang);
        pf_field_float(f, "wx", s->world[i].x);
        pf_field_float(f, "wy", s->world[i].y);
        pf_field_float(f, "ang", s->world[i].ang);
        pf_record_end(f);
    }
    fclose(f);
    return 0;
}

static int skel_import(const char* file, Skeleton* s) {
    FILE* f = fopen(file, "r");
    if (!f) return -1;

    int n = 0;
    if (pf_read_header(f, "skeleton", &n) != 0 || n > MAX_BONES) {
        fclose(f);
        return -1;
    }

    skel_init(s);
    PfRecord r;
    for (int k = 0; k < n; k++) {
        if (pf_read_record(f, &r) != 0 || strcmp(r.name, "bone") != 0) {
            fclose(f);
            return -1;
        }
        int idx, parent;
        float lx, ly, lang, wx, wy, ang;
        if (pf_get_int(&r, "index", &idx) || pf_get_int(&r, "parent", &parent) ||
            pf_get_float(&r, "lx", &lx) || pf_get_float(&r, "ly", &ly) ||
            pf_get_float(&r, "lang", &lang) || pf_get_float(&r, "wx", &wx) ||
            pf_get_float(&r, "wy", &wy) || pf_get_float(&r, "ang", &ang) ||
            idx != k) {
            fclose(f);
            return -1;
        }
        s->parent[k] = parent;
        s->local[k].x = lx; s->local[k].y = ly; s->local[k].ang = lang;
        s->world[k].x = wx; s->world[k].y = wy; s->world[k].ang = ang;
        s->count = k + 1;
    }
    fclose(f);
    return 0;
}

/* ---- legacy skeleton format (pre-unification): headerless CSV ----------
 *
 *   index,parent,lx,ly,lang,wx,wy,ang
 *
 * No magic/version/kind/count, fixed column order, no field names.  Kept so
 * files written before the shared format still load.
 */

static int skel_export_legacy(const Skeleton* s, const char* file) {
    FILE* f = fopen(file, "w");
    if (!f) return -1;
    for (int i = 0; i < s->count; i++) {
        fprintf(f, "%d,%d,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g\n",
                i, s->parent[i],
                (double)s->local[i].x, (double)s->local[i].y, (double)s->local[i].ang,
                (double)s->world[i].x, (double)s->world[i].y, (double)s->world[i].ang);
    }
    fclose(f);
    return 0;
}

static int skel_import_legacy(const char* file, Skeleton* s) {
    FILE* f = fopen(file, "r");
    if (!f) return -1;

    skel_init(s);
    char line[256];
    int n = 0, ok = 1;
    while (ok && fgets(line, sizeof line, f)) {
        char* p = line;
        while (*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n') p++;
        if (*p == '\0' || *p == '#') continue;

        int idx, parent;
        float lx, ly, lang, wx, wy, ang;
        if (sscanf(p, "%d,%d,%f,%f,%f,%f,%f,%f",
                   &idx, &parent, &lx, &ly, &lang, &wx, &wy, &ang) != 8 ||
            idx != n || n >= MAX_BONES) {
            ok = 0;
            break;
        }
        s->parent[n] = parent;
        s->local[n].x = lx; s->local[n].y = ly; s->local[n].ang = lang;
        s->world[n].x = wx; s->world[n].y = wy; s->world[n].ang = ang;
        s->count = ++n;
    }
    fclose(f);
    return (ok && n > 0) ? 0 : -1;
}

static int skel_same(const Skeleton* a, const Skeleton* b) {
    if (a->count != b->count) return 0;
    for (int i = 0; i < a->count; i++) {
        if (a->parent[i] != b->parent[i]) return 0;
        if (a->local[i].x != b->local[i].x || a->local[i].y != b->local[i].y ||
            a->local[i].ang != b->local[i].ang) return 0;
        if (a->world[i].x != b->world[i].x || a->world[i].y != b->world[i].y ||
            a->world[i].ang != b->world[i].ang) return 0;
    }
    return 1;
}

/* ------------------------------------------------------------ CLI modes */

static int skeleton_cli(int argc, char** argv) {
    const char* mode = argv[1];
    Skeleton src, back;

    if (strcmp(mode, "--pose-export") == 0 && argc == 3) {
        skel_pose_build(&src);
        return skel_export(&src, argv[2]) == 0 ? 0 : 1;
    }
    if (strcmp(mode, "--pose-import") == 0 && argc == 3) {
        if (skel_import(argv[2], &back) != 0) { printf("imported=0\n"); return 1; }
        printf("imported=1\ncount=%d\n", back.count);
        return 0;
    }
    if (strcmp(mode, "--pose-roundtrip") == 0 && argc == 3) {
        skel_pose_build(&src);
        int ok = skel_export(&src, argv[2]) == 0 &&
                 skel_import(argv[2], &back) == 0 &&
                 skel_same(&src, &back);
        printf("roundtrip=%d\n", ok ? 1 : 0);
        return ok ? 0 : 1;
    }
    if (strcmp(mode, "--legacy-export") == 0 && argc == 3) {
        skel_pose_build(&src);
        return skel_export_legacy(&src, argv[2]) == 0 ? 0 : 1;
    }
    if (strcmp(mode, "--legacy-check") == 0 && argc == 3) {
        skel_pose_build(&src);
        int ok = skel_import_legacy(argv[2], &back) == 0 && skel_same(&src, &back);
        printf("legacy=%d\n", ok ? 1 : 0);
        return ok ? 0 : 1;
    }

    fprintf(stderr, "usage: %s [--pose-export|--pose-import|--pose-roundtrip|"
                    "--legacy-export|--legacy-check] <file>\n", argv[0]);
    return 2;
}

int main(int argc, char** argv) {
    if (argc > 1) return skeleton_cli(argc, argv);

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
