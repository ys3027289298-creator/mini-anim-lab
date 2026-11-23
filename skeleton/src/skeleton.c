
/* skeleton: bone hierarchy, world transforms and skinning matrix assembly. */
#include <stdio.h>
#include <string.h>
#include <math.h>

#ifndef MAL_MODULE
#define MAL_MODULE "skeleton"
#endif
#define MAL_TIME_ITERS 600000L
#include "../../common/harness.h"

#define MAX_BONES 16

typedef struct { float x, y; float ang; } BoneLocal;
typedef struct { float x, y; float ang; } BoneWorld;

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

/* three-bone chain folded 90 degrees at the root, the pose used by the
   rotate_chain check; shared by the evidence paths below. */
MAL_UNUSED static void skel_build_pose(Skeleton* s) {
    skel_init(s);
    int root = skel_add(s, -1, 0.0f, 0.0f, 0.0f);
    int b1 = skel_add(s, root, 1.0f, 0.0f, 0.0f);
    skel_add(s, b1, 1.0f, 0.0f, 0.0f);
    s->local[root].ang = 1.5707963f;
    skel_update(s);
}

static void run_checks(void) {
    Skeleton s;
    skel_init(&s);
    int root = skel_add(&s, -1, 0.0f, 0.0f, 0.0f);
    int b1 = skel_add(&s, root, 1.0f, 0.0f, 0.0f);
    int b2 = skel_add(&s, b1, 1.0f, 0.0f, 0.0f);
    mal_check("add_root", root == 0);
    mal_check("add_child", b1 == 1 && b2 == 2);

    skel_update(&s);
    mal_check("world_chain", fabsf(s.world[2].x - 2.0f) < 1e-5f && fabsf(s.world[2].y) < 1e-5f);

    s.local[root].ang = 1.5707963f;   /* 90 deg */
    skel_update(&s);
    mal_check("rotate_chain", fabsf(s.world[2].y - 2.0f) < 1e-4f && fabsf(s.world[2].x) < 1e-4f);

    float x0, y0, x1, y1;
    bone_segment(&s, b1, &x0, &y0, &x1, &y1);
    mal_check("segment_dir", fabsf(x1 - x0) < 1e-4f && fabsf(y1 - y0 - 1.0f) < 1e-4f);

    mal_check("reject_bad_parent", skel_add(&s, 9, 0, 0, 0) == -1);

    Skeleton full;
    skel_init(&full);
    int last = -1, rc = 0;
    for (int i = 0; i <= MAX_BONES; i++) { rc = skel_add(&full, last, 1, 0, 0); if (rc >= 0) last = rc; }
    mal_check("capacity", rc == -1 && full.count == MAX_BONES);
}

int main(void) {
    run_checks();

#ifdef MAL_EVIDENCE_POSE
    if (mal_failures == 0) {
        Skeleton s;
        skel_build_pose(&s);
        mal_pose("bones", s.count);
        for (int i = 0; i < s.count; i++) {
            char k[24];
            snprintf(k, sizeof k, "w%d_x", i);   mal_pose(k, s.world[i].x);
            snprintf(k, sizeof k, "w%d_y", i);   mal_pose(k, s.world[i].y);
            snprintf(k, sizeof k, "w%d_ang", i); mal_pose(k, s.world[i].ang);
        }
        float x0, y0, x1, y1;
        bone_segment(&s, 1, &x0, &y0, &x1, &y1);
        mal_pose("seg_x0", x0); mal_pose("seg_y0", y0);
        mal_pose("seg_x1", x1); mal_pose("seg_y1", y1);
        mal_emit_pose();
    }
#endif

#ifdef MAL_EVIDENCE_IMAGE
    if (mal_failures == 0) {
        Skeleton s;
        MalImage im;
        skel_build_pose(&s);
        mal_img_clear(&im);
        static const unsigned char cols[3][3] = {
            { 80, 200, 255 }, { 255, 180, 90 }, { 120, 235, 150 }
        };
        for (int i = 0; i < s.count; i++) {
            float x0, y0, x1, y1;
            bone_segment(&s, i, &x0, &y0, &x1, &y1);
            mal_img_line(&im, 64.0 + x0 * 40.0, 110.0 - y0 * 40.0,
                         64.0 + x1 * 40.0, 110.0 - y1 * 40.0,
                         cols[i % 3][0], cols[i % 3][1], cols[i % 3][2]);
        }
        mal_emit_png(&im, MAL_MODULE);
    }
#endif

#ifdef MAL_EVIDENCE_TIMING
    if (mal_failures == 0) {
        Skeleton s;
        skel_build_pose(&s);
        double acc = 0.0;
        double t0 = mal_now();
        for (long it = 0; it < MAL_TIME_ITERS; it++) {
            s.local[0].ang = 1.5707963f * (float)(it % 7) / 7.0f;
            skel_update(&s);
            for (int i = 0; i < s.count; i++) {
                float x0, y0, x1, y1;
                bone_segment(&s, i, &x0, &y0, &x1, &y1);
                acc += x1 + y1 + 0.5f * (x0 + y0);
            }
        }
        double t1 = mal_now();
        if (!isfinite(acc)) { mal_check("time_finite", 0); }
        mal_emit_time(MAL_TIME_TIER, t1 - t0, MAL_TIME_ITERS);
    }
#endif

    printf("all_ok=%d\n", mal_failures == 0);
    return mal_failures == 0 ? 0 : 1;
}
