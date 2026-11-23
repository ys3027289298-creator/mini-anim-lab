
/* skeleton: bone hierarchy, world transforms and skinning matrix assembly. */
#include <stdio.h>
#include <string.h>
#include <math.h>

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

int main(void) {
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
