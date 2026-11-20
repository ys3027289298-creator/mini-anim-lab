
/* skeleton: bone hierarchy, world transforms and capsule-swept bone bounds.
 *
 * Bounds pipeline:
 *   1. every bone is treated as a capsule (unit segment swept by its radius)
 *   2. each capsule is sampled as a circumscribed polygon (contains the cap)
 *   3. per-bone hulls are merged bottom-up through the hierarchy
 *   4. the minimum-area oriented box is searched over hull edge directions
 *   5. the result is checked against a tight per-vertex AABB and a dense set
 *      of true capsule surface samples: missed samples must be exactly zero,
 *      area redundancy is bounded by the circumscribed-polygon overshoot.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#define MAX_BONES 16
#define CAPSULE_SEGMENTS 32                    /* vertices per circumscribed cap */
#define CAPSULE_POINTS (2 * CAPSULE_SEGMENTS)  /* both cap polygons */
#define MAX_HULL_POINTS (MAX_BONES * CAPSULE_POINTS)
/* a merge holds two hulls before re-hulling (parent accumulation + child) */
#define MERGE_CAP (2 * MAX_HULL_POINTS)
#define SKIN_SAMPLES 128                       /* dense surface samples per cap */
#define MAX_SKIN_POINTS (MAX_BONES * (2 * SKIN_SAMPLES + 4))
#define PI_F 3.14159265358979323846f

/* A circumscribed N-gon covers its disc with radius scaled by 1/cos(pi/N);
   a box built on that polygon can exceed the tight skin AABB by at most this
   factor squared. For N=32 the bound is ~0.97%. */
#define POLY_OVERSHOOT (1.0f / cosf(PI_F / (float)CAPSULE_SEGMENTS))
#define BOX_REDUNDANCY_BOUND (POLY_OVERSHOOT * POLY_OVERSHOOT - 1.0f)

/* validation error codes; any rejected call still returns -1 */
enum {
    SKEL_OK = 0,
    SKEL_ERR_CAPACITY = -1,
    SKEL_ERR_CYCLE = -2,
    SKEL_ERR_RADIUS = -3,
    SKEL_ERR_PARENT = -4,
    SKEL_ERR_INDEX = -5
};

typedef struct { float x, y; float ang; } BoneLocal;
typedef struct { float x, y, ang; } BoneWorld;
typedef struct { float x, y; } V2;
typedef struct { V2 a, b; float r; } Capsule;

typedef struct {
    int parent[MAX_BONES];
    float radius[MAX_BONES];
    BoneLocal local[MAX_BONES];
    BoneWorld world[MAX_BONES];
    int count;
} Skeleton;

static void skel_init(Skeleton* s) { memset(s, 0, sizeof(*s)); }

/* last rejection reason, cleared by the next successful call */
static int g_skel_error = SKEL_OK;
static int skel_error(void) { return g_skel_error; }
static void skel_clear_error(void) { g_skel_error = SKEL_OK; }
static int skel_fail(int code) { g_skel_error = code; return -1; }

static int skel_add(Skeleton* s, int parent, float x, float y, float ang) {
    if (s->count >= MAX_BONES) return skel_fail(SKEL_ERR_CAPACITY);
    if (parent < -1 || parent >= s->count) return skel_fail(SKEL_ERR_PARENT);
    int i = s->count++;
    s->parent[i] = parent;
    s->local[i].x = x; s->local[i].y = y; s->local[i].ang = ang;
    s->radius[i] = 0.0f;
    skel_clear_error();
    return i;
}

/* add a bone carrying an explicit capsule radius; negative/non-finite is an error */
static int skel_add_capsule(Skeleton* s, int parent, float x, float y, float ang, float radius) {
    if (!(radius >= 0.0f) || !isfinite(radius)) return skel_fail(SKEL_ERR_RADIUS);
    int i = skel_add(s, parent, x, y, ang);
    if (i < 0) return i;
    s->radius[i] = radius;
    return i;
}

static int skel_set_radius(Skeleton* s, int i, float radius) {
    if (i < 0 || i >= s->count) return skel_fail(SKEL_ERR_INDEX);
    if (!(radius >= 0.0f) || !isfinite(radius)) return skel_fail(SKEL_ERR_RADIUS);
    s->radius[i] = radius;
    skel_clear_error();
    return 0;
}

/* Reparent with a cycle check. Parents must precede their children so that a
   single in-order world-transform pass is valid. */
static int skel_set_parent(Skeleton* s, int i, int parent) {
    if (i < 0 || i >= s->count) return skel_fail(SKEL_ERR_INDEX);
    if (parent < -1 || parent >= s->count) return skel_fail(SKEL_ERR_PARENT);
    for (int p = parent, step = 0; p >= 0; p = s->parent[p], step++) {
        if (p == i || step > s->count) return skel_fail(SKEL_ERR_CYCLE);
    }
    if (parent >= i) return skel_fail(SKEL_ERR_PARENT);
    s->parent[i] = parent;
    skel_clear_error();
    return 0;
}

/* Full structural check; also catches corrupted/externally built hierarchies.
   Cycles are reported before the parent-index ordering rule. */
static int skel_validate(const Skeleton* s) {
    if (s->count < 0 || s->count > MAX_BONES) return SKEL_ERR_CAPACITY;
    for (int i = 0; i < s->count; i++) {
        if (!(s->radius[i] >= 0.0f) || !isfinite(s->radius[i])) return SKEL_ERR_RADIUS;
        int p = s->parent[i];
        if (p < -1 || p >= s->count) return SKEL_ERR_PARENT;
        for (int step = 0; p >= 0; step++) {
            if (p == i || step > s->count) return SKEL_ERR_CYCLE;
            p = s->parent[p];
            if (p < -1 || p >= s->count) return SKEL_ERR_PARENT;
        }
        if (s->parent[i] >= i) return SKEL_ERR_PARENT;
    }
    return SKEL_OK;
}

static void rotate(float x, float y, float a, float* ox, float* oy) {
    float c = cosf(a), s = sinf(a);
    *ox = x * c - y * s;
    *oy = x * s + y * c;
}

/* returns a validation error code, or SKEL_OK */
static int skel_update(Skeleton* s) {
    int rc = skel_validate(s);
    if (rc != SKEL_OK) return rc;
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
    return SKEL_OK;
}

/* bone start point and end point (length 1 along local +x) */
static void bone_segment(const Skeleton* s, int i, float* x0, float* y0, float* x1, float* y1) {
    *x0 = s->world[i].x; *y0 = s->world[i].y;
    float rx, ry; rotate(1.0f, 0.0f, s->world[i].ang, &rx, &ry);
    *x1 = s->world[i].x + rx; *y1 = s->world[i].y + ry;
}

static Capsule bone_capsule(const Skeleton* s, int i) {
    Capsule c;
    bone_segment(s, i, &c.a.x, &c.a.y, &c.b.x, &c.b.y);
    c.r = s->radius[i];
    return c;
}

/* ---------- 2D geometry helpers ---------- */

static float cross2(V2 a, V2 b) { return a.x * b.y - a.y * b.x; }

static float dot2(V2 a, V2 b) { return a.x * b.x + a.y * b.y; }

static V2 v2sub(V2 a, V2 b) { V2 r = { a.x - b.x, a.y - b.y }; return r; }

static int cmp_point_x(const void* pa, const void* pb) {
    const V2* a = (const V2*)pa;
    const V2* b = (const V2*)pb;
    if (a->x < b->x) return -1;
    if (a->x > b->x) return 1;
    if (a->y < b->y) return -1;
    if (a->y > b->y) return 1;
    return 0;
}

/* Andrew monotone chain on a caller-owned copy; writes the CCW hull back into
   pts and returns its size. Hull vertices from the sorted list are collected
   into scratch because the upper hull still reads the original order; the
   chain can hold up to n + (n-1) points before popping, hence 2n scratch. */
static int convex_hull(V2* pts, int n) {
    if (n < 3) return n;
    qsort(pts, (size_t)n, sizeof(V2), cmp_point_x);
    V2 out[MERGE_CAP];
    int k = 0;
    for (int i = 0; i < n; i++) {
        while (k >= 2 && cross2(v2sub(out[k - 1], out[k - 2]), v2sub(pts[i], out[k - 2])) <= 0.0f)
            k--;
        out[k++] = pts[i];
    }
    int lower = k + 1;
    for (int i = n - 2; i >= 0; i--) {
        while (k >= lower && cross2(v2sub(out[k - 1], out[k - 2]), v2sub(pts[i], out[k - 2])) <= 0.0f)
            k--;
        out[k++] = pts[i];
    }
    k -= 1; /* drop the repeated first point */
    memcpy(pts, out, sizeof(V2) * (size_t)k);
    return k;
}

/* smallest-area oriented box over a CCW hull; one side is collinear with a
   hull edge (rotating-calipers property), so edges give all candidate axes */
static int hull_min_box(const V2* hull, int n, V2 out[4]) {
    if (n <= 0) return 0;
    if (n == 1) {
        for (int i = 0; i < 4; i++) out[i] = hull[0];
        return 1;
    }
    int found = 0;
    float best_area = FLT_MAX;
    for (int i = 0; i < n; i++) {
        V2 e = v2sub(hull[(i + 1) % n], hull[i]);
        float el = sqrtf(e.x * e.x + e.y * e.y);
        if (el < 1e-9f) continue;
        V2 ax = { e.x / el, e.y / el };
        V2 ay = { -ax.y, ax.x };
        float minu = FLT_MAX, maxu = -FLT_MAX, minv = FLT_MAX, maxv = -FLT_MAX;
        for (int j = 0; j < n; j++) {
            float u = dot2(hull[j], ax), v = dot2(hull[j], ay);
            if (u < minu) minu = u;
            if (u > maxu) maxu = u;
            if (v < minv) minv = v;
            if (v > maxv) maxv = v;
        }
        float area = (maxu - minu) * (maxv - minv);
        if (area < best_area) {
            best_area = area;
            float cu = (minu + maxu) * 0.5f, cv = (minv + maxv) * 0.5f;
            V2 center = { ax.x * cu + ay.x * cv, ax.y * cu + ay.y * cv };
            float hu = (maxu - minu) * 0.5f, hv = (maxv - minv) * 0.5f;
            out[0] = (V2){ center.x - ax.x * hu - ay.x * hv, center.y - ax.y * hu - ay.y * hv };
            out[1] = (V2){ center.x + ax.x * hu - ay.x * hv, center.y + ax.y * hu - ay.y * hv };
            out[2] = (V2){ center.x + ax.x * hu + ay.x * hv, center.y + ax.y * hu + ay.y * hv };
            out[3] = (V2){ center.x - ax.x * hu + ay.x * hv, center.y - ax.y * hu + ay.y * hv };
            found = 1;
        }
    }
    return found;
}

/* axis-aligned box of a point set */
static void aabb_box(const V2* pts, int n, V2 out[4], float* area) {
    float minx = FLT_MAX, maxx = -FLT_MAX, miny = FLT_MAX, maxy = -FLT_MAX;
    for (int i = 0; i < n; i++) {
        if (pts[i].x < minx) minx = pts[i].x;
        if (pts[i].x > maxx) maxx = pts[i].x;
        if (pts[i].y < miny) miny = pts[i].y;
        if (pts[i].y > maxy) maxy = pts[i].y;
    }
    if (n <= 0) { minx = maxx = miny = maxy = 0.0f; }
    out[0] = (V2){ minx, miny }; out[1] = (V2){ maxx, miny };
    out[2] = (V2){ maxx, maxy }; out[3] = (V2){ minx, maxy };
    *area = (maxx - minx) * (maxy - miny);
}

static float quad_area(const V2 q[4]) {
    float a = 0.0f;
    for (int i = 0; i < 4; i++) a += cross2(q[i], q[(i + 1) % 4]);
    return fabsf(a) * 0.5f;
}

/* true-distance excess of p outside a box; zero means p is contained */
static float point_box_excess(V2 p, const V2 q[4]) {
    V2 u = v2sub(q[1], q[0]);
    V2 v = v2sub(q[3], q[0]);
    float lu = dot2(u, u), lv = dot2(v, v);
    float du = lu > 1e-18f ? dot2(v2sub(p, q[0]), u) / lu : 0.0f;
    float dv = lv > 1e-18f ? dot2(v2sub(p, q[0]), v) / lv : 0.0f;
    float eu = fmaxf(0.0f, fmaxf(-du, du - 1.0f));
    float ev = fmaxf(0.0f, fmaxf(-dv, dv - 1.0f));
    if (eu == 0.0f && ev == 0.0f) return 0.0f;
    /* eu/ev are normalized; scale back by the box half-extents */
    float ex = eu * 0.5f * sqrtf(lu);
    float ey = ev * 0.5f * sqrtf(lv);
    return sqrtf(ex * ex + ey * ey);
}

/* ---------- capsule sweep volumes ---------- */

/* Capsule = convex hull of its two end discs, so two circumscribed N-gons at
   the endpoints enclose the true capsule, rounded caps included. Vertex angles
   are aligned to the bone normal so the polygon sides are tangent to the body.
   Returns CAPSULE_POINTS for r>0, two endpoints for r==0. */
static int capsule_points(Capsule c, V2* out) {
    if (c.r <= 0.0f) {
        out[0] = c.a;
        out[1] = c.b;
        return 2;
    }
    V2 d = v2sub(c.b, c.a);
    float len = sqrtf(d.x * d.x + d.y * d.y);
    float base = (len < 1e-6f) ? 0.0f : atan2f(d.y, d.x) + 0.5f * PI_F;
    float rad = c.r * POLY_OVERSHOOT;
    for (int j = 0; j < CAPSULE_SEGMENTS; j++) {
        float t = base + 2.0f * PI_F * (float)j / (float)CAPSULE_SEGMENTS;
        out[j].x = c.a.x + rad * cosf(t);
        out[j].y = c.a.y + rad * sinf(t);
        out[CAPSULE_SEGMENTS + j].x = c.b.x + rad * cosf(t);
        out[CAPSULE_SEGMENTS + j].y = c.b.y + rad * sinf(t);
    }
    return CAPSULE_POINTS;
}

/* dense samples of the true capsule surface, standing in for skinned vertices */
static int capsule_skin_points(Capsule c, V2* out) {
    V2 d = v2sub(c.b, c.a);
    float len = sqrtf(d.x * d.x + d.y * d.y);
    int k = 0;
    for (int i = 0; i < 2; i++) {
        V2 p = (i == 0) ? c.a : c.b;
        for (int j = 0; j < SKIN_SAMPLES; j++) {
            float t = 2.0f * PI_F * (float)j / (float)SKIN_SAMPLES;
            out[k].x = p.x + c.r * cosf(t);
            out[k].y = p.y + c.r * sinf(t);
            k++;
        }
    }
    if (len > 1e-6f && c.r > 0.0f) {
        V2 u = { d.x / len, d.y / len }, n = { -u.y, u.x };
        for (int i = 0; i < 2; i++) {
            V2 p = (i == 0) ? c.a : c.b;
            out[k].x = p.x + n.x * c.r; out[k].y = p.y + n.y * c.r; k++;
            out[k].x = p.x - n.x * c.r; out[k].y = p.y - n.y * c.r; k++;
        }
    }
    return k;
}

/* Merge per-bone hulls bottom-up (children have larger indices, so descending
   index order finishes a subtree before it reaches its parent). Each merge
   re-hulls the union, keeping every intermediate hull tight. `hull` must hold
   MERGE_CAP points; on return *hull_n is the global hull size. */
static int skel_hull_layered(const Skeleton* s, V2* hull, int* hull_n) {
    int rc = skel_validate(s);
    if (rc != SKEL_OK) return rc;
    *hull_n = 0;
    if (s->count == 0) return SKEL_OK;
    V2* acc = (V2*)malloc(sizeof(V2) * (size_t)MERGE_CAP * MAX_BONES);
    if (!acc) return SKEL_ERR_CAPACITY;
    int sz[MAX_BONES];
    memset(sz, 0, sizeof(sz));
    for (int i = s->count - 1; i >= 0; i--) {
        V2* h = acc + (size_t)i * MERGE_CAP;
        V2 pts[CAPSULE_POINTS];
        int n = capsule_points(bone_capsule(s, i), pts);
        memcpy(h + sz[i], pts, sizeof(V2) * (size_t)n);
        sz[i] = convex_hull(h, sz[i] + n);

        int p = s->parent[i];
        if (p >= 0) {
            V2* ph = acc + (size_t)p * MERGE_CAP;
            memcpy(ph + sz[p], h, sizeof(V2) * (size_t)sz[i]);
            sz[p] = convex_hull(ph, sz[p] + sz[i]);
        } else {
            memcpy(hull + *hull_n, h, sizeof(V2) * (size_t)sz[i]);
            *hull_n = convex_hull(hull, *hull_n + sz[i]);
        }
    }
    free(acc);
    return SKEL_OK;
}

typedef struct {
    V2 oriented[4];
    float oriented_area;
    V2 vertex_box[4];
    float vertex_area;
} BoneBounds;

/* per-bone bounds: optimal box plus the per-vertex AABB reference */
static int bone_bounds(const Skeleton* s, int i, BoneBounds* out) {
    if (skel_validate(s) != SKEL_OK || i < 0 || i >= s->count) return -1;
    V2 pts[CAPSULE_POINTS];
    int n = capsule_points(bone_capsule(s, i), pts);
    V2 hull[CAPSULE_POINTS];
    memcpy(hull, pts, sizeof(V2) * (size_t)n);
    int hn = convex_hull(hull, n);
    if (!hull_min_box(hull, hn, out->oriented)) return -1;
    out->oriented_area = quad_area(out->oriented);
    aabb_box(pts, n, out->vertex_box, &out->vertex_area);
    return 0;
}

/* ---------- tests ---------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

/* count skinned surface samples that fall outside box */static int coverage_misses(const Skeleton* s, const V2 box[4], int* total_out) {
    int misses = 0, total = 0;
    V2 skin[MAX_SKIN_POINTS];
    for (int i = 0; i < s->count; i++) {
        int n = capsule_skin_points(bone_capsule(s, i), skin);
        if (total + n > MAX_SKIN_POINTS) n = MAX_SKIN_POINTS - total;
        for (int j = 0; j < n; j++) {
            total++;
            if (point_box_excess(skin[j], box) > 1e-3f) misses++;
        }
    }
    if (total_out) *total_out = total;
    return misses;
}

/* count surface samples of one bone that fall outside box */
static int bone_coverage_misses(const Skeleton* s, int bone, const V2 box[4], int* total_out) {
    int misses = 0;
    V2 skin[2 * SKIN_SAMPLES + 4];
    int n = capsule_skin_points(bone_capsule(s, bone), skin);
    for (int j = 0; j < n; j++)
        if (point_box_excess(skin[j], box) > 1e-3f) misses++;
    if (total_out) *total_out = n;
    return misses;
}

/* smallest capsule radius in the skeleton */
static float min_radius(const Skeleton* s) {
    float m = FLT_MAX;
    for (int i = 0; i < s->count; i++) if (s->radius[i] < m) m = s->radius[i];
    return m;
}

/* tight per-vertex AABB over all dense skin samples */
static float skin_vertex_aabb(const Skeleton* s, V2 box[4]) {
    V2 skin[MAX_SKIN_POINTS];
    int total = 0;
    for (int i = 0; i < s->count; i++) {
        int n = capsule_skin_points(bone_capsule(s, i), skin + total);
        total += n;
    }
    float area = 0.0f;
    aabb_box(skin, total, box, &area);
    return area;
}

/* helper: bounds of a whole skeleton via layered hull + optimal box */
static int skeleton_box(const Skeleton* s, V2 box[4], float* area, int* hull_n) {
    V2 hull[MERGE_CAP];
    int hn = 0;
    int rc = skel_hull_layered(s, hull, &hn);
    if (rc != SKEL_OK) return rc;
    if (!hull_min_box(hull, hn, box)) return -1;
    if (area) *area = quad_area(box);
    if (hull_n) *hull_n = hn;
    return SKEL_OK;
}

int main(void) {
    Skeleton s;
    skel_init(&s);
    int root = skel_add(&s, -1, 0.0f, 0.0f, 0.0f);
    int b1 = skel_add(&s, root, 1.0f, 0.0f, 0.0f);
    int b2 = skel_add(&s, b1, 1.0f, 0.0f, 0.0f);
    check("add_root", root == 0);
    check("add_child", b1 == 1 && b2 == 2);

    check("update_ok", skel_update(&s) == SKEL_OK);
    check("world_chain", fabsf(s.world[2].x - 2.0f) < 1e-5f && fabsf(s.world[2].y) < 1e-5f);

    s.local[root].ang = 1.5707963f;   /* 90 deg */
    check("update_ok_rot", skel_update(&s) == SKEL_OK);
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

    /* --- capsule radii --- */
    Skeleton r;
    skel_init(&r);
    int r0 = skel_add_capsule(&r, -1, 0, 0, 0.0f, 0.25f);
    check("radius_add", r0 == 0 && fabsf(r.radius[0] - 0.25f) < 1e-6f);
    check("radius_reject_negative", skel_add_capsule(&r, r0, 1, 0, 0, -0.1f) == -1);
    check("radius_negative_code", skel_error() == SKEL_ERR_RADIUS);
    check("radius_reject_nan", skel_set_radius(&r, r0, NAN) == -1 && skel_error() == SKEL_ERR_RADIUS);
    check("radius_set_ok", skel_set_radius(&r, r0, 0.25f) == 0 && fabsf(r.radius[0] - 0.25f) < 1e-6f);
    check("radius_reject_negative_set", skel_set_radius(&r, r0, -1.0f) == -1 && skel_error() == SKEL_ERR_RADIUS);
    check("radius_reject_index", skel_set_radius(&r, 5, 1.0f) == -1 && skel_error() == SKEL_ERR_INDEX);

    /* --- hierarchy cycle rejection --- */
    Skeleton cyc;
    skel_init(&cyc);
    int c0 = skel_add(&cyc, -1, 0, 0, 0);
    int c1 = skel_add(&cyc, c0, 1, 0, 0);
    int c2 = skel_add(&cyc, c1, 1, 0, 0);
    check("cycle_self_parent", skel_set_parent(&cyc, c1, c1) == -1 && skel_error() == SKEL_ERR_CYCLE);
    check("cycle_direct", skel_set_parent(&cyc, c1, c2) == -1 && skel_error() == SKEL_ERR_CYCLE);
    check("cycle_indirect", skel_set_parent(&cyc, c0, c2) == -1 && skel_error() == SKEL_ERR_CYCLE);
    check("cycle_preserves", cyc.parent[0] == -1 && cyc.parent[1] == c0 && cyc.parent[2] == c1);
    check("reparent_ok", skel_set_parent(&cyc, c2, c0) == 0 && cyc.parent[2] == c0);

    Skeleton corrupt;
    skel_init(&corrupt);
    corrupt.count = 2;
    corrupt.parent[0] = 1;
    corrupt.parent[1] = 0;
    check("validate_cycle", skel_validate(&corrupt) == SKEL_ERR_CYCLE);
    check("update_rejects_cycle", skel_update(&corrupt) == SKEL_ERR_CYCLE);
    V2 bad_hull[MERGE_CAP];
    check("hull_rejects_cycle", skel_hull_layered(&corrupt, bad_hull, &rc) == SKEL_ERR_CYCLE);
    corrupt.parent[0] = -1; corrupt.parent[1] = 0;
    corrupt.radius[1] = -0.5f;
    check("validate_radius", skel_validate(&corrupt) == SKEL_ERR_RADIUS);
    check("update_rejects_radius", skel_update(&corrupt) == SKEL_ERR_RADIUS);

    /* acyclic but ordered after its parent: single-pass update would read a
       stale parent transform, so validation rejects the ordering */
    Skeleton order;
    skel_init(&order);
    order.count = 3;
    order.parent[0] = 1; order.parent[1] = 2; order.parent[2] = -1;
    check("validate_order", skel_validate(&order) == SKEL_ERR_PARENT);
    check("set_parent_order", skel_set_parent(&order, 0, 2) == -1 && skel_error() == SKEL_ERR_PARENT);

    /* --- bent arm: coverage, optimality, redundancy --- */
    Skeleton arm;
    skel_init(&arm);
    int shoulder = skel_add_capsule(&arm, -1, 0.0f, 0.0f, 0.30f, 0.35f);
    int elbow = skel_add_capsule(&arm, shoulder, 1.0f, 0.0f, -0.90f, 0.28f);
    int wrist = skel_add_capsule(&arm, elbow, 1.0f, 0.0f, 1.40f, 0.20f);
    check("arm_chain", shoulder == 0 && elbow == 1 && wrist == 2);
    check("arm_update", skel_update(&arm) == SKEL_OK);

    V2 box[4];
    float oarea = 0.0f;
    int hn = 0;
    check("arm_box", skeleton_box(&arm, box, &oarea, &hn) == SKEL_OK && hn > 2);

    /* every true capsule surface point must be inside: zero misses allowed */
    int skin_total = 0, misses = coverage_misses(&arm, box, &skin_total);
    printf("skin_points=%d\nmisses=%d\n", skin_total, misses);
    check("misses_zero", misses == 0);

    /* every bone's own caps must be inside too */
    float minr = min_radius(&arm);
    printf("min_radius=%.6f\n", minr);
    int bone_misses = 0;
    for (int i = 0; i < arm.count; i++) bone_misses += bone_coverage_misses(&arm, i, box, NULL);
    check("per_bone_caps_covered", bone_misses == 0);

    /* optimal box is no larger than the hull AABB */
    V2 hull[MERGE_CAP];
    int hn2 = 0;
    skel_hull_layered(&arm, hull, &hn2);
    V2 hbox[4];
    float harea;
    aabb_box(hull, hn2, hbox, &harea);
    printf("oriented_area=%.6f\naabb_area=%.6f\n", oarea, harea);
    check("box_beats_aabb", oarea <= harea + 1e-4f);

    /* layered merging must match a single hull over all bone samples */
    {
        V2 all[MAX_HULL_POINTS];
        int k = 0;
        for (int i = 0; i < arm.count; i++) {
            V2 pts[CAPSULE_POINTS];
            int n = capsule_points(bone_capsule(&arm, i), pts);
            memcpy(all + k, pts, sizeof(V2) * (size_t)n);
            k += n;
        }
        int an = convex_hull(all, k);
        V2 abox[4];
        check("direct_hull_box", hull_min_box(all, an, abox) == 1);
        check("layered_matches_direct", fabsf(quad_area(abox) - oarea) <= 1e-4f * (1.0f + oarea));
    }

    /* redundancy against the tight per-vertex AABB must respect the bound */
    V2 tight[4];
    float tarea = skin_vertex_aabb(&arm, tight);
    float redundancy = (oarea - tarea) / tarea;
    printf("tight_vertex_area=%.6f\nredundancy=%.6f\nbound=%.6f\n", tarea, redundancy, BOX_REDUNDANCY_BOUND);
    check("redundancy_bound", redundancy <= BOX_REDUNDANCY_BOUND);
    check("tight_misses_zero", coverage_misses(&arm, tight, NULL) == 0);

    /* per-bone bounds are consistent with the budget */
    for (int i = 0; i < arm.count; i++) {
        BoneBounds bb;
        check("bone_bounds_ok", bone_bounds(&arm, i, &bb) == 0);
        check("bone_box_covers", bone_coverage_misses(&arm, i, bb.oriented, NULL) == 0);
        check("bone_box_within_bound", (bb.oriented_area - bb.vertex_area) / bb.vertex_area <= BOX_REDUNDANCY_BOUND);
        check("bone_box_not_larger", bb.oriented_area <= bb.vertex_area + 1e-4f);
    }

    /* zero-radius bone is a bare segment, still covered */
    Skeleton dot;
    skel_init(&dot);
    int d0 = skel_add_capsule(&dot, -1, 2.0f, -1.0f, 0.0f, 0.0f);
    check("zero_radius_add", d0 == 0 && skel_update(&dot) == SKEL_OK);
    BoneBounds bb;
    check("zero_radius_box", bone_bounds(&dot, d0, &bb) == 0);
    check("zero_radius_area", bb.oriented_area < 1e-6f);
    check("zero_radius_misses_zero", bone_coverage_misses(&dot, d0, bb.oriented, NULL) == 0);

    /* rigid rotation of the whole chain keeps the box area equal */
    Skeleton arm_rot = arm;
    arm_rot.local[shoulder].ang += 0.7f;
    skel_update(&arm_rot);
    V2 rbox[4];
    float rarea = 0.0f;
    check("rot_box", skeleton_box(&arm_rot, rbox, &rarea, NULL) == SKEL_OK);
    check("rotation_invariant", fabsf(rarea - oarea) <= 1e-3f * (1.0f + oarea));
    check("rotation_misses_zero", coverage_misses(&arm_rot, rbox, NULL) == 0);

    /* a fat cap must widen the bounds: stadium of length 1, radius 0.9 has a
       provably optimal enclosing-box area of (1+2r)*(2r) = 5.04 */
    Skeleton cap;
    skel_init(&cap);
    int cap0 = skel_add_capsule(&cap, -1, 0.0f, 0.0f, 0.0f, 0.9f);
    skel_update(&cap);
    V2 cbox[4];
    float carea = 0.0f;
    check("fat_cap_box", cap0 == 0 && skeleton_box(&cap, cbox, &carea, NULL) == SKEL_OK);
    check("fat_cap_misses_zero", coverage_misses(&cap, cbox, NULL) == 0);
    printf("fat_cap_area=%.6f\n", carea);
    check("fat_cap_area_bound", carea >= 5.04f - 1e-3f && carea <= 5.04f * (1.0f + BOX_REDUNDANCY_BOUND + 1e-4f));

    /* full MAX_BONES hierarchy: layered merge stays covering and bounded */
    Skeleton many;
    skel_init(&many);
    int prev = -1;
    for (int i = 0; i < MAX_BONES; i++)
        prev = skel_add_capsule(&many, prev, 0.6f, 0.0f, 0.35f, 0.15f);
    skel_update(&many);
    V2 mbox[4];
    float marea = 0.0f;
    int mhn = 0;
    check("many_bones_box", skeleton_box(&many, mbox, &marea, &mhn) == SKEL_OK && mhn > 2);
    check("many_bones_misses_zero", coverage_misses(&many, mbox, NULL) == 0);
    {
        V2 mtight[4];
        float mtarea = skin_vertex_aabb(&many, mtight);
        float mred = (marea - mtarea) / mtarea;
        printf("many_redundancy=%.6f\n", mred);
        /* the optimal box may beat the axis-aligned skin AABB; it may never
           exceed it by more than the polygon overshoot bound */
        check("many_bones_redundancy", mred <= BOX_REDUNDANCY_BOUND);
    }
    {
        /* branching: two children off one root, layered merge sees both */
        Skeleton branch;
        skel_init(&branch);
        int b0 = skel_add_capsule(&branch, -1, 0.0f, 0.0f, 0.6f, 0.2f);
        int bl = skel_add_capsule(&branch, b0, 1.0f, 0.0f, 0.9f, 0.2f);
        int br = skel_add_capsule(&branch, b0, 1.0f, 0.0f, -0.9f, 0.3f);
        skel_update(&branch);
        check("branch_chain", bl == 1 && br == 2);
        V2 bbox[4];
        check("branch_box", skeleton_box(&branch, bbox, &carea, NULL) == SKEL_OK);
        check("branch_misses_zero", coverage_misses(&branch, bbox, NULL) == 0);
    }

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
