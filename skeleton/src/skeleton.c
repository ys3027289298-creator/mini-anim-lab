
/* skeleton: bone hierarchy, world transforms and skinning matrix assembly.
 *
 * Two call paths, selected at build time by MAL_UNIFIED (default 1):
 *   1 -> every solve goes through the unified mal_solve() entry point
 *   0 -> the original static routines are called directly (legacy path)
 * The math below is shared by both paths, so their poses are identical.
 */
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "solver.h"
#include "skeleton.h"
#include "mal_err.h"

#ifndef MAL_UNIFIED
#define MAL_UNIFIED 1
#endif

#ifdef MAL_NO_MAIN
#define MAL_DISPATCH_FN
#else
#define MAL_DISPATCH_FN static
#endif

/* ---- original math, unchanged ---------------------------------------- */

static void skel_init(mal_skeleton* s) { memset(s, 0, sizeof(*s)); }

static int skel_add(mal_skeleton* s, int parent, float x, float y, float ang) {
    if (s->count >= MAL_SKEL_MAX_BONES) return -1;
    if (parent >= s->count) return -1;
    int i = s->count++;
    s->parent[i] = parent;
    s->local[i].x = x; s->local[i].y = y; s->local_ang[i] = ang;
    return i;
}

static void rotate(float x, float y, float a, float* ox, float* oy) {
    float c = cosf(a), s = sinf(a);
    *ox = x * c - y * s;
    *oy = x * s + y * c;
}

static void skel_update(mal_skeleton* s) {
    for (int i = 0; i < s->count; i++) {
        if (s->parent[i] < 0) {
            s->world[i].x = s->local[i].x;
            s->world[i].y = s->local[i].y;
            s->world_ang[i] = s->local_ang[i];
        } else {
            int p = s->parent[i];
            float rx, ry;
            rotate(s->local[i].x, s->local[i].y, s->world_ang[p], &rx, &ry);
            s->world[i].x = s->world[p].x + rx;
            s->world[i].y = s->world[p].y + ry;
            s->world_ang[i] = s->world_ang[p] + s->local_ang[i];
        }
    }
}

/* bone start point and end point (length 1 along local +x) */
#ifndef MAL_NO_MAIN
static void bone_segment(const mal_skeleton* s, int i,
                         float* x0, float* y0, float* x1, float* y1) {
    *x0 = s->world[i].x; *y0 = s->world[i].y;
    float rx, ry; rotate(1.0f, 0.0f, s->world_ang[i], &rx, &ry);
    *x1 = s->world[i].x + rx; *y1 = s->world[i].y + ry;
}
#endif

/* ---- unified-entry wrappers ------------------------------------------ */

void mal_skeleton_init(mal_skeleton* s) { skel_init(s); }
int mal_skeleton_capacity(void) { return MAL_SKEL_MAX_BONES; }

mal_status mal_skeleton_add(mal_skeleton* s, const mal_skel_add_req* req,
                            mal_skel_add_resp* resp, char* why, int whycap) {
    if (!s) return mal_fail(why, whycap, MAL_ERR_PARAM,
                            "skeleton: null instance");
    if (!req) return mal_fail(why, whycap, MAL_ERR_PARAM,
                              "skeleton: null request payload");
    if (!mal_finite_f(req->x) || !mal_finite_f(req->y) ||
        !mal_finite_f(req->ang))
        return mal_fail(why, whycap, MAL_ERR_PARAM,
                        "skeleton: non-finite local transform");
    if (req->parent < -1 || req->parent >= s->count)
        return mal_fail(why, whycap, MAL_ERR_PARAM,
                        "skeleton: parent id %d out of range [0,%d)",
                        req->parent, s->count);
    if (s->count >= MAL_SKEL_MAX_BONES)
        return mal_fail(why, whycap, MAL_ERR_STATE,
                        "skeleton: capacity %d reached", MAL_SKEL_MAX_BONES);
    int idx = skel_add(s, req->parent, req->x, req->y, req->ang);
    if (resp) resp->index = idx;
    return MAL_OK;
}

mal_status mal_skeleton_update(mal_skeleton* s, char* why, int whycap) {
    if (!s) return mal_fail(why, whycap, MAL_ERR_PARAM,
                            "skeleton: null instance");
    if (s->count <= 0)
        return mal_fail(why, whycap, MAL_ERR_STATE,
                        "skeleton: nothing to update, skeleton is empty");
    skel_update(s);
    return MAL_OK;
}

mal_status mal_skeleton_snapshot(const mal_skeleton* s, mal_skel_pose* out,
                                 char* why, int whycap) {
    if (!s) return mal_fail(why, whycap, MAL_ERR_PARAM,
                            "skeleton: null instance");
    if (!out) return mal_fail(why, whycap, MAL_ERR_PARAM,
                              "skeleton: null snapshot output");
    out->count = s->count;
    for (int i = 0; i < s->count; i++) {
        out->world[i] = s->world[i];
        out->ang[i] = s->world_ang[i];
    }
    return MAL_OK;
}

/* dispatch entry used by the unified mal_solve(). Non-static in a
 * consolidated build (MAL_NO_MAIN) so pose/ can register it. */
MAL_DISPATCH_FN mal_status skeleton_dispatch(int op, const void* req,
                                             void* resp,
                                             char* why, int whycap) {
    if (op == MAL_OP_SKELETON_UPDATE) {
        const mal_skel_update_req* q = (const mal_skel_update_req*)req;
        if (!q->s) return mal_fail(why, whycap, MAL_ERR_PARAM,
                                   "skeleton: null instance");
        return mal_skeleton_update(q->s, why, whycap);
    }
    if (op == MAL_OP_SKELETON_ADD) {
        const mal_skel_add_wrap* q = (const mal_skel_add_wrap*)req;
        if (!q->s) return mal_fail(why, whycap, MAL_ERR_PARAM,
                                   "skeleton: null instance");
        return mal_skeleton_add(q->s, &q->in,
                                (mal_skel_add_resp*)resp, why, whycap);
    }
    return mal_fail(why, whycap, MAL_ERR_UNSUPPORTED,
                    "skeleton: no op %d", op);
}

#ifndef MAL_NO_MAIN
const mal_dispatch_entry mal_dispatch_table[] = {
    { MAL_MODULE_SKELETON, skeleton_dispatch },
    { (mal_module)0, 0 }
};
#endif

/* ---- self test ------------------------------------------------------- */

#ifndef MAL_NO_MAIN
static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

/* unified routing helpers; the legacy path calls the statics directly */
#if MAL_UNIFIED
static int u_skel_add(mal_skeleton* s, int parent, float x, float y, float a) {
    mal_skel_add_wrap w; w.s = s; w.in.parent = parent;
    w.in.x = x; w.in.y = y; w.in.ang = a;
    mal_skel_add_resp out;
    mal_request r = { MAL_MODULE_SKELETON, MAL_OP_SKELETON_ADD, &w, &out };
    if (mal_solve(&r, 0, 0) != MAL_OK) return -1;
    return out.index;
}
static void u_skel_update(mal_skeleton* s) {
    mal_skel_update_req q = { s };
    mal_request r = { MAL_MODULE_SKELETON, MAL_OP_SKELETON_UPDATE, &q, 0 };
    (void)mal_solve(&r, 0, 0);
}
#define SKEL_ADD(s,p,x,y,a) u_skel_add((s),(p),(x),(y),(a))
#define SKEL_UPDATE(s)       u_skel_update(s)
#else
#define SKEL_ADD(s,p,x,y,a) skel_add((s),(p),(x),(y),(a))
#define SKEL_UPDATE(s)       skel_update(s)
#endif

int main(void) {
    mal_skeleton s;
    skel_init(&s);
    int root = SKEL_ADD(&s, -1, 0.0f, 0.0f, 0.0f);
    int b1 = SKEL_ADD(&s, root, 1.0f, 0.0f, 0.0f);
    int b2 = SKEL_ADD(&s, b1, 1.0f, 0.0f, 0.0f);
    check("add_root", root == 0);
    check("add_child", b1 == 1 && b2 == 2);

    SKEL_UPDATE(&s);
    check("world_chain", fabsf(s.world[2].x - 2.0f) < 1e-5f && fabsf(s.world[2].y) < 1e-5f);

    s.local_ang[root] = 1.5707963f;   /* 90 deg */
    SKEL_UPDATE(&s);
    check("rotate_chain", fabsf(s.world[2].y - 2.0f) < 1e-4f && fabsf(s.world[2].x) < 1e-4f);

    float x0, y0, x1, y1;
    bone_segment(&s, b1, &x0, &y0, &x1, &y1);
    check("segment_dir", fabsf(x1 - x0) < 1e-4f && fabsf(y1 - y0 - 1.0f) < 1e-4f);

    check("reject_bad_parent", SKEL_ADD(&s, 9, 0, 0, 0) == -1);

    mal_skeleton full;
    skel_init(&full);
    int last = -1, rc = 0;
    for (int i = 0; i <= MAL_SKEL_MAX_BONES; i++) { rc = SKEL_ADD(&full, last, 1, 0, 0); if (rc >= 0) last = rc; }
    check("capacity", rc == -1 && full.count == MAL_SKEL_MAX_BONES);

#if MAL_UNIFIED
    /* unified entry must name the module on illegal parameters */
    char why[128];
    mal_request bad = { MAL_MODULE_SKELETON, MAL_OP_SKELETON_ADD, 0, 0 };
    check("unified_null_payload",
          mal_solve(&bad, why, sizeof why) == MAL_ERR_PARAM &&
          strstr(why, "skeleton:") == why);

    mal_skel_add_wrap w; w.s = &s;
    w.in.parent = 99; w.in.x = 0; w.in.y = 0; w.in.ang = 0;
    bad.req = &w; bad.resp = 0;
    check("unified_bad_parent",
          mal_solve(&bad, why, sizeof why) == MAL_ERR_PARAM &&
          strstr(why, "skeleton:") == why);

    mal_request bogus = { (mal_module)99, 0, &w, 0 };
    check("unified_bad_module",
          mal_solve(&bogus, why, sizeof why) == MAL_ERR_PARAM &&
          strstr(why, "solver:") == why);

    mal_request unlinked = { MAL_MODULE_IK, 0, &w, 0 };
    check("unified_unlinked_module",
          mal_solve(&unlinked, why, sizeof why) == MAL_ERR_UNSUPPORTED &&
          strstr(why, "solver:") == why);
#endif

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
#endif
