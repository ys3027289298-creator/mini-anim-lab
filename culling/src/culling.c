
/* culling: hierarchical AABB culling against a 2D camera frustum.
 *
 * Pipeline:
 *   1. hc_validate      - hierarchy sanity, DFS cycle detection (loop -> error)
 *   2. hc_rebuild       - bottom-up merged AABB rebuild in DFS postorder
 *   3. hc_cull          - root-down frustum culling, per node
 *                          OUTSIDE (subtree pruned) / INSIDE (subtree passes)
 *                          / INTERSECT (descend)
 *   4. hc_reference     - per-bone brute-force coarse test, ground-truth set
 *   5. hc_compare       - missed/wrong cull counts with structural upper bounds
 *   6. HcCache          - previous-frame conclusions with coherence validation;
 *                          cache vs fresh-rebuild disagreement is reported
 *
 * Soundness invariant: while child AABBs stay contained in parent AABBs
 * (hc_check_containment), a node whose own box intersects the frustum can
 * never be hierarchy-culled, so wrong cull count is provably 0. Missed culls
 * (conservative over-cover from merged boxes) are bounded by the total size
 * of INTERSECT-classified subtrees.
 */
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>

#define HC_MAX_NODES  32
#define HC_MAX_PLANES 8
#define HC_EPS        1e-5f
#define HC_PI         3.14159265358979323846f

typedef enum {
    HC_OK = 0,
    HC_ERR_BAD_PARENT,   /* parent index out of range or self parent */
    HC_ERR_CYCLE,        /* loop in hierarchy */
    HC_ERR_CONTAINMENT,  /* rebuilt child box not inside parent box */
    HC_ERR_CACHE,        /* cached conclusion disagrees with fresh rebuild */
    HC_ERR_CAPACITY
} HcError;

static const char* hc_err_str(int e) {
    switch (e) {
        case HC_OK:               return "ok";
        case HC_ERR_BAD_PARENT:   return "bad_parent";
        case HC_ERR_CYCLE:        return "cycle";
        case HC_ERR_CONTAINMENT:  return "containment";
        case HC_ERR_CACHE:        return "cache_inconsistent";
        case HC_ERR_CAPACITY:     return "capacity";
        default:                  return "unknown";
    }
}

typedef struct { float x0, y0, x1, y1; } AABB2;          /* min corner, max corner */
typedef struct { float nx, ny, d; } Plane2;              /* inside: nx*x+ny*y+d >= 0 */
typedef struct { Plane2 pl[HC_MAX_PLANES]; int n; } Frustum2;

typedef struct {
    int parent[HC_MAX_NODES];
    AABB2 local[HC_MAX_NODES];      /* own geometry bound, e.g. one bone segment */
    AABB2 world[HC_MAX_NODES];      /* rebuilt merged bound */
    int order[HC_MAX_NODES];        /* DFS postorder: children before parents */
    int subtree_size[HC_MAX_NODES];
    int count;
} HcTree;

typedef struct {
    int valid;
    int fid;                        /* frustum identity; conclusions only
                                       comparable within the same view */
    AABB2 box[HC_MAX_NODES];
    int vis[HC_MAX_NODES];
} HcCache;

enum { HC_CULLED = 0, HC_VISIBLE = 1 };

typedef enum { BC_OUTSIDE = 0, BC_INTERSECT = 1, BC_INSIDE = 2 } BoundClass;

typedef struct {
    int ref_visible;        /* brute-force visible bones */
    int hier_visible;       /* hierarchy visible nodes */
    int missed_cull;        /* lou ti: reference culled, hierarchy kept visible */
    int wrong_cull;         /* wu ti: reference visible, hierarchy culled */
    int missed_bound;       /* structural upper bound on missed_cull */
    int wrong_bound;        /* upper bound on wrong_cull, provably 0 */
    int nodes_tested;       /* AABB/frustum classifications actually performed */
    int subtree_shortcuts;  /* nodes passed via an INSIDE ancestor without a test */
} HcStats;

/* ---------- basic geometry ---------- */

static AABB2 aabb_make(float x0, float y0, float x1, float y1) {
    AABB2 b = { x0, y0, x1, y1 };
    return b;
}

static AABB2 aabb_empty(void) {
    AABB2 b = { FLT_MAX, FLT_MAX, -FLT_MAX, -FLT_MAX };
    return b;
}

static void aabb_union(AABB2* a, const AABB2* b) {
    if (b->x0 < a->x0) a->x0 = b->x0;
    if (b->y0 < a->y0) a->y0 = b->y0;
    if (b->x1 > a->x1) a->x1 = b->x1;
    if (b->y1 > a->y1) a->y1 = b->y1;
}

static int aabb_contains(const AABB2* outer, const AABB2* inner) {
    return inner->x0 >= outer->x0 - HC_EPS &&
           inner->y0 >= outer->y0 - HC_EPS &&
           inner->x1 <= outer->x1 + HC_EPS &&
           inner->y1 <= outer->y1 + HC_EPS;
}

static int aabb_equal(const AABB2* a, const AABB2* b) {
    return fabsf(a->x0 - b->x0) <= HC_EPS && fabsf(a->y0 - b->y0) <= HC_EPS &&
           fabsf(a->x1 - b->x1) <= HC_EPS && fabsf(a->y1 - b->y1) <= HC_EPS;
}

static Plane2 plane_make(float nx, float ny, float px, float py) {
    Plane2 p = { nx, ny, -(nx * px + ny * py) };
    return p;
}

/* axis-aligned rectangular view region */
static void frustum_rect(Frustum2* f, float x0, float y0, float x1, float y1) {
    f->n = 4;
    f->pl[0] = plane_make(1.0f, 0.0f, x0, 0.0f);   /* x >= x0 */
    f->pl[1] = plane_make(-1.0f, 0.0f, x1, 0.0f);  /* x <= x1 */
    f->pl[2] = plane_make(0.0f, 1.0f, 0.0f, y0);   /* y >= y0 */
    f->pl[3] = plane_make(0.0f, -1.0f, 0.0f, y1);  /* y <= y1 */
}

/* genuine 2D camera frustum: near/far planes plus two FOV side rays */
static void frustum_camera(Frustum2* f, float ex, float ey,
                           float dir_ang, float half_fov,
                           float near, float far) {
    float dx = cosf(dir_ang), dy = sinf(dir_ang);
    float rr = dir_ang + half_fov, rl = dir_ang - half_fov;
    float rrx = cosf(rr), rry = sinf(rr);   /* right side ray */
    float rlx = cosf(rl), rly = sinf(rl);   /* left side ray */
    f->n = 4;
    f->pl[0] = plane_make(dx, dy, ex + dx * near, ey + dy * near);      /* near */
    f->pl[1] = plane_make(-dx, -dy, ex + dx * far, ey + dy * far);      /* far  */
    f->pl[2] = plane_make(rry, -rrx, ex, ey);   /* right side: left of ray  */
    f->pl[3] = plane_make(-rly, rlx, ex, ey);   /* left side: right of ray  */
}

/* AABB vs convex half-plane set: OUTSIDE / INTERSECT / fully INSIDE */
static BoundClass aabb_classify(const AABB2* b, const Frustum2* f) {
    float cx = 0.5f * (b->x0 + b->x1), cy = 0.5f * (b->y0 + b->y1);
    float hx = 0.5f * (b->x1 - b->x0), hy = 0.5f * (b->y1 - b->y0);
    int crossed = 0;
    for (int k = 0; k < f->n; k++) {
        float nx = f->pl[k].nx, ny = f->pl[k].ny;
        float dist = nx * cx + ny * cy + f->pl[k].d;
        float extent = fabsf(nx) * hx + fabsf(ny) * hy;
        if (dist < -extent - HC_EPS) return BC_OUTSIDE;
        if (dist < extent + HC_EPS) crossed = 1;
    }
    return crossed ? BC_INTERSECT : BC_INSIDE;
}

/* ---------- tree construction and validation ---------- */

static void hc_init(HcTree* t) { memset(t, 0, sizeof(*t)); }

static int hc_add(HcTree* t, int parent, AABB2 own) {
    if (t->count >= HC_MAX_NODES) return -HC_ERR_CAPACITY;
    if (parent != -1 && parent >= t->count) return -HC_ERR_BAD_PARENT;
    int i = t->count++;
    t->parent[i] = parent;
    t->local[i] = own;
    t->world[i] = aabb_empty();
    return i;
}

static int hc_dfs(HcTree* t, int i, int* color, int* pos) {
    if (color[i] == 1) return HC_ERR_CYCLE;
    if (color[i] == 2) return HC_OK;
    color[i] = 1;
    for (int j = 0; j < t->count; j++) {
        if (t->parent[j] == i) {
            int e = hc_dfs(t, j, color, pos);
            if (e != HC_OK) return e;
        }
    }
    color[i] = 2;
    t->order[(*pos)++] = i;   /* postorder: every child recorded first */
    return HC_OK;
}

/* rejects out-of-range parents, self loops and any cycle; fills postorder */
static int hc_validate(HcTree* t) {
    for (int i = 0; i < t->count; i++) {
        int p = t->parent[i];
        if (p != -1 && (p < 0 || p >= t->count || p == i))
            return HC_ERR_BAD_PARENT;
    }
    int color[HC_MAX_NODES] = { 0 };
    int pos = 0;
    for (int i = 0; i < t->count; i++) {
        if (t->parent[i] == -1) {
            int e = hc_dfs(t, i, color, &pos);
            if (e != HC_OK) return e;
        }
    }
    if (pos != t->count) return HC_ERR_CYCLE;  /* nodes only reachable via a loop */
    return HC_OK;
}

/* ---------- bottom-up rebuild ---------- */

static void hc_rebuild(HcTree* t) {
    for (int k = 0; k < t->count; k++) {
        int i = t->order[k];
        t->world[i] = t->local[i];
        t->subtree_size[i] = 1;
        for (int j = 0; j < t->count; j++) {
            if (t->parent[j] == i) {
                aabb_union(&t->world[i], &t->world[j]);
                t->subtree_size[i] += t->subtree_size[j];
            }
        }
    }
}

/* every rebuilt child box must stay inside its parent box */
static int hc_check_containment(const HcTree* t, int* bad_index) {
    for (int i = 0; i < t->count; i++) {
        int p = t->parent[i];
        if (p == -1) continue;
        if (!aabb_contains(&t->world[p], &t->world[i])) {
            if (bad_index) *bad_index = i;
            return HC_ERR_CONTAINMENT;
        }
    }
    return HC_OK;
}

/* ---------- hierarchical frustum cull ---------- */

static void hc_visit(const HcTree* t, int i, int inside_ancestor,
                     const Frustum2* f, int* vis, HcStats* st) {
    if (inside_ancestor) {
        vis[i] = HC_VISIBLE;
        st->subtree_shortcuts++;
        for (int j = 0; j < t->count; j++)
            if (t->parent[j] == i)
                hc_visit(t, j, 1, f, vis, st);
        return;
    }
    BoundClass cls = aabb_classify(&t->world[i], f);
    st->nodes_tested++;
    if (cls == BC_OUTSIDE) {
        /* subtree stays HC_CULLED: pruned here, no child tests at all */
        return;
    }
    vis[i] = HC_VISIBLE;
    if (cls == BC_INTERSECT) {
        /* everything below a partial box is the conservative region where a
         * missed cull can land; accumulated subtree sizes bound the count */
        st->missed_bound += t->subtree_size[i];
    }
    for (int j = 0; j < t->count; j++) {
        if (t->parent[j] == i)
            hc_visit(t, j, cls == BC_INSIDE, f, vis, st);
    }
}

static int hc_cull(const HcTree* t, const Frustum2* f, int* vis, HcStats* st) {
    memset(st, 0, sizeof(*st));
    for (int i = 0; i < t->count; i++) vis[i] = HC_CULLED;
    for (int i = 0; i < t->count; i++) {
        if (t->parent[i] == -1)
            hc_visit(t, i, 0, f, vis, st);
    }
    st->wrong_bound = 0;   /* containment makes wrong culls impossible */
    for (int i = 0; i < t->count; i++)
        if (vis[i] == HC_VISIBLE) st->hier_visible++;
    return HC_OK;
}

/* ---------- brute-force reference and set comparison ---------- */

static void hc_reference(const HcTree* t, const Frustum2* f, int* vis_ref) {
    for (int i = 0; i < t->count; i++)
        vis_ref[i] = aabb_classify(&t->local[i], f) == BC_OUTSIDE
                         ? HC_CULLED : HC_VISIBLE;
}

static void hc_compare(const HcTree* t, const int* vis_ref,
                       const int* vis_hier, HcStats* st) {
    for (int i = 0; i < t->count; i++) {
        if (vis_ref[i] == HC_VISIBLE) st->ref_visible++;
        if (vis_ref[i] == HC_CULLED && vis_hier[i] == HC_VISIBLE)
            st->missed_cull++;
        if (vis_ref[i] == HC_VISIBLE && vis_hier[i] == HC_CULLED)
            st->wrong_cull++;
    }
}

/* ---------- frame cache and coherence validation ---------- */

static void hc_cache_store(HcCache* c, int fid, const HcTree* t, const int* vis) {
    for (int i = 0; i < t->count; i++) {
        c->box[i] = t->world[i];
        c->vis[i] = vis[i];
    }
    c->fid = fid;
    c->valid = 1;
}

/* same frustum + unchanged cached box => cached conclusion must match the
 * fresh one. A view or box change invalidates the entry instead. */
static int hc_cache_coherence(const HcCache* c, int fid, const HcTree* t,
                              const int* fresh_vis, int* bad_index) {
    if (!c->valid || c->fid != fid) return HC_OK;
    for (int i = 0; i < t->count; i++) {
        if (aabb_equal(&c->box[i], &t->world[i]) && c->vis[i] != fresh_vis[i]) {
            if (bad_index) *bad_index = i;
            return HC_ERR_CACHE;
        }
    }
    return HC_OK;
}

/* one full frame: validate -> rebuild -> containment -> cull -> compare ->
 * cache coherence check -> cache refresh. A stale-cache warning is returned
 * while the fresh visibility set is still produced. */
static int hc_cull_frame(HcTree* t, int fid, const Frustum2* f, HcCache* cache,
                         int* vis, HcStats* st) {
    int e = hc_validate(t);
    if (e != HC_OK) return e;
    hc_rebuild(t);
    int bad;
    e = hc_check_containment(t, &bad);
    if (e != HC_OK) return e;
    hc_cull(t, f, vis, st);

    int ref[HC_MAX_NODES];
    hc_reference(t, f, ref);
    hc_compare(t, ref, vis, st);

    int warn = hc_cache_coherence(cache, fid, t, vis, &bad);
    hc_cache_store(cache, fid, t, vis);
    return warn;
}

/* ---------- self tests ---------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static HcTree build_demo_tree(void) {
    HcTree t;
    hc_init(&t);
    int r  = hc_add(&t, -1, aabb_make(-1, -1, 1, 1));     /* 0 root      */
    int a1 = hc_add(&t, r,  aabb_make(2, -1, 3, 1));      /* 1 chain A   */
    int a2 = hc_add(&t, a1, aabb_make(4, -1, 5, 1));      /* 2           */
    int b1 = hc_add(&t, r,  aabb_make(-1, 2, 1, 3));      /* 3 chain B   */
    int b2 = hc_add(&t, b1, aabb_make(-1, 4, 1, 5));      /* 4           */
    hc_add(&t, a2, aabb_make(6, -1, 7, 1));               /* 5 tip of A  */
    (void)b2;
    return t;
}

int main(void) {
    /* ---- hierarchy validation: bad parent and cycles ---- */
    HcTree t = build_demo_tree();
    check("validate_ok", hc_validate(&t) == HC_OK);

    HcTree bad = t;
    bad.parent[0] = 99;
    check("bad_parent_index", hc_validate(&bad) == HC_ERR_BAD_PARENT);

    bad = t;
    bad.parent[3] = 3;
    check("self_loop", hc_validate(&bad) == HC_ERR_BAD_PARENT);

    bad = t;
    bad.parent[1] = 2;   /* 1 <-> 2 loop */
    bad.parent[2] = 1;
    check("cycle_detected", hc_validate(&bad) == HC_ERR_CYCLE);
    check("err_str", strcmp(hc_err_str(HC_ERR_CYCLE), "cycle") == 0);
    Frustum2 guard;
    frustum_rect(&guard, 0, 0, 1, 1);
    HcCache dummy_cache;
    memset(&dummy_cache, 0, sizeof(dummy_cache));
    int guard_vis[HC_MAX_NODES];
    HcStats guard_st;
    check("cull_refuses_cycle",
          hc_cull_frame(&bad, 0, &guard, &dummy_cache, guard_vis, &guard_st)
          == HC_ERR_CYCLE);

    /* ---- bottom-up rebuild gives exact merged boxes ---- */
    hc_rebuild(&t);
    check("rebuild_leaf",
          aabb_equal(&t.world[5], &t.local[5]));
    AABB2 expect_a1 = aabb_make(2, -1, 7, 1);
    AABB2 expect_root = aabb_make(-1, -1, 7, 5);
    check("rebuild_merge_chain", aabb_equal(&t.world[1], &expect_a1));
    check("rebuild_merge_root", aabb_equal(&t.world[0], &expect_root));
    check("subtree_size", t.subtree_size[0] == 6 && t.subtree_size[1] == 3
                          && t.subtree_size[3] == 2);

    int bad_idx = -1;
    check("containment_ok", hc_check_containment(&t, &bad_idx) == HC_OK);
    AABB2 saved_w1 = t.world[1];
    t.world[1] = aabb_make(2.5f, -0.5f, 3, 0.5f);  /* simulate corrupt rebuild */
    check("containment_violation",
          hc_check_containment(&t, &bad_idx) == HC_ERR_CONTAINMENT &&
          bad_idx == 2);
    t.world[1] = saved_w1;

    /* ---- camera frustum plane orientation ---- */
    Frustum2 cam;
    frustum_camera(&cam, 0, 0, 0.0f, HC_PI / 4.0f, 1.0f, 10.0f);
    AABB2 pt_center = aabb_make(4.99f, -0.01f, 5.01f, 0.01f);
    AABB2 pt_side  = aabb_make(4.99f, 8.99f, 5.01f, 9.01f);
    AABB2 pt_behind = aabb_make(-1.01f, -0.01f, -0.99f, 0.01f);
    AABB2 pt_far = aabb_make(10.99f, -0.01f, 11.01f, 0.01f);
    check("camera_inside", aabb_classify(&pt_center, &cam) != BC_OUTSIDE);
    check("camera_fov_side", aabb_classify(&pt_side, &cam) == BC_OUTSIDE);
    check("camera_near", aabb_classify(&pt_behind, &cam) == BC_OUTSIDE);
    check("camera_far", aabb_classify(&pt_far, &cam) == BC_OUTSIDE);

    /* ---- frame 1: frustum covers everything, root passes fully inside ---- */
    HcCache cache;
    memset(&cache, 0, sizeof(cache));
    int vis[HC_MAX_NODES], ref[HC_MAX_NODES];
    HcStats st;
    Frustum2 wide;
    frustum_rect(&wide, -2, -2, 8, 6);
    int rc = hc_cull_frame(&t, 1, &wide, &cache, vis, &st);
    check("frame1_ok", rc == HC_OK);
    check("all_visible", st.hier_visible == 6 && st.ref_visible == 6);
    check("root_shortcut", st.nodes_tested == 1 && st.subtree_shortcuts == 5);
    check("no_disagreement", st.missed_cull == 0 && st.wrong_cull == 0);

    /* ---- frame 2: only chain A region visible; chain B pruned high up ---- */
    Frustum2 side;
    frustum_rect(&side, 2.5f, -2, 8, 2);
    rc = hc_cull_frame(&t, 2, &side, &cache, vis, &st);
    check("frame2_cache_coherent", rc == HC_OK);
    hc_reference(&t, &side, ref);
    int b_culled = vis[3] == HC_CULLED && vis[4] == HC_CULLED;
    check("branch_pruned", b_culled);
    check("fewer_tests_than_bones", st.nodes_tested < t.count);
    check("no_wrong_cull", st.wrong_cull == 0 && st.wrong_bound == 0);
    /* root's own box is outside while its merged box intersects: the classic
     * conservative missed cull, and it must respect the structural bound */
    check("missed_bounded",
          st.missed_cull >= 1 && st.missed_cull <= st.missed_bound);
    printf("stats side: ref=%d hier=%d missed=%d (bound %d) wrong=%d tests=%d\n",
           st.ref_visible, st.hier_visible, st.missed_cull, st.missed_bound,
           st.wrong_cull, st.nodes_tested);

    /* ---- frustum on node 4 only: parent 3 over-covers -> missed, never wrong */
    Frustum2 tip;
    frustum_rect(&tip, -2, 3.5f, 2, 6);
    rc = hc_cull_frame(&t, 3, &tip, &cache, vis, &st);
    check("frame3_cache_coherent", rc == HC_OK);
    check("tip_visible", vis[4] == HC_VISIBLE);
    check("overcover_missed_only",
          st.wrong_cull == 0 && st.missed_cull >= 1
          && st.missed_cull <= st.missed_bound);

    /* ---- empty space: root outside prunes the whole skeleton in one test ---- */
    Frustum2 nowhere;
    frustum_rect(&nowhere, 100, 100, 101, 101);
    rc = hc_cull_frame(&t, 4, &nowhere, &cache, vis, &st);
    check("empty_frame", rc == HC_OK && st.hier_visible == 0
                        && st.nodes_tested == 1 && st.wrong_cull == 0);

    /* ---- same view repeated: unchanged boxes and conclusions stay coherent ---- */
    Frustum2 side2;
    frustum_rect(&side2, 2.5f, -2, 8, 2);
    rc = hc_cull_frame(&t, 2, &side2, &cache, vis, &st);
    check("same_view_coherent", rc == HC_OK);
    hc_rebuild(&t);
    hc_cull(&t, &side2, vis, &st);
    int cache_bad = -1;
    check("cache_fresh_matches",
          hc_cache_coherence(&cache, 2, &t, vis, &cache_bad) == HC_OK);
    cache.vis[1] = HC_CULLED;   /* corrupt a cached conclusion, box unchanged */
    int coerce = hc_cache_coherence(&cache, 2, &t, vis, &cache_bad);
    check("cache_inconsistency_detected",
          coerce == HC_ERR_CACHE && cache_bad == 1);

    /* ---- moved bone: changed box invalidates the cache entry, no false alarm */
    t.local[5] = aabb_make(20, 20, 21, 21);
    hc_rebuild(&t);
    hc_cull(&t, &side2, vis, &st);
    check("move_invalidates_cache",
          hc_cache_coherence(&cache, 2, &t, vis, &cache_bad) == HC_OK);
    check("moved_bone_culled", vis[5] == HC_CULLED);

    /* ---- capacity guard ---- */
    HcTree full;
    hc_init(&full);
    int last = -1, cap_rc = 0;
    for (int i = 0; i <= HC_MAX_NODES; i++) {
        cap_rc = hc_add(&full, last, aabb_make(0, 0, 1, 1));
        if (cap_rc >= 0) last = cap_rc;
    }
    check("capacity", cap_rc == -HC_ERR_CAPACITY && full.count == HC_MAX_NODES);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
