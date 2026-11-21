
/* retarget: skeleton retargeting across renamed rigs.
   Mapping by hierarchy path + bone-length chain ratios, conflicts returned as
   candidate lists, out-of-tolerance ratio drift reported instead of forced.
   Post-mapping: bone-length recovery and mesh edge-length checks against a
   hand-authored mapping. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define RT_MAX_BONES 24
#define RT_MAX_PATH  128
#define RT_NAME      32
#define RT_TOKENS    8
#define RT_MAX_CAND  6
#define RT_MAX_VERTS 64
#define RT_MAX_EDGES 96

/* ---- tunable bounds (the upper bounds the pipeline is checked against) ---- */
#define RATIO_THRESH     0.10f   /* chain-ratio deviation above this is reported, not mapped */
#define RATIO_WARN       0.05f   /* mapped but flagged for review */
#define SCORE_MARGIN     0.15f   /* winner must beat runner-up by this much, else conflict */
#define RESID_BOUND      1e-5f   /* max relative length residual after recovery */
#define RIGID_RATIO_BND  (1.0f + 1e-4f)  /* same-bone edges must stay rigid */
#define EDGE_RATIO_BND   1.05f   /* overall max edge length ratio, mild pose */
#define DISTORT_BOUND    0.0f    /* auto vs manual edge distortion must be exactly zero */

typedef struct { float x, y; } V2;

static V2 vrot(float x, float y, float a) {
    float c = cosf(a), s = sinf(a);
    V2 r = { x * c - y * s, x * s + y * c };
    return r;
}

/* ------------------------------- rig model -------------------------------- */

typedef struct {
    char name[RT_NAME];
    int parent;
    float rest_ang;
    float len;          /* rest bone length, joint to joint */
    float cur_len;      /* working length, restored by rt_recover() */
    float local_ang;
    /* forward kinematics output */
    float wx, wy, wang;
    /* hierarchy path descriptor */
    int depth;
    int chain[RT_MAX_BONES];   /* root..self indices */
    int chain_n;
    float ratio[RT_MAX_BONES]; /* len[k] / chain total, aligned to chain */
    int leaf;
} RtBone;

typedef struct {
    RtBone b[RT_MAX_BONES];
    int count;
} Rig;

static void rig_init(Rig* r) { memset(r, 0, sizeof(*r)); }

static int rig_add(Rig* r, int parent, const char* name, float len, float rest_ang) {
    if (r->count >= RT_MAX_BONES) return -1;
    if (parent < -1 || parent >= r->count) return -1;
    if (len < 0.0f) return -1;
    int i = r->count++;
    RtBone* bn = &r->b[i];
    memset(bn, 0, sizeof(*bn));
    strncpy(bn->name, name, RT_NAME - 1);
    bn->parent = parent;
    bn->len = bn->cur_len = len;
    bn->rest_ang = bn->local_ang = rest_ang;
    if (parent < 0) {
        bn->depth = 0;
        bn->chain_n = 1;
        bn->chain[0] = i;
    } else {
        const RtBone* p = &r->b[parent];
        bn->depth = p->depth + 1;
        bn->chain_n = p->chain_n + 1;
        memcpy(bn->chain, p->chain, p->chain_n * sizeof(int));
        bn->chain[p->chain_n] = i;
        r->b[parent].leaf = 0;
    }
    bn->leaf = 1;
    float total = 0.0f;
    for (int k = 0; k < bn->chain_n; k++) total += r->b[bn->chain[k]].len;
    for (int k = 0; k < bn->chain_n; k++)
        bn->ratio[k] = total > 0.0f ? r->b[bn->chain[k]].len / total : 0.0f;
    return i;
}

static void rig_path(const Rig* r, int i, char* buf, int n) {
    const RtBone* bn = &r->b[i];
    buf[0] = '\0';
    for (int k = 0; k < bn->chain_n; k++) {
        if (k) strncat(buf, "/", n - (int)strlen(buf) - 1);
        strncat(buf, r->b[bn->chain[k]].name, n - (int)strlen(buf) - 1);
    }
}

/* local rotations absolute per bone; root sits at origin */
static void rig_pose(Rig* r, const float* ang) {
    for (int i = 0; i < r->count; i++) r->b[i].local_ang = ang[i];
}

static void rig_fk(Rig* r) {
    for (int i = 0; i < r->count; i++) {
        RtBone* bn = &r->b[i];
        if (bn->parent < 0) {
            bn->wx = 0.0f; bn->wy = 0.0f; bn->wang = bn->local_ang;
        } else {
            RtBone* p = &r->b[bn->parent];
            V2 off = vrot(bn->cur_len, 0.0f, p->wang);
            bn->wx = p->wx + off.x;
            bn->wy = p->wy + off.y;
            bn->wang = p->wang + bn->local_ang;
        }
    }
}

/* max relative deviation between two chains' length-ratio profiles */
static float chain_ratio_dev(const Rig* a, int ai, const Rig* b, int bi) {
    const RtBone* x = &a->b[ai];
    const RtBone* y = &b->b[bi];
    if (x->chain_n != y->chain_n) return 1.0f;
    float dev = 0.0f;
    for (int k = 0; k < x->chain_n; k++)
        dev = fmaxf(dev, fabsf(x->ratio[k] - y->ratio[k]));
    return dev;
}

/* --------------------------- name/path similarity -------------------------- */

static int tokenize(const char* s, char out[RT_TOKENS][16]) {
    int n = 0;
    char tmp[RT_NAME];
    strncpy(tmp, s, RT_NAME - 1);
    tmp[RT_NAME - 1] = '\0';
    char* p = strtok(tmp, "_");
    while (p && n < RT_TOKENS) {
        int k = 0;
        for (; p[k] && k < 15; k++) {
            char c = p[k];
            if (c >= 'A' && c <= 'Z') c = (char)(c - 'A' + 'a');
            out[n][k] = c;
        }
        out[n][k] = '\0';
        n++;
        p = strtok(NULL, "_");
    }
    return n;
}

static float name_sim(const char* a, const char* b) {
    char ta[RT_TOKENS][16], tb[RT_TOKENS][16];
    int na = tokenize(a, ta), nb = tokenize(b, tb);
    if (na == 0 || nb == 0) return 0.0f;
    int inter = 0;
    for (int i = 0; i < na; i++)
        for (int j = 0; j < nb; j++)
            if (strcmp(ta[i], tb[j]) == 0) inter++;
    int uni = na + nb - inter;
    return uni > 0 ? (float)inter / (float)uni : 0.0f;
}

/* -------------------------------- mapping --------------------------------- */

#define RT_MAP_PATH    0   /* exact hierarchy path */
#define RT_MAP_UNIQUE  1   /* unique best candidate */
#define RT_MAP_CONFLICT 2  /* several indistinguishable candidates */
#define RT_MAP_RATIO   3   /* ratio drift over threshold, refused */
#define RT_MAP_NONE    4   /* no structural candidate at all */

#define RT_ISS_CONFLICT 1
#define RT_ISS_RATIO   2
#define RT_ISS_WARN    3
#define RT_ISS_NOMATCH 4

typedef struct {
    int dst;
    int reason;
    float dev;
    int cand[RT_MAX_CAND];
    float cand_score[RT_MAX_CAND];
    int cand_n;
} RtMap;

typedef struct {
    int src;
    int code;
    float dev;
    int cand[RT_MAX_CAND];
    int cand_n;
} RtIssue;

typedef struct {
    RtMap m[RT_MAX_BONES];
    RtIssue issues[RT_MAX_BONES * 2];
    int issue_n;
} RtMapping;

typedef struct { int idx; float score; float dev; } Cand;

static void add_issue(RtMapping* out, int src, int code, float dev,
                      const int* cand, int cand_n) {
    RtIssue* is = &out->issues[out->issue_n++];
    is->src = src; is->code = code; is->dev = dev;
    is->cand_n = cand_n < RT_MAX_CAND ? cand_n : RT_MAX_CAND;
    for (int k = 0; k < is->cand_n; k++) is->cand[k] = cand[k];
}

static void build_mapping(const Rig* src, const Rig* tgt, RtMapping* out) {
    memset(out, 0, sizeof(*out));
    for (int i = 0; i < src->count; i++) {
        RtMap* mp = &out->m[i];
        mp->dst = -1;
        char sp[RT_MAX_PATH], tp[RT_MAX_PATH];
        rig_path(src, i, sp, sizeof(sp));

        /* exact hierarchy path wins outright */
        for (int j = 0; j < tgt->count; j++) {
            rig_path(tgt, j, tp, sizeof(tp));
            if (strcmp(sp, tp) == 0) {
                mp->dst = j;
                mp->reason = RT_MAP_PATH;
                mp->dev = chain_ratio_dev(src, i, tgt, j);
                break;
            }
        }
        if (mp->dst >= 0) continue;

        Cand all[RT_MAX_BONES];
        int all_n = 0, ok_n = 0;
        float best_blocked_dev = 0.0f;
        for (int j = 0; j < tgt->count; j++) {
            const RtBone* tb = &tgt->b[j];
            if (tb->depth != src->b[i].depth) continue;
            if (tb->chain_n != src->b[i].chain_n) continue;
            if (tb->leaf != src->b[i].leaf) continue;
            float dev = chain_ratio_dev(src, i, tgt, j);
            float score = 2.0f * name_sim(src->b[i].name, tb->name) - dev;
            all[all_n].idx = j;
            all[all_n].score = score;
            all[all_n].dev = dev;
            all_n++;
            best_blocked_dev = fmaxf(best_blocked_dev, dev);
            if (dev <= RATIO_THRESH) ok_n++;
        }
        /* insertion sort by score desc */
        for (int a = 1; a < all_n; a++) {
            Cand key = all[a];
            int b = a - 1;
            while (b >= 0 && all[b].score < key.score) { all[b + 1] = all[b]; b--; }
            all[b + 1] = key;
        }

        if (all_n == 0) {
            mp->reason = RT_MAP_NONE;
            add_issue(out, i, RT_ISS_NOMATCH, 0.0f, NULL, 0);
            continue;
        }
        if (ok_n == 0) {
            /* every structural match breaks the ratio budget: report, do not force */
            mp->reason = RT_MAP_RATIO;
            mp->dev = best_blocked_dev;
            int ids[RT_MAX_CAND], n = 0;
            for (int k = 0; k < all_n && n < RT_MAX_CAND; k++) ids[n++] = all[k].idx;
            add_issue(out, i, RT_ISS_RATIO, best_blocked_dev, ids, n);
            continue;
        }

        int n = 0;
        for (int k = 0; k < all_n && n < RT_MAX_CAND; k++) {
            if (all[k].dev > RATIO_THRESH) continue;
            mp->cand[n] = all[k].idx;
            mp->cand_score[n] = all[k].score;
            n++;
        }
        mp->cand_n = n;
        mp->dev = all[0].dev;

        int unique = n == 1 ||
                     (all[0].score - all[1].score >= SCORE_MARGIN && all[0].dev <= RATIO_THRESH);
        if (unique) {
            mp->dst = all[0].idx;
            mp->reason = RT_MAP_UNIQUE;
            if (mp->dev > RATIO_WARN)
                add_issue(out, i, RT_ISS_WARN, mp->dev, NULL, 0);
        } else {
            mp->reason = RT_MAP_CONFLICT;
            add_issue(out, i, RT_ISS_CONFLICT, mp->dev, mp->cand, mp->cand_n);
        }
    }
}

/* --------------------------- length recovery ------------------------------- */

/* restore rest lengths on the target rig, recompute FK, return the max
   relative residual between actual joint spacing and rest length */
static float rt_recover(Rig* tgt) {
    for (int i = 0; i < tgt->count; i++) tgt->b[i].cur_len = tgt->b[i].len;
    rig_fk(tgt);
    float resid = 0.0f;
    for (int i = 0; i < tgt->count; i++) {
        int p = tgt->b[i].parent;
        if (p < 0 || tgt->b[i].len <= 0.0f) continue;
        float dx = tgt->b[i].wx - tgt->b[p].wx;
        float dy = tgt->b[i].wy - tgt->b[p].wy;
        float actual = sqrtf(dx * dx + dy * dy);
        resid = fmaxf(resid, fabsf(actual - tgt->b[i].len) / tgt->b[i].len);
    }
    return resid;
}

/* ------------------------------ skin mesh ---------------------------------- */

typedef struct {
    int bone;
    float lx, ly;            /* local coords in bind frame */
    float bx, by, bang;      /* bind frame */
} RtVert;

typedef struct { int a, b; } RtEdge;

typedef struct {
    RtVert v[RT_MAX_VERTS];
    RtEdge e[RT_MAX_EDGES];
    int nv, ne;
} RtMesh;

static int mesh_add_vert(RtMesh* m, const Rig* bind, int bone, float wx, float wy) {
    if (m->nv >= RT_MAX_VERTS) return -1;
    const RtBone* bn = &bind->b[bone];
    float dx = wx - bn->wx, dy = wy - bn->wy;
    V2 loc = vrot(dx, dy, -bn->wang);
    RtVert* v = &m->v[m->nv];
    v->bone = bone;
    v->lx = loc.x; v->ly = loc.y;
    v->bx = bn->wx; v->by = bn->wy; v->bang = bn->wang;
    return m->nv++;
}

static int mesh_add_edge(RtMesh* m, int a, int b) {
    if (m->ne >= RT_MAX_EDGES) return -1;
    m->e[m->ne].a = a; m->e[m->ne].b = b;
    return m->ne++;
}

static void mesh_skin(const RtMesh* m, const Rig* posed, V2* out) {
    for (int i = 0; i < m->nv; i++) {
        const RtVert* v = &m->v[i];
        const RtBone* bn = &posed->b[v->bone];
        V2 d = vrot(v->lx, v->ly, bn->wang - v->bang);
        out[i].x = bn->wx + d.x;
        out[i].y = bn->wy + d.y;
    }
}

static float edge_len(const V2* p, const RtEdge* e) {
    float dx = p[e->a].x - p[e->b].x, dy = p[e->a].y - p[e->b].y;
    return sqrtf(dx * dx + dy * dy);
}

/* max length ratio of posed edges against the same mesh's bind edges */
static float edge_max_ratio(const RtMesh* m, const Rig* posed, const Rig* bind) {
    V2 pp[RT_MAX_VERTS], bp[RT_MAX_VERTS];
    mesh_skin(m, posed, pp);
    mesh_skin(m, bind, bp);
    float worst = 1.0f;
    for (int i = 0; i < m->ne; i++) {
        float d = edge_len(pp, &m->e[i]), ref = edge_len(bp, &m->e[i]);
        if (ref <= 0.0f) continue;
        worst = fmaxf(worst, fmaxf(d / ref, ref / d));
    }
    return worst;
}

/* max absolute edge-length difference between two posed skin results */
static float edge_distortion(const RtMesh* m, const Rig* ra, const Rig* rb) {
    V2 pa[RT_MAX_VERTS], pb[RT_MAX_VERTS];
    mesh_skin(m, ra, pa);
    mesh_skin(m, rb, pb);
    float worst = 0.0f;
    for (int i = 0; i < m->ne; i++)
        worst = fmaxf(worst, fabsf(edge_len(pa, &m->e[i]) - edge_len(pb, &m->e[i])));
    return worst;
}

/* per-edge max ratio restricted to edges whose two verts share a bone */
static float rigid_edge_ratio(const RtMesh* m, const Rig* posed, const Rig* bind) {
    V2 pp[RT_MAX_VERTS], bp[RT_MAX_VERTS];
    mesh_skin(m, posed, pp);
    mesh_skin(m, bind, bp);
    float worst = 1.0f;
    for (int i = 0; i < m->ne; i++) {
        const RtEdge* e = &m->e[i];
        if (m->v[e->a].bone != m->v[e->b].bone) continue;
        float d = edge_len(pp, e), ref = edge_len(bp, e);
        if (ref <= 0.0f) continue;
        worst = fmaxf(worst, fmaxf(d / ref, ref / d));
    }
    return worst;
}

/* --------------------------------- tests ----------------------------------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static const char* reason_name(int r) {
    static const char* names[] = { "path", "unique", "conflict", "ratio_block", "no_match" };
    return names[r];
}

/* old naming: Hips/LeftUpLeg/LeftLeg/LeftFoot, uniform unit proportions */
static Rig build_src_rig(void) {
    Rig r; rig_init(&r);
    int hip  = rig_add(&r, -1, "Hips",      0.50f, 0.0f);
    int thigh = rig_add(&r, hip, "LeftUpLeg", 2.00f, 0.0f);
    int shin  = rig_add(&r, thigh, "LeftLeg",  1.90f, 0.0f);
    rig_add(&r, shin, "LeftFoot", 0.60f, 0.0f);
    return r;
}

/* new naming, same topology, uniform 1.1x scale -> identical ratios */
static Rig build_clean_rig(void) {
    Rig r; rig_init(&r);
    int hip  = rig_add(&r, -1, "pelvis",  0.55f, 0.0f);
    int thigh = rig_add(&r, hip, "thigh_l", 2.20f, 0.0f);
    int shin  = rig_add(&r, thigh, "calf_l", 2.09f, 0.0f);
    rig_add(&r, shin, "foot_l", 0.66f, 0.0f);
    return r;
}

/* two indistinguishable legs -> conflict candidates */
static Rig build_conflict_rig(void) {
    Rig r; rig_init(&r);
    int hip = rig_add(&r, -1, "pelvis", 0.55f, 0.0f);
    int ta = rig_add(&r, hip, "legA_thigh", 2.20f, 0.0f);
    int sa = rig_add(&r, ta, "legA_shin", 2.09f, 0.0f);
    rig_add(&r, sa, "legA_foot", 0.66f, 0.0f);
    int tb = rig_add(&r, hip, "legB_thigh", 2.20f, 0.0f);
    int sb = rig_add(&r, tb, "legB_shin", 2.09f, 0.0f);
    rig_add(&r, sb, "legB_foot", 0.66f, 0.0f);
    return r;
}

/* same topology, unique names, but proportions break the ratio budget */
static Rig build_drift_rig(void) {
    Rig r; rig_init(&r);
    int hip = rig_add(&r, -1, "pelvis_drift", 0.55f, 0.0f);
    int thigh = rig_add(&r, hip, "femur_drift", 1.10f, 0.0f);
    int shin = rig_add(&r, thigh, "tibia_drift", 3.30f, 0.0f);
    rig_add(&r, shin, "foot_drift", 0.33f, 0.0f);
    return r;
}

/* mesh of per-bone vertex chains plus cross-bone edges across each joint.
   Vertices sit inset from the joint so cross edges have a non-zero bind
   length; a fixed sideways offset makes the ribbon non-degenerate. */
static RtMesh build_leg_mesh(const Rig* bind) {
    RtMesh m;
    memset(&m, 0, sizeof(m));
    const float side = 0.04f, inset = 0.20f;
    for (int i = 0; i < bind->count; i++) {
        const RtBone* bn = &bind->b[i];
        float a = fminf(inset, bn->len * 0.25f);
        V2 p0 = vrot(a, side, bn->wang);
        V2 p1 = vrot(bn->len - a, side, bn->wang);
        mesh_add_vert(&m, bind, i, bn->wx + p0.x, bn->wy + p0.y);
        mesh_add_vert(&m, bind, i, bn->wx + p1.x, bn->wy + p1.y);
    }
    for (int i = 0; i < bind->count; i++) {
        mesh_add_edge(&m, 2 * i, 2 * i + 1);               /* rigid, same bone */
        if (i > 0) mesh_add_edge(&m, 2 * (i - 1) + 1, 2 * i); /* cross the joint */
    }
    return m;
}

int main(void) {
    Rig src = build_src_rig();
    Rig clean = build_clean_rig();
    rig_fk(&src);
    rig_fk(&clean);

    /* ---- case 1: renamed but proportional rig maps fully and uniquely ---- */
    RtMapping mp;
    build_mapping(&src, &clean, &mp);
    int all_mapped = 1, all_unique = 1;
    float worst_dev = 0.0f;
    for (int i = 0; i < src.count; i++) {
        all_mapped = all_mapped && mp.m[i].dst == i;
        all_unique = all_unique &&
                     (mp.m[i].reason == RT_MAP_UNIQUE || mp.m[i].reason == RT_MAP_PATH);
        worst_dev = fmaxf(worst_dev, mp.m[i].dev);
    }
    check("clean_full_map", all_mapped);
    check("clean_no_conflict", all_unique);
    check("clean_ratio_zero", worst_dev < 1e-6f);
    check("clean_no_issues", mp.issue_n == 0);
    for (int i = 0; i < src.count; i++)
        printf("map[%s] -> %s reason=%s dev=%.5f\n",
               src.b[i].name, clean.b[mp.m[i].dst].name,
               reason_name(mp.m[i].reason), mp.m[i].dev);

    /* ---- case 2: two equal legs produce conflicts with candidates -------- */
    Rig conflict = build_conflict_rig();
    rig_fk(&conflict);
    RtMapping mc;
    build_mapping(&src, &conflict, &mc);
    int hip_ok = mc.m[0].dst == 0;
    int leg_conflicts = 1, two_cands = 1, conflict_issued = 0;
    for (int i = 1; i < src.count; i++) {
        leg_conflicts = leg_conflicts && mc.m[i].reason == RT_MAP_CONFLICT && mc.m[i].dst < 0;
        two_cands = two_cands && mc.m[i].cand_n == 2;
    }
    for (int k = 0; k < mc.issue_n; k++)
        if (mc.issues[k].code == RT_ISS_CONFLICT) conflict_issued++;
    check("conflict_root_mapped", hip_ok);
    check("conflict_legs_refused", leg_conflicts);
    check("conflict_candidates", two_cands);
    check("conflict_reported", conflict_issued == 3);
    printf("conflict thigh candidates: %d,%d scores=%.3f,%.3f\n",
           mc.m[1].cand[0], mc.m[1].cand[1],
           mc.m[1].cand_score[0], mc.m[1].cand_score[1]);

    /* ---- case 3: ratio drift over threshold is reported, not forced ------ */
    Rig drift = build_drift_rig();
    rig_fk(&drift);
    RtMapping md;
    build_mapping(&src, &drift, &md);
    int ratio_issues = 0, none_forced = 1;
    for (int k = 0; k < md.issue_n; k++)
        if (md.issues[k].code == RT_ISS_RATIO) ratio_issues++;
    for (int i = 0; i < src.count; i++)
        if (md.m[i].reason == RT_MAP_RATIO) none_forced = none_forced && md.m[i].dst < 0;
    check("ratio_drift_reported", ratio_issues >= 1);
    check("ratio_drift_not_forced", none_forced);
    for (int k = 0; k < md.issue_n; k++)
        if (md.issues[k].code == RT_ISS_RATIO)
            printf("ratio block on src bone %d: best dev=%.4f threshold=%.2f\n",
                   md.issues[k].src, md.issues[k].dev, RATIO_THRESH);

    /* ---- case 4: bone-length recovery after a bad scale pull ------------- */
    Rig pulled = build_clean_rig();
    rig_fk(&pulled);
    float pose_ang[RT_MAX_BONES] = { 0.10f, 0.30f, 0.50f, -0.20f };
    rig_pose(&pulled, pose_ang);
    pulled.b[1].cur_len *= 1.35f;   /* the "whole leg pulled sideways" failure */
    pulled.b[2].cur_len *= 0.70f;
    rig_fk(&pulled);
    float pre = 0.0f;
    for (int i = 1; i < pulled.count; i++) {
        int p = pulled.b[i].parent;
        float dx = pulled.b[i].wx - pulled.b[p].wx;
        float dy = pulled.b[i].wy - pulled.b[p].wy;
        pre = fmaxf(pre, fabsf(sqrtf(dx * dx + dy * dy) - pulled.b[i].len) / pulled.b[i].len);
    }
    float resid = rt_recover(&pulled);
    check("pull_detected_before_recover", pre > 0.30f);
    check("recover_residual_bound", resid <= RESID_BOUND);
    printf("length residual: pre=%.5f post=%.6f bound=%.0e\n", pre, resid, RESID_BOUND);

    /* ---- case 5: auto mapping vs hand mapping, edge distortion == 0 ------ */
    /* mild pose: largest inter-bone bend 0.15 rad so straight-at-bind cross
       edges stay inside the overall ratio bound */
    float src_pose[RT_MAX_BONES] = { 0.05f, 0.10f, 0.15f, -0.08f };
    rig_pose(&src, src_pose);
    rig_fk(&src);

    Rig auto_rig = build_clean_rig();
    Rig manual_rig = build_clean_rig();
    float mapped_ang[RT_MAX_BONES];
    for (int i = 0; i < src.count; i++) mapped_ang[mp.m[i].dst] = src_pose[i];
    rig_pose(&auto_rig, mapped_ang);
    rig_fk(&auto_rig);

    /* hand-authored mapping happens to be the identity permutation here */
    int hand[RT_MAX_BONES] = { 0, 1, 2, 3 };
    float hand_ang[RT_MAX_BONES];
    for (int i = 0; i < src.count; i++) hand_ang[hand[i]] = src_pose[i];
    rig_pose(&manual_rig, hand_ang);
    rig_fk(&manual_rig);

    Rig bind_rig = build_clean_rig();
    rig_fk(&bind_rig);
    RtMesh mesh = build_leg_mesh(&bind_rig);

    float distort = edge_distortion(&mesh, &auto_rig, &manual_rig);
    check("auto_vs_manual_zero_distortion", distort <= DISTORT_BOUND);
    printf("edge distortion auto vs manual = %.8f bound = %.1f\n", distort, DISTORT_BOUND);

    float rigid_r = rigid_edge_ratio(&mesh, &auto_rig, &bind_rig);
    float overall_r = edge_max_ratio(&mesh, &auto_rig, &bind_rig);
    check("rigid_edge_ratio_bound", rigid_r <= RIGID_RATIO_BND);
    check("overall_edge_ratio_bound", overall_r <= EDGE_RATIO_BND);
    printf("edge max ratio: rigid=%.6f (bound %.6f) overall=%.6f (bound %.2f)\n",
           rigid_r, RIGID_RATIO_BND, overall_r, EDGE_RATIO_BND);

    /* ---- case 6: a wrong mapping is caught by the same edge check -------- */
    Rig bad_rig = build_clean_rig();
    float bad_ang[RT_MAX_BONES];
    int perm[RT_MAX_BONES] = { 0, 2, 1, 3 };   /* thigh and calf swapped */
    for (int i = 0; i < src.count; i++) bad_ang[perm[i]] = src_pose[i];
    rig_pose(&bad_rig, bad_ang);
    rig_fk(&bad_rig);
    float bad_r = edge_max_ratio(&mesh, &bad_rig, &bind_rig);
    float bad_distort = edge_distortion(&mesh, &bad_rig, &manual_rig);
    /* the wrong permutation diverges from manual ground truth even though the
       pose is mild; any nonzero distortion is the detection signal */
    check("bad_mapping_detected", bad_distort > 1e-7f);
    printf("wrong-map edge ratio=%.6f distortion vs manual=%.6f\n", bad_r, bad_distort);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
