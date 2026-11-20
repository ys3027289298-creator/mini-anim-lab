
/* lod: screen-occupancy driven bone LOD.
 *
 * Bones are ranked by screen occupancy scaled by a configurable importance
 * weight; distal bones are merged into their nearest retained ancestor, up to
 * a per distance-tier merge level. Skinning weights are compacted onto the
 * retained bones and must keep summing to one. The LOD pose is compared
 * against the full-skeleton pose, with end-effector error accumulated per
 * camera-distance tier. When LOD is disabled the output pose is the full pose.
 */
#include <stdio.h>
#include <string.h>
#include <math.h>

#define MAX_BONES 16
#define MAX_VERTS 8
#define MAX_TIERS 4
#define LOD_EPS 1e-5f

enum {
    LOD_OK = 0,
    LOD_ERR_NULL = -1,
    LOD_ERR_PARAM = -2,
    LOD_ERR_PARENT = -3,
    LOD_ERR_SCREEN = -4,
    LOD_ERR_LEVELS = -5,
    LOD_ERR_WEIGHT = -6,
    LOD_ERR_TIER = -7,
    LOD_ERR_DIST = -8,
    LOD_ERR_SUM = -9
};

typedef struct { float x, y, ang; } Transform;

typedef struct {
    int parent[MAX_BONES];
    Transform local[MAX_BONES];
    Transform world[MAX_BONES];
    int count;
} LodSkeleton;

typedef struct {
    float maxDist;      /* upper bound of camera distance; last tier may be INFINITY */
    int mergeLevels;    /* leaf levels allowed to collapse; 1..root leaf distance */
    float maxError;     /* configurable end-effector error budget */
} LodTier;

typedef struct {
    int enabled;
    float importance[MAX_BONES];   /* configurable ranking weight per bone */
    LodTier tiers[MAX_TIERS];
    int tierCount;
} LodConfig;

typedef struct {
    int enabled;
    int tier;                      /* selected distance tier, -1 when disabled */
    int activeCount;               /* retained bones */
    int retained[MAX_BONES];       /* 1 if the bone stays active */
    int target[MAX_BONES];         /* merge target; retained bones target self */
    int rankOrder[MAX_BONES];      /* bones sorted by screen occupancy desc */
    Transform world[MAX_BONES];    /* LOD world pose */
    float weights[MAX_VERTS][MAX_BONES]; /* compacted skinning weights */
} LodLevel;

typedef struct {
    int samples;
    float maxError;
    float sumError;
    int overBudget;
} LodTierReport;

typedef struct {
    LodTierReport tiers[MAX_TIERS];
    int tierCount;
} LodReport;

static void rotate(float x, float y, float a, float* ox, float* oy) {
    float c = cosf(a), s = sinf(a);
    *ox = x * c - y * s;
    *oy = x * s + y * c;
}

/* forward kinematics; parents always precede children in index order */
static void fk(const LodSkeleton* sk, Transform* out) {
    for (int i = 0; i < sk->count; i++) {
        if (sk->parent[i] < 0) {
            out[i] = sk->local[i];
        } else {
            int p = sk->parent[i];
            float rx, ry;
            rotate(sk->local[i].x, sk->local[i].y, out[p].ang, &rx, &ry);
            out[i].x = out[p].x + rx;
            out[i].y = out[p].y + ry;
            out[i].ang = out[p].ang + sk->local[i].ang;
        }
    }
}

/* distance in edges from each bone to its nearest descendant leaf */
static void leaf_distance(const LodSkeleton* sk, int* leafDist) {
    int hasChild[MAX_BONES] = {0};
    for (int i = 0; i < sk->count; i++)
        if (sk->parent[i] >= 0) hasChild[sk->parent[i]] = 1;
    for (int i = sk->count - 1; i >= 0; i--) {
        if (!hasChild[i]) {
            leafDist[i] = 0;
            continue;
        }
        int best = MAX_BONES;
        for (int c = i + 1; c < sk->count; c++)
            if (sk->parent[c] == i && leafDist[c] < best) best = leafDist[c];
        leafDist[i] = best + 1;
    }
}

static int pick_tier(const LodConfig* cfg, float camDist) {
    for (int t = 0; t < cfg->tierCount; t++)
        if (camDist <= cfg->tiers[t].maxDist) return t;
    return -1;
}

static int lod_build(const LodSkeleton* sk, const float* screen,
                     float (*inWeights)[MAX_BONES], int vertCount,
                     const LodConfig* cfg, float camDist,
                     LodLevel* out) {
    if (!sk || !cfg || !out) return LOD_ERR_NULL;
    if (sk->count <= 0 || sk->count > MAX_BONES) return LOD_ERR_PARAM;
    if (vertCount < 0 || vertCount > MAX_VERTS) return LOD_ERR_PARAM;
    if (vertCount > 0 && !inWeights) return LOD_ERR_NULL;

    for (int i = 0; i < sk->count; i++) {
        int p = sk->parent[i];
        if (p >= i || (p < 0 && i != 0) || p < -1) return LOD_ERR_PARENT;
    }

    memset(out, 0, sizeof(*out));
    for (int i = 0; i < sk->count; i++) out->target[i] = i;

    Transform full[MAX_BONES];
    fk(sk, full);

    /* LOD disabled: identity mapping, pose identical to the full skeleton */
    if (!cfg->enabled) {
        out->enabled = 0;
        out->tier = -1;
        out->activeCount = sk->count;
        for (int i = 0; i < sk->count; i++) {
            out->retained[i] = 1;
            out->rankOrder[i] = i;
            out->world[i] = full[i];
        }
        for (int v = 0; v < vertCount; v++)
            memcpy(out->weights[v], inWeights[v], sizeof(float) * sk->count);
        return LOD_OK;
    }

    if (!screen) return LOD_ERR_NULL;
    for (int i = 0; i < sk->count; i++)
        if (!(screen[i] > 0.0f)) return LOD_ERR_SCREEN;
    for (int i = 0; i < sk->count; i++)
        if (!(cfg->importance[i] > 0.0f)) return LOD_ERR_PARAM;
    if (cfg->tierCount <= 0 || cfg->tierCount > MAX_TIERS) return LOD_ERR_TIER;

    int leafDist[MAX_BONES];
    leaf_distance(sk, leafDist);
    int rootLeafDist = leafDist[0];

    float prev = -1.0f;
    for (int t = 0; t < cfg->tierCount; t++) {
        const LodTier* tr = &cfg->tiers[t];
        if (!(tr->maxDist > prev)) return LOD_ERR_TIER;   /* strictly increasing; INF ok */
        if (tr->mergeLevels <= 0 || tr->mergeLevels > rootLeafDist) return LOD_ERR_LEVELS;
        if (!(tr->maxError >= 0.0f)) return LOD_ERR_TIER;
        prev = tr->maxDist;
    }

    int tier = pick_tier(cfg, camDist);
    if (tier < 0) return LOD_ERR_DIST;
    int levels = cfg->tiers[tier].mergeLevels;

    for (int v = 0; v < vertCount; v++) {
        float sum = 0.0f;
        for (int i = 0; i < sk->count; i++) sum += inWeights[v][i];
        if (fabsf(sum - 1.0f) > LOD_EPS) return LOD_ERR_WEIGHT;
    }

    /* rank by screen occupancy scaled by importance, descending; index tiebreak */
    float score[MAX_BONES];
    for (int i = 0; i < sk->count; i++) {
        score[i] = screen[i] * cfg->importance[i];
        out->rankOrder[i] = i;
    }
    for (int i = 1; i < sk->count; i++) {
        int key = out->rankOrder[i];
        int j = i - 1;
        while (j >= 0 && (score[out->rankOrder[j]] < score[key] ||
                          (score[out->rankOrder[j]] == score[key] && out->rankOrder[j] > key))) {
            out->rankOrder[j + 1] = out->rankOrder[j];
            j--;
        }
        out->rankOrder[j + 1] = key;
    }

    /* merge distal bones: leaf distance below the tier's merge level */
    for (int i = 0; i < sk->count; i++)
        out->retained[i] = (i == 0 || leafDist[i] >= levels) ? 1 : 0;

    for (int i = 1; i < sk->count; i++) {
        if (out->retained[i]) {
            out->target[i] = i;
        } else {
            int p = sk->parent[i];
            out->target[i] = out->retained[p] ? p : out->target[p];
        }
    }

    out->enabled = 1;
    out->tier = tier;
    for (int i = 0; i < sk->count; i++) {
        out->activeCount += out->retained[i];
        out->world[i] = out->retained[i] ? full[i] : full[out->target[i]];
    }

    /* compact skinning weights onto retained bones */
    for (int v = 0; v < vertCount; v++) {
        for (int i = 0; i < sk->count; i++)
            out->weights[v][out->target[i]] += inWeights[v][i];
        float sum = 0.0f;
        for (int i = 0; i < sk->count; i++) sum += out->weights[v][i];
        if (fabsf(sum - 1.0f) > LOD_EPS) return LOD_ERR_SUM;
    }
    return LOD_OK;
}

/* accumulate end-effector (leaf bone) position error into the tier buckets */
static void lod_compare(const LodSkeleton* sk, const LodLevel* lod,
                        const LodConfig* cfg, LodReport* rep) {
    if (!lod->enabled || lod->tier < 0) return;

    int hasChild[MAX_BONES] = {0};
    for (int i = 0; i < sk->count; i++)
        if (sk->parent[i] >= 0) hasChild[sk->parent[i]] = 1;

    LodTierReport* r = &rep->tiers[lod->tier];
    for (int i = 0; i < sk->count; i++) {
        if (hasChild[i]) continue;
        float dx = lod->world[i].x - sk->world[i].x;
        float dy = lod->world[i].y - sk->world[i].y;
        float err = hypotf(dx, dy);
        r->samples++;
        r->sumError += err;
        if (err > r->maxError) r->maxError = err;
        if (err > cfg->tiers[lod->tier].maxError) r->overBudget = 1;
    }
}

static void report_init(const LodConfig* cfg, LodReport* rep) {
    memset(rep, 0, sizeof(*rep));
    rep->tierCount = cfg->tierCount;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static LodSkeleton make_chain(void) {
    LodSkeleton sk;
    memset(&sk, 0, sizeof(sk));
    sk.count = 5;
    for (int i = 0; i < sk.count; i++) {
        sk.parent[i] = i == 0 ? -1 : i - 1;
        sk.local[i].x = 1.0f; sk.local[i].y = 0.0f; sk.local[i].ang = 0.0f;
    }
    sk.local[2].ang = 1.5707963f;  /* bend 90 deg at bone 2 */
    fk(&sk, sk.world);
    return sk;
}

static LodConfig base_config(int enabled) {
    LodConfig cfg;
    memset(&cfg, 0, sizeof(cfg));
    cfg.enabled = enabled;
    for (int i = 0; i < MAX_BONES; i++) cfg.importance[i] = 1.0f;
    cfg.tierCount = 3;
    cfg.tiers[0] = (LodTier){ 10.0f, 1, 1.0f };
    cfg.tiers[1] = (LodTier){ 30.0f, 2, 2.0f };
    cfg.tiers[2] = (LodTier){ INFINITY, 3, 2.5f };
    return cfg;
}

int main(void) {
    LodSkeleton sk = make_chain();
    float screen[5] = { 0.20f, 0.16f, 0.12f, 0.08f, 0.04f };
    float weights[2][MAX_BONES] = {
        { 0.5f, 0.3f, 0.2f, 0.0f, 0.0f },
        { 0.0f, 0.0f, 0.1f, 0.6f, 0.3f }
    };

    LodConfig off = base_config(0);
    LodLevel lod;
    int rc = lod_build(&sk, screen, weights, 2, &off, 5.0f, &lod);
    check("disabled_ok", rc == LOD_OK);
    int identical = 1, allRetained = 1, selfTarget = 1, weightsCopied = 1;
    for (int i = 0; i < sk.count; i++) {
        if (lod.world[i].x != sk.world[i].x || lod.world[i].y != sk.world[i].y ||
            lod.world[i].ang != sk.world[i].ang) identical = 0;
        if (!lod.retained[i]) allRetained = 0;
        if (lod.target[i] != i) selfTarget = 0;
    }
    for (int v = 0; v < 2; v++)
        for (int i = 0; i < sk.count; i++)
            if (lod.weights[v][i] != weights[v][i]) weightsCopied = 0;
    check("disabled_pose_identical", identical);
    check("disabled_all_retained", allRetained && lod.activeCount == 5);
    check("disabled_target_self", selfTarget);
    check("disabled_weights_copied", weightsCopied);

    LodConfig cfg = base_config(1);
    LodReport rep;
    report_init(&cfg, &rep);

    LodLevel nearLod;
    rc = lod_build(&sk, screen, weights, 2, &cfg, 5.0f, &nearLod);
    check("near_ok", rc == LOD_OK && nearLod.tier == 0);
    check("near_merge_leaf", nearLod.activeCount == 4 &&
          nearLod.retained[4] == 0 && nearLod.target[4] == 3);
    check("near_weight_redist",
          fabsf(nearLod.weights[1][3] - 0.9f) < LOD_EPS &&
          fabsf(nearLod.weights[1][4] - 0.0f) < LOD_EPS);
    int sumOne = 1;
    for (int v = 0; v < 2; v++) {
        float sum = 0.0f;
        for (int i = 0; i < sk.count; i++) sum += nearLod.weights[v][i];
        if (fabsf(sum - 1.0f) > LOD_EPS) sumOne = 0;
    }
    check("near_weights_sum_one", sumOne);
    lod_compare(&sk, &nearLod, &cfg, &rep);
    check("near_error_one_bone", fabsf(rep.tiers[0].maxError - 1.0f) < 1e-4f);

    int rankDesc = 1, permutation = 1;
    for (int i = 1; i < sk.count; i++) {
        float a = screen[nearLod.rankOrder[i - 1]], b = screen[nearLod.rankOrder[i]];
        if (a < b) rankDesc = 0;
    }
    int seen[MAX_BONES] = {0};
    for (int i = 0; i < sk.count; i++) {
        if (seen[nearLod.rankOrder[i]]) permutation = 0;
        seen[nearLod.rankOrder[i]] = 1;
    }
    check("rank_by_screen", rankDesc && nearLod.rankOrder[0] == 0 &&
          nearLod.rankOrder[4] == 4);
    check("rank_permutation", permutation);

    LodConfig weighted = cfg;
    weighted.importance[4] = 100.0f;
    LodLevel flipped;
    rc = lod_build(&sk, screen, weights, 2, &weighted, 5.0f, &flipped);
    check("rank_importance_flip", rc == LOD_OK && flipped.rankOrder[0] == 4);

    LodLevel midLod, farLod;
    rc = lod_build(&sk, screen, weights, 2, &cfg, 20.0f, &midLod);
    check("mid_ok", rc == LOD_OK && midLod.tier == 1);
    check("mid_merge_two_levels", midLod.activeCount == 3 &&
          midLod.target[3] == 2 && midLod.target[4] == 2);
    lod_compare(&sk, &midLod, &cfg, &rep);

    rc = lod_build(&sk, screen, weights, 2, &cfg, 100.0f, &farLod);
    check("far_ok", rc == LOD_OK && farLod.tier == 2);
    check("far_merge_three_levels", farLod.activeCount == 2 &&
          farLod.target[2] == 1 && farLod.target[4] == 1);
    check("retained_pose_unchanged",
          farLod.world[1].x == sk.world[1].x && farLod.world[1].y == sk.world[1].y &&
          farLod.world[1].ang == sk.world[1].ang);
    lod_compare(&sk, &farLod, &cfg, &rep);

    check("report_tier_counts", rep.tiers[0].samples == 1 &&
          rep.tiers[1].samples == 1 && rep.tiers[2].samples == 1);
    /* joint2 (2,0), joint4 (2,2): two-level err 2.0; joint1 (1,0): far err sqrt5 */
    check("report_mid_error", fabsf(rep.tiers[1].maxError - 2.0f) < 1e-4f);
    check("report_far_error", fabsf(rep.tiers[2].maxError - 2.2360680f) < 1e-4f);
    check("report_within_budget", rep.tiers[0].overBudget == 0 &&
          rep.tiers[1].overBudget == 0 && rep.tiers[2].overBudget == 0);

    LodConfig tight = cfg;
    tight.tiers[0].maxError = 0.01f;
    LodReport tightRep;
    report_init(&tight, &tightRep);
    lod_compare(&sk, &nearLod, &tight, &tightRep);
    check("over_budget_flag", tightRep.tiers[0].overBudget == 1);

    LodReport offRep;
    report_init(&off, &offRep);
    lod_compare(&sk, &lod, &off, &offRep);
    check("disabled_no_report", offRep.tiers[0].samples == 0);

    /* branched tree: merged chains must walk up to the nearest retained ancestor */
    LodSkeleton tree;
    memset(&tree, 0, sizeof(tree));
    tree.count = 6;
    int parent6[] = { -1, 0, 1, 0, 3, 4 };  /* 0->1->2 and 0->3->4->5 */
    for (int i = 0; i < 6; i++) {
        tree.parent[i] = parent6[i];
        tree.local[i].x = 1.0f;
    }
    fk(&tree, tree.world);
    float screen6[6] = { 0.30f, 0.20f, 0.10f, 0.20f, 0.12f, 0.05f };
    LodConfig treeCfg = base_config(1);
    treeCfg.tierCount = 1;
    treeCfg.tiers[0] = (LodTier){ INFINITY, 2, 1.0f };  /* root leaf distance is 2 */
    LodLevel treeLod;
    rc = lod_build(&tree, screen6, NULL, 0, &treeCfg, 1.0f, &treeLod);
    check("branch_ok", rc == LOD_OK && treeLod.activeCount == 2);
    check("branch_target_walk", treeLod.target[1] == 0 && treeLod.target[2] == 0 &&
          treeLod.target[4] == 3 && treeLod.target[5] == 3);

    /* rejection cases */
    float zeroScreen[5] = { 0.20f, 0.16f, 0.12f, 0.08f, 0.0f };
    check("reject_zero_screen",
          lod_build(&sk, zeroScreen, weights, 2, &cfg, 5.0f, &lod) == LOD_ERR_SCREEN);

    LodConfig badLevels = cfg;
    badLevels.tiers[0].mergeLevels = 0;
    check("reject_zero_levels",
          lod_build(&sk, screen, weights, 2, &badLevels, 5.0f, &lod) == LOD_ERR_LEVELS);
    badLevels = cfg;
    badLevels.tiers[0].mergeLevels = 5;   /* deeper than root leaf distance 4 */
    check("reject_levels_too_deep",
          lod_build(&sk, screen, weights, 2, &badLevels, 5.0f, &lod) == LOD_ERR_LEVELS);
    badLevels = treeCfg;
    badLevels.tiers[0].mergeLevels = 3;   /* would merge the branch root */
    check("reject_branch_root_merge",
          lod_build(&tree, screen6, NULL, 0, &badLevels, 1.0f, &lod) == LOD_ERR_LEVELS);

    float badWeights[2][MAX_BONES] = {
        { 0.5f, 0.3f, 0.1f, 0.0f, 0.0f },  /* sums to 0.9 */
        { 0.0f, 0.0f, 0.1f, 0.6f, 0.3f }
    };
    check("reject_bad_weight",
          lod_build(&sk, screen, badWeights, 2, &cfg, 5.0f, &lod) == LOD_ERR_WEIGHT);

    LodConfig finite = cfg;
    finite.tierCount = 1;
    finite.tiers[0] = (LodTier){ 10.0f, 1, 1.0f };
    check("reject_uncovered_distance",
          lod_build(&sk, screen, weights, 2, &finite, 99.0f, &lod) == LOD_ERR_DIST);

    LodSkeleton badParent = sk;
    badParent.parent[3] = 4;  /* parent after child */
    check("reject_bad_parent",
          lod_build(&badParent, screen, weights, 2, &cfg, 5.0f, &lod) == LOD_ERR_PARENT);

    check("reject_bad_vert_count",
          lod_build(&sk, screen, weights, MAX_VERTS + 1, &cfg, 5.0f, &lod) == LOD_ERR_PARAM);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
