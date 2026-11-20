
/* skin: greedy LBS influence-set pruning.
 *
 * Each vertex carries up to MAX_INF bone influences with raw, unnormalized
 * weights. We pick an influence set greedily by descending contribution
 * (weight * |displaced position|), renormalize the survivors, and at every
 * step publish an upper bound on the discarded weight and on the displacement
 * error against the exhaustive optimal set of the same size.
 *
 * Bound derivation (bind pose at origin, so bone displacement is p_i):
 *   W = sum w_i, q_i = w_i / W, full blend p = sum q_i p_i.
 *   For a kept subset S with discarded set D:
 *     p = rho_S * m_D + (1 - rho_S) * p_S,
 *     rho_S = (W - W_S) / W, m_D = (sum_D w_i p_i) / (W - W_S).
 *   With R = max |p_i| over positive-weight influences:
 *     |p_S - p| = rho_S |m_D - p| <= rho_S (|m_D| + |p|) <= 2 R rho_S.
 *   Greedy set G and optimal same-size set O (O minimizes |p_O - p|):
 *     |p_G - p_O| <= |p_G - p| + |p_O - p| <= 4 R rho_G.
 *   G_k grows with k, so rho_G (and the bound) is non-increasing and hits 0
 *   once every influence is kept. */
#include <stdio.h>
#include <string.h>
#include <math.h>

#define MAX_INF 8
#define EPS 1e-5f

typedef struct { float x, y; } V2;

static float vlen(V2 v) { return sqrtf(v.x * v.x + v.y * v.y); }
static float vdist(V2 a, V2 b) { return vlen((V2){ a.x - b.x, a.y - b.y }); }

typedef struct {
    V2 pos;      /* bone-transformed position (bind pose is the origin) */
    float w;     /* raw, unnormalized non-negative weight */
} Influence;

typedef struct {
    int idx;         /* selected influence index */
    float renorm;    /* renormalized weight inside the current set */
} Picked;

typedef struct {
    Picked pick;
    float discarded_w;     /* total raw weight cut at this step */
    float discarded_frac;  /* rho_k: discarded fraction of the full weight */
    float bound;           /* B_k = 4 R rho_k, displacement upper bound */
    V2 blend;              /* renormalized greedy blend for this set */
} GreedyStep;

static int inputs_valid(const Influence* in, int n, int max_inf) {
    if (in == NULL || n <= 0 || n > MAX_INF || max_inf <= 0) return 0;
    float total = 0.0f;
    for (int i = 0; i < n; i++) {
        if (!isfinite(in[i].w) || in[i].w < 0.0f) return 0;
        if (!isfinite(in[i].pos.x) || !isfinite(in[i].pos.y)) return 0;
        total += in[i].w;
    }
    if (!(total > 0.0f)) return 0;   /* all-zero weights cannot be normalized */
    return 1;
}

static V2 blend_mask(const Influence* in, int n, int mask) {
    float ws = 0.0f;
    V2 acc = { 0.0f, 0.0f };
    for (int i = 0; i < n; i++) {
        if (mask & (1 << i)) {
            ws += in[i].w;
            acc.x += in[i].w * in[i].pos.x;
            acc.y += in[i].w * in[i].pos.y;
        }
    }
    if (!(ws > 0.0f)) return (V2){ 0.0f, 0.0f };
    return (V2){ acc.x / ws, acc.y / ws };
}

/* contribution: how much this influence moves the vertex, raw-weight scaled */
static float contribution(const Influence* v) {
    return v->w * vlen(v->pos);
}

/* greedy ranking: contribution desc, then raw weight desc, then index asc */
static int ranks_before(const Influence* in, int a, int b) {
    float ca = contribution(&in[a]), cb = contribution(&in[b]);
    if (ca != cb) return ca > cb;
    if (in[a].w != in[b].w) return in[a].w > in[b].w;
    return a < b;
}

/* returns number of steps (>=1), or -1 on rejected input.
 * steps must hold min(max_inf, n) entries; full receives the full blend. */
static int greedy_select(const Influence* in, int n, int max_inf,
                         GreedyStep* steps, V2* full, float* out_R) {
    if (!inputs_valid(in, n, max_inf)) return -1;

    float total = 0.0f, R = 0.0f;
    for (int i = 0; i < n; i++) {
        total += in[i].w;
        if (in[i].w > 0.0f) {
            float d = vlen(in[i].pos);
            if (d > R) R = d;
        }
    }
    int full_mask = (1 << n) - 1;
    *full = blend_mask(in, n, full_mask);
    *out_R = R;

    int used[MAX_INF] = { 0 };
    int kept_mask = 0;
    int count = max_inf < n ? max_inf : n;
    for (int k = 0; k < count; k++) {
        int best = -1;
        for (int i = 0; i < n; i++) {
            if (used[i]) continue;
            if (best < 0 || ranks_before(in, i, best)) best = i;
        }
        used[best] = 1;
        kept_mask |= 1 << best;

        float kept_w = 0.0f;
        for (int i = 0; i < n; i++)
            if (kept_mask & (1 << i)) kept_w += in[i].w;

        GreedyStep* s = &steps[k];
        s->pick.idx = best;
        s->pick.renorm = kept_w > 0.0f ? in[best].w / kept_w : 0.0f;
        s->discarded_w = total - kept_w;
        s->discarded_frac = s->discarded_w / total;
        s->bound = 4.0f * R * s->discarded_frac;
        s->blend = blend_mask(in, n, kept_mask);
    }
    return count;
}

static float total_weight(const Influence* in, int n) {
    float t = 0.0f;
    for (int i = 0; i < n; i++) t += in[i].w;
    return t;
}

/* kept raw weights renormalize to one by construction; verify the weight
 * accounting and that the current pick stores its renormalized share.
 * earlier pick.renorm entries are stale once more bones join. */
static int renorm_sums_one(const Influence* in, const GreedyStep* steps,
                           int k, float total) {
    float kept = 0.0f;
    for (int j = 0; j <= k; j++) kept += in[steps[j].pick.idx].w;
    if (fabsf(kept - (total - steps[k].discarded_w)) > EPS) return 0;
    float expected = in[steps[k].pick.idx].w / kept;
    return fabsf(expected - steps[k].pick.renorm) <= EPS;
}

static int popcount(int x) {
    int c = 0;
    while (x) { c += x & 1; x >>= 1; }
    return c;
}

/* best size-k subset: smallest displacement from the full blend */
static V2 exhaustive_best(const Influence* in, int n, int k, V2 full) {
    int best_mask = -1;
    float best_d = 0.0f;
    for (int mask = 0; mask < (1 << n); mask++) {
        if (popcount(mask) != k) continue;
        float d = vdist(blend_mask(in, n, mask), full);
        if (best_mask < 0 || d < best_d) { best_mask = mask; best_d = d; }
    }
    return blend_mask(in, n, best_mask);
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

/* xorshift32: deterministic stress-test rng */
static unsigned rng_state = 0x12345678u;
static float randf(void) {
    unsigned x = rng_state;
    x ^= x << 13; x ^= x >> 17; x ^= x << 5;
    rng_state = x;
    return (float)(x & 0xffffff) / (float)0x1000000;
}

static Influence hand_case[MAX_INF];

static void build_hand_case(Influence* out) {
    out[0] = (Influence){ { 1.0f,  0.0f}, 5.0f };  /* contrib 5.00 */
    out[1] = (Influence){ { 0.0f,  2.0f}, 4.0f };  /* contrib 8.00 */
    out[2] = (Influence){ {-1.0f,  1.0f}, 3.0f };  /* contrib 4.24 */
    out[3] = (Influence){ { 2.0f, -2.0f}, 2.5f };  /* contrib 7.07 */
    out[4] = (Influence){ { 0.0f, -1.0f}, 2.0f };  /* contrib 2.00 */
    out[5] = (Influence){ { 1.5f,  0.0f}, 1.5f };  /* contrib 2.25 */
    out[6] = (Influence){ { 0.0f,  0.0f}, 3.0f };  /* contrib 0.00, heavy but static */
    out[7] = (Influence){ {-2.0f,  0.0f}, 0.1f };  /* contrib 0.20, light but moving */
}

int main(void) {
    /* ---- rejection paths ---- */
    build_hand_case(hand_case);
    GreedyStep steps[MAX_INF];
    V2 full; float R;

    check("reject_max_inf_zero",
          greedy_select(hand_case, MAX_INF, 0, steps, &full, &R) == -1);
    check("reject_max_inf_negative",
          greedy_select(hand_case, MAX_INF, -2, steps, &full, &R) == -1);

    Influence neg[MAX_INF];
    build_hand_case(neg);
    neg[3].w = -0.5f;
    check("reject_negative_weight",
          greedy_select(neg, MAX_INF, 4, steps, &full, &R) == -1);

    Influence zeros[MAX_INF];
    build_hand_case(zeros);
    for (int i = 0; i < MAX_INF; i++) zeros[i].w = 0.0f;
    check("reject_all_zero_weights",
          greedy_select(zeros, MAX_INF, 4, steps, &full, &R) == -1);

    /* ---- greedy run on the hand case ---- */
    int count = greedy_select(hand_case, MAX_INF, MAX_INF, steps, &full, &R);
    check("greedy_runs", count == MAX_INF);
    check("first_pick_top_contribution", steps[0].pick.idx == 1);
    check("contribution_beats_weight", steps[6].pick.idx == 7 &&
                                       steps[7].pick.idx == 6);

    int renorm_ok = 1, monotone_ok = 1, bound_ok = 1, single_ok = 1;
    int full_step_ok = 0;
    float hand_total = total_weight(hand_case, MAX_INF);
    printf("step discarded_w discarded_frac bound |pG-pOpt| |pG-full|\n");
    for (int k = 0; k < count; k++) {
        int set_size = k + 1;
        V2 opt = exhaustive_best(hand_case, MAX_INF, set_size, full);
        float err_opt = vdist(steps[k].blend, opt);
        float err_full = vdist(steps[k].blend, full);
        printf("%4d %11.4f %14.5f %6.4f %9.5f %9.5f\n",
               set_size, steps[k].discarded_w, steps[k].discarded_frac,
               steps[k].bound, err_opt, err_full);

        if (!renorm_sums_one(hand_case, steps, k, hand_total)) renorm_ok = 0;

        if (k > 0 && steps[k].bound > steps[k - 1].bound + EPS) monotone_ok = 0;
        if (err_opt > steps[k].bound + EPS) bound_ok = 0;
        if (err_full > 2.0f * R * steps[k].discarded_frac + EPS) single_ok = 0;

        if (set_size == MAX_INF &&
            steps[k].discarded_w == 0.0f && steps[k].bound <= EPS &&
            vdist(steps[k].blend, full) <= EPS) full_step_ok = 1;
    }
    check("kept_weights_renormalized", renorm_ok);
    check("bound_monotone_tightening", monotone_ok);
    check("error_vs_optimal_within_bound", bound_ok);
    check("single_set_inequality_2Rrho", single_ok);
    check("full_set_closes_the_bound", full_step_ok);

    /* ---- order independence: shuffle the input array ---- */
    Influence shuf[MAX_INF];
    int perm[MAX_INF];
    for (int i = 0; i < MAX_INF; i++) perm[i] = i;
    for (int i = MAX_INF - 1; i > 0; i--) {
        int j = (int)(randf() * (i + 1));
        int t = perm[i]; perm[i] = perm[j]; perm[j] = t;
    }
    for (int i = 0; i < MAX_INF; i++) shuf[i] = hand_case[perm[i]];
    GreedyStep s2[MAX_INF]; V2 full2; float R2;
    int c2 = greedy_select(shuf, MAX_INF, MAX_INF, s2, &full2, &R2);
    int order_ok = (c2 == count && vdist(full, full2) <= EPS && R == R2);
    for (int k = 0; k < count && order_ok; k++) {
        V2 a = steps[k].blend, b = s2[k].blend;
        if (vdist(a, b) > EPS ||
            fabsf(steps[k].bound - s2[k].bound) > EPS ||
            fabsf(steps[k].discarded_frac - s2[k].discarded_frac) > EPS)
            order_ok = 0;
    }
    check("greedy_order_independent", order_ok);

    /* ---- tie breaking: equal contribution, heavier raw weight wins ---- */
    Influence tie[4] = {
        { { 1.0f, 0.0f}, 2.0f },   /* contrib 2.0 */
        { { 0.0f, 2.0f}, 1.0f },   /* contrib 2.0 */
        { { 0.5f, 0.0f}, 1.0f },   /* contrib 0.5 */
        { { 0.0f, 0.0f}, 0.5f },   /* contrib 0.0 */
    };
    GreedyStep t1[4], t2s[4]; V2 f1, f2; float rr1, rr2;
    greedy_select(tie, 4, 4, t1, &f1, &rr1);
    Influence tie_swap[4] = { tie[1], tie[0], tie[2], tie[3] };
    greedy_select(tie_swap, 4, 4, t2s, &f2, &rr2);
    check("tie_prefers_heavier_weight", t1[0].pick.idx == 0);
    check("tie_break_order_independent",
          vdist(t1[0].blend, t2s[0].blend) <= EPS);

    /* ---- randomized stress: random weights/positions, incl. zero weights ---- */
    int stress_ok = 1, trials = 500;
    for (int tr = 0; tr < trials; tr++) {
        Influence s[MAX_INF];
        int n = 2 + (int)(randf() * (MAX_INF - 1));   /* 2..8 */
        for (int i = 0; i < n; i++) {
            s[i].pos.x = randf() * 4.0f - 2.0f;
            s[i].pos.y = randf() * 4.0f - 2.0f;
            s[i].w = (i % 5 == 0) ? 0.0f : randf() * 5.0f;
        }
        int max_inf = 1 + (int)(randf() * n);
        GreedyStep gs[MAX_INF]; V2 gf; float gR;
        int gc = greedy_select(s, n, max_inf, gs, &gf, &gR);
        if (gc <= 0) { stress_ok = 0; break; }

        float s_total = total_weight(s, n);
        float prev_bound = 1e30f;
        for (int k = 0; k < gc; k++) {
            V2 opt = exhaustive_best(s, n, k + 1, gf);
            if (vdist(gs[k].blend, opt) > gs[k].bound + EPS) { stress_ok = 0; }
            if (gs[k].bound > prev_bound + EPS) stress_ok = 0;
            prev_bound = gs[k].bound;
            if (!renorm_sums_one(s, gs, k, s_total)) stress_ok = 0;
            if (vdist(gs[k].blend, gf) >
                2.0f * gR * gs[k].discarded_frac + EPS) stress_ok = 0;
        }
    }
    check("stress_500_random_vertices", stress_ok);

    /* ---- stress on rejection: random negative injections ---- */
    int reject_ok = 1;
    for (int tr = 0; tr < 100; tr++) {
        Influence s[MAX_INF];
        for (int i = 0; i < MAX_INF; i++) {
            s[i].pos.x = randf(); s[i].pos.y = randf();
            s[i].w = randf();
        }
        s[(int)(randf() * MAX_INF)].w = -randf() - 0.01f;
        if (greedy_select(s, MAX_INF, 4, steps, &full, &R) != -1) reject_ok = 0;
    }
    check("stress_negative_weights_rejected", reject_ok);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
