
/* topo: dependency-ordered constraint solver on a general particle graph.
   Constraints are sorted by topological level (BFS depth from pinned
   roots) so corrections propagate root->leaf in a single pass; constraints
   closing dependency cycles are reported and iterated separately. */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define TP_MAX_P 32
#define TP_MAX_C 64

typedef struct { float x, y; float inv_mass; int pinned; int depth; } TPParticle;
typedef struct { int a, b; float rest; int level; int on_cycle; } TPConstraint;

typedef struct {
    TPParticle p[TP_MAX_P];
    int n;
    TPConstraint c[TP_MAX_C];
    int m;
    int order[TP_MAX_C]; /* constraint indices in topological solve order */
    int depth;           /* topological depth: max constraint level + 1 */
    int n_cycle;         /* number of constraints lying on a cycle */
} TPGraph;

static void tp_init(TPGraph* g) { memset(g, 0, sizeof(*g)); }

static int tp_add_particle(TPGraph* g, float x, float y, int pinned) {
    int i = g->n++;
    g->p[i].x = x; g->p[i].y = y;
    g->p[i].pinned = pinned;
    g->p[i].inv_mass = pinned ? 0.0f : 1.0f;
    g->p[i].depth = -1;
    return i;
}

static int tp_add_constraint(TPGraph* g, int a, int b, float rest) {
    int i = g->m++;
    g->c[i].a = a; g->c[i].b = b; g->c[i].rest = rest;
    g->c[i].level = 0; g->c[i].on_cycle = 0;
    return i;
}

static int uf_find(int* parent, int x) {
    while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; }
    return x;
}

/* build solve order, topological depth, and cycle report */
static void tp_analyze(TPGraph* g) {
    int uf[TP_MAX_P], head[TP_MAX_P];
    int fa_to[TP_MAX_C * 2], fa_ci[TP_MAX_C * 2], fa_next[TP_MAX_C * 2];
    int fe = 0;
    for (int i = 0; i < g->n; i++) { uf[i] = i; head[i] = -1; }

    /* union-find over constraints; a constraint whose endpoints are already
       connected closes a cycle: mark it and every forest edge on the path */
    g->n_cycle = 0;
    for (int ci = 0; ci < g->m; ci++) {
        int a = g->c[ci].a, b = g->c[ci].b;
        int ra = uf_find(uf, a), rb = uf_find(uf, b);
        if (ra != rb) {
            uf[ra] = rb;
            fa_to[fe] = b; fa_ci[fe] = ci; fa_next[fe] = head[a]; head[a] = fe++;
            fa_to[fe] = a; fa_ci[fe] = ci; fa_next[fe] = head[b]; head[b] = fe++;
        } else {
            int prev[TP_MAX_P], via[TP_MAX_P], q[TP_MAX_P];
            int qh = 0, qt = 0;
            g->c[ci].on_cycle = 1;
            for (int i = 0; i < g->n; i++) prev[i] = -1;
            q[qt++] = a; prev[a] = a;
            while (qh < qt) {
                int u = q[qh++];
                if (u == b) break;
                for (int e = head[u]; e != -1; e = fa_next[e]) {
                    int v = fa_to[e];
                    if (prev[v] == -1) { prev[v] = u; via[v] = fa_ci[e]; q[qt++] = v; }
                }
            }
            for (int v = b; v != a; v = prev[v]) g->c[via[v]].on_cycle = 1;
        }
    }
    for (int ci = 0; ci < g->m; ci++) g->n_cycle += g->c[ci].on_cycle;

    /* particle adjacency over all constraints */
    int ah[TP_MAX_P], at[TP_MAX_C * 2], an[TP_MAX_C * 2];
    int ae = 0;
    for (int i = 0; i < g->n; i++) ah[i] = -1;
    for (int ci = 0; ci < g->m; ci++) {
        int a = g->c[ci].a, b = g->c[ci].b;
        at[ae] = b; an[ae] = ah[a]; ah[a] = ae++;
        at[ae] = a; an[ae] = ah[b]; ah[b] = ae++;
    }

    /* BFS depth from pinned roots; unreached components root at their
       lowest-index particle so the solver stays total */
    int q[TP_MAX_P], qh = 0, qt = 0;
    for (int i = 0; i < g->n; i++) g->p[i].depth = -1;
    for (int i = 0; i < g->n; i++)
        if (g->p[i].pinned) { g->p[i].depth = 0; q[qt++] = i; }
    for (int s = 0; s < g->n; s++) {
        if (qh == qt && g->p[s].depth == -1) { g->p[s].depth = 0; q[qt++] = s; }
        while (qh < qt) {
            int u = q[qh++];
            for (int e = ah[u]; e != -1; e = an[e]) {
                int v = at[e];
                if (g->p[v].depth == -1) { g->p[v].depth = g->p[u].depth + 1; q[qt++] = v; }
            }
        }
    }

    /* constraint level = shallower endpoint depth; order by (level, index) */
    g->depth = 0;
    for (int ci = 0; ci < g->m; ci++) {
        int da = g->p[g->c[ci].a].depth, db = g->p[g->c[ci].b].depth;
        g->c[ci].level = da < db ? da : db;
        if (g->c[ci].level + 1 > g->depth) g->depth = g->c[ci].level + 1;
        g->order[ci] = ci;
    }
    for (int i = 1; i < g->m; i++) {
        int k = g->order[i], j = i - 1;
        while (j >= 0 && (g->c[g->order[j]].level > g->c[k].level ||
                          (g->c[g->order[j]].level == g->c[k].level && g->order[j] > k))) {
            g->order[j + 1] = g->order[j]; j--;
        }
        g->order[j + 1] = k;
    }
}

/* indices of constraints lying on cycles; returns count */
static int tp_cycle_constraints(const TPGraph* g, int* out) {
    int k = 0;
    for (int i = 0; i < g->m; i++)
        if (g->c[i].on_cycle) out[k++] = i;
    return k;
}

static void tp_project_w(TPGraph* g, TPConstraint* c, float wa, float wb) {
    TPParticle* a = &g->p[c->a];
    TPParticle* b = &g->p[c->b];
    float dx = b->x - a->x, dy = b->y - a->y;
    float d = sqrtf(dx * dx + dy * dy);
    if (d < 1e-6f) return;
    float diff = (d - c->rest) / d;
    float ws = wa + wb;
    if (ws <= 0) return;
    a->x += dx * diff * (wa / ws);
    a->y += dy * diff * (wa / ws);
    b->x -= dx * diff * (wb / ws);
    b->y -= dy * diff * (wb / ws);
}

static void tp_project(TPGraph* g, TPConstraint* c) {
    tp_project_w(g, c, g->p[c->a].inv_mass, g->p[c->b].inv_mass);
}

/* one iteration = one ordered pass, then cycle_iterations extra symmetric
   passes over the cyclic constraints only. In the ordered pass a
   constraint whose endpoints sit at different depths is projected
   toward the deeper particle only (follow the leader), so corrections
   propagate root->leaf within a single pass; same-depth (cycle)
   constraints use the symmetric split. */
static void tp_solve(TPGraph* g, const int* order, int iterations, int cycle_iterations) {
    for (int it = 0; it < iterations; it++) {
        for (int k = 0; k < g->m; k++) {
            TPConstraint* c = &g->c[order[k]];
            int da = g->p[c->a].depth, db = g->p[c->b].depth;
            if (da < db) tp_project_w(g, c, 0.0f, g->p[c->b].inv_mass);
            else if (db < da) tp_project_w(g, c, g->p[c->a].inv_mass, 0.0f);
            else tp_project(g, c);
        }
        for (int r = 0; r < cycle_iterations; r++)
            for (int k = 0; k < g->m; k++)
                if (g->c[k].on_cycle) tp_project(g, &g->c[k]);
    }
}

static float tp_residual(const TPGraph* g) {
    float worst = 0.0f;
    for (int ci = 0; ci < g->m; ci++) {
        const TPConstraint* c = &g->c[ci];
        float dx = g->p[c->b].x - g->p[c->a].x;
        float dy = g->p[c->b].y - g->p[c->a].y;
        float e = fabsf(sqrtf(dx * dx + dy * dy) - c->rest);
        if (e > worst) worst = e;
    }
    return worst;
}

static float tp_residual_sum(const TPGraph* g) {
    float sum = 0.0f;
    for (int ci = 0; ci < g->m; ci++) {
        const TPConstraint* c = &g->c[ci];
        float dx = g->p[c->b].x - g->p[c->a].x;
        float dy = g->p[c->b].y - g->p[c->a].y;
        sum += fabsf(sqrtf(dx * dx + dy * dy) - c->rest);
    }
    return sum;
}

static unsigned rng_state;
static unsigned rng_next(void) {
    rng_state ^= rng_state << 13; rng_state ^= rng_state >> 17; rng_state ^= rng_state << 5;
    return rng_state;
}
static void shuffled_order(const TPGraph* g, int* order, unsigned seed) {
    rng_state = seed ? seed : 1;
    for (int i = 0; i < g->m; i++) order[i] = i;
    for (int i = g->m - 1; i > 0; i--) {
        int j = (int)(rng_next() % (unsigned)(i + 1));
        int t = order[i]; order[i] = order[j]; order[j] = t;
    }
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

/* stretched-out chain, pinned at particle 0 */
static void build_chain(TPGraph* g, int n, float seg, float stretch) {
    tp_init(g);
    for (int i = 0; i < n; i++)
        tp_add_particle(g, i * seg * stretch, -0.3f * i, i == 0);
    for (int i = 0; i < n - 1; i++)
        tp_add_constraint(g, i, i + 1, seg);
    tp_analyze(g);
}

int main(void) {
    /* 1. topological order solves a chain in a single pass */
    TPGraph g;
    build_chain(&g, 12, 1.0f, 1.6f);
    printf("chain_topo_depth=%d\n", g.depth);
    check("chain_depth", g.depth == 11);
    check("chain_acyclic", g.n_cycle == 0);
    TPGraph t = g;
    tp_solve(&t, t.order, 1, 0);
    float topo_e1 = tp_residual(&t);
    printf("topo_residual_1pass=%.6f\n", topo_e1);
    check("topo_one_pass", topo_e1 < 1e-3f);

    /* 2. random order at the same iteration count: worse and unstable */
    float rmin = 1e9f, rmax = 0.0f;
    for (int s = 0; s < 8; s++) {
        TPGraph r = g;
        int ord[TP_MAX_C];
        shuffled_order(&r, ord, 1234 + s * 77);
        tp_solve(&r, ord, 1, 0);
        float e = tp_residual(&r);
        if (e < rmin) rmin = e;
        if (e > rmax) rmax = e;
    }
    printf("random_residual_min=%.6f max=%.6f\n", rmin, rmax);
    check("random_worse_same_iters", rmin > 0.01f && topo_e1 < rmin);
    check("random_unstable", rmax - rmin > 0.01f);
    TPGraph t1 = g, t2 = g;
    tp_solve(&t1, t1.order, 2, 0);
    tp_solve(&t2, t2.order, 2, 0);
    check("topo_deterministic", tp_residual(&t1) == tp_residual(&t2));

    /* 3. residual decreases monotonically with iterations (on a ring,
       where one pass is not enough) */
    TPGraph ring;
    tp_init(&ring);
    for (int i = 0; i < 6; i++) {
        float ang = 6.2831853f * (float)i / 6.0f;
        float rad = 2.0f + 0.3f * (float)(i % 2);
        tp_add_particle(&ring, cosf(ang) * rad, sinf(ang) * rad, i == 0);
    }
    for (int i = 0; i < 6; i++)
        tp_add_constraint(&ring, i, (i + 1) % 6, 1.0f);
    tp_analyze(&ring);
    float prev = 1e9f;
    int mono = 1;
    for (int it = 1; it <= 6; it++) {
        TPGraph r = ring;
        tp_solve(&r, r.order, it, 2);
        float e = tp_residual(&r);
        if (e > prev + 1e-6f) mono = 0;
        prev = e;
    }
    check("topo_monotonic", mono);

    /* 4. deeper chains: topological order stays converged in one pass,
       random order degrades as depth grows */
    int depths[] = {4, 8, 16};
    float avg[3];
    int topo_ok = 1;
    for (int k = 0; k < 3; k++) {
        TPGraph g2;
        build_chain(&g2, depths[k], 1.0f, 1.6f);
        TPGraph tt = g2;
        tp_solve(&tt, tt.order, 1, 0);
        if (tp_residual(&tt) > 1e-3f) topo_ok = 0;
        float sum = 0.0f;
        for (int s = 0; s < 8; s++) {
            TPGraph r = g2;
            int ord[TP_MAX_C];
            shuffled_order(&r, ord, 999 + s * 31);
            tp_solve(&r, ord, 1, 0);
            sum += tp_residual_sum(&r);
        }
        avg[k] = sum / 8.0f;
    }
    printf("random_avg_by_depth=%.4f,%.4f,%.4f\n", avg[0], avg[1], avg[2]);
    check("topo_depth_invariant", topo_ok);
    check("random_degrades_with_depth", avg[1] > avg[0] && avg[2] > avg[1]);

    /* 5. cycle report and iterative solve of the cyclic part */
    printf("ring_topo_depth=%d\n", ring.depth);
    int cyc[TP_MAX_C];
    int nc = tp_cycle_constraints(&ring, cyc);
    printf("cycle_constraints=");
    for (int i = 0; i < nc; i++) printf("%s%d", i ? "," : "", cyc[i]);
    printf("\n");
    check("cycle_detected", ring.n_cycle == 6 && nc == 6);
    TPGraph r0 = ring, r1 = ring;
    tp_solve(&r0, r0.order, 8, 0);
    tp_solve(&r1, r1.order, 8, 4);
    printf("ring_residual_no_cycle_iter=%.6f with=%.6f\n", tp_residual(&r0), tp_residual(&r1));
    check("cycle_iters_help", tp_residual(&r1) < tp_residual(&r0));
    check("cycle_converges", tp_residual(&r1) < 0.01f);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
