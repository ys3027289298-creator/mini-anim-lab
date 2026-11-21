
/* constraint: verlet particle chain with distance constraints, iteration
   count and damping; used as the toy model for cloth/rope style solvers.

   Also contains an alternating-projection solver that groups constraints
   (distance links) and collisions (circles) separately, rotates the group
   order every round and watches the residual decay.  A unified solve is used
   as the reference, and a non-decaying residual is reported together with
   the groups that are stuck. */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define MAX_P 16

typedef struct { float x, y, px, py; float inv_mass; int pinned; } Particle;

typedef struct {
    Particle p[MAX_P];
    int n;
    float rest[MAX_P];      /* rest length between i and i+1 */
    int iterations;
    float damping;
} Chain;

static void chain_init(Chain* c, int n, float seg, float x0, float y0) {
    memset(c, 0, sizeof(*c));
    c->n = n;
    c->iterations = 8;
    c->damping = 0.99f;
    for (int i = 0; i < n; i++) {
        c->p[i].x = x0 + seg * i; c->p[i].y = y0;
        c->p[i].px = c->p[i].x; c->p[i].py = c->p[i].y;
        c->p[i].inv_mass = 1.0f;
        c->p[i].pinned = 0;
        c->rest[i] = seg;
    }
    c->p[0].pinned = 1;
    c->p[0].inv_mass = 0.0f;
}

static void chain_step(Chain* c, float dt, float gx, float gy) {
    for (int i = 0; i < c->n; i++) {
        Particle* p = &c->p[i];
        if (p->pinned) continue;
        float vx = (p->x - p->px) * c->damping;
        float vy = (p->y - p->py) * c->damping;
        p->px = p->x; p->py = p->y;
        p->x += vx + gx * dt * dt;
        p->y += vy + gy * dt * dt;
    }
    for (int it = 0; it < c->iterations; it++) {
        for (int i = 0; i < c->n - 1; i++) {
            Particle* a = &c->p[i];
            Particle* b = &c->p[i + 1];
            float dx = b->x - a->x, dy = b->y - a->y;
            float d = sqrtf(dx * dx + dy * dy);
            if (d < 1e-6f) continue;
            float diff = (d - c->rest[i]) / d;
            float wa = a->inv_mass, wb = b->inv_mass;
            float wsum = wa + wb;
            if (wsum <= 0) continue;
            a->x += dx * diff * (wa / wsum);
            a->y += dy * diff * (wa / wsum);
            b->x -= dx * diff * (wb / wsum);
            b->y -= dy * diff * (wb / wsum);
        }
    }
}

static float chain_seg_error(const Chain* c) {
    float worst = 0.0f;
    for (int i = 0; i < c->n - 1; i++) {
        float dx = c->p[i + 1].x - c->p[i].x;
        float dy = c->p[i + 1].y - c->p[i].y;
        float d = sqrtf(dx * dx + dy * dy);
        float e = fabsf(d - c->rest[i]);
        if (e > worst) worst = e;
    }
    return worst;
}

/* ---- alternating projection: constraint group vs collision group --------
   Constraints (distance links) and collisions (circles) are projected as
   separate groups.  The alternating solver rotates the group order every
   round, records the residual of both groups, and decides convergence from
   the residual decay instead of only an absolute tolerance.
   A unified solver (one flat sweep over links and collisions, fixed order)
   provides the reference result the alternating one is compared against. */

#define MAX_ITERS 2048
#define MAX_CIRCLE 8
#define DECAY_EPS 1e-4f      /* relative residual drop over a window = "decaying" */
#define STALL_WINDOW 6       /* rounds compared; even, so both group orders
                                appear equally often on each side */

#define STUCK_NONE  0
#define STUCK_LINKS 1
#define STUCK_COLS  2

/* extended constraint: a distance link with a target distance that differs
   from rest, used to model the "constraint group" of a scene */
typedef struct { int a, b; float rest; } Link;

#define MAX_LINK 2 * MAX_P

typedef struct { float cx, cy, r; } Circle;

typedef struct {
    int n;
    Link links[MAX_LINK];
    int n_links;
    int pinned[MAX_P];
    Circle circ[MAX_CIRCLE];/* circle obstacles, up to MAX_P particles collide */
    int n_circ;
} Scene;

typedef struct {
    int rounds;
    float res;              /* combined residual after the last round */
    float res_links;
    float res_cols;
    int converged;          /* res <= tol */
    int stalled;            /* residual stopped decaying */
    int stall_round;        /* 1-based round the stall was detected on */
    int stuck;              /* STUCK_* bitmask of groups that stopped decaying */
    int rotated;
    float hist[MAX_ITERS];  /* combined residual per round */
    float hist_links[MAX_ITERS];
    float hist_cols[MAX_ITERS];
    int hist_cols_first[MAX_ITERS];
} SolveReport;

typedef enum { SOLVE_ALTERNATING, SOLVE_UNIFIED } SolveMode;

static void scene_init(Scene* s, int n) {
    memset(s, 0, sizeof(*s));
    s->n = n;
    for (int i = 0; i < n; i++) s->pinned[i] = 0;
    s->pinned[0] = 1;
}

static void scene_add_link(Scene* s, int a, int b, float rest) {
    if (s->n_links >= MAX_LINK || a < 0 || b < 0 || a >= s->n || b >= s->n || a == b)
        return;
    Link* l = &s->links[s->n_links++];
    l->a = a; l->b = b; l->rest = rest;
}

static void scene_add_circle(Scene* s, float cx, float cy, float r) {
    if (s->n_circ >= MAX_CIRCLE) return;
    s->circ[s->n_circ].cx = cx;
    s->circ[s->n_circ].cy = cy;
    s->circ[s->n_circ].r = r;
    s->n_circ++;
}

static void project_link(const Scene* s, int k, float* pos) {
    const Link* L = &s->links[k];
    int a = L->a, b = L->b;
    float wa = s->pinned[a] ? 0.0f : 1.0f;
    float wb = s->pinned[b] ? 0.0f : 1.0f;
    float wsum = wa + wb;
    if (wsum <= 0.0f) return;
    float dx = pos[2 * b] - pos[2 * a];
    float dy = pos[2 * b + 1] - pos[2 * a + 1];
    float d = sqrtf(dx * dx + dy * dy);
    if (d < 1e-6f) return;
    float diff = (d - L->rest) / d;
    pos[2 * a]     += dx * diff * (wa / wsum);
    pos[2 * a + 1] += dy * diff * (wa / wsum);
    pos[2 * b]     -= dx * diff * (wb / wsum);
    pos[2 * b + 1] -= dy * diff * (wb / wsum);
}

static void project_particle_collision(const Scene* s, int i, float* pos) {
    if (s->pinned[i]) return;
    for (int k = 0; k < s->n_circ; k++) {
        Circle c = s->circ[k];
        float dx = pos[2 * i] - c.cx;
        float dy = pos[2 * i + 1] - c.cy;
        float d = sqrtf(dx * dx + dy * dy);
        if (d < c.r) {
            if (d < 1e-6f) { dx = 1.0f; dy = 0.0f; d = 1.0f; }
            pos[2 * i]     = c.cx + dx / d * c.r;
            pos[2 * i + 1] = c.cy + dy / d * c.r;
        }
    }
}

static void project_links(const Scene* s, float* pos) {
    for (int k = 0; k < s->n_links; k++) project_link(s, k, pos);
}

static void project_cols(const Scene* s, float* pos) {
    for (int i = 0; i < s->n; i++) project_particle_collision(s, i, pos);
}

/* unified sweep: links and collisions in one flat, fixed order */
static void project_unified(const Scene* s, float* pos) {
    for (int k = 0; k < s->n_links; k++) {
        if (k < s->n) project_particle_collision(s, k, pos);
        project_link(s, k, pos);
    }
    for (int i = s->n_links; i < s->n; i++) project_particle_collision(s, i, pos);
}

static float scene_res_links(const Scene* s, const float* pos) {
    float worst = 0.0f;
    for (int k = 0; k < s->n_links; k++) {
        const Link* L = &s->links[k];
        float dx = pos[2 * L->b] - pos[2 * L->a];
        float dy = pos[2 * L->b + 1] - pos[2 * L->a + 1];
        float e = fabsf(sqrtf(dx * dx + dy * dy) - L->rest);
        if (e > worst) worst = e;
    }
    return worst;
}

static float scene_res_cols(const Scene* s, const float* pos) {
    float worst = 0.0f;
    for (int i = 0; i < s->n; i++) {
        if (s->pinned[i]) continue;
        for (int k = 0; k < s->n_circ; k++) {
            Circle c = s->circ[k];
            float dx = pos[2 * i] - c.cx;
            float dy = pos[2 * i + 1] - c.cy;
            float e = c.r - sqrtf(dx * dx + dy * dy);
            if (e > worst) worst = e;
        }
    }
    return worst;
}

/* A group is stuck when its residual is still above tol but has stopped
   decaying.  An alternating sweep leaves a different group clean every round
   (whichever ran last), so a naive adjacent-round comparison reads parity
   oscillation as progress.  Instead, sample each group only on the rounds
   where it ran first -- then its residual is the violation left over after
   the other group moved the particles, which is the quantity that stalling
   is supposed to detect. */
static int last_round_with_phase(const SolveReport* r, int cols_first) {
    /* a fixed order (or the unified sweep) has a single phase for both groups */
    if (!r->rotated) return r->rounds - 1;
    for (int i = r->rounds - 1; i >= 0; i--)
        if (r->hist_cols_first[i] == cols_first) return i;
    return -1;
}

static int group_stuck(const SolveReport* r, const float* hist, int cols_first,
                       float tol) {
    int i1 = last_round_with_phase(r, cols_first);
    int i0 = i1 - STALL_WINDOW;
    if (i0 < 0 || i1 < 0) return 0;
    float prev = hist[i0], now = hist[i1];
    if (now <= tol) return 0;
    float rel = prev > 1e-9f ? (prev - now) / prev : 0.0f;
    return rel < DECAY_EPS;
}

static int stuck_groups(const SolveReport* r, float tol) {
    int stuck = STUCK_NONE;
    if (group_stuck(r, r->hist_links, 0, tol)) stuck |= STUCK_LINKS;
    if (group_stuck(r, r->hist_cols, 1, tol)) stuck |= STUCK_COLS;
    return stuck;
}

static SolveReport solve_scene(const Scene* s, float* pos, SolveMode mode,
                               int order, int rotate, float tol,
                               int max_rounds, int log) {
    SolveReport r;
    memset(&r, 0, sizeof(r));
    r.rotated = rotate;
    if (max_rounds > MAX_ITERS) max_rounds = MAX_ITERS;

    /* measure convergence against the same group order every STALL_WINDOW
       rounds (STALL_WINDOW is even), which filters the parity oscillation */
    for (int round = 0; round < max_rounds; round++) {
        int cols_first = 0;
        if (mode == SOLVE_UNIFIED) {
            project_unified(s, pos);
        } else {
            cols_first = rotate ? ((round + order) & 1) : (order & 1);
            if (cols_first) { project_cols(s, pos); project_links(s, pos); }
            else            { project_links(s, pos); project_cols(s, pos); }
        }

        r.res_links = scene_res_links(s, pos);
        r.res_cols = scene_res_cols(s, pos);
        r.res = fmaxf(r.res_links, r.res_cols);
        r.hist[round] = r.res;
        r.hist_links[round] = r.res_links;
        r.hist_cols[round] = r.res_cols;
        r.hist_cols_first[round] = cols_first;
        r.rounds = round + 1;
        if (log)
            printf("round=%d %s res=%.6f links=%.6f cols=%.6f\n",
                   r.rounds,
                   mode == SOLVE_UNIFIED ? "unified" : (cols_first ? "cols_first" : "links_first"),
                   r.res, r.res_links, r.res_cols);

        if (tol > 0.0f && r.res <= tol) { r.converged = 1; break; }

        /* With a rotating order each group runs first every other round, so a
           group needs STALL_WINDOW/2 * 2 + 1 rounds to give two same-phase
           samples STALL_WINDOW apart.  Wait for that before judging. */
        if (tol > 0.0f && round > STALL_WINDOW) {
            float prev = r.hist[round - STALL_WINDOW];
            float rel = prev > 1e-9f ? (prev - r.res) / prev : 0.0f;
            int same_phase = rotate == 0 ||
                             r.hist_cols_first[round] == r.hist_cols_first[round - STALL_WINDOW];
            if (rel < DECAY_EPS && same_phase) {
                r.stalled = 1;
                r.stall_round = r.rounds;
                r.stuck = stuck_groups(&r, tol);
                break;
            }
        }
    }
    return r;
}

static float pos_diff(int n, const float* a, const float* b) {
    float worst = 0.0f;
    for (int i = 0; i < 2 * n; i++) {
        float d = fabsf(a[i] - b[i]);
        if (d > worst) worst = d;
    }
    return worst;
}

static const char* stuck_name(int stuck) {
    if (stuck == STUCK_LINKS) return "constraints";
    if (stuck == STUCK_COLS) return "collisions";
    if ((stuck & (STUCK_LINKS | STUCK_COLS)) == (STUCK_LINKS | STUCK_COLS))
        return "constraints+collisions";
    return "none";
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    Chain c;
    chain_init(&c, 8, 1.0f, 0.0f, 0.0f);
    check("init_len", c.n == 8 && fabsf(chain_seg_error(&c)) < 1e-5f);

    /* let it fall for a while */
    for (int i = 0; i < 200; i++) chain_step(&c, 0.016f, 0.0f, -9.8f);
    check("stays_in_stretch", chain_seg_error(&c) < 0.02f);
    check("pinned_anchor", fabsf(c.p[0].x) < 1e-6f && fabsf(c.p[0].y) < 1e-6f);
    check("falls_down", c.p[7].y < -1.0f);

    /* more iterations must not increase error */
    Chain a, b;
    chain_init(&a, 8, 1.0f, 0.0f, 0.0f);
    chain_init(&b, 8, 1.0f, 0.0f, 0.0f);
    a.iterations = 1; b.iterations = 16;
    for (int i = 0; i < 100; i++) { chain_step(&a, 0.016f, 0.0f, -9.8f); chain_step(&b, 0.016f, 0.0f, -9.8f); }
    check("more_iters_better", chain_seg_error(&b) <= chain_seg_error(&a) + 1e-6f);

    Chain cap;
    chain_init(&cap, 16, 0.5f, 0.0f, 0.0f);
    check("capacity_ok", cap.n == 16);

    /* ---- alternating projection: constraints vs collisions ------------- */
    {
        /* 5-particle chain pinned at both ends, draped over a circle that sits
           high between the anchors: links want the chain to hang away from the
           circle, collisions push it up over it.  The bounds keep the solution
           unique, so both solvers can be compared point by point. */
        Scene s;
        scene_init(&s, 5);
        scene_add_link(&s, 0, 1, 1.086f);
        scene_add_link(&s, 1, 2, 0.727f);
        scene_add_link(&s, 2, 3, 1.030f);
        scene_add_link(&s, 3, 4, 0.605f);
        scene_add_circle(&s, -0.131f, 0.492f, 0.741f);
        s.pinned[0] = 1;
        s.pinned[4] = 1;
        float start[2 * MAX_P] = {
            -1.5f, 0.0f,  0.484f, -0.277f,  1.835f, 0.335f,  2.527f, 0.598f,  1.5f, 0.0f
        };

        float alt[2 * MAX_P], uni[2 * MAX_P], alt2[2 * MAX_P];
        memcpy(alt, start, sizeof(alt));
        memcpy(uni, start, sizeof(uni));
        memcpy(alt2, start, sizeof(alt2));

        printf("\n-- alternating (links first, rotating) --\n");
        SolveReport ra = solve_scene(&s, alt, SOLVE_ALTERNATING, 0, 1, 1e-5f, 512, 0);
        printf("\n-- unified reference --\n");
        SolveReport ru = solve_scene(&s, uni, SOLVE_UNIFIED, 0, 0, 1e-5f, 512, 0);

        check("alt_converges", ra.converged == 1);
        check("uni_converges", ru.converged == 1);
        /* the alternating result is compared against the unified reference,
           and its residual must sit under a hard upper bound */
        float vd = pos_diff(s.n, alt, uni);
        float ub = 1e-4f;
        printf("alt_vs_unified_pos=%.7f residual_upper_bound=%.7f rounds_alt=%d rounds_uni=%d\n",
               vd, ub, ra.rounds, ru.rounds);
        check("alt_matches_unified", vd < 1e-3f);
        check("alt_residual_bounded", ra.res <= ub && ru.res <= ub);
        check("alt_no_false_stall", ra.stalled == 0);

        /* residual decays: strictly smaller than round one, by a wide margin */
        check("alt_residual_decays",
              ra.rounds >= 2 && ra.hist[ra.rounds - 1] < 0.05f * ra.hist[0]);

        /* group order changes the convergence speed: run the four orderings
           (rotating/fixed x links-first/cols-first) from the same start and
           report the spread in rounds */
        printf("\n-- alternating (cols first, rotating) --\n");
        SolveReport rb = solve_scene(&s, alt2, SOLVE_ALTERNATING, 1, 1, 1e-5f, 512, 0);
        printf("order_links_first_rounds=%d order_cols_first_rounds=%d gap=%d\n",
               ra.rounds, rb.rounds, ra.rounds - rb.rounds);
        check("alt_order_both_converge", ra.converged == 1 && rb.converged == 1);
        check("alt_order_same_solution", pos_diff(s.n, alt, alt2) < 1e-3f);

        float fixed_l[2 * MAX_P], fixed_c[2 * MAX_P];
        memcpy(fixed_l, start, sizeof(fixed_l));
        memcpy(fixed_c, start, sizeof(fixed_c));
        SolveReport rf  = solve_scene(&s, fixed_l, SOLVE_ALTERNATING, 0, 0, 1e-5f, 512, 0);
        SolveReport rfc = solve_scene(&s, fixed_c, SOLVE_ALTERNATING, 1, 0, 1e-5f, 512, 0);
        int cheapest = ra.rounds, priciest = ra.rounds;
        int all[4] = { ra.rounds, rb.rounds, rf.rounds, rfc.rounds };
        for (int i = 0; i < 4; i++) {
            if (all[i] < cheapest) cheapest = all[i];
            if (all[i] > priciest) priciest = all[i];
        }
        printf("order_table rounds: rotating_links_first=%d rotating_cols_first=%d "
               "fixed_links_first=%d fixed_cols_first=%d | fastest=%d slowest=%d spread=%d\n",
               ra.rounds, rb.rounds, rf.rounds, rfc.rounds, cheapest, priciest, priciest - cheapest);
        check("alt_order_speed_differs", priciest - cheapest >= 2);
        check("alt_fixed_order_converges", rf.converged == 1 && rfc.converged == 1);

        /* stuck case: the pinned anchor sits at the centre of a circle that
           swallows the loop, so the pinned particle can never leave it while
           the links pull it away -- the residual hits a plateau and the solver
           reports which groups stopped decaying */
        Scene st;
        scene_init(&st, 4);
        for (int i = 0; i < 4; i++) scene_add_link(&st, i, (i + 1) % 4, 1.0f);
        scene_add_circle(&st, 0.0f, 0.0f, 2.0f);
        float bad[2 * MAX_P];
        for (int i = 0; i < st.n; i++) { bad[2 * i] = 0.5f * i; bad[2 * i + 1] = 0.0f; }
        printf("\n-- stuck (pinned anchor inside circle) --\n");
        SolveReport rs = solve_scene(&st, bad, SOLVE_ALTERNATING, 0, 1, 1e-6f, 64, 0);
        printf("stalled=%d stall_round=%d stuck=%s res=%.6f\n",
               rs.stalled, rs.stall_round, stuck_name(rs.stuck), rs.res);
        check("stall_detected", rs.stalled == 1 && rs.converged == 0);

        /* every particle is trapped inside the circle and the links are far
           shorter than the circle: neither group can give ground, so both are
           reported as stuck on the same plateau */
        printf("stuck_mask=%d names=%s\n", rs.stuck, stuck_name(rs.stuck));
        check("stall_names_both_groups", rs.stuck == (STUCK_LINKS | STUCK_COLS));
        check("stall_residual_above_tol", rs.res > 1e-6f);

        /* single-group stall: the collision group is satisfiable but a link is
           held apart by a pinned pair, so only the constraints can be stuck.
           A two-particle link pinned at both ends to the wrong distance is
           infeasible for the links and has no collision at all. */
        Scene sl;
        scene_init(&sl, 2);
        scene_add_link(&sl, 0, 1, 1.0f);
        sl.pinned[1] = 1;
        float split[2 * MAX_P] = { 3.0f, 0.0f, 0.0f, 0.0f };
        printf("\n-- stuck (pinned pair held too far apart) --\n");
        SolveReport rsl = solve_scene(&sl, split, SOLVE_ALTERNATING, 0, 1, 1e-6f, 64, 0);
        printf("stalled=%d stuck=%s res=%.6f\n", rsl.stalled, stuck_name(rsl.stuck), rsl.res);
        check("stall_links_only", rsl.stalled == 1 && rsl.stuck == STUCK_LINKS);
    }

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
