/* cloth: spatial-block Verlet cloth with parallel block projection.

   The sheet is a grid of distance-constrained vertices.  Cells (the gaps
   between columns) are split into contiguous blocks; a seam column is one
   that lies on a block boundary and is therefore shared by the two blocks
   either side of it.

   Each block owns a contiguous range of columns.  A block projects every
   constraint that touches one of its owned vertices and keeps the result in
   its own correction buffer; a vertex is applied only by the block that owns
   it, reading only that block's buffer.  Seam constraints are therefore
   computed twice (once by each neighbour) from the same pre-iteration
   positions and summed into the two buffers identically - which is exactly
   the arbitration that keeps the seam closed.

   Each vertex sums its incident constraints in global constraint index order,
   so splitting the grid changes which thread computes a term, never the order
   of the floating point additions.  A blocked solve is bit-identical to the
   whole-cloth solve for every block count, and the seam mismatch is exactly
   zero rather than merely small.

   With seam sharing disabled (share_seams = 0) seam constraints are not
   projected at all, so each block solves its interior in isolation;
   cloth_mismatch() then reports the tear against the whole-cloth solve.

   make test runs the self-checks, the exact seam comparison, the tear
   detection and the block-count timing table. */
#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdatomic.h>
#include <sched.h>
#include <pthread.h>
#include <time.h>

#define MAX_W 96
#define MAX_H 96
#define MAX_V (MAX_W * MAX_H)
#define MAX_CONS ((MAX_W - 1) * MAX_H + MAX_W * (MAX_H - 1))
/* Seam constraints appear in two blocks' projection lists, so the lists
   together hold at most ncon + nseam <= 2 * MAX_CONS entries. */
#define MAX_BLOCK_CONS (2 * MAX_CONS)
#define MAX_BLOCKS 8
#define MAX_DEG 4          /* left/right/up/down */

typedef struct {
    float x, y;
    float px, py;
    float inv_mass;
    int pinned;
} CVertex;

typedef struct {
    int a, b;              /* vertex indices */
    float rest;
    int owner;             /* block the constraint belongs to, -1 if seam */
    int sa, sb;            /* slot in a's / b's incident-constraint list */
} CConstraint;

/* per-vertex correction slots, indexed by the constraint's slot in that
   vertex's incident-constraint list, and by vertex index.  Each block has its
   own buffer: a block writes every constraint it projects into its own
   buffer, and a vertex is applied only by its owner from that owner's
   buffer, so no two threads ever write the same cache line. */
typedef struct { float d[MAX_V][MAX_DEG][2]; } Buf;

typedef struct {
    int w, h;
    float spacing;
    CVertex v[MAX_V];

    CConstraint con[MAX_CONS];
    int ncon;
    int deg[MAX_V];                 /* constraints incident to each vertex */
    int vcon[MAX_V][MAX_DEG];       /* ... and their global indices */

    int iterations;
    float damping;
    float gravity_x, gravity_y;
    float dt;
    float relax;

    int nb;
    int cell0[MAX_BLOCKS], cell1[MAX_BLOCKS];   /* cell-column ranges */
    int oc0[MAX_BLOCKS], oc1[MAX_BLOCKS];       /* owned column ranges */
    int col_block[MAX_W];                       /* owner, -1 if on a seam */
    int seam_col[MAX_W];

    int block_cons[MAX_BLOCK_CONS];  /* constraints touching owned vertices */
    int block_cstart[MAX_BLOCKS + 1];

    Buf* block_corr;                /* one correction buffer per block */
} Cloth;

/* ------------------------------------------------------------------ setup */

static int idx(const Cloth* c, int x, int y) { return x * c->h + y; }

static int cloth_init(Cloth* c, int w, int h, float spacing) {
    if (w < 2 || h < 2 || w > MAX_W || h > MAX_H || spacing <= 0.0f) return -1;
    memset(c, 0, sizeof(*c));
    c->w = w; c->h = h;
    c->spacing = spacing;
    c->iterations = 8;
    c->damping = 0.99f;
    c->gravity_x = 0.0f;
    c->gravity_y = -9.8f;
    c->dt = 0.016f;
    c->relax = 1.0f;
    c->nb = 0;
    for (int x = 0; x < w; x++) {
        for (int y = 0; y < h; y++) {
            CVertex* p = &c->v[idx(c, x, y)];
            p->x = spacing * (float)x;
            p->y = spacing * (float)y;
            p->px = p->x;
            p->py = p->y;
            p->inv_mass = 1.0f;
            p->pinned = (y == h - 1);          /* top row hangs the sheet */
            if (p->pinned) p->inv_mass = 0.0f;
        }
    }
    return 0;
}

static void cloth_free(Cloth* c) {
    free(c->block_corr);
    c->block_corr = NULL;
}

/* add one distance constraint and register it in both endpoint lists, keeping
   each vertex's list in global constraint order */
static void add_con(Cloth* c, int a, int b, float rest, int owner) {
    CConstraint* k = &c->con[c->ncon];
    k->a = a; k->b = b;
    k->rest = rest;
    k->owner = owner;
    k->sa = c->deg[a];
    k->sb = c->deg[b];
    c->vcon[a][c->deg[a]++] = c->ncon;
    c->vcon[b][c->deg[b]++] = c->ncon;
    c->ncon++;
}

/* Cut the cells into nb contiguous ranges and rebuild the constraint tables.
   Ownership is contiguous: block k owns columns oc0[k]..oc1[k], so every
   column has exactly one owner and every vertex is applied once per
   iteration. */
static int cloth_build_blocks(Cloth* c, int nb) {
    int cells = c->w - 1;
    if (nb < 1 || nb > cells || nb > MAX_BLOCKS) return -1;
    c->nb = nb;

    int base = cells / nb, rem = cells % nb, cur = 0;
    for (int k = 0; k < nb; k++) {
        int n = base + (k < rem ? 1 : 0);
        c->cell0[k] = cur;
        c->cell1[k] = cur + n - 1;
        c->oc0[k] = cur;
        c->oc1[k] = cur + n - 1;
        cur += n;
    }
    c->oc1[nb - 1] = c->w - 1;

    memset(c->seam_col, 0, sizeof(c->seam_col));
    for (int k = 0; k + 1 < nb; k++) c->seam_col[c->cell1[k] + 1] = 1;

    for (int x = 0; x < c->w; x++) {
        c->col_block[x] = -1;
        if (c->seam_col[x]) continue;
        for (int k = 0; k < nb; k++) {
            if (x >= c->oc0[k] && x <= c->oc1[k]) { c->col_block[x] = k; break; }
        }
    }

    memset(c->deg, 0, sizeof(c->deg));
    c->ncon = 0;
    for (int x = 0; x < c->w - 1; x++) {
        int owner = (c->col_block[x] < 0 || c->col_block[x + 1] < 0)
                        ? -1 : c->col_block[x];
        for (int y = 0; y < c->h; y++)
            add_con(c, idx(c, x, y), idx(c, x + 1, y), c->spacing, owner);
    }
    for (int x = 0; x < c->w; x++) {
        int owner = c->col_block[x] < 0 ? -1 : c->col_block[x];
        for (int y = 0; y < c->h - 1; y++)
            add_con(c, idx(c, x, y), idx(c, x, y + 1), c->spacing, owner);
    }

    /* each block's projection list is every constraint with an endpoint in
       one of its owned columns, in global order */
    for (int k = 0; k < nb; k++) c->block_cstart[k] = 0;
    for (int k = 0; k < nb; k++) {
        int n = 0, x0 = c->oc0[k] * c->h, x1 = (c->oc1[k] + 1) * c->h;
        for (int i = 0; i < c->ncon; i++) {
            int a = c->con[i].a, b = c->con[i].b;
            if ((a >= x0 && a < x1) || (b >= x0 && b < x1))
                c->block_cons[c->block_cstart[k] + n++] = i;
        }
        c->block_cstart[k + 1] = c->block_cstart[k] + n;
    }
    return 0;
}

/* ---------------------------------------------------------------- buffers */

static int cloth_alloc_buffers(Cloth* c) {
    cloth_free(c);
    c->block_corr = calloc((size_t)c->nb, sizeof(Buf));
    return c->block_corr ? 0 : -1;
}

/* ----------------------------------------------------------------- solver */

static void integrate_vertex(CVertex* p, float dt, float gx, float gy, float damping) {
    if (p->pinned) return;
    float vx = (p->x - p->px) * damping;
    float vy = (p->y - p->py) * damping;
    p->px = p->x; p->py = p->y;
    p->x += vx + gx * dt * dt;
    p->y += vy + gy * dt * dt;
}

static void integrate_range(Cloth* c, int v0, int v1) {
    for (int v = v0; v < v1; v++)
        integrate_vertex(&c->v[v], c->dt, c->gravity_x, c->gravity_y, c->damping);
}

/* Compute one constraint's correction from the iteration-start positions.
   Seam constraints (owner < 0) are skipped when sharing is off, which is the
   tear scenario.  Both neighbours project a shared seam constraint from the
   same positions, so the two buffers receive identical contributions. */
static void project_con(Cloth* c, int ci, int active, Buf* corr) {
    const CConstraint* k = &c->con[ci];
    float* da = corr->d[k->a][k->sa];
    float* db = corr->d[k->b][k->sb];
    da[0] = da[1] = db[0] = db[1] = 0.0f;
    if (k->owner < 0 && !active) return;   /* seam dropped: tear case */

    const CVertex* a = &c->v[k->a];
    const CVertex* b = &c->v[k->b];
    float dx = b->x - a->x, dy = b->y - a->y;
    float dist = sqrtf(dx * dx + dy * dy);
    float wa = a->inv_mass, wb = b->inv_mass;
    float ws = wa + wb;
    if (dist <= 1e-6f || ws <= 0.0f) return;

    float diff = (dist - k->rest) / dist;
    float fa = diff * (wa / ws) * c->relax;
    float fb = diff * (wb / ws) * c->relax;
    da[0] = dx * fa;  da[1] = dy * fa;
    db[0] = -dx * fb; db[1] = -dy * fb;
}

/* Sum a vertex's incoming corrections in global constraint order. */
static void apply_vertex(Cloth* c, int v, Buf* corr) {
    CVertex* p = &c->v[v];
    if (p->pinned) return;
    float dx = 0.0f, dy = 0.0f;
    for (int j = 0; j < c->deg[v]; j++)
        { dx += corr->d[v][j][0]; dy += corr->d[v][j][1]; }
    p->x += dx;
    p->y += dy;
}

static void apply_range(Cloth* c, int v0, int v1, Buf* corr) {
    for (int v = v0; v < v1; v++) apply_vertex(c, v, corr);
}

/* whole-cloth reference: same arithmetic, one correction buffer, no partition */
static void cloth_solve_whole(Cloth* c, int steps) {
    static Buf corr;
    memset(&corr, 0, sizeof(corr));
    for (int s = 0; s < steps; s++) {
        integrate_range(c, 0, c->w * c->h);
        for (int it = 0; it < c->iterations; it++) {
            for (int i = 0; i < c->ncon; i++) project_con(c, i, 1, &corr);
            apply_range(c, 0, c->w * c->h, &corr);
        }
    }
}

/* Sense-reversing barrier.  A cloth does ~1k projections between barriers
   while pthread_barrier_wait costs a syscall per waiter; spinning keeps
   barrier cost far below the work it separates. */
typedef struct {
    _Atomic int count;
    _Atomic int gen;
    int n;
} Barrier;

static void barrier_init(Barrier* b, int n) {
    atomic_store(&b->count, 0);
    atomic_store(&b->gen, 0);
    b->n = n;
}

/* Spin with a short pause, then yield: a pure spin burns the sibling hardware
   thread for nothing, while yielding on every probe makes the barrier as
   expensive as a futex wait.  The bounded spin keeps latency in the
   submicrosecond range until a thread is genuinely descheduled. */
static void barrier_wait(Barrier* b) {
    int gen = atomic_load_explicit(&b->gen, memory_order_acquire);
    if (atomic_fetch_add_explicit(&b->count, 1, memory_order_acq_rel) + 1 == b->n) {
        atomic_store_explicit(&b->count, 0, memory_order_relaxed);
        atomic_fetch_add_explicit(&b->gen, 1, memory_order_release);
        return;
    }
    int spins = 0;
    while (atomic_load_explicit(&b->gen, memory_order_acquire) == gen) {
        if (++spins < 64) {
#if defined(__x86_64__) || defined(__i386__)
            __builtin_ia32_pause();
#elif defined(__aarch64__)
            __asm__ __volatile__("yield");
#endif
        } else {
            sched_yield();
        }
    }
}

typedef struct {
    Cloth* c;
    int block;
    int steps;
    int share;
    Barrier* bar;
} Worker;

/* A block integrates and applies its owned columns and projects every
   constraint incident to them.  It reads only iteration-start positions, so
   no cross-block ordering is needed between projection and application. */
static void* worker_main(void* arg) {
    Worker* wk = (Worker*)arg;
    Cloth* c = wk->c;
    int b = wk->block;
    int v0 = c->oc0[b] * c->h;
    int v1 = (c->oc1[b] + 1) * c->h;
    int cs = c->block_cstart[b], ce = c->block_cstart[b + 1];
    Buf* corr = &c->block_corr[b];

    memset(corr, 0, sizeof(Buf));
    for (int s = 0; s < wk->steps; s++) {
        integrate_range(c, v0, v1);
        for (int it = 0; it < c->iterations; it++) {
            barrier_wait(wk->bar);                  /* positions are stable */
            for (int i = cs; i < ce; i++)
                project_con(c, c->block_cons[i], wk->share, corr);
            barrier_wait(wk->bar);                  /* projections done */
            apply_range(c, v0, v1, corr);
        }
    }
    return NULL;
}

/* Parallel blocked solve.  Returns -1 if there are no blocks. */
static int cloth_solve_blocks(Cloth* c, int steps, int share_seams) {
    if (c->nb < 1) return -1;
    if (steps < 0) return -1;

    Barrier bar;
    barrier_init(&bar, c->nb);

    pthread_t th[MAX_BLOCKS];
    Worker wk[MAX_BLOCKS];
    int created = 0;
    for (int k = 0; k < c->nb; k++) {
        wk[k].c = c;
        wk[k].block = k;
        wk[k].steps = steps;
        wk[k].share = share_seams;
        wk[k].bar = &bar;
        if (pthread_create(&th[k], NULL, worker_main, &wk[k]) != 0) break;
        created++;
    }
    for (int k = 0; k < created; k++) pthread_join(th[k], NULL);
    return created == c->nb ? 0 : -1;
}

/* --------------------------------------------------------------- compare */

/* Largest per-vertex distance between two cloth states.  Reports how many
   vertices exceed tol and how many of those sit on a seam column, so a tear
   at a block boundary can be told apart from ordinary divergence. */
static float cloth_mismatch(const Cloth* a, const Cloth* b, float tol,
                            int* over, int* seam_over) {
    float worst = 0.0f;
    int cnt = 0, scnt = 0;
    for (int x = 0; x < a->w; x++) {
        for (int y = 0; y < a->h; y++) {
            int v = idx(a, x, y);
            float dx = a->v[v].x - b->v[v].x;
            float dy = a->v[v].y - b->v[v].y;
            float d = sqrtf(dx * dx + dy * dy);
            if (d > worst) worst = d;
            if (d > tol) {
                cnt++;
                if (a->seam_col[x]) scnt++;
            }
        }
    }
    if (over) *over = cnt;
    if (seam_over) *seam_over = scnt;
    return worst;
}

/* --------------------------------------------------------------- selftest */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static float max_stretch(const Cloth* c) {
    float worst = 0.0f;
    for (int i = 0; i < c->ncon; i++) {
        const CConstraint* k = &c->con[i];
        float dx = c->v[k->b].x - c->v[k->a].x;
        float dy = c->v[k->b].y - c->v[k->a].y;
        float e = fabsf(sqrtf(dx * dx + dy * dy) - k->rest);
        if (e > worst) worst = e;
    }
    return worst;
}

static double now_ms(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec * 1000.0 + (double)ts.tv_nsec / 1e6;
}


#define BENCH_STEPS 180
#define BENCH_REPS 5

int main(void) {
    static Cloth base, ref, blk;

    check("init", cloth_init(&base, 12, 12, 0.25f) == 0 &&
                  base.v[idx(&base, 11, 11)].pinned == 1 &&
                  base.v[idx(&base, 0, 0)].pinned == 0);
    check("reject_small", cloth_init(&blk, 1, 12, 0.25f) == -1);

    /* no blocks is an error */
    check("zero_blocks_rejected", cloth_build_blocks(&base, 0) == -1);
    check("negative_blocks_rejected", cloth_build_blocks(&base, -3) == -1);
    check("no_blocks_solve_rejected", base.nb == 0 &&
                                      cloth_solve_blocks(&base, 5, 1) == -1);
    check("build_one_block", cloth_build_blocks(&base, 1) == 0);
    check("buffer_alloc", cloth_alloc_buffers(&base) == 0);

    /* whole solve: sheet hangs, constraints stay close to rest length */
    memcpy(&ref, &base, sizeof(Cloth));
    cloth_solve_whole(&ref, 200);
    check("whole_stable", max_stretch(&ref) < 0.05f);
    float sag = ref.v[idx(&ref, 0, 11)].y - ref.v[idx(&ref, 0, 0)].y;
    check("whole_falls", sag > 1.0f);

    /* blocked solve must reproduce the whole solve exactly, for any shot */
    int nb_list[] = { 1, 2, 3, 4, 8 };
    char name[64];
    float worst_exact = 0.0f;
    for (int t = 0; t < (int)(sizeof(nb_list) / sizeof(nb_list[0])); t++) {
        int nb = nb_list[t];
        if (nb > base.w - 1) {
            snprintf(name, sizeof(name), "skip_blocks_%d", nb);
            check(name, cloth_build_blocks(&base, nb) == -1);
            continue;
        }
        check("build_blocks", cloth_build_blocks(&base, nb) == 0);
        check("buffer_alloc", cloth_alloc_buffers(&base) == 0);
        memcpy(&blk, &base, sizeof(Cloth));
        blk.block_corr = NULL;
        check("buffer_alloc_clone", cloth_alloc_buffers(&blk) == 0);
        check("solve_blocks", cloth_solve_blocks(&blk, 200, 1) == 0);
        int over = 0, seam_over = 0;
        float m = cloth_mismatch(&blk, &ref, 0.0f, &over, &seam_over);
        if (m > worst_exact) worst_exact = m;
        snprintf(name, sizeof(name), "blocks_%d_exact_match", nb);
        check(name, m == 0.0f && over == 0);
        snprintf(name, sizeof(name), "blocks_%d_seam_mismatch_zero", nb);
        check(name, seam_over == 0);
        snprintf(name, sizeof(name), "blocks_%d_stable", nb);
        check(name, max_stretch(&blk) < 0.05f);
        cloth_free(&blk);
    }
    printf("max_vertex_mismatch=%e\n", (double)worst_exact);
    (void)worst_exact;

    /* At the maximum grid size the per-block projection lists together hold
       more entries than there are constraints (shared seam constraints are
       listed twice); check the storage is sized for that. */
    check("max_grid_init", cloth_init(&base, MAX_W, MAX_H, 0.25f) == 0);
    check("max_grid_blocks", cloth_build_blocks(&base, MAX_BLOCKS) == 0);
    int max_lists = base.block_cstart[MAX_BLOCKS];
    check("list_storage_fits", max_lists <= MAX_BLOCK_CONS &&
                               max_lists > base.ncon);
    memcpy(&ref, &base, sizeof(Cloth));
    cloth_solve_whole(&ref, 20);
    check("max_grid_buffers", cloth_alloc_buffers(&base) == 0);
    memcpy(&blk, &base, sizeof(Cloth));
    blk.block_corr = NULL;
    check("max_grid_buffers_clone", cloth_alloc_buffers(&blk) == 0);
    check("max_grid_solve", cloth_solve_blocks(&blk, 20, 1) == 0);
    int mo = 0, ms = 0;
    check("max_grid_exact", cloth_mismatch(&blk, &ref, 0.0f, &mo, &ms) == 0.0f);
    cloth_free(&blk);

    /* dropped seam sharing must tear, and the compare must catch it.  This
       runs on its own small sheet so the tear numbers stay readable. */
    check("tear_init", cloth_init(&base, 16, 16, 0.25f) == 0);
    check("build_tear_case", cloth_build_blocks(&base, 4) == 0);
    memcpy(&ref, &base, sizeof(Cloth));
    cloth_solve_whole(&ref, 200);
    check("buffer_alloc", cloth_alloc_buffers(&base) == 0);
    memcpy(&blk, &base, sizeof(Cloth));
    blk.block_corr = NULL;
    check("buffer_alloc_clone", cloth_alloc_buffers(&blk) == 0);
    check("solve_torn", cloth_solve_blocks(&blk, 200, 0) == 0);
    int torn_over = 0, torn_seam = 0;
    float torn = cloth_mismatch(&blk, &ref, 1e-4f, &torn_over, &torn_seam);
    printf("tear_max=%e tear_vertices=%d tear_seam_vertices=%d\n",
           (double)torn, torn_over, torn_seam);
    check("tear_detected", torn > 1e-3f && torn_over > 0 && torn_seam > 0);
    /* torn vertices are concentrated on the seam columns, not scattered */
    check("tear_on_seams", torn_seam * 2 > torn_over && torn_seam >= 3 * 16 - 3);
    cloth_free(&blk);

    /* Block-count scaling.  The sheet is large enough that a block's own
       columns dominate the seam work and the barrier cost.  Each timing is
       the best of several runs: the median is a poor estimator on a shared
       machine where an unrelated task can steal a core mid-run. */
    check("bench_init", cloth_init(&base, 88, 88, 0.25f) == 0);

    static const int nb_bench[] = { 1, 2, 4, 8 };
    double times[4];
    for (int t = 0; t < 4; t++) {
        int nb = nb_bench[t];
        check("bench_build", cloth_build_blocks(&base, nb) == 0);
        memcpy(&ref, &base, sizeof(Cloth));
        cloth_solve_whole(&ref, BENCH_STEPS);
        double best = 1e30, exact = 0.0f;
        for (int r = 0; r < BENCH_REPS; r++) {
            memcpy(&blk, &base, sizeof(Cloth));
            blk.block_corr = NULL;
            check("bench_buffers_clone", cloth_alloc_buffers(&blk) == 0);
            double t0 = now_ms();
            int rc = cloth_solve_blocks(&blk, BENCH_STEPS, 1);
            double dt = now_ms() - t0;
            if (dt < best) best = dt;
            check("bench_solve", rc == 0);
            int over = 0, so = 0;
            float m = cloth_mismatch(&blk, &ref, 0.0f, &over, &so);
            if (m > exact) exact = m;
            cloth_free(&blk);
        }
        times[t] = best;
        check("bench_exact", exact == 0.0f);
        printf("timing blocks=%d ms=%.3f\n", nb, times[t]);
    }
    printf("speedup_2=%.2f speedup_4=%.2f speedup_8=%.2f\n",
           times[0] / times[1], times[0] / times[2], times[0] / times[3]);
    check("speedup_2", times[1] < times[0]);
    check("speedup_4", times[2] < times[1]);
    /* Past the physical core count the gain flattens (barrier cost and SMT
       sharing dominate), so the 8-block check only asks for a solid overall
       win rather than a strict step down from 4. */
    check("speedup_8", times[3] < times[0] / 2.0);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
