
/* equiv: substep/iteration equivalence analysis for the verlet chain.
   Per-frame compute budget B = substeps * iterations (projection sweeps).
   Compares substep-heavy vs iteration-heavy residual curves, estimates the
   optimal split by golden-section search, verifies it against exhaustive
   enumeration, and reports an upper bound on the split deviation.
   A zero (or negative) total budget is an error. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define MAX_P 16
#define FRAMES 120
#define DT 0.016f

typedef struct { float x, y, px, py; float inv_mass; int pinned; } Particle;

typedef struct {
    Particle p[MAX_P];
    int n;
    float rest[MAX_P];
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

/* residual = worst segment-length error after FRAMES frames of free fall;
   per frame: `substeps` steps of DT/substeps with `iters` sweeps each.
   If res_out is given, per-frame residuals are recorded (FRAMES entries). */
static float simulate(int substeps, int iters, float* res_out) {
    Chain c;
    chain_init(&c, 8, 1.0f, 0.0f, 0.0f);
    c.iterations = iters;
    float r = 0.0f;
    for (int f = 0; f < FRAMES; f++) {
        for (int s = 0; s < substeps; s++)
            chain_step(&c, DT / substeps, 0.0f, -9.8f);
        r = chain_seg_error(&c);
        if (res_out) res_out[f] = r;
    }
    return r;
}

typedef struct {
    float curve_a[FRAMES];  /* substep-heavy:  s=B, i=1   */
    float curve_b[FRAMES];  /* iteration-heavy: s=1, i=B  */
    float res_a, res_b;
    int split_est;          /* golden-section estimate of optimal s */
    int split_opt;          /* exhaustive optimum */
    int dev;                /* |split_est - split_opt|               */
    int dev_bound;          /* golden bracket width: dev <= bound    */
    float res_est, res_opt; /* residuals at est / opt split          */
    float subopt;           /* res_est - res_opt (>= 0)              */
    int cost_ok;            /* every scanned split has s*(B/s) <= B  */
    float* scan;            /* scan[s-1] = residual at split s, len B */
} Result;

/* golden-section search for the split s in [1, B] minimising the final
   residual; assumes unimodality, *bound receives the final bracket width
   which upper-bounds |estimate - true optimum| when that holds. */
static int golden_split(int B, int* bound) {
    int lo = 1, hi = B;
    int tol = B / 16;
    if (tol < 2) tol = 2;
    while (hi - lo > tol) {
        int span = hi - lo;
        int m1 = lo + (int)(span * 0.382f);
        int m2 = hi - (int)(span * 0.382f);
        if (m1 < lo + 1) m1 = lo + 1;
        if (m2 > hi - 1) m2 = hi - 1;
        if (m1 >= m2) break;
        float f1 = simulate(m1, B / m1, NULL);
        float f2 = simulate(m2, B / m2, NULL);
        if (f1 <= f2) hi = m2; else lo = m1;
    }
    int best = lo;
    float br = simulate(lo, B / lo, NULL);
    for (int s = lo + 1; s <= hi; s++) {
        float r = simulate(s, B / s, NULL);
        if (r < br) { br = r; best = s; }
    }
    *bound = hi - lo;
    return best;
}

/* returns 0 on success, -1 if the total budget is not positive */
static int analyze(int B, Result* r) {
    if (B <= 0) return -1;
    r->res_a = simulate(B, 1, r->curve_a);
    r->res_b = simulate(1, B, r->curve_b);
    r->split_est = golden_split(B, &r->dev_bound);
    r->scan = malloc(sizeof(float) * (size_t)B);
    if (!r->scan) return -1;
    r->cost_ok = 1;
    r->split_opt = 1;
    r->res_opt = 0.0f;
    for (int s = 1; s <= B; s++) {
        int iters = B / s;
        if (s * iters > B) r->cost_ok = 0;
        float res = simulate(s, iters, NULL);
        r->scan[s - 1] = res;
        if (s == 1 || res < r->res_opt) { r->res_opt = res; r->split_opt = s; }
    }
    r->res_est = r->scan[r->split_est - 1];
    r->dev = abs(r->split_est - r->split_opt);
    r->subopt = r->res_est - r->res_opt;
    return 0;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static void write_csvs(const char* dir, int B, const Result* r) {
    char path[512];
    snprintf(path, sizeof(path), "%s/equiv_curves.csv", dir);
    FILE* f = fopen(path, "w");
    if (f) {
        fprintf(f, "frame,substep_heavy,iter_heavy\n");
        for (int i = 0; i < FRAMES; i++)
            fprintf(f, "%d,%.6g,%.6g\n", i, r->curve_a[i], r->curve_b[i]);
        fclose(f);
    }
    snprintf(path, sizeof(path), "%s/equiv_scan.csv", dir);
    f = fopen(path, "w");
    if (f) {
        fprintf(f, "substeps,iters,final_res\n");
        for (int s = 1; s <= B; s++)
            fprintf(f, "%d,%d,%.6g\n", s, B / s, r->scan[s - 1]);
        fclose(f);
    }
}

int main(int argc, char** argv) {
    int B = argc > 1 ? atoi(argv[1]) : 64;
    const char* outdir = argc > 2 ? argv[2] : NULL;
    if (B <= 0) {
        fprintf(stderr, "error: total budget must be positive (got %d)\n", B);
        return 2;
    }

    Result r;
    memset(&r, 0, sizeof(r));
    if (analyze(B, &r) != 0) {
        fprintf(stderr, "error: analysis failed for budget %d\n", B);
        return 1;
    }

    printf("budget=%d\n", B);
    printf("frames=%d\n", FRAMES);
    printf("final_res_substep_heavy=%.6g\n", r.res_a);
    printf("final_res_iter_heavy=%.6g\n", r.res_b);
    printf("split_est=%d\n", r.split_est);
    printf("split_opt=%d\n", r.split_opt);
    printf("split_dev=%d\n", r.dev);
    printf("split_dev_bound=%d\n", r.dev_bound);
    printf("resid_at_est=%.6g\n", r.res_est);
    printf("resid_at_opt=%.6g\n", r.res_opt);
    printf("resid_subopt=%.6g\n", r.subopt);

    if (outdir) write_csvs(outdir, B, &r);

    check("zero_budget_errors", analyze(0, &(Result){0}) != 0);
    check("cost_within_budget", r.cost_ok);
    check("dev_within_bound", r.dev <= r.dev_bound);
    check("subopt_nonnegative", r.subopt >= -1e-6f);

    free(r.scan);
    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
