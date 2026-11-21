
/* constraint: verlet particle chain with distance constraints, iteration
   count and damping; used as the toy model for cloth/rope style solvers.

   Stability for mixed masses (light+heavy on one constraint):
   - projection is weighted by inverse mass (light moves more);
   - effective inverse-mass ratio is regularized (capped) so a light
     particle never absorbs ~100% of a correction in one projection
     (that teleport becomes a huge verlet velocity kick -> "explosion");
   - each projection's per-particle correction is clamped to
     max_correction, bounding the velocity a single solve can inject;
   - chain_warn_mass_ratio() flags mass ratios above MASS_RATIO_WARN. */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define MAX_P 16

#define MASS_RATIO_CAP  100.0f    /* effective inv-mass ratio cap in the solver */
#define MASS_RATIO_WARN 1000.0f   /* warn when a chain's mass ratio exceeds this */

typedef struct { float x, y, px, py; float inv_mass; int pinned; } Particle;

typedef struct {
    Particle p[MAX_P];
    int n;
    float rest[MAX_P];      /* rest length between i and i+1 */
    int iterations;
    float damping;
    float max_correction;   /* per-projection correction clamp (world units) */
    int warned;             /* mass-ratio warning already emitted */
} Chain;

static void chain_init(Chain* c, int n, float seg, float x0, float y0) {
    memset(c, 0, sizeof(*c));
    c->n = n;
    c->iterations = 8;
    c->damping = 0.99f;
    c->max_correction = 0.2f;
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

/* worst pairwise mass ratio among non-pinned particles (1.0 = uniform) */
static float chain_mass_ratio(const Chain* c) {
    float worst = 1.0f;
    for (int i = 0; i < c->n; i++) {
        float wi = c->p[i].inv_mass;
        if (wi <= 0.0f) continue;
        for (int j = i + 1; j < c->n; j++) {
            float wj = c->p[j].inv_mass;
            if (wj <= 0.0f) continue;
            float r = wi > wj ? wi / wj : wj / wi;
            if (r > worst) worst = r;
        }
    }
    return worst;
}

/* warn once per chain when the mass ratio exceeds MASS_RATIO_WARN;
   returns 1 if the ratio is extreme, 0 otherwise */
static int chain_warn_mass_ratio(Chain* c) {
    float r = chain_mass_ratio(c);
    if (r <= MASS_RATIO_WARN) return 0;
    if (!c->warned) {
        fprintf(stderr,
                "warning: chain mass ratio %.0f exceeds %.0f; "
                "solver regularizes effective ratio to %.0f, expect soft constraints\n",
                r, (float)MASS_RATIO_WARN, (float)MASS_RATIO_CAP);
        c->warned = 1;
    }
    return 1;
}

static void chain_step(Chain* c, float dt, float gx, float gy) {
    chain_warn_mass_ratio(c);
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
            /* regularize extreme mass ratios: cap the effective ratio so
               the lighter particle cannot absorb ~all of the correction */
            if (wa > 0.0f && wb > 0.0f) {
                if (wa > wb * MASS_RATIO_CAP) wa = wb * MASS_RATIO_CAP;
                else if (wb > wa * MASS_RATIO_CAP) wb = wa * MASS_RATIO_CAP;
            }
            float wsum = wa + wb;
            if (wsum <= 0) continue;
            /* inverse-mass weighted projection */
            float cax = dx * diff * (wa / wsum), cay = dy * diff * (wa / wsum);
            float cbx = dx * diff * (wb / wsum), cby = dy * diff * (wb / wsum);
            /* clamp each projection's per-particle correction */
            float mc = c->max_correction;
            if (mc > 0.0f) {
                float ma = sqrtf(cax * cax + cay * cay);
                if (ma > mc) { float s = mc / ma; cax *= s; cay *= s; }
                float mb = sqrtf(cbx * cbx + cby * cby);
                if (mb > mc) { float s = mc / mb; cbx *= s; cby *= s; }
            }
            a->x += cax; a->y += cay;
            b->x -= cbx; b->y -= cby;
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

/* analytic statics: a chain hanging from a pinned anchor under gravity
   (inextensible constraints) rests straight down at (x0, y0 - i*seg);
   returns the worst particle displacement from that reference */
static float chain_statics_error(const Chain* c, float x0, float y0, float seg) {
    float worst = 0.0f;
    for (int i = 0; i < c->n; i++) {
        float ex = x0, ey = y0 - seg * i;
        float e = hypotf(c->p[i].x - ex, c->p[i].y - ey);
        if (e > worst) worst = e;
    }
    return worst;
}

/* count frames in which any particle is non-finite or beyond a sanity
   bound; must always be zero */
static int chain_divergence_frames(Chain* c, int frames, float dt, float gy, float bound) {
    int bad = 0;
    for (int f = 0; f < frames; f++) {
        chain_step(c, dt, 0.0f, gy);
        for (int i = 0; i < c->n; i++) {
            if (!isfinite(c->p[i].x) || !isfinite(c->p[i].y) ||
                fabsf(c->p[i].x) > bound || fabsf(c->p[i].y) > bound) {
                bad++;
                break;
            }
        }
    }
    return bad;
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

    /* inverse-mass weighting: with masses 1 and 4 the light particle
       takes 4/5 of a single projection, the heavy one 1/5 */
    {
        Chain w;
        chain_init(&w, 2, 1.0f, 0.0f, 0.0f);
        w.p[0].pinned = 0; w.p[0].inv_mass = 1.0f;   /* mass 1 */
        w.p[1].inv_mass = 0.25f;                     /* mass 4 */
        w.p[1].x += 1.0f; w.p[1].px = w.p[1].x;    /* stretch to 2, at rest */
        w.iterations = 1; w.damping = 1.0f;
        w.max_correction = 0.0f;                     /* disable clamp for exact check */
        chain_step(&w, 0.016f, 0.0f, 0.0f);
        float ma = fabsf(w.p[0].x - 0.0f);
        float mb = fabsf(w.p[1].x - 2.0f);
        check("weighted_projection",
              fabsf(ma - 0.8f) < 1e-4f && fabsf(mb - 0.2f) < 1e-4f);
    }

    /* correction clamp: one projection may move a particle at most
       max_correction, however large the violation */
    {
        Chain w;
        chain_init(&w, 2, 1.0f, 0.0f, 0.0f);
        w.p[0].pinned = 0; w.p[0].inv_mass = 1.0f;
        w.p[1].x += 50.0f; w.p[1].px = w.p[1].x;   /* huge overstretch, at rest */
        w.iterations = 1; w.damping = 1.0f;
        chain_step(&w, 0.016f, 0.0f, 0.0f);
        check("correction_clamped",
              fabsf(w.p[0].x) <= w.max_correction + 1e-6f &&
              fabsf(w.p[1].x - 51.0f) <= w.max_correction + 1e-6f);
    }

    /* extreme mass ratios: divergence frames must be zero */
    {
        float ratios[] = { 1.0f, 1e2f, 1e4f, 1e6f };
        int ok = 1;
        for (int r = 0; r < 4; r++) {
            Chain h;
            chain_init(&h, 6, 1.0f, 0.0f, 0.0f);
            h.p[5].inv_mass = 1.0f / ratios[r];      /* heavy tail */
            int bad = chain_divergence_frames(&h, 600, 0.016f, -9.8f, 1e4f);
            printf("divergence_frames_ratio_%.0f=%d\n", ratios[r], bad);
            if (bad != 0) ok = 0;
        }
        check("no_divergence", ok);
    }

    /* extreme mass ratios: residual stretch stays bounded */
    {
        Chain h;
        chain_init(&h, 6, 1.0f, 0.0f, 0.0f);
        h.p[5].inv_mass = 1e-6f;
        for (int i = 0; i < 600; i++) chain_step(&h, 0.016f, 0.0f, -9.8f);
        float e = chain_seg_error(&h);
        printf("heavy_tail_1e6_stretch=%.4f bound=0.5\n", e);
        check("heavy_tail_stretch_bounded", e < 0.5f);
    }

    /* statics vs analytic reference: uniform chain hangs straight down,
       displacement error bounded by 0.02 */
    {
        Chain s;
        chain_init(&s, 8, 1.0f, 0.0f, 0.0f);
        for (int i = 0; i < 2000; i++) chain_step(&s, 0.016f, 0.0f, -9.8f);
        float e = chain_statics_error(&s, 0.0f, 0.0f, 1.0f);
        printf("statics_uniform_err=%.4f bound=0.02\n", e);
        check("statics_uniform_bound", e < 0.02f);
    }

    /* statics with a heavy tail (ratio 1e3): regularization softens the
       constraint, so the bound is looser but still small (0.2) */
    {
        Chain s;
        chain_init(&s, 6, 1.0f, 0.0f, 0.0f);
        s.p[5].inv_mass = 1e-3f;
        for (int i = 0; i < 4000; i++) chain_step(&s, 0.016f, 0.0f, -9.8f);
        float e = chain_statics_error(&s, 0.0f, 0.0f, 1.0f);
        printf("statics_heavy_1e3_err=%.4f bound=0.2\n", e);
        check("statics_heavy_bound", e < 0.2f);
    }

    /* mass-ratio warning: fires above threshold, quiet for uniform */
    {
        Chain w;
        chain_init(&w, 3, 1.0f, 0.0f, 0.0f);
        check("warn_quiet_uniform", chain_warn_mass_ratio(&w) == 0);
        w.p[2].inv_mass = 1e-5f;                     /* ratio 1e5 > 1e3 */
        check("warn_fires_extreme", chain_warn_mass_ratio(&w) == 1);
    }

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
