
/* collision: verlet particle chain against a 1d heightfield ground.
   per-segment analytic penetration depth, constraint pushback along the
   ground normal, coulomb-style clamp on tangential motion. */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define MAX_P 16
#define MAX_ITERS 32

/* ---------- heightfield ---------- */

typedef struct {
    const float* h;     /* n samples, spacing dx, first sample at x0 */
    int n;
    float x0, dx;
} HeightField;

static int hf_valid(const HeightField* hf) {
    if (!hf || !hf->h) return 0;
    if (hf->n < 2 || !(hf->dx > 0.0f)) return 0;
    if (!isfinite(hf->x0)) return 0;
    for (int i = 0; i < hf->n; i++)
        if (!isfinite(hf->h[i])) return 0;
    return 1;
}

static float hf_height(const HeightField* hf, float x) {
    float fx = (x - hf->x0) / hf->dx;
    if (fx <= 0.0f) return hf->h[0];
    if (fx >= (float)(hf->n - 1)) return hf->h[hf->n - 1];
    int i = (int)fx;
    float f = fx - (float)i;
    return hf->h[i] * (1.0f - f) + hf->h[i + 1] * f;
}

static void hf_normal(const HeightField* hf, float x, float* nx, float* ny) {
    float e = hf->dx * 0.5f;
    float slope = (hf_height(hf, x + e) - hf_height(hf, x - e)) / (2.0f * e);
    float len = sqrtf(slope * slope + 1.0f);
    *nx = -slope / len;
    *ny = 1.0f / len;
}

/* ---------- analytic segment penetration ---------- */

/* max over t in [0,1] of ground(x(t)) - y(t); both terms are piecewise
   linear in t so the max sits at t=0, t=1, or a heightfield knot crossing.
   returns 0 on success, -1 when the input is rejected (invalid field,
   non-finite endpoints, zero-length segment). */
static int seg_penetration(const HeightField* hf,
                           float ax, float ay, float bx, float by,
                           float* out_depth, float* out_t) {
    if (!hf_valid(hf)) return -1;
    if (!isfinite(ax) || !isfinite(ay) || !isfinite(bx) || !isfinite(by)) return -1;
    float sx = bx - ax, sy = by - ay;
    if (sx * sx + sy * sy < 1e-18f) return -1;

    float best = -1e30f, best_t = 0.0f;
    /* endpoints */
    for (int k = 0; k < 2; k++) {
        float t = (float)k;
        float x = ax + sx * t, y = ay + sy * t;
        float d = hf_height(hf, x) - y;
        if (d > best) { best = d; best_t = t; }
    }
    /* knot crossings of x(t) */
    if (fabsf(sx) > 1e-12f) {
        float lo = sx > 0.0f ? ax : bx;
        float hi = sx > 0.0f ? bx : ax;
        int k0 = (int)ceilf((lo - hf->x0) / hf->dx);
        int k1 = (int)floorf((hi - hf->x0) / hf->dx);
        if (k0 < 1) k0 = 1;
        if (k1 > hf->n - 1) k1 = hf->n - 1;
        for (int k = k0; k <= k1; k++) {
            float xk = hf->x0 + hf->dx * (float)k;
            float t = (xk - ax) / sx;
            if (t <= 0.0f || t >= 1.0f) continue;
            float y = ay + sy * t;
            float d = hf_height(hf, xk) - y;
            if (d > best) { best = d; best_t = t; }
        }
    }
    *out_depth = best;
    *out_t = best_t;
    return 0;
}

/* brute-force scan, used only to cross-check the analytic value */
static float seg_penetration_scan(const HeightField* hf,
                                  float ax, float ay, float bx, float by,
                                  int samples) {
    float best = -1e30f;
    for (int i = 0; i <= samples; i++) {
        float t = (float)i / (float)samples;
        float x = ax + (bx - ax) * t, y = ay + (by - ay) * t;
        float d = hf_height(hf, x) - y;
        if (d > best) best = d;
    }
    return best;
}

/* ---------- chain with ground contact ---------- */

typedef struct { float x, y, px, py; float inv_mass; int pinned; } Particle;

typedef struct {
    Particle p[MAX_P];
    int n;
    float rest[MAX_P];
    int iterations;
    float damping;
    float mu;           /* friction coefficient */
} ColChain;

static void col_chain_init(ColChain* c, int n, float seg, float x0, float y0) {
    memset(c, 0, sizeof(*c));
    c->n = n;
    c->iterations = 12;
    c->damping = 0.99f;
    c->mu = 0.6f;
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

/* clamp tangential displacement to mu * normal correction */
static float clamp_tangential(float t, float normal_corr, float mu) {
    float lim = mu * normal_corr;
    if (t > lim) return lim;
    if (t < -lim) return -lim;
    return t;
}

/* current worst positive penetration over all segments */
static float col_residual(const ColChain* c, const HeightField* hf) {
    float worst = 0.0f;
    for (int i = 0; i < c->n - 1; i++) {
        float d, t;
        if (seg_penetration(hf, c->p[i].x, c->p[i].y,
                            c->p[i + 1].x, c->p[i + 1].y, &d, &t) != 0)
            continue;
        if (d > worst) worst = d;
    }
    return worst;
}

/* one simulation step; if residuals is non-null, stores the max segment
   penetration after each solver iteration (up to max_res entries). */
static int col_chain_step(ColChain* c, const HeightField* hf,
                          float dt, float gx, float gy,
                          float* residuals, int max_res) {
    if (!hf_valid(hf)) return -1;
    float pre_x[MAX_P], pre_y[MAX_P];
    float nc[MAX_P], cnx[MAX_P], cny[MAX_P];
    for (int i = 0; i < c->n; i++) {
        pre_x[i] = c->p[i].x; pre_y[i] = c->p[i].y;
        nc[i] = 0.0f; cnx[i] = 0.0f; cny[i] = 1.0f;
    }

    for (int i = 0; i < c->n; i++) {
        Particle* p = &c->p[i];
        if (p->pinned) continue;
        float vx = (p->x - p->px) * c->damping;
        float vy = (p->y - p->py) * c->damping;
        p->px = p->x; p->py = p->y;
        p->x += vx + gx * dt * dt;
        p->y += vy + gy * dt * dt;
    }

    int iters = c->iterations;
    if (iters > MAX_ITERS) iters = MAX_ITERS;
    for (int it = 0; it < iters; it++) {
        /* distance constraints */
        for (int i = 0; i < c->n - 1; i++) {
            Particle* a = &c->p[i];
            Particle* b = &c->p[i + 1];
            float dx = b->x - a->x, dy = b->y - a->y;
            float d = sqrtf(dx * dx + dy * dy);
            if (d < 1e-6f) continue;
            float diff = (d - c->rest[i]) / d;
            float wa = a->inv_mass, wb = b->inv_mass;
            float wsum = wa + wb;
            if (wsum <= 0.0f) continue;
            a->x += dx * diff * (wa / wsum);
            a->y += dy * diff * (wa / wsum);
            b->x -= dx * diff * (wb / wsum);
            b->y -= dy * diff * (wb / wsum);
        }
        /* segment/ground pushback */
        for (int i = 0; i < c->n - 1; i++) {
            Particle* a = &c->p[i];
            Particle* b = &c->p[i + 1];
            float depth, t;
            if (seg_penetration(hf, a->x, a->y, b->x, b->y, &depth, &t) != 0)
                continue;
            if (depth <= 0.0f) continue;
            float cx = a->x + (b->x - a->x) * t;
            float nx, ny;
            hf_normal(hf, cx, &nx, &ny);
            /* move the contact point out by depth along the normal,
               split by barycentric-weighted inverse mass */
            float wa = a->inv_mass * (1.0f - t) * (1.0f - t);
            float wb = b->inv_mass * t * t;
            float ws = wa + wb;
            if (ws <= 0.0f) continue;
            float ca = depth * a->inv_mass * (1.0f - t) / ws;
            float cb = depth * b->inv_mass * t / ws;
            a->x += nx * ca; a->y += ny * ca;
            b->x += nx * cb; b->y += ny * cb;
            nc[i] += ca; nc[i + 1] += cb;
            cnx[i] = nx; cny[i] = ny;
            cnx[i + 1] = nx; cny[i + 1] = ny;
        }
        if (residuals && it < max_res)
            residuals[it] = col_residual(c, hf);
    }

    /* friction: tangential motion since step start limited by mu * normal
       correction accumulated this step */
    for (int i = 0; i < c->n; i++) {
        if (nc[i] <= 0.0f || c->p[i].pinned) continue;
        Particle* p = &c->p[i];
        float dx = p->x - pre_x[i], dy = p->y - pre_y[i];
        float dn = dx * cnx[i] + dy * cny[i];
        float tx = dx - cnx[i] * dn, ty = dy - cny[i] * dn;
        float tl = sqrtf(tx * tx + ty * ty);
        float lim = c->mu * nc[i];
        if (tl > lim && tl > 1e-12f) {
            float s = lim / tl;
            p->x = pre_x[i] + cnx[i] * dn + tx * s;
            p->y = pre_y[i] + cny[i] * dn + ty * s;
        }
    }
    return 0;
}

/* ---------- self test ---------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static unsigned rng_state = 12345u;
static float frand(float lo, float hi) {
    rng_state = rng_state * 1664525u + 1013904223u;
    float u = (float)(rng_state >> 8) / 16777216.0f;
    return lo + (hi - lo) * u;
}

int main(void) {
    float flat_h[] = { 0.0f, 0.0f };
    HeightField flat = { flat_h, 2, -100.0f, 200.0f };
    float bump_h[] = { 0.2f, 0.5f, 0.1f, 0.8f, 0.3f, 0.6f, 0.2f, 0.4f };
    HeightField bump = { bump_h, 8, 0.0f, 0.5f };

    /* rejection of invalid input */
    float bad_h[] = { 0.0f, NAN };
    HeightField bad = { bad_h, 2, 0.0f, 1.0f };
    check("rejects_nan_ground", !hf_valid(&bad));
    float inf_h[] = { 0.0f, INFINITY };
    HeightField inf = { inf_h, 2, 0.0f, 1.0f };
    check("rejects_inf_ground", !hf_valid(&inf));
    HeightField baddx = { flat_h, 2, 0.0f, 0.0f };
    check("rejects_zero_spacing", !hf_valid(&baddx));
    float d, t;
    check("rejects_zero_len_seg",
          seg_penetration(&flat, 1.0f, 1.0f, 1.0f, 1.0f, &d, &t) == -1);
    check("rejects_nan_seg",
          seg_penetration(&flat, 0.0f, NAN, 1.0f, 1.0f, &d, &t) == -1);
    check("accepts_valid",
          seg_penetration(&bump, 0.0f, 1.0f, 3.0f, 1.0f, &d, &t) == 0);

    /* analytic depth matches a dense scan, and bounds any coarse scan */
    int analytic_ok = 1, upper_ok = 1;
    for (int i = 0; i < 200; i++) {
        float ax = frand(-0.5f, 3.5f), ay = frand(-0.5f, 1.5f);
        float bx = frand(-0.5f, 3.5f), by = frand(-0.5f, 1.5f);
        float depth, tt;
        if (seg_penetration(&bump, ax, ay, bx, by, &depth, &tt) != 0) continue;
        float fine = seg_penetration_scan(&bump, ax, ay, bx, by, 20000);
        float coarse = seg_penetration_scan(&bump, ax, ay, bx, by, 64);
        if (depth < fine - 1e-3f || depth > fine + 1e-3f) analytic_ok = 0;
        if (depth < coarse - 1e-6f) upper_ok = 0;
    }
    check("analytic_matches_scan", analytic_ok);
    check("max_pen_is_upper_bound", upper_ok);

    /* friction clamp unit behavior */
    check("friction_clamps", fabsf(clamp_tangential(1.0f, 0.1f, 0.5f) - 0.05f) < 1e-6f);
    check("friction_zero_mu", clamp_tangential(1.0f, 0.1f, 0.0f) == 0.0f);
    check("friction_allows_small", fabsf(clamp_tangential(0.01f, 0.1f, 0.5f) - 0.01f) < 1e-6f);

    /* pushback residual converges within one step */
    ColChain c;
    col_chain_init(&c, 8, 1.0f, 0.0f, 0.5f);
    c.p[0].pinned = 0; c.p[0].inv_mass = 1.0f;   /* free rope */
    for (int i = 0; i < c.n; i++) c.p[i].y = -0.3f;  /* start underground */
    float res[MAX_ITERS] = { 0 };
    check("step_rejects_bad_ground", col_chain_step(&c, &bad, 0.016f, 0.0f, -9.8f, 0, 0) == -1);
    check("step_ok", col_chain_step(&c, &flat, 0.016f, 0.0f, -9.8f, res, MAX_ITERS) == 0);
    int monotone = 1;
    for (int i = 1; i < c.iterations; i++)
        if (res[i] > res[i - 1] + 1e-5f) monotone = 0;
    check("residual_nonincreasing", monotone);
    check("residual_converges", res[c.iterations - 1] < 1e-3f);
    check("residual_shrinks", res[c.iterations - 1] <= res[0] * 0.05f + 1e-6f);

    /* long run: rope rests on the ground, nothing sinks or tunnels */
    ColChain rope;
    col_chain_init(&rope, 10, 0.5f, -1.0f, 2.0f);
    for (int i = 0; i < 400; i++)
        col_chain_step(&rope, &flat, 0.016f, 0.0f, -9.8f, 0, 0);
    check("no_sink_after_settle", col_residual(&rope, &flat) < 1e-3f);
    float min_y = 1e30f;
    for (int i = 0; i < rope.n; i++) if (rope.p[i].y < min_y) min_y = rope.p[i].y;
    check("no_tunnel_below", min_y > -0.05f);

    /* rope on bumpy ground follows the terrain */
    ColChain rb;
    col_chain_init(&rb, 8, 0.5f, 0.2f, 2.0f);
    rb.p[0].pinned = 0; rb.p[0].inv_mass = 1.0f;
    for (int i = 0; i < 400; i++)
        col_chain_step(&rb, &bump, 0.016f, 0.0f, -9.8f, 0, 0);
    check("no_sink_on_terrain", col_residual(&rb, &bump) < 5e-3f);

    /* friction: sliding rope with mu=0 barely slides, high mu keeps sliding */
    ColChain s0, s1;
    col_chain_init(&s0, 6, 0.5f, 0.0f, 1.0f);
    col_chain_init(&s1, 6, 0.5f, 0.0f, 1.0f);
    s0.p[0].pinned = 0; s0.p[0].inv_mass = 1.0f;
    s1.p[0].pinned = 0; s1.p[0].inv_mass = 1.0f;
    s0.mu = 0.0f; s1.mu = 2.0f;
    for (int i = 0; i < s0.n; i++) {   /* initial horizontal velocity */
        s0.p[i].px -= 0.05f; s1.p[i].px -= 0.05f;
    }
    for (int i = 0; i < 300; i++) {
        col_chain_step(&s0, &flat, 0.016f, 0.0f, -9.8f, 0, 0);
        col_chain_step(&s1, &flat, 0.016f, 0.0f, -9.8f, 0, 0);
    }
    check("friction_limits_slide",
          s1.p[s1.n - 1].x > s0.p[s0.n - 1].x + 0.3f);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
