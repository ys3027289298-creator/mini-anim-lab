
/* constraint: cloth-style position-based distance constraints for a particle
   chain (rope).  Each segment is a bilateral distance constraint: the length
   error is split between the two endpoints in inverse-mass proportion, so a
   pinned particle (infinite mass) absorbs none of the correction and the
   remaining inertial particles share it exactly.  Convergence is driven by
   substeps (smaller integration step per solve) and Gauss-Seidel sweeps;
   every solve reports the worst remaining length error, and the best iterate
   is kept so a larger iteration budget can never make that residual worse. */
#include <stdio.h>
#include <string.h>
#include <math.h>

#define MAX_P 64
#define MAX_SEG (MAX_P - 1)
#define MIN_DIST 1e-12   /* coincident points: direction is undefined */
#define STIFF_MAX 1.0    /* over-relaxed constraints are clamped, not rejected */

typedef struct {
    double x, y;        /* position */
    double px, py;      /* previous position (Verlet velocity) */
    double mass;        /* strictly positive */
    double inv_mass;    /* 1/mass, or 0 when pinned */
    int pinned;
} Particle;

typedef struct {
    Particle p[MAX_P];
    double rest[MAX_SEG];       /* rest length of segment i (p[i]..p[i+1]) */
    double stiffness[MAX_SEG];  /* 0 disables the segment, 1 = rigid */
    int n;
    int iterations;             /* max forward+backward sweeps per substep */
    int substeps;
    double damping;
    double tolerance;           /* stop a solve below this residual (0 = off) */
    double last_residual;       /* best residual after the last solve/step */
    double peak_residual;       /* worst per-substep residual of the last step */
    int last_sweeps;            /* sweeps actually executed in the last step */
} Chain;

static int chain_init(Chain* c, int n, double seg, double x0, double y0) {
    if (!c || n < 2 || n > MAX_P || !(seg > 0.0) || !isfinite(seg)) return -1;
    memset(c, 0, sizeof(*c));
    c->n = n;
    c->iterations = 8;
    c->substeps = 1;
    c->damping = 0.99;
    for (int i = 0; i < n; i++) {
        Particle* p = &c->p[i];
        p->x = p->px = x0 + seg * i;
        p->y = p->py = y0;
        p->mass = 1.0;
        p->inv_mass = 1.0;
    }
    for (int i = 0; i < n - 1; i++) {
        c->rest[i] = seg;
        c->stiffness[i] = 1.0;
    }
    c->p[0].pinned = 1;
    c->p[0].inv_mass = 0.0;
    return 0;
}

/* zero, negative, NaN and infinite masses are invalid: a particle with no
   inertia cannot receive a finite share of the constraint impulse. */
static int chain_set_mass(Chain* c, int i, double m) {
    if (!c || i < 0 || i >= c->n) return -1;
    if (!(m > 0.0) || !isfinite(m)) return -1;
    c->p[i].mass = m;
    if (!c->p[i].pinned) c->p[i].inv_mass = 1.0 / m;
    return 0;
}

/* negative (and NaN/infinite) stiffness is invalid; >1 is clamped to rigid. */
static int chain_set_stiffness(Chain* c, int i, double k) {
    if (!c || i < 0 || i >= c->n - 1) return -1;
    if (!(k >= 0.0) || !isfinite(k)) return -1;
    c->stiffness[i] = k > STIFF_MAX ? STIFF_MAX : k;
    return 0;
}

static int chain_pin(Chain* c, int i, int on) {
    if (!c || i < 0 || i >= c->n) return -1;
    c->p[i].pinned = on ? 1 : 0;
    c->p[i].inv_mass = on ? 0.0 : 1.0 / c->p[i].mass;
    return 0;
}

/* worst absolute violation of any segment, in length units */
static double chain_residual(const Chain* c) {
    if (!c || c->n < 2) return 0.0;
    double worst = 0.0;
    for (int i = 0; i < c->n - 1; i++) {
        double dx = c->p[i + 1].x - c->p[i].x;
        double dy = c->p[i + 1].y - c->p[i].y;
        double e = fabs(hypot(dx, dy) - c->rest[i]);
        if (e > worst) worst = e;
    }
    return worst;
}

static double chain_length(const Chain* c) {
    double len = 0.0;
    for (int i = 0; i < c->n - 1; i++)
        len += hypot(c->p[i + 1].x - c->p[i].x, c->p[i + 1].y - c->p[i].y);
    return len;
}

/* One symmetric Gauss-Seidel sweep: forward then backward, which removes the
   directional bias the old forward-only update had. */
static void chain_relax(Chain* c) {
    if (!c || c->n < 2) return;
    int segs = c->n - 1;
    for (int pass = 0; pass < 2; pass++) {
        for (int q = 0; q < segs; q++) {
            int i = pass == 0 ? q : segs - 1 - q;
            Particle* a = &c->p[i];
            Particle* b = &c->p[i + 1];
            double wa = a->inv_mass, wb = b->inv_mass;
            double wsum = wa + wb;
            if (wsum <= 0.0) continue;               /* both pinned */
            double k = c->stiffness[i];
            if (k <= 0.0) continue;                  /* segment disabled */
            double dx = b->x - a->x, dy = b->y - a->y;
            double d = hypot(dx, dy);
            if (d < MIN_DIST) continue;
            /* inverse-mass split:  delta_a = -w_a/(w_a+w_b) * error */
            double diff = (d - c->rest[i]) / d * k;
            double ca = diff * (wa / wsum);
            double cb = diff * (wb / wsum);
            a->x += dx * ca;  a->y += dy * ca;
            b->x -= dx * cb;  b->y -= dy * cb;
        }
    }
}

/* Solve the current configuration with up to `iterations` sweeps.  The best
   iterate seen is restored at the end, so residual(2*budget) <= residual(budget)
   by construction even though Gauss-Seidel itself need not be monotone in the
   max norm for mixed masses and stiffnesses. */
static double chain_solve(Chain* c) {
    if (!c || c->n < 2) return 0.0;
    double best = chain_residual(c);
    double bx[MAX_P], by[MAX_P];
    for (int i = 0; i < c->n; i++) { bx[i] = c->p[i].x; by[i] = c->p[i].y; }
    int budget = c->iterations > 0 ? c->iterations : 1;
    int sweeps = 0;
    for (int it = 0; it < budget; it++) {
        if (c->tolerance > 0.0 && best <= c->tolerance) break;
        chain_relax(c);
        sweeps = it + 1;
        double r = chain_residual(c);
        if (r < best) {
            best = r;
            for (int i = 0; i < c->n; i++) { bx[i] = c->p[i].x; by[i] = c->p[i].y; }
        }
    }
    for (int i = 0; i < c->n; i++) { c->p[i].x = bx[i]; c->p[i].y = by[i]; }
    c->last_residual = best;
    c->last_sweeps = sweeps;
    return best;
}

/* Advance one frame: `substeps` Verlet integrations, each followed by a
   constraint solve.  Returns the largest residual of any substep. */
static double chain_step(Chain* c, double dt, double gx, double gy) {
    if (!c || c->n < 2) return 0.0;
    int subs = c->substeps > 0 ? c->substeps : 1;
    if (!(dt > 0.0) || !isfinite(dt)) return chain_residual(c);
    if (!isfinite(gx) || !isfinite(gy)) return chain_residual(c);
    double h = dt / (double)subs;
    double peak = 0.0;
    int sweeps = 0;
    for (int s = 0; s < subs; s++) {
        for (int i = 0; i < c->n; i++) {
            Particle* p = &c->p[i];
            if (p->pinned) continue;
            double vx = (p->x - p->px) * c->damping;
            double vy = (p->y - p->py) * c->damping;
            p->px = p->x;
            p->py = p->y;
            p->x += vx + gx * h * h;
            p->y += vy + gy * h * h;
        }
        double r = chain_solve(c);
        sweeps += c->last_sweeps;
        if (r > peak) peak = r;
    }
    c->peak_residual = peak;
    c->last_residual = chain_residual(c);
    c->last_sweeps = sweeps;
    return peak;
}

/* analytic catenary through two pins at the same height, span `span`, total
   arc length `len`.  a solves 2a sinh(span / (2a)) = len. */
static double catenary_a(double span, double len) {
    double lo = 1e-3, hi = span;
    while (2.0 * hi * sinh(span / (2.0 * hi)) > len) hi *= 2.0;
    for (int it = 0; it < 200; it++) {
        double mid = 0.5 * (lo + hi);
        double f = 2.0 * mid * sinh(span / (2.0 * mid));
        if (f > len) lo = mid; else hi = mid;
    }
    return 0.5 * (lo + hi);
}

/* node i of a discretised catenary with equal arc spacing h: it sits `s`
   arc length from the lowest point (which is at mid-span). */
static void catenary_node(int i, int n, double span, double len,
                          double* ox, double* oy) {
    double a = catenary_a(span, len);
    double h = len / (double)(n - 1);
    double s = (double)i * h - 0.5 * len;
    double u = s / a;
    *ox = 0.5 * span + a * asinh(u);
    *oy = a * sqrt(1.0 + u * u) - a * cosh(0.5 * span / a);
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

/* settle a pinned-pinned rope until it stops moving; returns frames used */
static int settle(Chain* c, double dt, double gy, double eps, int max_frames) {
    double before[MAX_P], before_y[MAX_P];
    for (int f = 0; f < max_frames; f++) {
        for (int i = 0; i < c->n; i++) { before[i] = c->p[i].x; before_y[i] = c->p[i].y; }
        chain_step(c, dt, 0.0, gy);
        double moved = 0.0;
        for (int i = 0; i < c->n; i++) {
            double d = hypot(c->p[i].x - before[i], c->p[i].y - before_y[i]);
            if (d > moved) moved = d;
        }
        if (moved < eps) return f + 1;
    }
    return max_frames;
}

int main(void) {
    Chain c;

    /* --- construction and baseline behaviour ----------------------------- */
    check("init", chain_init(&c, 8, 1.0, 0.0, 0.0) == 0);
    check("init_residual", chain_residual(&c) == 0.0);
    check("init_rejects_short", chain_init(&c, 1, 1.0, 0, 0) == -1);
    check("init_rejects_big", chain_init(&c, MAX_P + 1, 1.0, 0, 0) == -1);
    check("init_rejects_seg", chain_init(&c, 4, 0.0, 0, 0) == -1);

    check("init2", chain_init(&c, 8, 1.0, 0.0, 0.0) == 0);
    for (int i = 0; i < 200; i++) chain_step(&c, 0.016, 0.0, -9.8);
    check("stays_in_stretch", chain_residual(&c) < 0.02);
    check("pinned_anchor", fabs(c.p[0].x) < 1e-12 && fabs(c.p[0].y) < 1e-12);
    check("falls_down", c.p[7].y < -1.0);
    check("peak_reported", c.peak_residual >= c.last_residual - 1e-12);

    /* the rope must not creep longer frame after frame (the old rubber band) */
    double len_after = chain_length(&c) / 7.0;
    for (int i = 0; i < 2000; i++) chain_step(&c, 0.016, 0.0, -9.8);
    check("no_creep", chain_length(&c) / 7.0 < len_after + 1e-3);

    /* --- inverse-mass weighting ------------------------------------------ */
    Chain m;
    chain_init(&m, 2, 1.0, 0.0, 0.0);
    chain_pin(&m, 0, 0);
    check("mass_set", chain_set_mass(&m, 0, 1.0) == 0 && chain_set_mass(&m, 1, 3.0) == 0);
    m.p[1].x = m.p[1].px = 2.0;
    double com0 = (m.p[0].x * m.p[0].mass + m.p[1].x * m.p[1].mass) /
                  (m.p[0].mass + m.p[1].mass);
    chain_relax(&m);
    /* light particle moves 3x as far as the heavy one */
    check("inverse_mass_split",
          fabs(m.p[0].x - 0.75) < 1e-12 && fabs(m.p[1].x - 1.75) < 1e-12);
    check("com_preserved",
          fabs((m.p[0].x * m.p[0].mass + m.p[1].x * m.p[1].mass) /
               (m.p[0].mass + m.p[1].mass) - com0) < 1e-12);
    check("mass_resolved", chain_residual(&m) < 1e-12);

    /* equal-mass free chain: internal constraints conserve the centre of mass */
    Chain freec;
    chain_init(&freec, 5, 1.0, 0.0, 0.0);
    chain_pin(&freec, 0, 0);
    freec.p[4].x = 7.0; freec.p[4].px = 7.0;
    double comx = 0.0;
    for (int i = 0; i < 5; i++) comx += freec.p[i].x;
    chain_relax(&freec);
    double comx2 = 0.0;
    for (int i = 0; i < 5; i++) comx2 += freec.p[i].x;
    check("free_com", fabs(comx2 - comx) < 1e-12);

    /* --- input validation ------------------------------------------------- */
    Chain v;
    chain_init(&v, 4, 1.0, 0.0, 0.0);
    check("reject_zero_mass", chain_set_mass(&v, 1, 0.0) == -1);
    check("reject_neg_mass", chain_set_mass(&v, 1, -2.0) == -1);
    check("reject_nan_mass", chain_set_mass(&v, 1, NAN) == -1);
    check("reject_inf_mass", chain_set_mass(&v, 1, INFINITY) == -1);
    check("reject_bad_particle", chain_set_mass(&v, 4, 1.0) == -1);
    check("reject_neg_stiff", chain_set_stiffness(&v, 0, -0.5) == -1);
    check("reject_nan_stiff", chain_set_stiffness(&v, 0, NAN) == -1);
    check("reject_bad_segment", chain_set_stiffness(&v, 3, 1.0) == -1);
    check("clamp_stiff", chain_set_stiffness(&v, 0, 2.0) == 0 && v.stiffness[0] == 1.0);
    check("accept_stiff", chain_set_stiffness(&v, 0, 0.5) == 0 && v.stiffness[0] == 0.5);

    /* --- substeps control convergence ------------------------------------ */
    double peaks[4];
    int sub_vals[4] = { 1, 2, 4, 8 };
    Chain s;
    chain_init(&s, 16, 0.5, 0.0, 0.0);
    s.iterations = 2;
    s.damping = 1.0;
    for (int k = 0; k < 4; k++) {
        Chain t = s;
        t.substeps = sub_vals[k];
        peaks[k] = chain_step(&t, 0.05, 0.0, -9.8);
    }
    check("substeps_monotone",
          peaks[1] < peaks[0] && peaks[2] < peaks[1] && peaks[3] < peaks[2]);
    check("substeps_lower_peak", peaks[3] < 0.5 * peaks[0]);

    /* --- doubling sweeps must never raise the residual ------------------- */
    Chain base;
    chain_init(&base, 9, 1.0, 0.0, 0.0);
    for (int i = 0; i < 9; i++) {
        base.p[i].x = 1.3 * i;
        base.p[i].px = base.p[i].x;
        base.p[i].y = -0.2 * sin(0.7 * i);
        base.p[i].py = base.p[i].y;
    }
    base.iterations = 1;
    base.damping = 0.5;
    /* settle the state first so the measured residual is the solver's, not
       the initial guess's */
    for (int i = 0; i < 4000; i++) chain_step(&base, 0.01, 0.0, -9.8);
    double prev = INFINITY, first = 0.0, last = 0.0;
    int monotone = 1;
    for (int it = 1; it <= 512; it *= 2) {
        Chain t = base;
        t.iterations = it;
        double r = chain_solve(&t);
        if (it == 1) first = r;
        last = r;
        if (r > prev + 1e-15) monotone = 0;
        prev = r;
    }
    check("iters_monotone", monotone);
    check("iters_strict_gain", last < 0.1 * first);
    check("iters_residual_small", last < 1e-12);

    /* the best-iterate guard matters: raw Gauss-Seidel is not monotone in
       the max-norm for mixed masses and stiffnesses */
    Chain mesh;
    chain_init(&mesh, 12, 1.0, 0.0, 0.0);
    for (int i = 0; i < 12; i++) {
        mesh.p[i].x = 2.5 * i + 0.3 * (i & 1);
        mesh.p[i].px = mesh.p[i].x;
        mesh.p[i].y = -0.1 * sin(i);
        mesh.p[i].py = mesh.p[i].y;
        chain_set_mass(&mesh, i, 0.3 + 0.4 * i);
    }
    for (int i = 0; i < 11; i++) chain_set_stiffness(&mesh, i, 0.5 + 0.04 * i);
    int raw_bad = 0;
    {
        Chain g = mesh;
        double prev_raw = INFINITY;
        for (int it = 0; it < 64; it++) {
            chain_relax(&g);
            double rr = chain_residual(&g);
            if (rr > prev_raw + 1e-15) raw_bad = 1;
            prev_raw = rr;
        }
    }
    int guarded_bad = 0;
    {
        double prev_g = INFINITY;
        for (int it = 1; it <= 64; it *= 2) {
            Chain t = mesh;
            t.iterations = it;
            double rr = chain_solve(&t);
            if (rr > prev_g + 1e-15) guarded_bad = 1;
            prev_g = rr;
        }
    }
    printf("raw_gs_monotone=%d guarded_monotone=%d\n", !raw_bad, !guarded_bad);
    check("guard_restores_monotone", guarded_bad == 0);

    /* tolerance stops a solve early once the target is met */
    Chain tol = base;
    tol.iterations = 256;
    tol.tolerance = 1e-5;
    double tr = chain_solve(&tol);
    check("tolerance_met", tr <= 1e-5 && tol.last_sweeps < 256);

    /* --- catenary: discrete rope vs analytic static equilibrium ---------- */
    /* Two error sources bound the deviation, both second order:
         discrete: equal-chord nodes sit on the exact discrete catenary;
                   the chord offset from the smooth curve is <= h^2 kappa / 8,
                   at most a quarter of that bound is used here;
         temporal: the PBD fixed point is a projection of the Verlet step, so
                   it differs from the static equilibrium by O(dt^2). */
    const double span = 4.0, rope_len = 5.0, dt = 0.01, g = -9.8;
    int meshes[2] = { 21, 41 };
    for (int mi = 0; mi < 2; mi++) {
        int n = meshes[mi];
        double h = rope_len / (double)(n - 1);
        Chain r;
        chain_init(&r, n, h, 0.0, 0.0);
        /* start on a shallow parabola (shorter than rest: it has to stretch) */
        for (int i = 0; i < n; i++) {
            double u = (double)i / (double)(n - 1);
            r.p[i].x = r.p[i].px = span * u;
            r.p[i].y = r.p[i].py = -4.0 * 1.0 * u * (1.0 - u);
        }
        chain_pin(&r, n - 1, 1);
        r.iterations = 256;
        r.substeps = 1;
        r.damping = 0.5;
        r.tolerance = 0.0;
        int frames = settle(&r, dt, g, 1e-10, 40000);
        double err = 0.0;
        for (int i = 1; i < n - 1; i++) {
            double ax, ay;
            catenary_node(i, n, span, rope_len, &ax, &ay);
            double e = hypot(r.p[i].x - ax, r.p[i].y - ay);
            if (e > err) err = e;
        }
        double a = catenary_a(span, rope_len);
        double kappa = cosh(0.5 * span / a) / a;
        double bound = 0.25 * h * h * kappa + 2.0 * fabs(g) * dt * dt;
        printf("catenary_%d: frames=%d residual=%.3e err=%.6e bound=%.6e\n",
               n, frames, chain_residual(&r), err, bound);
        char name[64];
        snprintf(name, sizeof(name), "catenary_settles_%d", n);
        check(name, frames < 40000);
        snprintf(name, sizeof(name), "catenary_bound_%d", n);
        check(name, err <= bound);
        snprintf(name, sizeof(name), "catenary_second_order_%d", n);
        check(name, err < 0.5 * h * h * kappa + 3.0 * dt * dt);
    }

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
