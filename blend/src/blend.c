
/* blend: dual-parameter (fade duration x easing exponent) transition optimizer.
 *
 * Decision variables:
 *   T     fade duration, clamped to [tmin, tmax]
 *   alpha exponent of power easing w(u) = u^alpha, u in [0,1], clamped to [amin, amax]
 * Blended end-effector velocity: v(t) = (1-w(u))*vA(t) + w(u)*vB(t), u = t/T.
 * End displacement error is verified by numerical integration:
 *   e = | integral_0^T (vB(t) - v(t)) dt |  (reference B displacement minus blended one)
 * Feasibility requires e <= emax and peak blended acceleration <= acc_max.
 * Two solvers are provided and cross-checked in the self test:
 *   - exhaustive grid search (reference solution)
 *   - multi-start coordinate pattern search (fast solution)
 * Infeasible constraints are reported (SOLVE_INFEASIBLE), never replaced by defaults.
 */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define PROF_MAX 8
#define QUAD_N  128   /* Simpson substeps, must be even */
#define CHECK_N 512   /* independent trapezoid verification grid */

/* validation bounds advertised by the self test */
#define COST_BOUND 5e-3f   /* max allowed |J_fast - J_exhaustive|, normalized cost units */
#define ERR_BOUND  2e-3f   /* max allowed |e_fast - e_exhaustive|, meters */
#define QUAD_BOUND 1e-4f   /* max Simpson vs trapezoid discrepancy, meters */

typedef struct { float x, y; } V2;

/* piecewise-constant end-effector velocity profile of an action clip */
typedef struct {
    int n;
    float t[PROF_MAX];   /* segment start times, strictly increasing, t[0] = 0 */
    V2 v[PROF_MAX];      /* velocity held on [t[i], t[i+1]) */
} VelProfile;

typedef struct {
    VelProfile va, vb;
    float tmin, tmax;          /* duration bounds, seconds */
    float amin, amax;          /* easing exponent bounds, strictly positive */
    float emax;                /* max end displacement error, meters */
    float acc_max;             /* max blended acceleration, m/s^2 */
    float start_acc_max;       /* max acceleration at the switch instant, 0 = disabled */
    float w_dur, w_err, w_acc; /* cost weights, non-negative, positive sum */
} BlendProblem;

typedef struct {
    float T, alpha;
    V2 err_vec;     /* reference B displacement minus blended displacement */
    float err;      /* magnitude of err_vec, meters */
    float acc_peak; /* peak blended acceleration, m/s^2 */
    float start_acc;/* forward-difference acceleration at the switch instant */
    float cost;     /* normalized weighted cost */
    int feasible;
} Candidate;

typedef enum { SOLVE_OK = 0, SOLVE_INFEASIBLE = 1, SOLVE_BAD_INPUT = 2 } SolveStatus;

typedef struct {
    SolveStatus status;
    Candidate cand;        /* cand.T == -1 when status != SOLVE_OK */
    int evals;
    float min_violation;   /* diagnostic: smallest max-normalized constraint overshoot */
} SolveResult;

/* ---------- small vector / profile helpers ---------- */

static V2 vsub2(V2 a, V2 b) { V2 r = { a.x - b.x, a.y - b.y }; return r; }
static float vlen2(V2 a) { return hypotf(a.x, a.y); }

static V2 vel_at(const VelProfile* p, float t) {
    V2 zero = { 0, 0 };
    if (p->n <= 0) return zero;
    int seg = 0;
    for (int i = 1; i < p->n; i++) {
        if (t >= p->t[i]) seg = i;
        else break;
    }
    return p->v[seg];
}

static float clampf(float x, float lo, float hi) {
    return x < lo ? lo : (x > hi ? hi : x);
}

static float ease_weight(float u, float alpha) { return powf(u, alpha); }

static int problem_valid(const BlendProblem* p) {
    if (p->va.n <= 0 || p->vb.n <= 0) return 0;
    if (p->va.n > PROF_MAX || p->vb.n > PROF_MAX) return 0;
    if (p->va.t[0] != 0.0f || p->vb.t[0] != 0.0f) return 0;
    for (int i = 1; i < p->va.n; i++) if (p->va.t[i] <= p->va.t[i - 1]) return 0;
    for (int i = 1; i < p->vb.n; i++) if (p->vb.t[i] <= p->vb.t[i - 1]) return 0;
    if (!(p->tmin > 0.0f && p->tmax >= p->tmin)) return 0;
    if (!(p->amin > 0.0f && p->amax >= p->amin)) return 0;
    if (!(p->emax > 0.0f && p->acc_max > 0.0f)) return 0;
    if (p->start_acc_max < 0.0f) return 0;
    if (p->w_dur < 0 || p->w_err < 0 || p->w_acc < 0) return 0;
    if (!(p->w_dur + p->w_err + p->w_acc > 0)) return 0;
    return 1;
}

/* ---------- integration ---------- */

static float simpson(const float* f, float h) {
    float s = f[0] + f[QUAD_N];
    for (int i = 1; i < QUAD_N; i++) s += (i & 1) ? 4 * f[i] : 2 * f[i];
    return s * h / 3.0f;
}

/* independent cross-check: composite trapezoid on a finer grid */
static V2 integrate_error_trap(const BlendProblem* p, float T, float alpha, int n) {
    float h = T / n;
    V2 acc = { 0, 0 };
    for (int i = 0; i <= n; i++) {
        float t = i * h;
        float w = ease_weight(t / T, alpha);
        V2 d = vsub2(vel_at(&p->vb, t), vel_at(&p->va, t));
        float k = (1.0f - w) * (i == 0 || i == n ? 0.5f : 1.0f);
        acc.x += k * d.x;
        acc.y += k * d.y;
    }
    acc.x *= h; acc.y *= h;
    return acc;
}

/* ---------- candidate evaluation ---------- */

static Candidate evaluate(const BlendProblem* p, float T, float alpha) {
    Candidate c;
    memset(&c, 0, sizeof(c));
    c.T = T; c.alpha = alpha;

    float h = T / QUAD_N;
    V2 vblend[QUAD_N + 1];
    float fx[QUAD_N + 1], fy[QUAD_N + 1], bx[QUAD_N + 1], by[QUAD_N + 1];

    for (int i = 0; i <= QUAD_N; i++) {
        float t = i * h;
        float w = ease_weight((float)i / QUAD_N, alpha);
        V2 va = vel_at(&p->va, t);
        V2 vb = vel_at(&p->vb, t);
        vblend[i].x = (1 - w) * va.x + w * vb.x;
        vblend[i].y = (1 - w) * va.y + w * vb.y;
        fx[i] = vblend[i].x; fy[i] = vblend[i].y;
        bx[i] = vb.x;          by[i] = vb.y;
    }

    /* displacements over the fade window */
    V2 disp_blend = { simpson(fx, h), simpson(fy, h) };
    V2 disp_b     = { simpson(bx, h), simpson(by, h) };
    c.err_vec = vsub2(disp_b, disp_blend);
    c.err = vlen2(c.err_vec);

    /* peak acceleration from blended velocity differences */
    for (int i = 0; i <= QUAD_N; i++) {
        V2 a;
        if (i == 0) {
            a = vsub2(vblend[1], vblend[0]);
            a.x /= h; a.y /= h;
            c.start_acc = vlen2(a);   /* switch-instant forward difference */
        } else if (i == QUAD_N) {
            a = vsub2(vblend[QUAD_N], vblend[QUAD_N - 1]);
            a.x /= h; a.y /= h;
        } else {
            a = vsub2(vblend[i + 1], vblend[i - 1]);
            a.x /= 2 * h; a.y /= 2 * h;
        }
        float m = vlen2(a);
        if (m > c.acc_peak) c.acc_peak = m;
    }

    float wsum = p->w_dur + p->w_err + p->w_acc;
    c.cost = (p->w_dur * T / p->tmax
            + p->w_err * c.err / p->emax
            + p->w_acc * c.acc_peak / p->acc_max) / wsum;

    c.feasible =
        T >= p->tmin && T <= p->tmax &&
        alpha >= p->amin && alpha <= p->amax &&
        c.err <= p->emax + 1e-7f &&
        c.acc_peak <= p->acc_max + 1e-7f &&
        (p->start_acc_max == 0.0f || c.start_acc <= p->start_acc_max + 1e-7f);
    return c;
}

/* barrier used by the fast solver so it can cross infeasible regions */
static float penalized(const BlendProblem* p, const Candidate* c) {
    if (c->feasible) return c->cost;
    float oe = c->err > p->emax ? (c->err - p->emax) / p->emax : 0.0f;
    float oa = c->acc_peak > p->acc_max ? (c->acc_peak - p->acc_max) / p->acc_max : 0.0f;
    float os = (p->start_acc_max > 0.0f && c->start_acc > p->start_acc_max)
               ? (c->start_acc - p->start_acc_max) / p->start_acc_max : 0.0f;
    return 1e6f + 1e3f * (oe + oa + os);
}

static float violation_of(const BlendProblem* p, const Candidate* c) {
    float oe = c->err > p->emax ? (c->err - p->emax) / p->emax : 0.0f;
    float oa = c->acc_peak > p->acc_max ? (c->acc_peak - p->acc_max) / p->acc_max : 0.0f;
    float os = (p->start_acc_max > 0.0f && c->start_acc > p->start_acc_max)
               ? (c->start_acc - p->start_acc_max) / p->start_acc_max : 0.0f;
    return fmaxf(fmaxf(oe, oa), os);
}

static SolveResult bad_result(SolveStatus st) {
    SolveResult r;
    memset(&r, 0, sizeof(r));
    r.status = st;
    r.cand.T = -1.0f; r.cand.alpha = -1.0f;
    r.min_violation = INFINITY;
    return r;
}

/* ---------- constraint-ridge refinement (shared by both solvers) ----------
 * At a fixed T, feasibility in alpha forms short intervals (acceleration limits
 * cut the small-alpha side, the displacement-error budget cuts the large-alpha
 * side). A ternary contraction over T re-brackets those intervals at every
 * probe and follows the active constraint ridge to high precision. Without
 * this, derivative-free search stalls on the narrow curved feasible sliver.
 */

/* At a fixed T, feasibility in alpha is the union of up to a few disjoint
 * intervals: small alpha fails the start/peak acceleration limits, large
 * alpha fails the displacement error limit, and peak acceleration itself is
 * non-monotone in alpha, so the box midpoint need not be feasible.
 * Two passes: a coarse scan locates feasible runs, bisection brackets edges. */
static int feasible_alpha_intervals(const BlendProblem* p, float T,
                                    float lo_out[], float hi_out[],
                                    int max_seg, int* evals) {
    enum { SCAN = 24 };
    int feas[SCAN];
    float av[SCAN];
    for (int i = 0; i < SCAN; i++) {
        av[i] = p->amin + (p->amax - p->amin) * i / (SCAN - 1);
        Candidate c = evaluate(p, T, av[i]); (*evals)++;
        feas[i] = c.feasible;
    }

    int nseg = 0;
    int i = 0;
    while (i < SCAN && nseg < max_seg) {
        while (i < SCAN && !feas[i]) i++;
        if (i >= SCAN) break;
        int run_start = i;                       /* first feasible scan index */
        while (i < SCAN && feas[i]) i++;
        int run_end = i - 1;                     /* last feasible scan index */

        float seg_lo = (run_start == 0) ? p->amin : 0.0f;
        float seg_hi = (run_end == SCAN - 1) ? p->amax : 0.0f;

        if (seg_lo == 0.0f) {                    /* bracket lower edge */
            float l = av[run_start - 1], h = av[run_start];
            for (int k = 0; k < 30; k++) {
                float m = 0.5f * (l + h);
                Candidate cm = evaluate(p, T, m); (*evals)++;
                if (cm.feasible) h = m; else l = m;
            }
            seg_lo = h;
        }
        if (seg_hi == 0.0f) {                    /* bracket upper edge */
            float l = av[run_end], h = av[run_end + 1];
            for (int k = 0; k < 30; k++) {
                float m = 0.5f * (l + h);
                Candidate cm = evaluate(p, T, m); (*evals)++;
                if (cm.feasible) l = m; else h = m;
            }
            seg_hi = l;
        }
        lo_out[nseg] = seg_lo; hi_out[nseg] = seg_hi; nseg++;
    }
    return nseg;
}

/* best feasible candidate at fixed T; every feasible interval contributes its
 * two edges and two interior samples, covering interior minima and ridges */
static int best_at_T(const BlendProblem* p, float T, Candidate* out, int* evals) {
    float lo[4], hi[4];
    int nseg = feasible_alpha_intervals(p, T, lo, hi, 4, evals);
    if (nseg == 0) return 0;
    float best_cost = INFINITY;
    int found = 0;
    for (int s = 0; s < nseg; s++) {
        float as[4] = { lo[s], hi[s],
                        lo[s] + 0.33f * (hi[s] - lo[s]),
                        lo[s] + 0.67f * (hi[s] - lo[s]) };
        for (int k = 0; k < 4; k++) {
            Candidate c = evaluate(p, T, as[k]); (*evals)++;
            if (c.feasible && c.cost < best_cost) { best_cost = c.cost; *out = c; found = 1; }
        }
    }
    return found;
}

/* like best_at_T, but reuses the previous probe's alpha bracket.
 * Ridges move continuously with T, so a local bisection starting from the old
 * bracket is enough; a full 24-point rescan only runs if the bracket goes dry. */
static int best_at_T_near(const BlendProblem* p, float T,
                          float guess_lo, float guess_hi,
                          Candidate* out, int* evals) {
    float l = guess_lo, h = guess_hi;
    int iters = 0;
    while (iters < 24) {
        Candidate cl = evaluate(p, T, l); (*evals)++;
        Candidate ch = evaluate(p, T, h); (*evals)++;
        if (cl.feasible && ch.feasible) break;
        float pad = 0.25f * (h - l);
        if (!cl.feasible) { l = l - pad; if (l < p->amin) l = p->amin; }
        if (!ch.feasible) { h = h + pad; if (h > p->amax) h = p->amax; }
        iters++;
        if (l <= p->amin && h >= p->amax &&
            (!cl.feasible || !ch.feasible)) return best_at_T(p, T, out, evals);
    }
    /* tighten to the actual feasibility edges within the working bracket */
    float bl = l, bh = h;
    Candidate cbl = evaluate(p, T, bl); (*evals)++;
    if (!cbl.feasible) {
        float a = bl, b = bh;
        for (int k = 0; k < 20; k++) {
            float m = 0.5f * (a + b);
            Candidate cm = evaluate(p, T, m); (*evals)++;
            if (cm.feasible) b = m; else a = m;
        }
        bl = b;
    }
    float tl = l, th = h;
    Candidate cth = evaluate(p, T, th); (*evals)++;
    if (!cth.feasible) {
        float a = tl, b = th;
        for (int k = 0; k < 20; k++) {
            float m = 0.5f * (a + b);
            Candidate cm = evaluate(p, T, m); (*evals)++;
            if (cm.feasible) a = m; else b = m;
        }
        th = a;
    }
    if (th < bl) return 0;
    float as[4] = { bl, th, bl + 0.33f * (th - bl), bl + 0.67f * (th - bl) };
    float best_cost = INFINITY;
    int found = 0;
    for (int k = 0; k < 4; k++) {
        Candidate c = evaluate(p, T, as[k]); (*evals)++;
        if (c.feasible && c.cost < best_cost) { best_cost = c.cost; *out = c; found = 1; }
    }
    return found;
}

/* ternary contraction along T; alpha brackets are carried between probes.
 * 20 iterations shrink the T bracket by (2/3)^20 ~ 3e-4, which resolves the
 * shallow minima on the narrow feasible slivers seen in practice. */
static void ridge_refine(const BlendProblem* p, float T_lo, float T_hi,
                         float a_guess_lo, float a_guess_hi,
                         Candidate* best, float* min_violation, int* evals) {
    for (int k = 0; k < 20; k++) {
        float m1 = T_lo + (T_hi - T_lo) / 3.0f;
        float m2 = T_hi - (T_hi - T_lo) / 3.0f;
        Candidate c1, c2;
        int f1 = best_at_T_near(p, m1, a_guess_lo, a_guess_hi, &c1, evals);
        int f2 = best_at_T_near(p, m2, a_guess_lo, a_guess_hi, &c2, evals);
        if (f1) { a_guess_lo = c1.alpha - 0.05f; a_guess_hi = c1.alpha + 0.05f; }
        if (f2) { a_guess_lo = fminf(a_guess_lo, c2.alpha - 0.05f);
                  a_guess_hi = fmaxf(a_guess_hi, c2.alpha + 0.05f); }
        if (a_guess_lo < p->amin) a_guess_lo = p->amin;
        if (a_guess_hi > p->amax) a_guess_hi = p->amax;
        if (f1 && c1.cost < best->cost) *best = c1;
        if (f2 && c2.cost < best->cost) *best = c2;
        if (!f1 || (f2 && c2.cost < c1.cost)) T_lo = m1; else T_hi = m2;
        if (!f1 || !f2) {
            Candidate probes[2] = { c1, c2 };
            for (int q = 0; q < 2; q++) {
                float v = violation_of(p, &probes[q]);
                if (v < *min_violation) *min_violation = v;
            }
        }
    }
}

/* global entry: locate feasible T brackets on a coarse scan, refine each */
static Candidate ridge_solve(const BlendProblem* p, int scan_n,
                             float* min_violation, int* evals) {
    Candidate best; memset(&best, 0, sizeof(best));
    best.cost = INFINITY;

    float t_lo[32], t_hi[32], a_lo[32], a_hi[32];
    int nseg = 0;
    float open_lo = p->tmin;
    int was_feasible = 0;
    float last_a_lo = p->amin, last_a_hi = p->amax;
    for (int i = 0; i < scan_n; i++) {
        float T = p->tmin + (p->tmax - p->tmin) * i / (scan_n - 1);
        Candidate c;
        int ok = best_at_T(p, T, &c, evals);
        if (ok) {
            if (c.cost < best.cost) best = c;
            if (!was_feasible) {
                open_lo = (i > 0)
                    ? p->tmin + (p->tmax - p->tmin) * (i - 1) / (scan_n - 1)
                    : p->tmin;
            }
            last_a_lo = c.alpha - 0.1f; last_a_hi = c.alpha + 0.1f;
            if (last_a_lo < p->amin) last_a_lo = p->amin;
            if (last_a_hi > p->amax) last_a_hi = p->amax;
            was_feasible = 1;
        } else {
            if (was_feasible && nseg < 32) {
                t_lo[nseg] = open_lo; t_hi[nseg] = T;
                a_lo[nseg] = last_a_lo; a_hi[nseg] = last_a_hi; nseg++;
            }
            was_feasible = 0;
            Candidate raw = evaluate(p, T, 0.5f * (p->amin + p->amax));
            (*evals)++;
            float v = violation_of(p, &raw);
            if (v < *min_violation) *min_violation = v;
        }
    }
    if (was_feasible && nseg < 32) {
        t_lo[nseg] = open_lo; t_hi[nseg] = p->tmax;
        a_lo[nseg] = last_a_lo; a_hi[nseg] = last_a_hi; nseg++;
    }

    for (int s = 0; s < nseg; s++)
        ridge_refine(p, t_lo[s], t_hi[s], a_lo[s], a_hi[s],
                     &best, min_violation, evals);
    return best;
}

/* ---------- exhaustive reference solver ---------- */

static SolveResult solve_exhaustive(const BlendProblem* p, int nT, int na) {
    if (!problem_valid(p)) return bad_result(SOLVE_BAD_INPUT);
    SolveResult r;
    memset(&r, 0, sizeof(r));
    r.cand.T = -1.0f; r.cand.alpha = -1.0f;
    r.min_violation = INFINITY;

    Candidate best_feasible; memset(&best_feasible, 0, sizeof(best_feasible));
    best_feasible.cost = INFINITY;
    Candidate least_bad;   memset(&least_bad, 0, sizeof(least_bad));

    for (int i = 0; i < nT; i++) {
        float T = p->tmin + (p->tmax - p->tmin) * i / (nT - 1);
        for (int j = 0; j < na; j++) {
            float a = p->amin + (p->amax - p->amin) * j / (na - 1);
            Candidate c = evaluate(p, T, a);
            r.evals++;
            if (c.feasible) {
                if (c.cost < best_feasible.cost) best_feasible = c;
            } else {
                float v = violation_of(p, &c);
                if (v < r.min_violation) { r.min_violation = v; least_bad = c; }
            }
        }
    }

    /* polish the grid winner along constraint ridges so the reference is not
     * limited by grid quantization; grid + polish, not grid alone */
    if (best_feasible.cost != INFINITY) {
        Candidate polished = ridge_solve(p, 24, &r.min_violation, &r.evals);
        if (polished.cost < best_feasible.cost) best_feasible = polished;
    }

    if (best_feasible.cost != INFINITY) {
        r.status = SOLVE_OK;
        r.cand = best_feasible;
    } else {
        r.status = SOLVE_INFEASIBLE;
        r.cand = least_bad;      /* diagnostic only; T kept as sentinel below */
        r.cand.T = -1.0f;
        r.cand.alpha = -1.0f;
    }
    return r;
}

/* ---------- fast solver: multi-start coordinate pattern search ---------- */

static Candidate coord_search(const BlendProblem* p, float T0, float a0, int* evals) {
    float lo[2] = { p->tmin, p->amin };
    float hi[2] = { p->tmax, p->amax };
    float x[2] = { clampf(T0, lo[0], hi[0]), clampf(a0, lo[1], hi[1]) };
    float base_step[2] = { 0.1f * (hi[0] - lo[0]), 0.1f * (hi[1] - lo[1]) };
    float step[2] = { base_step[0], base_step[1] };

    Candidate cur = evaluate(p, x[0], x[1]);
    (*evals)++;
    float cur_pen = penalized(p, &cur);

    /* restarts re-open large diagonal steps so the search can keep sliding
     * along a curved constraint ridge after a fine-grained stall */
    for (int restart = 0; restart < 3; restart++) {
        if (restart > 0) { step[0] = base_step[0]; step[1] = base_step[1]; }
        int moved_this_phase = 0;
        for (int iter = 0; iter < 300; iter++) {
            int improved = 0;
            for (int k = 0; k < 2; k++) {
                for (int sgn = -1; sgn <= 1; sgn += 2) {
                    float y[2] = { x[0], x[1] };
                    y[k] = clampf(x[k] + sgn * step[k], lo[k], hi[k]);
                    if (fabsf(y[k] - x[k]) < 1e-12f) continue;
                    Candidate c = evaluate(p, y[0], y[1]);
                    (*evals)++;
                    float pen = penalized(p, &c);
                    if (pen < cur_pen - 1e-10f) {
                        cur = c; cur_pen = pen;
                        x[0] = y[0]; x[1] = y[1];
                        improved = 1; moved_this_phase = 1;
                        break;
                    }
                }
            }
            /* diagonal moves let the search slide along curved constraint ridges */
            if (!improved) {
                for (int s0 = -1; s0 <= 1 && !improved; s0 += 2)
                    for (int s1 = -1; s1 <= 1; s1 += 2) {
                        float y[2];
                        y[0] = clampf(x[0] + s0 * step[0], lo[0], hi[0]);
                        y[1] = clampf(x[1] + s1 * step[1], lo[1], hi[1]);
                        if (fabsf(y[0] - x[0]) < 1e-12f && fabsf(y[1] - x[1]) < 1e-12f) continue;
                        Candidate c = evaluate(p, y[0], y[1]);
                        (*evals)++;
                        float pen = penalized(p, &c);
                        if (pen < cur_pen - 1e-10f) {
                            cur = c; cur_pen = pen;
                            x[0] = y[0]; x[1] = y[1];
                            improved = 1; moved_this_phase = 1;
                            break;
                        }
                    }
            }
            if (!improved) {
                step[0] *= 0.5f; step[1] *= 0.5f;
                if (step[0] < 1e-6f * (hi[0] - lo[0]) &&
                    step[1] < 1e-6f * (hi[1] - lo[1])) break;
            }
        }
        if (!moved_this_phase) break;   /* a phase that never moved is a true local optimum */
    }
    return cur;
}

static SolveResult solve_fast(const BlendProblem* p) {
    if (!problem_valid(p)) return bad_result(SOLVE_BAD_INPUT);
    SolveResult r;
    memset(&r, 0, sizeof(r));
    r.cand.T = -1.0f; r.cand.alpha = -1.0f;
    r.min_violation = INFINITY;

    /* coarse 6x6 grid picks the multi-start seeds */
    float seed_T[6], seed_a[6];
    for (int i = 0; i < 6; i++) {
        seed_T[i] = p->tmin + (p->tmax - p->tmin) * i / 5.0f;
        seed_a[i] = p->amin + (p->amax - p->amin) * i / 5.0f;
    }
    Candidate grid_best; memset(&grid_best, 0, sizeof(grid_best));
    grid_best.cost = INFINITY;
    Candidate grid_least; memset(&grid_least, 0, sizeof(grid_least));
    float least_v = INFINITY;
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 6; j++) {
            Candidate c = evaluate(p, seed_T[i], seed_a[j]);
            r.evals++;
            if (c.feasible) {
                if (c.cost < grid_best.cost) grid_best = c;
            } else {
                float v = violation_of(p, &c);
                if (v < least_v) { least_v = v; grid_least = c; }
            }
        }

    Candidate best_feasible; memset(&best_feasible, 0, sizeof(best_feasible));
    best_feasible.cost = INFINITY;
    Candidate least_bad = grid_least;
    if (least_v != INFINITY) r.min_violation = least_v;

    float starts[4][2];
    int ns = 0;
    if (grid_best.cost != INFINITY) {
        starts[ns][0] = grid_best.T; starts[ns][1] = grid_best.alpha; ns++;
    }
    if (least_v != INFINITY) {
        starts[ns][0] = grid_least.T; starts[ns][1] = grid_least.alpha; ns++;
    }
    starts[ns][0] = 0.5f * (p->tmin + p->tmax);
    starts[ns][1] = 0.5f * (p->amin + p->amax); ns++;
    starts[ns][0] = p->tmin; starts[ns][1] = p->amin; ns++;

    for (int s = 0; s < ns; s++) {
        Candidate c = coord_search(p, starts[s][0], starts[s][1], &r.evals);
        if (c.feasible) {
            if (c.cost < best_feasible.cost) best_feasible = c;
        } else {
            float v = violation_of(p, &c);
            if (v < r.min_violation) { r.min_violation = v; least_bad = c; }
        }
    }

    /* ridge-aware sweep only when the switch-instant limit is active:
     * that is the constraint which bends the feasible set into a narrow
     * sliver where pattern search stalls. Without it, the box-plus-two-limits
     * geometry is wide enough for the multi-start search. */
    if (p->start_acc_max > 0.0f) {
        float ridge_violation = INFINITY;
        Candidate ridge = ridge_solve(p, 6, &ridge_violation, &r.evals);
        if (ridge.cost < best_feasible.cost) best_feasible = ridge;
        if (ridge_violation < r.min_violation) r.min_violation = ridge_violation;
    }

    if (best_feasible.cost != INFINITY) {
        r.status = SOLVE_OK;
        r.cand = best_feasible;
    } else {
        r.status = SOLVE_INFEASIBLE;
        r.cand = least_bad;
        r.cand.T = -1.0f;       /* sentinel: no answer, caller must not use it */
        r.cand.alpha = -1.0f;
    }
    return r;
}

/* ---------- self test ---------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static VelProfile profile_const(float vx, float vy) {
    VelProfile p;
    memset(&p, 0, sizeof(p));
    p.n = 1; p.t[0] = 0.0f;
    p.v[0].x = vx; p.v[0].y = vy;
    return p;
}

static void report(const char* tag, const SolveResult* r) {
    if (r->status == SOLVE_OK) {
        const Candidate* c = &r->cand;
        printf("[%s] status=OK T=%.4f alpha=%.4f J=%.5f err=%.5f acc=%.3f start_acc=%.3f evals=%d\n",
               tag, c->T, c->alpha, c->cost, c->err, c->acc_peak, c->start_acc, r->evals);
    } else if (r->status == SOLVE_INFEASIBLE) {
        printf("[%s] status=INFEASIBLE min_violation=%.5f evals=%d (no fallback)\n",
               tag, r->min_violation, r->evals);
    } else {
        printf("[%s] status=BAD_INPUT evals=%d\n", tag, r->evals);
    }
}

int main(void) {
    /* easing shape: fixed endpoints and strictly monotonic */
    check("ease_endpoints", fabsf(ease_weight(0.0f, 1.7f)) < 1e-7f &&
                            fabsf(ease_weight(1.0f, 1.7f) - 1.0f) < 1e-7f);
    int mono = 1;
    for (int i = 1; i <= 100; i++)
        if (!(ease_weight(i / 100.0f, 0.6f) >= ease_weight((i - 1) / 100.0f, 0.6f))) mono = 0;
    check("ease_monotonic", mono);

    /* blended velocity equals vA at fade start and vB at fade end */
    BlendProblem base;
    memset(&base, 0, sizeof(base));
    base.va = profile_const(1.0f, 0.0f);
    base.vb = profile_const(2.5f, 0.0f);
    base.tmin = 0.05f; base.tmax = 0.5f;
    base.amin = 0.5f;  base.amax = 3.0f;
    base.emax = 0.2f;  base.acc_max = 80.0f;
    base.w_dur = base.w_err = base.w_acc = 1.0f;
    {
        float T = 0.3f, alpha = 1.5f;
        Candidate c = evaluate(&base, T, alpha);
        /* endpoint behavior follows directly from w(0)=0, w(1)=1 */
        check("blend_endpoint_weight",
              fabsf(ease_weight(0.0f, alpha)) < 1e-7f &&
              fabsf(ease_weight(1.0f, alpha) - 1.0f) < 1e-7f);

        /* analytic integral for constant velocities:
         * e = |dv| * T * alpha / (alpha + 1) */
        float analytic = 1.5f * T * alpha / (alpha + 1.0f);
        printf("integral_simpson=%.6f integral_analytic=%.6f\n", c.err, analytic);
        check("integral_matches_analytic", fabsf(c.err - analytic) < 1e-5f);

        V2 et = integrate_error_trap(&base, T, alpha, CHECK_N);
        float etrap = vlen2(et);
        printf("integral_trap_fine=%.6f quad_discrepancy=%.2e bound=%.0e\n",
               etrap, fabsf(etrap - c.err), QUAD_BOUND);
        check("quadrature_crosscheck", fabsf(etrap - c.err) < QUAD_BOUND);
    }

    /* fast motion: large velocity difference, tight error budget */
    BlendProblem fast = base;
    fast.va = profile_const(0.0f, 0.0f);
    fast.vb = profile_const(4.0f, 0.0f);
    fast.emax = 0.15f;
    fast.acc_max = 60.0f;

    SolveResult fe = solve_exhaustive(&fast, 80, 80);
    SolveResult ff = solve_fast(&fast);
    report("fast/exhaustive", &fe);
    report("fast/fast", &ff);
    check("fast_both_ok", fe.status == SOLVE_OK && ff.status == SOLVE_OK);
    check("fast_feasible", ff.cand.err <= fast.emax + 1e-6f &&
                           ff.cand.acc_peak <= fast.acc_max + 1e-6f);
    float cost_dev = fabsf(ff.cand.cost - fe.cand.cost);
    float err_dev  = fabsf(ff.cand.err  - fe.cand.err);
    printf("fast cost_dev=%.6f bound=%.4f ok=%d\n", cost_dev, COST_BOUND, cost_dev <= COST_BOUND);
    printf("fast err_dev=%.6f bound=%.4f ok=%d\n", err_dev, ERR_BOUND, err_dev <= ERR_BOUND);
    check("fast_cost_bound", cost_dev <= COST_BOUND);
    check("fast_err_bound", err_dev <= ERR_BOUND);
    check("fast_fewer_evals", ff.evals < fe.evals / 4);

    /* slow motion: small velocity difference, strict acceleration limit.
     * acc <= 4 forces T >= ~0.15s, while the error budget keeps T bounded,
     * so the optimum is an interior trade-off rather than a box corner. */
    BlendProblem slow = base;
    slow.va = profile_const(0.2f, 0.0f);
    slow.vb = profile_const(0.8f, 0.0f);
    slow.emax = 0.1f;
    slow.acc_max = 4.0f;

    SolveResult se = solve_exhaustive(&slow, 80, 80);
    SolveResult sf = solve_fast(&slow);
    report("slow/exhaustive", &se);
    report("slow/fast", &sf);
    check("slow_both_ok", se.status == SOLVE_OK && sf.status == SOLVE_OK);
    check("slow_feasible", sf.cand.err <= slow.emax + 1e-6f &&
                           sf.cand.acc_peak <= slow.acc_max + 1e-6f);
    float scost_dev = fabsf(sf.cand.cost - se.cand.cost);
    float serr_dev  = fabsf(sf.cand.err  - se.cand.err);
    printf("slow cost_dev=%.6f bound=%.4f ok=%d\n", scost_dev, COST_BOUND, scost_dev <= COST_BOUND);
    printf("slow err_dev=%.6f bound=%.4f ok=%d\n", serr_dev, ERR_BOUND, serr_dev <= ERR_BOUND);
    check("slow_cost_bound", scost_dev <= COST_BOUND);
    check("slow_err_bound", serr_dev <= ERR_BOUND);

    /* the optimizer must adapt instead of returning one shared setting */
    check("params_adapt_to_motion",
          fabsf(ff.cand.T - sf.cand.T) > 1e-3f ||
          fabsf(ff.cand.alpha - sf.cand.alpha) > 0.05f);

    /* the previously hard-coded default (T=0.2s, linear easing) fails the fast case */
    Candidate dflt = evaluate(&fast, 0.2f, 1.0f);
    printf("hardcoded_default feasible=%d err=%.5f emax=%.4f acc=%.3f amax=%.2f\n",
           dflt.feasible, dflt.err, fast.emax, dflt.acc_peak, fast.acc_max);
    check("hardcoded_default_fails", dflt.feasible == 0);

    /* pop-sensitive switch: a bound on the switch-instant acceleration makes
     * linear easing infeasible; the optimum must move into the alpha > 1 region
     * and trade T against alpha, i.e. genuinely use both decision variables. */
    BlendProblem pop = base;
    pop.va = profile_const(0.0f, 0.0f);
    pop.vb = profile_const(5.0f, 0.0f);
    pop.tmin = 0.05f; pop.tmax = 0.5f;
    pop.amin = 0.5f;  pop.amax = 4.0f;
    pop.emax = 0.5f;
    pop.acc_max = 80.0f;
    pop.start_acc_max = 15.0f;
    Candidate lin = evaluate(&pop, 0.3f, 1.0f);
    printf("linear_easing start_acc=%.3f start_limit=%.2f feasible=%d\n",
           lin.start_acc, pop.start_acc_max, lin.feasible);
    check("linear_easing_violates_start", lin.feasible == 0);
    SolveResult pe = solve_exhaustive(&pop, 120, 120);
    SolveResult pf = solve_fast(&pop);
    report("pop/exhaustive", &pe);
    report("pop/fast", &pf);
    check("pop_both_ok", pe.status == SOLVE_OK && pf.status == SOLVE_OK);
    check("pop_uses_easing", pf.cand.alpha > 1.05f);
    check("pop_start_respected", pf.cand.start_acc <= pop.start_acc_max + 1e-6f);
    check("pop_all_respected", pf.cand.err <= pop.emax + 1e-6f &&
                               pf.cand.acc_peak <= pop.acc_max + 1e-6f);
    float pcost_dev = fabsf(pf.cand.cost - pe.cand.cost);
    float perr_dev  = fabsf(pf.cand.err  - pe.cand.err);
    printf("pop cost_dev=%.6f bound=%.4f ok=%d\n", pcost_dev, COST_BOUND, pcost_dev <= COST_BOUND);
    printf("pop err_dev=%.6f bound=%.4f ok=%d\n", perr_dev, ERR_BOUND, perr_dev <= ERR_BOUND);
    check("pop_cost_bound", pcost_dev <= COST_BOUND);
    check("pop_err_bound", perr_dev <= ERR_BOUND);
    check("pop_fewer_evals", pf.evals < pe.evals / 4);
    /* the same clip without the start constraint returns to near-linear easing:
     * alpha is driven by the constraint, not by a blind preference. */
    BlendProblem pop_free = pop;
    pop_free.start_acc_max = 0.0f;
    SolveResult pfree = solve_fast(&pop_free);
    report("pop/free", &pfree);
    check("easing_driven_by_constraint", pfree.status == SOLVE_OK &&
          pfree.cand.alpha < pf.cand.alpha - 0.15f);

    /* infeasible problem: error budget below the physical minimum.
     * smallest achievable error is at Tmin with the smallest alpha. */
    BlendProblem impossible = fast;
    impossible.emax = 0.01f;
    SolveResult ie = solve_exhaustive(&impossible, 80, 80);
    SolveResult iff = solve_fast(&impossible);
    report("impossible/exhaustive", &ie);
    report("impossible/fast", &iff);
    check("infeasible_reported",
          ie.status == SOLVE_INFEASIBLE && iff.status == SOLVE_INFEASIBLE);
    check("no_default_fallback",
          ie.cand.T < 0.0f && iff.cand.T < 0.0f &&
          ie.min_violation > 0.0f && iff.min_violation > 0.0f);

    /* invalid problem descriptions are rejected, not silently clamped */
    BlendProblem bad1 = fast; bad1.tmax = 0.01f;          /* tmax < tmin */
    BlendProblem bad2 = fast; bad2.amin = -1.0f;          /* non-positive exponent */
    BlendProblem bad3 = fast; bad3.emax = 0.0f;           /* non-positive budget */
    check("bad_input_duration", solve_fast(&bad1).status == SOLVE_BAD_INPUT);
    check("bad_input_alpha",    solve_fast(&bad2).status == SOLVE_BAD_INPUT);
    check("bad_input_budget",   solve_fast(&bad3).status == SOLVE_BAD_INPUT);

    /* time-varying (piecewise) velocity profile: both solvers must agree on status */
    VelProfile va_step;
    memset(&va_step, 0, sizeof(va_step));
    va_step.n = 2; va_step.t[0] = 0.0f; va_step.t[1] = 0.15f;
    va_step.v[0].x = 1.0f; va_step.v[0].y = 0.0f;
    va_step.v[1].x = 0.3f; va_step.v[1].y = 0.0f;
    BlendProblem varying = base;
    varying.va = va_step;
    varying.vb = profile_const(2.0f, 0.0f);
    SolveResult ve = solve_exhaustive(&varying, 64, 64);
    SolveResult vf = solve_fast(&varying);
    report("varying/exhaustive", &ve);
    report("varying/fast", &vf);
    check("varying_status_agrees", ve.status == vf.status);
    if (ve.status == SOLVE_OK) {
        float vcd = fabsf(vf.cand.cost - ve.cand.cost);
        float ved = fabsf(vf.cand.err - ve.cand.err);
        check("varying_cost_bound", vcd <= COST_BOUND);
        check("varying_err_bound", ved <= ERR_BOUND);
    }

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
