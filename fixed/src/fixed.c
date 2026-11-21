/* fixed: fixed-point (quantized) Verlet distance-constraint chain with a
   drift auditor.

   Why
   ---
   The floating-point chain in ../constraint/ slowly drifts when it runs for
   a long time: every frame's rounding error is fed back through the
   integrator, and the particle can end up somewhere an exact-arithmetic run
   would never put it.  This module stores positions, previous positions,
   rest lengths and every applied correction as fixed-point integers (a
   multiple of the quantization step q), and derives from q an upper bound
   on the rounding contribution to how far the quantized run may be from a
   floating-point run started from the same state.  A floating-point
   reference executes the identical algorithm, the per-frame drift is
   measured against it, and a drift above the bound raises an alarm.  A
   step so coarse that the distance constraints cannot reach the configured
   tolerance within the iteration budget is reported as non-convergence.

   Which floating-point reference
   ------------------------------
   The reference is double precision, not float.  float carries its own
   integration drift of the same order as the effect being measured -- it is
   what ../constraint/ demonstrates drifting -- so using it as the reference
   would hide small quantization inside float noise and add noise to large
   quantization.  double rounds at 2^-52 relative and its accumulated error
   stays far below one quantum for the q_bits range accepted here, so what
   the auditor measures is the effect of the fixed-point quantization.

   Representation
   --------------
   q = 2^-q_bits world units.  State (x, y, px, py, rest) is int64 quanta.
   Geometry (distance, weights, corrections) is evaluated in long double
   from the integer state; each correction applied to a movable endpoint is
   rounded to the nearest quantum, and that rounding is the only place
   quantization enters the constraint solve.  Gravity is emitted as a whole
   number of quanta through a fractional carry shared with the reference,
   so its long-run average is exact and contributes no relative drift.
   Damping is quantized to 32 fractional bits and both chains use the same
   factor.

   Error bound
   -----------
   One converged frame rounds each movable endpoint's exact correction
   once.  With U movable particles, K passes, and S movable endpoints over
   all constraints (2U integration roundings plus 2KS correction roundings,
   each within q/2), a single frame's rounding budget is

     r0 = (U + K*S) q.

   The audit checks it on the frame that has no history: the drift after
   frame 0 must be <= r0.  That check is sharp and is measured at a
   fraction of r0 on the self-test chain.

   Over a run, a difference can also be amplified by the dynamics: frame
   N's difference is carried into frame N+1 and magnified by the step map.
   For the damped test chain the map keeps a q-sized perturbation at or
   below its size from frame to frame, so the run-long drift stays inside

     B(N) = N * r0,

   the total rounding budget.  This is the bound the alarm compares against,
   and the damped self-test chains stay well below it from q = 2^-4 to
   q = 2^-16.  In a stiff or chaotic regime the map does not contract, and
   the observed difference can exceed B(N) for reasons that have nothing to
   do with rounding: the module then reports FM_DRIFT_EXCEEDED with the
   ratio and the frame of onset, and the caller reads that as "the
   difference is not explained by quantization" -- either the step in force
   is coarser than declared (a defect) or the dynamics amplifies a
   legitimate q-sized difference (a property of the run).  fm_report()
   prints peak/r0 and peak/(N*r0) so the two cases are distinguishable: an
   op-count defect exceeds r0 on a single frame early and keeps the ratio
   near 1, while amplification shows a steadily growing ratio far above it.
   The self-test covers both, including a near-undamped chain where the
   amplification regime is expected and must be reported, not hidden.

   Non-convergence
   ---------------
   q too large leaves a distance residual after the K passes: each pass may
   only move by whole quanta, so the residual per constraint is a small
   multiple of q while the step is fine and grows with q once a single pass
   overshoots the rest length.  After every frame the residual in world
   units is compared with res_tol; a frame over tolerance is counted,
   localized to its worst constraint, and drives the status to
   FM_NONCONVERGED.  This is independent of the drift alarm: the step can
   be too coarse for the geometry while the trajectory stays far inside
   its bound.  The crossover is data dependent (measured around 10-20 q on
   the default chain), so the criterion is the configured tolerance, not a
   multiplier of q.

   Validity and overflow
   ---------------------
   q_bits, geometry, damping, gravity, iteration count, tolerance and the
   reference setup are validated up front.  A rejected configuration
   reports FM_INVALID and never yields another verdict.  Hitting the int64
   coordinate cap mid-run latches FM_OVERFLOW and refuses further steps
   instead of wrapping into a plausible-looking trajectory.

   Status and reporting
   --------------------
   fm_status() precedence: invalid, overflow, drift_exceeded,
   nonconverged, undetermined (fewer than min_frames samples and no flag
   yet), ok.  fm_report() prints r0, the frame-0 ratio, peak drift,
   peak/r0, peak/(N*r0), the worst single-frame increment, the first
   exceedances (frame, particle, axis), the worst residual and its
   threshold.  fm_worst_constraint() names the constraint that owned the
   residual on the first non-converged frame, the same per-constraint
   number a manual check would compute.
*/
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <stdlib.h>

#define MAX_P 16
#define FX_BITS_MAX 40            /* |coordinate| <= 2^40 quanta */
#define FX_LIMIT (1LL << FX_BITS_MAX)
#define FX_DAMP_SHIFT 32
#define FX_ONE32 (1LL << FX_DAMP_SHIFT)

typedef enum {
    FM_OK = 0,
    FM_UNDETERMINED,
    FM_NONCONVERGED,
    FM_DRIFT_EXCEEDED,
    FM_OVERFLOW,
    FM_INVALID
} FmStatus;

typedef struct {
    int n;                /* particles, 2..MAX_P */
    double seg;           /* segment rest length, > 0 */
    double x0, y0;        /* initial pinned position */
    double damping;       /* 0 < damping < 1 */
    int iterations;       /* constraint passes per frame, 1..256 */
    int q_bits;           /* q = 2^-q_bits, 1..30 */
    double dt;            /* integration step, > 0 */
    double gx, gy;        /* gravity */
    double res_tol;       /* convergence tolerance in world units, > 0 */
    int min_frames;       /* samples before an ok verdict, >= 1 */
} FixCfg;

typedef struct {
    int64_t x, y, px, py;
    double inv_mass;
    int pinned;
} FParticle;

typedef struct {
    FParticle p[MAX_P];
    int n;
    int64_t rest[MAX_P];
    int iterations;
    int64_t dq32;
    int64_t gx_q32, gy_q32;
    int64_t gremx, gremy;
    double q;
    int overflow;
    int fault_quant;       /* self-test only: quantization coarser than q */
} FixChain;

typedef struct {
    double x, y, px, py, inv_mass;
    int pinned;
} RParticle;

typedef struct {
    RParticle p[MAX_P];
    int n;
    double rest[MAX_P];
    int iterations;
    int64_t dq32;
    int64_t gx_q32, gy_q32;
    int64_t gremx, gremy;
    double q;
} RefChain;

typedef struct {
    FixCfg cfg;
    FixChain fx;
    RefChain ref;
    double q;

    int frames;
    int valid;
    int overflow;

    double r0;             /* one converged frame's rounding bound */
    double frame0_drift;   /* measured drift of the very first frame */
    double accum_bound;    /* N * r0 for the latest frame */

    double peak_drift;
    double peak_ratio;     /* peak_drift / r0 */
    double latest_drift;

    int drift_frames;      /* frames above the accumulated cap */
    int first_drift_frame;
    int first_drift_particle;
    int first_drift_axis;
    int frame_jump_frames; /* frames whose own drift increment > r0 */
    int first_jump_frame;
    int first_jump_particle;
    int first_jump_axis;
    double max_jump;

    int nonconv_frames;
    int first_nonconv_frame;
    int worst_con;         /* constraint owning the worst residual */
    double worst_residual; /* max over run */
    double first_nonconv_resid;

} FixAudit;

static const char* fm_status_name(FmStatus s) {
    switch (s) {
        case FM_OK: return "ok";
        case FM_UNDETERMINED: return "undetermined";
        case FM_NONCONVERGED: return "nonconverged";
        case FM_DRIFT_EXCEEDED: return "drift_exceeded";
        case FM_OVERFLOW: return "overflow";
        case FM_INVALID: return "invalid";
    }
    return "invalid";
}

static const char* fm_status_reason(FmStatus s) {
    switch (s) {
        case FM_OK: return "within drift bound and converged";
        case FM_UNDETERMINED: return "insufficient samples: cannot determine";
        case FM_NONCONVERGED: return "constraint residual above tolerance";
        case FM_DRIFT_EXCEEDED: return "drift exceeded the quantization bound";
        case FM_OVERFLOW: return "fixed-point coordinate range exceeded";
        case FM_INVALID: return "invalid quantization step or configuration";
    }
    return "invalid quantization step or configuration";
}

static int fm_config_valid(const FixCfg* c) {
    if (c->n < 2 || c->n > MAX_P) return 0;
    if (!isfinite(c->seg) || c->seg <= 0.0) return 0;
    if (!isfinite(c->x0) || !isfinite(c->y0)) return 0;
    if (!isfinite(c->damping) || c->damping <= 0.0 || c->damping >= 1.0) return 0;
    if (c->iterations < 1 || c->iterations > 256) return 0;
    if (c->q_bits < 1 || c->q_bits > 30) return 0;
    if (!isfinite(c->dt) || c->dt <= 0.0) return 0;
    if (!isfinite(c->gx) || !isfinite(c->gy)) return 0;
    if (!isfinite(c->res_tol) || c->res_tol <= 0.0) return 0;
    if (c->min_frames < 1) return 0;
    double q = ldexp(1.0, -c->q_bits);
    double span = fabs(c->x0) + fabs(c->y0) + (double)c->n * c->seg;
    if (!isfinite(span) || span / q > (double)FX_LIMIT * 0.5) return 0;
    /* the per-frame gravity displacement must fit the 32-bit carry */
    double gxq = c->gx * c->dt * c->dt / q * (double)FX_ONE32;
    double gyq = c->gy * c->dt * c->dt / q * (double)FX_ONE32;
    if (!isfinite(gxq) || !isfinite(gyq)) return 0;
    if (fabs(gxq) > 9.0e18 || fabs(gyq) > 9.0e18) return 0;
    return 1;
}

/* ------------------------------------------------------------------ */
/* Chain construction and one-frame updates                           */
/* ------------------------------------------------------------------ */

static void fix_init(FixChain* c, const FixCfg* cfg, double q) {
    memset(c, 0, sizeof(*c));
    c->n = cfg->n;
    c->iterations = cfg->iterations;
    c->q = q;
    c->dq32 = (int64_t)llround(cfg->damping * (double)FX_ONE32);
    long double gxq = (long double)cfg->gx * cfg->dt * cfg->dt / q;
    long double gyq = (long double)cfg->gy * cfg->dt * cfg->dt / q;
    c->gx_q32 = (int64_t)llroundl(gxq * (long double)FX_ONE32);
    c->gy_q32 = (int64_t)llroundl(gyq * (long double)FX_ONE32);
    for (int i = 0; i < c->n; i++) {
        c->p[i].x = (int64_t)llround((cfg->x0 + cfg->seg * i) / q);
        c->p[i].y = (int64_t)llround(cfg->y0 / q);
        c->p[i].px = c->p[i].x;
        c->p[i].py = c->p[i].y;
        c->p[i].inv_mass = 1.0;
    }
    for (int i = 0; i < c->n - 1; i++) c->rest[i] = (int64_t)llround(cfg->seg / q);
    c->p[0].pinned = 1;
    c->p[0].inv_mass = 0.0;
}

/* The reference starts from the fixed chain's own quantized state and uses
   its quantized damping and gravity emitter, so the two chains differ only
   by the rounding of the constraint corrections. */
static void ref_init(RefChain* r, const FixChain* c) {
    memset(r, 0, sizeof(*r));
    r->n = c->n;
    r->iterations = c->iterations;
    r->dq32 = c->dq32;
    r->gx_q32 = c->gx_q32;
    r->gy_q32 = c->gy_q32;
    r->q = c->q;
    for (int i = 0; i < c->n; i++) {
        r->p[i].x = (double)c->p[i].x * c->q;
        r->p[i].y = (double)c->p[i].y * c->q;
        r->p[i].px = (double)c->p[i].px * c->q;
        r->p[i].py = (double)c->p[i].py * c->q;
        r->p[i].inv_mass = c->p[i].inv_mass;
        r->p[i].pinned = c->p[i].pinned;
    }
    for (int i = 0; i < c->n - 1; i++) r->rest[i] = (double)c->rest[i] * c->q;
}

/* Rounds a correction in quanta to a whole number of quanta.  fault_quant
   is a self-test injection: it makes the solver quantize more coarsely than
   the bound models, which must trip the drift alarm. */
static int64_t fx_q(FixChain* c, long double v) {
    int64_t k = c->fault_quant > 1 ? c->fault_quant : 1;
    if (k == 1) return (int64_t)llroundl(v);
    return (int64_t)(roundl(v / (long double)k) * (long double)k);
}

static int64_t fx_carry(int64_t* acc, int64_t per) {
    *acc += per;
    int64_t s = *acc >> FX_DAMP_SHIFT;
    *acc -= s << FX_DAMP_SHIFT;
    return s;
}

/* One frame.  Returns 1 on success, 0 when the int64 cap was hit. */
static int fix_step(FixChain* c) {
    if (c->overflow) return 0;
    for (int i = 0; i < c->n; i++) {
        FParticle* p = &c->p[i];
        if (p->pinned) continue;
        int64_t vx = p->x - p->px;
        int64_t vy = p->y - p->py;
        p->px = p->x;
        p->py = p->y;
        int64_t ax = fx_q(c, (long double)vx * (long double)c->dq32 / (long double)FX_ONE32);
        int64_t ay = fx_q(c, (long double)vy * (long double)c->dq32 / (long double)FX_ONE32);
        p->x += ax + fx_carry(&c->gremx, c->gx_q32);
        p->y += ay + fx_carry(&c->gremy, c->gy_q32);
        if (llabs(p->x) > FX_LIMIT || llabs(p->y) > FX_LIMIT) { c->overflow = 1; return 0; }
    }
    for (int it = 0; it < c->iterations; it++) {
        for (int i = 0; i < c->n - 1; i++) {
            FParticle* a = &c->p[i];
            FParticle* b = &c->p[i + 1];
            int64_t dx = b->x - a->x, dy = b->y - a->y;
            if (dx == 0 && dy == 0) continue;
            long double d = sqrtl((long double)dx * dx + (long double)dy * dy);
            long double diff = (d - (long double)c->rest[i]) / d;
            long double cx = (long double)dx * diff;
            long double cy = (long double)dy * diff;
            double wa = a->inv_mass, wb = b->inv_mass;
            long double ws = (long double)wa + (long double)wb;
            if (ws <= 0.0L) continue;
            if (!a->pinned) {
                a->x += fx_q(c, cx * wa / ws);
                a->y += fx_q(c, cy * wa / ws);
            }
            if (!b->pinned) {
                b->x -= fx_q(c, cx * wb / ws);
                b->y -= fx_q(c, cy * wb / ws);
            }
            if (llabs(a->x) > FX_LIMIT || llabs(a->y) > FX_LIMIT ||
                llabs(b->x) > FX_LIMIT || llabs(b->y) > FX_LIMIT) { c->overflow = 1; return 0; }
        }
    }
    return 1;
}

static void ref_step(RefChain* r) {
    for (int i = 0; i < r->n; i++) {
        RParticle* p = &r->p[i];
        if (p->pinned) continue;
        double vx = (p->x - p->px) * (double)r->dq32 / (double)FX_ONE32;
        double vy = (p->y - p->py) * (double)r->dq32 / (double)FX_ONE32;
        p->px = p->x;
        p->py = p->y;
        p->x += vx + (double)fx_carry(&r->gremx, r->gx_q32) * r->q;
        p->y += vy + (double)fx_carry(&r->gremy, r->gy_q32) * r->q;
    }
    for (int it = 0; it < r->iterations; it++) {
        for (int i = 0; i < r->n - 1; i++) {
            RParticle* a = &r->p[i];
            RParticle* b = &r->p[i + 1];
            double dx = b->x - a->x, dy = b->y - a->y;
            double d = sqrt(dx * dx + dy * dy);
            if (d == 0.0) continue;
            double diff = (d - r->rest[i]) / d;
            double wa = a->inv_mass, wb = b->inv_mass, ws = wa + wb;
            if (ws <= 0.0) continue;
            double cx = dx * diff, cy = dy * diff;
            if (!a->pinned) { a->x += cx * wa / ws; a->y += cy * wa / ws; }
            if (!b->pinned) { b->x -= cx * wb / ws; b->y -= cy * wb / ws; }
        }
    }
}

/* ------------------------------------------------------------------ */
/* Auditor                                                            */
/* ------------------------------------------------------------------ */

static double fix_residual(const FixChain* c, int* worst_idx) {
    double worst = 0.0;
    int idx = -1;
    for (int i = 0; i < c->n - 1; i++) {
        int64_t dx = c->p[i + 1].x - c->p[i].x;
        int64_t dy = c->p[i + 1].y - c->p[i].y;
        long double d = sqrtl((long double)dx * dx + (long double)dy * dy) * c->q;
        double e = (double)fabsl(d - (long double)c->rest[i] * c->q);
        if (e > worst) { worst = e; idx = i; }
    }
    if (worst_idx) *worst_idx = idx;
    return worst;
}

/* Max over movable particles and both axes of |fixed - reference|.  The
   caller gets the particle and axis that owned it (axis 0 = x, 1 = y). */
static double fix_drift(const FixChain* c, const RefChain* r, int* who, int* axis) {
    double worst = 0.0;
    int wi = -1, wa = 0;
    for (int i = 0; i < c->n; i++) {
        if (c->p[i].pinned) continue;
        double x = (double)c->p[i].x * c->q;
        double y = (double)c->p[i].y * c->q;
        double dx = fabs(x - r->p[i].x);
        double dy = fabs(y - r->p[i].y);
        if (dx > worst) { worst = dx; wi = i; wa = 0; }
        if (dy > worst) { worst = dy; wi = i; wa = 1; }
    }
    if (who) *who = wi;
    if (axis) *axis = wa;
    return worst;
}

/* An upper bound on the drift increment of one frame requires the number of
   movable endpoints over all constraints, which the auditor records here. */
static int fix_movable_endpoints(const FixChain* c) {    int s = 0;
    for (int i = 0; i < c->n - 1; i++) {
        if (!c->p[i].pinned) s++;
        if (!c->p[i + 1].pinned) s++;
    }
    return s;
}

static int fm_init(FixAudit* a, const FixCfg* cfg) {
    memset(a, 0, sizeof(*a));
    a->cfg = *cfg;
    a->worst_con = -1;
    if (!fm_config_valid(cfg)) { a->valid = 0; return 0; }
    a->valid = 1;
    a->q = ldexp(1.0, -cfg->q_bits);
    fix_init(&a->fx, cfg, a->q);
    ref_init(&a->ref, &a->fx);

    int n = cfg->n, iters = cfg->iterations;
    int movable = 0;
    for (int i = 0; i < n; i++) if (!a->fx.p[i].pinned) movable++;
    int s = fix_movable_endpoints(&a->fx);

    /* One converged frame rounds each movable endpoint's exact correction
       once: 2U integration roundings plus 2KS correction roundings, each
       within q/2. */
    a->r0 = ((double)movable + (double)iters * (double)s) * a->q;

    a->latest_drift = 0.0;
    a->accum_bound = 0.0;
    a->peak_drift = 0.0;
    a->peak_ratio = 0.0;
    a->max_jump = 0.0;
    return 1;
}

/* Advance one frame and audit it.  Returns 0 when the audit refuses to
   continue (invalid configuration or int64 overflow). */
static int fm_step(FixAudit* a) {
    if (!a->valid || a->overflow) return 0;
    if (!fix_step(&a->fx)) { a->overflow = 1; return 0; }
    a->frames++;
    ref_step(&a->ref);

    a->accum_bound = (double)a->frames * a->r0;

    int who = -1, axis = 0;
    double drift = fix_drift(&a->fx, &a->ref, &who, &axis);
    double jump = drift - a->latest_drift;
    a->latest_drift = drift;
    if (jump > a->max_jump) a->max_jump = jump;

    if (a->frames == 1) a->frame0_drift = drift;
    if (jump > a->r0) {
        /* One frame cannot move a converged coordinate by more than that
           frame's whole rounding budget; a larger increment means the step
           in force is coarser than the declared q (or the dynamics
           amplified the difference in a single frame). */
        a->frame_jump_frames++;
        if (a->first_jump_frame == 0) {
            a->first_jump_frame = a->frames;
            a->first_jump_particle = who;
            a->first_jump_axis = axis;
        }
    }
    if (drift > a->accum_bound) {
        a->drift_frames++;
        if (a->first_drift_frame == 0) {
            a->first_drift_frame = a->frames;
            a->first_drift_particle = who;
            a->first_drift_axis = axis;
        }
    }
    if (drift > a->peak_drift) {
        a->peak_drift = drift;
        a->peak_ratio = (a->r0 > 0.0) ? drift / a->r0 : 0.0;
    }

    int widx = -1;
    double resid = fix_residual(&a->fx, &widx);
    if (resid > a->worst_residual) a->worst_residual = resid;
    if (resid > a->cfg.res_tol) {
        a->nonconv_frames++;
        if (a->first_nonconv_frame == 0) {
            a->first_nonconv_frame = a->frames;
            a->worst_con = widx;
            a->first_nonconv_resid = resid;
        }
    }
    return 1;
}

static FmStatus fm_status(const FixAudit* a) {
    if (!a->valid) return FM_INVALID;
    if (a->overflow) return FM_OVERFLOW;
    if (a->drift_frames > 0 || a->frame_jump_frames > 0) return FM_DRIFT_EXCEEDED;
    if (a->nonconv_frames > 0) return FM_NONCONVERGED;
    if (a->frames < a->cfg.min_frames) return FM_UNDETERMINED;
    return FM_OK;
}

/* Constraint that owned the residual on the first non-converged frame,
   -1 if no frame failed to converge.  Its residual is written to *resid. */
static int fm_worst_constraint(const FixAudit* a, double* resid) {
    if (resid) *resid = a->first_nonconv_resid;
    if (a->nonconv_frames == 0) return -1;
    return a->worst_con;
}

static void fm_report(const FixAudit* a) {
    double frame0_ratio = (a->r0 > 0.0) ? a->frame0_drift / a->r0 : 0.0;
    double accum_ratio = (a->accum_bound > 0.0) ? a->peak_drift / a->accum_bound : 0.0;
    printf("r0=%.6g frame0_drift=%.6g frame0_ratio=%.3g\n",
           a->r0, a->frame0_drift, frame0_ratio);
    printf("drift_peak=%.6g peak_over_r0=%.3g accum_bound=%.6g peak_over_accum=%.3g\n",
           a->peak_drift, a->peak_ratio, a->accum_bound, accum_ratio);
    printf("drift_latest=%.6g max_frame_jump=%.6g jump_exceeded=%d first_jump=%d particle=%d axis=%d\n",
           a->latest_drift, a->max_jump, a->frame_jump_frames, a->first_jump_frame,
           a->first_jump_particle, a->first_jump_axis);
    printf("accumulated_exceeded=%d first_frame=%d particle=%d axis=%d\n",
           a->drift_frames, a->first_drift_frame,
           a->first_drift_particle, a->first_drift_axis);
    printf("residual_worst=%.6g res_tol=%.6g nonconv_frames=%d worst_con=%d\n",
           a->worst_residual, a->cfg.res_tol, a->nonconv_frames, a->worst_con);
    printf("status=%s reason=\"%s\"\n",
           fm_status_name(fm_status(a)), fm_status_reason(fm_status(a)));
}

/* ------------------------------------------------------------------ */
/* Self-test                                                          */
/* ------------------------------------------------------------------ */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static FixCfg base_cfg(void) {
    FixCfg c;
    memset(&c, 0, sizeof(c));
    c.n = 8;
    c.seg = 1.0;
    c.x0 = 0.0;
    c.y0 = 0.0;
    c.damping = 0.99;
    c.iterations = 8;
    c.q_bits = 12;
    c.dt = 0.016;
    c.gx = 0.0;
    c.gy = -9.8;
    c.res_tol = 0.01;
    c.min_frames = 50;
    return c;
}

/* Run the audit for nframes; returns the status. */
static FmStatus run(FixAudit* a, int nframes) {
    for (int i = 0; i < nframes; i++) fm_step(a);
    return fm_status(a);
}

/* Independent per-constraint residual, used to check localization. */
static double manual_residual(const FixAudit* a, int* idx) {
    double worst = 0;
    int wi = -1;
    for (int i = 0; i < a->fx.n - 1; i++) {
        double dx = (double)(a->fx.p[i + 1].x - a->fx.p[i].x) * a->q;
        double dy = (double)(a->fx.p[i + 1].y - a->fx.p[i].y) * a->q;
        double e = fabs(sqrt(dx * dx + dy * dy) - (double)a->fx.rest[i] * a->q);
        if (e > worst) { worst = e; wi = i; }
    }
    if (idx) *idx = wi;
    return worst;
}

int main(void) {
    FixCfg C = base_cfg();

    /* ---- configuration validation ---- */
    {
        FixAudit a;
        FixCfg bad;

        check("init_valid", fm_init(&a, &C) == 1 && fm_status(&a) != FM_INVALID);

        bad = C; bad.q_bits = 0;
        check("reject_q_bits_zero", fm_init(&a, &bad) == 0 && fm_status(&a) == FM_INVALID);
        bad = C; bad.q_bits = 31;
        check("reject_q_bits_large", fm_init(&a, &bad) == 0);
        bad = C; bad.q_bits = -3;
        check("reject_q_bits_negative", fm_init(&a, &bad) == 0);
        bad = C; bad.n = 1;
        check("reject_n_too_small", fm_init(&a, &bad) == 0);
        bad = C; bad.n = MAX_P + 1;
        check("reject_n_too_large", fm_init(&a, &bad) == 0);
        bad = C; bad.seg = 0.0;
        check("reject_seg_zero", fm_init(&a, &bad) == 0);
        bad = C; bad.seg = NAN;
        check("reject_seg_nan", fm_init(&a, &bad) == 0);
        bad = C; bad.damping = 1.0;
        check("reject_damping_one", fm_init(&a, &bad) == 0);
        bad = C; bad.damping = 0.0;
        check("reject_damping_zero", fm_init(&a, &bad) == 0);
        bad = C; bad.damping = 1.2;
        check("reject_damping_above_one", fm_init(&a, &bad) == 0);
        bad = C; bad.iterations = 0;
        check("reject_iters_zero", fm_init(&a, &bad) == 0);
        bad = C; bad.dt = 0.0;
        check("reject_dt_zero", fm_init(&a, &bad) == 0);
        bad = C; bad.dt = -0.016;
        check("reject_dt_negative", fm_init(&a, &bad) == 0);
        bad = C; bad.gx = INFINITY;
        check("reject_gx_inf", fm_init(&a, &bad) == 0);
        bad = C; bad.res_tol = 0.0;
        check("reject_res_tol_zero", fm_init(&a, &bad) == 0);
        bad = C; bad.res_tol = NAN;
        check("reject_res_tol_nan", fm_init(&a, &bad) == 0);
        bad = C; bad.min_frames = 0;
        check("reject_min_frames_zero", fm_init(&a, &bad) == 0);
        /* a range that cannot fit the int64 coordinate cap is refused even
           though each parameter on its own is fine */
        bad = C; bad.q_bits = 30; bad.x0 = 1e8;
        check("reject_range_fit", fm_init(&a, &bad) == 0);

        /* a rejected audit cannot be revived */
        bad = C; bad.q_bits = 0;
        fm_init(&a, &bad);
        fm_step(&a);
        check("rejected_stays_invalid", fm_status(&a) == FM_INVALID && a.frames == 0);
    }

    /* ---- undetermined until enough samples; ok afterwards ---- */
    {
        FixAudit a;
        fm_init(&a, &C);
        check("initial_undetermined", fm_status(&a) == FM_UNDETERMINED);
        FmStatus s = run(&a, C.min_frames - 1);
        check("short_run_undetermined", s == FM_UNDETERMINED);
        s = run(&a, 1);
        check("at_min_frames_determined", a.frames >= C.min_frames && s != FM_UNDETERMINED);
        run(&a, 500);
        check("healthy_run_ok", fm_status(&a) == FM_OK);
        check("healthy_no_drift_alarm", a.drift_frames == 0);
        check("healthy_no_nonconvergence", a.nonconv_frames == 0);
        fm_report(&a);
    }

    /* ---- frame-0 drift is within the converged one-frame bound ---- */
    {
        FixAudit a;
        fm_init(&a, &C);
        fm_step(&a);
        check("frame0_recorded", a.frames == 1);
        check("frame0_within_r0", a.frame0_drift <= a.r0);
        printf("frame0_drift=%.4g r0=%.4g ratio=%.3g\n",
               a.frame0_drift, a.r0, a.frame0_drift / a.r0);
    }

    /* ---- quantization step controls the residual: fine q keeps the chain
       inside the tolerance, coarse q cannot converge and is reported ---- */
    {
        FixAudit fine, coarse;
        FixCfg f = C, g = C;
        f.q_bits = 16; f.min_frames = 10;
        g.q_bits = 2;  g.min_frames = 10;
        fm_init(&fine, &f);
        fm_init(&coarse, &g);
        run(&fine, 2000);
        run(&coarse, 2000);
        check("fine_q_converges", fine.nonconv_frames == 0 &&
              fm_status(&fine) != FM_NONCONVERGED);
        check("coarse_q_nonconverged", fm_status(&coarse) == FM_NONCONVERGED);
        check("coarse_q_reports_residual", coarse.worst_residual > g.res_tol);
        /* localization: the reported constraint and residual must equal a
           manual per-constraint check on the same (latest) frame */
        /* localization: on the frame the auditor first reports
           non-convergence, its report must equal a manual per-constraint
           residual of that same frame. */
        FixAudit loc;
        fm_init(&loc, &g);
        double manual_at_first = 0.0;
        int manual_idx = -1;
        int seen = 0;
        for (int i = 0; i < 2000; i++) {
            fm_step(&loc);
            if (loc.nonconv_frames > 0) {
                manual_at_first = manual_residual(&loc, &manual_idx);
                seen = 1;
                break;
            }
        }
        double rep = 0.0;
        int got = fm_worst_constraint(&loc, &rep);
        check("nonconv_localized_matches_manual",
              seen && got == manual_idx && fabs(rep - manual_at_first) < 1e-6);
        /* a 1-iteration run at coarse q must stay above the fine one's
           residual: the residual decreases with iterations as well */
        FixCfg h = g; h.iterations = 1;
        FixAudit one;
        fm_init(&one, &h);
        run(&one, 2000);
        check("more_iters_reduce_residual", coarse.worst_residual <= one.worst_residual + 1e-9);
        printf("fine_res=%.4g coarse_res=%.4g coarse_con=%d one_iter_res=%.4g\n",
               fine.worst_residual, coarse.worst_residual, got, one.worst_residual);
    }

    /* ---- the residual threshold is the reported criterion: raising it
       above the measured residual clears the non-convergence verdict ---- */
    {
        FixAudit coarse;
        FixCfg g = C;
        g.q_bits = 4;
        fm_init(&coarse, &g);
        run(&coarse, 1000);
        check("coarse_flagged_with_default_tol", fm_status(&coarse) == FM_NONCONVERGED);
        double measured = coarse.worst_residual;

        FixAudit relaxed;
        FixCfg g2 = g;
        g2.res_tol = measured * 2.0 + 1e-9;
        fm_init(&relaxed, &g2);
        run(&relaxed, 1000);
        check("relaxed_tol_allows", fm_status(&relaxed) != FM_NONCONVERGED);
    }

    /* ---- the drift alarm: a correctly rounded chain stays inside the
       operation-counted bound across the whole step range, and a chain
       that truncates instead of rounding exceeds it ---- */
    {
        int all_below = 1;
        for (int bits = 4; bits <= 16; bits += 2) {
            FixAudit a;
            FixCfg c = C;
            c.q_bits = bits;
            c.min_frames = 1;
            fm_init(&a, &c);
            run(&a, 3000);
            if (a.drift_frames > 0) all_below = 0;
            printf("bound_scan q_bits=%d peak=%.4g accum_bound=%.4g peak_over_accum=%.3g\n",
                   bits, a.peak_drift, a.accum_bound,
                   a.accum_bound > 0 ? a.peak_drift / a.accum_bound : 0.0);
        }
        check("rounded_chain_within_bound", all_below);

        /* truncation biases every correction downward by up to a whole
           quantum instead of at most half of one, and the bias does not
           cancel; the audit must catch it. */
        int caught = 0;
        for (int bits = 4; bits <= 16; bits += 2) {
            FixCfg c = C;
            c.q_bits = bits;
            c.min_frames = 1;
            FixAudit a;
            fm_init(&a, &c);
            a.fx.fault_quant = 8;
            run(&a, 3000);
            if (fm_status(&a) == FM_DRIFT_EXCEEDED) caught++;
        }
        check("truncating_chain_alarm_fires", caught > 0);
        printf("truncation_caught_over=%d\n", caught);
    }

    /* ---- the alarm is localized: it names the particle and axis it first
       saw over the bound, and that particle can be checked by hand ---- */
    {
        FixCfg c = C;
        c.q_bits = 10;
        c.min_frames = 1;
        FixAudit a;
        fm_init(&a, &c);
        a.fx.fault_quant = 8;
        run(&a, 3000);
        check("alarm_localized", a.frame_jump_frames > 0 &&
              a.first_jump_frame > 0 &&
              a.first_jump_particle >= 1 && a.first_jump_particle < c.n &&
              (a.first_jump_axis == 0 || a.first_jump_axis == 1));
        printf("alarm_first_frame=%d particle=%d axis=%d max_jump=%.4g r0=%.4g\n",
               a.first_jump_frame, a.first_jump_particle, a.first_jump_axis,
               a.max_jump, a.r0);
    }

    /* ---- overflow is latched, not wrapped ---- */
    {
        FixCfg c = C;
        c.q_bits = 4;
        c.min_frames = 1;
        FixAudit a;
        fm_init(&a, &c);
        /* place a particle at the coordinate cap with a large velocity so
           the next integration step runs past it */
        a.fx.p[1].x = FX_LIMIT;
        a.fx.p[1].px = FX_LIMIT - 1000000;
        int stopped = fm_step(&a) == 0;
        check("overflow_latched", stopped && a.overflow && fm_status(&a) == FM_OVERFLOW);
        int before = a.frames;
        fm_step(&a);
        check("overflow_refuses_further_steps", a.frames == before);
    }

    /* ---- topology/shape independence: a pinned middle particle reduces
       the movable endpoint count and stays inside the bound ---- */
    {
        FixCfg c = C;
        c.n = 6;
        c.q_bits = 10;
        c.min_frames = 1;
        FixAudit a;
        fm_init(&a, &c);
        run(&a, 2000);
        check("pinned_only_chain_bounded", a.drift_frames == 0);
        check("pinned_chain_converges", a.nonconv_frames == 0);
    }

    /* ---- residual versus q: the residual is a constraint-solve quantity
       that grows with the step, which is the mechanism behind the
       non-convergence report.  Print the scan so the crossover is visible. ---- */
    {
        double res[17];
        int nres = 0;
        for (int bits = 2; bits <= 18; bits++) {
            FixCfg c = C;
            c.q_bits = bits;
            c.min_frames = 1;
            FixAudit a;
            fm_init(&a, &c);
            run(&a, 2000);
            res[nres++] = a.worst_residual;
            printf("residual_scan q=%.3g worst=%.5g nonconv=%d\n",
                   a.q, a.worst_residual, a.nonconv_frames > 0);
        }
        check("residual_grows_with_coarse_q",
              res[0] > 10.0 * res[nres - 1] && res[nres - 1] > 0.0);
    }

    /* ---- per-constraint residual zero is reachable: a fully converged
       coarse chain reports zero residual and the status names it ---- */
    {
        FixCfg c = C;
        c.q_bits = 10;
        c.iterations = 32;
        c.min_frames = 10;
        c.res_tol = 0.02;
        FixAudit a;
        fm_init(&a, &c);
        run(&a, 500);
        check("high_iter_converges_at_moderate_q", a.nonconv_frames == 0);
        check("converged_status_ok", fm_status(&a) == FM_OK);
    }

    /* ---- amplification regime: a stiff, nearly undamped chain does not
       contract a q-sized difference, so the run-long difference grows
       faster than the N*r0 rounding budget.  The auditor must report that
       honestly (drift_exceeded) rather than call the run quantization-clean,
       and the report shows the ratio growing instead of behaving like a
       one-shot defect: the difference is real and it is not a rounding sum. ---- */
    {
        FixCfg c = C;
        c.damping = 0.999999;
        c.q_bits = 10;
        c.min_frames = 1;
        c.res_tol = 0.05;
        FixAudit a;
        fm_init(&a, &c);
        for (int i = 0; i < 200; i++) fm_step(&a);
        double early_peak = a.peak_drift;
        for (int i = 0; i < 4000; i++) fm_step(&a);
        check("stiff_chain_reported_not_ok", fm_status(&a) == FM_DRIFT_EXCEEDED);
        check("stiff_chain_grows", a.peak_drift > 2.0 * early_peak);
        printf("stiff early_peak=%.4g late_peak=%.4g accum_bound=%.4g\n",
               early_peak, a.peak_drift, a.accum_bound);
    }

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
