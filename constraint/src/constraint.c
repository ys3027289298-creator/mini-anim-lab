
/* constraint: verlet particle chain with distance constraints, iteration
   count and damping; used as the toy model for cloth/rope style solvers.

   Includes a per-frame conservation auditor: total length and total
   mechanical energy are compared against their theoretical values every
   frame, violations are localized to the specific distance constraint,
   and a run shorter than the configured sample minimum is reported as
   undetermined rather than as a pass. */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define MAX_P 16
#define MAX_FRAMES 2048

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

/* ------------------------------------------------------------------ */
/* Conservation auditor                                               */
/* ------------------------------------------------------------------ */
/*
   Theory
   ------
   length:  L_ref = sum(rest_i), constant.  The distance constraints are
            what should keep L(t) at L_ref at every frame.
   energy:  E(t) = K + U, with K = sum 1/2 m |(x - px)/dt|^2 and
            U = -sum m (gx*x + gy*y).  The theoretical value is E at the
            first observed frame, i.e. whatever energy the run starts
            with is the baseline.  Under damping <= 1 the system may only
            lose energy, so the default policy flags only *gains* beyond
            en_tol (the "floats up on its own" signature).
            CONS_ENERGY_ABSOLUTE instead flags any deviation from the
            baseline in either direction (useful to catch over-damping
            too).  Note: the semi-implicit Verlet integrator has bounded
            O(dt) oscillation even when undamped, so pick en_tol above
            that floor for nondissipative runs.

   Localization
   ------------
   Constraint i is reported when |len_i - rest_i| > len_tol, exactly the
   criterion a manual per-constraint check uses, so the reported set and
   the manual set agree without gaps.  Energy is a global quantity: when
   energy is violated but no single constraint exceeds len_tol, the
   auditor records a global violation (constraint -1) instead of blaming
   a segment.

   Undetermined
   ------------
   If fewer than cfg.min_frames frames were observed and nothing has been
   seen to violate, the status is CONS_UNDETERMINED -- never CONS_OK.
   An already-observed violation is reported even in a short run, because
   withholding it would be a miss relative to a manual check.  Invalid
   thresholds/config are rejected up front with CONS_INVALID.
*/

typedef enum {
    CONS_OK = 0,
    CONS_UNDETERMINED,
    CONS_VIOLATION,
    CONS_INVALID
} ConsStatus;

typedef enum {
    CONS_ENERGY_NONINCREASE = 0,  /* flag only E(t) > E(0) + en_tol */
    CONS_ENERGY_ABSOLUTE = 1      /* flag |E(t) - E(0)| > en_tol */
} ConsEnergyPolicy;

typedef struct {
    float len_tol;        /* per segment and total-length tolerance, > 0 */
    float en_tol;         /* energy tolerance, > 0 */
    int min_frames;       /* minimum samples needed before judging, >= 1 */
    int energy_policy;    /* ConsEnergyPolicy */
    float dt;             /* integration step used by the chain, > 0 */
    float gx, gy;         /* gravity used by the chain */
} ConsCfg;

typedef struct {
    float length, kinetic, potential, energy;
    float worst_residual;     /* max_i |len_i - rest_i| */
    int worst_constraint;     /* argmax_i, -1 when n_con == 0 */
    int n_hits;               /* constraints with residual > len_tol */
    unsigned hits_mask;       /* bit i set when constraint i is violated */
    float resid[MAX_P];       /* per-constraint |len_i - rest_i| */
} ConsFrame;

typedef struct {
    ConsCfg cfg;
    int n_con;
    float rest[MAX_P];
    float theory_length;
    float theory_energy;
    int frames;
    int invalid;
    int capped;               /* exceeded MAX_FRAMES of stored detail */
    ConsFrame fr[MAX_FRAMES];
    int length_bad_frames;
    int energy_bad_frames;
    int energy_global_frames; /* energy bad, no constraint localized */
    int hit_count[MAX_P];
    int first_hit_frame[MAX_P];
    float peak_resid[MAX_P];  /* worst residual seen per constraint */
} ConsCheck;

static const char* cons_status_name(ConsStatus s) {
    switch (s) {
        case CONS_OK: return "ok";
        case CONS_UNDETERMINED: return "undetermined";
        case CONS_VIOLATION: return "violation";
        case CONS_INVALID: return "invalid";
    }
    return "invalid";
}

static const char* cons_status_reason(ConsStatus s) {
    switch (s) {
        case CONS_OK: return "conserved within tolerance";
        case CONS_UNDETERMINED: return "insufficient samples: cannot determine";
        case CONS_VIOLATION: return "conservation violated";
        case CONS_INVALID: return "invalid threshold or config";
    }
    return "invalid threshold or config";
}

static int cons_config_valid(const ConsCfg* cfg, int n_con) {
    if (!isfinite(cfg->len_tol) || cfg->len_tol <= 0.0f) return 0;
    if (!isfinite(cfg->en_tol) || cfg->en_tol <= 0.0f) return 0;
    if (!isfinite(cfg->dt) || cfg->dt <= 0.0f) return 0;
    if (!isfinite(cfg->gx) || !isfinite(cfg->gy)) return 0;
    if (cfg->min_frames < 1) return 0;
    if (cfg->energy_policy != CONS_ENERGY_NONINCREASE &&
        cfg->energy_policy != CONS_ENERGY_ABSOLUTE) return 0;
    if (n_con < 1 || n_con > MAX_P - 1) return 0;
    return 1;
}

static int cons_energy_bad(const ConsCfg* cfg, float energy, float theory) {
    float dev = energy - theory;
    if (cfg->energy_policy == CONS_ENERGY_ABSOLUTE) return fabsf(dev) > cfg->en_tol;
    return dev > cfg->en_tol;
}

static void cons_measure(const Chain* c, const ConsCheck* chk, ConsFrame* f) {
    const ConsCfg* cfg = &chk->cfg;
    memset(f, 0, sizeof(*f));
    f->worst_constraint = -1;

    float L = 0.0f, worst = 0.0f;
    for (int i = 0; i < chk->n_con; i++) {
        float dx = c->p[i + 1].x - c->p[i].x;
        float dy = c->p[i + 1].y - c->p[i].y;
        float d = sqrtf(dx * dx + dy * dy);
        L += d;
        float r = fabsf(d - chk->rest[i]);
        f->resid[i] = r;
        if (r > cfg->len_tol) {
            f->n_hits++;
            f->hits_mask |= 1u << i;
        }
        if (r > worst) { worst = r; f->worst_constraint = i; }
    }
    f->length = L;
    f->worst_residual = worst;

    float K = 0.0f, U = 0.0f;
    for (int i = 0; i < c->n; i++) {
        const Particle* p = &c->p[i];
        float m = (p->inv_mass > 0.0f) ? 1.0f / p->inv_mass : 0.0f;
        float vx = (p->x - p->px) / cfg->dt;
        float vy = (p->y - p->py) / cfg->dt;
        K += 0.5f * m * (vx * vx + vy * vy);
        U += -m * (cfg->gx * p->x + cfg->gy * p->y);
    }
    f->kinetic = K;
    f->potential = U;
    f->energy = K + U;
}

static void cons_commit(ConsCheck* chk, const ConsFrame* f) {
    int idx = chk->frames;
    if (idx < MAX_FRAMES) chk->fr[idx] = *f;
    else chk->capped = 1;

    if (fabsf(f->length - chk->theory_length) > chk->cfg.len_tol)
        chk->length_bad_frames++;
    if (cons_energy_bad(&chk->cfg, f->energy, chk->theory_energy)) {
        chk->energy_bad_frames++;
        if (f->n_hits == 0) chk->energy_global_frames++;
    }
    for (int i = 0; i < chk->n_con; i++) {
        if (f->resid[i] > chk->peak_resid[i]) chk->peak_resid[i] = f->resid[i];
        if (f->hits_mask & (1u << i)) {
            chk->hit_count[i]++;
            if (chk->first_hit_frame[i] < 0) chk->first_hit_frame[i] = idx;
        }
    }
    chk->frames++;
}

/* Returns 1 when the config is accepted, 0 when rejected (chk->invalid). */
static int cons_init(ConsCheck* chk, const Chain* c, const ConsCfg* cfg) {
    memset(chk, 0, sizeof(*chk));
    chk->cfg = *cfg;
    chk->n_con = c->n - 1;
    for (int i = 0; i < MAX_P; i++) chk->first_hit_frame[i] = -1;

    int ok = cons_config_valid(cfg, chk->n_con);
    if (ok) {
        chk->theory_length = 0.0f;
        for (int i = 0; i < chk->n_con; i++) {
            if (!isfinite(c->rest[i]) || c->rest[i] <= 0.0f) { ok = 0; break; }
            chk->rest[i] = c->rest[i];
            chk->theory_length += c->rest[i];
        }
    }
    if (!ok) { chk->invalid = 1; return 0; }

    ConsFrame f0;
    cons_measure(c, chk, &f0);
    chk->theory_energy = f0.energy;
    cons_commit(chk, &f0);
    return 1;
}

static void cons_observe(ConsCheck* chk, const Chain* c) {
    if (chk->invalid || chk->frames == 0) return;
    /* The topology is fixed at init; a resized chain is a config error,
       not something to silently average over. */
    if (c->n - 1 != chk->n_con) { chk->invalid = 1; return; }
    ConsFrame f;
    cons_measure(c, chk, &f);
    cons_commit(chk, &f);
}

static ConsStatus cons_status(const ConsCheck* chk) {
    if (chk->invalid) return CONS_INVALID;
    /* An observed violation is already a definite finding: more samples
       cannot undo it, and withholding it would disagree with (and miss
       relative to) a manual per-constraint check. */
    if (chk->length_bad_frames > 0 || chk->energy_bad_frames > 0) return CONS_VIOLATION;
    for (int i = 0; i < chk->n_con; i++)
        if (chk->hit_count[i] > 0) return CONS_VIOLATION;
    /* Nothing observed yet, but too few frames to call it conserved. */
    if (chk->frames < chk->cfg.min_frames) return CONS_UNDETERMINED;
    return CONS_OK;
}

/* Worst localized constraint across the run; its peak residual is written
   to *residual when provided.  Returns -1 if no constraint ever exceeded
   len_tol (e.g. an energy-only violation). */
static int cons_worst_constraint(const ConsCheck* chk, float* residual) {
    int idx = -1;
    float worst = 0.0f;
    for (int i = 0; i < chk->n_con; i++) {
        if (chk->hit_count[i] > 0 && (idx < 0 || chk->peak_resid[i] > worst)) {
            idx = i;
            worst = chk->peak_resid[i];
        }
    }
    if (residual) *residual = worst;
    return idx;
}

/* Manual per-constraint check, independent of the auditor. */
static int manual_hits(const Chain* c, float tol, int* out) {
    int n = 0;
    for (int i = 0; i < c->n - 1; i++) {
        float dx = c->p[i + 1].x - c->p[i].x;
        float dy = c->p[i + 1].y - c->p[i].y;
        float d = sqrtf(dx * dx + dy * dy);
        if (fabsf(d - c->rest[i]) > tol) out[n++] = i;
    }
    return n;
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

    /* ---------------- conservation auditor ---------------- */
    const float DT = 0.016f, GX = 0.0f, GY = -9.8f;
    ConsCfg cfg;
    memset(&cfg, 0, sizeof(cfg));
    cfg.len_tol = 0.05f;
    cfg.en_tol = 1e-2f;
    cfg.min_frames = 50;
    cfg.energy_policy = CONS_ENERGY_NONINCREASE;
    cfg.dt = DT; cfg.gx = GX; cfg.gy = GY;

    /* healthy run: length stays at sum(rest), energy never rises, and the
       first frame establishes the theory energy baseline. */
    {
        Chain h;
        chain_init(&h, 8, 1.0f, 0.0f, 0.0f);
        ConsCheck chk;
        check("cons_accepts_valid_cfg", cons_init(&chk, &h, &cfg) == 1);
        check("cons_theory_length", fabsf(chk.theory_length - 7.0f) < 1e-5f);
        check("cons_theory_energy_from_frame0",
              fabsf(chk.theory_energy - chk.fr[0].energy) < 1e-6f);
        check("cons_initial_frame_recorded", chk.frames == 1);
    }

    /* per-frame recording and healthy verdict */
    {
        Chain h;
        chain_init(&h, 8, 1.0f, 0.0f, 0.0f);
        ConsCheck chk;
        cons_init(&chk, &h, &cfg);
        int hits = 0;
        for (int i = 0; i < 200; i++) {
            chain_step(&h, DT, GX, GY);
            cons_observe(&chk, &h);
            if (chk.fr[chk.frames - 1].n_hits) hits++;
        }
        check("cons_records_every_frame", chk.frames == 201);
        check("cons_healthy_ok", cons_status(&chk) == CONS_OK);
        check("cons_healthy_no_constraint_hits", hits == 0 && chk.length_bad_frames == 0);
        check("cons_healthy_no_energy_gain", chk.energy_bad_frames == 0);
        /* theory values actually compared: length close to 7 every frame */
        float maxdev = 0.0f;
        for (int i = 0; i < chk.frames; i++) {
            float d = fabsf(chk.fr[i].length - chk.theory_length);
            if (d > maxdev) maxdev = d;
        }
        check("cons_length_within_tol", maxdev <= cfg.len_tol);
        printf("cons_healthy_max_len_dev=%.6g\n", maxdev);
    }

    /* insufficient sampling is explicitly undetermined, never "ok". */
    {
        Chain h;
        chain_init(&h, 8, 1.0f, 0.0f, 0.0f);
        ConsCfg few = cfg;
        few.min_frames = 10;
        ConsCheck chk;
        cons_init(&chk, &h, &few);
        for (int i = 0; i < 4; i++) { chain_step(&h, DT, GX, GY); cons_observe(&chk, &h); }
        ConsStatus s = cons_status(&chk);
        check("cons_insufficient_undetermined", s == CONS_UNDETERMINED);
        check("cons_insufficient_not_ok", s != CONS_OK);
        printf("cons_insufficient_reason=\"%s\"\n", cons_status_reason(s));
    }

    /* exactly min_frames samples is enough to judge. */
    {
        Chain h;
        chain_init(&h, 8, 1.0f, 0.0f, 0.0f);
        ConsCfg bound = cfg;
        bound.min_frames = 5;
        ConsCheck chk;
        cons_init(&chk, &h, &bound);       /* frame 0 */
        for (int i = 0; i < 4; i++) { chain_step(&h, DT, GX, GY); cons_observe(&chk, &h); }
        check("cons_frames_equal_min_determined",
              chk.frames == 5 && cons_status(&chk) != CONS_UNDETERMINED);
    }

    /* sample count alone does not mask an observed violation: a short run
       with a broken segment is still reported (no misses vs manual check),
       while the same short run with nothing wrong stays undetermined. */
    {
        Chain h;
        chain_init(&h, 8, 1.0f, 0.0f, 0.0f);
        h.p[4].x += 5.0f;
        ConsCfg few = cfg;
        few.min_frames = 10;
        ConsCheck chk;
        cons_init(&chk, &h, &few);
        h.p[4].x += 5.0f;                  /* keep it visibly broken */
        cons_observe(&chk, &h);
        check("cons_short_bad_run_reported_not_masked",
              chk.frames < few.min_frames && cons_status(&chk) == CONS_VIOLATION);
    }

    /* localization: a deliberately corrupted segment is reported, and the
       reported set equals the manual per-constraint set. */
    {
        Chain h;
        chain_init(&h, 8, 1.0f, 0.0f, 0.0f);
        ConsCheck chk;
        cons_init(&chk, &h, &cfg);
        h.p[4].x += 3.0f;                  /* stretches constraints 3 and 4 */
        cons_observe(&chk, &h);
        const ConsFrame* f = &chk.fr[chk.frames - 1];
        int mi[MAX_P];
        int mn = manual_hits(&h, cfg.len_tol, mi);
        int same = (mn == f->n_hits);
        for (int k = 0; k < mn; k++) same = same && ((f->hits_mask >> mi[k]) & 1u);
        check("cons_localize_matches_manual", same && mn == 2);
        check("cons_localize_worst_segment",
              f->worst_constraint == 3 || f->worst_constraint == 4);
        check("cons_localize_status_violation", cons_status(&chk) == CONS_VIOLATION);
        printf("cons_localized_hits=%d worst=%d\n", f->n_hits, f->worst_constraint);
    }

    /* energy can be violated while all segments look fine ("floats up"):
       the auditor must flag it, and must not invent a segment to blame. */
    {
        Chain h;
        chain_init(&h, 8, 1.0f, 0.0f, 0.0f);
        h.damping = 1.005f;                /* injects energy every step */
        ConsCfg ecfg = cfg;
        ecfg.min_frames = 50;
        ecfg.en_tol = 1e-2f;
        ConsCheck chk;
        cons_init(&chk, &h, &ecfg);
        for (int i = 0; i < 400; i++) { chain_step(&h, DT, GX, GY); cons_observe(&chk, &h); }
        ConsStatus s = cons_status(&chk);
        check("cons_floatup_detected", s == CONS_VIOLATION && chk.energy_bad_frames > 0);
        check("cons_floatup_global_not_misattributed",
              chk.energy_global_frames > 0);
        check("cons_floatup_blames_no_constraint",
              cons_worst_constraint(&chk, NULL) == -1);
        printf("cons_floatup_status=%s energy_bad=%d global=%d max_seg=%.4g\n",
               cons_status_name(s), chk.energy_bad_frames, chk.energy_global_frames,
               chk.energy_global_frames ? chk.fr[chk.frames - 1].worst_residual : 0.0f);
    }

    /* absolute energy policy also flags energy *loss* (e.g. violent
       damping), which the non-increase policy deliberately allows. */
    {
        Chain h;
        chain_init(&h, 8, 1.0f, 0.0f, 0.0f);
        h.damping = 0.5f;
        ConsCfg acfg = cfg;
        acfg.min_frames = 50;
        acfg.energy_policy = CONS_ENERGY_ABSOLUTE;
        ConsCheck chk;
        cons_init(&chk, &h, &acfg);
        for (int i = 0; i < 200; i++) { chain_step(&h, DT, GX, GY); cons_observe(&chk, &h); }
        check("cons_absolute_flags_energy_loss", cons_status(&chk) == CONS_VIOLATION);

        Chain g;
        chain_init(&g, 8, 1.0f, 0.0f, 0.0f);
        g.damping = 0.5f;
        ConsCfg ncfg = cfg;
        ncfg.min_frames = 50;
        ncfg.energy_policy = CONS_ENERGY_NONINCREASE;
        ConsCheck nchk;
        cons_init(&nchk, &g, &ncfg);
        for (int i = 0; i < 200; i++) { chain_step(&g, DT, GX, GY); cons_observe(&nchk, &g); }
        check("cons_nonincrease_allows_energy_loss", cons_status(&nchk) == CONS_OK);
    }

    /* zero misses against a manual per-constraint check across a long,
       deliberately unstable run. */
    {
        Chain h;
        chain_init(&h, 8, 1.0f, 0.0f, 0.0f);
        h.damping = 1.05f;
        ConsCfg acfg = cfg;
        acfg.min_frames = 50;
        ConsCheck chk;
        cons_init(&chk, &h, &acfg);
        int misses = 0, falsepos = 0, localized = 0;
        for (int i = 0; i < 300; i++) {
            chain_step(&h, DT, GX, GY);
            cons_observe(&chk, &h);
            const ConsFrame* f = &chk.fr[chk.frames - 1];
            int mi[MAX_P];
            int mn = manual_hits(&h, acfg.len_tol, mi);
            int want[MAX_P] = { 0 };
            for (int k = 0; k < mn; k++) {
                want[mi[k]] = 1;
                if (!((f->hits_mask >> mi[k]) & 1u)) misses++;
            }
            for (int k = 0; k < chk.n_con; k++) {
                int got = (int)((f->hits_mask >> k) & 1u);
                if (got && !want[k]) falsepos++;
                if (got) localized++;
            }
        }
        check("cons_zero_misses_vs_manual", misses == 0);
        check("cons_no_false_positive_vs_manual", falsepos == 0);
        check("cons_unstable_run_has_localized_hits", localized > 0);
        printf("cons_audit_misses=%d falsepos=%d hits=%d frames=%d\n",
               misses, falsepos, localized, chk.frames);
    }

    /* invalid thresholds and configs are rejected outright. */
    {
        Chain h;
        chain_init(&h, 8, 1.0f, 0.0f, 0.0f);
        ConsCheck chk;
        ConsCfg bad = cfg;

        bad = cfg; bad.len_tol = 0.0f;
        check("cons_reject_len_tol_zero", cons_init(&chk, &h, &bad) == 0 &&
              cons_status(&chk) == CONS_INVALID);
        bad = cfg; bad.len_tol = -1.0f;
        check("cons_reject_len_tol_negative", cons_init(&chk, &h, &bad) == 0);
        bad = cfg; bad.len_tol = NAN;
        check("cons_reject_len_tol_nan", cons_init(&chk, &h, &bad) == 0);
        bad = cfg; bad.len_tol = INFINITY;
        check("cons_reject_len_tol_inf", cons_init(&chk, &h, &bad) == 0);
        bad = cfg; bad.en_tol = 0.0f;
        check("cons_reject_en_tol_zero", cons_init(&chk, &h, &bad) == 0);
        bad = cfg; bad.en_tol = NAN;
        check("cons_reject_en_tol_nan", cons_init(&chk, &h, &bad) == 0);
        bad = cfg; bad.en_tol = INFINITY;
        check("cons_reject_en_tol_inf", cons_init(&chk, &h, &bad) == 0);
        bad = cfg; bad.min_frames = 0;
        check("cons_reject_min_frames_zero", cons_init(&chk, &h, &bad) == 0);
        bad = cfg; bad.dt = 0.0f;
        check("cons_reject_dt_zero", cons_init(&chk, &h, &bad) == 0);
        bad = cfg; bad.dt = -0.016f;
        check("cons_reject_dt_negative", cons_init(&chk, &h, &bad) == 0);
        bad = cfg; bad.energy_policy = 42;
        check("cons_reject_bad_policy", cons_init(&chk, &h, &bad) == 0);
        /* observe after rejection must not revive the check */
        cons_observe(&chk, &h);
        check("cons_rejected_stays_invalid",
              cons_status(&chk) == CONS_INVALID && chk.frames == 0);
        printf("cons_invalid_reason=\"%s\"\n", cons_status_reason(CONS_INVALID));
    }

    /* topology changes after init are rejected, not averaged in. */
    {
        Chain h;
        chain_init(&h, 8, 1.0f, 0.0f, 0.0f);
        ConsCheck chk;
        cons_init(&chk, &h, &cfg);
        h.n = 7;                            /* shrink behind the auditor */
        cons_observe(&chk, &h);
        check("cons_reject_topology_change",
              cons_status(&chk) == CONS_INVALID && chk.frames == 1);
    }

    /* detail buffer overflow must not lose the verdict: counters keep
       accumulating past MAX_FRAMES and status stays correct. */
    {
        Chain h;
        chain_init(&h, 8, 1.0f, 0.0f, 0.0f);
        ConsCfg lc = cfg;
        lc.min_frames = 10;
        ConsCheck chk;
        cons_init(&chk, &h, &lc);
        h.p[4].x += 5.0f;
        cons_observe(&chk, &h);             /* one clearly bad frame */
        for (int i = 0; i < MAX_FRAMES + 100; i++) cons_observe(&chk, &h);
        check("cons_overflow_keeps_counting",
              chk.capped && chk.frames == MAX_FRAMES + 102);
        check("cons_overflow_keeps_verdict",
              cons_status(&chk) == CONS_VIOLATION && chk.hit_count[3] > 0);
    }

    /* randomized differential: auditor's localized set must equal the
       manual per-constraint set on every frame, with zero misses. */
    {
        unsigned seed = 12345u;
        int total_misses = 0, total_fp = 0, frames = 0, trips = 0;
        for (int trial = 0; trial < 40; trial++) {
            Chain h;
            chain_init(&h, 4 + (int)(seed % 12), 0.25f + (float)(seed % 7) * 0.25f, 0.0f, 0.0f);
            seed = seed * 1103515245u + 12345u;
            h.damping = 0.9f + (float)(seed % 21) * 0.01f;   /* 0.90..1.10 */
            ConsCfg rc = cfg;
            rc.min_frames = 1;
            ConsCheck chk;
            if (!cons_init(&chk, &h, &rc)) continue;
            for (int f = 0; f < 120; f++) {
                chain_step(&h, DT, GX, GY);
                if ((seed & 63u) == 0) {     /* occasional impulse */
                    h.p[1 + (int)(seed % (unsigned)(h.n - 1))].x += 2.0f;
                }
                seed = seed * 1103515245u + 12345u;
                cons_observe(&chk, &h);
                const ConsFrame* fr = &chk.fr[chk.frames - 1];
                int mi[MAX_P];
                int mn = manual_hits(&h, rc.len_tol, mi);
                int want[MAX_P] = { 0 };
                for (int k = 0; k < mn; k++) {
                    want[mi[k]] = 1;
                    if (!((fr->hits_mask >> mi[k]) & 1u)) total_misses++;
                }
                for (int k = 0; k < chk.n_con; k++) {
                    int got = (int)((fr->hits_mask >> k) & 1u);
                    if (got && !want[k]) total_fp++;
                }
                if (mn) trips++;
                frames++;
            }
        }
        check("cons_diff_no_misses", total_misses == 0);
        check("cons_diff_no_false_positive", total_fp == 0);
        check("cons_diff_exercised", frames > 1000 && trips > 0);
        printf("cons_diff_frames=%d constrained_trips=%d misses=%d fp=%d\n",
               frames, trips, total_misses, total_fp);
    }

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
