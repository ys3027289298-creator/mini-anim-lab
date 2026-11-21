
/* sleep: rigid-body sleep detection from mean kinetic energy and positional
   variance, with hysteresis, an AABB that is maintained while a body sleeps,
   and frame-level scoring against hand-labelled rest intervals.

   The public-ish sleep_* functions are the library; main() is the self-test,
   matching the other modules in this repo (no header, run `make test`).

   Sign conventions and units
   --------------------------
   energy[i]   : 0.5 * inv_mass * |v|^2, i.e. per-body kinetic energy, in
                 mass * distance^2 / time^2 (J when mass is kg).
   variance[i] : sum over x,y,z of the per-axis sample variance of the last
                 SL_NSTEP positions, in distance^2, i.e. the trace of the 3x3
                 position covariance.  It is a sum rather than a mean so that
                 a single displaced sample does not get diluted by axes that
                 happen to be still: a 1D +/-d two-state jitter reads as d^2
                 and a 3D jitter of the same amplitude reads as 3d^2.
   Both quantities are averaged over bodies, so thresholds do not have to be
   re-tuned when the body count changes.  The total (sum) is still reported
   as energy_total for diagnostics.

   Hysteresis: a body enters sleep only when BOTH quantities sit below the
   "enter" thresholds, and leaves only when EITHER quantity exceeds the
   "exit" thresholds, with enter < exit.  This is the standard schmitt
   trigger; without the gap a body parked on a threshold changes state every
   tick and the sleep count flickers.  Negative, NaN or infinite thresholds
   are rejected at configuration time and again on every tick.
*/
#include <stdio.h>
#include <string.h>
#include <math.h>

#define SL_NSTEP    4       /* velocity-free position history for variance */
#define SL_MAXBODY  8
#define SL_MAXFRAME 256

typedef struct { float x, y, z; } SlVec3;

typedef struct { SlVec3 lo, hi; int valid; } SlAabb;

typedef struct {
    SlVec3 pos, vel;
    SlVec3 half;            /* AABB half extents in body space */
    float inv_mass;         /* > 0; a fixed/infinite-mass body is rejected */
} SlBody;

typedef struct {
    float e_enter;          /* kinetic energy must be BELOW this to sleep */
    float e_exit;           /* kinetic energy ABOVE this wakes the body */
    float var_enter;        /* positional variance BELOW this to sleep */
    float var_exit;         /* positional variance ABOVE this wakes it */
    float hold_time;        /* seconds both conditions must hold, > 0 */
    float fixed_dt;         /* seconds per frame, > 0 */
    int   substeps;         /* physics substeps per frame, >= 1 */
} SlSleepConfig;

typedef struct { SlVec3 pos[SL_MAXBODY]; } SlHistFrame;

typedef struct {
    int body_count;
    SlSleepConfig cfg;
    int hold_ticks;         /* hold_time expressed in whole ticks */
    int candidate[SL_MAXBODY];   /* consecutive quiet ticks per body */
    unsigned char sleeping[SL_MAXBODY];
    unsigned char initialized;
    int ticks;              /* history entries, saturates at SL_NSTEP */
    SlHistFrame hist[SL_NSTEP];
    SlAabb aabb[SL_MAXBODY];   /* per-body swept box, valid while asleep */
    SlAabb bounds;             /* union of every body, valid while any sleeps */
    int wake_count;         /* sleep-state transitions (either direction) */
    int valid;
} SlState;

typedef struct {
    float energy[SL_MAXBODY];
    float variance[SL_MAXBODY];
    float energy_mean;      /* mean over bodies, this is what is thresholded */
    float variance_mean;
    float energy_total;     /* sum over bodies, reported for diagnostics */
} SlMeasure;

typedef struct {
    int valid;
    float energy_mean;
    float variance_mean;
    float energy_total;         /* 0.5 * sum(m * v^2) over the rig */
    unsigned char sleeping[SL_MAXBODY];
    SlAabb bbox[SL_MAXBODY];   /* per-body swept AABB while asleep */
    SlAabb bounds;             /* union of all bodies, valid while any sleeps */
    int wake_count;
} SlTick;

typedef struct { int start; int end; } SlInterval;  /* [start, end) frames */

typedef struct {
    float fa_frac;          /* allowed fraction of labelled-asleep frames awake */
    float miss_frac;        /* allowed fraction of labelled-awake frames asleep */
    float lat_frac;         /* allowed fraction of wake events over the cap */
    int   max_latency;      /* frames, >= 1 */
} SlEvalConfig;

typedef struct {
    int ok;                 /* all measured counts within their bounds */
    int fa_frames;          /* truth=asleep, prediction=awake */
    int truth_asleep;
    int fa_bound;
    int miss_frames;        /* truth=awake, prediction=asleep */
    int truth_awake;
    int miss_bound;
    int wake_runs;          /* truth wake events inside the window */
    int wake_exceeded;      /* those whose detection latency > max_latency */
    int wake_bound;
    int worst_latency;      /* frames, 0 when no wake event */
} SlEval;

/* ------------------------------------------------------------------ */
/* helpers                                                             */
/* ------------------------------------------------------------------ */

static SlVec3 v3(float x, float y, float z) { SlVec3 v = { x, y, z }; return v; }

/* Seconds represented by one sleep tick.  A tick advances the physics by
   fixed_dt / substeps of simulated time, so hold_time and wake latency are
   independent of how many substeps the caller runs per frame. */
static float sleep_timescale(const SlSleepConfig* cfg) {
    float t = cfg->fixed_dt;
    if (cfg->substeps > 0) t /= (float)cfg->substeps;
    if (!(t > 0.0f) || !isfinite(t)) t = 1.0f;   /* defensive; cfg is checked */
    return t;
}

/* Number of quiet ticks a body must accumulate.  The deadline is rounded up
   to a whole tick, so 0.5 s at 60 Hz means 30 ticks rather than 31.25.  It
   is kept as an integer because accumulating fractional seconds in float
   drifts: 0.016 * 32 is 0.51199997, which compares below 0.512 and adds a
   spurious extra frame. */
static int sleep_hold_ticks(const SlSleepConfig* cfg) {
    float ts = sleep_timescale(cfg);
    float ticks = cfg->hold_time / ts;
    if (!(ticks >= 1.0f)) ticks = 1.0f;
    if (ticks > 1.0e7f) ticks = 1.0e7f;
    return (int)ceilf(ticks);
}

static int finite3(SlVec3 v) {
    return isfinite(v.x) && isfinite(v.y) && isfinite(v.z);
}

/* Thresholds must be strictly positive and ordered; a negative threshold is
   meaningless for a non-negative measure, and enter >= exit reproduces the
   flicker the hysteresis exists to remove. */
static int sleep_config_valid(const SlSleepConfig* cfg) {
    if (!cfg) return 0;
    if (!(cfg->e_enter > 0.0f) || !isfinite(cfg->e_enter)) return 0;
    if (!(cfg->e_exit > 0.0f) || !isfinite(cfg->e_exit)) return 0;
    if (!(cfg->var_enter > 0.0f) || !isfinite(cfg->var_enter)) return 0;
    if (!(cfg->var_exit > 0.0f) || !isfinite(cfg->var_exit)) return 0;
    if (!(cfg->e_exit > cfg->e_enter)) return 0;
    if (!(cfg->var_exit > cfg->var_enter)) return 0;
    if (!(cfg->hold_time > 0.0f) || !isfinite(cfg->hold_time)) return 0;
    if (!(cfg->fixed_dt > 0.0f) || !isfinite(cfg->fixed_dt)) return 0;
    if (cfg->substeps < 1) return 0;
    return 1;
}

static int sleep_body_valid(const SlBody* b) {
    if (!(b->inv_mass > 0.0f) || !isfinite(b->inv_mass)) return 0;
    if (!finite3(b->pos) || !finite3(b->vel) || !finite3(b->half)) return 0;
    if (b->half.x < 0.0f || b->half.y < 0.0f || b->half.z < 0.0f) return 0;
    return 1;
}

/* Defaults are derived in README.md; keep the two comments in sync.
   e_enter = 0.5 * m * v_quiet^2 with v_quiet = 4.5e-2 u/s and m = 1:
           0.5 * 1 * 0.045^2 = 1.0e-3.
   var_enter = d^2 with d = 1.0e-2 u positional jitter -> 1.0e-4.
   Both exits are 2x their enter thresholds (r = 2 sits mid-range in
   [1.5, 4]: r >= 1.5 keeps resting-contact noise from waking the body,
   r <= 4 keeps a slow push from waiting more than ~2 frames to wake).
   hold_time = 0.5 s is 30 ticks at 60 Hz, enough for the 4-sample variance
   to fill and to ride out a single contact-resolution spike. */
static SlSleepConfig sleep_config_default(void) {
    SlSleepConfig c;
    memset(&c, 0, sizeof(c));
    c.e_enter = 1.0e-3f;
    c.e_exit = 2.0e-3f;
    c.var_enter = 1.0e-4f;
    c.var_exit = 2.0e-4f;
    c.hold_time = 0.5f;
    c.fixed_dt = 1.0f / 60.0f;
    c.substeps = 1;
    return c;
}

static SlEvalConfig sleep_eval_config_default(void) {
    SlEvalConfig c;
    c.fa_frac = 0.05f;
    c.miss_frac = 0.05f;
    c.lat_frac = 0.10f;
    c.max_latency = 6;
    return c;
}

/* ------------------------------------------------------------------ */
/* measurement                                                         */
/* ------------------------------------------------------------------ */

/* Per-body kinetic energy and positional variance.  Welford's algorithm is
   used instead of sum-of-squares because positions can be far from the
   origin (a character at x = 1000 with 1 cm jitter overflows the precision
   of E[x^2] - E[x]^2 in float). */
static SlMeasure sleep_measure(const SlBody* b, int n,
                               const SlHistFrame* hist, int count) {
    SlMeasure m;
    memset(&m, 0, sizeof(m));
    if (count < 1) count = 1;
    float inv_count = 1.0f / (float)count;

    for (int i = 0; i < n; i++) {
        float vx = b[i].vel.x, vy = b[i].vel.y, vz = b[i].vel.z;
        m.energy[i] = 0.5f * (vx * vx + vy * vy + vz * vz) / b[i].inv_mass;

        float mean[3] = { 0.0f, 0.0f, 0.0f };
        float m2[3] = { 0.0f, 0.0f, 0.0f };
        for (int k = 0; k < count; k++) {
            float x[3] = { hist[k].pos[i].x, hist[k].pos[i].y, hist[k].pos[i].z };
            float inv_k = 1.0f / (float)(k + 1);
            for (int a = 0; a < 3; a++) {
                float d = x[a] - mean[a];
                mean[a] += d * inv_k;
                m2[a] += d * (x[a] - mean[a]);
            }
        }
        float var = 0.0f;
        for (int a = 0; a < 3; a++) {
            float v = m2[a] * inv_count;
            if (v > 0.0f) var += v;     /* guard against a tiny negative */
        }
        m.variance[i] = var;

        m.energy_mean += m.energy[i];
        m.variance_mean += m.variance[i];
        m.energy_total += m.energy[i];
    }
    m.energy_mean /= (float)n;
    m.variance_mean /= (float)n;
    return m;
}

/* ------------------------------------------------------------------ */
/* AABB maintenance                                                    */
/* ------------------------------------------------------------------ */

static void aabb_grow_point(SlAabb* a, SlVec3 p) {
    if (!a->valid) {
        a->lo = p;
        a->hi = p;
        a->valid = 1;
        return;
    }
    if (p.x < a->lo.x) a->lo.x = p.x;
    if (p.y < a->lo.y) a->lo.y = p.y;
    if (p.z < a->lo.z) a->lo.z = p.z;
    if (p.x > a->hi.x) a->hi.x = p.x;
    if (p.y > a->hi.y) a->hi.y = p.y;
    if (p.z > a->hi.z) a->hi.z = p.z;
}

static void aabb_union(SlAabb* a, const SlAabb* b) {
    if (!b->valid) return;
    aabb_grow_point(a, b->lo);
    aabb_grow_point(a, b->hi);
}

/* ------------------------------------------------------------------ */
/* detector                                                            */
/* ------------------------------------------------------------------ */

static int sleep_init(SlState* st, int body_count, const SlSleepConfig* cfg) {
    if (!st || !cfg) return 0;
    if (body_count < 1 || body_count > SL_MAXBODY) return 0;
    memset(st, 0, sizeof(*st));
    if (!sleep_config_valid(cfg)) return 0;
    st->body_count = body_count;
    st->cfg = *cfg;
    st->hold_ticks = sleep_hold_ticks(cfg);
    st->valid = 1;
    return 1;
}

/* Force a body awake.  Use this when an external event (a contact, an
   impulse, a scripted animation) must not wait for the next measure; wake
   latency in SlEval is measured against this kind of immediate wake. */
static int sleep_wake(SlState* st, int i) {
    if (!st || !st->valid) return 0;
    if (i < 0 || i >= st->body_count) return 0;
    if (st->sleeping[i]) {
        st->sleeping[i] = 0;
        st->wake_count++;
    }
    st->candidate[i] = 0;
    st->aabb[i].valid = 0;
    return 1;
}

/* One detector tick.  bodies[i] must stay at index i between calls; state is
   per body, not per pointer.  Returns 0 on bad input (NaN position, invalid
   config, wrong body count), leaving *out zeroed. */
static int sleep_tick(const SlBody* bodies, int n, SlState* st, SlTick* out) {
    if (!bodies || !st || !out) return 0;
    memset(out, 0, sizeof(*out));
    if (!st->valid || n != st->body_count || n < 1 || n > SL_MAXBODY) return 0;
    if (!sleep_config_valid(&st->cfg)) return 0;
    for (int i = 0; i < n; i++)
        if (!sleep_body_valid(&bodies[i])) return 0;

    /* history: index 0 is the current frame, SL_NSTEP-1 the oldest */
    for (int k = SL_NSTEP - 1; k > 0; k--) st->hist[k] = st->hist[k - 1];
    for (int i = 0; i < n; i++) st->hist[0].pos[i] = bodies[i].pos;
    if (st->ticks < SL_NSTEP) st->ticks++;

    SlMeasure m = sleep_measure(bodies, n, st->hist, st->ticks);
    if (!isfinite(m.energy_mean) || !isfinite(m.variance_mean)) return 0;

    /* State transition.  The first tick only seeds the history; it must not
       report every body as asleep merely because no variance exists yet. */
    for (int i = 0; i < n; i++) {
        if (!st->initialized) break;
        int before = st->sleeping[i];
        if (st->sleeping[i]) {
            /* while asleep, hold until EITHER measure leaves the band */
            if (m.energy[i] > st->cfg.e_exit || m.variance[i] > st->cfg.var_exit) {
                st->sleeping[i] = 0;
                st->candidate[i] = 0;
            }
        } else {
            /* while awake, require BOTH measures inside the enter band for
               hold_ticks consecutive ticks; any excursion restarts the clock */
            if (m.energy[i] < st->cfg.e_enter && m.variance[i] < st->cfg.var_enter) {
                if (st->candidate[i] < st->hold_ticks) st->candidate[i]++;
                if (st->candidate[i] >= st->hold_ticks) {
                    st->sleeping[i] = 1;
                    st->candidate[i] = 0;
                }
            } else {
                st->candidate[i] = 0;
            }
        }
        if (st->sleeping[i] != before) st->wake_count++;
    }
    st->initialized = 1;

    /* Physical bounds: every body's current world AABB.  Each sleeping body
       grows its own stored box by its own current extent, so a contact proxy
       keeps covering the pose even if the integrator nudges a "sleeping"
       body; the union box covers the whole rig for a single broad-phase
       entry.  A body's box is dropped as soon as it wakes. */
    SlAabb body_now[SL_MAXBODY];
    for (int i = 0; i < n; i++) {
        body_now[i].valid = 1;
        body_now[i].lo = v3(bodies[i].pos.x - bodies[i].half.x,
                            bodies[i].pos.y - bodies[i].half.y,
                            bodies[i].pos.z - bodies[i].half.z);
        body_now[i].hi = v3(bodies[i].pos.x + bodies[i].half.x,
                            bodies[i].pos.y + bodies[i].half.y,
                            bodies[i].pos.z + bodies[i].half.z);
    }
    SlAabb bounds;
    memset(&bounds, 0, sizeof(bounds));
    for (int i = 0; i < n; i++) {
        if (st->sleeping[i]) {
            if (!st->aabb[i].valid) st->aabb[i] = body_now[i];
            aabb_union(&st->aabb[i], &body_now[i]);
            out->bbox[i] = st->aabb[i];
            aabb_union(&bounds, &st->aabb[i]);
        } else {
            st->aabb[i].valid = 0;      /* stale once awake; rebuilt on sleep */
            out->bbox[i].valid = 0;
            aabb_union(&bounds, &body_now[i]);
        }
        out->sleeping[i] = st->sleeping[i];
    }
    out->bounds = bounds;

    out->energy_mean = m.energy_mean;
    out->variance_mean = m.variance_mean;
    out->energy_total = m.energy_total;
    out->wake_count = st->wake_count;
    out->valid = 1;
    return 1;
}

/* ------------------------------------------------------------------ */
/* scoring against hand-labelled intervals                             */
/* ------------------------------------------------------------------ */

/* Frames are clipped to [0, nframes); an interval with start >= end and an
   interval entirely outside the window contribute nothing.  Overlapping
   intervals are unioned, so a frame is labelled asleep at most once. */
static void sl_clip_interval(const SlInterval* in, int nframes, int* s, int* e) {
    int a = in->start, b = in->end;
    if (a < 0) a = 0;
    if (a > nframes) a = nframes;
    if (b < 0) b = 0;
    if (b > nframes) b = nframes;
    if (b < a) b = a;
    *s = a; *e = b;
}

static int sleep_evaluate(const unsigned char* pred, int nframes,
                          const SlInterval* truth, int tcount,
                          const SlEvalConfig* cfg, SlEval* ev) {
    if (!ev) return 0;
    memset(ev, 0, sizeof(*ev));
    if (!pred || !cfg) return 0;
    if (tcount < 0) return 0;
    if (tcount > 0 && !truth) return 0;   /* no labels is legal; null ones are not */
    if (nframes < 1 || nframes > SL_MAXFRAME) return 0;
    if (!isfinite(cfg->fa_frac) || !isfinite(cfg->miss_frac) ||
        !isfinite(cfg->lat_frac)) return 0;
    if (cfg->fa_frac < 0.0f || cfg->fa_frac > 1.0f) return 0;
    if (cfg->miss_frac < 0.0f || cfg->miss_frac > 1.0f) return 0;
    if (cfg->lat_frac < 0.0f || cfg->lat_frac > 1.0f) return 0;
    if (cfg->max_latency < 1) return 0;

    unsigned char is_asleep[SL_MAXFRAME];
    memset(is_asleep, 0, sizeof(is_asleep));
    for (int k = 0; k < tcount; k++) {
        int s, e;
        sl_clip_interval(&truth[k], nframes, &s, &e);
        for (int f = s; f < e; f++) is_asleep[f] = 1;
    }
    for (int f = 0; f < nframes; f++) {
        if (is_asleep[f]) {
            ev->truth_asleep++;
            if (!pred[f]) ev->fa_frames++;
        } else {
            ev->truth_awake++;
            if (pred[f]) ev->miss_frames++;
        }
    }

    /* Wake latency is measured from the truth boundary: the first frame that
       the artist labelled awake after a labelled rest.  If the detector is
       still asleep there, count frames until it reports awake.  g == nframes
       means it did not wake inside the window and yields nframes - f, a
       deliberately large value that the cap will flag. */
    for (int f = 1; f < nframes; f++) {
        if (!is_asleep[f - 1] || is_asleep[f]) continue;
        int g = f;
        while (g < nframes && pred[g]) g++;   /* advance while still asleep */
        int lat = g - f;
        ev->wake_runs++;
        if (lat > cfg->max_latency) ev->wake_exceeded++;
        if (lat > ev->worst_latency) ev->worst_latency = lat;
    }

    /* Upper bounds: ceil so a single frame is tolerated when a positive
       fraction is requested; a fraction of 0 means "no deviation allowed". */
    ev->fa_bound = (int)ceilf(cfg->fa_frac * (float)ev->truth_asleep - 1e-6f);
    ev->miss_bound = (int)ceilf(cfg->miss_frac * (float)ev->truth_awake - 1e-6f);
    ev->wake_bound = (int)ceilf(cfg->lat_frac * (float)ev->wake_runs - 1e-6f);

    ev->ok = (ev->fa_frames <= ev->fa_bound) &&
             (ev->miss_frames <= ev->miss_bound) &&
             (ev->wake_exceeded <= ev->wake_bound);
    return 1;
}

/* ------------------------------------------------------------------ */
/* self-test                                                           */
/* ------------------------------------------------------------------ */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static SlBody make_body(float x, float y, float z, float hx, float hy) {
    SlBody b;
    memset(&b, 0, sizeof(b));
    b.pos = v3(x, y, z);
    b.half = v3(hx, hy, 0.0f);
    b.inv_mass = 1.0f;
    return b;
}

/* thresholds chosen so the test arithmetic is obvious:
   enter at E=10, V=100; exit at E=20, V=200. */
static SlSleepConfig test_cfg(void) {
    SlSleepConfig c = sleep_config_default();
    c.e_enter = 10.0f;
    c.e_exit = 20.0f;
    c.var_enter = 100.0f;
    c.var_exit = 200.0f;
    c.hold_time = 0.5f;
    c.fixed_dt = 0.016f;
    c.substeps = 1;
    return c;
}

int main(void) {
    /* ---- config validation, including negative thresholds ---- */
    {
        SlSleepConfig c = sleep_config_default();
        check("cfg_default_valid", sleep_config_valid(&c) == 1);

        SlSleepConfig bad = c;
        bad.e_enter = -1.0f;
        check("cfg_negative_energy_enter", sleep_config_valid(&bad) == 0);
        bad = c; bad.e_exit = -1.0f;
        check("cfg_negative_energy_exit", sleep_config_valid(&bad) == 0);
        bad = c; bad.var_enter = -1.0f;
        check("cfg_negative_variance_enter", sleep_config_valid(&bad) == 0);
        bad = c; bad.var_exit = -1.0f;
        check("cfg_negative_variance_exit", sleep_config_valid(&bad) == 0);
        bad = c; bad.e_exit = bad.e_enter;
        check("cfg_hysteresis_required", sleep_config_valid(&bad) == 0);
        bad = c; bad.var_exit = bad.var_enter;
        check("cfg_variance_hysteresis_required", sleep_config_valid(&bad) == 0);
        bad = c; bad.hold_time = 0.0f;
        check("cfg_zero_hold", sleep_config_valid(&bad) == 0);
        bad = c; bad.fixed_dt = -0.016f;
        check("cfg_negative_dt", sleep_config_valid(&bad) == 0);
        bad = c; bad.substeps = 0;
        check("cfg_zero_substeps", sleep_config_valid(&bad) == 0);
        bad = c; bad.e_enter = NAN;
        check("cfg_nan_energy", sleep_config_valid(&bad) == 0);
        bad = c; bad.e_exit = INFINITY;
        check("cfg_inf_energy", sleep_config_valid(&bad) == 0);

        SlState st;
        SlSleepConfig neg = c;
        neg.e_enter = -1.0f;
        check("init_rejects_negative", sleep_init(&st, 1, &neg) == 0);
    }

    /* ---- at rest => sleeps within the hold window ---- */
    {
        SlSleepConfig cfg = test_cfg();
        SlState st;
        check("init_ok", sleep_init(&st, 1, &cfg) == 1);
        SlBody b = make_body(0, 0, 0, 1, 1);
        SlTick t;
        int slept_at = -1;
        int all_ticks_ok = 1;
        for (int f = 0; f < 80; f++) {
            if (!sleep_tick(&b, 1, &st, &t)) all_ticks_ok = 0;
            if (t.sleeping[0] && slept_at < 0) slept_at = f;
        }
        check("tick_ok", all_ticks_ok == 1);
        check("rest_sleeps", slept_at >= 0);
        check("rest_hold_ticks", st.hold_ticks == 32);
        /* tick 0 only seeds the history, so the deadline lands exactly on
           the hold tick rather than one past it */
        check("rest_hold_30", slept_at == st.hold_ticks);
        check("rest_awake_before_hold", slept_at != 0);
        check("rest_energy_zero", fabsf(t.energy_total) <= 0.0f &&
                                  fabsf(t.energy_mean) <= 0.0f);
        check("rest_bbox_valid", t.bbox[0].valid == 1);
        check("rest_bbox_lo", fabsf(t.bbox[0].lo.x + 1.0f) < 1e-6f &&
                               fabsf(t.bbox[0].lo.y + 1.0f) < 1e-6f);
        check("rest_bbox_hi", fabsf(t.bbox[0].hi.x - 1.0f) < 1e-6f &&
                              fabsf(t.bbox[0].hi.y - 1.0f) < 1e-6f);
    }

    /* ---- moving => never sleeps, energy gate ---- */
    {
        SlSleepConfig cfg = test_cfg();
        SlState st;
        sleep_init(&st, 1, &cfg);
        SlBody b = make_body(0, 0, 0, 1, 1);
        SlTick t;
        int any_sleep = 0;
        for (int f = 0; f < 90; f++) {
            b.vel = v3(5.0f, 0, 0);         /* E = 0.5*25 = 12.5 > 10 */
            b.pos.x = 0.2f * sinf(0.3f * (float)f);
            sleep_tick(&b, 1, &st, &t);
            if (t.sleeping[0]) any_sleep = 1;
        }
        check("moving_stays_awake", any_sleep == 0);
    }

    /* ---- still but jittering => never sleeps, variance gate ---- */
    {
        SlSleepConfig cfg = test_cfg();
        SlState st;
        sleep_init(&st, 1, &cfg);
        SlBody b = make_body(0, 0, 0, 1, 1);
        SlTick t;
        int any_sleep = 0;
        for (int f = 0; f < 90; f++) {
            b.vel = v3(0, 0, 0);
            b.pos.x = (f & 1) ? 20.0f : -20.0f;  /* V -> 200 > 100 */
            sleep_tick(&b, 1, &st, &t);
            if (t.sleeping[0]) any_sleep = 1;
        }
        check("jitter_blocks_sleep", any_sleep == 0);
    }

    /* ---- slow creep below both thresholds still sleeps ---- */
    {
        SlSleepConfig cfg = test_cfg();
        SlState st;
        sleep_init(&st, 1, &cfg);
        SlBody b = make_body(0, 0, 0, 1, 1);
        SlTick t;
        int slept_at = -1;
        for (int f = 0; f < 90; f++) {
            b.vel = v3(1.0f, 0, 0);              /* E = 0.5 < 10 */
            b.pos.x = sinf(0.2f * (float)f);     /* V ~ 0.5 < 100 */
            sleep_tick(&b, 1, &st, &t);
            if (t.sleeping[0] && slept_at < 0) slept_at = f;
        }
        check("slow_creep_sleeps", slept_at >= 0);
    }

    /* ---- hysteresis: band (enter, exit] neither wakes nor re-sleeps ---- */
    {
        SlSleepConfig cfg = test_cfg();
        SlState st;
        sleep_init(&st, 1, &cfg);
        SlBody b = make_body(0, 0, 0, 1, 1);
        SlTick t;

        int slept = 0;
        for (int f = 0; f < 60; f++) {
            sleep_tick(&b, 1, &st, &t);
            if (t.sleeping[0]) slept = 1;
        }
        check("hysteresis_reached_sleep", slept == 1);

        /* E = 15 sits between enter (10) and exit (20): must not wake */
        int woke = 0;
        for (int f = 0; f < 60; f++) {
            b.vel = v3(sqrtf(30.0f), 0, 0);
            sleep_tick(&b, 1, &st, &t);
            if (!t.sleeping[0]) woke = 1;
        }
        check("hysteresis_holds_through_band", woke == 0);

        /* E = 25 is above exit: must wake on the first such tick */
        b.vel = v3(sqrtf(50.0f), 0, 0);
        sleep_tick(&b, 1, &st, &t);
        check("hysteresis_wakes_above_exit", t.sleeping[0] == 0);

        /* a body that never got to sleep must not sleep at E = 15 either */
        SlState st2;
        sleep_init(&st2, 1, &cfg);
        SlBody b2 = make_body(0, 0, 0, 1, 1);
        int any_sleep = 0;
        for (int f = 0; f < 90; f++) {
            b2.vel = v3(sqrtf(30.0f), 0, 0);
            sleep_tick(&b2, 1, &st2, &t);
            if (t.sleeping[0]) any_sleep = 1;
        }
        check("hysteresis_blocks_flicker", any_sleep == 0);

        /* variance half of the hysteresis: a small +/-5 wobble stays inside
           the exit band, +/-22 (axis variance 484) exceeds it.  Positions sit
           at x = 1500 to prove the variance is origin-independent. */
        SlState st3;
        sleep_init(&st3, 1, &cfg);
        SlBody b3 = make_body(0, 0, 0, 1, 1);
        int slept3 = 0;
        for (int f = 0; f < 60; f++) {
            b3.pos.x = 1500.0f;
            sleep_tick(&b3, 1, &st3, &t);
            if (t.sleeping[0]) slept3 = 1;
        }
        check("var_hysteresis_reached_sleep", slept3 == 1);
        int woke3 = 0;
        for (int f = 0; f < 12; f++) {
            b3.pos.x = 1500.0f + ((f & 1) ? 5.0f : -5.0f);   /* V = 25 < 200 */
            sleep_tick(&b3, 1, &st3, &t);
            if (!t.sleeping[0]) woke3 = 1;
        }
        check("var_hysteresis_holds", woke3 == 0);
        for (int f = 0; f < 12; f++) {
            b3.pos.x = 1500.0f + ((f & 1) ? 22.0f : -22.0f); /* V = 484 > 200 */
            sleep_tick(&b3, 1, &st3, &t);
            if (!t.sleeping[0]) woke3 = 1;
        }
        check("var_hysteresis_wakes", woke3 == 1);
    }

    /* ---- AABB grows while asleep and resets on wake ---- */
    {
        SlSleepConfig cfg = test_cfg();
        SlState st;
        sleep_init(&st, 1, &cfg);
        SlBody b = make_body(0, 0, 0, 1, 2);
        SlTick t;
        for (int f = 0; f < 60; f++) sleep_tick(&b, 1, &st, &t);
        check("bbox_awake_before_sleep", t.bbox[0].valid == 1); /* now asleep */
        check("bbox_union_covers_rig", t.bounds.valid == 1 &&
              fabsf(t.bounds.hi.y - 2.0f) < 1e-6f);
        /* nudge the "sleeping" body far to the right; the box must follow */
        b.pos = v3(10.0f, 0, 0);
        sleep_tick(&b, 1, &st, &t);
        check("bbox_tracks_sleeping_body", t.bbox[0].valid == 1 &&
              fabsf(t.bbox[0].hi.x - 11.0f) < 1e-6f);
        check("bbox_keeps_history", fabsf(t.bbox[0].lo.x + 1.0f) < 1e-6f);
        check("bbox_union_follows", fabsf(t.bounds.lo.x + 1.0f) < 1e-6f &&
              fabsf(t.bounds.hi.x - 11.0f) < 1e-6f);
        /* waking clears validity; the next sleep rebuilds it */
        b.vel = v3(sqrtf(50.0f), 0, 0);
        sleep_tick(&b, 1, &st, &t);
        check("bbox_cleared_on_wake", t.sleeping[0] == 0 && t.bbox[0].valid == 0);
    }

    /* ---- variance window: one displaced sample out of N ---- */
    {
        SlSleepConfig cfg = test_cfg();
        SlState st;
        sleep_init(&st, 1, &cfg);
        SlBody b = make_body(0, 0, 0, 1, 1);
        SlTick t;
        for (int f = 0; f < SL_NSTEP; f++) sleep_tick(&b, 1, &st, &t);
        b.pos.x = 0.1f;                  /* exactly one sample displaced */
        sleep_tick(&b, 1, &st, &t);
        /* a single outlier gives (N-1)/N^2 * a^2 = 3/16 * 0.01 on one axis;
           the other two axes are still, so the trace equals that value */
        float expect = (3.0f / 16.0f) * 0.1f * 0.1f;
        check("variance_single_spike", fabsf(t.variance_mean - expect) < 1e-7f);
        check("variance_short_window_caps", st.ticks == SL_NSTEP);
    }

    /* ---- external wake entry point ---- */
    {
        SlSleepConfig cfg = test_cfg();
        SlState st;
        sleep_init(&st, 2, &cfg);
        SlBody b[2] = { make_body(0, 0, 0, 1, 1), make_body(3, 0, 0, 1, 1) };
        SlTick t;
        for (int f = 0; f < 60; f++) sleep_tick(b, 2, &st, &t);
        check("wake_precondition_asleep", t.sleeping[0] == 1 && t.sleeping[1] == 1);
        int wakes_before = t.wake_count;

        check("wake_ok", sleep_wake(&st, 1) == 1);
        check("wake_external_immediate", st.sleeping[1] == 0 && st.sleeping[0] == 1);
        check("wake_counts_transition", st.wake_count == wakes_before + 1);
        check("wake_clears_bbox", st.aabb[1].valid == 0);

        /* waking an awake body is a no-op, not an error */
        check("wake_idempotent", sleep_wake(&st, 1) == 1 &&
              st.wake_count == wakes_before + 1);
        check("wake_rejects_bad_index", sleep_wake(&st, 2) == 0 &&
              sleep_wake(&st, -1) == 0);

        /* a woken body must re-earn its hold window */
        int reslept = 0;
        for (int f = 0; f < 60; f++) {
            sleep_tick(b, 2, &st, &t);
            if (t.sleeping[1]) reslept = 1;
        }
        check("wake_reacquires_sleep", reslept == 1);
    }

    /* ---- body input rejection ---- */
    {
        SlSleepConfig cfg = test_cfg();
        SlState st;
        sleep_init(&st, 2, &cfg);
        SlBody b[2] = { make_body(0, 0, 0, 1, 1), make_body(5, 0, 0, 1, 1) };
        SlTick t;

        /* means are body-count independent; totals are not */
        b[0].vel = v3(sqrtf(2.0f * 12.5f), 0, 0);   /* E0 = 12.5 */
        b[1].vel = v3(sqrtf(2.0f * 7.5f), 0, 0);    /* E1 = 7.5  */
        sleep_tick(b, 2, &st, &t);
        check("energy_mean_two_bodies", fabsf(t.energy_mean - 10.0f) < 1e-4f);
        check("energy_total_two_bodies", fabsf(t.energy_total - 20.0f) < 1e-4f);

        SlBody nanb = make_body(0, 0, 0, 1, 1);
        nanb.vel = v3(NAN, 0, 0);
        check("tick_rejects_nan_vel", sleep_tick(&nanb, 1, &st, &t) == 0);
        nanb = make_body(0, 0, 0, 1, 1);
        nanb.pos = v3(INFINITY, 0, 0);
        check("tick_rejects_inf_pos", sleep_tick(&nanb, 1, &st, &t) == 0);
        nanb = make_body(0, 0, 0, 1, 1);
        nanb.inv_mass = 0.0f;
        check("tick_rejects_zero_mass", sleep_tick(&nanb, 1, &st, &t) == 0);
        nanb = make_body(0, 0, 0, -1, 1);
        check("tick_rejects_negative_half", sleep_tick(&nanb, 1, &st, &t) == 0);
        SlBody good = make_body(0, 0, 0, 1, 1);
        check("tick_rejects_wrong_count", sleep_tick(&good, 1, &st, &t) == 0);
    }

    /* ---- simulated wake latency: energy held in band for 3 frames ---- */
    {
        SlSleepConfig cfg = test_cfg();
        SlState st;
        sleep_init(&st, 1, &cfg);
        SlBody b = make_body(0, 0, 0, 1, 1);
        SlTick t;
        unsigned char pred[80];
        int s0 = -1, wake = -1;
        int W = 50;

        for (int f = 0; f < 80; f++) {
            if (f >= W && f <= W + 2) b.vel = v3(sqrtf(30.0f), 0, 0); /* E=15 */
            else if (f > W + 2)      b.vel = v3(sqrtf(50.0f), 0, 0); /* E=25 */
            else                     b.vel = v3(0, 0, 0);

            sleep_tick(&b, 1, &st, &t);
            pred[f] = t.sleeping[0];
            if (s0 < 0 && t.sleeping[0]) s0 = f;
            if (wake < 0 && s0 >= 0 && f > s0 && !t.sleeping[0]) wake = f;
        }
        check("lat_slept_before_event", s0 >= 0 && s0 < W);
        check("lat_woke_three_frames_late", wake == W + 3);

        SlInterval truth[1] = { { s0, W } };
        SlEval ev;
        SlEvalConfig ec = sleep_eval_config_default();
        ec.fa_frac = 0.05f;
        ec.miss_frac = 0.20f;
        ec.lat_frac = 0.0f;
        ec.max_latency = 3;
        check("lat_eval_ok", sleep_evaluate(pred, 80, truth, 1, &ec, &ev) == 1);
        check("lat_no_false_alarm", ev.fa_frames == 0);
        check("lat_reports_worst", ev.worst_latency == 3);
        check("lat_at_cap_passes", ev.wake_exceeded == 0 && ev.ok == 1);

        ec.max_latency = 2;
        sleep_evaluate(pred, 80, truth, 1, &ec, &ev);
        check("lat_over_cap_flagged", ev.wake_exceeded == 1 && ev.ok == 0);

        ec.max_latency = 2;
        ec.lat_frac = 1.0f;
        sleep_evaluate(pred, 80, truth, 1, &ec, &ev);
        check("lat_over_cap_reported_not_failed", ev.wake_exceeded == 1 &&
              ev.ok == 1);

        /* never-wake latency: hold energy inside the hysteresis band after
           the artist's wake frame, so the detector stays asleep */
        SlState stn;
        sleep_init(&stn, 1, &cfg);
        SlBody bn = make_body(0, 0, 0, 1, 1);
        unsigned char predn[80];
        int s0n = -1;
        for (int f = 0; f < 80; f++) {
            bn.vel = (f >= W) ? v3(sqrtf(30.0f), 0, 0) : v3(0, 0, 0);
            sleep_tick(&bn, 1, &stn, &t);
            predn[f] = t.sleeping[0];
            if (s0n < 0 && t.sleeping[0]) s0n = f;
        }
        SlInterval truthn[1] = { { s0n, W } };
        ec = sleep_eval_config_default();
        ec.lat_frac = 0.0f;
        ec.max_latency = 2;
        ec.fa_frac = 0.0f;
        ec.miss_frac = 1.0f;      /* not the failure under test */
        sleep_evaluate(predn, 80, truthn, 1, &ec, &ev);
        check("lat_never_wakes_large", ev.worst_latency == 80 - W &&
              ev.wake_exceeded == 1 && ev.ok == 0);
    }

    /* ---- synthetic scoring: exact counts and upper bounds ---- */
    {
        enum { NF = 20 };
        unsigned char pred[NF];
        /* pred: asleep 0..7, awake 8..19.  Truth: asleep 0..4 and 7..8. */
        memset(pred, 1, 8);
        memset(pred + 8, 0, NF - 8);

        SlInterval truth[] = {
            { 0, 5 },       /* frames 0..4 asleep */
            { 7, 9 },       /* frames 7..8 asleep */
            { -3, -1 },     /* clips to empty */
            { 20, 25 },     /* clips to empty (start == nframes) */
            { 12, 10 },     /* malformed, contributes nothing */
        };
        SlEvalConfig ec = sleep_eval_config_default();
        ec.fa_frac = 0.20f;     /* bound = ceil(0.2 * 7) = 2 */
        ec.miss_frac = 0.20f;   /* bound = ceil(0.2 * 13) = 3 */
        ec.lat_frac = 0.50f;    /* bound = ceil(0.5 * 2) = 1 */
        ec.max_latency = 2;

        SlEval ev;
        check("eval_ok", sleep_evaluate(pred, NF, truth, 5, &ec, &ev) == 1);
        check("eval_truth_counts", ev.truth_asleep == 7 && ev.truth_awake == 13);
        check("eval_false_alarm_count", ev.fa_frames == 1 && ev.fa_bound == 2);
        check("eval_miss_count", ev.miss_frames == 2 && ev.miss_bound == 3);
        check("eval_wake_runs", ev.wake_runs == 2);
        check("eval_worst_latency", ev.worst_latency == 3);
        check("eval_wake_within_bound", ev.wake_exceeded == 1 && ev.wake_bound == 1);
        check("eval_passes", ev.ok == 1);

        /* tightening the limits must flip ok, and bounds are reported */
        ec.lat_frac = 0.0f;
        sleep_evaluate(pred, NF, truth, 5, &ec, &ev);
        check("eval_latency_fails", ev.wake_bound == 0 && ev.wake_exceeded == 1 &&
              ev.ok == 0);
        ec.lat_frac = 0.50f;
        ec.fa_frac = 0.0f;
        sleep_evaluate(pred, NF, truth, 5, &ec, &ev);
        check("eval_false_alarm_fails", ev.fa_bound == 0 && ev.fa_frames == 1 &&
              ev.ok == 0);
        ec.fa_frac = 0.20f;
        ec.miss_frac = 0.0f;
        sleep_evaluate(pred, NF, truth, 5, &ec, &ev);
        check("eval_miss_fails", ev.miss_bound == 0 && ev.miss_frames == 2 &&
              ev.ok == 0);

        /* detector that never sleeps: alarm bound trips, latency is zero */
        unsigned char never[NF];
        memset(never, 0, sizeof(never));
        ec = sleep_eval_config_default();
        ec.fa_frac = 0.20f; ec.miss_frac = 0.20f; ec.lat_frac = 0.50f;
        ec.max_latency = 2;
        sleep_evaluate(never, NF, truth, 5, &ec, &ev);
        check("eval_never_sleeps_flagged", ev.fa_frames == 7 && ev.ok == 0);
        check("eval_never_sleeps_latency_zero", ev.worst_latency == 0 &&
              ev.wake_runs == 2);

        /* detector that never wakes: miss bound trips, latency is large */
        unsigned char never_wake[NF];
        memset(never_wake, 1, sizeof(never_wake));
        sleep_evaluate(never_wake, NF, truth, 5, &ec, &ev);
        check("eval_never_wakes_flagged", ev.miss_frames == 13 && ev.ok == 0);
        check("eval_never_wakes_latency_large", ev.worst_latency == 15 &&
              ev.wake_exceeded == 2);

        /* invalid scoring configuration is rejected */
        ec.fa_frac = 1.5f;
        check("eval_rejects_bad_frac", sleep_evaluate(pred, NF, truth, 5, &ec, &ev) == 0);
        ec = sleep_eval_config_default();
        ec.max_latency = 0;
        check("eval_rejects_zero_cap", sleep_evaluate(pred, NF, truth, 5, &ec, &ev) == 0);
        SlEvalConfig dflt = sleep_eval_config_default();
        check("eval_rejects_zero_frames",
              sleep_evaluate(pred, 0, truth, 5, &dflt, &ev) == 0);
        check("eval_rejects_no_pred",
              sleep_evaluate(NULL, NF, truth, 5, &dflt, &ev) == 0);
        check("eval_rejects_null_truth_with_count",
              sleep_evaluate(pred, NF, NULL, 5, &dflt, &ev) == 0);

        /* zero labels is legal: every predicted sleep is a miss */
        SlEvalConfig loose = sleep_eval_config_default();
        loose.miss_frac = 1.0f;
        check("eval_zero_labels_ok",
              sleep_evaluate(pred, NF, NULL, 0, &loose, &ev) == 1 &&
              ev.truth_asleep == 0 && ev.miss_frames == 8);

        /* frame count at the cap is accepted, one past it is rejected */
        static unsigned char cap_pred[SL_MAXFRAME];
        memset(cap_pred, 0, sizeof(cap_pred));
        check("eval_max_frames_ok",
              sleep_evaluate(cap_pred, SL_MAXFRAME, NULL, 0, &loose, &ev) == 1);
        check("eval_over_max_frames_rejected",
              sleep_evaluate(cap_pred, SL_MAXFRAME + 1, NULL, 0, &loose, &ev) == 0);

        /* intervals clip at both ends and overlapping ones are unioned */
        unsigned char on[NF];
        memset(on, 1, sizeof(on));
        SlInterval clipped[] = { { -5, 1 }, { -3, -1 }, { 18, 99 }, { 18, 25 } };
        check("eval_clip_both_ends",
              sleep_evaluate(on, NF, clipped, 4, &loose, &ev) == 1 &&
              ev.truth_asleep == 3 && ev.fa_frames == 0);
    }

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
