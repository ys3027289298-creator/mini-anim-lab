
/* lod: constraint level-of-detail for a ragdoll-style distance-constraint
   chain.

   Reduction. Keep every `stride`-th particle as an *active* node. Active
   nodes are integrated and tied by tethers whose target length is the rope
   length of the full edges they replace. Skipped particles are reconstructed
   each frame by rest-length interpolation along the active chain. Cost:
       nodes = ceil((n-1)/stride) + 1,  cost ~ (nodes-1) * iterations
   The iteration budget is part of the tier, so a tier cuts both constraint
   count and solve iterations.

   Cull tier. stride 0 keeps every particle simulated but drops all tethers:
   the chain free-falls. It is a legal tier, reported explicitly as culled
   and *unbounded*. It is never silently treated as equivalent to a solved
   tier.

   Switch continuity.
     - lod_switch() never touches particle state. At the switch instant the
       pose and Verlet previous positions are bit-identical: no teleport.
     - a tether that survives a switch keeps its goal; a new tether starts at
       its current endpoint distance, so the first solve after the switch
       applies no correction. This is what keeps a fine<->coarse switch (and
       returning from a long cull) continuous.
     - particles that become active keep a valid Verlet previous position
       (their last reconstructed position), so their velocity does not reset.
     - a particle that becomes inactive starts a short blend toward the
       interpolated chord, beginning at weight 0 so the switch frame itself
       moves nothing. A particle already inactive follows the moving chord
       exactly, so reconstruction adds no motion lag.
     - a switch sets `engaging`, and lod_step then clamps each active
       particle's per-frame constraint displacement to ENGAGE_RATE*dt until
       the frame's correction fits comfortably under the cap. That is the
       hard guarantee behind "no jump": at most ENGAGE_RATE*dt of
       constraint-induced motion per frame, on every frame. A large
       correction (re-closing edges after a long cull) is spread over as
       many frames as it needs instead of popping in one. The clamp is
       released with headroom, so releasing it is not a jump either.
       `last_correction` exposes the number per frame.
     - band_select() has hysteresis, so a camera parked on a band boundary
       does not flip tiers every frame.

   Error reporting. `bound` is a declared, benchmark-scoped upper bound on
   the pose error of a band against the full solver, checked every run:
     - settle_error: after a gravity settle, max particle distance between
       the reduced pose and the full pose. This is the well-posed steady-
       state number; it is what `settle_bound` constrains.
     - sway_error: the same comparison while the anchor swings. Reported as
       data only: a swinging chain is chaotic, so long free runs diverge
       regardless of approximation quality and cannot carry a tight bound.
     - step_deviation: worst per-frame particle distance between one reduced
       step and one full step, both from the identical state. This is the
       local error a switch can inject and does not accumulate; it is what
       `step_bound` constrains. */

#include <stdio.h>
#include <math.h>
#include <string.h>

#define MAX_P 16
#define MAX_CONS (MAX_P - 1)
#define MAX_TIERS 4

typedef struct { float x, y, px, py; float inv_mass; int pinned; } Particle;

typedef struct {
    Particle p[MAX_P];
    int n;
    float rest[MAX_P];        /* rest length of full edge i -> i+1 */
    float cum[MAX_P];         /* cumulative rest length at particle i */
    float damping;
} Chain;

typedef struct { int a, b; float goal; float target; } Con;

typedef struct {
    Con cons[MAX_CONS];
    int n_cons;
    int iterations;

    int active[MAX_P];
    int n_active;
    int stride;               /* 0 = cull tier */

    float blend[MAX_P];       /* inactive particles: 0 = still at its
                                 pre-switch position, 1 = fully on the chord */
    float rec_px[MAX_P];      /* last frame's reconstructed position, for
                                 inactive particles. Kept separately from the
                                 Particle's Verlet history so a particle's
                                 velocity survives being dropped and later
                                 reactivated. */
    float rec_py[MAX_P];

    int engaging;             /* 1 while a switch's correction is still being
                                 eased in. lod_step then clamps each active
                                 particle's per-frame constraint displacement
                                 to engage_rate*dt, so a tier change never
                                 moves a particle further than that on any
                                 frame while engaging. Cleared once the
                                 unclamped correction fits under half the cap,
                                 which leaves headroom so releasing the clamp
                                 is itself not a step. */
    float last_correction;    /* largest constraint-induced particle
                                 displacement in the most recent step. While
                                 engaging it is <= engage_rate*dt by
                                 construction; the self-test asserts that. */
} Lod;

/* A distance band: tier parameters, the camera rule that selects it, and the
   declared pose-error bounds. A negative bound means "unbounded" and is
   reported as such. */
typedef struct {
    int stride;
    int iterations;
    float max_dist;           /* eligible while camera distance <= this */
    float min_cover;          /* ... and screen coverage >= this */
    float settle_bound;       /* declared max pose error vs full, settled */
    float step_bound;         /* declared max per-frame step deviation */
} Band;

static const char* band_name(int i) {
    static const char* names[MAX_TIERS] = { "near", "mid", "far", "cull" };
    return (i >= 0 && i < MAX_TIERS) ? names[i] : "?";
}

/* ---- chain / solver ---------------------------------------------------- */

static void chain_init(Chain* c, int n, float seg, float x0, float y0) {
    memset(c, 0, sizeof(*c));
    c->n = n;
    c->damping = 0.99f;
    for (int i = 0; i < n; i++) {
        c->p[i].x = x0 + seg * i; c->p[i].y = y0;
        c->p[i].px = c->p[i].x; c->p[i].py = c->p[i].y;
        c->p[i].inv_mass = 1.0f;
        c->p[i].pinned = 0;
        c->rest[i] = seg;
        c->cum[i] = seg * (float)i;
    }
    c->p[0].pinned = 1;
    c->p[0].inv_mass = 0.0f;
}

static float pdist(const Particle* a, const Particle* b) {
    float dx = a->x - b->x, dy = a->y - b->y;
    return sqrtf(dx * dx + dy * dy);
}

/* Solve one distance correction. */
static void solve_constraint(Particle* a, Particle* b, float goal) {
    float dx = b->x - a->x, dy = b->y - a->y;
    float d = sqrtf(dx * dx + dy * dy);
    if (d < 1e-6f) return;
    float diff = (d - goal) / d;
    float wa = a->inv_mass, wb = b->inv_mass;
    float wsum = wa + wb;
    if (wsum <= 0) return;
    a->x += dx * diff * (wa / wsum);
    a->y += dy * diff * (wa / wsum);
    b->x -= dx * diff * (wb / wsum);
    b->y -= dy * diff * (wb / wsum);
}

/* ---- tier construction ------------------------------------------------- */

static int tier_active_count(int n, int stride) {
    if (stride <= 1) return n;
    int k = 0;
    for (int i = 0; i < n; i += stride) k++;
    if ((n - 1) % stride != 0) k++;      /* final particle always active */
    return k;
}

static Con* find_con(Lod* lod, int a, int b) {
    for (int i = 0; i < lod->n_cons; i++)
        if (lod->cons[i].a == a && lod->cons[i].b == b) return &lod->cons[i];
    return NULL;
}

/* Canonical full solver: every edge constrained at its rest length, no goal
   relaxation, no reconstruction. Used as the reference and as the control
   branch of a step-deviation probe. */
static void lod_full(Lod* lod, const Chain* c, int iterations) {
    memset(lod, 0, sizeof *lod);
    lod->stride = 1;
    lod->iterations = iterations;
    lod->n_active = c->n;
    lod->engaging = 0;
    for (int i = 0; i < c->n; i++) {
        lod->active[i] = 1;
        lod->blend[i] = 1.0f;
        if (i < c->n - 1) {
            lod->cons[lod->n_cons].a = i;
            lod->cons[lod->n_cons].b = i + 1;
            lod->cons[lod->n_cons].goal = c->rest[i];
            lod->cons[lod->n_cons].target = c->rest[i];
            lod->n_cons++;
        }
    }
}

/* Change tier without touching pose state. Surviving tethers keep their
   goal; new tethers start at the current endpoint distance. Particles that
   become inactive start a blend at weight 0; particles already inactive keep
   following the chord. Every switch starts `engaging`, so the first frames
   cap how far any particle may be moved by the constraint solve. */
static void lod_switch(const Chain* c, Lod* lod, int stride, int iterations) {
    Lod old = *lod;
    int have_old = old.n_active > 0;

    memset(lod, 0, sizeof *lod);
    lod->stride = stride;
    lod->iterations = iterations;
    /* Engage only when there is a previous tier to stay continuous with.
       A fresh setup has no prior pose, so there is nothing to pop from and
       the first frame solves unmodified. */
    lod->engaging = have_old ? 1 : 0;

    if (stride <= 0) {
        /* cull: simulate everything, tie nothing */
        for (int i = 0; i < c->n; i++) lod->active[i] = 1;
    } else if (stride == 1) {
        for (int i = 0; i < c->n; i++) lod->active[i] = 1;
    } else {
        for (int i = 0; i < c->n; i += stride) lod->active[i] = 1;
        lod->active[c->n - 1] = 1;
    }
    lod->n_active = 0;
    for (int i = 0; i < c->n; i++) lod->n_active += lod->active[i];

    if (stride > 0) {
        int prev = -1;
        for (int i = 0; i < c->n; i++) {
            if (!lod->active[i]) continue;
            if (prev >= 0) {
                Con* k = &lod->cons[lod->n_cons++];
                k->a = prev;
                k->b = i;
                k->target = c->cum[i] - c->cum[prev];
                Con* live = have_old ? find_con(&old, prev, i) : NULL;
                k->goal = live ? live->goal : pdist(&c->p[prev], &c->p[i]);
            }
            prev = i;
        }
    }

    for (int i = 0; i < c->n; i++) {
        if (lod->active[i]) {
            lod->blend[i] = 1.0f;
            lod->rec_px[i] = c->p[i].px;
            lod->rec_py[i] = c->p[i].py;
        } else if (have_old && old.active[i]) {
            /* just dropped: ease onto the chord, but keep the velocity the
               particle had while active (its Verlet previous position) */
            lod->blend[i] = 0.0f;
            lod->rec_px[i] = c->p[i].px;
            lod->rec_py[i] = c->p[i].py;
        } else if (have_old) {
            lod->blend[i] = old.blend[i];
            lod->rec_px[i] = old.rec_px[i];
            lod->rec_py[i] = old.rec_py[i];
        } else {
            lod->blend[i] = 0.0f;
            lod->rec_px[i] = c->p[i].x;
            lod->rec_py[i] = c->p[i].y;
        }
    }
}

/* Move every tether goal toward its target by at most rate*dt. */
static void lod_relax(Lod* lod, float rate, float dt) {
    float lim = rate * dt;
    for (int i = 0; i < lod->n_cons; i++) {
        Con* k = &lod->cons[i];
        float d = k->target - k->goal;
        if (d > lim) d = lim;
        if (d < -lim) d = -lim;
        k->goal += d;
    }
}

/* Reconstruct inactive particles from the active chain by rest-length
   interpolation, weighted by blend. The weight is applied before it is
   advanced, so the switch frame itself moves nothing. Verlet previous
   positions are maintained for inactive particles too, so reactivating them
   on a finer tier preserves velocity. */
static void lod_reconstruct(Chain* c, Lod* lod, float blend_rate, float dt) {
    float step = blend_rate * dt;
    if (step > 1.0f) step = 1.0f;
    int prev = -1;
    for (int i = 0; i < c->n; i++) {
        if (lod->active[i]) { prev = i; continue; }
        int next = i + 1;
        while (next < c->n && !lod->active[next]) next++;
        if (prev < 0 || next >= c->n) continue;
        float span = c->cum[next] - c->cum[prev];
        float t = span > 0.0f ? (c->cum[i] - c->cum[prev]) / span : 0.0f;
        float cx = c->p[prev].x + (c->p[next].x - c->p[prev].x) * t;
        float cy = c->p[prev].y + (c->p[next].y - c->p[prev].y) * t;
        float w = lod->blend[i];
        c->p[i].px = lod->rec_px[i];   /* velocity from the last real frame */
        c->p[i].py = lod->rec_py[i];
        c->p[i].x += (cx - c->p[i].x) * w;
        c->p[i].y += (cy - c->p[i].y) * w;
        lod->rec_px[i] = c->p[i].x;
        lod->rec_py[i] = c->p[i].y;
        lod->blend[i] = w + step < 1.0f ? w + step : 1.0f;
    }
}

static void lod_step(Chain* c, Lod* lod, float dt, float gx, float gy,
                     float relax_rate, float blend_rate, float engage_rate) {
    for (int i = 0; i < c->n; i++) {
        Particle* p = &c->p[i];
        if (!lod->active[i] || p->pinned) continue;
        float vx = (p->x - p->px) * c->damping;
        float vy = (p->y - p->py) * c->damping;
        p->px = p->x; p->py = p->y;
        p->x += vx + gx * dt * dt;
        p->y += vy + gy * dt * dt;
    }
    lod_relax(lod, relax_rate, dt);

    float bx[MAX_P], by[MAX_P];
    for (int i = 0; i < c->n; i++) { bx[i] = c->p[i].x; by[i] = c->p[i].y; }
    for (int it = 0; it < lod->iterations; it++)
        for (int i = 0; i < lod->n_cons; i++)
            solve_constraint(&c->p[lod->cons[i].a], &c->p[lod->cons[i].b],
                             lod->cons[i].goal);

    float worst = 0.0f;
    for (int i = 0; i < c->n; i++) {
        if (!lod->active[i]) continue;
        float dx = c->p[i].x - bx[i], dy = c->p[i].y - by[i];
        float m = sqrtf(dx * dx + dy * dy);
        if (m > worst) worst = m;
    }
    if (lod->engaging) {
        /* Clamp this frame's constraint correction so no active particle
           moves further than engage_rate*dt because of it. The switch frame
           itself starts from a zero correction, so a switch can never pop.
           Engagement ends once the *unclamped* correction fits under half
           the cap, leaving headroom so releasing the clamp does not produce
           a step of its own. */
        float cap = engage_rate * dt;
        if (worst > 0.5f * cap && worst > 1e-9f) {
            float s = worst > cap ? cap / worst : 1.0f;
            for (int i = 0; i < c->n; i++) {
                if (!lod->active[i]) continue;
                c->p[i].x = bx[i] + (c->p[i].x - bx[i]) * s;
                c->p[i].y = by[i] + (c->p[i].y - by[i]) * s;
            }
            worst *= s;
        } else {
            lod->engaging = 0;
        }
    }
    lod->last_correction = worst;
    lod_reconstruct(c, lod, blend_rate, dt);
}

/* Worst full-edge stretch at the current pose. Diagnostic only. */
static float chain_seg_error(const Chain* c) {
    float worst = 0.0f;
    for (int i = 0; i < c->n - 1; i++) {
        float e = fabsf(pdist(&c->p[i], &c->p[i + 1]) - c->rest[i]);
        if (e > worst) worst = e;
    }
    return worst;
}

/* ---- distance / screen selection --------------------------------------- */

/* Fraction of viewport height covered by a world-space length at distance
   `dist` with vertical half-FOV `half_fov`. Clamped to [0,1]. */
static float screen_coverage(float world_len, float dist, float half_fov) {
    if (dist < 1e-4f) return 1.0f;
    float half_h = dist * tanf(half_fov);
    if (half_h <= 1e-6f) return 1.0f;
    float frac = world_len / (2.0f * half_h);
    if (frac < 0.0f) frac = 0.0f;
    if (frac > 1.0f) frac = 1.0f;
    return frac;
}

/* World units per screen pixel at `dist`, for a viewport `view_px` high and
   vertical half-FOV `half_fov`. Used to express a band's error bound in
   pixels at the distance where it is selected. */
static float world_per_pixel(float dist, float half_fov, float view_px) {
    float half_h = dist * tanf(half_fov);
    return 2.0f * half_h / view_px;
}

/* Finest band whose distance and coverage conditions both hold, else the
   coarsest band (the cull band if it is last). `current` is the band in use
   now, or -1 to select from scratch. Promotion needs the finer band's
   condition to hold with a margin; demotion needs the current band's own
   condition to be violated. The gap between the two is hysteresis, so a
   camera sitting on a boundary does not flip tiers every frame. */
static int band_select(const Band* bands, int n, float dist, float cover,
                       int current, float hysteresis) {
    int want = n - 1;
    for (int i = 0; i < n; i++)
        if (dist <= bands[i].max_dist && cover >= bands[i].min_cover) {
            want = i;
            break;
        }
    if (current < 0 || current >= n) return want;
    if (want < current) {
        /* promote only once clearly inside the finer band's range */
        float lim_d = bands[want].max_dist * (1.0f - hysteresis);
        float lim_c = bands[want].min_cover * (1.0f + hysteresis);
        int clear = dist <= lim_d && cover >= lim_c;
        return clear ? want : current;
    }
    if (want > current) {
        /* demote only once clearly outside the current band's range */
        float lim_d = bands[current].max_dist * (1.0f + hysteresis);
        float lim_c = bands[current].min_cover * (1.0f - hysteresis);
        int outside = dist > lim_d || cover < lim_c;
        return outside ? want : current;
    }
    return want;
}

/* ---- metrics ----------------------------------------------------------- */

static float pose_error_mean(const Chain* got, const Chain* ref) {
    if (got->n != ref->n) return 1e30f;
    float sum = 0.0f;
    for (int i = 0; i < got->n; i++) sum += pdist(&got->p[i], &ref->p[i]);
    return sum / (float)got->n;
}

static float pose_error_max(const Chain* got, const Chain* ref) {
    if (got->n != ref->n) return 1e30f;
    float worst = 0.0f;
    for (int i = 0; i < got->n; i++) {
        float d = pdist(&got->p[i], &ref->p[i]);
        if (d > worst) worst = d;
    }
    return worst;
}

static float pos_delta(const Chain* a, const Chain* b) {
    float worst = 0.0f;
    for (int i = 0; i < a->n; i++) {
        float d = pdist(&a->p[i], &b->p[i]);
        if (d > worst) worst = d;
    }
    return worst;
}

/* One frame of the reduced tier vs one frame of the full solver, both
   started from `state`. Local, non-accumulating error. */
static float step_deviation(const Chain* state, const Lod* reduced,
                            float dt, float gx, float gy,
                            float relax_rate, float blend_rate,
                            float engage_rate, int full_iters) {
    Chain a = *state, b = *state;
    Lod la = *reduced, lb;
    lod_full(&lb, &b, full_iters);
    lod_step(&a, &la, dt, gx, gy, relax_rate, blend_rate, engage_rate);
    lod_step(&b, &lb, dt, gx, gy, relax_rate, blend_rate, engage_rate);
    return pos_delta(&a, &b);
}

/* Worst per-edge stretch of the reduced chain; shows how well the kept
   constraints hold the shape. Diagnostic only. */
static float lod_seg_error(const Chain* c, const Lod* lod) {
    float worst = 0.0f;
    for (int i = 0; i < lod->n_cons; i++) {
        float e = fabsf(pdist(&c->p[lod->cons[i].a], &c->p[lod->cons[i].b]) -
                        lod->cons[i].goal);
        if (e > worst) worst = e;
    }
    return worst;
}

/* How far each skipped particle of a given pose sits from the chord between
   its surrounding active nodes: the geometry-only floor of a band's error.
   Reported to separate it from solver drift. */
static float interp_error(const Chain* ref, int stride) {
    if (stride <= 1) return 0.0f;
    int active[MAX_P] = { 0 };
    for (int i = 0; i < ref->n; i += stride) active[i] = 1;
    active[ref->n - 1] = 1;
    float worst = 0.0f;
    int prev = -1;
    for (int i = 0; i < ref->n; i++) {
        if (active[i]) { prev = i; continue; }
        int next = i + 1;
        while (next < ref->n && !active[next]) next++;
        if (prev < 0 || next >= ref->n) continue;
        float span = ref->cum[next] - ref->cum[prev];
        float t = span > 0.0f ? (ref->cum[i] - ref->cum[prev]) / span : 0.0f;
        float lx = ref->p[prev].x + (ref->p[next].x - ref->p[prev].x) * t;
        float ly = ref->p[prev].y + (ref->p[next].y - ref->p[prev].y) * t;
        float d = sqrtf((ref->p[i].x - lx) * (ref->p[i].x - lx) +
                        (ref->p[i].y - ly) * (ref->p[i].y - ly));
        if (d > worst) worst = d;
    }
    return worst;
}

/* ---- self test --------------------------------------------------------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

#define N_PARTICLES 10
#define SEG_LEN 1.0f
#define DT 0.016f
#define GRAV_Y (-9.8f)
#define RELAX_RATE 4.0f     /* tether goal change, world units/second */
#define BLEND_RATE 12.0f    /* inactive-particle blend, per second */
#define ENGAGE_RATE 4.5f    /* switch correction cap, world units/second */
#define FULL_ITERS 16
#define FRAMES 240
#define HALF_FOV 0.5f
#define VIEW_PX 1080.0f
#define HYS 0.10f           /* band-selection hysteresis, fraction of the edge */
#define FULL_RATES RELAX_RATE, BLEND_RATE, ENGAGE_RATE

/* Vertical chain under gravity; the anchor is at the origin and pinned. */
static void init_hang(Chain* c) {
    chain_init(c, N_PARTICLES, SEG_LEN, 0.0f, 0.0f);
    for (int i = 0; i < N_PARTICLES; i++) {
        c->p[i].x = 0.0f;
        c->p[i].y = -SEG_LEN * (float)i;
        c->p[i].px = c->p[i].x;
        c->p[i].py = c->p[i].y;
    }
}

static void drive_anchor(Chain* c, int frame) {
    float t = (float)frame * DT;
    c->p[0].x = 0.6f * sinf(t * 1.2f);
    c->p[0].y = 0.0f;
}

int main(void) {
    const Band bands[MAX_TIERS] = {
        /* stride iters  max_dist min_cover  settle_bound step_bound */
        {  1, FULL_ITERS, 6.0f,  0.30f,       0.000001f,  0.000001f },
        {  2,  8,         20.0f, 0.10f,       0.010000f,  0.010000f },
        {  3,  4,         60.0f, 0.02f,       0.020000f,  0.020000f },
        {  0,  0,         1e9f,  0.00f,      -1.0f,      -1.0f }  /* cull */
    };
    check("band_count", MAX_TIERS == 4);

    /* ---- settle benchmark: declared bounds checked here ---- */
    Chain ref;
    init_hang(&ref);
    Lod ref_lod;
    lod_full(&ref_lod, &ref, FULL_ITERS);
    for (int f = 0; f < FRAMES; f++)
        lod_step(&ref, &ref_lod, DT, 0, GRAV_Y, FULL_RATES);
    check("ref_stays_connected", chain_seg_error(&ref) < 0.05f);

    float settle_err[MAX_TIERS];
    for (int t = 0; t < MAX_TIERS; t++) {
        Chain c;
        init_hang(&c);
        Lod lod;
        memset(&lod, 0, sizeof lod);
        lod_switch(&c, &lod, bands[t].stride, bands[t].iterations);
        for (int f = 0; f < FRAMES; f++)
            lod_step(&c, &lod, DT, 0, GRAV_Y, FULL_RATES);

        float settle_max = pose_error_max(&c, &ref);
        float settle_mean = pose_error_mean(&c, &ref);
        settle_err[t] = settle_max;

        char name[64];
        snprintf(name, sizeof name, "band_%s_nodes", band_name(t));
        check(name, lod.n_active == tier_active_count(N_PARTICLES, bands[t].stride));
        snprintf(name, sizeof name, "band_%s_cons", band_name(t));
        check(name, lod.n_cons == (bands[t].stride ? lod.n_active - 1 : 0));

        if (bands[t].settle_bound < 0.0f) {
            /* unbounded band: not silently passed; it is reported as such */
            printf("band_%s stride=%d iters=%d nodes=%d cons=%d cost=%d "
                   "settle_mean=%.6f settle_max=%.6f bound=unbounded\n",
                   band_name(t), bands[t].stride, bands[t].iterations,
                   lod.n_active, lod.n_cons, lod.n_cons * bands[t].iterations,
                   settle_mean, settle_max);
            check("band_cull_marked_unbounded", bands[t].stride == 0);
        } else {
            int bound_ok = settle_max <= bands[t].settle_bound;
            snprintf(name, sizeof name, "band_%s_settle_within_bound", band_name(t));
            check(name, bound_ok);
            printf("band_%s stride=%d iters=%d nodes=%d cons=%d cost=%d "
                   "settle_mean=%.6f settle_max=%.6f bound=%.6f within=%d\n",
                   band_name(t), bands[t].stride, bands[t].iterations,
                   lod.n_active, lod.n_cons, lod.n_cons * bands[t].iterations,
                   settle_mean, settle_max, bands[t].settle_bound, bound_ok);
        }
    }
    /* settle error must grow with the reduction level, and each band's
       declared bound must sit above its own measured error. This is the
       per-distance-bound report the task asks for: a number per band plus
       the guard that the number is actually respected. */
    check("settle_error_monotone",
          settle_err[0] <= settle_err[1] && settle_err[1] <= settle_err[2]);
    for (int t = 0; t < MAX_TIERS; t++) {
        if (bands[t].settle_bound < 0.0f) continue;
        printf("report_band_%s distance<=%.1f settle_err=%.6f bound=%.6f "
               "bound_px=%.3f\n",
               band_name(t), bands[t].max_dist, settle_err[t],
               bands[t].settle_bound,
               bands[t].settle_bound /
                   world_per_pixel(bands[t].max_dist, HALF_FOV, VIEW_PX));
    }
    /* the near band is the full solver: it must reproduce the reference
       exactly (float equality; memcmp would trip on signed zero) */
    {
        Chain c;
        init_hang(&c);
        Lod lod;
        memset(&lod, 0, sizeof lod);
        lod_switch(&c, &lod, bands[0].stride, bands[0].iterations);
        for (int f = 0; f < FRAMES; f++)
            lod_step(&c, &lod, DT, 0, GRAV_Y, FULL_RATES);
        int same = lod.n_cons == ref_lod.n_cons && c.n == ref.n;
        for (int i = 0; same && i < c.n; i++)
            same = c.p[i].x == ref.p[i].x && c.p[i].y == ref.p[i].y &&
                   c.p[i].px == ref.p[i].px && c.p[i].py == ref.p[i].py;
        check("full_tier_bit_exact", same);
    }

    /* cull tier must be reported as culled and unbounded */
    {
        Chain c;
        init_hang(&c);
        Lod lod;
        memset(&lod, 0, sizeof lod);
        lod_switch(&c, &lod, 0, 0);
        for (int f = 0; f < FRAMES; f++)
            lod_step(&c, &lod, DT, 0, GRAV_Y, FULL_RATES);
        float e = pose_error_max(&c, &ref);
        printf("band_cull culled=1 cons=0 bound=inf "
               "freefall_y=%.3f err_max=%.6f explicitly_unbounded=1\n",
               c.p[N_PARTICLES - 1].y, e);
        check("zero_tier_active", lod.n_active == N_PARTICLES);
        check("zero_tier_no_constraints", lod.n_cons == 0);
        check("zero_tier_freefall", c.p[N_PARTICLES - 1].y < -20.0f);
        check("zero_tier_unbounded", e > 1.0f);
    }

    /* ---- sway benchmark: report pose drift and step deviation ---- */
    {
        Chain sref;
        init_hang(&sref);
        Lod sref_lod;
        lod_full(&sref_lod, &sref, FULL_ITERS);

        Chain sc[MAX_TIERS];
        Lod sl[MAX_TIERS];
        float dev_max[MAX_TIERS], err_mean[MAX_TIERS], err_max[MAX_TIERS];
        for (int t = 0; t < MAX_TIERS; t++) {
            init_hang(&sc[t]);
            memset(&sl[t], 0, sizeof sl[t]);
            lod_switch(&sc[t], &sl[t], bands[t].stride, bands[t].iterations);
            dev_max[t] = err_mean[t] = err_max[t] = 0.0f;
        }

        for (int f = 0; f < FRAMES; f++) {
            drive_anchor(&sref, f);
            lod_step(&sref, &sref_lod, DT, 0, GRAV_Y, FULL_RATES);
            for (int t = 0; t < MAX_TIERS; t++) {
                drive_anchor(&sc[t], f);
                float dev = step_deviation(&sc[t], &sl[t], DT, 0, GRAV_Y,
                                           FULL_RATES, FULL_ITERS);
                lod_step(&sc[t], &sl[t], DT, 0, GRAV_Y, FULL_RATES);
                if (f > FRAMES / 2) {
                    float e = pose_error_max(&sc[t], &sref);
                    err_mean[t] += e;
                    if (e > err_max[t]) err_max[t] = e;
                    if (dev > dev_max[t]) dev_max[t] = dev;
                }
            }
        }
        for (int t = 0; t < MAX_TIERS; t++) {
            err_mean[t] /= (float)(FRAMES / 2);
            char name[64];
            if (bands[t].step_bound < 0.0f) {
                printf("sway_%s err_mean=%.6f err_max=%.6f "
                       "step_dev_max=%.6f step_bound=unbounded\n",
                       band_name(t), err_mean[t], err_max[t], dev_max[t]);
                snprintf(name, sizeof name, "sway_%s_zero_cons",
                         band_name(t));
                check(name, bands[t].stride == 0 && sl[t].n_cons == 0);
            } else {
                int ok = dev_max[t] <= bands[t].step_bound;
                printf("sway_%s err_mean=%.6f err_max=%.6f "
                       "step_dev_max=%.6f step_bound=%.6f within=%d interp=%.6f\n",
                       band_name(t), err_mean[t], err_max[t], dev_max[t],
                       bands[t].step_bound, ok, interp_error(&sref, bands[t].stride));
                snprintf(name, sizeof name, "sway_%s_step_within_bound", band_name(t));
                check(name, ok);
            }
        }
        /* step deviation is monotone in the reduction level */
        check("step_dev_monotone", dev_max[0] <= dev_max[1] && dev_max[1] <= dev_max[2]);
        /* the swinging reference actually hangs: the scenario is live */
        check("sway_reference_hangs", sref.p[N_PARTICLES - 1].y < -5.0f);
    }

/* ---- selection by distance and screen coverage ---- */
    {
        int cur = -1;
        cur = band_select(bands, 4, 3.0f, 0.80f, cur, HYS);
        check("select_near", cur == 0);
        cur = band_select(bands, 4, 12.0f, 0.20f, cur, HYS);
        check("select_mid", cur == 1);
        cur = band_select(bands, 4, 40.0f, 0.05f, cur, HYS);
        check("select_far", cur == 2);
        cur = band_select(bands, 4, 500.0f, 0.00f, cur, HYS);
        check("select_cull", cur == 3);
    }
    check("select_coverage_can_demote",
          band_select(bands, 4, 3.0f, 0.01f, -1, HYS) != 0);

    /* hysteresis: a camera parked on the near/mid boundary must not flip.
       At the boundary distance the fine band's condition still holds; just
       past it, the coarse band's condition also holds. */
    {
        float edge = bands[0].max_dist;
        float just_past = edge * 1.01f;
        float just_before = edge * 0.99f;
        int a = band_select(bands, 4, just_past, 0.5f, 0, HYS);
        check("hysteresis_holds_mid", a == 0);
        int b = band_select(bands, 4, just_past, 0.5f, 0, 0.0f);
        check("no_hysteresis_flips", b == 1);
        int c = band_select(bands, 4, just_before, 0.5f, 1, HYS);
        check("hysteresis_holds_near", c == 1);
        check("hysteresis_clear_promotes",
              band_select(bands, 4, 2.0f, 0.5f, 1, HYS) == 0);
        check("hysteresis_far_demotes",
              band_select(bands, 4, 30.0f, 0.5f, 0, HYS) == 2);
    }

    check("coverage_monotonic",
          screen_coverage(9.0f, 5.0f, HALF_FOV) > screen_coverage(9.0f, 50.0f, HALF_FOV));
    check("coverage_clamped", screen_coverage(1e6f, 1.0f, 1.4f) <= 1.0f);
    {
        /* band error bounds expressed in pixels at the far edge of each band */
        float px_near = world_per_pixel(bands[0].max_dist, HALF_FOV, VIEW_PX);
        float px_mid = world_per_pixel(bands[1].max_dist, HALF_FOV, VIEW_PX);
        float px_far = world_per_pixel(bands[2].max_dist, HALF_FOV, VIEW_PX);
        printf("bound_px near=%.3f mid=%.3f far=%.3f\n",
               bands[0].settle_bound / px_near,
               bands[1].settle_bound / px_mid,
               bands[2].settle_bound / px_far);
        check("bounds_subpixel_far",
              bands[2].settle_bound / px_far < 1.0f);
    }

    /* ---- switch continuity: fine <-> coarse <-> cull ---- */
    {
        Chain c;
        init_hang(&c);
        Lod lod;
        memset(&lod, 0, sizeof lod);
        lod_switch(&c, &lod, 1, FULL_ITERS);
        for (int f = 0; f < 120; f++)
            lod_step(&c, &lod, DT, 0, GRAV_Y, FULL_RATES);

        /* the switch call itself must not move anything */
        Chain snap = c;
        lod_switch(&c, &lod, 2, 8);
        check("switch_positions_untouched", pos_delta(&c, &snap) == 0.0f);
        check("switch_velocity_untouched",
              memcmp(c.p, snap.p, sizeof c.p) == 0);
        check("switch_engaging", lod.engaging == 1);

        /* First frame after the switch: the constraint correction is capped
           at ENGAGE_RATE*dt, so the particle cannot move further than that
           because of the switch on any single frame. */
        float cap = ENGAGE_RATE * DT;
        lod_step(&c, &lod, DT, 0, GRAV_Y, FULL_RATES);
        printf("switch_fine_to_mid_frame1_correction=%.6f cap=%.6f\n",
               lod.last_correction, cap);
        check("switch_fine_to_mid_continuous", lod.last_correction <= cap + 1e-6f);

        /* through the whole engagement: no frame may exceed the cap */
        float worst = 0.0f;
        for (int f = 0; f < 120; f++) {
            lod_step(&c, &lod, DT, 0, GRAV_Y, FULL_RATES);
            if (lod.last_correction > worst) worst = lod.last_correction;
        }
        printf("switch_fine_to_mid_engage_max_correction=%.6f\n", worst);
        check("switch_fine_to_mid_no_later_pop", worst <= cap + 1e-6f);
        check("switch_engagement_finishes", lod.engaging == 0);

        /* after engagement the mid tier must sit close to the full pose:
           the cap may only delay convergence, not prevent it */
        Chain full_ref;
        init_hang(&full_ref);
        Lod full_ref_lod;
        lod_full(&full_ref_lod, &full_ref, FULL_ITERS);
        for (int f = 0; f < 120; f++)
            lod_step(&full_ref, &full_ref_lod, DT, 0, GRAV_Y, FULL_RATES);
        float settled = pose_error_max(&c, &full_ref);
        printf("switch_fine_to_mid_settled_err=%.6f\n", settled);
        check("switch_fine_to_mid_converges", settled < 0.05f);

        /* back to fine: same cap guarantee */
        Chain pre = c;
        lod_switch(&c, &lod, 1, FULL_ITERS);
        check("switch_back_positions_untouched", pos_delta(&c, &pre) == 0.0f);
        lod_step(&c, &lod, DT, 0, GRAV_Y, FULL_RATES);
        printf("switch_mid_to_fine_frame1_correction=%.6f\n", lod.last_correction);
        check("switch_mid_to_fine_continuous", lod.last_correction <= cap + 1e-6f);

        /* into cull: exact at the instant, free fall after. No constraints,
           so re-engaging later reports a zero correction on the switch
           frame itself. */
        snap = c;
        lod_switch(&c, &lod, 0, 0);
        check("switch_to_cull_positions_untouched", pos_delta(&c, &snap) == 0.0f);
        lod_step(&c, &lod, DT, 0, GRAV_Y, FULL_RATES);
        check("switch_to_cull_freefall_step",
              c.p[N_PARTICLES - 1].y < snap.p[N_PARTICLES - 1].y);

        /* stay culled long enough that every edge has stretched far beyond
           its rest length, then return to full detail */
        for (int f = 0; f < 60; f++)
            lod_step(&c, &lod, DT, 0, GRAV_Y, FULL_RATES);
        float cull_seg = chain_seg_error(&c);
        snap = c;
        lod_switch(&c, &lod, 1, FULL_ITERS);
        check("switch_from_cull_positions_untouched", pos_delta(&c, &snap) == 0.0f);
        float cull_worst = 0.0f;
        int cull_frames = 0;
        for (int f = 0; f < 240; f++) {
            lod_step(&c, &lod, DT, 0, GRAV_Y, FULL_RATES);
            if (lod.last_correction > cull_worst) cull_worst = lod.last_correction;
            cull_frames++;
            if (!lod.engaging) break;
        }
        printf("cull_seg_error_before_return=%.6f\n", cull_seg);
        printf("switch_cull_to_fine_max_correction=%.6f cap=%.6f frames=%d\n",
               cull_worst, cap, cull_frames);
        check("switch_from_cull_continuous", cull_worst <= cap + 1e-6f);
        /* and it must actually re-close the stretched edges. Convergence
           after a big stretch is not instant: it takes as long as the cap
           requires, so allow it to settle before asserting. */
        for (int f = 0; f < 240; f++)
            lod_step(&c, &lod, DT, 0, GRAV_Y, FULL_RATES);
        printf("switch_cull_to_fine_seg_error_after=%.6f\n", chain_seg_error(&c));
        check("switch_from_cull_recloses", chain_seg_error(&c) < 0.05f);
    }

    /* goals relax toward targets at a bounded rate, no overshoot */
    {
        Chain c;
        init_hang(&c);
        Lod lod;
        memset(&lod, 0, sizeof lod);
        lod_switch(&c, &lod, 1, FULL_ITERS);
        for (int f = 0; f < 30; f++)
            lod_step(&c, &lod, DT, 0, GRAV_Y, FULL_RATES);
        lod_switch(&c, &lod, 3, 4);
        float g0 = lod.cons[1].goal, target = lod.cons[1].target;
        int overshoot = 0, out_of_rate = 0;
        for (int f = 0; f < 120; f++) {
            float prev = lod.cons[1].goal;
            lod_step(&c, &lod, DT, 0, GRAV_Y, FULL_RATES);
            if (fabsf(lod.cons[1].goal - prev) > RELAX_RATE * DT + 1e-5f)
                out_of_rate = 1;
            if (fabsf(lod.cons[1].goal - target) > fabsf(g0 - target) + 1e-4f)
                overshoot = 1;
        }
        printf("relax_goal_start=%.6f relax_goal_end=%.6f target=%.6f\n",
               g0, lod.cons[1].goal, target);
        check("relax_rate_limited", !out_of_rate);
        check("relax_no_overshoot", !overshoot);
        check("relax_reaches_target", fabsf(lod.cons[1].goal - target) < 1e-4f);
        check("relax_kept_tethers_hold", lod_seg_error(&c, &lod) < 0.05f);
    }

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
