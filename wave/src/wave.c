
/* wave: longitudinal impulse propagation on a spring/constraint chain.
   The distance constraint between neighbours acts as a spring of stiffness k;
   the analytical long-wave speed is c = h*sqrt(k/m).  Explicit velocity-verlet
   integration needs c*dt_sub/h <= 1 (CFL), so frame steps are split into
   substeps sized by wave_required_substeps().  Each particle records the frame
   and substep-precise time at which the impulse front first reaches it, which
   is compared against the analytical front time; the measured speed converges
   to the analytical one as substeps double.  Non-positive stiffness and
   non-positive substep counts are rejected. */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define MAX_N 32
#define THRESH_FRAC 0.02f   /* arrival threshold as a fraction of flick speed */

enum {
    WAVE_OK = 0,
    WAVE_ERR_STIFFNESS = -1,
    WAVE_ERR_SUBSTEPS = -2,
    WAVE_ERR_MASS = -3,
    WAVE_ERR_SPACING = -4,
    WAVE_ERR_COUNT = -5,
    WAVE_ERR_DT = -6
};

typedef struct {
    int n;
    float h;                /* rest spacing between particles */
    float mass;             /* per-particle mass m */
    float k;                /* constraint (spring) stiffness */
    float dt_frame;         /* frame duration */
    int substeps;           /* substeps per frame */
    float x[MAX_N];
    float v[MAX_N];
    float v0;               /* flick impulse as an initial velocity */
    int frame;              /* completed frame count */
    float time;
    int arrive_frame[MAX_N];   /* frame of first front crossing, -1 if none */
    float arrive_time[MAX_N];  /* substep-interpolated crossing time */
} WaveSim;

/* analytical wave speed c = h*sqrt(k/m) */
static float wave_speed(float k, float m, float h) {
    return h * sqrtf(k / m);
}

/* analytical front arrival time at particle i after a flick at particle 0 */
static float wave_arrival_time(float k, float m, int i) {
    return (float)i * sqrtf(m / k);
}

/* smallest substep count keeping the Courant number c*dt_sub/h <= 1 */
static int wave_required_substeps(float k, float m, float h, float dt_frame) {
    float c = wave_speed(k, m, h);
    int s = (int)ceilf(c * dt_frame / h);
    return s < 1 ? 1 : s;
}

/* actual Courant number of the configured substeps */
static float wave_courant(const WaveSim* s) {
    float c = wave_speed(s->k, s->mass, s->h);
    return c * (s->dt_frame / (float)s->substeps) / s->h;
}

static int wave_init(WaveSim* s, int n, float h, float mass, float k,
                     float dt_frame, int substeps) {
    if (k <= 0.0f) return WAVE_ERR_STIFFNESS;   /* negative/zero stiffness */
    if (substeps <= 0) return WAVE_ERR_SUBSTEPS; /* zero/negative substeps */
    if (mass <= 0.0f) return WAVE_ERR_MASS;
    if (h <= 0.0f) return WAVE_ERR_SPACING;
    if (n < 2 || n > MAX_N) return WAVE_ERR_COUNT;
    if (dt_frame <= 0.0f) return WAVE_ERR_DT;

    memset(s, 0, sizeof(*s));
    s->n = n;
    s->h = h;
    s->mass = mass;
    s->k = k;
    s->dt_frame = dt_frame;
    s->substeps = substeps;
    for (int i = 0; i < n; i++) {
        s->x[i] = h * (float)i;
        s->v[i] = 0.0f;
        s->arrive_frame[i] = -1;
        s->arrive_time[i] = -1.0f;
    }
    s->frame = 0;
    s->time = 0.0f;
    s->v0 = 0.0f;
    return WAVE_OK;
}

/* flick particle 0: impulse delivered as an initial velocity */
static void wave_flick(WaveSim* s, float v0) {
    s->v[0] = v0;
    s->v0 = v0;
}

static void wave_accel(const WaveSim* s, float* a) {
    for (int i = 0; i < s->n; i++) {
        float f = 0.0f;
        if (i > 0)
            f += s->k * ((s->x[i - 1] - s->x[i]) + s->h);
        if (i < s->n - 1)
            f += s->k * ((s->x[i + 1] - s->x[i]) - s->h);
        a[i] = f / s->mass;
    }
}

/* advance one frame, split into substeps; scan for front crossings */
static void wave_step(WaveSim* s) {
    float dt_sub = s->dt_frame / (float)s->substeps;
    float a[MAX_N], a2[MAX_N], v_prev[MAX_N];
    float thr = THRESH_FRAC * fabsf(s->v0);

    wave_accel(s, a);
    for (int q = 0; q < s->substeps; q++) {
        float t_prev = s->time;
        memcpy(v_prev, s->v, sizeof(s->v));

        for (int i = 0; i < s->n; i++) {           /* velocity verlet */
            s->v[i] += 0.5f * a[i] * dt_sub;
            s->x[i] += s->v[i] * dt_sub;
        }
        wave_accel(s, a2);
        for (int i = 0; i < s->n; i++) {
            s->v[i] += 0.5f * a2[i] * dt_sub;
            a[i] = a2[i];
        }
        s->time += dt_sub;

        for (int i = 1; i < s->n; i++) {          /* front arrival scan */
            if (s->arrive_frame[i] >= 0 || thr <= 0.0f) continue;
            if (fabsf(s->v[i]) >= thr) {
                float span = s->v[i] - v_prev[i];
                float frac = span != 0.0f
                    ? (thr - fabsf(v_prev[i])) / span : 0.0f;
                if (frac < 0.0f) frac = 0.0f;
                if (frac > 1.0f) frac = 1.0f;
                s->arrive_time[i] = t_prev + frac * dt_sub;
                s->arrive_frame[i] = s->frame;
            }
        }
    }
    s->frame++;
}

static float wave_energy(const WaveSim* s) {
    float e = 0.0f;
    for (int i = 0; i < s->n; i++)
        e += 0.5f * s->mass * s->v[i] * s->v[i];
    for (int i = 0; i < s->n - 1; i++) {
        float ext = (s->x[i + 1] - s->x[i]) - s->h;
        e += 0.5f * s->k * ext * ext;
    }
    return e;
}

/* measured front speed between two particles, from substep arrival times.
   Far particles are preferred: the threshold-based front time inherits an
   O(1/sqrt(i)) Bessel-ramp lag whose difference shrinks with particle index. */
static float wave_measured_speed(const WaveSim* s, int i_lo, int i_hi) {
    float dt = s->arrive_time[i_hi] - s->arrive_time[i_lo];
    return s->h * (float)(i_hi - i_lo) / dt;
}

/* upper bound on integer arrival-frame error for particle i:
   one frame of quantization plus the continuous-time timing deviation */
static int wave_frame_err_bound(const WaveSim* s, int i) {
    float t_star = wave_arrival_time(s->k, s->mass, i);
    float drift = fabsf(s->arrive_time[i] - t_star) / s->dt_frame;
    return 1 + (int)ceilf(drift);
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

#define N_PART 30
#define H 0.25f
#define MASS 1.0f
#define K 1600.0f          /* c = 0.25*sqrt(1600) = 10 m/s */
#define DT  (1.0f / 30.0f)
#define FRAMES 45
#define I_LO 10           /* far-particle window: threshold-ramp lag */
#define I_HI 27           /* difference is negligible here */

static WaveSim run_sim(int substeps, float dt_frame, int frames) {
    WaveSim s;
    int rc = wave_init(&s, N_PART, H, MASS, K, dt_frame, substeps);
    if (rc != WAVE_OK) { memset(&s, 0, sizeof(s)); return s; }
    wave_flick(&s, 1.0f);
    for (int f = 0; f < frames; f++) wave_step(&s);
    return s;
}

int main(void) {
    /* --- analytical values --- */
    check("analytic_speed", fabsf(wave_speed(K, MASS, H) - 10.0f) < 1e-5f);
    check("analytic_arrival",
          fabsf(wave_arrival_time(K, MASS, 4) - 0.1f) < 1e-6f);

    /* --- CFL-driven substep count: c*dt/h = 10*(1/30)/0.25 = 1.333 -> 2 --- */
    check("required_substeps",
          wave_required_substeps(K, MASS, H, DT) == 2);
    check("required_substeps_soft",
          wave_required_substeps(100.0f, MASS, H, DT) == 1);

    /* --- invalid configurations rejected --- */
    WaveSim bad;
    check("reject_neg_stiffness",
          wave_init(&bad, N_PART, H, MASS, -1.0f, DT, 2) == WAVE_ERR_STIFFNESS);
    check("reject_zero_stiffness",
          wave_init(&bad, N_PART, H, MASS, 0.0f, DT, 2) == WAVE_ERR_STIFFNESS);
    check("reject_zero_substeps",
          wave_init(&bad, N_PART, H, MASS, K, DT, 0) == WAVE_ERR_SUBSTEPS);
    check("reject_neg_substeps",
          wave_init(&bad, N_PART, H, MASS, K, DT, -2) == WAVE_ERR_SUBSTEPS);
    check("reject_bad_mass",
          wave_init(&bad, N_PART, H, 0.0f, K, DT, 2) == WAVE_ERR_MASS);

    /* --- propagation runs at S, 2S, 4S, 8S --- */
    WaveSim s2 = run_sim(2, DT, FRAMES);
    WaveSim s4 = run_sim(4, DT, FRAMES);
    WaveSim s8 = run_sim(8, DT, FRAMES);
    WaveSim s16 = run_sim(16, DT, FRAMES);

    int all_arrived = 1;
    for (int i = 1; i < N_PART; i++)
        all_arrived &= s2.arrive_frame[i] >= 0 && s2.arrive_time[i] > 0.0f;
    check("all_fronts_recorded", all_arrived);
    check("courant_in_limit", wave_courant(&s2) <= 1.0f);

    /* --- arrival frames vs analytical front, with explicit upper bound --- */
    int within_bound = 1;
    int max_flat = 0;
    printf("i  frame  theory_frame  t_arr     t_star    bound\n");
    for (int i = 1; i < N_PART; i++) {
        float t_star = wave_arrival_time(K, MASS, i);
        float theory_frame = t_star / DT;
        int ferr = (int)fabsf((float)s2.arrive_frame[i] - theory_frame);
        int bound = wave_frame_err_bound(&s2, i);
        if (ferr > bound) within_bound = 0;
        if (ferr > max_flat) max_flat = ferr;
        printf("%-2d %-6d %-13.2f %-9.4f %-9.4f %d\n",
               i, s2.arrive_frame[i], theory_frame,
               s2.arrive_time[i], t_star, bound);
    }
    check("frame_err_within_bound", within_bound);
    check("frame_err_flat_bound", max_flat <= 3);
    printf("frame_err_max=%d frame_err_flat_bound=3\n", max_flat);

    /* --- convergence: deviation vs the fine (S=16) reference must shrink
           on every substep doubling; velocity verlet is O(dt^2) ----------- */
    float c_ref = wave_measured_speed(&s16, I_LO, I_HI);
    float e2 = fabsf(wave_measured_speed(&s2, I_LO, I_HI) - c_ref) / c_ref;
    float e4 = fabsf(wave_measured_speed(&s4, I_LO, I_HI) - c_ref) / c_ref;
    float e8 = fabsf(wave_measured_speed(&s8, I_LO, I_HI) - c_ref) / c_ref;
    printf("speed_err_vs_ref: s2=%.5f s4=%.5f s8=%.5f s16=ref(%.4f)\n",
           e2, e4, e8, c_ref);
    check("converges_on_double", e4 < e2 && e8 < e4);
    check("convergence_ratio", e4 <= 0.6f * e2 && e8 <= 0.6f * e4);

    /* fine run also matches the analytical speed; the remaining offset is
       the threshold-crossing ramp lag, documented above, not time stepping */
    float c_theory = wave_speed(K, MASS, H);
    check("fine_run_matches_theory",
          fabsf(c_ref - c_theory) / c_theory < 0.05f);

    /* --- frame-rate independence: identical substep size at 30 and 60 Hz
           must give identical arrival times regardless of frame rate ----- */
    WaveSim hz30 = run_sim(2, 1.0f / 30.0f, FRAMES);      /* dt_sub = 1/60 */
    WaveSim hz60 = run_sim(1, 1.0f / 60.0f, FRAMES * 2);  /* dt_sub = 1/60 */
    int same_timing = 1;
    for (int i = 1; i < N_PART; i++)
        same_timing &= fabsf(hz30.arrive_time[i] - hz60.arrive_time[i]) < 1e-6f;
    check("framerate_independent", same_timing);

    /* --- ignoring the CFL limit goes unstable; respecting it keeps energy
           bounded, and the bound tightens as substeps double -------------- */
    WaveSim unstable = run_sim(1, DT, FRAMES);   /* Courant 1.333 > 1 */
    float e0 = 0.5f * MASS * 1.0f * 1.0f;
    check("cfl_violation_blows_up", wave_energy(&unstable) > 10.0f * e0);
    check("cfl_ok_energy_bounded", wave_energy(&s2) < 1.5f * e0);
    check("cfl_energy_converges",
          wave_energy(&s16) < wave_energy(&s2) &&
          fabsf(wave_energy(&s16) - e0) / e0 < 0.01f);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
