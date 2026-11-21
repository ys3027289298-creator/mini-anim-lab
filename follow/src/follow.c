
/* follow: spring-damper camera follow with target-velocity feedforward.
   Mass-normalized model (mass = 1), g(t) is the target:

       e  = x - g                       (signed lag)
       x'' = -wn^2 * e - 2*zeta*wn * (x' - g')

   The damper acts on *relative* velocity, so a target moving at constant
   speed is tracked with zero steady-state lag (feedforward).
   Critical damping: c = 2*sqrt(k) = 2*wn with k = wn^2.
   Overshoot fraction Os maps to zeta = sqrt(L^2 / (pi^2 + L^2)), L = ln Os.
   Configs that would yield zeta <= 0 (wn <= 0, Os >= 1, ...) are rejected. */
#include <stdio.h>
#include <math.h>
#include <string.h>

static float pi_f(void) { return acosf(-1.0f); }

static float follow1d_analytic_overshoot(float zeta);

/* ---------------- 1-D follower ---------------- */

typedef struct {
    float x, v;       /* follower position / velocity */
    float wn;         /* natural angular frequency, rad/s */
    float zeta;       /* damping ratio; 1 = critical, >1 overdamped */
    float k;          /* stiffness, always wn^2 */
    float c;          /* damping coefficient, always 2*zeta*wn */
    float overshoot;  /* configured overshoot fraction, 0..1 (0 => critical) */
} Follow1D;

static int follow1d_init_zeta(Follow1D* f, float wn, float zeta) {
    /* zeta must be strictly positive: zero is undamped (never settles),
       negative pumps energy into the system instead of removing it. */
    if (isnan(wn) || isinf(wn) || wn <= 0.0f) return -1;
    if (isnan(zeta) || isinf(zeta) || zeta <= 0.0f) return -1;
    memset(f, 0, sizeof(*f));
    f->wn = wn;
    f->zeta = zeta;
    f->k = wn * wn;
    f->c = 2.0f * zeta * wn;
    f->overshoot = follow1d_analytic_overshoot(zeta);
    return 0;
}

static int follow1d_init_critical(Follow1D* f, float wn) {
    if (isnan(wn) || isinf(wn) || wn <= 0.0f) return -1;
    memset(f, 0, sizeof(*f));
    f->wn = wn;
    f->zeta = 1.0f;
    f->k = wn * wn;
    f->c = 2.0f * wn;          /* critical coefficient: c^2 = 4k */
    f->overshoot = 0.0f;
    return 0;
}

static int follow1d_init_overshoot(Follow1D* f, float wn, float overshoot) {
    if (isnan(wn) || isinf(wn) || wn <= 0.0f) return -1;
    /* Os in (0,1) maps to 0 < zeta < 1; Os = 0 means critical damping.
       Os >= 1 would require zeta = 0 (pure oscillation), so reject it. */
    if (isnan(overshoot) || isinf(overshoot) ||
        overshoot < 0.0f || overshoot >= 1.0f) return -1;
    float zeta;
    if (overshoot <= 0.0f) {
        zeta = 1.0f;
    } else {
        float L = logf(overshoot);
        zeta = sqrtf((L * L) / (pi_f() * pi_f() + L * L));
    }
    return follow1d_init_zeta(f, wn, zeta) == 0
               ? (f->overshoot = overshoot, 0)
               : -1;
}

static float follow1d_step(Follow1D* f, float dt,
                           float target, float target_vel) {
    if (dt <= 0.0f) return f->x;
    /* substep so that wn*h <= 0.05: semi-implicit euler then tracks the
       analytic flow closely even with long frame times (e.g. 1/30 s). */
    int n = (int)ceilf(f->wn * dt / 0.05f);
    if (n < 1) n = 1;
    float h = dt / (float)n;
    for (int i = 0; i < n; i++) {
        float lag = f->x - target;                 /* target held over frame */
        float rel_v = f->v - target_vel;           /* feedforward lives here */
        float a = -f->k * lag - f->c * rel_v;
        f->v += a * h;                             /* semi-implicit euler: */
        f->x += f->v * h;                          /* new v advances x */
    }
    return f->x;
}

/* ---------------- analytic reference (unit step, zero IC) ---------------- */

static float follow1d_analytic_step(float zeta, float wn, float t) {
    if (t <= 0.0f) return 0.0f;
    if (zeta < 1.0f) {                             /* underdamped */
        float root = sqrtf(1.0f - zeta * zeta);
        float wd = wn * root;
        float phi = atanf(root / zeta);
        return 1.0f - expf(-zeta * wn * t) / root * sinf(wd * t + phi);
    }
    if (zeta == 1.0f) {                            /* critical */
        float s = wn * t;
        return 1.0f - (1.0f + s) * expf(-s);
    }
    float r = sqrtf(zeta * zeta - 1.0f);           /* overdamped */
    float a = zeta - r, b = zeta + r;
    return 1.0f + (b * expf(-a * wn * t) - a * expf(-b * wn * t)) / (a - b);
}

/* first peak height of the unit step response; zero for zeta >= 1 */
static float follow1d_analytic_overshoot(float zeta) {
    if (!(zeta > 0.0f) || zeta >= 1.0f) return 0.0f;
    return expf(-pi_f() * zeta / sqrtf(1.0f - zeta * zeta));
}

/* upper bound on the 2% settling time: after this instant the analytic
   step response is guaranteed to stay within +-2% of the target. */
static float follow1d_settling_bound_2pct(float zeta, float wn) {
    if (zeta < 1.0f) {
        /* |e| <= exp(-zeta*wn*t)/sqrt(1-zeta^2) */
        float root = sqrtf(1.0f - zeta * zeta);
        return -logf(0.02f * root) / (zeta * wn);
    }
    if (zeta == 1.0f) {
        /* |e| = (1+s)e^-s, s = wn*t; solve envelope == 0.02 */
        float lo = 0.0f, hi = 1.0f;
        while ((1.0f + hi) * expf(-hi) > 0.02f) hi *= 2.0f;
        for (int i = 0; i < 60; i++) {
            float mid = 0.5f * (lo + hi);
            if ((1.0f + mid) * expf(-mid) > 0.02f) lo = mid;
            else hi = mid;
        }
        return 0.5f * (lo + hi) / wn;
    }
    /* |e| <= b/(b-a) * exp(-a*wn*t), a = zeta-sqrt(zeta^2-1) (slow pole) */
    float r = sqrtf(zeta * zeta - 1.0f);
    float a = zeta - r, b = zeta + r;
    return -logf(0.02f * (b - a) / b) / (a * wn);
}

/* ---------------- 3-D wrapper: one independent follower per axis ---------------- */

typedef struct { Follow1D axis[3]; } Follow3D;

static int follow3d_init_overshoot(Follow3D* f, float wn, float overshoot) {
    int rc = 0;
    for (int i = 0; i < 3; i++)
        if (follow1d_init_overshoot(&f->axis[i], wn, overshoot) != 0) rc = -1;
    return rc;
}

static int follow3d_init_critical(Follow3D* f, float wn) {
    int rc = 0;
    for (int i = 0; i < 3; i++)
        if (follow1d_init_critical(&f->axis[i], wn) != 0) rc = -1;
    return rc;
}

static void follow3d_step(Follow3D* f, float dt,
                          const float target[3], const float target_vel[3],
                          float out[3]) {
    for (int i = 0; i < 3; i++)
        out[i] = follow1d_step(&f->axis[i], dt, target[i], target_vel[i]);
}

/* ---------------- self tests ---------------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

typedef struct {
    int mode;             /* 0 critical, 1 overshoot, 2 raw zeta */
    float wn, os, zeta;
} StepCfg;

typedef struct {
    float max_overshoot;  /* max(x - 1) over the run */
    float max_abs_err;    /* numeric vs analytic step response */
    float settling_t;     /* last t with |x - 1| > 0.02 (0 if never) */
    float final_x;
    float zeta;
} StepResult;

/* the bound is derived from the analytic response envelope, so verify it
   against the analytic flow itself over a full second past the bound. */
static int analytic_bound_holds(float zeta, float wn) {
    float bd = follow1d_settling_bound_2pct(zeta, wn);
    for (float t = bd; t <= bd + 1.0f + 1e-9f; t += 0.01f) {
        float e = fabsf(follow1d_analytic_step(zeta, wn, t) - 1.0f);
        if (e > 0.02f + 1e-5f) return 0;
    }
    return 1;
}

static StepResult simulate_step(const StepCfg* cfg, float dt, float T) {
    Follow1D f;
    if (cfg->mode == 0) follow1d_init_critical(&f, cfg->wn);
    else if (cfg->mode == 1) follow1d_init_overshoot(&f, cfg->wn, cfg->os);
    else follow1d_init_zeta(&f, cfg->wn, cfg->zeta);

    StepResult r;
    memset(&r, 0, sizeof(r));
    r.zeta = f.zeta;
    int steps = (int)(T / dt);
    for (int i = 1; i <= steps; i++) {
        float t = i * dt;
        float x = follow1d_step(&f, dt, 1.0f, 0.0f);
        float an = follow1d_analytic_step(f.zeta, f.wn, t);
        float err = fabsf(x - an);
        if (err > r.max_abs_err) r.max_abs_err = err;
        if (x - 1.0f > r.max_overshoot) r.max_overshoot = x - 1.0f;
        if (fabsf(x - 1.0f) > 0.02f) r.settling_t = t;
    }
    r.final_x = f.x;
    return r;
}

int main(void) {
    /* --- parameter relations --- */
    Follow1D a;
    follow1d_init_critical(&a, 10.0f);
    check("critical_relation",
          a.zeta == 1.0f && fabsf(a.c * a.c - 4.0f * a.k) < 1e-4f &&
          fabsf(a.k - a.wn * a.wn) < 1e-4f);

    Follow1D b;
    check("init_overshoot_ok", follow1d_init_overshoot(&b, 10.0f, 0.05f) == 0);
    check("os_roundtrip",
          fabsf(follow1d_analytic_overshoot(b.zeta) - 0.05f) < 1e-4f);
    check("zeta_derived_in_range", b.zeta > 0.0f && b.zeta < 1.0f);
    check("coeff_relation",
          fabsf(b.c - 2.0f * b.zeta * b.wn) < 1e-5f &&
          fabsf(b.k - b.wn * b.wn) < 1e-4f);

    Follow1D lo;
    follow1d_init_overshoot(&lo, 10.0f, 0.30f);
    check("more_overshoot_less_zeta", lo.zeta < b.zeta && lo.zeta > 0.0f);

    Follow1D d0;
    check("os_zero_is_critical",
          follow1d_init_overshoot(&d0, 10.0f, 0.0f) == 0 && d0.zeta == 1.0f);

    /* --- reject configs that give zeta <= 0 or are otherwise invalid --- */
    Follow1D tmp;
    check("rej_wn_zero", follow1d_init_critical(&tmp, 0.0f) == -1);
    check("rej_wn_neg", follow1d_init_critical(&tmp, -5.0f) == -1);
    check("rej_wn_nan", follow1d_init_critical(&tmp, NAN) == -1);
    check("rej_os_one", follow1d_init_overshoot(&tmp, 10.0f, 1.0f) == -1);
    check("rej_os_big", follow1d_init_overshoot(&tmp, 10.0f, 1.5f) == -1);
    check("rej_os_neg", follow1d_init_overshoot(&tmp, 10.0f, -0.1f) == -1);
    check("rej_zeta_zero", follow1d_init_zeta(&tmp, 10.0f, 0.0f) == -1);
    check("rej_zeta_neg", follow1d_init_zeta(&tmp, 10.0f, -0.5f) == -1);
    check("rej_zeta_nan", follow1d_init_zeta(&tmp, 10.0f, NAN) == -1);
    check("rej_zeta_inf", follow1d_init_zeta(&tmp, 10.0f, INFINITY) == -1);
    check("overdamped_accepted",
          follow1d_init_zeta(&tmp, 10.0f, 2.0f) == 0 && tmp.zeta == 2.0f);

    Follow1D keep;
    follow1d_init_critical(&keep, 10.0f);
    int rc = follow1d_init_overshoot(&keep, -1.0f, 0.05f);
    check("failed_init_preserves", rc == -1 && keep.wn == 10.0f &&
                                   keep.zeta == 1.0f);

    Follow1D frozen;
    follow1d_init_critical(&frozen, 10.0f);
    float x0 = frozen.x;
    check("nonpositive_dt_noop",
          follow1d_step(&frozen, 0.0f, 1.0f, 0.0f) == x0 &&
          follow1d_step(&frozen, -1.0f, 1.0f, 0.0f) == x0);

    /* --- numeric step response vs analytic solution, three regimes --- */
    const StepCfg under_cfg = {1, 10.0f, 0.05f, 0.0f};   /* zeta ~ 0.69 */
    const StepCfg crit_cfg  = {0, 10.0f, 0.0f, 0.0f};
    const StepCfg over_cfg  = {2, 10.0f, 0.0f, 2.0f};

    StepResult ru = simulate_step(&under_cfg, 0.001f, 3.0f);
    StepResult rcrit = simulate_step(&crit_cfg, 0.001f, 3.0f);
    StepResult rov = simulate_step(&over_cfg, 0.001f, 3.0f);
    printf("under: zeta=%.4f max_err=%.6f peak_os=%.4f settle=%.4f\n",
           ru.zeta, ru.max_abs_err, ru.max_overshoot, ru.settling_t);
    printf("crit:  zeta=1.0000 max_err=%.6f peak_os=%.4f settle=%.4f\n",
           rcrit.max_abs_err, rcrit.max_overshoot, rcrit.settling_t);
    printf("over:  zeta=2.0000 max_err=%.6f peak_os=%.4f settle=%.4f\n",
           rov.max_abs_err, rov.max_overshoot, rov.settling_t);

    check("under_matches_analytic", ru.max_abs_err < 0.005f);
    check("crit_matches_analytic", rcrit.max_abs_err < 0.005f);
    check("over_matches_analytic", rov.max_abs_err < 0.005f);

    /* overshoot upper bounds: configured value (analytic) and measured */
    float os_bound = follow1d_analytic_overshoot(ru.zeta);
    check("under_os_within_config",
          ru.max_overshoot <= 0.05f + 0.002f &&
          ru.max_overshoot <= os_bound + 0.002f &&
          ru.max_overshoot > 0.03f);
    check("crit_no_overshoot", rcrit.max_overshoot <= 0.002f);
    check("over_no_overshoot", rov.max_overshoot <= 0.002f);

    /* settling time: analytic flow must obey the envelope upper bound;
       the numeric flow gets a small slack for integrator + grid error. */
    float bd_u = follow1d_settling_bound_2pct(ru.zeta, 10.0f);
    float bd_c = follow1d_settling_bound_2pct(1.0f, 10.0f);
    float bd_o = follow1d_settling_bound_2pct(2.0f, 10.0f);
    printf("settling bounds 2%%: under=%.4f crit=%.4f over=%.4f\n",
           bd_u, bd_c, bd_o);
    check("under_analytic_bound", analytic_bound_holds(ru.zeta, 10.0f));
    check("crit_analytic_bound", analytic_bound_holds(1.0f, 10.0f));
    check("over_analytic_bound", analytic_bound_holds(2.0f, 10.0f));
    check("under_settle_bound", ru.settling_t <= bd_u + 3.0f * 0.001f);
    check("crit_settle_bound", rcrit.settling_t <= bd_c + 3.0f * 0.001f);
    check("over_settle_bound", rov.settling_t <= bd_o + 3.0f * 0.001f);
    check("all_reach_target",
          fabsf(ru.final_x - 1.0f) < 0.001f &&
          fabsf(rcrit.final_x - 1.0f) < 0.001f &&
          fabsf(rov.final_x - 1.0f) < 0.001f);

    /* --- velocity feedforward: constant-velocity ramp target --- */
    Follow1D ramp;
    follow1d_init_critical(&ramp, 8.0f);
    const float v = 3.0f;
    float max_lag = 0.0f;
    for (int i = 0; i < 4000; i++) {
        float t = (i + 1) * 0.001f;
        follow1d_step(&ramp, 0.001f, v * t, v);   /* g' fed in */
        if (t >= 1.0f) {
            float lag = fabsf(ramp.x - v * t);
            if (lag > max_lag) max_lag = lag;
        }
    }
    printf("ramp feedforward max |lag| after t>=1s: %.6f\n", max_lag);
    check("feedforward_zero_lag", max_lag < 0.005f);

    Follow1D noff;
    follow1d_init_critical(&noff, 8.0f);
    for (int i = 0; i < 4000; i++) {
        float t = (i + 1) * 0.001f;
        follow1d_step(&noff, 0.001f, v * t, 0.0f);  /* g' withheld */
    }
    float noff_lag = fabsf(noff.x - v * 4.0f);
    printf("ramp no-feedforward lag at t=4s: %.6f (theory 2*zeta*v/wn=%.4f)\n",
           noff_lag, 2.0f * v / 8.0f);
    check("no_feedforward_shows_lag", noff_lag > 0.5f && noff_lag < 0.9f);

    /* --- long frame times must stay stable thanks to substepping --- */
    Follow1D coarse;
    follow1d_init_overshoot(&coarse, 20.0f, 0.05f);
    float worst = 0.0f;
    for (int i = 0; i < 120; i++) {                 /* 4 s at 30 fps */
        float x = follow1d_step(&coarse, 1.0f / 30.0f, 1.0f, 0.0f);
        if (fabsf(x) > worst) worst = fabsf(x);
    }
    check("coarse_dt_bounded", worst < 1.5f);
    check("coarse_dt_settles", fabsf(coarse.x - 1.0f) < 0.02f);

    /* a single 250 ms hitch must not explode or fling the follower */
    Follow1D hitch;
    follow1d_init_critical(&hitch, 12.0f);
    follow1d_step(&hitch, 0.25f, 1.0f, 0.0f);
    check("hitch_bounded", isfinite(hitch.x) && isfinite(hitch.v) &&
                           fabsf(hitch.x) < 1.5f && fabsf(hitch.v) < 20.0f);
    for (int i = 0; i < 300; i++) follow1d_step(&hitch, 0.01f, 1.0f, 0.0f);
    check("hitch_recovers", fabsf(hitch.x - 1.0f) < 0.02f);

    /* --- 3-D wrapper --- */
    Follow3D f3;
    check("f3_init_ok", follow3d_init_critical(&f3, 10.0f) == 0);
    Follow3D f3bad;
    check("f3_rejects_bad", follow3d_init_overshoot(&f3bad, -1.0f, 0.05f) == -1);
    float tgt[3] = {1.0f, -2.0f, 0.5f};
    float zero3[3] = {0.0f, 0.0f, 0.0f};
    float out[3];
    for (int i = 0; i < 200; i++)
        follow3d_step(&f3, 0.01f, tgt, zero3, out);
    check("f3_tracks",
          fabsf(out[0] - tgt[0]) < 0.02f &&
          fabsf(out[1] - tgt[1]) < 0.02f &&
          fabsf(out[2] - tgt[2]) < 0.02f);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
