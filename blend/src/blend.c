
/* blend: animation-to-ragdoll transition blending.
 *
 * On entering ragdoll the verlet chain would otherwise collapse in one
 * frame. This module keeps the physics chain integrated (distance
 * constraints still solved every step) while the emitted pose is blended
 * from a frozen animation pose toward the live physics pose:
 *
 *   e_i(t) = 1 - alpha(t) * (1 - w[i])
 *   out_i  = e_i * anim_i + (1 - e_i) * physics_i
 *
 * w[i] is the per-bone animation weight at the end of the blend, alpha is
 * a smoothstep ramp over a duration derived from impact strength (stronger
 * impact -> shorter blend). All entry points validate their inputs:
 * out-of-range or non-finite weights and non-positive durations are
 * rejected with an error code.
 */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define MAX_B 16
#define ITERATIONS 8
#define DAMPING 0.99f
#define DT 0.016f
#define PI 3.14159265f
#define TAU 6.2831853f

enum {
    BLEND_OK = 0,
    BLEND_ERR_NULL = -1,
    BLEND_ERR_COUNT = -2,
    BLEND_ERR_WEIGHT = -3,
    BLEND_ERR_DURATION = -4,
    BLEND_ERR_STATE = -5
};

typedef struct { float x, y, ang; } Pose;

typedef struct { float x, y, px, py, inv_mass; int pinned; } Particle;

typedef struct {
    Particle p[MAX_B];
    float rest[MAX_B];   /* rest length between i and i+1 */
    int n;
} PhysChain;

typedef struct {
    int n;
    PhysChain phys;
    Pose anim[MAX_B];
    Pose out[MAX_B];
    float weight[MAX_B];   /* per-bone animation influence at end of blend */
    float duration;
    float elapsed;
    int active;
    /* smoothness bookkeeping on the end effector */
    int have_prev;
    float prev_ee_x, prev_ee_y;
    float prev_vel_x, prev_vel_y;
    float prev_alpha;
    float max_accel;       /* largest frame-to-frame velocity change */
    int reversals;         /* per-axis velocity sign flips past deadzone */
    int alpha_regressed;
} Blend;

/* ---------- validation ---------- */

static int blend_validate_weights(const float* w, int n) {
    if (!w || n <= 0 || n > MAX_B) return BLEND_ERR_WEIGHT;
    for (int i = 0; i < n; i++) {
        if (!isfinite(w[i]) || w[i] < 0.0f || w[i] > 1.0f)
            return BLEND_ERR_WEIGHT;
    }
    return BLEND_OK;
}

/* ---------- impact -> duration ---------- */

#define BLEND_MIN_DUR 0.15f
#define BLEND_MAX_DUR 1.20f

/* stronger impact gives a shorter blend, clamped to [MIN, MAX].
 * negative or non-finite impact is rejected. */
static int blend_duration_from_impact(float impact, float* duration) {
    if (!duration) return BLEND_ERR_NULL;
    if (!isfinite(impact) || impact < 0.0f) return BLEND_ERR_DURATION;
    float d = BLEND_MAX_DUR - impact;
    if (d < BLEND_MIN_DUR) d = BLEND_MIN_DUR;
    if (d > BLEND_MAX_DUR) d = BLEND_MAX_DUR;
    *duration = d;
    return BLEND_OK;
}

/* ---------- physics chain (verlet + distance constraints) ---------- */

static void phys_step(PhysChain* c, float dt, float gx, float gy) {
    for (int i = 0; i < c->n; i++) {
        Particle* p = &c->p[i];
        if (p->pinned) continue;
        float vx = (p->x - p->px) * DAMPING;
        float vy = (p->y - p->py) * DAMPING;
        p->px = p->x; p->py = p->y;
        p->x += vx + gx * dt * dt;
        p->y += vy + gy * dt * dt;
    }
    for (int it = 0; it < ITERATIONS; it++) {
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
    }
}

static float phys_stretch(const PhysChain* c) {
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

static float phys_bone_angle(const PhysChain* c, int i) {
    if (i < c->n - 1)
        return atan2f(c->p[i + 1].y - c->p[i].y,
                      c->p[i + 1].x - c->p[i].x);
    return atan2f(c->p[i].y - c->p[i - 1].y,
                  c->p[i].x - c->p[i - 1].x);
}

/* ---------- angle helpers ---------- */

static float wrap_angle(float a) {
    while (a > PI) a -= TAU;
    while (a < -PI) a += TAU;
    return a;
}

static float lerp_angle(float a, float b, float t) {
    return wrap_angle(a + wrap_angle(b - a) * t);
}

/* ---------- blend controller ---------- */

static int blend_start(Blend* b, const Pose* anim, const Pose* phys,
                       const float* weight, int n, float duration) {
    if (!b || !anim || !phys || !weight) return BLEND_ERR_NULL;
    if (n <= 0 || n > MAX_B) return BLEND_ERR_COUNT;
    int rc = blend_validate_weights(weight, n);
    if (rc != BLEND_OK) return rc;
    if (!isfinite(duration) || duration <= 0.0f) return BLEND_ERR_DURATION;

    memset(b, 0, sizeof(*b));
    b->n = n;
    b->duration = duration;
    for (int i = 0; i < n; i++) {
        b->anim[i] = anim[i];
        b->weight[i] = weight[i];
        b->phys.p[i].x = phys[i].x;
        b->phys.p[i].y = phys[i].y;
        b->phys.p[i].px = phys[i].x;
        b->phys.p[i].py = phys[i].y;
        b->phys.p[i].pinned = (i == 0);
        b->phys.p[i].inv_mass = (i == 0) ? 0.0f : 1.0f;
        b->out[i] = anim[i];
    }
    for (int i = 0; i < n - 1; i++) {
        float dx = phys[i + 1].x - phys[i].x;
        float dy = phys[i + 1].y - phys[i].y;
        b->phys.rest[i] = sqrtf(dx * dx + dy * dy);
    }
    b->phys.n = n;
    b->active = 1;
    return BLEND_OK;
}

/* normalized blend factor: 0 at trigger, 1 at end of transition */
static float blend_alpha(const Blend* b) {
    float u = b->elapsed / b->duration;
    if (u <= 0.0f) return 0.0f;
    if (u >= 1.0f) return 1.0f;
    return u * u * (3.0f - 2.0f * u);
}

/* advances one frame: physics integrated and solved, then output blended.
 * returns BLEND_ERR_STATE once the transition has completed. */
static int blend_step(Blend* b, float dt, float gx, float gy) {
    if (!b) return BLEND_ERR_NULL;
    if (!b->active) return BLEND_ERR_STATE;
    if (!isfinite(dt) || dt <= 0.0f) return BLEND_ERR_DURATION;

    /* constraints keep being solved for the whole transition */
    phys_step(&b->phys, dt, gx, gy);

    b->elapsed += dt;
    float alpha = blend_alpha(b);
    if (alpha < b->prev_alpha - 1e-7f) b->alpha_regressed = 1;
    b->prev_alpha = alpha;

    for (int i = 0; i < b->n; i++) {
        float e = 1.0f - alpha * (1.0f - b->weight[i]);
        b->out[i].x = e * b->anim[i].x + (1.0f - e) * b->phys.p[i].x;
        b->out[i].y = e * b->anim[i].y + (1.0f - e) * b->phys.p[i].y;
        /* e is the animation coefficient, so lerp toward physics by 1-e */
        b->out[i].ang = lerp_angle(b->anim[i].ang,
                                   phys_bone_angle(&b->phys, i), 1.0f - e);
    }

    /* end-effector smoothness metrics */
    Pose* ee = &b->out[b->n - 1];
    if (b->have_prev) {
        float vx = ee->x - b->prev_ee_x;
        float vy = ee->y - b->prev_ee_y;
        float ax = vx - b->prev_vel_x;
        float ay = vy - b->prev_vel_y;
        float acc = sqrtf(ax * ax + ay * ay);
        if (acc > b->max_accel) b->max_accel = acc;
        const float eps = 1e-4f;
        if (vx * b->prev_vel_x < 0.0f &&
            fabsf(vx) > eps && fabsf(b->prev_vel_x) > eps)
            b->reversals++;
        if (vy * b->prev_vel_y < 0.0f &&
            fabsf(vy) > eps && fabsf(b->prev_vel_y) > eps)
            b->reversals++;
        b->prev_vel_x = vx;
        b->prev_vel_y = vy;
    } else {
        b->have_prev = 1;
    }
    b->prev_ee_x = ee->x;
    b->prev_ee_y = ee->y;

    if (b->elapsed >= b->duration) b->active = 0;
    return BLEND_OK;
}

/* ---------- hand-tuned artist reference transition ---------- */

typedef struct { float u; float x, y, ang; } RefKey;

/* manually authored keyframes (normalized time): rounded artist hand-tune
 * of the same sag transition, kept independent of the procedural math.
 * The short pinned chain mostly translates downward, so the end segment
 * only rotates a few hundredths of a radian. */
static const RefKey ref_keys[] = {
    { 0.00f, 7.00f,  0.00f,  0.00f },
    { 0.25f, 7.00f, -0.01f, -0.00f },
    { 0.50f, 6.98f, -0.15f, -0.01f },
    { 0.75f, 6.87f, -0.55f, -0.01f },
    { 1.00f, 6.64f, -1.15f, -0.02f }
};
#define REF_N ((int)(sizeof(ref_keys) / sizeof(ref_keys[0])))

static void ref_sample(float u, Pose* out) {
    if (u <= ref_keys[0].u) { *out = (Pose){ ref_keys[0].x, ref_keys[0].y, ref_keys[0].ang }; return; }
    if (u >= ref_keys[REF_N - 1].u) {
        const RefKey* k = &ref_keys[REF_N - 1];
        *out = (Pose){ k->x, k->y, k->ang };
        return;
    }
    for (int i = 1; i < REF_N; i++) {
        if (u <= ref_keys[i].u) {
            const RefKey* a = &ref_keys[i - 1];
            const RefKey* c = &ref_keys[i];
            float t = (u - a->u) / (c->u - a->u);
            out->x = a->x + (c->x - a->x) * t;
            out->y = a->y + (c->y - a->y) * t;
            out->ang = lerp_angle(a->ang, c->ang, t);
            return;
        }
    }
}

/* ---------- self test ---------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

#define N_BONES 8
#define EE_BOUND_POS 0.06f
#define EE_BOUND_ANG 0.03f
#define JITTER_ACCEL 0.02f

int main(void) {
    /* frozen animation pose: horizontal chain, root pinned at origin */
    Pose anim[N_BONES];
    Pose phys0[N_BONES];
    float weight[N_BONES];
    for (int i = 0; i < N_BONES; i++) {
        anim[i] = (Pose){ (float)i, 0.0f, 0.0f };
        phys0[i] = anim[i];
        weight[i] = 1.0f - 0.7f * (float)i / (N_BONES - 1); /* 1.0 .. 0.3 */
    }

    /* validation: weights must stay inside [0,1] */
    float bad_lo[N_BONES], bad_hi[N_BONES], bad_nan[N_BONES];
    for (int i = 0; i < N_BONES; i++) {
        bad_lo[i] = bad_hi[i] = bad_nan[i] = 0.5f;
    }
    bad_lo[3] = -0.01f;
    bad_hi[3] = 1.01f;
    bad_nan[3] = NAN;
    Blend tmp;
    float dummy_dur = 0.0f;
    check("reject_weight_low",
          blend_start(&tmp, anim, phys0, bad_lo, N_BONES, 0.6f) == BLEND_ERR_WEIGHT);
    check("reject_weight_high",
          blend_start(&tmp, anim, phys0, bad_hi, N_BONES, 0.6f) == BLEND_ERR_WEIGHT);
    check("reject_weight_nan",
          blend_start(&tmp, anim, phys0, bad_nan, N_BONES, 0.6f) == BLEND_ERR_WEIGHT);

    /* validation: duration must be positive and finite */
    check("reject_duration_zero",
          blend_start(&tmp, anim, phys0, weight, N_BONES, 0.0f) == BLEND_ERR_DURATION);
    check("reject_duration_negative",
          blend_start(&tmp, anim, phys0, weight, N_BONES, -0.5f) == BLEND_ERR_DURATION);
    check("reject_duration_nan",
          blend_start(&tmp, anim, phys0, weight, N_BONES, NAN) == BLEND_ERR_DURATION);
    check("reject_null",
          blend_start(NULL, anim, phys0, weight, N_BONES, 0.6f) == BLEND_ERR_NULL);
    check("reject_bad_count",
          blend_start(&tmp, anim, phys0, weight, 0, 0.6f) == BLEND_ERR_COUNT);
    check("reject_bad_impact",
          blend_duration_from_impact(-1.0f, &dummy_dur) == BLEND_ERR_DURATION);

    /* impact mapping: stronger impact -> shorter duration, clamped */
    float d0, d1, d2;
    blend_duration_from_impact(0.0f, &d0);
    blend_duration_from_impact(0.6f, &d1);
    blend_duration_from_impact(100.0f, &d2);
    check("impact_clamp_max", fabsf(d0 - BLEND_MAX_DUR) < 1e-6f);
    check("impact_shorter", d1 < d0);
    check("impact_clamp_min", fabsf(d2 - BLEND_MIN_DUR) < 1e-6f);

    /* run the transition */
    Blend b;
    float duration;
    blend_duration_from_impact(0.6f, &duration);
    int rc = blend_start(&b, anim, phys0, weight, N_BONES, duration);
    check("start_ok", rc == BLEND_OK && b.active);

    /* first emitted frame must still sit on the animation pose: no collapse */
    blend_step(&b, DT, 0.0f, -9.8f);
    float first_dx = b.out[N_BONES - 1].x - anim[N_BONES - 1].x;
    float first_dy = b.out[N_BONES - 1].y - anim[N_BONES - 1].y;
    check("no_instant_collapse",
          sqrtf(first_dx * first_dx + first_dy * first_dy) < 0.05f);

    /* compare every frame against the hand-tuned reference while the
       physics chain keeps solving constraints */
    float max_pos_err = 0.0f, max_ang_err = 0.0f, worst_stretch = 0.0f;
    int frames = 1;
    while (b.active) {
        rc = blend_step(&b, DT, 0.0f, -9.8f);
        if (rc != BLEND_OK) break;
        frames++;
        float u = b.elapsed / b.duration;
        if (u > 1.0f) u = 1.0f;
        Pose ref = { 0.0f, 0.0f, 0.0f };
        ref_sample(u, &ref);
        Pose* ee = &b.out[N_BONES - 1];
        float ex = ee->x - ref.x, ey = ee->y - ref.y;
        float pe = sqrtf(ex * ex + ey * ey);
        float ae = fabsf(wrap_angle(ee->ang - ref.ang));
        if (pe > max_pos_err) max_pos_err = pe;
        if (ae > max_ang_err) max_ang_err = ae;
        float st = phys_stretch(&b.phys);
        if (st > worst_stretch) worst_stretch = st;
    }
    printf("ee_pos_err_max=%.4f (bound %.2f)\n", max_pos_err, EE_BOUND_POS);
    printf("ee_ang_err_max=%.4f (bound %.2f)\n", max_ang_err, EE_BOUND_ANG);
    check("ee_pos_within_bound", max_pos_err <= EE_BOUND_POS);
    check("ee_ang_within_bound", max_ang_err <= EE_BOUND_ANG);
    check("constraints_hold_during_blend", worst_stretch < 0.02f);

    /* no jitter: monotone ramp, no velocity reversals, bounded accel */
    printf("max_ee_accel=%.5f (bound %.2f), reversals=%d, frames=%d\n",
           b.max_accel, JITTER_ACCEL, b.reversals, frames);
    check("alpha_monotone", b.alpha_regressed == 0);
    check("no_velocity_reversal", b.reversals == 0);
    check("accel_bounded", b.max_accel <= JITTER_ACCEL);

    /* end pose: root weight 1 keeps it at the animation pose; the end
       bone's animation influence equals its configured weight */
    check("root_held", fabsf(b.out[0].x - anim[0].x) < 1e-5f &&
                       fabsf(b.out[0].y - anim[0].y) < 1e-5f);
    float end_e = weight[N_BONES - 1];
    float expect_x = end_e * anim[N_BONES - 1].x +
                     (1.0f - end_e) * b.phys.p[N_BONES - 1].x;
    check("end_bone_weight", fabsf(b.out[N_BONES - 1].x - expect_x) < 1e-5f);

    /* stepping after completion is rejected */
    check("step_after_finish",
          blend_step(&b, DT, 0.0f, -9.8f) == BLEND_ERR_STATE);

    /* edge case: every bone fully animated locks the output to the anim
       pose even though the physics chain keeps integrating underneath */
    float w_anim[N_BONES], w_phys[N_BONES];
    for (int i = 0; i < N_BONES; i++) { w_anim[i] = 1.0f; w_phys[i] = 0.0f; }
    Blend locked;
    rc = blend_start(&locked, anim, phys0, w_anim, N_BONES, 0.6f);
    check("locked_start", rc == BLEND_OK);
    while (locked.active) blend_step(&locked, DT, 0.0f, -9.8f);
    int locked_ok = 1;
    for (int i = 0; i < N_BONES; i++) {
        if (fabsf(locked.out[i].x - anim[i].x) > 1e-5f ||
            fabsf(locked.out[i].y - anim[i].y) > 1e-5f ||
            fabsf(wrap_angle(locked.out[i].ang - anim[i].ang)) > 1e-5f)
            locked_ok = 0;
    }
    check("full_anim_weight_locks_pose", locked_ok);

    /* edge case: zero animation weight and no gravity means the output
       tracks a stationary chain, so the jitter metrics must read exactly
       zero (a genuine no-drive scenario cannot wobble) */
    Blend free_chain;
    rc = blend_start(&free_chain, anim, phys0, w_phys, N_BONES, 0.6f);
    check("free_start", rc == BLEND_OK);
    while (free_chain.active) blend_step(&free_chain, DT, 0.0f, 0.0f);
    check("no_drive_no_jitter",
          free_chain.max_accel < 1e-6f && free_chain.reversals == 0);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
