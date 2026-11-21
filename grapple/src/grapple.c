
/* grapple: aim-ray hook selection, rope-length constraint and swing tension.
   The rope is a one-sided (inextensible) distance constraint: it can pull the
   body toward the anchor but never pushes it, so slack rope applies no force. */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define MAX_HOOKS 16
#define MIN_ROPE 1e-3f      /* shorter than this is treated as a zero rope */
#define ROPE_EPS 1e-4f      /* allowed numerical overshoot of the rope length */

typedef struct { float x, y; } V2;

static V2 v2(float x, float y) { V2 r = { x, y }; return r; }
static V2 vadd(V2 a, V2 b) { V2 r = { a.x + b.x, a.y + b.y }; return r; }
static V2 vsub(V2 a, V2 b) { V2 r = { a.x - b.x, a.y - b.y }; return r; }
static V2 vmul(V2 a, float s) { V2 r = { a.x * s, a.y * s }; return r; }
static float vdot(V2 a, V2 b) { return a.x * b.x + a.y * b.y; }
static float vlen(V2 a) { return sqrtf(a.x * a.x + a.y * a.y); }
static V2 vnorm(V2 a) {
    float l = vlen(a);
    V2 r = { 0.0f, 0.0f };
    if (l > 0.0f) { r.x = a.x / l; r.y = a.y / l; }
    return r;
}
static int vfinite(V2 a) { return isfinite(a.x) && isfinite(a.y); }

typedef enum {
    GRAPPLE_OK = 0,
    GRAPPLE_MISS = -1,          /* aim ray selected no hook */
    GRAPPLE_BAD_ANCHOR = -2,    /* anchor is not finite */
    GRAPPLE_ZERO_ROPE = -3      /* hook coincides with the body */
} GrappleResult;

typedef struct {
    V2 hooks[MAX_HOOKS];
    int hook_count;
    V2 pos, vel, gravity;
    V2 anchor;
    float rope;         /* current rope length, 0 when released */
    float tension;      /* last computed tension, NAN when slack/released */
    int attached;
    int taut;
} Grapple;

static void grapple_init(Grapple* g, V2 pos, V2 vel, V2 gravity) {
    memset(g, 0, sizeof(*g));
    g->pos = pos;
    g->vel = vel;
    g->gravity = gravity;
    g->tension = NAN;
}

static int grapple_add_hook(Grapple* g, V2 hook) {
    if (g->hook_count >= MAX_HOOKS) return -1;
    if (!vfinite(hook)) return -1;
    g->hooks[g->hook_count] = hook;
    return g->hook_count++;
}

static int grapple_is_attached(const Grapple* g) { return g->attached; }
static int grapple_is_taut(const Grapple* g) { return g->attached && g->taut; }
static float grapple_tension(const Grapple* g) { return g->attached ? g->tension : NAN; }

/* pick the closest hook along the aim ray within range and angular tolerance */
static int grapple_pick(const Grapple* g, V2 origin, V2 dir, float max_range, float tolerance) {
    if (!vfinite(origin) || !vfinite(dir)) return -1;
    if (!(max_range > 0.0f) || !isfinite(max_range)) return -1;
    if (!(tolerance >= 0.0f) || !isfinite(tolerance)) return -1;
    float dl = vlen(dir);
    if (!(dl > 1e-6f)) return -1;
    V2 u = vmul(dir, 1.0f / dl);

    int best = -1;
    float best_s = 0.0f, best_perp = 0.0f;
    for (int i = 0; i < g->hook_count; i++) {
        V2 rel = vsub(g->hooks[i], origin);
        float s = vdot(rel, u);
        if (s < 0.0f || s > max_range) continue;         /* behind or out of range */
        float perp = vlen(vsub(rel, vmul(u, s)));
        if (perp > tolerance) continue;                  /* outside aim cone */
        if (best < 0 || s < best_s - 1e-6f ||
            (fabsf(s - best_s) <= 1e-6f && perp < best_perp)) {
            best = i;
            best_s = s;
            best_perp = perp;
        }
    }
    return best;
}

static GrappleResult grapple_attach(Grapple* g, V2 hook) {
    if (!vfinite(hook) || !vfinite(g->pos)) return GRAPPLE_BAD_ANCHOR;
    float d = vlen(vsub(hook, g->pos));
    if (!(d > MIN_ROPE)) return GRAPPLE_ZERO_ROPE;       /* zero / degenerate rope */
    g->anchor = hook;
    g->rope = d;
    g->attached = 1;
    g->taut = 0;
    g->tension = NAN;
    return GRAPPLE_OK;
}

static GrappleResult grapple_fire(Grapple* g, V2 origin, V2 dir, float max_range, float tolerance) {
    int i = grapple_pick(g, origin, dir, max_range, tolerance);
    if (i < 0) return GRAPPLE_MISS;
    return grapple_attach(g, g->hooks[i]);
}

static void grapple_release(Grapple* g) {
    g->attached = 0;
    g->taut = 0;
    g->rope = 0.0f;
    g->tension = NAN;
}

/* radial overshoot beyond the rope; the solver must keep this at zero */
static float grapple_rope_violation(const Grapple* g) {
    if (!g->attached) return 0.0f;
    float over = vlen(vsub(g->pos, g->anchor)) - g->rope;
    return over > 0.0f ? over : 0.0f;
}

/* tension of a taut rope: centripetal term from tangential speed plus the
   radial component of gravity. n points from the anchor to the body. */
static float tension_at(V2 anchor, float rope, V2 pos, V2 vel, V2 gravity) {
    if (!(rope > 0.0f)) return NAN;
    V2 r = vsub(pos, anchor);
    float d = vlen(r);
    if (!(d > 1e-6f)) return NAN;
    V2 n = vmul(r, 1.0f / d);
    float vr = vdot(vel, n);
    float vt2 = vdot(vel, vel) - vr * vr;
    if (vt2 < 0.0f) vt2 = 0.0f;
    return vt2 / rope + vdot(gravity, n);
}

static void grapple_step(Grapple* g, float dt) {
    if (!(dt > 0.0f)) return;

    g->vel = vadd(g->vel, vmul(g->gravity, dt));
    g->pos = vadd(g->pos, vmul(g->vel, dt));
    if (!g->attached) { g->taut = 0; g->tension = NAN; return; }

    V2 r = vsub(g->pos, g->anchor);
    float d = vlen(r);
    if (!(d > 1e-6f)) { g->taut = 0; g->tension = NAN; return; }  /* body at anchor */
    V2 n = vmul(r, 1.0f / d);
    float vr = vdot(g->vel, n);
    float t = tension_at(g->anchor, g->rope, g->pos, g->vel, g->gravity);

    /* Stretched past the rope, or taut (tension > 0): snap onto the circle and
       cancel radial velocity. Tension <= 0 inside the circle means slack, and a
       rope that can only pull applies no force there. */
    if (d > g->rope || t > 0.0f) {
        g->pos = vadd(g->anchor, vmul(n, g->rope));
        g->vel = vsub(g->vel, vmul(n, vr));
        g->taut = 1;
        g->tension = t;
    } else {
        g->taut = 0;
        g->tension = NAN;
    }
}

/* authored swing curve: small-angle pendulum, L = 4, anchor (0, 6), sampled
   every 0.125 s from t = 0 to t = 2. Used as the hand-tuned reference. */
static const float swing_keys[17][3] = {
    { 0.000f,  0.399334f, 2.019983f },
    { 0.125f,  0.391739f, 2.019229f },
    { 0.250f,  0.369237f, 2.017078f },
    { 0.375f,  0.332665f, 2.013857f },
    { 0.500f,  0.283389f, 2.010051f },
    { 0.625f,  0.223264f, 2.006236f },
    { 0.750f,  0.154571f, 2.002988f },
    { 0.875f,  0.079935f, 2.000799f },
    { 1.000f,  0.002219f, 2.000001f },
    { 1.125f, -0.075581f, 2.000714f },
    { 1.250f, -0.150471f, 2.002831f },
    { 1.375f, -0.219574f, 2.006031f },
    { 1.500f, -0.280249f, 2.009830f },
    { 1.625f, -0.330194f, 2.013652f },
    { 1.750f, -0.367529f, 2.016920f },
    { 1.875f, -0.390856f, 2.019142f },
    { 2.000f, -0.399309f, 2.019981f },
};

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    /* ---- aim-ray hook selection ---- */
    Grapple aim;
    grapple_init(&aim, v2(0, -1), v2(0, 0), v2(0, -9.8f));
    check("hook_add_a", grapple_add_hook(&aim, v2(0, 5)) == 0);
    check("hook_add_b", grapple_add_hook(&aim, v2(10, 0)) == 1);
    check("hook_add_c", grapple_add_hook(&aim, v2(3, 3)) == 2);
    check("hook_reject_nan", grapple_add_hook(&aim, v2(NAN, 0)) == -1);

    check("aim_x_picks_b", grapple_pick(&aim, v2(0, 0), v2(1, 0), 12.0f, 0.5f) == 1);
    check("aim_up_picks_a", grapple_pick(&aim, v2(0, 0), v2(0, 1), 12.0f, 0.5f) == 0);
    check("aim_diag_picks_c", grapple_pick(&aim, v2(0, 0), v2(1, 1), 5.0f, 0.5f) == 2);
    check("aim_cone_miss", grapple_pick(&aim, v2(0, 0), v2(0.98f, 0.2f), 12.0f, 0.5f) == -1);
    check("aim_range_reject", grapple_pick(&aim, v2(0, 0), v2(1, 1), 4.0f, 0.5f) == -1);
    check("aim_behind_reject", grapple_pick(&aim, v2(0, 0), v2(-1, 0), 12.0f, 0.5f) == -1);
    check("aim_zero_dir", grapple_pick(&aim, v2(0, 0), v2(0, 0), 12.0f, 0.5f) == -1);
    check("aim_nan_origin", grapple_pick(&aim, v2(NAN, 0), v2(1, 0), 12.0f, 0.5f) == -1);

    Grapple fire;
    grapple_init(&fire, v2(0, 0), v2(0, 0), v2(0, -9.8f));
    grapple_add_hook(&fire, v2(10, 0));
    check("fire_ok", grapple_fire(&fire, v2(0, 0), v2(1, 0), 12.0f, 0.5f) == GRAPPLE_OK);
    check("fire_attached", grapple_is_attached(&fire) == 1);
    check("fire_rope_len", fabsf(fire.rope - 10.0f) < 1e-5f);
    check("fire_miss_keeps_state",
          grapple_fire(&fire, v2(0, 0), v2(-1, 0), 12.0f, 0.5f) == GRAPPLE_MISS &&
          grapple_is_attached(&fire) == 1);

    /* ---- invalid anchor / zero rope rejection ---- */
    Grapple bad;
    grapple_init(&bad, v2(1, 2), v2(0, 0), v2(0, -9.8f));
    check("attach_ok", grapple_attach(&bad, v2(0, 6)) == GRAPPLE_OK);
    check("attach_rope_len", fabsf(bad.rope - sqrtf(17.0f)) < 1e-5f);
    check("attach_nan_reject", grapple_attach(&bad, v2(NAN, 6)) == GRAPPLE_BAD_ANCHOR);
    check("attach_inf_reject", grapple_attach(&bad, v2(INFINITY, 6)) == GRAPPLE_BAD_ANCHOR);
    check("attach_kept_after_reject", grapple_is_attached(&bad) == 1 && fabsf(bad.rope - sqrtf(17.0f)) < 1e-5f);

    Grapple zero;
    grapple_init(&zero, v2(2, 3), v2(0, 0), v2(0, -9.8f));
    check("attach_zero_rope_reject", grapple_attach(&zero, v2(2, 3)) == GRAPPLE_ZERO_ROPE);
    check("attach_near_zero_reject", grapple_attach(&zero, v2(2, 3.0005f)) == GRAPPLE_ZERO_ROPE);
    check("zero_not_attached", grapple_is_attached(&zero) == 0);

    Grapple zfire;
    grapple_init(&zfire, v2(0, 0), v2(0, 0), v2(0, -9.8f));
    grapple_add_hook(&zfire, v2(0, 0));
    check("fire_zero_rope_reject",
          grapple_fire(&zfire, v2(0, 0), v2(1, 0), 12.0f, 0.5f) == GRAPPLE_ZERO_ROPE &&
          grapple_is_attached(&zfire) == 0);

    /* ---- tension from tangential speed ---- */
    V2 grav = v2(0, -9.8f);
    float t1 = tension_at(v2(0, 0), 2.0f, v2(0, -2), v2(3, 0), grav);
    check("tension_value", fabsf(t1 - 14.3f) < 1e-4f);
    float t2 = tension_at(v2(0, 0), 2.0f, v2(0, -2), v2(3, 1.5f), grav);
    check("tension_ignores_radial_vel", fabsf(t2 - t1) < 1e-5f);
    float t3 = tension_at(v2(0, 0), 2.0f, v2(0, -2), v2(6, 0), grav);
    check("tension_scales_v_squared", fabsf((t3 - 9.8f) - 4.0f * (t1 - 9.8f)) < 1e-5f);
    check("tension_slack_is_nan", !isfinite(tension_at(v2(0, 0), 0.0f, v2(0, -2), v2(3, 0), grav)));

    /* ---- rope is one-sided: no force when slack, no push ---- */
    Grapple slack;
    grapple_init(&slack, v2(0, 0), v2(0, 0), v2(0, 0));
    grapple_attach(&slack, v2(0, 0));   /* rejected? pos == hook, rope 0 */
    check("slack_setup", !grapple_is_attached(&slack));
    grapple_init(&slack, v2(0, -1), v2(0, 0), v2(0, 0));
    check("slack_attach", grapple_attach(&slack, v2(0, 0)) == GRAPPLE_OK);
    slack.pos = v2(0, -0.5f);           /* inside the circle, zero velocity */
    grapple_step(&slack, 0.1f);
    check("slack_no_push", fabsf(slack.pos.y + 0.5f) < 1e-6f && !grapple_is_taut(&slack));
    slack.vel = v2(0, 1);               /* inward: rope must not stop it */
    grapple_step(&slack, 0.1f);
    check("slack_allows_inward", fabsf(slack.pos.y + 0.4f) < 1e-5f && !grapple_is_taut(&slack));

    Grapple top;
    grapple_init(&top, v2(0, 4), v2(1, 0), v2(0, -9.8f));
    check("top_attach", grapple_attach(&top, v2(0, 0)) == GRAPPLE_OK);
    grapple_step(&top, 1.0f / 240.0f);
    check("top_slack_no_push", !grapple_is_taut(&top) &&
          vlen(vsub(top.pos, top.anchor)) < top.rope + 1e-6f);

    /* ---- rope length cap with outward velocity ---- */
    Grapple cap;
    grapple_init(&cap, v2(0, -4), v2(0, -20), v2(0, -9.8f));
    check("cap_attach", grapple_attach(&cap, v2(0, 0)) == GRAPPLE_OK);
    grapple_step(&cap, 1.0f / 60.0f);
    V2 cap_n = vnorm(vsub(cap.pos, cap.anchor));
    check("cap_distance", grapple_rope_violation(&cap) <= ROPE_EPS);
    check("cap_radial_zero", fabsf(vdot(cap.vel, cap_n)) < 1e-5f);
    check("cap_taut", grapple_is_taut(&cap) && grapple_tension(&cap) > 0.0f);

    /* ---- authored swing trajectory comparison ---- */
    Grapple sw;
    grapple_init(&sw, v2(swing_keys[0][1], swing_keys[0][2]), v2(0, 0), v2(0, -9.8f));
    check("swing_attach", grapple_attach(&sw, v2(0, 6)) == GRAPPLE_OK);
    check("swing_rope", fabsf(sw.rope - 4.0f) < 1e-4f);
    float max_dev = 0.0f;
    int violations = 0;
    float min_t = INFINITY, max_t = -INFINITY;
    for (int step = 1; step <= 480; step++) {
        grapple_step(&sw, 1.0f / 240.0f);
        if (grapple_rope_violation(&sw) > ROPE_EPS) violations++;
        if (grapple_is_taut(&sw)) {
            float t = grapple_tension(&sw);
            if (t < min_t) min_t = t;
            if (t > max_t) max_t = t;
        }
        if (step % 30 == 0) {
            int k = step / 30;
            float dev = vlen(vsub(sw.pos, v2(swing_keys[k][1], swing_keys[k][2])));
            if (dev > max_dev) max_dev = dev;
        }
    }
    printf("swing_max_dev=%.6f\n", max_dev);
    check("swing_dev_under_bound", max_dev <= 0.02f);
    check("swing_rope_violations_zero", violations == 0);
    check("swing_tension_positive", min_t > 0.0f && isfinite(min_t) && max_t > min_t);

    /* ---- release / re-attach ---- */
    Grapple mode;
    grapple_init(&mode, v2(0, -2), v2(5, 0), v2(0, -9.8f));
    check("mode_attach", grapple_attach(&mode, v2(0, 0)) == GRAPPLE_OK);
    for (int i = 0; i < 120; i++) grapple_step(&mode, 1.0f / 240.0f);
    check("mode_taut", grapple_is_taut(&mode) && grapple_tension(&mode) > 0.0f);
    grapple_release(&mode);
    check("mode_released", !grapple_is_attached(&mode) && !isfinite(grapple_tension(&mode)));
    V2 v_release = mode.vel;
    for (int i = 0; i < 120; i++) grapple_step(&mode, 1.0f / 240.0f);
    check("mode_free_vx", fabsf(mode.vel.x - v_release.x) < 1e-4f);
    check("mode_free_gravity", fabsf((mode.vel.y - v_release.y) + 9.8f * 0.5f) < 1e-4f);
    check("mode_no_rope_after_release", grapple_rope_violation(&mode) == 0.0f);
    check("mode_reattach", grapple_attach(&mode, v2(10, 10)) == GRAPPLE_OK &&
          fabsf(mode.rope - vlen(vsub(v2(10, 10), mode.pos))) < 1e-5f);

    /* ---- long stress run: violations must stay at zero ---- */
    Grapple hard;
    grapple_init(&hard, v2(1, -3), v2(12, 6), v2(0, -9.8f));
    check("stress_attach", grapple_attach(&hard, v2(0, 5)) == GRAPPLE_OK);
    int hard_viol = 0;
    for (int i = 0; i < 2000; i++) {
        grapple_step(&hard, 1.0f / 240.0f);
        if (grapple_rope_violation(&hard) > ROPE_EPS) hard_viol++;
        if (i % 400 == 399) hard.vel = vadd(hard.vel, v2(3, -3));
    }
    check("stress_violations_zero", hard_viol == 0);
    check("stress_finite", vfinite(hard.pos) && vfinite(hard.vel));

    /* ---- fixed capacity ---- */
    Grapple full;
    grapple_init(&full, v2(0, 0), v2(0, 0), v2(0, 0));
    int rc = 0;
    for (int i = 0; i <= MAX_HOOKS; i++) rc = grapple_add_hook(&full, v2((float)i, 1));
    check("capacity", rc == -1 && full.hook_count == MAX_HOOKS);

    /* ---- non-positive dt is ignored ---- */
    Grapple guard;
    grapple_init(&guard, v2(1, 1), v2(2, 2), v2(0, -9.8f));
    grapple_attach(&guard, v2(0, 0));
    V2 p0 = guard.pos;
    grapple_step(&guard, 0.0f);
    grapple_step(&guard, -0.1f);
    check("dt_guard", fabsf(guard.pos.x - p0.x) < 1e-9f && fabsf(guard.pos.y - p0.y) < 1e-9f);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
