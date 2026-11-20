
/* ik: analytic two-bone IK with pole vector and joint angle limits.
   The raw solver refuses unreachable targets; the policy layer below clamps
   the request onto the reachable annulus, walks there smoothly instead of
   snapping, and exposes the decision through IkStatus. */
#include <stdio.h>
#include <string.h>
#include <math.h>

typedef struct { float x, y; } V2;

#define IK_PI 3.14159265358979f

/* shared with the raw solver's own guards so both give the same verdict */
#define IK_REACH_EPS 1e-5f
#define IK_TARGET_EPS 1e-6f

/* clamped targets travel at most this many units per second */
#define IK_DEFAULT_SMOOTH_SPEED 8.0f
/* used when the caller passes a non-positive or absurd dt */
#define IK_NOMINAL_DT (1.0f / 60.0f)
#define IK_MAX_DT 0.1f

static float vlen(V2 v) { return sqrtf(v.x * v.x + v.y * v.y); }
static V2 vsub(V2 a, V2 b) { V2 r = { a.x - b.x, a.y - b.y }; return r; }
static V2 vadd(V2 a, V2 b) { V2 r = { a.x + b.x, a.y + b.y }; return r; }
static V2 vscale(V2 v, float s) { V2 r = { v.x * s, v.y * s }; return r; }
static V2 vnorm(V2 v) { float l = vlen(v); V2 r = { 0, 0 }; if (l > 0) { r.x = v.x / l; r.y = v.y / l; } return r; }

typedef struct {
    int ok;
    V2 joint;
    float a1;   /* angle of first bone */
    float a2;   /* angle of second bone */
} IkResult;

/* where the requested target sits relative to the reachable annulus */
typedef enum {
    IK_REACHABLE = 0,   /* |l1-l2| <= dist <= l1+l2 */
    IK_TOO_CLOSE = 1,   /* inside the inner hole */
    IK_UNREACHABLE = 2  /* beyond the outer radius */
} IkReach;

/* why a frame produced no pose */
typedef enum {
    IK_ERR_NONE = 0,
    IK_ERR_DEGENERATE_LIMB = 1,    /* a bone has zero (or negative) length */
    IK_ERR_DEGENERATE_TARGET = 2   /* target coincides with the root */
} IkError;

/* queryable outcome of one solve */
typedef struct {
    int ok;             /* 1 when a pose was produced */
    int status;         /* IkReach verdict for the request */
    int err;            /* IkError; takes precedence over ok/status */
    int clamped;        /* request was outside the annulus and got clamped */
    int smoothing;      /* clamped target is still transitioning */
    float t;            /* progress 0..1 of the current transition, 1 = settled */
    float dist;         /* ||request - root|| */
    V2 request;         /* target exactly as the caller asked for it */
    V2 applied;         /* target the solver actually used this frame */
} IkStatus;

/* target smoothing state: the applied target approaches the clamped target
   along the annulus, never moving more than speed*dt per step */
typedef struct {
    V2 pos;             /* last target actually applied */
    float speed;        /* units per second */
    float start_gap;    /* gap at the start of the current transition */
    int active;
    float t;
    int seeded;         /* 0 until the first solve establishes a pose */
} IkSmoother;

typedef struct {
    V2 root;
    float l1, l2;
    V2 pole;
    IkSmoother smooth;
    IkStatus status;
} IkChain;

/* ---- raw solver: unchanged behaviour ---------------------------------- */

/* solve for target with bone lengths l1,l2; pole decides bend direction */
static IkResult solve_two_bone(V2 root, V2 target, float l1, float l2, V2 pole) {
    IkResult r;
    memset(&r, 0, sizeof(r));
    V2 d = vsub(target, root);
    float dist = vlen(d);
    if (l1 <= 0 || l2 <= 0) return r;
    if (dist > l1 + l2 + 1e-5f) return r;                 /* unreachable */
    if (dist < fabsf(l1 - l2) - 1e-5f) return r;          /* too close */
    if (dist < 1e-6f) return r;

    float cos_a = (l1 * l1 + dist * dist - l2 * l2) / (2 * l1 * dist);
    if (cos_a > 1.0f) cos_a = 1.0f;
    if (cos_a < -1.0f) cos_a = -1.0f;
    float a = acosf(cos_a);

    V2 dir = vnorm(d);
    float side = (dir.x * pole.y - dir.y * pole.x) >= 0 ? 1.0f : -1.0f;
    float ang = atan2f(dir.y, dir.x) + side * a;

    r.joint.x = root.x + l1 * cosf(ang);
    r.joint.y = root.y + l1 * sinf(ang);
    r.a1 = ang;
    V2 d2 = vsub(target, r.joint);
    r.a2 = atan2f(d2.y, d2.x);
    r.ok = 1;
    return r;
}

/* ---- reachability ------------------------------------------------------ */

/* same thresholds as solve_two_bone, but reported instead of swallowed */
static int reach_level(V2 root, V2 target, float l1, float l2) {
    if (l1 <= 0.0f || l2 <= 0.0f) return IK_UNREACHABLE;
    float dist = vlen(vsub(target, root));
    if (dist > l1 + l2 + IK_REACH_EPS) return IK_UNREACHABLE;
    if (dist < fabsf(l1 - l2) - IK_REACH_EPS) return IK_TOO_CLOSE;
    return IK_REACHABLE;
}

/* nearest point on the reachable annulus, keeping the request's direction */
static V2 clamp_to_ring(V2 root, V2 target, float l1, float l2) {
    float lo = fabsf(l1 - l2);
    float hi = l1 + l2;
    V2 d = vsub(target, root);
    float dist = vlen(d);
    V2 dir = vnorm(d);
    if (dir.x == 0.0f && dir.y == 0.0f) dir = (V2){ 1, 0 };
    float r = dist < lo ? lo : hi;
    if (dist <= hi && dist >= lo) r = dist;   /* already on the annulus */
    return vadd(root, vscale(dir, r));
}

static void ik_chain_init(IkChain* c, V2 root, float l1, float l2, V2 pole) {
    memset(c, 0, sizeof(*c));
    c->root = root;
    c->l1 = l1;
    c->l2 = l2;
    c->pole = pole;
    c->smooth.speed = IK_DEFAULT_SMOOTH_SPEED;
    c->smooth.pos = root;
    c->smooth.t = 1.0f;
    c->status.err = IK_ERR_NONE;
    c->status.t = 1.0f;
}

/* one smoothing step along the annulus, keeping every intermediate point
   exactly on/inside the reachable annulus */
static V2 ik_smooth_step(IkSmoother* s, V2 root, float lo, float hi, V2 desired, float dt) {
    if (!s->active) {
        float gap = vlen(vsub(desired, s->pos));
        s->start_gap = gap;
        s->t = 0.0f;
        s->active = 1;
        if (gap <= IK_TARGET_EPS) {
            s->pos = desired;
            s->active = 0;
            s->t = 1.0f;
            return s->pos;
        }
    }

    float max_step = s->speed * dt;
    float remaining = vlen(vsub(desired, s->pos));
    if (remaining <= max_step) {
        s->pos = desired;
        s->active = 0;
        s->t = 1.0f;
        return s->pos;
    }

    /* polar interpolation: radius stays inside [lo,hi] by construction, so
       the smoothed path can never cut through the unreachable inner hole */
    V2 pc = vsub(s->pos, root);
    V2 pd = vsub(desired, root);
    float rc = vlen(pc), rd = vlen(pd);
    float ac = atan2f(pc.y, pc.x), ad = atan2f(pd.y, pd.x);
    float da = ad - ac;
    if (da > IK_PI) da -= 2.0f * IK_PI;
    if (da < -IK_PI) da += 2.0f * IK_PI;
    float rref = rc > rd ? rc : rd;
    if (rref < IK_TARGET_EPS) rref = IK_TARGET_EPS;
    float max_ang = max_step / rref;
    float sa = da > max_ang ? max_ang : (da < -max_ang ? -max_ang : da);
    float dr = rd - rc;
    float sr = dr > max_step ? max_step : (dr < -max_step ? -max_step : dr);

    /* radial and angular budgets are capped separately; fold them so the
       chord actually travelled this frame never exceeds speed*dt */
    float rmax = rc > rd ? rc : rd;
    float budget = fabsf(sr) + rmax * fabsf(sa);
    if (budget > max_step) {
        float k = max_step / budget;
        sr *= k;
        sa *= k;
    }

    float nr = rc + sr;
    if (nr < lo) nr = lo;
    if (nr > hi) nr = hi;
    s->pos = (V2){ root.x + nr * cosf(ac + sa), root.y + nr * sinf(ac + sa) };

    remaining = vlen(vsub(desired, s->pos));
    s->t = (s->start_gap > IK_TARGET_EPS) ? 1.0f - remaining / s->start_gap : 1.0f;
    if (s->t < 0.0f) s->t = 0.0f;
    if (s->t > 1.0f) s->t = 1.0f;
    return s->pos;
}

/* query the outcome of the last solve (by value; never NULL) */
static IkStatus ik_chain_status(const IkChain* c) {
    IkStatus empty;
    memset(&empty, 0, sizeof(empty));
    empty.t = 1.0f;
    return c ? c->status : empty;
}

/* query the last applied target, or 0 when nothing has been solved yet */
static int ik_chain_target(const IkChain* c, V2* out) {
    if (!c || !out || !c->smooth.seeded) return 0;
    *out = c->smooth.pos;
    return 1;
}

/* speed of clamped-target transitions, in units per second */
static void ik_chain_set_smooth_speed(IkChain* c, float speed) {
    if (c && speed > 0.0f) c->smooth.speed = speed;
}

/* solve one frame: validate, classify, clamp, smooth, query. */
static IkResult ik_chain_solve(IkChain* c, V2 request, float dt) {
    IkResult none;
    memset(&none, 0, sizeof(none));
    memset(&c->status, 0, sizeof(c->status));
    c->status.request = request;
    c->status.applied = request;
    c->status.t = 1.0f;

    if (!(c->l1 > 0.0f) || !(c->l2 > 0.0f)) {
        c->status.err = IK_ERR_DEGENERATE_LIMB;
        c->status.status = IK_UNREACHABLE;
        c->status.dist = vlen(vsub(request, c->root));
        return none;
    }

    float dist = vlen(vsub(request, c->root));
    c->status.dist = dist;
    c->status.status = reach_level(c->root, request, c->l1, c->l2);

    if (dist < IK_TARGET_EPS) {
        c->status.err = IK_ERR_DEGENERATE_TARGET;
        return none;
    }

    if (!(dt > 0.0f)) dt = IK_NOMINAL_DT;   /* also rejects NaN */
    if (dt > IK_MAX_DT) dt = IK_MAX_DT;

    if (c->status.status == IK_REACHABLE) {
        /* reachable targets go straight to the solver, bit for bit as before */
        c->status.applied = request;
        c->status.clamped = 0;
        c->status.smoothing = 0;
        c->status.t = 1.0f;
        c->smooth.pos = request;
        c->smooth.active = 0;
        c->smooth.t = 1.0f;
        c->smooth.seeded = 1;
        IkResult r = solve_two_bone(c->root, request, c->l1, c->l2, c->pole);
        c->status.ok = r.ok;
        return r;
    }

    /* outside the annulus: clamp radiply, then let the smoother walk there */
    float lo = fabsf(c->l1 - c->l2);
    float hi = c->l1 + c->l2;
    V2 desired = clamp_to_ring(c->root, request, c->l1, c->l2);
    c->status.clamped = 1;

    if (!c->smooth.seeded) {
        /* no previous pose to transition from */
        c->smooth.pos = desired;
        c->smooth.seeded = 1;
        c->smooth.active = 0;
        c->smooth.t = 1.0f;
    } else {
        ik_smooth_step(&c->smooth, c->root, lo, hi, desired, dt);
    }

    c->status.applied = c->smooth.pos;
    c->status.smoothing = c->smooth.active;
    c->status.t = c->smooth.active ? c->smooth.t : 1.0f;

    IkResult r = solve_two_bone(c->root, c->status.applied, c->l1, c->l2, c->pole);
    c->status.ok = r.ok;
    return r;
}

/* ---- self-test --------------------------------------------------------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static int veq(V2 a, V2 b) { return a.x == b.x && a.y == b.y; }

int main(void) {
    V2 root = { 0, 0 }, target = { 1, 0 }, pole = { 0, 1 };
    IkResult r = solve_two_bone(root, target, 1.0f, 1.0f, pole);
    check("reachable", r.ok == 1);
    check("joint_on_circle", fabsf(vlen(vsub(r.joint, root)) - 1.0f) < 1e-4f);
    check("end_matches", fabsf(vlen(vsub(r.joint, target)) - 1.0f) < 1e-4f);
    check("pole_up", r.joint.y > 0.0f);

    IkResult r2 = solve_two_bone(root, target, 1.0f, 1.0f, (V2){ 0, -1 });
    check("pole_down", r2.joint.y < 0.0f);

    IkResult far = solve_two_bone(root, (V2){ 5, 0 }, 1.0f, 1.0f, pole);
    check("unreachable", far.ok == 0);

    IkResult close = solve_two_bone(root, (V2){ 0.1f, 0 }, 1.0f, 0.2f, pole);
    check("too_close", close.ok == 0);

    IkResult zero = solve_two_bone(root, target, 0.0f, 1.0f, pole);
    check("zero_len", zero.ok == 0);

    /* analytic classification agrees with the raw solver's acceptance */
    check("reach_class_far", reach_level(root, (V2){ 5, 0 }, 1.0f, 1.0f) == IK_UNREACHABLE);
    check("reach_class_close", reach_level(root, (V2){ 0.1f, 0 }, 1.0f, 0.2f) == IK_TOO_CLOSE);
    check("reach_class_in", reach_level(root, (V2){ 1, 0 }, 1.0f, 1.0f) == IK_REACHABLE);

    int class_agrees = 1;
    for (int i = 1; i < 400; i++) {
        float d = 0.01f + 0.01f * (float)i;
        V2 t = { d, 0.3f };
        IkResult rr = solve_two_bone(root, t, 1.0f, 1.5f, pole);
        int lvl = reach_level(root, t, 1.0f, 1.5f);
        if ((rr.ok == 1) != (lvl == IK_REACHABLE)) { class_agrees = 0; break; }
    }
    check("class_matches_solver", class_agrees);

    /* clamping lands exactly on the annulus along the request direction */
    V2 c1 = clamp_to_ring(root, (V2){ 5, 0 }, 1.0f, 1.0f);
    check("clamp_far_on_ring", fabsf(vlen(vsub(c1, root)) - 2.0f) < 1e-6f && fabsf(c1.x - 2.0f) < 1e-6f);
    V2 c2 = clamp_to_ring(root, (V2){ -0.1f, 0 }, 1.0f, 0.2f);
    check("clamp_close_on_ring", fabsf(vlen(vsub(c2, root)) - 0.8f) < 1e-6f && c2.x < 0.0f);

    /* unreachable request: pose is produced, status is queryable */
    IkChain ch;
    ik_chain_init(&ch, root, 1.0f, 1.0f, pole);
    IkResult cr = ik_chain_solve(&ch, (V2){ 5, 0 }, 1.0f / 60.0f);
    check("chain_unreachable_solves", cr.ok == 1);
    check("chain_status_unreachable",
          ch.status.status == IK_UNREACHABLE && ch.status.clamped == 1 &&
          ch.status.err == IK_ERR_NONE && ch.status.ok == 1);
    check("chain_applied_on_ring", fabsf(vlen(vsub(ch.status.applied, root)) - 2.0f) < 1e-6f);
    check("chain_request_preserved", ch.status.request.x == 5.0f && ch.status.request.y == 0.0f);
    check("chain_applied_differs", !veq(ch.status.applied, ch.status.request));
    check("chain_status_dist", fabsf(ch.status.dist - 5.0f) < 1e-6f);

    /* degenerate limb and degenerate target both report errors */
    IkChain zb;
    ik_chain_init(&zb, root, 0.0f, 1.0f, pole);
    IkResult zbr = ik_chain_solve(&zb, target, 1.0f / 60.0f);
    check("err_zero_bone", zbr.ok == 0 && zb.status.err == IK_ERR_DEGENERATE_LIMB && zb.status.ok == 0);
    ik_chain_init(&zb, root, 1.0f, -0.5f, pole);
    zbr = ik_chain_solve(&zb, target, 1.0f / 60.0f);
    check("err_negative_bone", zbr.ok == 0 && zb.status.err == IK_ERR_DEGENERATE_LIMB);

    IkChain zt;
    ik_chain_init(&zt, root, 1.0f, 1.0f, pole);
    IkResult ztr = ik_chain_solve(&zt, root, 1.0f / 60.0f);
    check("err_zero_target", ztr.ok == 0 && zt.status.err == IK_ERR_DEGENERATE_TARGET && zt.status.ok == 0);
    ik_chain_init(&zt, root, 1.0f, 0.2f, pole);
    ztr = ik_chain_solve(&zt, root, 1.0f / 60.0f);
    check("err_zero_target_small_inner", ztr.ok == 0 && zt.status.err == IK_ERR_DEGENERATE_TARGET);

    /* reachable poses are bit-identical to the raw solver, as before */
    IkChain ri;
    ik_chain_init(&ri, root, 1.0f, 1.0f, pole);
    IkResult a = ik_chain_solve(&ri, target, 1.0f / 60.0f);
    IkResult b = solve_two_bone(root, target, 1.0f, 1.0f, pole);
    check("reachable_bitwise",
          a.ok == b.ok && a.joint.x == b.joint.x && a.joint.y == b.joint.y &&
          a.a1 == b.a1 && a.a2 == b.a2);
    check("reachable_applied_is_request",
          veq(ri.status.applied, target) && ri.status.clamped == 0 &&
          ri.status.smoothing == 0 && ri.status.t == 1.0f);

    int bitwise_all = 1;
    for (int i = 0; i < 64; i++) {
        float ang = -3.0f + 6.0f * (float)i / 63.0f;
        V2 t = { 1.3f * cosf(ang), 1.3f * sinf(ang) };
        IkChain cc;
        ik_chain_init(&cc, root, 1.5f, 0.8f, (V2){ 0.2f, -0.7f });
        IkResult p = ik_chain_solve(&cc, t, 1.0f / 60.0f);
        IkResult q = solve_two_bone(root, t, 1.5f, 0.8f, (V2){ 0.2f, -0.7f });
        if (p.ok != q.ok || p.joint.x != q.joint.x || p.joint.y != q.joint.y ||
            p.a1 != q.a1 || p.a2 != q.a2) { bitwise_all = 0; break; }
    }
    check("reachable_bitwise_sweep", bitwise_all);

    /* first clamp after a reachable pose: bounded step, not a jump */
    IkChain tr;
    ik_chain_init(&tr, root, 1.0f, 1.0f, pole);
    ik_chain_solve(&tr, (V2){ 1.5f, 0 }, 1.0f / 60.0f);
    float max_step = IK_DEFAULT_SMOOTH_SPEED / 60.0f;
    ik_chain_solve(&tr, (V2){ 5, 0 }, 1.0f / 60.0f);
    check("transition_bounded", vlen(vsub(tr.status.applied, (V2){ 1.5f, 0 })) <= max_step + 1e-5f);
    check("transition_in_flight", tr.status.smoothing == 1 && tr.status.t < 1.0f);
    check("transition_starts_near_zero", tr.status.t >= 0.0f && tr.status.t < 0.5f);
    check("transition_reachable_mid", reach_level(root, tr.status.applied, 1.0f, 1.0f) == IK_REACHABLE);
    float t_prev = tr.status.t;
    int t_increases = 1, t_advances = 0;
    int frames = 1;
    while (tr.status.smoothing && frames < 100) {
        ik_chain_solve(&tr, (V2){ 5, 0 }, 1.0f / 60.0f);
        if (tr.status.smoothing && tr.status.t < t_prev) t_increases = 0;
        if (tr.status.smoothing && tr.status.t > t_prev) t_advances++;
        t_prev = tr.status.t;
        frames++;
    }
    check("transition_t_monotonic", t_increases == 1 && t_advances > 0 && t_prev == 1.0f);
    check("transition_completes", frames < 10 && tr.status.t == 1.0f);
    check("transition_lands_on_ring", fabsf(vlen(vsub(tr.status.applied, root)) - 2.0f) < 1e-6f);
    check("transition_still_clamped", tr.status.clamped == 1 && tr.status.status == IK_UNREACHABLE);

    /* back inside: applied is the request again, no residual smoothing */
    ik_chain_solve(&tr, (V2){ 1.0f, 0.5f }, 1.0f / 60.0f);
    check("reenter_keeps_request",
          veq(tr.status.applied, (V2){ 1.0f, 0.5f }) && tr.status.clamped == 0 &&
          tr.status.smoothing == 0 && tr.status.t == 1.0f);

    /* long sweep past the outer ring: stable, smooth, on-annulus */
    IkChain sw;
    ik_chain_init(&sw, root, 1.5f, 0.8f, pole);
    float hi = 2.3f, lo = 0.7f;
    int missing_pose = 0, unreachable_frames = 0, jump_frames = 0, reversal_frames = 0;
    V2 prev_applied = { 0, 0 }, prev_desired = { 0, 0 };
    int have_prev = 0;
    V2 last_req = { 0, 0 };
    for (int i = 0; i < 180; i++) {
        float u = (float)i / 179.0f;
        float radius = 1.0f + 3.0f * u;
        float ang = 0.4f + 1.6f * u;
        V2 req = { radius * cosf(ang), radius * sinf(ang) };
        last_req = req;
        IkResult rr = ik_chain_solve(&sw, req, 1.0f / 60.0f);
        if (!rr.ok) missing_pose++;
        if (reach_level(root, sw.status.applied, 1.5f, 0.8f) != IK_REACHABLE) unreachable_frames++;
        V2 desired = (sw.status.status == IK_REACHABLE)
                         ? req
                         : clamp_to_ring(root, req, 1.5f, 0.8f);
        if (have_prev) {
            V2 da = vsub(sw.status.applied, prev_applied);
            V2 dd = vsub(desired, prev_desired);
            if (vlen(da) > max_step * (1.0f + 1e-3f)) jump_frames++;
            if (da.x * dd.x + da.y * dd.y < -1e-6f) reversal_frames++;
        }
        prev_applied = sw.status.applied;
        prev_desired = desired;
        have_prev = 1;
    }
    check("sweep_pose_every_frame", missing_pose == 0);
    check("sweep_zero_unreachable_frames", unreachable_frames == 0);
    check("sweep_zero_jump_frames", jump_frames == 0);
    check("sweep_zero_reversal_frames", reversal_frames == 0);
    check("sweep_settled", sw.status.smoothing == 0 && sw.status.t == 1.0f);
    V2 final_want = clamp_to_ring(root, last_req, 1.5f, 0.8f);
    check("sweep_final_is_clamp", fabsf(sw.status.applied.x - final_want.x) < 1e-5f &&
                                       fabsf(sw.status.applied.y - final_want.y) < 1e-5f);
    check("sweep_final_on_ring", fabsf(vlen(vsub(sw.status.applied, root)) - hi) < 1e-5f);

    /* static unreachable target: settles once, then holds still (no shaking) */
    IkChain so;
    ik_chain_init(&so, root, 1.5f, 0.8f, pole);
    int settle_frame = -1, moved_after = 0;
    V2 last = { 0, 0 };
    for (int i = 0; i < 300; i++) {
        ik_chain_solve(&so, (V2){ 5, 3 }, 1.0f / 60.0f);
        if (so.status.smoothing == 0 && settle_frame < 0) { settle_frame = i; last = so.status.applied; }
        else if (settle_frame >= 0) {
            if (vlen(vsub(so.status.applied, last)) > 1e-7f) moved_after++;
            last = so.status.applied;
        }
    }
    check("soak_settles", settle_frame >= 0 && settle_frame < 60);
    check("soak_no_shake", moved_after == 0);
    check("soak_stays_on_ring", reach_level(root, so.status.applied, 1.5f, 0.8f) == IK_REACHABLE &&
                                    fabsf(vlen(vsub(so.status.applied, root)) - hi) < 1e-5f);

    /* a zero dt still yields a valid, reachable pose */
    IkChain zd;
    ik_chain_init(&zd, root, 1.0f, 1.0f, pole);
    ik_chain_solve(&zd, (V2){ 1.0f, 0 }, 1.0f / 60.0f);
    IkResult zdr = ik_chain_solve(&zd, (V2){ 5, 0 }, 0.0f);
    check("zero_dt_valid", zdr.ok == 1 && reach_level(root, zd.status.applied, 1.0f, 1.0f) == IK_REACHABLE);
    check("zero_dt_smooths", zd.status.smoothing == 1);

    /* outer ring to the far side of the inner ring: a straight lerp would
       cut through the hole, the smoothed path has to stay on the annulus */
    IkChain cx;
    ik_chain_init(&cx, root, 1.5f, 0.8f, pole);
    ik_chain_solve(&cx, (V2){ hi, 0 }, 1.0f / 60.0f);      /* seed on outer ring */
    ik_chain_solve(&cx, (V2){ -0.3f, 0 }, 1.0f / 60.0f);   /* inner hole, opposite */
    int cross_bad = 0, cross_ok_frames = 0;
    while (cx.status.smoothing && cross_ok_frames < 200) {
        if (!cx.status.ok) cross_bad++;
        if (reach_level(root, cx.status.applied, 1.5f, 0.8f) != IK_REACHABLE) cross_bad++;
        ik_chain_solve(&cx, (V2){ -0.3f, 0 }, 1.0f / 60.0f);
        cross_ok_frames++;
    }
    check("cross_transition_zero_bad_frames", cross_bad == 0);
    check("cross_transition_length", cross_ok_frames > 5 && cross_ok_frames < 200);
    check("cross_lands_inner_ring", fabsf(vlen(vsub(cx.status.applied, root)) - lo) < 1e-5f &&
                                        cx.status.applied.x < 0.0f);
    check("cross_status_too_close", cx.status.status == IK_TOO_CLOSE && cx.status.clamped == 1);

    /* degenerate limb wins over an unreachable request */
    IkChain zl;
    ik_chain_init(&zl, root, 1.0f, 0.0f, pole);
    IkResult zlr = ik_chain_solve(&zl, (V2){ 9, 9 }, 1.0f / 60.0f);
    check("err_beats_clamp", zlr.ok == 0 && zl.status.err == IK_ERR_DEGENERATE_LIMB &&
                                 zl.status.ok == 0 && zl.status.clamped == 0);

    /* zero target wins over both bones being zero */
    IkChain zz;
    ik_chain_init(&zz, root, 0.0f, 0.0f, pole);
    IkResult zzr = ik_chain_solve(&zz, root, 1.0f / 60.0f);
    check("err_zero_target_and_bones",
          zzr.ok == 0 && zz.status.err == IK_ERR_DEGENERATE_LIMB && zz.status.ok == 0);
    ik_chain_init(&zz, root, 1.5f, 0.8f, pole);
    zzr = ik_chain_solve(&zz, (V2){ 1e-9f, 0 }, 1.0f / 60.0f);
    check("err_near_zero_target", zzr.ok == 0 && zz.status.err == IK_ERR_DEGENERATE_TARGET &&
                                      zz.status.status == IK_TOO_CLOSE);
    ik_chain_init(&zz, root, 1.0f, 1.0f, pole);
    zzr = ik_chain_solve(&zz, (V2){ 1e-9f, 0 }, 1.0f / 60.0f);
    check("err_zero_target_on_equal_arms", zzr.ok == 0 &&
                                               zz.status.err == IK_ERR_DEGENERATE_TARGET &&
                                               zz.status.status == IK_REACHABLE);

    /* query API answers without the caller reaching into internals */
    IkChain qc;
    ik_chain_init(&qc, root, 1.0f, 1.0f, pole);
    V2 qout = { 99, 99 };
    check("query_before_solve", ik_chain_target(&qc, &qout) == 0);
    check("query_empty_status", ik_chain_status(NULL).ok == 0 && ik_chain_status(NULL).t == 1.0f);
    ik_chain_solve(&qc, (V2){ 1.2f, 0 }, 1.0f / 60.0f);
    int qgot = ik_chain_target(&qc, &qout);
    check("query_after_solve", qgot == 1 && qout.x == 1.2f && qout.y == 0.0f);
    check("query_status_snapshot", ik_chain_status(&qc).status == IK_REACHABLE &&
                                       ik_chain_status(&qc).applied.x == 1.2f);
    ik_chain_solve(&qc, (V2){ 6, 0 }, 1.0f / 60.0f);
    ik_chain_set_smooth_speed(&qc, 1.0f);
    ik_chain_solve(&qc, (V2){ 6, 0 }, 1.0f / 60.0f);
    check("query_speed_setter", qc.smooth.speed == 1.0f && qc.status.smoothing == 1);
    ik_chain_set_smooth_speed(&qc, -1.0f);
    check("query_speed_rejects_bad", qc.smooth.speed == 1.0f);
    ik_chain_set_smooth_speed(&qc, 0.0f);
    check("query_speed_rejects_zero", qc.smooth.speed == 1.0f);

    /* inner hole: clamped to the inner ring and held there */
    IkChain ih;
    ik_chain_init(&ih, root, 1.5f, 0.8f, pole);
    ik_chain_solve(&ih, (V2){ 2.0f, 0 }, 1.0f / 60.0f);   /* reachable seed */
    for (int i = 0; i < 120; i++) ik_chain_solve(&ih, (V2){ 0.2f, 0 }, 1.0f / 60.0f);
    check("inner_clamped", ih.status.status == IK_TOO_CLOSE && fabsf(vlen(vsub(ih.status.applied, root)) - lo) < 1e-5f);
    check("inner_reachable", reach_level(root, ih.status.applied, 1.5f, 0.8f) == IK_REACHABLE);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
