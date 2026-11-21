
/* constraint: verlet particle chain with distance constraints, iteration
   count and damping; used as the toy model for cloth/rope style solvers.
   Also carries conical joint limits: the child segment of an interior
   particle is constrained to a half-angle around the parent segment. */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define MAX_P 16

typedef struct { float x, y; } V2;

typedef struct { float x, y, px, py; float inv_mass; int pinned; } Particle;

typedef struct {
    Particle p[MAX_P];
    int n;
    float rest[MAX_P];      /* rest length between i and i+1 */
    int iterations;
    float damping;
    float cone[MAX_P];      /* half-angle at particle i, i in [1,n-2]; <0 disables */
} Chain;

/* status codes for the cone projection */
enum {
    CONE_OK = 0,
    CONE_ERR_NEG_ANGLE = 1,   /* negative cone half-angle */
    CONE_ERR_ZERO_AXIS = 2,   /* degenerate cone axis */
    CONE_ERR_ZERO_DIR = 3,    /* degenerate direction */
    CONE_ERR_BAD_JOINT = 4    /* joint index outside [1, n-2] */
};

static float v2_len(V2 v) { return sqrtf(v.x * v.x + v.y * v.y); }
static float v2_dot(V2 a, V2 b) { return a.x * b.x + a.y * b.y; }
static float v2_cross(V2 a, V2 b) { return a.x * b.y - a.y * b.x; }

static float wrap_pi(float a) {
    const float pi = 3.14159265358979323846f;
    const float tau = 6.28318530717958647692f;
    while (a > pi) a -= tau;
    while (a < -pi) a += tau;
    return a;
}

static void chain_init(Chain* c, int n, float seg, float x0, float y0) {
    memset(c, 0, sizeof(*c));
    c->n = n;
    c->iterations = 8;
    c->damping = 0.99f;
    for (int i = 0; i < MAX_P; i++) c->cone[i] = -1.0f;
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

static void chain_solve_distances(Chain* c) {
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
    chain_solve_distances(c);
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

/* Project dir into the cone of half-angle max_angle about axis.

   The direction length is preserved; only its angle changes.  A negative
   half-angle or a zero-length axis/direction is a hard error, not a silent
   no-op.  removed (optional) receives the angle that had to be removed. */
static int cone_project(V2* dir, V2 axis, float max_angle, float* removed) {
    if (removed) *removed = 0.0f;
    if (max_angle < 0.0f) return CONE_ERR_NEG_ANGLE;
    float al = v2_len(axis);
    if (al < 1e-6f) return CONE_ERR_ZERO_AXIS;
    float dl = v2_len(*dir);
    if (dl < 1e-6f) return CONE_ERR_ZERO_DIR;

    V2 an = { axis.x / al, axis.y / al };
    float rel = atan2f(v2_cross(an, *dir), v2_dot(an, *dir));  /* [-pi, pi] */
    float excess = fabsf(rel) - max_angle;
    if (excess <= 0.0f) return CONE_OK;

    float base = atan2f(an.y, an.x);
    float ang = base + (rel < 0.0f ? -max_angle : max_angle);
    dir->x = dl * cosf(ang);
    dir->y = dl * sinf(ang);
    if (removed) *removed = excess;
    return CONE_OK;
}

/* Rigidly rotate particles pivot+1..n-1 about the pivot, including their
   previous positions so the Verlet velocity is rotated consistently. */
static void chain_rotate_from(Chain* c, int pivot, float delta) {
    Particle* o = &c->p[pivot];
    float cs = cosf(delta), sn = sinf(delta);
    for (int j = pivot + 1; j < c->n; j++) {
        Particle* q = &c->p[j];
        float dx = q->x - o->x, dy = q->y - o->y;
        q->x = o->x + dx * cs - dy * sn;
        q->y = o->y + dx * sn + dy * cs;
        dx = q->px - o->px; dy = q->py - o->py;
        q->px = o->px + dx * cs - dy * sn;
        q->py = o->py + dx * sn + dy * cs;
    }
}

/* Enable a cone limit at interior particle joint (1..n-2). */
static int chain_set_cone(Chain* c, int joint, float max_angle) {
    if (joint < 1 || joint + 1 >= c->n) return CONE_ERR_BAD_JOINT;
    if (max_angle < 0.0f) return CONE_ERR_NEG_ANGLE;
    c->cone[joint] = max_angle;
    return CONE_OK;
}

/* Enforce one joint's cone.  The cone axis is the parent segment direction;
   the child sub-chain is rotated rigidly about the joint, so every segment
   length is preserved.  out_drift is the largest displacement of a moved
   particle, bounded by removed_angle * max_radius. */
static int chain_cone_limit(Chain* c, int joint, float* out_removed, float* out_drift) {
    if (out_removed) *out_removed = 0.0f;
    if (out_drift) *out_drift = 0.0f;
    if (joint < 1 || joint + 1 >= c->n) return CONE_ERR_BAD_JOINT;
    float max_angle = c->cone[joint];
    if (max_angle < 0.0f) return CONE_OK;            /* no limit configured */

    Particle* a = &c->p[joint - 1];
    Particle* b = &c->p[joint];
    Particle* d = &c->p[joint + 1];
    V2 axis = { b->x - a->x, b->y - a->y };
    V2 old = { d->x - b->x, d->y - b->y };
    V2 dir = old;

    float removed = 0.0f;
    int rc = cone_project(&dir, axis, max_angle, &removed);
    if (rc != CONE_OK) return rc;
    if (removed <= 0.0f) return CONE_OK;

    float delta = wrap_pi(atan2f(dir.y, dir.x) - atan2f(old.y, old.x));
    float radius = 0.0f;
    for (int j = joint + 1; j < c->n; j++) {
        float dx = c->p[j].x - b->x, dy = c->p[j].y - b->y;
        float r = sqrtf(dx * dx + dy * dy);
        if (r > radius) radius = r;
    }
    chain_rotate_from(c, joint, delta);

    if (out_removed) *out_removed = removed;
    if (out_drift) *out_drift = 2.0f * radius * fabsf(sinf(0.5f * delta));
    return CONE_OK;
}

/* Remaining angle outside the cone at joint, zero when inside. */
static float chain_cone_excess(const Chain* c, int joint) {
    if (joint < 1 || joint + 1 >= c->n) return 0.0f;
    if (c->cone[joint] < 0.0f) return 0.0f;
    const Particle* a = &c->p[joint - 1];
    const Particle* b = &c->p[joint];
    const Particle* d = &c->p[joint + 1];
    V2 axis = { b->x - a->x, b->y - a->y };
    V2 dir = { d->x - b->x, d->y - b->y };
    float al = v2_len(axis), dl = v2_len(dir);
    if (al < 1e-6f || dl < 1e-6f) return 0.0f;
    float rel = atan2f(v2_cross(axis, dir), v2_dot(axis, dir));
    float excess = fabsf(rel) - c->cone[joint];
    return excess > 0.0f ? excess : 0.0f;
}

/* Reference implementation for tests: clamp the child's world angle to the
   parent's world angle +- max_angle.  Independent of cone_project's
   vector path, but must produce the same pose. */
static int chain_cone_limit_manual(Chain* c, int joint, float* out_removed) {
    if (out_removed) *out_removed = 0.0f;
    if (joint < 1 || joint + 1 >= c->n) return CONE_ERR_BAD_JOINT;
    float max_angle = c->cone[joint];
    if (max_angle < 0.0f) return CONE_OK;

    const Particle* a = &c->p[joint - 1];
    const Particle* b = &c->p[joint];
    const Particle* d = &c->p[joint + 1];
    V2 axis = { b->x - a->x, b->y - a->y };
    V2 dir = { d->x - b->x, d->y - b->y };
    if (v2_len(axis) < 1e-6f) return CONE_ERR_ZERO_AXIS;
    if (v2_len(dir) < 1e-6f) return CONE_ERR_ZERO_DIR;

    float rel = wrap_pi(atan2f(dir.y, dir.x) - atan2f(axis.y, axis.x));
    float clamped = rel;
    if (clamped > max_angle) clamped = max_angle;
    if (clamped < -max_angle) clamped = -max_angle;
    float delta = clamped - rel;
    if (fabsf(delta) < 1e-9f) return CONE_OK;

    chain_rotate_from(c, joint, delta);
    if (out_removed) *out_removed = fabsf(delta);
    return CONE_OK;
}

static int chain_project_cone_limits(Chain* c, float* out_removed, float* out_drift) {
    float removed_total = 0.0f, drift_total = 0.0f;
    int rc = CONE_OK;
    for (int j = 1; j + 1 < c->n; j++) {
        if (c->cone[j] < 0.0f) continue;
        float removed = 0.0f, drift = 0.0f;
        rc = chain_cone_limit(c, j, &removed, &drift);
        if (rc != CONE_OK) break;
        removed_total += removed;
        drift_total += drift;
    }
    if (out_removed) *out_removed = removed_total;
    if (out_drift) *out_drift = drift_total;
    return rc;
}

/* Required order: project the pose into the cones, then re-solve the chain.
   Projection is a rigid rotation of the child sub-chain about the joint, so
   segment lengths are unchanged and the distance pass does not reintroduce
   any cone excess. */
static int chain_apply_cone_limits(Chain* c, float* out_removed, float* out_drift) {
    chain_solve_distances(c);                       /* consistent lengths first */
    int rc = chain_project_cone_limits(c, out_removed, out_drift);
    if (rc != CONE_OK) return rc;
    chain_solve_distances(c);                       /* no-op on a valid pose */
    return CONE_OK;
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

    /* --- cone projection primitive --- */
    {
        V2 axis = { 1.0f, 0.0f };
        V2 out = { cosf(1.2f), sinf(1.2f) };
        float removed = 0.0f;
        int rc = cone_project(&out, axis, 0.5f, &removed);
        check("cone_project_ok", rc == CONE_OK);
        check("cone_project_angle", fabsf(atan2f(out.y, out.x) - 0.5f) < 1e-5f);
        check("cone_project_removed", fabsf(removed - 0.7f) < 1e-5f);
        check("cone_project_keeps_len", fabsf(v2_len(out) - 1.0f) < 1e-5f);

        V2 in = { cosf(0.3f), sinf(0.3f) };
        removed = -1.0f;
        rc = cone_project(&in, axis, 0.5f, &removed);
        check("cone_inside_untouched",
              rc == CONE_OK && removed == 0.0f &&
              fabsf(atan2f(in.y, in.x) - 0.3f) < 1e-6f);

        V2 below = { cosf(-2.0f), sinf(-2.0f) };
        removed = 0.0f;
        rc = cone_project(&below, axis, 0.4f, &removed);
        check("cone_project_negative_side",
              rc == CONE_OK &&
              fabsf(atan2f(below.y, below.x) + 0.4f) < 1e-5f &&
              fabsf(removed - 1.6f) < 1e-5f);

        V2 tmp = { 1.0f, 0.0f };
        check("cone_negative_angle_error",
              cone_project(&tmp, axis, -0.1f, NULL) == CONE_ERR_NEG_ANGLE);
        V2 zero_axis = { 0.0f, 0.0f };
        check("cone_zero_axis_error",
              cone_project(&tmp, zero_axis, 0.5f, NULL) == CONE_ERR_ZERO_AXIS);
        V2 zero_dir = { 0.0f, 0.0f };
        check("cone_zero_dir_error",
              cone_project(&zero_dir, axis, 0.5f, NULL) == CONE_ERR_ZERO_DIR);
    }

    /* --- cone limits on a folded ragdoll chain --- */
    Chain lim;
    chain_init(&lim, 4, 1.0f, 0.0f, 0.0f);
    lim.p[1].x = 1.0f; lim.p[1].y = 0.0f;
    lim.p[2].x = 1.0f + cosf(-2.0f); lim.p[2].y = sinf(-2.0f);
    lim.p[3].x = lim.p[2].x + cosf(-0.4f); lim.p[3].y = lim.p[2].y + sinf(-0.4f);
    for (int i = 0; i < lim.n; i++) { lim.p[i].px = lim.p[i].x; lim.p[i].py = lim.p[i].y; }

    check("cone_setup_lengths", chain_seg_error(&lim) < 1e-5f);
    check("cone_set", chain_set_cone(&lim, 1, 0.5f) == CONE_OK &&
                      chain_set_cone(&lim, 2, 0.4f) == CONE_OK);
    check("cone_violated_before", chain_cone_excess(&lim, 1) > 1.4f &&
                                   chain_cone_excess(&lim, 2) > 1.1f);

    float radius1 = 0.0f;
    for (int i = 2; i < lim.n; i++) {
        float dx = lim.p[i].x - lim.p[1].x, dy = lim.p[i].y - lim.p[1].y;
        float r = sqrtf(dx * dx + dy * dy);
        if (r > radius1) radius1 = r;
    }
    Chain auto_c = lim;
    float rem1 = 0.0f, drift1 = 0.0f;
    int rc1 = chain_cone_limit(&auto_c, 1, &rem1, &drift1);
    check("cone_project_j1", rc1 == CONE_OK && fabsf(rem1 - 1.5f) < 1e-4f);
    check("cone_drift_bound_j1", drift1 <= rem1 * radius1 + 1e-5f);

    float radius2 = 0.0f;
    for (int i = 3; i < auto_c.n; i++) {
        float dx = auto_c.p[i].x - auto_c.p[2].x, dy = auto_c.p[i].y - auto_c.p[2].y;
        float r = sqrtf(dx * dx + dy * dy);
        if (r > radius2) radius2 = r;
    }
    float rem2 = 0.0f, drift2 = 0.0f;
    int rc2 = chain_cone_limit(&auto_c, 2, &rem2, &drift2);
    check("cone_project_j2", rc2 == CONE_OK && fabsf(rem2 - 1.2f) < 1e-4f);
    check("cone_drift_bound_j2", drift2 <= rem2 * radius2 + 1e-5f);

    check("cone_excess_zero_after", chain_cone_excess(&auto_c, 1) < 1e-5f &&
                                     chain_cone_excess(&auto_c, 2) < 1e-5f);
    check("cone_lengths_preserved", chain_seg_error(&auto_c) < 1e-5f);

    /* manual clamp must land on the same pose */
    Chain man = lim;
    float mrem = 0.0f;
    int mrc1 = chain_cone_limit_manual(&man, 1, &mrem);
    int mrc2 = chain_cone_limit_manual(&man, 2, &mrem);
    check("cone_manual_ok", mrc1 == CONE_OK && mrc2 == CONE_OK);
    float pose_diff = 0.0f;
    for (int i = 0; i < man.n; i++) {
        float dx = man.p[i].x - auto_c.p[i].x;
        float dy = man.p[i].y - auto_c.p[i].y;
        float d = sqrtf(dx * dx + dy * dy);
        if (d > pose_diff) pose_diff = d;
    }
    check("cone_matches_manual", pose_diff < 1e-5f);

    /* project first, then re-solve the distance chain */
    Chain solved = lim;
    float rem_total = 0.0f, drift_total = 0.0f;
    int src = chain_apply_cone_limits(&solved, &rem_total, &drift_total);
    check("cone_then_distances_ok", src == CONE_OK);
    check("cone_post_chain_excess_zero",
          chain_cone_excess(&solved, 1) < 1e-5f &&
          chain_cone_excess(&solved, 2) < 1e-5f);
    check("cone_post_chain_lengths", chain_seg_error(&solved) < 1e-5f);
    check("cone_total_removed", fabsf(rem_total - (rem1 + rem2)) < 1e-4f);

    float net = 0.0f;
    for (int i = 0; i < solved.n; i++) {
        float dx = solved.p[i].x - lim.p[i].x;
        float dy = solved.p[i].y - lim.p[i].y;
        float d = sqrtf(dx * dx + dy * dy);
        if (d > net) net = d;
    }
    check("cone_drift_net_bounded", net <= drift_total + 1e-4f);
    check("cone_drift_analytic_bound",
          drift_total <= rem1 * radius1 + rem2 * radius2 + 1e-4f);

    check("cone_set_negative_rejected",
          chain_set_cone(&lim, 1, -0.2f) == CONE_ERR_NEG_ANGLE);
    check("cone_set_bad_joint_rejected",
          chain_set_cone(&lim, 0, 0.5f) == CONE_ERR_BAD_JOINT);

    Chain degen;
    chain_init(&degen, 3, 1.0f, 0.0f, 0.0f);
    degen.p[1].x = degen.p[0].x; degen.p[1].y = degen.p[0].y;
    degen.p[1].px = degen.p[1].x; degen.p[1].py = degen.p[1].y;
    chain_set_cone(&degen, 1, 0.5f);
    float rr = 0.0f, dd = 0.0f;
    check("cone_zero_axis_chain_error",
          chain_cone_limit(&degen, 1, &rr, &dd) == CONE_ERR_ZERO_AXIS);

    Chain zdir;
    chain_init(&zdir, 3, 1.0f, 0.0f, 0.0f);
    zdir.p[2].x = zdir.p[1].x; zdir.p[2].y = zdir.p[1].y;
    zdir.p[2].px = zdir.p[2].x; zdir.p[2].py = zdir.p[2].y;
    chain_set_cone(&zdir, 1, 0.5f);
    check("cone_zero_dir_chain_error",
          chain_cone_limit(&zdir, 1, &rr, &dd) == CONE_ERR_ZERO_DIR);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
