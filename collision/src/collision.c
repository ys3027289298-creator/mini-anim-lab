
/* collision: swept (continuous) collision for the verlet particle chain.

   A naive solver only projects the *current* position out of an obstacle.
   At high swing speed a particle can be on the far side of a circle after
   one frame: the projection pushes it back through, so visually the rope
   passes through and springs back (tunnelling).

   The fix here:
     1. continuous detection against the segment from the previous-frame
        position to the current position (seg_circle_hit, entry time only);
     2. on a hit, roll the particle back to the nearest surface point at the
        first contact time and move its verlet anchor (px,py) with it, so the
        rebuilt velocity never gains a spurious component;
     3. cap the residual penetration of any single step (pen_cap);
     4. zero restitution: the inward normal velocity is removed but never
        reversed, so a rollback cannot produce a rebound jump.
*/
#include <stdio.h>
#include <math.h>
#include <string.h>

#define MAX_P 16
#define EPS 1e-7f

typedef struct { float x, y; } V2;

static V2 v2(float x, float y) { return (V2){ x, y }; }
static V2 vadd(V2 a, V2 b) { return v2(a.x + b.x, a.y + b.y); }
static V2 vsub(V2 a, V2 b) { return v2(a.x - b.x, a.y - b.y); }
static V2 vmul(V2 a, float s) { return v2(a.x * s, a.y * s); }
static float vdot(V2 a, V2 b) { return a.x * b.x + a.y * b.y; }
static float vlen(V2 a) { return sqrtf(vdot(a, a)); }

typedef struct { float x, y, px, py, inv_mass; int pinned; } Particle;

typedef struct {
    V2 c;
    float r;
} Circle;

typedef struct {
    Particle p[MAX_P];
    int n;
    float rest[MAX_P];      /* rest length between i and i+1 */
    int iterations;
    float damping;
    Circle obs;
    int has_obs;
    float pen_cap;          /* upper bound on residual per-step penetration */
    /* diagnostics from the last rope_step() call */
    float worst_pen;        /* largest residual penetration after resolution */
    float worst_raw_pen;    /* deepest sweep penetration before resolution */
    int contacts;           /* particle contacts resolved this step */
    int reversed;           /* count of post-fix normal-velocity reversals */
    int cap_clamps;         /* times raw penetration exceeded the cap */
} Rope;

static void rope_init(Rope* rp, int n, float seg, float x0, float y0) {
    memset(rp, 0, sizeof(*rp));
    rp->n = n;
    rp->iterations = 8;
    rp->damping = 0.99f;
    rp->pen_cap = 1e-3f;
    for (int i = 0; i < n; i++) {
        rp->p[i].x = x0 + seg * (float)i;
        rp->p[i].y = y0;
        rp->p[i].px = rp->p[i].x;
        rp->p[i].py = rp->p[i].y;
        rp->p[i].inv_mass = 1.0f;
        rp->p[i].pinned = 0;
        rp->rest[i] = seg;
    }
    rp->p[0].pinned = 1;
    rp->p[0].inv_mass = 0.0f;
}

static void rope_set_obstacle(Rope* rp, V2 c, float r) {
    rp->obs.c = c;
    rp->obs.r = r;
    rp->has_obs = 1;
}

/* Segment a->b swept against circle o. Returns 1 on an *entry* hit with
   t in [0,1], contact point q and outward unit normal n. Entry only: a
   segment starting inside the circle reports no hit, its penetration is
   handled separately by the positional pass (deep penetration fallback). */
static int seg_circle_hit(V2 a, V2 b, Circle o, float* t_out,
                          V2* q_out, V2* n_out) {
    V2 d = vsub(b, a);
    V2 f = vsub(a, o.c);
    float aa = vdot(d, d);
    if (aa < 1e-12f) {
        /* zero-length sweep: no entry time exists, penetration (if any) is
           handled by the caller's deep-penetration fallback */
        return 0;
    }
    float bb = 2.0f * vdot(f, d);
    float cc = vdot(f, f) - o.r * o.r;
    if (cc <= 0.0f) return 0;                    /* starts on/inside: no entry */
    float disc = bb * bb - 4.0f * aa * cc;
    if (disc < 0.0f) return 0;
    float t = (-bb - sqrtf(disc)) / (2.0f * aa);
    if (t < 0.0f || t > 1.0f) return 0;
    V2 q = vadd(a, vmul(d, t));
    V2 n = vsub(q, o.c);
    float nl = vlen(n);
    if (nl < EPS) return 0;
    n = vmul(n, 1.0f / nl);
    if (t_out) *t_out = t;
    if (q_out) *q_out = q;
    if (n_out) *n_out = n;
    return 1;
}

static float point_pen(V2 x, Circle o) {
    float d = vlen(vsub(x, o.c));
    float pen = o.r - d;
    return pen > 0.0f ? pen : 0.0f;
}

/* Deepest penetration reached anywhere along the swept segment a->b:
   the distance from the circle centre to the closest point on the segment.
   Unlike point_pen(b) this also catches a tunnel that ends outside. */
static float sweep_pen(V2 a, V2 b, Circle o) {
    V2 d = vsub(b, a);
    float dd = vdot(d, d);
    float t = 0.0f;
    if (dd > 1e-12f) {
        t = vdot(vsub(o.c, a), d) / dd;
        if (t < 0.0f) t = 0.0f;
        if (t > 1.0f) t = 1.0f;
    }
    V2 closest = vadd(a, vmul(d, t));
    float pen = o.r - vlen(vsub(o.c, closest));
    return pen > 0.0f ? pen : 0.0f;
}

static void solve_distance(Rope* rp) {
    for (int it = 0; it < rp->iterations; it++) {
        for (int i = 0; i < rp->n - 1; i++) {
            Particle* a = &rp->p[i];
            Particle* b = &rp->p[i + 1];
            float dx = b->x - a->x, dy = b->y - a->y;
            float d = sqrtf(dx * dx + dy * dy);
            if (d < 1e-6f) continue;
            float diff = (d - rp->rest[i]) / d;
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

/* Positional collision pass for one particle. frame is its position at the
   start of the step, so the swept segment (frame -> now) is stable across
   constraint iterations. Sets touched when this particle made contact this
   pass, incoming when that contact was made while moving inward, and stores
   the contact-time outward normal in contact_n (which must not be
   recomputed from the final position: the particle can wrap past the
   centre and invert the geometric normal). */
static void resolve_particle(Rope* rp, int i, V2 frame,
                             int* touched, int* incoming, V2* contact_n) {
    Particle* pt = &rp->p[i];
    if (pt->pinned) return;
    V2 now = v2(pt->x, pt->y);
    /* raw = deepest penetration the unconstrained step would have reached */
    float raw = sweep_pen(frame, now, rp->obs);
    if (raw > rp->worst_raw_pen) rp->worst_raw_pen = raw;

    V2 q, n;
    float t;
    if (seg_circle_hit(frame, now, rp->obs, &t, &q, &n)) {
        /* roll back to the nearest surface point at first contact; the
           verlet anchor travels by the same delta so its rebuilt velocity
           keeps only the tangential component. */
        (void)t;
        V2 target = q;
        pt->px += target.x - now.x;
        pt->py += target.y - now.y;
        pt->x = target.x;
        pt->y = target.y;
        *touched = 1;
        *incoming = 1;                       /* an entry hit is inward by construction */
        *contact_n = n;
        return;
    }

    /* no entry hit this segment: handle residual / deep penetration */
    if (point_pen(now, rp->obs) > 0.0f) {
        V2 out = vsub(now, rp->obs.c);
        float d = vlen(out);
        if (d < EPS) return;
        n = vmul(out, 1.0f / d);
        /* never allow more than pen_cap of penetration in one step */
        float pen = point_pen(now, rp->obs);
        float keep = pen > rp->pen_cap ? rp->pen_cap : pen;
        if (pen > rp->pen_cap) rp->cap_clamps++;
        V2 target = vadd(rp->obs.c, vmul(n, rp->obs.r - keep));
        V2 delta = vsub(target, now);
        /* anchor follows the push-out: rebuilt velocity loses its outward
           normal component instead of turning into a rebound jump */
        pt->px += delta.x;
        pt->py += delta.y;
        pt->x = target.x;
        pt->y = target.y;
        *touched = 1;
        *contact_n = n;
        if (vdot(vsub(now, frame), n) < -EPS) *incoming = 1;
    }
}

/* Zero-restitution velocity fixup against the contact-time normal. Remove
   the inward normal component of the verlet velocity for every particle
   that made contact this step; the normal component is clamped at zero,
   never negated. A reversal jump is only meaningful for a particle that
   was actually moving inward at contact. */
static void fix_contact_velocities(Rope* rp, const int* touched,
                                   const int* incoming,
                                   const V2* contact_n) {
    for (int i = 0; i < rp->n; i++) {
        if (!touched[i] || rp->p[i].pinned) continue;
        Particle* pt = &rp->p[i];
        V2 v = v2(pt->x - pt->px, pt->y - pt->py);
        V2 n = contact_n[i];
        if (vlen(n) < EPS) continue;
        float vn = vdot(v, n);
        if (vn < -EPS) {
            /* remove the inward component: shift the anchor so the rebuilt
               velocity's normal part is exactly zero, never reversed */
            pt->px += vn * n.x;
            pt->py += vn * n.y;
        } else if (vn > EPS && incoming[i]) {
            /* inward-moving contact ended with outward velocity: rebound */
            rp->reversed++;
        }
    }
}

/* One verlet step. Returns -1 when dt is not strictly positive. */
static int rope_step(Rope* rp, float dt, float gx, float gy) {
    if (!(dt > 0.0f)) return -1;

    V2 frame[MAX_P];
    for (int i = 0; i < rp->n; i++) {
        Particle* pt = &rp->p[i];
        frame[i] = v2(pt->x, pt->y);
        if (pt->pinned) continue;
        float vx = (pt->x - pt->px) * rp->damping;
        float vy = (pt->y - pt->py) * rp->damping;
        pt->px = pt->x;
        pt->py = pt->y;
        pt->x += vx + gx * dt * dt;
        pt->y += vy + gy * dt * dt;
    }

    rp->worst_pen = 0.0f;
    rp->worst_raw_pen = 0.0f;
    rp->contacts = 0;
    rp->reversed = 0;
    rp->cap_clamps = 0;

    int touched[MAX_P] = { 0 };
    int incoming[MAX_P] = { 0 };
    V2 contact_n[MAX_P];
    memset(contact_n, 0, sizeof(contact_n));
    if (rp->has_obs) {
        for (int it = 0; it < rp->iterations; it++) {
            solve_distance(rp);
            for (int i = 0; i < rp->n; i++)
                resolve_particle(rp, i, frame[i], &touched[i],
                                 &incoming[i], &contact_n[i]);
        }
        fix_contact_velocities(rp, touched, incoming, contact_n);
        for (int i = 0; i < rp->n; i++) {
            if (touched[i]) rp->contacts++;
            float pen = point_pen(v2(rp->p[i].x, rp->p[i].y), rp->obs);
            if (pen > rp->worst_pen) rp->worst_pen = pen;
        }
    } else {
        solve_distance(rp);
    }
    return 0;
}

/* Naive baseline used to reproduce the bug: integrate and project the
   current position out of the circle only. A fast particle ending up on
   the far side is projected back through the obstacle. */
static int rope_step_naive(Rope* rp, float dt, float gx, float gy) {
    if (!(dt > 0.0f)) return -1;
    for (int i = 0; i < rp->n; i++) {
        Particle* pt = &rp->p[i];
        if (pt->pinned) continue;
        float vx = (pt->x - pt->px) * rp->damping;
        float vy = (pt->y - pt->py) * rp->damping;
        pt->px = pt->x;
        pt->py = pt->y;
        pt->x += vx + gx * dt * dt;
        pt->y += vy + gy * dt * dt;
    }
    for (int it = 0; it < rp->iterations; it++) {
        solve_distance(rp);
        if (!rp->has_obs) continue;
        for (int i = 0; i < rp->n; i++) {
            Particle* pt = &rp->p[i];
            if (pt->pinned) continue;
            V2 out = vsub(v2(pt->x, pt->y), rp->obs.c);
            float d = vlen(out);
            if (d >= rp->obs.r || d < EPS) continue;
            V2 n = vmul(out, 1.0f / d);
            V2 target = vadd(rp->obs.c, vmul(n, rp->obs.r));
            pt->x = target.x;
            pt->y = target.y;
        }
    }
    rp->worst_pen = 0.0f;
    for (int i = 0; i < rp->n; i++) {
        float pen = point_pen(v2(rp->p[i].x, rp->p[i].y), rp->obs);
        if (pen > rp->worst_pen) rp->worst_pen = pen;
    }
    return 0;
}

static float chain_seg_error(const Rope* rp) {
    float worst = 0.0f;
    for (int i = 0; i < rp->n - 1; i++) {
        float dx = rp->p[i + 1].x - rp->p[i].x;
        float dy = rp->p[i + 1].y - rp->p[i].y;
        float d = sqrtf(dx * dx + dy * dy);
        float e = fabsf(d - rp->rest[i]);
        if (e > worst) worst = e;
    }
    return worst;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

/* one free particle (anchor pin unused here) launched straight through
   a circle; the step length dwarfs the diameter, a guaranteed tunnel. */
static Rope make_cannon(float speed, float dt) {
    Rope rp;
    rope_init(&rp, 1, 1.0f, -4.0f, 0.0f);
    rp.p[0].pinned = 0;
    rp.p[0].inv_mass = 1.0f;
    rp.p[0].px = rp.p[0].x - speed * dt;   /* physical velocity via anchor */
    rp.p[0].py = rp.p[0].y;
    rp.damping = 1.0f;
    rope_set_obstacle(&rp, v2(0.0f, 0.0f), 1.0f);
    return rp;
}

int main(void) {
    /* ---- unit: swept segment / circle entry hit ---- */
    Circle unit = { v2(0.0f, 0.0f), 1.0f };
    V2 q, n;
    float t;
    check("ccd_hit", seg_circle_hit(v2(-4, 0), v2(4, 0), unit, &t, &q, &n)
          && fabsf(t - 3.0f / 8.0f) < 1e-5f
          && fabsf(q.x + 1.0f) < 1e-5f && fabsf(q.y) < 1e-5f
          && fabsf(n.x + 1.0f) < 1e-5f);
    check("ccd_miss", !seg_circle_hit(v2(-4, 3), v2(4, 3), unit, &t, &q, &n));
    check("ccd_inside_no_entry",
          !seg_circle_hit(v2(0, 0), v2(4, 0), unit, &t, &q, &n));
    check("ccd_glancing", seg_circle_hit(v2(-4, 0.5f), v2(4, 0.5f),
                                         unit, &t, &q, &n)
          && fabsf(vlen(vsub(q, unit.c)) - 1.0f) < 1e-5f);

    /* ---- zero / negative time steps must be rejected ---- */
    Rope rz;
    rope_init(&rz, 4, 1.0f, 0.0f, 0.0f);
    check("dt_zero_error", rope_step(&rz, 0.0f, 0.0f, -9.8f) == -1);
    check("dt_neg_error", rope_step(&rz, -0.016f, 0.0f, -9.8f) == -1);
    check("dt_naive_zero_error",
          rope_step_naive(&rz, 0.0f, 0.0f, -9.8f) == -1);
    check("dt_positive_ok", rope_step(&rz, 0.016f, 0.0f, -9.8f) == 0);

    /* ---- reproduce the bug with the naive solver ---- */
    Rope naive = make_cannon(64.0f, 0.1f);
    int tunnel_frames = 0;
    for (int k = 0; k < 20; k++) {
        float before_x = naive.p[0].x;
        rope_step_naive(&naive, 0.1f, 0.0f, 0.0f);
        /* frame started outside on the left, ended on the far side without
           the solver ever holding the surface: a pass-through event */
        if (before_x < -unit.r && naive.p[0].x > unit.r) tunnel_frames++;
    }
    check("naive_tunnels", tunnel_frames > 0);

    /* ---- CCD fix: same cannon cannot pass through ---- */
    Rope ccd = make_cannon(64.0f, 0.1f);
    float ccd_worst_pen = 0.0f, ccd_worst_raw = 0.0f;
    int ccd_reversed = 0, ccd_any_contact = 0, ccd_cap_violations = 0;
    for (int k = 0; k < 20; k++) {
        rope_step(&ccd, 0.1f, 0.0f, 0.0f);
        /* diagnostics describe this step only: accumulate on the frame
           where contact actually happens, before the next step resets them */
        if (ccd.contacts > 0) {
            if (ccd.worst_raw_pen > ccd_worst_raw) ccd_worst_raw = ccd.worst_raw_pen;
            ccd_reversed += ccd.reversed;
            ccd_any_contact += ccd.contacts;
        }
        if (ccd.worst_pen > ccd_worst_pen) ccd_worst_pen = ccd.worst_pen;
        if (ccd.worst_pen > ccd.pen_cap + 1e-6f) ccd_cap_violations++;
    }
    check("ccd_no_cap_violation", ccd_cap_violations == 0);
    check("ccd_stops_at_surface", ccd.p[0].x <= unit.r + 1e-4f);
    check("ccd_pen_under_cap", ccd_worst_pen <= ccd.pen_cap + 1e-6f);
    check("ccd_residual_near_zero", ccd_worst_pen < 1e-4f);
    check("ccd_raw_was_deep", ccd_worst_raw > 0.5f);
    check("ccd_made_contact", ccd_any_contact > 0);
    check("ccd_no_velocity_reversal", ccd_reversed == 0);
    /* head-on rollback leaves no normal velocity at all */
    check("ccd_headon_normal_zero",
          fabsf((ccd.p[0].x - ccd.p[0].px)) < 1e-4f);

    /* ---- glancing hit preserves tangential velocity, no reversal ---- */
    Rope gl;
    rope_init(&gl, 1, 1.0f, -4.0f, 0.5f);
    gl.p[0].pinned = 0;
    gl.p[0].inv_mass = 1.0f;
    gl.damping = 1.0f;
    rope_set_obstacle(&gl, v2(0.0f, 0.0f), 1.0f);
    gl.p[0].px = gl.p[0].x - 64.0f * 0.1f;
    gl.p[0].py = gl.p[0].y;
    rope_step(&gl, 0.1f, 0.0f, 0.0f);
    {
        V2 v = vsub(v2(gl.p[0].x, gl.p[0].y), v2(gl.p[0].px, gl.p[0].py));
        V2 m = vsub(v2(gl.p[0].x, gl.p[0].y), unit.c);
        float ml = vlen(m);
        V2 nn = vmul(m, 1.0f / ml);
        float vn = vdot(v, nn);
        V2 tangent = vsub(v, vmul(nn, vn));
        /* contact normal is (-0.866, 0.5), tangent is (0.5, 0.866); the
           (6.4, 0) incoming velocity keeps 6.4 * 0.5 = 3.2 tangential */
        check("glance_touches_surface", fabsf(ml - 1.0f) < 1e-4f);
        check("glance_no_reversal", vn <= 1e-4f);
        check("glance_tangent_kept", fabsf(vlen(tangent) - 3.2f) < 1e-3f);
    }

    /* ---- full rope, violent downward swing through a peg ---- */
    Rope rope;
    rope_init(&rope, 6, 1.0f, 0.0f, 4.0f);
    rope_set_obstacle(&rope, v2(3.0f, 1.0f), 0.8f);
    for (int i = 1; i < rope.n; i++) {
        rope.p[i].px = rope.p[i].x;
        rope.p[i].py = rope.p[i].y + 300.0f * 0.016f;   /* hard downward kick */
    }
    float rope_worst_pen = 0.0f, rope_worst_raw = 0.0f;
    int rope_reversed = 0, rope_contacts = 0, cap_violations = 0;
    for (int k = 0; k < 120; k++) {
        rope_step(&rope, 0.016f, 0.0f, -9.8f);
        if (rope.worst_pen > rope_worst_pen) rope_worst_pen = rope.worst_pen;
        if (rope.contacts > 0 && rope.worst_raw_pen > rope_worst_raw)
            rope_worst_raw = rope.worst_raw_pen;
        rope_reversed += rope.reversed;
        rope_contacts += rope.contacts;
        if (rope.worst_pen > rope.pen_cap + 1e-6f) cap_violations++;
    }
    check("rope_made_contact", rope_contacts > 0);
    check("rope_pen_under_cap", rope_worst_pen <= rope.pen_cap + 1e-6f);
    check("rope_no_cap_violation", cap_violations == 0);
    check("rope_raw_was_deep", rope_worst_raw > 0.3f);
    check("rope_no_velocity_reversal", rope_reversed == 0);
    check("rope_stays_in_stretch", chain_seg_error(&rope) < 0.05f);
    check("rope_anchor_pinned",
          fabsf(rope.p[0].x) < 1e-6f && fabsf(rope.p[0].y - 4.0f) < 1e-6f);

    /* ---- reference comparison: coarse CCD vs fine substeps ---- */
    Rope fine = make_cannon(64.0f, 0.1f);
    fine.damping = 1.0f;
    Rope coarse = make_cannon(64.0f, 0.1f);
    coarse.damping = 1.0f;
    const int SUB = 64;
    float fine_worst_pen = 0.0f, coarse_worst = 0.0f;
    for (int k = 0; k < 20; k++) {
        for (int s = 0; s < SUB; s++)
            rope_step(&fine, 0.1f / (float)SUB, 0.0f, 0.0f);
        rope_step(&coarse, 0.1f, 0.0f, 0.0f);
        if (fine.worst_pen > fine_worst_pen) fine_worst_pen = fine.worst_pen;
        if (coarse.worst_pen > coarse_worst) coarse_worst = coarse.worst_pen;
    }
    check("fine_ref_no_pen", fine_worst_pen < 1e-5f);
    check("coarse_pen_matches_ref_order", coarse_worst < 1e-4f);
    check("coarse_tracks_fine",
          fabsf(coarse.p[0].x - fine.p[0].x) < 0.1f
          && fabsf(coarse.p[0].y - fine.p[0].y) < 0.1f);

    /* ---- no obstacle: ordinary verlet behaviour is unchanged ---- */
    Rope fall;
    rope_init(&fall, 8, 1.0f, 0.0f, 0.0f);
    for (int k = 0; k < 200; k++) rope_step(&fall, 0.016f, 0.0f, -9.8f);
    check("no_obs_stretch", chain_seg_error(&fall) < 0.02f);
    check("no_obs_pinned",
          fabsf(fall.p[0].x) < 1e-6f && fabsf(fall.p[0].y) < 1e-6f);
    check("no_obs_falls", fall.p[7].y < -1.0f);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
