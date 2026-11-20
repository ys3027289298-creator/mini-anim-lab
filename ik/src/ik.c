
/* ik: analytic two-bone IK with pole vector and joint angle limits. */
#include <stdio.h>
#include <math.h>
#include <string.h>

typedef struct { float x, y; } V2;

static float vlen(V2 v) { return sqrtf(v.x * v.x + v.y * v.y); }
static V2 vsub(V2 a, V2 b) { V2 r = { a.x - b.x, a.y - b.y }; return r; }
static V2 vnorm(V2 v) { float l = vlen(v); V2 r = { 0, 0 }; if (l > 0) { r.x = v.x / l; r.y = v.y / l; } return r; }

typedef struct {
    int ok;
    int err;    /* non-zero: invalid argument (e.g. zero pole tolerance) */
    V2 joint;
    float a1;   /* angle of first bone */
    float a2;   /* angle of second bone */
} IkResult;

/* joint position for a given bend side; a = angle between aim dir and bone 1 */
static V2 joint_for_side(V2 root, V2 dir, float l1, float a, float side) {
    float ang = atan2f(dir.y, dir.x) + side * a;
    V2 j = { root.x + l1 * cosf(ang), root.y + l1 * sinf(ang) };
    return j;
}

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

    r.joint = joint_for_side(root, dir, l1, a, side);
    r.a1 = atan2f(r.joint.y - root.y, r.joint.x - root.x);
    V2 d2 = vsub(target, r.joint);
    r.a2 = atan2f(d2.y, d2.x);
    r.ok = 1;
    return r;
}

/* per-chain state for time-continuous solving */
typedef struct {
    float side;     /* current bend side, +1/-1; 0 = uninitialized */
    float blend;    /* 1 = settled on current side, <1 = transitioning */
} IkChainState;

/* Time-continuous solve: keeps the previous frame's bend side, only flips
   when the pole cross product leaves a hysteresis band (|cross| > pole_eps),
   and eases the joint across when a flip does happen.
   pole_eps must be > 0: a zero tolerance would re-enable per-frame flipping,
   so it is reported as an error. blend_rate is in 1/seconds. */
static IkResult solve_two_bone_anim(V2 root, V2 target, float l1, float l2, V2 pole,
                                    float pole_eps, float blend_rate, float dt,
                                    IkChainState* st) {
    IkResult r;
    memset(&r, 0, sizeof(r));
    if (pole_eps <= 0.0f) { r.err = 1; return r; }        /* zero tolerance: error */
    V2 d = vsub(target, root);
    float dist = vlen(d);
    if (l1 <= 0 || l2 <= 0) return r;
    if (dist > l1 + l2 + 1e-5f) return r;
    if (dist < fabsf(l1 - l2) - 1e-5f) return r;
    if (dist < 1e-6f) return r;

    float cos_a = (l1 * l1 + dist * dist - l2 * l2) / (2 * l1 * dist);
    if (cos_a > 1.0f) cos_a = 1.0f;
    if (cos_a < -1.0f) cos_a = -1.0f;
    float a = acosf(cos_a);

    V2 dir = vnorm(d);
    float cross = dir.x * pole.y - dir.y * pole.x;

    /* hysteresis: outside the dead zone the pole side wins; inside, keep previous */
    if (st->side == 0.0f) {
        st->side = cross >= 0.0f ? 1.0f : -1.0f;
        st->blend = 1.0f;
    } else if (fabsf(cross) > pole_eps) {
        float want = cross > 0.0f ? 1.0f : -1.0f;
        if (want != st->side) {
            st->side = want;
            st->blend = 0.0f;                             /* start transition */
        }
    }

    V2 jcur = joint_for_side(root, dir, l1, a, st->side);
    if (st->blend < 1.0f) {
        st->blend += blend_rate * dt;
        if (st->blend > 1.0f) st->blend = 1.0f;
        float t = st->blend * st->blend * (3.0f - 2.0f * st->blend);  /* smoothstep */
        V2 jold = joint_for_side(root, dir, l1, a, -st->side);
        r.joint.x = jold.x + (jcur.x - jold.x) * t;
        r.joint.y = jold.y + (jcur.y - jold.y) * t;
    } else {
        r.joint = jcur;
    }

    r.a1 = atan2f(r.joint.y - root.y, r.joint.x - root.x);
    V2 d2 = vsub(target, r.joint);
    r.a2 = atan2f(d2.y, d2.x);
    r.ok = 1;
    return r;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

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

    /* zero pole tolerance must be reported as an error */
    IkChainState st0 = { 0, 0 };
    IkResult ztol = solve_two_bone_anim(root, target, 1.0f, 1.0f, pole, 0.0f, 4.0f, 1.0f / 60.0f, &st0);
    check("zero_pole_eps_err", ztol.ok == 0 && ztol.err == 1);

    /* sweep the target across the pole line: no frame may flip the joint
       discontinuously, and settled frames must match the analytic solution */
    {
        IkChainState st = { 0, 0 };
        const int frames = 240;
        const float dt = 1.0f / 60.0f;
        const float pole_eps = 0.05f, blend_rate = 4.0f;
        int flip_frames = 0, analytic_bad = 0, solved = 1;
        V2 prev = { 0, 0 };
        float prev_target_y = 0.0f;
        for (int i = 0; i < frames; i++) {
            float ty = -0.6f + 1.2f * (float)i / (float)(frames - 1);  /* crosses y=0 */
            V2 t = { 1.0f, ty };
            IkResult r = solve_two_bone_anim(root, t, 1.0f, 1.0f, pole,
                                             pole_eps, blend_rate, dt, &st);
            if (!r.ok) { solved = 0; break; }
            if (i > 0) {
                float jump = vlen(vsub(r.joint, prev));
                float target_step = fabsf(ty - prev_target_y);
                /* joint may track the target, but never leap to the other side */
                if (jump > target_step + 0.1f) flip_frames++;
            }
            /* settled (not blending) -> must equal analytic joint for current side */
            if (st.blend == 1.0f) {
                IkResult ref = solve_two_bone(root, t, 1.0f, 1.0f,
                                              st.side > 0 ? (V2){ 0, 1 } : (V2){ 0, -1 });
                if (!ref.ok || vlen(vsub(r.joint, ref.joint)) > 1e-4f) analytic_bad++;
            }
            prev = r.joint;
            prev_target_y = ty;
        }
        check("sweep_solved", solved);
        check("flip_frames_zero", flip_frames == 0);
        check("matches_analytic", analytic_bad == 0);
    }

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
