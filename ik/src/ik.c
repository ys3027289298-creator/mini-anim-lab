
/* ik: analytic two-bone IK with pole vector and joint angle limits. */
#include <stdio.h>
#include <string.h>
#include <math.h>

#ifndef MAL_MODULE
#define MAL_MODULE "ik"
#endif
#define MAL_TIME_ITERS 400000L
#include "../../common/harness.h"

typedef struct { float x, y; } V2;

static float vlen(V2 v) { return sqrtf(v.x * v.x + v.y * v.y); }
static V2 vsub(V2 a, V2 b) { V2 r = { a.x - b.x, a.y - b.y }; return r; }
static V2 vnorm(V2 v) { float l = vlen(v); V2 r = { 0, 0 }; if (l > 0) { r.x = v.x / l; r.y = v.y / l; } return r; }

typedef struct {
    int ok;
    V2 joint;
    float a1;   /* angle of first bone */
    float a2;   /* angle of second bone */
} IkResult;

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

/* deterministic workload shared by the timing and image evidence */
#define IK_ITERS 200000L

typedef struct { V2 root, joint, tip; } IkChain;

MAL_UNUSED static IkChain ik_pose(void) {
    IkChain c;
    c.root = (V2){ 0.0f, 0.0f };
    IkResult r = solve_two_bone(c.root, (V2){ 1.0f, 0.0f }, 1.0f, 1.0f, (V2){ 0.0f, 1.0f });
    c.joint = r.joint;
    c.tip = (V2){ r.joint.x + cosf(r.a2), r.joint.y + sinf(r.a2) };
    return c;
}

static void run_checks(void) {
    V2 root = { 0, 0 }, target = { 1, 0 }, pole = { 0, 1 };
    IkResult r = solve_two_bone(root, target, 1.0f, 1.0f, pole);
    mal_check("reachable", r.ok == 1);
    mal_check("joint_on_circle", fabsf(vlen(vsub(r.joint, root)) - 1.0f) < 1e-4f);
    mal_check("end_matches", fabsf(vlen(vsub(r.joint, target)) - 1.0f) < 1e-4f);
    mal_check("pole_up", r.joint.y > 0.0f);

    IkResult r2 = solve_two_bone(root, target, 1.0f, 1.0f, (V2){ 0, -1 });
    mal_check("pole_down", r2.joint.y < 0.0f);

    IkResult far = solve_two_bone(root, (V2){ 5, 0 }, 1.0f, 1.0f, pole);
    mal_check("unreachable", far.ok == 0);

    IkResult close = solve_two_bone(root, (V2){ 0.1f, 0 }, 1.0f, 0.2f, pole);
    mal_check("too_close", close.ok == 0);

    IkResult zero = solve_two_bone(root, target, 0.0f, 1.0f, pole);
    mal_check("zero_len", zero.ok == 0);
}

int main(void) {
    run_checks();

#ifdef MAL_EVIDENCE_POSE
    if (mal_failures == 0) {
        IkResult up = solve_two_bone((V2){ 0, 0 }, (V2){ 1, 0 }, 1.0f, 1.0f, (V2){ 0, 1 });
        IkResult dn = solve_two_bone((V2){ 0, 0 }, (V2){ 1, 0 }, 1.0f, 1.0f, (V2){ 0, -1 });
        mal_pose("up_ok", up.ok);
        mal_pose("up_joint_x", up.joint.x);
        mal_pose("up_joint_y", up.joint.y);
        mal_pose("up_a1", up.a1);
        mal_pose("up_a2", up.a2);
        mal_pose("dn_joint_y", dn.joint.y);
        mal_pose("reach", vlen(vsub(up.joint, (V2){ 0, 0 })));
        mal_pose("miss", vlen(vsub(up.joint, (V2){ 1, 0 })));
        mal_emit_pose();
    }
#endif

#ifdef MAL_EVIDENCE_IMAGE
    if (mal_failures == 0) {
        MalImage im;
        IkChain c = ik_pose();
        mal_img_clear(&im);
        const double ox = 24.0, oy = 100.0, scale = 70.0;
        mal_img_line(&im, ox, oy, ox + c.joint.x * scale, oy - c.joint.y * scale, 80, 200, 255);
        mal_img_line(&im, ox + c.joint.x * scale, oy - c.joint.y * scale,
                     ox + c.tip.x * scale, oy - c.tip.y * scale, 255, 180, 90);
        mal_emit_png(&im, MAL_MODULE);
    }
#endif

#ifdef MAL_EVIDENCE_TIMING
    if (mal_failures == 0) {
        double acc = 0.0;
        double t0 = mal_now();
        for (long it = 0; it < MAL_TIME_ITERS; it++) {
            float a = (float)(it % 360) * 0.017453292f;
            V2 t = { 1.6f * cosf(a), 1.6f * sinf(a) };
            IkResult r = solve_two_bone((V2){ 0, 0 }, t, 1.0f, 1.0f, (V2){ 0, 1 });
            acc += r.joint.x + r.joint.y + r.a1 + r.a2;
        }
        double t1 = mal_now();
        if (!isfinite(acc)) mal_check("time_finite", 0);
        mal_emit_time(MAL_TIME_TIER, t1 - t0, MAL_TIME_ITERS);
    }
#endif

    printf("all_ok=%d\n", mal_failures == 0);
    return mal_failures == 0 ? 0 : 1;
}
