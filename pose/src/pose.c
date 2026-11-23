/*
 * pose: cross-module consistency check for the unified solver.
 *
 * One binary links all four modules. Each module source is compiled with
 * MAL_NO_MAIN and MAL_UNIFIED=1, so they contribute their wrappers and
 * per-module dispatch entries but no standalone main(); this file supplies
 * the merged dispatch table consumed by mal_solve().
 *
 * Scenario: an IK chain and a Verlet distance-constraint chain are asked
 * for the same end-effector target, and a skeleton is posed from the IK
 * angles. The three poses must agree; a curve drives the moving target.
 */
#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "solver.h"
#include "skeleton.h"
#include "ik.h"
#include "curve.h"
#include "constraint.h"

/* Merged dispatch table: every module linked into this binary. */
mal_status skeleton_dispatch(int, const void*, void*, char*, int);
mal_status ik_dispatch(int, const void*, void*, char*, int);
mal_status curve_dispatch(int, const void*, void*, char*, int);
mal_status chain_dispatch(int, const void*, void*, char*, int);

const mal_dispatch_entry mal_dispatch_table[] = {
    { MAL_MODULE_SKELETON,   skeleton_dispatch },
    { MAL_MODULE_IK,         ik_dispatch },
    { MAL_MODULE_CURVE,      curve_dispatch },
    { MAL_MODULE_CONSTRAINT, chain_dispatch },
    { (mal_module)0, 0 }
};

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static float vdist(mal_v2 a, mal_v2 b) {
    float dx = a.x - b.x, dy = a.y - b.y;
    return sqrtf(dx * dx + dy * dy);
}

/* Every call in this file goes through the unified entry. */
static mal_status solve(const mal_request* r) {
    char why[160];
    mal_status st = mal_solve(r, why, sizeof why);
    if (st != MAL_OK) printf("unexpected_error=%s\n", why);
    return st;
}

static double now_s(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

/* One full pose round: curve target -> IK -> skeleton -> constraint chain.
 *
 * The analytic IK bends toward `pole`. A Verlet distance chain alone has
 * no bend preference (the straight configuration is a valid solution), so
 * the chain is driven by the same pole: its middle particle gets a spring
 * pull toward the IK elbow while the end particle is pinned to the target.
 * With equal rest lengths the distance constraints then converge onto the
 * IK/skeleton pose. That convergence is what the cross-module check
 * asserts.
 */
typedef struct {
    mal_v2 ik_joint;      /* elbow from analytic IK */
    mal_v2 skel_joint;    /* elbow from skeleton world transforms */
    mal_v2 chain_joint;   /* middle particle from distance constraints */
    mal_v2 end;           /* end effector position */
} pose_row;

static int pose_once(mal_v2 root, mal_v2 target, mal_v2 pole,
                     int chain_iters, pose_row* row) {
    /* 1. IK through the unified entry. */
    mal_ik_req ikq = { root, target, pole, 1.0f, 1.0f };
    mal_ik_resp ikr;
    mal_request rik = { MAL_MODULE_IK, MAL_OP_IK_SOLVE, &ikq, &ikr };
    if (solve(&rik) != MAL_OK || !ikr.ok) return 0;
    row->ik_joint = ikr.joint;

    /* 2. Skeleton posed from the IK world angles. Two unit bones:
     * root at `root` rotated by a1; child offset 1 along root +x with
     * relative angle (a2 - a1). Bone i's world position is its start
     * joint, matching the IK root/joint/end points. */
    mal_skeleton sk;
    mal_skeleton_init(&sk);
    mal_skel_add_wrap aw; aw.s = &sk;
    aw.in.parent = -1; aw.in.x = root.x; aw.in.y = root.y; aw.in.ang = ikr.a1;
    mal_request ra = { MAL_MODULE_SKELETON, MAL_OP_SKELETON_ADD, &aw, 0 };
    if (solve(&ra) != MAL_OK) return 0;
    aw.in.parent = 0; aw.in.x = 1.0f; aw.in.y = 0.0f;
    aw.in.ang = ikr.a2 - ikr.a1;
    if (solve(&ra) != MAL_OK) return 0;
    mal_skel_update_req uq = { &sk };
    mal_request ru = { MAL_MODULE_SKELETON, MAL_OP_SKELETON_UPDATE, &uq, 0 };
    if (solve(&ru) != MAL_OK) return 0;
    row->skel_joint = sk.world[1];
    /* end = joint + unit vector along the second bone's world angle */
    row->end.x = sk.world[1].x + cosf(sk.world_ang[1]);
    row->end.y = sk.world[1].y + sinf(sk.world_ang[1]);

    /* 3. Constraint chain of the same two bone lengths, same pole drive. */
    mal_chain ch;
    mal_chain_init_wrap iw; iw.ch = &ch;
    iw.in.n = 3; iw.in.seg = 1.0f; iw.in.x0 = root.x; iw.in.y0 = root.y;
    mal_request ri = { MAL_MODULE_CONSTRAINT, MAL_OP_CHAIN_INIT, &iw, 0 };
    if (solve(&ri) != MAL_OK) return 0;
    mal_chain_set_iterations(&ch, chain_iters);

    /* Bend force on the middle particle points at the IK elbow (the pole
     * solution); the end particle is pinned to the target and the anchor
     * is already pinned by chain init. Gravity is zero, so only the drive
     * plus distance constraints shape the chain. */
    const float dt = 0.016f;
    const float k = 30.0f;
    ch.p[2].x = target.x;
    ch.p[2].y = target.y;
    ch.p[2].px = target.x;
    ch.p[2].py = target.y;
    ch.p[2].pinned = 1;
    ch.p[2].inv_mass = 0.0f;
    mal_chain_step_wrap sw; sw.ch = &ch;
    sw.in.dt = dt; sw.in.gx = 0.0f; sw.in.gy = 0.0f;
    mal_request rs = { MAL_MODULE_CONSTRAINT, MAL_OP_CHAIN_STEP, &sw, 0 };
    for (int i = 0; i < 600; i++) {
        float fx = k * (ikr.joint.x - ch.p[1].x);
        float fy = k * (ikr.joint.y - ch.p[1].y);
        ch.p[1].x += fx * dt * dt;
        ch.p[1].y += fy * dt * dt;
        if (solve(&rs) != MAL_OK) return 0;
    }
    mal_v2 mid;
    mal_chain_get(&ch, 1, &mid, 0, 0);
    row->chain_joint = mid;
    return 1;
}

int main(void) {
    mal_v2 root = { 0.0f, 0.0f };
    mal_v2 pole = { 0.0f, 1.0f };

    /* Curve drives the target x from 0.6 to 1.4 (both reachable). */
    mal_curve cv;
    mal_curve_init(&cv);
    mal_curve_add_wrap cw; cw.c = &cv;
    cw.in.t = 0.0f; cw.in.v = 0.6f; cw.in.tan_in = 0.0f; cw.in.tan_out = 1.0f;
    mal_request rc = { MAL_MODULE_CURVE, MAL_OP_CURVE_ADD, &cw, 0 };
    check("curve_key0", solve(&rc) == MAL_OK);
    cw.in.t = 1.0f; cw.in.v = 1.4f; cw.in.tan_in = 1.0f; cw.in.tan_out = 0.0f;
    check("curve_key1", solve(&rc) == MAL_OK);

    pose_row rows[5];
    int consistent = 1;
    for (int k = 0; k < 5; k++) {
        float t = k / 4.0f;
        mal_curve_eval_wrap ew; ew.c = &cv; ew.in.t = t;
        mal_curve_eval_resp ev;
        mal_request re = { MAL_MODULE_CURVE, MAL_OP_CURVE_EVAL, &ew, &ev };
        if (solve(&re) != MAL_OK) { consistent = 0; break; }
        mal_v2 target = { ev.v, 0.0f };
        if (!pose_once(root, target, pole, 8, &rows[k])) {
            consistent = 0; break;
        }
        float e_ik_skel = vdist(rows[k].ik_joint, rows[k].skel_joint);
        float e_ik_chain = vdist(rows[k].ik_joint, rows[k].chain_joint);
        float e_end = vdist(rows[k].end, target);
        if (e_ik_skel > 1e-3f || e_ik_chain > 0.05f || e_end > 1e-3f)
            consistent = 0;
        printf("pose_%d target=(%.4f,0) ik=(%.4f,%.4f) skel=(%.4f,%.4f) "
               "chain=(%.4f,%.4f) d_ik_skel=%.5f d_ik_chain=%.5f d_end=%.5f\n",
               k, ev.v,
               rows[k].ik_joint.x, rows[k].ik_joint.y,
               rows[k].skel_joint.x, rows[k].skel_joint.y,
               rows[k].chain_joint.x, rows[k].chain_joint.y,
               vdist(rows[k].ik_joint, rows[k].skel_joint),
               vdist(rows[k].ik_joint, rows[k].chain_joint),
               vdist(rows[k].end, target));
    }
    check("poses_consistent", consistent);

    /* Illegal parameters through the unified entry must name the module. */
    char why[160];
    mal_ik_req bad = { root, root, pole, -1.0f, 1.0f };
    mal_ik_resp badout;
    mal_request rb = { MAL_MODULE_IK, MAL_OP_IK_SOLVE, &bad, &badout };
    check("ik_param_error",
          mal_solve(&rb, why, sizeof why) == MAL_ERR_PARAM &&
          strncmp(why, "ik:", 3) == 0);

    mal_chain_step_wrap sw2; sw2.ch = 0; sw2.in.dt = 0.016f;
    sw2.in.gx = 0.0f; sw2.in.gy = 0.0f;
    mal_request rch = { MAL_MODULE_CONSTRAINT, MAL_OP_CHAIN_STEP, &sw2, 0 };
    check("constraint_param_error",
          mal_solve(&rch, why, sizeof why) == MAL_ERR_PARAM &&
          strncmp(why, "constraint:", 11) == 0);

    check("solver_null_request",
          mal_solve(0, why, sizeof why) == MAL_ERR_PARAM &&
          strncmp(why, "solver:", 7) == 0);

    /* ---- timing, two tiers ------------------------------------------- */
    pose_row tmp;
    mal_v2 tgt = { 1.0f, 0.0f };

    double t0 = now_s();
    pose_once(root, tgt, pole, 8, &tmp);          /* tier 1: one round */
    double t1 = now_s();

    int N = 2000;
    double t2 = now_s();
    for (int i = 0; i < N; i++)
        pose_once(root, tgt, pole, 8, &tmp);      /* tier 2: bulk rounds */
    double t3 = now_s();
    printf("time_one_round_us=%.3f\n", (t1 - t0) * 1e6);
    printf("time_%d_rounds_ms=%.3f\n", N, (t3 - t2) * 1e3);
    printf("time_per_round_us=%.3f\n", (t3 - t2) * 1e6 / N);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
