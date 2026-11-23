/*
 * invalid: every illegal-parameter path of the unified entry.
 *
 * For each case the program prints:
 *   <module>.<case> status=<n> msg="<prefix>: ..."
 * Every case must return MAL_ERR_PARAM (or MAL_ERR_UNSUPPORTED for unknown
 * ops/modules) and the message must start with the owning module name,
 * proving the entry points at the module that rejected the input.
 */
#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "solver.h"
#include "skeleton.h"
#include "ik.h"
#include "curve.h"
#include "constraint.h"

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

static void expect(const char* tag, mal_status want, const char* prefix,
                   const mal_request* r) {
    char why[192] = {0};
    mal_status got = mal_solve(r, why, sizeof why);
    int ok = got == want && strncmp(why, prefix, strlen(prefix)) == 0;
    printf("%s status=%d msg=\"%s\" ok=%d\n", tag, (int)got, why, ok);
    if (!ok) failures++;
}

int main(void) {
    /* ---- cross-cutting solver checks ---- */
    expect("solver.null_request", MAL_ERR_PARAM, "solver:", 0);

    mal_request badmod = { (mal_module)42, 0, (const void*)1, 0 };
    expect("solver.bad_module_id", MAL_ERR_PARAM, "solver:", &badmod);

    mal_skeleton sk;
    mal_skeleton_init(&sk);
    mal_skel_update_req uq = { &sk };
    mal_request nullpayload = { MAL_MODULE_SKELETON, 0, 0, 0 };
    expect("solver.null_payload", MAL_ERR_PARAM, "skeleton:", &nullpayload);

    /* ---- skeleton ---- */
    mal_skel_add_wrap saw; saw.s = &sk;
    saw.in.parent = 99; saw.in.x = 0; saw.in.y = 0; saw.in.ang = 0;
    mal_request rsk = { MAL_MODULE_SKELETON, MAL_OP_SKELETON_ADD, &saw, 0 };
    expect("skeleton.bad_parent", MAL_ERR_PARAM, "skeleton:", &rsk);

    saw.in.parent = -1; saw.in.x = NAN;
    expect("skeleton.non_finite", MAL_ERR_PARAM, "skeleton:", &rsk);
    saw.in.x = 0;

    /* fill to capacity, then add must fail with STATE, message still named */
    mal_skeleton full;
    mal_skeleton_init(&full);
    int last = -1;
    for (int i = 0; i < MAL_SKEL_MAX_BONES; i++) {
        mal_skel_add_wrap fw; fw.s = &full;
        fw.in.parent = last; fw.in.x = 1; fw.in.y = 0; fw.in.ang = 0;
        mal_request fr = { MAL_MODULE_SKELETON, MAL_OP_SKELETON_ADD, &fw, 0 };
        char why[128];
        if (mal_solve(&fr, why, sizeof why) != MAL_OK) { printf("setup_failure=%s\n", why); return 2; }
        last = i;
    }
    mal_skel_add_wrap fw; fw.s = &full;
    fw.in.parent = 0; fw.in.x = 1; fw.in.y = 0; fw.in.ang = 0;
    mal_request fr = { MAL_MODULE_SKELETON, MAL_OP_SKELETON_ADD, &fw, 0 };
    expect("skeleton.capacity", MAL_ERR_STATE, "skeleton:", &fr);

    mal_skeleton empty;
    mal_skeleton_init(&empty);
    mal_skel_update_req euq = { &empty };
    mal_request reu = { MAL_MODULE_SKELETON, MAL_OP_SKELETON_UPDATE, &euq, 0 };
    expect("skeleton.empty_update", MAL_ERR_STATE, "skeleton:", &reu);

    mal_request skop = { MAL_MODULE_SKELETON, 77, &uq, 0 };
    expect("skeleton.bad_op", MAL_ERR_UNSUPPORTED, "skeleton:", &skop);

    /* ---- ik ---- */
    mal_v2 root = {0,0}, tgt = {1,0}, pole = {0,1};
    mal_ik_req ikbad = { root, tgt, pole, -1.0f, 1.0f };
    mal_ik_resp ikout;
    mal_request rik = { MAL_MODULE_IK, MAL_OP_IK_SOLVE, &ikbad, &ikout };
    expect("ik.bad_length", MAL_ERR_PARAM, "ik:", &rik);

    ikbad.l1 = 1.0f; ikbad.target.x = INFINITY;
    expect("ik.non_finite", MAL_ERR_PARAM, "ik:", &rik);
    ikbad.target.x = 1.0f;

    mal_request riknull = { MAL_MODULE_IK, MAL_OP_IK_SOLVE, &ikbad, 0 };
    expect("ik.null_resp", MAL_ERR_PARAM, "ik:", &riknull);

    mal_request rikop = { MAL_MODULE_IK, 5, &ikbad, &ikout };
    expect("ik.bad_op", MAL_ERR_UNSUPPORTED, "ik:", &rikop);

    /* ---- curve ---- */
    mal_curve cv;
    mal_curve_init(&cv);
    mal_curve_add_wrap cw; cw.c = &cv;
    cw.in.t = 0; cw.in.v = 0; cw.in.tan_in = 0; cw.in.tan_out = 1;
    mal_request rc = { MAL_MODULE_CURVE, MAL_OP_CURVE_ADD, &cw, 0 };
    char why[128];
    if (mal_solve(&rc, why, sizeof why) != MAL_OK) { printf("setup_failure=%s\n", why); return 2; }
    cw.in.t = 0;   /* not strictly increasing */
    expect("curve.unsorted", MAL_ERR_PARAM, "curve:", &rc);

    cw.in.t = 1; cw.in.v = NAN;
    expect("curve.non_finite", MAL_ERR_PARAM, "curve:", &rc);
    cw.in.v = 1;

    mal_curve_eval_wrap ew; ew.c = 0; ew.in.t = 0.5f;
    mal_request rev = { MAL_MODULE_CURVE, MAL_OP_CURVE_EVAL, &ew, 0 };
    expect("curve.null_instance", MAL_ERR_PARAM, "curve:", &rev);

    mal_curve blank;
    mal_curve_init(&blank);
    mal_curve_eval_wrap bew; bew.c = &blank; bew.in.t = 0.5f;
    mal_request rbe = { MAL_MODULE_CURVE, MAL_OP_CURVE_EVAL, &bew, 0 };
    expect("curve.empty_eval", MAL_ERR_STATE, "curve:", &rbe);

    /* ---- constraint ---- */
    mal_chain ch;
    mal_chain_init_wrap iw; iw.ch = &ch;
    iw.in.n = 1; iw.in.seg = 1; iw.in.x0 = 0; iw.in.y0 = 0;
    mal_request ri = { MAL_MODULE_CONSTRAINT, MAL_OP_CHAIN_INIT, &iw, 0 };
    expect("constraint.bad_count", MAL_ERR_PARAM, "constraint:", &ri);

    iw.in.n = 8; iw.in.seg = 0;
    expect("constraint.bad_seg", MAL_ERR_PARAM, "constraint:", &ri);

    iw.in.seg = 1; iw.in.x0 = INFINITY;
    expect("constraint.non_finite_anchor", MAL_ERR_PARAM, "constraint:", &ri);
    iw.in.x0 = 0;
    if (mal_solve(&ri, why, sizeof why) != MAL_OK) { printf("setup_failure=%s\n", why); return 2; }

    mal_chain_step_wrap sw; sw.ch = &ch;
    sw.in.dt = 0; sw.in.gx = 0; sw.in.gy = -9.8f;
    mal_request rs = { MAL_MODULE_CONSTRAINT, MAL_OP_CHAIN_STEP, &sw, 0 };
    expect("constraint.bad_dt", MAL_ERR_PARAM, "constraint:", &rs);

    sw.in.dt = 0.016f; sw.in.gy = NAN;
    expect("constraint.non_finite_gravity", MAL_ERR_PARAM, "constraint:", &rs);
    sw.in.gy = -9.8f;

    mal_chain_step_wrap nsw; nsw.ch = 0;
    nsw.in.dt = 0.016f; nsw.in.gx = 0; nsw.in.gy = 0;
    mal_request rns = { MAL_MODULE_CONSTRAINT, MAL_OP_CHAIN_STEP, &nsw, 0 };
    expect("constraint.null_instance", MAL_ERR_PARAM, "constraint:", &rns);

    mal_request chop = { MAL_MODULE_CONSTRAINT, 42, &sw, 0 };
    expect("constraint.bad_op", MAL_ERR_UNSUPPORTED, "constraint:", &chop);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
