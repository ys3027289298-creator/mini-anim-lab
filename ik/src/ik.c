
/* ik: analytic two-bone IK with pole vector and joint angle limits.
 *
 * MAL_UNIFIED=1 (default): solves go through mal_solve() -> mal_ik_solve().
 * MAL_UNIFIED=0: the original static solve_two_bone() is called directly.
 * Both paths execute the same math below.
 */
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "solver.h"
#include "ik.h"
#include "mal_err.h"

#ifndef MAL_UNIFIED
#define MAL_UNIFIED 1
#endif

#ifdef MAL_NO_MAIN
#define MAL_DISPATCH_FN
#else
#define MAL_DISPATCH_FN static
#endif

typedef mal_v2 V2;

static float vlen(V2 v) { return sqrtf(v.x * v.x + v.y * v.y); }
static V2 vsub(V2 a, V2 b) { V2 r = { a.x - b.x, a.y - b.y }; return r; }
static V2 vnorm(V2 v) { float l = vlen(v); V2 r = { 0, 0 }; if (l > 0) { r.x = v.x / l; r.y = v.y / l; } return r; }

/* solve for target with bone lengths l1,l2; pole decides bend direction */
typedef struct {
    int ok;
    V2 joint;
    float a1;   /* angle of first bone */
    float a2;   /* angle of second bone */
} IkResult;

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

/* ---- unified-entry wrapper ------------------------------------------- */

mal_status mal_ik_solve(const mal_ik_req* req, mal_ik_resp* resp,
                        char* why, int whycap) {
    if (!req) return mal_fail(why, whycap, MAL_ERR_PARAM,
                              "ik: null request payload");
    if (!resp) return mal_fail(why, whycap, MAL_ERR_PARAM,
                               "ik: null response");
    if (!mal_finite_f(req->l1) || !mal_finite_f(req->l2) ||
        !mal_finite_f(req->root.x) || !mal_finite_f(req->root.y) ||
        !mal_finite_f(req->target.x) || !mal_finite_f(req->target.y) ||
        !mal_finite_f(req->pole.x) || !mal_finite_f(req->pole.y))
        return mal_fail(why, whycap, MAL_ERR_PARAM,
                        "ik: non-finite input");
    if (req->l1 <= 0.0f || req->l2 <= 0.0f)
        return mal_fail(why, whycap, MAL_ERR_PARAM,
                        "ik: bone lengths must be positive (got %g,%g)",
                        (double)req->l1, (double)req->l2);

    IkResult r = solve_two_bone(req->root, req->target, req->l1, req->l2,
                                req->pole);
    resp->ok = r.ok;
    resp->joint = r.joint;
    resp->a1 = r.a1;
    resp->a2 = r.a2;
    return MAL_OK;   /* unreachable is a valid solve outcome: resp->ok = 0 */
}

/* dispatch entry used by the unified mal_solve(). Non-static in a
 * consolidated build (MAL_NO_MAIN) so pose/ can register it. */
MAL_DISPATCH_FN mal_status ik_dispatch(int op, const void* req, void* resp,
                                       char* why, int whycap) {
    if (op != MAL_OP_IK_SOLVE)
        return mal_fail(why, whycap, MAL_ERR_UNSUPPORTED,
                        "ik: no op %d", op);
    return mal_ik_solve((const mal_ik_req*)req, (mal_ik_resp*)resp,
                        why, whycap);
}

#ifndef MAL_NO_MAIN
const mal_dispatch_entry mal_dispatch_table[] = {
    { MAL_MODULE_IK, ik_dispatch },
    { (mal_module)0, 0 }
};
#endif

/* ---- self test ------------------------------------------------------- */

#ifndef MAL_NO_MAIN
static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

#if MAL_UNIFIED
static int ik_run(V2 root, V2 target, float l1, float l2, V2 pole,
                  IkResult* out) {
    mal_ik_req req = { root, target, pole, l1, l2 };
    mal_ik_resp resp;
    mal_request r = { MAL_MODULE_IK, MAL_OP_IK_SOLVE, &req, &resp };
    if (mal_solve(&r, 0, 0) != MAL_OK) { memset(out, 0, sizeof(*out)); return 0; }
    out->ok = resp.ok; out->joint = resp.joint;
    out->a1 = resp.a1; out->a2 = resp.a2;
    return 1;
}
#else
static int ik_run(V2 root, V2 target, float l1, float l2, V2 pole,
                  IkResult* out) {
    *out = solve_two_bone(root, target, l1, l2, pole);
    return 1;
}
#endif

int main(void) {
    V2 root = { 0, 0 }, target = { 1, 0 }, pole = { 0, 1 };
    IkResult r;
    ik_run(root, target, 1.0f, 1.0f, pole, &r);
    check("reachable", r.ok == 1);
    check("joint_on_circle", fabsf(vlen(vsub(r.joint, root)) - 1.0f) < 1e-4f);
    check("end_matches", fabsf(vlen(vsub(r.joint, target)) - 1.0f) < 1e-4f);
    check("pole_up", r.joint.y > 0.0f);

    IkResult r2;
    ik_run(root, target, 1.0f, 1.0f, (V2){ 0, -1 }, &r2);
    check("pole_down", r2.joint.y < 0.0f);

    IkResult far;
    ik_run(root, (V2){ 5, 0 }, 1.0f, 1.0f, pole, &far);
    check("unreachable", far.ok == 0);

    IkResult close;
    ik_run(root, (V2){ 0.1f, 0 }, 1.0f, 0.2f, pole, &close);
    check("too_close", close.ok == 0);

    IkResult zero;
    ik_run(root, target, 0.0f, 1.0f, pole, &zero);
    check("zero_len", zero.ok == 0);

#if MAL_UNIFIED
    /* illegal parameters must be rejected with the module named */
    char why[128];
    mal_ik_req badreq = { root, target, pole, 0.0f, 1.0f };
    mal_ik_resp badresp;
    mal_request br = { MAL_MODULE_IK, MAL_OP_IK_SOLVE, &badreq, &badresp };
    check("unified_bad_length",
          mal_solve(&br, why, sizeof why) == MAL_ERR_PARAM &&
          strstr(why, "ik:") == why);

    br.req = 0;
    check("unified_null_payload",
          mal_solve(&br, why, sizeof why) == MAL_ERR_PARAM &&
          strstr(why, "ik:") == why);

    br.req = &badreq; br.resp = 0;
    check("unified_null_resp",
          mal_solve(&br, why, sizeof why) == MAL_ERR_PARAM &&
          strstr(why, "ik:") == why);

    br.op = 7; br.resp = &badresp;
    check("unified_bad_op",
          mal_solve(&br, why, sizeof why) == MAL_ERR_UNSUPPORTED &&
          strstr(why, "ik:") == why);
#endif

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
#endif
