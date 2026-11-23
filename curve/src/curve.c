
/* curve: keyframe sampling with cubic hermite tangents and clamping.
 *
 * MAL_UNIFIED=1 (default): add/eval go through mal_solve().
 * MAL_UNIFIED=0: the original static curve_add/curve_eval are called
 * directly. The math below is shared by both paths.
 */
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "solver.h"
#include "curve.h"
#include "mal_err.h"

#ifndef MAL_UNIFIED
#define MAL_UNIFIED 1
#endif

#ifdef MAL_NO_MAIN
#define MAL_DISPATCH_FN
#else
#define MAL_DISPATCH_FN static
#endif

static void curve_init(mal_curve* c) { memset(c, 0, sizeof(*c)); }

static int curve_add(mal_curve* c, float t, float v, float tin, float tout) {
    if (c->count >= MAL_CURVE_MAX_KEYS) return -1;
    if (c->count > 0 && t <= c->keys[c->count - 1].t) return -1;  /* strictly increasing */
    int idx = c->count++;
    c->keys[idx].t = t; c->keys[idx].v = v;
    c->keys[idx].tan_in = tin; c->keys[idx].tan_out = tout;
    return idx;
}

static float curve_eval(const mal_curve* c, float t) {
    if (c->count == 0) return 0.0f;
    if (t <= c->keys[0].t) return c->keys[0].v;
    if (t >= c->keys[c->count - 1].t) return c->keys[c->count - 1].v;
    for (int i = 1; i < c->count; i++) {
        const mal_curve_key* a = &c->keys[i - 1];
        const mal_curve_key* b = &c->keys[i];
        if (t <= b->t) {
            float dt = b->t - a->t;
            float u = (t - a->t) / dt;
            float u2 = u * u, u3 = u2 * u;
            float h00 = 2 * u3 - 3 * u2 + 1;
            float h10 = u3 - 2 * u2 + u;
            float h01 = -2 * u3 + 3 * u2;
            float h11 = u3 - u2;
            return h00 * a->v + h10 * dt * a->tan_out + h01 * b->v + h11 * dt * b->tan_in;
        }
    }
    return c->keys[c->count - 1].v;
}

/* derivative for continuity checks */
#ifndef MAL_NO_MAIN
static float curve_deriv(const mal_curve* c, float t) {
    float h = 1e-3f;
    float a = curve_eval(c, t - h);
    float b = curve_eval(c, t + h);
    return (b - a) / (2 * h);
}
#endif

/* ---- unified-entry wrappers ------------------------------------------ */

void mal_curve_init(mal_curve* c) { curve_init(c); }
int mal_curve_capacity(void) { return MAL_CURVE_MAX_KEYS; }

mal_status mal_curve_add(mal_curve* c, const mal_curve_add_req* req,
                         mal_curve_add_resp* resp, char* why, int whycap) {
    if (!c) return mal_fail(why, whycap, MAL_ERR_PARAM,
                            "curve: null instance");
    if (!req) return mal_fail(why, whycap, MAL_ERR_PARAM,
                              "curve: null request payload");
    if (!mal_finite_f(req->t) || !mal_finite_f(req->v) ||
        !mal_finite_f(req->tan_in) || !mal_finite_f(req->tan_out))
        return mal_fail(why, whycap, MAL_ERR_PARAM,
                        "curve: non-finite keyframe value");
    if (c->count >= MAL_CURVE_MAX_KEYS)
        return mal_fail(why, whycap, MAL_ERR_STATE,
                        "curve: capacity %d reached", MAL_CURVE_MAX_KEYS);
    if (c->count > 0 && req->t <= c->keys[c->count - 1].t)
        return mal_fail(why, whycap, MAL_ERR_PARAM,
                        "curve: keyframe time %g must be strictly increasing",
                        (double)req->t);
    int idx = curve_add(c, req->t, req->v, req->tan_in, req->tan_out);
    if (resp) resp->index = idx;
    return MAL_OK;
}

mal_status mal_curve_eval(mal_curve* c, const mal_curve_eval_req* req,
                          mal_curve_eval_resp* resp, char* why, int whycap) {
    if (!c) return mal_fail(why, whycap, MAL_ERR_PARAM,
                            "curve: null instance");
    if (!req) return mal_fail(why, whycap, MAL_ERR_PARAM,
                              "curve: null request payload");
    if (!mal_finite_f(req->t))
        return mal_fail(why, whycap, MAL_ERR_PARAM,
                        "curve: non-finite sample time");
    if (c->count == 0)
        return mal_fail(why, whycap, MAL_ERR_STATE,
                        "curve: cannot eval an empty curve");
    float v = curve_eval(c, req->t);
    if (resp) resp->v = v;
    return MAL_OK;
}

/* dispatch entry used by the unified mal_solve(). Non-static in a
 * consolidated build (MAL_NO_MAIN) so pose/ can register it. */
MAL_DISPATCH_FN mal_status curve_dispatch(int op, const void* req,
                                          void* resp,
                                          char* why, int whycap) {
    if (op == MAL_OP_CURVE_EVAL) {
        const mal_curve_eval_wrap* q = (const mal_curve_eval_wrap*)req;
        if (!q->c) return mal_fail(why, whycap, MAL_ERR_PARAM,
                                   "curve: null instance");
        return mal_curve_eval(q->c, &q->in,
                              (mal_curve_eval_resp*)resp, why, whycap);
    }
    if (op == MAL_OP_CURVE_ADD) {
        const mal_curve_add_wrap* q = (const mal_curve_add_wrap*)req;
        if (!q->c) return mal_fail(why, whycap, MAL_ERR_PARAM,
                                   "curve: null instance");
        return mal_curve_add(q->c, &q->in,
                             (mal_curve_add_resp*)resp, why, whycap);
    }
    return mal_fail(why, whycap, MAL_ERR_UNSUPPORTED,
                    "curve: no op %d", op);
}

#ifndef MAL_NO_MAIN
const mal_dispatch_entry mal_dispatch_table[] = {
    { MAL_MODULE_CURVE, curve_dispatch },
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
static int u_curve_add(mal_curve* c, float t, float v, float tin, float tout) {
    mal_curve_add_wrap w; w.c = c; w.in.t = t; w.in.v = v;
    w.in.tan_in = tin; w.in.tan_out = tout;
    mal_curve_add_resp out;
    mal_request r = { MAL_MODULE_CURVE, MAL_OP_CURVE_ADD, &w, &out };
    if (mal_solve(&r, 0, 0) != MAL_OK) return -1;
    return out.index;
}
static float u_curve_eval(mal_curve* c, float t) {
    mal_curve_eval_wrap w; w.c = c; w.in.t = t;
    mal_curve_eval_resp out;
    mal_request r = { MAL_MODULE_CURVE, MAL_OP_CURVE_EVAL, &w, &out };
    if (mal_solve(&r, 0, 0) != MAL_OK) return 0.0f;
    return out.v;
}
#define CURVE_ADD(c,t,v,ti,to) u_curve_add((c),(t),(v),(ti),(to))
#define CURVE_EVAL(c,t)        u_curve_eval((c),(t))
#else
#define CURVE_ADD(c,t,v,ti,to) curve_add((c),(t),(v),(ti),(to))
#define CURVE_EVAL(c,t)        curve_eval((c),(t))
#endif

int main(void) {
    mal_curve c;
    curve_init(&c);
    check("add_first", CURVE_ADD(&c, 0.0f, 0.0f, 0.0f, 1.0f) == 0);
    check("add_second", CURVE_ADD(&c, 1.0f, 1.0f, 1.0f, 0.0f) == 1);
    check("clamp_left", fabsf(CURVE_EVAL(&c, -1.0f)) < 1e-6f);
    check("clamp_right", fabsf(CURVE_EVAL(&c, 2.0f) - 1.0f) < 1e-6f);
    check("mid_monotonic", CURVE_EVAL(&c, 0.5f) > 0.0f && CURVE_EVAL(&c, 0.5f) < 1.0f);
    check("endpoints_exact", fabsf(CURVE_EVAL(&c, 0.0f)) < 1e-6f && fabsf(CURVE_EVAL(&c, 1.0f) - 1.0f) < 1e-6f);

    mal_curve flat;
    curve_init(&flat);
    CURVE_ADD(&flat, 0.0f, 0.5f, 0.0f, 0.0f);
    CURVE_ADD(&flat, 1.0f, 0.5f, 0.0f, 0.0f);
    check("flat_curve", fabsf(CURVE_EVAL(&flat, 0.5f) - 0.5f) < 1e-5f);
    check("flat_deriv_zero", fabsf(curve_deriv(&flat, 0.5f)) < 1e-3f);

    check("reject_unsorted", CURVE_ADD(&flat, 0.5f, 1.0f, 0, 0) == -1);

    mal_curve full;
    curve_init(&full);
    int rc = 0;
    for (int i = 0; i <= MAL_CURVE_MAX_KEYS; i++) rc = CURVE_ADD(&full, (float)i, (float)i, 0, 0);
    check("capacity", rc == -1 && full.count == MAL_CURVE_MAX_KEYS);

#if MAL_UNIFIED
    /* illegal parameters must be rejected with the module named */
    char why[128];
    mal_curve_add_wrap w; w.c = &flat;
    w.in.t = 0.5f; w.in.v = 1.0f; w.in.tan_in = 0; w.in.tan_out = 0;
    mal_request bad = { MAL_MODULE_CURVE, MAL_OP_CURVE_ADD, &w, 0 };
    check("unified_unsorted",
          mal_solve(&bad, why, sizeof why) == MAL_ERR_PARAM &&
          strstr(why, "curve:") == why);

    mal_curve_eval_wrap ew; ew.c = 0; ew.in.t = 0.0f;
    bad.req = &ew; bad.op = MAL_OP_CURVE_EVAL;
    check("unified_null_instance",
          mal_solve(&bad, why, sizeof why) == MAL_ERR_PARAM &&
          strstr(why, "curve:") == why);
#endif

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
#endif
