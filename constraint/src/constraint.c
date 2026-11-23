
/* constraint: verlet particle chain with distance constraints, iteration
   count and damping; used as the toy model for cloth/rope style solvers.
 *
 * MAL_UNIFIED=1 (default): init/step go through mal_solve().
 * MAL_UNIFIED=0: the original static chain_init/chain_step are called
 * directly. The math below is shared by both paths.
 */
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "solver.h"
#include "constraint.h"
#include "mal_err.h"

#ifndef MAL_UNIFIED
#define MAL_UNIFIED 1
#endif

#ifdef MAL_NO_MAIN
#define MAL_DISPATCH_FN
#else
#define MAL_DISPATCH_FN static
#endif

static void chain_init(mal_chain* c, int n, float seg, float x0, float y0) {
    memset(c, 0, sizeof(*c));
    c->n = n;
    c->iterations = 8;
    c->damping = 0.99f;
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

static void chain_step(mal_chain* c, float dt, float gx, float gy) {
    for (int i = 0; i < c->n; i++) {
        mal_particle* p = &c->p[i];
        if (p->pinned) continue;
        float vx = (p->x - p->px) * c->damping;
        float vy = (p->y - p->py) * c->damping;
        p->px = p->x; p->py = p->y;
        p->x += vx + gx * dt * dt;
        p->y += vy + gy * dt * dt;
    }
    for (int it = 0; it < c->iterations; it++) {
        for (int i = 0; i < c->n - 1; i++) {
            mal_particle* a = &c->p[i];
            mal_particle* b = &c->p[i + 1];
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

static float chain_seg_error(const mal_chain* c) {
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

/* ---- unified-entry wrappers ------------------------------------------ */

void mal_chain_get(const mal_chain* c, int i, mal_v2* pos, float* inv_mass,
                   int* pinned) {
    if (pos) { pos->x = c->p[i].x; pos->y = c->p[i].y; }
    if (inv_mass) *inv_mass = c->p[i].inv_mass;
    if (pinned) *pinned = c->p[i].pinned;
}

int mal_chain_count(const mal_chain* c) { return c->n; }
float mal_chain_rest(const mal_chain* c, int i) { return c->rest[i]; }
void mal_chain_set_iterations(mal_chain* c, int iterations) {
    c->iterations = iterations;
}
float mal_chain_seg_error(const mal_chain* c) { return chain_seg_error(c); }

mal_status mal_chain_init(mal_chain* c, const mal_chain_init_req* req,
                          char* why, int whycap) {
    if (!c) return mal_fail(why, whycap, MAL_ERR_PARAM,
                            "constraint: null instance");
    if (!req) return mal_fail(why, whycap, MAL_ERR_PARAM,
                              "constraint: null request payload");
    if (req->n < 2 || req->n > MAL_CHAIN_MAX)
        return mal_fail(why, whycap, MAL_ERR_PARAM,
                        "constraint: particle count %d out of range [2,%d]",
                        req->n, MAL_CHAIN_MAX);
    if (!mal_finite_f(req->seg) || req->seg <= 0.0f)
        return mal_fail(why, whycap, MAL_ERR_PARAM,
                        "constraint: segment length must be positive (got %g)",
                        (double)req->seg);
    if (!mal_finite_f(req->x0) || !mal_finite_f(req->y0))
        return mal_fail(why, whycap, MAL_ERR_PARAM,
                        "constraint: non-finite anchor position");
    chain_init(c, req->n, req->seg, req->x0, req->y0);
    return MAL_OK;
}

mal_status mal_chain_step(mal_chain* c, const mal_chain_step_req* req,
                          char* why, int whycap) {
    if (!c) return mal_fail(why, whycap, MAL_ERR_PARAM,
                            "constraint: null instance");
    if (!req) return mal_fail(why, whycap, MAL_ERR_PARAM,
                              "constraint: null request payload");
    if (c->n < 2)
        return mal_fail(why, whycap, MAL_ERR_STATE,
                        "constraint: chain not initialized");
    if (!mal_finite_f(req->dt) || req->dt <= 0.0f)
        return mal_fail(why, whycap, MAL_ERR_PARAM,
                        "constraint: dt must be positive (got %g)",
                        (double)req->dt);
    if (!mal_finite_f(req->gx) || !mal_finite_f(req->gy))
        return mal_fail(why, whycap, MAL_ERR_PARAM,
                        "constraint: non-finite gravity");
    if (c->iterations < 0)
        return mal_fail(why, whycap, MAL_ERR_PARAM,
                        "constraint: iterations must be >= 0 (got %d)",
                        c->iterations);
    chain_step(c, req->dt, req->gx, req->gy);
    return MAL_OK;
}

/* dispatch entry used by the unified mal_solve(). Non-static in a
 * consolidated build (MAL_NO_MAIN) so pose/ can register it. */
MAL_DISPATCH_FN mal_status chain_dispatch(int op, const void* req,
                                          void* resp,
                                          char* why, int whycap) {
    (void)resp;
    if (op == MAL_OP_CHAIN_STEP) {
        const mal_chain_step_wrap* q = (const mal_chain_step_wrap*)req;
        if (!q->ch) return mal_fail(why, whycap, MAL_ERR_PARAM,
                                    "constraint: null instance");
        return mal_chain_step(q->ch, &q->in, why, whycap);
    }
    if (op == MAL_OP_CHAIN_INIT) {
        const mal_chain_init_wrap* q = (const mal_chain_init_wrap*)req;
        if (!q->ch) return mal_fail(why, whycap, MAL_ERR_PARAM,
                                    "constraint: null instance");
        return mal_chain_init(q->ch, &q->in, why, whycap);
    }
    return mal_fail(why, whycap, MAL_ERR_UNSUPPORTED,
                    "constraint: no op %d", op);
}

#ifndef MAL_NO_MAIN
const mal_dispatch_entry mal_dispatch_table[] = {
    { MAL_MODULE_CONSTRAINT, chain_dispatch },
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
static void u_chain_init(mal_chain* c, int n, float seg, float x0, float y0) {
    mal_chain_init_wrap w; w.ch = c;
    w.in.n = n; w.in.seg = seg; w.in.x0 = x0; w.in.y0 = y0;
    mal_request r = { MAL_MODULE_CONSTRAINT, MAL_OP_CHAIN_INIT, &w, 0 };
    (void)mal_solve(&r, 0, 0);
}
static void u_chain_step(mal_chain* c, float dt, float gx, float gy) {
    mal_chain_step_wrap w; w.ch = c; w.in.dt = dt; w.in.gx = gx; w.in.gy = gy;
    mal_request r = { MAL_MODULE_CONSTRAINT, MAL_OP_CHAIN_STEP, &w, 0 };
    (void)mal_solve(&r, 0, 0);
}
#define CHAIN_INIT(c,n,seg,x0,y0) u_chain_init((c),(n),(seg),(x0),(y0))
#define CHAIN_STEP(c,dt,gx,gy)    u_chain_step((c),(dt),(gx),(gy))
#else
#define CHAIN_INIT(c,n,seg,x0,y0) chain_init((c),(n),(seg),(x0),(y0))
#define CHAIN_STEP(c,dt,gx,gy)    chain_step((c),(dt),(gx),(gy))
#endif

int main(void) {
    mal_chain c;
    CHAIN_INIT(&c, 8, 1.0f, 0.0f, 0.0f);
    check("init_len", c.n == 8 && fabsf(chain_seg_error(&c)) < 1e-5f);

    /* let it fall for a while */
    for (int i = 0; i < 200; i++) CHAIN_STEP(&c, 0.016f, 0.0f, -9.8f);
    check("stays_in_stretch", chain_seg_error(&c) < 0.02f);
    check("pinned_anchor", fabsf(c.p[0].x) < 1e-6f && fabsf(c.p[0].y) < 1e-6f);
    check("falls_down", c.p[7].y < -1.0f);

    /* more iterations must not increase error */
    mal_chain a, b;
    CHAIN_INIT(&a, 8, 1.0f, 0.0f, 0.0f);
    CHAIN_INIT(&b, 8, 1.0f, 0.0f, 0.0f);
    a.iterations = 1; b.iterations = 16;
    for (int i = 0; i < 100; i++) { CHAIN_STEP(&a, 0.016f, 0.0f, -9.8f); CHAIN_STEP(&b, 0.016f, 0.0f, -9.8f); }
    check("more_iters_better", chain_seg_error(&b) <= chain_seg_error(&a) + 1e-6f);

    mal_chain cap;
    CHAIN_INIT(&cap, 16, 0.5f, 0.0f, 0.0f);
    check("capacity_ok", cap.n == 16);

#if MAL_UNIFIED
    /* illegal parameters must be rejected with the module named */
    char why[128];
    mal_chain_init_wrap w; w.ch = &c;
    w.in.n = 0; w.in.seg = 1.0f; w.in.x0 = 0.0f; w.in.y0 = 0.0f;
    mal_request bad = { MAL_MODULE_CONSTRAINT, MAL_OP_CHAIN_INIT, &w, 0 };
    check("unified_bad_count",
          mal_solve(&bad, why, sizeof why) == MAL_ERR_PARAM &&
          strstr(why, "constraint:") == why);

    w.in.n = 8; w.in.seg = 0.0f;
    check("unified_bad_seg",
          mal_solve(&bad, why, sizeof why) == MAL_ERR_PARAM &&
          strstr(why, "constraint:") == why);

    mal_chain_step_wrap sw; sw.ch = &c;
    sw.in.dt = 0.0f; sw.in.gx = 0.0f; sw.in.gy = -9.8f;
    bad.req = &sw; bad.op = MAL_OP_CHAIN_STEP;
    check("unified_bad_dt",
          mal_solve(&bad, why, sizeof why) == MAL_ERR_PARAM &&
          strstr(why, "constraint:") == why);
#endif

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
#endif
