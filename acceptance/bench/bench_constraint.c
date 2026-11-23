/*
 * bench_constraint: chain-pose snapshot and two-tier timing.
 * Compares the unified mal_solve() stepping path against the direct
 * wrapper over the identical Verlet math, and times one step plus a bulk
 * batch. The two paths must produce byte-identical particle positions.
 */
#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "solver.h"
#include "constraint.h"

mal_status chain_dispatch(int, const void*, void*,
                                          char*, int);
const mal_dispatch_entry mal_dispatch_table[] = {
    { MAL_MODULE_CONSTRAINT, chain_dispatch }, { (mal_module)0, 0 }
};

static double now_s(void) {
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

static void dump(const char* tag, const mal_chain* c) {
    for (int i = 0; i < c->n; i++)
        printf("%s p%d pos=(%.7f,%.7f)\n", tag, i, c->p[i].x, c->p[i].y);
    printf("%s seg_error=%.7f\n", tag, mal_chain_seg_error(c));
}

int main(void) {
    mal_chain cu, cd;
    mal_chain_init_wrap iu; iu.ch = &cu;
    iu.in.n = 8; iu.in.seg = 1.0f; iu.in.x0 = 0; iu.in.y0 = 0;
    mal_request ri = { MAL_MODULE_CONSTRAINT, MAL_OP_CHAIN_INIT, &iu, 0 };
    if (mal_solve(&ri, 0, 0) != MAL_OK) return 2;
    mal_chain_init_wrap id_; id_.ch = &cd; id_.in = iu.in;
    if (mal_chain_init(&cd, &id_.in, 0, 0) != MAL_OK) return 2;

    mal_chain_step_wrap sw; sw.ch = &cu;
    sw.in.dt = 0.016f; sw.in.gx = 0.0f; sw.in.gy = -9.8f;
    mal_request rs = { MAL_MODULE_CONSTRAINT, MAL_OP_CHAIN_STEP, &sw, 0 };

    double t0 = now_s();
    (void)mal_solve(&rs, 0, 0);
    double t1 = now_s();

    int N = 10000;
    double t2 = now_s();
    for (int i = 0; i < N; i++) (void)mal_solve(&rs, 0, 0);
    double t3 = now_s();

    /* replay the same N+1 steps directly to prove byte-identical poses */
    for (int i = 0; i < N + 1; i++)
        (void)mal_chain_step(&cd, &sw.in, 0, 0);

    dump("unified", &cu);
    dump("direct ", &cd);
    int same = memcmp(cu.p, cd.p, sizeof(cu.p)) == 0;
    printf("poses_byte_identical=%d\n", same);

    printf("time_one_step_us=%.3f\n", (t1 - t0) * 1e6);
    printf("time_%d_steps_ms=%.3f\n", N, (t3 - t2) * 1e3);
    printf("time_per_step_us=%.4f\n", (t3 - t2) * 1e6 / N);
    return same ? 0 : 1;
}
