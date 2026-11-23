/*
 * bench_ik: joint-pose snapshot and two-tier timing for analytic IK.
 * Compares the unified mal_solve() path against the direct wrapper (which
 * runs the identical solve_two_bone math) and times one solve plus a bulk
 * batch.
 */
#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "solver.h"
#include "ik.h"

mal_status ik_dispatch(int, const void*, void*, char*, int);
const mal_dispatch_entry mal_dispatch_table[] = {
    { MAL_MODULE_IK, ik_dispatch }, { (mal_module)0, 0 }
};

static double now_s(void) {
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

static void dump(const char* tag, const mal_ik_resp* r) {
    printf("%s ok=%d joint=(%.7f,%.7f) a1=%.7f a2=%.7f\n",
           tag, r->ok, r->joint.x, r->joint.y, r->a1, r->a2);
}

int main(void) {
    mal_v2 root = {0,0}, pole = {0,1};
    mal_ik_req req = { root, (mal_v2){1.2f, 0.4f}, pole, 1.0f, 1.0f };
    mal_ik_resp ru_out, rd_out;
    mal_request r = { MAL_MODULE_IK, MAL_OP_IK_SOLVE, &req, &ru_out };

    double t0 = now_s();
    mal_status st = mal_solve(&r, 0, 0);
    double t1 = now_s();
    (void)st;
    dump("unified", &ru_out);

    int N = 100000;
    double t2 = now_s();
    for (int i = 0; i < N; i++) (void)mal_solve(&r, 0, 0);
    double t3 = now_s();

    (void)mal_ik_solve(&req, &rd_out, 0, 0);
    dump("direct ", &rd_out);

    printf("time_one_solve_us=%.3f\n", (t1 - t0) * 1e6);
    printf("time_%d_solves_ms=%.3f\n", N, (t3 - t2) * 1e3);
    printf("time_per_solve_us=%.4f\n", (t3 - t2) * 1e6 / N);
    return 0;
}
