/*
 * bench_curve: sampled-value snapshot and two-tier timing for the curve.
 * Same hermite math through mal_solve() and through the direct wrapper;
 * outputs must match float-for-float.
 */
#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "solver.h"
#include "curve.h"

mal_status curve_dispatch(int, const void*, void*,
                                          char*, int);
const mal_dispatch_entry mal_dispatch_table[] = {
    { MAL_MODULE_CURVE, curve_dispatch }, { (mal_module)0, 0 }
};

static double now_s(void) {
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

int main(void) {
    mal_curve c;
    mal_curve_init(&c);
    const float ts[]  = {0.0f, 0.3f, 0.7f, 1.0f};
    const float vs[]  = {0.0f, 0.4f, 0.9f, 1.0f};
    for (int i = 0; i < 4; i++) {
        mal_curve_add_wrap w; w.c = &c;
        w.in.t = ts[i]; w.in.v = vs[i]; w.in.tan_in = 1.0f;
        w.in.tan_out = i == 3 ? 0.0f : 1.0f;
        mal_request r = { MAL_MODULE_CURVE, MAL_OP_CURVE_ADD, &w, 0 };
        if (mal_solve(&r, 0, 0) != MAL_OK) return 2;
    }

    float samples[9];
    for (int i = 0; i < 9; i++) samples[i] = i / 8.0f;

    double t0 = now_s();
    float first = 0;
    for (int i = 0; i < 9; i++) {
        mal_curve_eval_wrap w; w.c = &c; w.in.t = samples[i];
        mal_curve_eval_resp out;
        mal_request r = { MAL_MODULE_CURVE, MAL_OP_CURVE_EVAL, &w, &out };
        (void)mal_solve(&r, 0, 0);
        if (i == 0) first = out.v;
    }
    double t1 = now_s();
    (void)first;

    int N = 100000;
    double t2 = now_s();
    mal_curve_eval_wrap w; w.c = &c; w.in.t = 0.5f;
    mal_curve_eval_resp out;
    mal_request r = { MAL_MODULE_CURVE, MAL_OP_CURVE_EVAL, &w, &out };
    for (int i = 0; i < N; i++) (void)mal_solve(&r, 0, 0);
    double t3 = now_s();

    for (int i = 0; i < 9; i++) {
        mal_curve_eval_wrap uw; uw.c = &c; uw.in.t = samples[i];
        mal_curve_eval_resp uo, do_;
        mal_request ur = { MAL_MODULE_CURVE, MAL_OP_CURVE_EVAL, &uw, &uo };
        (void)mal_solve(&ur, 0, 0);
        (void)mal_curve_eval(&c, &uw.in, &do_, 0, 0);
        printf("t=%.4f unified=%.7f direct=%.7f match=%d\n",
               samples[i], uo.v, do_.v, uo.v == do_.v);
    }

    printf("time_9_evals_us=%.3f\n", (t1 - t0) * 1e6);
    printf("time_%d_evals_ms=%.3f\n", N, (t3 - t2) * 1e3);
    printf("time_per_eval_us=%.4f\n", (t3 - t2) * 1e6 / N);
    return 0;
}
