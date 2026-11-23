/*
 * bench_skeleton: pose snapshot and two-tier timing for the skeleton.
 *
 * Prints the same world pose under both call paths:
 *   legacy  -> direct static routines (MAL_UNIFIED=0 behavior)
 *   unified -> mal_solve() dispatch
 * followed by one-round and 10000-round timing. The math is shared, so
 * the two pose blocks must match to the last float.
 */
#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "solver.h"
#include "skeleton.h"

mal_status skeleton_dispatch(int, const void*, void*,
                                             char*, int);
const mal_dispatch_entry mal_dispatch_table[] = {
    { MAL_MODULE_SKELETON, skeleton_dispatch }, { (mal_module)0, 0 }
};

static double now_s(void) {
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

static void build(mal_skeleton* s) {
    mal_skeleton_init(s);
    mal_skel_add_wrap w; w.s = s;
    w.in.parent = -1; w.in.x = 0; w.in.y = 0; w.in.ang = 0;
    mal_request r0 = { MAL_MODULE_SKELETON, MAL_OP_SKELETON_ADD, &w, 0 };
    (void)mal_solve(&r0, 0, 0);
    w.in.parent = 0; w.in.x = 1; w.in.y = 0; w.in.ang = 0.6f;
    (void)mal_solve(&r0, 0, 0);
    w.in.parent = 1; w.in.x = 1; w.in.y = 0; w.in.ang = -0.9f;
    (void)mal_solve(&r0, 0, 0);
    w.in.parent = 2; w.in.x = 1; w.in.y = 0; w.in.ang = 0.3f;
    (void)mal_solve(&r0, 0, 0);
}

static void dump_pose(const char* tag, const mal_skeleton* s) {
    for (int i = 0; i < s->count; i++)
        printf("%s bone%d pos=(%.7f,%.7f) ang=%.7f\n",
               tag, i, s->world[i].x, s->world[i].y, s->world_ang[i]);
}

int main(void) {
    mal_skeleton s;
    build(&s);

    /* unified path */
    mal_skel_update_req uq = { &s };
    mal_request ru = { MAL_MODULE_SKELETON, MAL_OP_SKELETON_UPDATE, &uq, 0 };
    double t0 = now_s();
    (void)mal_solve(&ru, 0, 0);
    double t1 = now_s();
    dump_pose("unified", &s);

    int N = 10000;
    double t2 = now_s();
    for (int i = 0; i < N; i++) (void)mal_solve(&ru, 0, 0);
    double t3 = now_s();

    /* direct-wrapper path (same math, no dispatch indirection) */
    (void)mal_skeleton_update(&s, 0, 0);
    dump_pose("direct ", &s);

    printf("time_one_update_us=%.3f\n", (t1 - t0) * 1e6);
    printf("time_%d_updates_ms=%.3f\n", N, (t3 - t2) * 1e3);
    printf("time_per_update_us=%.4f\n", (t3 - t2) * 1e6 / N);
    return 0;
}
