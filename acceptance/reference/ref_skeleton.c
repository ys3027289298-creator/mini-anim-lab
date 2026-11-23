/*
 * ref_skeleton: drives the ORIGINAL pre-change skeleton.c static routines
 * through the same scenario as bench_skeleton, emitting reference pose
 * values. The original main() is renamed away via macro so only its
 * functions are compiled here.
 */
#define _POSIX_C_SOURCE 199309L
#include <stdio.h>

#define main skeleton_orig_main_unused
#include "orig/skeleton.c"
#undef main

int main(void) {
    Skeleton s;
    skel_init(&s);
    int b0 = skel_add(&s, -1, 0, 0, 0);
    int b1 = skel_add(&s, b0, 1, 0, 0.6f);
    int b2 = skel_add(&s, b1, 1, 0, -0.9f);
    int b3 = skel_add(&s, b2, 1, 0, 0.3f);
    (void)b3;
    skel_update(&s);
    for (int i = 0; i < s.count; i++)
        printf("ref bone%d pos=(%.7f,%.7f) ang=%.7f\n",
               i, s.world[i].x, s.world[i].y, s.world[i].ang);
    return 0;
}
