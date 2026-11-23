/* ref_constraint: original pre-change Verlet chain, same init and the
 * same 10001-step replay count as bench_constraint (one timed step plus
 * the 10000 bulk steps). */
#define _POSIX_C_SOURCE 199309L
#include <stdio.h>

#define main chain_orig_main_unused
#include "orig/constraint.c"
#undef main

int main(void) {
    Chain c;
    chain_init(&c, 8, 1.0f, 0.0f, 0.0f);
    for (int i = 0; i < 10001; i++)
        chain_step(&c, 0.016f, 0.0f, -9.8f);
    for (int i = 0; i < c.n; i++)
        printf("ref p%d pos=(%.7f,%.7f)\n", i, c.p[i].x, c.p[i].y);
    printf("ref seg_error=%.7f\n", chain_seg_error(&c));
    return 0;
}
