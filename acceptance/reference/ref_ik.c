/* ref_ik: original pre-change IK math, same scenario as bench_ik. */
#define _POSIX_C_SOURCE 199309L
#include <stdio.h>

#define main ik_orig_main_unused
#include "orig/ik.c"
#undef main

int main(void) {
    V2 root = {0, 0}, pole = {0, 1};
    V2 target = {1.2f, 0.4f};
    IkResult r = solve_two_bone(root, target, 1.0f, 1.0f, pole);
    printf("ref ok=%d joint=(%.7f,%.7f) a1=%.7f a2=%.7f\n",
           r.ok, r.joint.x, r.joint.y, r.a1, r.a2);
    return 0;
}
