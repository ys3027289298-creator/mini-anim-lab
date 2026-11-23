/* ref_curve: original pre-change curve math, same keys/samples as
 * bench_curve. */
#define _POSIX_C_SOURCE 199309L
#include <stdio.h>

#define main curve_orig_main_unused
#include "orig/curve.c"
#undef main

int main(void) {
    Curve c;
    curve_init(&c);
    float ts[] = {0.0f, 0.3f, 0.7f, 1.0f};
    float vs[] = {0.0f, 0.4f, 0.9f, 1.0f};
    for (int i = 0; i < 4; i++)
        curve_add(&c, ts[i], vs[i], 1.0f, i == 3 ? 0.0f : 1.0f);
    for (int i = 0; i < 9; i++) {
        float t = i / 8.0f;
        printf("ref t=%.4f v=%.7f\n", t, curve_eval(&c, t));
    }
    return 0;
}
