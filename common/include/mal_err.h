#ifndef MAL_ERR_H
#define MAL_ERR_H

#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include "solver.h"

/*
 * Write a module-prefixed error message (e.g. "ik: bone length must be
 * positive"). why/whycap may be NULL/0: callers that only want the status
 * code are allowed to skip the message.
 */
static inline mal_status mal_fail(char* why, int whycap, mal_status st,
                                  const char* fmt, ...) {
    if (why && whycap > 0) {
        va_list ap;
        va_start(ap, fmt);
        vsnprintf(why, (size_t)whycap, fmt, ap);
        va_end(ap);
    }
    return st;
}

static inline int mal_finite_f(float v) {
    /* isfinite is a macro that works on float without promoting issues. */
    return isfinite(v);
}

#endif
