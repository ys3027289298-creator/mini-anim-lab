#ifndef CURVE_PUB_H
#define CURVE_PUB_H

#include "solver.h"

#define MAL_CURVE_MAX_KEYS 16

typedef struct { float t, v, tan_in, tan_out; } mal_curve_key;

struct Curve {
    mal_curve_key keys[MAL_CURVE_MAX_KEYS];
    int count;
};

struct mal_curve_add_req {
    float t, v, tan_in, tan_out;
};
struct mal_curve_add_resp {
    int index;
};

struct mal_curve_eval_req {
    float t;
};
struct mal_curve_eval_resp {
    float v;
};

/* Stateful ops bundle the instance with their inputs. */
struct mal_curve_eval_wrap { mal_curve* c; struct mal_curve_eval_req in; };
struct mal_curve_add_wrap { mal_curve* c; struct mal_curve_add_req in; };
void mal_curve_init(mal_curve* c);
int mal_curve_capacity(void);
mal_status mal_curve_add(mal_curve* c, const mal_curve_add_req* req,
                         mal_curve_add_resp* resp, char* why, int whycap);
mal_status mal_curve_eval(mal_curve* c, const mal_curve_eval_req* req,
                          mal_curve_eval_resp* resp, char* why, int whycap);

#endif
