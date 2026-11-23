#ifndef CONSTRAINT_PUB_H
#define CONSTRAINT_PUB_H

#include "solver.h"
#include "mal_vec.h"

#define MAL_CHAIN_MAX 16

typedef struct {
    float x, y, px, py, inv_mass;
    int pinned;
} mal_particle;

struct Chain {
    mal_particle p[MAL_CHAIN_MAX];
    int n;
    float rest[MAL_CHAIN_MAX];   /* rest length between i and i+1 */
    int iterations;
    float damping;
};

struct mal_chain_init_req {
    int n;
    float seg, x0, y0;
};

struct mal_chain_step_req {
    float dt, gx, gy;
};

/* Stateful ops bundle the instance with their inputs. */
struct mal_chain_init_wrap { mal_chain* ch; struct mal_chain_init_req in; };
struct mal_chain_step_wrap { mal_chain* ch; struct mal_chain_step_req in; };

void mal_chain_get(const mal_chain* c, int i, mal_v2* pos, float* inv_mass,
                   int* pinned);
int mal_chain_count(const mal_chain* c);
float mal_chain_rest(const mal_chain* c, int i);
void mal_chain_set_iterations(mal_chain* c, int iterations);
float mal_chain_seg_error(const mal_chain* c);

mal_status mal_chain_init(mal_chain* c, const mal_chain_init_req* req,
                          char* why, int whycap);
mal_status mal_chain_step(mal_chain* c, const mal_chain_step_req* req,
                          char* why, int whycap);

#endif
