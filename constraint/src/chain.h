/* chain.h: verlet particle chain with pbd distance constraints.
   Solver-only API, split out of the original single-file self-test so the
   conservation auditor (../audit) links the exact same solver code. */
#ifndef MINI_CHAIN_H
#define MINI_CHAIN_H

#define MAX_P 16

typedef struct { float x, y, px, py; float inv_mass; int pinned; } Particle;

/* Optional read-only observation hooks. Both default to NULL. The solver
   result never depends on whether a hook is attached: a hook only observes
   values the solver already computed and must not modify the chain. */
typedef void (*ChainProjectHook)(void* user, int iter, int ci,
                                 float dxa, float dya,
                                 float dxb, float dyb);
typedef void (*ChainStageHook)(void* user);

typedef struct {
    ChainProjectHook project;   /* called after every constraint projection */
    ChainStageHook integrated;  /* called once after the verlet integration pass */
    void* user;
} ChainObserver;

typedef struct {
    Particle p[MAX_P];
    int n;
    float rest[MAX_P];      /* rest length between i and i+1 */
    int iterations;
    float damping;
    ChainObserver obs;      /* zeroed by chain_init => no observation */
} Chain;

void chain_init(Chain* c, int n, float seg, float x0, float y0);
void chain_step(Chain* c, float dt, float gx, float gy);
float chain_seg_error(const Chain* c);

#endif
