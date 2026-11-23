
/* constraint: verlet particle chain with distance constraints, iteration
   count and damping; used as the toy model for cloth/rope style solvers. */
#include <stdio.h>
#include <math.h>
#include <string.h>

#ifndef MAL_MODULE
#define MAL_MODULE "constraint"
#endif
#define MAL_TIME_ITERS 40000L
#include "../../common/harness.h"

#define MAX_P 16

typedef struct { float x, y, px, py; float inv_mass; int pinned; } Particle;

typedef struct {
    Particle p[MAX_P];
    int n;
    float rest[MAX_P];      /* rest length between i and i+1 */
    int iterations;
    float damping;
} Chain;

static void chain_init(Chain* c, int n, float seg, float x0, float y0) {
    memset(c, 0, sizeof(*c));
    c->n = n;
    c->iterations = 8;
    c->damping = 0.99f;
    for (int i = 0; i < n; i++) {
        c->p[i].x = x0 + seg * i; c->p[i].y = y0;
        c->p[i].px = c->p[i].x; c->p[i].py = c->p[i].y;
        c->p[i].inv_mass = 1.0f;
        c->p[i].pinned = 0;
        c->rest[i] = seg;
    }
    c->p[0].pinned = 1;
    c->p[0].inv_mass = 0.0f;
}

static void chain_step(Chain* c, float dt, float gx, float gy) {
    for (int i = 0; i < c->n; i++) {
        Particle* p = &c->p[i];
        if (p->pinned) continue;
        float vx = (p->x - p->px) * c->damping;
        float vy = (p->y - p->py) * c->damping;
        p->px = p->x; p->py = p->y;
        p->x += vx + gx * dt * dt;
        p->y += vy + gy * dt * dt;
    }
    for (int it = 0; it < c->iterations; it++) {
        for (int i = 0; i < c->n - 1; i++) {
            Particle* a = &c->p[i];
            Particle* b = &c->p[i + 1];
            float dx = b->x - a->x, dy = b->y - a->y;
            float d = sqrtf(dx * dx + dy * dy);
            if (d < 1e-6f) continue;
            float diff = (d - c->rest[i]) / d;
            float wa = a->inv_mass, wb = b->inv_mass;
            float wsum = wa + wb;
            if (wsum <= 0) continue;
            a->x += dx * diff * (wa / wsum);
            a->y += dy * diff * (wa / wsum);
            b->x -= dx * diff * (wb / wsum);
            b->y -= dy * diff * (wb / wsum);
        }
    }
}

static float chain_seg_error(const Chain* c) {
    float worst = 0.0f;
    for (int i = 0; i < c->n - 1; i++) {
        float dx = c->p[i + 1].x - c->p[i].x;
        float dy = c->p[i + 1].y - c->p[i].y;
        float d = sqrtf(dx * dx + dy * dy);
        float e = fabsf(d - c->rest[i]);
        if (e > worst) worst = e;
    }
    return worst;
}

/* the settled chain from the falls_down check, reused by the evidence paths */
MAL_UNUSED static void constraint_settle(Chain* c) {
    chain_init(c, 8, 1.0f, 0.0f, 0.0f);
    for (int i = 0; i < 200; i++) chain_step(c, 0.016f, 0.0f, -9.8f);
}

static void run_checks(void) {
    Chain c;
    chain_init(&c, 8, 1.0f, 0.0f, 0.0f);
    mal_check("init_len", c.n == 8 && fabsf(chain_seg_error(&c)) < 1e-5f);

    /* let it fall for a while */
    for (int i = 0; i < 200; i++) chain_step(&c, 0.016f, 0.0f, -9.8f);
    mal_check("stays_in_stretch", chain_seg_error(&c) < 0.02f);
    mal_check("pinned_anchor", fabsf(c.p[0].x) < 1e-6f && fabsf(c.p[0].y) < 1e-6f);
    mal_check("falls_down", c.p[7].y < -1.0f);

    /* more iterations must not increase error */
    Chain a, b;
    chain_init(&a, 8, 1.0f, 0.0f, 0.0f);
    chain_init(&b, 8, 1.0f, 0.0f, 0.0f);
    a.iterations = 1; b.iterations = 16;
    for (int i = 0; i < 100; i++) { chain_step(&a, 0.016f, 0.0f, -9.8f); chain_step(&b, 0.016f, 0.0f, -9.8f); }
    mal_check("more_iters_better", chain_seg_error(&b) <= chain_seg_error(&a) + 1e-6f);

    Chain cap;
    chain_init(&cap, 16, 0.5f, 0.0f, 0.0f);
    mal_check("capacity_ok", cap.n == 16);
}

int main(void) {
    run_checks();

#ifdef MAL_EVIDENCE_POSE
    if (mal_failures == 0) {
        Chain c;
        constraint_settle(&c);
        mal_pose("n", c.n);
        for (int i = 0; i < c.n; i++) {
            char k[16];
            snprintf(k, sizeof k, "p%d_x", i); mal_pose(k, c.p[i].x);
            snprintf(k, sizeof k, "p%d_y", i); mal_pose(k, c.p[i].y);
        }
        mal_pose("seg_error", chain_seg_error(&c));
        mal_emit_pose();
    }
#endif

#ifdef MAL_EVIDENCE_IMAGE
    if (mal_failures == 0) {
        MalImage im;
        Chain c;
        constraint_settle(&c);
        mal_img_clear(&im);
        const double scale = 24.0, ox = 64.0, oy = 20.0;
        for (int i = 0; i < c.n - 1; i++)
            mal_img_line(&im, ox + c.p[i].x * scale, oy - c.p[i].y * scale,
                         ox + c.p[i + 1].x * scale, oy - c.p[i + 1].y * scale,
                         120, 235, 150);
        mal_emit_png(&im, MAL_MODULE);
    }
#endif

#ifdef MAL_EVIDENCE_TIMING
    if (mal_failures == 0) {
        Chain c;
        chain_init(&c, 8, 1.0f, 0.0f, 0.0f);
        double acc = 0.0;
        const long settle = MAL_TIME_ITERS / 200 > 0 ? MAL_TIME_ITERS / 200 : 1;
        double t0 = mal_now();
        for (long it = 0; it < MAL_TIME_ITERS; it++) {
            chain_step(&c, 0.016f, 0.0f, -9.8f);
            acc += chain_seg_error(&c) + c.p[7].y;
            if ((it % settle) == 0) chain_init(&c, 8, 1.0f, 0.0f, 0.0f);
        }
        double t1 = mal_now();
        if (!isfinite(acc)) mal_check("time_finite", 0);
        mal_emit_time(MAL_TIME_TIER, t1 - t0, MAL_TIME_ITERS);
    }
#endif

    printf("all_ok=%d\n", mal_failures == 0);
    return mal_failures == 0 ? 0 : 1;
}
