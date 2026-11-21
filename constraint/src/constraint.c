
/* constraint: verlet particle chain with distance constraints, iteration
   count and damping; used as the toy model for cloth/rope style solvers. */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define MAX_P 16

typedef struct { float x, y, px, py; float inv_mass; int pinned; } Particle;

typedef struct {
    Particle p[MAX_P];
    int n;
    float rest[MAX_P];      /* rest length between i and i+1 */
    int iterations;
    float damping;
} Chain;

/* frame-stepping policy: split a (possibly large) frame dt into substeps,
   clamp per-substep displacement, and enforce a solver iteration floor. */
typedef struct {
    float target_dt;    /* substep size we aim for (s) */
    int   max_substeps; /* hard cap, avoids spiral of death; <=0 means no cap */
    float max_disp;     /* max per-substep displacement per particle; <=0 disables */
    int   min_iterations; /* solver iteration floor */
} StepPolicy;

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

static void chain_integrate(Chain* c, float dt, float gx, float gy, float max_disp) {
    /* damping is defined per 1/60 s; scale by dt so behavior is
       framerate-independent (pow keeps substep and reference runs consistent) */
    float damp = powf(c->damping, dt * 60.0f);
    for (int i = 0; i < c->n; i++) {
        Particle* p = &c->p[i];
        if (p->pinned) continue;
        float mx = (p->x - p->px) * damp + gx * dt * dt;
        float my = (p->y - p->py) * damp + gy * dt * dt;
        if (max_disp > 0.0f) {
            float m = sqrtf(mx * mx + my * my);
            if (m > max_disp) {
                float s = max_disp / m;
                mx *= s; my *= s;
            }
        }
        p->px = p->x; p->py = p->y;
        p->x += mx; p->y += my;
    }
}

static void chain_solve(Chain* c, int iterations) {
    for (int it = 0; it < iterations; it++) {
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

static void chain_step(Chain* c, float dt, float gx, float gy) {
    chain_integrate(c, dt, gx, gy, 0.0f);
    chain_solve(c, c->iterations);
}

/* substeps needed to cover frame_dt with steps of about target_dt.
   returns 0 on invalid input (non-positive dt) — callers must treat
   a zero substep count as an error, never silently skip simulation. */
static int chain_substeps_for(float frame_dt, float target_dt, int max_substeps) {
    if (frame_dt <= 0.0f || target_dt <= 0.0f) return 0;
    int n = (int)ceilf(frame_dt / target_dt);
    if (n < 1) n = 1;
    if (max_substeps > 0 && n > max_substeps) n = max_substeps;
    return n;
}

/* advance one render frame: adaptive substepping + displacement clamp +
   iteration floor. returns 0 on success, -1 if the substep count is zero. */
static int chain_step_frame(Chain* c, float frame_dt, float gx, float gy,
                            const StepPolicy* pol) {
    int sub = chain_substeps_for(frame_dt, pol->target_dt, pol->max_substeps);
    if (sub <= 0) {
        fprintf(stderr, "chain_step_frame: zero substeps (frame_dt=%g target_dt=%g)\n",
                (double)frame_dt, (double)pol->target_dt);
        return -1;
    }
    int iters = c->iterations < pol->min_iterations ? pol->min_iterations
                                                    : c->iterations;
    float dt = frame_dt / (float)sub;
    for (int s = 0; s < sub; s++) {
        chain_integrate(c, dt, gx, gy, pol->max_disp);
        chain_solve(c, iters);
    }
    return 0;
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

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    Chain c;
    chain_init(&c, 8, 1.0f, 0.0f, 0.0f);
    check("init_len", c.n == 8 && fabsf(chain_seg_error(&c)) < 1e-5f);

    /* let it fall for a while */
    for (int i = 0; i < 200; i++) chain_step(&c, 0.016f, 0.0f, -9.8f);
    check("stays_in_stretch", chain_seg_error(&c) < 0.02f);
    check("pinned_anchor", fabsf(c.p[0].x) < 1e-6f && fabsf(c.p[0].y) < 1e-6f);
    check("falls_down", c.p[7].y < -1.0f);

    /* more iterations must not increase error */
    Chain a, b;
    chain_init(&a, 8, 1.0f, 0.0f, 0.0f);
    chain_init(&b, 8, 1.0f, 0.0f, 0.0f);
    a.iterations = 1; b.iterations = 16;
    for (int i = 0; i < 100; i++) { chain_step(&a, 0.016f, 0.0f, -9.8f); chain_step(&b, 0.016f, 0.0f, -9.8f); }
    check("more_iters_better", chain_seg_error(&b) <= chain_seg_error(&a) + 1e-6f);

    Chain cap;
    chain_init(&cap, 16, 0.5f, 0.0f, 0.0f);
    check("capacity_ok", cap.n == 16);

    /* --- frame-stepping policy: substeps, clamp, iteration floor --- */
    StepPolicy pol = { 1.0f / 120.0f, 64, 0.05f, 4 };

    check("substeps_adaptive", chain_substeps_for(0.05f, pol.target_dt, pol.max_substeps) == 6);
    check("substeps_min_one", chain_substeps_for(0.001f, pol.target_dt, pol.max_substeps) == 1);
    check("substeps_capped", chain_substeps_for(10.0f, pol.target_dt, pol.max_substeps) == pol.max_substeps);
    check("substeps_zero_on_bad_dt", chain_substeps_for(0.0f, pol.target_dt, pol.max_substeps) == 0
                                  && chain_substeps_for(0.016f, 0.0f, pol.max_substeps) == 0);

    Chain zc;
    chain_init(&zc, 8, 1.0f, 0.0f, 0.0f);
    check("frame_step_zero_substeps_err", chain_step_frame(&zc, 0.0f, 0.0f, -9.8f, &pol) == -1);
    check("frame_step_ok", chain_step_frame(&zc, 0.016f, 0.0f, -9.8f, &pol) == 0);

    /* displacement clamp: huge implied velocity must not explode the step */
    Chain cl;
    chain_init(&cl, 8, 1.0f, 0.0f, 0.0f);
    cl.p[7].px = cl.p[7].x + 5.0f;
    {
        float bx = cl.p[7].x, by = cl.p[7].y;
        chain_integrate(&cl, 1.0f / 120.0f, 0.0f, -9.8f, pol.max_disp);
        float dx = cl.p[7].x - bx, dy = cl.p[7].y - by;
        check("max_disp_clamped", sqrtf(dx * dx + dy * dy) <= pol.max_disp + 1e-4f);
    }

    /* iteration floor: iterations=0 must still converge via the floor */
    Chain mi;
    chain_init(&mi, 8, 1.0f, 0.0f, 0.0f);
    mi.iterations = 0;
    for (int i = 0; i < 120; i++) chain_step_frame(&mi, 1.0f / 60.0f, 0.0f, -9.8f, &pol);
    check("min_iter_floor", chain_seg_error(&mi) < 0.02f);

    /* low frame rate vs high frame rate reference. ref runs at 240 Hz,
       low at 20 Hz. Pose error is only bounded over the first 0.5 s:
       past that the swinging chain is chaotic and any two step sizes
       diverge regardless of solver quality. Jitter (second difference
       of the end position) is bounded over the whole run. */
    {
        Chain ref, low, naive;
        chain_init(&ref, 8, 1.0f, 0.0f, 0.0f);
        chain_init(&low, 8, 1.0f, 0.0f, 0.0f);
        chain_init(&naive, 8, 1.0f, 0.0f, 0.0f);
        const float ref_dt = 1.0f / 240.0f, low_dt = 1.0f / 20.0f;
        const int frames = 60;         /* 3 s total */
        const int pose_frames = 10;    /* pose compare window: first 0.5 s */
        float max_pose_err = 0.0f, max_naive_err = 0.0f, jitter = 0.0f, naive_jitter = 0.0f;
        float ex1 = 0, ey1 = 0, ex2 = 0, ey2 = 0; /* end pos, 1 and 2 frames back */
        float nx1 = 0, ny1 = 0, nx2 = 0, ny2 = 0; /* same for the naive chain */
        int have = 0;
        for (int f = 0; f < frames; f++) {
            for (int k = 0; k < 12; k++) chain_step(&ref, ref_dt, 0.0f, -9.8f);
            chain_step_frame(&low, low_dt, 0.0f, -9.8f, &pol);
            chain_step(&naive, low_dt, 0.0f, -9.8f);
            Particle* e = &low.p[low.n - 1];
            Particle* r = &ref.p[ref.n - 1];
            Particle* v = &naive.p[naive.n - 1];
            if (f < pose_frames) {
                float perr = sqrtf((e->x - r->x) * (e->x - r->x) + (e->y - r->y) * (e->y - r->y));
                float nerr = sqrtf((v->x - r->x) * (v->x - r->x) + (v->y - r->y) * (v->y - r->y));
                if (perr > max_pose_err) max_pose_err = perr;
                if (nerr > max_naive_err) max_naive_err = nerr;
            }
            if (have >= 2) {
                float jx = e->x - 2 * ex1 + ex2, jy = e->y - 2 * ey1 + ey2;
                float j = sqrtf(jx * jx + jy * jy);
                if (j > jitter) jitter = j;
                float kx = v->x - 2 * nx1 + nx2, ky = v->y - 2 * ny1 + ny2;
                float kj = sqrtf(kx * kx + ky * ky);
                if (kj > naive_jitter) naive_jitter = kj;
            }
            ex2 = ex1; ey2 = ey1; ex1 = e->x; ey1 = e->y;
            nx2 = nx1; ny2 = ny1; nx1 = v->x; ny1 = v->y;
            if (have < 2) have++;
        }
        printf("max_pose_err=%.4f max_naive_err=%.4f jitter=%.4f naive_jitter=%.4f\n",
               (double)max_pose_err, (double)max_naive_err, (double)jitter, (double)naive_jitter);
        check("lowfps_pose_err_bound", max_pose_err < 0.02f);
        check("lowfps_jitter_bound", jitter < 0.05f);
        check("substeps_beat_naive", max_pose_err <= max_naive_err + 1e-6f);
        check("substeps_calmer_than_naive", jitter <= naive_jitter + 1e-6f);
    }

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
