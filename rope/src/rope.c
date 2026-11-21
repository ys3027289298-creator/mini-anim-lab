
/* rope: verlet particle chain with unilateral (pull-only) distance constraints.
   Tension is back-derived from accumulated PBD corrections; a segment whose
   tension crosses the limit fractures and degrades to a free end. */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define MAX_P 16

typedef struct { float x, y, px, py; float inv_mass; int pinned; } RopeParticle;

typedef struct {
    RopeParticle p[MAX_P];
    int n;
    float rest[MAX_P];      /* rest length of segment i (between p[i], p[i+1]) */
    int broken[MAX_P];      /* fractured segment: constraint never enforced */
    float tension[MAX_P];   /* tension estimate from the most recent step */
    float tension_limit;    /* fracture threshold, INFINITY = unbreakable */
    int iterations;
    float damping;
    int step_count;
    int break_seg;          /* segment fractured on the latest step, -1 if none */
    int break_step;         /* step index of the latest fracture, -1 */
} Rope;

static void rope_init(Rope* r, int n, float seg, float x0, float y0) {
    memset(r, 0, sizeof(*r));
    r->n = n;
    r->iterations = 32;
    r->damping = 0.99f;
    r->tension_limit = INFINITY;
    r->break_seg = -1;
    r->break_step = -1;
    for (int i = 0; i < n; i++) {
        r->p[i].x = x0;
        r->p[i].y = y0 - seg * (float)i;          /* hangs straight down */
        r->p[i].px = r->p[i].x;
        r->p[i].py = r->p[i].y;
        r->p[i].inv_mass = 1.0f;                  /* unit mass: m = 1 */
    }
    r->p[0].pinned = 1;
    r->p[0].inv_mass = 0.0f;
    for (int i = 0; i < n - 1; i++) r->rest[i] = seg;
}

/* non-positive limits (including NaN) are rejected; the old limit is kept */
static int rope_set_limit(Rope* r, float tmax) {
    if (!(tmax > 0.0f)) return -1;
    r->tension_limit = tmax;
    return 0;
}

static int rope_broken(const Rope* r, int seg) {
    return seg >= 0 && seg < r->n - 1 && r->broken[seg];
}

static float rope_tension(const Rope* r, int seg) {
    if (seg < 0 || seg >= r->n - 1) return 0.0f;
    return r->tension[seg];
}

/* worst residual stretch over the segments still in one piece */
static float rope_seg_error(const Rope* r) {
    float worst = 0.0f;
    for (int i = 0; i < r->n - 1; i++) {
        if (r->broken[i]) continue;
        float dx = r->p[i + 1].x - r->p[i].x;
        float dy = r->p[i + 1].y - r->p[i].y;
        float d = sqrtf(dx * dx + dy * dy);
        float e = fabsf(d - r->rest[i]);
        if (e > worst) worst = e;
    }
    return worst;
}

static void rope_step(Rope* r, float dt, float gx, float gy) {
    if (dt <= 0.0f) return;
    for (int i = 0; i < r->n; i++) {
        RopeParticle* p = &r->p[i];
        if (p->pinned) continue;
        float vx = (p->x - p->px) * r->damping;
        float vy = (p->y - p->py) * r->damping;
        p->px = p->x; p->py = p->y;
        p->x += vx + gx * dt * dt;
        p->y += vy + gy * dt * dt;
    }

    float tent[MAX_P] = {0};
    for (int it = 0; it < r->iterations; it++) {
        for (int i = 0; i < r->n - 1; i++) {
            if (r->broken[i]) continue;
            RopeParticle* a = &r->p[i];
            RopeParticle* b = &r->p[i + 1];
            float dx = b->x - a->x, dy = b->y - a->y;
            float d = sqrtf(dx * dx + dy * dy);
            if (d < 1e-6f) continue;
            float stretch = d - r->rest[i];
            if (stretch <= 0.0f) continue;                 /* slack: rope never pushes */
            float wa = a->inv_mass, wb = b->inv_mass;
            float wsum = wa + wb;
            if (wsum <= 0.0f) continue;
            float k = stretch / d;
            a->x += dx * k * (wa / wsum);
            a->y += dy * k * (wa / wsum);
            b->x -= dx * k * (wb / wsum);
            b->y -= dy * k * (wb / wsum);
            /* PBD correction equals impulse displacement w_b*T*dt^2, hence
               T = stretch / (wsum * dt^2); sum iterations for the step. */
            tent[i] += stretch / (wsum * dt * dt);
        }
    }
    for (int i = 0; i < r->n - 1; i++) r->tension[i] = tent[i];

    r->break_seg = -1;
    /* only the single most overloaded segment gives way this step */
    int worst = -1;
    float worst_ratio = 1.0f;
    for (int i = 0; i < r->n - 1; i++) {
        if (r->broken[i]) continue;
        float ratio = r->tension[i] / r->tension_limit;
        if (ratio > worst_ratio) { worst_ratio = ratio; worst = i; }
    }
    r->step_count++;
    if (worst >= 0) {
        /* flag flip only: positions are untouched, so the break is continuous */
        r->broken[worst] = 1;
        r->break_seg = worst;
        r->break_step = r->step_count;
    }
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    /* --- basic hanging rope ------------------------------------------------ */
    Rope c;
    rope_init(&c, 8, 1.0f, 0.0f, 0.0f);
    check("init_len", c.n == 8 && rope_seg_error(&c) < 1e-5f);
    for (int i = 0; i < 200; i++) rope_step(&c, 0.016f, 0.0f, -9.8f);
    check("stays_in_stretch", rope_seg_error(&c) < 0.02f);
    check("pinned_anchor", c.p[0].x == 0.0f && c.p[0].y == 0.0f);
    check("falls_down", c.p[7].y < -1.0f);
    check("tension_nonneg", rope_tension(&c, 0) >= 0.0f && rope_tension(&c, 6) >= 0.0f);
    check("no_break_default", c.break_seg == -1 && c.break_step == -1);

    /* --- zero / negative / NaN limits must be rejected --------------------- */
    Rope v;
    rope_init(&v, 4, 1.0f, 0.0f, 0.0f);
    check("reject_zero", rope_set_limit(&v, 0.0f) == -1);
    check("reject_neg", rope_set_limit(&v, -10.0f) == -1);
    check("reject_nan", rope_set_limit(&v, NAN) == -1);
    check("limit_untouched", isinf(v.tension_limit));
    check("accept_positive", rope_set_limit(&v, 1.0f) == 0 && v.tension_limit == 1.0f);

    /* --- analytic statics: vertical chain, unit masses, settled ----------- */
    enum { SN = 6 };
    Rope s;
    rope_init(&s, SN, 1.0f, 0.0f, 0.0f);
    for (int i = 0; i < 3000; i++) rope_step(&s, 0.01f, 0.0f, -9.8f);
    check("settled_rest", fabsf(s.p[SN - 1].y - s.p[SN - 1].py) < 1e-3f);
    check("settled_vertical", fabsf(s.p[SN - 1].x) < 1e-5f);

    float max_err = 0.0f, max_bound = 0.0f;
    int identity_ok = 1, triangle_ok = 1;
    for (int i = 0; i < SN - 1; i++) {
        float t_static = (float)(SN - 1 - i) * 9.8f;       /* weight below seg i */
        float err = s.tension[i] - t_static;
        float sum_r = 0.0f, bound = 0.0f;
        for (int j = i + 1; j < SN; j++) {
            float up = s.tension[j - 1];
            float down = (j < SN - 1) ? s.tension[j] : 0.0f;
            float rj = up - down - 9.8f;                  /* force residual on p[j], m=1 */
            sum_r += rj;
            bound += fabsf(rj);
        }
        if (fabsf(err - sum_r) > 1e-3f) identity_ok = 0;  /* telescoping identity */
        if (fabsf(err) > bound + 1e-4f) triangle_ok = 0;  /* guaranteed upper bound */
        if (fabsf(err) > max_err) max_err = fabsf(err);
        if (bound > max_bound) max_bound = bound;
    }
    check("residual_identity", identity_ok);
    check("err_under_bound", triangle_ok);
    check("tension_matches_statics", max_err < 0.02f * (SN - 1) * 9.8f);
    printf("threshold_bound=%.5f measured_max_err=%.5f\n", max_bound, max_err);

    /* --- fracture: limit between top-link (4g) and next link (3g) --------- */
    enum { BN = 5 };
    Rope live, ref;
    rope_init(&live, BN, 1.0f, 0.0f, 0.0f);
    rope_init(&ref, BN, 1.0f, 0.0f, 0.0f);
    for (int i = 0; i < 3000; i++) {
        rope_step(&live, 0.01f, 0.0f, -9.8f);
        rope_step(&ref, 0.01f, 0.0f, -9.8f);
    }
    rope_set_limit(&live, 3.5f * 9.8f);                    /* 34.3: only 4g = 39.2 exceeds */
    int broke_step = -1;
    for (int i = 0; i < 2000 && broke_step < 0; i++) {
        rope_step(&live, 0.01f, 0.0f, -9.8f);
        rope_step(&ref, 0.01f, 0.0f, -9.8f);
        if (live.break_seg >= 0) broke_step = live.step_count;
    }
    check("break_happened", broke_step > 0);
    check("break_top_segment", live.break_seg == 0 && rope_broken(&live, 0));
    check("only_one_break", live.break_step == broke_step);

    /* continuity: at the break step both ropes must still be bit-identical,
       since fracturing only flips a flag after the solve */
    int continuous = 1;
    for (int i = 0; i < BN; i++) {
        if (live.p[i].x != ref.p[i].x || live.p[i].y != ref.p[i].y ||
            live.p[i].px != ref.p[i].px || live.p[i].py != ref.p[i].py)
            continuous = 0;
    }
    check("break_position_continuous", continuous);

    /* free end: internal constraint forces cancel in pairs over the detached
       sub-chain, so its center of mass follows the exact unconstrained verlet
       step even while its own segment constraints keep firing */
    float sx = 0, spx = 0, sy = 0, spy = 0;
    for (int i = 1; i < BN; i++) {
        sx += live.p[i].x;   spx += live.p[i].px;
        sy += live.p[i].y;   spy += live.p[i].py;
    }
    float m = (float)(BN - 1);
    float e_sx = sx + (sx - spx) * live.damping;
    float e_sy = sy + (sy - spy) * live.damping + (-9.8f) * 0.01f * 0.01f * m;
    rope_step(&live, 0.01f, 0.0f, -9.8f);
    float asx = 0, asy = 0;
    for (int i = 1; i < BN; i++) { asx += live.p[i].x; asy += live.p[i].y; }
    check("free_end_unconstrained", fabsf(asx - e_sx) < 1e-5f && fabsf(asy - e_sy) < 1e-5f);

    for (int i = 0; i < 60; i++) rope_step(&live, 0.01f, 0.0f, -9.8f);
    float d01 = hypotf(live.p[1].x - live.p[0].x, live.p[1].y - live.p[0].y);
    check("broken_link_opens", d01 > live.rest[0] + 0.5f);
    float inner_worst = 0.0f;
    for (int i = 1; i < BN - 1; i++) {
        float d = hypotf(live.p[i + 1].x - live.p[i].x, live.p[i + 1].y - live.p[i].y);
        if (fabsf(d - live.rest[i]) > inner_worst) inner_worst = fabsf(d - live.rest[i]);
    }
    check("free_fall_together", inner_worst < 0.02f);      /* lower chain stays intact */
    check("stays_broken", rope_broken(&live, 0));
    check("anchor_unmoved", live.p[0].x == 0.0f && live.p[0].y == 0.0f);
    check("ref_never_broke", ref.break_seg == -1);

    /* --- slack rope: compression must never read as tension --------------- */
    Rope slack;
    rope_init(&slack, 4, 1.0f, 0.0f, 0.0f);
    rope_set_limit(&slack, 0.01f);                         /* tiny positive limit */
    for (int i = 0; i < slack.n; i++) {                    /* squash segments */
        slack.p[i].y = -0.1f * (float)i;
        slack.p[i].py = slack.p[i].y;
    }
    for (int i = 0; i < 50; i++) rope_step(&slack, 0.01f, 0.0f, 0.0f);
    int slack_ok = 1;
    for (int i = 0; i < slack.n - 1; i++)
        if (rope_tension(&slack, i) != 0.0f || rope_broken(&slack, i)) slack_ok = 0;
    check("slack_no_tension", slack_ok);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
