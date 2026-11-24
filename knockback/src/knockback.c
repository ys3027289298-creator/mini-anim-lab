
/* knockback: impulse-driven hit reaction.
   an attack applies an impulse J (its strength); the receiver's velocity
   jumps by J/m and is then bled off by drag, while position is integrated
   from that velocity. the slide therefore shows the whole force phase
   instead of teleporting the body by a fixed offset. */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define DT          (1.0f / 60.0f)
#define REST_SPEED  0.05f   /* under this speed the body counts as stopped */
#define MAX_STEPS   6000

typedef struct {
    float mass;
    float k;          /* designer drag rate, 1/s: dv = -k v dt */
    float x;
    float v;
    float peak;       /* largest |x| reached */
} Knockback;

static void kb_init(Knockback* kb, float mass, float k) {
    memset(kb, 0, sizeof(*kb));
    kb->mass = mass;
    kb->k = k;
}

/* momentum delivered by the hit: v += J/m, position does not jump */
static void kb_hit(Knockback* kb, float impulse) {
    kb->v += impulse / kb->mass;
}

/* correct drag: F = -m k v => dv/dt = -k v, same decay for every mass.
   semi-implicit euler on v, then x is the integral of v. */
static void kb_step(Knockback* kb, float dt) {
    kb->v -= kb->k * kb->v * dt;
    if (fabsf(kb->v) < REST_SPEED) kb->v = 0.0f;
    kb->x += kb->v * dt;
    if (fabsf(kb->x) > kb->peak) kb->peak = fabsf(kb->x);
}

/* broken drag kept to exercise the diagnostic: F = -c v with one fixed raw
   coefficient, so dv/dt = -(c/m) v. doubling the mass halves the decay —
   the damping term is inversely scaled by the mass ratio. */
typedef struct {
    Knockback kb;
    float c_raw;
} BuggyKnockback;

static void buggy_step(BuggyKnockback* b, float dt) {
    Knockback* kb = &b->kb;
    kb->v -= (b->c_raw / kb->mass) * kb->v * dt;
    if (fabsf(kb->v) < REST_SPEED) kb->v = 0.0f;
    kb->x += kb->v * dt;
    if (fabsf(kb->x) > kb->peak) kb->peak = fabsf(kb->x);
}

typedef struct {
    float ratio;     /* J/m */
    float mass;
    float v0;        /* speed right after the impulse */
    float peak;      /* peak displacement */
    float rest_t;    /* time until |v| < REST_SPEED */
    float coef;      /* drag coefficient c in F = -c v */
    float lambda;    /* effective per-second decay rate c/m */
} RunResult;

static RunResult run_correct(float J, float m, float k, float dt) {
    Knockback kb;
    kb_init(&kb, m, k);
    kb_hit(&kb, J);
    RunResult r;
    memset(&r, 0, sizeof(r));
    r.ratio = J / m;
    r.mass = m;
    r.v0 = kb.v;
    r.coef = m * k;
    r.lambda = k;
    int n = 0;
    for (; n < MAX_STEPS && kb.v != 0.0f; n++) kb_step(&kb, dt);
    r.peak = kb.peak;
    r.rest_t = n * dt;
    return r;
}

static RunResult run_buggy(float J, float m, float c, float dt) {
    BuggyKnockback b;
    memset(&b, 0, sizeof(b));
    b.kb.mass = m;
    b.c_raw = c;
    b.kb.v = J / m;
    RunResult r;
    memset(&r, 0, sizeof(r));
    r.ratio = J / m;
    r.mass = m;
    r.v0 = b.kb.v;
    r.coef = c;
    r.lambda = c / m;
    int n = 0;
    for (; n < MAX_STEPS && b.kb.v != 0.0f; n++) buggy_step(&b, dt);
    r.peak = b.kb.peak;
    r.rest_t = n * dt;
    return r;
}

static void print_table(const char* title, const RunResult* r, int n) {
    printf("%s\n", title);
    printf("%8s %6s %8s %9s %8s %9s %9s\n",
           "ratio", "mass", "v0", "peak_x", "rest_s", "coef_c", "rate_c/m");
    for (int i = 0; i < n; i++) {
        printf("%8.2f %6.2f %8.3f %9.3f %8.3f %9.3f %9.3f\n",
               r[i].ratio, r[i].mass, r[i].v0, r[i].peak,
               r[i].rest_t, r[i].coef, r[i].lambda);
    }
}

/* the task's rule: with groups ordered by descending ratio, if the biggest
   J/m group stops first, damping has been inversely scaled by mass. */
static int damping_inverted(const RunResult* r, int n) {
    int shortest = 0;
    for (int i = 1; i < n; i++)
        if (r[i].rest_t < r[shortest].rest_t) shortest = i;
    return shortest == 0;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    const float J = 10.0f;
    const float k = 2.0f;                 /* drag rate, 1/s */
    const float masses[3] = {1.0f, 2.0f, 4.0f};  /* ratios 10, 5, 2.5 */

    RunResult good[3], bad[3];
    for (int i = 0; i < 3; i++) {
        good[i] = run_correct(J, masses[i], k, DT);
        bad[i] = run_buggy(J, masses[i], k, DT);  /* same numeric coef, raw */
    }

    /* impulse gives v0 = J/m, nothing moves before the hit */
    Knockback kb;
    kb_init(&kb, 2.0f, k);
    check("no_impulse_idle", kb.x == 0.0f && kb.v == 0.0f);
    kb_hit(&kb, J);
    check("impulse_v0", fabsf(kb.v - J / kb.mass) < 1e-6f);

    /* velocity keeps decaying and never reverses; position is its integral */
    float prev_v = kb.v, prev_x = kb.x;
    int monotone = 1;
    for (int i = 0; i < 30; i++) {
        kb_step(&kb, DT);
        monotone = monotone && kb.v >= 0.0f && kb.v <= prev_v && kb.x >= prev_x;
        prev_v = kb.v; prev_x = kb.x;
    }
    check("velocity_decays", monotone);

    /* continuous solution: total slide is the integral v0/k (euler within 5%) */
    check("slide_is_integral",
          fabsf(good[0].peak - good[0].v0 / k) / (good[0].v0 / k) < 0.05f);

    /* same impulse: a heavier body must not slide as far */
    check("heavier_shorter_slide", good[2].peak < good[1].peak &&
                                   good[1].peak < good[0].peak);

    /* diagnostic: correct drag never has the biggest-ratio group stop first;
       the raw-coefficient bug always does */
    check("diagnostic_clean", damping_inverted(good, 3) == 0);
    check("diagnostic_fires_on_bug", damping_inverted(bad, 3) == 1);

    /* bug signature: with F = -c v the total slide J/c is mass independent */
    float spread = fabsf(bad[0].peak - bad[2].peak) / bad[0].peak;
    check("bug_peak_mass_independent", spread < 0.05f);

    print_table("\n-- correct drag F = -m*k*v (k=2.00/s) --", good, 3);
    print_table("\n-- buggy drag F = -c*v, fixed c (damping inversely scaled) --",
                bad, 3);
    printf("\ninverted_damping_correct=%d\n", damping_inverted(good, 3));
    printf("inverted_damping_buggy=%d\n", damping_inverted(bad, 3));

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
