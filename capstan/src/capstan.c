
/* capstan: rope wrapped around a post with anchor friction. tension ratio
   follows the wrap angle (Euler-Eytelwein, T_load/T_hold = e^(mu*theta));
   the numeric model integrates tangential friction over discrete contact
   points and is checked against the analytic ratio with an error bound.
   slip accumulates while overloaded and the anchor releases only once the
   accumulated slip passes a threshold. */
#include <stdio.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define MAX_WRAP (16.0f * (float)M_PI)  /* more than 8 full turns is rejected */
#define MAX_SEG 4096

typedef struct {
    float mu;          /* friction coefficient, must be >= 0 */
    float wrap_angle;  /* total wrap angle in radians, (0, MAX_WRAP] */
    int segments;      /* contact point count along the wrap */
} Capstan;

/* returns 0 on success, -1 when parameters are rejected */
static int capstan_init(Capstan* c, float mu, float wrap_angle, int segments) {
    if (!isfinite(mu) || mu < 0.0f) return -1;
    if (!isfinite(wrap_angle) || wrap_angle <= 0.0f || wrap_angle > MAX_WRAP) return -1;
    if (segments < 1 || segments > MAX_SEG) return -1;
    c->mu = mu;
    c->wrap_angle = wrap_angle;
    c->segments = segments;
    return 0;
}

/* analytic Euler-Eytelwein tension ratio e^(mu*theta) */
static float capstan_ratio_analytic(const Capstan* c) {
    return expf(c->mu * c->wrap_angle);
}

/* numeric ratio: at each contact point the normal load is T*dtheta and the
   tangential friction adds at most mu*T*dtheta, so T *= (1 + mu*dtheta) */
static float capstan_ratio_numeric(const Capstan* c) {
    float dtheta = c->wrap_angle / (float)c->segments;
    float t = 1.0f;
    for (int i = 0; i < c->segments; i++)
        t *= 1.0f + c->mu * dtheta;
    return t;
}

/* a priori upper bound on |numeric - analytic| for the forward discretization:
   e^x - (1 + x/n)^n <= e^x * x^2 / (2n) with x = mu*theta */
static float capstan_error_bound(const Capstan* c) {
    float x = c->mu * c->wrap_angle;
    float n = (float)c->segments;
    return expf(x) * x * x / (2.0f * n);
}

typedef struct {
    Capstan c;
    float slip;            /* accumulated slip while overloaded */
    float slip_threshold;  /* anchor releases only beyond this */
    int released;
} Anchor;

static int anchor_init(Anchor* a, float mu, float wrap_angle, int segments,
                       float slip_threshold) {
    if (capstan_init(&a->c, mu, wrap_angle, segments) != 0) return -1;
    if (!isfinite(slip_threshold) || slip_threshold <= 0.0f) return -1;
    a->slip = 0.0f;
    a->slip_threshold = slip_threshold;
    a->released = 0;
    return 0;
}

/* one step under load; returns 1 once the anchor has released.
   while the demanded ratio t_load/t_hold exceeds the analytic capacity the
   rope slips proportionally to the overload; below capacity it sticks. */
static int anchor_load(Anchor* a, float t_load, float t_hold, float dt) {
    if (a->released) return 1;
    if (t_hold <= 0.0f || dt <= 0.0f) return 0;
    float capacity = capstan_ratio_analytic(&a->c);
    float overload = t_load / t_hold - capacity;
    if (overload > 0.0f) {
        a->slip += overload * dt;
        if (a->slip > a->slip_threshold) a->released = 1;
    }
    return a->released;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    Capstan c;
    check("init_ok", capstan_init(&c, 0.5f, (float)M_PI, 64) == 0);

    /* invalid parameters are rejected */
    check("reject_neg_mu", capstan_init(&c, -0.1f, 1.0f, 64) != 0);
    check("reject_nan_mu", capstan_init(&c, NAN, 1.0f, 64) != 0);
    check("reject_zero_angle", capstan_init(&c, 0.5f, 0.0f, 64) != 0);
    check("reject_neg_angle", capstan_init(&c, 0.5f, -1.0f, 64) != 0);
    check("reject_huge_angle", capstan_init(&c, 0.5f, 17.0f * (float)M_PI, 64) != 0);
    check("reject_no_segments", capstan_init(&c, 0.5f, 1.0f, 0) != 0);

    /* known value: mu=0.5, theta=pi -> e^(pi/2) ~ 4.8105 */
    capstan_init(&c, 0.5f, (float)M_PI, 256);
    float analytic = capstan_ratio_analytic(&c);
    check("analytic_value", fabsf(analytic - 4.8105f) < 1e-3f);

    /* numeric ratio tracks analytic within the a priori bound */
    float numeric = capstan_ratio_numeric(&c);
    float err = fabsf(numeric - analytic);
    float bound = capstan_error_bound(&c);
    check("within_bound", err <= bound);
    check("bound_positive", bound > 0.0f);

    /* finer contact discretization shrinks the error */
    Capstan coarse, fine;
    capstan_init(&coarse, 0.5f, (float)M_PI, 8);
    capstan_init(&fine, 0.5f, (float)M_PI, 512);
    float err_coarse = fabsf(capstan_ratio_numeric(&coarse) - capstan_ratio_analytic(&coarse));
    float err_fine = fabsf(capstan_ratio_numeric(&fine) - capstan_ratio_analytic(&fine));
    check("refinement_helps", err_fine < err_coarse);

    /* zero friction means no mechanical advantage */
    capstan_init(&c, 0.0f, (float)M_PI, 64);
    check("zero_mu_ratio_one", fabsf(capstan_ratio_numeric(&c) - 1.0f) < 1e-6f);

    /* anchor: sticks within capacity, slips under overload, releases only
       past the slip threshold */
    Anchor a;
    check("anchor_init_ok", anchor_init(&a, 0.5f, (float)M_PI, 64, 0.5f) == 0);
    check("anchor_reject_neg_mu", anchor_init(&a, -0.5f, (float)M_PI, 64, 0.5f) != 0);
    check("anchor_reject_bad_angle", anchor_init(&a, 0.5f, -1.0f, 64, 0.5f) != 0);
    check("anchor_reject_bad_thresh", anchor_init(&a, 0.5f, (float)M_PI, 64, 0.0f) != 0);

    anchor_init(&a, 0.5f, (float)M_PI, 64, 0.5f);
    float cap = capstan_ratio_analytic(&a.c);
    for (int i = 0; i < 100; i++) anchor_load(&a, cap * 0.9f, 1.0f, 0.016f);
    check("sticks_within_capacity", !a.released && a.slip < 1e-6f);

    /* brief overload: slips a little but must not release yet */
    anchor_init(&a, 0.5f, (float)M_PI, 64, 0.5f);
    for (int i = 0; i < 5; i++) anchor_load(&a, cap * 1.1f, 1.0f, 0.016f);
    check("brief_overload_holds", !a.released && a.slip > 0.0f);

    /* sustained overload: slip accumulates past the threshold and releases */
    anchor_init(&a, 0.5f, (float)M_PI, 64, 0.5f);
    int released = 0;
    for (int i = 0; i < 2000 && !released; i++)
        released = anchor_load(&a, cap * 1.5f, 1.0f, 0.016f);
    check("sustained_overload_releases", released);
    check("slip_passed_threshold", a.slip > a.slip_threshold);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
