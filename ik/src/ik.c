
/* ik: analytic two-bone IK + FABRIK iterative solver for chains of any length. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define IK_MAX_BONES 16
#define IK_EPS 1e-6f

typedef struct { float x, y; } V2;

static float vlen(V2 v) { return sqrtf(v.x * v.x + v.y * v.y); }
static V2 vsub(V2 a, V2 b) { V2 r = { a.x - b.x, a.y - b.y }; return r; }
static V2 vadd(V2 a, V2 b) { V2 r = { a.x + b.x, a.y + b.y }; return r; }
static V2 vscale(V2 v, float s) { V2 r = { v.x * s, v.y * s }; return r; }
static V2 vnorm(V2 v) { float l = vlen(v); V2 r = { 0, 0 }; if (l > 0) { r.x = v.x / l; r.y = v.y / l; } return r; }
static float vdist(V2 a, V2 b) { return vlen(vsub(a, b)); }

typedef struct {
    int ok;
    V2 joint;
    float a1;   /* angle of first bone */
    float a2;   /* angle of second bone */
} IkResult;

/* solve for target with bone lengths l1,l2; pole decides bend direction */
static IkResult solve_two_bone(V2 root, V2 target, float l1, float l2, V2 pole) {
    IkResult r;
    memset(&r, 0, sizeof(r));
    V2 d = vsub(target, root);
    float dist = vlen(d);
    if (l1 <= 0 || l2 <= 0) return r;
    if (dist > l1 + l2 + 1e-5f) return r;                 /* unreachable */
    if (dist < fabsf(l1 - l2) - 1e-5f) return r;          /* too close */
    if (dist < 1e-6f) return r;

    float cos_a = (l1 * l1 + dist * dist - l2 * l2) / (2 * l1 * dist);
    if (cos_a > 1.0f) cos_a = 1.0f;
    if (cos_a < -1.0f) cos_a = -1.0f;
    float a = acosf(cos_a);

    V2 dir = vnorm(d);
    float side = (dir.x * pole.y - dir.y * pole.x) >= 0 ? 1.0f : -1.0f;
    float ang = atan2f(dir.y, dir.x) + side * a;

    r.joint.x = root.x + l1 * cosf(ang);
    r.joint.y = root.y + l1 * sinf(ang);
    r.a1 = ang;
    V2 d2 = vsub(target, r.joint);
    r.a2 = atan2f(d2.y, d2.x);
    r.ok = 1;
    return r;
}

/* ------------------------------------------------------------------------- */
/* iterative solver                                                          */
/* ------------------------------------------------------------------------- */

typedef enum {
    IK_OK = 0,
    IK_ERR_BAD_CHAIN = 1,          /* count out of range or non-finite input */
    IK_ERR_BAD_LENGTH = 2,         /* bone length <= 0 or non-finite */
    IK_ERR_BAD_TOLERANCE = 3,      /* tolerance <= 0 or non-finite */
    IK_ERR_BAD_ITERATIONS = 4,     /* max_iterations <= 0 */
    IK_ERR_ITERATION_DISABLED = 5, /* iterative off and chain longer than 2 */
    IK_ERR_NOT_CONVERGED = 6,      /* hit max_iterations above tolerance */
    IK_ERR_NO_SOLUTION = 7         /* analytic path has no solution */
} IkStatus;

typedef struct {
    int iterative;        /* 0 = legacy path: analytic for 2 bones, no-op otherwise */
    int max_iterations;   /* > 0 when iterative is on */
    float tolerance;      /* > 0; convergence when max joint displacement <= this */
} IkOptions;

typedef struct {
    int bone_count;
    V2 root;
    float length[IK_MAX_BONES];
    V2 pole;                    /* bend direction for the analytic path */
    V2 pose[IK_MAX_BONES + 1];  /* seed pose for iterative; pose[0] is replaced by root */
} IkChain;

typedef struct {
    IkStatus status;
    int iterations_used;
    float residual;                 /* last iteration's max joint displacement */
    int used_analytic;
    V2 joints[IK_MAX_BONES + 1];    /* bone_count + 1 points, root first */
} IkSolveResult;

static IkOptions ik_default_options(void) {
    IkOptions o;
    o.iterative = 1;
    o.max_iterations = 64;
    o.tolerance = 1e-5f;
    return o;
}

static int finite2(V2 v) { return isfinite(v.x) && isfinite(v.y); }

static IkStatus validate_chain(const IkChain* c, int need_tolerance, const IkOptions* opt) {
    if (!c || c->bone_count < 1 || c->bone_count > IK_MAX_BONES) return IK_ERR_BAD_CHAIN;
    if (!finite2(c->root)) return IK_ERR_BAD_CHAIN;
    for (int i = 0; i < c->bone_count; i++) {
        if (!isfinite(c->length[i]) || c->length[i] <= 0.0f) return IK_ERR_BAD_LENGTH;
    }
    if (need_tolerance) {
        if (!opt) return IK_ERR_BAD_CHAIN;
        if (!isfinite(opt->tolerance) || opt->tolerance <= 0.0f) return IK_ERR_BAD_TOLERANCE;
        if (opt->max_iterations <= 0) return IK_ERR_BAD_ITERATIONS;
    }
    return IK_OK;
}

/* place a at distance len from b, along b->a; fallback_dir is a unit vector used
   when the two points coincide. */
static void place_from(V2* a, V2 b, float len, V2 fallback_dir) {
    V2 d = vsub(*a, b);
    float r = vlen(d);
    V2 u = r > IK_EPS ? vscale(d, 1.0f / r) : vnorm(fallback_dir);
    *a = vadd(b, vscale(u, len));
}

/* FABRIK: alternate a forward pass (pin the end effector to the target) and a
   backward pass (pin the root).  Convergence is measured by the largest joint
   displacement between the start and the end of one full iteration. */
static IkSolveResult fabrik_solve(const IkChain* chain, V2 target,
                                  const IkOptions* opt,
                                  float* history, int history_capacity) {
    IkSolveResult out;
    memset(&out, 0, sizeof(out));
    out.status = validate_chain(chain, 1, opt);
    if (out.status != IK_OK) return out;
    if (!finite2(target)) { out.status = IK_ERR_BAD_CHAIN; return out; }

    int n = chain->bone_count;
    V2 p[IK_MAX_BONES + 1];
    V2 seed_dir[IK_MAX_BONES];
    for (int i = 0; i <= n; i++) p[i] = chain->pose[i];
    p[0] = chain->root;
    for (int i = 0; i <= n; i++) {
        if (!finite2(p[i])) { out.status = IK_ERR_BAD_CHAIN; return out; }
    }
    for (int i = 0; i < n; i++) {
        seed_dir[i] = vnorm(vsub(p[i + 1], p[i]));
        if (vlen(seed_dir[i]) < 0.5f) { seed_dir[i].x = 1.0f; seed_dir[i].y = 0.0f; }
    }

    int recorded = 0, converged = 0;
    float residual = 0.0f;
    for (int iter = 1; iter <= opt->max_iterations; iter++) {
        V2 before[IK_MAX_BONES + 1];
        for (int i = 0; i <= n; i++) before[i] = p[i];

        /* forward: end effector at target, walk back to the root */
        p[n] = target;
        for (int i = n - 1; i >= 0; i--) {
            place_from(&p[i], p[i + 1], chain->length[i],
                       vscale(seed_dir[i], -1.0f));
        }

        /* backward: root fixed, walk out to the end */
        p[0] = chain->root;
        for (int i = 1; i <= n; i++) {
            place_from(&p[i], p[i - 1], chain->length[i - 1], seed_dir[i - 1]);
        }

        float max_disp = 0.0f;
        for (int i = 1; i <= n; i++) {
            float d = vdist(p[i], before[i]);
            if (d > max_disp) max_disp = d;
        }
        residual = max_disp;
        if (history && recorded < history_capacity) history[recorded++] = max_disp;

        out.iterations_used = iter;
        if (max_disp <= opt->tolerance) { converged = 1; break; }
    }

    for (int i = 0; i <= n; i++) out.joints[i] = p[i];
    out.residual = residual;
    out.status = converged ? IK_OK : IK_ERR_NOT_CONVERGED;
    return out;
}

/* Unified entry point.
   iterative off: two bones follow the legacy analytic solver (identical pose);
   longer chains are left untouched and reported as IK_ERR_ITERATION_DISABLED. */
static IkSolveResult ik_solve(const IkChain* chain, V2 target, IkOptions opt) {
    IkSolveResult out;
    memset(&out, 0, sizeof(out));

    IkStatus st = validate_chain(chain, opt.iterative != 0, &opt);
    if (st != IK_OK) { out.status = st; return out; }

    if (!opt.iterative) {
        if (chain->bone_count == 2) {
            IkResult r = solve_two_bone(chain->root, target, chain->length[0],
                                        chain->length[1], chain->pole);
            if (!r.ok) { out.status = IK_ERR_NO_SOLUTION; return out; }
            out.status = IK_OK;
            out.used_analytic = 1;
            out.iterations_used = 0;
            out.joints[0] = chain->root;
            out.joints[1] = r.joint;
            out.joints[2] = target;
            return out;
        }
        for (int i = 0; i <= chain->bone_count; i++) out.joints[i] = chain->pose[i];
        out.status = IK_ERR_ITERATION_DISABLED;
        return out;
    }

    return fabrik_solve(chain, target, &opt, NULL, 0);
}

/* ------------------------------------------------------------------------- */
/* self tests                                                                */
/* ------------------------------------------------------------------------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

#define LEN_TOL 1e-4f
/* For an exactly reachable two-bone target, the converged position error is
   bounded by this multiple of the configured tolerance (measured across the
   direction x reach grid: 4.98x .. 5.61x). */
#define TWO_BONE_ERR_RATIO 10.0f

static IkChain make_chain(int n, V2 root, const float* lengths) {
    IkChain c;
    memset(&c, 0, sizeof(c));
    c.bone_count = n;
    c.root = root;
    c.pole = (V2){ 0.0f, 1.0f };
    for (int i = 0; i < n; i++) c.length[i] = lengths[i];
    c.pose[0] = root;
    for (int i = 1; i <= n; i++) c.pose[i] = vadd(c.pose[i - 1], (V2){ 1.0f, 0.0f });
    return c;
}

static int lengths_preserved(const IkChain* c, const IkSolveResult* r) {
    for (int i = 0; i < c->bone_count; i++) {
        if (fabsf(vdist(r->joints[i], r->joints[i + 1]) - c->length[i]) > LEN_TOL) return 0;
    }
    return 1;
}

static void test_analytic(void) {
    V2 root = { 0, 0 }, target = { 1, 0 }, pole = { 0, 1 };
    IkResult r = solve_two_bone(root, target, 1.0f, 1.0f, pole);
    check("reachable", r.ok == 1);
    check("joint_on_circle", fabsf(vlen(vsub(r.joint, root)) - 1.0f) < 1e-4f);
    check("end_matches", fabsf(vlen(vsub(r.joint, target)) - 1.0f) < 1e-4f);
    check("pole_up", r.joint.y > 0.0f);

    IkResult r2 = solve_two_bone(root, target, 1.0f, 1.0f, (V2){ 0, -1 });
    check("pole_down", r2.joint.y < 0.0f);

    IkResult far = solve_two_bone(root, (V2){ 5, 0 }, 1.0f, 1.0f, pole);
    check("unreachable", far.ok == 0);

    IkResult close = solve_two_bone(root, (V2){ 0.1f, 0 }, 1.0f, 0.2f, pole);
    check("too_close", close.ok == 0);

    IkResult zero = solve_two_bone(root, target, 0.0f, 1.0f, pole);
    check("zero_len", zero.ok == 0);
}

static void test_fabrik_chains(void) {
    float l3[3] = { 1.0f, 1.0f, 1.0f };
    IkChain c3 = make_chain(3, (V2){ 0, 0 }, l3);
    V2 t3 = { 1.2f, 1.4f };
    IkOptions opt = ik_default_options();
    opt.tolerance = 1e-6f;
    opt.max_iterations = 128;

    IkSolveResult r3 = ik_solve(&c3, t3, opt);
    check("fabrik3_converged", r3.status == IK_OK);
    check("fabrik3_end", vdist(r3.joints[3], t3) < 1e-4f);
    check("fabrik3_lengths", lengths_preserved(&c3, &r3));
    check("fabrik3_root", vdist(r3.joints[0], c3.root) < 1e-7f);

    float l4[4] = { 1.0f, 1.0f, 1.0f, 1.0f };
    IkChain c4 = make_chain(4, (V2){ 0, 0 }, l4);
    V2 t4 = { 0.8f, 2.6f };
    IkSolveResult r4 = ik_solve(&c4, t4, opt);
    check("fabrik4_converged", r4.status == IK_OK);
    check("fabrik4_end", vdist(r4.joints[4], t4) < 1e-4f);
    check("fabrik4_lengths", lengths_preserved(&c4, &r4));

    /* unreachable: chain straightens toward the target */
    V2 far = { 10.0f, 0.0f };
    IkSolveResult rf = ik_solve(&c3, far, opt);
    check("fabrik_unreachable_converged", rf.status == IK_OK);
    check("fabrik_unreachable_end", fabsf(vdist(rf.joints[3], c3.root) - 3.0f) < 1e-4f);
    check("fabrik_unreachable_dir", fabsf(rf.joints[3].y) < 1e-4f);
}

static void test_rejections(void) {
    float l2[2] = { 1.0f, 1.0f };
    IkChain c = make_chain(2, (V2){ 0, 0 }, l2);
    V2 target = { 1.0f, 0.5f };
    IkOptions opt = ik_default_options();
    opt.tolerance = 1e-6f;

    IkOptions zero_tol = opt; zero_tol.tolerance = 0.0f;
    check("reject_tol_zero", ik_solve(&c, target, zero_tol).status == IK_ERR_BAD_TOLERANCE);

    IkOptions neg_tol = opt; neg_tol.tolerance = -1e-3f;
    check("reject_tol_negative", ik_solve(&c, target, neg_tol).status == IK_ERR_BAD_TOLERANCE);

    c.length[1] = -1.0f;
    check("reject_len_negative", ik_solve(&c, target, opt).status == IK_ERR_BAD_LENGTH);
    c.length[1] = 0.0f;
    check("reject_len_zero", ik_solve(&c, target, opt).status == IK_ERR_BAD_LENGTH);
    c.length[1] = 1.0f;

    IkOptions zero_iter = opt; zero_iter.max_iterations = 0;
    check("reject_iter_zero", ik_solve(&c, target, zero_iter).status == IK_ERR_BAD_ITERATIONS);

    float l3[3] = { 1.0f, -0.5f, 1.0f };
    IkChain c3 = make_chain(3, (V2){ 0, 0 }, l3);
    check("reject_len_negative_3", ik_solve(&c3, target, opt).status == IK_ERR_BAD_LENGTH);
}

static void test_iteration_switch(void) {
    float l3[3] = { 1.0f, 1.0f, 1.0f };
    IkChain c3 = make_chain(3, (V2){ 0, 0 }, l3);
    V2 target = { 1.0f, 1.0f };
    IkOptions off = ik_default_options();
    off.iterative = 0;

    IkSolveResult r = ik_solve(&c3, target, off);
    int same = 1;
    for (int i = 0; i <= 3; i++) {
        if (r.joints[i].x != c3.pose[i].x || r.joints[i].y != c3.pose[i].y) same = 0;
    }
    check("off3_status", r.status == IK_ERR_ITERATION_DISABLED);
    check("off3_pose_unchanged", same);

    /* two bones, iteration off: bit-identical to the analytic solver */
    float l2[2] = { 1.0f, 0.8f };
    IkChain c2 = make_chain(2, (V2){ 0, 0 }, l2);
    c2.pole = (V2){ 0.0f, 1.0f };
    V2 t2 = { 0.7f, 1.1f };
    IkResult an = solve_two_bone(c2.root, t2, l2[0], l2[1], c2.pole);
    IkSolveResult ra = ik_solve(&c2, t2, off);
    check("off2_status", ra.status == IK_OK && ra.used_analytic == 1);
    check("off2_pose_identical",
          an.ok &&
          ra.joints[0].x == c2.root.x && ra.joints[0].y == c2.root.y &&
          ra.joints[1].x == an.joint.x && ra.joints[1].y == an.joint.y &&
          ra.joints[2].x == t2.x && ra.joints[2].y == t2.y);
}

static void test_iteration_limit(void) {
    float l4[4] = { 1.0f, 1.0f, 1.0f, 1.0f };
    IkChain c = make_chain(4, (V2){ 0, 0 }, l4);
    V2 target = { 0.5f, 2.0f };

    IkOptions one = ik_default_options();
    one.tolerance = 1e-8f;
    one.max_iterations = 1;
    IkSolveResult r1 = ik_solve(&c, target, one);
    check("iter_limit_status", r1.status == IK_ERR_NOT_CONVERGED);
    check("iter_limit_used", r1.iterations_used == 1);

    IkOptions many = ik_default_options();
    many.tolerance = 1e-8f;
    many.max_iterations = 200;
    IkSolveResult r2 = ik_solve(&c, target, many);
    check("iter_many_converged", r2.status == IK_OK);
    check("iter_many_within_cap", r2.iterations_used > 1 && r2.iterations_used <= 200);
}

/* Compare iterative FABRIK against the analytic two-bone solver over a grid of
   directions and distances.  The reported error is the largest positional
   difference for either the joint or the end effector. */
static float two_bone_grid(const IkOptions* opt, int* cases, int* converged,
                           int* worst_it, int* other_errors) {
    const float l1 = 1.0f, l2 = 0.8f, total = l1 + l2;
    const float fractions[] = { 0.30f, 0.50f, 0.70f, 0.90f, 0.98f };
    const int dirs = 16;
    const float pi = 3.14159265358979f;
    float max_err = 0.0f;
    *cases = 0;
    *converged = 0;
    *worst_it = 0;
    *other_errors = 0;

    for (int d = 0; d < dirs; d++) {
        float th = 2.0f * pi * (float)d / (float)dirs;
        V2 dir = { cosf(th), sinf(th) };
        V2 pole = { -dir.y, dir.x };   /* positive bend branch */
        for (unsigned k = 0; k < sizeof(fractions) / sizeof(fractions[0]); k++) {
            float dist = fractions[k] * total;
            V2 target = vscale(dir, dist);

            IkResult an = solve_two_bone((V2){ 0, 0 }, target, l1, l2, pole);
            if (!an.ok) continue;

            /* seed the elbow on the same branch as the analytic solution */
            float seed_ang = th + 0.9f;
            IkChain c = make_chain(2, (V2){ 0, 0 }, (float[]){ l1, l2 });
            c.pose[1] = (V2){ l1 * cosf(seed_ang), l1 * sinf(seed_ang) };
            c.pose[2] = vadd(c.pose[1], (V2){ l2 * cosf(seed_ang), l2 * sinf(seed_ang) });

            IkSolveResult fab = ik_solve(&c, target, *opt);
            (*cases)++;
            if (fab.status == IK_OK) {
                (*converged)++;
                if (fab.iterations_used > *worst_it) *worst_it = fab.iterations_used;
                V2 an_end = vadd(an.joint, (V2){ l2 * cosf(an.a2), l2 * sinf(an.a2) });
                float e_joint = vdist(fab.joints[1], an.joint);
                float e_end = vdist(fab.joints[2], an_end);
                float e = e_joint > e_end ? e_joint : e_end;
                if (e > max_err) max_err = e;
            } else if (fab.status != IK_ERR_NOT_CONVERGED) {
                (*other_errors)++;
            }
        }
    }
    return max_err;
}

static void test_two_bone_against_analytic(void) {
    IkOptions opt = ik_default_options();
    opt.tolerance = 1e-6f;
    opt.max_iterations = 128;

    int cases = 0, converged = 0, worst_it = 0, other = 0;
    float max_err = two_bone_grid(&opt, &cases, &converged, &worst_it, &other);
    float err_bound = TWO_BONE_ERR_RATIO * opt.tolerance;

    printf("two_bone_cases=%d\n", cases);
    printf("two_bone_converged=%d\n", converged);
    printf("two_bone_max_iterations=%d\n", worst_it);
    printf("two_bone_err_max=%.9g\n", (double)max_err);
    printf("two_bone_err_bound=%.9g\n", (double)err_bound);
    check("two_bone_all_converged", cases > 0 && converged == cases);
    check("two_bone_no_unexpected_errors", other == 0);
    check("two_bone_err_within_bound", max_err <= err_bound);
}

/* Error and iteration count as a function of the tolerance.
   Float32 coordinates put a practical floor near 1e-6; below that the residual
   may stall above the tolerance in the near-straight configuration.  Such a run
   must report IK_ERR_NOT_CONVERGED, never a silently inaccurate OK. */
#define TOL_FLOOR 1e-6f
static void two_bone_tolerance_table(void) {
    const float tols[] = { 1e-3f, 1e-4f, 1e-5f, 1e-6f, 1e-7f };
    int ok = 1, floor_ok = 1;
    printf("tol,max_iterations,converged,cases,max_err,err_over_tol\n");
    for (unsigned i = 0; i < sizeof(tols) / sizeof(tols[0]); i++) {
        IkOptions opt = ik_default_options();
        opt.tolerance = tols[i];
        opt.max_iterations = 4096;
        int cases = 0, converged = 0, worst_it = 0, other = 0;
        float max_err = two_bone_grid(&opt, &cases, &converged, &worst_it, &other);
        printf("%.0e,%d,%d,%d,%.9g,%.2f\n", (double)tols[i], worst_it,
               converged, cases, (double)max_err, (double)(max_err / tols[i]));
        if (other != 0) ok = 0;
        if (tols[i] >= TOL_FLOOR) {
            if (max_err > TWO_BONE_ERR_RATIO * tols[i] || converged != cases) ok = 0;
        } else {
            /* below the float32 floor: converged runs must still respect the
               bound; failed runs must have reported NOT_CONVERGED. */
            if (max_err > TWO_BONE_ERR_RATIO * tols[i]) floor_ok = 0;
        }
    }
    check("two_bone_tolerance_table_bound", ok);
    check("two_bone_below_floor_safe", floor_ok);
}

/* Residual history for the iteration-vs-residual curve. */
static int curve_case(IkSolveResult* out, float* hist, int cap, IkOptions opt) {
    float l4[4] = { 1.0f, 1.0f, 1.0f, 1.0f };
    IkChain c = make_chain(4, (V2){ 0, 0 }, l4);
    c.pose[1] = (V2){ 0.7f, 0.7f };
    c.pose[2] = (V2){ 1.4f, 1.2f };
    c.pose[3] = (V2){ 2.1f, 1.5f };
    c.pose[4] = (V2){ 2.8f, 1.6f };

    *out = fabrik_solve(&c, (V2){ 1.0f, 2.4f }, &opt, hist, cap);

    int n = out->iterations_used;
    int decreasing = n >= 2 && hist[n - 1] < hist[0];
    int last_ok = n > 0 && hist[n - 1] <= opt.tolerance;
    check("curve_converged", out->status == IK_OK);
    check("curve_residual_decreasing", decreasing);
    check("curve_final_below_tolerance", last_ok);
    check("curve_history_len", n <= cap);
    return n;
}

int main(int argc, char** argv) {
    int want_curve = 0;
    float curve_tol = 1e-7f;
    int curve_iters = 200;
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--curve") == 0) want_curve = 1;
        else if (strcmp(argv[i], "--curve-tol") == 0 && i + 1 < argc) curve_tol = strtof(argv[++i], NULL);
        else if (strcmp(argv[i], "--curve-iters") == 0 && i + 1 < argc) curve_iters = (int)strtol(argv[++i], NULL, 10);
    }

    if (want_curve) {
        if (!(curve_tol > 0.0f) || curve_iters <= 0) {
            fprintf(stderr, "bad --curve-tol (needs > 0) or --curve-iters (needs > 0)\n");
            return 2;
        }
        float hist[4096];
        int cap = (int)(sizeof(hist) / sizeof(hist[0]));
        if (curve_iters > cap) curve_iters = cap;
        IkOptions opt = ik_default_options();
        opt.tolerance = curve_tol;
        opt.max_iterations = curve_iters;
        IkSolveResult r;
        int n = curve_case(&r, hist, cap, opt);
        printf("iteration,residual\n");
        for (int i = 0; i < n; i++) printf("%d,%.9g\n", i + 1, (double)hist[i]);
        fprintf(stderr, "curve_converged=%d iterations=%d residual=%.9g\n",
                r.status == IK_OK, r.iterations_used, (double)r.residual);
        return r.status == IK_OK && failures == 0 ? 0 : 1;
    }

    test_analytic();
    test_fabrik_chains();
    test_rejections();
    test_iteration_switch();
    test_iteration_limit();
    test_two_bone_against_analytic();
    two_bone_tolerance_table();

    float hist[256];
    IkSolveResult cr;
    IkOptions curve_opt = ik_default_options();
    curve_opt.tolerance = 1e-7f;
    curve_opt.max_iterations = 200;
    curve_case(&cr, hist, (int)(sizeof(hist) / sizeof(hist[0])), curve_opt);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
