
/* order: pipeline order comparison between skin->constraint (A) and
   constraint->skin (B). defines both orders, measures condition numbers
   (skinning map and constraint jacobian) and per-sweep error growth,
   gives a per-scenario selection basis (error when there is none),
   bounds the joint position difference between the two orders, and
   reports when a condition number crosses the instability threshold. */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define MAX_V 8
#define MAX_B 8
#define MAX_IT 32

/* float has ~7 digits; above KAPPA_WARN the solve loses meaningful
   precision, above KAPPA_UNSTABLE (~1/sqrt(eps)^2 territory) the
   linearized solve is dominated by roundoff and is reported unstable. */
#define KAPPA_WARN 1e3f
#define KAPPA_UNSTABLE 1e6f

/* ---------- rig and skinning ---------- */

typedef struct { float x, y, ang; } Pose;   /* absolute 2d transform */

typedef struct {
    Pose bind[MAX_B];
    Pose pose[MAX_B];
    int n;
} Rig;

typedef struct {
    float bx, by;       /* bind position */
    float w[MAX_B];     /* bone weights, not necessarily normalized */
} Vert;

static void rot(float x, float y, float a, float* ox, float* oy) {
    float c = cosf(a), s = sinf(a);
    *ox = x * c - y * s; *oy = x * s + y * c;
}

/* skin a bind point: sum_b w_b * (T_pose_b * T_bind_b^-1) * p */
static void skin_point(const Rig* r, const Vert* v, float* ox, float* oy) {
    float ax = 0.0f, ay = 0.0f;
    for (int b = 0; b < r->n; b++) {
        float w = v->w[b];
        if (w == 0.0f) continue;
        float rel = r->pose[b].ang - r->bind[b].ang;
        float dx = v->bx - r->bind[b].x, dy = v->by - r->bind[b].y;
        float rx, ry;
        rot(dx, dy, rel, &rx, &ry);
        ax += w * (r->pose[b].x + rx);
        ay += w * (r->pose[b].y + ry);
    }
    *ox = ax; *oy = ay;
}

/* the skinning jacobian of one vertex is sum_b w_b R(rel_b), which in 2d
   is always a similarity; both singular values equal this scale. */
static float skin_scale(const Rig* r, const Vert* v) {
    float re = 0.0f, im = 0.0f;
    for (int b = 0; b < r->n; b++) {
        float a = r->pose[b].ang - r->bind[b].ang;
        re += v->w[b] * cosf(a);
        im += v->w[b] * sinf(a);
    }
    return sqrtf(re * re + im * im);
}

/* global skinning condition number: max/min vertex scale. unnormalized
   weights or phase cancellation between bone rotations degrade it. */
static float skin_cond(const Rig* r, const Vert* vs, int n) {
    float lo = INFINITY, hi = 0.0f;
    for (int i = 0; i < n; i++) {
        float s = skin_scale(r, &vs[i]);
        if (s < lo) lo = s;
        if (s > hi) hi = s;
    }
    if (lo <= 1e-12f) return INFINITY;
    return hi / lo;
}

static float skin_max_scale(const Rig* r, const Vert* vs, int n) {
    float hi = 0.0f;
    for (int i = 0; i < n; i++) {
        float s = skin_scale(r, &vs[i]);
        if (s > hi) hi = s;
    }
    return hi;
}

/* ---------- distance constraint chain ---------- */

typedef struct {
    float x[MAX_V], y[MAX_V];
    float rest[MAX_V];
    float inv_mass[MAX_V];
    int n;
} Chain;

typedef struct {
    float resid[MAX_IT];    /* max |C_i| after each sweep: error growth */
    int iters;
    float disp;             /* max per-point displacement over the solve */
    float final_resid;
} SolveStats;

static float proj_sweep(Chain* c) {
    for (int i = 0; i < c->n - 1; i++) {
        float dx = c->x[i + 1] - c->x[i], dy = c->y[i + 1] - c->y[i];
        float d = sqrtf(dx * dx + dy * dy);
        if (d < 1e-6f) continue;
        float diff = (d - c->rest[i]) / d;
        float wa = c->inv_mass[i], wb = c->inv_mass[i + 1];
        float ws = wa + wb;
        if (ws <= 0.0f) continue;
        c->x[i]     += dx * diff * (wa / ws);
        c->y[i]     += dy * diff * (wa / ws);
        c->x[i + 1] -= dx * diff * (wb / ws);
        c->y[i + 1] -= dy * diff * (wb / ws);
    }
    float worst = 0.0f;
    for (int i = 0; i < c->n - 1; i++) {
        float dx = c->x[i + 1] - c->x[i], dy = c->y[i + 1] - c->y[i];
        float e = fabsf(sqrtf(dx * dx + dy * dy) - c->rest[i]);
        if (e > worst) worst = e;
    }
    return worst;
}

static SolveStats chain_solve(Chain* c, int max_it) {
    SolveStats st;
    memset(&st, 0, sizeof(st));
    float sx[MAX_V], sy[MAX_V];
    memcpy(sx, c->x, sizeof(sx));
    memcpy(sy, c->y, sizeof(sy));
    for (int it = 0; it < max_it && it < MAX_IT; it++) {
        st.resid[it] = proj_sweep(c);
        st.iters++;
        if (st.resid[it] < 1e-6f) break;
    }
    st.final_resid = st.resid[st.iters - 1];
    float disp = 0.0f;
    for (int i = 0; i < c->n; i++) {
        float d = hypotf(c->x[i] - sx[i], c->y[i] - sy[i]);
        if (d > disp) disp = d;
    }
    st.disp = disp;
    return st;
}

/* worst per-sweep residual ratio; >1 means the error grows */
static float error_growth(const SolveStats* st) {
    float worst = 0.0f;
    for (int i = 1; i < st->iters; i++) {
        if (st->resid[i - 1] > 1e-12f) {
            float r = st->resid[i] / st->resid[i - 1];
            if (r > worst) worst = r;
        }
    }
    return worst;
}

/* ---------- constraint jacobian condition number ---------- */

/* jacobi eigenvalue iteration for a symmetric m x m matrix */
static void jacobi_eigen(float a[MAX_V][MAX_V], int m, float* lmin, float* lmax) {
    for (int sweep = 0; sweep < 64; sweep++) {
        float off = 0.0f;
        for (int p = 0; p < m; p++)
            for (int q = p + 1; q < m; q++) off += a[p][q] * a[p][q];
        if (off < 1e-12f) break;
        for (int p = 0; p < m - 1; p++) {
            for (int q = p + 1; q < m; q++) {
                float apq = a[p][q];
                if (fabsf(apq) < 1e-12f) continue;
                float theta = 0.5f * atan2f(2.0f * apq, a[q][q] - a[p][p]);
                float c = cosf(theta), s = sinf(theta);
                for (int k = 0; k < m; k++) {
                    float akp = a[k][p], akq = a[k][q];
                    a[k][p] = c * akp - s * akq;
                    a[k][q] = s * akp + c * akq;
                }
                for (int k = 0; k < m; k++) {
                    float apk = a[p][k], aqk = a[q][k];
                    a[p][k] = c * apk - s * aqk;
                    a[q][k] = s * apk + c * aqk;
                }
            }
        }
    }
    *lmin = *lmax = a[0][0];
    for (int i = 1; i < m; i++) {
        if (a[i][i] < *lmin) *lmin = a[i][i];
        if (a[i][i] > *lmax) *lmax = a[i][i];
    }
}

/* condition number of the constraint jacobian via G = J J^T (tridiagonal
   for a chain). stretches toward a straight line -> kappa grows like
   O(n^2); a degenerate (zero-length) segment is infinitely ill-conditioned. */
static float chain_cond(const Chain* c) {
    int m = c->n - 1;
    if (m < 1) return 1.0f;
    float ux[MAX_V], uy[MAX_V];
    for (int i = 0; i < m; i++) {
        float dx = c->x[i + 1] - c->x[i], dy = c->y[i + 1] - c->y[i];
        float d = sqrtf(dx * dx + dy * dy);
        if (d < 1e-6f) return INFINITY;
        ux[i] = dx / d; uy[i] = dy / d;
    }
    float g[MAX_V][MAX_V];
    memset(g, 0, sizeof(g));
    for (int i = 0; i < m; i++) g[i][i] = 2.0f;
    for (int i = 0; i + 1 < m; i++) {
        float dot = -(ux[i] * ux[i + 1] + uy[i] * uy[i + 1]);
        g[i][i + 1] = dot;
        g[i + 1][i] = dot;
    }
    float lmin, lmax;
    jacobi_eigen(g, m, &lmin, &lmax);
    if (lmin <= 1e-12f) return INFINITY;
    return lmax / lmin;
}

/* 0 ok, 1 warn, 2 unstable; prints the crossed threshold */
static int kappa_report(const char* tag, float k) {
    if (!(k < KAPPA_UNSTABLE)) {    /* catches inf and nan */
        printf("kappa_unstable %s k=%g threshold=%g\n", tag, (double)k, (double)KAPPA_UNSTABLE);
        return 2;
    }
    if (!(k < KAPPA_WARN)) {
        printf("kappa_warn %s k=%g threshold=%g\n", tag, (double)k, (double)KAPPA_WARN);
        return 1;
    }
    return 0;
}

/* ---------- the two orders ---------- */

/* order A: skin to world space, then solve constraints there */
static void run_order_a(const Rig* r, const Vert* vs, int n, const float* rest,
                        int iters, Chain* out, SolveStats* st) {
    Chain c;
    memset(&c, 0, sizeof(c));
    c.n = n;
    for (int i = 0; i < n; i++) {
        skin_point(r, &vs[i], &c.x[i], &c.y[i]);
        c.rest[i] = rest[i];
        c.inv_mass[i] = 1.0f;
    }
    *st = chain_solve(&c, iters);
    *out = c;
}

/* order B: solve constraints in bind space, then skin the result */
static void run_order_b(const Rig* r, const Vert* vs, int n, const float* rest,
                        int iters, Chain* out, SolveStats* st) {
    Chain c;
    memset(&c, 0, sizeof(c));
    c.n = n;
    for (int i = 0; i < n; i++) {
        c.x[i] = vs[i].bx;
        c.y[i] = vs[i].by;
        c.rest[i] = rest[i];
        c.inv_mass[i] = 1.0f;
    }
    *st = chain_solve(&c, iters);
    for (int i = 0; i < n; i++)
        skin_point(r, &vs[i], &out->x[i], &out->y[i]);
    out->n = n;
}

/* a-posteriori upper bound on the per-joint difference between orders:
   |A - B| <= |A - S(v)| + |S(v) - S(P_b(v))|
          <= (dispA + kc*residA) + smax*(dispB + kc*residB),
   where kc*resid bounds the remaining distance to the exact projection. */
static float order_diff_bound(const SolveStats* a, const SolveStats* b,
                              float smax, float kc) {
    return a->disp + kc * a->final_resid
         + smax * (b->disp + kc * b->final_resid);
}

/* ---------- order selection ---------- */

enum {
    ORDER_ERR_UNSTABLE = -2,
    ORDER_ERR_NO_BASIS = -1,
    ORDER_SKIN_FIRST = 1,        /* skin -> constraint */
    ORDER_CONSTRAINT_FIRST = 2   /* constraint -> skin */
};

typedef struct {
    int known_world_anchor;     /* is it known whether constraints hold in world space? */
    int world_anchor;           /* constraints reference world-space anchors/colliders */
    int known_weights;          /* is weight normalization known? */
    int weights_normalized;
    float kappa_skin;           /* <=0 when unknown */
    float kappa_constraint;     /* <=0 when unknown */
} OrderScenario;

static int order_select(const OrderScenario* sc, const char** why) {
    float ks = sc->kappa_skin > 0.0f ? sc->kappa_skin : 1.0f;
    float kc = sc->kappa_constraint > 0.0f ? sc->kappa_constraint : 1.0f;
    if (!(ks < KAPPA_UNSTABLE) || !(kc < KAPPA_UNSTABLE)) {
        *why = "condition number beyond instability threshold";
        return ORDER_ERR_UNSTABLE;
    }
    if (sc->known_world_anchor && sc->world_anchor) {
        *why = "constraints hold in world space: solve after skinning";
        return ORDER_SKIN_FIRST;
    }
    if (sc->known_weights && !sc->weights_normalized) {
        *why = "unnormalized weights distort world-space rest lengths: solve in bind space";
        return ORDER_CONSTRAINT_FIRST;
    }
    if (!sc->known_world_anchor && !sc->known_weights &&
        sc->kappa_skin <= 0.0f && sc->kappa_constraint <= 0.0f) {
        *why = "no scenario facts: no basis for order selection";
        return ORDER_ERR_NO_BASIS;
    }
    *why = "well-conditioned and no world-space facts: orders equivalent within solver tolerance, skin first";
    return ORDER_SKIN_FIRST;
}

/* ---------- tests ---------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static Vert mkvert(float bx, float by, float w0, float w1) {
    Vert v;
    memset(&v, 0, sizeof(v));
    v.bx = bx; v.by = by;
    v.w[0] = w0; v.w[1] = w1;
    return v;
}

static float joint_diff(const Chain* a, const Chain* b) {
    float worst = 0.0f;
    for (int i = 0; i < a->n; i++) {
        float d = hypotf(a->x[i] - b->x[i], a->y[i] - b->y[i]);
        if (d > worst) worst = d;
    }
    return worst;
}

static int near(float a, float b, float tol) { return fabsf(a - b) < tol; }

int main(void) {
    /* hand-computed scenario: bone0 fixed at origin, bone1 bound at (1,0);
       v0=(0,0) on bone0, v1=(2,0) on bone1, rest length 2.
       frame 1: bone1 rotated 90 deg -> skinned v1=(1,1), dist sqrt(2).
       order A pulls both joints to distance 2 along (1,1)/sqrt(2);
       order B sees rest satisfied in bind space and does nothing. */
    Rig rig;
    memset(&rig, 0, sizeof(rig));
    rig.n = 2;
    rig.bind[0] = (Pose){ 0.0f, 0.0f, 0.0f };
    rig.pose[0] = (Pose){ 0.0f, 0.0f, 0.0f };
    rig.bind[1] = (Pose){ 1.0f, 0.0f, 0.0f };

    Vert vs[2];
    vs[0] = mkvert(0.0f, 0.0f, 1.0f, 0.0f);
    vs[1] = mkvert(2.0f, 0.0f, 0.0f, 1.0f);
    float rest[1] = { 2.0f };

    Chain ca, cb;
    SolveStats sa, sb;

    /* frame 1, hand-computed: A joints at (0.5-+1/sqrt(2)), B at (0,0),(1,1),
       per-joint difference 1 - 1/sqrt(2) = 0.29289322 */
    rig.pose[1] = (Pose){ 1.0f, 0.0f, 1.5707963f };
    run_order_a(&rig, vs, 2, rest, 8, &ca, &sa);
    run_order_b(&rig, vs, 2, rest, 8, &cb, &sb);
    check("f1_order_a_j0", near(ca.x[0], -0.20710678f, 5e-4f) && near(ca.y[0], -0.20710678f, 5e-4f));
    check("f1_order_a_j1", near(ca.x[1], 1.20710678f, 5e-4f) && near(ca.y[1], 1.20710678f, 5e-4f));
    check("f1_order_b", near(cb.x[0], 0.0f, 5e-4f) && near(cb.x[1], 1.0f, 5e-4f) && near(cb.y[1], 1.0f, 5e-4f));
    float d1 = joint_diff(&ca, &cb);
    check("f1_diff_hand", near(d1, 0.29289322f, 5e-4f));
    float kc1 = chain_cond(&ca);
    float bound1 = order_diff_bound(&sa, &sb, skin_max_scale(&rig, vs, 2), kc1);
    check("f1_diff_bound", d1 <= bound1 + 1e-5f);

    /* frame 2, hand-computed: bone1 at 60 deg -> skinned v1=(1.5, sqrt(3)/2),
       dist sqrt(3); A joints at (3/4-+sqrt(3)/2, sqrt(3)/4-+1/2)... i.e.
       j0=(-0.11602540,-0.06698730), j1=(1.61602540,0.93301270);
       difference sqrt(7/4 - sqrt(3)) = 0.1339746 */
    rig.pose[1] = (Pose){ 1.0f, 0.0f, 1.0471975f };
    run_order_a(&rig, vs, 2, rest, 8, &ca, &sa);
    run_order_b(&rig, vs, 2, rest, 8, &cb, &sb);
    check("f2_order_a_j0", near(ca.x[0], -0.11602540f, 5e-4f) && near(ca.y[0], -0.06698730f, 5e-4f));
    check("f2_order_a_j1", near(ca.x[1], 1.61602540f, 5e-4f) && near(ca.y[1], 0.93301270f, 5e-4f));
    check("f2_order_b", near(cb.x[1], 1.5f, 5e-4f) && near(cb.y[1], 0.86602540f, 5e-4f));
    float d2 = joint_diff(&ca, &cb);
    check("f2_diff_hand", near(d2, 0.1339746f, 5e-4f));
    float bound2 = order_diff_bound(&sa, &sb, skin_max_scale(&rig, vs, 2), chain_cond(&ca));
    check("f2_diff_bound", d2 <= bound2 + 1e-5f);

    /* the two frames differ visually: orders genuinely do not commute */
    check("orders_differ", d1 > 0.1f && d2 > 0.05f);

    /* skinning condition numbers */
    check("kappa_skin_rigid", near(skin_cond(&rig, vs, 2), 1.0f, 1e-4f));
    Vert half[2];
    half[0] = mkvert(0.0f, 0.0f, 1.0f, 0.0f);
    half[1] = mkvert(2.0f, 0.0f, 0.0f, 0.5f);   /* weights sum to 0.5 */
    check("kappa_skin_unnormalized", near(skin_cond(&rig, half, 2), 2.0f, 1e-3f));

    /* phase cancellation: rotations 0 and 180 deg with equal weights
       shrink one vertex's scale toward zero, so the max/min ratio
       across vertices blows past the instability threshold */
    Rig cancel;
    memset(&cancel, 0, sizeof(cancel));
    cancel.n = 2;
    cancel.pose[1].ang = 3.14159265f;
    Vert cv[2];
    cv[0] = mkvert(0.0f, 0.0f, 0.5f, 0.5f);   /* scales cancel */
    cv[1] = mkvert(0.0f, 0.0f, 1.0f, 0.0f);   /* rigid, scale 1 */
    float kc_cancel = skin_cond(&cancel, cv, 2);
    check("kappa_skin_cancellation", kc_cancel > KAPPA_UNSTABLE);
    check("kappa_report_unstable", kappa_report("skin_cancel", kc_cancel) == 2);
    check("kappa_report_ok", kappa_report("skin_rigid", 1.0f) == 0);

    /* constraint jacobian conditioning: straight 8-point chain has
       kappa = (1-cos(7pi/8))/(1-cos(pi/8)) ~= 25.27; a 45 deg zigzag
       gives orthogonal consecutive rows -> kappa = 1 */
    Chain straight;
    memset(&straight, 0, sizeof(straight));
    straight.n = 8;
    for (int i = 0; i < 8; i++) { straight.x[i] = (float)i; straight.rest[i] = 1.0f; straight.inv_mass[i] = 1.0f; }
    float k_straight = chain_cond(&straight);
    printf("kappa_chain_straight=%.2f\n", (double)k_straight);
    check("kappa_chain_straight", k_straight > 20.0f && k_straight < 30.0f);

    Chain zig;
    memset(&zig, 0, sizeof(zig));
    zig.n = 8;
    for (int i = 0; i < 8; i++) { zig.x[i] = (float)i; zig.y[i] = (float)(i % 2); zig.rest[i] = 1.0f; zig.inv_mass[i] = 1.0f; }
    float k_zig = chain_cond(&zig);
    printf("kappa_chain_zigzag=%.2f\n", (double)k_zig);
    check("kappa_chain_zigzag", k_zig < k_straight && near(k_zig, 1.0f, 0.5f));

    /* degenerate segment -> infinite condition number -> reported */
    Chain degen;
    memset(&degen, 0, sizeof(degen));
    degen.n = 3;
    degen.rest[0] = 1.0f; degen.rest[1] = 1.0f;
    check("kappa_chain_degenerate", chain_cond(&degen) == INFINITY);
    check("kappa_report_degenerate", kappa_report("chain_degen", chain_cond(&degen)) == 2);

    /* error growth: projection sweeps must not amplify the residual */
    Chain stretched;
    memset(&stretched, 0, sizeof(stretched));
    stretched.n = 6;
    for (int i = 0; i < 6; i++) { stretched.x[i] = 1.5f * i; stretched.rest[i] = 1.0f; stretched.inv_mass[i] = 1.0f; }
    SolveStats sg = chain_solve(&stretched, 16);
    float growth = error_growth(&sg);
    printf("error_growth_stretched=%.3f\n", (double)growth);
    check("error_not_growing", growth <= 1.0f + 1e-6f);
    check("error_decays", sg.final_resid < sg.resid[0]);

    /* selection basis */
    const char* why = 0;
    OrderScenario empty;
    memset(&empty, 0, sizeof(empty));
    check("select_no_basis", order_select(&empty, &why) == ORDER_ERR_NO_BASIS);

    OrderScenario anchor;
    memset(&anchor, 0, sizeof(anchor));
    anchor.known_world_anchor = 1;
    anchor.world_anchor = 1;
    check("select_world_anchor", order_select(&anchor, &why) == ORDER_SKIN_FIRST);

    OrderScenario unnorm;
    memset(&unnorm, 0, sizeof(unnorm));
    unnorm.known_weights = 1;
    unnorm.weights_normalized = 0;
    check("select_unnormalized", order_select(&unnorm, &why) == ORDER_CONSTRAINT_FIRST);

    OrderScenario wild;
    memset(&wild, 0, sizeof(wild));
    wild.known_world_anchor = 1;
    wild.world_anchor = 1;
    wild.kappa_skin = INFINITY;
    check("select_unstable", order_select(&wild, &why) == ORDER_ERR_UNSTABLE);

    OrderScenario calm;
    memset(&calm, 0, sizeof(calm));
    calm.known_weights = 1;
    calm.weights_normalized = 1;
    calm.kappa_skin = 1.0f;
    calm.kappa_constraint = 4.0f;
    check("select_default", order_select(&calm, &why) == ORDER_SKIN_FIRST);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
