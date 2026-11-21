
/* mass: non-uniform linear density along a straight rope (thick at one end,
   thin at the other). Trapezoidal lumped masses give total mass, center of
   mass and the 2x2 inertia tensor; per-node inverse masses feed constraint
   solvers so the heavy end barely moves and the rope stops flinging.
   Discretization error bounds (vs. analytic linear-density values) and an
   angular velocity drift bound under rotation are provided. */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define MAX_N 32

typedef struct { float x, y; } P2;

typedef struct {
    P2 a, b;            /* endpoints */
    float rho[MAX_N];   /* linear density at each node (mass per length) */
    int n;
    float len;
} Rope;

static void rope_init(Rope* r, float x0, float y0, float x1, float y1) {
    memset(r, 0, sizeof(*r));
    r->a.x = x0; r->a.y = y0;
    r->b.x = x1; r->b.y = y1;
    float dx = x1 - x0, dy = y1 - y0;
    r->len = sqrtf(dx * dx + dy * dy);
}

/* evenly spaced nodes; density must be positive everywhere */
static int rope_set_density(Rope* r, int n, const float* rho) {
    if (n < 2 || n > MAX_N) return -1;
    if (r->len < 1e-6f) return -1;
    for (int i = 0; i < n; i++)
        if (!(rho[i] > 0.0f)) return -1;   /* rejects zero, negative, NaN */
    r->n = n;
    memcpy(r->rho, rho, sizeof(float) * (size_t)n);
    return 0;
}

static P2 rope_node(const Rope* r, int i) {
    float t = (float)i / (float)(r->n - 1);
    P2 p = { r->a.x + (r->b.x - r->a.x) * t, r->a.y + (r->b.y - r->a.y) * t };
    return p;
}

/* trapezoidal lumped masses; returns total mass */
static float rope_masses(const Rope* r, float* m) {
    float h = r->len / (float)(r->n - 1);
    float total = 0.0f;
    for (int i = 0; i < r->n; i++) {
        float w = (i == 0 || i == r->n - 1) ? 0.5f : 1.0f;
        m[i] = r->rho[i] * h * w;
        total += m[i];
    }
    return total;
}

static float rope_com(const Rope* r, P2* out) {
    float m[MAX_N];
    float total = rope_masses(r, m);
    float cx = 0.0f, cy = 0.0f;
    for (int i = 0; i < r->n; i++) {
        P2 p = rope_node(r, i);
        cx += m[i] * p.x;
        cy += m[i] * p.y;
    }
    out->x = cx / total;
    out->y = cy / total;
    return total;
}

/* 2x2 inertia tensor about the center of mass: i[0]=Ixx i[1]=Iyy i[2]=Ixy;
   returns the trace, i.e. the scalar moment about the out-of-plane axis */
static float rope_inertia(const Rope* r, float i[3]) {
    float m[MAX_N];
    rope_masses(r, m);
    P2 c;
    rope_com(r, &c);
    float ixx = 0.0f, iyy = 0.0f, ixy = 0.0f;
    for (int k = 0; k < r->n; k++) {
        P2 p = rope_node(r, k);
        float dx = p.x - c.x, dy = p.y - c.y;
        ixx += m[k] * dy * dy;
        iyy += m[k] * dx * dx;
        ixy -= m[k] * dx * dy;
    }
    i[0] = ixx; i[1] = iyy; i[2] = ixy;
    return ixx + iyy;
}

/* inverse masses for constraint weighting: heavy end moves less */
static float rope_inv_masses(const Rope* r, float* out) {
    float m[MAX_N];
    float total = rope_masses(r, m);
    for (int i = 0; i < r->n; i++) out[i] = 1.0f / m[i];
    return total;
}

/* trapezoid rule error bound: |E| <= L h^2 / 12 * max|f''|.
   Bounds below assume density varies linearly between endpoint densities. */

/* first moment int rho*s ds: f'' = 2k for rho(s) = rho0 + k s */
static float mass_first_moment_bound(float len, int n, float k) {
    double h = (double)len / (double)(n - 1);
    return (float)((double)len * h * h * fabs((double)k) / 6.0);
}

/* second moment int rho*s^2 ds: f'' = 2 rho0 + 6 k s, max |f''| at an endpoint */
static float mass_second_moment_bound(float len, int n, float rho0, float k) {
    double h = (double)len / (double)(n - 1);
    double fa = fabs(2.0 * (double)rho0);
    double fb = fabs(2.0 * (double)rho0 + 6.0 * (double)k * (double)len);
    double f2 = fa > fb ? fa : fb;
    return (float)((double)len * h * h * f2 / 12.0);
}

/* center-of-mass error bound along the rope axis (mass is exact for linear
   density, so the com error is the first-moment error over M) */
static float mass_com_bound(float len, int n, float rho0, float rho1) {
    double k = ((double)rho1 - (double)rho0) / (double)len;
    double M = (double)len * ((double)rho0 + (double)rho1) / 2.0;
    return (float)((double)mass_first_moment_bound(len, n, (float)k) / M);
}

/* inertia (about com) error bound: I = mu2 - M xbar^2 with M exact, so
   dI <= E_mu2 + 2 M xbar E_x + M E_x^2 */
static float mass_inertia_bound(float len, int n, float rho0, float rho1) {
    double k = ((double)rho1 - (double)rho0) / (double)len;
    double M = (double)len * ((double)rho0 + (double)rho1) / 2.0;
    double xbar = (double)len * ((double)rho0 + 2.0 * (double)rho1) /
                  (3.0 * ((double)rho0 + (double)rho1));
    double ex = (double)mass_com_bound(len, n, rho0, rho1);
    double e2 = (double)mass_second_moment_bound(len, n, rho0, (float)k);
    return (float)(e2 + 2.0 * M * xbar * ex + M * ex * ex);
}

/* angular velocity drift bound: with angular momentum L = I w conserved and
   inertia known to within delta_i, |dw| <= |w| delta_i / (I - delta_i);
   returns -1 if delta_i >= inertia (inertia too uncertain to bound) */
static float mass_omega_drift_bound(float inertia, float delta_i, float omega) {
    if (!(delta_i >= 0.0f) || delta_i >= inertia) return -1.0f;
    return fabsf(omega) * delta_i / (inertia - delta_i);
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    /* thick at s=0 (rho=2), thin at s=L (rho=0.5), L=4, along the x axis */
    const int n = 17;
    const float L = 4.0f, r0 = 2.0f, r1 = 0.5f;
    float rho[MAX_N];
    for (int i = 0; i < n; i++)
        rho[i] = r0 + (r1 - r0) * (float)i / (float)(n - 1);

    Rope r;
    rope_init(&r, 0.0f, 0.0f, L, 0.0f);
    check("set_ok", rope_set_density(&r, n, rho) == 0);

    float bad_zero[4] = { 1.0f, 0.0f, 1.0f, 1.0f };
    float bad_neg[4] = { 1.0f, -0.5f, 1.0f, 1.0f };
    check("reject_zero", rope_set_density(&r, 4, bad_zero) == -1);
    check("reject_negative", rope_set_density(&r, 4, bad_neg) == -1);
    check("reject_few_nodes", rope_set_density(&r, 1, rho) == -1);
    check("reject_too_many", rope_set_density(&r, MAX_N + 1, rho) == -1);
    Rope degenerate;
    rope_init(&degenerate, 1.0f, 1.0f, 1.0f, 1.0f);
    check("reject_zero_len", rope_set_density(&degenerate, 4, rho) == -1);
    check("state_kept_after_reject", r.n == n);

    /* analytic values for linear density on a straight rope */
    float M_ana = L * (r0 + r1) * 0.5f;
    float x_ana = L * (r0 + 2.0f * r1) / (3.0f * (r0 + r1));
    float mu2_ana = L * L * L * (r0 + 3.0f * r1) / 12.0f;
    float I_ana = mu2_ana - M_ana * x_ana * x_ana;

    float m[MAX_N];
    float M_num = rope_masses(&r, m);
    check("mass_matches_analytic", fabsf(M_num - M_ana) < 1e-4f);

    float sum = 0.0f;
    for (int i = 0; i < n; i++) sum += m[i];
    check("mass_conserved", fabsf(sum - M_num) < 1e-5f);

    P2 c;
    rope_com(&r, &c);
    float com_bound = mass_com_bound(L, n, r0, r1);
    check("com_toward_thick_end", c.x < L * 0.5f);
    check("com_on_axis", fabsf(c.y) < 1e-6f);
    check("com_within_bound", fabsf(c.x - x_ana) <= com_bound + 1e-5f);

    float t[3];
    float I_num = rope_inertia(&r, t);
    float i_bound = mass_inertia_bound(L, n, r0, r1);
    check("inertia_within_bound", fabsf(I_num - I_ana) <= i_bound + 1e-4f);
    check("tensor_axis_aligned", fabsf(t[0]) < 1e-6f && fabsf(t[2]) < 1e-6f);
    check("tensor_trace_is_rod", fabsf(t[1] - I_num) < 1e-5f);

    /* tilted rope: same densities, direction (0.6, 0.8), length 5 */
    Rope tilted;
    rope_init(&tilted, 0.0f, 0.0f, 3.0f, 4.0f);
    check("tilted_set_ok", rope_set_density(&tilted, n, rho) == 0);
    float tt[3];
    float I_tilt = rope_inertia(&tilted, tt);
    float L5 = 5.0f;
    float M5 = L5 * (r0 + r1) * 0.5f;
    float x5 = L5 * (r0 + 2.0f * r1) / (3.0f * (r0 + r1));
    float I5_ana = L5 * L5 * L5 * (r0 + 3.0f * r1) / 12.0f - M5 * x5 * x5;
    check("tilted_trace_within_bound",
          fabsf(I_tilt - I5_ana) <= mass_inertia_bound(L5, n, r0, r1) + 1e-4f);
    /* a thin rod has no inertia about its own axis: I * u == 0 */
    float ux = 0.6f, uy = 0.8f;
    float rx = tt[0] * ux + tt[2] * uy;
    float ry = tt[2] * ux + tt[1] * uy;
    check("tensor_kills_axis", fabsf(rx) < 1e-3f && fabsf(ry) < 1e-3f);
    check("tensor_psd", tt[0] >= -1e-6f && tt[1] >= -1e-6f &&
          tt[0] * tt[1] - tt[2] * tt[2] >= -1e-4f);

    /* uniform density sanity: com at midpoint, I = M L^2 / 12 */
    Rope uni;
    rope_init(&uni, 0.0f, 0.0f, 2.0f, 0.0f);
    float one[9];
    for (int i = 0; i < 9; i++) one[i] = 1.0f;
    rope_set_density(&uni, 9, one);
    P2 uc;
    float Mu = rope_com(&uni, &uc);
    float ut[3];
    float Iu = rope_inertia(&uni, ut);
    check("uniform_midpoint", fabsf(uc.x - 1.0f) < 1e-5f);
    check("uniform_inertia", fabsf(Iu - Mu * 4.0f / 12.0f) <=
          mass_inertia_bound(2.0f, 9, 1.0f, 1.0f) + 1e-5f);

    /* correction weights: inverse mass splits a constraint correction so
       the thin (light) end moves more than the thick (heavy) end */
    float inv[MAX_N];
    rope_inv_masses(&r, inv);
    check("heavy_end_stiffer", inv[0] < inv[n - 1]);
    float wsum = inv[0] + inv[n - 1];
    check("thin_end_moves_more", inv[n - 1] / wsum > 0.5f);
    check("split_matches_mass_ratio",
          fabsf(inv[n - 1] / inv[0] - m[0] / m[n - 1]) < 1e-4f);

    /* angular velocity drift: L_ang = I w conserved, inertia off by <= dI */
    float omega = 3.0f;
    float L_ang = I_ana * omega;
    float omega_num = L_ang / I_num;
    float drift_bound = mass_omega_drift_bound(I_ana, i_bound, omega);
    check("drift_bound_positive", drift_bound > 0.0f);
    check("drift_within_bound", fabsf(omega_num - omega) <= drift_bound + 1e-4f);
    check("drift_bound_invalid", mass_omega_drift_bound(1.0f, 1.5f, 3.0f) < 0.0f);

    /* refining the mesh must shrink the bounds */
    check("refinement_shrinks_bounds",
          mass_inertia_bound(L, 33, r0, r1) < i_bound &&
          mass_com_bound(L, 33, r0, r1) < com_bound);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
