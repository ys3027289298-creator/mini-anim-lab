
/* grasp: coupled five-finger grasp solve.
 *
 * Each finger is a two-bone chain solved analytically as a triangle (law of
 * cosines); the elbow takes a fixed, explicitly chosen side of the root-target
 * line, so a finger is a pure function of its contact angle. A single finger
 * therefore needs no iteration and the solver never falls into the mirrored
 * elbow branch by floating-point noise.
 *
 * The coupling is what needs solving: all five contact angles are adjusted
 * together so that
 *   - every finger stays inside its reach annulus,
 *   - every finger pair keeps capsule clearance (distance between bone spines),
 *   - adjacent fingers share a bending bound,
 *   - each contact stays close to its nominal angle.
 * That joint problem is an augmented Lagrangian around a damped least-squares
 * inner loop. When the tolerance cannot be met the solver reports the residual
 * and flags the result invalid instead of returning a near-miss pose.
 *
 * A second, closed-form step distributes a required contact wrench over the
 * fingers in proportion to per-finger weights, and reports the wrench residual
 * when proportional forces cannot realise it.
 */
#include <stdio.h>
#include <string.h>
#include <math.h>

#define NF      5           /* fingers */
#define NPAIR   10          /* unordered finger pairs */
#define NADJ    4           /* adjacent finger pairs */
#define NG      (2 * NF + NPAIR + NADJ + NF)   /* 29 constraints */
#define NV      NF          /* decision variables: contact angles */

typedef struct { double x, y; } V2;

static V2 v2(double x, double y) { V2 r = { x, y }; return r; }
static V2 vsub(V2 a, V2 b) { return v2(a.x - b.x, a.y - b.y); }
static V2 vadd(V2 a, V2 b) { return v2(a.x + b.x, a.y + b.y); }
static V2 vmul(V2 a, double s) { return v2(a.x * s, a.y * s); }
static double vlen(V2 a) { return sqrt(a.x * a.x + a.y * a.y); }
static double vdot(V2 a, V2 b) { return a.x * b.x + a.y * b.y; }
static double vcross(V2 a, V2 b) { return a.x * b.y - a.y * b.x; }
static double clampd(double x, double lo, double hi) { return x < lo ? lo : (x > hi ? hi : x); }

/* ---------- geometry ---------- */

/* Closest distance from point p to segment ab. */
static double point_seg(V2 p, V2 a, V2 b) {
    V2 ab = vsub(b, a), ap = vsub(p, a);
    double l2 = vdot(ab, ab);
    if (l2 < 1e-14) return vlen(ap);
    double t = clampd(vdot(ap, ab) / l2, 0.0, 1.0);
    return vlen(vsub(p, vadd(a, vmul(ab, t))));
}

/* Proper interior crossing test (shared endpoints are not a crossing). */
static int seg_cross(V2 a, V2 b, V2 c, V2 d) {
    double d1 = vcross(vsub(d, c), vsub(a, c));
    double d2 = vcross(vsub(d, c), vsub(b, c));
    double d3 = vcross(vsub(b, a), vsub(c, a));
    double d4 = vcross(vsub(b, a), vsub(d, a));
    return ((d1 > 0) != (d2 > 0)) && ((d3 > 0) != (d4 > 0));
}

/* Distance between two segments; a proper crossing is reported as 0. */
static double seg_seg(V2 a, V2 b, V2 c, V2 d) {
    if (seg_cross(a, b, c, d)) return 0.0;
    double m = point_seg(a, c, d);
    double v = point_seg(b, c, d); if (v < m) m = v;
    v = point_seg(c, a, b);        if (v < m) m = v;
    v = point_seg(d, a, b);        if (v < m) m = v;
    return m;
}

/* One finger: root, its two segments [P,J] and [J,T]. */
typedef struct {
    V2 P, J, T;
    double bend;        /* interior bend: pi - angle between the two bones */
    double dist;        /* root-to-target distance */
    double dmin, dmax;  /* reach annulus */
    int ok;
} Finger;

/* Analytic two-bone solve with an explicit elbow side (+1 / -1). The joint is
 * placed at angle(root->target) + side * acos(...) from the root, which is
 * exactly one of the two triangle solutions and never depends on noise. */
static Finger solve_finger_side(V2 P, V2 T, double l1, double l2, int side) {
    Finger f;
    memset(&f, 0, sizeof(f));
    f.P = P; f.T = T; f.dmin = fabs(l1 - l2); f.dmax = l1 + l2;
    if (l1 <= 0 || l2 <= 0) return f;
    V2 d = vsub(T, P);
    double dist = vlen(d);
    f.dist = dist;
    if (dist < 1e-9) return f;
    if (dist > l1 + l2 + 1e-9 || dist < fabs(l1 - l2) - 1e-9) return f;

    double ca = clampd((l1 * l1 + dist * dist - l2 * l2) / (2 * l1 * dist), -1.0, 1.0);
    double a = acos(ca);
    double ang = atan2(d.y, d.x) + (side >= 0 ? a : -a);
    f.J = v2(P.x + l1 * cos(ang), P.y + l1 * sin(ang));

    V2 w = vsub(P, f.J), v = vsub(T, f.J);
    double lw = vlen(w), lv = vlen(v);
    f.bend = 3.14159265358979323846 - acos(clampd(vdot(w, v) / (lw * lv), -1.0, 1.0));
    f.ok = 1;
    return f;
}

/* Which elbow side points away from the grasp axis. Fingers left of centre
 * (root x < 0) take the branch whose elbow is further left; right side mirrors
 * it; the centre finger takes +1. At the symmetric start both branches sit at
 * the same radius, so an explicit rule is required rather than a cross product
 * that is exactly zero there. */
static int outboard_side(V2 P, V2 T, double l1, double l2) {
    Finger a = solve_finger_side(P, T, l1, l2, +1);
    Finger b = solve_finger_side(P, T, l1, l2, -1);
    if (!a.ok || !b.ok) return +1;
    double rx = P.x;
    if (fabs(rx) < 1e-9) return +1;
    double da = (a.J.x - P.x) * rx;
    double db = (b.J.x - P.x) * rx;
    return da >= db ? +1 : -1;
}

/* Clearance between two fingers: closest approach of their capsule spines. */
static double finger_clearance(const Finger* a, const Finger* b) {
    V2 as[2][2] = { { a->P, a->J }, { a->J, a->T } };
    V2 bs[2][2] = { { b->P, b->J }, { b->J, b->T } };
    double m = 1e300;
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++) {
            double d = seg_seg(as[i][0], as[i][1], bs[j][0], bs[j][1]);
            if (d < m) m = d;
        }
    return m;
}

/* ---------- grasp description ---------- */

typedef struct {
    double db;        /* spread of the finger roots on the palm */
    double H;         /* roots sit at y = -H */
    double R;         /* object radius */
    double l1, l2;    /* bone lengths */
    double r;         /* finger capsule radius; capsules penetrate below 2r */
    double dc0;       /* nominal contact angle spacing */
    double margin;    /* extra clearance the solver aims for */
    double max_shift; /* allowed contact-angle deviation from nominal */
    double max_bend;  /* allowed bend difference between adjacent fingers */
    int    max_outer; /* augmented-Lagrangian outer iterations */
    double tol;       /* acceptance threshold on the constraint residual */
    double mu0, alpha, mu_max;
} Spec;

static Spec spec_default(void) {
    Spec s;
    s.db = 0.40; s.H = 2.0; s.R = 0.90;
    s.l1 = 0.85; s.l2 = 0.85; s.r = 0.08;
    s.dc0 = 0.20; s.margin = 0.006;
    s.max_shift = 0.60; s.max_bend = 0.30;
    s.max_outer = 150; s.tol = 1e-7;
    s.mu0 = 1.0; s.alpha = 4.0; s.mu_max = 1e10;
    return s;
}

static V2 spec_root(const Spec* s, int i) { return v2((i - 2) * s->db, -s->H); }
static double spec_nominal_phi(const Spec* s, int i) { return (i - 2) * s->dc0; }
static V2 spec_contact(const Spec* s, double phi) {
    return v2(s->R * sin(phi), -s->R * cos(phi));
}
static double spec_clear_target(const Spec* s) { return 2 * s->r + s->margin; }

static const int PAIRS[NPAIR][2] = {
    {0,1},{1,2},{2,3},{3,4},          /* adjacent */
    {0,2},{0,3},{0,4},{1,3},{1,4},{2,4}  /* skip-one */
};
static const int ADJ[NADJ][2] = { {0,1},{1,2},{2,3},{3,4} };

static void poses_for(const Spec* s, const double* phi, Finger* f) {
    for (int i = 0; i < NF; i++) {
        V2 P = spec_root(s, i), T = spec_contact(s, phi[i]);
        f[i] = solve_finger_side(P, T, s->l1, s->l2, outboard_side(P, T, s->l1, s->l2));
    }
}

/* Constraint vector: g[k] >= 0 means the k-th constraint holds. */
static void eval_constraints(const Spec* s, const double* phi, double* g) {
    Finger f[NF];
    poses_for(s, phi, f);
    int k = 0;
    for (int i = 0; i < NF; i++) {
        if (!f[i].ok) { g[k++] = -1e3; g[k++] = -1e3; continue; }
        g[k++] = f[i].dmax - f[i].dist - 1e-4;
        g[k++] = f[i].dist - f[i].dmin - 1e-4;
    }
    double target = spec_clear_target(s);
    for (int p = 0; p < NPAIR; p++) {
        int i = PAIRS[p][0], j = PAIRS[p][1];
        g[k++] = (f[i].ok && f[j].ok) ? finger_clearance(&f[i], &f[j]) - target : -1e3;
    }
    for (int p = 0; p < NADJ; p++) {
        int i = ADJ[p][0], j = ADJ[p][1];
        g[k++] = (f[i].ok && f[j].ok)
                 ? s->max_bend - fabs(f[i].bend - f[j].bend) : -1e3;
    }
    for (int i = 0; i < NF; i++)
        g[k++] = s->max_shift - fabs(phi[i] - spec_nominal_phi(s, i));
}

static double max_violation(const double* g, int n) {
    double worst = 0.0;
    for (int i = 0; i < n; i++) if (-g[i] > worst) worst = -g[i];
    return worst;
}

/* ---------- inner damped least squares (5 x 5) ---------- */

/* Solve A x = b by Gaussian elimination with partial pivoting. */
static int solve_linear(int n, double A[NV][NV], const double* b, double* x) {
    double M[NV][NV + 1];
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) M[i][j] = A[i][j];
        M[i][n] = b[i];
    }
    for (int c = 0; c < n; c++) {
        int piv = c;
        for (int r = c + 1; r < n; r++)
            if (fabs(M[r][c]) > fabs(M[piv][c])) piv = r;
        if (fabs(M[piv][c]) < 1e-18) return 0;
        if (piv != c)
            for (int j = c; j <= n; j++) { double t = M[c][j]; M[c][j] = M[piv][j]; M[piv][j] = t; }
        double d = M[c][c];
        for (int j = c; j <= n; j++) M[c][j] /= d;
        for (int r = 0; r < n; r++) {
            if (r == c) continue;
            double f = M[r][c];
            if (f == 0.0) continue;
            for (int j = c; j <= n; j++) M[r][j] -= f * M[c][j];
        }
    }
    for (int i = 0; i < n; i++) x[i] = M[i][n];
    return 1;
}

/* Augmented-Lagrangian residual for fixed multipliers and penalty weight. */
static void al_residual(const Spec* s, const double* phi,
                        const double* lambda, double mu, double* res) {
    for (int i = 0; i < NF; i++) res[i] = phi[i] - spec_nominal_phi(s, i);
    double g[NG];
    eval_constraints(s, phi, g);
    double scale = sqrt(2.0 / mu);
    for (int k = 0; k < NG; k++) {
        double v = lambda[k] - mu * g[k];
        res[NF + k] = v > 0 ? scale * v : 0.0;
    }
}

static double sumsq(const double* r, int n) {
    double s = 0.0;
    for (int i = 0; i < n; i++) s += r[i] * r[i];
    return s;
}

/* Levenberg-Marquardt on the augmented-Lagrangian residual, with a finite
 * difference Jacobian. Refines the decision variables in place. */
static void lm_inner(const Spec* s, double* x, const double* lambda, double mu,
                     int maxit) {
    const int nr = NF + NG;
    double r[NF + NG], rp[NF + NG], rm[NF + NG];
    al_residual(s, x, lambda, mu, r);
    double cost = 0.5 * sumsq(r, nr);
    double lam = 1e-3;
    const double h = 1e-6;

    for (int it = 0; it < maxit; it++) {
        double J[nr][NV], A[NV][NV], g[NV];
        double xp[NV], xm[NV];
        for (int k = 0; k < NV; k++) {
            memcpy(xp, x, sizeof(double) * NV);
            memcpy(xm, x, sizeof(double) * NV);
            xp[k] += h; xm[k] -= h;
            al_residual(s, xp, lambda, mu, rp);
            al_residual(s, xm, lambda, mu, rm);
            for (int i = 0; i < nr; i++) J[i][k] = (rp[i] - rm[i]) / (2 * h);
        }
        memset(A, 0, sizeof(A));
        memset(g, 0, sizeof(g));
        for (int i = 0; i < nr; i++)
            for (int a = 0; a < NV; a++) {
                g[a] += J[i][a] * r[i];
                for (int b = 0; b < NV; b++) A[a][b] += J[i][a] * J[i][b];
            }
        double dg[NV];
        for (int a = 0; a < NV; a++) dg[a] = A[a][a] > 1e-8 ? A[a][a] : 1e-8;

        int accepted = 0;
        double dx[NV] = { 0 };
        for (int attempt = 0; attempt < 15; attempt++) {
            double M[NV][NV], rhs[NV];
            for (int a = 0; a < NV; a++) {
                for (int b = 0; b < NV; b++) M[a][b] = A[a][b] + (a == b ? lam * dg[a] : 0.0);
                rhs[a] = -g[a];
            }
            double step[NV];
            if (!solve_linear(NV, M, rhs, step)) { lam *= 10; continue; }
            double nrm = 0.0;
            for (int a = 0; a < NV; a++) if (fabs(step[a]) > nrm) nrm = fabs(step[a]);
            if (nrm > 0.3) for (int a = 0; a < NV; a++) step[a] *= 0.3 / nrm;

            double xn[NV];
            for (int a = 0; a < NV; a++) xn[a] = x[a] + step[a];
            al_residual(s, xn, lambda, mu, rp);
            double cn = 0.5 * sumsq(rp, nr);
            if (cn < cost - 1e-16) {
                memcpy(x, xn, sizeof(double) * NV);
                memcpy(r, rp, sizeof(double) * nr);
                cost = cn;
                lam = lam * 0.3 > 1e-12 ? lam * 0.3 : 1e-12;
                memcpy(dx, step, sizeof(double) * NV);
                accepted = 1;
                break;
            }
            lam *= 4;
        }
        if (!accepted) break;
        double dmax = 0.0;
        for (int a = 0; a < NV; a++) if (fabs(dx[a]) > dmax) dmax = fabs(dx[a]);
        if (dmax < 1e-10) break;
    }
}

/* ---------- coupled solve ---------- */

typedef struct {
    int ok;             /* converged within tol: the pose is valid */
    double phi[NF];     /* solved contact angles */
    double residual;    /* worst constraint violation at the returned point */
    int iters;          /* outer iterations used */
    Finger finger[NF];
} SolveResult;

/* Per-finger independent solve: every contact on the ray from the object
 * centre through its own root. This is the naive baseline; the five fingers
 * are computed without reference to each other. */
static void independent_start(const Spec* s, double* phi) {
    for (int i = 0; i < NF; i++)
        phi[i] = atan2((i - 2) * s->db, s->H);
}

static void solve_coupled(const Spec* s, const double* x0, SolveResult* out) {
    double x[NV], lambda[NG], g[NG];
    memcpy(x, x0, sizeof(double) * NV);
    memset(lambda, 0, sizeof(lambda));
    double mu = s->mu0;
    int outer = 0;

    for (outer = 1; outer <= s->max_outer; outer++) {
        lm_inner(s, x, lambda, mu, 120);
        eval_constraints(s, x, g);
        double resid = max_violation(g, NG);
        if (resid <= s->tol) break;
        for (int k = 0; k < NG; k++) {
            double v = lambda[k] - mu * g[k];
            lambda[k] = v > 0 ? v : 0.0;
        }
        mu = mu * s->alpha < s->mu_max ? mu * s->alpha : s->mu_max;
    }

    eval_constraints(s, x, g);
    out->residual = max_violation(g, NG);
    out->iters = outer > s->max_outer ? s->max_outer : outer;
    out->ok = out->residual <= s->tol;
    memcpy(out->phi, x, sizeof(double) * NV);
    poses_for(s, x, out->finger);
}

/* Number of finger pairs whose capsules overlap (spine distance < 2r). */
static int penetration_count(const Spec* s, const double* phi) {
    Finger f[NF];
    poses_for(s, phi, f);
    int count = 0;
    for (int p = 0; p < NPAIR; p++) {
        int i = PAIRS[p][0], j = PAIRS[p][1];
        if (!f[i].ok || !f[j].ok) { count++; continue; }
        if (finger_clearance(&f[i], &f[j]) < 2 * s->r - 1e-9) count++;
    }
    return count;
}

/* ---------- grasp-matrix force allocation ---------- */

typedef struct { V2 p; V2 n; } Contact;   /* contact point and inward normal */

typedef struct {
    double f[NF];       /* contact force per finger (compression positive) */
    double lambda;      /* common ratio f_i / w_i */
    double resid[3];    /* contact wrench delivered minus required */
    double resid_norm;
    double sum;         /* sum of contact forces */
    int ok;             /* wrench met within tol, and lambda >= 0 */
} ForceResult;

/* Required contact wrench w_ext (fx, fy, torque) is distributed over the
 * contacts in proportion to weights: f_i = lambda * w_i. The single scale
 * lambda is the least-squares fit to the wrench balance G f = w_ext, so the
 * allocation is proportional by construction; if those proportional forces
 * cannot realise the wrench, the leftover is reported, not hidden. */
static ForceResult allocate_force(const Contact* c, int n, const double* w,
                                  const double w_ext[3], double tol) {
    ForceResult out;
    memset(&out, 0, sizeof(out));
    if (n <= 0 || n > NF) return out;

    /* Grasp matrix columns: [n_x, n_y, p x n]. */
    double J[3][NF];
    for (int k = 0; k < n; k++) {
        J[0][k] = c[k].n.x;
        J[1][k] = c[k].n.y;
        J[2][k] = c[k].p.x * c[k].n.y - c[k].p.y * c[k].n.x;
    }
    double Jw[3];
    for (int a = 0; a < 3; a++) {
        Jw[a] = 0.0;
        for (int k = 0; k < n; k++) Jw[a] += J[a][k] * w[k];
    }
    double den = Jw[0] * Jw[0] + Jw[1] * Jw[1] + Jw[2] * Jw[2];
    double num = Jw[0] * w_ext[0] + Jw[1] * w_ext[1] + Jw[2] * w_ext[2];
    out.lambda = den > 0 ? num / den : 0.0;

    for (int k = 0; k < n; k++) {
        out.f[k] = out.lambda * w[k];
        out.sum += out.f[k];
    }
    for (int a = 0; a < 3; a++) {
        double v = 0.0;
        for (int k = 0; k < n; k++) v += J[a][k] * out.f[k];
        out.resid[a] = v - w_ext[a];
    }
    out.resid_norm = sqrt(out.resid[0] * out.resid[0] +
                          out.resid[1] * out.resid[1] +
                          out.resid[2] * out.resid[2]);
    out.ok = (out.resid_norm <= tol) && (out.lambda >= 0.0);
    return out;
}

/* ---------- self tests ---------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    Spec s = spec_default();

    /* A single finger is an analytic triangle: the joint lies on both bone
     * circles, and the two sides are mirror images of each other. */
    Finger f = solve_finger_side(v2(0, 0), v2(1.0, 0), 1.0, 1.0, +1);
    check("finger_reachable", f.ok == 1);
    check("joint_on_bone1", fabs(vlen(vsub(f.J, f.P)) - 1.0) < 1e-9);
    check("joint_on_bone2", fabs(vlen(vsub(f.J, f.T)) - 1.0) < 1e-9);
    check("finger_side_up", f.J.y > 0.0);
    Finger fm = solve_finger_side(v2(0, 0), v2(1.0, 0), 1.0, 1.0, -1);
    check("finger_side_mirror", fm.J.y < 0.0 && fabs(fm.J.x - f.J.x) < 1e-9);
    Finger far = solve_finger_side(v2(0, 0), v2(3.0, 0), 1.0, 1.0, +1);
    check("finger_unreachable", far.ok == 0);

    /* The elbow side rule is deterministic and mirror-symmetric. */
    check("side_rule_left", outboard_side(v2(-0.8, -2.0), spec_contact(&s, -0.4), s.l1, s.l2) == +1);
    check("side_rule_right", outboard_side(v2( 0.8, -2.0), spec_contact(&s,  0.4), s.l1, s.l2) == -1);

    /* Independent per-finger solve: every contact on its own root ray. */
    double indep[NV];
    independent_start(&s, indep);
    int indep_pen = penetration_count(&s, indep);
    check("independent_penetrates", indep_pen > 0);

    /* Coupled solve from that same starting point. */
    SolveResult r;
    solve_coupled(&s, indep, &r);
    check("coupled_converges", r.ok == 1);
    check("coupled_residual_small", r.residual <= s.tol);
    check("coupled_zero_penetration", penetration_count(&s, r.phi) == 0);
    check("coupled_beats_independent", penetration_count(&s, r.phi) < indep_pen);

    /* The solved pose honours its own shared constraints. */
    int reach_ok = 1, bend_ok = 1, shift_ok = 1;
    for (int i = 0; i < NF; i++) {
        if (!r.finger[i].ok) reach_ok = 0;
        if (fabs(r.phi[i] - spec_nominal_phi(&s, i)) > s.max_shift + 1e-9) shift_ok = 0;
    }
    for (int p = 0; p < NADJ; p++)
        if (fabs(r.finger[ADJ[p][0]].bend - r.finger[ADJ[p][1]].bend) > s.max_bend + 1e-9)
            bend_ok = 0;
    check("coupled_reachable", reach_ok);
    check("coupled_bend_shared", bend_ok);
    check("coupled_shift_bounded", shift_ok);

    /* Infeasible specifications must report a residual, not a pose. */
    Spec tight_bend = s; tight_bend.max_bend = 0.05;
    SolveResult rb;
    solve_coupled(&tight_bend, indep, &rb);
    check("infeasible_bend_not_ok", rb.ok == 0);
    check("infeasible_bend_reports_residual", rb.residual > s.tol);

    Spec tight_shift = s; tight_shift.max_shift = 0.02;
    SolveResult rs;
    solve_coupled(&tight_shift, indep, &rs);
    check("infeasible_shift_not_ok", rs.ok == 0);
    check("infeasible_shift_reports_residual", rs.residual > s.tol);

    /* Force allocation on a simplified flat set of contacts (all normals the
     * same). Here a vertical load is exactly realisable, so conservation and
     * proportionality can be checked cleanly. */
    Contact c[NF];
    double pos[NF] = { -0.6, -0.3, 0.0, 0.3, 0.6 };
    for (int i = 0; i < NF; i++) {
        c[i].p = v2(pos[i], 0.0);
        c[i].n = v2(0.0, 1.0);     /* fingers push up into the object */
    }
    double w_eq[NF] = { 1, 1, 1, 1, 1 };
    double load[3] = { 0.0, 10.0, 0.0 };   /* contacts must supply +10 up */
    ForceResult f1 = allocate_force(c, NF, w_eq, load, 1e-9);
    check("force_balanced", f1.ok == 1);
    check("force_sum_conserved", fabs(f1.sum - 10.0) < 1e-9);
    int prop_eq = 1;
    for (int i = 0; i < NF; i++)
        if (fabs(f1.f[i] - 2.0) > 1e-9) prop_eq = 0;
    check("force_equal_split", prop_eq);

    /* Redistributing weights keeps the delivered total fixed. */
    double w_sym[NF] = { 1, 2, 3, 2, 1 };
    ForceResult f2 = allocate_force(c, NF, w_sym, load, 1e-9);
    int prop_ok = 1;
    for (int i = 1; i < NF; i++)
        if (fabs(f2.f[i] / w_sym[i] - f2.f[0] / w_sym[0]) > 1e-9) prop_ok = 0;
    check("force_proportional", prop_ok);
    check("force_redistribute_conserved", fabs(f2.sum - 10.0) < 1e-9);
    check("force_redistribute_balanced", f2.ok == 1);
    check("force_weighted_shares", fabs(f2.f[2] / f2.f[0] - 3.0) < 1e-9);

    /* Weights that break the torque balance cannot be realised by
     * proportional forces: the residual is reported and the result is not ok. */
    double w_bad[NF] = { 1, 1, 1, 1, 3 };
    ForceResult f3 = allocate_force(c, NF, w_bad, load, 1e-9);
    check("force_unbalanced_not_ok", f3.ok == 0);
    check("force_unbalanced_reports_residual", f3.resid_norm > 1e-6);

    /* A wrench no proportional allocation can supply is likewise flagged. */
    double lateral[3] = { 5.0, 0.0, 0.0 };   /* all normals vertical */
    ForceResult f4 = allocate_force(c, NF, w_eq, lateral, 1e-9);
    check("force_lateral_not_ok", f4.ok == 0);
    check("force_lateral_reports_residual", f4.resid_norm > 1e-6);

    /* Contacts taken from the coupled pose itself (circular object). The
     * normals are radial and not parallel, so a single proportional scale
     * cannot satisfy both force and wrench exactly: the solver must report the
     * residual rather than claim success. */
    Contact rc[NF];
    double dir_load[3] = { 0.0, 10.0, 0.0 };
    for (int i = 0; i < NF; i++) {
        V2 T = r.finger[i].T;
        double l = vlen(T);
        rc[i].p = T;
        rc[i].n = l > 0 ? v2(-T.x / l, -T.y / l) : v2(0, 1);   /* inward */
    }
    ForceResult f5 = allocate_force(rc, NF, w_eq, dir_load, 1e-9);
    double torque = 0.0;
    for (int i = 0; i < NF; i++) torque += f5.f[i] * (rc[i].p.x * rc[i].n.y - rc[i].p.y * rc[i].n.x);
    /* p x n vanishes for radial normals; the residual that matters is the
     * force mismatch under one shared proportional scale. */
    check("force_real_contacts_flagged", f5.ok == 0 && f5.resid_norm > 1e-3);

    printf("independent_penetrations=%d\n", indep_pen);
    printf("coupled_penetrations=%d\n", penetration_count(&s, r.phi));
    printf("coupled_residual=%.3e\n", r.residual);
    printf("coupled_iters=%d\n", r.iters);
    printf("infeasible_bend_residual=%.3e\n", rb.residual);
    printf("infeasible_shift_residual=%.3e\n", rs.residual);
    printf("force_sum=%.6f\n", f2.sum);
    printf("force_ratio=%.6f\n", f2.f[2] / f2.f[0]);
    printf("force_unbalanced_residual=%.6f\n", f3.resid_norm);
    printf("force_real_contact_torque=%.6f\n", torque);
    printf("force_real_contact_residual=%.6f\n", f5.resid_norm);
    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
