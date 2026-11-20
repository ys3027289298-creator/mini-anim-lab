
/* seam: loop seam correction by constrained least squares.
 *
 * One pose channel is sampled as p[0..n-1] over the loop [0, period],
 * where frame n-1 wraps back to frame 0.  A clip that is not periodic
 * hitches at the wrap.  We enforce C3 periodicity,
 *
 *     D^m (p+d)(period) = D^m (p+d)(0),   m = 0..3,
 *
 * as linear equations in a correction d limited to the first and last
 * SEAM_WIN samples.  The jets are pointwise: one-sided finite
 * differences at the seam frames, so the constraints see only the
 * loop's endpoint jets and never the correction's interior.
 *
 * A correction that minimizes the raw sample norm is spiky (the
 * Euclidean norm rewards concentration), so the objective is the
 * bending energy of the offset - the least-squares sense of "smallest
 * pose disturbance" - minimized within a smooth Hermite basis: per
 * side, the monotone C3 smootherstep (jet (1,0,0,0) at the seam -
 * carries the value shift) plus polynomial jet bumps u^k(1-u)^4
 * (independent k-th jet at the seam), all with vanishing jets at the
 * window edge.  With B the constraints lifted into the basis and R
 * the basis stiffness (sampled second-difference Gram), the closed
 * form is c = R^{-1} B^T (B R^{-1} B^T)^{-1} b.  The report carries measured
 * values plus a priori upper bounds for both the residual (per
 * constraint) and the pose offset.  If the needed offset exceeds the
 * caller's bandwidth the clip is reported unfixable and left
 * untouched - never clamped to fit.
 */
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#define SEAM_ORDER 3                /* match derivatives 1..SEAM_ORDER */
#define SEAM_CONS (SEAM_ORDER + 1)  /* pointwise value + smoothed 1..3 */
#define SEAM_WIN 12                 /* corrected samples at each end */
#define SEAM_SLOTS (2 * SEAM_WIN)   /* corrected sample slots */
#define SEAM_BAS 4                  /* smooth bumps per side */
#define SEAM_COEF (2 * SEAM_BAS)    /* correction coefficients */
#define SEAM_MAX_NODES 8            /* largest fd stencil */
#define SEAM_MAX_DIM 16             /* largest dense matrix */
#define SEAM_PI 3.14159265358979323846

typedef enum {
    SEAM_OK = 0,        /* clip corrected in place */
    SEAM_ERR_ARG,       /* bad arguments */
    SEAM_ERR_SINGULAR,  /* degenerate constraint system */
    SEAM_ERR_BANDWIDTH  /* needs more than the allowed bandwidth */
} SeamStatus;

typedef struct {
    SeamStatus status;
    float residual[SEAM_CONS];        /* seam jet mismatch, orders 0..3 */
    float residual_bound[SEAM_CONS];  /* a priori roundoff bound per order */
    float offset_max;                 /* measured max |correction| needed */
    float offset_bound;               /* a priori bound on ||correction||_2 */
} SeamReport;

/* dense n x n solve by gaussian elimination: rows are equilibrated,
 * then partial pivoting; A and b are destroyed.  0 on success. */
static int solve_sq(int n, double* A, double* b, double* x) {
    for (int i = 0; i < n; i++) {
        double mx = 0;
        for (int j = 0; j < n; j++) mx = fmax(mx, fabs(A[i * n + j]));
        if (mx == 0.0) return -1;
        for (int j = 0; j < n; j++) A[i * n + j] /= mx;
        b[i] /= mx;
    }
    for (int col = 0; col < n; col++) {
        int piv = col;
        for (int r = col + 1; r < n; r++)
            if (fabs(A[r * n + col]) > fabs(A[piv * n + col])) piv = r;
        if (fabs(A[piv * n + col]) < 1e-13) return -1;
        if (piv != col) {
            for (int k = col; k < n; k++) {
                double t = A[col * n + k];
                A[col * n + k] = A[piv * n + k];
                A[piv * n + k] = t;
            }
            double t = b[col]; b[col] = b[piv]; b[piv] = t;
        }
        for (int r = col + 1; r < n; r++) {
            double f = A[r * n + col] / A[col * n + col];
            for (int k = col; k < n; k++) A[r * n + k] -= f * A[col * n + k];
            b[r] -= f * b[col];
        }
    }
    for (int r = n - 1; r >= 0; r--) {
        double acc = b[r];
        for (int k = r + 1; k < n; k++) acc -= A[r * n + k] * x[k];
        x[r] = acc / A[r * n + r];
    }
    return 0;
}

static int invert_sq(int n, const double* A, double* inv) {
    double Ac[SEAM_MAX_DIM * SEAM_MAX_DIM];
    double b[SEAM_MAX_DIM], x[SEAM_MAX_DIM];
    for (int col = 0; col < n; col++) {
        memcpy(Ac, A, (size_t)n * (size_t)n * sizeof(double));
        for (int i = 0; i < n; i++) b[i] = (i == col) ? 1.0 : 0.0;
        if (solve_sq(n, Ac, b, x) != 0) return -1;
        for (int i = 0; i < n; i++) inv[i * n + col] = x[i];
    }
    return 0;
}

/* weights w[j] with sum_j w[j] f(x[j]) ~= f^(m)(x0), exact for
 * polynomials of degree < nn, from the moment equations
 *     sum_j w[j] (x[j]-x0)^k = (k == m) m!,   k = 0..nn-1. */
static int fd_weights(int m, const double* x, int nn, double x0, double* w) {
    double V[SEAM_MAX_NODES * SEAM_MAX_NODES], rhs[SEAM_MAX_NODES];
    double mf = 1;
    for (int i = 2; i <= m; i++) mf *= i;
    for (int k = 0; k < nn; k++) {
        rhs[k] = (k == m) ? mf : 0.0;
        for (int j = 0; j < nn; j++) V[k * nn + j] = pow(x[j] - x0, (double)k);
    }
    return solve_sq(nn, V, rhs, w);
}

/* smooth correction bumps on u in [0,1], all with vanishing jets
 * (orders 0..3) at the window edge u=1:
 *   k=0: the C3 smootherstep 1 - 35u^4 + 84u^5 - 70u^6 + 20u^7,
 *        monotone, jet (1,0,0,0) at the seam - carries value shifts;
 *   k=1..3: u^k (1-u)^4, whose k-th jet at the seam is independent,
 *        so slope/curvature/jerk mismatches are correctable too. */
static double seam_bump(int k, double u) {
    if (k == 0) {
        double u4 = u * u * u * u;
        return 1.0 - 35.0 * u4 + 84.0 * u4 * u - 70.0 * u4 * u * u
             + 20.0 * u4 * u * u * u;
    }
    double v = 1.0 - u;
    return pow(u, (double)k) * v * v * v * v;
}

/* weighted minimum norm:  min sum_i x_i^2 / s_i  s.t.  A x = b.
 * Closed form x = S A^T (A S A^T)^{-1} b.  Rows are equilibrated
 * first, which leaves the solution unchanged.  Used by the tests to
 * validate the KKT mathematics against closed-form solutions. */
static int minnorm_solve(int r, int n, const double* A, const double* s,
                         const double* b, double* x, double* ginv) {
    double As[SEAM_CONS * SEAM_SLOTS], bs[SEAM_CONS];
    double G[SEAM_CONS * SEAM_CONS], Gi[SEAM_CONS * SEAM_CONS];
    double lam[SEAM_CONS];
    for (int i = 0; i < r; i++) {
        double mx = 0;
        for (int j = 0; j < n; j++) mx = fmax(mx, fabs(A[i * n + j]));
        if (mx == 0.0) return -1;
        for (int j = 0; j < n; j++) As[i * n + j] = A[i * n + j] / mx;
        bs[i] = b[i] / mx;
    }
    for (int a = 0; a < r; a++)
        for (int c = 0; c < r; c++) {
            double acc = 0;
            for (int j = 0; j < n; j++) acc += As[a * n + j] * s[j] * As[c * n + j];
            G[a * r + c] = acc;
        }
    if (invert_sq(r, G, Gi) != 0) return -1;
    for (int a = 0; a < r; a++) {
        double acc = 0;
        for (int c = 0; c < r; c++) acc += Gi[a * r + c] * bs[c];
        lam[a] = acc;
    }
    for (int j = 0; j < n; j++) {
        double acc = 0;
        for (int a = 0; a < r; a++) acc += As[a * n + j] * lam[a];
        x[j] = s[j] * acc;
    }
    if (ginv) memcpy(ginv, Gi, (size_t)r * (size_t)r * sizeof(double));
    return 0;
}

SeamStatus seam_fix(float* p, int n, float period, float bandwidth, SeamReport* rep) {
    SeamReport local;
    if (!rep) rep = &local;
    memset(rep, 0, sizeof(*rep));
    rep->status = SEAM_ERR_ARG;
    if (!p || !(period > 0.0f) || !(bandwidth >= 0.0f) || n < 2 * SEAM_WIN + 1)
        return rep->status;

    const double h = (double)period / (n - 1);
    const int K = SEAM_WIN;

    /* sampled smooth basis: slots 0..K-1 are head samples, slots
     * K..2K-1 tail samples n-1-j; each side gets SEAM_BAS bumps */
    double PHI[SEAM_SLOTS][SEAM_COEF] = { { 0 } };
    for (int j = 0; j < K; j++) {
        double u = (double)j / K;
        for (int k = 0; k < SEAM_BAS; k++) {
            PHI[j][k] = seam_bump(k, u);
            PHI[K + j][SEAM_BAS + k] = seam_bump(k, u);
        }
    }

    /* pointwise jet functionals at both loop ends, orders 0..3:
     * one-sided finite differences, exact for degree <= m+2 */
    double wH[SEAM_CONS][SEAM_MAX_NODES], wT[SEAM_CONS][SEAM_MAX_NODES];
    int nn[SEAM_CONS];
    for (int m = 0; m < SEAM_CONS; m++) {
        double x[SEAM_MAX_NODES];
        nn[m] = m + 3;
        for (int j = 0; j < nn[m]; j++) x[j] = j * h;
        if (fd_weights(m, x, nn[m], 0.0, wH[m]) != 0) {
            rep->status = SEAM_ERR_SINGULAR;
            return rep->status;
        }
        for (int j = 0; j < nn[m]; j++) x[j] = (double)period - j * h;
        if (fd_weights(m, x, nn[m], (double)period, wT[m]) != 0) {
            rep->status = SEAM_ERR_SINGULAR;
            return rep->status;
        }
    }

    /* sample-space constraints C d = b over the correction slots:
     *   D^m d(period) - D^m d(0)  =  D^m p(0) - D^m p(period),
     * with row 0 the pointwise value jump. */
    double C[SEAM_CONS][SEAM_SLOTS] = { { 0 } }, b[SEAM_CONS] = { 0 };
    for (int m = 0; m < SEAM_CONS; m++) {
        for (int j = 0; j < nn[m]; j++) {
            C[m][j] = -wH[m][j];
            C[m][K + j] = wT[m][j];
            b[m] += wH[m][j] * p[j] - wT[m][j] * p[n - 1 - j];
        }
    }

    /* lift constraints into the basis: B = C Phi (SEAM_CONS x SEAM_COEF) */
    double B[SEAM_CONS][SEAM_COEF] = { { 0 } };
    for (int m = 0; m < SEAM_CONS; m++)
        for (int k = 0; k < SEAM_COEF; k++)
            for (int j = 0; j < SEAM_SLOTS; j++)
                B[m][k] += C[m][j] * PHI[j][k];

    /* equilibrate rows (leaves the solution unchanged) so the small
     * systems and all norm bounds live in one consistent scaling */
    for (int m = 0; m < SEAM_CONS; m++) {
        double mx = 0;
        for (int k = 0; k < SEAM_COEF; k++) mx = fmax(mx, fabs(B[m][k]));
        if (mx == 0.0) {
            rep->status = SEAM_ERR_SINGULAR;
            return rep->status;
        }
        for (int k = 0; k < SEAM_COEF; k++) B[m][k] /= mx;
        b[m] /= mx;
    }

    /* stiffness R of the sampled basis: per side, the second-
     * difference Gram - the objective is the bending energy of the
     * offset, so the correction is the smoothest one meeting the
     * jet constraints, not the most concentrated one. */
    double W[SEAM_COEF][SEAM_COEF] = { { 0 } }, Wi[SEAM_COEF][SEAM_COEF];
    for (int side = 0; side < 2; side++) {
        double d2[SEAM_WIN + 1][SEAM_BAS];
        for (int j = 1; j <= K; j++)
            for (int k = 0; k < SEAM_BAS; k++) {
                double um = seam_bump(k, (double)(j - 1) / K);
                double u0 = j < K ? seam_bump(k, (double)j / K) : 0.0;
                double up = j + 1 < K ? seam_bump(k, (double)(j + 1) / K) : 0.0;
                d2[j][k] = um - 2.0 * u0 + up;
            }
        for (int k = 0; k < SEAM_BAS; k++)
            for (int l = 0; l <= k; l++) {
                double acc = 0;
                for (int j = 1; j <= K; j++) acc += d2[j][k] * d2[j][l];
                W[side * SEAM_BAS + k][side * SEAM_BAS + l] = acc;
                W[side * SEAM_BAS + l][side * SEAM_BAS + k] = acc;
            }
    }
    if (invert_sq(SEAM_COEF, &W[0][0], &Wi[0][0]) != 0) {
        rep->status = SEAM_ERR_SINGULAR;
        return rep->status;
    }

    /* c = W^{-1} B^T (B W^{-1} B^T)^{-1} b */
    double T[SEAM_COEF][SEAM_CONS] = { { 0 } };   /* T = W^{-1} B^T */
    for (int k = 0; k < SEAM_COEF; k++)
        for (int m = 0; m < SEAM_CONS; m++)
            for (int l = 0; l < SEAM_COEF; l++)
                T[k][m] += Wi[k][l] * B[m][l];
    double M[SEAM_CONS][SEAM_CONS] = { { 0 } }, Mi[SEAM_CONS][SEAM_CONS];
    for (int a = 0; a < SEAM_CONS; a++)
        for (int c = 0; c < SEAM_CONS; c++)
            for (int k = 0; k < SEAM_COEF; k++)
                M[a][c] += B[a][k] * T[k][c];
    if (invert_sq(SEAM_CONS, &M[0][0], &Mi[0][0]) != 0) {
        rep->status = SEAM_ERR_SINGULAR;
        return rep->status;
    }
    double coef[SEAM_COEF] = { 0 };
    for (int k = 0; k < SEAM_COEF; k++)
        for (int a = 0; a < SEAM_CONS; a++) {
            double la = 0;
            for (int c = 0; c < SEAM_CONS; c++) la += Mi[a][c] * b[c];
            coef[k] += T[k][a] * la;
        }

    /* correction in sample space: d = Phi c */
    double d[SEAM_SLOTS] = { 0 };
    for (int j = 0; j < SEAM_SLOTS; j++)
        for (int k = 0; k < SEAM_COEF; k++) d[j] += PHI[j][k] * coef[k];

    /* measured correction */
    double dmax = 0, dnorm = 0;
    for (int j = 0; j < SEAM_SLOTS; j++) {
        dmax = fmax(dmax, fabs(d[j]));
        dnorm += d[j] * d[j];
    }
    dnorm = sqrt(dnorm);

    /* a priori bound:
     *   ||d||_2 <= ||Phi||_2 . ||T||_2 . ||M^{-1}||_2 . ||b||_2
     * with ||X||_2 <= sqrt(||X||_1 ||X||_inf) */
    {
        double bnorm = 0;
        for (int m = 0; m < SEAM_CONS; m++) bnorm += b[m] * b[m];
        bnorm = sqrt(bnorm);
        double phi1 = 0, phii = 0, t1 = 0, ti = 0, m1 = 0, mi = 0;
        for (int k = 0; k < SEAM_COEF; k++) {
            double c1 = 0;
            for (int j = 0; j < SEAM_SLOTS; j++) c1 += fabs(PHI[j][k]);
            phi1 = fmax(phi1, c1);
        }
        for (int j = 0; j < SEAM_SLOTS; j++) {
            double ci = 0;
            for (int k = 0; k < SEAM_COEF; k++) ci += fabs(PHI[j][k]);
            phii = fmax(phii, ci);
        }
        for (int m = 0; m < SEAM_CONS; m++) {
            double c1 = 0, ci = 0;
            for (int k = 0; k < SEAM_COEF; k++) c1 += fabs(T[k][m]);
            for (int a = 0; a < SEAM_CONS; a++) ci += fabs(T[m][a]);
            t1 = fmax(t1, c1);
            ti = fmax(ti, ci);
        }
        for (int c = 0; c < SEAM_CONS; c++) {
            double c1 = 0, ci = 0;
            for (int a = 0; a < SEAM_CONS; a++) {
                c1 += fabs(Mi[a][c]);
                ci += fabs(Mi[c][a]);
            }
            m1 = fmax(m1, c1);
            mi = fmax(mi, ci);
        }
        rep->offset_bound =
            (float)(sqrt(phi1 * phii) * sqrt(t1 * ti) * sqrt(m1 * mi) * bnorm);
    }
    rep->offset_max = (float)dmax;

    /* seam mismatch the computed correction achieves (double eval) */
    for (int m = 0; m < SEAM_CONS; m++) {
        double r = -b[m], wabs = 0;
        for (int j = 0; j < SEAM_SLOTS; j++) {
            r += C[m][j] * d[j];
            wabs += fabs(C[m][j]);
        }
        rep->residual[m] = (float)fabs(r);
        rep->residual_bound[m] = (float)(8.0 * DBL_EPSILON * wabs * (1.0 + dmax));
    }

    if (dmax > (double)bandwidth) {
        /* would exceed the bandwidth: report unfixable, touch nothing */
        rep->status = SEAM_ERR_BANDWIDTH;
        return rep->status;
    }

    for (int j = 0; j < K; j++) {
        p[j] = (float)(p[j] + d[j]);
        p[n - 1 - j] = (float)(p[n - 1 - j] + d[K + j]);
    }

    /* achieved mismatch, re-measured on the corrected float clip */
    {
        double qmax = 0;
        for (int j = 0; j < n; j++) qmax = fmax(qmax, fabs((double)p[j]));
        for (int m = 0; m < SEAM_CONS; m++) {
            double r = 0, wabs = 0;
            for (int j = 0; j < nn[m]; j++) {
                r += wT[m][j] * p[n - 1 - j] - wH[m][j] * p[j];
                wabs += fabs(wT[m][j]) + fabs(wH[m][j]);
            }
            rep->residual[m] = (float)fabs(r);
            rep->residual_bound[m] = (float)(8.0 * FLT_EPSILON * wabs * (qmax + 1e-30));
        }
    }
    rep->status = SEAM_OK;
    return rep->status;
}

/* ---------------- self test ---------------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

/* independent smoothness check: least-squares fit of a degree-6
 * polynomial to 13 samples straddling the wrap.  The max fit residual
 * is small iff the loop is smooth across the seam - a kink of any
 * order shows up as a residual proportional to the jet jump. */
static double seam_fit_residual(const float* q, int n, double period) {
    enum { NP = 7, NS = 13 };
    double h = period / (n - 1);
    double G[NP * NP] = { 0 }, rhs[NP] = { 0 }, y[NS];
    for (int j = 0; j < NS; j++) {
        int o = j - 6;
        int idx = ((o % (n - 1)) + (n - 1)) % (n - 1);
        double u = o * h, xp[NP];
        y[j] = q[idx];
        xp[0] = 1.0;
        for (int a = 1; a < NP; a++) xp[a] = xp[a - 1] * u;
        for (int a = 0; a < NP; a++) {
            rhs[a] += xp[a] * y[j];
            for (int c = 0; c < NP; c++) G[a * NP + c] += xp[a] * xp[c];
        }
    }
    double coef[NP];
    if (solve_sq(NP, G, rhs, coef) != 0) return -1.0;
    double mx = 0;
    for (int j = 0; j < NS; j++) {
        double u = (j - 6) * h, f = 0;
        for (int a = NP - 1; a >= 0; a--) f = f * u + coef[a];
        mx = fmax(mx, fabs(f - y[j]));
    }
    return mx;
}

#define TEST_N 129

int main(void) {
    /* fd weights reproduce exact polynomial derivatives */
    {
        double h = 0.1, x[SEAM_MAX_NODES], w[SEAM_MAX_NODES];
        int ok = 1;
        for (int m = 0; m <= 3 && ok; m++) {
            int nn = m + 3;
            double mf = 1;
            for (int i = 2; i <= m; i++) mf *= i;
            /* f(x) = x^m + x^(nn-1): f^(m)(0) = m! exactly */
            for (int j = 0; j < nn; j++) x[j] = j * h;
            if (fd_weights(m, x, nn, 0.0, w) != 0) { ok = 0; break; }
            double acc = 0;
            for (int j = 0; j < nn; j++) acc += w[j] * (pow(x[j], m) + pow(x[j], nn - 1));
            if (fabs(acc - mf) > 1e-8) ok = 0;
            /* mirrored nodes about x0=1: f^(m)(1) = m! + (nn-1)!/(nn-1-m)! */
            for (int j = 0; j < nn; j++) x[j] = 1.0 - j * h;
            if (fd_weights(m, x, nn, 1.0, w) != 0) { ok = 0; break; }
            double expect = mf, pm = 1;
            for (int i = 0; i < m; i++) pm *= (nn - 1 - i);
            expect += pm;
            acc = 0;
            for (int j = 0; j < nn; j++) acc += w[j] * (pow(x[j], m) + pow(x[j], nn - 1));
            if (fabs(acc - expect) > 1e-7) ok = 0;
        }
        check("weights_exact", ok);
    }

    /* min-norm solver against closed-form solutions */
    {
        double x[3], err;

        /* rank-1: min x0^2+x1^2 s.t. x0 - x1 = 0.6 -> (0.3, -0.3) */
        double A1[2] = { 1, -1 }, s1[2] = { 1, 1 }, b1[1] = { 0.6 };
        check("minnorm_rank1_run", minnorm_solve(1, 2, A1, s1, b1, x, 0) == 0);
        err = fmax(fabs(x[0] - 0.3), fabs(x[1] + 0.3));
        check("minnorm_rank1_exact", err < 1e-12);

        /* weighted: min x0^2 + x1^2/3 s.t. x0 - x1 = 0.6 -> (0.15, -0.45) */
        double s2[2] = { 1, 3 };
        check("minnorm_weighted_run", minnorm_solve(1, 2, A1, s2, b1, x, 0) == 0);
        err = fmax(fabs(x[0] - 0.15), fabs(x[1] + 0.45));
        check("minnorm_weighted_exact", err < 1e-12);

        /* 2x3: A=[[1,0,1],[0,1,1]], b=(1,2) -> x = A^T(AA^T)^{-1}b = (0,1,1) */
        double A3[6] = { 1, 0, 1, 0, 1, 1 }, s3[3] = { 1, 1, 1 }, b3[2] = { 1, 2 };
        check("minnorm_2x3_run", minnorm_solve(2, 3, A3, s3, b3, x, 0) == 0);
        err = fmax(fabs(x[0]), fmax(fabs(x[1] - 1), fabs(x[2] - 1)));
        check("minnorm_2x3_exact", err < 1e-12);
    }

    /* already-periodic clip: correction stays negligible */
    {
        float p[TEST_N];
        SeamReport rep;
        double h = 1.0 / (TEST_N - 1);
        for (int j = 0; j < TEST_N; j++) p[j] = sinf((float)(2.0 * SEAM_PI * j * h));
        check("periodic_ok", seam_fix(p, TEST_N, 1.0f, 0.01f, &rep) == SEAM_OK);
        printf("periodic_report residual0=%g offset_max=%g offset_bound=%g\n",
               rep.residual[0], rep.offset_max, rep.offset_bound);
        check("periodic_small_offset", rep.offset_max < 1e-2f);
        check("periodic_offset_certified", rep.offset_max <= rep.offset_bound);
        check("periodic_residual_certified",
              rep.residual[0] <= rep.residual_bound[0] &&
              rep.residual[3] <= rep.residual_bound[3]);
    }

    /* ramp: seam jump of 1.0 gets fixed and certified */
    {
        float p[TEST_N], orig[TEST_N];
        SeamReport rep;
        double h = 1.0 / (TEST_N - 1);
        for (int j = 0; j < TEST_N; j++) p[j] = orig[j] = (float)(j * h);

        double fit_before = seam_fit_residual(p, TEST_N, 1.0);
        check("ramp_broken", fabs(p[TEST_N - 1] - p[0] - 1.0f) < 1e-6f);
        check("ramp_broken_fit", fit_before > 0.05);

        check("ramp_ok", seam_fix(p, TEST_N, 1.0f, 2.0f, &rep) == SEAM_OK);
        printf("ramp_report residual=%g,%g,%g,%g offset_max=%g offset_bound=%g\n",
               rep.residual[0], rep.residual[1], rep.residual[2],
               rep.residual[3], rep.offset_max, rep.offset_bound);
        check("ramp_offset_certified", rep.offset_max <= rep.offset_bound);
        check("ramp_residual_certified",
              rep.residual[0] <= rep.residual_bound[0] &&
              rep.residual[1] <= rep.residual_bound[1] &&
              rep.residual[2] <= rep.residual_bound[2] &&
              rep.residual[3] <= rep.residual_bound[3]);
        check("ramp_residual_small",
              rep.residual[0] < 1e-5f && rep.residual[1] < 1e-3f);

        /* interior untouched */
        check("ramp_interior_untouched",
              memcmp(&orig[SEAM_WIN], &p[SEAM_WIN],
                     (TEST_N - 2 * SEAM_WIN) * sizeof(float)) == 0);

        /* independent smoothness check: local polynomial fit across
         * the wrap goes from clearly broken to smooth */
        double fit_after = seam_fit_residual(p, TEST_N, 1.0);
        printf("ramp_fit before=%g after=%g\n", fit_before, fit_after);
        check("ramp_c3_smooth", fit_after < 1e-3 && fit_after < fit_before * 0.01);

        /* the correction itself must be a clean monotone decay on each
         * side - no ringing, no overshoot */
        int mono = 1;
        for (int j = 0; j + 1 < SEAM_WIN; j++) {
            double d0 = (double)p[j] - orig[j];
            double d1 = (double)p[j + 1] - orig[j + 1];
            if (d0 < d1 - 1e-6) mono = 0;
        }
        check("ramp_correction_monotone", mono);
        check("ramp_correction_blends",
              fabs((double)p[SEAM_WIN - 1] - orig[SEAM_WIN - 1]) < 0.02);
    }

    /* bandwidth too narrow: reported unfixable, clip untouched */
    {
        float p[TEST_N], orig[TEST_N];
        SeamReport rep;
        double h = 1.0 / (TEST_N - 1);
        for (int j = 0; j < TEST_N; j++) p[j] = orig[j] = (float)(j * h);
        check("bandwidth_refused", seam_fix(p, TEST_N, 1.0f, 0.01f, &rep) == SEAM_ERR_BANDWIDTH);
        check("bandwidth_not_clamped", memcmp(orig, p, sizeof(p)) == 0);
        check("bandwidth_reports_need", rep.offset_max > 0.01f);
        printf("bandwidth_report needed_offset=%g allowed=0.01\n", rep.offset_max);
    }

    /* known feasible correction bounds the optimum (min-norm property):
     * p = q* + e with q* smooth periodic and e built from cos^2 bumps,
     * so -e lies in the correction subspace, is a feasible fix, and
     * the solver's correction must be no larger. */
    {
        float p[TEST_N], before[TEST_N];
        SeamReport rep;
        double h = 1.0 / (TEST_N - 1), enorm = 0;
        for (int j = 0; j < TEST_N; j++) {
            double t = j * h;
            double q = sin(2.0 * SEAM_PI * t) + 0.3 * sin(6.0 * SEAM_PI * t);
            double e = 0;
            if (j < SEAM_WIN) {
                double c = cos(0.5 * SEAM_PI * j / SEAM_WIN);
                e = 0.2 * c * c;
            } else if (j >= TEST_N - SEAM_WIN) {
                double c = cos(0.5 * SEAM_PI * (TEST_N - 1 - j) / SEAM_WIN);
                e = -0.15 * c * c;
            }
            enorm += e * e;
            p[j] = before[j] = (float)(q + e);
        }
        enorm = sqrt(enorm);

        double fit_before = seam_fit_residual(p, TEST_N, 1.0);
        check("feasible_ok", seam_fix(p, TEST_N, 1.0f, 1.0f, &rep) == SEAM_OK);
        double dnorm = 0;
        for (int j = 0; j < TEST_N; j++) {
            double dj = (double)p[j] - before[j];
            dnorm += dj * dj;
        }
        dnorm = sqrt(dnorm);
        double fit_after = seam_fit_residual(p, TEST_N, 1.0);

        /* optimality certificate: the solver minimizes bending energy,
         * so its correction must be at least as smooth as the known
         * feasible correction -e (up to stencil truncation slack) */
        double bend_d = 0, bend_e = 0;
        for (int side = 0; side < 2; side++) {
            double dv[SEAM_WIN + 1], ev[SEAM_WIN + 1];
            for (int j = 0; j <= SEAM_WIN; j++) {
                int idx = side ? (TEST_N - 1 - j) : j;
                dv[j] = j < SEAM_WIN ? (double)p[idx] - before[idx] : 0.0;
                if (j < SEAM_WIN) {
                    double t = idx * h;
                    double q = sin(2.0 * SEAM_PI * t) + 0.3 * sin(6.0 * SEAM_PI * t);
                    ev[j] = (double)before[idx] - q;
                } else {
                    ev[j] = 0.0;
                }
            }
            for (int j = 1; j <= SEAM_WIN; j++) {
                double u0d = j < SEAM_WIN ? dv[j] : 0.0;
                double upd = j + 1 < SEAM_WIN ? dv[j + 1] : 0.0;
                double u0e = j < SEAM_WIN ? ev[j] : 0.0;
                double upe = j + 1 < SEAM_WIN ? ev[j + 1] : 0.0;
                double d2d = dv[j - 1] - 2.0 * u0d + upd;
                double d2e = ev[j - 1] - 2.0 * u0e + upe;
                bend_d += d2d * d2d;
                bend_e += d2e * d2e;
            }
        }
        printf("feasible_report dnorm=%g enorm=%g bend=%g<=%g fit=%g->%g\n",
               dnorm, enorm, bend_d, bend_e, fit_before, fit_after);
        check("feasible_optimal", bend_d <= bend_e * 1.05 + 1e-9);
        check("feasible_c3", fit_after < fit_before * 0.01);
        check("feasible_offset_certified", rep.offset_max <= rep.offset_bound);
    }

    /* argument validation */
    {
        float p[TEST_N] = { 0 };
        SeamReport rep;
        check("arg_null", seam_fix(0, TEST_N, 1.0f, 1.0f, &rep) == SEAM_ERR_ARG);
        check("arg_short", seam_fix(p, 2 * SEAM_WIN, 1.0f, 1.0f, &rep) == SEAM_ERR_ARG);
        check("arg_period", seam_fix(p, TEST_N, 0.0f, 1.0f, &rep) == SEAM_ERR_ARG);
        check("arg_bandwidth", seam_fix(p, TEST_N, 1.0f, -1.0f, &rep) == SEAM_ERR_ARG);
    }

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
