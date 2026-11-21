
/* unified: hard and soft constraints solved in one stiffness system.
 *
 * Soft constraints are stiffness terms.  Their stiffness matrix is block
 * diagonal by construction, and the soft solve is done block by block.
 * Hard constraints are linear projections C d = b, imposed exactly, once.
 *
 * The unified solution of
 *     minimise 1/2 d^T K d - f^T d   subject to   C d = b
 * is
 *     d* = d_s - K^-1 C^T S^-1 (C d_s - b),   S = C K^-1 C^T
 * where d_s = K^-1 f is the soft-only solve.  The correction is driven by
 * the *residual* the soft solve leaves on the hard manifold: the hard
 * constraint never pulls a soft block away from its own equilibrium.
 *
 * The naive alternative (solve separately, then add)
 *     d_n = d_s + K^-1 C^T S^-1 b
 * adds the hard displacement as if the soft load were zero.  It is in
 * equilibrium (K d_n - f + C^T psi_n = 0) but violates the hard constraint
 * by exactly C d_s.  The two solutions differ by
 *     delta = d* - d_n = -K^-1 C^T S^-1 (C d_s)
 *     ||delta||_K^2 = (C d_s)^T S^-1 (C d_s) <= ||C d_s||^2 / lambda_min(S)
 * which is the bound reported here.  A cruder bound follows from
 *     ||delta||_K <= ||C d_s|| * sigma_max(C) * sqrt(kappa(K)) / lambda_min(S)
 * and shows why an ill-conditioned K inflates any gap.
 */
#include <stdio.h>
#include <string.h>
#include <math.h>

#define MAX_DOF   24
#define MAX_SPR   64
#define MAX_HARD  6
#define MAX_BLOCK 8

typedef struct { int a, b; double k; } Spring;
typedef struct { double c[MAX_DOF]; double target; } Hard;

typedef struct {
    int ndof;
    Spring spr[MAX_SPR];
    int nspr;
    double load[MAX_DOF];
    int fixed[MAX_DOF];
    Hard hard[MAX_HARD];
    int nhard;
    int block[MAX_BLOCK][MAX_DOF];
    int bsize[MAX_BLOCK];
    int nblock;
} Model;

typedef struct {
    int ok;
    const char* why;

    double d[MAX_DOF];            /* unified solution */
    double d_naive[MAX_DOF];      /* separately solved, then added */
    double d_soft[MAX_DOF];       /* soft-only solve */
    double psi[MAX_HARD];         /* hard multipliers C d = b */

    double hard_residual;         /* max |C d - b|, unified */
    double hard_residual_naive;   /* max |C d_n - b| */
    double balance;               /* max |K d - f + C^T psi| */
    double block_delta;           /* max |block solve - dense solve| */

    double cond_k;                /* lambda_max(K) / lambda_min(K) */
    double cond_s;                /* lambda_max(S) / lambda_min(S) */
    double diff_k;                /* ||d* - d_n||_K, actual */
    double bound;                 /* ||C d_s|| / sqrt(lambda_min(S)) */
    double bound_cond;            /* conservative bound using kappa(K) */
} Solution;

static void model_init(Model* m, int ndof) {
    memset(m, 0, sizeof(*m));
    m->ndof = ndof;
}

static int model_spring(Model* m, int a, int b, double k) {
    if (m->nspr >= MAX_SPR) return -1;
    if (a < 0 || a >= m->ndof || b < 0 || b >= m->ndof) return -1;
    Spring* s = &m->spr[m->nspr];
    s->a = a; s->b = b; s->k = k;
    return m->nspr++;
}

static int model_hard(Model* m, const double* c, double target) {
    if (m->nhard >= MAX_HARD) return -1;
    Hard* h = &m->hard[m->nhard];
    memset(h->c, 0, sizeof(h->c));
    for (int i = 0; i < m->ndof; i++) h->c[i] = c[i];
    h->target = target;
    return m->nhard++;
}

static int model_block(Model* m, const int* idx, int n) {
    if (m->nblock >= MAX_BLOCK || n <= 0 || n > MAX_DOF) return -1;
    memcpy(m->block[m->nblock], idx, (size_t)n * sizeof(int));
    m->bsize[m->nblock] = n;
    return m->nblock++;
}

/* ---- dense linear algebra ------------------------------------------- */

/* in-place Cholesky, lower triangle incl. diagonal; 0 on success */
static int chol_factor(double A[MAX_DOF][MAX_DOF], int n) {
    double maxd = 0.0;
    for (int i = 0; i < n; i++) if (A[i][i] > maxd) maxd = A[i][i];
    double tol = 1e-12 * maxd;
    for (int j = 0; j < n; j++) {
        double d = A[j][j];
        for (int k = 0; k < j; k++) d -= A[j][k] * A[j][k];
        if (!(d > tol)) return -1;
        double s = sqrt(d);
        A[j][j] = s;
        for (int i = j + 1; i < n; i++) {
            double v = A[i][j];
            for (int k = 0; k < j; k++) v -= A[i][k] * A[j][k];
            A[i][j] = v / s;
        }
    }
    return 0;
}

static void chol_solve(const double L[MAX_DOF][MAX_DOF], int n,
                       const double* b, double* x) {
    for (int i = 0; i < n; i++) {
        double v = b[i];
        for (int k = 0; k < i; k++) v -= L[i][k] * x[k];
        x[i] = v / L[i][i];
    }
    for (int i = n - 1; i >= 0; i--) {
        double v = x[i];
        for (int k = i + 1; k < n; k++) v -= L[k][i] * x[k];
        x[i] = v / L[i][i];
    }
}

/* symmetric eigenvalues by cyclic Jacobi (A is copied, not modified) */
static void sym_eig(const double A[MAX_DOF][MAX_DOF], int n, double* ev) {
    double a[MAX_DOF][MAX_DOF];
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) a[i][j] = A[i][j];
    for (int sweep = 0; sweep < 80; sweep++) {
        double off = 0.0;
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++) off += a[i][j] * a[i][j];
        if (off < 1e-30) break;
        for (int p = 0; p < n - 1; p++) {
            for (int q = p + 1; q < n; q++) {
                double apq = a[p][q];
                if (fabs(apq) < 1e-18) continue;
                double tau = (a[q][q] - a[p][p]) / (2.0 * apq);
                double t = (tau >= 0.0 ? 1.0 : -1.0) /
                           (fabs(tau) + sqrt(tau * tau + 1.0));
                double c = 1.0 / sqrt(t * t + 1.0);
                double s = t * c;
                for (int k = 0; k < n; k++) {
                    double akp = a[k][p], akq = a[k][q];
                    a[k][p] = c * akp - s * akq;
                    a[k][q] = s * akp + c * akq;
                }
                for (int k = 0; k < n; k++) {
                    double apk = a[p][k], aqk = a[q][k];
                    a[p][k] = c * apk - s * aqk;
                    a[q][k] = s * apk + c * aqk;
                }
            }
        }
    }
    for (int i = 0; i < n; i++) ev[i] = a[i][i];
}

static double dot(const double* a, const double* b, int n) {
    double s = 0.0;
    for (int i = 0; i < n; i++) s += a[i] * b[i];
    return s;
}

/* ---- unified solve --------------------------------------------------- */

static int solve_unified(const Model* m, Solution* sol) {
    memset(sol, 0, sizeof(*sol));
    sol->ok = 0;

    int fmap[MAX_DOF];
    int nf = 0;
    for (int i = 0; i < m->ndof; i++) fmap[i] = m->fixed[i] ? -1 : nf++;
    if (nf == 0) { sol->why = "no free dofs"; return -1; }

    double K[MAX_DOF][MAX_DOF];
    double f[MAX_DOF];
    memset(K, 0, sizeof(K));
    memset(f, 0, sizeof(f));
    for (int i = 0; i < m->ndof; i++)
        if (fmap[i] >= 0) f[fmap[i]] = m->load[i];

    for (int i = 0; i < m->nspr; i++) {
        const Spring* sp = &m->spr[i];
        if (sp->k < 0.0) { sol->why = "negative stiffness"; return -1; }
        int ia = fmap[sp->a], ib = fmap[sp->b];
        if (sp->a == sp->b) {                    /* spring to ground */
            if (ia >= 0) K[ia][ia] += sp->k;
            continue;
        }
        if (ia >= 0) K[ia][ia] += sp->k;
        if (ib >= 0) K[ib][ib] += sp->k;
        if (ia >= 0 && ib >= 0) {
            K[ia][ib] -= sp->k;
            K[ib][ia] -= sp->k;
        }
    }

    double ev[MAX_DOF];
    sym_eig(K, nf, ev);
    double lmax = ev[0], lmin = ev[0];
    for (int i = 1; i < nf; i++) {
        if (ev[i] > lmax) lmax = ev[i];
        if (ev[i] < lmin) lmin = ev[i];
    }
    if (!(lmin > 0.0)) { sol->why = "stiffness not positive definite"; return -1; }
    sol->cond_k = lmax / lmin;

    double L[MAX_DOF][MAX_DOF];
    memcpy(L, K, sizeof(K));
    if (chol_factor(L, nf) != 0) {
        sol->why = "stiffness not positive definite";
        return -1;
    }

    double ds[MAX_DOF];
    double ds_dense[MAX_DOF];
    chol_solve(L, nf, f, ds_dense);

    if (m->nblock > 0) {
        int bof[MAX_DOF];
        for (int i = 0; i < nf; i++) bof[i] = -1;
        for (int b = 0; b < m->nblock; b++) {
            if (m->bsize[b] <= 0 || m->bsize[b] > nf) {
                sol->why = "invalid block size"; return -1;
            }
            for (int k = 0; k < m->bsize[b]; k++) {
                int dof = m->block[b][k];
                if (dof < 0 || dof >= m->ndof || fmap[dof] < 0) {
                    sol->why = "block contains fixed dof"; return -1;
                }
                int g = fmap[dof];
                if (bof[g] >= 0) { sol->why = "block partition overlaps"; return -1; }
                bof[g] = b;
            }
        }
        for (int i = 0; i < nf; i++)
            if (bof[i] < 0) { sol->why = "blocks do not cover all free dofs"; return -1; }
        for (int i = 0; i < m->nspr; i++) {
            const Spring* sp = &m->spr[i];
            int ia = fmap[sp->a], ib = fmap[sp->b];
            if (ia >= 0 && ib >= 0 && ia != ib && bof[ia] != bof[ib]) {
                sol->why = "blocks cut a stiffness coupling"; return -1;
            }
        }

        for (int b = 0; b < m->nblock; b++) {
            int n = m->bsize[b];
            int g[MAX_DOF];
            double A[MAX_DOF][MAX_DOF];
            double rhs[MAX_DOF], out[MAX_DOF];
            for (int i = 0; i < n; i++) g[i] = fmap[m->block[b][i]];
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) A[i][j] = K[g[i]][g[j]];
                rhs[i] = f[g[i]];
            }
            if (chol_factor(A, n) != 0) {
                sol->why = "soft block not positive definite"; return -1;
            }
            chol_solve(A, n, rhs, out);
            for (int i = 0; i < n; i++) ds[g[i]] = out[i];
        }

        double worst = 0.0;
        for (int i = 0; i < nf; i++) {
            double e = fabs(ds[i] - ds_dense[i]);
            if (e > worst) worst = e;
        }
        sol->block_delta = worst;
        if (worst > 1e-8) {
            sol->why = "block solve disagrees with dense solve";
            return -1;
        }
    } else {
        for (int i = 0; i < nf; i++) ds[i] = ds_dense[i];
    }

    /* hard constraints over free dofs */
    int h = m->nhard;
    double C[MAX_DOF][MAX_DOF];
    memset(C, 0, sizeof(C));
    for (int i = 0; i < h; i++) {
        for (int k = 0; k < m->ndof; k++) {
            if (m->hard[i].c[k] == 0.0) continue;
            if (fmap[k] < 0) continue;   /* fixed dofs have zero displacement */
            C[i][fmap[k]] = m->hard[i].c[k];
        }
        if (dot(C[i], C[i], nf) < 1e-30) {
            sol->why = "hard constraint has no free dof"; return -1;
        }
    }

    for (int i = 0; i < nf; i++) sol->d[i] = ds[i];

    if (h > 0) {
        /* S = C K^-1 C^T */
        double KC[MAX_DOF][MAX_DOF];     /* K^-1 C^T, columns */
        double S[MAX_DOF][MAX_DOF];
        double tmp[MAX_DOF];
        for (int j = 0; j < h; j++) {
            chol_solve(L, nf, C[j], tmp);
            for (int k = 0; k < nf; k++) KC[k][j] = tmp[k];
        }
        for (int i = 0; i < h; i++) {
            for (int j = 0; j < h; j++) {
                double v = 0.0;
                for (int k = 0; k < nf; k++) v += C[i][k] * KC[k][j];
                S[i][j] = v;
            }
        }

        double evs[MAX_DOF];
        sym_eig(S, h, evs);
        double smax = evs[0], smin = evs[0];
        for (int i = 1; i < h; i++) {
            if (evs[i] > smax) smax = evs[i];
            if (evs[i] < smin) smin = evs[i];
        }
        if (!(smin > 0.0) || smin <= 1e-12 * smax) {
            sol->why = "hard constraints linearly dependent"; return -1;
        }
        sol->cond_s = smax / smin;

        double CC[MAX_DOF][MAX_DOF];
        for (int i = 0; i < h; i++)
            for (int j = 0; j < h; j++)
                CC[i][j] = dot(C[i], C[j], nf);
        double evc[MAX_DOF];
        sym_eig(CC, h, evc);
        double sigmax_c = 0.0;
        for (int i = 0; i < h; i++) if (evc[i] > sigmax_c) sigmax_c = evc[i];
        sigmax_c = sqrt(sigmax_c);

        double Sf[MAX_DOF][MAX_DOF];
        memcpy(Sf, S, sizeof(S));
        if (chol_factor(Sf, h) != 0) {
            sol->why = "hard constraint system not factorable"; return -1;
        }

        double r[MAX_DOF], psi[MAX_DOF], g[MAX_DOF], corr[MAX_DOF];
        for (int i = 0; i < h; i++)
            r[i] = dot(C[i], ds, nf) - m->hard[i].target;
        chol_solve(Sf, h, r, psi);
        for (int i = 0; i < h; i++) sol->psi[i] = psi[i];

        memset(g, 0, sizeof(g));
        for (int i = 0; i < h; i++)
            for (int k = 0; k < nf; k++) g[k] += C[i][k] * psi[i];
        chol_solve(L, nf, g, corr);
        for (int k = 0; k < nf; k++) sol->d[k] = ds[k] - corr[k];

        /* naive: same hard projection but with b, soft load ignored */
        double rb[MAX_DOF], pn[MAX_DOF], gn[MAX_DOF], cn[MAX_DOF];
        for (int i = 0; i < h; i++) rb[i] = m->hard[i].target;
        chol_solve(Sf, h, rb, pn);
        memset(gn, 0, sizeof(gn));
        for (int i = 0; i < h; i++)
            for (int k = 0; k < nf; k++) gn[k] += C[i][k] * pn[i];
        chol_solve(L, nf, gn, cn);
        /* d_s + K^-1 C^T S^-1 b: hard displacement added on top of the
           soft equilibrium, as if the soft load contributed nothing */
        for (int k = 0; k < nf; k++) sol->d_naive[k] = ds[k] + cn[k];

        /* residuals */
        double worst = 0.0, worst_n = 0.0;
        for (int i = 0; i < h; i++) {
            double res = fabs(dot(C[i], sol->d, nf) - m->hard[i].target);
            double res_n = fabs(dot(C[i], sol->d_naive, nf) - m->hard[i].target);
            if (res > worst) worst = res;
            if (res_n > worst_n) worst_n = res_n;
        }
        sol->hard_residual = worst;
        sol->hard_residual_naive = worst_n;

        double Kd[MAX_DOF], bnd = 0.0;
        for (int k = 0; k < nf; k++) {
            double v = 0.0;
            for (int j = 0; j < nf; j++) v += K[k][j] * sol->d[j];
            Kd[k] = v - f[k] + g[k];
            if (fabs(Kd[k]) > bnd) bnd = fabs(Kd[k]);
        }
        sol->balance = bnd;

        double dlt[MAX_DOF], q = 0.0;
        for (int k = 0; k < nf; k++) dlt[k] = sol->d[k] - sol->d_naive[k];
        for (int k = 0; k < nf; k++) {
            double v = 0.0;
            for (int j = 0; j < nf; j++) v += K[k][j] * dlt[j];
            q += dlt[k] * v;
        }
        sol->diff_k = sqrt(q > 0.0 ? q : 0.0);

        /* ||delta||_K^2 = (C ds)^T S^-1 (C ds) <= ||C ds||^2 / lambda_min(S) */
        double r2[MAX_DOF];
        for (int i = 0; i < h; i++) r2[i] = dot(C[i], ds, nf);
        double rn = sqrt(dot(r2, r2, h));
        sol->bound = rn / sqrt(smin);
        sol->bound_cond = rn * sigmax_c * sqrt(sol->cond_k) / smin;
    } else {
        for (int k = 0; k < nf; k++) sol->d_naive[k] = ds[k];
        sol->hard_residual = 0.0;
        sol->hard_residual_naive = 0.0;
        sol->bound = 0.0;
        sol->bound_cond = 0.0;

        double Kd[MAX_DOF], bnd = 0.0;
        for (int k = 0; k < nf; k++) {
            double v = 0.0;
            for (int j = 0; j < nf; j++) v += K[k][j] * sol->d[j];
            Kd[k] = v - f[k];
            if (fabs(Kd[k]) > bnd) bnd = fabs(Kd[k]);
        }
        sol->balance = bnd;
    }

    for (int i = 0; i < nf; i++) {
        sol->d_soft[i] = ds[i];
    }
    /* expand back to model dof numbering (fixed dofs stay at zero) */
    {
        double dfull[MAX_DOF], nfull[MAX_DOF], sfull[MAX_DOF];
        memset(dfull, 0, sizeof(dfull));
        memset(nfull, 0, sizeof(nfull));
        memset(sfull, 0, sizeof(sfull));
        for (int i = 0; i < m->ndof; i++) {
            if (fmap[i] < 0) continue;
            dfull[i] = sol->d[fmap[i]];
            nfull[i] = sol->d_naive[fmap[i]];
            sfull[i] = ds[fmap[i]];
        }
        memcpy(sol->d, dfull, sizeof(dfull));
        memcpy(sol->d_naive, nfull, sizeof(nfull));
        memcpy(sol->d_soft, sfull, sizeof(sfull));
    }

    sol->ok = 1;
    return 0;
}

/* ---- self tests ------------------------------------------------------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static int near(double a, double b, double tol) { return fabs(a - b) <= tol; }

/* two weakly coupled soft blocks pulled by two hard constraints; this is
   where "separately solved, then added" drags the soft block off. */
static void build_mixed(Model* m) {
    model_init(m, 5);
    model_spring(m, 0, 0, 1000.0);     /* anchor */
    model_spring(m, 0, 1, 100.0);
    model_spring(m, 1, 2, 100.0);
    model_spring(m, 3, 3, 100.0);      /* anchor */
    model_spring(m, 3, 4, 10.0);
    m->load[0] = 1.0; m->load[1] = -2.0; m->load[2] = 3.0;
    m->load[3] = 4.0; m->load[4] = -5.0;
    int b0[] = { 0, 1, 2 }, b1[] = { 3, 4 };
    model_block(m, b0, 3);
    model_block(m, b1, 2);
    double c0[] = { 1, 0, 0, 0, 0 };
    double c1[] = { 0, 1, 0, 0, -1 };
    model_hard(m, c0, 0.4);
    model_hard(m, c1, 0.1);
}

static void report(const char* title, const Model* m) {
    Solution s;
    if (solve_unified(m, &s) != 0) {
        printf("%s: rejected: %s\n", title, s.why);
        return;
    }
    printf("== %s ==\n", title);
    printf("dof  soft-only   unified     naive\n");
    for (int i = 0; i < m->ndof; i++)
        printf("%3d  %9.6f  %9.6f  %9.6f\n",
               i, s.d_soft[i], s.d[i], s.d_naive[i]);
    printf("hard_residual unified=%.3e naive=%.3e\n",
           s.hard_residual, s.hard_residual_naive);
    printf("balance=%.3e cond(K)=%.3f cond(S)=%.3f\n",
           s.balance, s.cond_k, s.cond_s);
    printf("diff_K=%.6f bound=%.6f bound_cond=%.3f\n",
           s.diff_k, s.bound, s.bound_cond);
}

int main(void) {
    Model m;
    Solution s;

    build_mixed(&m);
    report("mixed soft/hard", &m);

    check("mixed_ok", solve_unified(&m, &s) == 0 && s.ok);
    check("blocks_match_dense", near(s.block_delta, 0.0, 1e-9));
    check("unified_satisfies_hard", s.hard_residual < 1e-9);
    check("naive_gets_dragged", near(s.hard_residual_naive, 0.522, 1e-3));
    check("equilibrium", s.balance < 1e-8);
    check("bound_holds", s.diff_k <= s.bound + 1e-9);
    check("cond_reported", s.cond_k > 1.0 && s.cond_s > 1.0);
    check("soft_solve_value", near(s.d_soft[0], 0.002, 1e-9) &&
                              near(s.d_soft[4], -0.51, 1e-9));
    check("unified_value", near(s.d[0], 0.4, 1e-9) &&
                           near(s.d[1], 0.3416666667, 1e-9) &&
                           near(s.d[4], 0.2416666667, 1e-9));
    check("naive_value", near(s.d_naive[0], 0.402, 1e-9) &&
                         near(s.d_naive[4], -0.235, 1e-9));
    check("gap_nonzero", s.diff_k > 1.0);

    /* negative stiffness is refused before any solve */
    Model neg;
    build_mixed(&neg);
    neg.spr[1].k = -100.0;
    check("reject_negative_k", solve_unified(&neg, &s) != 0 && !s.ok);

    /* hard constraints that are scalar multiples leave S singular */
    Model dep;
    build_mixed(&dep);
    for (int k = 0; k < 5; k++) dep.hard[1].c[k] = 2.0 * dep.hard[0].c[k];
    dep.hard[1].target = 0.8;
    check("reject_dependent_hard", solve_unified(&dep, &s) != 0 && !s.ok);

    /* a block split that cuts a spring is refused */
    Model cut;
    build_mixed(&cut);
    cut.block[0][2] = 3;         /* block {0,1,3} cuts spring 1-2 */
    check("reject_block_cut", solve_unified(&cut, &s) != 0 && !s.ok);

    /* a floating dof makes K singular */
    Model flo;
    model_init(&flo, 2);
    model_spring(&flo, 0, 0, 10.0);
    check("reject_singular_k", solve_unified(&flo, &s) != 0 && !s.ok);

    /* one dof: the soft solve already meets the hard target, so the naive
       answer still adds the full hard displacement and overshoots */
    Model one;
    model_init(&one, 1);
    model_spring(&one, 0, 0, 100.0);
    one.load[0] = 50.0;
    double c1[] = { 1 };
    model_hard(&one, c1, 0.5);
    check("one_ok", solve_unified(&one, &s) == 0 && s.ok);
    check("one_unified", near(s.d[0], 0.5, 1e-12) && s.hard_residual < 1e-12);
    check("one_naive", near(s.d_naive[0], 1.0, 1e-12) &&
                       near(s.hard_residual_naive, 0.5, 1e-12));
    /* ||delta||_K = |C ds| / sqrt(S); here S is tight, so bound == gap */
    check("one_bound_holds", near(s.diff_k, 5.0, 1e-9) &&
                             near(s.bound, 5.0, 1e-9));
    /* fixed dof: excluded from the system, stays at zero */
    Model fix;
    model_init(&fix, 2);
    fix.fixed[0] = 1;
    model_spring(&fix, 0, 0, 1000.0);
    model_spring(&fix, 0, 1, 100.0);
    fix.load[1] = 10.0;
    check("fixed_ok", solve_unified(&fix, &s) == 0 && s.ok);
    check("fixed_value", near(s.d[0], 0.0, 1e-15) && near(s.d[1], 0.1, 1e-12));
    check("fixed_cond", near(s.cond_k, 1.0, 1e-9));

    /* no hard constraints: unified is the soft solve, gap is zero */
    Model soft;
    model_init(&soft, 1);
    model_spring(&soft, 0, 0, 100.0);
    soft.load[0] = 50.0;
    check("soft_only", solve_unified(&soft, &s) == 0 && near(s.d[0], 0.5, 1e-12) &&
                       near(s.diff_k, 0.0, 1e-15) && near(s.bound, 0.0, 1e-15));

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
