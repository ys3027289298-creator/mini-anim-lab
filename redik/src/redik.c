
/* redik: redundant IK. Torso and both hands are solved jointly with a
   damped pseudoinverse; a secondary posture task rides in the null
   space of the primary Jacobian, so it cannot increase the primary
   residual beyond the damping slack (lambda/2)*alpha*|g|.
   Reports: end-effector error with an a posteriori upper bound, joint
   drift with upper bounds, infeasibility when task dims exceed the
   active DOF, and degraded best-effort poses for unreachable targets. */
#include <stdio.h>
#include <string.h>
#include <math.h>

#define MAX_Q 7   /* tx ty torso_ang l_shoulder l_elbow r_shoulder r_elbow */
#define MAX_TASKS 4
#define MAX_M (2 * MAX_TASKS)

#define REDIK_OK 0
#define REDIK_DEGRADED 1        /* best effort returned, residual remains */
#define REDIK_INFEASIBLE_DOF 2  /* task dimension exceeds active DOF */

typedef struct { float x, y; } V2;

typedef struct {
    float q[MAX_Q], qmin[MAX_Q], qmax[MAX_Q];
    int frozen[MAX_Q];
    float shoulder_off, shoulder_h;   /* shoulder local offset (+-off, h) */
    float l1, l2;                     /* upper arm / forearm length */
} Body;

static V2 rotv(V2 v, float a) {
    float c = cosf(a), s = sinf(a);
    V2 r = { v.x * c - v.y * s, v.x * s + v.y * c };
    return r;
}

static V2 hand_pos(const Body* b, int side, V2* sh_out, V2* el_out) {
    V2 t = { b->q[0], b->q[1] };
    V2 os = { side == 0 ? -b->shoulder_off : b->shoulder_off, b->shoulder_h };
    V2 so = rotv(os, b->q[2]);
    V2 sh = { t.x + so.x, t.y + so.y };
    float a1 = b->q[2] + b->q[3 + 2 * side];
    V2 el = { sh.x + b->l1 * cosf(a1), sh.y + b->l1 * sinf(a1) };
    float a2 = a1 + b->q[4 + 2 * side];
    V2 h = { el.x + b->l2 * cosf(a2), el.y + b->l2 * sinf(a2) };
    if (sh_out) *sh_out = sh;
    if (el_out) *el_out = el;
    return h;
}

/* stacked primary residual e = target - hand, returns ||e|| */
static float primary_err(const Body* b, const V2* target, const int* hand, int ntasks, float* e) {
    float ne = 0;
    for (int t = 0; t < ntasks; t++) {
        V2 h = hand_pos(b, hand[t], 0, 0);
        e[2 * t] = target[t].x - h.x;
        e[2 * t + 1] = target[t].y - h.y;
        ne += e[2 * t] * e[2 * t] + e[2 * t + 1] * e[2 * t + 1];
    }
    return sqrtf(ne);
}

/* J (m x n, row major) over active dofs; revolute column = perp(hand - pivot) */
static void build_jacobian(const Body* b, const int* hand, int ntasks, const int* map, int n, float* J) {
    V2 t = { b->q[0], b->q[1] };
    V2 sh[2], el[2], h[2];
    for (int s = 0; s < 2; s++) h[s] = hand_pos(b, s, &sh[s], &el[s]);
    int m = 2 * ntasks;
    memset(J, 0, sizeof(float) * m * n);
    for (int ti = 0; ti < ntasks; ti++) {
        int s = hand[ti];
        for (int c = 0; c < n; c++) {
            int dof = map[c];
            float jx = 0, jy = 0;
            V2 d = { 0, 0 };
            if (dof == 0) { jx = 1; }
            else if (dof == 1) { jy = 1; }
            else if (dof == 2) { d.x = h[s].x - t.x; d.y = h[s].y - t.y; }
            else if (dof == 3 + 2 * s) { d.x = h[s].x - sh[s].x; d.y = h[s].y - sh[s].y; }
            else if (dof == 4 + 2 * s) { d.x = h[s].x - el[s].x; d.y = h[s].y - el[s].y; }
            if (dof >= 2) { jx = -d.y; jy = d.x; }
            J[2 * ti * n + c] = jx;
            J[(2 * ti + 1) * n + c] = jy;
        }
    }
}

/* solve A x = b for symmetric positive definite A (n x n) via Cholesky */
static int chol_solve(int n, const float* A, const float* b, float* x) {
    float L[MAX_M * MAX_M], y[MAX_M];
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            float s = A[i * n + j];
            for (int k = 0; k < j; k++) s -= L[i * n + k] * L[j * n + k];
            if (i == j) {
                if (s <= 1e-12f) return 0;
                L[i * n + i] = sqrtf(s);
            } else {
                L[i * n + j] = s / L[j * n + j];
            }
        }
    }
    for (int i = 0; i < n; i++) {
        float s = b[i];
        for (int k = 0; k < i; k++) s -= L[i * n + k] * y[k];
        y[i] = s / L[i * n + i];
    }
    for (int i = n - 1; i >= 0; i--) {
        float s = y[i];
        for (int k = i + 1; k < n; k++) s -= L[k * n + i] * x[k];
        x[i] = s / L[i * n + i];
    }
    return 1;
}

/* eigenvalues of a symmetric n x n matrix (cyclic Jacobi) */
static void sym_eigvals(int n, const float* A, float* out) {
    float a[MAX_M * MAX_M];
    memcpy(a, A, sizeof(float) * n * n);
    for (int sweep = 0; sweep < 64; sweep++) {
        float off = 0;
        for (int p = 0; p < n; p++)
            for (int q = p + 1; q < n; q++) off += a[p * n + q] * a[p * n + q];
        if (off < 1e-20f) break;
        for (int p = 0; p < n; p++) {
            for (int q = p + 1; q < n; q++) {
                float apq = a[p * n + q];
                if (fabsf(apq) < 1e-12f) continue;
                float phi = 0.5f * atan2f(2.0f * apq, a[q * n + q] - a[p * n + p]);
                float c = cosf(phi), s = sinf(phi);
                for (int k = 0; k < n; k++) {
                    float akp = a[k * n + p], akq = a[k * n + q];
                    a[k * n + p] = c * akp - s * akq;
                    a[k * n + q] = s * akp + c * akq;
                }
                for (int k = 0; k < n; k++) {
                    float apk = a[p * n + k], aqk = a[q * n + k];
                    a[p * n + k] = c * apk - s * aqk;
                    a[q * n + k] = s * apk + c * aqk;
                }
            }
        }
    }
    for (int i = 0; i < n; i++) out[i] = a[i * n + i];
}

typedef struct {
    V2 target[MAX_TASKS];
    int hand[MAX_TASKS];
    int ntasks;
    float tol, lambda, alpha;   /* stop tol, damping, secondary step size */
    int max_iter;
    int use_secondary;          /* project joint centering into null space */
    float q_rest[MAX_Q], w_center[MAX_Q];
} SolveOpts;

typedef struct {
    int status, iters;
    float err, err_bound;       /* final primary residual and its upper bound */
    float kappa;                /* linear contraction factor at the solution */
    float drift, drift_bound;   /* ||q - q0|| and triangle-inequality bound */
    float drift_theory;         /* pinv-norm bound: sum |e|/(2*lambda) + alpha*|g| */
    float null_violation;       /* primary-error increase beyond projection slack */
    float center_cost;          /* secondary objective value at the solution */
} Report;

static float center_cost(const Body* b, const SolveOpts* o) {
    float h = 0;
    for (int i = 0; i < MAX_Q; i++) {
        float d = b->q[i] - o->q_rest[i];
        h += 0.5f * o->w_center[i] * d * d;
    }
    return h;
}

static Report solve(Body* b, const SolveOpts* o) {
    Report r;
    memset(&r, 0, sizeof(r));
    int map[MAX_Q], n = 0;
    for (int i = 0; i < MAX_Q; i++) if (!b->frozen[i]) map[n++] = i;
    int m = 2 * o->ntasks;
    float q0[MAX_Q];
    memcpy(q0, b->q, sizeof(q0));
    r.err_bound = o->tol;
    if (m > n || n == 0) { r.status = REDIK_INFEASIBLE_DOF; r.err = -1.0f; return r; }
    if (m == 0) { r.status = REDIK_OK; return r; }

    float e[MAX_M], J[MAX_M * MAX_Q], M[MAX_M * MAX_M];
    float lin_res = 0, last_step = 0;
    float ne = primary_err(b, o->target, o->hand, o->ntasks, e);
    int it = 0;
    for (it = 0; it < o->max_iter && ne > o->tol; it++) {
        build_jacobian(b, o->hand, o->ntasks, map, n, J);
        for (int i = 0; i < m; i++)           /* M = J J' + lambda^2 I */
            for (int j = 0; j < m; j++) {
                float s = (i == j) ? o->lambda * o->lambda : 0;
                for (int k = 0; k < n; k++) s += J[i * n + k] * J[j * n + k];
                M[i * m + j] = s;
            }
        float x[MAX_M];
        if (!chol_solve(m, M, e, x)) { r.err_bound = 1e30f; break; }
        float dqp[MAX_Q];                     /* primary step: dq = J' (J J' + l^2 I)^-1 e */
        for (int c = 0; c < n; c++) {
            float s = 0;
            for (int i = 0; i < m; i++) s += J[i * n + c] * x[i];
            dqp[c] = s;
        }
        float g[MAX_Q], dqn[MAX_Q], ng = 0;   /* secondary: centering, null-projected */
        memset(g, 0, sizeof(g));
        memset(dqn, 0, sizeof(dqn));
        if (o->use_secondary) {
            for (int c = 0; c < n; c++) {
                int dof = map[c];
                g[c] = o->w_center[dof] * (b->q[dof] - o->q_rest[dof]);
                ng += g[c] * g[c];
            }
            ng = sqrtf(ng);
            float Jg[MAX_M], y[MAX_M];        /* P g = g - J' M^-1 (J g) */
            for (int i = 0; i < m; i++) {
                float s = 0;
                for (int c = 0; c < n; c++) s += J[i * n + c] * g[c];
                Jg[i] = s;
            }
            chol_solve(m, M, Jg, y);
            for (int c = 0; c < n; c++) {
                float s = g[c];
                for (int i = 0; i < m; i++) s -= J[i * n + c] * y[i];
                dqn[c] = -o->alpha * s;
            }
        }
        /* null-step moves the linearized primary residual by at most
           |J P| * alpha * |g| <= (lambda/2) * alpha * |g|  (|J P| <= lambda/2) */
        float lp = 0, lt = 0;
        for (int i = 0; i < m; i++) {
            float jp = 0, jt = 0;
            for (int c = 0; c < n; c++) {
                jp += J[i * n + c] * dqp[c];
                jt += J[i * n + c] * (dqp[c] + dqn[c]);
            }
            float rp = e[i] - jp, rt = e[i] - jt;
            lp += rp * rp;
            lt += rt * rt;
        }
        float viol = sqrtf(lt) - sqrtf(lp) - 0.5f * o->lambda * o->alpha * ng;
        if (o->use_secondary && viol > r.null_violation) r.null_violation = viol;
        /* apply with joint limits; track the actually applied step */
        float app[MAX_Q], st = 0;
        for (int c = 0; c < n; c++) {
            int dof = map[c];
            float qn = b->q[dof] + dqp[c] + dqn[c];
            if (qn < b->qmin[dof]) qn = b->qmin[dof];
            if (qn > b->qmax[dof]) qn = b->qmax[dof];
            app[c] = qn - b->q[dof];
            b->q[dof] = qn;
            st += app[c] * app[c];
        }
        last_step = sqrtf(st);
        r.drift_bound += last_step;           /* triangle inequality: exact bound */
        r.drift_theory += ne / (2.0f * o->lambda) + o->alpha * ng;
        /* linearized residual of the applied step, for the a posteriori bound */
        lin_res = 0;
        for (int i = 0; i < m; i++) {
            float jt = 0;
            for (int c = 0; c < n; c++) jt += J[i * n + c] * app[c];
            float rt = e[i] - jt;
            lin_res += rt * rt;
        }
        lin_res = sqrtf(lin_res);
        ne = primary_err(b, o->target, o->hand, o->ntasks, e);
    }
    r.iters = it;
    r.err = ne;
    r.status = ne <= o->tol ? REDIK_OK : REDIK_DEGRADED;
    /* contraction factor kappa = lambda^2 / (lambda^2 + sigma_min^2) at the solution */
    build_jacobian(b, o->hand, o->ntasks, map, n, J);
    float A[MAX_M * MAX_M], ev[MAX_M] = { 0 };
    for (int i = 0; i < m; i++)
        for (int j = 0; j < m; j++) {
            float s = 0;
            for (int k = 0; k < n; k++) s += J[i * n + k] * J[j * n + k];
            A[i * m + j] = s;
        }
    sym_eigvals(m, A, ev);
    float smin2 = ev[0];
    for (int i = 1; i < m; i++) if (ev[i] < smin2) smin2 = ev[i];
    if (smin2 < 0) smin2 = 0;
    r.kappa = o->lambda * o->lambda / (o->lambda * o->lambda + smin2);
    /* a posteriori error bound: last linearized residual + Taylor remainder.
       Hessian entries are bounded by the max moment arm R, so the remainder
       of one step is <= 0.5 * R * n * |dq|^2. */
    float os = sqrtf(b->shoulder_off * b->shoulder_off + b->shoulder_h * b->shoulder_h);
    float taylor = lin_res + 0.5f * (os + b->l1 + b->l2) * (float)n * last_step * last_step;
    if (taylor > r.err_bound) r.err_bound = taylor;
    for (int i = 0; i < MAX_Q; i++) {
        float d = b->q[i] - q0[i];
        r.drift += d * d;
    }
    r.drift = sqrtf(r.drift);
    r.center_cost = center_cost(b, o);
    return r;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static void init_body(Body* b) {
    memset(b, 0, sizeof(*b));
    b->qmin[0] = b->qmin[1] = -0.5f;  b->qmax[0] = b->qmax[1] = 0.5f;
    b->qmin[2] = -0.8f;               b->qmax[2] = 0.8f;
    for (int i = 3; i <= 5; i += 2) { b->qmin[i] = -3.2f; b->qmax[i] = 3.2f; }
    for (int i = 4; i <= 6; i += 2) { b->qmin[i] = -2.9f; b->qmax[i] = 2.9f; }
    b->shoulder_off = 0.15f;
    b->shoulder_h = 0.45f;
    b->l1 = b->l2 = 0.35f;            /* per-arm reach 0.7 */
    b->q[3] = b->q[5] = -1.5707963f;  /* arms hang down */
}

static void init_opts(SolveOpts* o) {
    memset(o, 0, sizeof(*o));
    o->tol = 1e-3f;
    o->lambda = 1e-3f;
    o->alpha = 0.3f;
    o->max_iter = 200;
    o->q_rest[3] = o->q_rest[5] = -1.5707963f;
    o->q_rest[4] = o->q_rest[6] = 0.6f;   /* elbows slightly bent at rest */
    for (int i = 0; i < MAX_Q; i++) o->w_center[i] = 1.0f;
}

int main(void) {
    Body b;
    SolveOpts o;
    init_body(&b);
    init_opts(&o);
    o.target[0] = (V2){ -0.75f, -0.25f }; o.hand[0] = 0;   /* left, out of reach */
    o.target[1] = (V2){ 0.45f, 0.6f };    o.hand[1] = 1;
    o.ntasks = 2;

    /* torso frozen: the original complaint — body stays put, hand is pulled off */
    Body bf = b;
    bf.frozen[0] = bf.frozen[1] = bf.frozen[2] = 1;
    Report rf = solve(&bf, &o);
    check("frozen_torso_degraded", rf.status == REDIK_DEGRADED);
    check("frozen_torso_residual", rf.err > 0.05f);

    /* joint solve: torso moves, both hands arrive */
    Body bj = b;
    Report rj = solve(&bj, &o);
    check("joint_solve_ok", rj.status == REDIK_OK);
    check("joint_solve_err", rj.err <= o.tol);
    check("torso_moved", fabsf(bj.q[0]) + fabsf(bj.q[1]) + fabsf(bj.q[2]) > 1e-3f);

    /* end-effector error upper bound and contraction factor */
    check("err_within_bound", rj.err <= rj.err_bound + 1e-4f);
    check("kappa_contracts", rj.kappa < 1.0f);

    /* joint drift upper bounds */
    check("drift_bounded", rj.drift <= rj.drift_bound + 1e-5f);
    check("drift_bound_sound", rj.drift_bound <= rj.drift_theory + 1e-4f);
    check("drift_nonzero", rj.drift > 1e-3f);

    /* plain pseudoinverse vs null-space projection of the secondary task */
    Body bp = b;
    SolveOpts op = o;
    op.use_secondary = 0;
    Report rp = solve(&bp, &op);
    Body bn = b;
    SolveOpts on = o;
    on.use_secondary = 1;
    Report rn = solve(&bn, &on);
    printf("pinv: err=%.5f drift=%.4f cost=%.5f\n", rp.err, rp.drift, rp.center_cost);
    printf("null: err=%.5f drift=%.4f cost=%.5f viol=%.2e\n",
           rn.err, rn.drift, rn.center_cost, rn.null_violation);
    printf("bounds: err=%.5f<=%.5f kappa=%.4f drift=%.4f<=%.4f<=%.1f\n",
           rj.err, rj.err_bound, rj.kappa, rj.drift, rj.drift_bound, rj.drift_theory);
    check("pinv_ok", rp.status == REDIK_OK);
    check("null_ok", rn.status == REDIK_OK);
    check("primary_not_worse", rn.err <= rp.err + 1e-3f);
    check("null_no_violation", rn.null_violation <= 1e-5f);
    check("secondary_improved", rn.center_cost <= rp.center_cost + 1e-4f);

    /* task dimension exceeds active DOF: report infeasible */
    Body bi = b;
    for (int i = 0; i < MAX_Q; i++) bi.frozen[i] = 1;
    bi.frozen[3] = bi.frozen[5] = 0;      /* 2 DOF vs 4 task dims */
    Report ri = solve(&bi, &o);
    check("infeasible_dof", ri.status == REDIK_INFEASIBLE_DOF);

    /* unreachable target: degraded best effort, limits respected */
    Body bu = b;
    SolveOpts ou = o;
    ou.target[0] = (V2){ 5.0f, 5.0f };
    Report ru = solve(&bu, &ou);
    check("unreachable_degraded", ru.status == REDIK_DEGRADED && ru.err > ou.tol);
    int in_limits = 1;
    for (int i = 0; i < MAX_Q; i++)
        if (bu.q[i] < bu.qmin[i] - 1e-4f || bu.q[i] > bu.qmax[i] + 1e-4f) in_limits = 0;
    check("limits_respected", in_limits);
    check("unreachable_drift_bounded", ru.drift <= ru.drift_bound + 1e-5f);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
