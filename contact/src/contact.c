
/* contact: static friction-contact solver for a single 2d rigid body
   resting on a set of contact points (the ragdoll-lying-on-the-ground
   case). builds the per-contact complementarity conditions
   (fn >= 0, an >= 0, fn*an = 0, |ft| <= mu*fn), solves the friction cone
   problem with projected gauss-seidel, and reports a convergence
   residual. zero-length normals and penetrating contacts are rejected;
   when no static solution exists or the iteration does not converge,
   the conflicting contact pair is reported. */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define CT_MAX_CONTACTS 32
#define CT_NORM_EPS 1e-6f   /* normals shorter than this are rejected */
#define CT_PEN_TOL  1e-4f   /* gap below -CT_PEN_TOL counts as penetration */
#define CT_GAP_TOL  1e-4f   /* gap above this means the contact is open */
#define CT_EQ_TOL   1e-2f   /* equilibrium residual tolerance, force units */

enum {
    CT_OK = 0,
    CT_ERR_ZERO_NORMAL,    /* a contact normal has (near) zero length */
    CT_ERR_PENETRATION,    /* a contact gap is negative beyond tolerance */
    CT_ERR_NO_SOLUTION,    /* no cone-feasible force set balances the body */
    CT_ERR_NOT_CONVERGED   /* projected iteration stalled */
};

typedef struct { float x, y; } ct_vec2;

typedef struct {
    ct_vec2 point;   /* world position of the contact */
    ct_vec2 normal;  /* unit normal, surface -> body; zero normals rejected */
    float mu;        /* friction coefficient */
    float gap;       /* signed distance to the surface; must be >= 0 */
} ct_contact;

typedef struct {
    float mass, inertia;
    ct_vec2 com;
    ct_vec2 gravity;
    ct_vec2 ext_force;  /* extra applied force, at the com */
    float ext_torque;
    ct_contact contacts[CT_MAX_CONTACTS];
    int n_contacts;
    int max_iterations;
    float tolerance;    /* convergence residual tolerance */
} ct_problem;

typedef struct {
    float fn[CT_MAX_CONTACTS];
    float ft[CT_MAX_CONTACTS];
    float residual;     /* final fixed-point residual of the iteration */
    float penetration;  /* max penetration depth over contacts; must be 0 */
    float eq_residual;  /* |sum f + f_ext| and torque imbalance, max norm */
    int iterations;
    int converged;
    int error;          /* CT_OK or one of CT_ERR_* */
    int error_contact;  /* contact that caused a rejection, else -1 */
    int conflict_a, conflict_b; /* conflicting pair on failure, else -1 */
} ct_result;

static float dot2(ct_vec2 a, ct_vec2 b) { return a.x * b.x + a.y * b.y; }
static float cross2(ct_vec2 a, ct_vec2 b) { return a.x * b.y - a.y * b.x; }

/* project (zn, zt) onto the friction cone { n >= 0, |t| <= mu n } */
static void cone_project(float zn, float zt, float mu, float* pn, float* pt) {
    float at = fabsf(zt);
    if (at <= mu * zn) { *pn = zn; *pt = zt; return; }   /* inside */
    if (mu * at <= -zn) { *pn = 0.0f; *pt = 0.0f; return; } /* polar cone */
    float s = zt >= 0.0f ? 1.0f : -1.0f;
    float lam = (zn + mu * at) / (1.0f + mu * mu);       /* onto the rim */
    if (lam <= 0.0f) { *pn = 0.0f; *pt = 0.0f; return; }
    *pn = lam; *pt = s * mu * lam;
}

typedef struct {
    int n;
    int active[CT_MAX_CONTACTS];
    float mu[CT_MAX_CONTACTS];
    float wnn[CT_MAX_CONTACTS][CT_MAX_CONTACTS];
    float wnt[CT_MAX_CONTACTS][CT_MAX_CONTACTS];
    float wtn[CT_MAX_CONTACTS][CT_MAX_CONTACTS];
    float wtt[CT_MAX_CONTACTS][CT_MAX_CONTACTS];
    float an[CT_MAX_CONTACTS], at[CT_MAX_CONTACTS]; /* contact accelerations */
    float fn[CT_MAX_CONTACTS], ft[CT_MAX_CONTACTS];
} ct_workspace;

/* fixed-point residual: max over contacts of |f - proj_cone(f - rho a)|;
   zero iff the complementarity conditions hold exactly */
static float compute_residual(ct_workspace* w, float* cr) {
    float worst = 0.0f;
    for (int i = 0; i < w->n; i++) {
        cr[i] = 0.0f;
        if (!w->active[i]) continue;
        float rho = 1.0f / (w->wnn[i][i] + w->wtt[i][i]);
        float zn = w->fn[i] - rho * w->an[i];
        float zt = w->ft[i] - rho * w->at[i];
        float pn, pt;
        cone_project(zn, zt, w->mu[i], &pn, &pt);
        float dn = w->fn[i] - pn, dt = w->ft[i] - pt;
        cr[i] = sqrtf(dn * dn + dt * dt);
        if (cr[i] > worst) worst = cr[i];
    }
    return worst;
}

/* phase a sweep: normal forces only, friction frozen. the return value is
   the fixed-point residual |fn - max(0, fn - rho an)| of the sweep */
static float normal_sweep(ct_workspace* w) {
    float worst = 0.0f;
    for (int i = 0; i < w->n; i++) {
        if (!w->active[i]) continue;
        float rho = 1.0f / w->wnn[i][i];
        float nfn = w->fn[i] - rho * w->an[i];
        if (nfn < 0.0f) nfn = 0.0f;
        float d = fabsf(nfn - w->fn[i]);
        if (d > worst) worst = d;
        float dn = nfn - w->fn[i];
        w->fn[i] = nfn;
        for (int k = 0; k < w->n; k++) {
            w->an[k] += w->wnn[k][i] * dn;
            w->at[k] += w->wtn[k][i] * dn;
        }
    }
    return worst;
}

/* phase b sweep: friction only, normal forces frozen, box-clamped to the
   cone rim |ft| <= mu fn */
static float friction_sweep(ct_workspace* w) {
    float worst = 0.0f;
    for (int i = 0; i < w->n; i++) {
        if (!w->active[i]) continue;
        float rho = 1.0f / w->wtt[i][i];
        float nft = w->ft[i] - rho * w->at[i];
        float lim = w->mu[i] * w->fn[i];
        if (nft > lim) nft = lim;
        if (nft < -lim) nft = -lim;
        float d = fabsf(nft - w->ft[i]);
        if (d > worst) worst = d;
        float dt = nft - w->ft[i];
        w->ft[i] = nft;
        for (int k = 0; k < w->n; k++) {
            w->an[k] += w->wnt[k][i] * dt;
            w->at[k] += w->wtt[k][i] * dt;
        }
    }
    return worst;
}

/* phase c: one projected gauss-seidel update of contact i: solve the 2x2
   block with the other forces frozen, then keep the result inside the cone */
static void block_update(ct_workspace* w, int i) {
    float a_n = w->an[i] - w->wnn[i][i] * w->fn[i] - w->wnt[i][i] * w->ft[i];
    float a_t = w->at[i] - w->wtn[i][i] * w->fn[i] - w->wtt[i][i] * w->ft[i];
    float wnn = w->wnn[i][i], wnt = w->wnt[i][i];
    float wtn = w->wtn[i][i], wtt = w->wtt[i][i];
    float mu = w->mu[i];
    float nfn = 0.0f, nft = 0.0f;
    int done = 0;

    /* stick: exact 2x2 solve, valid if it lands inside the cone */
    float det = wnn * wtt - wnt * wtn;
    if (det > 1e-12f) {
        float cfn = (wnt * a_t - wtt * a_n) / det;
        float cft = (wtn * a_n - wnn * a_t) / det;
        if (cfn >= 0.0f && fabsf(cft) <= mu * cfn) {
            nfn = cfn; nft = cft; done = 1;
        }
    }
    /* separation: contact opens, no force needed */
    if (!done && a_n >= 0.0f) { nfn = 0.0f; nft = 0.0f; done = 1; }
    /* sliding: ft pinned to the cone rim, solve the normal equation */
    if (!done) {
        for (int s = -1; s <= 1 && !done; s += 2) {
            float denom = wnn + s * mu * wnt;
            if (fabsf(denom) < 1e-12f) continue;
            float f = -a_n / denom;
            float ares = a_t + (wtn + s * mu * wtt) * f;
            if (f > 0.0f && ares * s <= 0.0f) {
                nfn = f; nft = s * mu * f; done = 1;
            }
        }
    }
    /* fallback: clamped scalar update */
    if (!done) {
        nfn = fmaxf(0.0f, -a_n / wnn);
        float lim = mu * nfn;
        nft = -a_t / wtt;
        if (nft > lim) nft = lim;
        if (nft < -lim) nft = -lim;
    }

    float dn = nfn - w->fn[i], dt = nft - w->ft[i];
    w->fn[i] = nfn; w->ft[i] = nft;
    for (int k = 0; k < w->n; k++) {
        w->an[k] += w->wnn[k][i] * dn + w->wnt[k][i] * dt;
        w->at[k] += w->wtn[k][i] * dn + w->wtt[k][i] * dt;
    }
}

/* the two active contacts with the largest residuals */
static void pick_conflict(const float* cr, const int* active, int n,
                          int* a, int* b) {
    int bi = -1, bj = -1;
    for (int i = 0; i < n; i++) {
        if (!active[i]) continue;
        if (bi < 0 || cr[i] > cr[bi]) { bj = bi; bi = i; }
        else if (bj < 0 || cr[i] > cr[bj]) { bj = i; }
    }
    *a = bi; *b = bj;
}

void ct_solve(const ct_problem* prob, ct_result* res) {
    memset(res, 0, sizeof(*res));
    res->error = CT_OK;
    res->error_contact = -1;
    res->conflict_a = -1;
    res->conflict_b = -1;
    int n = prob->n_contacts;

    static ct_workspace w;
    w.n = n;
    ct_vec2 nrm[CT_MAX_CONTACTS], tgt[CT_MAX_CONTACTS], r[CT_MAX_CONTACTS];
    float cn[CT_MAX_CONTACTS], ctl[CT_MAX_CONTACTS];

    /* validate: reject zero normals, reject penetration */
    for (int i = 0; i < n; i++) {
        ct_vec2 nn = prob->contacts[i].normal;
        float len = sqrtf(dot2(nn, nn));
        if (len < CT_NORM_EPS) {
            res->error = CT_ERR_ZERO_NORMAL;
            res->error_contact = i;
            return;
        }
        nrm[i].x = nn.x / len; nrm[i].y = nn.y / len;
        tgt[i].x = -nrm[i].y; tgt[i].y = nrm[i].x;
        r[i].x = prob->contacts[i].point.x - prob->com.x;
        r[i].y = prob->contacts[i].point.y - prob->com.y;
        cn[i] = cross2(r[i], nrm[i]);
        ctl[i] = cross2(r[i], tgt[i]);
        w.mu[i] = prob->contacts[i].mu;
        w.active[i] = prob->contacts[i].gap <= CT_GAP_TOL;
        float pen = -prob->contacts[i].gap;
        if (pen > res->penetration) res->penetration = pen;
    }
    if (res->penetration > CT_PEN_TOL) {
        res->error = CT_ERR_PENETRATION;
        for (int i = 0; i < n; i++) {
            if (-prob->contacts[i].gap > CT_PEN_TOL) {
                res->error_contact = i;
                break;
            }
        }
        return;
    }

    /* delassus operator blocks W = J M^-1 J^T and the free acceleration */
    float inv_m = 1.0f / prob->mass;
    float inv_i = 1.0f / prob->inertia;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            w.wnn[i][j] = dot2(nrm[i], nrm[j]) * inv_m + cn[i] * cn[j] * inv_i;
            w.wnt[i][j] = dot2(nrm[i], tgt[j]) * inv_m + cn[i] * ctl[j] * inv_i;
            w.wtn[i][j] = dot2(tgt[i], nrm[j]) * inv_m + ctl[i] * cn[j] * inv_i;
            w.wtt[i][j] = dot2(tgt[i], tgt[j]) * inv_m + ctl[i] * ctl[j] * inv_i;
        }
    }
    ct_vec2 gacc = {
        prob->gravity.x + prob->ext_force.x * inv_m,
        prob->gravity.y + prob->ext_force.y * inv_m
    };
    float alpha = prob->ext_torque * inv_i;
    for (int i = 0; i < n; i++) {
        w.an[i] = dot2(nrm[i], gacc) + cn[i] * alpha;
        w.at[i] = dot2(tgt[i], gacc) + ctl[i] * alpha;
        w.fn[i] = 0.0f;
        w.ft[i] = 0.0f;
    }

    /* analytic feasibility pre-check when all contacts share one plane:
       the load line must fall inside the support interval and friction
       must be able to hold the tangential load */
    int common = n > 0;
    for (int i = 1; i < n && common; i++)
        if (fabsf(dot2(nrm[i], nrm[0]) - 1.0f) > 1e-3f) common = 0;
    if (common) {
        float h = dot2(r[0], nrm[0]);
        for (int i = 1; i < n && common; i++)
            if (fabsf(dot2(r[i], nrm[0]) - h) > 1e-3f) common = 0;
    }
    if (common) {
        ct_vec2 gtot = {
            prob->mass * prob->gravity.x + prob->ext_force.x,
            prob->mass * prob->gravity.y + prob->ext_force.y
        };
        float gn = dot2(gtot, nrm[0]);
        float gt = dot2(gtot, tgt[0]);
        float need_n = -gn;                 /* total normal force required */
        float mumin = w.mu[0];
        float smin = 0.0f, smax = 0.0f, h = dot2(r[0], nrm[0]);
        int imin = -1, imax = -1, nact = 0;
        for (int i = 0; i < n; i++) {
            if (!w.active[i]) continue;
            float s = dot2(r[i], tgt[0]);   /* com projects to s = 0 */
            if (imin < 0 || s < smin) { smin = s; imin = i; }
            if (imax < 0 || s > smax) { smax = s; imax = i; }
            if (w.mu[i] < mumin) mumin = w.mu[i];
            nact++;
        }
        float q = prob->ext_torque - h * gt; /* moment the normals must give */
        float eps = 1e-4f * (1.0f + fabsf(need_n));
        int bad = nact == 0
            || gn >= -1e-6f                             /* no support at all */
            || fabsf(gt) > mumin * need_n + eps         /* friction too weak */
            || q < need_n * smin - eps                  /* load line outside */
            || q > need_n * smax + eps;                 /* the support span  */
        if (bad) {
            res->error = CT_ERR_NO_SOLUTION;
            res->conflict_a = imin;
            res->conflict_b = imax;
            return;
        }
    }

    /* projected iteration on the complementarity problem, in three
       phases: normal forces first, then friction, then the coupled
       problem. the force distribution is indeterminate, and this order
       lands on the low-friction solution instead of the cone rim */
    float cr[CT_MAX_CONTACTS];
    int it = 0;
    while (it < prob->max_iterations) {
        float d = normal_sweep(&w);
        it++;
        if (d <= prob->tolerance) break;
    }
    while (it < prob->max_iterations) {
        float d = friction_sweep(&w);
        it++;
        if (d <= prob->tolerance) break;
    }
    res->residual = compute_residual(&w, cr);
    while (res->residual > prob->tolerance && it < prob->max_iterations) {
        for (int i = 0; i < n; i++)
            if (w.active[i]) block_update(&w, i);
        res->residual = compute_residual(&w, cr);
        it++;
    }
    res->iterations = it;
    res->converged = res->residual <= prob->tolerance;

    /* equilibrium residual: contact forces must cancel the applied load */
    ct_vec2 gtot = {
        prob->mass * prob->gravity.x + prob->ext_force.x,
        prob->mass * prob->gravity.y + prob->ext_force.y
    };
    float sfx = gtot.x, sfy = gtot.y, sm = prob->ext_torque;
    for (int i = 0; i < n; i++) {
        ct_vec2 f = {
            nrm[i].x * w.fn[i] + tgt[i].x * w.ft[i],
            nrm[i].y * w.fn[i] + tgt[i].y * w.ft[i]
        };
        sfx += f.x; sfy += f.y;
        sm += cross2(r[i], f);
    }
    res->eq_residual = fmaxf(fmaxf(fabsf(sfx), fabsf(sfy)), fabsf(sm));

    for (int i = 0; i < n; i++) {
        res->fn[i] = w.fn[i];
        res->ft[i] = w.ft[i];
    }

    if (!res->converged) {
        res->error = CT_ERR_NOT_CONVERGED;
        pick_conflict(cr, w.active, n, &res->conflict_a, &res->conflict_b);
    } else if (res->eq_residual > CT_EQ_TOL) {
        /* complementarity solved but the body still is not balanced:
           the contact set cannot hold the load (degenerate normals) */
        res->converged = 0;
        res->error = CT_ERR_NO_SOLUTION;
        pick_conflict(cr, w.active, n, &res->conflict_a, &res->conflict_b);
    }
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

/* uniform rod lying flat on two contact points */
static ct_problem flat_rod(void) {
    ct_problem p;
    memset(&p, 0, sizeof(p));
    p.mass = 2.0f;
    p.inertia = 0.5f;
    p.com.x = 0.0f; p.com.y = 0.5f;
    p.gravity.x = 0.0f; p.gravity.y = -9.8f;
    p.max_iterations = 128;
    p.tolerance = 1e-6f;
    p.n_contacts = 2;
    p.contacts[0].point.x = -0.5f; p.contacts[0].point.y = 0.0f;
    p.contacts[0].normal.x = 0.0f; p.contacts[0].normal.y = 1.0f;
    p.contacts[0].mu = 0.5f; p.contacts[0].gap = 0.0f;
    p.contacts[1].point.x = 0.5f; p.contacts[1].point.y = 0.0f;
    p.contacts[1].normal.x = 0.0f; p.contacts[1].normal.y = 1.0f;
    p.contacts[1].mu = 0.5f; p.contacts[1].gap = 0.0f;
    return p;
}

int main(void) {
    /* flat rest: converges, residual reported, zero penetration */
    ct_problem p = flat_rod();
    ct_result r;
    ct_solve(&p, &r);
    printf("residual=%.3g iterations=%d eq_residual=%.3g\n",
           (double)r.residual, r.iterations, (double)r.eq_residual);
    check("flat_rest_converges", r.error == CT_OK && r.converged);
    check("residual_reported", r.residual <= p.tolerance);
    check("zero_penetration", r.penetration == 0.0f);
    check("equilibrium_small", r.eq_residual < 1e-3f);

    /* analytic statics, com shifted off center: the normal forces are
       determinate, fn = mg * (0.4, 0.6). reduced equilibrium system
       A [fn0;fn1] = [mg; 0] with A = [[1,1],[-0.6,0.4]], |A^-1|_inf = 1.6,
       so |f - f_exact|_inf <= 1.6 * |e|_inf where e is the residual of
       the solved forces on those two equilibrium rows (the friction
       torque -0.5*(ft0+ft1) enters the moment row's right-hand side) */
    ct_problem pa = flat_rod();
    pa.com.x = 0.1f;
    ct_result ra;
    ct_solve(&pa, &ra);
    float mg = 2.0f * 9.8f;
    float e_y = ra.fn[0] + ra.fn[1] - mg;
    float e_m = -0.6f * ra.fn[0] + 0.4f * ra.fn[1]
              - 0.5f * (ra.ft[0] + ra.ft[1]);
    float eres = fmaxf(fabsf(e_y), fabsf(e_m));
    float bound = 1.6f * eres;
    float err = fmaxf(fabsf(ra.fn[0] - mg * 0.4f),
                      fabsf(ra.fn[1] - mg * 0.6f));
    printf("analytic_err=%.3g bound=%.3g\n", (double)err, (double)bound);
    check("matches_analytic", err < 1e-3f);
    check("error_bound_valid", err <= bound + 1e-5f);
    check("no_friction_at_rest", fabsf(r.ft[0]) < 1e-3f && fabsf(r.ft[1]) < 1e-3f);

    /* sideways load: friction must balance it and stay inside the cone.
       the tangent convention is t = (-n.y, n.x) = (-1,0) here, so the
       friction forces must sum to +4 to cancel the +4 x load */
    ct_problem pf = flat_rod();
    pf.ext_force.x = 4.0f;
    ct_result rf;
    ct_solve(&pf, &rf);
    check("friction_balances", rf.error == CT_OK
        && fabsf(rf.ft[0] + rf.ft[1] - 4.0f) < 1e-2f);
    check("inside_cone", fabsf(rf.ft[0]) <= 0.5f * rf.fn[0] + 1e-3f
        && fabsf(rf.ft[1]) <= 0.5f * rf.fn[1] + 1e-3f);

    /* an open contact (gap > 0) must stay unloaded */
    ct_problem ph = flat_rod();
    ph.n_contacts = 3;
    ph.contacts[2].point.x = 0.0f; ph.contacts[2].point.y = 0.2f;
    ph.contacts[2].normal.x = 0.0f; ph.contacts[2].normal.y = 1.0f;
    ph.contacts[2].mu = 0.5f; ph.contacts[2].gap = 0.2f;
    ct_result rh;
    ct_solve(&ph, &rh);
    check("hovering_unloaded", rh.error == CT_OK && rh.fn[2] == 0.0f);

    /* zero normal is rejected and named */
    ct_problem pz = flat_rod();
    pz.contacts[1].normal.x = 0.0f; pz.contacts[1].normal.y = 0.0f;
    ct_result rz;
    ct_solve(&pz, &rz);
    check("zero_normal_rejected", rz.error == CT_ERR_ZERO_NORMAL
        && rz.error_contact == 1);

    /* penetrating contact is rejected: penetration must be zero */
    ct_problem pp = flat_rod();
    pp.contacts[1].gap = -0.01f;
    ct_result rp;
    ct_solve(&pp, &rp);
    check("penetration_rejected", rp.error == CT_ERR_PENETRATION
        && rp.error_contact == 1 && rp.penetration > 0.0f);

    /* com outside the support span: no static solution, conflict pair
       must be the two contacts bounding the span */
    ct_problem po = flat_rod();
    po.com.x = 2.0f;
    ct_result ro;
    ct_solve(&po, &ro);
    check("no_solution_detected", ro.error == CT_ERR_NO_SOLUTION);
    check("conflict_pair_reported", ro.conflict_a >= 0 && ro.conflict_b >= 0
        && ro.conflict_a != ro.conflict_b);

    /* sideways normals cannot support weight: must fail with a named
       conflict pair, one way or the other */
    ct_problem px = flat_rod();
    px.contacts[0].normal.x = 1.0f; px.contacts[0].normal.y = 0.0f;
    px.contacts[1].normal.x = -1.0f; px.contacts[1].normal.y = 0.0f;
    ct_result rx;
    ct_solve(&px, &rx);
    check("unsupported_load_fails", rx.error == CT_ERR_NO_SOLUTION
        || rx.error == CT_ERR_NOT_CONVERGED);
    check("unsupported_conflict_named", rx.conflict_a >= 0
        && rx.conflict_b >= 0 && rx.conflict_a != rx.conflict_b);

    /* iteration budget exhausted: not converged, conflict pair named */
    ct_problem ps = flat_rod();
    ps.max_iterations = 0;
    ct_result rs;
    ct_solve(&ps, &rs);
    check("not_converged_reported", rs.error == CT_ERR_NOT_CONVERGED
        && !rs.converged);
    check("stall_conflict_named", rs.conflict_a >= 0 && rs.conflict_b >= 0
        && rs.conflict_a != rs.conflict_b);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
