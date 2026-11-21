
/* hysteresis: rate-independent material hysteresis for a stretchy rope.
   Loading follows a piecewise-linear backbone (stiffness varies with the
   stretch ratio). Unloading branches off with the initial tangent and
   rejoins a residual-strain shifted copy of the backbone; the residual
   strain is solved so the loop area equals a configurable fraction of the
   work done during loading. Load/unload samples are recorded on two curves
   and can be compared against a hand tensile test, with an upper bound on
   the return-point error. Invalid parameters and non-positive stiffness
   are rejected with error codes. */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define MAX_KEYS   16
#define MAX_BRANCH 32
#define MAX_SAMP   128

enum {
    HY_OK = 0,
    HY_ERR_TOO_FEW     = -1,
    HY_ERR_ORDER       = -2,
    HY_ERR_ORIGIN      = -3,
    HY_ERR_STIFFNESS   = -4,
    HY_ERR_DISSIPATION = -5,
    HY_ERR_STRAIN      = -6,
    HY_ERR_CAPACITY    = -7,
    HY_ERR_FINITE      = -8
};

typedef struct { float eps; float sig; } HyPoint;

typedef struct {
    HyPoint key[MAX_KEYS];      /* loading backbone, strain ascending */
    int n;
    float e0;                   /* initial tangent stiffness */
    float dissip;               /* target loop area fraction d in [0,1) */
    float eps_max;              /* largest strain reached in this loop */
    float eps_prev;
    float turn;                 /* turn strain of the active branch */
    float resid;                /* residual strain where unload hits zero */
    HyPoint br[MAX_BRANCH];     /* active unload/reload branch, descending */
    int bn;
    int phase;                  /* 0 = loading, 1 = unloading */
    int has_prev;
} HMat;

static const char* hy_err_name(int rc) {
    switch (rc) {
        case HY_OK:                return "ok";
        case HY_ERR_TOO_FEW:       return "too_few_keys";
        case HY_ERR_ORDER:         return "strain_not_increasing";
        case HY_ERR_ORIGIN:        return "first_key_not_origin";
        case HY_ERR_STIFFNESS:     return "non_positive_stiffness";
        case HY_ERR_DISSIPATION:   return "dissipation_out_of_range";
        case HY_ERR_STRAIN:        return "negative_or_nan_strain";
        case HY_ERR_CAPACITY:      return "capacity_exceeded";
        case HY_ERR_FINITE:        return "non_finite_value";
        default:                   return "unknown";
    }
}

static int hy_init(HMat* m, const HyPoint* keys, int n, float dissip) {
    if (!isfinite(dissip) || dissip < 0.0f || dissip >= 1.0f)
        return HY_ERR_DISSIPATION;
    if (n < 2 || n > MAX_KEYS) return HY_ERR_TOO_FEW;
    for (int i = 0; i < n; i++)
        if (!isfinite(keys[i].eps) || !isfinite(keys[i].sig))
            return HY_ERR_FINITE;
    if (keys[0].eps != 0.0f || keys[0].sig != 0.0f)
        return HY_ERR_ORIGIN;
    for (int i = 1; i < n; i++) {
        if (!(keys[i].eps > keys[i - 1].eps)) return HY_ERR_ORDER;
        float slope = (keys[i].sig - keys[i - 1].sig) /
                      (keys[i].eps - keys[i - 1].eps);
        if (!isfinite(slope) || !(slope > 0.0f))
            return HY_ERR_STIFFNESS;   /* zero or negative stiffness */
    }
    memset(m, 0, sizeof(*m));
    memcpy(m->key, keys, (size_t)n * sizeof(HyPoint));
    m->n = n;
    m->e0 = (keys[1].sig - keys[0].sig) / (keys[1].eps - keys[0].eps);
    m->dissip = dissip;
    return HY_OK;
}

/* stress on the loading backbone; extrapolates with the last slope */
static float hy_backbone(const HMat* m, float eps) {
    if (eps <= 0.0f) return 0.0f;
    if (eps >= m->key[m->n - 1].eps) {
        const HyPoint* a = &m->key[m->n - 2];
        const HyPoint* b = &m->key[m->n - 1];
        float k = (b->sig - a->sig) / (b->eps - a->eps);
        return b->sig + k * (eps - b->eps);
    }
    for (int i = 1; i < m->n; i++) {
        if (eps <= m->key[i].eps) {
            const HyPoint* a = &m->key[i - 1];
            const HyPoint* b = &m->key[i];
            float u = (eps - a->eps) / (b->eps - a->eps);
            return a->sig + u * (b->sig - a->sig);
        }
    }
    return 0.0f;
}

/* work stored under the backbone up to strain eps (trapezoidal) */
static float hy_work(const HMat* m, float eps) {
    if (eps <= 0.0f) return 0.0f;
    float w = 0.0f;
    int i = 1;
    for (; i < m->n && m->key[i].eps <= eps; i++)
        w += 0.5f * (m->key[i].sig + m->key[i - 1].sig) *
             (m->key[i].eps - m->key[i - 1].eps);
    if (eps > m->key[m->n - 1].eps) {
        w += 0.5f * (m->key[m->n - 1].sig + hy_backbone(m, eps)) *
             (eps - m->key[m->n - 1].eps);
    } else if (i < m->n) {
        w += 0.5f * (m->key[i - 1].sig + hy_backbone(m, eps)) *
             (eps - m->key[i - 1].eps);
    }
    return w;
}

/* Unloading is an affine re-mapping of the loading backbone onto the
   interval [er, et] (a scaled Masing rule):
     sig_un(eps) = F( (eps - er) * et / (et - er) )
   It leaves the turn point (et, F(et)) continuously, hits zero stress at
   the residual strain er, and its enclosed loop area works out exactly to
     W_load - W_unload = d * W(et),  er = d * et
   so the configured dissipation fraction is met in closed form with no
   iteration. Evaluating through the backbone (rather than a chord through
   branch keypoints) keeps the mapped curve and its loop area exact. */
static float hy_eval_branch(const HMat* m, float eps) {
    if (m->bn == 0) return 0.0f;
    if (m->resid <= 0.0f) {                  /* d = 0: perfectly elastic */
        if (eps >= m->turn) return hy_backbone(m, m->turn);
        return hy_backbone(m, eps);
    }
    if (m->turn <= m->resid) return 0.0f;
    if (eps >= m->turn) return hy_backbone(m, m->turn);
    if (eps <= m->resid) return 0.0f;
    float u = (eps - m->resid) * m->turn / (m->turn - m->resid);
    return hy_backbone(m, u);
}

/* build the unload branch at a loading reversal. The residual strain
   follows directly from the dissipation fraction: er = d * et */
static int hy_build_branch(HMat* m, float et) {
    float sigt = hy_backbone(m, et);
    float er = m->dissip * et;
    m->turn = et;
    m->resid = er;
    m->phase = 1;
    m->bn = 0;
    if (sigt <= 0.0f || et <= er) {          /* reversal at rest: zero branch */
        m->br[m->bn++] = (HyPoint){0.0f, 0.0f};
        return HY_OK;
    }
    float span = et - er;
    /* image of every backbone keypoint under the affine map, descending */
    m->br[m->bn++] = (HyPoint){et, sigt};
    for (int i = m->n - 2; i >= 1; i--) {
        float u = m->key[i].eps;
        if (u >= et) continue;
        if (m->bn >= MAX_BRANCH) return HY_ERR_CAPACITY;
        float e = er + span * (u / et);
        m->br[m->bn++] = (HyPoint){e, hy_backbone(m, u)};
    }
    if (m->bn >= MAX_BRANCH) return HY_ERR_CAPACITY;
    m->br[m->bn++] = (HyPoint){er, 0.0f};
    return HY_OK;
}

/* constitutive query: tracks load/unload direction and remembers the loop */
static int hy_stress(HMat* m, float eps, float* sig) {
    if (!isfinite(eps) || eps < 0.0f) return HY_ERR_STRAIN;
    if (!m->has_prev) {
        m->has_prev = 1;
        m->eps_prev = eps;
        m->eps_max = eps;
        *sig = hy_backbone(m, eps);
        return HY_OK;
    }
    int rc = HY_OK;
    if (eps > m->eps_prev) {                                  /* loading */
        if (m->phase == 1 && eps < m->turn) {
            *sig = hy_eval_branch(m, eps);                    /* reloading */
        } else {
            if (m->phase == 1) { m->phase = 0; m->bn = 0; }   /* rejoined */
            *sig = hy_backbone(m, eps);
        }
    } else if (eps < m->eps_prev) {                           /* unloading */
        if (m->phase != 1 || m->eps_max > m->turn) {
            rc = hy_build_branch(m, m->eps_max);
            if (rc != HY_OK) return rc;
        }
        *sig = hy_eval_branch(m, eps);
    } else {
        *sig = m->phase == 1 ? hy_eval_branch(m, eps) : hy_backbone(m, eps);
    }
    m->eps_prev = eps;
    if (eps > m->eps_max) m->eps_max = eps;
    return rc;
}

/* ---- recorded load/unload curves ---- */

typedef struct {
    HyPoint load[MAX_SAMP]; int nl;
    HyPoint unload[MAX_SAMP]; int nu;
} Recorder;

static void rec_reset(Recorder* r) { memset(r, 0, sizeof(*r)); }

static int rec_step(Recorder* r, const HMat* m, float eps, float sig) {
    if (m->phase == 0) {
        if (r->nl >= MAX_SAMP) return HY_ERR_CAPACITY;
        r->load[r->nl++] = (HyPoint){eps, sig};
    } else {
        if (r->nu >= MAX_SAMP) return HY_ERR_CAPACITY;
        r->unload[r->nu++] = (HyPoint){eps, sig};
    }
    return HY_OK;
}

static float curve_area(const HyPoint* p, int n) {
    float a = 0.0f;
    for (int i = 1; i < n; i++)
        a += 0.5f * (p[i].sig + p[i - 1].sig) *
             fabsf(p[i].eps - p[i - 1].eps);
    return a;
}

static int run_cycle(HMat* m, Recorder* r, float emax, int steps) {
    rec_reset(r);
    for (int i = 0; i <= steps; i++) {      /* loading including the turn */
        float eps = emax * (float)i / (float)steps;
        float sig;
        int rc = hy_stress(m, eps, &sig);
        if (rc != HY_OK) return rc;
        rc = rec_step(r, m, eps, sig);
        if (rc != HY_OK) return rc;
    }
    /* open the unload branch at the turn; the turn is recorded again as
       the first point of the unload curve so both legs span [0, et] */
    int rc = hy_build_branch(m, emax);
    if (rc != HY_OK) return rc;
    {
        float tsig = hy_eval_branch(m, emax);
        rc = rec_step(r, m, emax, tsig);
        if (rc != HY_OK) return rc;
    }
    for (int i = 1; i <= steps; i++) {
        float eps = emax * (1.0f - (float)i / (float)steps);
        float sig = hy_eval_branch(m, eps);
        rc = rec_step(r, m, eps, sig);
        if (rc != HY_OK) return rc;
    }
    return HY_OK;
}

/* ---- comparison with a hand tensile test ---- */

typedef struct { float eps; float sig; int phase; } MeasPoint; /* 0 load, 1 unload */

typedef struct {
    float max_abs;       /* worst stress deviation, load and unload */
    float max_rel;
    float resid_bound;   /* upper bound on return-point strain error */
} CompareReport;

static float hy_min_backbone_slope(const HMat* m, float upto) {
    float kmin = 1e30f;
    for (int i = 1; i < m->n && m->key[i - 1].eps < upto; i++) {
        float k = (m->key[i].sig - m->key[i - 1].sig) /
                  (m->key[i].eps - m->key[i - 1].eps);
        if (k < kmin) kmin = k;
    }
    return kmin;
}

/* the branch is the backbone re-parameterized by u = et*(eps-er)/(et-er),
   so its slope is the backbone slope times the constant map ratio */
static float hy_min_unload_slope(const HMat* m) {
    if (m->bn == 0 || m->turn <= m->resid)
        return hy_min_backbone_slope(m, m->turn);
    float ratio = m->turn / (m->turn - m->resid);
    return hy_min_backbone_slope(m, m->turn) * ratio;
}

static int hy_compare(const HMat* m, const MeasPoint* s, int n,
                      CompareReport* rep) {
    float max_abs = 0.0f, max_rel = 0.0f, max_un = 0.0f;
    for (int i = 0; i < n; i++) {
        if (!isfinite(s[i].eps) || !isfinite(s[i].sig) || s[i].eps < 0.0f)
            return HY_ERR_FINITE;
        float model = s[i].phase == 0 ? hy_backbone(m, s[i].eps)
                                      : hy_eval_branch(m, s[i].eps);
        float dev = fabsf(model - s[i].sig);
        if (dev > max_abs) max_abs = dev;
        float rel = dev / (fabsf(model) + 1e-6f);
        if (rel > max_rel) max_rel = rel;
        if (s[i].phase == 1 && dev > max_un) max_un = dev;
    }
    /* a stress deviation of ds shifts the zero crossing by at most ds/k */
    float kmin = hy_min_unload_slope(m);
    rep->max_abs = max_abs;
    rep->max_rel = max_rel;
    rep->resid_bound = (max_un > 0.0f ? max_un : max_abs) / kmin;
    return HY_OK;
}

/* linear interpolation of the strain where the unload curve crosses zero */
static float zero_crossing(const HyPoint* p, int n) {
    for (int i = 1; i < n; i++) {
        if (p[i - 1].sig > 0.0f && p[i].sig <= 0.0f) {
            float u = p[i - 1].sig / (p[i - 1].sig - p[i].sig);
            return p[i - 1].eps + u * (p[i].eps - p[i - 1].eps);
        }
    }
    return 0.0f;
}

/* ---- minimal verlet rope driver: anchor at 0, one free end ---- */

typedef struct {
    HMat* mat;
    float rest;
    float x, px;
    float inv_mass;
    float damping;
} Rope;

static void rope_init(Rope* r, HMat* m, float rest, float mass) {
    r->mat = m;
    r->rest = rest;
    r->x = r->px = rest;
    r->inv_mass = 1.0f / mass;
    r->damping = 0.6f;
}

static void rope_stretch(Rope* r, float length) {
    r->x = r->px = length;   /* grab and hold: zero velocity */
}

static void rope_step(Rope* r, float dt) {
    float v = (r->x - r->px) * r->damping;
    r->px = r->x;
    r->x += v;
    float eps = (r->x - r->rest) / r->rest;
    float sig = 0.0f;
    if (eps > 0.0f) hy_stress(r->mat, eps, &sig);
    else { hy_stress(r->mat, 0.0f, &sig); sig = 0.0f; }  /* slack: no tension */
    r->x += -sig * r->inv_mass * dt * dt;
}

/* ---- self tests ---- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    /* rope-like stiffening backbone: stiffness grows with stretch ratio */
    HyPoint keys[] = {
        {0.0f, 0.0f}, {0.2f, 0.3f}, {0.4f, 0.8f}, {0.6f, 1.5f}
    };
    int nk = (int)(sizeof(keys) / sizeof(keys[0]));

    HMat m;
    check("init_ok", hy_init(&m, keys, nk, 0.3f) == HY_OK);
    check("stiffness_by_strain",
          fabsf(hy_backbone(&m, 0.1f) - 0.15f) < 1e-5f &&   /* k = 1.5 */
          fabsf(hy_backbone(&m, 0.5f) - 1.15f) < 1e-5f &&   /* k = 3.5 */
          fabsf(hy_backbone(&m, 0.8f) - 2.2f)  < 1e-5f);    /* extrapolate */

    /* parameter validation: out-of-range dissipation and bad stiffness */
    HMat bad;
    int rc = hy_init(&bad, keys, nk, -0.1f);
    printf("reject_negative_d=%s\n", hy_err_name(rc));
    check("reject_negative_d", rc == HY_ERR_DISSIPATION);
    rc = hy_init(&bad, keys, nk, 1.0f);
    check("reject_full_dissipation", rc == HY_ERR_DISSIPATION);

    HyPoint flat[] = {{0.0f, 0.0f}, {0.2f, 0.3f}, {0.4f, 0.3f}};
    rc = hy_init(&bad, flat, 3, 0.3f);
    printf("reject_zero_stiffness=%s\n", hy_err_name(rc));
    check("reject_zero_stiffness", rc == HY_ERR_STIFFNESS);

    HyPoint neg[] = {{0.0f, 0.0f}, {0.2f, 0.3f}, {0.4f, 0.1f}};
    rc = hy_init(&bad, neg, 3, 0.3f);
    printf("reject_negative_stiffness=%s\n", hy_err_name(rc));
    check("reject_negative_stiffness", rc == HY_ERR_STIFFNESS);

    HyPoint unsorted[] = {{0.0f, 0.0f}, {0.4f, 0.8f}, {0.2f, 0.3f}};
    check("reject_unsorted", hy_init(&bad, unsorted, 3, 0.3f) == HY_ERR_ORDER);
    HyPoint origin[] = {{0.1f, 0.0f}, {0.2f, 0.3f}};
    check("reject_bad_origin", hy_init(&bad, origin, 2, 0.3f) == HY_ERR_ORIGIN);
    check("reject_too_few", hy_init(&bad, keys, 1, 0.3f) == HY_ERR_TOO_FEW);
    HyPoint nanpt[] = {{0.0f, 0.0f}, {0.2f, NAN}};
    check("reject_non_finite", hy_init(&bad, nanpt, 2, 0.3f) == HY_ERR_FINITE);

    HMat em;
    check("elastic_init", hy_init(&em, keys, nk, 0.0f) == HY_OK);
    float sig;
    check("reject_negative_strain", hy_stress(&em, -0.01f, &sig) == HY_ERR_STRAIN);

    /* elastic material: zero loop area, full return to rest */
    HMat em2;
    hy_init(&em2, keys, nk, 0.0f);
    Recorder rec;
    check("elastic_cycle", run_cycle(&em2, &rec, 0.5f, 24) == HY_OK);
    float e_area = curve_area(rec.load, rec.nl) - curve_area(rec.unload, rec.nu);
    check("elastic_zero_area", fabsf(e_area) < 1e-5f);
    check("elastic_full_return", fabsf(rec.unload[rec.nu - 1].sig) < 1e-6f &&
                                 fabsf(rec.unload[rec.nu - 1].eps) < 1e-6f);

    /* hysteretic material: recorded loop area matches the configured target */
    HMat hm;
    check("hyst_init", hy_init(&hm, keys, nk, 0.3f) == HY_OK);
    check("hyst_cycle", run_cycle(&hm, &rec, 0.5f, 24) == HY_OK);
    float w_load = curve_area(rec.load, rec.nl);
    float w_un = curve_area(rec.unload, rec.nu);
    float loop = w_load - w_un;
    float target = 0.3f * hy_work(&hm, 0.5f);
    printf("loop_area=%.5f target=%.5f work=%.5f residual=%.5f\n",
           loop, target, hy_work(&hm, 0.5f), hm.resid);

    /* exact identity of the affine model. The branch maps [er, et] onto
       [0, et], so its work is W(et) scaled by the interval width ratio,
       and the loop area is exactly d * W(et) */
    float w_analytic = hy_work(&hm, hm.turn) *
                       (hm.turn - hm.resid) / hm.turn;
    check("unload_work_identity",
          fabsf(w_analytic - (1.0f - hm.dissip) * hy_work(&hm, hm.turn)) < 1e-6f);
    check("loop_area_exact",
          fabsf(hm.dissip * hy_work(&hm, hm.turn) - target) < 1e-6f);

    /* densely sampled cycle: trapezoidal integration then reproduces the
       configured loop area to within integration error */
    Recorder fine;
    HMat hf;
    hy_init(&hf, keys, nk, 0.3f);
    check("fine_cycle", run_cycle(&hf, &fine, 0.5f, 100) == HY_OK);
    float fine_loop = curve_area(fine.load, fine.nl) -
                      curve_area(fine.unload, fine.nu);
    printf("fine_loop_area=%.5f target=%.5f\n", fine_loop, target);
    check("loop_area_target", fabsf(fine_loop - target) < 0.06f * target);
    check("area_positive", loop > 0.0f);

    /* unloading is softer than loading at the same strain, and leaves a
       residual strain: the rope does not spring all the way back */
    float s_load = hy_backbone(&hm, 0.25f);
    float s_un = hy_eval_branch(&hm, 0.25f);
    check("unload_softer", s_un < s_load && s_un > 0.0f);
    check("residual_strain", hm.resid > 0.05f &&
          fabsf(hy_eval_branch(&hm, hm.resid)) < 1e-5f &&
          fabsf(hy_eval_branch(&hm, 0.0f)) < 1e-5f);

    /* the unload branch leaves the turn point with matching stress; inside
       the branch its slope is the backbone slope scaled by the map ratio */
    float h = 1e-3f;
    float e_in = 0.35f;
    float k_branch = (hy_eval_branch(&hm, e_in + h) -
                      hy_eval_branch(&hm, e_in - h)) / (2 * h);
    float ratio = 0.5f / (0.5f - hm.resid);
    float u_in = (e_in - hm.resid) * 0.5f / (0.5f - hm.resid);
    float k_back = (hy_backbone(&hm, u_in + h) -
                    hy_backbone(&hm, u_in - h)) / (2 * h);
    check("unload_turn_continuous",
          fabsf(hy_eval_branch(&hm, 0.5f) - hy_backbone(&hm, 0.5f)) < 1e-5f &&
          fabsf(k_branch - ratio * k_back) < 0.02f * k_branch);
    /* residual strain equals d * turn strain in closed form */
    check("residual_formula", fabsf(hm.resid - 0.3f * 0.5f) < 1e-6f);

    /* reload memory: strain returns along the same branch and rejoins at
       the turn point without a stress jump */
    float sig_at_turn;
    hy_stress(&hm, 0.5f, &sig_at_turn);      /* reload back to the turn */
    check("reload_rejoins", fabsf(sig_at_turn - hy_backbone(&hm, 0.5f)) < 1e-5f);
    hy_stress(&hm, 0.3f, &sig);              /* unload again: same branch */
    check("branch_remembered", fabsf(sig - hy_eval_branch(&hm, 0.3f)) < 1e-5f &&
                               fabsf(hm.turn - 0.5f) < 1e-6f);

    /* hand tensile test: shift every measured unload stress down by a
       constant offset so its zero crossing lands above the model residual;
       the reported error bound must cover the crossing displacement */
    MeasPoint samp[2 * MAX_SAMP];
    int ns = 0;
    const float delta = 0.02f;
    for (int i = 0; i < rec.nl && ns < 2 * MAX_SAMP; i++)
        samp[ns++] = (MeasPoint){rec.load[i].eps, rec.load[i].sig, 0};
    for (int i = 0; i < rec.nu && ns < 2 * MAX_SAMP; i++) {
        float s = rec.unload[i].sig - delta;
        samp[ns++] = (MeasPoint){rec.unload[i].eps, s < 0.0f ? 0.0f : s, 1};
    }
    CompareReport rep;
    memset(&rep, 0, sizeof(rep));
    check("compare_ok", hy_compare(&hm, samp, ns, &rep) == HY_OK);
    printf("compare_max_abs=%.5f resid_bound=%.5f\n",
           rep.max_abs, rep.resid_bound);
    check("compare_finds_deviation", fabsf(rep.max_abs - delta) < 1e-6f);

    HyPoint meas_un[MAX_SAMP];
    for (int i = 0; i < rec.nu; i++) {
        float s = rec.unload[i].sig - delta;
        meas_un[i] = (HyPoint){rec.unload[i].eps, s < 0.0f ? 0.0f : s};
    }
    float meas_resid = zero_crossing(meas_un, rec.nu);
    float resid_err = fabsf(meas_resid - hm.resid);
    printf("measured_residual=%.5f model_residual=%.5f err=%.5f\n",
           meas_resid, hm.resid, resid_err);
    check("residual_bound_holds", resid_err <= rep.resid_bound + 1e-6f);
    check("bound_is_finite", isfinite(rep.resid_bound) && rep.resid_bound > 0.0f);

    /* recorder capacity: each leg takes steps+1 slots; the load leg fills
       the buffer first, then the unload turn point overflows it. A fresh
       material keeps this test independent of the loop state above. */
    HMat cm;
    hy_init(&cm, keys, nk, 0.3f);
    Recorder full;
    rc = run_cycle(&cm, &full, 0.5f, MAX_SAMP);
    check("recorder_capacity", rc == HY_ERR_CAPACITY &&
                               full.nl == MAX_SAMP && full.nu == 0);

    /* end-to-end: stretch a verlet rope, let go, measure where it settles.
       fresh materials here — the ones above carry loop state */
    HMat rm, re;
    hy_init(&rm, keys, nk, 0.3f);
    hy_init(&re, keys, nk, 0.0f);
    Rope rope, erope;
    rope_init(&rope, &rm, 1.0f, 1.0f);
    rope_init(&erope, &re, 1.0f, 1.0f);
    rope_stretch(&rope, 1.5f);    /* stretch ratio 1.5, eps = 0.5 */
    rope_stretch(&erope, 1.5f);
    for (int i = 0; i < 30000; i++) { rope_step(&rope, 0.01f);
                                      rope_step(&erope, 0.01f); }
    float rope_eps = (rope.x - rope.rest) / rope.rest;
    float erope_eps = (erope.x - erope.rest) / erope.rest;
    printf("settled_hysteretic_eps=%.5f (residual %.5f)\n",
           rope_eps, rm.resid);
    printf("settled_elastic_eps=%.5f\n", erope_eps);
    check("rope_keeps_residual", fabsf(rope_eps - rm.resid) < 0.02f &&
                                 rope_eps > 0.05f);
    /* the elastic rope comes back far closer to rest than the hysteretic
       one and keeps shrinking, i.e. there is no residual-strain floor */
    check("rope_elastic_returns", fabsf(erope_eps) < 0.02f &&
                                  erope_eps < rope_eps - 0.05f);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
