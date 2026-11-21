
/* curve: keyframe sampling with cubic hermite tangents and clamping, plus a
   C2 periodic (cyclic) cubic spline used to continue loops outside the keys.
   The periodic spline solves for one second derivative per key under the
   cyclic constraint, so a loop's value, slope and curvature are continuous
   across the seam; extrapolation reduces t modulo the period instead of
   continuing the last cubic, which is what removes the runaway acceleration. */
#include <stdio.h>
#include <string.h>
#include <math.h>

#define MAX_KEYS 16

typedef struct { float t; float v; float tan_in; float tan_out; } Key;

typedef struct {
    Key keys[MAX_KEYS];
    int count;
} Curve;

static void curve_init(Curve* c) { memset(c, 0, sizeof(*c)); }

static int curve_add(Curve* c, float t, float v, float tin, float tout) {
    if (c->count >= MAX_KEYS) return -1;
    if (c->count > 0 && t <= c->keys[c->count - 1].t) return -1;  /* strictly increasing */
    Key* k = &c->keys[c->count++];
    k->t = t; k->v = v; k->tan_in = tin; k->tan_out = tout;
    return c->count - 1;
}

static float curve_eval(const Curve* c, float t) {
    if (c->count == 0) return 0.0f;
    if (t <= c->keys[0].t) return c->keys[0].v;
    if (t >= c->keys[c->count - 1].t) return c->keys[c->count - 1].v;
    for (int i = 1; i < c->count; i++) {
        const Key* a = &c->keys[i - 1];
        const Key* b = &c->keys[i];
        if (t <= b->t) {
            float dt = b->t - a->t;
            float u = (t - a->t) / dt;
            float u2 = u * u, u3 = u2 * u;
            float h00 = 2 * u3 - 3 * u2 + 1;
            float h10 = u3 - 2 * u2 + u;
            float h01 = -2 * u3 + 3 * u2;
            float h11 = u3 - u2;
            return h00 * a->v + h10 * dt * a->tan_out + h01 * b->v + h11 * dt * b->tan_in;
        }
    }
    return c->keys[c->count - 1].v;
}

/* derivative for continuity checks */
static float curve_deriv(const Curve* c, float t) {
    float h = 1e-3f;
    float a = curve_eval(c, t - h);
    float b = curve_eval(c, t + h);
    return (b - a) / (2 * h);
}

/* ------------------------------------------------------------------ *
 * Periodic C2 cubic spline.
 *
 * The keys describe exactly one period [t0, t0 + period).  Key i is followed
 * by key i+1 for i < n-1, and the last key wraps to key 0 over a closing
 * segment of length (t0 + period) - t[n-1].  We solve, under that cyclic
 * coupling, for the second derivative M[i] at each key.  Because there is a
 * single M per key, the cubics on either side of a key share value, slope and
 * curvature there: C2 at every knot, including the wrap seam.
 *
 * Evaluation always reduces t into the base period, so the value at
 * t + k*period equals the value at t for every integer k.  There is no growing
 * terminal cubic and therefore no high-frequency blow-up past the last key.
 * ------------------------------------------------------------------ */

typedef enum {
    PC_OK = 0,
    PC_ERR_PERIOD = -1,   /* period is zero, negative or not finite        */
    PC_ERR_FEW_KEYS = -2, /* fewer than three distinct keys                 */
    PC_ERR_SPACING = -3,  /* non-finite times/values, non-increasing times,
                             or key span not strictly less than the period  */
    PC_ERR_CAPACITY = -4, /* more keys than MAX_PC_KEYS                     */
    PC_ERR_GAIN = -5      /* requested gain cap exceeded by the closure     */
} PCStatus;

#define MAX_PC_KEYS MAX_KEYS
#define PC_MAX_MODES (MAX_PC_KEYS / 2 + 1)

typedef struct {
    float t[MAX_PC_KEYS];
    float v[MAX_PC_KEYS];
    float m[MAX_PC_KEYS];   /* second derivative at each key */
    int n;
    float period;
} PeriodicCurve;

static const char* pc_status_str(PCStatus s) {
    switch (s) {
        case PC_OK:           return "ok";
        case PC_ERR_PERIOD:   return "period must be finite and > 0";
        case PC_ERR_FEW_KEYS: return "need at least 3 keys";
        case PC_ERR_SPACING:  return "key times must strictly increase and span < period";
        case PC_ERR_CAPACITY: return "too many keys";
        case PC_ERR_GAIN:     return "closure-segment gain exceeds requested cap";
    }
    return "unknown";
}

/* interval lengths; h[n-1] is the closing segment back to key 0 */
static void pc_gaps(const PeriodicCurve* pc, double h[MAX_PC_KEYS]) {
    int n = pc->n;
    for (int i = 0; i < n; i++) {
        double next = (i + 1 < n) ? (double)pc->t[i + 1]
                                  : (double)pc->t[0] + (double)pc->period;
        h[i] = next - (double)pc->t[i];
    }
}

/* Solve the cyclic system A M = r for the second derivatives, given values v.
   A is symmetric and strictly diagonally dominant (all h > 0), so Gaussian
   elimination without pivoting is stable; n <= MAX_PC_KEYS keeps it cheap. */
static PCStatus pc_solve(PeriodicCurve* pc, const float* v) {
    int n = pc->n;
    double A[MAX_PC_KEYS][MAX_PC_KEYS];
    double r[MAX_PC_KEYS], m[MAX_PC_KEYS], h[MAX_PC_KEYS];
    pc_gaps(pc, h);

    for (int i = 0; i < n; i++) {
        int im = (i + n - 1) % n;
        int ip = (i + 1) % n;
        for (int j = 0; j < n; j++) A[i][j] = 0.0;
        A[i][im] += h[im];
        A[i][i]  += 2.0 * (h[im] + h[i]);
        A[i][ip] += h[i];
        r[i] = 6.0 * (((double)v[ip] - v[i]) / h[i]
                    - ((double)v[i] - v[im]) / h[im]);
    }
    for (int col = 0; col < n; col++) {
        double piv = A[col][col];
        if (fabs(piv) < 1e-12) return PC_ERR_SPACING;
        for (int row = col + 1; row < n; row++) {
            double f = A[row][col] / piv;
            if (f == 0.0) continue;
            for (int j = col; j < n; j++) A[row][j] -= f * A[col][j];
            r[row] -= f * r[col];
        }
    }
    for (int i = n - 1; i >= 0; i--) {
        double s = r[i];
        for (int j = i + 1; j < n; j++) s -= A[i][j] * m[j];
        m[i] = s / A[i][i];
    }
    for (int i = 0; i < n; i++) pc->m[i] = (float)m[i];
    return PC_OK;
}

/* reduce t into [t0, t0 + period) */
static double pc_wrap(const PeriodicCurve* pc, double t) {
    double t0 = pc->t[0], P = pc->period;
    double x = fmod(t - t0, P);
    if (x < 0.0) x += P;
    if (x >= P) x = 0.0;               /* guard fmod rounding up to P */
    return t0 + x;
}

/* index i with t[i] <= x < t[i+1] (cyclically); x must already be wrapped */
static int pc_segment(const PeriodicCurve* pc, double x) {
    int n = pc->n, i = 0;
    while (i < n - 1 && x >= (double)pc->t[i + 1]) i++;
    return i;
}

/* Evaluate segment i on x in [t[i], t[i+1]] with no wrapping or searching.
   order: 0 = value, 1 = first derivative, 2 = second derivative.
   Used both for arithmetic (pc_eval) and for exact one-sided seam limits. */
static double pc_seg_eval(const PeriodicCurve* pc, int i, double x, int order) {
    int n = pc->n, j = (i + 1) % n;
    double ti = pc->t[i];
    double tj = (i + 1 < n) ? (double)pc->t[i + 1]
                            : (double)pc->t[0] + (double)pc->period;
    double hi = tj - ti;
    double a = tj - x;                 /* distance to the right knot */
    double b = x - ti;                 /* distance to the left knot  */
    double Mi = pc->m[i], Mj = pc->m[j];
    double vi = pc->v[i], vj = pc->v[j];
    if (order == 0)
        return Mi * a * a * a / (6.0 * hi) + Mj * b * b * b / (6.0 * hi)
             + (vi - Mi * hi * hi / 6.0) * a / hi
             + (vj - Mj * hi * hi / 6.0) * b / hi;
    if (order == 1)
        return -Mi * a * a / (2.0 * hi) + Mj * b * b / (2.0 * hi)
             + (vj - vi) / hi - hi * (Mj - Mi) / 6.0;
    return (Mi * a + Mj * b) / hi;     /* order == 2 */
}

/* Evaluate the wrapped spline at t. */
static float pc_eval(const PeriodicCurve* pc, float t, int order) {
    double x = pc_wrap(pc, (double)t);
    int i = pc_segment(pc, x);
    return (float)pc_seg_eval(pc, i, x, order);
}

/* Exact one-sided limits at every knot, including the wrap seam.  At knot r the
   left segment is (r-1) if r>0 else the closing segment n-1 evaluated at x=t0+P;
   the right segment is r if r<n else segment 0 evaluated at x=t0.  Returns the
   largest |left - right| for the requested order across all knots. */
static double pc_knot_jump(const PeriodicCurve* pc, int order) {
    int n = pc->n;
    double worst = 0.0;
    for (int r = 0; r < n; r++) {
        double left, right;
        if (r == 0) {
            left  = pc_seg_eval(pc, n - 1, (double)pc->t[0] + (double)pc->period, order);
            right = pc_seg_eval(pc, 0, (double)pc->t[0], order);
        } else {
            left  = pc_seg_eval(pc, r - 1, (double)pc->t[r], order);
            right = pc_seg_eval(pc, r, (double)pc->t[r], order);
        }
        double d = fabs(left - right);
        if (d > worst) worst = d;
    }
    return worst;
}

/* ------------------------------------------------------------------ *
 * Gain and gain curve.
 *
 * Mode k of the loop is v[i] = cos(2*pi*k*(t[i]-t0)/period), for k = 0..n/2,
 * the discrete frequencies the keys can represent.  Feeding that through the
 * spline and dividing the RMS of the reconstruction by the RMS of the ideal
 * continuous cosine gives the gain of mode k on a given window.  On the base
 * period and on any period translated from it the gain is ~1 (the loop is
 * reproduced); on the synthesised closing segment it can exceed 1, and that is
 * what the cap bounds.  Sampling farther out never changes the gain because
 * evaluation wraps.
 * ------------------------------------------------------------------ */

static double pc_mode_gain(const PeriodicCurve* pc, int k, double lo, double hi,
                           int samples) {
    int n = pc->n;
    double t0 = pc->t[0], P = pc->period;
    double w = 2.0 * 3.14159265358979323846 * (double)k / P;

    PeriodicCurve tmp = *pc;
    for (int i = 0; i < n; i++) tmp.v[i] = (float)cos(w * ((double)pc->t[i] - t0));
    if (pc_solve(&tmp, tmp.v) != PC_OK) return 1e30;

    double out_sq = 0.0, ref_sq = 0.0;
    for (int s = 0; s < samples; s++) {
        double t = lo + (hi - lo) * (double)s / (double)(samples - 1);
        double y = (double)pc_eval(&tmp, (float)t, 0);
        double ref = cos(w * (t - t0));
        out_sq += y * y;
        ref_sq += ref * ref;
    }
    if (ref_sq <= 1e-30) return 1.0;
    return sqrt(out_sq / ref_sq);
}

typedef struct { int k; double base_gain, closure_gain, far_gain; } PCGainRow;

/* Fill rows[k] for k = 0..n/2: gain over the base period, over the closing
   segment, and over the period one period later.  Returns the row count. */
static int pc_gain_curve(const PeriodicCurve* pc, int samples,
                         PCGainRow rows[PC_MAX_MODES]) {
    int n = pc->n, modes = n / 2 + 1;
    double t0 = pc->t[0], P = pc->period;
    for (int k = 0; k < modes; k++) {
        rows[k].k = k;
        rows[k].base_gain    = pc_mode_gain(pc, k, t0, t0 + P, samples);
        rows[k].closure_gain = pc_mode_gain(pc, k, (double)pc->t[n - 1], t0 + P, samples);
        rows[k].far_gain     = pc_mode_gain(pc, k, t0 + P, t0 + 2.0 * P, samples);
    }
    return modes;
}

static double pc_max_of(const PCGainRow rows[PC_MAX_MODES], int modes,
                        int which /* 0 base, 1 closure, 2 far */) {
    double worst = 0.0;
    for (int k = 0; k < modes; k++) {
        double g = which == 0 ? rows[k].base_gain
                 : which == 1 ? rows[k].closure_gain : rows[k].far_gain;
        if (g > worst) worst = g;
    }
    return worst;
}

/* Build the periodic continuation.  gain_cap <= 0 disables the gain check;
   otherwise init fails with PC_ERR_GAIN when the closing segment amplifies any
   resolvable mode beyond gain_cap.  In that case *out_gain (if non-NULL) gets
   the offending gain so the caller can report it. */
static PCStatus periodic_init(PeriodicCurve* pc, const float* t, const float* v,
                              int n, float period, float gain_cap, double* out_gain) {
    if (out_gain) *out_gain = 0.0;
    memset(pc, 0, sizeof(*pc));
    if (!isfinite(period) || period <= 0.0f) return PC_ERR_PERIOD;
    if (n < 3) return PC_ERR_FEW_KEYS;
    if (n > MAX_PC_KEYS) return PC_ERR_CAPACITY;
    for (int i = 0; i < n; i++) {
        if (!isfinite(t[i]) || !isfinite(v[i])) return PC_ERR_SPACING;
        if (i > 0 && !(t[i] > t[i - 1])) return PC_ERR_SPACING;
    }
    if (!(((double)t[n - 1] - (double)t[0]) < (double)period)) return PC_ERR_SPACING;

    pc->n = n;
    pc->period = period;
    for (int i = 0; i < n; i++) { pc->t[i] = t[i]; pc->v[i] = v[i]; }

    PCStatus st = pc_solve(pc, pc->v);
    if (st != PC_OK) return st;

    if (gain_cap > 0.0f) {
        PCGainRow rows[PC_MAX_MODES];
        int modes = pc_gain_curve(pc, 257, rows);
        double g = pc_max_of(rows, modes, 1);
        if (out_gain) *out_gain = g;
        if (!(g <= (double)gain_cap)) return PC_ERR_GAIN;
    }
    return PC_OK;
}

/* ------------------------------------------------------------------ *
 * Whole-segment rebuild reference, used only to quantify the failure the
 * periodic spline avoids.  It fits one open natural cubic spline to the same
 * keys (free ends: M[0] = M[n-1] = 0) and continues the final cubic past the
 * last key, which is the usual "rebuild and extrapolate" behaviour.
 * ------------------------------------------------------------------ */
static void ref_natural_M(const float* t, const float* v, int n, double* M) {
    double h[MAX_PC_KEYS], a[MAX_PC_KEYS], b[MAX_PC_KEYS], c[MAX_PC_KEYS], d[MAX_PC_KEYS];
    for (int i = 0; i < n - 1; i++) h[i] = (double)t[i + 1] - (double)t[i];
    M[0] = M[n - 1] = 0.0;
    if (n <= 2) return;
    for (int i = 1; i < n - 1; i++) {
        a[i] = h[i - 1];
        b[i] = 2.0 * (h[i - 1] + h[i]);
        c[i] = h[i];
        d[i] = 6.0 * (((double)v[i + 1] - v[i]) / h[i]
                    - ((double)v[i] - v[i - 1]) / h[i - 1]);
    }
    for (int i = 2; i < n - 1; i++) {     /* Thomas on interior unknowns 1..n-2 */
        double f = a[i] / b[i - 1];
        b[i] -= f * c[i - 1];
        d[i] -= f * d[i - 1];
    }
    M[n - 2] = d[n - 2] / b[n - 2];
    for (int i = n - 3; i >= 1; i--) M[i] = (d[i] - c[i] * M[i + 1]) / b[i];
}

static int ref_seg(const float* t, int n, double x) {
    if (x <= (double)t[n - 1]) {
        int i = 0;
        while (i < n - 2 && x > (double)t[i + 1]) i++;
        return i;
    }
    return n - 2;                         /* continue the terminal cubic */
}

static double ref_value(const float* t, const float* v, const double* M, int n, double x) {
    int i = ref_seg(t, n, x);
    double h = (double)t[i + 1] - (double)t[i];
    double a = (double)t[i + 1] - x, b = x - (double)t[i];
    return M[i] * a * a * a / (6.0 * h) + M[i + 1] * b * b * b / (6.0 * h)
         + ((double)v[i] - M[i] * h * h / 6.0) * a / h
         + ((double)v[i + 1] - M[i + 1] * h * h / 6.0) * b / h;
}

static double ref_d2(const float* t, const double* M, int n, double x) {
    int i = ref_seg(t, n, x);
    double h = (double)t[i + 1] - (double)t[i];
    return (M[i] * ((double)t[i + 1] - x) + M[i + 1] * (x - (double)t[i])) / h;
}

/* worst-mode gain of the rebuilt curve on [lo, hi]; same definition as
   pc_mode_gain so the two numbers are directly comparable */
static double ref_gain_on(const float* t, int n, float period,
                          double lo, double hi, int samples) {
    double w0 = 2.0 * 3.14159265358979323846 / (double)period;
    double worst = 0.0;
    for (int k = 0; k <= n / 2; k++) {
        float vv[MAX_PC_KEYS];
        for (int i = 0; i < n; i++)
            vv[i] = (float)cos(w0 * (double)k * ((double)t[i] - (double)t[0]));
        double M[MAX_PC_KEYS];
        ref_natural_M(t, vv, n, M);
        double out_sq = 0.0, ref_sq = 0.0;
        for (int s = 0; s < samples; s++) {
            double x = lo + (hi - lo) * (double)s / (double)(samples - 1);
            double y = ref_value(t, vv, M, n, x);
            double r = cos(w0 * (double)k * (x - (double)t[0]));
            out_sq += y * y;
            ref_sq += r * r;
        }
        if (ref_sq > 1e-30) {
            double g = sqrt(out_sq / ref_sq);
            if (g > worst) worst = g;
        }
    }
    return worst;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

/* one period of a loop, deliberately returning far from where it started */
static const float LOOP_T[8] = { 0, 1, 2, 3, 4, 5, 6, 7 };
static const float LOOP_V[8] = { 0.0f, 1.0f, 0.5f, -0.5f, 0.0f, 0.8f, 0.2f, -0.3f };
#define LOOP_N 8
#define LOOP_P 8.0f

static void print_gain_curve(const PeriodicCurve* pc) {
    PCGainRow rows[PC_MAX_MODES];
    int modes = pc_gain_curve(pc, 401, rows);
    printf("# gain curve (k, base, closure, far)\n");
    for (int k = 0; k < modes; k++)
        printf("gain k=%d base=%.6f closure=%.6f far=%.6f\n",
               rows[k].k, rows[k].base_gain, rows[k].closure_gain, rows[k].far_gain);
}

int main(void) {
    Curve c;
    curve_init(&c);
    check("add_first", curve_add(&c, 0.0f, 0.0f, 0.0f, 1.0f) == 0);
    check("add_second", curve_add(&c, 1.0f, 1.0f, 1.0f, 0.0f) == 1);
    check("clamp_left", fabsf(curve_eval(&c, -1.0f)) < 1e-6f);
    check("clamp_right", fabsf(curve_eval(&c, 2.0f) - 1.0f) < 1e-6f);
    check("mid_monotonic", curve_eval(&c, 0.5f) > 0.0f && curve_eval(&c, 0.5f) < 1.0f);
    check("endpoints_exact", fabsf(curve_eval(&c, 0.0f)) < 1e-6f && fabsf(curve_eval(&c, 1.0f) - 1.0f) < 1e-6f);

    Curve flat;
    curve_init(&flat);
    curve_add(&flat, 0.0f, 0.5f, 0.0f, 0.0f);
    curve_add(&flat, 1.0f, 0.5f, 0.0f, 0.0f);
    check("flat_curve", fabsf(curve_eval(&flat, 0.5f) - 0.5f) < 1e-5f);
    check("flat_deriv_zero", fabsf(curve_deriv(&flat, 0.5f)) < 1e-3f);

    check("reject_unsorted", curve_add(&flat, 0.5f, 1.0f, 0, 0) == -1);

    Curve full;
    curve_init(&full);
    int rc = 0;
    for (int i = 0; i <= MAX_KEYS; i++) rc = curve_add(&full, (float)i, (float)i, 0, 0);
    check("capacity", rc == -1 && full.count == MAX_KEYS);

    /* ---------------- periodic C2 continuation ---------------- */

    PeriodicCurve pc;
    PCStatus st = periodic_init(&pc, LOOP_T, LOOP_V, LOOP_N, LOOP_P, 0.0f, NULL);
    check("pc_init_ok", st == PC_OK && pc.n == LOOP_N);
    if (st != PC_OK) {
        printf("pc_init_reason=%s\n", pc_status_str(st));
        printf("all_ok=0\n");
        return 1;
    }

    float key_err = 0.0f;
    for (int i = 0; i < LOOP_N; i++) {
        float e = fabsf(pc_eval(&pc, LOOP_T[i], 0) - LOOP_V[i]);
        if (e > key_err) key_err = e;
    }
    check("pc_interpolates_keys", key_err < 1e-4f);

    float period_err = 0.0f;
    for (int s = 0; s <= 64; s++) {
        float t = LOOP_P * (float)s / 64.0f + 0.013f;
        float y0 = pc_eval(&pc, t, 0);
        float e1 = fabsf(pc_eval(&pc, t + LOOP_P, 0) - y0);
        float e2 = fabsf(pc_eval(&pc, t + 13.0f * LOOP_P, 0) - y0);
        if (e1 > period_err) period_err = e1;
        if (e2 > period_err) period_err = e2;
    }
    check("pc_periodic_exact", period_err < 1e-4f);

    /* exact one-sided limits at every knot, including the wrap seam */
    double c0 = pc_knot_jump(&pc, 0);
    double c1 = pc_knot_jump(&pc, 1);
    double c2 = pc_knot_jump(&pc, 2);
    printf("pc_seam_c0=%.3e pc_seam_c1=%.3e pc_seam_c2=%.3e\n", c0, c1, c2);
    check("pc_c0", c0 < 1e-5);
    check("pc_c1", c1 < 1e-4);
    check("pc_c2", c2 < 1e-3);

    print_gain_curve(&pc);

    PCGainRow rows[PC_MAX_MODES];
    int modes = pc_gain_curve(&pc, 401, rows);
    double base_gain = pc_max_of(rows, modes, 0);
    double close_gain = pc_max_of(rows, modes, 1);
    double far_gain = pc_max_of(rows, modes, 2);
    printf("pc_base_gain=%.6f pc_closure_gain=%.6f pc_far_gain=%.6f\n",
           base_gain, close_gain, far_gain);
    check("pc_base_gain_bounded", base_gain < 1.01);
    check("pc_closure_gain_bounded", close_gain < 1.5);
    check("pc_far_equals_base", fabs(far_gain - base_gain) < 1e-3 && far_gain < 1.01);

    /* ---------------- whole-segment rebuild for comparison ---------------- */

    double Mn[8];
    ref_natural_M(LOOP_T, LOOP_V, LOOP_N, Mn);
    double rb_seam = fabs(ref_d2(LOOP_T, Mn, LOOP_N, (double)LOOP_P) - Mn[0]);
    double rb_far = ref_gain_on(LOOP_T, LOOP_N, LOOP_P, (double)LOOP_P, 2.0 * (double)LOOP_P, 401);
    double rb_far5 = ref_gain_on(LOOP_T, LOOP_N, LOOP_P, 5.0 * (double)LOOP_P, 6.0 * (double)LOOP_P, 401);
    printf("rebuild_seam_d2_jump=%.6f rebuild_far_gain=%.4f rebuild_far5_gain=%.1f\n",
           rb_seam, rb_far, rb_far5);
    check("rebuild_seam_jump_present", rb_seam > 0.1);
    check("rebuild_far_gain_blows_up", rb_far > 10.0);
    check("rebuild_gain_grows_with_distance", rb_far5 > rb_far);
    check("periodic_seam_beats_rebuild", c2 < rb_seam);
    check("periodic_far_gain_beats_rebuild", far_gain * 10.0 < rb_far);

    /* ---------------- gain cap as a limiter ----------------
       A loop whose keys stop well short of the period (a long closing gap)
       amplifies its top resolved mode over the synthesised segment; the cap
       must accept a loose bound, reject a tight one, and report the measured
       gain so the caller can say why. */
    static const float GAP_T[5] = { 0.0f, 1.0f, 2.0f, 3.0f, 3.5f };
    static const float GAP_V[5] = { 0.0f, 1.0f, -1.0f, 1.0f, -1.0f };
    PeriodicCurve capped;
    double gap_closure = 0.0;
    check("pc_gap_init", periodic_init(&capped, GAP_T, GAP_V, 5, 8.0f, 0.0f, NULL) == PC_OK);
    {
        PCGainRow grows[PC_MAX_MODES];
        int gm = pc_gain_curve(&capped, 257, grows);
        gap_closure = pc_max_of(grows, gm, 1);
    }
    printf("pc_gap_closure_gain=%.4f\n", gap_closure);
    check("pc_gap_closure_amplifies", gap_closure > 1.5);
    check("pc_gain_cap_accept", periodic_init(&capped, GAP_T, GAP_V, 5,
                                              8.0f, (float)gap_closure + 0.1f, NULL) == PC_OK);
    double rejected_gain = 0.0;
    PCStatus cap_st = periodic_init(&capped, GAP_T, GAP_V, 5,
                                    8.0f, (float)gap_closure - 0.1f, &rejected_gain);
    check("pc_gain_cap_reject", cap_st == PC_ERR_GAIN && rejected_gain > 0.0);
    printf("pc_gain_cap_reason=%s measured=%.4f\n", pc_status_str(cap_st), rejected_gain);

    /* ---------------- rejection paths (each prints its reason) ---------------- */

    PeriodicCurve rej;
    double dummy;
    PCStatus r1 = periodic_init(&rej, LOOP_T, LOOP_V, LOOP_N, 0.0f, 0, &dummy);
    check("reject_period_zero", r1 == PC_ERR_PERIOD);
    printf("reject_period_zero_reason=%s\n", pc_status_str(r1));

    PCStatus r2 = periodic_init(&rej, LOOP_T, LOOP_V, LOOP_N, -2.0f, 0, &dummy);
    check("reject_period_negative", r2 == PC_ERR_PERIOD);
    printf("reject_period_negative_reason=%s\n", pc_status_str(r2));

    PCStatus r3 = periodic_init(&rej, LOOP_T, LOOP_V, LOOP_N, NAN, 0, &dummy);
    check("reject_period_nan", r3 == PC_ERR_PERIOD);
    printf("reject_period_nan_reason=%s\n", pc_status_str(r3));

    PCStatus r4 = periodic_init(&rej, LOOP_T, LOOP_V, 2, LOOP_P, 0, &dummy);
    check("reject_two_keys", r4 == PC_ERR_FEW_KEYS);
    printf("reject_two_keys_reason=%s\n", pc_status_str(r4));

    PCStatus r5 = periodic_init(&rej, LOOP_T, LOOP_V, LOOP_N, 7.0f, 0, &dummy);
    check("reject_span_equals_period", r5 == PC_ERR_SPACING);
    printf("reject_span_reason=%s\n", pc_status_str(r5));

    float dup_t[8] = { 0, 1, 2, 2, 4, 5, 6, 7 };
    PCStatus r6 = periodic_init(&rej, dup_t, LOOP_V, LOOP_N, LOOP_P, 0, &dummy);
    check("reject_duplicate_time", r6 == PC_ERR_SPACING);

    float many_t[MAX_KEYS + 1], many_v[MAX_KEYS + 1];
    for (int i = 0; i <= MAX_KEYS; i++) { many_t[i] = (float)i; many_v[i] = (float)i; }
    PCStatus r7 = periodic_init(&rej, many_t, many_v, MAX_KEYS + 1, (float)MAX_KEYS, 0, &dummy);
    check("reject_capacity", r7 == PC_ERR_CAPACITY);
    printf("reject_capacity_reason=%s\n", pc_status_str(r7));

    /* ---------------- non-uniform spacing still C2 and periodic ------------ */

    float nu_t[5] = { 0.0f, 0.7f, 2.9f, 3.0f, 6.1f };
    float nu_v[5] = { 1.0f, -2.0f, 0.5f, 3.0f, 2.0f };
    PeriodicCurve npc;
    check("pc_nonuniform_init", periodic_init(&npc, nu_t, nu_v, 5, 7.3f, 0, &dummy) == PC_OK);
    double nu_c2 = pc_knot_jump(&npc, 2);
    float nu_per = 0.0f;
    for (int s = 0; s <= 64; s++) {
        float t = 7.3f * (float)s / 64.0f;
        float e = fabsf(pc_eval(&npc, t + 7.3f, 0) - pc_eval(&npc, t, 0));
        if (e > nu_per) nu_per = e;
    }
    printf("pc_nonuniform_c2=%.3e\n", nu_c2);
    check("pc_nonuniform_c2", nu_c2 < 1e-3);
    check("pc_nonuniform_periodic", nu_per < 1e-4);

    /* ---------------- measured upper bounds, swept over configurations ----
       Establish the advertised bounds empirically rather than from one lucky
       loop:
         * well-conditioned loops (no near-coincident keys): C2 to roundoff and
           exact periodicity, since |f(t+P)-f(t)| is then just the input's float
           resolution;
         * loops with one clustered pair: still C2 to roundoff, still gain-neutral
           far out, and the base-window gain stays under a finite constant.
       Note: with a near-coincident pair the slope is large, so comparing f(t)
       and f(t+P) as floats is limited by the ulp of t+P, not by the spline; the
       periodicity assertion is therefore made only on the well-conditioned set. */
    {
        double sweep_c2 = 0.0, sweep_period_err = 0.0;
        double sweep_base = 0.0, sweep_far_over_base = 1e9;
        unsigned seed = 0x9e3779b9u;
        for (int trial = 0; trial < 400; trial++) {
            int clustered = trial >= 200;
            seed = seed * 1664525u + 1013904223u;
            int n = 3 + (int)(seed % 14u);
            seed = seed * 1664525u + 1013904223u;
            float period = 2.0f + (float)(seed % 60u);
            double gaps[MAX_PC_KEYS], sum = 0.0;
            for (int i = 0; i < n; i++) {
                seed = seed * 1664525u + 1013904223u;
                gaps[i] = 0.4 + (double)(seed % 1000u) / 500.0;
                sum += gaps[i];
            }
            if (clustered) {
                seed = seed * 1664525u + 1013904223u;
                int cluster = (int)(seed % (unsigned)n);
                seed = seed * 1664525u + 1013904223u;
                gaps[cluster] = pow(10.0, -(2.0 + (double)(seed % 4u)));
            }
            sum = 0.0;
            for (int i = 0; i < n; i++) sum += gaps[i];
            float ts[16], vs[16];
            double acc = 0.0, scale = (double)period / sum * 0.999;
            for (int i = 0; i < n; i++) {
                ts[i] = (float)acc;
                seed = seed * 1664525u + 1013904223u;
                vs[i] = (float)((double)(seed % 2000u) / 1000.0 - 1.0);
                acc += gaps[i] * scale;
            }
            PeriodicCurve spc;
            if (periodic_init(&spc, ts, vs, n, period, 0.0f, NULL) != PC_OK) continue;
            double sc2 = pc_knot_jump(&spc, 2);
            if (sc2 > sweep_c2) sweep_c2 = sc2;
            if (!clustered) {
                for (int s = 0; s <= 16; s++) {
                    float t0s = period * (float)s / 16.0f + 1e-3f;
                    double e = fabs(pc_eval(&spc, t0s + period, 0) - pc_eval(&spc, t0s, 0));
                    if (e > sweep_period_err) sweep_period_err = e;
                }
            }
            PCGainRow srows[PC_MAX_MODES];
            int sm = pc_gain_curve(&spc, 129, srows);
            double b = pc_max_of(srows, sm, 0);
            double f = pc_max_of(srows, sm, 2);
            if (b > sweep_base) sweep_base = b;
            if (f / b < sweep_far_over_base) sweep_far_over_base = f / b;
        }
        printf("sweep_c2_max=%.3e sweep_period_err_max=%.3e sweep_base_gain_max=%.4f "
               "sweep_far_over_base_min=%.6f\n",
               sweep_c2, sweep_period_err, sweep_base, sweep_far_over_base);
        check("sweep_c2_upper_bound", sweep_c2 < 1e-3);
        check("sweep_periodic_upper_bound", sweep_period_err < 1e-4);
        check("sweep_base_gain_finite_bound", sweep_base < 3.0);
        check("sweep_far_gain_matches_base",
              fabs(sweep_far_over_base - 1.0) < 1e-3);
    }


    /* ---------------- far window does not accelerate ---------------- */

    float base_dd = 0.0f, far_dd = 0.0f;
    for (int win = 0; win < 2; win++) {
        float t0w = win == 0 ? 8.0f * LOOP_P : 0.0f;
        float prev = pc_eval(&pc, t0w, 0);
        float cur = pc_eval(&pc, t0w + 0.125f, 0);
        for (int s = 2; s <= 128; s++) {
            float t = t0w + 0.125f * (float)s;
            float next = pc_eval(&pc, t, 0);
            float dd = fabsf(next - 2.0f * cur + prev);
            if (win == 0 && dd > far_dd) far_dd = dd;
            if (win == 1 && dd > base_dd) base_dd = dd;
            prev = cur; cur = next;
        }
    }
    printf("pc_base_2nd_diff=%.6f pc_far_2nd_diff=%.6f\n", base_dd, far_dd);
    check("pc_no_acceleration_far", far_dd < 1.5f * base_dd + 1e-3f);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
