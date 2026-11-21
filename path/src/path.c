
/* path: polyline arc-length table, constant-speed sampling, curvature and
   closest point queries.

   Second half: arc-length reparameterisation of a sampled smooth curve.
   A cumulative chord table is built by uniform or curvature-adaptive
   sampling, distances are inverted to (segment, ratio) by binary search,
   and the result is checked against an analytic arc length: position-error
   bound, speed-fluctuation bound and the log-log convergence slope. */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define MAX_PTS 32

#define ARC_MAX 1024
#define ARC_PI 3.14159265358979323846
#define ARC_EPS 1e-9

/* return codes */
#define ARC_OK 0
#define ARC_E_ARG (-1)      /* bad argument / non-finite input */
#define ARC_E_ZERO (-2)     /* zero-length segment */
#define ARC_E_RANGE (-3)    /* distance outside [0, total] */
#define ARC_E_CAP (-4)      /* node capacity exhausted */
#define ARC_E_STATE (-5)    /* table not built / inconsistent */

typedef struct { float x, y; } P2;

typedef struct {
    P2 p[MAX_PTS];
    float cum[MAX_PTS];    /* cumulative arc length at each point */
    int n;
    float total;
} Path;

static void path_init(Path* pa) { memset(pa, 0, sizeof(*pa)); }

static int path_add(Path* pa, float x, float y) {
    if (pa->n >= MAX_PTS) return -1;
    if (!isfinite(x) || !isfinite(y)) return -1;
    if (pa->n > 0) {
        P2 prev = pa->p[pa->n - 1];
        float dx = x - prev.x, dy = y - prev.y;
        float d = sqrtf(dx * dx + dy * dy);
        if (d < 1e-6f) return -1;               /* zero-length segment rejected */
        pa->cum[pa->n] = pa->cum[pa->n - 1] + d;
        pa->total = pa->cum[pa->n];
    }
    pa->p[pa->n].x = x;
    pa->p[pa->n].y = y;
    pa->n++;
    return pa->n - 1;
}

/* constant speed sample: s in [0, total] */
static int path_sample(const Path* pa, float s, P2* out, float* seg_t) {
    if (pa->n < 2) return -1;
    if (!isfinite(s) || s < 0.0f || s > pa->total) return -1;
    for (int i = 1; i < pa->n; i++) {
        if (s <= pa->cum[i]) {
            float seg = pa->cum[i] - pa->cum[i - 1];
            if (seg <= 0.0f) return -1;         /* zero-length segment: reject */
            float t = (s - pa->cum[i - 1]) / seg;
            out->x = pa->p[i - 1].x + (pa->p[i].x - pa->p[i - 1].x) * t;
            out->y = pa->p[i - 1].y + (pa->p[i].y - pa->p[i - 1].y) * t;
            if (seg_t) *seg_t = t;
            return i - 1;
        }
    }
    return -1;
}

/* closest point on path to q, returns distance and arc length */
static float path_closest(const Path* pa, P2 q, float* out_s) {
    float best = 1e30f, best_s = 0.0f;
    for (int i = 1; i < pa->n; i++) {
        P2 a = pa->p[i - 1], b = pa->p[i];
        float vx = b.x - a.x, vy = b.y - a.y;
        float wx = q.x - a.x, wy = q.y - a.y;
        float vv = vx * vx + vy * vy;
        float t = vv > 0.0f ? (wx * vx + wy * vy) / vv : 0.0f;
        if (t < 0.0f) t = 0.0f;
        if (t > 1.0f) t = 1.0f;
        float px = a.x + vx * t, py = a.y + vy * t;
        float dx = q.x - px, dy = q.y - py;
        float d = sqrtf(dx * dx + dy * dy);
        if (d < best) {
            best = d;
            best_s = pa->cum[i - 1] + (pa->cum[i] - pa->cum[i - 1]) * t;
        }
    }
    if (out_s) *out_s = best_s;
    return best;
}

/* signed curvature at segment i (turn angle over length) */
static float path_curvature(const Path* pa, int i) {
    if (i < 1 || i >= pa->n - 1) return 0.0f;
    P2 a = pa->p[i - 1], b = pa->p[i], c = pa->p[i + 1];
    float abx = b.x - a.x, aby = b.y - a.y;
    float bcx = c.x - b.x, bcy = c.y - b.y;
    float cross = abx * bcy - aby * bcx;
    float la = sqrtf(abx * abx + aby * aby);
    float lb = sqrtf(bcx * bcx + bcy * bcy);
    float lc = sqrtf((c.x - a.x) * (c.x - a.x) + (c.y - a.y) * (c.y - a.y));
    if (la * lb * lc == 0.0f) return 0.0f;
    return 2.0f * cross / (la * lb * lc);
}

/* ------------------------------------------------------------------ */
/* arc-length reparameterisation of a sampled smooth curve            */

typedef struct { double x, y, s, u; } ArcNode;

typedef struct {
    ArcNode node[ARC_MAX];
    int n;
    double total;       /* sum of chord lengths */
    double max_sag;     /* largest sampled deviation from a chord */
    double max_turn;    /* largest sampled turning angle over a segment */
} ArcPath;

/* parametric curve point; u is the raw parameter, not arc length */
typedef void (*ArcEval)(double u, double* x, double* y, void* ctx);
/* analytic arc length from the first node parameter u0 to u */
typedef double (*ArcLen)(double u, void* ctx);
/* exact point at true arc-length distance s from the curve start */
typedef void (*ArcRefAt)(double s, double* x, double* y, void* ctx);

static double arc_clamp(double v, double lo, double hi) {
    return v < lo ? lo : (v > hi ? hi : v);
}

static double arc_angle_between(double ax, double ay, double bx, double by) {
    double la = hypot(ax, ay), lb = hypot(bx, by);
    if (la < ARC_EPS || lb < ARC_EPS) return 0.0;
    double c = (ax * bx + ay * by) / (la * lb);
    return acos(arc_clamp(c, -1.0, 1.0));
}

static double arc_point_chord_dist(double px, double py,
                                   double ax, double ay, double bx, double by) {
    double vx = bx - ax, vy = by - ay;
    double vv = vx * vx + vy * vy;
    double t = vv > 0.0 ? ((px - ax) * vx + (py - ay) * vy) / vv : 0.0;
    t = arc_clamp(t, 0.0, 1.0);
    double dx = px - (ax + vx * t), dy = py - (ay + vy * t);
    return hypot(dx, dy);
}

/* append a node, rejecting zero-length segments and non-finite points */
static int arc_push(ArcPath* ap, double x, double y, double u) {
    if (ap->n >= ARC_MAX) return ARC_E_CAP;
    if (!isfinite(x) || !isfinite(y) || !isfinite(u)) return ARC_E_ARG;
    if (ap->n > 0) {
        const ArcNode* prev = &ap->node[ap->n - 1];
        double d = hypot(x - prev->x, y - prev->y);
        if (!(d > ARC_EPS)) return ARC_E_ZERO;  /* zero-length segment rejected */
        double s = prev->s + d;
        if (!isfinite(s)) return ARC_E_ARG;
        ap->node[ap->n].s = s;
        ap->total = s;
    }
    ap->node[ap->n].x = x;
    ap->node[ap->n].y = y;
    ap->node[ap->n].u = u;
    ap->n++;
    return ARC_OK;
}

static void arc_reset(ArcPath* ap) {
    ap->n = 0;
    ap->total = 0.0;
    ap->max_sag = 0.0;
    ap->max_turn = 0.0;
}

/* measure max deviation from the chord and total turning on one segment
   using four sub-chords; cheap, but resolves curvature the nodes miss */
static void arc_measure(const ArcEval f, void* ctx, const ArcNode* a,
                        const ArcNode* b, double* sag_out, double* turn_out) {
    double qx[5], qy[5];
    qx[0] = a->x; qy[0] = a->y;
    qx[4] = b->x; qy[4] = b->y;
    for (int k = 1; k < 4; k++) {
        double u = a->u + (b->u - a->u) * (k * 0.25);
        f(u, &qx[k], &qy[k], ctx);
    }
    double sag = 0.0, turn = 0.0;
    for (int k = 1; k < 4; k++) {
        double d = arc_point_chord_dist(qx[k], qy[k], a->x, a->y, b->x, b->y);
        if (d > sag) sag = d;
    }
    for (int k = 1; k < 4; k++) {
        turn += arc_angle_between(qx[k] - qx[k - 1], qy[k] - qy[k - 1],
                                  qx[k + 1] - qx[k], qy[k + 1] - qy[k]);
    }
    *sag_out = sag;
    *turn_out = turn;
}

static void arc_metrics(ArcPath* ap, const ArcEval f, void* ctx) {
    ap->max_sag = 0.0;
    ap->max_turn = 0.0;
    for (int i = 1; i < ap->n; i++) {
        double sag, turn;
        arc_measure(f, ctx, &ap->node[i - 1], &ap->node[i], &sag, &turn);
        if (sag > ap->max_sag) ap->max_sag = sag;
        if (turn > ap->max_turn) ap->max_turn = turn;
    }
}

/* sample u0..u1 at nseg uniform parameter steps */
static int arc_build_uniform(ArcPath* ap, const ArcEval f, void* ctx,
                             double u0, double u1, int nseg) {
    if (!ap || !f) return ARC_E_ARG;
    if (nseg < 1 || nseg > ARC_MAX - 1) return ARC_E_ARG;
    if (!isfinite(u0) || !isfinite(u1) || !(u1 > u0)) return ARC_E_ARG;
    arc_reset(ap);
    for (int i = 0; i <= nseg; i++) {
        double u = (i == nseg) ? u1 : u0 + (u1 - u0) * i / nseg;
        double x, y;
        f(u, &x, &y, ctx);
        int rc = arc_push(ap, x, y, u);
        if (rc != ARC_OK) return rc;
    }
    arc_metrics(ap, f, ctx);
    return ARC_OK;
}

/* split every segment whose midpoint sagitta exceeds tol, so node spacing
   scales like 1/sqrt(curvature): denser where the curve bends harder */
static int arc_build_adaptive(ArcPath* ap, const ArcEval f, void* ctx,
                              double u0, double u1, double tol, int max_nodes) {
    if (!ap || !f) return ARC_E_ARG;
    if (!isfinite(u0) || !isfinite(u1) || !(u1 > u0)) return ARC_E_ARG;
    if (!(tol > 0.0) || !isfinite(tol)) return ARC_E_ARG;
    if (max_nodes < 2) return ARC_E_ARG;
    if (max_nodes > ARC_MAX) max_nodes = ARC_MAX;

    double u[ARC_MAX];
    int n = 2;
    u[0] = u0; u[1] = u1;
    for (int iter = 0; iter < 40; iter++) {
        int split = 0;
        for (int i = 0; i < n - 1; i++) {
            double a = u[i], b = u[i + 1];
            double ux, uy, vx, vy, mx, my;
            f(a, &ux, &uy, ctx);
            f(b, &vx, &vy, ctx);
            f(0.5 * (a + b), &mx, &my, ctx);
            double sag = arc_point_chord_dist(mx, my,
                                              ux, uy, vx, vy);
            if (sag > tol) {
                if (n + 1 > max_nodes) return ARC_E_CAP;
                memmove(&u[i + 2], &u[i + 1], (size_t)(n - i - 1) * sizeof(double));
                u[i + 1] = 0.5 * (a + b);
                n++;
                split = 1;
                i++;                    /* let the next pass finish the halves */
            }
        }
        if (!split) break;
    }

    arc_reset(ap);
    for (int i = 0; i < n; i++) {
        double x, y;
        f(u[i], &x, &y, ctx);
        int rc = arc_push(ap, x, y, u[i]);
        if (rc != ARC_OK) return rc;
    }
    arc_metrics(ap, f, ctx);
    return ARC_OK;
}

/* distance -> (segment index, intra-segment ratio); binary search on s */
static int arc_lookup(const ArcPath* ap, double s,
                      double* x, double* y, int* seg, double* ratio) {
    if (!ap || ap->n < 2) return ARC_E_STATE;
    if (!isfinite(s)) return ARC_E_ARG;
    if (s < 0.0 || s > ap->total) return ARC_E_RANGE;
    int lo = 1, hi = ap->n - 1;
    while (lo < hi) {
        int mid = lo + (hi - lo) / 2;
        if (ap->node[mid].s < s) lo = mid + 1;
        else hi = mid;
    }
    int i = lo;
    double h = ap->node[i].s - ap->node[i - 1].s;
    if (!(h > 0.0)) return ARC_E_ZERO;          /* zero-length segment: reject */
    double t = arc_clamp((s - ap->node[i - 1].s) / h, 0.0, 1.0);
    if (x) *x = ap->node[i - 1].x + (ap->node[i].x - ap->node[i - 1].x) * t;
    if (y) *y = ap->node[i - 1].y + (ap->node[i].y - ap->node[i - 1].y) * t;
    if (seg) *seg = i - 1;
    if (ratio) *ratio = t;
    return ARC_OK;
}

/* structural check: strict, finite, strictly increasing length table */
static int arc_validate(const ArcPath* ap) {
    if (!ap || ap->n < 2) return ARC_E_STATE;
    if (!isfinite(ap->total) || !(ap->total > 0.0)) return ARC_E_STATE;
    double acc = 0.0;
    for (int i = 1; i < ap->n; i++) {
        const ArcNode* a = &ap->node[i - 1];
        const ArcNode* b = &ap->node[i];
        double d = hypot(b->x - a->x, b->y - a->y);
        if (!(d > ARC_EPS)) return ARC_E_ZERO;
        if (!(b->s > a->s)) return ARC_E_STATE;
        acc += d;
        double tol = 1e-9 * (1.0 + fabs(acc));
        if (fabs(acc - b->s) > tol) return ARC_E_STATE;
    }
    if (fabs(acc - ap->total) > 1e-9 * (1.0 + fabs(acc))) return ARC_E_STATE;
    return ARC_OK;
}

/* relative error of the cumulative table against analytic arc length */
static double arc_length_rel_error(const ArcPath* ap, const ArcLen L, void* ctx) {
    if (!ap || ap->n < 2 || !L) return NAN;
    double a = L(ap->node[0].u, ctx);
    double b = L(ap->node[ap->n - 1].u, ctx);
    double total = b - a;
    if (!(total > 0.0)) return NAN;
    return fabs(ap->total - total) / total;
}

/* measured speed fluctuation: each table segment is traversed at constant
   speed, so the spread of (true arc / chord) is the speed error */
static double arc_speed_fluctuation(const ArcPath* ap, const ArcLen L, void* ctx) {
    double worst = 0.0;
    for (int i = 1; i < ap->n; i++) {
        double chord = ap->node[i].s - ap->node[i - 1].s;
        double truth = L(ap->node[i].u, ctx) - L(ap->node[i - 1].u, ctx);
        if (!(chord > 0.0) || !(truth >= 0.0)) return NAN;
        double dev = fabs(truth / chord - 1.0);
        if (dev > worst) worst = dev;
    }
    return worst;
}

/* table-only upper bound on speed fluctuation: on a chord, chord/arc >=
   cos(turn/2), so the deviation is at most 1 - cos(turn_max/2). Conservative
   by roughly 3x versus the sharp sinc bound, but converges at the same order. */
static double arc_speed_bound(const ArcPath* ap) {
    double th = arc_clamp(ap->max_turn, 0.0, ARC_PI);
    return 1.0 - cos(th / 2.0);
}

/* upper bound on position error: sampled chord sagitta (2x for the point
   between sub-samples) plus the accumulated chord-vs-arc drift */
static double arc_position_bound(const ArcPath* ap, double analytic_total) {
    double drift = analytic_total - ap->total;
    if (!(drift > 0.0)) drift = 0.0;
    return 2.0 * ap->max_sag + drift;
}

/* sampled max |table point - analytic point at the same distance s| */
static double arc_max_position_error(const ArcPath* ap, const ArcRefAt ref,
                                     void* ctx, int nsamples) {
    if (!ap || !ref || nsamples < 1) return NAN;
    double worst = 0.0;
    for (int k = 0; k <= nsamples; k++) {
        double s = ap->total * k / nsamples;
        double x, y, rx, ry;
        if (arc_lookup(ap, s, &x, &y, NULL, NULL) != ARC_OK) return NAN;
        ref(s, &rx, &ry, ctx);
        double d = hypot(x - rx, y - ry);
        if (d > worst) worst = d;
    }
    return worst;
}

/* least-squares slope of log(err) against log(n) */
static double arc_fit_slope(const double* n, const double* err, int count) {
    if (!n || !err || count < 2) return NAN;
    double sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0;
    for (int i = 0; i < count; i++) {
        double x = log((double)n[i]);
        double e = err[i] > 0.0 && isfinite(err[i]) ? err[i] : 1e-300;
        double y = log(e);
        sx += x; sy += y; sxx += x * x; sxy += x * y;
    }
    double den = count * sxx - sx * sx;
    if (fabs(den) < 1e-300) return NAN;
    return (count * sxy - sx * sy) / den;
}

/* ------------------------------------------------------------------ */
/* analytic test curves                                                */

typedef struct { double r; } CircleCtx;

static void circle_eval(double u, double* x, double* y, void* ctx) {
    double r = ((CircleCtx*)ctx)->r;
    *x = r * cos(u);
    *y = r * sin(u);
}
static double circle_len(double u, void* ctx) { return ((CircleCtx*)ctx)->r * u; }
static void circle_ref(double s, double* x, double* y, void* ctx) {
    double r = ((CircleCtx*)ctx)->r;
    *x = r * cos(s / r);
    *y = r * sin(s / r);
}

static double par_S(double x) {
    return 0.5 * (x * sqrt(1.0 + x * x) + asinh(x));
}

typedef struct { double x0, x1; } ParCtx;

static void par_eval(double x, double* ox, double* oy, void* ctx) {
    (void)ctx;
    *ox = x;
    *oy = 0.5 * x * x;
}
static double par_len(double x, void* ctx) {
    ParCtx* c = (ParCtx*)ctx;
    return par_S(x) - par_S(c->x0);
}
static void par_ref(double s, double* ox, double* oy, void* ctx) {
    ParCtx* c = (ParCtx*)ctx;
    double target = par_S(c->x0) + s;
    double lo = c->x0, hi = c->x1;
    for (int i = 0; i < 100; i++) {
        double mid = 0.5 * (lo + hi);
        if (par_S(mid) < target) lo = mid; else hi = mid;
    }
    double x = 0.5 * (lo + hi);
    *ox = x;
    *oy = 0.5 * x * x;
}

/* ------------------------------------------------------------------ */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    Path pa;
    path_init(&pa);
    check("add0", path_add(&pa, 0, 0) == 0);
    check("add1", path_add(&pa, 3, 0) == 1);
    check("add2", path_add(&pa, 3, 4) == 2);
    check("total_len", fabsf(pa.total - 7.0f) < 1e-5f);

    P2 s;
    check("sample_mid", path_sample(&pa, 1.5f, &s, NULL) == 0 && fabsf(s.x - 1.5f) < 1e-5f);
    check("sample_corner", path_sample(&pa, 3.0f, &s, NULL) >= 0 && fabsf(s.x - 3.0f) < 1e-5f && fabsf(s.y) < 1e-5f);
    check("sample_end", path_sample(&pa, 7.0f, &s, NULL) >= 0 && fabsf(s.y - 4.0f) < 1e-4f);
    check("sample_out", path_sample(&pa, 8.0f, &s, NULL) == -1);

    check("dup_rejected", path_add(&pa, 3, 4) == -1);
    check("nan_rejected", path_add(&pa, NAN, 1.0f) == -1);

    float best_s = 0;
    float d = path_closest(&pa, (P2){ 1.5f, 2.0f }, &best_s);
    check("closest_dist", fabsf(d - 1.5f) < 1e-4f);
    check("closest_s", fabsf(best_s - 5.0f) < 1e-4f);

    Path straight;
    path_init(&straight);
    path_add(&straight, 0, 0); path_add(&straight, 1, 0); path_add(&straight, 2, 0);
    check("curvature_straight", fabsf(path_curvature(&straight, 1)) < 1e-5f);

    Path bend;
    path_init(&bend);
    path_add(&bend, 0, 0); path_add(&bend, 1, 0); path_add(&bend, 1, 1);
    check("curvature_turn", path_curvature(&bend, 1) > 0.0f);

    /* equal distance steps must stay equal-length on straight runs; the
       corner only turns the direction, it must not change the speed */
    {
        float ds = 0.05f, worst = 0.0f;
        int straight_steps = 0;
        for (float q = 0.0f; q + ds <= pa.total + 1e-6f; q += ds) {
            P2 a, b;
            int sa = path_sample(&pa, q, &a, NULL);
            int sb = path_sample(&pa, q + ds, &b, NULL);
            if (sa != sb || sa < 0) continue;
            straight_steps++;
            float e = fabsf(hypotf(b.x - a.x, b.y - a.y) - ds);
            if (e > worst) worst = e;
        }
        check("uniform_speed_straight", straight_steps > 50 && worst < 1e-4f);
    }

    /* ---- analytic circle ---- */
    CircleCtx cc = { 1.0 };
    ArcPath circ;
    int rc = arc_build_uniform(&circ, circle_eval, &cc, 0.0, 2.0 * ARC_PI, 64);
    check("circle_build", rc == ARC_OK && arc_validate(&circ) == ARC_OK);
    double circ_total = 2.0 * ARC_PI * cc.r;
    double circ_len_err = arc_length_rel_error(&circ, circle_len, &cc);
    double circ_pos = arc_max_position_error(&circ, circle_ref, &cc, 64 * 64);
    double circ_pos_b = arc_position_bound(&circ, circ_total);
    double circ_spd = arc_speed_fluctuation(&circ, circle_len, &cc);
    double circ_spd_b = arc_speed_bound(&circ);
    printf("circle_len_rel=%.6e\n", circ_len_err);
    printf("circle_pos=%.6e pos_bound=%.6e\n", circ_pos, circ_pos_b);
    printf("circle_speed=%.6e speed_bound=%.6e\n", circ_spd, circ_spd_b);
    check("circle_pos_within_bound", circ_pos <= circ_pos_b);
    check("circle_speed_within_bound", circ_spd <= circ_spd_b);
    check("circle_len_second_order", circ_len_err < 1e-3);

    double x, y, ratio;
    int seg;
    check("lookup_ok", arc_lookup(&circ, circ_total * 0.5, &x, &y, &seg, &ratio) == ARC_OK);
    check("lookup_ratio", ratio >= 0.0 && ratio <= 1.0 && seg >= 0 && seg < circ.n - 1);
    check("lookup_neg_rejected", arc_lookup(&circ, -1e-3, &x, &y, &seg, &ratio) == ARC_E_RANGE);
    check("lookup_past_rejected", arc_lookup(&circ, circ_total + 1e-3, &x, &y, &seg, &ratio) == ARC_E_RANGE);
    check("lookup_nan_rejected", arc_lookup(&circ, NAN, &x, &y, &seg, &ratio) == ARC_E_ARG);
    check("lookup_end", arc_lookup(&circ, circ.total, &x, &y, &seg, &ratio) == ARC_OK && fabs(ratio - 1.0) < 1e-12);
    check("lookup_past_table_total", arc_lookup(&circ, circ.total + 1e-9, &x, &y, &seg, &ratio) == ARC_E_RANGE);

    /* zero-length segment rejection: a constant curve has no arc length */
    {
        /* constant evaluator via circle with r = 0 */
        CircleCtx zero = { 0.0 };
        ArcPath degen;
        check("zero_segment_rejected",
              arc_build_uniform(&degen, circle_eval, &zero, 0.0, 1.0, 8) == ARC_E_ZERO);
        check("zero_adaptive_rejected",
              arc_build_adaptive(&degen, circle_eval, &zero, 0.0, 1.0, 1e-3, 256) == ARC_E_ZERO);
        check("degenerate_lookup_state", arc_lookup(&degen, 0.0, &x, &y, &seg, &ratio) != ARC_OK);
    }

    /* ---- parabola: non-uniform speed, curvature varies along the path ---- */
    ParCtx pc = { -1.0, 1.0 };
    double par_total = par_S(1.0) - par_S(-1.0);

    /* naive parameter sampling (equal du) has speed spread sqrt(2) on this
       path; arc-length reparameterisation must flatten it */
    double naive_min = 1e30, naive_max = 0.0;
    for (int i = 0; i <= 256; i++) {
        double uu = -1.0 + 2.0 * i / 256;
        double sp = sqrt(1.0 + uu * uu);           /* |dP/du| */
        if (sp < naive_min) naive_min = sp;
        if (sp > naive_max) naive_max = sp;
    }
    printf("naive_speed_spread=%.6f\n", naive_max / naive_min);

    enum { kn = 6 };
    double dens[kn], e_len[kn], e_pos[kn], e_spd[kn];
    int densities[kn] = { 16, 32, 64, 128, 256, 512 };
    int all_ok = 1, bounds_ok = 1, monotone = 1;
    for (int k = 0; k < kn; k++) {
        int nseg = densities[k];
        ArcPath ap;
        int brc = arc_build_uniform(&ap, par_eval, &pc, pc.x0, pc.x1, nseg);
        if (brc != ARC_OK || arc_validate(&ap) != ARC_OK) { all_ok = 0; continue; }
        dens[k] = nseg;
        e_len[k] = arc_length_rel_error(&ap, par_len, &pc);
        e_pos[k] = arc_max_position_error(&ap, par_ref, &pc, 32 * nseg);
        e_spd[k] = arc_speed_fluctuation(&ap, par_len, &pc);
        double pb = arc_position_bound(&ap, par_total);
        double sb = arc_speed_bound(&ap);
        if (!(e_pos[k] <= pb) || !(e_spd[k] <= sb)) bounds_ok = 0;
        if (k > 0 && !(e_pos[k] < e_pos[k - 1])) monotone = 0;
        if (nseg == 64) {
            printf("parabola_N64_len_rel=%.6e\n", e_len[k]);
            printf("parabola_N64_pos=%.6e pos_bound=%.6e\n", e_pos[k], pb);
            printf("parabola_N64_speed=%.6e speed_bound=%.6e\n", e_spd[k], sb);
        }
    }
    check("parabola_builds", all_ok);
    check("parabola_bounds_hold", bounds_ok);
    check("parabola_error_decreases", monotone);

    double slope_len = arc_fit_slope(dens, e_len, kn);
    double slope_pos = arc_fit_slope(dens, e_pos, kn);
    double slope_spd = arc_fit_slope(dens, e_spd, kn);
    const double theory = -2.0;
    const double slope_tol = 0.15;
    printf("slope_len=%.5f theory=%.1f dev=%+.5f\n", slope_len, theory, slope_len - theory);
    printf("slope_pos=%.5f theory=%.1f dev=%+.5f\n", slope_pos, theory, slope_pos - theory);
    printf("slope_spd=%.5f theory=%.1f dev=%+.5f\n", slope_spd, theory, slope_spd - theory);
    check("slope_len_matches_theory", fabs(slope_len - theory) <= slope_tol);
    check("slope_pos_matches_theory", fabs(slope_pos - theory) <= slope_tol);
    check("slope_spd_matches_theory", fabs(slope_spd - theory) <= slope_tol);

    /* ---- curvature-adaptive refinement ---- */
    {
        ArcPath ap;
        double tol = 1e-3;
        int arc = arc_build_adaptive(&ap, par_eval, &pc, pc.x0, pc.x1, tol, 1024);
        check("adaptive_build", arc == ARC_OK && arc_validate(&ap) == ARC_OK);
        check("adaptive_sag_respected", ap.max_sag <= tol + 1e-12);
        /* node density per unit arc near the nose (x=0) must exceed the
           density near the ends where curvature is lower */
        double nodes_center = 0.0, nodes_end = 0.0, arc_center, arc_end;
        for (int i = 1; i < ap.n; i++) {
            double um = 0.5 * (ap.node[i - 1].u + ap.node[i].u);
            double ds = ap.node[i].s - ap.node[i - 1].s;
            if (fabs(um) < 0.2) nodes_center += ds;
            if (fabs(um) > 0.8) nodes_end += ds;
        }
        (void)nodes_center; (void)nodes_end;
        /* compare node counts rather than lengths, per unit analytic arc */
        int cnt_center = 0, cnt_end = 0;
        for (int i = 1; i < ap.n; i++) {
            double um = 0.5 * (ap.node[i - 1].u + ap.node[i].u);
            if (fabs(um) < 0.2) cnt_center++;
            if (fabs(um) > 0.8) cnt_end++;
        }
        arc_center = par_S(0.2) - par_S(-0.2);
        arc_end = (par_S(1.0) - par_S(0.8)) + (par_S(-0.8) - par_S(-1.0));
        double dens_ratio = (cnt_center / arc_center) / (cnt_end / arc_end);
        printf("adaptive_nodes=%d dense_ratio=%.4f\n", ap.n, dens_ratio);
        check("adaptive_denser_at_high_curvature", dens_ratio > 1.15);
        double pb = arc_position_bound(&ap, par_total);
        double pe = arc_max_position_error(&ap, par_ref, &pc, 64 * 1024);
        double se = arc_speed_fluctuation(&ap, par_len, &pc);
        printf("adaptive_pos=%.6e pos_bound=%.6e\n", pe, pb);
        printf("adaptive_speed=%.6e speed_bound=%.6e\n", se, arc_speed_bound(&ap));
        check("adaptive_position_bound", pe <= pb);
        check("adaptive_speed_bound", se <= arc_speed_bound(&ap));

        /* tolerance sweep: error vs actual node count keeps the same slope */
        enum { kt = 5 };
        double td[kt], tp[kt];
        double tols[kt] = { 8e-3, 4e-3, 2e-3, 1e-3, 5e-4 };
        int adapt_ok = 1;
        for (int k = 0; k < kt; k++) {
            ArcPath ta;
            int trc = arc_build_adaptive(&ta, par_eval, &pc, pc.x0, pc.x1, tols[k], 1024);
            if (trc != ARC_OK) { adapt_ok = 0; td[k] = 1; tp[k] = 1; continue; }
            td[k] = ta.n - 1;
            tp[k] = arc_max_position_error(&ta, par_ref, &pc, 64 * (ta.n - 1));
        }
        check("adaptive_tolerance_builds", adapt_ok);
        double adapt_slope = arc_fit_slope(td, tp, kt);
        printf("adaptive_slope_pos=%.5f theory=%.1f dev=%+.5f\n",
               adapt_slope, theory, adapt_slope - theory);
        check("adaptive_slope_matches_theory", fabs(adapt_slope - theory) <= 0.25);
    }

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
