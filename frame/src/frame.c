
/* frame: parallel transport frame (PTF) along a polyline, with closure
   correction for closed paths.

   Transport rule: at every vertex the frame is rotated by the wrapped tangent
   difference angle, so headings are unwrapped accumulations (the naive
   per-vertex atan2 frame flips by pi at reversals; accumulating wrapped turns
   removes that flip). For a closed path the full-loop accumulation gives the
   holonomy, which is cross-checked against the analytic winding (the turning
   number theorem: total turning of a closed plane curve is 2pi * integer).
   The residual gap to the nearest integer turn is spread uniformly over the
   n transport steps and the frames are rebuilt, so the orientation error after
   one full loop is zero. A zero tangent edge, or a jump (half-turn flip, or an
   oversized per-step spread) after correction, is rejected. */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define MAX_FRAME 64
#define EPS 1e-6f          /* zero-vector guard */
#define CLOSE_EPS 1e-5f    /* numerical zero for the closed orientation error */
#define TWO_PI 6.28318530717958647692f
#define PI 3.14159265358979323846f
/* uniform spread larger than this per step is treated as a jump */
#define SPREAD_JUMP_LIM (0.5f * PI)
/* transported winding may disagree with the analytic one by at most this */
#define WIND_MATCH_TOL 0.25f

typedef struct { float x, y; } V2;

typedef struct {
    V2 p[MAX_FRAME];
    int n;
} Poly;

typedef struct {
    V2 tan;        /* unit tangent */
    V2 normal;     /* unit normal, tangent rotated 90 degrees ccw */
    float heading; /* transported frame angle (radians, unwrapped) */
} Frame2;

typedef struct {
    Frame2 f[MAX_FRAME];
    float delta[MAX_FRAME]; /* turn step before frame i; delta[0] is the seam */
    float holonomy;         /* full-loop accumulated turn (radians) */
    float residual;         /* gap to the nearest integer turn */
    float winding;          /* holonomy / 2pi */
    int n;
    int closed;
} FramePath;

static void poly_init(Poly* po) { memset(po, 0, sizeof(*po)); }

static int poly_add(Poly* po, float x, float y) {
    if (po->n >= MAX_FRAME) return -1;
    po->p[po->n].x = x;
    po->p[po->n].y = y;
    po->n++;
    return po->n - 1;
}

/* wrap angle into (-pi, pi] */
static float wrap_pi(float a) {
    while (a > PI) a -= TWO_PI;
    while (a <= -PI) a += TWO_PI;
    return a;
}

static V2 v2_unit(float x, float y, int* ok) {
    V2 r = { 0.0f, 0.0f };
    float d = sqrtf(x * x + y * y);
    if (d < EPS) { *ok = 0; return r; }   /* zero tangent rejected */
    *ok = 1;
    r.x = x / d;
    r.y = y / d;
    return r;
}

static Frame2 make_frame(float heading) {
    Frame2 fr;
    fr.heading = heading;
    fr.tan.x = cosf(heading);
    fr.tan.y = sinf(heading);
    fr.normal.x = -fr.tan.y;
    fr.normal.y = fr.tan.x;
    return fr;
}

/* parallel transport: rotate the frame by each wrapped tangent difference.
   A closed path of n points has n edges, the last one closing back onto p0;
   the turn across that closing vertex is stored as delta[0] (the seam). */
static int frame_transport(const Poly* po, int closed, FramePath* out) {
    memset(out, 0, sizeof(*out));
    int n = po->n;
    int edges = closed ? n : n - 1;
    if (n < 2 || edges < 1 || edges > MAX_FRAME) return -1;

    V2 t[MAX_FRAME];
    int ok;
    for (int i = 0; i < edges; i++) {
        V2 a = po->p[i];
        V2 b = closed ? po->p[(i + 1) % n] : po->p[i + 1];
        t[i] = v2_unit(b.x - a.x, b.y - a.y, &ok);
        if (!ok) return -2;                 /* zero tangent */
    }

    out->n = edges;
    out->closed = closed;
    out->f[0] = make_frame(atan2f(t[0].y, t[0].x));
    out->delta[0] = closed
        ? wrap_pi(atan2f(t[0].y, t[0].x) - atan2f(t[n - 1].y, t[n - 1].x))
        : 0.0f;
    for (int i = 1; i < edges; i++) {
        float da = wrap_pi(atan2f(t[i].y, t[i].x)
                           - atan2f(t[i - 1].y, t[i - 1].x));
        out->delta[i] = da;
        out->f[i] = make_frame(out->f[i - 1].heading + da);
    }
    return 0;
}

/* analytic winding of a closed polyline: total signed tangent turning / 2pi.
   The turning number theorem guarantees this is an integer. Independent code
   path from the transport accumulation above. */
static float analytic_winding(const Poly* po) {
    int n = po->n;
    float theta[MAX_FRAME];
    int ok;
    for (int i = 0; i < n; i++) {
        V2 a = po->p[i];
        V2 b = po->p[(i + 1) % n];
        V2 t = v2_unit(b.x - a.x, b.y - a.y, &ok);
        if (!ok) return 0.0f;
        theta[i] = atan2f(t.y, t.x);
    }
    float total = 0.0f;
    for (int i = 0; i < n; i++)
        total += wrap_pi(theta[(i + 1) % n] - theta[i]);
    return total / TWO_PI;
}

/* independent integer winding of point q by ray crossing */
static int point_winding(const Poly* po, V2 q) {
    int w = 0;
    for (int i = 0; i < po->n; i++) {
        V2 a = po->p[i], b = po->p[(i + 1) % po->n];
        if (a.y <= q.y) {
            if (b.y > q.y && (b.x - a.x) * (q.y - a.y)
                           - (q.x - a.x) * (b.y - a.y) > 0.0f) w++;
        } else {
            if (b.y <= q.y && (b.x - a.x) * (q.y - a.y)
                            - (q.x - a.x) * (b.y - a.y) < 0.0f) w--;
        }
    }
    return w;
}

/* spread the closure residual uniformly over the n steps and rebuild.
   Returns 0 on success, -5 if the spread introduces a jump, -6 if the
   closed orientation error is not zero. */
static int frame_apply_closure(FramePath* fp) {
    int n = fp->n;
    if (n < 2) return -1;

    float hol = 0.0f;
    for (int i = 0; i < n; i++) hol += fp->delta[i];
    fp->holonomy = hol;
    fp->winding = hol / TWO_PI;
    float w_int = roundf(fp->winding);
    fp->residual = hol - w_int * TWO_PI;    /* the "few degrees" gap */

    float step = fp->residual / n;
    if (fabsf(step) > SPREAD_JUMP_LIM) return -5;   /* jump after spread */

    /* corrected steps: steps 1..n-1 each pay an equal share of the residual;
       the seam step 0 absorbs the exact remainder so the full-loop turn is
       constructively 2pi*w and the closing orientation error is zero. */
    float dcorr[MAX_FRAME];
    float interior = 0.0f;
    for (int i = 1; i < n; i++) {
        dcorr[i] = fp->delta[i] - step;
        if (fabsf(dcorr[i]) >= PI - EPS) return -5; /* half-turn flip */
        interior += dcorr[i];
    }
    dcorr[0] = w_int * TWO_PI - interior;
    if (fabsf(dcorr[0]) >= PI - EPS) return -5;
    if (fabsf(dcorr[0] - fp->delta[0]) > SPREAD_JUMP_LIM) return -5;

    /* rebuild frames, keeping the seam frame f[0] fixed */
    for (int i = 0; i < n; i++) fp->delta[i] = dcorr[i];
    float h = fp->f[0].heading;
    fp->f[0] = make_frame(h);
    for (int i = 1; i < n; i++) {
        h += dcorr[i];
        fp->f[i] = make_frame(h);
    }

    /* full loop: corrected seam step must land back on f[0] exactly */
    float back = h + dcorr[0];
    if (fabsf(wrap_pi(back - fp->f[0].heading)) > CLOSE_EPS) return -6;
    /* corrected winding is now an integer turn */
    fp->winding = w_int;
    fp->residual = 0.0f;
    return 0;
}

/* holonomy implied by the raw transport, before any correction */
static float winding_pre_holonomy(const FramePath* fp) {
    float hol = 0.0f;
    for (int i = 0; i < fp->n; i++) hol += fp->delta[i];
    return hol / TWO_PI;
}

/* full closed-path pipeline: analytic cross-check, then uniform spread.
   -3: analytic winding is not an integer (degenerate geometry)
   -4: transported holonomy disagrees with the analytic winding
   -5/-6: propagated from frame_apply_closure */
static int frame_close(FramePath* fp, const Poly* po) {
    if (!fp->closed || fp->n != po->n) return -1;

    float wa = analytic_winding(po);
    if (fabsf(wa - roundf(wa)) > CLOSE_EPS) return -3;
    if (fabsf(winding_pre_holonomy(fp) - wa) > WIND_MATCH_TOL) return -4;
    return frame_apply_closure(fp);
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    /* unit square: four 90 degree turns, analytic winding +1 */
    Poly sq;
    poly_init(&sq);
    poly_add(&sq, 0, 0); poly_add(&sq, 1, 0);
    poly_add(&sq, 1, 1); poly_add(&sq, 0, 1);
    FramePath fs;
    check("square_transport", frame_transport(&sq, 1, &fs) == 0);
    check("square_frames", fs.n == 4 && fs.closed == 1);
    check("square_holonomy_pre", fabsf(winding_pre_holonomy(&fs) - 1.0f) < CLOSE_EPS);
    check("square_winding_analytic", fabsf(analytic_winding(&sq) - 1.0f) < CLOSE_EPS);
    check("square_winding_ray", point_winding(&sq, (V2){ 0.5f, 0.5f }) == 1);
    check("square_close", frame_close(&fs, &sq) == 0);
    check("square_closure_zero",
          fabsf(wrap_pi(fs.f[3].heading + fs.delta[0] - fs.f[0].heading)) < CLOSE_EPS);
    check("square_winding_int", fabsf(fs.winding - 1.0f) < EPS);
    check("square_residual_zero", fabsf(fs.residual) < EPS);
    printf("square_holonomy=%.6f\n", fs.holonomy);
    check("square_holonomy_2pi", fabsf(fs.holonomy - TWO_PI) < CLOSE_EPS);
    check("square_ortho",
          fabsf(fs.f[0].tan.x * fs.f[0].normal.x +
                fs.f[0].tan.y * fs.f[0].normal.y) < EPS &&
          fabsf(fs.f[0].normal.x + fs.f[0].tan.y) < EPS);
    check("square_unit_tan",
          fabsf(fs.f[2].tan.x * fs.f[2].tan.x +
                fs.f[2].tan.y * fs.f[2].tan.y - 1.0f) < EPS);

    /* sampled circle: the closure error must vanish after the spread */
    Poly cir;
    poly_init(&cir);
    for (int i = 0; i < 32; i++) {
        float a = TWO_PI * i / 32.0f;
        poly_add(&cir, cosf(a), sinf(a));
    }
    FramePath fc;
    check("circle_transport", frame_transport(&cir, 1, &fc) == 0);
    check("circle_winding_analytic", fabsf(analytic_winding(&cir) - 1.0f) < CLOSE_EPS);
    check("circle_winding_ray", point_winding(&cir, (V2){ 0.0f, 0.0f }) == 1);
    check("circle_close", frame_close(&fc, &cir) == 0);
    check("circle_closure_zero",
          fabsf(wrap_pi(fc.f[31].heading + fc.delta[0] - fc.f[0].heading)) < CLOSE_EPS);
    check("circle_winding_int", fabsf(fc.winding - 1.0f) < EPS);
    for (int i = 0; i < fc.n; i++) {
        float tl = sqrtf(fc.f[i].tan.x * fc.f[i].tan.x +
                         fc.f[i].tan.y * fc.f[i].tan.y);
        check("circle_unit_tan", fabsf(tl - 1.0f) < EPS);
        float ta = atan2f(fc.f[i].tan.y, fc.f[i].tan.x);
        check("circle_frame_aligned", fabsf(wrap_pi(fc.f[i].heading - ta)) < CLOSE_EPS);
    }

    /* zero tangent: a duplicated consecutive point must be rejected */
    Poly bad;
    poly_init(&bad);
    poly_add(&bad, 0, 0); poly_add(&bad, 1, 0); poly_add(&bad, 1, 0);
    poly_add(&bad, 0, 1);
    FramePath fb;
    check("zero_tan_rejected", frame_transport(&bad, 1, &fb) == -2);

    /* open path: transport only, no closure correction */
    Poly op;
    poly_init(&op);
    poly_add(&op, 0, 0); poly_add(&op, 1, 0); poly_add(&op, 1, 1);
    FramePath fo;
    check("open_transport", frame_transport(&op, 0, &fo) == 0 && fo.n == 2);
    check("open_turn", fabsf(wrap_pi(fo.f[1].heading - fo.f[0].heading)
                             - 0.5f * PI) < CLOSE_EPS);

    /* straight open path: no turning at all, no flip */
    Poly st;
    poly_init(&st);
    poly_add(&st, 0, 0); poly_add(&st, 1, 0); poly_add(&st, 2, 0);
    FramePath fst;
    frame_transport(&st, 0, &fst);
    check("straight_no_turn",
          fabsf(wrap_pi(fst.f[1].heading - fst.f[0].heading)) < EPS);

    /* injected integration drift: hand-built frames whose holonomy carries a
       0.05 rad residual; uniform spread must close the loop with zero error */
    FramePath drift;
    memset(&drift, 0, sizeof(drift));
    drift.n = 8;
    drift.closed = 1;
    for (int i = 0; i < 8; i++) drift.delta[i] = TWO_PI / 8.0f;
    drift.delta[3] += 0.05f;                       /* simulate accumulated drift */
    drift.f[0] = make_frame(0.0f);
    float hh = 0.0f;
    for (int i = 1; i < 8; i++) { hh += drift.delta[i]; drift.f[i] = make_frame(hh); }
    check("drift_detected", fabsf(winding_pre_holonomy(&drift) - 1.0f) > 1e-3f);
    check("drift_close", frame_apply_closure(&drift) == 0);
    check("drift_closure_zero",
          fabsf(wrap_pi(drift.f[7].heading + drift.delta[0] - drift.f[0].heading))
          < CLOSE_EPS);
    check("drift_winding_int", fabsf(drift.winding - 1.0f) < EPS);
    check("drift_residual_zero", fabsf(drift.residual) < EPS);

    /* jump rejection: a frame flipped by pi gives a half-turn after spread */
    FramePath flip;
    memset(&flip, 0, sizeof(flip));
    flip.n = 4;
    flip.closed = 1;
    flip.delta[0] = 0.0f; flip.delta[1] = 0.0f;
    flip.delta[2] = PI; flip.delta[3] = PI;        /* two flips still total 2pi */
    flip.f[0] = make_frame(0.0f);
    float fh = 0.0f;
    for (int i = 1; i < 4; i++) { fh += flip.delta[i]; flip.f[i] = make_frame(fh); }
    check("flip_rejected", frame_apply_closure(&flip) == -5);

    /* winding mismatch: tampered transport vs the real square geometry */
    FramePath tamper;
    frame_transport(&sq, 1, &tamper);
    tamper.delta[1] += 2.0f;                    /* 0.318 turn disagreement */
    check("winding_mismatch_rejected", frame_close(&tamper, &sq) == -4);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
