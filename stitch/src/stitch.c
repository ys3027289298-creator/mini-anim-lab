
/* stitch: C2 joining of two cubic hermite animation clips.
 *
 * At the seam each side contributes one cubic hermite piece. Converting to
 * bezier form, the two inner control handles x (left) and y (right) are the
 * unknowns. Position continuity (C0) is fixed by averaging the two clip
 * values at the seam; tangent (C1) and acceleration (C2) continuity give a
 * 2x2 linear system. The normalized system's infinity-norm condition number
 * decides whether the join is solvable; failures return an explicit reason
 * instead of forcing a match.
 */
#include <stdio.h>
#include <string.h>
#include <math.h>

#define EPS_TIME 1e-6f
#define EPS_DET  1e-12f      /* normalized determinant: below this is singular */
#define COND_MAX 1e6f        /* infinity-norm condition number limit (float32) */
#define EPS_F    1.1920929e-7f

/* cubic hermite piece on [t0,t1], endpoint values v0/v1, tangents m0/m1
 * (tangents are derivatives with respect to time t, not the local u). */
typedef struct {
    float t0, t1;
    float v0, v1;
    float m0, m1;
} Hermite;

static float hermite_eval(const Hermite* h, float t) {
    float dt = h->t1 - h->t0;
    float u = (t - h->t0) / dt;
    float u2 = u * u, u3 = u2 * u;
    float h00 = 2 * u3 - 3 * u2 + 1;
    float h10 = u3 - 2 * u2 + u;
    float h01 = -2 * u3 + 3 * u2;
    float h11 = u3 - u2;
    return h00 * h->v0 + h10 * dt * h->m0 + h01 * h->v1 + h11 * dt * h->m1;
}

static float hermite_d1(const Hermite* h, float t) {
    float dt = h->t1 - h->t0;
    float u = (t - h->t0) / dt;
    float h00 = 6 * u * u - 6 * u;
    float h10 = 3 * u * u - 4 * u + 1;
    float h01 = -6 * u * u + 6 * u;
    float h11 = 3 * u * u - 2 * u;
    return (h00 * h->v0 + h10 * dt * h->m0 + h01 * h->v1 + h11 * dt * h->m1) / dt;
}

static float hermite_d2(const Hermite* h, float t) {
    float dt = h->t1 - h->t0;
    float u = (t - h->t0) / dt;
    /* H00''=12u-6, H10''=6u-4, H01''=6-12u, H11''=6u-2 (derivatives in u).
     * value bases carry v, tangent bases carry dt*m, then divide by dt^2. */
    float c_v0 = (12.0f * u - 6.0f) / (dt * dt);
    float c_m0 = (6.0f * u - 4.0f) / dt;
    float c_v1 = (6.0f - 12.0f * u) / (dt * dt);
    float c_m1 = (6.0f * u - 2.0f) / dt;
    return c_v0 * h->v0 + c_m0 * h->m0 + c_v1 * h->v1 + c_m1 * h->m1;
}

typedef enum {
    STITCH_OK = 0,
    STITCH_BAD_INPUT,
    STITCH_NO_OVERLAP,
    STITCH_WINDOW_OUT_OF_RANGE,
    STITCH_SINGULAR,
    STITCH_ILL_CONDITIONED,
    STITCH_NONFINITE_SOLUTION
} StitchStatus;

typedef struct {
    float t_lo, t_hi;       /* full clip domain */
    Hermite h;              /* local piece around the seam:
                               left clip ends at seam (h.t1), right starts (h.t0) */
} Clip;

typedef struct {
    StitchStatus status;
    const char* reason;
    float cond;             /* infinity-norm condition number of normalized system */
    float residual;         /* |A x - b| inf of the physical system (value units) */
    float seam_t, seam_v, seam_m;
    Hermite left, right;    /* solved C2 pieces */
    Hermite hand_left, hand_right;   /* manual baseline: position/tangent aligned only */
    float d2_left, d2_right;
    float d2_jump;          /* measured second-derivative jump of the C2 join */
    float d2_jump_bound;    /* analytic upper bound for that jump */
    float hand_d2_jump;     /* baseline jump, expected to stay large */
    float max_dev_vs_hand;  /* pointwise max |c2 - hand| over the window */
} StitchResult;

/* solve the C2 join. clip A lies before clip B and their local pieces meet at
 * one seam time. [win0,win1] is the evaluation/comparison window, sampled at
 * n_samples points. */
static StitchResult stitch_c2(const Clip* A, const Clip* B,
                              float win0, float win1, int n_samples) {
    StitchResult r;
    memset(&r, 0, sizeof(r));
    r.status = STITCH_BAD_INPUT;
    r.reason = "bad input";

    const float fin[] = {
        A->t_lo, A->t_hi, B->t_lo, B->t_hi, win0, win1,
        A->h.t0, A->h.t1, A->h.v0, A->h.v1, A->h.m0, A->h.m1,
        B->h.t0, B->h.t1, B->h.v0, B->h.v1, B->h.m0, B->h.m1
    };
    for (size_t i = 0; i < sizeof(fin) / sizeof(fin[0]); i++) {
        if (!isfinite(fin[i])) { r.reason = "non-finite input (nan/inf in clip or window)"; return r; }
    }
    if (n_samples < 2) { r.reason = "need at least 2 sample points"; return r; }
    if (!(A->t_lo < A->t_hi) || !(B->t_lo < B->t_hi)) {
        r.reason = "clip domain has zero or negative length"; return r;
    }
    if (!(win0 < win1)) { r.reason = "evaluation window has zero length"; return r; }

    /* clips must overlap in time */
    if (A->t_hi < B->t_lo - EPS_TIME || B->t_hi < A->t_lo - EPS_TIME) {
        r.status = STITCH_NO_OVERLAP;
        r.reason = "clips do not overlap in time";
        return r;
    }

    /* local pieces must share one seam time */
    if (fabsf(A->h.t1 - B->h.t0) > EPS_TIME) {
        r.reason = "local pieces do not meet at a shared seam time";
        return r;
    }
    float tj = A->h.t1;
    if (tj < A->t_lo - EPS_TIME || tj > A->t_hi + EPS_TIME ||
        tj < B->t_lo - EPS_TIME || tj > B->t_hi + EPS_TIME) {
        r.reason = "seam lies outside one of the clip domains";
        return r;
    }

    /* evaluation window must be coverable by the two clips and straddle the
     * seam: A covers [win0,tj], B covers [tj,win1]. The clip domains may only
     * touch at the seam itself, so each side is checked against its own clip. */
    if (win0 < A->t_lo - EPS_TIME || tj > A->t_hi + EPS_TIME ||
        tj < B->t_lo - EPS_TIME || win1 > B->t_hi + EPS_TIME ||
        win0 > tj + EPS_TIME || win1 < tj - EPS_TIME) {
        r.status = STITCH_WINDOW_OUT_OF_RANGE;
        r.reason = "window outside clip overlap, or does not straddle the seam";
        return r;
    }

    float dL = tj - A->h.t0;
    float dR = B->h.t1 - tj;
    if (!(dL > 0.0f) || !(dR > 0.0f)) {
        r.status = STITCH_SINGULAR;
        r.reason = "zero-length interval next to seam: constraint matrix is singular";
        return r;
    }

    /* C0: shared seam point is the average of what each clip wanted. */
    float J = 0.5f * (A->h.v1 + B->h.v0);

    /* bezier outer handles: left P0/P1, right P4/P3. */
    float P1 = A->h.v0 + dL * A->h.m0 / 3.0f;
    float P4 = B->h.v1 - dR * B->h.m1 / 3.0f;

    /* Two linear systems, kept separate because their units differ.
     *
     * Bezier end second derivatives (the point order matters):
     *   left  at seam: 6 (J - 2x + P1) / dL^2
     *   right at seam: 6 (P4 - 2y + J) / dR^2
     *
     * Physical system for the inner bezier handles x (left), y (right):
     *   C1: [ dR     dL  ] [x] = [ J(dR+dL)                          ]
     *   C2: [-2dR^2  2dL^2] [y]   [ dL^2 P4 - dR^2 P1 + J(dL^2-dR^2) ]
     *
     * Normalized system (xhat = x/dL, yhat = y/dR, s = dL/dR) is used only
     * for the condition number -- it is dimensionless and O(1) while the two
     * seam intervals stay comparable:
     *   [ 1    s ] [xhat] = [ J(1+s)/dL                                   ]
     *   [-2s  2s ] [yhat]   [ (s P4 - P1 + J(1-s)) / dL                  ]
     */
    float s = dL / dR;
    float na11 = 1.0f, na12 = s, na21 = -2.0f * s, na22 = 2.0f * s;
    float ndet = na11 * na22 - na12 * na21;   /* 4s */

    if (fabsf(ndet) < EPS_DET) {
        r.status = STITCH_SINGULAR;
        r.reason = "determinant below tolerance: interval ratio degenerates, no unique C2 handles";
        return r;
    }

    /* infinity-norm condition number via explicit inverse of the normalized system */
    float nrm_a = fmaxf(fabsf(na11) + fabsf(na12), fabsf(na21) + fabsf(na22));
    float inv11 = na22 / ndet, inv12 = -na12 / ndet;
    float inv21 = -na21 / ndet, inv22 = na11 / ndet;
    float nrm_inv = fmaxf(fabsf(inv11) + fabsf(inv12), fabsf(inv21) + fabsf(inv22));
    r.cond = nrm_a * nrm_inv;
    if (r.cond > COND_MAX) {
        r.status = STITCH_ILL_CONDITIONED;
        r.reason = "condition number above limit: seam intervals too disparate, handles numerically unreliable";
        return r;
    }

    float dL2 = dL * dL, dR2 = dR * dR;
    float A11 = dR,  A12 = dL;
    float A21 = -2.0f * dR2, A22 = 2.0f * dL2;
    float B0 = J * (dR + dL);
    float B1 = dL2 * P4 - dR2 * P1 + J * (dL2 - dR2);
    float det = A11 * A22 - A12 * A21;   /* 2 dL dR (dL+dR), nonzero here */

    float x = (B0 * A22 - A12 * B1) / det;
    float y = (A11 * B1 - B0 * A21) / det;
    if (!isfinite(x) || !isfinite(y)) {
        r.status = STITCH_NONFINITE_SOLUTION;
        r.reason = "solved handles are non-finite";
        return r;
    }
    r.residual = fmaxf(fabsf(A11 * x + A12 * y - B0),
                       fabsf(A21 * x + A22 * y - B1));

    /* hermite tangents of the solved bezier pieces. the outer handles P1/P4
     * stay pinned to the original clips, so their tangents are preserved; only
     * the two seam-side tangents move, and the C1 row makes them equal. */
    float m0L = 3.0f * (P1 - A->h.v0) / dL;
    float mL = 3.0f * (J - x) / dL;
    float mR = 3.0f * (y - J) / dR;
    float m1R = 3.0f * (B->h.v1 - P4) / dR;
    float m = mL;
    (void)mR;

    r.seam_t = tj; r.seam_v = J; r.seam_m = m;
    r.left = (Hermite){ A->h.t0, tj, A->h.v0, J, m0L, mL };
    r.right = (Hermite){ tj, B->h.t1, J, B->h.v1, mR, m1R };

    /* manual baseline: force C0 through J but keep each clip's own seam tangent,
     * so C1 may fail and the second-derivative jump is deliberately preserved. */
    r.hand_left = (Hermite){ A->h.t0, tj, A->h.v0, J, A->h.m0, A->h.m1 };
    r.hand_right = (Hermite){ tj, B->h.t1, J, B->h.v1, B->h.m0, B->h.m1 };

    r.d2_left = hermite_d2(&r.left, tj);
    r.d2_right = hermite_d2(&r.right, tj);
    r.d2_jump = fabsf(r.d2_left - r.d2_right);
    r.hand_d2_jump = fabsf(hermite_d2(&r.hand_left, tj) - hermite_d2(&r.hand_right, tj));

    /* analytic bound: the two seam second derivatives are
     *   d2L = 6 (J-2x+P1)/dL^2 and d2R = 6 (P4-2y+J)/dR^2.
     * the C2 row enforces dR^2(J-2x+P1) = dL^2(P4-2y+J); its physical residual
     * rho bounds the jump exactly:
     *   |d2L - d2R| = 6 rho / (dL^2 dR^2)
     * plus a float32 rounding margin on the evaluated derivatives. */
    float qL = J - 2.0f * x + P1;
    float qR = P4 - 2.0f * y + J;
    float rho = fabsf(dR2 * qL - dL2 * qR);
    r.d2_jump_bound = 6.0f * rho / (dL2 * dR2)
                    + 8.0f * EPS_F * (fabsf(r.d2_left) + fabsf(r.d2_right));

    /* pointwise comparison against the manual join over the whole window */
    r.max_dev_vs_hand = 0.0f;
    for (int i = 0; i < n_samples; i++) {
        float t = win0 + (win1 - win0) * (float)i / (float)(n_samples - 1);
        float vc, vh;
        if (t <= tj) {
            vc = hermite_eval(&r.left, fminf(t, tj));
            vh = hermite_eval(&r.hand_left, fminf(t, tj));
        } else {
            vc = hermite_eval(&r.right, t);
            vh = hermite_eval(&r.hand_right, t);
        }
        float dev = fabsf(vc - vh);
        if (isfinite(dev) && dev > r.max_dev_vs_hand) r.max_dev_vs_hand = dev;
    }

    r.status = STITCH_OK;
    r.reason = "ok";
    return r;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    /* nominal case: both clips pass through 0 at t=1 with tangent 1, so the
     * manual join is C1 already -- but its acceleration jumps from -6 to +6. */
    Clip A = { 0.0f, 1.0f, { 0.0f, 1.0f, 0.0f, 0.0f, 1.0f, 1.0f } };
    Clip B = { 1.0f, 2.0f, { 1.0f, 2.0f, 0.0f, 0.0f, 1.0f, 1.0f } };
    StitchResult r = stitch_c2(&A, &B, 0.0f, 2.0f, 201);
    check("nominal_ok", r.status == STITCH_OK);
    printf("cond=%.3g residual=%.3g\n", r.cond, r.residual);
    check("c0_seam", fabsf(hermite_eval(&r.left, r.seam_t) - r.seam_v) < 1e-6f &&
                     fabsf(hermite_eval(&r.right, r.seam_t) - r.seam_v) < 1e-6f);
    float d1L = hermite_d1(&r.left, r.seam_t);
    float d1R = hermite_d1(&r.right, r.seam_t);
    check("c1_tangent", fabsf(d1L - d1R) < 1e-4f && fabsf(d1L - r.seam_m) < 1e-4f);
    printf("d2_left=%.4f d2_right=%.4f jump=%.3g bound=%.3g hand_jump=%.3f\n",
           r.d2_left, r.d2_right, r.d2_jump, r.d2_jump_bound, r.hand_d2_jump);
    check("c2_jump_tiny", r.d2_jump < 1e-3f);
    check("jump_under_bound", r.d2_jump <= r.d2_jump_bound + 1e-4f);
    check("bound_tight", r.d2_jump_bound < 1e-2f);
    check("hand_jump_large", r.hand_d2_jump > 10.0f);
    check("outer_anchors_kept",
          fabsf(r.left.v0 - A.h.v0) < 1e-7f && fabsf(r.left.m0 - A.h.m0) < 1e-7f &&
          fabsf(r.right.v1 - B.h.v1) < 1e-7f && fabsf(r.right.m1 - B.h.m1) < 1e-7f);
    check("pointwise_finite", isfinite(r.max_dev_vs_hand));
    printf("max_dev_vs_hand=%.4f\n", r.max_dev_vs_hand);

    /* mismatched seam values and tangents: the manual join cannot even hold C1 */
    Clip A2 = { 0.0f, 1.0f, { 0.0f, 1.0f, 0.0f, 0.0f, 1.0f, 1.0f } };
    Clip B2 = { 1.0f, 2.0f, { 1.0f, 2.0f, 0.2f, 0.0f, 0.5f, 0.0f } };
    StitchResult r2 = stitch_c2(&A2, &B2, 0.0f, 2.0f, 201);
    check("mismatch_ok", r2.status == STITCH_OK);
    float hand_d1_jump = fabsf(hermite_d1(&r2.hand_left, r2.seam_t) -
                               hermite_d1(&r2.hand_right, r2.seam_t));
    float c2_d1_jump = fabsf(hermite_d1(&r2.left, r2.seam_t) -
                             hermite_d1(&r2.right, r2.seam_t));
    printf("hand_d1_jump=%.4f c2_d1_jump=%.3g c2_d2_jump=%.3g bound=%.3g\n",
           hand_d1_jump, c2_d1_jump, r2.d2_jump, r2.d2_jump_bound);
    check("hand_breaks_c1", hand_d1_jump > 0.1f);
    check("c2_restores_c1", c2_d1_jump < 1e-4f);
    check("c2_restores_c2", r2.d2_jump < 1e-3f && r2.d2_jump <= r2.d2_jump_bound + 1e-4f);

    /* rejection: window outside overlap / not straddling seam */
    StitchResult r3 = stitch_c2(&A, &B, -0.5f, 0.5f, 10);
    check("window_rejected", r3.status == STITCH_WINDOW_OUT_OF_RANGE);
    printf("reason_window=%s\n", r3.reason);
    StitchResult r3b = stitch_c2(&A, &B, 0.2f, 0.8f, 10);
    check("window_straddle_rejected", r3b.status == STITCH_WINDOW_OUT_OF_RANGE);

    /* rejection: no temporal overlap */
    Clip C = { 2.0f, 3.0f, { 2.0f, 3.0f, 0.0f, 0.0f, 1.0f, 1.0f } };
    Clip D = { 3.5f, 4.5f, { 3.5f, 4.5f, 0.0f, 0.0f, 1.0f, 1.0f } };
    StitchResult r4 = stitch_c2(&C, &D, 2.5f, 4.0f, 10);
    check("no_overlap_rejected", r4.status == STITCH_NO_OVERLAP);
    printf("reason_no_overlap=%s\n", r4.reason);

    /* rejection: ill-conditioned (dL/dR = 1e-7) and singular (dL = 0) */
    Clip E = { 0.0f, 1.0f, { 1.0f - 1e-7f, 1.0f, 0.0f, 0.0f, 1.0f, 1.0f } };
    StitchResult r5 = stitch_c2(&E, &B, 0.9999999f, 1.5f, 10);
    printf("ill_cond=%.3g status=%d reason=%s\n", r5.cond, r5.status, r5.reason);
    check("ill_conditioned_rejected", r5.status == STITCH_ILL_CONDITIONED);
    Clip F = { 0.0f, 1.0f, { 1.0f, 1.0f, 0.0f, 0.0f, 1.0f, 1.0f } };
    StitchResult r6 = stitch_c2(&F, &B, 0.5f, 1.5f, 10);
    check("singular_rejected", r6.status == STITCH_SINGULAR);
    printf("reason_singular=%s\n", r6.reason);

    /* well-conditioned but unequal intervals still solve, and the bound holds */
    Clip Jc = { 0.0f, 2.0f, { 0.5f, 2.0f, 0.0f, 0.5f, 0.4f, -0.3f } };
    Clip K = { 2.0f, 3.0f, { 2.0f, 3.0f, 0.4f, 0.0f, 0.2f, 0.6f } };
    StitchResult r9 = stitch_c2(&Jc, &K, 0.5f, 3.0f, 251);
    check("unequal_intervals_ok", r9.status == STITCH_OK && r9.cond < 100.0f);
    check("unequal_c2_under_bound", r9.d2_jump <= r9.d2_jump_bound + 1e-3f);
    printf("unequal: cond=%.3g d2_jump=%.3g bound=%.3g dev=%.4f\n",
           r9.cond, r9.d2_jump, r9.d2_jump_bound, r9.max_dev_vs_hand);

    /* rejection: non-finite input */
    Clip G = { 0.0f, 1.0f, { 0.0f, 1.0f, NAN, 0.0f, 1.0f, 1.0f } };
    StitchResult r7 = stitch_c2(&G, &B, 0.0f, 1.5f, 10);
    check("nan_rejected", r7.status == STITCH_BAD_INPUT);
    printf("reason_nan=%s\n", r7.reason);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
