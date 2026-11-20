
/* foot: foot-contact solver for sloped terrain.
   Pipeline per frame: sole contact from foot facing -> pin to ground
   height -> back-solve ankle -> two-bone IK -> ankle twist from sole
   facing. Illegal ground queries and poses the bones cannot reach are
   rejected with explicit codes. */
#include <stdio.h>
#include <string.h>
#include <math.h>

#define MAX_GROUND 32
#define MAX_FRAMES 32
#define EPS 1e-6f

typedef struct { float x, y; } V2;

static float vlen(V2 v) { return sqrtf(v.x * v.x + v.y * v.y); }
static V2 vsub(V2 a, V2 b) { V2 r = { a.x - b.x, a.y - b.y }; return r; }
static V2 vnorm(V2 v) { float l = vlen(v); V2 r = { 0, 0 }; if (l > 0) { r.x = v.x / l; r.y = v.y / l; } return r; }
static int vfinite(V2 v) { return isfinite(v.x) && isfinite(v.y); }

static float wrap_pi(float a) {
    const float PI = 3.14159265f;
    const float TWO_PI = 6.2831853f;
    while (a > PI) a -= TWO_PI;
    while (a < -PI) a += TWO_PI;
    return a;
}

/* ---------------- ground: piecewise-linear 2D heightfield ---------------- */

typedef struct {
    float xs[MAX_GROUND];
    float hs[MAX_GROUND];   /* ground height y at each sampled x */
    int count;
} Ground;

/* validates: 2..MAX_GROUND samples, strictly increasing x, all finite */
static int ground_init(Ground* g, const float* xs, const float* hs, int n) {
    if (n < 2 || n > MAX_GROUND) return 0;
    for (int i = 0; i < n; i++) {
        if (!isfinite(xs[i]) || !isfinite(hs[i])) return 0;
        if (i > 0 && xs[i] <= xs[i - 1]) return 0;
    }
    memcpy(g->xs, xs, (size_t)n * sizeof(float));
    memcpy(g->hs, hs, (size_t)n * sizeof(float));
    g->count = n;
    return 1;
}

/* height h and slope (dy/dx) at x; returns 0 when x leaves the sampled range */
static int ground_sample(const Ground* g, float x, float* h, float* slope) {
    if (!isfinite(x) || g->count < 2) return 0;
    if (x < g->xs[0] || x > g->xs[g->count - 1]) return 0;
    int i = 0;
    while (i < g->count - 1 && x > g->xs[i + 1]) i++;
    float dx = g->xs[i + 1] - g->xs[i];
    float t = (x - g->xs[i]) / dx;
    *h = g->hs[i] + t * (g->hs[i + 1] - g->hs[i]);
    *slope = (g->hs[i + 1] - g->hs[i]) / dx;
    return 1;
}

/* sole angle that sits flush with the terrain at x; 0 where illegal */
static float ground_slope_angle(const Ground* g, float x) {
    float h, s;
    if (!ground_sample(g, x, &h, &s)) return 0.0f;
    return atanf(s);
}

/* ---------------- two-bone IK (mirrors the ik module) ---------------- */

typedef struct {
    int ok;
    V2 joint;
    float a1;   /* angle of first bone */
    float a2;   /* angle of second bone */
} IkResult;

static IkResult solve_two_bone(V2 root, V2 target, float l1, float l2, V2 pole) {
    IkResult r;
    memset(&r, 0, sizeof(r));
    V2 d = vsub(target, root);
    float dist = vlen(d);
    if (l1 <= 0 || l2 <= 0) return r;
    if (dist > l1 + l2 + 1e-5f) return r;                 /* unreachable */
    if (dist < fabsf(l1 - l2) - 1e-5f) return r;          /* too close */
    if (dist < 1e-6f) return r;

    float cos_a = (l1 * l1 + dist * dist - l2 * l2) / (2 * l1 * dist);
    if (cos_a > 1.0f) cos_a = 1.0f;
    if (cos_a < -1.0f) cos_a = -1.0f;
    float a = acosf(cos_a);

    V2 dir = vnorm(d);
    float side = (dir.x * pole.y - dir.y * pole.x) >= 0 ? 1.0f : -1.0f;
    float ang = atan2f(dir.y, dir.x) + side * a;

    r.joint.x = root.x + l1 * cosf(ang);
    r.joint.y = root.y + l1 * sinf(ang);
    r.a1 = ang;
    V2 d2 = vsub(target, r.joint);
    r.a2 = atan2f(d2.y, d2.x);
    r.ok = 1;
    return r;
}

/* ---------------- foot-contact solver ---------------- */

enum {
    FC_OK = 0,
    FC_BAD_INPUT = 1,        /* non-finite frame, non-positive bone length, ... */
    FC_GROUND_ILLEGAL = 2,   /* contact x outside ground samples / non-finite */
    FC_OUT_OF_REACH = 3      /* corrected ankle is beyond what the bones reach */
};

typedef struct {
    V2 hip;
    V2 ankle;        /* animation-desired ankle */
    float foot_ang;  /* desired sole world angle; use slope angle for flush plant */
} FootFrame;

typedef struct {
    V2 knee;
    V2 ankle;        /* corrected ankle the IK actually reaches */
    V2 contact;      /* corrected sole contact, pinned to the ground */
    float foot_ang;
    float shin_ang;  /* world angle of knee -> ankle */
    float twist;     /* local ankle rotation: wrap_pi(foot_ang - shin_ang) */
    float ground_err;/* |contact y - ground h| after correction */
    float fix_err;   /* contact y move: the raw float/sink against the ground */
} FootSolve;

static int foot_solve(const Ground* g, float l1, float l2, float foot_len,
                      V2 pole, const FootFrame* f, FootSolve* out) {
    if (l1 <= 0 || l2 <= 0 || foot_len < 0 || !isfinite(foot_len)) return FC_BAD_INPUT;
    if (!vfinite(f->hip) || !vfinite(f->ankle) || !isfinite(f->foot_ang) || !vfinite(pole))
        return FC_BAD_INPUT;

    /* 1. desired sole contact: ankle plus the foot along its facing */
    V2 dir = { cosf(f->foot_ang), sinf(f->foot_ang) };
    V2 contact0 = { f->ankle.x + foot_len * dir.x,
                    f->ankle.y + foot_len * dir.y };

    /* 2. ground height correction; out-of-range / non-finite is rejected */
    float h, slope;
    if (!ground_sample(g, contact0.x, &h, &slope)) return FC_GROUND_ILLEGAL;
    V2 contact = { contact0.x, h };

    /* 3. back-solve the ankle keeping the same facing */
    V2 ankle = { contact.x - foot_len * dir.x,
                 contact.y - foot_len * dir.y };

    /* 4. two-bone IK for the planted ankle; short bones are rejected */
    IkResult ik = solve_two_bone(f->hip, ankle, l1, l2, pole);
    if (!ik.ok) return FC_OUT_OF_REACH;

    /* 5. shin/ankle twist that makes the sole hold the requested facing */
    float shin_ang = atan2f(ankle.y - ik.joint.y, ankle.x - ik.joint.x);

    memset(out, 0, sizeof(*out));
    out->knee = ik.joint;
    out->ankle = ankle;
    out->contact = contact;
    out->foot_ang = f->foot_ang;
    out->shin_ang = shin_ang;
    out->twist = wrap_pi(f->foot_ang - shin_ang);
    out->ground_err = fabsf(contact.y - h);
    out->fix_err = fabsf(contact.y - contact0.y);
    return FC_OK;
}

/* -------- comparison against per-frame hand-placed feet -------- */

typedef struct {
    float max_fix_err;       /* solver contact vs manual contact (y move) */
    float max_manual_ground; /* manual contact vs ground: hand-placement error */
    float max_solver_ground; /* solver contact vs ground after pinning */
    int solved;
    int rejected;
} FootReport;

/* Runs the solver over hand-placed frames. Because correction only rewrites
   the contact y (x and facing are kept), each frame's fix equals its manual
   ground error, so max_manual_ground is an upper bound on max_fix_err.
   Returns 0 when any frame was rejected. */
static int foot_compare(const Ground* g, float l1, float l2, float foot_len,
                        V2 pole, const FootFrame* frames, int n, FootReport* rep) {
    memset(rep, 0, sizeof(*rep));
    if (n <= 0 || n > MAX_FRAMES) return 0;
    for (int i = 0; i < n; i++) {
        const FootFrame* f = &frames[i];
        V2 dir = { cosf(f->foot_ang), sinf(f->foot_ang) };
        V2 manual_contact = { f->ankle.x + foot_len * dir.x,
                              f->ankle.y + foot_len * dir.y };
        float h, slope;
        int legal = ground_sample(g, manual_contact.x, &h, &slope);

        FootSolve s;
        int rc = foot_solve(g, l1, l2, foot_len, pole, f, &s);
        if (rc != FC_OK || !legal) { rep->rejected++; continue; }

        float manual_ground = fabsf(manual_contact.y - h);
        float fix = fabsf(s.contact.y - manual_contact.y);
        if (manual_ground > rep->max_manual_ground) rep->max_manual_ground = manual_ground;
        if (fix > rep->max_fix_err) rep->max_fix_err = fix;
        if (s.ground_err > rep->max_solver_ground) rep->max_solver_ground = s.ground_err;
        rep->solved++;
    }
    return rep->rejected == 0 ? 1 : 0;
}

/* ---------------- self tests ---------------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    /* downhill: y = -0.4 x over [0, 10] */
    Ground slope_g;
    float sx[2] = { 0.0f, 10.0f };
    float sh[2] = { 0.0f, -4.0f };
    check("ground_init_valid", ground_init(&slope_g, sx, sh, 2) == 1);

    Ground bad;
    float one[1] = { 0.0f };
    check("reject_too_few", ground_init(&bad, one, one, 1) == 0);
    float ux[2] = { 1.0f, 1.0f }, uh[2] = { 0.0f, 1.0f };
    check("reject_unsorted", ground_init(&bad, ux, uh, 2) == 0);
    float nx[2] = { 0.0f, 1.0f }, nh[2] = { 0.0f, NAN };
    check("reject_non_finite", ground_init(&bad, nx, nh, 2) == 0);
    float big_x[MAX_GROUND + 1], big_h[MAX_GROUND + 1];
    for (int i = 0; i <= MAX_GROUND; i++) { big_x[i] = (float)i; big_h[i] = 0.0f; }
    check("reject_capacity", ground_init(&bad, big_x, big_h, MAX_GROUND + 1) == 0);

    float h, s;
    check("sample_height", ground_sample(&slope_g, 2.2f, &h, &s) == 1
          && fabsf(h + 0.88f) < 1e-5f);
    check("sample_slope", fabsf(s + 0.4f) < 1e-6f);
    check("sample_out_left", ground_sample(&slope_g, -0.1f, &h, &s) == 0);
    check("sample_out_right", ground_sample(&slope_g, 10.1f, &h, &s) == 0);
    check("sample_nan", ground_sample(&slope_g, NAN, &h, &s) == 0);

    /* flat ground for the floating-foot case */
    Ground flat;
    float fx[2] = { -2.0f, 8.0f }, fh[2] = { 0.0f, 0.0f };
    ground_init(&flat, fx, fh, 2);

    const float L1 = 1.0f, L2 = 1.0f, FL = 0.2f;
    V2 pole = { 1.0f, 0.0f };

    /* floating foot on flat ground: contact at y=-0.5 gets pulled to 0 */
    FootFrame fl = { { 0.0f, 1.68f }, { 0.0f, -0.5f }, 0.0f };
    FootSolve fs;
    check("flat_solved", foot_solve(&flat, L1, L2, FL, pole, &fl, &fs) == FC_OK);
    check("flat_pinned", fabsf(fs.contact.y) < EPS && fs.ground_err < EPS);
    check("flat_fix_amount", fabsf(fs.fix_err - 0.5f) < 1e-5f);

    /* downhill hand-placed frames: pelvis follows the slope, feet stay at the
       flat-ground height y=-1.6 -> the foot sinks early and floats late */
    FootFrame frames[3];
    for (int i = 0; i < 3; i++) {
        float x = 2.0f + 2.0f * i;                 /* 2, 4, 6 */
        frames[i].hip.x = x;
        frames[i].hip.y = -0.4f * x + 1.6f;
        frames[i].ankle.x = x;
        frames[i].ankle.y = -1.6f;
        frames[i].foot_ang = 0.0f;
    }
    FootSolve a, b, c;
    check("down_sink_solved", foot_solve(&slope_g, L1, L2, FL, pole, &frames[0], &a) == FC_OK);
    check("down_mid_solved", foot_solve(&slope_g, L1, L2, FL, pole, &frames[1], &b) == FC_OK);
    check("down_float_solved", foot_solve(&slope_g, L1, L2, FL, pole, &frames[2], &c) == FC_OK);

    check("all_on_ground", a.ground_err < EPS && b.ground_err < EPS && c.ground_err < EPS);
    check("sink_amount", fabsf(a.fix_err - 0.72f) < 1e-5f);   /* h(2.2) = -0.88 */
    check("float_amount", fabsf(c.fix_err - 0.88f) < 1e-5f);  /* h(6.2) = -2.48 */
    check("contact_x_kept", fabsf(a.contact.x - 2.2f) < EPS
          && fabsf(c.contact.x - 6.2f) < EPS);

    /* pole vector bends the knee forward of the hip->ankle line */
    check("knee_pole_side", a.knee.x > frames[0].hip.x);

    /* twist: sole facing minus shin world angle, wrapped to [-pi, pi] */
    check("twist_identity",
          fabsf(wrap_pi(a.twist + a.shin_ang - a.foot_ang)) < 1e-5f &&
          fabsf(wrap_pi(c.twist + c.shin_ang - c.foot_ang)) < 1e-5f);
    check("twist_in_range", fabsf(a.twist) <= 3.1415927f && fabsf(c.twist) <= 3.1415927f);

    /* wrap across +/-pi: foot facing backward while the shin points down */
    FootFrame back = { { 0.0f, 0.5f }, { 0.0f, -1.6f }, 3.14159265f };
    FootSolve bs;
    check("wrap_case_solved", foot_solve(&flat, L1, L2, FL, pole, &back, &bs) == FC_OK);
    check("wrap_pi_branch", fabsf(bs.twist) <= 3.1415927f
          && fabsf(wrap_pi(bs.twist + bs.shin_ang - bs.foot_ang)) < 1e-5f);

    /* sole flush with the slope: ankle (4,-1.6) lies on y=-0.4x, so the
       contact already sits on the terrain and no correction is needed */
    FootFrame flush = { { 4.0f, 0.0f }, { 4.0f, -1.6f }, ground_slope_angle(&slope_g, 4.0f) };
    FootSolve fls;
    check("slope_angle", fabsf(flush.foot_ang - atanf(-0.4f)) < 1e-6f);
    check("flush_solved", foot_solve(&slope_g, L1, L2, FL, pole, &flush, &fls) == FC_OK);
    check("flush_zero_fix", fls.fix_err < 1e-5f && fls.ground_err < EPS);
    check("flush_twist", fabsf(wrap_pi(fls.twist + fls.shin_ang - fls.foot_ang)) < 1e-5f);

    /* rejection: contact x outside the sampled range */
    FootFrame over = { { 9.0f, -2.0f }, { 9.95f, -1.6f }, 0.0f };
    check("reject_ground_right", foot_solve(&slope_g, L1, L2, FL, pole, &over, &fs) == FC_GROUND_ILLEGAL);
    FootFrame under = { { 0.0f, 1.0f }, { -0.3f, -0.5f }, 0.0f };
    check("reject_ground_left", foot_solve(&slope_g, L1, L2, FL, pole, &under, &fs) == FC_GROUND_ILLEGAL);

    /* rejection: corrected ankle beyond bone reach */
    FootFrame far = { { 5.0f, 10.0f }, { 5.0f, -1.6f }, 0.0f };
    check("reject_out_of_reach", foot_solve(&slope_g, L1, L2, FL, pole, &far, &fs) == FC_OUT_OF_REACH);

    /* rejection: malformed input */
    check("reject_bone_len", foot_solve(&flat, 0.0f, L2, FL, pole, &fl, &fs) == FC_BAD_INPUT);
    FootFrame nanf = { { 0.0f, 1.0f }, { NAN, -0.5f }, 0.0f };
    check("reject_nan_frame", foot_solve(&flat, L1, L2, FL, pole, &nanf, &fs) == FC_BAD_INPUT);
    check("reject_foot_len", foot_solve(&flat, L1, L2, -0.1f, pole, &fl, &fs) == FC_BAD_INPUT);

    /* comparison vs manual frames: manual ground error bounds the contact move,
       and every solved contact ends on the ground */
    FootReport rep;
    check("compare_all_solved", foot_compare(&slope_g, L1, L2, FL, pole, frames, 3, &rep) == 1);
    check("compare_solver_on_ground", rep.max_solver_ground < EPS);
    check("compare_bound", rep.max_fix_err <= rep.max_manual_ground + 1e-5f
          && rep.max_fix_err > 0.0f);
    check("compare_bound_value", fabsf(rep.max_manual_ground - 0.88f) < 1e-5f
          && fabsf(rep.max_fix_err - 0.88f) < 1e-5f);

    FootFrame with_bad[4];
    memcpy(with_bad, frames, sizeof(frames));
    with_bad[3] = over;
    FootReport rep2;
    check("compare_counts_reject", foot_compare(&slope_g, L1, L2, FL, pole, with_bad, 4, &rep2) == 0
          && rep2.solved == 3 && rep2.rejected == 1);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
