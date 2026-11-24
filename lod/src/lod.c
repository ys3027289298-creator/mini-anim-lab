/* lod: skeleton LOD handoff with a cost-gated, phase-aligned transition.

   Symptom: after the LOD swap the distal limb joints pop.  The near and far
   skeletons play separately authored clips - the far clip starts later and
   bends less - and nothing bridges them.  The swap fires the instant the
   camera distance crosses a threshold, so the far clip's joints land
   wherever they happen to be, and the next frame they snap onto the far
   clip's motion.

   Fix, in three parts:
     1. the distance band [DIST_NEAR, DIST_FAR] is the switching interval.
        Entering it arms a swap in the entry direction (+1 receding, -1
        approaching).  The swap is committed by a cost threshold while still
        inside the band; a cheap frame (similar clips, calm motion) fires it,
        a busy frame waits.  If the band is left before one is found, the
        swap is forced and flagged, so the LOD never gets stuck.
     2. the transition window (0.2 s) starts on the trigger frame and runs
        along the entry direction.  The incoming clip's matched key points
        are phase-aligned and committed at that same frame, so the window
        opens at the best overlap the two clips can reach.
     3. near and far each record their own trigger time, handoff error and
        worst displayed step inside their window.

   Timing is a pure function of the frame index (frame * DT): there is no
   accumulated clock, no randomness and no wall-clock input, so repeated
   runs trigger on exactly the same frame.  The measured deviation is 0
   frames, inside the half-frame budget. */
#include <stdio.h>
#include <string.h>
#include <math.h>

#define TAU 6.283185307179586f
#define DT (1.0f / 60.0f)
#define NF 360                      /* 6 s at 60 fps */

#define NKEY 3                      /* matched key points per rig */
#define MAX_TRIG 8
#define ALIGN_SUB 4                 /* alignment step DT/4 = 1/240 s */
#define ALIGN_SPAN 0.40f            /* +-0.4 s search around the base offset */
#define TRANSITION_FRAMES 12        /* 0.2 s window */
#define COST_THRESHOLD 0.85f        /* in units of natural per-frame motion */
#define DIST_FAR 5.40f
#define DIST_NEAR 3.60f
#define HALF_FRAME (0.5f * DT)

typedef struct { float x, y; } V2;

/* ------------------------------------------------------------------ rigs */

typedef struct {
    int bones;
    float len;
    int key_joint[NKEY];            /* joint index of each matched key point */
} Rig;

/* near: 6 x 0.6 units; far: 3 x 1.2 units.  Same reach, different joint
   count and segment size, so the matched key points sit on different joints
   and each rig gets its own recorded error. */
static const Rig RIGS[2] = {
    { 6, 0.60f, { 2, 4, 6 } },
    { 3, 1.20f, { 1, 2, 3 } },
};

static float bone_angle_near(int i, float t) {
    return 0.55f * sinf(TAU * 0.90f * t + 0.35f * (float)i)
         + 0.12f * sinf(TAU * 1.60f * t + 0.80f * (float)i);
}

/* far clip: even sub-bone phases (coarser bend), 0.35 s late, 7% lower
   amplitude.  The lag is what key-point alignment removes; the amplitude
   mismatch is the residual that no time shift can take out. */
#define FAR_LAG 0.35f
#define FAR_AMP 0.93f

static float bone_angle_far(int j, float t) {
    float tt = t - FAR_LAG;
    int i = 2 * j;
    return FAR_AMP * (0.55f * sinf(TAU * 0.90f * tt + 0.35f * (float)i)
                    + 0.12f * sinf(TAU * 1.60f * tt + 0.80f * (float)i));
}

static void rig_keys(int s, float clock, V2 out[NKEY]) {
    const Rig* r = &RIGS[s];
    V2 j[8];
    memset(j, 0, sizeof(j));
    float ang = 0.0f;
    for (int i = 0; i < r->bones; i++) {
        ang += s ? bone_angle_far(i, clock) : bone_angle_near(i, clock);
        j[i + 1].x = j[i].x + r->len * cosf(ang);
        j[i + 1].y = j[i].y + r->len * sinf(ang);
    }
    for (int m = 0; m < NKEY; m++) out[m] = j[r->key_joint[m]];
}

static float d2(V2 a, V2 b) {
    float dx = a.x - b.x, dy = a.y - b.y;
    return dx * dx + dy * dy;
}

static float max_gap(const V2 a[NKEY], const V2 b[NKEY]) {
    float w = 0.0f;
    for (int m = 0; m < NKEY; m++) {
        float e = sqrtf(d2(a[m], b[m]));
        if (e > w) w = e;
    }
    return w;
}

/* worst movement of any key point in one frame of rig s */
static float natural_step(int s, float clock) {
    V2 a[NKEY], b[NKEY];
    rig_keys(s, clock, a);
    rig_keys(s, clock - DT, b);
    return max_gap(a, b);
}

/* -------------------------------------------------------------- scenario */

static float smooth01(float x) {
    if (x <= 0.0f) return 0.0f;
    if (x >= 1.0f) return 1.0f;
    return x * x * (3.0f - 2.0f * x);
}

/* camera pulls back 2.5 -> 7.0, holds, then comes back. */
static float camera_distance(float t) {
    if (t <= 1.0f) return 2.5f;
    if (t <= 3.0f) return 2.5f + 4.5f * smooth01((t - 1.0f) / 2.0f);
    if (t <= 4.0f) return 7.0f;
    if (t <= 6.0f) return 7.0f - 4.5f * smooth01((t - 4.0f) / 2.0f);
    return 2.5f;
}

/* alignment residual: position at the junction plus a quarter weight on the
   position one frame later (that is the direction term).  Minimax over the
   key points, so no single key is left behind. */
static float align_resid(int out_rig, int in_rig, float t, const float off[2],
                         float cand) {
    V2 o0[NKEY], o1[NKEY], i0[NKEY], i1[NKEY];
    rig_keys(out_rig, t + off[out_rig], o0);
    rig_keys(out_rig, t + DT + off[out_rig], o1);
    rig_keys(in_rig, t + cand, i0);
    rig_keys(in_rig, t + DT + cand, i1);
    float worst = 0.0f;
    for (int m = 0; m < NKEY; m++) {
        float e = sqrtf(d2(o0[m], i0[m])) + 0.25f * sqrtf(d2(o1[m], i1[m]));
        if (e > worst) worst = e;
    }
    return worst;
}

/* Best offset for the incoming clip.  Deterministic scan, earliest wins on
   ties.  Committed only if no key point gets worse than at the current
   offset (Pareto guard), so "aligned" is never a regression. */
static float align_offset(int out_rig, int in_rig, float t, const float off[2]) {
    int steps = (int)(ALIGN_SPAN / DT * (float)ALIGN_SUB);
    float base = off[in_rig];
    float best = base, best_err = 1e30f;
    for (int k = -steps; k <= steps; k++) {
        float cand = base + (float)k * DT / (float)ALIGN_SUB;
        float e = align_resid(out_rig, in_rig, t, off, cand);
        if (e < best_err - 1e-9f) { best_err = e; best = cand; }
    }
    if (best == base) return base;

    V2 o[NKEY], raw[NKEY], al[NKEY];
    rig_keys(out_rig, t + off[out_rig], o);
    rig_keys(in_rig, t + base, raw);
    rig_keys(in_rig, t + best, al);
    for (int m = 0; m < NKEY; m++)
        if (sqrtf(d2(o[m], al[m])) > sqrtf(d2(o[m], raw[m])) + 1e-6f)
            return base;
    return best;
}

/* Cost of handing off on this frame, in units of natural per-frame motion:
   the matched-key gap that alignment can actually achieve, divided by the
   faster clip's own step.  A busy frame is expensive even when the clips are
   close, so the gate waits for a quiet one. */
static float handoff_cost(int out_rig, int in_rig, float t, const float off[2]) {
    float cand = align_offset(out_rig, in_rig, t, off);
    V2 o[NKEY], i[NKEY];
    rig_keys(out_rig, t + off[out_rig], o);
    rig_keys(in_rig, t + cand, i);
    float nat = natural_step(out_rig, t + off[out_rig]);
    float nf = natural_step(in_rig, t + cand);
    if (nf > nat) nat = nf;
    if (nat < 1e-4f) nat = 1e-4f;
    return max_gap(o, i) / nat;
}

/* ------------------------------------------------------------------ run */

typedef struct {
    int out_joint, in_joint;
    V2 out_pt, in_raw, in_aligned;
    float res_raw, res_aligned;     /* per-key position residual */
} KeyMark;

typedef struct {
    int frame, from, to;
    int gated;                      /* 1 cost threshold, 0 band edge */
    float time, dist, cost, entry;  /* entry +1 receding, -1 approaching */
    int waited;                     /* frames armed before firing */
    float align_base, align_off;
    float gap_raw, gap_aligned;
    float jump_at, window_peak;
    float nat_peak_window, window_excess;
    KeyMark key[NKEY];
} Trigger;

typedef struct {
    int frame;
    float dist, cost;
    int state, window;
} FrameInfo;

typedef struct {
    FrameInfo f[NF];
    float jump[NF];
    Trigger trig[MAX_TRIG];
    int ntrig;
    /* per skeleton size: 0 = near, 1 = far */
    int attach_frame[2], detach_frame[2];
    float attach_time[2], detach_time[2];
    float attach_err[2];            /* aligned handoff residual it absorbs */
    float detach_err[2];            /* same residual, left behind */
    float attach_step[2];           /* displayed step on its attach frame */
    float window_step[2];           /* largest displayed step in its window */
    float natural_peak[2];
    float jump_peak;
    int forced;
} Run;

typedef enum { POLICY_BAND, POLICY_CROSS } Policy;

static void simulate(Policy policy, float threshold, int window, int use_align,
                     Run* R) {
    memset(R, 0, sizeof(*R));
    for (int s = 0; s < 2; s++) {
        R->attach_frame[s] = R->detach_frame[s] = -1;
        R->attach_err[s] = R->detach_err[s] = -1.0f;
        R->attach_step[s] = R->window_step[s] = -1.0f;
    }

    int state = 0;
    int blending = 0, from = 0, to = 0, start = 0;
    int pending = -1, pending_since = 0;
    float off[2] = { 0.0f, 0.0f };
    V2 prev[NKEY];
    int have_prev = 0;

    for (int n = 0; n < NF; n++) {
        float t = (float)n * DT;
        float dist = camera_distance(t);

        R->f[n].frame = n;
        R->f[n].dist = dist;
        R->f[n].state = state;
        R->f[n].window = blending;

        /* which rig does this frame want? */
        int want = state;
        if (policy == POLICY_BAND) {
            if (state == 0 && dist > DIST_NEAR) want = 1;   /* arm receding */
            else if (state == 1 && dist < DIST_FAR) want = 0; /* arm approaching */
        } else {
            if (dist >= DIST_FAR) want = 1;   /* the old behaviour */
            else if (dist <= DIST_NEAR) want = 0;
        }
        int in_band = (dist >= DIST_NEAR && dist <= DIST_FAR);

        float cost = 0.0f;
        if (!blending && want != state) {
            if (pending != want) { pending = want; pending_since = n; }
            cost = handoff_cost(state, want, t, off);

            int gate_ok = 0, must = 0;
            if (policy == POLICY_BAND && in_band) {
                gate_ok = (threshold >= 1e8f) ? 1 : (cost <= threshold);
            } else {
                must = 1;               /* band edge or old crossing policy */
            }

            if (gate_ok || must) {
                from = state;
                to = want;
                start = n;
                blending = 1;

                float base = off[to];
                V2 out_pre[NKEY], in_pre[NKEY];
                rig_keys(from, t + off[from], out_pre);
                rig_keys(to, t + base, in_pre);

                Trigger* g = &R->trig[R->ntrig++];
                g->frame = n;
                g->from = from;
                g->to = to;
                g->gated = gate_ok && !must;
                g->time = t;
                g->dist = dist;
                g->cost = cost;
                g->entry = (to == 1) ? 1.0f : -1.0f;
                g->waited = n - pending_since;
                g->align_base = base;
                g->gap_raw = max_gap(out_pre, in_pre);
                g->nat_peak_window = 0.0f;
                g->window_excess = 0.0f;

                float cand = use_align ? align_offset(from, to, t, off) : base;
                off[to] = cand;
                g->align_off = cand;

                V2 out_al[NKEY], in_al[NKEY];
                rig_keys(from, t + off[from], out_al);
                rig_keys(to, t + cand, in_al);
                g->gap_aligned = max_gap(out_al, in_al);

                for (int m = 0; m < NKEY; m++) {
                    KeyMark* k = &g->key[m];
                    k->out_joint = RIGS[from].key_joint[m];
                    k->in_joint = RIGS[to].key_joint[m];
                    k->out_pt = out_pre[m];
                    k->in_raw = in_pre[m];
                    k->in_aligned = in_al[m];
                    k->res_raw = sqrtf(d2(out_pre[m], in_pre[m]));
                    k->res_aligned = sqrtf(d2(out_al[m], in_al[m]));
                }

                if (must) R->forced = 1;
                pending = -1;

                R->attach_frame[to] = n;
                R->attach_time[to] = t;
                R->detach_frame[from] = n;
                R->detach_time[from] = t;
                R->attach_err[to] = g->gap_aligned;
                R->detach_err[from] = g->gap_aligned;
            }
        }
        R->f[n].cost = cost;

        V2 disp[NKEY];
        if (blending) {
            float u = (window > 0) ? (float)(n - start) / (float)window : 1.0f;
            if (u > 1.0f) u = 1.0f;
            float w = u * u * (3.0f - 2.0f * u);
            V2 a[NKEY], b[NKEY];
            rig_keys(from, t + off[from], a);
            rig_keys(to, t + off[to], b);
            for (int m = 0; m < NKEY; m++) {
                disp[m].x = a[m].x + (b[m].x - a[m].x) * w;
                disp[m].y = a[m].y + (b[m].y - a[m].y) * w;
            }
            if (n - start >= window) { state = to; blending = 0; }
        } else {
            rig_keys(state, t + off[state], disp);
        }

        float jump = have_prev ? max_gap(disp, prev) : 0.0f;
        R->jump[n] = jump;
        if (jump > R->jump_peak) R->jump_peak = jump;

        if (blending) {
            float nat0 = natural_step(0, t + off[0]);
            float nat1 = natural_step(1, t + off[1]);
            float natmax = nat0 > nat1 ? nat0 : nat1;
            Trigger* cur = &R->trig[R->ntrig - 1];
            float ex = jump - natmax;
            if (ex > cur->window_excess) cur->window_excess = ex;
            if (natmax > cur->nat_peak_window) cur->nat_peak_window = natmax;
            if (jump > R->window_step[to]) R->window_step[to] = jump;
            if (n == R->attach_frame[to]) R->attach_step[to] = jump;
            if (nat0 > R->natural_peak[0]) R->natural_peak[0] = nat0;
            if (nat1 > R->natural_peak[1]) R->natural_peak[1] = nat1;
        } else if (n == R->attach_frame[state]) {
            R->attach_step[state] = jump;
        }

        for (int m = 0; m < NKEY; m++) prev[m] = disp[m];
        have_prev = 1;
    }

    for (int i = 0; i < R->ntrig; i++) {
        int n0 = R->trig[i].frame;
        R->trig[i].jump_at = R->jump[n0];
        float peak = R->jump[n0];
        for (int n = n0; n <= n0 + window && n < NF; n++)
            if (R->jump[n] > peak) peak = R->jump[n];
        R->trig[i].window_peak = peak;
    }
}

/* ------------------------------------------------------------- reporting */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static void report_run(const char* tag, const Run* R) {
    printf("%s_ntrig=%d\n", tag, R->ntrig);
    for (int i = 0; i < R->ntrig; i++) {
        const Trigger* g = &R->trig[i];
        printf("%s_trig%d_frame=%d\n", tag, i, g->frame);
        printf("%s_trig%d_time=%.6f\n", tag, i, (double)g->time);
        printf("%s_trig%d_from=%d\n", tag, i, g->from);
        printf("%s_trig%d_to=%d\n", tag, i, g->to);
        printf("%s_trig%d_entry=%+.0f\n", tag, i, (double)g->entry);
        printf("%s_trig%d_dist=%.4f\n", tag, i, (double)g->dist);
        printf("%s_trig%d_cost=%.6f\n", tag, i, (double)g->cost);
        printf("%s_trig%d_gated=%d\n", tag, i, g->gated);
        printf("%s_trig%d_waited=%d\n", tag, i, g->waited);
        printf("%s_trig%d_align_base=%.6f\n", tag, i, (double)g->align_base);
        printf("%s_trig%d_align_off=%.6f\n", tag, i, (double)g->align_off);
        printf("%s_trig%d_gap_raw=%.6f\n", tag, i, (double)g->gap_raw);
        printf("%s_trig%d_gap_aligned=%.6f\n", tag, i, (double)g->gap_aligned);
        printf("%s_trig%d_jump_at=%.6f\n", tag, i, (double)g->jump_at);
        printf("%s_trig%d_jump_window=%.6f\n", tag, i, (double)g->window_peak);
        printf("%s_trig%d_nat_window=%.6f\n", tag, i, (double)g->nat_peak_window);
        printf("%s_trig%d_window_excess=%.6f\n", tag, i, (double)g->window_excess);
        printf("%s_trig%d_window_sec=%.6f\n", tag, i,
               (double)(TRANSITION_FRAMES * DT));
        for (int m = 0; m < NKEY; m++) {
            const KeyMark* k = &g->key[m];
            printf("%s_trig%d_key%d_out_joint=%d\n", tag, i, m, k->out_joint);
            printf("%s_trig%d_key%d_in_joint=%d\n", tag, i, m, k->in_joint);
            printf("%s_trig%d_key%d_res_raw=%.6f\n", tag, i, m,
                   (double)k->res_raw);
            printf("%s_trig%d_key%d_res_aligned=%.6f\n", tag, i, m,
                   (double)k->res_aligned);
            printf("%s_trig%d_key%d_out_pt=%.4f,%.4f\n", tag, i, m,
                   (double)k->out_pt.x, (double)k->out_pt.y);
            printf("%s_trig%d_key%d_in_pt=%.4f,%.4f\n", tag, i, m,
                   (double)k->in_aligned.x, (double)k->in_aligned.y);
        }
    }
    printf("%s_forced=%d\n", tag, R->forced);
    for (int s = 0; s < 2; s++) {
        const char* rig = s ? "far" : "near";
        if (R->attach_frame[s] >= 0) {
            printf("%s_%s_on_frame=%d\n", tag, rig, R->attach_frame[s]);
            printf("%s_%s_on_time=%.6f\n", tag, rig, (double)R->attach_time[s]);
        }
        if (R->detach_frame[s] >= 0) {
            printf("%s_%s_off_frame=%d\n", tag, rig, R->detach_frame[s]);
            printf("%s_%s_off_time=%.6f\n", tag, rig, (double)R->detach_time[s]);
        }
        printf("%s_%s_handoff_err=%.6f\n", tag, rig, (double)R->attach_err[s]);
        printf("%s_%s_handoff_err_out=%.6f\n", tag, rig, (double)R->detach_err[s]);
        printf("%s_%s_attach_step=%.6f\n", tag, rig, (double)R->attach_step[s]);
        printf("%s_%s_window_step=%.6f\n", tag, rig, (double)R->window_step[s]);
        printf("%s_%s_natural_peak=%.6f\n", tag, rig, (double)R->natural_peak[s]);
    }
    printf("%s_jump_peak=%.6f\n", tag, (double)R->jump_peak);
}

int main(void) {
    Run fixed, cross, ungated, noalign, forced, r2, r3;
    simulate(POLICY_BAND, COST_THRESHOLD, TRANSITION_FRAMES, 1, &fixed);
    simulate(POLICY_CROSS, COST_THRESHOLD, 0, 0, &cross);   /* the bug */
    simulate(POLICY_BAND, 1e9f, TRANSITION_FRAMES, 1, &ungated);
    simulate(POLICY_BAND, COST_THRESHOLD, TRANSITION_FRAMES, 0, &noalign);
    simulate(POLICY_BAND, 0.0f, TRANSITION_FRAMES, 1, &forced);

    report_run("fixed", &fixed);
    report_run("cross", &cross);
    report_run("ungated", &ungated);
    report_run("noalign", &noalign);
    report_run("forced", &forced);

    printf("dist_far=%.4f\n", (double)DIST_FAR);
    printf("dist_near=%.4f\n", (double)DIST_NEAR);
    printf("threshold=%.4f\n", (double)COST_THRESHOLD);
    printf("window_sec=%.4f\n", (double)(TRANSITION_FRAMES * DT));
    printf("dt=%.6f\n", (double)DT);
    printf("half_frame=%.6f\n", (double)HALF_FRAME);
    printf("far_lag=%.4f\n", (double)FAR_LAG);

    check("rig_reach_equal",
          fabsf(RIGS[0].bones * RIGS[0].len -
                RIGS[1].bones * RIGS[1].len) < 1e-5f);
    check("rig_joint_count_diff", RIGS[0].bones != RIGS[1].bones);

    /* (1) the cost threshold takes over inside the switching interval */
    check("two_swaps", fixed.ntrig == 2);
    check("both_gated", fixed.trig[0].gated == 1 && fixed.trig[1].gated == 1);
    check("not_forced", fixed.forced == 0);
    check("inside_band",
          fixed.trig[0].dist >= DIST_NEAR && fixed.trig[0].dist <= DIST_FAR &&
          fixed.trig[1].dist >= DIST_NEAR && fixed.trig[1].dist <= DIST_FAR);
    check("cost_below_threshold",
          fixed.trig[0].cost <= COST_THRESHOLD &&
          fixed.trig[1].cost <= COST_THRESHOLD);
    check("cross_policy_ungated",
          cross.trig[0].gated == 0 && cross.trig[1].gated == 0);
    check("gate_delays_outbound",
          fixed.trig[0].frame > ungated.trig[0].frame &&
          fixed.trig[0].cost <= ungated.trig[0].cost);
    check("gate_delays_inbound",
          fixed.trig[1].frame > ungated.trig[1].frame &&
          fixed.trig[1].cost <= ungated.trig[1].cost);
    check("entry_direction",
          fixed.trig[0].entry > 0.0f && fixed.trig[0].to == 1 &&
          fixed.trig[1].entry < 0.0f && fixed.trig[1].to == 0);
    check("gate_waits_for_cheap_frame",
          fixed.trig[0].waited >= 1 && fixed.trig[1].waited >= 1);
    /* safety valve: an unreachable threshold swaps at the band edge */
    check("valve_fires_when_no_cheap_frame",
          forced.ntrig == 2 && forced.forced == 1 && forced.trig[0].gated == 0);
    check("valve_fires_at_band_edge",
          fabsf(forced.trig[0].dist - DIST_FAR) < 0.05f &&
          fabsf(forced.trig[1].dist - DIST_NEAR) < 0.05f);

    /* (2) transition window along the entry direction, key points aligned */
    check("window_length", fabsf(TRANSITION_FRAMES * DT - 0.2f) < 1e-6f);
    check("window_starts_at_trigger",
          fixed.attach_frame[1] == fixed.trig[0].frame &&
          fixed.attach_frame[0] == fixed.trig[1].frame);
    for (int i = 0; i < fixed.ntrig; i++)
        for (int m = 0; m < NKEY; m++)
            check("key_aligned_no_worse",
                  fixed.trig[i].key[m].res_aligned <=
                  fixed.trig[i].key[m].res_raw + 1e-6f);
    check("align_reduces_gap",
          fixed.trig[0].gap_aligned < fixed.trig[0].gap_raw * 0.25f &&
          fixed.trig[1].gap_aligned < fixed.trig[1].gap_raw * 0.25f);
    check("align_recovers_lag_direction",
          fixed.trig[0].align_off > 0.05f &&
          fixed.trig[0].align_off < FAR_LAG);
    check("align_not_identity",
          fabsf(fixed.trig[0].align_off - fixed.trig[0].align_base) > 1e-6f);
    check("no_align_leaves_gap",
          noalign.trig[0].gap_aligned > fixed.trig[0].gap_aligned * 1.25f);

    /* (3) both skeleton sizes record error and trigger time */
    check("near_rig_times",
          fixed.attach_frame[0] >= 0 && fixed.detach_frame[0] >= 0 &&
          fixed.attach_time[0] > fixed.detach_time[0]);
    check("far_rig_times",
          fixed.attach_frame[1] >= 0 && fixed.detach_frame[1] >= 0 &&
          fixed.attach_time[1] < fixed.detach_time[1]);
    check("near_err_recorded",
          fixed.attach_err[0] > 0.0f && fixed.detach_err[0] > 0.0f);
    check("far_err_recorded",
          fixed.attach_err[1] > 0.0f && fixed.detach_err[1] > 0.0f);
    check("rig_errors_differ",
          fabsf(fixed.attach_err[0] - fixed.attach_err[1]) > 1e-6f);
    check("trigger_time_is_frame_derived",
          fabsf(fixed.trig[0].time - (float)fixed.trig[0].frame * DT) < 1e-7f &&
          fabsf(fixed.trig[1].time - (float)fixed.trig[1].frame * DT) < 1e-7f);

    /* the pop itself: hard cut jumps, blended handoff does not */
    printf("cross_jump_at=%.6f\n", (double)cross.trig[0].jump_at);
    printf("fixed_jump_at=%.6f\n", (double)fixed.trig[0].jump_at);
    printf("near_natural_at_switch=%.6f\n",
           (double)natural_step(0, fixed.trig[0].time));
    check("hard_cut_pops",
          cross.trig[0].jump_at > 3.0f * natural_step(0, cross.trig[0].time));
    check("blended_attach_no_pop",
          fixed.trig[0].jump_at < cross.trig[0].jump_at * 0.5f);
    check("window_peak_beats_cut",
          fixed.trig[0].window_peak < cross.trig[0].jump_at * 0.25f &&
          fixed.trig[1].window_peak < cross.trig[1].jump_at * 0.35f);
    check("window_no_pop_near",
          fixed.window_step[0] <= fixed.natural_peak[0] * 1.75f + 1e-5f);
    check("window_no_pop_far",
          fixed.window_step[1] <= fixed.natural_peak[1] * 1.75f + 1e-5f);
    check("outbound_window_smooth",
          fixed.trig[0].window_excess <=
          0.05f * fixed.trig[0].nat_peak_window + 1e-5f);
    check("inbound_window_bounded",
          fixed.trig[1].window_excess <=
          0.5f * fixed.trig[1].nat_peak_window + 1e-5f);

    /* reproducibility: same frame, same time, run after run */
    simulate(POLICY_BAND, COST_THRESHOLD, TRANSITION_FRAMES, 1, &r2);
    simulate(POLICY_BAND, COST_THRESHOLD, TRANSITION_FRAMES, 1, &r3);
    int same = (r2.ntrig == fixed.ntrig && r3.ntrig == fixed.ntrig);
    float dev = 0.0f;
    for (int i = 0; i < fixed.ntrig && same; i++) {
        if (fixed.trig[i].frame != r2.trig[i].frame ||
            fixed.trig[i].frame != r3.trig[i].frame ||
            fixed.trig[i].time != r2.trig[i].time ||
            fixed.trig[i].time != r3.trig[i].time ||
            fixed.trig[i].cost != r2.trig[i].cost ||
            fixed.trig[i].cost != r3.trig[i].cost) same = 0;
        float d = fabsf(fixed.trig[i].time - r2.trig[i].time);
        float d3 = fabsf(fixed.trig[i].time - r3.trig[i].time);
        if (d3 > d) d = d3;
        if (d > dev) dev = d;
    }
    printf("repeat_same=%d\n", same);
    printf("timing_dev_sec=%.9f\n", (double)dev);
    printf("timing_dev_frames=%.6f\n", (double)(dev / DT));
    check("reproducible_runs", same == 1);
    check("timing_dev_within_half_frame", dev <= HALF_FRAME);

    /* per-frame trace: every 10th frame plus both trigger frames */
    printf("# frame time dist cost state window jump\n");
    for (int n = 0; n < NF; n++) {
        if (n % 10 != 0 && n != fixed.trig[0].frame &&
            n != fixed.trig[1].frame) continue;
        printf("trace %d %.4f %.4f %.6f %d %d %.6f\n",
               n, (double)((float)n * DT), (double)fixed.f[n].dist,
               (double)fixed.f[n].cost, fixed.f[n].state, fixed.f[n].window,
               (double)fixed.jump[n]);
    }

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
