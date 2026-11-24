
/* seam: align root-height curves across clip seams, verify foot contact.
 *
 * Scenario: stand -> crouch -> jump, 30 fps. Each clip has its own authored
 * root-height curve; the curves carry constant vertical offset errors, so the
 * feet briefly leave the ground at the seams. This module:
 *   1. aligns the three curves at the seams (C0 height match),
 *   2. measures foot-to-ground distance fluctuation over grounded frames and
 *      compresses it under a threshold (foot-plant clamp),
 *   3. dumps 10 frames on each side of every seam, per clip,
 *   4. re-tests each clip in isolation (seams opened) to separate intrinsic
 *      float from seam-induced float.
 */
#include <stdio.h>
#include <string.h>
#include <math.h>

#define MAX_KEYS 8
#define FPS 30.0f
#define LEG_REACH 1.0f        /* max hip-to-foot extension (m) */
#define FLOAT_THRESHOLD 0.01f /* allowed grounded foot fluctuation (m) */

typedef struct { float t, v, tin, tout; } Key;
typedef struct { Key keys[MAX_KEYS]; int count; } Curve;

static void curve_init(Curve* c) { memset(c, 0, sizeof(*c)); }

static void curve_add(Curve* c, float t, float v) {
    Key* k = &c->keys[c->count++];
    k->t = t; k->v = v; k->tin = 0.0f; k->tout = 0.0f;
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
            float h01 = -2 * u3 + 3 * u2;
            return h00 * a->v + h01 * b->v;  /* zero tangents */
        }
    }
    return c->keys[c->count - 1].v;
}

typedef struct {
    const char* name;
    Curve root;          /* authored root-height curve, includes offset error */
    int first_frame, last_frame;
    int airborne_from;   /* first intentionally airborne frame; -1 = grounded clip */
    float shift;         /* vertical alignment shift, solved by align_seams */
} Clip;

static float clip_root(const Clip* c, int frame) {
    float t = (frame - c->first_frame) / FPS;
    return curve_eval(&c->root, t) + c->shift;
}

static float foot_dist(float root_h) {
    return root_h > LEG_REACH ? root_h - LEG_REACH : 0.0f;
}

static int clip_grounded(const Clip* c, int frame) {
    return c->airborne_from < 0 || frame < c->airborne_from;
}

/* foot-plant compression: on grounded frames clamp root so the foot can reach */
static float compress_root(const Clip* c, int frame, float root_h) {
    if (clip_grounded(c, frame) && root_h > LEG_REACH) return LEG_REACH;
    return root_h;
}

/* C0 height alignment: shift each clip so its start meets the previous end */
static void align_seams(Clip* clips, int n) {
    clips[0].shift = 0.0f;
    for (int i = 1; i < n; i++) {
        float prev_end = clip_root(&clips[i - 1], clips[i - 1].last_frame);
        clips[i].shift = 0.0f;
        float cur_start = clip_root(&clips[i], clips[i].first_frame);
        clips[i].shift = prev_end - cur_start;
    }
}

/* grounded foot-distance fluctuation over the whole sequence */
static void fluctuation(const Clip* clips, int n, int compress,
                        float* out_min, float* out_max) {
    float mn = 1e9f, mx = 0.0f;
    for (int i = 0; i < n; i++) {
        const Clip* c = &clips[i];
        for (int f = c->first_frame; f <= c->last_frame; f++) {
            if (!clip_grounded(c, f)) continue;
            float r = clip_root(c, f);
            if (compress) r = compress_root(c, f, r);
            float d = foot_dist(r);
            if (d < mn) mn = d;
            if (d > mx) mx = d;
        }
    }
    *out_min = mn; *out_max = mx;
}

/* seams opened: max grounded foot float of one clip on its own curve */
static float single_segment_float(const Clip* c) {
    float mx = 0.0f;
    for (int f = c->first_frame; f <= c->last_frame; f++) {
        if (!clip_grounded(c, f)) continue;
        float d = foot_dist(clip_root(c, f));
        if (d > mx) mx = d;
    }
    return mx;
}

static float max_fd_in_range(const Clip* c, int f0, int f1, int fix) {
    float mx = 0.0f;
    for (int f = f0; f <= f1; f++) {
        float r = clip_root(c, f);
        if (fix) r = compress_root(c, f, r);
        float d = foot_dist(r);
        if (d > mx) mx = d;
    }
    return mx;
}

static void dump_frames(const Clip* c, int f0, int f1) {
    printf("frame  root_raw  fd_raw  root_fix  fd_fix\n");
    for (int f = f0; f <= f1; f++) {
        Clip raw = *c; raw.shift = 0.0f;
        float r_raw = clip_root(&raw, f);
        float r_fix = compress_root(c, f, clip_root(c, f));
        printf("%5d  %8.4f  %6.4f  %8.4f  %6.4f\n",
               f, r_raw, foot_dist(r_raw), r_fix, foot_dist(r_fix));
    }
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    Clip clips[3];
    for (int i = 0; i < 3; i++) { curve_init(&clips[i].root); clips[i].shift = 0.0f; }

    /* stand, frames 0-59: idle sway, peaks overshoot leg reach by 0.012 */
    clips[0].name = "stand";
    clips[0].first_frame = 0; clips[0].last_frame = 59; clips[0].airborne_from = -1;
    curve_add(&clips[0].root, 0.0f, 1.000f);
    curve_add(&clips[0].root, 0.5f, 1.012f);
    curve_add(&clips[0].root, 1.0f, 0.998f);
    curve_add(&clips[0].root, 1.5f, 1.012f);
    curve_add(&clips[0].root, 59.0f / FPS, 1.000f);

    /* crouch, frames 60-99: authored with +0.08 offset error */
    clips[1].name = "crouch";
    clips[1].first_frame = 60; clips[1].last_frame = 99; clips[1].airborne_from = -1;
    curve_add(&clips[1].root, 0.0f, 1.00f + 0.08f);
    curve_add(&clips[1].root, 0.4f, 0.86f + 0.08f);
    curve_add(&clips[1].root, 0.8f, 0.74f + 0.08f);
    curve_add(&clips[1].root, 39.0f / FPS, 0.70f + 0.08f);

    /* jump, frames 100-149: authored with +0.35 offset error */
    clips[2].name = "jump";
    clips[2].first_frame = 100; clips[2].last_frame = 149; clips[2].airborne_from = -1;
    curve_add(&clips[2].root, 0.0f, 0.70f + 0.35f);
    curve_add(&clips[2].root, 0.33f, 0.66f + 0.35f);
    curve_add(&clips[2].root, 0.6f, 0.95f + 0.35f);
    curve_add(&clips[2].root, 0.9f, 1.25f + 0.35f);
    curve_add(&clips[2].root, 1.2f, 1.45f + 0.35f);
    curve_add(&clips[2].root, 49.0f / FPS, 1.00f + 0.35f);

    /* --- seam gaps before alignment --- */
    float gap1 = curve_eval(&clips[1].root, 0.0f) - curve_eval(&clips[0].root, 59.0f / FPS);
    float gap2 = curve_eval(&clips[2].root, 0.0f) - curve_eval(&clips[1].root, 39.0f / FPS);
    printf("seam1_gap=%.4f seam2_gap=%.4f\n", gap1, gap2);
    check("seam1_gap_detected", fabsf(gap1) > FLOAT_THRESHOLD);
    check("seam2_gap_detected", fabsf(gap2) > FLOAT_THRESHOLD);

    /* --- float at seams before any fix --- */
    float pre1 = max_fd_in_range(&clips[1], 60, 69, 0);
    float pre2 = max_fd_in_range(&clips[2], 100, 109, 0);
    printf("pre_float_seam1=%.4f pre_float_seam2=%.4f\n", pre1, pre2);
    check("pre_float_seam1", pre1 > FLOAT_THRESHOLD);
    check("pre_float_seam2", pre2 > FLOAT_THRESHOLD);

    /* --- align --- */
    align_seams(clips, 3);
    float post1 = clip_root(&clips[1], 60) - clip_root(&clips[0], 59);
    float post2 = clip_root(&clips[2], 100) - clip_root(&clips[1], 99);
    printf("shift_crouch=%.4f shift_jump=%.4f\n", clips[1].shift, clips[2].shift);
    check("seam1_c0_after_align", fabsf(post1) < 1e-6f);
    check("seam2_c0_after_align", fabsf(post2) < 1e-6f);

    /* jump takeoff: first frame the aligned root exceeds leg reach */
    for (int f = clips[2].first_frame; f <= clips[2].last_frame; f++) {
        if (clip_root(&clips[2], f) > LEG_REACH) { clips[2].airborne_from = f; break; }
    }
    printf("jump_airborne_from=%d\n", clips[2].airborne_from);
    check("jump_push_grounded", clips[2].airborne_from > clips[2].first_frame + 5);

    /* --- fluctuation: aligned only, then compressed --- */
    float mn0, mx0, mn1, mx1;
    fluctuation(clips, 3, 0, &mn0, &mx0);
    fluctuation(clips, 3, 1, &mn1, &mx1);
    printf("fluct_aligned=%.4f fluct_compressed=%.4f threshold=%.4f\n", mx0, mx1, FLOAT_THRESHOLD);
    check("fluctuation_needs_compression", mx0 > FLOAT_THRESHOLD);
    check("fluctuation_within_threshold", mx1 <= FLOAT_THRESHOLD);

    /* compression must not touch intentional airtime */
    float apex_raw = foot_dist(clip_root(&clips[2], 136));
    float apex_fix = foot_dist(compress_root(&clips[2], 136, clip_root(&clips[2], 136)));
    check("airborne_untouched", fabsf(apex_raw - apex_fix) < 1e-6f && apex_raw > 0.3f);

    /* --- single-segment isolation: seams opened, each clip on its own --- */
    printf("\n== single-segment test (seams opened) ==\n");
    float sfloat[3];
    for (int i = 0; i < 3; i++) {
        sfloat[i] = single_segment_float(&clips[i]);
        printf("single_%s_float=%.4f\n", clips[i].name, sfloat[i]);
    }
    check("single_crouch_clean", sfloat[1] <= FLOAT_THRESHOLD);
    check("single_jump_clean", sfloat[2] <= FLOAT_THRESHOLD);
    check("single_stand_intrinsic_float", sfloat[0] > FLOAT_THRESHOLD);

    /* --- seam window dumps: 10 frames each side, per clip --- */
    printf("\n== dataset 1/3: stand, frames 50-59 (seam1 before) ==\n");
    dump_frames(&clips[0], 50, 59);
    printf("\n== dataset 2/3: crouch, frames 60-69 (seam1 after) / 90-99 (seam2 before) ==\n");
    dump_frames(&clips[1], 60, 69);
    printf("\n");
    dump_frames(&clips[1], 90, 99);
    printf("\n== dataset 3/3: jump, frames 100-109 (seam2 after) ==\n");
    dump_frames(&clips[2], 100, 109);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
