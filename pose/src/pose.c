
/* pose: key pose interpolation. Each frame orients the bend plane from a
 * reference vector derived from the two end poses' joint directions
 * (never hardcoded), then solves two-bone IK inside that plane.
 * Reports knee-direction angle vs ground per frame and flags flips. */
#include <stdio.h>
#include <string.h>
#include <math.h>

typedef struct { float x, y; } V2;

static float vlen(V2 v) { return sqrtf(v.x * v.x + v.y * v.y); }
static float vdot(V2 a, V2 b) { return a.x * b.x + a.y * b.y; }
static V2 vsub(V2 a, V2 b) { V2 r = { a.x - b.x, a.y - b.y }; return r; }
static V2 vadd(V2 a, V2 b) { V2 r = { a.x + b.x, a.y + b.y }; return r; }
static V2 vscale(V2 v, float s) { V2 r = { v.x * s, v.y * s }; return r; }
static V2 vnorm(V2 v) { float l = vlen(v); V2 r = { 0, 0 }; if (l > 0) { r.x = v.x / l; r.y = v.y / l; } return r; }
static V2 vlerp(V2 a, V2 b, float t) { return vadd(vscale(a, 1.0f - t), vscale(b, t)); }

typedef struct {
    int ok;
    V2 joint;
} IkResult;

/* solve for target with bone lengths l1,l2; pole decides bend direction */
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
    r.ok = 1;
    return r;
}

typedef struct {
    V2 hip;
    V2 ankle;
    V2 knee;
} Pose;

/* bend direction: knee offset from the hip->ankle axis, normalized.
 * This is the "joint facing" used to orient the trajectory plane. */
static V2 bend_dir(V2 hip, V2 ankle, V2 knee) {
    V2 axis = vsub(ankle, hip);
    float len2 = vdot(axis, axis);
    if (len2 < 1e-12f) return vnorm(vsub(knee, hip));
    float u = vdot(vsub(knee, hip), axis) / len2;
    V2 on_axis = vadd(hip, vscale(axis, u));
    return vnorm(vsub(knee, on_axis));
}

/* knee direction angle vs ground (x axis), degrees */
static float ground_angle(V2 dir) {
    return atan2f(dir.y, dir.x) * 180.0f / 3.14159265f;
}

#define FRAMES_PER_SEG 20
#define NUM_SEGS 3
#define TOTAL_FRAMES (FRAMES_PER_SEG * NUM_SEGS)
#define FLIP_THRESHOLD 90.0f   /* deg; knee direction jump beyond this = flip */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    const float l1 = 0.45f, l2 = 0.45f;
    const V2 ankle = { 0.0f, 0.0f };
    const V2 key_pole = { 1.0f, 0.0f };   /* authoring hint for key poses only */

    /* stand, half-squat, full-squat, rise: hip heights, foot planted */
    const float key_hip_x[4] = { 0.00f, 0.02f, 0.05f, 0.00f };
    const float key_hip_y[4] = { 0.88f, 0.62f, 0.40f, 0.88f };

    Pose keys[4];
    for (int i = 0; i < 4; i++) {
        keys[i].hip.x = key_hip_x[i];
        keys[i].hip.y = key_hip_y[i];
        keys[i].ankle = ankle;
        IkResult r = solve_two_bone(keys[i].hip, ankle, l1, l2, key_pole);
        if (!r.ok) { printf("key_pose_%d_unsolvable=1\n", i); return 1; }
        keys[i].knee = r.joint;
    }

    float angles[TOTAL_FRAMES];
    V2 refs[TOTAL_FRAMES];
    V2 actuals[TOTAL_FRAMES];
    int flips[TOTAL_FRAMES];
    int frame = 0, flip_count = 0, solve_fail = 0;

    for (int seg = 0; seg < NUM_SEGS; seg++) {
        const Pose* a = &keys[seg];
        const Pose* b = &keys[seg + 1];
        /* orientation rule: reference vector interpolated between the two
         * end poses' joint directions, re-normalized per frame */
        V2 ref_a = bend_dir(a->hip, a->ankle, a->knee);
        V2 ref_b = bend_dir(b->hip, b->ankle, b->knee);

        for (int i = 1; i <= FRAMES_PER_SEG; i++) {
            float t = (float)i / FRAMES_PER_SEG;
            V2 hip = vlerp(a->hip, b->hip, t);
            V2 ref = vnorm(vlerp(ref_a, ref_b, t));
            if (vlen(ref) < 1e-6f) ref = ref_a;   /* degenerate: keep segment start */

            IkResult r = solve_two_bone(hip, ankle, l1, l2, ref);
            if (!r.ok) { solve_fail++; r.joint = hip; }

            V2 actual = bend_dir(hip, ankle, r.joint);
            float ang = ground_angle(actual);

            int flip = 0;
            if (frame > 0) {
                float d = fabsf(ang - angles[frame - 1]);
                if (d > 180.0f) d = 360.0f - d;
                flip = d > FLIP_THRESHOLD;
                if (flip) flip_count++;
            }

            angles[frame] = ang;
            refs[frame] = ref;
            actuals[frame] = actual;
            flips[frame] = flip;
            printf("frame=%02d seg=%d t=%.2f angle=%8.3f%s\n",
                   frame + 1, seg, t, ang, flip ? " FLIP" : "");
            frame++;
        }
    }

    printf("flip_frames=%d\n", flip_count);
    if (flip_count > 2) {
        /* reference vector suspect: show chosen ref vs actual facing */
        for (int f = 0; f < TOTAL_FRAMES; f++) {
            if (!flips[f]) continue;
            printf("flip_frame=%02d ref=(%.4f,%.4f) actual=(%.4f,%.4f)\n",
                   f + 1, refs[f].x, refs[f].y, actuals[f].x, actuals[f].y);
        }
    }

    check("all_frames_solved", solve_fail == 0);
    check("flips_within_tolerance", flip_count <= 2);
    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
