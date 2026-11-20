
/* skeleton: bone hierarchy, world transforms and skinning matrix assembly.
   Optional unified rotation representation: rotations are accumulated as
   quaternions (normalized every frame) and composed down the parent chain,
   instead of summing ever-growing scalar angles. Off by default; when off,
   poses come from the original angle-sum path unchanged. */
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>

#define MAX_BONES 16

typedef struct { float x, y; float ang; } BoneLocal;
typedef struct { float x, y, ang; } BoneWorld;
typedef struct { float w, x, y, z; } Quat;

typedef struct {
    int parent[MAX_BONES];
    BoneLocal local[MAX_BONES];
    BoneWorld world[MAX_BONES];
    Quat local_q[MAX_BONES];   /* unified mode: local rotation state */
    Quat world_q[MAX_BONES];   /* unified mode: accumulated world rotation */
    int unified;               /* 0 = legacy angle-sum path (default) */
    int quat_error;            /* sticky: zero-norm quaternion encountered */
    int count;
} Skeleton;

static void skel_init(Skeleton* s) { memset(s, 0, sizeof(*s)); }

static Quat quat_from_angle_z(float a) {
    Quat q = { cosf(a * 0.5f), 0.0f, 0.0f, sinf(a * 0.5f) };
    return q;
}

static Quat quat_mul(Quat a, Quat b) {
    Quat r = {
        a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z,
        a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y,
        a.w * b.y - a.x * b.z + a.y * b.w + a.z * b.x,
        a.w * b.z + a.x * b.y - a.y * b.x + a.z * b.w
    };
    return r;
}

/* normalize in place; fails (returns -1) on a zero-norm quaternion */
static int quat_normalize(Quat* q) {
    float n = q->w * q->w + q->x * q->x + q->y * q->y + q->z * q->z;
    if (!(n > 0.0f)) return -1;    /* zero, denormal-underflowed or NaN */
    float inv = 1.0f / sqrtf(n);
    q->w *= inv; q->x *= inv; q->y *= inv; q->z *= inv;
    return 0;
}

static float quat_norm(Quat q) {
    return sqrtf(q.w * q.w + q.x * q.x + q.y * q.y + q.z * q.z);
}

/* rotate a vector in the xy plane by q (assumed unit) */
static void quat_rotate(Quat q, float x, float y, float* ox, float* oy) {
    float tx = -2.0f * q.z * y;              /* t = 2 * cross(q.xyz, v) */
    float ty =  2.0f * q.z * x;
    float tz =  2.0f * (q.x * y - q.y * x);
    *ox = x + q.w * tx + q.y * tz - q.z * ty;
    *oy = y + q.w * ty + q.z * tx - q.x * tz;
}

/* rotation angle about z represented by q */
static float quat_to_angle_z(Quat q) {
    return atan2f(2.0f * (q.w * q.z + q.x * q.y),
                  1.0f - 2.0f * (q.y * q.y + q.z * q.z));
}

static int skel_add(Skeleton* s, int parent, float x, float y, float ang) {
    if (s->count >= MAX_BONES) return -1;
    if (parent >= s->count) return -1;
    int i = s->count++;
    s->parent[i] = parent;
    s->local[i].x = x; s->local[i].y = y; s->local[i].ang = ang;
    s->local_q[i] = quat_from_angle_z(ang);
    return i;
}

/* enable/disable the unified quaternion path; enabling syncs the
   quaternion state from the current local angles */
static void skel_enable_unified(Skeleton* s, int on) {
    s->unified = on;
    s->quat_error = 0;
    if (on) {
        for (int i = 0; i < s->count; i++)
            s->local_q[i] = quat_from_angle_z(s->local[i].ang);
    }
}

/* absolute local rotation (e.g. from a curve); keeps both representations
   in sync so the flag can be toggled at any time */
static void skel_set_local_angle(Skeleton* s, int i, float ang) {
    s->local[i].ang = ang;
    s->local_q[i] = quat_from_angle_z(ang);
}

/* incremental local rotation (e.g. a spinning bone). Legacy mode grows the
   scalar angle; unified mode accumulates a normalized quaternion instead,
   so the state stays bounded no matter how long the run lasts. */
static int skel_rotate_local(Skeleton* s, int i, float dang) {
    if (!s->unified) {
        s->local[i].ang += dang;
        return 0;
    }
    Quat dq = quat_from_angle_z(dang);
    Quat q = quat_mul(s->local_q[i], dq);
    if (quat_normalize(&q) != 0) { s->quat_error = 1; return -1; }
    s->local_q[i] = q;
    s->local[i].ang = quat_to_angle_z(q);
    return 0;
}

static void rotate(float x, float y, float a, float* ox, float* oy) {
    float c = cosf(a), s = sinf(a);
    *ox = x * c - y * s;
    *oy = x * s + y * c;
}

/* original angle-sum path, kept byte-for-byte as the reference behavior */
static void skel_update_legacy(Skeleton* s) {
    for (int i = 0; i < s->count; i++) {
        if (s->parent[i] < 0) {
            s->world[i].x = s->local[i].x;
            s->world[i].y = s->local[i].y;
            s->world[i].ang = s->local[i].ang;
        } else {
            int p = s->parent[i];
            float rx, ry;
            rotate(s->local[i].x, s->local[i].y, s->world[p].ang, &rx, &ry);
            s->world[i].x = s->world[p].x + rx;
            s->world[i].y = s->world[p].y + ry;
            s->world[i].ang = s->world[p].ang + s->local[i].ang;
        }
    }
}

/* unified path: compose the accumulated parent world quaternion with the
   local quaternion (parent-chain compensation), normalize every frame,
   and derive positions by rotating offsets with the quaternion directly */
static int skel_update_unified(Skeleton* s) {
    for (int i = 0; i < s->count; i++) {
        Quat qw;
        float px, py;
        if (s->parent[i] < 0) {
            qw = s->local_q[i];
            px = s->local[i].x;
            py = s->local[i].y;
        } else {
            int p = s->parent[i];
            float rx, ry;
            qw = quat_mul(s->world_q[p], s->local_q[i]);
            quat_rotate(s->world_q[p], s->local[i].x, s->local[i].y, &rx, &ry);
            px = s->world[p].x + rx;
            py = s->world[p].y + ry;
        }
        if (quat_normalize(&qw) != 0) { s->quat_error = 1; return -1; }
        s->world_q[i] = qw;
        s->world[i].x = px;
        s->world[i].y = py;
        s->world[i].ang = quat_to_angle_z(qw);
    }
    return 0;
}

static int skel_update(Skeleton* s) {
    if (!s->unified) { skel_update_legacy(s); return 0; }
    return skel_update_unified(s);
}

/* bone start point and end point (length 1 along local +x) */
static void bone_segment(const Skeleton* s, int i, float* x0, float* y0, float* x1, float* y1) {
    *x0 = s->world[i].x; *y0 = s->world[i].y;
    float rx, ry; rotate(1.0f, 0.0f, s->world[i].ang, &rx, &ry);
    *x1 = s->world[i].x + rx; *y1 = s->world[i].y + ry;
}

/* conservative worst-case bound on end-effector position error of the
   unified path: each quaternion compose+normalize rounds the represented
   angle by at most ~4*FLT_EPSILON, and rounding adds linearly in the worst
   case (over accumulated local updates and down the chain) */
static double unified_tip_bound(double total_bone_len, long accum_ops) {
    return total_bone_len * (double)accum_ops * 4.0 * (double)FLT_EPSILON;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    Skeleton s;
    skel_init(&s);
    int root = skel_add(&s, -1, 0.0f, 0.0f, 0.0f);
    int b1 = skel_add(&s, root, 1.0f, 0.0f, 0.0f);
    int b2 = skel_add(&s, b1, 1.0f, 0.0f, 0.0f);
    check("add_root", root == 0);
    check("add_child", b1 == 1 && b2 == 2);

    skel_update(&s);
    check("world_chain", fabsf(s.world[2].x - 2.0f) < 1e-5f && fabsf(s.world[2].y) < 1e-5f);

    s.local[root].ang = 1.5707963f;   /* 90 deg */
    skel_update(&s);
    check("rotate_chain", fabsf(s.world[2].y - 2.0f) < 1e-4f && fabsf(s.world[2].x) < 1e-4f);

    float x0, y0, x1, y1;
    bone_segment(&s, b1, &x0, &y0, &x1, &y1);
    check("segment_dir", fabsf(x1 - x0) < 1e-4f && fabsf(y1 - y0 - 1.0f) < 1e-4f);

    check("reject_bad_parent", skel_add(&s, 9, 0, 0, 0) == -1);

    Skeleton full;
    skel_init(&full);
    int last = -1, rc = 0;
    for (int i = 0; i <= MAX_BONES; i++) { rc = skel_add(&full, last, 1, 0, 0); if (rc >= 0) last = rc; }
    check("capacity", rc == -1 && full.count == MAX_BONES);

    /* --- unified rotation representation --- */

    check("unified_off_by_default", s.unified == 0);

    /* legacy path byte-identical whether reached directly or via dispatcher */
    Skeleton la, lb;
    skel_init(&la); skel_init(&lb);
    skel_add(&la, -1, 0.5f, -0.25f, 0.3f); skel_add(&la, 0, 1.0f, 0.2f, -0.7f);
    skel_add(&lb, -1, 0.5f, -0.25f, 0.3f); skel_add(&lb, 0, 1.0f, 0.2f, -0.7f);
    skel_update_legacy(&la);
    skel_update(&lb);   /* unified == 0 -> dispatches to legacy */
    check("legacy_path_identical",
          memcmp(la.world, lb.world, sizeof(la.world)) == 0);

    /* static pose: unified matches legacy */
    Skeleton ua, ub;
    skel_init(&ua); skel_init(&ub);
    for (int i = 0; i < 2; i++) {
        Skeleton* t = i ? &ub : &ua;
        skel_add(t, -1, 0.0f, 0.0f, 0.4f);
        skel_add(t, 0, 1.0f, 0.0f, -0.9f);
        skel_add(t, 1, 1.0f, 0.0f, 1.1f);
    }
    skel_enable_unified(&ub, 1);
    skel_update(&ua);
    skel_update(&ub);
    check("unified_matches_legacy",
          fabsf(ua.world[2].x - ub.world[2].x) < 1e-5f &&
          fabsf(ua.world[2].y - ub.world[2].y) < 1e-5f &&
          fabsf(ua.world[2].ang - ub.world[2].ang) < 1e-4f);

    /* zero-norm quaternion must be reported, not silently used */
    Quat qz = { 0.0f, 0.0f, 0.0f, 0.0f };
    check("quat_normalize_zero", quat_normalize(&qz) == -1);

    Skeleton zs;
    skel_init(&zs);
    skel_add(&zs, -1, 0.0f, 0.0f, 0.0f);
    skel_add(&zs, 0, 1.0f, 0.0f, 0.0f);
    skel_enable_unified(&zs, 1);
    memset(&zs.local_q[1], 0, sizeof(Quat));   /* corrupt the state */
    check("zero_norm_error", skel_update(&zs) == -1 && zs.quat_error == 1);

    /* drift: 10 minutes at 60 fps. Root spins by a fixed increment every
       frame (angles grow without bound in legacy mode); children oscillate.
       Reference pose: matrices rebuilt from scratch each frame in double
       precision from the commanded angles. */
    enum { FRAMES = 10 * 60 * 60, N_BONES = 4 };
    const float total_len = (float)N_BONES;    /* unit bones */
    const float drift_threshold = 1e-3f;

    Skeleton uni, leg;
    skel_init(&uni); skel_init(&leg);
    for (int i = 0; i < N_BONES; i++) {
        skel_add(&uni, i - 1, i ? 1.0f : 0.0f, 0.0f, 0.0f);
        skel_add(&leg, i - 1, i ? 1.0f : 0.0f, 0.0f, 0.0f);
    }
    skel_enable_unified(&uni, 1);
    skel_set_local_angle(&uni, 3, 0.25f);
    leg.local[3].ang = 0.25f;

    double cmd[N_BONES] = { 0.0, 0.0, 0.0, 0.25 };  /* commanded local angles */
    double max_err_uni = 0.0, max_err_leg = 0.0, max_qdev = 0.0;
    int update_ok = 1;

    for (int f = 0; f < FRAMES; f++) {
        double t = (double)f / 60.0;
        cmd[0] += 0.01;                          /* spin: grows to 360 rad */
        cmd[1] = 0.6 * sin(t * 2.1);
        cmd[2] = 0.4 * cos(t * 1.3);

        if (skel_rotate_local(&uni, 0, 0.01f) != 0) update_ok = 0;
        skel_set_local_angle(&uni, 1, (float)cmd[1]);
        skel_set_local_angle(&uni, 2, (float)cmd[2]);
        leg.local[0].ang += 0.01f;
        leg.local[1].ang = (float)cmd[1];
        leg.local[2].ang = (float)cmd[2];
        if (skel_update(&uni) != 0) update_ok = 0;
        skel_update(&leg);

        /* reference: rebuild the chain matrices from scratch, in double */
        double wx[N_BONES], wy[N_BONES], wa[N_BONES];
        for (int i = 0; i < N_BONES; i++) {
            double lx = i ? 1.0 : 0.0;
            if (i == 0) { wx[i] = lx; wy[i] = 0.0; wa[i] = cmd[i]; }
            else {
                wa[i] = wa[i - 1] + cmd[i];
                wx[i] = wx[i - 1] + lx * cos(wa[i - 1]);
                wy[i] = wy[i - 1] + lx * sin(wa[i - 1]);
            }
        }
        double rtx = wx[N_BONES - 1] + cos(wa[N_BONES - 1]);
        double rty = wy[N_BONES - 1] + sin(wa[N_BONES - 1]);

        float ux0, uy0, ux1, uy1, lx0, ly0, lx1, ly1;
        bone_segment(&uni, N_BONES - 1, &ux0, &uy0, &ux1, &uy1);
        bone_segment(&leg, N_BONES - 1, &lx0, &ly0, &lx1, &ly1);
        double eu = hypot((double)ux1 - rtx, (double)uy1 - rty);
        double el = hypot((double)lx1 - rtx, (double)ly1 - rty);
        if (eu > max_err_uni) max_err_uni = eu;
        if (el > max_err_leg) max_err_leg = el;

        for (int i = 0; i < N_BONES; i++) {
            double dev = fabs((double)quat_norm(uni.world_q[i]) - 1.0);
            if (dev > max_qdev) max_qdev = dev;
        }
    }

    double bound = unified_tip_bound(total_len, FRAMES + N_BONES);
    printf("drift_max_unified=%.3g\n", max_err_uni);
    printf("drift_max_legacy=%.3g\n", max_err_leg);
    printf("drift_bound_worst_case=%.3g\n", bound);
    printf("quat_norm_max_dev=%.3g\n", max_qdev);

    check("ten_min_update_ok", update_ok == 1 && uni.quat_error == 0);
    check("quat_stays_unit", max_qdev < 1e-4);
    check("drift_under_threshold", max_err_uni < drift_threshold);
    check("drift_within_bound", max_err_uni < bound);
    check("unified_beats_legacy", max_err_uni < max_err_leg);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
