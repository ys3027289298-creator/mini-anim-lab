
/* retarget: mirror a pose across an automatically fitted symmetry plane.
 *
 * A pose is a tree of bone frames.  Each bone stores a local offset in its
 * parent's animated frame and a local rotation about z.  Mirroring uses the
 * reflection matrix
 *
 *     M = I - 2 n n^T                  (n = unit normal of the symmetry plane)
 *
 *   positions:   p' = M (p - o) + o    (o = a point on the plane)
 *   rotations:   R' = M R M            (conjugation; in 2D, a' = -a)
 *
 * Translation flips only the component along n -- the "specified axis" -- and
 * leaves the in-plane component untouched.  The M(p-o)+o form (rather than
 * M p) is the origin fix-up that makes a second mirror return the original
 * position; the a -> -a from the conjugation is the sign correction.
 *
 * Composition identity: for a reflection,
 *
 *     R(-A) M = M R(A)
 *
 * so if every local offset is reflected by M and every local rotation is
 * negated, the composed world origins are exact reflections, whatever the
 * plane angle.  mirror_pose does exactly that and is therefore an exact
 * involution.  The stored angle describes a bone *frame* (a proper rotation),
 * so the reflected frame is represented by -A; that is the reflection
 * composed with a half turn about n.  It is not the same number as the
 * reflected world direction (2*phi - A), and only the frame representation can
 * drive child offsets.  End effectors are consequently compared as bone
 * origins, which do reflect exactly.
 *
 * The plane is fitted by orthogonal least squares over the bones flagged as
 * sagittal (root, spine, neck, head): the normal is the minor eigenvector of
 * their location scatter, and the residual is the RMS perpendicular distance
 * of those points from the fitted line.  Limb bones are excluded because
 * their rest spread is exactly what mirroring is supposed to cancel.  If the
 * residual exceeds the caller's tolerance the plane is still returned, but
 * flagged MIRROR_WARN so the host can report it.
 *
 * mirror_retarget copies the conjugated local rotations onto the paired bones
 * (left limb -> right limb) and keeps the target skeleton's own rest offsets.
 * On a structurally asymmetric skeleton the right hand origin lands off the
 * mirrored left hand origin; the error is bounded by the sum of per-segment
 * length differences along the limb.
 */
#include <stdio.h>
#include <string.h>
#include <math.h>

#define MAX_BONES 16
#define PI 3.14159265358979f

typedef struct { float x, y; } V2;
typedef struct { V2 pos; float ang; } BonePose;
typedef struct { float x, y, ang; } World;

typedef struct {
    int parent[MAX_BONES];
    int sagittal[MAX_BONES];   /* 1 = lies on the symmetry plane by convention */
    int pair[MAX_BONES];       /* mirror partner bone (self when sagittal) */
    BonePose rest[MAX_BONES];  /* target rest pose, already mirrored */
    int count;
} Skeleton;

typedef struct { V2 origin; V2 normal; } Plane;

#define MIRROR_OK 0
#define MIRROR_WARN 1          /* fitted, but residual above tolerance */
#define MIRROR_NO_DATA 2

typedef struct {
    Plane plane;
    float residual;            /* RMS perpendicular distance, plane points */
    int n_points;
    int status;
} MirrorFit;

static V2 vmake(float x, float y) { V2 r = { x, y }; return r; }
static V2 vadd(V2 a, V2 b) { return vmake(a.x + b.x, a.y + b.y); }
static V2 vsub(V2 a, V2 b) { return vmake(a.x - b.x, a.y - b.y); }
static V2 vmul(V2 a, float s) { return vmake(a.x * s, a.y * s); }
static float vdot(V2 a, V2 b) { return a.x * b.x + a.y * b.y; }
static float vlen(V2 a) { return sqrtf(vdot(a, a)); }
static V2 vrot(V2 v, float a) {
    float c = cosf(a), s = sinf(a);
    return vmake(v.x * c - v.y * s, v.x * s + v.y * c);
}

static float wrap_pi(float a) {
    a = fmodf(a + PI, 2.0f * PI);
    if (a <= 0.0f) a += 2.0f * PI;
    return a - PI;
}

/* p' = M(p - o) + o; applying twice is the identity */
static V2 plane_reflect(const Plane* p, V2 v) {
    float d = vdot(vsub(v, p->origin), p->normal);
    return vsub(v, vmul(p->normal, 2.0f * d));
}

/* R' = M R M for M = I - 2 n n^T.  In 2D the result is rotation by -a for any
 * plane angle; computed through the matrix so the mapping is explicit. */
static float conj_angle(float a, const Plane* p) {
    float n0 = p->normal.x, n1 = p->normal.y;
    float m00 = 1.0f - 2.0f * n0 * n0, m01 = -2.0f * n0 * n1;
    float m10 = m01,                   m11 = 1.0f - 2.0f * n1 * n1;
    float c = cosf(a), s = sinf(a);
    float a00 = m00 * c + m01 * s, a01 = -m00 * s + m01 * c;
    float a10 = m10 * c + m11 * s, a11 = -m10 * s + m11 * c;
    float b00 = a00 * m00 + a01 * m10;
    float b10 = a10 * m00 + a11 * m10;
    return atan2f(b10, b00);
}

/* reflect a relative offset: translations have no origin, so only the M part
 * applies.  Rest offsets are stored in the parent's animated frame (rest local
 * rotations are zero), which is what makes d' = M d compose into an exact
 * world-space reflection. */
static V2 reflect_offset(const Plane* p, V2 d) {
    return vsub(d, vmul(p->normal, 2.0f * vdot(d, p->normal)));
}

static void pose_world(const Skeleton* s, const BonePose* pose, World* w) {
    for (int i = 0; i < s->count; i++) {
        if (s->parent[i] < 0) {
            w[i].x = pose[i].pos.x; w[i].y = pose[i].pos.y; w[i].ang = pose[i].ang;
        } else {
            int q = s->parent[i];
            V2 d = vrot(pose[i].pos, w[q].ang);
            w[i].x = w[q].x + d.x;
            w[i].y = w[q].y + d.y;
            w[i].ang = w[q].ang + pose[i].ang;
        }
    }
}

/* Orthogonal least squares line through the sagittal bones of the rest pose.
 * Normal = minor eigenvector of the location scatter; residual = RMS distance
 * of the points from the fitted line. */
static void fit_symmetry_axis(const Skeleton* s, float tol, MirrorFit* out) {
    World w[MAX_BONES];
    pose_world(s, s->rest, w);

    V2 pts[MAX_BONES];
    int n = 0;
    V2 mean = vmake(0, 0);
    for (int i = 0; i < s->count; i++) {
        if (!s->sagittal[i]) continue;
        pts[n++] = vmake(w[i].x, w[i].y);
        mean = vadd(mean, pts[n - 1]);
    }
    out->n_points = n;
    out->plane.origin = vmake(0, 0);
    out->plane.normal = vmake(1, 0);
    out->residual = 0.0f;
    out->status = MIRROR_NO_DATA;
    if (n < 2) return;
    mean = vmul(mean, 1.0f / (float)n);

    float sxx = 0, sxy = 0, syy = 0;
    for (int i = 0; i < n; i++) {
        V2 d = vsub(pts[i], mean);
        sxx += d.x * d.x; sxy += d.x * d.y; syy += d.y * d.y;
    }
    if (sxx + syy < 1e-12f) return;

    float theta = 0.5f * atan2f(2.0f * sxy, sxx - syy);   /* major axis */
    V2 normal = vmake(-sinf(theta), cosf(theta));          /* minor axis */
    if (normal.x < -1e-6f) normal = vmul(normal, -1.0f);   /* deterministic side */
    else if (normal.x <= 1e-6f && normal.y < 0) normal = vmul(normal, -1.0f);

    float sum = 0;
    for (int i = 0; i < n; i++) {
        float d = vdot(vsub(pts[i], mean), normal);
        sum += d * d;
    }
    out->plane.origin = mean;
    out->plane.normal = normal;
    out->residual = sqrtf(sum / (float)n);
    out->status = out->residual > tol ? MIRROR_WARN : MIRROR_OK;
}

/* geometric mirror on the same bone indices; exact involution */
static void mirror_pose(const Skeleton* s, const BonePose* in, const Plane* p, BonePose* out) {
    for (int i = 0; i < s->count; i++) {
        int q = s->parent[i];
        if (q < 0) {
            out[i].pos = plane_reflect(p, in[i].pos);
        } else {
            out[i].pos = reflect_offset(p, in[i].pos);
        }
        out[i].ang = conj_angle(in[i].ang, p);
    }
}

/* retarget: mirrored local rotations onto the paired bones; the target keeps
 * its own (pre-mirrored) rest offsets.  Root translation is mirrored across
 * the plane, child offsets come from the target's rest pose. */
static void mirror_retarget(const Skeleton* s, const BonePose* in, const Plane* p, BonePose* out) {
    for (int i = 0; i < s->count; i++) {
        int j = s->pair[i];
        if (j < 0 || j >= s->count) j = i;
        out[j] = s->rest[j];
        out[j].ang = conj_angle(in[i].ang, p);
    }
    int root = -1;
    for (int i = 0; i < s->count; i++) if (s->parent[i] < 0) { root = i; break; }
    if (root >= 0) out[s->pair[root]].pos = plane_reflect(p, in[root].pos);
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static int add_bone(Skeleton* s, int parent, float x, float y) {
    if (s->count >= MAX_BONES) return -1;
    int i = s->count++;
    s->parent[i] = parent;
    s->rest[i].pos = vmake(x, y);
    s->rest[i].ang = 0.0f;
    s->sagittal[i] = 0;
    s->pair[i] = i;
    return i;
}

/* front-view humanoid; +x is the character's left, +y up.  head_dx != 0
 * pulls the sagittal points off the y axis.  The right arm has different
 * segment lengths, so the skeleton is not symmetric in size. */
static void build_body(Skeleton* s, float head_dx) {
    memset(s, 0, sizeof(*s));
    add_bone(s, -1, 0.0f, 0.0f);          /* 0 pelvis  */
    add_bone(s, 0, 0.0f, 0.6f);           /* 1 spine   */
    add_bone(s, 1, 0.0f, 0.6f);           /* 2 chest   */
    add_bone(s, 2, 0.0f, 0.35f);          /* 3 neck    */
    add_bone(s, 3, head_dx, 0.3f);        /* 4 head    */
    add_bone(s, 2, 0.22f, 0.5f);          /* 5 l_shoulder */
    add_bone(s, 5, 0.45f, -0.05f);        /* 6 l_elbow */
    add_bone(s, 6, 0.42f, 0.0f);          /* 7 l_hand  */
    add_bone(s, 2, -0.22f, 0.5f);         /* 8 r_shoulder */
    add_bone(s, 8, -0.40f, -0.05f);       /* 9 r_elbow  (shorter) */
    add_bone(s, 9, -0.38f, 0.0f);         /* 10 r_hand  (shorter) */
    for (int i = 0; i <= 4; i++) s->sagittal[i] = 1;
    s->pair[5] = 8;  s->pair[8] = 5;
    s->pair[6] = 9;  s->pair[9] = 6;
    s->pair[7] = 10; s->pair[10] = 7;
}

static void build_pose(const Skeleton* s, BonePose* p) {
    for (int i = 0; i < s->count; i++) p[i] = s->rest[i];
    p[0].pos = vmake(0.05f, 0.02f);
    p[0].ang = 0.03f;
    p[1].ang = 0.10f;
    p[2].ang = 0.05f;
    p[3].ang = -0.08f;
    p[4].ang = 0.04f;
    p[5].ang = 0.60f;
    p[6].ang = 0.45f;
    p[7].ang = 0.20f;
    p[8].ang = -0.15f;
    p[9].ang = 0.30f;
    p[10].ang = -0.10f;
}

/* hand-authored mirror, independent of the fitted plane: negate every local
 * rotation and flip only the root translation's x (valid once the plane is
 * checked to be the vertical axis). */
static void manual_mirror(const Skeleton* s, const BonePose* in, BonePose* out) {
    for (int i = 0; i < s->count; i++) {
        int j = s->pair[i];
        out[j] = s->rest[j];
        out[j].ang = -in[i].ang;
    }
    int root = -1;
    for (int i = 0; i < s->count; i++) if (s->parent[i] < 0) { root = i; break; }
    if (root >= 0) out[root].pos = vmake(-in[root].pos.x, in[root].pos.y);
}

static float max_pose_diff(const Skeleton* s, const BonePose* a, const BonePose* b) {
    float m = 0;
    for (int i = 0; i < s->count; i++) {
        float d = vlen(vsub(a[i].pos, b[i].pos));
        float da = fabsf(wrap_pi(a[i].ang - b[i].ang));
        if (d > m) m = d;
        if (da > m) m = da;
    }
    return m;
}

/* bound on |mirror_f(T) - mirror_t(T)| for fitted plane f vs true plane t:
 *   2|dn|(p + |T - o_t|) + 2|do|,  p = distance of T from the true plane */
static float mirror_error_bound(const Plane* f, const Plane* t, V2 T) {
    float dn = vlen(vsub(f->normal, t->normal));
    float doff = vlen(vsub(f->origin, t->origin));
    V2 rel = vsub(T, t->origin);
    float p = fabsf(vdot(rel, t->normal));
    return 2.0f * dn * (p + vlen(rel)) + 2.0f * doff;
}

int main(void) {
    Skeleton body;
    build_body(&body, 0.0f);
    BonePose pose[MAX_BONES], mirrored[MAX_BONES], twice[MAX_BONES];
    BonePose manual[MAX_BONES], retargeted[MAX_BONES];
    build_pose(&body, pose);

    /* --- symmetry axis fit ------------------------------------------- */
    MirrorFit fit;
    fit_symmetry_axis(&body, 0.01f, &fit);
    check("fit_ok", fit.status == MIRROR_OK);
    check("fit_on_plane", fit.n_points == 5 && fit.residual < 1e-5f);
    check("axis_is_x", fit.plane.normal.x > 0.9999f && fabsf(fit.plane.normal.y) < 1e-4f);

    /* translation flips only the plane-normal axis (here x) */
    V2 t = plane_reflect(&fit.plane, vmake(0.3f, 0.7f));
    check("translation_only_x", fabsf(t.x + 0.3f) < 1e-5f && fabsf(t.y - 0.7f) < 1e-5f);

    /* R' = M R M is the sign correction a -> -a for any 2D plane */
    check("rotation_conjugation", fabsf(wrap_pi(conj_angle(0.7f, &fit.plane) + 0.7f)) < 1e-5f);

    /* --- geometric mirror -------------------------------------------- */
    mirror_pose(&body, pose, &fit.plane, mirrored);

    World sw[MAX_BONES], mw[MAX_BONES];
    pose_world(&body, pose, sw);
    pose_world(&body, mirrored, mw);
    float werr = 0.0f, aerr = 0.0f;
    for (int i = 0; i < body.count; i++) {
        V2 r = plane_reflect(&fit.plane, vmake(sw[i].x, sw[i].y));
        float d = vlen(vsub(vmake(mw[i].x, mw[i].y), r));
        float da = fabsf(wrap_pi(mw[i].ang + sw[i].ang));
        if (d > werr) werr = d;
        if (da > aerr) aerr = da;
    }
    check("mirror_world_position_exact", werr < 1e-5f);
    check("mirror_world_angle_negates", aerr < 1e-5f);
    float angle_only = 0.0f;
    for (int i = 0; i < body.count; i++) {
        float da = fabsf(wrap_pi(mirrored[i].ang + pose[i].ang));
        if (da > angle_only) angle_only = da;
    }
    check("local_angles_negate", angle_only < 1e-5f);
    check("root_translation_x_only",
          fabsf(mirrored[0].pos.x + pose[0].pos.x) < 1e-5f &&
          fabsf(mirrored[0].pos.y - pose[0].pos.y) < 1e-5f);

    mirror_pose(&body, mirrored, &fit.plane, twice);
    check("mirror_twice_identity", max_pose_diff(&body, pose, twice) < 1e-5f);

    /* --- retarget to the paired limb, vs hand-authored mirror --------- */
    mirror_retarget(&body, pose, &fit.plane, retargeted);
    manual_mirror(&body, pose, manual);
    check("retarget_matches_manual", max_pose_diff(&body, retargeted, manual) < 1e-5f);

    /* end effector: right hand origin vs mirrored left hand origin.  The
     * right arm segments are shorter, so the error is bounded by the sum of
     * the per-segment length differences. */
    World rw[MAX_BONES];
    pose_world(&body, retargeted, rw);
    V2 tgt = vmake(rw[10].x, rw[10].y);
    V2 ref = plane_reflect(&fit.plane, vmake(sw[7].x, sw[7].y));
    float ee_err = vlen(vsub(tgt, ref));
    float ee_bound = 0.0f;
    int chain_src[3] = { 5, 6, 7 }, chain_dst[3] = { 8, 9, 10 };
    for (int k = 0; k < 3; k++)
        ee_bound += fabsf(vlen(body.rest[chain_src[k]].pos) - vlen(body.rest[chain_dst[k]].pos));
    check("retarget_endpoint_within_bound", ee_err <= ee_bound + 1e-5f);
    printf("endpoint_bound=%.6f endpoint_error=%.6f\n", ee_bound, ee_err);

    /* --- asymmetric skeleton: residual above tolerance must be reported - */
    Skeleton skew;
    build_body(&skew, 0.04f);
    MirrorFit strict, loose;
    fit_symmetry_axis(&skew, 0.01f, &strict);
    fit_symmetry_axis(&skew, 0.03f, &loose);
    check("residual_warn_above_tol", strict.status == MIRROR_WARN);
    check("residual_ok_below_tol", loose.status == MIRROR_OK);
    check("residual_reported", strict.residual > 0.01f && strict.residual == loose.residual);
    if (strict.status == MIRROR_WARN)
        printf("report=symmetry_residual residual=%.6f tol=%.6f status=warn\n",
               strict.residual, 0.01f);

    /* fitted-plane mirror vs true-plane (x = 0) mirror, with a proved bound */
    Plane truth;
    truth.origin = vmake(0.0f, strict.plane.origin.y);
    truth.normal = vmake(1.0f, 0.0f);
    BonePose sp[MAX_BONES], sf[MAX_BONES], st[MAX_BONES];
    build_pose(&skew, sp);
    mirror_pose(&skew, sp, &strict.plane, sf);
    mirror_pose(&skew, sp, &truth, st);
    World sw2[MAX_BONES], sfw[MAX_BONES], stw[MAX_BONES];
    pose_world(&skew, sp, sw2);
    pose_world(&skew, sf, sfw);
    pose_world(&skew, st, stw);
    V2 T = vmake(sw2[7].x, sw2[7].y);
    float axis_err = vlen(vsub(vmake(sfw[7].x, sfw[7].y), vmake(stw[7].x, stw[7].y)));
    float axis_bound = mirror_error_bound(&strict.plane, &truth, T);
    check("fitted_axis_endpoint_within_bound", axis_err <= axis_bound + 1e-5f);
    printf("axis_bound=%.6f axis_error=%.6f\n", axis_bound, axis_err);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
