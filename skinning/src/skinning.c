
/* skinning: linear blend skinning (LBS) vs dual quaternion skinning (DQS).
 *
 * DQS pipeline per vertex:
 *   1. each bone rigid transform (R,t) is composed into a dual quaternion
 *        real part  r = R
 *        dual part  d = 0.5 * t_q * r,  t_q = (0, tx, ty, tz)
 *   1b. skin transform S = W * B^-1 (bind inverse, rigid for both matrices
 *       and dual quaternions)
 *   2. weight-blend the skin dual quaternions (with antipodality correction)
 *   3. normalize the blended dual quaternion
 *   4. position: p' = r p r* + 2 * (d r*).xyz
 *      normal : rotation part r only (translation never touches normals)
 *   5. a zero dual part or a non-positive weight sum returns an error code
 *
 * LBS blends 4x4 matrices; positions use the blended matrix, normals use the
 * inverse-transpose direction (implemented via the cofactor matrix so a
 * near-degenerate blend cannot divide by a zero determinant).
 *
 * The two skinners run on the same mesh and the same poses; per joint they
 * report relative volume error and mean normal angle error.
 */
#include <stdio.h>
#include <string.h>
#include <math.h>

#define MAX_BONES 16
#define BONES 3                 /* bone chain length used by the test mesh */
#define SEG_PER_BONE 4
#define RADIAL 8
#define RINGS (BONES * SEG_PER_BONE + 1)
#define MAX_VERTS (RINGS * RADIAL)
#define MAX_TRIS ((RINGS - 1) * RADIAL * 2)
#define TUBE_R 0.25f
#define BONE_LEN 1.0f
#define PI_F 3.14159265358979323846f
#define EPS 1e-8f

typedef struct { float x, y, z; } V3;
typedef struct { float w, x, y, z; } Quat;               /* w is scalar part */
typedef struct { Quat r, d; } DQ;                         /* real and dual */
typedef struct { float m[4][4]; } M4;

enum { DQ_OK = 0, DQ_ERR_ZERO_DUAL = 1, DQ_ERR_ZERO_WEIGHT = 2 };

/* ---------- vec3 ---------- */

static V3 vadd(V3 a, V3 b) { return (V3){ a.x + b.x, a.y + b.y, a.z + b.z }; }
static V3 vsub(V3 a, V3 b) { return (V3){ a.x - b.x, a.y - b.y, a.z - b.z }; }
static V3 vscale(V3 a, float s) { return (V3){ a.x * s, a.y * s, a.z * s }; }
static float vdot(V3 a, V3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static V3 vcross(V3 a, V3 b) {
    return (V3){ a.y * b.z - a.z * b.y,
                 a.z * b.x - a.x * b.z,
                 a.x * b.y - a.y * b.x };
}
static float vlen(V3 a) { return sqrtf(vdot(a, a)); }
static V3 vnormalize(V3 a) {
    float l = vlen(a);
    return l > EPS ? vscale(a, 1.0f / l) : (V3){ 0, 0, 0 };
}

/* ---------- quaternion ---------- */

static Quat q_mul(Quat a, Quat b) {
    return (Quat){
        a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z,
        a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y,
        a.w * b.y - a.x * b.z + a.y * b.w + a.z * b.x,
        a.w * b.z + a.x * b.y - a.y * b.x + a.z * b.w
    };
}
static Quat q_conj(Quat q) { return (Quat){ q.w, -q.x, -q.y, -q.z }; }
static float q_dot(Quat a, Quat b) {
    return a.w * b.w + a.x * b.x + a.y * b.y + a.z * b.z;
}
static float q_norm(Quat a) { return sqrtf(q_dot(a, a)); }
static Quat q_scale(Quat a, float s) {
    return (Quat){ a.w * s, a.x * s, a.y * s, a.z * s };
}
static Quat q_axis(V3 axis, float rad) {
    V3 u = vnormalize(axis);
    float h = rad * 0.5f, s = sinf(h);
    return (Quat){ cosf(h), u.x * s, u.y * s, u.z * s };
}
static V3 q_rotate(Quat q, V3 v) {
    Quat p = { 0, v.x, v.y, v.z };
    Quat r = q_mul(q_mul(q, p), q_conj(q));
    return (V3){ r.x, r.y, r.z };
}

/* ---------- dual quaternion ---------- */

/* compose rigid transform (R,t): r = R, d = 0.5 * t_q * r */
static DQ dq_rigid(Quat r, V3 t) {
    Quat tq = { 0, t.x, t.y, t.z };
    DQ q;
    q.r = r;
    q.d = q_scale(q_mul(tq, r), 0.5f);
    return q;
}

/* dual quaternion product: (r1 + e d1)(r2 + e d2) = r1 r2 + e(d1 r2 + r1 d2) */
static DQ dq_mul(DQ a, DQ b) {
    DQ q;
    q.r = q_mul(a.r, b.r);
    q.d = q_mul(a.d, b.r);
    Quat tmp = q_mul(a.r, b.d);
    q.d = (Quat){ q.d.w + tmp.w, q.d.x + tmp.x, q.d.y + tmp.y, q.d.z + tmp.z };
    return q;
}

/* inverse of a unit rigid dual quaternion (R,t): (R^T, -R^T t) */
static DQ dq_inverse_rigid(const DQ* q) {
    Quat ri = q_conj(q->r);
    Quat tq = q_mul(q->d, q_conj(q->r));
    V3 t = { 2.0f * tq.x, 2.0f * tq.y, 2.0f * tq.z };
    return dq_rigid(ri, vscale(t, -1.0f));
}

/* a dual quaternion with a zero dual part carries no translation and is
 * rejected as degenerate by every skinning entry point */
static int dq_check(const DQ* q) {
    float dl = q->d.w * q->d.w + q->d.x * q->d.x +
               q->d.y * q->d.y + q->d.z * q->d.z;
    return dl > EPS ? DQ_OK : DQ_ERR_ZERO_DUAL;
}

/* weight-blend then normalize; antipodality fix keeps all real parts in the
 * same hemisphere before accumulation */
static int dq_blend(const DQ* dqs, const float* w, int n, DQ* out) {
    float wsum = 0.0f;
    for (int i = 0; i < n; i++) wsum += w[i];
    if (!(wsum > EPS)) return DQ_ERR_ZERO_WEIGHT;

    DQ acc;
    memset(&acc, 0, sizeof(acc));
    for (int i = 0; i < n; i++) {
        int rc = dq_check(&dqs[i]);
        if (rc != DQ_OK) return rc;
        DQ q = dqs[i];
        if (i > 0 && q_dot(q.r, dqs[0].r) < 0.0f) {
            q.r = q_scale(q.r, -1.0f);
            q.d = q_scale(q.d, -1.0f);
        }
        acc.r = (Quat){ acc.r.w + w[i] * q.r.w, acc.r.x + w[i] * q.r.x,
                        acc.r.y + w[i] * q.r.y, acc.r.z + w[i] * q.r.z };
        acc.d = (Quat){ acc.d.w + w[i] * q.d.w, acc.d.x + w[i] * q.d.x,
                        acc.d.y + w[i] * q.d.y, acc.d.z + w[i] * q.d.z };
    }
    float nl = q_norm(acc.r);
    if (!(nl > EPS)) return DQ_ERR_ZERO_DUAL;
    out->r = q_scale(acc.r, 1.0f / nl);
    out->d = q_scale(acc.d, 1.0f / nl);
    return DQ_OK;
}

/* point transform: p' = r p r* + 2 (d r*).xyz */
static V3 dq_point(const DQ* q, V3 p) {
    Quat tq = q_mul(q->d, q_conj(q->r));
    V3 t = { 2.0f * tq.x, 2.0f * tq.y, 2.0f * tq.z };
    return vadd(q_rotate(q->r, p), t);
}

/* normals take the rotation part only */
static V3 dq_normal(const DQ* q, V3 n) {
    return vnormalize(q_rotate(q->r, n));
}

/* ---------- 4x4 matrix path (LBS) ---------- */

static M4 m_rigid(Quat q, V3 t) {
    float x = q.x, y = q.y, z = q.z, w = q.w;
    M4 m;
    memset(&m, 0, sizeof(m));
    m.m[0][0] = 1 - 2 * (y * y + z * z);
    m.m[0][1] = 2 * (x * y - z * w);
    m.m[0][2] = 2 * (x * z + y * w);
    m.m[1][0] = 2 * (x * y + z * w);
    m.m[1][1] = 1 - 2 * (x * x + z * z);
    m.m[1][2] = 2 * (y * z - x * w);
    m.m[2][0] = 2 * (x * z - y * w);
    m.m[2][1] = 2 * (y * z + x * w);
    m.m[2][2] = 1 - 2 * (x * x + y * y);
    m.m[0][3] = t.x; m.m[1][3] = t.y; m.m[2][3] = t.z;
    m.m[3][3] = 1.0f;
    return m;
}

static M4 m_weighted(const M4* mats, const float* w, int n) {
    M4 acc;
    memset(&acc, 0, sizeof(acc));
    for (int i = 0; i < n; i++)
        for (int r = 0; r < 4; r++)
            for (int c = 0; c < 4; c++)
                acc.m[r][c] += w[i] * mats[i].m[r][c];
    return acc;
}

static V3 m_point(const M4* m, V3 p) {
    return (V3){ m->m[0][0] * p.x + m->m[0][1] * p.y + m->m[0][2] * p.z + m->m[0][3],
                 m->m[1][0] * p.x + m->m[1][1] * p.y + m->m[1][2] * p.z + m->m[1][3],
                 m->m[2][0] * p.x + m->m[2][1] * p.y + m->m[2][2] * p.z + m->m[2][3] };
}

/* inverse of a rigid 4x4 [R t; 0 1]: [R^T -R^T t; 0 1] */
static M4 m_inverse_rigid(const M4* m) {
    M4 inv;
    memset(&inv, 0, sizeof(inv));
    for (int r = 0; r < 3; r++)
        for (int c = 0; c < 3; c++)
            inv.m[r][c] = m->m[c][r];
    for (int r = 0; r < 3; r++) {
        float s = 0.0f;
        for (int c = 0; c < 3; c++) s += inv.m[r][c] * m->m[c][3];
        inv.m[r][3] = -s;
    }
    inv.m[3][3] = 1.0f;
    return inv;
}

/* inverse-transpose direction without division: C = cof(A) = det(A) A^-T,
 * so normalize(C n) equals the inverse-transpose normal for any nonzero det */
static V3 m_normal(const M4* m, V3 n) {
    const float (*a)[4] = m->m;
    float C[3][3];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            int rows[2], cols[2], nr = 0, nc = 0;
            for (int k = 0; k < 3; k++) {
                if (k != i) rows[nr++] = k;
                if (k != j) cols[nc++] = k;
            }
            float minor = a[rows[0]][cols[0]] * a[rows[1]][cols[1]]
                        - a[rows[0]][cols[1]] * a[rows[1]][cols[0]];
            C[i][j] = ((i + j) & 1) ? -minor : minor;
        }
    }
    V3 r = { C[0][0] * n.x + C[0][1] * n.y + C[0][2] * n.z,
             C[1][0] * n.x + C[1][1] * n.y + C[1][2] * n.z,
             C[2][0] * n.x + C[2][1] * n.y + C[2][2] * n.z };
    return vnormalize(r);
}

/* ---------- test mesh: a tube around the bone chain ---------- */

typedef struct { V3 pos, nrm; float w[MAX_BONES]; } Vert;
typedef struct { int a, b, c; } Tri;
typedef struct { Vert v[MAX_VERTS]; Tri t[MAX_TRIS]; int nv, nt; } Mesh;
typedef struct { V3 pos[MAX_VERTS], nrm[MAX_VERTS]; int nv; } Skin;

static int ring_vid(int ring, int k) { return ring * RADIAL + (k % RADIAL); }

static void mesh_build(Mesh* mesh, V3 base) {
    memset(mesh, 0, sizeof(*mesh));
    for (int i = 0; i < RINGS; i++) {
        float s = (float)i / SEG_PER_BONE;          /* bone-space parameter */
        int b0 = (int)s;
        float f = s - b0;
        if (b0 >= BONES - 1) { b0 = BONES - 1; f = 0.0f; }   /* last bone: full weight */
        for (int k = 0; k < RADIAL; k++) {
            float th = 2.0f * PI_F * k / RADIAL;
            Vert* vt = &mesh->v[mesh->nv];
            vt->pos = (V3){ base.x + i * BONE_LEN / SEG_PER_BONE,
                            base.y + TUBE_R * cosf(th),
                            base.z + TUBE_R * sinf(th) };
            vt->nrm = (V3){ 0, cosf(th), sinf(th) };
            vt->w[b0] = 1.0f - f;
            if (f > 0.0f) vt->w[b0 + 1] = f;
            mesh->nv++;
        }
    }
    for (int i = 0; i + 1 < RINGS; i++)
        for (int k = 0; k < RADIAL; k++) {
            int a = ring_vid(i, k), b = ring_vid(i, k + 1);
            int c = ring_vid(i + 1, k), d = ring_vid(i + 1, k + 1);
            /* winding chosen so the geometric normal points away from the axis */
            mesh->t[mesh->nt++] = (Tri){ a, b, c };
            mesh->t[mesh->nt++] = (Tri){ b, d, c };
        }
}

static void skin_lbs(const Mesh* mesh, const M4* mats, Skin* out) {
    out->nv = mesh->nv;
    for (int i = 0; i < mesh->nv; i++) {
        M4 m = m_weighted(mats, mesh->v[i].w, BONES);
        out->pos[i] = m_point(&m, mesh->v[i].pos);
        out->nrm[i] = m_normal(&m, mesh->v[i].nrm);
    }
}

static int skin_dqs(const Mesh* mesh, const DQ* dqs, Skin* out) {
    out->nv = 0;
    for (int i = 0; i < mesh->nv; i++) {
        DQ blended;
        int rc = dq_blend(dqs, mesh->v[i].w, BONES, &blended);
        if (rc != DQ_OK) return rc;
        out->pos[i] = dq_point(&blended, mesh->v[i].pos);
        out->nrm[i] = dq_normal(&blended, mesh->v[i].nrm);
    }
    out->nv = mesh->nv;
    return DQ_OK;
}

/* ---------- per-joint metrics ---------- */
/* joint j owns the slab made of rings jK-1 .. jK+1 (one segment each side),
 * closed with cap fans so the signed-tetrahedron mesh volume is well defined */

static V3 ring_center(const V3* pos, int ring) {
    V3 c = { 0, 0, 0 };
    for (int k = 0; k < RADIAL; k++) c = vadd(c, pos[ring_vid(ring, k)]);
    return vscale(c, 1.0f / RADIAL);
}

static float tet_vol(V3 a, V3 b, V3 c) { return vdot(a, vcross(b, c)) / 6.0f; }

/* add one triangle's tetrahedron volume, flipping the winding when its
 * geometric normal faces inward so the whole closed slab is outward-consistent
 * even after extreme bends */
static void add_oriented(float* vol, V3 a, V3 b, V3 c, V3 outward) {
    V3 n = vcross(vsub(b, a), vsub(c, a));
    if (vdot(n, outward) < 0.0f) { V3 t = b; b = c; c = t; }
    *vol += tet_vol(a, b, c);
}

static float slab_volume(const V3* pos, int joint) {
    int r0 = joint * SEG_PER_BONE - 1, r1 = joint * SEG_PER_BONE,
        r2 = joint * SEG_PER_BONE + 1;
    V3 c0 = ring_center(pos, r0), c1 = ring_center(pos, r1),
       c2 = ring_center(pos, r2);
    float vol = 0.0f;
    for (int k = 0; k < RADIAL; k++) {
        int kn = (k + 1) % RADIAL;
        /* side triangles over the two segments; outward = away from ring center */
        V3 o1 = vsub(pos[ring_vid(r1, k)], c1);
        add_oriented(&vol, pos[ring_vid(r0, k)], pos[ring_vid(r1, kn)],
                     pos[ring_vid(r1, k)], o1);
        add_oriented(&vol, pos[ring_vid(r0, k)], pos[ring_vid(r0, kn)],
                     pos[ring_vid(r1, kn)], o1);
        V3 o2 = vsub(pos[ring_vid(r2, k)], c2);
        add_oriented(&vol, pos[ring_vid(r1, k)], pos[ring_vid(r2, kn)],
                     pos[ring_vid(r2, k)], o2);
        add_oriented(&vol, pos[ring_vid(r1, k)], pos[ring_vid(r1, kn)],
                     pos[ring_vid(r2, kn)], o2);
        /* cap fans; outward = slab axis direction at each end */
        V3 axis0 = vnormalize(vsub(c1, c0)), axis2 = vnormalize(vsub(c2, c1));
        add_oriented(&vol, c0, pos[ring_vid(r0, k)],
                     pos[ring_vid(r0, kn)], vscale(axis0, -1.0f));
        add_oriented(&vol, c2, pos[ring_vid(r2, k)],
                     pos[ring_vid(r2, kn)], axis2);
    }
    return fabsf(vol);
}

/* ground-truth normal of the deformed surface at each vertex: area-weighted
 * average of the normals of the incident triangles, recomputed from the
 * deformed positions. This measures whether a skinner's output normal matches
 * the geometry it actually produced — the candy-wrapper artifact shows up here
 * as an output normal pointing away from the real surface */
static void geo_normals(const Mesh* mesh, const V3* pos, V3* out) {
    V3 acc[MAX_VERTS];
    memset(acc, 0, sizeof(acc));
    for (int t = 0; t < mesh->nt; t++) {
        Tri f = mesh->t[t];
        V3 n = vcross(vsub(pos[f.b], pos[f.a]), vsub(pos[f.c], pos[f.a]));
        acc[f.a] = vadd(acc[f.a], n);
        acc[f.b] = vadd(acc[f.b], n);
        acc[f.c] = vadd(acc[f.c], n);
    }
    for (int i = 0; i < mesh->nv; i++) out[i] = vnormalize(acc[i]);
}

/* mean angle (deg) between the skinner's output normals and the recomputed
 * surface normals over the joint's three-ring slab */
static float slab_normal_err_joint(const Skin* s, const V3* geo, int joint) {
    float sum = 0.0f;
    int cnt = 0;
    for (int ring = joint * SEG_PER_BONE - 1;
         ring <= joint * SEG_PER_BONE + 1; ring++) {
        for (int k = 0; k < RADIAL; k++) {
            int id = ring_vid(ring, k);
            float c = vdot(geo[id], s->nrm[id]);
            if (c > 1.0f) c = 1.0f;
            if (c < -1.0f) c = -1.0f;
            sum += acosf(c) * 180.0f / PI_F;
            cnt++;
        }
    }
    return sum / (float)cnt;
}

/* ---------- skeleton forward kinematics ---------- */

typedef struct { int bone; V3 axis; float deg; } PoseOp;
typedef struct { const char* name; PoseOp op[4]; int nop; } Pose;

/* forward kinematics: world transforms, then skin transforms S = W * B^-1.
 * When bind_inv is NULL (bind-pose construction pass) the world transforms
 * themselves are returned */
static void fk(const Pose* pose, V3 base,
               M4* mats, DQ* dqs, Quat* skin_r,
               const M4* bind_inv_m, const DQ* bind_inv_dq) {
    Quat world_r[MAX_BONES];
    V3 world_t[MAX_BONES];
    Quat local[MAX_BONES];
    for (int i = 0; i < BONES; i++) local[i] = (Quat){ 1, 0, 0, 0 };
    for (int i = 0; i < pose->nop; i++) {
        PoseOp o = pose->op[i];
        Quat q = q_axis(o.axis, o.deg * PI_F / 180.0f);
        local[o.bone] = q_mul(local[o.bone], q);
    }
    for (int i = 0; i < BONES; i++) {
        if (i == 0) {
            world_r[i] = local[i];
            world_t[i] = base;
        } else {
            world_r[i] = q_mul(world_r[i - 1], local[i]);
            V3 off = q_rotate(world_r[i - 1], (V3){ BONE_LEN, 0, 0 });
            world_t[i] = vadd(world_t[i - 1], off);
        }
        M4 wm = m_rigid(world_r[i], world_t[i]);
        DQ wq = dq_rigid(world_r[i], world_t[i]);
        if (bind_inv_m == NULL) {
            mats[i] = wm;
            dqs[i] = wq;
        } else {
            /* S = W B^-1 */
            M4 sm;
            memset(&sm, 0, sizeof(sm));
            for (int r = 0; r < 4; r++)
                for (int c = 0; c < 4; c++)
                    for (int k = 0; k < 4; k++)
                        sm.m[r][c] += wm.m[r][k] * bind_inv_m[i].m[k][c];
            mats[i] = sm;
            dqs[i] = dq_mul(wq, bind_inv_dq[i]);
        }
        skin_r[i] = world_r[i];   /* bind rotations are identity */
    }
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    /* skeleton base is deliberately off-origin so every bone has a nonzero
     * translation and therefore a nonzero dual part */
    V3 base = { 1.0f, 1.0f, 0.0f };

    Mesh mesh;
    mesh_build(&mesh, base);

    /* every pose first rotates the root: then every skin transform carries a
     * nonzero translation, so every per-bone dual quaternion has a nonzero
     * dual part (a zero dual part is defined as an error and must never feed
     * the blender) */
    Pose poses[] = {
        { "bend120_j1",       {{ 0, { 0, 0, 1 }, 10.0f },
                                { 1, { 0, 0, 1 }, 120.0f }}, 2 },
        { "bend-135_j1",      {{ 0, { 0, 0, 1 }, 10.0f },
                                { 1, { 0, 0, 1 }, -135.0f }}, 2 },
        { "bend150_j2",       {{ 0, { 0, 0, 1 }, 10.0f },
                                { 2, { 0, 0, 1 }, 150.0f }}, 2 },
        { "twist150_j1",      {{ 0, { 0, 0, 1 }, 10.0f },
                                { 1, { 1, 0, 0 }, 150.0f }}, 2 },   /* candy wrapper */
        { "twist180_j2",      {{ 0, { 0, 0, 1 }, 10.0f },
                                { 2, { 1, 0, 0 }, 180.0f }}, 2 },
        { "bend+twist_j1",    {{ 0, { 0, 0, 1 }, 10.0f },
                                { 1, { 0, 0, 1 }, 120.0f },
                                { 1, { 1, 0, 0 }, 120.0f }}, 3 },
        { "fold90_j1j2",      {{ 0, { 0, 0, 1 }, 10.0f },
                                { 1, { 0, 0, 1 }, 90.0f },
                                { 2, { 0, 0, 1 }, 90.0f }}, 3 },
        { "counter140_j1j2",  {{ 0, { 0, 0, 1 }, 10.0f },
                                { 1, { 0, 0, 1 }, 140.0f },
                                { 2, { 0, 0, 1 }, -140.0f }}, 3 },
    };
    int nposes = (int)(sizeof(poses) / sizeof(poses[0]));

    /* bind pose: identity local rotations. Bind inverse turns world transforms
     * into skin transforms S = W B^-1, so bind-pose skinning reproduces the mesh */
    Pose bind = { "bind", {{ 0 }}, 0 };
    M4 bind_m[MAX_BONES]; DQ bind_dq[MAX_BONES]; Quat bind_r[MAX_BONES];
    M4 bind_inv_m[MAX_BONES]; DQ bind_inv_dq[MAX_BONES];
    fk(&bind, base, bind_m, bind_dq, bind_r, NULL, NULL);
    for (int i = 0; i < BONES; i++) {
        bind_inv_m[i] = m_inverse_rigid(&bind_m[i]);
        bind_inv_dq[i] = dq_inverse_rigid(&bind_dq[i]);
    }

    /* inverse round-trip must be identity on a sample point */
    V3 sample = { 1.7f, 0.4f, -0.6f };
    float inv_gap = 0.0f;
    for (int i = 0; i < BONES; i++) {
        M4 prod;
        memset(&prod, 0, sizeof(prod));
        for (int r = 0; r < 4; r++)
            for (int c = 0; c < 4; c++)
                for (int k = 0; k < 4; k++)
                    prod.m[r][c] += bind_m[i].m[r][k] * bind_inv_m[i].m[k][c];
        inv_gap += vlen(vsub(m_point(&prod, sample), sample));
        DQ rt = dq_mul(bind_dq[i], bind_inv_dq[i]);
        inv_gap += vlen(vsub(dq_point(&rt, sample), sample));
    }
    check("bind_inverse_roundtrip", inv_gap < 1e-5f);

    /* at bind pose every skin transform is identity: the matrix must be I and
     * the dual quaternion the unit identity (r=1, d=0). LBS through identity
     * matrices reproduces the mesh; DQS at bind pose is the identity too, but
     * its zero dual part correctly triggers the zero-dual guard instead of
     * entering the blender, which is the mandated error path */
    M4 id_m[MAX_BONES]; DQ id_dq[MAX_BONES]; Quat id_r[MAX_BONES];
    fk(&bind, base, id_m, id_dq, id_r, bind_inv_m, bind_inv_dq);
    int bind_is_identity = 1;
    for (int i = 0; i < BONES; i++) {
        for (int r = 0; r < 4; r++)
            for (int c = 0; c < 4; c++) {
                float want = (r == c) ? 1.0f : 0.0f;
                if (fabsf(id_m[i].m[r][c] - want) > 1e-6f) bind_is_identity = 0;
            }
        if (fabsf(id_dq[i].r.w - 1.0f) > 1e-6f ||
            id_dq[i].r.x * id_dq[i].r.x + id_dq[i].r.y * id_dq[i].r.y +
            id_dq[i].r.z * id_dq[i].r.z > 1e-12f ||
            id_dq[i].d.w * id_dq[i].d.w + id_dq[i].d.x * id_dq[i].d.x +
            id_dq[i].d.y * id_dq[i].d.y + id_dq[i].d.z * id_dq[i].d.z > 1e-12f)
            bind_is_identity = 0;
    }
    check("bind_skin_is_identity", bind_is_identity);

    Skin bind_lbs;
    skin_lbs(&mesh, id_m, &bind_lbs);
    float bind_repro_lbs = 0.0f;
    for (int i = 0; i < mesh.nv; i++)
        bind_repro_lbs += vlen(vsub(bind_lbs.pos[i], mesh.v[i].pos));
    check("bind_repro_lbs", bind_repro_lbs < 1e-3f);
    Skin rejected;
    check("dqs_bind_rejects_zero_dual", skin_dqs(&mesh, id_dq, &rejected) == DQ_ERR_ZERO_DUAL);

    /* on single-bone vertices the DQ and matrix skin transforms must agree;
     * use a bent pose so the path through W B^-1 is actually exercised */
    M4 sm0[MAX_BONES]; DQ sdq0[MAX_BONES]; Quat sr0[MAX_BONES];
    fk(&poses[0], base, sm0, sdq0, sr0, bind_inv_m, bind_inv_dq);
    float dq_m_gap = 0.0f;
    for (int i = 0; i < mesh.nv; i++) {
        int single = -1;
        for (int b = 0; b < BONES; b++) if (mesh.v[i].w[b] == 1.0f) single = b;
        if (single < 0) continue;
        dq_m_gap += vlen(vsub(dq_point(&sdq0[single], mesh.v[i].pos),
                              m_point(&sm0[single], mesh.v[i].pos)));
    }
    check("dq_matches_matrix", dq_m_gap < 1e-4f);

    /* blended dual quaternion stays unit: |r| = 1 and r.d = 0 */
    DQ probe;
    float probe_w[MAX_BONES] = { 0.1f, 0.6f, 0.3f };
    check("dq_blend_unit", dq_blend(sdq0, probe_w, BONES, &probe) == DQ_OK &&
          fabsf(q_norm(probe.r) - 1.0f) < 1e-5f &&
          fabsf(q_dot(probe.r, probe.d)) < 1e-5f);

    /* error paths: zero dual part, zero weight sum */
    DQ bad = { { 1, 0, 0, 0 }, { 0, 0, 0, 0 } };
    check("reject_zero_dual", dq_check(&bad) == DQ_ERR_ZERO_DUAL);
    DQ tmp_in[1] = { bad };
    float w_one[1] = { 1.0f };
    DQ tmp_out;
    check("blend_reject_zero_dual", dq_blend(tmp_in, w_one, 1, &tmp_out) == DQ_ERR_ZERO_DUAL);
    float w_zero[MAX_BONES] = { 0 };
    check("blend_reject_zero_weight", dq_blend(sdq0, w_zero, BONES, &tmp_out) == DQ_ERR_ZERO_WEIGHT);

    /* bind-pose slab volumes, identical for both methods */
    float vbind[BONES];
    for (int j = 1; j < BONES; j++) vbind[j] = slab_volume(bind_lbs.pos, j);

    /* rigid pose (root rotates only): the whole mesh undergoes one rigid
     * transform, so every closed slab must keep its exact volume under both
     * methods — this validates the volume measurement itself */
    Pose rigid_pose = { "rigid_root", {{ 0, { 0, 0, 1 }, 47.0f }}, 1 };
    M4 rm[MAX_BONES]; DQ rdq[MAX_BONES]; Quat rr[MAX_BONES];
    fk(&rigid_pose, base, rm, rdq, rr, bind_inv_m, bind_inv_dq);
    Skin rsl, rsd;
    skin_lbs(&mesh, rm, &rsl);
    int rrc = skin_dqs(&mesh, rdq, &rsd);
    int rigid_vol_ok = (rrc == DQ_OK && rsl.nv == rsd.nv);
    for (int j = 1; j < BONES && rigid_vol_ok; j++) {
        float vl = slab_volume(rsl.pos, j), vd = slab_volume(rsd.pos, j);
        if (fabsf(vl - vbind[j]) / vbind[j] > 1e-4f ||
            fabsf(vd - vbind[j]) / vbind[j] > 1e-4f)
            rigid_vol_ok = 0;
    }
    check("rigid_pose_volume_exact", rigid_vol_ok);

    printf("\nper-joint errors over %d extreme poses "
           "(vol: relative volume change, nrm: mean deg of output normal\n"
           "                    vs surface normal recomputed from deformed positions)\n", nposes);
    printf("%-16s j  vol_lbs   vol_dqs   nrm_lbs_deg nrm_dqs_deg\n", "pose");

    int count_match = 1, dqs_wins_twist = 1, dqs_vol_tight = 1;
    double sum_vol_lbs = 0, sum_vol_dqs = 0, sum_nrm_lbs = 0, sum_nrm_dqs = 0;
    int zones = 0;

    for (int p = 0; p < nposes; p++) {
        M4 mats[MAX_BONES]; DQ dqs[MAX_BONES]; Quat skin_r[MAX_BONES];
        fk(&poses[p], base, mats, dqs, skin_r, bind_inv_m, bind_inv_dq);
        Skin sl, sd;
        skin_lbs(&mesh, mats, &sl);
        int rc = skin_dqs(&mesh, dqs, &sd);
        if (rc != DQ_OK) { check("dqs_skin_run", 0); return 1; }

        /* each skinner is scored against the geometry of its own deformed mesh */
        V3 geo_lbs[MAX_VERTS], geo_dqs[MAX_VERTS];
        geo_normals(&mesh, sl.pos, geo_lbs);
        geo_normals(&mesh, sd.pos, geo_dqs);

        /* vertex counts must agree for the same pose */
        if (!(sl.nv == mesh.nv && sd.nv == mesh.nv && sl.nv == sd.nv))
            count_match = 0;

        for (int j = 1; j < BONES; j++) {
            float vl = slab_volume(sl.pos, j);
            float vd = slab_volume(sd.pos, j);
            float el = fabsf(vl - vbind[j]) / vbind[j];
            float ed = fabsf(vd - vbind[j]) / vbind[j];
            float nl = slab_normal_err_joint(&sl, geo_lbs, j);
            float nd = slab_normal_err_joint(&sd, geo_dqs, j);
            printf("%-16s %d  %.5f  %.5f  %8.3f   %8.3f\n",
                   poses[p].name, j, el, ed, nl, nd);
            sum_vol_lbs += el; sum_vol_dqs += ed;
            sum_nrm_lbs += nl; sum_nrm_dqs += nd;
            zones++;
            if (p == 3 && j == 1) {       /* twist150_j1, pinched joint */
                /* DQS blends rigid motions: no candy-wrapper collapse, only a
                 * small residual from rotation interpolation itself */
                if (!(ed < 0.1f && ed * 2.0f < el)) dqs_vol_tight = 0;
                if (!(el > ed)) dqs_wins_twist = 0;
            }
        }
    }

    /* single-bone vertices are rigid: LBS and DQS must coincide */
    int rigid_agree = 1;
    for (int p = 0; p < nposes && rigid_agree; p++) {
        M4 mats[MAX_BONES]; DQ dqs[MAX_BONES]; Quat skin_r[MAX_BONES];
        fk(&poses[p], base, mats, dqs, skin_r, bind_inv_m, bind_inv_dq);
        Skin sl, sd;
        skin_lbs(&mesh, mats, &sl);
        if (skin_dqs(&mesh, dqs, &sd) != DQ_OK) { rigid_agree = 0; break; }
        for (int i = 0; i < mesh.nv; i++) {
            int single = 0;
            for (int b = 0; b < BONES; b++) if (mesh.v[i].w[b] == 1.0f) single = 1;
            if (!single) continue;
            if (vlen(vsub(sl.pos[i], sd.pos[i])) > 1e-4f ||
                vlen(vsub(sl.nrm[i], sd.nrm[i])) > 1e-4f) {
                rigid_agree = 0;
                break;
            }
        }
    }

    printf("\nzone means: vol_lbs=%.5f vol_dqs=%.5f nrm_lbs=%.3fdeg nrm_dqs=%.3fdeg\n",
           sum_vol_lbs / zones, sum_vol_dqs / zones,
           sum_nrm_lbs / zones, sum_nrm_dqs / zones);

    check("vertex_count_match", count_match);
    check("rigid_verts_lbs_eq_dqs", rigid_agree);
    check("dqs_twist_volume_preserved", dqs_vol_tight);
    check("lbs_twist_pinched", dqs_wins_twist);
    check("dqs_mean_vol_below_lbs", sum_vol_dqs < sum_vol_lbs);
    check("dqs_mean_nrm_small", sum_nrm_dqs / zones < 15.0f);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
