
/* lod: skeletal LOD switching coupled with frustum culling.
 *
 * Two pipelines live in this file:
 *   linked=1  - weld merged LOD tips onto the full-LOD end effectors, rebuild
 *               cull bounds from the active LOD each frame, segment-vs-frustum
 *               test with a radius margin, parent/child visibility coherence,
 *               hysteresis at LOD switch distances.
 *   linked=0  - legacy pipeline untouched: merged tips keep authored chain
 *               length along the first sub-bone direction (no weld), culling
 *               is a point test on the two segment endpoints.
 *
 * Modes:
 *   ./lod                      self-test (see main()), prints name=0/1
 *   ./lod sim <outdir> <0|1>   deterministic headless simulation, writes
 *                              PPM frames, contact sheet, triptychs, CSVs.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_JOINTS 24
#define MAX_CHAINS 8
#define MAX_RENDER 32
#define MAX_SWITCHES 32

#define SCR_W 320
#define SCR_H 240
#define FOCAL 300.0f
#define FOV_HALF 0.26179939f          /* 15 degrees, half of a 30 deg fov */
#define ASPECT (4.0f / 3.0f)
#define CULL_MARGIN 0.12f             /* world units, covers radius + blend */
#define HYST 0.8f                     /* LOD distance hysteresis */
#define D_NEAR 9.0f                   /* below -> LOD0 */
#define D_FAR 13.0f                   /* above -> LOD2, in between LOD1 */

/* ------------------------------------------------------------------ math */

typedef struct { float x, y; } Vec2;

static Vec2 v2(float x, float y) { Vec2 r = { x, y }; return r; }
static Vec2 vadd(Vec2 a, Vec2 b) { return v2(a.x + b.x, a.y + b.y); }
static Vec2 vsub(Vec2 a, Vec2 b) { return v2(a.x - b.x, a.y - b.y); }
static float vlen(Vec2 a) { return sqrtf(a.x * a.x + a.y * a.y); }
static Vec2 vscale(Vec2 a, float k) { return v2(a.x * k, a.y * k); }
static Vec2 vnorm(Vec2 a) {
    float l = vlen(a);
    return l > 1e-9f ? v2(a.x / l, a.y / l) : v2(0.0f, 0.0f);
}
static float vdot(Vec2 a, Vec2 b) { return a.x * b.x + a.y * b.y; }

static Vec2 rotate2(Vec2 p, float a) {
    float c = cosf(a), s = sinf(a);
    return v2(p.x * c - p.y * s, p.x * s + p.y * c);
}

/* closest point on segment [a,b] to q */
static Vec2 seg_closest(Vec2 a, Vec2 b, Vec2 q) {
    Vec2 ab = vsub(b, a);
    float l2 = vdot(ab, ab);
    if (l2 < 1e-12f) return a;
    float t = (vdot(vsub(q, a), ab) / l2);
    if (t < 0.0f) t = 0.0f;
    if (t > 1.0f) t = 1.0f;
    return vadd(a, vscale(ab, t));
}

/* ------------------------------------------------------------- skeleton */

typedef struct {
    int parent;
    Vec2 local;                     /* offset in parent local space */
    float radius;                   /* render thickness of incoming bone */
} JointDef;

typedef struct {
    int count;
    JointDef def[MAX_JOINTS];
    Vec2 world[MAX_JOINTS];
    float ang[MAX_JOINTS];          /* current local angles (animated) */
    const char* name[MAX_JOINTS];
} Skeleton;

enum {
    J_ROOT = 0, J_SPINE0, J_CHEST, J_NECK, J_HEAD,
    J_LSHO, J_LELB, J_LWRI,
    J_RSHO, J_RELB, J_RWRI,
    J_LHIP, J_LKNE, J_LANK,
    J_RHIP, J_RKNE, J_RANK,
    JOINT_N
};

static void skel_build(Skeleton* s) {
    static const JointDef def[JOINT_N] = {
        /* parent, local offset, radius */
        { -1,     { 0.0f, 0.0f}, 0.00f }, /* root */
        { J_ROOT, { 0.0f, 1.2f}, 0.22f }, /* spine0 */
        { J_SPINE0,{0.0f, 1.3f}, 0.24f }, /* chest */
        { J_CHEST,{ 0.0f, 0.9f}, 0.16f }, /* neck */
        { J_NECK, { 0.0f, 0.7f}, 0.30f }, /* head */
        { J_CHEST,{ 0.9f, 0.1f}, 0.16f }, /* l shoulder */
        { J_LSHO, { 1.4f, 0.0f}, 0.13f }, /* l elbow */
        { J_LELB, { 1.3f, 0.0f}, 0.11f }, /* l wrist */
        { J_CHEST,{-0.9f, 0.1f}, 0.16f }, /* r shoulder */
        { J_RSHO, {-1.4f, 0.0f}, 0.13f }, /* r elbow */
        { J_RELB, {-1.3f, 0.0f}, 0.11f }, /* r wrist */
        { J_ROOT, { 0.5f,-0.3f}, 0.20f }, /* l hip */
        { J_LHIP, { 1.4f, 0.0f}, 0.15f }, /* l knee */
        { J_LKNE, { 1.3f, 0.0f}, 0.12f }, /* l ankle */
        { J_ROOT, {-0.5f,-0.3f}, 0.20f }, /* r hip */
        { J_RHIP, {-1.4f, 0.0f}, 0.15f }, /* r knee */
        { J_RKNE, {-1.3f, 0.0f}, 0.12f }, /* r ankle */
    };
    static const char* names[JOINT_N] = {
        "root", "spine0", "chest", "neck", "head",
        "lsho", "lelb", "lwri", "rsho", "relb", "rwri",
        "lhip", "lkne", "lank", "rhip", "rkne", "rank"
    };
    s->count = JOINT_N;
    memcpy(s->def, def, sizeof(def));
    memcpy(s->name, names, sizeof(names));
    memset(s->ang, 0, sizeof(s->ang));
}

/* deterministic procedural pose: local angle per joint at time t */
static void skel_pose(Skeleton* s, float t) {
    for (int i = 0; i < s->count; i++) s->ang[i] = 0.0f;
    s->ang[J_ROOT]   = 0.05f * sinf(2.0f * 3.14159265f * t / 7.0f);
    s->ang[J_SPINE0] = 0.04f * sinf(2.0f * 3.14159265f * t / 5.0f);
    s->ang[J_CHEST]  = 0.05f * sinf(2.0f * 3.14159265f * t / 5.0f + 1.0f);
    s->ang[J_LSHO]   = -0.60f + 0.50f * sinf(2.0f * 3.14159265f * t / 2.2f);
    s->ang[J_LELB]   =  0.90f + 0.40f * sinf(2.0f * 3.14159265f * t / 1.7f + 0.7f);
    s->ang[J_RSHO]   =  0.60f - 0.50f * sinf(2.0f * 3.14159265f * t / 2.2f);
    s->ang[J_RELB]   = -0.90f - 0.40f * sinf(2.0f * 3.14159265f * t / 1.7f + 0.7f);
    s->ang[J_LHIP]   =  0.25f * sinf(2.0f * 3.14159265f * t / 2.4f);
    s->ang[J_LKNE]   =  0.35f + 0.25f * sinf(2.0f * 3.14159265f * t / 1.2f + 1.0f);
    s->ang[J_RHIP]   = -0.25f * sinf(2.0f * 3.14159265f * t / 2.4f);
    s->ang[J_RKNE]   = -0.35f - 0.25f * sinf(2.0f * 3.14159265f * t / 1.2f + 1.0f);
}

static void skel_fk(Skeleton* s) {
    for (int i = 0; i < s->count; i++) {
        Vec2 loc = rotate2(s->def[i].local, s->ang[i]);
        if (s->def[i].parent < 0) {
            s->world[i] = loc;
        } else {
            int p = s->def[i].parent;
            s->world[i] = vadd(s->world[p], loc);
        }
    }
}

/* --------------------------------------------------------- LOD topology */

typedef struct {
    const char* name;
    int anchor;                     /* first joint, stays fixed */
    int joints[5];                  /* chain joints including anchor + tip */
    int n;
    int tip;
} Chain;

typedef struct {
    Chain c[MAX_CHAINS];
    int n;
} ChainTable;

static void chains_build(ChainTable* t) {
    t->n = 0;
    Chain* c;
    c = &t->c[t->n++];
    c->name = "torso"; c->anchor = J_ROOT; c->tip = J_CHEST;
    c->n = 3; c->joints[0] = J_ROOT; c->joints[1] = J_SPINE0; c->joints[2] = J_CHEST;
    c = &t->c[t->n++];
    c->name = "head"; c->anchor = J_CHEST; c->tip = J_HEAD;
    c->n = 3; c->joints[0] = J_CHEST; c->joints[1] = J_NECK; c->joints[2] = J_HEAD;
    c = &t->c[t->n++];
    c->name = "l_arm"; c->anchor = J_LSHO; c->tip = J_LWRI;
    c->n = 3; c->joints[0] = J_LSHO; c->joints[1] = J_LELB; c->joints[2] = J_LWRI;
    c = &t->c[t->n++];
    c->name = "r_arm"; c->anchor = J_RSHO; c->tip = J_RWRI;
    c->n = 3; c->joints[0] = J_RSHO; c->joints[1] = J_RELB; c->joints[2] = J_RWRI;
    c = &t->c[t->n++];
    c->name = "l_leg"; c->anchor = J_LHIP; c->tip = J_LANK;
    c->n = 3; c->joints[0] = J_LHIP; c->joints[1] = J_LKNE; c->joints[2] = J_LANK;
    c = &t->c[t->n++];
    c->name = "r_leg"; c->anchor = J_RHIP; c->tip = J_RANK;
    c->n = 3; c->joints[0] = J_RHIP; c->joints[1] = J_RKNE; c->joints[2] = J_RANK;
}

/* authored straight-line length of a chain (sum of local offset lengths) */
static float seq_length(const Skeleton* s, const int* seq, int n) {
    float sum = 0.0f;
    for (int k = 1; k < n; k++)
        sum += vlen(s->def[seq[k]].local);
    return sum;
}

typedef struct {
    Vec2 p0, p1;
    float radius;
    int cover;                      /* conservative rasterization (merged) */
    int chain;                      /* index into ChainTable, -1 if n/a */
    int tip_joint;                  /* LOD0 joint this segment ends at */
} RenderBone;

typedef struct {
    RenderBone b[MAX_RENDER];
    int n;
} RenderList;

/* full LOD: every authored bone */
static void build_lod0(const Skeleton* s, RenderList* r) {
    r->n = 0;
    for (int i = 1; i < s->count; i++) {
        int k = r->n++;
        int p = s->def[i].parent;
        r->b[k].p0 = s->world[p];
        r->b[k].p1 = s->world[i];
        r->b[k].radius = s->def[i].radius;
        r->b[k].cover = 0;
        r->b[k].chain = -1;
        r->b[k].tip_joint = i;
    }
}

/* merged LOD: one segment per chain.
 * linked  : tip welded to the current LOD0 end-effector position, and the
 *           capsule radius is inflated by the chain sagitta (max distance of
 *           an intermediate joint from the anchor-tip chord) so the straight
 *           segment fully covers the bent LOD0 silhouette - no holes.
 * !linked : legacy tip = anchor + authored chain length along the first
 *           sub-bone's current direction (drops the bend), no inflation. */
static void build_merged(const Skeleton* s, const ChainTable* t,
                         RenderList* r, int lod, int linked) {
    r->n = 0;
    for (int ci = 0; ci < t->n; ci++) {
        const Chain* c = &t->c[ci];
        /* LOD2 folds head into the torso segment */
        if (lod == 2 && (c->tip == J_HEAD)) continue;

        int seq[6] = {0}, seq_n = 0;
        for (int j = 0; j < c->n; j++) seq[seq_n++] = c->joints[j];
        if (lod == 2 && c->tip == J_CHEST) {
            seq[seq_n++] = J_NECK;
            seq[seq_n++] = J_HEAD;
        }

        int k = r->n++;
        int anchor_j = seq[0], tip_j = seq[seq_n - 1];
        Vec2 anchor = s->world[anchor_j];
        Vec2 tip = s->world[tip_j];
        float total = seq_length(s, seq, seq_n);
        float rad = 0.0f;
        for (int j = 1; j < seq_n; j++) {
            float rj = s->def[seq[j]].radius;
            if (rj > rad) rad = rj;
        }

        /* sagitta of the current pose against the anchor-tip chord.
         * each intermediate joint carries its incoming bone's radius, so the
         * envelope of the bent LOD0 chain bulges by (point distance + that
         * radius); the merged capsule must cover the whole envelope or its
         * straight silhouette opens a hole beside the bend. */
        float sag = 0.0f;
        if (linked) {
            for (int j = 1; j < seq_n - 1; j++) {
                Vec2 q = seg_closest(anchor, tip, s->world[seq[j]]);
                float d = vlen(vsub(s->world[seq[j]], q))
                        + s->def[seq[j]].radius;
                if (d > sag) sag = d;
            }
        }

        Vec2 end;
        if (linked) {
            end = tip;
            rad += sag + 1e-4f;
        } else {
            int first = seq[1];
            Vec2 dir = vnorm(vsub(s->world[first], anchor));
            end = vadd(anchor, vscale(dir, total));
        }
        r->b[k].p0 = anchor;
        r->b[k].p1 = end;
        r->b[k].radius = rad;
        r->b[k].cover = linked;
        r->b[k].chain = ci;
        r->b[k].tip_joint = tip_j;
    }
}

static void build_render(const Skeleton* s, const ChainTable* t,
                         int lod, int linked, RenderList* r) {
    if (lod == 0) build_lod0(s, r);
    else build_merged(s, t, r, lod, linked);
}

/* ---------------------------------------------------------------- camera */

typedef struct {
    Vec2 center;
    float half_w, half_h;
    float dist;
} Frustum;

static Frustum frustum_make(Vec2 center, float dist) {
    Frustum f;
    f.center = center;
    f.dist = dist;
    f.half_h = tanf(FOV_HALF) * dist;
    f.half_w = f.half_h * ASPECT;
    return f;
}

static int point_inside(const Frustum* f, Vec2 p) {
    return fabsf(p.x - f->center.x) <= f->half_w &&
           fabsf(p.y - f->center.y) <= f->half_h;
}

/* segment vs frustum AABB expanded by segment radius + cull margin */
static int segment_inside(const Frustum* f, Vec2 a, Vec2 b, float radius) {
    float lim = radius + CULL_MARGIN;
    Vec2 q = seg_closest(a, b, f->center);
    return fabsf(q.x - f->center.x) <= f->half_w + lim &&
           fabsf(q.y - f->center.y) <= f->half_h + lim;
}

/* anchor (attachment end) inside the frustum grown by the cull margin */
static int anchor_inside(const Frustum* f, Vec2 p) {
    return fabsf(p.x - f->center.x) <= f->half_w + CULL_MARGIN &&
           fabsf(p.y - f->center.y) <= f->half_h + CULL_MARGIN;
}

static void cull_list(const RenderList* r, const Frustum* f,
                      int linked, int visible[MAX_RENDER]) {
    for (int i = 0; i < r->n; i++) {
        const RenderBone* b = &r->b[i];
        if (linked) {
            /* keep the bone if the segment crosses the expanded frustum, or
             * if its attachment end is still on screen. The anchor rule is
             * what keeps a limb welded to the body (no detached gap); a bone
             * fully outside with an off-screen anchor is genuinely culled */
            visible[i] = segment_inside(f, b->p0, b->p1, b->radius) ||
                         anchor_inside(f, b->p0);
        } else {
            visible[i] = point_inside(f, b->p0) || point_inside(f, b->p1);
        }
    }
}

/* ------------------------------------------------------------ rasterizer */

typedef struct {
    int w, h;
    unsigned char* px;              /* RGB, row major */
} Framebuf;

static void fb_init(Framebuf* fb, int w, int h) {
    fb->w = w; fb->h = h;
    fb->px = (unsigned char*)calloc((size_t)w * h * 3, 1);
}
static void fb_free(Framebuf* fb) { free(fb->px); fb->px = NULL; }

static void fb_clear(Framebuf* fb) {
    static const unsigned char bg[3] = { 24, 24, 28 };
    for (int i = 0; i < fb->w * fb->h; i++)
        memcpy(fb->px + i * 3, bg, 3);
}

static void fb_set(Framebuf* fb, int x, int y,
                   unsigned char r, unsigned char g, unsigned char b) {
    if (x < 0 || y < 0 || x >= fb->w || y >= fb->h) return;
    unsigned char* p = fb->px + (y * fb->w + x) * 3;
    p[0] = r; p[1] = g; p[2] = b;
}

static void fb_rect(Framebuf* fb, int x0, int y0, int x1, int y1,
                    unsigned char r, unsigned char g, unsigned char b) {
    for (int y = y0; y <= y1; y++)
        for (int x = x0; x <= x1; x++) fb_set(fb, x, y, r, g, b);
}

typedef struct {
    Vec2 center;
    float scale;                    /* pixels per world unit */
} Camera;

static void draw_capsule(Framebuf* fb, const Camera* cam,
                         Vec2 a, Vec2 b, float radius, int cover,
                         unsigned char cr, unsigned char cg, unsigned char cb) {
    /* sub-pixel endpoint positions: projection is kept in continuous
     * coordinates so a welded merged capsule and the LOD0 sub-bones share
     * exactly the same geometry, no integer-endpoint drift between LODs */
    double ax = fb->w * 0.5 + cam->scale * (a.x - cam->center.x);
    double ay = fb->h * 0.5 - cam->scale * (a.y - cam->center.y);
    double bx = fb->w * 0.5 + cam->scale * (b.x - cam->center.x);
    double by = fb->h * 0.5 - cam->scale * (b.y - cam->center.y);    /* cover: round the pixel radius up so merged capsules never under-cover
     * the finer LOD0 silhouette at integer rasterization boundaries */
    double rr = cover ? ceil(radius * cam->scale)
                      : floor(radius * cam->scale);
    double r2 = rr * rr;
    double dx = bx - ax, dy = by - ay;
    double l2 = dx * dx + dy * dy;
    double x0d = fmin(ax, bx) - rr - 1.0, x1d = fmax(ax, bx) + rr + 1.0;
    double y0d = fmin(ay, by) - rr - 1.0, y1d = fmax(ay, by) + rr + 1.0;
    int x0 = (int)x0d, x1 = (int)x1d, y0 = (int)y0d, y1 = (int)y1d;
    double ring2 = (rr - 2.0) * (rr - 2.0);
    for (int y = y0; y <= y1; y++) {
        for (int x = x0; x <= x1; x++) {
            if (x < 0 || y < 0 || x >= fb->w || y >= fb->h) continue;
            double px = x + 0.5, py = y + 0.5;
            double t = 0.0;
            if (l2 > 1e-12) {
                t = ((px - ax) * dx + (py - ay) * dy) / l2;
                if (t < 0.0) t = 0.0;
                if (t > 1.0) t = 1.0;
            }
            double qx = ax + t * dx, qy = ay + t * dy;
            double ex = px - qx, ey = py - qy;
            double d2 = ex * ex + ey * ey;
            if (d2 <= r2) {
                int ring = d2 > ring2;
                fb_set(fb, x, y,
                       ring ? cr / 2 : cr,
                       ring ? cg / 2 : cg,
                       ring ? cb / 2 : cb);
            }
        }
    }
}

static void draw_frustum(Framebuf* fb, const Camera* cam, const Frustum* f) {
    int x0 = (int)(fb->w * 0.5f - cam->scale * f->half_w);
    int x1 = (int)(fb->w * 0.5f + cam->scale * f->half_w);
    int y0 = (int)(fb->h * 0.5f - cam->scale * f->half_h);
    int y1 = (int)(fb->h * 0.5f + cam->scale * f->half_h);
    for (int x = x0; x <= x1; x++) {
        fb_set(fb, x, y0, 120, 120, 130);
        fb_set(fb, x, y1, 120, 120, 130);
    }
    for (int y = y0; y <= y1; y++) {
        fb_set(fb, x0, y, 120, 120, 130);
        fb_set(fb, x1, y, 120, 120, 130);
    }
}

/* palette per chain index; LOD0 bones get a neutral tint */
static void chain_color(int chain, unsigned char* r,
                        unsigned char* g, unsigned char* b) {
    static const unsigned char pal[][3] = {
        {210, 170, 120},   /* torso */
        {230, 220, 180},   /* head */
        {120, 200, 230},   /* l arm */
        {90, 160, 220},    /* r arm */
        {160, 220, 140},   /* l leg */
        {120, 190, 110},   /* r leg */
    };
    if (chain < 0) { *r = 180; *g = 180; *b = 190; return; }
    *r = pal[chain % 6][0]; *g = pal[chain % 6][1]; *b = pal[chain % 6][2];
}

static void draw_badge(Framebuf* fb, int lod, int mode) {
    /* mode: 0 legacy (red strip), 1 linked (green), 2 ideal (blue) */
    static const unsigned char col[3][3] = {
        {200, 70, 70}, {70, 200, 110}, {90, 130, 230}
    };
    fb_rect(fb, fb->w - 10, 2, fb->w - 3, 5,
            col[mode][0], col[mode][1], col[mode][2]);
    for (int i = 0; i < 3; i++) {
        int on = i == lod;
        fb_rect(fb, 4 + i * 8, 4, 9 + i * 8, 9,
                on ? 235 : 70, on ? 235 : 70, on ? 235 : 80);
    }
}

static void render_scene(Framebuf* fb, const RenderList* r,
                         const int visible[MAX_RENDER],
                         const Frustum* f, int lod, int mode,
                         int all_visible) {
    Camera cam;
    cam.center = f->center;
    cam.scale = (fb->w == SCR_W ? FOCAL : FOCAL * 0.5f) / f->dist;
    fb_clear(fb);
    draw_frustum(fb, &cam, f);
    for (int i = 0; i < r->n; i++) {
        if (!all_visible && !visible[i]) continue;
        unsigned char cr, cg, cb;
        chain_color(r->b[i].chain, &cr, &cg, &cb);
        draw_capsule(fb, &cam, r->b[i].p0, r->b[i].p1, r->b[i].radius,
                     r->b[i].cover, cr, cg, cb);
    }
    draw_badge(fb, lod, mode);
}

static int write_ppm(const char* path, const Framebuf* fb) {
    FILE* fp = fopen(path, "wb");
    if (!fp) { perror(path); return -1; }
    fprintf(fp, "P6\n%d %d\n255\n", fb->w, fb->h);
    fwrite(fb->px, 1, (size_t)fb->w * fb->h * 3, fp);
    fclose(fp);
    return 0;
}

/* triptych: legacy | linked | ideal, side by side with dark separators */
static int write_triptych(const char* path,
                          const Framebuf* a, const Framebuf* b,
                          const Framebuf* c) {
    int gap = 6;
    Framebuf out;
    fb_init(&out, a->w * 3 + gap * 2, a->h);
    fb_clear(&out);
    for (int y = 0; y < a->h; y++) {
        memcpy(out.px + (y * out.w) * 3,
               a->px + y * a->w * 3, (size_t)a->w * 3);
        memcpy(out.px + (y * out.w + a->w + gap) * 3,
               b->px + y * a->w * 3, (size_t)a->w * 3);
        memcpy(out.px + (y * out.w + 2 * (a->w + gap)) * 3,
               c->px + y * a->w * 3, (size_t)a->w * 3);
    }
    int rc = write_ppm(path, &out);
    fb_free(&out);
    return rc;
}

/* ------------------------------------------------------ LOD with hysteresis */

static int lod_target(float dist) {
    if (dist < D_NEAR) return 0;
    if (dist > D_FAR) return 2;
    return 1;
}

/* advances LOD state one frame; returns 1 if the level changed */
static int lod_step(int* cur, float dist) {
    int want = lod_target(dist);
    int before = *cur;
    while (*cur != want) {
        int up = want > *cur;
        float boundary = up ? (*cur == 0 ? D_NEAR : D_FAR)
                            : (*cur == 2 ? D_FAR : D_NEAR);
        if (up) {
            if (dist <= boundary + HYST) break;
            (*cur)++;
        } else {
            if (dist >= boundary - HYST) break;
            (*cur)--;
        }
    }
    return *cur != before;
}

/* ----------------------------------------------------------- simulation */

typedef struct {
    int frame;
    float t, dist, camx;
    int lod_from, lod_to;
    int cov_before, cov_after, missing_after;
    float err_legacy, err_linked;
} SwitchRecord;

typedef struct {
    int frames;
    int lod;
    SwitchRecord sw[MAX_SWITCHES];
    int nsw;
    int worst_linked_missing;      /* vs culled LOD0 truth, linked run */
    int worst_legacy_missing;      /* vs culled LOD0 truth, legacy run */
    int legacy_hole_frames;
} Sim;

static float sim_distance(float t) {
    return 11.0f + 5.0f * sinf(2.0f * 3.14159265f * t / 8.0f);
}
static float sim_camx(float t) {
    return 2.6f * sinf(2.0f * 3.14159265f * t / 5.3f);
}

/* max end-effector error of a merged LOD against LOD0 joint tips */
static void tip_errors(const Skeleton* s, const ChainTable* t, int lod,
                       float* max_legacy, float* max_linked) {
    RenderList rl, ri;
    build_render(s, t, lod, 0, &rl);
    build_render(s, t, lod, 1, &ri);
    float el = 0.0f, ew = 0.0f;
    for (int i = 0; i < rl.n; i++) {
        Vec2 ref = s->world[rl.b[i].tip_joint];
        float dl = vlen(vsub(rl.b[i].p1, ref));
        float dw = vlen(vsub(ri.b[i].p1, ref));
        if (dl > el) el = dl;
        if (dw > ew) ew = dw;
    }
    *max_legacy = el;
    *max_linked = ew;
}

static int coverage(const Framebuf* fb) {
    int n = 0;
    for (int i = 0; i < fb->w * fb->h; i++) {
        const unsigned char* p = fb->px + i * 3;
        if (!(p[0] == 24 && p[1] == 24 && p[2] == 28)) n++;
    }
    return n;
}

/* pixels painted by ideal but missing from the tested framebuffer */
static int missing_pixels(const Framebuf* ideal, const Framebuf* test) {
    int n = 0;
    for (int i = 0; i < ideal->w * ideal->h; i++) {
        const unsigned char* a = ideal->px + i * 3;
        const unsigned char* b = test->px + i * 3;
        int aChar = !(a[0] == 24 && a[1] == 24 && a[2] == 28);
        int bChar = !(b[0] == 24 && b[1] == 24 && b[2] == 28);
        if (aChar && !bChar) n++;
    }
    return n;
}

/* run the full simulation. when outdir is NULL nothing is written, used by
 * the in-memory self-test. linked selects the active pipeline. */
static int simulate(int linked, int nframes, const char* outdir, Sim* sim) {
    Skeleton sk;
    ChainTable ct;
    skel_build(&sk);
    chains_build(&ct);
    sim->frames = nframes;
    /* the camera starts at 11 units, which selects LOD1; seed the state
     * there so the first frame is not recorded as a spurious switch */
    sim->lod = 1;
    sim->nsw = 0;

    FILE* fmetrics = NULL;
    FILE* ferrors = NULL;
    FILE* fswitch = NULL;
    if (outdir) {
        char path[512];
        snprintf(path, sizeof(path), "%s/frame_metrics.csv", outdir);
        fmetrics = fopen(path, "w");
        snprintf(path, sizeof(path), "%s/error_table.csv", outdir);
        ferrors = fopen(path, "w");
        snprintf(path, sizeof(path), "%s/switch_log.csv", outdir);
        fswitch = fopen(path, "w");
        if (!fmetrics || !ferrors || !fswitch) return -1;
        fprintf(fmetrics,
                "frame,t,dist,camx,lod,visible,culled,coverage,truth_px,"
                "missing_vs_truth,legacy_missing_px,max_tip_err_legacy,"
                "max_tip_err_linked\n");
        fprintf(ferrors,
                "frame,t,active_lod,table_lod,chain,ref_x,ref_y,"
                "legacy_x,legacy_y,linked_x,linked_y,err_legacy,err_linked\n");
        fprintf(fswitch,
                "frame,t,dist,camx,lod_from,lod_to,direction,"
                "coverage_before,coverage_after,missing_after_truth,"
                "max_err_legacy,max_err_linked,hysteresis\n");
    }

    /* fb        : the active pipeline's output for this run
     * fb_truth  : LOD0 full skeleton through the SAME linked culling (the
     *             no-hole reference: what the character actually looks like)
     * fb_legacy / fb_linked : both merged-LOD pipelines, for triptychs */
    Framebuf fb, fb_truth, fb_legacy, fb_linked;
    fb_init(&fb, SCR_W, SCR_H);
    fb_init(&fb_truth, SCR_W, SCR_H);
    fb_init(&fb_legacy, SCR_W, SCR_H);
    fb_init(&fb_linked, SCR_W, SCR_H);

    /* contact sheet: 6 x 4 half-resolution thumbs */
    Framebuf sheet;
    int thumb_w = SCR_W / 2, thumb_h = SCR_H / 2;
    int cols = 6, rows = 4;
    if (outdir) fb_init(&sheet, thumb_w * cols, thumb_h * rows);
    int sheet_slots = cols * rows;

    sim->worst_linked_missing = 0;
    sim->worst_legacy_missing = 0;
    sim->legacy_hole_frames = 0;
    int prev_cov = -1;
    float dt = 1.0f / 30.0f;

    for (int frame = 0; frame < nframes; frame++) {
        float t = frame * dt;
        float dist = sim_distance(t);
        Vec2 cam = v2(sim_camx(t), 0.0f);
        Frustum fr = frustum_make(cam, dist);

        skel_pose(&sk, t);
        skel_fk(&sk);

        int lod_before = sim->lod;
        int switched = lod_step(&sim->lod, dist);

        /* active-LOD render lists for both pipelines */
        RenderList r_active, r_legacy, r_linked, r_truth;
        build_render(&sk, &ct, sim->lod, linked, &r_active);
        build_render(&sk, &ct, sim->lod, 0, &r_legacy);
        build_render(&sk, &ct, sim->lod, 1, &r_linked);
        build_render(&sk, &ct, 0, 1, &r_truth);

        int vis[MAX_RENDER], vis_legacy[MAX_RENDER];
        int vis_linked[MAX_RENDER], vis_truth[MAX_RENDER];
        cull_list(&r_active, &fr, linked, vis);
        cull_list(&r_legacy, &fr, 0, vis_legacy);
        cull_list(&r_linked, &fr, 1, vis_linked);
        cull_list(&r_truth, &fr, 1, vis_truth);

        int nvis = 0;
        for (int i = 0; i < r_active.n; i++) if (vis[i]) nvis++;

        render_scene(&fb, &r_active, vis, &fr, sim->lod,
                     linked ? 1 : 0, 0);
        render_scene(&fb_truth, &r_truth, vis_truth, &fr, 0, 2, 0);
        render_scene(&fb_legacy, &r_legacy, vis_legacy, &fr, sim->lod, 0, 0);
        render_scene(&fb_linked, &r_linked, vis_linked, &fr, sim->lod, 1, 0);

        int cov_truth = coverage(&fb_truth);
        int miss_linked = missing_pixels(&fb_truth, &fb_linked);
        int miss_legacy = missing_pixels(&fb_truth, &fb_legacy);
        int miss_active = linked ? miss_linked : miss_legacy;
        if (miss_linked > sim->worst_linked_missing)
            sim->worst_linked_missing = miss_linked;
        if (miss_legacy > sim->worst_legacy_missing)
            sim->worst_legacy_missing = miss_legacy;
        if (miss_legacy > 0) sim->legacy_hole_frames++;

        float el1, ew1, el2, ew2;
        tip_errors(&sk, &ct, 1, &el1, &ew1);
        tip_errors(&sk, &ct, 2, &el2, &ew2);
        float el = el1 > el2 ? el1 : el2;
        float ew = ew1 > ew2 ? ew1 : ew2;

        if (fmetrics) {
            fprintf(fmetrics,
                    "%d,%.6f,%.6f,%.6f,%d,%d,%d,%d,%d,%d,%d,%.6f,%.6f\n",
                    frame, t, dist, cam.x, sim->lod, nvis,
                    r_active.n - nvis, coverage(&fb), cov_truth,
                    miss_active, miss_legacy, el, ew);
        }

        /* error table: every 12th frame and every switch frame */
        if (ferrors && (frame % 12 == 0 || switched)) {
            for (int table_lod = 1; table_lod <= 2; table_lod++) {
                RenderList rl, ri;
                build_render(&sk, &ct, table_lod, 0, &rl);
                build_render(&sk, &ct, table_lod, 1, &ri);
                for (int i = 0; i < rl.n; i++) {
                    Vec2 ref = sk.world[rl.b[i].tip_joint];
                    float dl = vlen(vsub(rl.b[i].p1, ref));
                    float dw = vlen(vsub(ri.b[i].p1, ref));
                    const char* nm = rl.b[i].chain >= 0
                        ? ct.c[rl.b[i].chain].name : "bone";
                    fprintf(ferrors,
                            "%d,%.6f,%d,%d,%s,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,"
                            "%.6f,%.6f\n",
                            frame, t, sim->lod, table_lod, nm,
                            ref.x, ref.y,
                            rl.b[i].p1.x, rl.b[i].p1.y,
                            ri.b[i].p1.x, ri.b[i].p1.y, dl, dw);
                }
            }
        }

        if (switched && sim->nsw < MAX_SWITCHES) {
            SwitchRecord* rec = &sim->sw[sim->nsw++];
            rec->frame = frame; rec->t = t; rec->dist = dist; rec->camx = cam.x;
            rec->lod_from = lod_before; rec->lod_to = sim->lod;
            rec->cov_before = prev_cov;
            rec->cov_after = coverage(&fb);
            rec->missing_after = miss_active;
            rec->err_legacy = el; rec->err_linked = ew;
            if (fswitch) {
                fprintf(fswitch,
                        "%d,%.6f,%.6f,%.6f,%d,%d,%s,%d,%d,%d,%.6f,%.6f,%.2f\n",
                        frame, t, dist, cam.x, lod_before, sim->lod,
                        sim->lod > lod_before ? "up" : "down",
                        rec->cov_before, rec->cov_after, rec->missing_after,
                        el, ew, HYST);
            }
            if (outdir) {
                char path[512];
                snprintf(path, sizeof(path),
                         "%s/triptych_%02d_f%04d.ppm",
                         outdir, sim->nsw - 1, frame);
                write_triptych(path, &fb_legacy, &fb_linked, &fb_truth);
            }
        }
        prev_cov = coverage(&fb);

        if (outdir) {
            char path[512];
            if (frame % 12 == 0 || switched) {
                snprintf(path, sizeof(path),
                         "%s/frame_%04d.ppm", outdir, frame);
                write_ppm(path, &fb);
            }
            int slot = (int)((long)frame * sheet_slots / nframes);
            if (slot < sheet_slots) {
                Framebuf thumb;
                fb_init(&thumb, thumb_w, thumb_h);
                render_scene(&thumb, &r_active, vis, &fr, sim->lod,
                             linked ? 1 : 0, 0);
                int ox = (slot % cols) * thumb_w;
                int oy = (slot / cols) * thumb_h;
                for (int y = 0; y < thumb_h; y++)
                    memcpy(sheet.px + ((oy + y) * sheet.w + ox) * 3,
                           thumb.px + y * thumb_w * 3,
                           (size_t)thumb_w * 3);
                fb_free(&thumb);
            }
        }
    }

    if (outdir) {
        char path[512];
        snprintf(path, sizeof(path), "%s/contact_sheet.ppm", outdir);
        write_ppm(path, &sheet);
        fb_free(&sheet);

        snprintf(path, sizeof(path), "%s/run_meta.txt", outdir);
        FILE* meta = fopen(path, "w");
        if (meta) {
            fprintf(meta,
                    "linked=%d\nframes=%d\nfov_deg=30\nwidth=%d height=%d\n"
                    "cull_margin=%.3f\nhysteresis=%.2f\n"
                    "thresholds=%.1f,%.1f\n"
                    "worst_linked_missing_px=%d\n"
                    "worst_legacy_missing_px=%d\nlegacy_hole_frames=%d\n"
                    "switches=%d\n",
                    linked, nframes, SCR_W, SCR_H, CULL_MARGIN, HYST,
                    D_NEAR, D_FAR, sim->worst_linked_missing,
                    sim->worst_legacy_missing, sim->legacy_hole_frames,
                    sim->nsw);
            fclose(meta);
        }
    }

    if (fmetrics) fclose(fmetrics);
    if (ferrors) fclose(ferrors);
    if (fswitch) fclose(fswitch);
    fb_free(&fb); fb_free(&fb_truth);
    fb_free(&fb_legacy); fb_free(&fb_linked);

    sim->frames = nframes;
    return 0;
}

/* ------------------------------------------------------------------ test */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

#ifndef LOD_NO_MAIN
int main(int argc, char** argv) {
    if (argc >= 2 && strcmp(argv[1], "sim") == 0) {
        if (argc < 4) {
            fprintf(stderr, "usage: %s sim <outdir> <linked 0|1> [frames]\n",
                    argv[0]);
            return 2;
        }
        const char* outdir = argv[2];
        int linked = atoi(argv[3]);
        int nframes = argc >= 5 ? atoi(argv[4]) : 240;
        Sim sim;
        int rc = simulate(linked, nframes, outdir, &sim);
        if (rc < 0) return 1;
        printf("sim_done linked=%d frames=%d switches=%d\n",
               linked, nframes, sim.nsw);
        return 0;
    }

    /* ---- self-test ---- */
    Skeleton sk;
    ChainTable ct;
    skel_build(&sk);
    chains_build(&ct);

    /* forward kinematics on a known straight pose (all joint angles zero) */
    memset(sk.ang, 0, sizeof(sk.ang));
    skel_fk(&sk);
    check("fk_vertical",
          fabsf(sk.world[J_CHEST].y - 2.5f) < 1e-5f &&
          fabsf(sk.world[J_CHEST].x) < 1e-5f);

    /* welded merged tips coincide with LOD0 end effectors */
    RenderList rw, r0;
    build_render(&sk, &ct, 1, 1, &rw);
    build_render(&sk, &ct, 0, 1, &r0);
    int weld_ok = 1;
    for (int i = 0; i < rw.n; i++)
        weld_ok = weld_ok &&
                  vlen(vsub(rw.b[i].p1, sk.world[rw.b[i].tip_joint])) == 0.0f;
    check("weld_tip_exact", weld_ok);

    /* a bent chain: legacy tip must land short of the real end effector */
    sk.ang[J_LSHO] = 0.0f;
    sk.ang[J_LELB] = 1.5707963f;              /* 90 degree forearm */
    skel_fk(&sk);
    RenderList rl;
    build_render(&sk, &ct, 1, 0, &rl);
    int arm = -1;
    for (int i = 0; i < rl.n; i++)
        if (rl.b[i].tip_joint == J_LWRI) arm = i;
    float legacy_err = vlen(vsub(rl.b[arm].p1, sk.world[J_LWRI]));
    check("legacy_tip_short_bent", legacy_err > 0.5f);
    build_render(&sk, &ct, 1, 1, &rw);
    for (int i = 0; i < rw.n; i++)
        if (rw.b[i].tip_joint == J_LWRI)
            check("weld_tip_zero_bent",
                  vlen(vsub(rw.b[i].p1, sk.world[J_LWRI])) == 0.0f);

    /* segment crossing the frustum with both endpoints outside */
    Frustum f = frustum_make(v2(0.0f, 0.0f), 10.0f);
    Vec2 a = v2(-f.half_w - 2.0f, 0.0f);
    Vec2 b = v2( f.half_w + 2.0f, 0.0f);
    int vis_legacy = point_inside(&f, a) || point_inside(&f, b);
    int vis_linked = segment_inside(&f, a, b, 0.15f);
    check("legacy_culls_crossing", !vis_legacy);
    check("linked_keeps_crossing", vis_linked);

    /* anchor coherence: a segment reaching off-screen but attached on-screen
     * is kept (no detached hole); a segment fully off-screen is culled */
    RenderList rc = {0};
    rc.n = 2;
    /* anchor just outside the raw frustum but inside the cull margin, tip
     * far outside: linkage keeps it, raw endpoint test drops it */
    rc.b[0].p0 = v2(f.half_w + CULL_MARGIN * 0.5f, 0.0f);
    rc.b[0].p1 = v2(f.half_w + 3.0f, 0.0f);
    rc.b[0].radius = 0.15f;
    rc.b[1].p0 = v2(f.half_w + 5.0f, 0.0f);            /* anchor outside */
    rc.b[1].p1 = v2(f.half_w + 7.0f, 0.0f);
    rc.b[1].radius = 0.15f;
    int vis[MAX_RENDER];
    cull_list(&rc, &f, 1, vis);
    check("anchor_coherence_kept", vis[0] == 1);
    check("offscreen_bone_culled", vis[1] == 0);
    cull_list(&rc, &f, 0, vis);
    check("legacy_drops_anchored", vis[0] == 0);

    /* hysteresis: hovering just past a threshold must not flicker */
    int lod = 0;
    lod_step(&lod, D_NEAR + 0.5f);          /* inside hysteresis band */
    check("hyst_hold_near", lod == 0);
    lod_step(&lod, D_NEAR + HYST + 0.1f);   /* past the band */
    check("hyst_switch", lod == 1);
    lod_step(&lod, D_NEAR + 0.1f);          /* back inside band */
    check("hyst_no_flick_back", lod == 1);

    /* deterministic in-memory full runs: linked has zero missing pixels
     * against the culled LOD0 truth, while the legacy pipeline demonstrably
     * opens holes */
    Sim sim_l, sim_o;
    check("sim_linked_rc", simulate(1, 240, NULL, &sim_l) == 0);
    check("sim_legacy_rc", simulate(0, 240, NULL, &sim_o) == 0);
    check("linked_no_holes",
          sim_l.worst_linked_missing == 0 &&
          sim_l.worst_legacy_missing > 0);
    check("legacy_reproduces_holes", sim_o.legacy_hole_frames > 0);
    check("switches_recorded",
          sim_l.nsw >= 4 && sim_l.nsw == sim_o.nsw);

    /* switch moments: no missing pixels and zero welded tip error */
    int switch_clean = 1;
    for (int i = 0; i < sim_l.nsw; i++)
        switch_clean = switch_clean &&
                       sim_l.sw[i].missing_after == 0 &&
                       sim_l.sw[i].err_linked == 0.0f;
    check("switch_moments_clean", switch_clean);

    /* pixel determinism: render the same frame twice, memcmp */
    skel_pose(&sk, 1.234f);
    skel_fk(&sk);
    Frustum f2 = frustum_make(v2(0.3f, 0.0f), 11.0f);
    RenderList ra, rb;
    build_render(&sk, &ct, 1, 0, &ra);
    build_render(&sk, &ct, 1, 0, &rb);
    int va[MAX_RENDER], vb[MAX_RENDER];
    cull_list(&ra, &f2, 0, va);
    cull_list(&rb, &f2, 0, vb);
    Framebuf fa, fb2;
    fb_init(&fa, SCR_W, SCR_H);
    fb_init(&fb2, SCR_W, SCR_H);
    render_scene(&fa, &ra, va, &f2, 1, 0, 0);
    render_scene(&fb2, &rb, vb, &f2, 1, 0, 0);
    check("render_deterministic",
          memcmp(fa.px, fb2.px, (size_t)SCR_W * SCR_H * 3) == 0);
    fb_free(&fa); fb_free(&fb2);

    /* linkage off must be pixel-identical to the all-legacy pipeline:
     * every merged list built with linked=0 has cover==0 and the legacy
     * endpoint cull path, so no new raster behavior can leak in */
    int legacy_pure = 1;
    for (int lodx = 1; lodx <= 2 && legacy_pure; lodx++) {
        RenderList rp;
        build_render(&sk, &ct, lodx, 0, &rp);
        for (int i = 0; i < rp.n; i++)
            if (rp.b[i].cover) legacy_pure = 0;
    }
    check("legacy_path_untouched", legacy_pure);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
#endif /* LOD_NO_MAIN */
