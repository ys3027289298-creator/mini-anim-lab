
/* solver: integrated two-bone IK + verlet constraint pipeline.
 *
 * Scene (2D): a pinned shoulder, a two-bone arm (elbow, hand) and a
 * circular torso. The IK target travels through the torso, so IK wants
 * the hand inside the body while the body-collision constraint wants it
 * outside. The two subsystems disagree, and the frame solve order decides
 * who wins:
 *
 *   fixed  (the fix): every frame runs
 *                       integrate -> IK -> distance/collision constraints
 *                     so the collision projection is the last word and the
 *                     hand always finishes on the torso surface.
 *   legacy (buggy):  the IK/constraint order is unspecified and flips per
 *                     frame; on constraint-then-IK frames IK overwrites the
 *                     collision result and the hand ends inside the torso.
 *
 * The analytic IK is ported from ik/, the verlet distance chain from
 * constraint/. Everything is deterministic: fixed dt, no RNG.
 */
#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>

typedef struct { float x, y; } V2;

/* ---- scene ------------------------------------------------------------- */

#define BODY_X 0.6f
#define BODY_Y 1.5f
#define BODY_R 0.9f                 /* torso radius */

#define SH_X (-0.5f)
#define SH_Y 1.5f
#define BONE_L1 1.0f                /* shoulder -> elbow */
#define BONE_L2 1.0f                /* elbow -> hand */

#define R_HAND 0.12f                /* hand particle/mask radius */
#define R_ELBOW 0.10f               /* elbow particle/mask radius */

/* Raster safety margin: the collision limit exceeds the mask radius by
 * more than half a high-tier pixel diagonal, so a particle resting on the
 * limit cannot share a pixel with the body mask at either tier. */
#define SCALE_HI 120
#define SCALE_LO 40
#define CELL_HI (1.0f / SCALE_HI)
#define MARGIN_HAND  (R_HAND  + 0.71f * CELL_HI + 0.001f)
#define MARGIN_ELBOW (R_ELBOW + 0.71f * CELL_HI + 0.001f)

#define VIEW_X0 (-1.5f)
#define VIEW_Y0 (-0.5f)
#define VIEW_W 4.5f
#define VIEW_H 4.0f

#define N_PART 3                    /* shoulder, elbow, hand */
#define MAX_FRAMES 256

/* target loop: centre (0.4, 1.5), radius 0.6 -> always inside the torso,
 * always reachable from the shoulder (distance 0.3 .. 1.5, reach 2.0). */
#define TGT_CX 0.4f
#define TGT_CY 1.5f
#define TGT_R 0.6f

/* ---- analytic two-bone IK, ported from ik/ik.c ------------------------- */

typedef struct {
    int ok;
    V2 joint;
    float a1;
    float a2;
} IkResult;

static float vlen(V2 v) { return sqrtf(v.x * v.x + v.y * v.y); }
static V2 vsub(V2 a, V2 b) { V2 r = { a.x - b.x, a.y - b.y }; return r; }
static V2 vnorm(V2 v) {
    float l = vlen(v);
    V2 r = { 0, 0 };
    if (l > 0) { r.x = v.x / l; r.y = v.y / l; }
    return r;
}

static IkResult ik_two_bone(V2 root, V2 target, float l1, float l2, V2 pole) {
    IkResult r;
    memset(&r, 0, sizeof(r));
    V2 d = vsub(target, root);
    float dist = vlen(d);
    if (l1 <= 0 || l2 <= 0) return r;
    if (dist > l1 + l2 + 1e-5f) return r;
    if (dist < fabsf(l1 - l2) - 1e-5f) return r;
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

/* ---- verlet particle chain, ported from constraint/constraint.c -------- */

typedef struct {
    V2 pos[N_PART];
    V2 prev[N_PART];
    float inv_mass[N_PART];
    int pinned[N_PART];
    float rest[N_PART - 1];
    int iterations;
    float damping;
    float gx, gy;
} ArmChain;

static void chain_init(ArmChain* c, int iterations) {
    memset(c, 0, sizeof(*c));
    c->iterations = iterations;
    c->damping = 0.95f;
    c->gx = 0.0f;
    c->gy = -2.0f;
    V2 p[N_PART] = {
        { SH_X, SH_Y },
        { SH_X + BONE_L1, SH_Y },
        { SH_X + BONE_L1 + BONE_L2, SH_Y },
    };
    for (int i = 0; i < N_PART; i++) {
        c->pos[i] = p[i];
        c->prev[i] = p[i];
        c->inv_mass[i] = 1.0f;
        c->pinned[i] = 0;
    }
    c->rest[0] = BONE_L1;
    c->rest[1] = BONE_L2;
    c->pinned[0] = 1;             /* shoulder anchored */
    c->inv_mass[0] = 0.0f;
}

static void chain_integrate(ArmChain* c, float dt) {
    for (int i = 0; i < N_PART; i++) {
        if (c->pinned[i]) continue;
        float vx = (c->pos[i].x - c->prev[i].x) * c->damping;
        float vy = (c->pos[i].y - c->prev[i].y) * c->damping;
        c->prev[i] = c->pos[i];
        c->pos[i].x += vx + c->gx * dt * dt;
        c->pos[i].y += vy + c->gy * dt * dt;
    }
}

/* IK stage: analytic solve pins elbow and hand positions directly;
 * previous positions are left alone so verlet keeps its velocity. */
static int chain_apply_ik(ArmChain* c, V2 target, V2 pole) {
    IkResult r = ik_two_bone(c->pos[0], target, BONE_L1, BONE_L2, pole);
    if (!r.ok) return 0;
    c->pos[1] = r.joint;
    V2 hand = { target.x, target.y };
    c->pos[2] = hand;
    return 1;
}

/* constraint stage: distance constraints followed by torso collision
 * projection, so collision is the last operation of every iteration. */
static void chain_apply_constraints(ArmChain* c) {
    static const float margin[N_PART] = { 0, MARGIN_ELBOW, MARGIN_HAND };
    V2 body = { BODY_X, BODY_Y };
    for (int it = 0; it < c->iterations; it++) {
        for (int i = 0; i < N_PART - 1; i++) {
            V2* a = &c->pos[i];
            V2* b = &c->pos[i + 1];
            float dx = b->x - a->x, dy = b->y - a->y;
            float d = sqrtf(dx * dx + dy * dy);
            if (d < 1e-6f) continue;
            float diff = (d - c->rest[i]) / d;
            float wa = c->inv_mass[i], wb = c->inv_mass[i + 1];
            float wsum = wa + wb;
            if (wsum <= 0) continue;
            a->x += dx * diff * (wa / wsum);
            a->y += dy * diff * (wa / wsum);
            b->x -= dx * diff * (wb / wsum);
            b->y -= dy * diff * (wb / wsum);
        }
        for (int i = 1; i < N_PART; i++) {
            V2* p = &c->pos[i];
            float dx = p->x - body.x, dy = p->y - body.y;
            float d = sqrtf(dx * dx + dy * dy);
            float lim = BODY_R + margin[i];
            if (d >= lim) continue;
            float nx, ny;
            if (d > 1e-9f) { nx = dx / d; ny = dy / d; }
            else { nx = 1.0f; ny = 0.0f; }
            p->x = body.x + nx * lim;
            p->y = body.y + ny * lim;
        }
    }
}

/* ---- per-frame logging ------------------------------------------------- */

typedef struct {
    int frame;
    char first;           /* 'I' = IK first (fixed order), 'C' = constraints first */
    int ik_ok;
    V2 target;
    V2 shoulder, elbow, hand;
    int hand_pen_lo, arm_pen_lo;
    int hand_pen_hi, arm_pen_hi;
    float clearance;     /* signed hand distance from torso surface */
} FrameLog;

static V2 frame_target(int frame, int n_frames) {
    float theta = 6.2831853f * (float)frame / (float)n_frames;
    V2 t = { TGT_CX + TGT_R * cosf(theta), TGT_CY + TGT_R * sinf(theta) };
    return t;
}

static void run_pipeline(int fixed_order, int iterations, int n_frames,
                         FrameLog* logs) {
    ArmChain c;
    chain_init(&c, iterations);
    V2 pole = { 0.0f, 1.0f };
    float dt = 2.0f / (float)n_frames;       /* one full target loop in 2 s */

    for (int f = 0; f < n_frames; f++) {
        V2 target = frame_target(f, n_frames);
        chain_integrate(&c, dt);

        /* fixed policy: IK always first. legacy policy: order flips per
         * frame, modelling an engine with unspecified subsystem order. */
        int ik_first = fixed_order || (f & 1);
        logs[f].frame = f;
        logs[f].first = ik_first ? 'I' : 'C';
        if (ik_first) {
            logs[f].ik_ok = chain_apply_ik(&c, target, pole);
            chain_apply_constraints(&c);
        } else {
            chain_apply_constraints(&c);
            logs[f].ik_ok = chain_apply_ik(&c, target, pole);
        }

        logs[f].target = target;
        logs[f].shoulder = c.pos[0];
        logs[f].elbow = c.pos[1];
        logs[f].hand = c.pos[2];
        V2 dh = vsub(c.pos[2], (V2){ BODY_X, BODY_Y });
        logs[f].clearance = vlen(dh) - BODY_R;
    }
}

/* ---- rasterization and clipping pixel counts --------------------------- */

typedef struct { int hand_pen, arm_pen; } PenCount;

static PenCount raster(const ArmChain* c, int scale, unsigned char* rgb) {
    int w = (int)(VIEW_W * scale);
    int h = (int)(VIEW_H * scale);
    PenCount pen = { 0, 0 };
    V2 body = { BODY_X, BODY_Y };

    for (int py = 0; py < h; py++) {
        float wy = VIEW_Y0 + ((float)py + 0.5f) / (float)scale;
        for (int px = 0; px < w; px++) {
            float wx = VIEW_X0 + ((float)px + 0.5f) / (float)scale;
            V2 q = { wx, wy };
            int body_m = vlen(vsub(q, body)) <= BODY_R;
            int hand_m = vlen(vsub(q, c->pos[2])) <= R_HAND;
            int elbow_m = vlen(vsub(q, c->pos[1])) <= R_ELBOW;
            int arm_m = hand_m || elbow_m;
            if (hand_m && body_m) pen.hand_pen++;
            if (arm_m && body_m) pen.arm_pen++;
            if (rgb) {
                unsigned char* o = &rgb[3 * (py * w + px)];
                unsigned char rr = 255, gg = 255, bb = 255;
                if (arm_m && body_m) { rr = 230; gg = 0; bb = 0; }       /* clipping: red */
                else if (body_m) { rr = 170; gg = 170; bb = 170; }       /* torso: grey */
                else if (hand_m) { rr = 0; gg = 190; bb = 0; }           /* hand: green */
                else if (elbow_m) { rr = 70; gg = 110; bb = 230; }       /* elbow: blue */
                o[0] = rr; o[1] = gg; o[2] = bb;
            }
        }
    }

    if (rgb) {
        /* bone lines: visualization only, never part of the pixel metric */
        for (int b = 0; b < N_PART - 1; b++) {
            int x0 = (int)((c->pos[b].x - VIEW_X0) * scale);
            int y0 = h - 1 - (int)((c->pos[b].y - VIEW_Y0) * scale);
            int x1 = (int)((c->pos[b + 1].x - VIEW_X0) * scale);
            int y1 = h - 1 - (int)((c->pos[b + 1].y - VIEW_Y0) * scale);
            int dx = abs(x1 - x0), dy = -abs(y1 - y0);
            int sx = x0 < x1 ? 1 : -1, sy = y0 < y1 ? 1 : -1;
            int err = dx + dy, x = x0, y = y0;
            for (;;) {
                if (x >= 0 && x < w && y >= 0 && y < h) {
                    unsigned char* o = &rgb[3 * (y * w + x)];
                    if (o[0] == 255 && o[1] == 255 && o[2] == 255) {
                        o[0] = 30; o[1] = 30; o[2] = 90;
                    }
                }
                if (x == x1 && y == y1) break;
                int e2 = 2 * err;
                if (e2 >= dy) { err += dy; x += sx; }
                if (e2 <= dx) { err += dx; y += sy; }
            }
        }
    }
    return pen;
}

/* replay a logged pipeline through the rasterizer at one tier */
static PenCount log_raster(const FrameLog* logs, int n_frames, int scale,
                           unsigned char* rgb) {
    ArmChain c;
    PenCount total = { 0, 0 };
    for (int f = 0; f < n_frames; f++) {
        c.pos[0] = logs[f].shoulder;
        c.pos[1] = logs[f].elbow;
        c.pos[2] = logs[f].hand;
        PenCount p = raster(&c, scale, rgb);
        total.hand_pen += p.hand_pen;
        total.arm_pen += p.arm_pen;
    }
    return total;
}

static void fill_pen_columns(FrameLog* logs, int n_frames) {
    for (int f = 0; f < n_frames; f++) {
        ArmChain c;
        c.pos[0] = logs[f].shoulder;
        c.pos[1] = logs[f].elbow;
        c.pos[2] = logs[f].hand;
        PenCount lo = raster(&c, SCALE_LO, NULL);
        PenCount hi = raster(&c, SCALE_HI, NULL);
        logs[f].hand_pen_lo = lo.hand_pen;
        logs[f].arm_pen_lo = lo.arm_pen;
        logs[f].hand_pen_hi = hi.hand_pen;
        logs[f].arm_pen_hi = hi.arm_pen;
    }
}

/* ---- output helpers ---------------------------------------------------- */

static void mkdir_quiet(const char* path) {
    mkdir(path, 0755);                 /* EEXIST etc. intentionally ignored */
}

static void write_pose_csv(const char* path, const FrameLog* logs, int n) {
    FILE* fp = fopen(path, "w");
    if (!fp) { perror(path); exit(2); }
    fprintf(fp, "frame,first_stage,target_x,target_y,shoulder_x,shoulder_y,"
                "elbow_x,elbow_y,hand_x,hand_y,hand_body_dist\n");
    for (int i = 0; i < n; i++) {
        const FrameLog* L = &logs[i];
        fprintf(fp, "%d,%c,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
                L->frame, L->first, L->target.x, L->target.y,
                L->shoulder.x, L->shoulder.y, L->elbow.x, L->elbow.y,
                L->hand.x, L->hand.y, L->clearance);
    }
    fclose(fp);
}

static void write_pen_csv(const char* path, const FrameLog* logs, int n) {
    FILE* fp = fopen(path, "w");
    if (!fp) { perror(path); exit(2); }
    fprintf(fp, "frame,hand_pen_lo,arm_pen_lo,hand_pen_hi,arm_pen_hi,"
                "hand_clearance\n");
    for (int i = 0; i < n; i++) {
        const FrameLog* L = &logs[i];
        fprintf(fp, "%d,%d,%d,%d,%d,%.6f\n", L->frame,
                L->hand_pen_lo, L->arm_pen_lo,
                L->hand_pen_hi, L->arm_pen_hi, L->clearance);
    }
    fclose(fp);
}

static void write_ppm(const char* path, int scale, const FrameLog* one) {
    int w = (int)(VIEW_W * scale);
    int h = (int)(VIEW_H * scale);
    unsigned char* rgb = malloc((size_t)w * h * 3);
    if (!rgb) { perror("malloc"); exit(2); }
    ArmChain c;
    c.pos[0] = one->shoulder;
    c.pos[1] = one->elbow;
    c.pos[2] = one->hand;
    raster(&c, scale, rgb);
    FILE* fp = fopen(path, "wb");
    if (!fp) { perror(path); exit(2); }
    fprintf(fp, "P6\n%d %d\n255\n", w, h);
    fwrite(rgb, 3, (size_t)w * h, fp);
    fclose(fp);
    free(rgb);
}

/* ---- timing ------------------------------------------------------------ */

static double now_ms(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec * 1000.0 + ts.tv_nsec / 1.0e6;
}

#define TIMING_REPEATS 12

static double avg_n(const double* a) {
    double s = 0;
    for (int i = 0; i < TIMING_REPEATS; i++) s += a[i];
    return s / TIMING_REPEATS;
}

static double min_n(const double* a) {
    double m = a[0];
    for (int i = 1; i < TIMING_REPEATS; i++) if (a[i] < m) m = a[i];
    return m;
}

static void time_policy(FILE* fp, const char* mode, int fixed_order,
                        int iterations, int n_frames) {
    double sim[TIMING_REPEATS], ras_lo[TIMING_REPEATS], ras_hi[TIMING_REPEATS];
    size_t bytes = sizeof(FrameLog) * (size_t)n_frames;
    FrameLog* logs = malloc(bytes);
    if (!logs) { perror("malloc"); exit(2); }

    for (int r = 0; r < TIMING_REPEATS; r++) {
        double t0 = now_ms();
        run_pipeline(fixed_order, iterations, n_frames, logs);
        double t1 = now_ms();

        PenCount lo = log_raster(logs, n_frames, SCALE_LO, NULL);
        double t2 = now_ms();
        PenCount hi = log_raster(logs, n_frames, SCALE_HI, NULL);
        double t3 = now_ms();

        /* keep the counts live so the rasterizer cannot be optimized away */
        if (lo.hand_pen < 0 || hi.arm_pen < 0) fprintf(stderr, "impossible\n");

        sim[r] = (t1 - t0) / n_frames;
        ras_lo[r] = (t2 - t1) / n_frames;
        ras_hi[r] = (t3 - t2) / n_frames;
    }

    fprintf(fp, "%s,K%d,lo,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f\n", mode, iterations,
            avg_n(sim), min_n(sim), avg_n(ras_lo), min_n(ras_lo),
            avg_n(sim) + avg_n(ras_lo), min_n(sim) + min_n(ras_lo));
    fprintf(fp, "%s,K%d,hi,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f\n", mode, iterations,
            avg_n(sim), min_n(sim), avg_n(ras_hi), min_n(ras_hi),
            avg_n(sim) + avg_n(ras_hi), min_n(sim) + min_n(ras_hi));
    free(logs);
}

/* ---- acceptance run ---------------------------------------------------- */

#define N_FRAMES 40

static int accept_run(const char* outdir) {
    char buf[512];
    mkdir_quiet(outdir);
    snprintf(buf, sizeof(buf), "%s/frames", outdir); mkdir_quiet(buf);
    snprintf(buf, sizeof(buf), "%s/frames/lo", outdir); mkdir_quiet(buf);
    snprintf(buf, sizeof(buf), "%s/frames/hi", outdir); mkdir_quiet(buf);

    FrameLog* fixed = calloc(N_FRAMES, sizeof(FrameLog));
    FrameLog* legacy = calloc(N_FRAMES, sizeof(FrameLog));
    if (!fixed || !legacy) { perror("calloc"); return 2; }

    /* primary pose/penetration artefacts use the middle iteration档 K=8 */
    run_pipeline(1, 8, N_FRAMES, fixed);
    run_pipeline(0, 8, N_FRAMES, legacy);
    fill_pen_columns(fixed, N_FRAMES);
    fill_pen_columns(legacy, N_FRAMES);

    snprintf(buf, sizeof(buf), "%s/pose_fixed.csv", outdir);
    write_pose_csv(buf, fixed, N_FRAMES);
    snprintf(buf, sizeof(buf), "%s/pose_legacy.csv", outdir);
    write_pose_csv(buf, legacy, N_FRAMES);
    snprintf(buf, sizeof(buf), "%s/penetration_fixed.csv", outdir);
    write_pen_csv(buf, fixed, N_FRAMES);
    snprintf(buf, sizeof(buf), "%s/penetration_legacy.csv", outdir);
    write_pen_csv(buf, legacy, N_FRAMES);

    /* PPMs: every frame at the low tier, key frames at the high tier */
    static const int keys[] = { 0, 10, 20, 30, 39 };
    for (int f = 0; f < N_FRAMES; f++) {
        snprintf(buf, sizeof(buf), "%s/frames/lo/fixed_f%02d.ppm", outdir, f);
        write_ppm(buf, SCALE_LO, &fixed[f]);
        snprintf(buf, sizeof(buf), "%s/frames/lo/legacy_f%02d.ppm", outdir, f);
        write_ppm(buf, SCALE_LO, &legacy[f]);
    }
    for (unsigned k = 0; k < sizeof(keys) / sizeof(keys[0]); k++) {
        int f = keys[k];
        snprintf(buf, sizeof(buf), "%s/frames/hi/fixed_f%02d.ppm", outdir, f);
        write_ppm(buf, SCALE_HI, &fixed[f]);
        snprintf(buf, sizeof(buf), "%s/frames/hi/legacy_f%02d.ppm", outdir, f);
        write_ppm(buf, SCALE_HI, &legacy[f]);
    }

    snprintf(buf, sizeof(buf), "%s/timing.csv", outdir);
    FILE* tf = fopen(buf, "w");
    if (!tf) { perror(buf); return 2; }
    fprintf(tf, "mode,quality,raster_tier,sim_avg_ms,sim_min_ms,"
                "raster_avg_ms,raster_min_ms,total_avg_ms,total_min_ms\n");
    time_policy(tf, "fixed", 1, 4, N_FRAMES);
    time_policy(tf, "fixed", 1, 16, N_FRAMES);
    time_policy(tf, "legacy", 0, 4, N_FRAMES);
    time_policy(tf, "legacy", 0, 16, N_FRAMES);
    fclose(tf);

    int sums[4][2] = { 0 };
    for (int f = 0; f < N_FRAMES; f++) {
        sums[0][0] += fixed[f].hand_pen_lo;  sums[0][1] += fixed[f].arm_pen_lo;
        sums[1][0] += fixed[f].hand_pen_hi;  sums[1][1] += fixed[f].arm_pen_hi;
        sums[2][0] += legacy[f].hand_pen_lo; sums[2][1] += legacy[f].arm_pen_lo;
        sums[3][0] += legacy[f].hand_pen_hi; sums[3][1] += legacy[f].arm_pen_hi;
    }
    printf("acceptance artefacts written to %s\n", outdir);
    printf("fixed  clipping pixels: lo hand=%d arm=%d | hi hand=%d arm=%d\n",
           sums[0][0], sums[0][1], sums[1][0], sums[1][1]);
    printf("legacy clipping pixels: lo hand=%d arm=%d | hi hand=%d arm=%d\n",
           sums[2][0], sums[2][1], sums[3][0], sums[3][1]);

    free(fixed);
    free(legacy);
    return 0;
}

/* ---- built-in self-test (no files written) ----------------------------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static int sum_pen(const FrameLog* logs, int n, int hi, int hand) {
    int s = 0;
    for (int i = 0; i < n; i++) {
        if (hi) s += hand ? logs[i].hand_pen_hi : logs[i].arm_pen_hi;
        else s += hand ? logs[i].hand_pen_lo : logs[i].arm_pen_lo;
    }
    return s;
}

static int selftest(void) {
    enum { N = N_FRAMES };
    FrameLog* a = calloc(N, sizeof(FrameLog));
    FrameLog* b = calloc(N, sizeof(FrameLog));
    FrameLog* l = calloc(N, sizeof(FrameLog));

    run_pipeline(1, 8, N, a);
    run_pipeline(1, 8, N, b);               /* determinism rerun */
    run_pipeline(0, 8, N, l);
    fill_pen_columns(a, N);
    fill_pen_columns(b, N);
    fill_pen_columns(l, N);

    check("deterministic_fixed", memcmp(a, b, sizeof(FrameLog) * N) == 0);

    int all_ik = 1, fixed_all_ikfirst = 1, legacy_has_both = 0;
    int legacy_c = 0;
    for (int i = 0; i < N; i++) {
        all_ik = all_ik && a[i].ik_ok && l[i].ik_ok;
        fixed_all_ikfirst = fixed_all_ikfirst && a[i].first == 'I';
        if (l[i].first == 'C') legacy_c++;
    }
    legacy_has_both = legacy_c > 0 && legacy_c < N;
    check("ik_solves_every_frame", all_ik);
    check("fixed_order_is_ik_first", fixed_all_ikfirst);
    check("legacy_order_flips", legacy_has_both);

    check("fixed_zero_hand_pen_lo", sum_pen(a, N, 0, 1) == 0);
    check("fixed_zero_arm_pen_lo", sum_pen(a, N, 0, 0) == 0);
    check("fixed_zero_hand_pen_hi", sum_pen(a, N, 1, 1) == 0);
    check("fixed_zero_arm_pen_hi", sum_pen(a, N, 1, 0) == 0);

    check("legacy_reproduces_pen_lo", sum_pen(l, N, 0, 1) > 0);
    check("legacy_reproduces_pen_hi", sum_pen(l, N, 1, 1) > 0);

    float min_clear_fixed = 1e30f, min_clear_legacy = 1e30f;
    for (int i = 0; i < N; i++) {
        if (a[i].clearance < min_clear_fixed) min_clear_fixed = a[i].clearance;
        if (l[i].clearance < min_clear_legacy) min_clear_legacy = l[i].clearance;
    }
    check("fixed_hand_outside_margin",
          min_clear_fixed >= MARGIN_HAND - 1e-6f);
    check("legacy_hand_enters_body", min_clear_legacy < 0.0f);

    /* zero clipping must hold at every iteration档 (K4 and K16) and both
     * raster tiers: the margin is derived from the pixel size, not tuned per
     * iteration count, so refinement cannot reintroduce a clipping pixel. */
    FrameLog* q4 = calloc(N, sizeof(FrameLog));
    FrameLog* q16 = calloc(N, sizeof(FrameLog));
    run_pipeline(1, 4, N, q4);
    run_pipeline(1, 16, N, q16);
    fill_pen_columns(q4, N);
    fill_pen_columns(q16, N);
    check("fixed_zero_pen_K4",
          sum_pen(q4, N, 0, 1) == 0 && sum_pen(q4, N, 0, 0) == 0 &&
          sum_pen(q4, N, 1, 1) == 0 && sum_pen(q4, N, 1, 0) == 0);
    check("fixed_zero_pen_K16",
          sum_pen(q16, N, 0, 1) == 0 && sum_pen(q16, N, 0, 0) == 0 &&
          sum_pen(q16, N, 1, 1) == 0 && sum_pen(q16, N, 1, 0) == 0);

    /* turning the fixed order off must restore the current poses exactly on
     * every IK-first frame: same IK primitive, same constraints, only the
     * per-frame schedule differs. */
    int revert_exact = 1;
    for (int i = 0; i < N; i++) {
        if (l[i].first != 'I') continue;
        revert_exact = revert_exact &&
            l[i].target.x == a[i].target.x && l[i].target.y == a[i].target.y &&
            l[i].elbow.x == a[i].elbow.x && l[i].elbow.y == a[i].elbow.y &&
            l[i].hand.x == a[i].hand.x && l[i].hand.y == a[i].hand.y &&
            l[i].clearance == a[i].clearance;
    }
    check("legacy_ik_frames_match_fixed", revert_exact);

    /* legacy C frames pin the hand onto the target, which is the mechanism
     * behind their penetration - assert it rather than trusting the plot. */
    int c_pinned = 1;
    for (int i = 0; i < N; i++) {
        if (l[i].first != 'C') continue;
        c_pinned = c_pinned &&
            fabsf(l[i].hand.x - l[i].target.x) < 1e-6f &&
            fabsf(l[i].hand.y - l[i].target.y) < 1e-6f;
    }
    check("legacy_C_frames_pin_target", c_pinned);

    printf("all_ok=%d\n", failures == 0);
    free(a); free(b); free(l); free(q4); free(q16);
    return failures == 0 ? 0 : 1;
}

int main(int argc, char** argv) {
    if (argc == 1) return selftest();
    if (argc == 3 && strcmp(argv[1], "--accept") == 0)
        return accept_run(argv[2]);
    fprintf(stderr, "usage: %s [--accept OUTDIR]\n", argv[0]);
    return 2;
}
