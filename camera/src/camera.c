
/* camera: joint accel/angular-velocity stress evaluation feeding an endpoint
   braking damper; motion-blur frame capture across three sensitivity
   presets; sudden-stop micro bounce with settle-time comparison of
   auto-damping vs manual counter-input modes. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define PI 3.14159265358979f
#define FPS 60.0f
#define DT (1.0f / FPS)
#define HALF_FRAME (0.5f * DT)

#define PX_PER_UNIT 40.0f    /* screen px per world unit */
#define FOCAL_PX 220.0f      /* screen px per radian of yaw */
#define IMG_W 320
#define IMG_H 90
#define SUB 12               /* exposure sub-samples for motion blur */
#define T_EXP (1.0f / 120.0f)

#define A_REF 60.0f          /* stress normalizers: linear accel, angular vel */
#define W_REF 6.0f

#define SETTLE_PX 0.5f       /* bounce settled below this screen offset */
#define TGT_POS 8.0f
#define TGT_YAW 1.1f

#define MOUNT_WN 16.0f       /* mount resonance rad/s (~2.5 Hz wobble) */
#define MOUNT_ZETA 0.06f     /* underdamped: this is what bounces */
#define MOUNT_COUP 0.5f      /* inertial coupling of base accel into head */
#define AUTO_MOUNT_C 22.0f   /* auto mode: extra mount damping near endpoint */

#define USER_DELAY 0.10f     /* manual mode: reaction latency before bracing */
#define USER_G 100.0f        /* brace stiffness */
#define USER_GV 14.0f        /* brace viscous damping */
#define USER_UMAX 40.0f      /* brace force limit */

/* hard end-stop: the base slams to rest within one frame when it reaches
   the endpoint still carrying speed; the decel spike excites the mount */
#define STOP_REM 0.2f
#define STOP_VEL 4.0f

static float clampf(float v, float lo, float hi) {
    return v < lo ? lo : v > hi ? hi : v;
}

/* static test scene: stripes + sharp edges along world x, in [0,1] */
static float scene(float wx) {
    float v = 0.5f
        + 0.22f * sinf(2.0f * PI * wx / 0.35f)
        + 0.16f * sinf(2.0f * PI * wx / 0.11f + 1.3f)
        + 0.25f * tanhf((wx - 2.0f) * 60.0f)
        - 0.25f * tanhf((wx - 5.5f) * 60.0f)
        + 0.20f * tanhf((wx - 7.2f) * 60.0f);
    return clampf(v, 0.0f, 1.0f);
}

/* joint evaluation of linear acceleration and angular velocity */
static float joint_stress(float acc, float yaw_vel) {
    return fabsf(acc) / A_REF + fabsf(yaw_vel) / W_REF;
}

typedef struct {
    const char* name;
    float sens;              /* input sensitivity: scales accel authority */
    float kp, kd;            /* chase controller gains */
    float a_max, w_max;      /* accel clamps (before sens), linear / yaw */
    float brake_radius;      /* endpoint distance where the brake ramps in */
    float brake_yaw_radius;
    float brake_c0;          /* base brake damping */
    float brake_k;           /* stress gain on the brake */
    int brake_on;
} Tune;

/* brake damping coefficient for a given joint stress */
static float brake_gain(const Tune* tu, float stress) {
    return tu->brake_c0 * (1.0f + tu->brake_k * stress);
}

typedef struct {
    float pos, vel, acc;
    float yaw, yvel, yacc;
    float stress;
} Cam;

static void cam_step(Cam* c, const Tune* tu, float dt) {
    float amax = tu->a_max * tu->sens;
    float wmax = tu->w_max * tu->sens;
    float acc = tu->kp * (TGT_POS - c->pos) - tu->kd * c->vel;
    float yacc = tu->kp * (TGT_YAW - c->yaw) - tu->kd * c->yvel;
    acc = clampf(acc, -amax, amax);
    yacc = clampf(yacc, -wmax, wmax);

    float stress = joint_stress(acc, c->yvel);

    /* endpoint braking damper, ramps in over the last stretch */
    if (tu->brake_on) {
        float ramp = 1.0f - clampf(fabsf(TGT_POS - c->pos) / tu->brake_radius, 0.0f, 1.0f);
        float yramp = 1.0f - clampf(fabsf(TGT_YAW - c->yaw) / tu->brake_yaw_radius, 0.0f, 1.0f);
        float cg = brake_gain(tu, stress);
        acc -= cg * ramp * c->vel;
        yacc -= cg * yramp * c->yvel;
    }
    c->stress = stress;
    c->acc = acc;
    c->yacc = yacc;
    c->vel += acc * dt;
    c->pos += c->vel * dt;
    c->yvel += yacc * dt;
    c->yaw += c->yvel * dt;
}

/* camera head on an underdamped mount: sudden base stops excite it */
typedef struct {
    float off, vel;          /* head offset from base, world units */
    float extra_c;           /* mode-dependent extra damping */
} Mount;

static void mount_step(Mount* m, float base_acc, float user_u, float dt) {
    float c = 2.0f * MOUNT_ZETA * MOUNT_WN + m->extra_c;
    float a = -MOUNT_WN * MOUNT_WN * m->off - c * m->vel
            - MOUNT_COUP * base_acc + user_u;
    m->vel += a * dt;
    m->off += m->vel * dt;
}

/* render one frame with motion blur integrated over the exposure */
static void render_frame(float* img, int w, int h,
                         float pos, float yaw, float vel, float yvel) {
    for (int y = 0; y < h; y++) {
        float vgain = 0.9f + 0.1f * sinf((float)y / h * PI);
        for (int x = 0; x < w; x++) {
            float acc = 0.0f;
            for (int j = 0; j < SUB; j++) {
                float tau = ((float)j + 0.5f) / SUB * T_EXP;
                float p = pos + vel * tau;
                float yw = yaw + yvel * tau;
                float wx = p + ((float)x - w * 0.5f + yw * FOCAL_PX) / PX_PER_UNIT;
                acc += scene(wx);
            }
            img[y * w + x] = acc / SUB * vgain;
        }
    }
}

/* sharpness: mean squared horizontal gradient of the center row */
static float sharpness(const float* img, int w, int h) {
    const float* row = img + (h / 2) * w;
    double s = 0.0;
    for (int x = 0; x < w - 1; x++) {
        float d = row[x + 1] - row[x];
        s += d * d;
    }
    return (float)(s / (w - 1));
}

static int write_pgm(const char* path, const float* img, int w, int h) {
    FILE* f = fopen(path, "wb");
    if (!f) { printf("warn: cannot write %s\n", path); return -1; }
    fprintf(f, "P5\n%d %d\n255\n", w, h);
    for (int i = 0; i < w * h; i++)
        fputc((int)(clampf(img[i], 0.0f, 1.0f) * 255.0f + 0.5f), f);
    fclose(f);
    return 0;
}

static void band_blit(float* dst, int dw, int dy, const float* src, int w, int h) {
    for (int y = 0; y < h; y++)
        memcpy(dst + (dy + y) * dw, src + y * w, (size_t)w * sizeof(float));
}

/* ---- violent-move capture -------------------------------------------- */

typedef struct {
    float pos, yaw, vel, yvel;
    float blur_px, t;
    float peak_stress;
} Capture;

static Capture run_capture(const Tune* tu) {
    Cam c; memset(&c, 0, sizeof(c));
    Capture best; memset(&best, 0, sizeof(best));
    for (int i = 0; i < (int)(1.2f / DT); i++) {
        cam_step(&c, tu, DT);
        float blur = (fabsf(c.vel) * PX_PER_UNIT + fabsf(c.yvel) * FOCAL_PX) * T_EXP;
        if (c.stress > best.peak_stress) best.peak_stress = c.stress;
        if (blur > best.blur_px) {
            best.blur_px = blur;
            best.pos = c.pos; best.yaw = c.yaw;
            best.vel = c.vel; best.yvel = c.yvel;
            best.t = i * DT;
        }
    }
    return best;
}

static float run_overshoot(const Tune* tu) {
    Cam c; memset(&c, 0, sizeof(c));
    float over = 0.0f;
    for (int i = 0; i < (int)(3.0f / DT); i++) {
        cam_step(&c, tu, DT);
        if (c.pos - TGT_POS > over) over = c.pos - TGT_POS;
    }
    return over;
}

/* ---- sudden-stop bounce ---------------------------------------------- */

enum { MODE_RAW, MODE_AUTO, MODE_MANUAL };

typedef struct {
    float t_stop, t_settle;
    float off_min, off_max;  /* screen px excursions after the stop */
    int crossings;           /* sign changes of the offset after the stop */
    int settled;
    /* capture of the worst bounce instant */
    float cap_pos, cap_yaw, cap_off;
} BounceResult;

static BounceResult run_bounce(const Tune* tu, int mode, float* trace, int max_n) {
    Cam c; memset(&c, 0, sizeof(c));
    Mount m; memset(&m, 0, sizeof(m));
    BounceResult r; memset(&r, 0, sizeof(r));
    r.t_stop = -1.0f;

    Tune tt = *tu;
    if (mode != MODE_AUTO) tt.brake_on = 0;

    float u = 0.0f;
    float prev_off = 0.0f, last_over = -1.0f;
    int ti = 0, have_prev = 0, stopped = 0;

    for (int i = 0; i < (int)(4.0f / DT); i++) {
        float t = i * DT;
        cam_step(&c, &tt, DT);

        /* auto mode: mount damping kicks in near the endpoint */
        if (mode == MODE_AUTO)
            m.extra_c = fabsf(TGT_POS - c.pos) < tu->brake_radius ? AUTO_MOUNT_C : 0.0f;

        /* sudden stop: base hits the end-stop and halts within a frame */
        if (!stopped && fabsf(TGT_POS - c.pos) < STOP_REM && fabsf(c.vel) < STOP_VEL) {
            c.acc = -c.vel / DT;       /* decel spike feeds the mount below */
            c.vel = 0.0f;
            c.yvel = 0.0f;
            stopped = 1;
            r.t_stop = t;
        }

        /* manual mode: after a reaction delay the user braces, applying
           force-limited viscous damping + stiffness against the bounce */
        if (mode == MODE_MANUAL && r.t_stop >= 0.0f && t >= r.t_stop + USER_DELAY)
            u = clampf(-USER_GV * m.vel - USER_G * m.off, -USER_UMAX, USER_UMAX);

        mount_step(&m, c.acc, u, DT);

        if (r.t_stop >= 0.0f) {
            float off_px = m.off * PX_PER_UNIT;
            if (trace && ti < max_n) trace[ti++] = off_px;
            if (off_px > r.off_max) r.off_max = off_px;
            if (off_px < r.off_min) r.off_min = off_px;
            if (fabsf(off_px) >= fabsf(r.cap_off * PX_PER_UNIT)) {
                r.cap_off = m.off;
                r.cap_pos = c.pos; r.cap_yaw = c.yaw;
            }
            if (have_prev && fabsf(off_px) > 0.1f &&
                ((off_px > 0.0f) != (prev_off > 0.0f)))
                r.crossings++;
            prev_off = off_px;
            have_prev = 1;
            if (fabsf(off_px) >= SETTLE_PX) last_over = t;
        }
    }
    if (r.t_stop >= 0.0f) {
        r.t_settle = last_over < 0.0f ? 0.0f : last_over - r.t_stop;
        r.settled = 1;
    }
    return r;
}

/* plot two bounce traces (auto bright, manual gray) into an image */
static void plot_traces(float* img, int w, int h,
                        const float* a, const float* b, int n, float ylim) {
    memset(img, 0, (size_t)w * h * sizeof(float));
    for (int x = 0; x < w; x++) {
        img[(h / 2) * w + x] = 0.25f;                       /* zero axis */
        int th = (int)(SETTLE_PX / ylim * (h / 2 - 2));
        img[(h / 2 - th) * w + x] = 0.12f;                  /* settle band */
        img[(h / 2 + th) * w + x] = 0.12f;
    }
    for (int x = 0; x < w; x++) {
        int i = x * n / w;
        if (i >= n) break;
        float ya = clampf(a[i] / ylim, -1.0f, 1.0f);
        float yb = clampf(b[i] / ylim, -1.0f, 1.0f);
        int ra = h / 2 - (int)(ya * (h / 2 - 2));
        int rb = h / 2 - (int)(yb * (h / 2 - 2));
        img[ra * w + x] = 1.0f;
        if (img[rb * w + x] < 0.6f) img[rb * w + x] = 0.6f;
    }
}

/* ---- self test -------------------------------------------------------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    Tune base = {
        "mid", 1.0f,
        22.0f, 3.0f,          /* kp, kd */
        55.0f, 30.0f,         /* a_max, w_max */
        1.6f, 0.5f,           /* brake radii */
        6.0f, 1.5f,           /* brake c0, stress gain */
        1
    };
    Tune low = base;  low.name = "low";  low.sens = 0.55f;
    Tune mid = base;
    Tune high = base; high.name = "high"; high.sens = 1.7f;

    /* joint stress is active during a violent move, brake scales with it */
    Capture cap_mid = run_capture(&mid);
    check("stress_joint_active", cap_mid.peak_stress > 1.2f);
    check("brake_scales_with_stress",
          brake_gain(&mid, 2.0f) > brake_gain(&mid, 0.5f));

    /* endpoint braking damps overshoot */
    Tune free_t = mid; free_t.brake_on = 0;
    float over_on = run_overshoot(&mid);
    float over_free = run_overshoot(&free_t);
    check("brake_reduces_overshoot",
          over_on < 0.25f * over_free && over_on < 0.5f);
    printf("overshoot: brake=%.3f free=%.3f\n", over_on, over_free);

    /* three sensitivities: capture the blurriest frame of each, compare */
    Tune tunes[3] = { low, mid, high };
    float frames[3][IMG_H * IMG_W];
    float sharp[3], blur[3];
    for (int k = 0; k < 3; k++) {
        Capture cap = run_capture(&tunes[k]);
        render_frame(frames[k], IMG_W, IMG_H, cap.pos, cap.yaw, cap.vel, cap.yvel);
        sharp[k] = sharpness(frames[k], IMG_W, IMG_H);
        blur[k] = cap.blur_px;
        char path[64];
        snprintf(path, sizeof(path), "out/frame_%s.pgm", tunes[k].name);
        write_pgm(path, frames[k], IMG_W, IMG_H);
    }
    static float ref[IMG_H * IMG_W];
    render_frame(ref, IMG_W, IMG_H, 4.0f, 0.5f, 0.0f, 0.0f);
    printf("blur_px: low=%.2f mid=%.2f high=%.2f\n", blur[0], blur[1], blur[2]);
    printf("sharpness: low=%.5f mid=%.5f high=%.5f static=%.5f\n",
           sharp[0], sharp[1], sharp[2], sharpness(ref, IMG_W, IMG_H));
    check("sharpness_low_gt_mid", sharp[0] > sharp[1]);
    check("sharpness_mid_gt_high", sharp[1] > sharp[2]);

    /* stacked overlay: 3 frames + 3 edge maps + blend of the three */
    {
        int eh = 40, sep = 4;
        int H = 3 * (IMG_H + eh + sep) + IMG_H;
        float* cmp = calloc((size_t)IMG_W * H, sizeof(float));
        int dy = 0;
        for (int k = 0; k < 3; k++) {
            band_blit(cmp, IMG_W, dy, frames[k], IMG_W, IMG_H);
            dy += IMG_H;
            for (int y = 0; y < eh; y++)
                for (int x = 0; x < IMG_W - 1; x++) {
                    float g = fabsf(frames[k][IMG_H * IMG_W / 2 + x + 1]
                                  - frames[k][IMG_H * IMG_W / 2 + x]);
                    cmp[(dy + y) * IMG_W + x] = clampf(g * 4.0f, 0.0f, 1.0f);
                }
            dy += eh + sep;
        }
        for (int i = 0; i < IMG_H * IMG_W; i++)
            cmp[dy * IMG_W + i] = (frames[0][i] + frames[1][i] + frames[2][i]) / 3.0f;
        write_pgm("out/sharpness_compare.pgm", cmp, IMG_W, H);
        free(cmp);
    }

    /* sudden-stop micro bounce: raw shows it, auto and manual counteract it */
    float trace_a[90], trace_m[90];
    BounceResult raw = run_bounce(&mid, MODE_RAW, NULL, 0);
    BounceResult au = run_bounce(&mid, MODE_AUTO, trace_a, 90);
    BounceResult ma = run_bounce(&mid, MODE_MANUAL, trace_m, 90);
    printf("bounce: raw_max=%.2fpx raw_min=%.2fpx crossings=%d\n",
           raw.off_max, raw.off_min, raw.crossings);
    printf("settle: raw=%.0fms auto=%.0fms manual=%.0fms gap=%.0fms half_frame=%.1fms\n",
           raw.t_settle * 1000.0f,
           au.t_settle * 1000.0f, ma.t_settle * 1000.0f,
           fabsf(au.t_settle - ma.t_settle) * 1000.0f, HALF_FRAME * 1000.0f);
    check("bounce_back_occurs", raw.crossings >= 2 && raw.off_max > SETTLE_PX);
    check("auto_settles", au.settled);
    check("manual_settles", ma.settled);
    check("auto_faster_than_manual", au.t_settle < ma.t_settle);
    /* spec: a mode gap shorter than half a frame counts as failure */
    check("settle_gap_ge_half_frame",
          fabsf(au.t_settle - ma.t_settle) >= HALF_FRAME);

    /* bounce capture: worst bounce instant vs settled frame */
    {
        static float peak[IMG_H * IMG_W], calm[IMG_H * IMG_W];
        render_frame(peak, IMG_W, IMG_H,
                     raw.cap_pos + raw.cap_off, raw.cap_yaw, 0.0f, 0.0f);
        render_frame(calm, IMG_W, IMG_H, TGT_POS, TGT_YAW, 0.0f, 0.0f);
        static float both[2 * IMG_H * IMG_W];
        band_blit(both, IMG_W, 0, peak, IMG_W, IMG_H);
        band_blit(both, IMG_W, IMG_H, calm, IMG_W, IMG_H);
        write_pgm("out/bounce_capture.pgm", both, IMG_W, 2 * IMG_H);
    }
    {
        static float plot[180 * IMG_W];
        plot_traces(plot, IMG_W, 180, trace_a, trace_m, 90, 3.0f);
        write_pgm("out/bounce_trace.pgm", plot, IMG_W, 180);
        FILE* f = fopen("out/bounce_trace.csv", "w");
        if (f) {
            fprintf(f, "t_ms,auto_px,manual_px\n");
            for (int i = 0; i < 90; i++)
                fprintf(f, "%.1f,%.3f,%.3f\n", i * DT * 1000.0f,
                        trace_a[i], trace_m[i]);
            fclose(f);
        }
    }

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
