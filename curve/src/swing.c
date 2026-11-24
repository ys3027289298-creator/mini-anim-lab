
/* swing: record arm-swing angle sequences at slow-walk / fast-walk / run
   cadences, smooth the extrema tangents, and count angular-velocity
   jumps against the curvature bound of the smoothed curve.
   usage: ./swing [out_dir]   (csvs go to out_dir, default ../artifacts) */
#define CURVE_NO_MAIN
#include "curve.c"

#define FPS 60.0f
#define DUR 4.0f
#define AMP 0.5f          /* swing amplitude, rad */
#define DUTY_DOWN 0.575f  /* peak->trough fraction of the period */
#define SPIKE_T 2.5f      /* one authored glitch: a peak overshoots here */
#define SPIKE_V 0.85f
#define NFRAMES ((int)(DUR * FPS) + 1)

typedef struct { const char* name; float freq; } Cadence;
static const Cadence CAD[] = {
    { "slow_walk", 0.9f },
    { "fast_walk", 1.5f },
    { "run",       2.4f },
};
#define NCAD (int)(sizeof(CAD) / sizeof(CAD[0]))

/* curvature budget for a smoothed swing: the steepest half-phase is the
   short (1-DUTY_DOWN)*T leg, whose zero-tangent hermite bound is
   12*A/dt^2; 10% headroom. */
static float kappa_max_for(float freq) {
    float T = 1.0f / freq;
    float dt_min = (1.0f - DUTY_DOWN) * T;
    return 1.1f * 12.0f * AMP / (dt_min * dt_min);
}

/* kinked swing curve: keys at swing extrema, tan_in = 0 but tan_out =
   full outgoing secant -- the tangent discontinuity that makes the arm
   hitch one frame at the top of the swing. */
static void build_swing(Curve* c, float freq) {
    float T = 1.0f / freq;
    float ts[MAX_KEYS], vs[MAX_KEYS];
    int n = 0;
    int k_spike = (int)(SPIKE_T / T + 0.5f);
    for (int k = -1; k * T < DUR + T && n < MAX_KEYS - 2; k++) {
        ts[n] = k * T;
        vs[n] = (k == k_spike) ? SPIKE_V : AMP;
        n++;
        ts[n] = k * T + DUTY_DOWN * T;
        vs[n] = -AMP;
        n++;
    }
    curve_init(c);
    for (int i = 0; i < n; i++) {
        float tout = 0.0f;
        if (i + 1 < n) tout = (vs[i + 1] - vs[i]) / (ts[i + 1] - ts[i]);
        curve_add(c, ts[i], vs[i], 0.0f, tout);
    }
}

static void bake(const Curve* c, float* th) {
    for (int i = 0; i < NFRAMES; i++) th[i] = curve_eval(c, i / FPS);
}

/* jump: per-frame omega change beyond what curvature kappa_max can
   produce in one frame (|d_omega| <= kappa*dt for any curve in budget).
   consecutive violating frames merge into one cluster. */
static int count_jumps(const float* th, float kappa_max, int* clusters) {
    float thresh = kappa_max / FPS;
    int jumps = 0, cl = 0, open = 0;
    for (int i = 1; i < NFRAMES - 1; i++) {
        float o0 = (th[i] - th[i - 1]) * FPS;
        float o1 = (th[i + 1] - th[i]) * FPS;
        if (fabsf(o1 - o0) > thresh) {
            jumps++;
            if (!open) { cl++; open = 1; }
        } else {
            open = 0;
        }
    }
    *clusters = cl;
    return jumps;
}

static const char* bound_str(char* buf, float b) {
    if (isinf(b)) return "inf";
    sprintf(buf, "%.1f", b);
    return buf;
}

int main(int argc, char** argv) {
    const char* out_dir = argc > 1 ? argv[1] : "../artifacts";
    float w_fixed = 1.0f / CAD[0].freq / 6.0f;  /* calibrated on slow walk, then frozen */

    printf("swing analysis: fps=%.0f amp=%.2frad duty=%.3f/%.3f spike@%.1fs\n",
           FPS, AMP, DUTY_DOWN, 1.0f - DUTY_DOWN, SPIKE_T);
    printf("fixed window w=%.4fs (=T_slow/6), scaled window w=T/6 per cadence\n\n",
           w_fixed);
    printf("%-10s %5s %8s | %-22s | %-22s | %-22s\n",
           "cadence", "f", "kappa_max",
           "raw: bound,clusters", "fixed-w: clusters,special", "scaled-w: clusters,special");

    int fast_fixed_clusters = 0, fast_expected = 0;
    for (int ci = 0; ci < NCAD; ci++) {
        float f = CAD[ci].freq;
        float T = 1.0f / f;
        float kmax = kappa_max_for(f);

        Curve raw, fixed, scaled;
        build_swing(&raw, f);
        fixed = raw;
        scaled = raw;
        int sk_fixed = 0, sk_scaled = 0;
        int sp_fixed = curve_smooth_extrema(&fixed, w_fixed, kmax, &sk_fixed);
        int sp_scaled = curve_smooth_extrema(&scaled, T / 6.0f, kmax, &sk_scaled);

        float th_raw[NFRAMES], th_fixed[NFRAMES], th_scaled[NFRAMES];
        bake(&raw, th_raw);
        bake(&fixed, th_fixed);
        bake(&scaled, th_scaled);

        int cl_raw, cl_fixed, cl_scaled;
        int j_raw = count_jumps(th_raw, kmax, &cl_raw);
        int j_fixed = count_jumps(th_fixed, kmax, &cl_fixed);
        int j_scaled = count_jumps(th_scaled, kmax, &cl_scaled);

        char b_raw[16], b_fixed[16], b_scaled[16];
        char col_raw[32], col_fixed[32], col_scaled[32];
        sprintf(col_raw, "%s, %d", bound_str(b_raw, curve_curvature_bound(&raw)), cl_raw);
        sprintf(col_fixed, "%s, %d, %d", bound_str(b_fixed, curve_curvature_bound(&fixed)),
                cl_fixed, sp_fixed);
        sprintf(col_scaled, "%s, %d, %d", bound_str(b_scaled, curve_curvature_bound(&scaled)),
                cl_scaled, sp_scaled);
        printf("%-10s %5.2f %8.1f | %-22s | %-22s | %-22s\n",
               CAD[ci].name, f, kmax, col_raw, col_fixed, col_scaled);
        printf("%-10s %5s %8s | jumps=%-16d | jumps=%-3d skipped=%-9d | jumps=%-3d skipped=%-9d\n",
               "", "", "", j_raw, j_fixed, sk_fixed, j_scaled, sk_scaled);

        if (ci == 1) { fast_fixed_clusters = cl_fixed; fast_expected = sp_scaled; }

        char path[256];
        sprintf(path, "%s/swing_%s.csv", out_dir, CAD[ci].name);
        FILE* fp = fopen(path, "w");
        if (fp) {
            fprintf(fp, "t,th_raw,th_fixed,th_scaled,om_raw,om_fixed,om_scaled\n");
            for (int i = 0; i < NFRAMES; i++) {
                float om_raw = 0, om_fixed = 0, om_scaled = 0;
                if (i < NFRAMES - 1) {
                    om_raw = (th_raw[i + 1] - th_raw[i]) * FPS;
                    om_fixed = (th_fixed[i + 1] - th_fixed[i]) * FPS;
                    om_scaled = (th_scaled[i + 1] - th_scaled[i]) * FPS;
                }
                fprintf(fp, "%.4f,%.6f,%.6f,%.6f,%.4f,%.4f,%.4f\n",
                        i / FPS, th_raw[i], th_fixed[i], th_scaled[i],
                        om_raw, om_fixed, om_scaled);
            }
            fclose(fp);
            printf("  recorded %s (%d frames)\n", path, NFRAMES);
        } else {
            printf("  warn: cannot write %s\n", path);
        }
    }

    printf("\nrelation: a tangent kink makes the curvature bound infinite and shows up as\n");
    printf("one omega-jump cluster per extremum, so raw clusters ~= #extrema. smoothing sets\n");
    printf("tan_in == tan_out; extrema that would exceed kappa_max are clamped and counted\n");
    printf("as special (handled separately). after smoothing clusters == 0 iff the bound is\n");
    printf("<= kappa_max everywhere; any remaining cluster localizes an infeasible extremum.\n");
    printf("\ndiagnostic: fast_walk fixed-window clusters = %d, expected ~= 0 after smoothing\n",
           fast_fixed_clusters);
    printf("window scaled with cadence: %s\n",
           fast_fixed_clusters > fast_expected + 1
               ? "NO -- fixed window overlaps neighboring extrema at higher cadence, they get skipped and keep their kinks"
               : "yes");
    return 0;
}
