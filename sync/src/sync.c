
/* sync: resampling-ratio estimation from paired event timestamps.
 *
 * Model: y_i = a + b*x_i + e_i, where x is the animation-event clock and
 * y is the audio clock.  The slope b is the effective resampling ratio;
 * drift between the two clocks accumulates as (b - 1) * elapsed time.
 *
 * Pipeline: Theil-Sen robust initial fit -> scale-based outlier rejection
 * with refit -> OLS on the inliers -> Student-t confidence interval ->
 * residual stationarity checks.  Estimates that fail the guards are
 * rejected with a reason code rather than returned as if trustworthy.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define SYNC_MAX_PTS 1024
#define SYNC_MIN_PTS 8            /* fewer pairs than this: refuse */
#define SYNC_MAX_OUTLIER_FRAC 0.25
#define SYNC_REJECT_K 3.5         /* residual flag threshold, in robust sigmas */
#define SYNC_MAX_REFITS 3
#define SYNC_STAT_Z 4.0           /* stationarity guard, in standard errors */
#define SYNC_TINY 1e-12
#define SYNC_PI 3.14159265358979323846

typedef struct { double x, y; } Pt;

typedef struct {
    int    rejected;             /* 1 => do not use ratio/ci_* */
    const char* reject_reason;   /* static string, NULL when accepted */
    double ratio;                /* estimated slope b */
    double offset;               /* estimated intercept a */
    double ci_lo, ci_hi;         /* two-sided 95% interval on b */
    double ci_halfwidth;
    double se;                   /* OLS standard error of the slope */
    double residual_sigma;       /* inlier residual standard deviation */
    double robust_scale;         /* robust residual scale used for flagging */
    double hard_bound;           /* worst-case |bhat - b| under |e| <= E,
                                    -1 when no noise bound was supplied */
    int    n, inlier_count, outlier_count;
    int    stationarity_ok;
    int    refits;
} SyncResult;

/* ---------- small stats helpers ---------- */

static int cmp_double(const void* pa, const void* pb) {
    double a = *(const double*)pa, b = *(const double*)pb;
    return (a > b) - (a < b);
}

static int cmp_pt_x(const void* pa, const void* pb) {
    double a = ((const Pt*)pa)->x, b = ((const Pt*)pb)->x;
    return (a > b) - (a < b);
}

static double median_sorted(const double* v, int n) {
    if (n <= 0) return 0.0;
    int h = n / 2;
    return (n & 1) ? v[h] : 0.5 * (v[h - 1] + v[h]);
}

/* Ordinary least squares.  Returns 0 on success, -1 on degenerate x. */
static int ols(const Pt* p, const char* mask, int n,
               double* b, double* a, double* sigma, double* se) {
    double sx = 0.0, sy = 0.0;
    int m = 0;
    for (int i = 0; i < n; i++) {
        if (mask[i]) { sx += p[i].x; sy += p[i].y; m++; }
    }
    if (m < 3) return -1;
    double mx = sx / m, my = sy / m;
    double sxx = 0.0, sxy = 0.0, sse = 0.0;
    for (int i = 0; i < n; i++) {
        if (!mask[i]) continue;
        double dx = p[i].x - mx;
        sxx += dx * dx;
        sxy += dx * (p[i].y - my);
    }
    if (sxx <= SYNC_TINY) return -1;
    double bb = sxy / sxx;
    double aa = my - bb * mx;
    for (int i = 0; i < n; i++) {
        if (!mask[i]) continue;
        double r = p[i].y - (aa + bb * p[i].x);
        sse += r * r;
    }
    *b = bb;
    *a = aa;
    *sigma = sqrt(sse / (m - 2));
    *se = *sigma / sqrt(sxx);
    return 0;
}

/* Two-sided 95% Student-t critical values, df 1..30, normal beyond. */
static double t975(int df) {
    static const double t[30] = {
        12.706, 4.303, 3.182, 2.776, 2.571, 2.447, 2.365, 2.306,
        2.262,  2.228, 2.201, 2.179, 2.160, 2.145, 2.131, 2.120,
        2.110,  2.101, 2.093, 2.086, 2.080, 2.074, 2.069, 2.064,
        2.060,  2.056, 2.052, 2.048, 2.045, 2.042
    };
    if (df < 1) return 0.0;
    if (df <= 30) return t[df - 1];
    return 1.96;
}

/* Robust residual scale that stays inlier-driven while outliers are the
 * majority: take the median of the lower half of |r| (the 25th percentile
 * of |r|) and undo the half-normal CDF at 0.25 (= 0.3186 sigma). */
static double robust_scale(const double* r, int n, double* scratch) {
    for (int i = 0; i < n; i++) scratch[i] = fabs(r[i]);
    qsort(scratch, (size_t)n, sizeof(double), cmp_double);
    int half = n / 2;
    if (half < 1) return 0.0;
    double q25 = median_sorted(scratch, half);
    return q25 / 0.3186;
}

/* Theil-Sen: median of all pairwise slopes, median intercept afterwards.
 * Tolerates up to ~29% gross contamination before its own breakdown. */
static int theil_sen(const Pt* p, int n, double* b, double* a,
                     double* scratch) {
    size_t k = 0, cap = (size_t)n * (n - 1) / 2;
    if (cap == 0) return -1;
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            double dx = p[j].x - p[i].x;
            if (fabs(dx) > SYNC_TINY) scratch[k++] = (p[j].y - p[i].y) / dx;
        }
    }
    if (k == 0) return -1;
    qsort(scratch, k, sizeof(double), cmp_double);
    double bb = median_sorted(scratch, (int)k);
    for (int i = 0; i < n; i++) scratch[i] = p[i].y - bb * p[i].x;
    qsort(scratch, (size_t)n, sizeof(double), cmp_double);
    *b = bb;
    *a = median_sorted(scratch, n);
    return 0;
}

/* Convert an F(d1,d2) statistic to an approximate normal z, so the same
 * SYNC_STAT_Z threshold applies to every guard.  Scale F onto a chi-square
 * variable and apply the Wilson-Hilferty cube-root normalizing transform.
 * Accurate enough at the tail we care about (d1 = 2, d2 >= ~30). */
static double fstat_to_z(double f, int d1, int d2) {
    double chi = d1 * f * (d2 > 100 ? 1.0 : (double)d2 / (d2 - 2.0));
    double wh = pow(chi / d1, 1.0 / 3.0);
    double mean = 1.0 - 2.0 / (9.0 * d1);
    double sd = sqrt(2.0 / (9.0 * d1));
    return (wh - mean) / sd;
}

/* Residuals must look stationary: a drifting ratio shows up as a slope or
 * ramp in the residuals even when every individual residual is small. */
static int residuals_stationary(const Pt* p, const char* mask, int n,
                                double b, double a, double sigma) {
    int m = 0;
    for (int i = 0; i < n; i++) if (mask[i]) m++;
    if (m < 4) return 0;

    double* r = (double*)malloc((size_t)m * sizeof(double));
    if (!r) return 0;
    int k = 0;
    for (int i = 0; i < n; i++)
        if (mask[i]) r[k++] = p[i].y - (a + b * p[i].x);  /* p sorted by x */

    int half = m / 2;
    int ok = 1;
    if (sigma <= SYNC_TINY) { free(r); return 1; }  /* exact fit: stationary */

    /* (1) half-window mean drift, Student-t sized */
    double mean1 = 0.0, mean2 = 0.0;
    for (int i = 0; i < half; i++) mean1 += r[i];
    for (int i = half; i < m; i++) mean2 += r[i];
    mean1 /= half; mean2 /= (m - half);
    double drift_z = fabs(mean1 - mean2) /
                     (sigma * sqrt(1.0 / half + 1.0 / (m - half)));
    if (drift_z > SYNC_STAT_Z) ok = 0;

    /* (2) Chow break test: does one line describe the whole span, or do the
     * two halves want different lines?  A slope change at mid-stream leaves
     * a V-shaped residual whose half-window slopes have opposite signs; a
     * curvature statistic loses most of that signal, but the F-test for the
     * two extra half-specific parameters does not. */
    double sse_full = 0.0;
    for (int i = 0; i < m; i++) sse_full += r[i] * r[i];

    double sse_halves = 0.0;
    int degenerate = 0;
    int idx;
    for (int pass = 0; pass < 2; pass++) {
        int k0 = pass == 0 ? 0 : half;
        int k1 = pass == 0 ? half : m;
        int hn = k1 - k0;
        if (hn < 2) { degenerate = 1; break; }
        double hmx = 0.0, ymean = 0.0;
        idx = 0;
        for (int i = 0; i < n; i++) {
            if (!mask[i]) continue;
            if (idx >= k0 && idx < k1) { hmx += p[i].x; ymean += p[i].y; }
            idx++;
        }
        hmx /= hn; ymean /= hn;
        double hsxx = 0.0, hsxy = 0.0;
        idx = 0;
        for (int i = 0; i < n; i++) {
            if (!mask[i]) continue;
            if (idx >= k0 && idx < k1) {
                double dx = p[i].x - hmx;
                hsxx += dx * dx;
                hsxy += dx * (p[i].y - ymean);
            }
            idx++;
        }
        if (hsxx <= SYNC_TINY) { degenerate = 1; break; }
        double hb = hsxy / hsxx;
        double ha = ymean - hb * hmx;
        idx = 0;
        for (int i = 0; i < n; i++) {
            if (!mask[i]) continue;
            if (idx >= k0 && idx < k1) {
                double er = p[i].y - (ha + hb * p[i].x);
                sse_halves += er * er;
            }
            idx++;
        }
    }
    if (!degenerate && sse_full > sse_halves + SYNC_TINY && m > 4) {
        double fstat = ((sse_full - sse_halves) / 2.0) /
                       (sse_halves / (m - 4));
        if (fstat_to_z(fstat, 2, m - 4) > SYNC_STAT_Z) ok = 0;
    }

    /* (3) cumulative residual walk must stay within random-walk size */
    double cum = 0.0, maxcum = 0.0;
    for (int i = 0; i < m; i++) {
        cum += r[i];
        if (fabs(cum) > maxcum) maxcum = fabs(cum);
    }
    if (sigma > 0.0 && maxcum > SYNC_STAT_Z * sigma * sqrt((double)m)) ok = 0;

    free(r);
    return ok;
}

/* Deterministic worst-case slope error for OLS under bounded inlier noise
 * |e_i| <= E:  bhat - b = sum w_i e_i with sum w_i = 0, hence the bound
 * E * sum |w_i| = E * sum|x_i - xbar| / Sxx. */
static double slope_error_bound(const Pt* p, const char* mask, int n, double E) {
    if (E <= 0.0) return -1.0;
    double sx = 0.0;
    int m = 0;
    for (int i = 0; i < n; i++) if (mask[i]) { sx += p[i].x; m++; }
    if (m < 3) return -1.0;
    double mx = sx / m, sxx = 0.0, sa = 0.0;
    for (int i = 0; i < n; i++) {
        if (!mask[i]) continue;
        double dx = p[i].x - mx;
        sxx += dx * dx;
        sa += fabs(dx);
    }
    if (sxx <= SYNC_TINY) return -1.0;
    return E * sa / sxx;
}

static void sync_reject(SyncResult* out, const char* reason) {
    out->rejected = 1;
    out->reject_reason = reason;
}

/* Main entry.  noise_bound E is optional (pass <= 0 to skip the hard bound).
 * x/y need not be sorted; they are copied internally. */
static void sync_estimate(const double* x, const double* y, int n,
                          double noise_bound, SyncResult* out) {
    memset(out, 0, sizeof(*out));
    out->hard_bound = -1.0;
    out->n = n;

    if (n < SYNC_MIN_PTS) { sync_reject(out, "too_few_points"); return; }

    Pt* p = (Pt*)malloc((size_t)n * sizeof(Pt));
    char* mask = (char*)malloc((size_t)n);
    double* resid = (double*)malloc((size_t)n * sizeof(double));
    double* scratch = (double*)malloc((size_t)n * (n - 1) / 2 * sizeof(double));
    if (!p || !mask || !resid || !scratch) {
        sync_reject(out, "alloc_failed");
        goto done;
    }
    for (int i = 0; i < n; i++) { p[i].x = x[i]; p[i].y = y[i]; }
    qsort(p, (size_t)n, sizeof(Pt), cmp_pt_x);
    if (p[n - 1].x - p[0].x <= SYNC_TINY) {
        sync_reject(out, "no_x_span");
        goto done;
    }

    double b, a;
    if (theil_sen(p, n, &b, &a, scratch) != 0) {
        sync_reject(out, "degenerate_fit");
        goto done;
    }

    /* Iterative flag -> refit.  The lower-half scale cannot be inflated by
     * a majority of gross outliers, so a contaminated majority gets flagged
     * and trips the fraction guard instead of being silently fitted. */
    for (int iter = 0; iter <= SYNC_MAX_REFITS; iter++) {
        for (int i = 0; i < n; i++)
            resid[i] = p[i].y - (a + b * p[i].x);
        double s = robust_scale(resid, n, scratch);
        out->robust_scale = s;

        int outliers = 0;
        if (s > SYNC_TINY) {
            for (int i = 0; i < n; i++)
                mask[i] = (fabs(resid[i]) <= SYNC_REJECT_K * s);
        } else {
            /* zero robust scale: the robust fit already passes through every
             * point, so nothing there is gross enough to flag */
            for (int i = 0; i < n; i++) mask[i] = 1;
        }
        for (int i = 0; i < n; i++) if (!mask[i]) outliers++;

        if ((double)outliers / n > SYNC_MAX_OUTLIER_FRAC) {
            sync_reject(out, "outlier_fraction_too_high");
            out->ratio = b; out->offset = a;
            out->outlier_count = outliers;
            out->inlier_count = n - outliers;
            goto done;
        }

        double nb, na, nsigma, nse;
        if (ols(p, mask, n, &nb, &na, &nsigma, &nse) != 0) {
            sync_reject(out, "degenerate_fit");
            goto done;
        }
        int stable = (fabs(nb - b) <= SYNC_TINY * (1.0 + fabs(b)));
        b = nb; a = na;
        out->refits = iter + 1;
        if (stable) break;
    }

    int inliers = 0;
    for (int i = 0; i < n; i++) if (mask[i]) inliers++;
    out->inlier_count = inliers;
    out->outlier_count = n - inliers;

    double sigma, se;
    if (ols(p, mask, n, &b, &a, &sigma, &se) != 0 || se <= 0.0) {
        sync_reject(out, "degenerate_fit");
        goto done;
    }
    out->ratio = b;
    out->offset = a;
    out->residual_sigma = sigma;
    out->se = se;

    out->stationarity_ok = residuals_stationary(p, mask, n, b, a, sigma);
    if (!out->stationarity_ok) {
        sync_reject(out, "residuals_nonstationary");
        goto done;
    }

    double tc = t975(inliers - 2);
    out->ci_halfwidth = tc * se;
    out->ci_lo = b - out->ci_halfwidth;
    out->ci_hi = b + out->ci_halfwidth;
    out->hard_bound = slope_error_bound(p, mask, n, noise_bound);

done:
    free(p); free(mask); free(resid); free(scratch);
}

/* ---------- deterministic PRNG for self-tests ---------- */

static unsigned long long rng_state;

static double rand_uniform(void) {
    unsigned long long x = rng_state;
    x ^= x << 13; x ^= x >> 7; x ^= x << 17;
    rng_state = x;
    return (double)(x >> 11) / (double)(1ULL << 53);  /* [0,1) */
}

static double rand_gaussian(void) {
    double u1 = rand_uniform(), u2 = rand_uniform();
    if (u1 < 1e-12) u1 = 1e-12;
    return sqrt(-2.0 * log(u1)) * cos(2.0 * SYNC_PI * u2);
}

/* ---------- self tests ---------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

#define MC_TRIALS 200
#define N_PTS 60
#define DT 0.0166667              /* ~60 Hz event cadence, seconds */
#define B_TRUE 1.0005             /* 500 ppm resampling ratio */
#define A_TRUE (-0.002)
#define SIGMA 0.0005              /* 0.5 ms gaussian timestamp jitter */

static void make_pair(double* x, double* y, int i, double extra) {
    x[i] = i * DT;
    y[i] = A_TRUE + B_TRUE * x[i] + SIGMA * rand_gaussian() + extra;
}

int main(void) {
    rng_state = 0x123456789abcdef0ULL;
    double x[N_PTS], y[N_PTS];
    SyncResult r;

    /* guard: too few points */
    for (int i = 0; i < 4; i++) { x[i] = i * DT; y[i] = x[i]; }
    sync_estimate(x, y, 4, 0.0, &r);
    check("too_few_points", r.rejected && strcmp(r.reject_reason, "too_few_points") == 0);

    /* guard: no x spread */
    for (int i = 0; i < 20; i++) { x[i] = 1.0; y[i] = 1.0 + 0.001 * i; }
    sync_estimate(x, y, 20, 0.0, &r);
    check("no_x_span", r.rejected && strcmp(r.reject_reason, "no_x_span") == 0);

    /* exact line: essentially machine-precision recovery */
    for (int i = 0; i < N_PTS; i++) {
        x[i] = i * DT;
        y[i] = A_TRUE + B_TRUE * x[i];
    }
    sync_estimate(x, y, N_PTS, 0.0, &r);
    check("exact_line_accepted", !r.rejected);
    check("exact_line_slope", !r.rejected && fabs(r.ratio - B_TRUE) < 1e-12);

    /* Monte Carlo: clean gaussian noise, CI coverage and stationary residuals */
    int covered = 0, accepted = 0, stationary = 0;
    double sum_err = 0.0, sum_hw = 0.0, max_err = 0.0;
    for (int t = 0; t < MC_TRIALS; t++) {
        for (int i = 0; i < N_PTS; i++) make_pair(x, y, i, 0.0);
        sync_estimate(x, y, N_PTS, 0.0, &r);
        if (!r.rejected) {
            accepted++;
            if (r.ci_lo <= B_TRUE && B_TRUE <= r.ci_hi) covered++;
            if (r.stationarity_ok) stationary++;
            double err = fabs(r.ratio - B_TRUE);
            sum_err += err;
            sum_hw += r.ci_halfwidth;
            if (err > max_err) max_err = err;
        }
    }
    check("clean_all_accepted", accepted == MC_TRIALS);
    check("clean_coverage_90pct", covered >= (int)(0.90 * MC_TRIALS));
    check("clean_all_stationary", stationary == MC_TRIALS);
    check("ci_hw_dominates_err", sum_hw > sum_err);
    printf("clean_coverage=%d/%d max_err_ppm=%.2f mean_hw_ppm=%.2f\n",
           covered, MC_TRIALS, max_err / B_TRUE * 1e6,
           sum_hw / MC_TRIALS / B_TRUE * 1e6);

    /* hard deterministic bound under uniform bounded noise */
    const double E = 0.001;       /* +/- 1 ms */
    int bound_ok = 0;
    double worst_err = 0.0, worst_bound = 0.0;
    for (int t = 0; t < 100; t++) {
        for (int i = 0; i < N_PTS; i++) {
            x[i] = i * DT;
            y[i] = A_TRUE + B_TRUE * x[i] + (2.0 * rand_uniform() - 1.0) * E;
        }
        sync_estimate(x, y, N_PTS, E, &r);
        if (!r.rejected && r.hard_bound > 0.0) {
            double err = fabs(r.ratio - B_TRUE);
            if (err > worst_err) worst_err = err;
            if (r.hard_bound > worst_bound) worst_bound = r.hard_bound;
            if (err <= r.hard_bound * (1.0 + 1e-12)) bound_ok++;
        }
    }
    check("hard_bound_holds", bound_ok == 100);
    printf("hard_bound_ppm worst_err=%.2f worst_bound=%.2f\n",
           worst_err / B_TRUE * 1e6, worst_bound / B_TRUE * 1e6);

    /* stationarity false-positive rate over 500 clean trials */
    int fp = 0;
    for (int t = 0; t < 500; t++) {
        for (int i = 0; i < N_PTS; i++) make_pair(x, y, i, 0.0);
        sync_estimate(x, y, N_PTS, 0.0, &r);
        if (r.rejected) fp++;
    }
    check("stationarity_false_positives", fp <= 5);  /* ~1% budget, z>=4 */

    /* gradual drift: ratio ramps linearly over the whole span, no jump */
    for (int i = 0; i < N_PTS; i++) {
        double bi = B_TRUE + 0.006 * (double)i / (N_PTS - 1);
        x[i] = i * DT;
        y[i] = A_TRUE + bi * x[i] + SIGMA * rand_gaussian();
    }
    sync_estimate(x, y, N_PTS, 0.0, &r);
    check("gradual_drift_rejected", r.rejected &&
          strcmp(r.reject_reason, "residuals_nonstationary") == 0);

    /* 10% gross outliers: rejected points removed, ratio still covered */
    for (int i = 0; i < N_PTS; i++) {
        double gross = (i % 10 == 3) ? ((i / 10 % 2) ? 0.020 : -0.020) : 0.0;
        make_pair(x, y, i, gross);
    }
    sync_estimate(x, y, N_PTS, 0.0, &r);
    check("gross10_accepted", !r.rejected);
    check("gross10_outliers_flagged", !r.rejected && r.outlier_count >= 4 && r.outlier_count <= 8);
    check("gross10_covered", !r.rejected && r.ci_lo <= B_TRUE && B_TRUE <= r.ci_hi);
    check("gross10_stationary", !r.rejected && r.stationarity_ok);
    if (!r.rejected)
        printf("gross10_err_ppm=%.2f outliers=%d\n",
               fabs(r.ratio - B_TRUE) / B_TRUE * 1e6, r.outlier_count);

    /* 60% contamination: refuse rather than fit through the garbage */
    for (int i = 0; i < N_PTS; i++) {
        double gross = (i % 10 < 6) ? 0.020 : 0.0;
        make_pair(x, y, i, gross);
    }
    sync_estimate(x, y, N_PTS, 0.0, &r);
    check("gross60_rejected", r.rejected &&
          strcmp(r.reject_reason, "outlier_fraction_too_high") == 0);

    /* ratio that changes mid-stream: residuals ramp, estimate must be refused */
    for (int i = 0; i < N_PTS; i++) {
        double bi = (i < N_PTS / 2) ? B_TRUE : 1.004;
        x[i] = i * DT;
        y[i] = A_TRUE + bi * x[i] + SIGMA * rand_gaussian();
    }
    sync_estimate(x, y, N_PTS, 0.0, &r);
    check("drifting_ratio_rejected", r.rejected &&
          strcmp(r.reject_reason, "residuals_nonstationary") == 0);

    /* CI sanity on a clean fit */
    for (int i = 0; i < N_PTS; i++) make_pair(x, y, i, 0.0);
    sync_estimate(x, y, N_PTS, 0.0, &r);
    check("ci_ordering", !r.rejected &&
          r.ci_lo < r.ratio && r.ratio < r.ci_hi && r.ci_halfwidth > 0.0);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
