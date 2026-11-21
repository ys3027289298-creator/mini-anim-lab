
/* warp: monotone time-warp alignment between two sampled animation tracks
   (e.g. position/rotation/scale curves whose frequencies fight each other).

   - dynamic program over the frame-pair grid
   - transition cost = path-length-weighted (position difference + velocity
     difference) + bend penalty |local slope - 1| + smoothness penalty on the
     slope change between adjacent steps. Length weighting means a long jump
     pays for every frame it skips, so the DP cannot cheapen itself by
     skipping cells; the smoothness term breaks ties on periodic signals.
   - path is strictly monotone; every step's physical-time slope is bounded
   - grid over a cell budget: over-budget calls do NOT run, they return a
     decimation plan (and Sakoe-Chiba band info) for the caller to apply
   - resample track B onto A's frame grid through the resulting mapping
   - compare against hand-authored anchor pairs and report the sync error
     upper bound

   Single-threaded: the DP workspace is file-scope static storage. */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define WARP_MAX_CHAN   6
#define WARP_MAX_N      256
#define WARP_MAX_STEP   8                  /* max frames a single DP step jumps */
#define WARP_BUDGET     4096               /* max DP cells solved in one call */
#define WARP_MAX_ANCHOR 64
#define WARP_INF        1e30f
#define WARP_EPS        1e-6f
#define WARP_PI         3.14159265358979323846f

typedef struct {
    float t[WARP_MAX_N];                   /* physical time, strictly increasing */
    float v[WARP_MAX_CHAN][WARP_MAX_N];
    int wrap_chan[WARP_MAX_CHAN];          /* 1 = circular channel (angles) */
    int n, nchan;
} Track;

typedef struct {
    float w_pos;                           /* weight on value difference */
    float w_vel;                           /* weight on velocity difference */
    float w_bend;                          /* weight on |slope - 1| bending */
    float w_smooth;                        /* weight on |delta slope| between steps */
    float slope_min;                       /* physical-time slope lower bound */
    float slope_max;                       /* physical-time slope upper bound */
} WarpConfig;

typedef enum {
    WARP_OK = 0,
    WARP_ERR_INPUT = -1,
    WARP_ERR_BUDGET = -2,                  /* caller must decimate, see plan */
    WARP_ERR_INFEASIBLE = -3               /* slope bounds cannot be met */
} WarpStatus;

/* reported only on WARP_ERR_BUDGET: how to downsample before retrying. */
typedef struct {
    int dec_a, dec_b;                      /* keep every kth frame */
    int n_a_after, n_b_after;
    int cells_after;                       /* decimated grid size, <= budget */
    int band_r;                            /* alternative: Sakoe-Chiba half-width */
    int band_cells;                        /* cells the band alternative would cost */
} DecimPlan;

typedef struct {
    int ia[WARP_MAX_N * 2];                /* indices into track A */
    int ib[WARP_MAX_N * 2];                /* aligned indices into track B */
    int len;
    float total_cost;
    float total_weight;                    /* summed transition path lengths */
    float mean_cell_cost;                  /* total_cost / total_weight */
    float max_slope;                       /* steepest step on the path */
    float max_slope_change;                /* biggest |delta slope| between steps */
} WarpPath;

typedef struct {
    float max_abs_err;                     /* sync error upper bound (seconds) */
    float mean_err;
    float rms_err;
    float max_gap;                         /* largest manual anchor spacing */
} SyncReport;

/* ---- DP workspace (file-scope, module is not reentrant) ---- */
static float g_ta[WARP_MAX_N], g_tb[WARP_MAX_N];
static float g_vel_a[WARP_MAX_CHAN][WARP_MAX_N];
static float g_vel_b[WARP_MAX_CHAN][WARP_MAX_N];
static float g_scale[WARP_MAX_CHAN];
static float g_dp[WARP_MAX_N * WARP_MAX_N];
static short g_back[WARP_MAX_N * WARP_MAX_N]; /* -1 unreachable, 0 start, >0 packed */
static float g_prev_slope[WARP_MAX_N * WARP_MAX_N]; /* slope used to reach each cell */

static void track_init(Track* tr) { memset(tr, 0, sizeof(*tr)); }

static void track_set_wrap(Track* tr, int chan, int wrap) {
    if (chan >= 0 && chan < WARP_MAX_CHAN) tr->wrap_chan[chan] = wrap ? 1 : 0;
}

static int track_add(Track* tr, float t, const float* v, int nchan) {
    if (tr->n >= WARP_MAX_N || nchan < 1 || nchan > WARP_MAX_CHAN) return -1;
    if (tr->n > 0 && t <= tr->t[tr->n - 1]) return -1;     /* strictly increasing */
    tr->t[tr->n] = t;
    for (int c = 0; c < nchan; c++) tr->v[c][tr->n] = v[c];
    tr->nchan = nchan;
    tr->n++;
    return tr->n - 1;
}

static WarpConfig warp_config_default(void) {
    WarpConfig cfg = { 1.0f, 0.5f, 0.1f, 0.05f, 0.5f, 2.0f };
    return cfg;
}

static int track_valid(const Track* tr) {
    if (!tr || tr->n < 2 || tr->n > WARP_MAX_N) return 0;
    if (tr->nchan < 1 || tr->nchan > WARP_MAX_CHAN) return 0;
    for (int i = 1; i < tr->n; i++)
        if (!(tr->t[i] > tr->t[i - 1])) return 0;
    return 1;
}

/* shortest signed angular difference, result in [-pi, pi] */
static float wrap_diff(float d) { return atan2f(sinf(d), cosf(d)); }

static float chan_diff(const Track* a, const Track* b, int c, int i, int j) {
    float d = a->v[c][i] - b->v[c][j];
    if (a->wrap_chan[c]) d = wrap_diff(d);                   /* period 2*pi */
    float s = g_scale[c];
    return d / s;
}

/* retained frame count when keeping every kth frame (index 0,k,2k,...),
   plus the last source frame when the sequence does not already end there */
static int retained(int n, int k) {
    int c = (n + k - 1) / k;
    if ((n - 1) % k != 0) c++;
    return c;
}

static DecimPlan make_plan(int na, int nb) {
    DecimPlan plan;
    memset(&plan, 0, sizeof(plan));
    int best_da = na, best_db = nb, best_cells = na * nb, best_score = na + nb + 1;
    for (int da = 1; da <= na; da++) {
        for (int db = 1; db <= nb; db++) {
            int ra = retained(na, da), rb = retained(nb, db);
            if (ra < 2 || rb < 2) continue;
            int cells = ra * rb;
            if (cells > WARP_BUDGET) continue;
            int score = da + db;                              /* discard as little as possible */
            if (score < best_score ||
                (score == best_score && cells < best_cells) ||
                (score == best_score && cells == best_cells && da < best_da)) {
                best_score = score; best_cells = cells;
                best_da = da; best_db = db;
            }
        }
    }
    plan.dec_a = best_da; plan.dec_b = best_db;
    plan.n_a_after = retained(na, best_da);
    plan.n_b_after = retained(nb, best_db);
    plan.cells_after = best_cells;
    int r = (WARP_BUDGET / na - 1) / 2;                      /* Sakoe-Chiba alternative */
    if (r < 0) r = 0;
    if (r > nb - 1) r = nb - 1;
    plan.band_r = r;
    plan.band_cells = na * (nb < 2 * r + 1 ? nb : 2 * r + 1);
    return plan;
}

/* copy every kth frame; the last source frame is always kept */
static int warp_decimate(const Track* src, int k, Track* dst) {
    if (!track_valid(src) || k < 1) return WARP_ERR_INPUT;
    track_init(dst);
    dst->nchan = src->nchan;
    memcpy(dst->wrap_chan, src->wrap_chan, sizeof(dst->wrap_chan));
    for (int i = 0; i < src->n; i += k) {
        float vals[WARP_MAX_CHAN];
        for (int c = 0; c < src->nchan; c++) vals[c] = src->v[c][i];
        if (track_add(dst, src->t[i], vals, src->nchan) < 0) return WARP_ERR_INPUT;
    }
    int last = src->n - 1;
    if (dst->t[dst->n - 1] != src->t[last]) {
        float vals[WARP_MAX_CHAN];
        for (int c = 0; c < src->nchan; c++) vals[c] = src->v[c][last];
        if (track_add(dst, src->t[last], vals, src->nchan) < 0)
            return WARP_ERR_INPUT;
    }
    return WARP_OK;
}

static float cell_cost(const Track* a, const Track* b, const WarpConfig* cfg,
                       int na, int nb, int i, int j) {
    (void)na; (void)nb;
    float pos = 0.0f, vel = 0.0f;
    for (int c = 0; c < a->nchan; c++) {
        float d = chan_diff(a, b, c, i, j);
        pos += d * d;
        float dv = (g_vel_a[c][i] - g_vel_b[c][j]) / g_scale[c];
        vel += dv * dv;
    }
    pos /= a->nchan;
    vel /= a->nchan;
    /* per-cell mismatch; the transition loop multiplies this by the step's
       path length so skipping frames can never reduce accumulated cost */
    return cfg->w_pos * pos + cfg->w_vel * vel;
}

/*
 * Solve the bounded-slope monotone alignment.
 * Path direction: (i, j) means A frame i is played against B frame j.
 * Returns WARP_ERR_BUDGET (plan filled) instead of allocating/running an
 * oversized grid, WARP_ERR_INFEASIBLE if the slope bounds block every path.
 */
static WarpStatus warp_estimate(const Track* a, const Track* b,
                                const WarpConfig* cfg, WarpPath* path,
                                DecimPlan* plan) {
    if (!track_valid(a) || !track_valid(b) || !cfg || !path) return WARP_ERR_INPUT;
    if (a->nchan != b->nchan) return WARP_ERR_INPUT;
    if (!(cfg->slope_min >= 0.0f && cfg->slope_max >= cfg->slope_min))
        return WARP_ERR_INPUT;
    if (cfg->w_pos < 0.0f || cfg->w_vel < 0.0f || cfg->w_bend < 0.0f ||
        cfg->w_smooth < 0.0f)
        return WARP_ERR_INPUT;

    int na = a->n, nb = b->n;
    int cells = na * nb;
    if (cells > WARP_BUDGET) {
        if (plan) *plan = make_plan(na, nb);
        return WARP_ERR_BUDGET;
    }

    float dur_a = a->t[na - 1] - a->t[0];
    float dur_b = b->t[nb - 1] - b->t[0];
    if (!(dur_a > 0.0f && dur_b > 0.0f)) return WARP_ERR_INPUT;

    /* endpoint average slope must itself respect the bounds */
    float end_slope = dur_b / dur_a;
    if (end_slope < cfg->slope_min - WARP_EPS || end_slope > cfg->slope_max + WARP_EPS)
        return WARP_ERR_INFEASIBLE;

    /* normalized time in [0,1] */
    for (int i = 0; i < na; i++) g_ta[i] = (a->t[i] - a->t[0]) / dur_a;
    for (int j = 0; j < nb; j++) g_tb[j] = (b->t[j] - b->t[0]) / dur_b;

    /* per-channel normalization scale; wrap channels normalize by pi */
    for (int c = 0; c < a->nchan; c++) {
        if (a->wrap_chan[c]) { g_scale[c] = WARP_PI; continue; }
        float lo = a->v[c][0], hi = lo;
        for (int i = 0; i < na; i++) {
            lo = fminf(lo, a->v[c][i]); hi = fmaxf(hi, a->v[c][i]);
        }
        for (int j = 0; j < nb; j++) {
            lo = fminf(lo, b->v[c][j]); hi = fmaxf(hi, b->v[c][j]);
        }
        g_scale[c] = (hi - lo) > 1e-8f ? hi - lo : 1.0f;
    }

    /* velocities d(normalized value)/d(normalized time) */
    for (int c = 0; c < a->nchan; c++) {
        for (int i = 0; i < na; i++) {
            int p = i > 0 ? i - 1 : i, q = i < na - 1 ? i + 1 : i;
            float d = a->v[c][q] - a->v[c][p];
            if (a->wrap_chan[c]) d = wrap_diff(d);
            g_vel_a[c][i] = (d / g_scale[c]) / (g_ta[q] - g_ta[p] + WARP_EPS);
        }
        for (int j = 0; j < nb; j++) {
            int p = j > 0 ? j - 1 : j, q = j < nb - 1 ? j + 1 : j;
            float d = b->v[c][q] - b->v[c][p];
            if (b->wrap_chan[c]) d = wrap_diff(d);
            g_vel_b[c][j] = (d / g_scale[c]) / (g_tb[q] - g_tb[p] + WARP_EPS);
        }
    }

    for (int k = 0; k < cells; k++) {
        g_dp[k] = WARP_INF; g_back[k] = -1; g_prev_slope[k] = -1.0f;
    }

    for (int i = 0; i < na; i++) {
        for (int j = 0; j < nb; j++) {
            float here = cell_cost(a, b, cfg, na, nb, i, j);
            if (i == 0 && j == 0) {
                g_dp[0] = here;
                g_back[0] = 0;
                g_prev_slope[0] = 1.0f;
                continue;
            }
            float best = WARP_INF;
            short best_code = -1;
            float best_prev_slope = 1.0f;
            int kmax = i > j ? i : j;
            if (kmax > WARP_MAX_STEP) kmax = WARP_MAX_STEP;
            for (int di = 0; di <= kmax; di++) {
                for (int dj = 0; dj <= kmax; dj++) {
                    if (di == 0 && dj == 0) continue;
                    int pi = i - di, pj = j - dj;
                    if (pi < 0 || pj < 0) continue;
                    int pidx = pi * nb + pj;
                    if (g_back[pidx] < 0) continue;          /* unreachable predecessor */
                    float slope;
                    if (di == 0) slope = WARP_INF;
                    else if (dj == 0) slope = 0.0f;
                    else slope = (dur_b * (g_tb[j] - g_tb[pj]))
                               / (dur_a * (g_ta[i] - g_ta[pi]));
                    if (slope < cfg->slope_min - WARP_EPS) continue;
                    if (slope > cfg->slope_max + WARP_EPS) continue;
                    /* length-normalized DTW: a long jump pays for every
                       frame it skips. Bend = local slope vs identity 1:1;
                       smooth = slope change vs the predecessor step, which
                       breaks ties on periodic signals and rejects jitter. */
                    float w = 0.5f * (float)(di + dj);
                    float bend = fabsf(slope - 1.0f);
                    float smooth = fabsf(slope - g_prev_slope[pidx]);
                    float cand = g_dp[pidx] + w * here
                               + cfg->w_bend * w * bend
                               + cfg->w_smooth * smooth;
                    if (cand < best) {
                        best = cand;
                        best_code = (short)(di * 32 + dj + 1);
                        best_prev_slope = slope;
                    }
                }
            }
            g_dp[i * nb + j] = best;
            g_back[i * nb + j] = best_code;
            g_prev_slope[i * nb + j] = best_prev_slope;
        }
    }

    int end = (na - 1) * nb + (nb - 1);
    if (g_back[end] < 0) return WARP_ERR_INFEASIBLE;

    /* backtrack into reversed scratch, then flip */
    int ri[WARP_MAX_N * 2], rj[WARP_MAX_N * 2];
    int len = 0, i = na - 1, j = nb - 1;
    while (!(i == 0 && j == 0)) {
        if (len >= WARP_MAX_N * 2) return WARP_ERR_INFEASIBLE;
        ri[len] = i; rj[len] = j; len++;
        short code = g_back[i * nb + j];
        int di = (code - 1) / 32, dj = (code - 1) % 32;
        i -= di; j -= dj;
    }
    ri[len] = 0; rj[len] = 0; len++;
    if (len > WARP_MAX_N * 2) return WARP_ERR_INFEASIBLE;

    path->len = len;
    float prev_slope = -1.0f, max_change = 0.0f, max_slope = 0.0f, total_w = 0.0f;
    for (int k = 0; k < len; k++) {
        path->ia[k] = ri[len - 1 - k];
        path->ib[k] = rj[len - 1 - k];
    }
    for (int k = 1; k < len; k++) {
        int di = path->ia[k] - path->ia[k - 1];
        int dj = path->ib[k] - path->ib[k - 1];
        total_w += 0.5f * (float)(di + dj);
        float slope = di == 0 ? WARP_INF
                    : (dur_b * (g_tb[path->ib[k]] - g_tb[path->ib[k - 1]]))
                    / (dur_a * (g_ta[path->ia[k]] - g_ta[path->ia[k - 1]]));
        if (dj == 0) slope = 0.0f;
        if (slope > max_slope) max_slope = slope;
        if (prev_slope >= 0.0f) {
            float ch = fabsf(slope - prev_slope);
            if (ch > max_change) max_change = ch;
        }
        prev_slope = slope;
    }
    path->total_cost = g_dp[end];
    path->total_weight = total_w;
    path->mean_cell_cost = total_w > 0.0f ? g_dp[end] / total_w : g_dp[end];
    path->max_slope = max_slope;
    path->max_slope_change = max_change;
    return WARP_OK;
}

/* B time that the mapping assigns to A time ta (piecewise linear) */
static float warp_map_time(const Track* a, const Track* b,
                           const WarpPath* p, float ta) {
    if (ta <= a->t[p->ia[0]]) return b->t[p->ib[0]];
    if (ta >= a->t[p->ia[p->len - 1]]) return b->t[p->ib[p->len - 1]];
    for (int k = 1; k < p->len; k++) {
        float t0 = a->t[p->ia[k - 1]], t1 = a->t[p->ia[k]];
        if (ta <= t1) {
            float u = t1 > t0 ? (ta - t0) / (t1 - t0) : 0.0f;
            return b->t[p->ib[k - 1]] + u * (b->t[p->ib[k]] - b->t[p->ib[k - 1]]);
        }
    }
    return b->t[p->ib[p->len - 1]];
}

/* resample B onto A's frame grid through the warp; wrap channels interpolate
   along the shortest angular arc */
static WarpStatus warp_resample(const Track* a, const Track* b,
                                const WarpPath* p, Track* out) {
    if (!track_valid(a) || !track_valid(b) || !p || p->len < 2 || !out)
        return WARP_ERR_INPUT;
    track_init(out);
    out->nchan = a->nchan;
    memcpy(out->wrap_chan, a->wrap_chan, sizeof(out->wrap_chan));
    int seg = 1;
    for (int i = 0; i < a->n; i++) {
        float ta = a->t[i];
        while (seg < p->len - 1 && ta > a->t[p->ia[seg]]) seg++;
        float t0 = a->t[p->ia[seg - 1]], t1 = a->t[p->ia[seg]];
        float u = t1 > t0 ? (ta - t0) / (t1 - t0) : 0.0f;
        if (u < 0.0f) u = 0.0f; else if (u > 1.0f) u = 1.0f;
        int j0 = p->ib[seg - 1], j1 = p->ib[seg];
        float vals[WARP_MAX_CHAN];
        for (int c = 0; c < a->nchan; c++) {
            float d = b->v[c][j1] - b->v[c][j0];
            if (a->wrap_chan[c]) d = wrap_diff(d);
            vals[c] = b->v[c][j0] + u * d;
        }
        if (track_add(out, ta, vals, a->nchan) < 0) return WARP_ERR_INPUT;
    }
    return WARP_OK;
}

/*
 * Compare the solved mapping with hand-authored anchor pairs (man_ta[k]
 * aligns with man_tb[k], sorted by man_ta). max_abs_err is the sync error
 * upper bound over the anchor set; both mappings are piecewise linear, so
 * between anchors the deviation is bounded by the endpoint values plus the
 * frame spacing, and max_gap reports the widest interval the bound is
 * interpolated across.
 */
static WarpStatus warp_manual_compare(const Track* a, const Track* b,
                                      const WarpPath* p,
                                      const float* man_ta, const float* man_tb,
                                      int m, SyncReport* rep) {
    if (!track_valid(a) || !track_valid(b) || !p || !man_ta || !man_tb || !rep)
        return WARP_ERR_INPUT;
    if (m < 1 || m > WARP_MAX_ANCHOR) return WARP_ERR_INPUT;
    for (int k = 1; k < m; k++)
        if (!(man_ta[k] > man_ta[k - 1])) return WARP_ERR_INPUT;
    memset(rep, 0, sizeof(*rep));
    float sum = 0.0f, sum2 = 0.0f, max_abs = 0.0f, max_gap = 0.0f;
    for (int k = 0; k < m; k++) {
        float mapped = warp_map_time(a, b, p, man_ta[k]);
        float e = mapped - man_tb[k];
        float ae = fabsf(e);
        if (ae > max_abs) max_abs = ae;
        sum += ae; sum2 += e * e;
        if (k > 0) {
            float gap = man_ta[k] - man_ta[k - 1];
            if (gap > max_gap) max_gap = gap;
        }
    }
    rep->max_abs_err = max_abs;
    rep->mean_err = sum / m;
    rep->rms_err = sqrtf(sum2 / m);
    rep->max_gap = max_gap;
    return WARP_OK;
}

/* ---------------- self test ---------------- */
static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static float ang_wrap(float x) {
    float v = fmodf(x, 2.0f * WARP_PI);
    if (v < 0.0f) v += 2.0f * WARP_PI;
    return v;
}

/* ground-truth warp used by the fixtures: tA = w(tB) */
static float true_warp(float u) { return u + 0.08f * sinf(2.0f * WARP_PI * u); }

/* animation-like multi-shape signal: ease-in-out position, oscillating
   rotation (circular channel), monotone ease scale. Mixed shapes break the
   value ambiguity a single pure sinusoid would suffer from. */
static void anim_signal(float t, float v[3]) {
    v[0] = 0.5f * (1.0f - cosf(WARP_PI * t));        /* position: ease 0..1 */
    v[1] = ang_wrap(4.0f * WARP_PI * t);             /* rotation: circular */
    v[2] = t * t;                                    /* scale: monotone */
}

static void build_tracks(Track* a, Track* b, int na, int nb) {
    track_init(a); track_init(b);
    track_set_wrap(a, 1, 1);
    track_set_wrap(b, 1, 1);
    for (int i = 0; i < na; i++) {
        float t = (float)i / (na - 1);
        float v[3];
        anim_signal(t, v);
        track_add(a, t, v, 3);
    }
    for (int j = 0; j < nb; j++) {
        float u = (float)j / (nb - 1);
        float w = true_warp(u);
        float v[3];
        anim_signal(w, v);
        track_add(b, u, v, 3);
    }
}

int main(void) {
    Track a, b;
    build_tracks(&a, &b, 41, 33);
    WarpConfig cfg = warp_config_default();
    WarpPath path;
    DecimPlan plan;

    WarpStatus st = warp_estimate(&a, &b, &cfg, &path, &plan);
    check("solve_ok", st == WARP_OK);
    check("anchors_corners",
          path.len > 2 && path.ia[0] == 0 && path.ib[0] == 0 &&
          path.ia[path.len - 1] == a.n - 1 && path.ib[path.len - 1] == b.n - 1);

    int mono = 1, slope_ok = 1;
    for (int k = 1; k < path.len; k++) {
        if (!(path.ia[k] > path.ia[k - 1] && path.ib[k] > path.ib[k - 1])) mono = 0;
        float dt_a = a.t[path.ia[k]] - a.t[path.ia[k - 1]];
        float dt_b = b.t[path.ib[k]] - b.t[path.ib[k - 1]];
        float s = dt_b / dt_a;
        if (s < cfg.slope_min - 1e-4f || s > cfg.slope_max + 1e-4f) slope_ok = 0;
    }
    check("strict_monotone", mono);
    check("slope_bounded", slope_ok);
    check("bend_stats_finite",
          isfinite(path.max_slope) && isfinite(path.max_slope_change) &&
          path.max_slope <= cfg.slope_max + 1e-4f && path.max_slope_change >= 0.0f);
    check("cost_positive", path.total_cost > 0.0f && path.mean_cell_cost > 0.0f);

    /* manual alignment: ground-truth anchor pairs tA = w(tB) */
    float man_ta[WARP_MAX_ANCHOR], man_tb[WARP_MAX_ANCHOR];
    int m = 11;
    for (int k = 0; k < m; k++) {
        float u = (float)k / (m - 1);
        man_ta[k] = true_warp(u);
        man_tb[k] = u;
    }
    SyncReport rep;
    check("manual_compare_ok",
          warp_manual_compare(&a, &b, &path, man_ta, man_tb, m, &rep) == WARP_OK);
    /* bound is a few B-frame spacings: on a discrete grid the piecewise
       linear mapping cannot beat its own time resolution */
    check("sync_error_bound", rep.max_abs_err < 2.5f / (b.n - 1) + 0.01f);
    check("sync_report_order",
          rep.max_abs_err >= rep.rms_err - 1e-6f && rep.rms_err >= rep.mean_err - 1e-6f);
    check("max_gap_reported", rep.max_gap > 0.0f);

    /* anchors sitting exactly on solved nodes must give zero error */
    float z_ta[8] = { 0 }, z_tb[8] = { 0 };
    int zm = 0;
    for (int k = 0; k < path.len && zm < 8; k += 5) {
        z_ta[zm] = a.t[path.ia[k]];
        z_tb[zm] = b.t[path.ib[k]];
        zm++;
    }
    SyncReport zrep;
    warp_manual_compare(&a, &b, &path, z_ta, z_tb, zm, &zrep);
    check("anchor_on_node_zero", zrep.max_abs_err < 1e-6f);

    /* resample B back onto A's grid: it should reproduce A. Channel 1 is
       circular, so its error is taken along the shortest angular arc. */
    Track aligned;
    check("resample_ok", warp_resample(&a, &b, &path, &aligned) == WARP_OK);
    check("resample_grid", aligned.n == a.n && aligned.nchan == a.nchan);
    float se_pos = 0.0f, se_ang = 0.0f, se_scl = 0.0f;
    for (int i = 0; i < a.n; i++) {
        float dp = aligned.v[0][i] - a.v[0][i];
        se_pos += dp * dp;
        float da = wrap_diff(aligned.v[1][i] - a.v[1][i]);
        se_ang += da * da;
        float ds = aligned.v[2][i] - a.v[2][i];
        se_scl += ds * ds;
    }
    check("resample_pos_rms_ok", sqrtf(se_pos / a.n) < 0.08f);
    check("resample_angle_rms_ok", sqrtf(se_ang / a.n) < 0.45f);
    check("resample_scale_rms_ok", sqrtf(se_scl / a.n) < 0.08f);

    /* mapping recovered in the time domain must stay within ~2 B-frame
       spacings of the ground-truth warp */
    float map_max = 0.0f;
    for (int i = 0; i < a.n; i++) {
        float ta = a.t[i];
        float mapped_tb = warp_map_time(&a, &b, &path, ta);
        float true_tb = 0.0f;
        for (int jj = 0; jj < b.n - 1; jj++) {
            float u0 = (float)jj / (b.n - 1), u1 = (float)(jj + 1) / (b.n - 1);
            float w0 = true_warp(u0), w1 = true_warp(u1);
            if (ta <= w1 || jj == b.n - 2) {
                float uu = w1 > w0 ? (ta - w0) / (w1 - w0) : 0.0f;
                true_tb = u0 + uu * (u1 - u0);
                break;
            }
        }
        float e = fabsf(mapped_tb - true_tb);
        if (e > map_max) map_max = e;
    }
    check("mapping_time_accurate", map_max < 2.0f / (b.n - 1) + 0.02f);

    /* budget: oversized grid must refuse and report a decimation strategy */
    Track big_a, big_b;
    build_tracks(&big_a, &big_b, 80, 80);
    DecimPlan bp;
    WarpStatus bst = warp_estimate(&big_a, &big_b, &cfg, &path, &bp);
    check("budget_refused", bst == WARP_ERR_BUDGET);
    check("plan_within_budget",
          bp.cells_after <= WARP_BUDGET && bp.dec_a >= 1 && bp.dec_b >= 1 &&
          bp.n_a_after >= 2 && bp.n_b_after >= 2 && bp.band_cells <= WARP_BUDGET);

    Track da, db;
    check("decimate_a", warp_decimate(&big_a, bp.dec_a, &da) == WARP_OK &&
                        da.n == bp.n_a_after);
    check("decimate_b", warp_decimate(&big_b, bp.dec_b, &db) == WARP_OK &&
                        db.n == bp.n_b_after);
    check("retry_after_decimate", warp_estimate(&da, &db, &cfg, &path, NULL) == WARP_OK);

    Track dec10;
    warp_decimate(&a, 3, &dec10);
    check("decimate_keeps_last",
          dec10.n == retained(a.n, 3) &&
          dec10.t[dec10.n - 1] == a.t[a.n - 1]);

    /* slope window above the endpoint average slope (1.0): infeasible */
    WarpConfig tight = cfg;
    tight.slope_min = 1.5f;
    tight.slope_max = 3.0f;
    check("infeasible_reported",
          warp_estimate(&a, &b, &tight, &path, NULL) == WARP_ERR_INFEASIBLE);

    /* input validation */
    Track bad;
    track_init(&bad);
    float v0[2] = { 0.0f, 0.0f };
    track_add(&bad, 0.0f, v0, 2);
    check("reject_short_track", warp_estimate(&bad, &b, &cfg, &path, NULL) == WARP_ERR_INPUT);
    check("reject_null", warp_estimate(NULL, &b, &cfg, &path, NULL) == WARP_ERR_INPUT);
    WarpConfig bad_cfg = cfg;
    bad_cfg.slope_max = 0.1f;
    check("reject_bad_bounds",
          warp_estimate(&a, &b, &bad_cfg, &path, NULL) == WARP_ERR_INPUT);
    Track one_chan;
    track_init(&one_chan);
    for (int i = 0; i < 10; i++) {
        float one[2] = { (float)i, 0.0f };
        track_add(&one_chan, (float)i, one, 1);
    }
    check("reject_nchan_mismatch",
          warp_estimate(&a, &one_chan, &cfg, &path, NULL) == WARP_ERR_INPUT);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
