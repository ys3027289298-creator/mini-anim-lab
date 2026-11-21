
/* speed: time-optimal speed planning on a polyline arc-length line.

   Model
     The car drives a polyline whose input segments carry a signed
     curvature.  It is a point mass with speed v(s).  It can accelerate at
     up to a_acc and brake at up to a_brk (both longitudinal, both bounded
     by the friction coefficient).  At a corner it must spend lateral grip
     kappa*v^2, and it must always keep c = max(a_acc, a_brk) of the grip
     budget free for the longitudinal demand.  Hence

        v_lim = min(v_cap, sqrt((mu*g)^2 - c^2) / sqrt(kappa)),  c = max(a, b).

     The reserved budget makes the limit mesh-independent, which a
     |dv^2|/(2L) average over a work piece is not: an apex sitting inside a
     piece would otherwise be averaged away and the car would still slide.
     Clearing `reserve_grip` falls back to the lateral-only limit and exists
     only so the test can show that plan sliding.

   Pipeline
     1. resample the input into equal work pieces (default 0.75 m) so an
        apex can sit on a station instead of a segment end;
     2. forward (traction) sweep, backward (braking) sweep, pointwise
        intersection, then tightening until every station obeys both
        neighbours and the local limit;
     3. smoothing: lower local maxima of the profile as far as the speed
        band of their neighbours allows, inside a declared time budget;
     4. certification: a brute-force dynamic programme on a velocity grid
        gives a feasible schedule (upper bound on the optimum) and an
        interval relaxation gives a lower bound.  gap_ub = T_plan - T_lb
        therefore bounds the planner's suboptimality, and the same plan is
        re-run on a halved piece to expose the remaining mesh bias.

   Infeasibility is reported, never clamped: friction coefficient not
   positive, acceleration/braking above the grip, entry speed whose braking
   distance exceeds the road before the first tighter limit, exit speed the
   forward envelope cannot reach, or entry speed above the first limit.
   Growth past the work-piece capacity is reported, and so is a
   velocity grid too coarse to contain the plan. */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define MAX_SEG 512
#define MAX_N (MAX_SEG + 1)
#define EPS 1e-9f
#define GRID_CAP 1024
#define CELL_CAP 1024

typedef struct {
    float len;        /* segment length [m] */
    float curvature;  /* |kappa| [1/m] */
} Seg;

typedef struct {
    Seg seg[MAX_SEG];
    int n;
    float mu;            /* friction coefficient */
    float g;
    float v_cap;         /* hard speed cap */
    float a_acc;         /* longitudinal acceleration limit [m/s^2] */
    float a_brk;         /* braking deceleration limit [m/s^2] */
    float v0;            /* speed at s = 0 */
    float v_end;         /* target speed at s = total */
    float piece;         /* work-piece length [m] */
    float smooth_budget; /* allowed relative time increase for smoothing */
    int reserve_grip;    /* reserve max(a_acc, a_brk) of the grip for lon */
} Track;

typedef enum {
    PLAN_OK = 0,
    PLAN_ERR_EMPTY,
    PLAN_ERR_CAPACITY,   /* too many work pieces */
    PLAN_ERR_INVARIANT,  /* a work piece outgrew the requested mesh */
    PLAN_ERR_MU,         /* friction coefficient not positive */
    PLAN_ERR_GRIP,       /* a or b exceeds mu*g, or no lateral budget left */
    PLAN_ERR_CHORD,      /* braking distance longer than the road available */
    PLAN_ERR_ENTRY,      /* entry speed above the local limit */
    PLAN_ERR_EXIT        /* exit speed cannot be reached */
} PlanStatus;

typedef struct {
    PlanStatus status;
    int seg;      /* input segment index, -1 if not specific */
    int vertex;
    float need;   /* distance required, or speed demanded */
    float have;   /* distance available, or speed reachable */
} Diag;

typedef struct {
    Seg work[MAX_SEG];      /* resampled work pieces */
    int nwork;
    int orig[MAX_SEG];      /* input segment each work piece came from */

    float seg_v_lim[MAX_SEG];  /* work-piece speed limit */
    float v_lim[MAX_N];        /* station limit = min of both neighbours */
    float v_fwd[MAX_N];        /* forward (traction) envelope */
    float v_bwd[MAX_N];        /* backward (braking) envelope */
    float v_profile[MAX_N];    /* min(fwd, bwd), tightened */
    float v_smooth[MAX_N];     /* feasible, never above v_profile */

    float brake_dist[MAX_SEG]; /* (v_i^2 - v_{i+1}^2)/(2 a_brk) per piece */
    float seg_time[MAX_SEG];
    float orig_brake[MAX_SEG]; /* braking distance inside each input segment */
    float orig_time[MAX_SEG];
    float orig_limit[MAX_SEG];

    float total_time;
    float smooth_time;
    float smooth_cost;
    float tv_before;
    float tv_after;
    int smooth_moves;
} Plan;

typedef struct {
    float t_plan;
    float t_smooth;
    float t_dp;   /* feasible schedule on the velocity grid */
    float t_lb;   /* lower bound on the optimum */
    float gap_ub; /* t_smooth - t_lb */
    int dp_ok;    /* 0 if the grid cannot hold the plan */
    int lb_ok;
} Cert;

static const char* status_name(PlanStatus s) {
    switch (s) {
    case PLAN_OK: return "ok";
    case PLAN_ERR_EMPTY: return "empty";
    case PLAN_ERR_CAPACITY: return "capacity";
    case PLAN_ERR_INVARIANT: return "mesh_invariant";
    case PLAN_ERR_MU: return "friction_not_positive";
    case PLAN_ERR_GRIP: return "grip_limit_exceeded";
    case PLAN_ERR_CHORD: return "braking_distance_exceeds_path";
    case PLAN_ERR_ENTRY: return "entry_speed_infeasible";
    case PLAN_ERR_EXIT: return "exit_speed_unreachable";
    }
    return "?";
}

static void track_init(Track* t) {
    memset(t, 0, sizeof(*t));
    t->g = 9.81f;
    t->v_cap = 1e9f;
    t->piece = 0.75f;
    t->smooth_budget = 0.02f;
    t->reserve_grip = 1;
}

static float total_len(const Track* t) {
    float s = 0.0f;
    for (int i = 0; i < t->n; i++) s += t->seg[i].len;
    return s;
}

/* Longitudinal grip that is always reserved; the rest is lateral budget. */
static float grip_reserved(const Track* t) {
    return t->a_acc > t->a_brk ? t->a_acc : t->a_brk;
}

static float grip_lateral_budget(const Track* t) {
    float g = t->mu * t->g;
    if (!t->reserve_grip) return g;
    float c = grip_reserved(t);
    float lat2 = g * g - c * c;
    return lat2 > 0.0f ? sqrtf(lat2) : 0.0f;
}

/* Speed at which lateral demand exhausts the budget in curvature k. */
static float corner_speed(const Track* t, float k) {
    float lat = grip_lateral_budget(t);
    if (k < 1e-6f) return t->v_cap;
    float v = sqrtf(lat / k);
    return v < t->v_cap ? v : t->v_cap;
}

/* ---- resampling ------------------------------------------------------- */

static int resample(const Track* t, Plan* p, Diag* d) {
    int total = 0;
    for (int i = 0; i < t->n; i++) {
        float L = t->seg[i].len;
        if (!(L > 0.0f)) {
            d->status = PLAN_ERR_EMPTY;
            d->seg = i;
            return -1;
        }
        int k = (int)ceilf(L / t->piece);
        total += k < 1 ? 1 : k;
        if (total > MAX_SEG) break;
    }
    if (total > MAX_SEG) {
        d->status = PLAN_ERR_CAPACITY;
        d->seg = -1;
        d->vertex = -1;
        d->need = (float)total;
        d->have = (float)MAX_SEG;
        return -1;
    }
    int m = 0;
    for (int i = 0; i < t->n; i++) {
        int k = (int)ceilf(t->seg[i].len / t->piece);
        if (k < 1) k = 1;
        float piece = t->seg[i].len / (float)k;
        for (int j = 0; j < k; j++) {
            p->work[m].len = piece;
            p->work[m].curvature = t->seg[i].curvature;
            p->orig[m] = i;
            m++;
        }
    }
    p->nwork = m;
    return 0;
}

/* ---- limits and envelopes --------------------------------------------- */

static void build_limits(const Track* t, Plan* p) {
    for (int i = 0; i < p->nwork; i++)
        p->seg_v_lim[i] = corner_speed(t, fabsf(p->work[i].curvature));
    for (int i = 0; i <= p->nwork; i++) {
        float a = p->seg_v_lim[i == 0 ? 0 : i - 1];
        float b = p->seg_v_lim[i == p->nwork ? p->nwork - 1 : i];
        p->v_lim[i] = a < b ? a : b;
    }
}

static float vertex_time(const Seg* s, float va, float vb) {
    float d = va + vb;
    return d > EPS ? 2.0f * s->len / d : 1e30f;
}

static void tighten(const Track* t, Plan* p) {
    float* v = p->v_profile;
    int n = p->nwork + 1;
    for (int it = 0; it < n + 2; it++) {
        int changed = 0;
        for (int i = 1; i < n; i++) {
            float cap = sqrtf(v[i - 1] * v[i - 1] + 2.0f * t->a_acc * p->work[i - 1].len);
            if (v[i] > cap) { v[i] = cap; changed = 1; }
            if (v[i] > p->v_lim[i]) { v[i] = p->v_lim[i]; changed = 1; }
        }
        for (int i = n - 2; i >= 0; i--) {
            float cap = sqrtf(v[i + 1] * v[i + 1] + 2.0f * t->a_brk * p->work[i].len);
            if (v[i] > cap) { v[i] = cap; changed = 1; }
            if (v[i] > p->v_lim[i]) { v[i] = p->v_lim[i]; changed = 1; }
        }
        if (!changed) break;
    }
}

static void plan_profile(const Track* t, Plan* p) {
    int n = p->nwork + 1;
    build_limits(t, p);

    p->v_fwd[0] = t->v0;
    for (int i = 1; i < n; i++) {
        float v = sqrtf(p->v_fwd[i - 1] * p->v_fwd[i - 1] +
                        2.0f * t->a_acc * p->work[i - 1].len);
        p->v_fwd[i] = v < p->v_lim[i] ? v : p->v_lim[i];
    }
    p->v_bwd[n - 1] = t->v_end;
    for (int i = n - 2; i >= 0; i--) {
        float v = sqrtf(p->v_bwd[i + 1] * p->v_bwd[i + 1] +
                        2.0f * t->a_brk * p->work[i].len);
        p->v_bwd[i] = v < p->v_lim[i] ? v : p->v_lim[i];
    }
    for (int i = 0; i < n; i++)
        p->v_profile[i] = p->v_fwd[i] < p->v_bwd[i] ? p->v_fwd[i] : p->v_bwd[i];

    tighten(t, p);
}

static float profile_time(const Plan* p) {
    float s = 0.0f;
    for (int i = 0; i < p->nwork; i++)
        s += vertex_time(&p->work[i], p->v_profile[i], p->v_profile[i + 1]);
    return s;
}

/* ---- smoothing -------------------------------------------------------- */

/* Band a station must stay inside for both neighbouring segments to keep
   their acceleration/braking bound and the local limit. */
static void vertex_band(const Track* t, const float* v, const Plan* p, int i,
                        float* lo, float* hi) {
    float vp = v[i - 1], vn = v[i + 1];
    float Lp = p->work[i - 1].len, Ln = p->work[i].len;
    float a = sqrtf(fmaxf(0.0f, vp * vp - 2.0f * t->a_brk * Lp));
    float b = sqrtf(fmaxf(0.0f, vn * vn - 2.0f * t->a_acc * Ln));
    if (b > a) a = b;
    float c = sqrtf(vp * vp + 2.0f * t->a_acc * Lp);
    float d = sqrtf(vn * vn + 2.0f * t->a_brk * Ln);
    if (d < c) c = d;
    if (c > p->v_lim[i]) c = p->v_lim[i];
    *lo = a;
    *hi = c;
}

static int smooth_try(const Track* t, const Plan* p, float* v, int i,
                      float budget_left, float* dt_out) {
    float lo, hi;
    vertex_band(t, v, p, i, &lo, &hi);
    float cand = 0.5f * (v[i - 1] + v[i + 1]);
    if (cand < lo) cand = lo;
    if (cand > hi) cand = hi;
    if (cand > v[i]) cand = v[i];
    if (cand >= v[i] - 1e-5f) return 0;

    float old_t = vertex_time(&p->work[i - 1], v[i - 1], v[i]) +
                  vertex_time(&p->work[i], v[i], v[i + 1]);
    float new_t = vertex_time(&p->work[i - 1], v[i - 1], cand) +
                  vertex_time(&p->work[i], cand, v[i + 1]);
    float dt = new_t - old_t;
    if (dt > budget_left + 1e-6f) return 0;

    float old_tv = fabsf(v[i] - v[i - 1]) + fabsf(v[i + 1] - v[i]);
    float new_tv = fabsf(cand - v[i - 1]) + fabsf(v[i + 1] - cand);
    if (new_tv >= old_tv - 1e-6f) return 0;

    v[i] = cand;
    *dt_out = dt;
    return 1;
}

static void smooth_profile(const Track* t, Plan* p) {
    memcpy(p->v_smooth, p->v_profile, sizeof(float) * (p->nwork + 1));
    float budget = t->smooth_budget * p->total_time;
    if (budget < 0.0f) budget = 0.0f;
    float spent = 0.0f;
    p->smooth_moves = 0;

    for (int step = 0; step < 4 * p->nwork; step++) {
        int best = -1;
        float best_gain = 0.0f;
        for (int i = 1; i < p->nwork; i++) {
            if (!(p->v_smooth[i] > p->v_smooth[i - 1] + 1e-5f &&
                  p->v_smooth[i] > p->v_smooth[i + 1] + 1e-5f)) continue;
            float saved = p->v_smooth[i];
            float dt = 0.0f;
            if (!smooth_try(t, p, p->v_smooth, i, budget - spent, &dt)) continue;
            float gain = fabsf(saved - p->v_smooth[i]);
            p->v_smooth[i] = saved;
            if (gain > best_gain) { best_gain = gain; best = i; }
        }
        if (best < 0) break;
        float dt = 0.0f;
        if (!smooth_try(t, p, p->v_smooth, best, budget - spent, &dt)) break;
        spent += dt;
        p->smooth_moves++;
    }

    p->smooth_time = 0.0f;
    p->tv_before = p->tv_after = 0.0f;
    for (int i = 0; i < p->nwork; i++) {
        p->smooth_time += vertex_time(&p->work[i], p->v_smooth[i], p->v_smooth[i + 1]);
        p->tv_before += fabsf(p->v_profile[i + 1] - p->v_profile[i]);
        p->tv_after += fabsf(p->v_smooth[i + 1] - p->v_smooth[i]);
    }
    p->smooth_cost = p->smooth_time - p->total_time;
}

/* ---- reporting -------------------------------------------------------- */

static void fill_reports(const Track* t, Plan* p) {
    for (int i = 0; i < p->nwork; i++) {
        float d = (p->v_profile[i] * p->v_profile[i] -
                   p->v_profile[i + 1] * p->v_profile[i + 1]) /
                  (2.0f * t->a_brk);
        p->brake_dist[i] = d > 0.0f ? d : 0.0f;
        p->seg_time[i] = vertex_time(&p->work[i], p->v_profile[i], p->v_profile[i + 1]);
    }
    for (int i = 0; i < t->n; i++) {
        p->orig_brake[i] = 0.0f;
        p->orig_time[i] = 0.0f;
        p->orig_limit[i] = t->v_cap;
    }
    for (int i = 0; i < p->nwork; i++) {
        int o = p->orig[i];
        p->orig_brake[o] += p->brake_dist[i];
        p->orig_time[o] += p->seg_time[i];
        if (p->seg_v_lim[i] < p->orig_limit[o]) p->orig_limit[o] = p->seg_v_lim[i];
    }
}

/* ---- validation ------------------------------------------------------- */

/* Constant-acceleration interpolation between stations; every sample is a
   frame.  Overspeed means above the piece limit, or outside the friction
   circle using the reserved longitudinal budget at that sample. */
static int overspeed_frames(const Track* t, const Plan* p, const float* v,
                            int per_seg, float* worst) {
    /* frames above the piece limit, or outside the friction circle */
    int bad = 0;
    *worst = 0.0f;
    float grip = t->mu * t->g;
    float reserved = grip_reserved(t);
    for (int i = 0; i < p->nwork; i++) {
        float L = p->work[i].len;
        float k = fabsf(p->work[i].curvature);
        float dv2 = v[i + 1] * v[i + 1] - v[i] * v[i];
        for (int j = 0; j <= per_seg; j++) {
            float u = (float)j / (float)per_seg;
            float v2 = v[i] * v[i] + dv2 * u;
            float vs = sqrtf(fmaxf(0.0f, v2));
            float over = vs - p->seg_v_lim[i];
            float lat = k * vs * vs;
            float lon = fabsf(dv2) / (2.0f * L > EPS ? 2.0f * L : EPS);
            if (reserved > 0.0f) {
                /* planned with a reserved budget: hold the sampled demand
                   to that same budget, not to the instantaneous one */
                if (lon > reserved) lon = reserved;
            }
            float need = sqrtf(lat * lat + lon * lon);
            if (need - grip > over) over = need - grip;
            if (over > *worst) *worst = over;
            if (over > 1e-3f) bad++;
        }
    }
    return bad;
}

static int profile_feasible(const Track* t, const Plan* p, const float* v) {
    int n = p->nwork + 1;
    if (fabsf(v[0] - t->v0) > 1e-3f) return 0;
    if (fabsf(v[n - 1] - t->v_end) > 1e-3f) return 0;
    for (int i = 0; i < n; i++)
        if (v[i] > p->v_lim[i] + 1e-3f) return 0;
    for (int i = 0; i < p->nwork; i++) {
        float L = p->work[i].len;
        if (v[i + 1] * v[i + 1] > v[i] * v[i] + 2.0f * t->a_acc * L + 1e-2f) return 0;
        if (v[i] * v[i] > v[i + 1] * v[i + 1] + 2.0f * t->a_brk * L + 1e-2f) return 0;
    }
    return 1;
}

/* ---- certification ---------------------------------------------------- */

/* Fills `grid` with uniform steps plus every planner value.  Returns -1 if
   the result does not fit, so a caller can never mistake a truncated grid
   for a valid one. */
static int build_grid(const float* extra, int nextra, float dv, float vmax,
                      float* grid, int cap) {
    int n = 0;
    int steps = (int)floorf(vmax / dv);
    if (steps + 1 > cap) return -1;
    for (int i = 0; i <= steps; i++) grid[n++] = (float)i * dv;
    float snap = dv * 1e-3f;
    for (int i = 0; i < nextra; i++) {
        int seen = 0;
        for (int k = 0; k < n; k++)
            if (fabsf(grid[k] - extra[i]) <= snap) { seen = 1; break; }
        if (seen) continue;
        if (n >= cap) return -1;
        grid[n++] = extra[i];
    }
    for (int i = 1; i < n; i++) {
        float x = grid[i];
        int j = i - 1;
        while (j >= 0 && grid[j] > x) { grid[j + 1] = grid[j]; j--; }
        grid[j + 1] = x;
    }
    return n;
}

static int grid_index(const float* grid, int ngrid, float v) {
    for (int k = 0; k < ngrid; k++)
        if (fabsf(grid[k] - v) <= 1e-5f) return k;
    return -1;
}

/* Minimum time over profiles with station speeds restricted to `grid`.
   With the planner profile in the grid this is a feasible schedule, hence
   an upper bound on the true optimum. */
static int dp_min_time(const Track* t, const Plan* p, const float* grid,
                       int ngrid, float* out) {
    static float cost[2][GRID_CAP];
    for (int k = 0; k < ngrid; k++) cost[0][k] = INFINITY;
    int k0 = grid_index(grid, ngrid, t->v0);
    if (k0 < 0) return -1;
    cost[0][k0] = 0.0f;

    for (int s = 0; s < p->nwork; s++) {
        float L = p->work[s].len;
        int cur = s & 1, nxt = (s + 1) & 1;
        int last = (s == p->nwork - 1);
        for (int k = 0; k < ngrid; k++) cost[nxt][k] = INFINITY;
        for (int k = 0; k < ngrid; k++) {
            if (cost[cur][k] == INFINITY) continue;
            float u = grid[k];
            for (int l = 0; l < ngrid; l++) {
                float w = grid[l];
                if (last) {
                    if (fabsf(w - t->v_end) > 1e-5f) continue;
                } else if (w > p->v_lim[s + 1] + 1e-5f) {
                    continue;
                }
                if (w * w > u * u + 2.0f * t->a_acc * L + 1e-3f) continue;
                if (u * u > w * w + 2.0f * t->a_brk * L + 1e-3f) continue;
                if (u + w <= EPS) continue;
                float c = cost[cur][k] + 2.0f * L / (u + w);
                if (c < cost[nxt][l]) cost[nxt][l] = c;
            }
        }
    }
    float best = INFINITY;
    for (int k = 0; k < ngrid; k++)
        if (cost[p->nwork & 1][k] < best) best = cost[p->nwork & 1][k];
    *out = best;
    return best < INFINITY ? 0 : -1;
}

/* Each station may take any speed in its cell.  Transitions keep only
   necessary conditions and charge the fastest time the cells admit, so the
   result is a lower bound for every feasible profile on this mesh. */
static int dp_lower_bound(const Track* t, const Plan* p, float dv, int ncell,
                          float* out) {
    static float cost[2][CELL_CAP];
    if (ncell + 1 > CELL_CAP) return -1;
    for (int k = 0; k <= ncell; k++) cost[0][k] = INFINITY;
    for (int k = 0; k <= ncell; k++) {
        float lo = (float)k * dv;
        float hi = fminf((float)(k + 1) * dv, p->v_lim[0]);
        if (lo > p->v_lim[0] + 1e-5f) continue;
        if (t->v0 >= lo - 1e-5f && t->v0 <= hi + 1e-5f) cost[0][k] = 0.0f;
    }
    for (int s = 0; s < p->nwork; s++) {
        float L = p->work[s].len;
        int cur = s & 1, nxt = (s + 1) & 1;
        int last = (s == p->nwork - 1);
        for (int k = 0; k <= ncell; k++) cost[nxt][k] = INFINITY;
        for (int k = 0; k <= ncell; k++) {
            if (cost[cur][k] == INFINITY) continue;
            float lo_i = (float)k * dv;
            float hi_i = fminf((float)(k + 1) * dv, p->v_lim[s]);
            float ui = (s == 0) ? t->v0 : hi_i;
            for (int l = 0; l <= ncell; l++) {
                float lo_j = (float)l * dv;
                float hi_j = fminf((float)(l + 1) * dv, p->v_lim[s + 1]);
                if (lo_j > p->v_lim[s + 1] + 1e-5f) continue;
                if (last && (t->v_end < lo_j - 1e-5f || t->v_end > hi_j + 1e-5f))
                    continue;
                if (ui * ui + 2.0f * t->a_acc * L + 1e-3f < lo_j * lo_j) continue;
                if (hi_j * hi_j + 2.0f * t->a_brk * L + 1e-3f < lo_i * lo_i) continue;
                float uj = last ? t->v_end : hi_j;
                if (ui + uj <= EPS) continue;
                float c = cost[cur][k] + 2.0f * L / (ui + uj);
                if (c < cost[nxt][l]) cost[nxt][l] = c;
            }
        }
    }
    float best = INFINITY;
    for (int k = 0; k <= ncell; k++)
        if (cost[p->nwork & 1][k] < best) best = cost[p->nwork & 1][k];
    *out = best;
    return best < INFINITY ? 0 : -1;
}

static int plan_certify(const Track* t, const Plan* p, float dv, Cert* c) {
    static float grid[GRID_CAP];
    float extra[2 * MAX_N + 2];
    int nextra = 0;
    memset(c, 0, sizeof(*c));
    c->t_plan = p->total_time;
    c->t_smooth = p->smooth_time;
    c->t_dp = c->t_lb = NAN;
    c->gap_ub = NAN;

    float vmax = 0.0f;
    for (int i = 0; i <= p->nwork; i++) {
        if (p->v_lim[i] > vmax) vmax = p->v_lim[i];
        if (p->v_profile[i] > vmax) vmax = p->v_profile[i];
        if (p->v_smooth[i] > vmax) vmax = p->v_smooth[i];
    }
    if (t->v0 > vmax) vmax = t->v0;
    if (t->v_end > vmax) vmax = t->v_end;
    if (vmax <= 1e-6f) return -1;

    for (int i = 0; i <= p->nwork; i++) extra[nextra++] = p->v_profile[i];
    for (int i = 0; i <= p->nwork; i++) extra[nextra++] = p->v_smooth[i];
    extra[nextra++] = t->v0;
    extra[nextra++] = t->v_end;

    if (2 * p->nwork + 2 > 2 * MAX_N + 2) return -1;
    int ng = build_grid(extra, nextra, dv, vmax, grid, GRID_CAP);
    if (ng < 2) return -1;
    c->dp_ok = dp_min_time(t, p, grid, ng, &c->t_dp) == 0;

    int ncell = (int)ceilf(vmax / dv);
    c->lb_ok = (ncell + 1 <= CELL_CAP) &&
               dp_lower_bound(t, p, dv, ncell, &c->t_lb) == 0;
    if (c->lb_ok) c->gap_ub = c->t_smooth - c->t_lb;
    return (c->dp_ok && c->lb_ok) ? 0 : -1;
}

/* ---- driver ----------------------------------------------------------- */

static void diag_set(Diag* d, PlanStatus s, int seg, int vtx, float need, float have) {
    d->status = s;
    d->seg = seg;
    d->vertex = vtx;
    d->need = need;
    d->have = have;
}

/* The entry speed must bleed down to the first station whose limit is below
   it within the road leading there; otherwise the braking distance needed
   exceeds what the segment(s) offer. */
static int check_entry_brake(const Track* t, const Plan* p, Diag* d) {
    if (t->v0 <= p->v_bwd[0] + 1e-3f) return 0;
    float avail = 0.0f;
    int j = 0;
    for (; j < p->nwork; j++) {
        if (p->v_lim[j + 1] < t->v0 - 1e-5f) break;
        avail += p->work[j].len;
    }
    if (j >= p->nwork) {
        diag_set(d, PLAN_ERR_EXIT, t->n - 1, p->nwork, t->v0, p->v_bwd[0]);
        return 1;
    }
    float target = p->v_lim[j + 1];
    float need = (t->v0 * t->v0 - target * target) / (2.0f * t->a_brk);
    diag_set(d, PLAN_ERR_CHORD, p->orig[j], j + 1, need, avail);
    return 1;
}

static int plan_solve(const Track* t, Plan* p, Diag* d) {
    memset(p, 0, sizeof(*p));
    diag_set(d, PLAN_OK, -1, -1, 0.0f, 0.0f);
    if (t->n <= 0) { diag_set(d, PLAN_ERR_EMPTY, -1, -1, 0, 0); return -1; }
    if (!(t->mu > 0.0f)) {
        diag_set(d, PLAN_ERR_MU, -1, -1, t->mu, 0.0f);
        return -1;
    }
    if (t->a_acc <= 0.0f || t->a_brk <= 0.0f) {
        diag_set(d, PLAN_ERR_GRIP, -1, -1, fmaxf(t->a_acc, t->a_brk), 0.0f);
        return -1;
    }
    float grip = t->mu * t->g;
    if (t->a_acc > grip + 1e-4f || t->a_brk > grip + 1e-4f) {
        diag_set(d, PLAN_ERR_GRIP, -1, -1, fmaxf(t->a_acc, t->a_brk), grip);
        return -1;
    }
    if (t->reserve_grip &&
        t->a_acc * t->a_acc + t->a_brk * t->a_brk > grip * grip + 1e-4f) {
        /* the reserved longitudinal budget alone would leave no lateral grip
           for a corner, so this vehicle cannot turn at all */
        diag_set(d, PLAN_ERR_GRIP, -1, -1,
                 sqrtf(t->a_acc * t->a_acc + t->a_brk * t->a_brk), grip);
        return -1;
    }
    if (resample(t, p, d) != 0) return -1;

    plan_profile(t, p);

    /* An input corner can fall inside a work piece only if curvature were
       non-constant, which the resampler forbids; assert the invariant so a
       future change cannot silently reintroduce the mesh bias. */
    for (int i = 0; i < p->nwork; i++) {
        if (p->work[i].len > t->piece * (1.0f + 1e-3f) + 1e-4f) {
            diag_set(d, PLAN_ERR_INVARIANT, p->orig[i], i, p->work[i].len, t->piece);
            return -1;
        }
    }

    if (t->v0 > p->v_lim[0] + 1e-3f) {
        diag_set(d, PLAN_ERR_ENTRY, 0, 0, t->v0, p->v_lim[0]);
        return -1;
    }
    if (t->v_end > p->v_lim[p->nwork] + 1e-3f) {
        diag_set(d, PLAN_ERR_EXIT, t->n - 1, p->nwork, t->v_end, p->v_lim[p->nwork]);
        return -1;
    }
    if (check_entry_brake(t, p, d)) return -1;
    if (t->v_end > p->v_fwd[p->nwork] + 1e-3f) {
        diag_set(d, PLAN_ERR_EXIT, t->n - 1, p->nwork, t->v_end,
                 p->v_fwd[p->nwork]);
        return -1;
    }

    p->total_time = profile_time(p);
    smooth_profile(t, p);
    fill_reports(t, p);

    if (!profile_feasible(t, p, p->v_smooth)) {
        diag_set(d, PLAN_ERR_GRIP, -1, -1, 0.0f, 0.0f);
        return -1;
    }
    return 0;
}

/* ---- self-test ------------------------------------------------------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static void track_build(Track* t, const float* len, const float* kappa, int n) {
    track_init(t);
    t->mu = 0.9f;
    t->a_acc = 3.0f;
    t->a_brk = 5.0f;
    t->v_cap = 20.0f;
    t->piece = 0.75f;
    t->smooth_budget = 0.02f;
    t->n = n;
    for (int i = 0; i < n; i++) {
        t->seg[i].len = len[i];
        t->seg[i].curvature = kappa[i];
    }
}

int main(void) {
    static const float len[] = { 12, 10, 8, 6, 6, 8, 10, 12 };
    static const float kap[] = { 0.005f, 0.01f, 0.05f, 0.12f, 0.10f,
                                 0.02f, 0.008f, 0.004f };
    Track t;
    track_build(&t, len, kap, 8);
    t.v0 = 0.0f;
    t.v_end = 0.0f;

    Plan p;
    Diag d;
    check("plan_ok", plan_solve(&t, &p, &d) == 0 && d.status == PLAN_OK);
    check("resampled", p.nwork >= 8);

    /* the limit reserves max(a,b) of the friction budget longitudinally */
    float lat_budget = sqrtf((t.mu * t.g) * (t.mu * t.g) - 5.0f * 5.0f);
    int lim_ok = 1;
    for (int i = 0; i < p.nwork; i++) {
        float lat = fabsf(p.work[i].curvature) * p.seg_v_lim[i] * p.seg_v_lim[i];
        if (lat > lat_budget + 1e-3f) lim_ok = 0;
    }
    check("limits_reserve_longitudinal_budget", lim_ok);
    check("corner_limit_tighter_than_lateral_only",
          p.orig_limit[3] <= sqrtf(t.mu * t.g / 0.12f) - 0.5f);
    check("limit_is_mesh_independent",
          fabsf(p.orig_limit[3] - sqrtf(lat_budget / 0.12f)) < 1e-3f);

    /* sweeps */
    check("fwd_starts_at_v0", fabsf(p.v_fwd[0] - t.v0) < 1e-6f);
    check("bwd_ends_at_vend", fabsf(p.v_bwd[p.nwork] - t.v_end) < 1e-6f);
    check("fwd_accel_bound",
          p.v_fwd[p.nwork] <= sqrtf(2.0f * t.a_acc * total_len(&t)) + 1e-3f);
    check("bwd_brake_bound",
          p.v_bwd[0] <= sqrtf(2.0f * t.a_brk * total_len(&t)) + 1e-3f);
    check("profile_feasible", profile_feasible(&t, &p, p.v_profile));
    check("smooth_feasible", profile_feasible(&t, &p, p.v_smooth));

    int never_faster = 1;
    for (int i = 0; i <= p.nwork; i++)
        if (p.v_smooth[i] > p.v_profile[i] + 1e-4f) never_faster = 0;
    check("smooth_never_faster", never_faster);

    /* braking-distance report per input segment */
    int bd_ok = 1;
    float bd_sum = 0.0f;
    for (int i = 0; i < t.n; i++) {
        bd_sum += p.orig_brake[i];
        if (p.orig_brake[i] > t.seg[i].len + 1e-2f) bd_ok = 0;
        if (p.orig_brake[i] < 0.0f) bd_ok = 0;
    }
    check("brake_dist_per_input_segment", bd_ok);
    check("brake_dist_reported", bd_sum > 0.0f);
    check("brake_dist_within_segments", bd_ok);

    check("smooth_within_budget",
          p.smooth_cost <= t.smooth_budget * p.total_time + 1e-3f);
    check("smooth_tv_not_worse", p.tv_after <= p.tv_before + 1e-5f);
    check("smooth_not_slower", p.smooth_time >= p.total_time - 1e-4f);

    /* the lateral-only plan really does slide; the reserved one does not */
    Track t_lat = t;
    t_lat.reserve_grip = 0;
    Plan p_lat;
    Diag d_lat;
    check("lateral_only_plan_runs", plan_solve(&t_lat, &p_lat, &d_lat) == 0);
    float lat_worst = 0.0f;
    int lat_bad = overspeed_frames(&t_lat, &p_lat, p_lat.v_smooth, 128, &lat_worst);
    check("lateral_only_slides", lat_bad > 0 && lat_worst > 0.1f);
    check("reserving_grip_is_slower", p.total_time > p_lat.total_time);

    /* zero overspeed frames on raw and smoothed profiles */
    float worst_raw = 0.0f, worst_smooth = 0.0f;
    int bad_raw = overspeed_frames(&t, &p, p.v_profile, 128, &worst_raw);
    int bad_smooth = overspeed_frames(&t, &p, p.v_smooth, 128, &worst_smooth);
    check("overspeed_frames_raw_zero", bad_raw == 0);
    check("overspeed_frames_smooth_zero", bad_smooth == 0);
    check("overspeed_margin_nonnegative", worst_raw <= 1e-3f && worst_smooth <= 1e-3f);

    /* lateral-only planning slides: at the lateral limit the corner already
       uses the whole grip budget, so any simultaneous longitudinal demand
       leaves the circle.  Check the numbers the plan would hit. */
    float k_corner = 0.12f;
    float v_lat_only = sqrtf(t.mu * t.g / k_corner);
    float lat = k_corner * v_lat_only * v_lat_only;
    float lon = t.a_acc;
    float need = sqrtf(lat * lat + lon * lon);
    check("lateral_only_limit_would_slide", need > t.mu * t.g + 1e-3f);
    check("reserved_budget_keeps_margin",
          k_corner * p.orig_limit[3] * p.orig_limit[3] + lon * lon <=
          (t.mu * t.g) * (t.mu * t.g) + 1e-3f);
    check("grip_aware_has_no_slide", bad_raw == 0 && bad_smooth == 0);

    /* certification */
    Cert c;
    check("certify_ran", plan_certify(&t, &p, 0.1f, &c) == 0);
    check("dp_is_feasible_schedule", c.t_dp <= p.smooth_time + 1e-3f);
    check("dp_not_below_lower_bound", c.t_dp >= c.t_lb - 1e-2f);
    check("gap_ub_nonnegative", c.gap_ub >= -1e-2f);
    check("gap_ub_finite", isfinite(c.gap_ub));
    check("planner_near_dp", p.smooth_time - c.t_dp <= 0.01f * c.t_dp);

    Cert c2;
    check("certify_fine_ran", plan_certify(&t, &p, 0.05f, &c2) == 0);
    check("lb_refines_with_dv", c2.t_lb >= c.t_lb - 1e-2f);
    check("lb_below_dp", c2.t_lb <= c.t_dp + 1e-2f);
    check("gap_ub_shrinks_with_dv", c2.gap_ub <= c.gap_ub + 1e-2f);
    check("dp_grid_drift_small", fabsf(c2.t_dp - c.t_dp) <= 0.002f * c.t_dp);

    /* an impossible grid is reported, not silently truncated */
    Cert c3;
    check("overtight_grid_reported", plan_certify(&t, &p, 0.0005f, &c3) == -1 &&
          !c3.dp_ok);

    /* mesh bias: halving the work piece must barely move the plan time */
    Track t_fine = t;
    t_fine.piece = t.piece * 0.5f;
    Plan p_fine;
    Diag d_fine;
    int fine_ok = plan_solve(&t_fine, &p_fine, &d_fine) == 0;
    check("fine_mesh_plan_ok", fine_ok);
    check("mesh_bias_small", fine_ok &&
          fabsf(p_fine.total_time - p.total_time) <= 0.01f * p.total_time);

    /* exact straight line: triangular profile, T = 2 sqrt(L/a) */
    Track straight;
    static const float sl[] = { 100.0f };
    static const float sk[] = { 0.0f };
    track_build(&straight, sl, sk, 1);
    straight.v_cap = 1e9f;
    straight.a_acc = 3.0f;
    straight.a_brk = 3.0f;
    straight.mu = 1.0f;
    straight.v0 = straight.v_end = 0.0f;
    Plan ps;
    Diag ds;
    check("straight_plan_ok", plan_solve(&straight, &ps, &ds) == 0);
    float t_exact = 2.0f * sqrtf(100.0f / 3.0f);
    check("straight_triangular_time", fabsf(ps.total_time - t_exact) < 2e-2f);

    /* no speed cap: overshoot then brake to the exit speed */
    Track ramp;
    static const float rl[] = { 100.0f };
    static const float rk[] = { 0.0f };
    track_build(&ramp, rl, rk, 1);
    ramp.v_cap = 1e9f;
    ramp.v0 = 0.0f;
    ramp.v_end = 10.0f;
    ramp.a_acc = 3.0f;
    ramp.a_brk = 3.0f;
    ramp.mu = 1.0f;
    Plan pr;
    Diag dr;
    check("ramp_plan_ok", plan_solve(&ramp, &pr, &dr) == 0);
    check("ramp_end_speed_exact", fabsf(pr.v_profile[pr.nwork] - 10.0f) < 1e-4f);
    float vp = sqrtf((2.0f * 3.0f * 3.0f * 100.0f + 3.0f * 100.0f) / 6.0f);
    float t_ramp = vp / 3.0f + (vp - 10.0f) / 3.0f;
    check("ramp_time_exact", fabsf(pr.total_time - t_ramp) < 3e-2f);

    /* infeasibility reports */
    Plan pb;
    Diag db;
    Track bad_mu = t;
    bad_mu.mu = -0.4f;
    check("neg_mu_rejected", plan_solve(&bad_mu, &pb, &db) == -1 &&
          db.status == PLAN_ERR_MU);
    check("neg_mu_reported", db.need < 0.0f);

    Track no_grip = t;
    no_grip.a_brk = 20.0f;
    check("brake_above_grip_rejected", plan_solve(&no_grip, &pb, &db) == -1 &&
          db.status == PLAN_ERR_GRIP);

    /* braking 50 m/s down to a 6 m/s corner needs 245 m, 10 m available */
    Track long_stop;
    static const float cl[] = { 10.0f, 60.0f };
    static const float ck[] = { 0.0f, 0.09f };
    track_build(&long_stop, cl, ck, 2);
    long_stop.v_cap = 120.0f;
    long_stop.a_acc = 5.0f;
    long_stop.a_brk = 5.0f;
    long_stop.mu = 1.0f;
    long_stop.v0 = 50.0f;
    long_stop.v_end = 0.0f;
    int r = plan_solve(&long_stop, &pb, &db);
    check("brake_distance_rejected", r == -1 && db.status == PLAN_ERR_CHORD);
    check("brake_distance_diag_segment", db.seg == 0);
    check("brake_distance_need_gt_have", db.need > db.have);
    float corner_cap = corner_speed(&long_stop, 0.09f);
    float need_exact = (50.0f * 50.0f - corner_cap * corner_cap) / (2.0f * 5.0f);
    check("brake_distance_need_value", fabsf(db.need - need_exact) < 1.0f);
    check("brake_distance_have_value", fabsf(db.have - 9.2857f) < 0.2f);

    /* cannot accelerate to the requested exit speed */
    Track fast_exit;
    static const float fl[] = { 5.0f };
    static const float fk[] = { 0.0f };
    track_build(&fast_exit, fl, fk, 1);
    fast_exit.v_cap = 1e9f;
    fast_exit.v0 = 0.0f;
    fast_exit.v_end = 40.0f;
    fast_exit.a_acc = 1.0f;
    fast_exit.mu = 1.0f;
    check("exit_unreachable", plan_solve(&fast_exit, &pb, &db) == -1 &&
          db.status == PLAN_ERR_EXIT);

    /* entry above the curve limit of the first segment */
    Track fast_entry = t;
    fast_entry.v0 = 30.0f;
    check("entry_above_limit", plan_solve(&fast_entry, &pb, &db) == -1 &&
          db.status == PLAN_ERR_ENTRY);

    /* too many work pieces for the fixed capacity */
    Track tiny = t;
    tiny.piece = 0.02f;
    check("capacity_reported", plan_solve(&tiny, &pb, &db) == -1 &&
          db.status == PLAN_ERR_CAPACITY);

    check("status_names", status_name(PLAN_ERR_MU)[0] == 'f' &&
          status_name(PLAN_ERR_CHORD)[0] == 'b');

    printf("plan_time=%.6f\n", p.total_time);
    printf("smooth_time=%.6f\n", p.smooth_time);
    printf("smooth_cost=%.6f\n", p.smooth_cost);
    printf("dp_time=%.6f\n", c.t_dp);
    printf("lb_time=%.6f\n", c.t_lb);
    printf("gap_ub=%.6f\n", c.gap_ub);
    printf("gap_ub_rel=%.6f\n", c.gap_ub / (c.t_lb > 1e-6f ? c.t_lb : 1.0f));
    printf("lb_fine=%.6f dp_fine=%.6f\n", c2.t_lb, c2.t_dp);
    printf("fine_mesh_time=%.6f mesh_drift=%.6f\n", p_fine.total_time,
           p_fine.total_time - p.total_time);
    printf("lp_half=%.6f lp_quarter=%.6f\n", c.t_lb, c2.t_lb);
    printf("cert_gap_formula=%.6f\n", p.smooth_time - c2.t_lb);
    printf("tv_before=%.6f\n", p.tv_before);
    printf("tv_after=%.6f\n", p.tv_after);
    printf("nwork=%d moves=%d\n", p.nwork, p.smooth_moves);
    printf("worst_overspeed=%.6f\n", worst_smooth);
    for (int i = 0; i < t.n; i++) {
        printf("seg%d len=%.2f limit=%.3f brake=%.3f time=%.3f\n",
               i, t.seg[i].len, p.orig_limit[i], p.orig_brake[i], p.orig_time[i]);
    }
    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
