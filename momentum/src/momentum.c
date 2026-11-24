
/* momentum: hard-stop slide accounting and inverse-mass momentum compensation.
 *
 * A runner enters a fixed 2.0 m stop zone at speed v and brakes hard. Each of
 * the four ground contacts (left/right heel and ball) must supply
 * a_req = v^2 / (2D) of deceleration, but friction supplies at most mu*g. The
 * deficit
 *
 *     a_slip = max(0, a_req - mu*g)
 *
 * is what the contact cannot hold, so it slides forward while engaged.
 * Engagement is staggered (heel before ball, left before right), and a contact
 * that engages late slides for less time:
 *
 *     T      = v / a_req
 *     slip_i = 0.5 * a_slip * max(0, T - delay_i)^2
 *     S      = sum_i slip_i
 *
 * The pelvis travels D = 2.0 m across the stop, so the trigger threshold is
 * half that travel, thr = D/2 = 1.0 m. When S > thr, momentum compensation
 * triggers and the excess E = S - thr is handed back to the joints inversely
 * proportional to their mass (the same inverse-mass weighting constraint/ uses):
 *
 *     quota_j = E * (1/m_j) / sum_k (1/m_k),   sum_j quota_j = E
 *
 * so the residual slide after compensation is exactly
 *
 *     residual = S - sum_j quota_j = min(S, thr).
 *
 * Who carries it. The inverse-mass formula always hands the most to the
 * lightest joints (hands, forearms, feet). But slide happens where the body
 * touches the ground, and only joints on a contact chain (the legs, plus the
 * pelvis that roots both) can deliver a correction there. The module also
 * reports that grounded share as a sensitivity view:
 *
 *     demand_leg   = quota_pelvis/2 + sum(quota on that leg)
 *     grounded_leg = min(demand_leg, slip on that leg's contacts)
 *     residual_grounded = S - sum_leg grounded_leg
 *
 * It is not the requested residual; it answers "how much of the allocated
 * correction can act on the sliding contacts".
 *
 * Everything is deterministic: closed form, fixed inputs, no RNG.
 */
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>

#define N_JOINTS 15
#define N_CONTACTS 4
#define N_RUNS 5
#define N_SKELETONS 2

#define GRAV 9.807
#define FRICTION 0.62
#define STOP_DIST 2.0          /* pelvis travel across the hard stop, m */
#define RUNUP_DIST 5.0         /* fixed acceleration lane, m */

#define CHAIN_NONE 0
#define CHAIN_PELVIS 1         /* roots both legs */
#define CHAIN_LEFT 2
#define CHAIN_RIGHT 3

typedef struct {
    const char* name;
    double mass[N_JOINTS];
    int chain[N_JOINTS];
} SkeletonDef;

typedef struct {
    const char* name;
    int leg;                   /* 0 = left, 1 = right */
    double delay;              /* engagement delay after brake onset, s */
    double x;                  /* planted x in the stop zone, m */
    double y;
} ContactDef;

typedef struct {
    double speed;
} RunSpec;

typedef struct {
    double speed;
    double a_req;              /* deceleration the stop demands, m/s^2 */
    double stop_time;
    double threshold;          /* half the pelvis travel across the stop */
    double a_slip;             /* deceleration the contacts cannot hold */
    double slip[N_CONTACTS];
    double total;              /* S */
    double excess;             /* E */
    int triggered;
} RunResult;

typedef struct {
    double weight[N_JOINTS];   /* (1/m) normalised to 1 */
    double quota[N_JOINTS];    /* inverse-mass allocation of E */
    double quota_sum;
    double residual;           /* S - sum(quota) = min(S, thr) when triggered */
    double demand_leg[2];
    double available_leg[2];
    double grounded_leg[2];
    double grounded[N_JOINTS]; /* quota scaled by the per-leg grounding ratio */
    double grounded_total;
    double residual_grounded;
    double leg_quota_share;    /* fraction of E offered to leg chains */
} CompResult;

/* Joint order is fixed; CSVs and figures rely on it. */
static const char* JOINT_NAME[N_JOINTS] = {
    "pelvis", "spine", "head",
    "upperarm_l", "forearm_l", "hand_l",
    "upperarm_r", "forearm_r", "hand_r",
    "thigh_l", "shin_l", "foot_l",
    "thigh_r", "shin_r", "foot_r"
};

/* Two 70 kg skeletons, same hierarchy, different mass split.
   A: heavy torso, light limbs.  B: heavy limbs, lighter torso. */
static const SkeletonDef SKELETONS[N_SKELETONS] = {
    { "A",
      { 12.0, 13.8, 5.0,
        1.8, 1.1, 0.45,
        1.8, 1.1, 0.45,
        9.5, 5.0, 1.75,
        9.5, 5.0, 1.75 },
      { 1, 0, 0, 0, 0, 0, 0, 0, 0,
        2, 2, 2, 3, 3, 3 } },
    { "B",
      { 13.0, 17.5, 6.5,
        3.5, 1.8, 1.1,
        3.5, 1.8, 1.1,
        6.0, 3.0, 1.1,
        6.0, 3.0, 1.1 },
      { 1, 0, 0, 0, 0, 0, 0, 0, 0,
        2, 2, 2, 3, 3, 3 } }
};

/* Footprints at brake onset: left foot slightly ahead, right slightly behind.
   y is the lateral offset inside the lane. */
static const ContactDef CONTACTS[N_CONTACTS] = {
    { "heel_l", 0, 0.00, RUNUP_DIST + 0.01, +0.13 },
    { "ball_l", 0, 0.03, RUNUP_DIST + 0.19, +0.13 },
    { "heel_r", 1, 0.06, RUNUP_DIST - 0.19, -0.13 },
    { "ball_r", 1, 0.09, RUNUP_DIST - 0.01, -0.13 }
};

/* Five controlled trials, same speeds for both skeletons so the only free
   variable is the mass distribution. */
static const RunSpec RUNS[N_RUNS] = {
    { 5.0 }, { 5.2 }, { 5.6 }, { 6.0 }, { 6.4 }
};

/* ---- simulation ------------------------------------------------------- */

static void simulate(const RunSpec* rs, RunResult* out) {
    memset(out, 0, sizeof(*out));
    const double v = rs->speed;
    out->speed = v;
    out->a_req = v * v / (2.0 * STOP_DIST);
    out->stop_time = v / out->a_req;
    out->threshold = 0.5 * STOP_DIST;

    out->a_slip = out->a_req - FRICTION * GRAV;
    if (out->a_slip < 0.0) out->a_slip = 0.0;

    double total = 0.0;
    for (int i = 0; i < N_CONTACTS; i++) {
        double engaged = out->stop_time - CONTACTS[i].delay;
        double s = engaged > 0.0 ? 0.5 * out->a_slip * engaged * engaged : 0.0;
        out->slip[i] = s;
        total += s;
    }
    out->total = total;
    out->triggered = total > out->threshold;
    out->excess = out->triggered ? total - out->threshold : 0.0;
}

/* ---- compensation ----------------------------------------------------- */

static void inverse_mass_weights(const double* mass, double* w) {
    double sum = 0.0;
    for (int j = 0; j < N_JOINTS; j++) {
        w[j] = 1.0 / mass[j];
        sum += w[j];
    }
    for (int j = 0; j < N_JOINTS; j++) w[j] /= sum;
}

static int leg_of(const int* chain, int j) {
    return chain[j] == CHAIN_LEFT ? 0 : chain[j] == CHAIN_RIGHT ? 1 : -1;
}

static void compensate(const SkeletonDef* sk, const RunResult* rr, CompResult* out) {
    memset(out, 0, sizeof(*out));
    inverse_mass_weights(sk->mass, out->weight);

    for (int j = 0; j < N_JOINTS; j++) {
        out->quota[j] = rr->excess * out->weight[j];
        out->quota_sum += out->quota[j];
    }
    out->residual = rr->total - out->quota_sum;

    for (int c = 0; c < N_CONTACTS; c++)
        out->available_leg[CONTACTS[c].leg] += rr->slip[c];

    /* The pelvis roots both legs, so its quota is offered to each at half;
       this keeps demand_leg[0] + demand_leg[1] <= E. */
    for (int j = 0; j < N_JOINTS; j++) {
        if (sk->chain[j] == CHAIN_PELVIS) {
            out->demand_leg[0] += 0.5 * out->quota[j];
            out->demand_leg[1] += 0.5 * out->quota[j];
        } else {
            int l = leg_of(sk->chain, j);
            if (l >= 0) out->demand_leg[l] += out->quota[j];
        }
    }
    out->leg_quota_share = rr->excess > 0.0
                         ? (out->demand_leg[0] + out->demand_leg[1]) / rr->excess
                         : 0.0;

    double ratio[2] = { 0.0, 0.0 };
    for (int l = 0; l < 2; l++) {
        out->grounded_leg[l] = out->demand_leg[l] < out->available_leg[l]
                             ? out->demand_leg[l] : out->available_leg[l];
        if (out->demand_leg[l] > 0.0)
            ratio[l] = out->grounded_leg[l] / out->demand_leg[l];
        out->grounded_total += out->grounded_leg[l];
    }
    for (int j = 0; j < N_JOINTS; j++) {
        if (sk->chain[j] == CHAIN_PELVIS)
            out->grounded[j] = 0.5 * out->quota[j] * (ratio[0] + ratio[1]);
        else {
            int l = leg_of(sk->chain, j);
            if (l >= 0) out->grounded[j] = out->quota[j] * ratio[l];
        }
    }
    out->residual_grounded = rr->total - out->grounded_total;
}

/* ---- path ------------------------------------------------------------- */

typedef struct { double t, x, v; int phase; } PathSample; /* 0 accel, 1 brake */

static double path_total_time(double speed) {
    double a_acc = speed * speed / (2.0 * RUNUP_DIST);
    double a_req = speed * speed / (2.0 * STOP_DIST);
    return speed / a_acc + speed / a_req;
}

/* Sample the pelvis along the fixed lane: accelerate from rest across
   RUNUP_DIST, then brake across STOP_DIST. per_phase samples per phase, with
   the phase boundary shared. */
static int path_samples(double speed, PathSample* out, int per_phase) {
    double a_acc = speed * speed / (2.0 * RUNUP_DIST);
    double t_acc = speed / a_acc;
    double a_req = speed * speed / (2.0 * STOP_DIST);
    double t_brk = speed / a_req;
    int n = 0;
    for (int i = 0; i < per_phase; i++) {
        double t = t_acc * i / (double)(per_phase - 1);
        out[n].t = t;
        out[n].x = 0.5 * a_acc * t * t;
        out[n].v = a_acc * t;
        out[n].phase = 0;
        n++;
    }
    for (int i = 1; i < per_phase; i++) {
        double t = t_brk * i / (double)(per_phase - 1);
        out[n].t = t_acc + t;
        out[n].x = RUNUP_DIST + speed * t - 0.5 * a_req * t * t;
        out[n].v = speed - a_req * t;
        out[n].phase = 1;
        n++;
    }
    return n;
}

/* ---- self-test -------------------------------------------------------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static int near(double a, double b, double eps) { return fabs(a - b) <= eps; }

static int path_x_monotonic(const PathSample* ps, int n) {
    for (int i = 1; i < n; i++) if (ps[i].x < ps[i - 1].x) return 0;
    return 1;
}

/* ---- CSV export ------------------------------------------------------- */

static void write_accept(const char* dir) {
    char path[512];
    mkdir(dir, 0755);
    FILE* f;

    RunResult rr[N_RUNS];
    for (int r = 0; r < N_RUNS; r++) simulate(&RUNS[r], &rr[r]);

    snprintf(path, sizeof(path), "%s/runs.csv", dir);
    f = fopen(path, "w");
    if (!f) { perror(path); return; }
    fprintf(f, "skeleton,run,trial_label,speed,a_req,stop_time,threshold,"
               "slip_heel_l,slip_ball_l,slip_heel_r,slip_ball_r,"
               "total_slide,triggered,excess,quota_sum,residual,residual_ratio,"
               "leg_quota_share,grounded,residual_grounded\n");
    for (int s = 0; s < N_SKELETONS; s++) {
        for (int r = 0; r < N_RUNS; r++) {
            CompResult cr;
            compensate(&SKELETONS[s], &rr[r], &cr);
            fprintf(f, "%s,%d,%s%d,%.3f,%.6f,%.6f,%.6f,"
                       "%.6f,%.6f,%.6f,%.6f,%.6f,%d,%.6f,%.6f,%.6f,%.6f,"
                       "%.6f,%.6f,%.6f\n",
                    SKELETONS[s].name, r + 1, SKELETONS[s].name, r + 1,
                    rr[r].speed, rr[r].a_req, rr[r].stop_time, rr[r].threshold,
                    rr[r].slip[0], rr[r].slip[1], rr[r].slip[2], rr[r].slip[3],
                    rr[r].total, rr[r].triggered, rr[r].excess,
                    cr.quota_sum, cr.residual,
                    rr[r].total > 0.0 ? cr.residual / rr[r].total : 0.0,
                    cr.leg_quota_share, cr.grounded_total, cr.residual_grounded);
        }
    }
    fclose(f);

    snprintf(path, sizeof(path), "%s/joints.csv", dir);
    f = fopen(path, "w");
    if (!f) { perror(path); return; }
    fprintf(f, "skeleton,joint,mass,inverse_mass_weight,chain,"
               "quota_total,quota_share,grounded_total,grounded_share\n");
    for (int s = 0; s < N_SKELETONS; s++) {
        double quota[N_JOINTS] = { 0 }, grounded[N_JOINTS] = { 0 };
        double qsum = 0.0, gsum = 0.0;
        for (int r = 0; r < N_RUNS; r++) {
            CompResult cr;
            compensate(&SKELETONS[s], &rr[r], &cr);
            for (int j = 0; j < N_JOINTS; j++) {
                quota[j] += cr.quota[j];
                grounded[j] += cr.grounded[j];
            }
        }
        for (int j = 0; j < N_JOINTS; j++) { qsum += quota[j]; gsum += grounded[j]; }
        double w[N_JOINTS];
        inverse_mass_weights(SKELETONS[s].mass, w);
        for (int j = 0; j < N_JOINTS; j++) {
            const char* chain = SKELETONS[s].chain[j] == CHAIN_NONE ? "none"
                              : SKELETONS[s].chain[j] == CHAIN_PELVIS ? "pelvis"
                              : SKELETONS[s].chain[j] == CHAIN_LEFT ? "left" : "right";
            fprintf(f, "%s,%s,%.4f,%.8f,%s,%.6f,%.4f,%.6f,%.4f\n",
                    SKELETONS[s].name, JOINT_NAME[j], SKELETONS[s].mass[j], w[j],
                    chain, quota[j], qsum > 0 ? quota[j] / qsum * 100.0 : 0.0,
                    grounded[j], gsum > 0 ? grounded[j] / gsum * 100.0 : 0.0);
        }
    }
    fclose(f);

    snprintf(path, sizeof(path), "%s/contacts.csv", dir);
    f = fopen(path, "w");
    if (!f) { perror(path); return; }
    fprintf(f, "skeleton,run,contact,leg,delay,x_plant,y,slip,grounded,residual\n");
    for (int s = 0; s < N_SKELETONS; s++) {
        for (int r = 0; r < N_RUNS; r++) {
            CompResult cr;
            compensate(&SKELETONS[s], &rr[r], &cr);
            for (int c = 0; c < N_CONTACTS; c++) {
                int l = CONTACTS[c].leg;
                double ratio = cr.available_leg[l] > 0.0
                             ? cr.grounded_leg[l] / cr.available_leg[l] : 0.0;
                double g = rr[r].slip[c] * ratio;
                fprintf(f, "%s,%d,%s,%s,%.2f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
                        SKELETONS[s].name, r + 1, CONTACTS[c].name,
                        l == 0 ? "left" : "right", CONTACTS[c].delay,
                        CONTACTS[c].x, CONTACTS[c].y,
                        rr[r].slip[c], g, rr[r].slip[c] - g);
            }
        }
    }
    fclose(f);

    snprintf(path, sizeof(path), "%s/path.csv", dir);
    f = fopen(path, "w");
    if (!f) { perror(path); return; }
    fprintf(f, "run,speed,t,x,v,phase\n");
    for (int r = 0; r < N_RUNS; r++) {
        PathSample ps[128];
        int n = path_samples(RUNS[r].speed, ps, 61);
        for (int i = 0; i < n; i++)
            fprintf(f, "%d,%.3f,%.6f,%.6f,%.6f,%d\n", r + 1, RUNS[r].speed,
                    ps[i].t, ps[i].x, ps[i].v, ps[i].phase);
    }
    fclose(f);

    printf("wrote %s/{runs,joints,contacts,path}.csv\n", dir);
}

int main(int argc, char** argv) {
    if (argc >= 3 && strcmp(argv[1], "--accept") == 0) {
        write_accept(argv[2]);
        return 0;
    }

    /* mass tables both total 70 kg, but split differently */
    double ma = 0.0, mb = 0.0;
    for (int j = 0; j < N_JOINTS; j++) {
        ma += SKELETONS[0].mass[j];
        mb += SKELETONS[1].mass[j];
    }
    check("mass_a_70", near(ma, 70.0, 1e-9));
    check("mass_b_70", near(mb, 70.0, 1e-9));
    check("mass_split_differs", SKELETONS[0].mass[5] != SKELETONS[1].mass[5]);

    double w[N_JOINTS];
    inverse_mass_weights(SKELETONS[0].mass, w);
    double wsum = 0.0;
    for (int j = 0; j < N_JOINTS; j++) wsum += w[j];
    check("weights_sum_1", near(wsum, 1.0, 1e-12));
    check("inverse_mass_order", w[5] > w[11] && w[11] > w[9] && w[9] > w[0]);
    check("lightest_dominates", w[5] > 0.2);

    RunResult rr[N_RUNS];
    for (int r = 0; r < N_RUNS; r++) simulate(&RUNS[r], &rr[r]);
    check("threshold_half_stop", near(rr[0].threshold, STOP_DIST / 2.0, 1e-12));
    check("slip_monotonic_speed", rr[0].total < rr[1].total &&
          rr[1].total < rr[2].total && rr[2].total < rr[3].total &&
          rr[3].total < rr[4].total);
    check("phase_delay_reduces_slip", rr[4].slip[0] > rr[4].slip[1] &&
          rr[4].slip[1] > rr[4].slip[2] && rr[4].slip[2] > rr[4].slip[3]);
    check("trigger_first_two_off", rr[0].triggered == 0 && rr[1].triggered == 0);
    check("trigger_rest_on", rr[2].triggered == 1 && rr[3].triggered == 1 &&
          rr[4].triggered == 1);
    check("excess_is_s_minus_thr", near(rr[4].excess, rr[4].total - rr[4].threshold, 1e-12));
    check("no_excess_when_off", near(rr[0].excess, 0.0, 0.0));
    check("stop_time_v_over_a", near(rr[4].stop_time, rr[4].speed / rr[4].a_req, 1e-12));
    check("a_req_stop_dist", near(rr[4].a_req, rr[4].speed * rr[4].speed / (2.0 * STOP_DIST), 1e-12));

    /* the rule: allocation is exactly E, in inverse-mass proportion */
    CompResult cr;
    compensate(&SKELETONS[0], &rr[4], &cr);
    check("quota_sum_excess", near(cr.quota_sum, rr[4].excess, 1e-9));
    check("residual_is_threshold", near(cr.residual, rr[4].threshold, 1e-9));
    check("quota_inverse_mass", cr.quota[5] > cr.quota[11] && cr.quota[11] > cr.quota[9]);
    check("quota_ratio_matches_mass",
          near(cr.quota[5] / cr.quota[9], SKELETONS[0].mass[9] / SKELETONS[0].mass[5], 1e-9));
    check("demand_matches_leg_quota",
          near(cr.demand_leg[0], 0.5 * cr.quota[0] + cr.quota[9] + cr.quota[10] + cr.quota[11], 1e-12) &&
          near(cr.demand_leg[1], 0.5 * cr.quota[0] + cr.quota[12] + cr.quota[13] + cr.quota[14], 1e-12));
    check("grounded_le_demand", cr.grounded_total <= cr.quota_sum + 1e-12);
    check("grounded_le_available",
          cr.grounded_leg[0] <= cr.available_leg[0] + 1e-12 &&
          cr.grounded_leg[1] <= cr.available_leg[1] + 1e-12);
    check("residual_grounded_ge_residual", cr.residual_grounded >= cr.residual - 1e-12);
    check("residual_grounded_le_total", cr.residual_grounded <= rr[4].total + 1e-12);
    check("offchain_never_grounded", cr.grounded[1] == 0.0 && cr.grounded[5] == 0.0);
    check("pelvis_grounded_half_split",
          near(cr.grounded[0], 0.5 * cr.quota[0] *
               (cr.grounded_leg[0] / cr.demand_leg[0] +
                cr.grounded_leg[1] / cr.demand_leg[1]), 1e-9));

    /* mass distribution changes who carries the correction */
    CompResult cra, crb;
    compensate(&SKELETONS[0], &rr[4], &cra);
    compensate(&SKELETONS[1], &rr[4], &crb);
    check("heavier_hands_less_quota", crb.quota[5] < cra.quota[5]);
    check("heavier_legs_ground_more", crb.grounded_total > cra.grounded_total);
    check("residual_grounded_less_for_b", crb.residual_grounded < cra.residual_grounded);

    /* cap binds when the quota exceeds what contacts can absorb */
    RunResult tiny;
    memset(&tiny, 0, sizeof(tiny));
    tiny.total = 10.0;
    tiny.threshold = 1.0;
    tiny.excess = 9.0;
    tiny.triggered = 1;
    for (int c = 0; c < N_CONTACTS; c++) tiny.slip[c] = 0.05;
    CompResult cap;
    compensate(&SKELETONS[0], &tiny, &cap);
    check("cap_binds", near(cap.grounded_leg[0], 0.10, 1e-12) &&
                       near(cap.grounded_leg[1], 0.10, 1e-12));
    check("cap_residual_grounded", near(cap.residual_grounded, 9.8, 1e-12));

    /* untriggered runs compensate nothing */
    CompResult off;
    compensate(&SKELETONS[1], &rr[0], &off);
    check("off_quota_zero", near(off.quota[0], 0.0, 0.0) && near(off.quota[5], 0.0, 0.0));
    check("off_residual_is_total", near(off.residual, rr[0].total, 1e-12));
    check("off_grounded_zero", near(off.grounded_total, 0.0, 0.0));

    /* path: fixed lane, ends at rest on the stop line */
    PathSample ps[128];
    int np = path_samples(6.4, ps, 61);
    check("path_sample_count", np == 121);
    check("path_starts_at_rest", near(ps[0].x, 0.0, 1e-12) && near(ps[0].v, 0.0, 1e-12));
    check("path_accel_end", near(ps[60].x, RUNUP_DIST, 1e-9) && near(ps[60].v, 6.4, 1e-9));
    check("path_stop_line", near(ps[np - 1].x, RUNUP_DIST + STOP_DIST, 1e-9));
    check("path_ends_at_rest", near(ps[np - 1].v, 0.0, 1e-9));
    check("path_monotonic_x", path_x_monotonic(ps, np));
    check("path_total_time", near(path_total_time(6.4), ps[np - 1].t, 1e-12));

    /* determinism */
    RunResult again;
    simulate(&RUNS[4], &again);
    check("deterministic_runs", memcmp(&again, &rr[4], sizeof(RunResult)) == 0);
    CompResult cr2;
    compensate(&SKELETONS[0], &rr[4], &cr2);
    check("deterministic_comp", memcmp(&cr2, &cr, sizeof(CompResult)) == 0);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
