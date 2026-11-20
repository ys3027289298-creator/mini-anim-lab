
/* retarget: root motion extraction, skeleton-ratio scaling and stride alignment.
 *
 * A clip is a frame list; each frame carries the world root (hips) position and
 * both foot positions. Extraction splits every frame into a root track and a
 * pose track (feet relative to the root). Retargeting then:
 *
 *   - scales the root track and the pose track by target_leg_length /
 *     source_leg_length, so a longer or shorter rig keeps its footfalls on the
 *     ground instead of sinking through it;
 *   - measures the source contact phase (first left-foot contact, in cycles) and
 *     resamples the clip shifted so contact happens at cfg.target_contact_phase.
 *     The resample wraps over one gait period and advances the root by the
 *     stride vector of each wrapped cycle, so planted feet stay planted;
 *   - rejects a zero / non-finite skeleton proportion with RET_ERR_SCALE, a
 *     phase difference wider than half a cycle with RET_ERR_PHASE, and plant
 *     slip beyond cfg.slide_limit with RET_ERR_SLIDE. Output is not written on
 *     any error.
 *
 * retarget_compare() measures the retargeted footfalls against a hand-retargeted
 * clip so the landing points can be checked directly.
 *
 * With cfg.enabled == 0 nothing is retargeted: the clip is copied verbatim.
 */
#include <stdio.h>
#include <string.h>
#include <math.h>

#define MAX_FRAMES 128
#define MAX_CONTACTS (MAX_FRAMES / 2 + 1)
#define CONTACT_EPS 0.02f
#define PI 3.14159265f

#define RET_OK 0
#define RET_ERR_INPUT (-1)
#define RET_ERR_SCALE (-2)
#define RET_ERR_PHASE (-3)
#define RET_ERR_SLIDE (-4)

typedef struct { float x, y; } V2;
typedef struct { V2 root; V2 foot[2]; } Frame;
typedef struct { Frame f[MAX_FRAMES]; int n; } Clip;

/* extracted motion: root track plus root-relative pose track */
typedef struct { V2 root[MAX_FRAMES]; V2 pose[2][MAX_FRAMES]; int n; } Motion;

typedef struct { int first, last; } Contact;

typedef struct {
    int enabled;
    float source_leg_length;
    float target_leg_length;
    float target_contact_phase;   /* cycles in [0,1) */
    float phase_tolerance;        /* cycles; <= 0 means half a cycle */
    float slide_limit;            /* target-scale world units; <= 0 disables */
} RetargetCfg;

typedef struct {
    float scale;
    int period;                   /* source gait period, frames */
    float source_phase;           /* left contact phase at frame 0, cycles */
    float phase_offset;           /* target minus source, cycles */
    int shift_frames;             /* source-time shift used for resampling */
    float slide;                  /* worst planted-foot drift */
    float landing_error;          /* filled by retarget_compare */
    int contacts_compared;
} RetargetStats;

/* ---------------------------------------------------------------- vector helpers */

static V2 v2(float x, float y) { V2 r = { x, y }; return r; }
static V2 vadd(V2 a, V2 b) { return v2(a.x + b.x, a.y + b.y); }
static V2 vsub(V2 a, V2 b) { return v2(a.x - b.x, a.y - b.y); }
static V2 vscale(V2 a, float s) { return v2(a.x * s, a.y * s); }
static float vlen(V2 v) { return sqrtf(v.x * v.x + v.y * v.y); }
static int finite2(V2 v) { return isfinite(v.x) && isfinite(v.y); }

static int floor_div(int a, int b) {
    int q = a / b;
    return (a % b != 0 && ((a < 0) != (b < 0))) ? q - 1 : q;
}

static int mod_pos(int a, int b) {
    int r = a % b;
    return r < 0 ? r + b : r;
}

/* ---------------------------------------------------------------- extraction */

static int motion_extract(const Clip* c, Motion* m) {
    if (!c || !m || c->n <= 0 || c->n > MAX_FRAMES) return RET_ERR_INPUT;
    m->n = c->n;
    for (int i = 0; i < c->n; i++) {
        if (!finite2(c->f[i].root) || !finite2(c->f[i].foot[0]) || !finite2(c->f[i].foot[1]))
            return RET_ERR_INPUT;
        m->root[i] = c->f[i].root;
        m->pose[0][i] = vsub(c->f[i].foot[0], c->f[i].root);
        m->pose[1][i] = vsub(c->f[i].foot[1], c->f[i].root);
    }
    return RET_OK;
}

/* reassemble: root translation and pose each go through the limb ratio */
static void motion_compose(const Motion* m, float scale, Clip* out) {
    out->n = m->n;
    for (int i = 0; i < m->n; i++) {
        V2 root = vscale(m->root[i], scale);
        out->f[i].root = root;
        out->f[i].foot[0] = vadd(root, vscale(m->pose[0][i], scale));
        out->f[i].foot[1] = vadd(root, vscale(m->pose[1][i], scale));
    }
}

/* ---------------------------------------------------------------- contacts */

static float min_foot_y(const Clip* c, int leg) {
    float best = INFINITY;
    for (int i = 0; i < c->n; i++)
        if (c->f[i].foot[leg].y < best) best = c->f[i].foot[leg].y;
    return best;
}

static int foot_contacts(const Clip* c, int leg, float eps, Contact* out, int max) {
    float floor_y = min_foot_y(c, leg);
    int count = 0, inside = 0;
    for (int i = 0; i < c->n; i++) {
        int down = c->f[i].foot[leg].y <= floor_y + eps;
        if (down && !inside) {
            if (count >= max) return count;
            out[count].first = i;
            out[count].last = i;
            inside = 1;
        } else if (down) {
            out[count].last = i;
        } else if (inside) {
            count++;
            inside = 0;
        }
    }
    if (inside) count++;
    return count;
}

static float max_slide(const Clip* c, float eps) {
    float worst = 0.0f;
    for (int leg = 0; leg < 2; leg++) {
        Contact iv[MAX_CONTACTS];
        int n = foot_contacts(c, leg, eps, iv, MAX_CONTACTS);
        for (int k = 0; k < n; k++) {
            V2 planted = c->f[iv[k].first].foot[leg];
            for (int i = iv[k].first; i <= iv[k].last; i++) {
                float d = vlen(vsub(c->f[i].foot[leg], planted));
                if (d > worst) worst = d;
            }
        }
    }
    return worst;
}

/* ---------------------------------------------------------------- retarget */

static int retarget_run(const Clip* src, const RetargetCfg* cfg, Clip* out, RetargetStats* st) {
    if (!src || !cfg || !out || !st) return RET_ERR_INPUT;
    memset(st, 0, sizeof(*st));
    st->scale = 1.0f;

    if (!cfg->enabled) {                      /* playback stays identical */
        if (src->n <= 0 || src->n > MAX_FRAMES) return RET_ERR_INPUT;
        *out = *src;
        return RET_OK;
    }
    if (src->n <= 0 || src->n > MAX_FRAMES) return RET_ERR_INPUT;

    /* a zero or unusable skeleton proportion is an error, never a silent divide */
    if (!isfinite(cfg->source_leg_length) || !isfinite(cfg->target_leg_length) ||
        cfg->source_leg_length <= 0.0f || cfg->target_leg_length <= 0.0f)
        return RET_ERR_SCALE;
    float scale = cfg->target_leg_length / cfg->source_leg_length;
    if (!isfinite(scale) || scale <= 0.0f) return RET_ERR_SCALE;

    Motion m;
    int rc = motion_extract(src, &m);
    if (rc != RET_OK) return rc;

    Contact iv[MAX_CONTACTS];
    int contacts = foot_contacts(src, 0, CONTACT_EPS, iv, MAX_CONTACTS);
    if (contacts < 2) return RET_ERR_INPUT;
    int period = iv[1].first - iv[0].first;
    if (period <= 0 || period >= src->n) return RET_ERR_INPUT;

    float source_phase = (float)mod_pos(iv[0].first, period) / (float)period;
    float offset = cfg->target_contact_phase - source_phase;
    float tol = cfg->phase_tolerance > 0.0f ? cfg->phase_tolerance : 0.5f;
    if (fabsf(offset) > tol) return RET_ERR_PHASE;

    int shift = (int)lroundf(-offset * (float)period);   /* source-time shift */

    /* one period of pose, with the root advanced by the stride each wrap */
    Motion tiled;
    tiled.n = src->n;
    V2 stride = vsub(m.root[iv[1].first], m.root[iv[0].first]);
    for (int i = 0; i < src->n; i++) {
        int t = i + shift;
        int cycle = floor_div(t, period);
        int phase = mod_pos(t, period);
        V2 root = vadd(m.root[phase], vscale(stride, (float)cycle));
        tiled.root[i] = root;
        tiled.pose[0][i] = m.pose[0][phase];
        tiled.pose[1][i] = m.pose[1][phase];
    }

    Clip tmp;
    motion_compose(&tiled, scale, &tmp);

    float slide = max_slide(&tmp, CONTACT_EPS * scale);
    if (cfg->slide_limit > 0.0f && slide > cfg->slide_limit) return RET_ERR_SLIDE;

    *out = tmp;
    st->scale = scale;
    st->period = period;
    st->source_phase = source_phase;
    st->phase_offset = offset;
    st->shift_frames = shift;
    st->slide = slide;
    return RET_OK;
}

/* landing points of the retargeted clip against a hand-retargeted reference */
static int retarget_compare(const Clip* out, const Clip* manual, float eps,
                            RetargetStats* st) {
    if (!out || !manual || !st) return RET_ERR_INPUT;
    float worst = 0.0f;
    int compared = 0;
    for (int leg = 0; leg < 2; leg++) {
        Contact a[MAX_CONTACTS], b[MAX_CONTACTS];
        int na = foot_contacts(out, leg, eps, a, MAX_CONTACTS);
        int nb = foot_contacts(manual, leg, eps, b, MAX_CONTACTS);
        for (int k = 0; k < na && k < nb; k++) {
            V2 p = out->f[a[k].first].foot[leg];
            float best = INFINITY;
            for (int j = 0; j < nb; j++) {
                float d = vlen(vsub(p, manual->f[b[j].first].foot[leg]));
                if (d < best) best = d;
            }
            if (best > worst) worst = best;
            compared++;
        }
    }
    st->landing_error = worst;
    st->contacts_compared = compared;
    return RET_OK;
}

/* ---------------------------------------------------------------- self test */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static void put(Clip* c, int i, V2 root, V2 lf, V2 rf) {
    c->f[i].root = root;
    c->f[i].foot[0] = lf;
    c->f[i].foot[1] = rf;
}

/* 24-frame walk, 12-frame period and 3 units of travel per period; drift is the
 * planted-foot slip per frame (0 = clean contact). */
static Clip walk_clip(float drift) {
    Clip c;
    memset(&c, 0, sizeof(c));
    c.n = 24;
    for (int i = 0; i < 24; i++) {
        int k = i / 12, ph = i % 12;
        V2 l, r;
        if (ph <= 5) {
            l = v2(3.0f * k + drift * (float)ph, 0.0f);
        } else {
            float u = (float)(ph - 5) / 6.5f;
            l = v2(3.0f * k + 3.0f * u, 0.4f * sinf(PI * u));
        }
        if (ph >= 6) {
            r = v2(3.0f * k + 1.5f, 0.0f);
        } else {
            float u = (float)(ph + 1) / 6.5f;
            r = v2(3.0f * (k - 1) + 1.5f + 3.0f * u, 0.4f * sinf(PI * u));
        }
        put(&c, i, v2(0.25f * (float)i, 0.0f), l, r);
    }
    return c;
}

static Clip scale_clip(const Clip* s, float k) {
    Clip c = *s;
    for (int i = 0; i < c.n; i++) {
        c.f[i].root = vscale(c.f[i].root, k);
        c.f[i].foot[0] = vscale(c.f[i].foot[0], k);
        c.f[i].foot[1] = vscale(c.f[i].foot[1], k);
    }
    return c;
}

int main(void) {
    Clip src = walk_clip(0.0f);

    RetargetCfg cfg;
    memset(&cfg, 0, sizeof(cfg));
    cfg.enabled = 1;
    cfg.source_leg_length = 1.0f;
    cfg.target_leg_length = 1.5f;
    cfg.target_contact_phase = 0.0f;
    cfg.slide_limit = 0.05f;

    /* root/pose split reconstructs the world foot exactly */
    Motion m;
    int rc = motion_extract(&src, &m);
    check("extract_ok", rc == RET_OK && m.n == src.n);
    check("extract_roundtrip",
          fabsf(m.root[7].x + m.pose[0][7].x - src.f[7].foot[0].x) < 1e-5f &&
          fabsf(m.root[7].y + m.pose[0][7].y - src.f[7].foot[0].y) < 1e-5f);

    Contact iv[MAX_CONTACTS];
    int lc = foot_contacts(&src, 0, CONTACT_EPS, iv, MAX_CONTACTS);
    check("contacts_left", lc == 2 && iv[0].first == 0 && iv[1].first == 12);
    check("contacts_right", foot_contacts(&src, 1, CONTACT_EPS, iv, MAX_CONTACTS) == 2);

    /* nominal retarget: target rig 1.5x, contact already aligned */
    Clip out;
    RetargetStats st;
    memset(&out, 0, sizeof(out));
    rc = retarget_run(&src, &cfg, &out, &st);
    check("run_ok", rc == RET_OK);
    check("scale_ratio", fabsf(st.scale - 1.5f) < 1e-5f);
    check("scaled_root", fabsf(out.f[8].root.x - src.f[8].root.x * 1.5f) < 1e-4f);
    check("scaled_foot", fabsf(out.f[2].foot[0].x - src.f[2].foot[0].x * 1.5f) < 1e-4f);
    check("phase_aligned", fabsf(st.phase_offset) < 1e-5f && st.shift_frames == 0);
    check("period_found", st.period == 12);
    check("slide_zero", st.slide < 1e-5f);
    check("feet_grounded", fabsf(min_foot_y(&out, 0)) < 1e-5f &&
                            fabsf(min_foot_y(&out, 1)) < 1e-5f);

    /* hand-retargeted landings match exactly */
    Clip manual = scale_clip(&src, 1.5f);
    RetargetStats cs = st;
    rc = retarget_compare(&out, &manual, CONTACT_EPS * 1.5f, &cs);
    check("compare_ok", rc == RET_OK && cs.contacts_compared == 4);
    check("landing_match", cs.landing_error < 1e-4f);

    Clip manual_bad = manual;
    manual_bad.f[0].foot[0].x += 0.1f;
    RetargetStats csb = st;
    retarget_compare(&out, &manual_bad, CONTACT_EPS * 1.5f, &csb);
    check("landing_error_seen", fabsf(csb.landing_error - 0.1f) < 1e-3f);

    /* quarter-cycle phase request moves contact to frame 3 and keeps plants still */
    RetargetCfg shifted = cfg;
    shifted.target_contact_phase = 0.25f;
    Clip out_shift;
    RetargetStats sts;
    memset(&out_shift, 0, sizeof(out_shift));
    rc = retarget_run(&src, &shifted, &out_shift, &sts);
    int sc = foot_contacts(&out_shift, 0, CONTACT_EPS * 1.5f, iv, MAX_CONTACTS);
    check("phase_shift_ok", rc == RET_OK && sts.shift_frames == -3);
    check("phase_shift_moved_contact", sc >= 2 && iv[0].first == 3);
    check("phase_shift_no_slide", sts.slide < 1e-5f);

    /* exactly half a cycle is still aligned, one step past it is not */
    RetargetCfg half = cfg;
    half.target_contact_phase = 0.5f;
    Clip out_half;
    RetargetStats sth;
    memset(&out_half, 0, sizeof(out_half));
    rc = retarget_run(&src, &half, &out_half, &sth);
    check("half_cycle_ok", rc == RET_OK && sth.shift_frames == -6);

    RetargetCfg bad_phase = cfg;
    bad_phase.target_contact_phase = 0.6f;
    Clip untouched = walk_clip(0.0f);
    Clip sentinel = walk_clip(0.0f);
    rc = retarget_run(&src, &bad_phase, &untouched, &st);
    check("phase_error", rc == RET_ERR_PHASE);
    check("error_leaves_output", memcmp(&untouched, &sentinel, sizeof(Clip)) == 0);

    /* negative offset direction: source contact phase is 0.25 and target asks 0 */
    Clip shifted_src = src;
    for (int i = 0; i < shifted_src.n; i++)
        put(&shifted_src, i, src.f[(i + 21) % src.n].root,
            src.f[(i + 21) % src.n].foot[0], src.f[(i + 21) % src.n].foot[1]);
    int sp = foot_contacts(&shifted_src, 0, CONTACT_EPS, iv, MAX_CONTACTS);
    check("shifted_source_phase", sp >= 2 && iv[0].first == 3 && iv[1].first == 15);
    RetargetCfg neg = cfg;
    neg.target_contact_phase = 0.0f;
    neg.phase_tolerance = 0.2f;
    check("phase_error_negative",
          retarget_run(&shifted_src, &neg, &untouched, &st) == RET_ERR_PHASE);
    /* symmetric: +0.5 is also outside a 0.25 tolerance */
    RetargetCfg pos = cfg;
    pos.target_contact_phase = 0.75f;
    pos.phase_tolerance = 0.25f;
    check("phase_error_beyond_half_cycle",
          retarget_run(&shifted_src, &pos, &untouched, &st) == RET_ERR_PHASE);

    /* zero skeleton proportion on either side is refused */
    RetargetCfg zero_src = cfg;
    zero_src.source_leg_length = 0.0f;
    check("zero_source_scale", retarget_run(&src, &zero_src, &untouched, &st) == RET_ERR_SCALE);
    RetargetCfg zero_tgt = cfg;
    zero_tgt.target_leg_length = 0.0f;
    check("zero_target_scale", retarget_run(&src, &zero_tgt, &untouched, &st) == RET_ERR_SCALE);
    RetargetCfg neg_src = cfg;
    neg_src.source_leg_length = -1.0f;
    check("negative_source_scale",
          retarget_run(&src, &neg_src, &untouched, &st) == RET_ERR_SCALE);
    RetargetCfg nan_tgt = cfg;
    nan_tgt.target_leg_length = NAN;
    check("nan_target_scale", retarget_run(&src, &nan_tgt, &untouched, &st) == RET_ERR_SCALE);

    /* planted-foot slip is bounded */
    Clip slipping = walk_clip(0.1f);
    check("slide_error", retarget_run(&slipping, &cfg, &untouched, &st) == RET_ERR_SLIDE);
    RetargetCfg loose = cfg;
    loose.slide_limit = 10.0f;
    Clip slipped;
    memset(&slipped, 0, sizeof(slipped));
    rc = retarget_run(&slipping, &loose, &slipped, &st);
    check("slide_measured", rc == RET_OK && st.slide > 0.5f);

    /* disabled retargeting is a byte-for-byte copy of the source */
    RetargetCfg off = cfg;
    off.enabled = 0;
    Clip passthrough;
    memset(&passthrough, 0, sizeof(passthrough));
    rc = retarget_run(&src, &off, &passthrough, &st);
    check("disabled_ok", rc == RET_OK);
    check("disabled_identity", memcmp(&passthrough, &src, sizeof(Clip)) == 0);

    /* disabled playback ignores retarget settings entirely */
    RetargetCfg off_bad = off;
    off_bad.source_leg_length = 0.0f;
    off_bad.target_contact_phase = 0.9f;
    memset(&passthrough, 0, sizeof(passthrough));
    rc = retarget_run(&src, &off_bad, &passthrough, &st);
    check("disabled_ignores_settings", rc == RET_OK &&
          memcmp(&passthrough, &src, sizeof(Clip)) == 0);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
