/* blend: transition buffer for collapse frames.
   Low curve stiffness lets a pose switch collapse and infect neighboring
   frames. This module:
     1. flags abnormal frames with a per-frame physical step limit
     2. derives isolation spacing S = ceil(E / v_max) from node offset E
     3. rebuilds the infected zone as a transition buffer anchored on the
        original frames outside the zone
     4. emits a position offset table marking abnormal frames
     5. exports collapsed poses next to hand-keyed reference poses */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define MAX_FRAMES 64
#define MAX_NODES 8

typedef struct { float x, y; } V2;

typedef struct {
    V2 pos[MAX_FRAMES][MAX_NODES];
    int n_frames;
    int n_nodes;
} Track;

typedef struct {
    float offset[MAX_FRAMES];   /* max node offset vs reference, per frame */
    int abnormal[MAX_FRAMES];
    int n;
} OffsetTable;

static float vdist(V2 a, V2 b) {
    float dx = a.x - b.x, dy = a.y - b.y;
    return sqrtf(dx * dx + dy * dy);
}

/* max node displacement between frame f-1 and frame f */
static float frame_step(const Track* t, int f) {
    float worst = 0.0f;
    for (int n = 0; n < t->n_nodes; n++) {
        float d = vdist(t->pos[f][n], t->pos[f - 1][n]);
        if (d > worst) worst = d;
    }
    return worst;
}

static float track_max_step(const Track* t) {
    float worst = 0.0f;
    for (int f = 1; f < t->n_frames; f++) {
        float s = frame_step(t, f);
        if (s > worst) worst = s;
    }
    return worst;
}

/* isolation spacing vs node offset degree:
   absorbing offset E linearly over S frames per side gives a per-frame
   correction of E/S; the physical limit requires E/S <= v_max, hence
       S = ceil(E / v_max)
   and the residual offset k frames away from the collapse center is
       r(k) = E * max(0, 1 - k/S)                                        */
static int isolation_spacing(float E, float v_max) {
    if (E <= 0.0f || v_max <= 0.0f) return 0;
    return (int)ceilf(E / v_max - 1e-6f);
}

static float residual_offset(float E, int S, int k) {
    if (k < 0) k = -k;
    if (S <= 0 || k >= S) return 0.0f;
    return E * (1.0f - (float)k / (float)S);
}

/* position offset table: per-frame max node offset vs reference,
   frames above threshold are marked abnormal */
static void offset_table_build(OffsetTable* ot, const Track* cur, const Track* ref, float threshold) {
    memset(ot, 0, sizeof(*ot));
    ot->n = cur->n_frames;
    for (int f = 0; f < cur->n_frames; f++) {
        float worst = 0.0f;
        for (int n = 0; n < cur->n_nodes; n++) {
            float d = vdist(cur->pos[f][n], ref->pos[f][n]);
            if (d > worst) worst = d;
        }
        ot->offset[f] = worst;
        ot->abnormal[f] = worst > threshold;
    }
}

static float smoothstep(float u) { return u * u * (3.0f - 2.0f * u); }

/* rebuild frames [lo..hi] as a smooth transition between the original
   anchor poses at lo-1 and hi+1 (frames the collapse never reached) */
static void build_transition_buffer(Track* t, int lo, int hi) {
    if (lo < 1) lo = 1;
    if (hi > t->n_frames - 2) hi = t->n_frames - 2;
    float span = (float)((hi + 1) - (lo - 1));
    for (int f = lo; f <= hi; f++) {
        float u = (float)(f - (lo - 1)) / span;
        float w = smoothstep(u);
        for (int n = 0; n < t->n_nodes; n++) {
            V2 a = t->pos[lo - 1][n];
            V2 b = t->pos[hi + 1][n];
            t->pos[f][n].x = a.x + (b.x - a.x) * w;
            t->pos[f][n].y = a.y + (b.y - a.y) * w;
        }
    }
}

/* cross-frame physical constraint: clamp per-frame node displacement to
   v_max, walking forward from frame 0 so a clamp never re-violates
   earlier frames */
static void apply_physical_limits(Track* t, float v_max) {
    for (int f = 1; f < t->n_frames; f++) {
        for (int n = 0; n < t->n_nodes; n++) {
            V2* p = &t->pos[f][n];
            float dx = p->x - t->pos[f - 1][n].x;
            float dy = p->y - t->pos[f - 1][n].y;
            float d = sqrtf(dx * dx + dy * dy);
            if (d > v_max && d > 0.0f) {
                float s = v_max / d;
                p->x = t->pos[f - 1][n].x + dx * s;
                p->y = t->pos[f - 1][n].y + dy * s;
            }
        }
    }
}

/* export collapsed poses next to hand-keyed reference poses, one block
   per abnormal frame, with per-node deviation */
static int export_collapsed(const char* path, const Track* cur, const Track* ref, const OffsetTable* ot) {
    FILE* fp = fopen(path, "w");
    if (!fp) return -1;
    fprintf(fp, "# collapsed pose vs hand-keyed reference\n");
    for (int f = 0; f < ot->n; f++) {
        if (!ot->abnormal[f]) continue;
        fprintf(fp, "frame %d offset %.4f\n", f, ot->offset[f]);
        for (int n = 0; n < cur->n_nodes; n++) {
            float d = vdist(cur->pos[f][n], ref->pos[f][n]);
            fprintf(fp, "node %d cur %.4f %.4f ref %.4f %.4f dev %.4f\n",
                    n, cur->pos[f][n].x, cur->pos[f][n].y,
                    ref->pos[f][n].x, ref->pos[f][n].y, d);
        }
    }
    fclose(fp);
    return 0;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

#define NF 41
#define NN 3
#define SWITCH_F 20
#define INFECT_R 3
#define COLLAPSE_E 0.8f
#define V_MAX 0.1f
#define OFFSET_THRESH 0.05f

/* hand-keyed reference: root fixed, two links rotating smoothly */
static void build_reference(Track* t) {
    t->n_frames = NF;
    t->n_nodes = NN;
    for (int f = 0; f < NF; f++) {
        float a = 0.02f * f;
        float b = 0.03f * f + 0.5f;
        t->pos[f][0].x = 0.0f; t->pos[f][0].y = 0.0f;
        t->pos[f][1].x = cosf(a); t->pos[f][1].y = sinf(a);
        t->pos[f][2].x = cosf(a) + cosf(b); t->pos[f][2].y = sinf(a) + sinf(b);
    }
}

/* low-stiffness collapse: nodes pulled toward the root, strongest at the
   switch frame, bleeding linearly into INFECT_R neighbors on each side */
static void inject_collapse(Track* t, const Track* ref) {
    *t = *ref;
    for (int f = SWITCH_F - INFECT_R; f <= SWITCH_F + INFECT_R; f++) {
        float w = 1.0f - (float)abs(f - SWITCH_F) / (INFECT_R + 1.0f);
        for (int n = 0; n < NN; n++) {
            float dx = -ref->pos[f][n].x, dy = -ref->pos[f][n].y;
            float l = sqrtf(dx * dx + dy * dy);
            if (l < 1e-6f) continue;    /* root stays */
            t->pos[f][n].x += dx / l * COLLAPSE_E * w;
            t->pos[f][n].y += dy / l * COLLAPSE_E * w;
        }
    }
}

int main(void) {
    Track ref, cur, collapsed;
    build_reference(&ref);
    inject_collapse(&cur, &ref);
    collapsed = cur;

    /* reference itself respects the physical limit */
    check("reference_within_limits", track_max_step(&ref) <= V_MAX);

    /* collapse breaks the physical limit around the switch frame */
    check("collapse_breaks_limit",
          frame_step(&cur, SWITCH_F) > V_MAX || frame_step(&cur, SWITCH_F + 1) > V_MAX);

    /* offset table marks exactly the infected frames as abnormal */
    OffsetTable ot;
    offset_table_build(&ot, &cur, &ref, OFFSET_THRESH);
    int exact = 1;
    for (int f = 0; f < NF; f++) {
        int expect = abs(f - SWITCH_F) <= INFECT_R;
        if (ot.abnormal[f] != expect) exact = 0;
    }
    check("offset_table_marks_infected", exact);
    check("offset_peak_at_switch", fabsf(ot.offset[SWITCH_F] - COLLAPSE_E) < 1e-4f);

    /* spacing derivation: S = ceil(E / v_max) = ceil(0.8 / 0.1) = 8 */
    float E = ot.offset[SWITCH_F];
    int S = isolation_spacing(E, V_MAX);
    check("spacing_formula", S == 8);
    check("spacing_covers_infection", S > INFECT_R);
    check("correction_within_limit", E / S <= V_MAX + 1e-6f);
    check("residual_decays",
          fabsf(residual_offset(E, S, 0) - E) < 1e-6f &&
          fabsf(residual_offset(E, S, S)) < 1e-6f &&
          fabsf(residual_offset(E, S, 4) - E * 0.5f) < 1e-4f);

    /* buffer anchors sit on clean original frames outside the collapse */
    check("anchors_clean",
          !ot.abnormal[SWITCH_F - S - 1] && !ot.abnormal[SWITCH_F + S + 1]);

    /* rebuild the infected zone from original anchors, then clamp */
    build_transition_buffer(&cur, SWITCH_F - S, SWITCH_F + S);
    apply_physical_limits(&cur, V_MAX);
    check("limits_hold_after_fix", track_max_step(&cur) <= V_MAX + 1e-5f);
    check("anchors_untouched",
          vdist(cur.pos[SWITCH_F - S - 1][2], ref.pos[SWITCH_F - S - 1][2]) < 1e-6f &&
          vdist(cur.pos[SWITCH_F + S + 1][2], ref.pos[SWITCH_F + S + 1][2]) < 1e-6f);

    OffsetTable fixed;
    offset_table_build(&fixed, &cur, &ref, OFFSET_THRESH);
    float worst = 0.0f;
    for (int f = 0; f < NF; f++)
        if (fixed.offset[f] > worst) worst = fixed.offset[f];
    printf("worst_deviation_after_fix=%.4f\n", worst);
    check("buffer_tracks_reference", worst < 0.25f);

    /* export collapsed poses next to the hand-keyed reference frames */
    const char* out = "blend_export.txt";
    check("export_written", export_collapsed(out, &collapsed, &ref, &ot) == 0);
    int flagged = 0;
    for (int f = 0; f < ot.n; f++) {
        if (!ot.abnormal[f]) continue;
        flagged++;
        printf("cmp frame %d offset %.4f\n", f, ot.offset[f]);
    }
    FILE* fp = fopen(out, "r");
    int lines = 0, ch;
    while ((ch = fgetc(fp)) != EOF) if (ch == '\n') lines++;
    fclose(fp);
    check("export_line_count", lines == 1 + flagged * (1 + NN));

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
