
/* curve: keyframe sampling with cubic hermite tangents and clamping. */
#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "../../shared/pose_export.h"

#define MAX_KEYS 16
#define POSE_SAMPLES 21

typedef struct { float t; float v; float tan_in; float tan_out; } Key;

typedef struct {
    Key keys[MAX_KEYS];
    int count;
} Curve;

static void curve_init(Curve* c) { memset(c, 0, sizeof(*c)); }

static int curve_add(Curve* c, float t, float v, float tin, float tout) {
    if (c->count >= MAX_KEYS) return -1;
    if (c->count > 0 && t <= c->keys[c->count - 1].t) return -1;  /* strictly increasing */
    Key* k = &c->keys[c->count++];
    k->t = t; k->v = v; k->tan_in = tin; k->tan_out = tout;
    return c->count - 1;
}

static float curve_eval(const Curve* c, float t) {
    if (c->count == 0) return 0.0f;
    if (t <= c->keys[0].t) return c->keys[0].v;
    if (t >= c->keys[c->count - 1].t) return c->keys[c->count - 1].v;
    for (int i = 1; i < c->count; i++) {
        const Key* a = &c->keys[i - 1];
        const Key* b = &c->keys[i];
        if (t <= b->t) {
            float dt = b->t - a->t;
            float u = (t - a->t) / dt;
            float u2 = u * u, u3 = u2 * u;
            float h00 = 2 * u3 - 3 * u2 + 1;
            float h10 = u3 - 2 * u2 + u;
            float h01 = -2 * u3 + 3 * u2;
            float h11 = u3 - u2;
            return h00 * a->v + h10 * dt * a->tan_out + h01 * b->v + h11 * dt * b->tan_in;
        }
    }
    return c->keys[c->count - 1].v;
}

/* derivative for continuity checks */
static float curve_deriv(const Curve* c, float t) {
    float h = 1e-3f;
    float a = curve_eval(c, t - h);
    float b = curve_eval(c, t + h);
    return (b - a) / (2 * h);
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

/* ---------- unified pose export (new format, see shared/pose_export.h) ----- */

/* canonical pose: hermite curve 0->1 with tangents 1->0, sampled on [0,1] */
static void curve_pose_build(Curve* c) {
    curve_init(c);
    curve_add(c, 0.0f, 0.0f, 0.0f, 1.0f);
    curve_add(c, 1.0f, 1.0f, 1.0f, 0.0f);
}

static int curve_pose_export(const char* outdir, int argc, char** argv) {
    if (pe_mkdir_p(outdir) != 0) return 2;

    float vals[POSE_SAMPLES];
    /* timing tier 1: pose numeric computation, 1000 repeats */
    Curve c;
    long t0 = pe_now_us();
    for (int rep = 0; rep < 1000; rep++) {
        curve_pose_build(&c);
        for (int i = 0; i < POSE_SAMPLES; i++)
            vals[i] = curve_eval(&c, (float)i / (POSE_SAMPLES - 1));
        if (vals[0] + vals[POSE_SAMPLES - 1] > 1e30f) return 2;
    }
    long t1 = pe_now_us();
    long compute_us = t1 - t0;

    /* timing tier 2: image rendering, single pass */
    PeCanvas cv;
    pe_canvas_init(&cv, -0.1f, -0.1f, 1.1f, 1.1f);
    long t2 = pe_now_us();
    for (int i = 0; i < POSE_SAMPLES - 1; i++) {
        float ta = (float)i / (POSE_SAMPLES - 1);
        float tb = (float)(i + 1) / (POSE_SAMPLES - 1);
        pe_line(&cv, ta, vals[i], tb, vals[i + 1], 40, 90, 180);
    }
    for (int i = 0; i < c.count; i++)
        pe_circle(&cv, c.keys[i].t, c.keys[i].v, 0.03f, 200, 60, 60);
    long t3 = pe_now_us();
    long render_us = t3 - t2;

    char ppm_path[600], json_path[600];
    snprintf(ppm_path, sizeof(ppm_path), "%s/pose.ppm", outdir);
    snprintf(json_path, sizeof(json_path), "%s/pose.json", outdir);
    if (pe_canvas_save(&cv, ppm_path) != 0) return 2;

    char hash[65];
    if (pe_sha256_file(ppm_path, hash) != 0) return 2;

    FILE* f = fopen(json_path, "wb");
    if (!f) return 2;
    fprintf(f, "{\n");
    fprintf(f, "  \"module\": \"curve\",\n");
    fprintf(f, "  \"format_version\": %d,\n", POSE_FORMAT_VERSION);
    fprintf(f, "  \"pose\": {\n");
    fprintf(f, "    \"key_count\": %d,\n", c.count);
    fprintf(f, "    \"sample_count\": %d,\n", POSE_SAMPLES);
    fprintf(f, "    \"samples\": [\n");
    for (int i = 0; i < POSE_SAMPLES; i++) {
        float t = (float)i / (POSE_SAMPLES - 1);
        fprintf(f, "      {\"t\": %.6f, \"v\": %.6f}%s\n", t, vals[i],
                i == POSE_SAMPLES - 1 ? "" : ",");
    }
    fprintf(f, "    ]\n");
    fprintf(f, "  },\n");
    fprintf(f, "  \"image\": {\n");
    fprintf(f, "    \"file\": \"pose.ppm\",\n");
    fprintf(f, "    \"format\": \"ppm-p6\",\n");
    fprintf(f, "    \"width\": %d,\n", POSE_CANVAS_W);
    fprintf(f, "    \"height\": %d,\n", POSE_CANVAS_H);
    fprintf(f, "    \"sha256\": \"%s\"\n", hash);
    fprintf(f, "  },\n");
    fprintf(f, "  \"timing\": {\n");
    fprintf(f, "    \"compute_pose_us\": %ld,\n", compute_us);
    fprintf(f, "    \"compute_repeats\": 1000,\n");
    fprintf(f, "    \"render_image_us\": %ld\n", render_us);
    fprintf(f, "  },\n");
    fprintf(f, "  \"reproduce\": ");
    fputc('[', f);
    for (int i = 0; i < argc; i++) {
        pe_write_json_string(f, argv[i]);
        fputs(i == argc - 1 ? "" : ", ", f);
    }
    fputs("],\n", f);
    fprintf(f, "  \"reproduce_hint\": \"run the command above from the repository root\"\n");
    fprintf(f, "}\n");
    fclose(f);
    return 0;
}

int main(int argc, char** argv) {
    if (strcmp(argc > 1 ? argv[1] : "", "--pose-export") == 0) {
        if (argc != 3) {
            fprintf(stderr, "usage: %s --pose-export <outdir>\n", argv[0]);
            return 2;
        }
        return curve_pose_export(argv[2], argc, argv);
    }

    Curve c;
    curve_init(&c);
    check("add_first", curve_add(&c, 0.0f, 0.0f, 0.0f, 1.0f) == 0);
    check("add_second", curve_add(&c, 1.0f, 1.0f, 1.0f, 0.0f) == 1);
    check("clamp_left", fabsf(curve_eval(&c, -1.0f)) < 1e-6f);
    check("clamp_right", fabsf(curve_eval(&c, 2.0f) - 1.0f) < 1e-6f);
    check("mid_monotonic", curve_eval(&c, 0.5f) > 0.0f && curve_eval(&c, 0.5f) < 1.0f);
    check("endpoints_exact", fabsf(curve_eval(&c, 0.0f)) < 1e-6f && fabsf(curve_eval(&c, 1.0f) - 1.0f) < 1e-6f);

    Curve flat;
    curve_init(&flat);
    curve_add(&flat, 0.0f, 0.5f, 0.0f, 0.0f);
    curve_add(&flat, 1.0f, 0.5f, 0.0f, 0.0f);
    check("flat_curve", fabsf(curve_eval(&flat, 0.5f) - 0.5f) < 1e-5f);
    check("flat_deriv_zero", fabsf(curve_deriv(&flat, 0.5f)) < 1e-3f);

    check("reject_unsorted", curve_add(&flat, 0.5f, 1.0f, 0, 0) == -1);

    Curve full;
    curve_init(&full);
    int rc = 0;
    for (int i = 0; i <= MAX_KEYS; i++) rc = curve_add(&full, (float)i, (float)i, 0, 0);
    check("capacity", rc == -1 && full.count == MAX_KEYS);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
