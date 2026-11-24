
/* curve: keyframe sampling with cubic hermite tangents and clamping. */
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "../../shared/pose_format.h"

#define MAX_KEYS 16

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

/* ===================== shared pose format (new) ==========================
 *
 * One "key" record per keyframe, carrying time, value and both Hermite
 * tangents:
 *
 *   key t=0 v=0 tin=0 tout=1
 *
 * Reading the file back rebuilds the Curve through curve_add(), so the
 * roundtrip checks the internal key array (and rejects the same reordered
 * input the in-memory API rejects).
 */

static Curve curve_pose_build(void) {
    Curve c;
    curve_init(&c);
    curve_add(&c, 0.0f, 0.0f, 0.0f, 1.0f);
    curve_add(&c, 0.5f, 0.75f, 0.4f, 0.2f);
    curve_add(&c, 1.0f, 1.0f, 1.0f, 0.0f);
    return c;
}

static int curve_export(const Curve* c, const char* file) {
    FILE* f = fopen(file, "w");
    if (!f) return -1;
    pf_write_header(f, "curve", c->count);
    for (int i = 0; i < c->count; i++) {
        pf_record_begin(f, "key");
        pf_field_float(f, "t", c->keys[i].t);
        pf_field_float(f, "v", c->keys[i].v);
        pf_field_float(f, "tin", c->keys[i].tan_in);
        pf_field_float(f, "tout", c->keys[i].tan_out);
        pf_record_end(f);
    }
    fclose(f);
    return 0;
}

static int curve_import(const char* file, Curve* c) {
    FILE* f = fopen(file, "r");
    if (!f) return -1;

    int n = 0;
    if (pf_read_header(f, "curve", &n) != 0 || n > MAX_KEYS) {
        fclose(f);
        return -1;
    }

    curve_init(c);
    PfRecord r;
    for (int k = 0; k < n; k++) {
        if (pf_read_record(f, &r) != 0 || strcmp(r.name, "key") != 0) {
            fclose(f);
            return -1;
        }
        float t, v, tin, tout;
        if (pf_get_float(&r, "t", &t) || pf_get_float(&r, "v", &v) ||
            pf_get_float(&r, "tin", &tin) || pf_get_float(&r, "tout", &tout) ||
            curve_add(c, t, v, tin, tout) < 0) {
            fclose(f);
            return -1;
        }
    }
    fclose(f);
    return 0;
}

/* ---- legacy curve format (pre-unification): one headerless CSV line -----
 *
 *   t,v,tan_in,tan_out
 */

static int curve_export_legacy(const Curve* c, const char* file) {
    FILE* f = fopen(file, "w");
    if (!f) return -1;
    for (int i = 0; i < c->count; i++)
        fprintf(f, "%.9g,%.9g,%.9g,%.9g\n",
                (double)c->keys[i].t, (double)c->keys[i].v,
                (double)c->keys[i].tan_in, (double)c->keys[i].tan_out);
    fclose(f);
    return 0;
}

static int curve_import_legacy(const char* file, Curve* c) {
    FILE* f = fopen(file, "r");
    if (!f) return -1;

    curve_init(c);
    char line[256];
    int ok = 1;
    while (ok && fgets(line, sizeof line, f)) {
        char* p = line;
        while (*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n') p++;
        if (*p == '\0' || *p == '#') continue;

        float t, v, tin, tout;
        if (sscanf(p, "%f,%f,%f,%f", &t, &v, &tin, &tout) != 4 ||
            curve_add(c, t, v, tin, tout) < 0) {
            ok = 0;
        }
    }
    fclose(f);
    return (ok && c->count > 0) ? 0 : -1;
}

static int curve_same(const Curve* a, const Curve* b) {
    if (a->count != b->count) return 0;
    for (int i = 0; i < a->count; i++) {
        if (a->keys[i].t != b->keys[i].t || a->keys[i].v != b->keys[i].v ||
            a->keys[i].tan_in != b->keys[i].tan_in ||
            a->keys[i].tan_out != b->keys[i].tan_out) return 0;
    }
    return 1;
}

/* ------------------------------------------------------------ CLI modes */

static int curve_cli(int argc, char** argv) {
    const char* mode = argv[1];
    Curve src = curve_pose_build(), back;

    if (strcmp(mode, "--pose-export") == 0 && argc == 3)
        return curve_export(&src, argv[2]) == 0 ? 0 : 1;
    if (strcmp(mode, "--pose-import") == 0 && argc == 3) {
        if (curve_import(argv[2], &back) != 0) { printf("imported=0\n"); return 1; }
        printf("imported=1\ncount=%d\n", back.count);
        return 0;
    }
    if (strcmp(mode, "--pose-roundtrip") == 0 && argc == 3) {
        int ok = curve_export(&src, argv[2]) == 0 &&
                 curve_import(argv[2], &back) == 0 &&
                 curve_same(&src, &back);
        printf("roundtrip=%d\n", ok ? 1 : 0);
        return ok ? 0 : 1;
    }
    if (strcmp(mode, "--legacy-export") == 0 && argc == 3)
        return curve_export_legacy(&src, argv[2]) == 0 ? 0 : 1;
    if (strcmp(mode, "--legacy-check") == 0 && argc == 3) {
        int ok = curve_import_legacy(argv[2], &back) == 0 && curve_same(&src, &back);
        printf("legacy=%d\n", ok ? 1 : 0);
        return ok ? 0 : 1;
    }

    fprintf(stderr, "usage: %s [--pose-export|--pose-import|--pose-roundtrip|"
                    "--legacy-export|--legacy-check] <file>\n", argv[0]);
    return 2;
}

int main(int argc, char** argv) {
    if (argc > 1) return curve_cli(argc, argv);

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
