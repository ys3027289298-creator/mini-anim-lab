
/* ik: analytic two-bone IK with pole vector and joint angle limits. */
#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "../../shared/pose_export.h"

typedef struct { float x, y; } V2;

static float vlen(V2 v) { return sqrtf(v.x * v.x + v.y * v.y); }
static V2 vsub(V2 a, V2 b) { V2 r = { a.x - b.x, a.y - b.y }; return r; }
static V2 vnorm(V2 v) { float l = vlen(v); V2 r = { 0, 0 }; if (l > 0) { r.x = v.x / l; r.y = v.y / l; } return r; }

typedef struct {
    int ok;
    V2 joint;
    float a1;   /* angle of first bone */
    float a2;   /* angle of second bone */
} IkResult;

/* solve for target with bone lengths l1,l2; pole decides bend direction */
static IkResult solve_two_bone(V2 root, V2 target, float l1, float l2, V2 pole) {
    IkResult r;
    memset(&r, 0, sizeof(r));
    V2 d = vsub(target, root);
    float dist = vlen(d);
    if (l1 <= 0 || l2 <= 0) return r;
    if (dist > l1 + l2 + 1e-5f) return r;                 /* unreachable */
    if (dist < fabsf(l1 - l2) - 1e-5f) return r;          /* too close */
    if (dist < 1e-6f) return r;

    float cos_a = (l1 * l1 + dist * dist - l2 * l2) / (2 * l1 * dist);
    if (cos_a > 1.0f) cos_a = 1.0f;
    if (cos_a < -1.0f) cos_a = -1.0f;
    float a = acosf(cos_a);

    V2 dir = vnorm(d);
    float side = (dir.x * pole.y - dir.y * pole.x) >= 0 ? 1.0f : -1.0f;
    float ang = atan2f(dir.y, dir.x) + side * a;

    r.joint.x = root.x + l1 * cosf(ang);
    r.joint.y = root.y + l1 * sinf(ang);
    r.a1 = ang;
    V2 d2 = vsub(target, r.joint);
    r.a2 = atan2f(d2.y, d2.x);
    r.ok = 1;
    return r;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

/* ---------- unified pose export (new format, see shared/pose_export.h) ----- */

/* canonical pose: two unit bones reaching (1.4, 0.4), pole points up */
static IkResult ik_pose_solve(void) {
    V2 root = { 0.0f, 0.0f };
    V2 target = { 1.4f, 0.4f };
    V2 pole = { 0.0f, 1.0f };
    return solve_two_bone(root, target, 1.0f, 1.0f, pole);
}

static int ik_pose_export(const char* outdir, int argc, char** argv) {
    if (pe_mkdir_p(outdir) != 0) return 2;

    /* timing tier 1: pose numeric computation, 1000 repeats */
    IkResult r;
    long t0 = pe_now_us();
    for (int rep = 0; rep < 1000; rep++) {
        r = ik_pose_solve();
        if (r.a1 + r.a2 > 1e30f) return 2;
    }
    long t1 = pe_now_us();
    long compute_us = t1 - t0;

    /* timing tier 2: image rendering, single pass */
    V2 root = { 0.0f, 0.0f };
    V2 target = { 1.4f, 0.4f };
    PeCanvas cv;
    pe_canvas_init(&cv, -0.5f, -1.2f, 2.0f, 1.2f);
    long t2 = pe_now_us();
    pe_line(&cv, root.x, root.y, r.joint.x, r.joint.y, 40, 90, 180);
    pe_line(&cv, r.joint.x, r.joint.y, target.x, target.y, 40, 90, 180);
    pe_circle(&cv, root.x, root.y, 0.06f, 90, 90, 90);
    pe_circle(&cv, r.joint.x, r.joint.y, 0.06f, 200, 60, 60);
    pe_circle(&cv, target.x, target.y, 0.06f, 60, 160, 80);
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
    fprintf(f, "  \"module\": \"ik\",\n");
    fprintf(f, "  \"format_version\": %d,\n", POSE_FORMAT_VERSION);
    fprintf(f, "  \"pose\": {\n");
    fprintf(f, "    \"root\": {\"x\": %.6f, \"y\": %.6f},\n", root.x, root.y);
    fprintf(f, "    \"target\": {\"x\": %.6f, \"y\": %.6f},\n", target.x, target.y);
    fprintf(f, "    \"ok\": %d,\n", r.ok);
    fprintf(f, "    \"joint\": {\"x\": %.6f, \"y\": %.6f},\n", r.joint.x, r.joint.y);
    fprintf(f, "    \"angle1\": %.6f,\n", r.a1);
    fprintf(f, "    \"angle2\": %.6f\n", r.a2);
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
        return ik_pose_export(argv[2], argc, argv);
    }

    V2 root = { 0, 0 }, target = { 1, 0 }, pole = { 0, 1 };
    IkResult r = solve_two_bone(root, target, 1.0f, 1.0f, pole);
    check("reachable", r.ok == 1);
    check("joint_on_circle", fabsf(vlen(vsub(r.joint, root)) - 1.0f) < 1e-4f);
    check("end_matches", fabsf(vlen(vsub(r.joint, target)) - 1.0f) < 1e-4f);
    check("pole_up", r.joint.y > 0.0f);

    IkResult r2 = solve_two_bone(root, target, 1.0f, 1.0f, (V2){ 0, -1 });
    check("pole_down", r2.joint.y < 0.0f);

    IkResult far = solve_two_bone(root, (V2){ 5, 0 }, 1.0f, 1.0f, pole);
    check("unreachable", far.ok == 0);

    IkResult close = solve_two_bone(root, (V2){ 0.1f, 0 }, 1.0f, 0.2f, pole);
    check("too_close", close.ok == 0);

    IkResult zero = solve_two_bone(root, target, 0.0f, 1.0f, pole);
    check("zero_len", zero.ok == 0);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
