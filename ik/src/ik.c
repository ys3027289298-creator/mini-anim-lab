
/* ik: analytic two-bone IK with pole vector and joint angle limits. */
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "../../shared/pose_format.h"

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

/* ===================== shared pose format (new) ==========================
 *
 * One "ik" record holds the solve inputs (root, end-effector target, pole,
 * bone lengths) and the result (joint position and both joint angles):
 *
 *   ik root_x=0 root_y=0 target_x=1.4 target_y=0.4 pole_x=0 pole_y=1
 *      l1=1 l2=1 ok=1 joint_x=.. joint_y=.. a1=.. a2=..
 *
 * Reading it back restores the IkResult bit-for-bit, so the roundtrip checks
 * the solver output rather than a re-solve that could drift.
 */

typedef struct {
    V2 root, target, pole;
    float l1, l2;
    IkResult result;
} IkPose;

static IkPose ik_pose_build(void) {
    IkPose p;
    p.root = (V2){ 0.0f, 0.0f };
    p.target = (V2){ 1.4f, 0.4f };
    p.pole = (V2){ 0.0f, 1.0f };
    p.l1 = 1.0f;
    p.l2 = 1.0f;
    p.result = solve_two_bone(p.root, p.target, p.l1, p.l2, p.pole);
    return p;
}

static int ik_export(const IkPose* p, const char* file) {
    FILE* f = fopen(file, "w");
    if (!f) return -1;
    pf_write_header(f, "ik", 1);
    pf_record_begin(f, "ik");
    pf_field_float(f, "root_x", p->root.x);
    pf_field_float(f, "root_y", p->root.y);
    pf_field_float(f, "target_x", p->target.x);
    pf_field_float(f, "target_y", p->target.y);
    pf_field_float(f, "pole_x", p->pole.x);
    pf_field_float(f, "pole_y", p->pole.y);
    pf_field_float(f, "l1", p->l1);
    pf_field_float(f, "l2", p->l2);
    pf_field_int(f, "ok", p->result.ok);
    pf_field_float(f, "joint_x", p->result.joint.x);
    pf_field_float(f, "joint_y", p->result.joint.y);
    pf_field_float(f, "a1", p->result.a1);
    pf_field_float(f, "a2", p->result.a2);
    pf_record_end(f);
    fclose(f);
    return 0;
}

static int ik_import(const char* file, IkPose* p) {
    FILE* f = fopen(file, "r");
    if (!f) return -1;

    int n = 0;
    if (pf_read_header(f, "ik", &n) != 0 || n != 1) {
        fclose(f);
        return -1;
    }

    PfRecord r;
    if (pf_read_record(f, &r) != 0 || strcmp(r.name, "ik") != 0) {
        fclose(f);
        return -1;
    }

    int ok = 0;
    int bad = pf_get_float(&r, "root_x", &p->root.x) ||
              pf_get_float(&r, "root_y", &p->root.y) ||
              pf_get_float(&r, "target_x", &p->target.x) ||
              pf_get_float(&r, "target_y", &p->target.y) ||
              pf_get_float(&r, "pole_x", &p->pole.x) ||
              pf_get_float(&r, "pole_y", &p->pole.y) ||
              pf_get_float(&r, "l1", &p->l1) ||
              pf_get_float(&r, "l2", &p->l2) ||
              pf_get_int(&r, "ok", &ok) ||
              pf_get_float(&r, "joint_x", &p->result.joint.x) ||
              pf_get_float(&r, "joint_y", &p->result.joint.y) ||
              pf_get_float(&r, "a1", &p->result.a1) ||
              pf_get_float(&r, "a2", &p->result.a2);
    p->result.ok = ok;
    fclose(f);
    return bad ? -1 : 0;
}

/* ---- legacy ik format (pre-unification): one headerless CSV line --------
 *
 *   root_x,root_y,target_x,target_y,pole_x,pole_y,l1,l2,ok,joint_x,joint_y,a1,a2
 */

static int ik_export_legacy(const IkPose* p, const char* file) {
    FILE* f = fopen(file, "w");
    if (!f) return -1;
    fprintf(f, "%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%d,%.9g,%.9g,%.9g,%.9g\n",
            (double)p->root.x, (double)p->root.y,
            (double)p->target.x, (double)p->target.y,
            (double)p->pole.x, (double)p->pole.y,
            (double)p->l1, (double)p->l2, p->result.ok,
            (double)p->result.joint.x, (double)p->result.joint.y,
            (double)p->result.a1, (double)p->result.a2);
    fclose(f);
    return 0;
}

static int ik_import_legacy(const char* file, IkPose* p) {
    FILE* f = fopen(file, "r");
    if (!f) return -1;

    char line[512];
    int ok = 0;
    int got = 0;
    while (fgets(line, sizeof line, f)) {
        char* q = line;
        while (*q == ' ' || *q == '\t' || *q == '\r' || *q == '\n') q++;
        if (*q == '\0' || *q == '#') continue;
        got = sscanf(q, "%f,%f,%f,%f,%f,%f,%f,%f,%d,%f,%f,%f,%f",
                     &p->root.x, &p->root.y, &p->target.x, &p->target.y,
                     &p->pole.x, &p->pole.y, &p->l1, &p->l2, &ok,
                     &p->result.joint.x, &p->result.joint.y,
                     &p->result.a1, &p->result.a2);
        break;
    }
    p->result.ok = ok;
    fclose(f);
    return got == 13 ? 0 : -1;
}

static int ik_pose_same(const IkPose* a, const IkPose* b) {
    return a->root.x == b->root.x && a->root.y == b->root.y &&
           a->target.x == b->target.x && a->target.y == b->target.y &&
           a->pole.x == b->pole.x && a->pole.y == b->pole.y &&
           a->l1 == b->l1 && a->l2 == b->l2 &&
           a->result.ok == b->result.ok &&
           a->result.joint.x == b->result.joint.x &&
           a->result.joint.y == b->result.joint.y &&
           a->result.a1 == b->result.a1 && a->result.a2 == b->result.a2;
}

/* ------------------------------------------------------------ CLI modes */

static int ik_cli(int argc, char** argv) {
    const char* mode = argv[1];
    IkPose src = ik_pose_build(), back;

    if (strcmp(mode, "--pose-export") == 0 && argc == 3)
        return ik_export(&src, argv[2]) == 0 ? 0 : 1;
    if (strcmp(mode, "--pose-import") == 0 && argc == 3) {
        if (ik_import(argv[2], &back) != 0) { printf("imported=0\n"); return 1; }
        printf("imported=1\nok=%d\n", back.result.ok);
        return 0;
    }
    if (strcmp(mode, "--pose-roundtrip") == 0 && argc == 3) {
        int ok = ik_export(&src, argv[2]) == 0 &&
                 ik_import(argv[2], &back) == 0 &&
                 ik_pose_same(&src, &back);
        printf("roundtrip=%d\n", ok ? 1 : 0);
        return ok ? 0 : 1;
    }
    if (strcmp(mode, "--legacy-export") == 0 && argc == 3)
        return ik_export_legacy(&src, argv[2]) == 0 ? 0 : 1;
    if (strcmp(mode, "--legacy-check") == 0 && argc == 3) {
        int ok = ik_import_legacy(argv[2], &back) == 0 && ik_pose_same(&src, &back);
        printf("legacy=%d\n", ok ? 1 : 0);
        return ok ? 0 : 1;
    }

    fprintf(stderr, "usage: %s [--pose-export|--pose-import|--pose-roundtrip|"
                    "--legacy-export|--legacy-check] <file>\n", argv[0]);
    return 2;
}

int main(int argc, char** argv) {
    if (argc > 1) return ik_cli(argc, argv);

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
