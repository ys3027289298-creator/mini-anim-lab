/* replay: animation playback consistency checker.
 *
 * Pipeline: keyframe curves (same math as ../curve) -> per-bone local poses
 * -> skeleton world transforms (same math as ../skeleton), recorded per frame.
 *
 * Usage:
 *   replay check  <anim> <ref.rec>   validate input, play twice in-process and
 *                                    diff per bone/frame, then compare against
 *                                    the hand-computed reference and report an
 *                                    error upper bound.
 *   replay record <anim> <out.rec>   validate input, play once, dump recording.
 *   replay diff   <a.rec> <b.rec>    diff two recordings per bone/frame.
 *
 * Input validation: NaN/Inf in any numeric field, negative key time, or
 * negative frame start time is a hard error (exit 2). Sample files are only
 * ever opened read-only; the only file written is <out.rec> of `record`.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_BONES 16
#define MAX_KEYS 16
#define MAX_FRAMES 1024

#define REF_TOL 1e-4f   /* accepted upper bound vs hand-computed reference */

static void die(const char* msg, const char* file, int line) {
    fprintf(stderr, "error: %s (%s:%d)\n", msg, file, line);
    exit(2);
}

/* ---- curve (same evaluation as ../curve/src/curve.c) ---- */

typedef struct { float t, v, tan_in, tan_out; } Key;
typedef struct { Key keys[MAX_KEYS]; int count; } Curve;

static int curve_add(Curve* c, float t, float v, float tin, float tout) {
    if (c->count >= MAX_KEYS) return -1;
    if (c->count > 0 && t <= c->keys[c->count - 1].t) return -1;
    Key* k = &c->keys[c->count++];
    k->t = t; k->v = v; k->tan_in = tin; k->tan_out = tout;
    return 0;
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

/* ---- skeleton (same evaluation as ../skeleton/src/skeleton.c) ---- */

typedef struct { float x, y, ang; } Pose;

typedef struct {
    int parent[MAX_BONES];
    Pose local[MAX_BONES];
    Pose world[MAX_BONES];
    int count;
} Skeleton;

static void rotate(float x, float y, float a, float* ox, float* oy) {
    float c = cosf(a), s = sinf(a);
    *ox = x * c - y * s;
    *oy = x * s + y * c;
}

static void skel_update(Skeleton* s) {
    for (int i = 0; i < s->count; i++) {
        if (s->parent[i] < 0) {
            s->world[i] = s->local[i];
        } else {
            int p = s->parent[i];
            float rx, ry;
            rotate(s->local[i].x, s->local[i].y, s->world[p].ang, &rx, &ry);
            s->world[i].x = s->world[p].x + rx;
            s->world[i].y = s->world[p].y + ry;
            s->world[i].ang = s->world[p].ang + s->local[i].ang;
        }
    }
}

/* ---- animation: skeleton rest pose + per-bone channel curves ---- */

typedef struct {
    int parent[MAX_BONES];
    Pose rest[MAX_BONES];
    Curve cx[MAX_BONES], cy[MAX_BONES], ca[MAX_BONES];
    int bones;
    int frames;
    float dt, t0;
} Anim;

static void anim_init(Anim* a) { memset(a, 0, sizeof(*a)); }

static float chan_eval(const Curve* c, float rest, float t) {
    return c->count > 0 ? curve_eval(c, t) : rest;
}

/* ---- recording: world pose of every bone at every frame ---- */

typedef struct {
    int frames, bones;
    Pose p[MAX_FRAMES][MAX_BONES];
} Rec;

static void play(const Anim* a, Rec* r) {
    r->frames = a->frames;
    r->bones = a->bones;
    Skeleton s;
    memset(&s, 0, sizeof(s));
    s.count = a->bones;
    for (int b = 0; b < a->bones; b++) s.parent[b] = a->parent[b];
    for (int f = 0; f < a->frames; f++) {
        float t = a->t0 + f * a->dt;
        for (int b = 0; b < a->bones; b++) {
            s.local[b].x = chan_eval(&a->cx[b], a->rest[b].x, t);
            s.local[b].y = chan_eval(&a->cy[b], a->rest[b].y, t);
            s.local[b].ang = chan_eval(&a->ca[b], a->rest[b].ang, t);
        }
        skel_update(&s);
        for (int b = 0; b < a->bones; b++) r->p[f][b] = s.world[b];
    }
}

/* ---- parsing ---- */

static float parse_float(const char* tok, const char* file, int line) {
    if (!tok) die("missing numeric field", file, line);
    char* end = NULL;
    float v = strtof(tok, &end);
    if (end == tok || *end != '\0') die("malformed number", file, line);
    if (isnan(v)) die("NaN in input", file, line);
    if (isinf(v)) die("Inf in input", file, line);
    return v;
}

static int parse_int(const char* tok, const char* file, int line) {
    if (!tok) die("missing integer field", file, line);
    char* end = NULL;
    long v = strtol(tok, &end, 10);
    if (end == tok || *end != '\0') die("malformed integer", file, line);
    return (int)v;
}

static void load_anim(const char* path, Anim* a) {
    FILE* fp = fopen(path, "r");   /* samples are read-only for us */
    if (!fp) { fprintf(stderr, "error: cannot open %s\n", path); exit(2); }
    anim_init(a);
    char buf[256];
    int line = 0, have_frames = 0;
    while (fgets(buf, sizeof(buf), fp)) {
        line++;
        char* tok[8];
        int n = 0;
        for (char* p = strtok(buf, " \t\r\n"); p && n < 8; p = strtok(NULL, " \t\r\n"))
            tok[n++] = p;
        if (n == 0 || tok[0][0] == '#') continue;
        if (strcmp(tok[0], "bone") == 0) {
            if (n != 5) die("bone: want 4 fields", path, line);
            if (a->bones >= MAX_BONES) die("too many bones", path, line);
            int parent = parse_int(tok[1], path, line);
            if (parent >= a->bones || parent < -1) die("bad parent index", path, line);
            int b = a->bones++;
            a->parent[b] = parent;
            a->rest[b].x = parse_float(tok[2], path, line);
            a->rest[b].y = parse_float(tok[3], path, line);
            a->rest[b].ang = parse_float(tok[4], path, line);
        } else if (strcmp(tok[0], "key") == 0) {
            if (n != 7) die("key: want 6 fields", path, line);
            int b = parse_int(tok[1], path, line);
            if (b < 0 || b >= a->bones) die("key: bad bone index", path, line);
            float t = parse_float(tok[3], path, line);
            if (t < 0.0f) die("negative key time", path, line);
            float v = parse_float(tok[4], path, line);
            float tin = parse_float(tok[5], path, line);
            float tout = parse_float(tok[6], path, line);
            Curve* c;
            if (strcmp(tok[2], "x") == 0) c = &a->cx[b];
            else if (strcmp(tok[2], "y") == 0) c = &a->cy[b];
            else if (strcmp(tok[2], "ang") == 0) c = &a->ca[b];
            else die("key: unknown channel", path, line);
            if (curve_add(c, t, v, tin, tout) != 0)
                die("key: curve full or times not increasing", path, line);
        } else if (strcmp(tok[0], "frames") == 0) {
            if (n != 4) die("frames: want 3 fields", path, line);
            a->frames = parse_int(tok[1], path, line);
            a->dt = parse_float(tok[2], path, line);
            a->t0 = parse_float(tok[3], path, line);
            if (a->frames <= 0 || a->frames > MAX_FRAMES) die("frames: bad count", path, line);
            if (a->dt <= 0.0f) die("frames: dt must be positive", path, line);
            if (a->t0 < 0.0f) die("frames: negative start time", path, line);
            have_frames = 1;
        } else {
            die("unknown directive", path, line);
        }
    }
    fclose(fp);
    if (a->bones == 0) die("no bones defined", path, line);
    if (!have_frames) die("no frames directive", path, line);
}

/* ---- recording file I/O (same format as the hand-computed reference) ---- */

static void write_rec(const char* path, const Rec* r) {
    FILE* fp = fopen(path, "w");
    if (!fp) { fprintf(stderr, "error: cannot write %s\n", path); exit(2); }
    fprintf(fp, "rec %d %d\n", r->frames, r->bones);
    for (int f = 0; f < r->frames; f++)
        for (int b = 0; b < r->bones; b++)
            fprintf(fp, "pose %d %d %.9g %.9g %.9g\n", f, b,
                    r->p[f][b].x, r->p[f][b].y, r->p[f][b].ang);
    fclose(fp);
}

static void load_rec(const char* path, Rec* r) {
    FILE* fp = fopen(path, "r");
    if (!fp) { fprintf(stderr, "error: cannot open %s\n", path); exit(2); }
    memset(r, 0, sizeof(*r));
    char buf[256];
    int line = 0, have_header = 0;
    while (fgets(buf, sizeof(buf), fp)) {
        line++;
        char* tok[8];
        int n = 0;
        for (char* p = strtok(buf, " \t\r\n"); p && n < 8; p = strtok(NULL, " \t\r\n"))
            tok[n++] = p;
        if (n == 0 || tok[0][0] == '#') continue;
        if (strcmp(tok[0], "rec") == 0) {
            if (n != 3) die("rec: want 2 fields", path, line);
            r->frames = parse_int(tok[1], path, line);
            r->bones = parse_int(tok[2], path, line);
            if (r->frames <= 0 || r->frames > MAX_FRAMES ||
                r->bones <= 0 || r->bones > MAX_BONES)
                die("rec: bad dimensions", path, line);
            have_header = 1;
        } else if (strcmp(tok[0], "pose") == 0) {
            if (!have_header) die("pose before rec header", path, line);
            if (n != 6) die("pose: want 5 fields", path, line);
            int f = parse_int(tok[1], path, line);
            int b = parse_int(tok[2], path, line);
            if (f < 0 || f >= r->frames || b < 0 || b >= r->bones)
                die("pose: index out of range", path, line);
            r->p[f][b].x = parse_float(tok[3], path, line);
            r->p[f][b].y = parse_float(tok[4], path, line);
            r->p[f][b].ang = parse_float(tok[5], path, line);
        } else {
            die("unknown directive", path, line);
        }
    }
    fclose(fp);
    if (!have_header) die("no rec header", path, line);
}

/* ---- diff: localize every mismatch to (frame, bone, channel) ---- */

static const char* CHAN[3] = { "x", "y", "ang" };

static float pose_component(const Pose* p, int c) {
    return c == 0 ? p->x : c == 1 ? p->y : p->ang;
}

/* Returns number of differing components; prints each mismatch and a
 * summary with the max abs error and where it occurs. */
static int diff_recs(const char* na, const Rec* a, const char* nb, const Rec* b,
                     float tol, int bitwise) {
    if (a->frames != b->frames || a->bones != b->bones) {
        printf("dim_mismatch %s=%dx%d %s=%dx%d\n", na, a->frames, a->bones,
               nb, b->frames, b->bones);
        return -1;
    }
    int diffs = 0;
    float max_err = 0.0f;
    int max_f = -1, max_b = -1, max_c = -1;
    for (int f = 0; f < a->frames; f++) {
        for (int bn = 0; bn < a->bones; bn++) {
            for (int c = 0; c < 3; c++) {
                float va = pose_component(&a->p[f][bn], c);
                float vb = pose_component(&b->p[f][bn], c);
                float err = fabsf(va - vb);
                int bad = bitwise ? (memcmp(&va, &vb, sizeof(float)) != 0)
                                  : (err > tol);
                if (bad) {
                    diffs++;
                    printf("diff frame=%d bone=%d chan=%s %s=%.9g %s=%.9g err=%.9g\n",
                           f, bn, CHAN[c], na, va, nb, vb, err);
                }
                if (err > max_err) { max_err = err; max_f = f; max_b = bn; max_c = c; }
            }
        }
    }
    printf("max_abs_err=%.9g at frame=%d bone=%d chan=%s\n",
           max_err, max_f, max_b, max_c >= 0 ? CHAN[max_c] : "-");
    printf("diff_count=%d\n", diffs);
    return diffs;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static void usage(void) {
    fprintf(stderr,
        "usage:\n"
        "  replay check  <anim> <ref.rec>\n"
        "  replay record <anim> <out.rec>\n"
        "  replay diff   <a.rec> <b.rec>\n");
    exit(2);
}

int main(int argc, char** argv) {
    if (argc < 2) usage();

    if (strcmp(argv[1], "check") == 0 && argc == 4) {
        Anim a;
        load_anim(argv[2], &a);
        printf("input_ok=1 bones=%d frames=%d\n", a.bones, a.frames);

        /* determinism: same input, played twice */
        static Rec run1, run2;
        play(&a, &run1);
        play(&a, &run2);
        int d = diff_recs("run1", &run1, "run2", &run2, 0.0f, 1);
        check("replay_deterministic", d == 0);

        /* accuracy: run1 vs hand-computed reference, error upper bound */
        static Rec ref;
        load_rec(argv[3], &ref);
        d = diff_recs("replay", &run1, "ref", &ref, REF_TOL, 0);
        check("within_ref_tol", d == 0);

        printf("all_ok=%d\n", failures == 0);
        return failures == 0 ? 0 : 1;
    }

    if (strcmp(argv[1], "record") == 0 && argc == 4) {
        Anim a;
        load_anim(argv[2], &a);
        static Rec r;
        play(&a, &r);
        write_rec(argv[3], &r);
        printf("recorded frames=%d bones=%d -> %s\n", r.frames, r.bones, argv[3]);
        return 0;
    }

    if (strcmp(argv[1], "diff") == 0 && argc == 4) {
        static Rec a, b;
        load_rec(argv[2], &a);
        load_rec(argv[3], &b);
        int d = diff_recs("a", &a, "b", &b, 0.0f, 1);
        printf("identical=%d\n", d == 0);
        return d == 0 ? 0 : 1;
    }

    usage();
    return 2;
}
