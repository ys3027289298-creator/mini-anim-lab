#define _POSIX_C_SOURCE 200809L

/* constraint: verlet particle chain with distance constraints, iteration
   count and damping; used as the toy model for cloth/rope style solvers.

   Debug export (off unless dbg_open() succeeds):
     residuals.csv  one row per constraint per frame: length, rest, residual
     worst.csv      one row per frame: index and value of the worst constraint
     heatmap.ppm    residual heat map, frame on x, constraint on y; the colour
                    ramp is green (0) -> red (DBG_RAMP_MAX), clamped above.
   dbg_verify() reads the csv and the ppm back and reports an error when the
   painted colour does not decode to the recorded residual. The solver never
   reads from DebugExport, so export on/off gives bit-identical particles. */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>

#define MAX_P 16

typedef struct { float x, y, px, py; float inv_mass; int pinned; } Particle;

typedef struct {
    Particle p[MAX_P];
    int n;
    float rest[MAX_P];      /* rest length between i and i+1 */
    int iterations;
    float damping;
} Chain;

static void chain_init(Chain* c, int n, float seg, float x0, float y0) {
    memset(c, 0, sizeof(*c));
    c->n = n;
    c->iterations = 8;
    c->damping = 0.99f;
    for (int i = 0; i < n; i++) {
        c->p[i].x = x0 + seg * i; c->p[i].y = y0;
        c->p[i].px = c->p[i].x; c->p[i].py = c->p[i].y;
        c->p[i].inv_mass = 1.0f;
        c->p[i].pinned = 0;
        c->rest[i] = seg;
    }
    c->p[0].pinned = 1;
    c->p[0].inv_mass = 0.0f;
}

static void chain_step(Chain* c, float dt, float gx, float gy) {
    for (int i = 0; i < c->n; i++) {
        Particle* p = &c->p[i];
        if (p->pinned) continue;
        float vx = (p->x - p->px) * c->damping;
        float vy = (p->y - p->py) * c->damping;
        p->px = p->x; p->py = p->y;
        p->x += vx + gx * dt * dt;
        p->y += vy + gy * dt * dt;
    }
    for (int it = 0; it < c->iterations; it++) {
        for (int i = 0; i < c->n - 1; i++) {
            Particle* a = &c->p[i];
            Particle* b = &c->p[i + 1];
            float dx = b->x - a->x, dy = b->y - a->y;
            float d = sqrtf(dx * dx + dy * dy);
            if (d < 1e-6f) continue;
            float diff = (d - c->rest[i]) / d;
            float wa = a->inv_mass, wb = b->inv_mass;
            float wsum = wa + wb;
            if (wsum <= 0) continue;
            a->x += dx * diff * (wa / wsum);
            a->y += dy * diff * (wa / wsum);
            b->x -= dx * diff * (wb / wsum);
            b->y -= dy * diff * (wb / wsum);
        }
    }
}

/* ---- debug export ------------------------------------------------------- */

#define DBG_CELL 16          /* pixels per (frame, constraint) cell */
#define DBG_RAMP_MAX 0.05f   /* residual that maps to the hottest colour */
#define DBG_DUMP_FRAMES 1024
#define DBG_DUMP_COLS 32

typedef struct {
    int on;
    char dir[512];
    FILE* residuals;
    FILE* worst;
    int frames;
    int cols;
    unsigned char* pixels;   /* 3 bytes per (frame, constraint) cell */
    size_t cap;
    int error;
} DebugExport;

/* Residual of constraint i: |d - rest|. chain_seg_error() uses the same
   helper, so the aggregate and the per-constraint numbers cannot drift. */
static float chain_residual(const Chain* c, int i, float* len, float* err) {
    float dx = c->p[i + 1].x - c->p[i].x;
    float dy = c->p[i + 1].y - c->p[i].y;
    float d = sqrtf(dx * dx + dy * dy);
    float e = fabsf(d - c->rest[i]);
    if (len) *len = d;
    if (err) *err = e;
    return e;
}

static float chain_seg_error(const Chain* c) {
    float worst = 0.0f;
    for (int i = 0; i < c->n - 1; i++) {
        float e = chain_residual(c, i, NULL, NULL);
        if (e > worst) worst = e;
    }
    return worst;
}

/* green (0) -> yellow -> red (>= DBG_RAMP_MAX).
   The colour is quantised to a level L in [0,255] and the green channel is
   derived from L, so red alone recovers L exactly on read-back. */
static int dbg_level(float e) {
    float t = e / DBG_RAMP_MAX;
    if (!(t >= 0.0f)) t = 0.0f;          /* also catches NaN */
    if (t > 1.0f) t = 1.0f;
    return (int)(t * 255.0f + 0.5f);
}

static void dbg_colour(float e, int* r, int* g, int* b) {
    int l = dbg_level(e);
    *r = l;
    *g = (int)(180.0f * (float)(255 - l) / 255.0f + 0.5f);
    *b = 0;
}

static int dbg_join(char* out, size_t n, const char* dir, const char* name) {
    int k = snprintf(out, n, "%s/%s", dir, name);
    return (k > 0 && (size_t)k < n) ? 0 : -1;
}

static int dbg_mkdir_p(const char* dir) {
    char buf[1024];
    size_t n = strlen(dir);
    if (n == 0 || n >= sizeof(buf)) return -1;
    memcpy(buf, dir, n + 1);
    for (size_t i = 1; i < n; i++) {
        if (buf[i] != '/') continue;
        buf[i] = '\0';
        if (mkdir(buf, 0755) != 0) {
            struct stat st;
            if (stat(buf, &st) != 0 || !S_ISDIR(st.st_mode)) return -1;
        }
        buf[i] = '/';
    }
    if (mkdir(buf, 0755) != 0) {
        struct stat st;
        if (stat(buf, &st) != 0 || !S_ISDIR(st.st_mode)) return -1;
    }
    return 0;
}

static int dbg_open(DebugExport* d, const char* dir) {
    memset(d, 0, sizeof(*d));
    if (dbg_mkdir_p(dir) != 0) {
        fprintf(stderr, "debug: cannot create %s\n", dir);
        return -1;
    }
    snprintf(d->dir, sizeof(d->dir), "%s", dir);

    char path[1024];
    if (dbg_join(path, sizeof(path), dir, "residuals.csv") != 0) return -1;
    d->residuals = fopen(path, "w");
    if (!d->residuals) { fprintf(stderr, "debug: cannot write %s\n", path); return -1; }
    fprintf(d->residuals, "frame,constraint,particle_a,particle_b,length,rest,residual\n");

    if (dbg_join(path, sizeof(path), dir, "worst.csv") != 0) goto fail;
    d->worst = fopen(path, "w");
    if (!d->worst) { fprintf(stderr, "debug: cannot write %s\n", path); goto fail; }
    fprintf(d->worst, "frame,constraint,residual\n");

    d->on = 1;
    return 0;

fail:
    fclose(d->residuals);
    d->residuals = NULL;
    return -1;
}

/* Called once per frame after chain_step(). Read-only with respect to the
   chain: enabling export cannot perturb the simulation. */
static void dbg_record(DebugExport* d, const Chain* c) {
    if (!d->on || d->error) return;
    int cols = c->n - 1;
    if (cols < 1) return;
    if (d->cols == 0) d->cols = cols;
    if (d->cols != cols) {
        fprintf(stderr, "debug: chain length changed mid-export (%d -> %d)\n", d->cols, cols);
        d->error = 1;
        return;
    }

    size_t need = (size_t)(d->frames + 1) * (size_t)d->cols * 3u;
    if (need > d->cap) {
        size_t cap = d->cap ? d->cap : (size_t)d->cols * 3u * 64u;
        while (cap < need) cap *= 2;
        unsigned char* p = realloc(d->pixels, cap);
        if (!p) { fprintf(stderr, "debug: out of memory\n"); d->error = 1; return; }
        d->pixels = p;
        d->cap = cap;
    }

    int f = d->frames;
    int worst_i = 0;
    float worst_e = 0.0f;
    for (int i = 0; i < cols; i++) {
        float len = 0.0f, e = 0.0f;
        chain_residual(c, i, &len, &e);
        if (i == 0 || e > worst_e) { worst_e = e; worst_i = i; }
        if (fprintf(d->residuals, "%d,%d,%d,%d,%.9g,%.9g,%.9g\n",
                    f, i, i, i + 1, (double)len, (double)c->rest[i], (double)e) < 0) {
            d->error = 1;
        }
        int r, g, b;
        dbg_colour(e, &r, &g, &b);
        unsigned char* px = &d->pixels[((size_t)f * (size_t)d->cols + (size_t)i) * 3u];
        px[0] = (unsigned char)r;
        px[1] = (unsigned char)g;
        px[2] = (unsigned char)b;
    }
    if (fprintf(d->worst, "%d,%d,%.9g\n", f, worst_i, (double)worst_e) < 0) d->error = 1;
    d->frames++;
}

static int dbg_write_heatmap(const DebugExport* d) {
    char path[1024];
    if (dbg_join(path, sizeof(path), d->dir, "heatmap.ppm") != 0) return -1;
    FILE* f = fopen(path, "wb");
    if (!f) { fprintf(stderr, "debug: cannot write %s\n", path); return -1; }

    int w = d->frames * DBG_CELL;
    int h = d->cols * DBG_CELL;
    fprintf(f, "P6\n%d %d\n255\n", w, h);
    unsigned char* line = malloc((size_t)w * 3u);
    if (!line) { fclose(f); return -1; }

    for (int ci = 0; ci < d->cols; ci++) {
        for (int py = 0; py < DBG_CELL; py++) {
            for (int fr = 0; fr < d->frames; fr++) {
                const unsigned char* px = &d->pixels[((size_t)fr * (size_t)d->cols + (size_t)ci) * 3u];
                for (int k = 0; k < DBG_CELL; k++) {
                    memcpy(line + ((size_t)fr * DBG_CELL + (size_t)k) * 3u, px, 3);
                }
            }
            if (fwrite(line, 3u, (size_t)w, f) != (size_t)w) {
                free(line); fclose(f);
                fprintf(stderr, "debug: short write to %s\n", path);
                return -1;
            }
        }
    }
    free(line);
    if (fclose(f) != 0) return -1;
    return 0;
}

static int dbg_close(DebugExport* d) {
    if (!d->on && !d->residuals && !d->worst) return 0;
    int rc = d->error ? -1 : 0;
    if (d->residuals && fclose(d->residuals) != 0) rc = -1;
    if (d->worst && fclose(d->worst) != 0) rc = -1;
    d->residuals = NULL;
    d->worst = NULL;
    if (!d->error && dbg_write_heatmap(d) != 0) rc = -1;
    free(d->pixels);
    d->pixels = NULL;
    d->on = 0;
    return rc;
}

/* ---- reading the export back for verification --------------------------- */

typedef struct {
    int frames;
    int cols;
    float residual[DBG_DUMP_FRAMES][DBG_DUMP_COLS];
    float length[DBG_DUMP_FRAMES][DBG_DUMP_COLS];
    float rest[DBG_DUMP_FRAMES][DBG_DUMP_COLS];
    int pa[DBG_DUMP_FRAMES][DBG_DUMP_COLS];
    int pb[DBG_DUMP_FRAMES][DBG_DUMP_COLS];
    int filled[DBG_DUMP_FRAMES][DBG_DUMP_COLS];
    int worst_i[DBG_DUMP_FRAMES];
    float worst_e[DBG_DUMP_FRAMES];
    int worst_filled[DBG_DUMP_FRAMES];
} DbgDump;

static int dbg_read_residuals(const char* path, DbgDump* d) {
    FILE* f = fopen(path, "r");
    if (!f) return -1;
    char line[512];
    if (!fgets(line, sizeof(line), f)) { fclose(f); return -1; }
    int maxf = -1, maxc = -1;
    while (fgets(line, sizeof(line), f)) {
        int fr, ci, a, b;
        float len, rest, e;
        if (sscanf(line, "%d,%d,%d,%d,%f,%f,%f", &fr, &ci, &a, &b, &len, &rest, &e) != 7) {
            fclose(f);
            return -1;
        }
        if (fr < 0 || fr >= DBG_DUMP_FRAMES || ci < 0 || ci >= DBG_DUMP_COLS) {
            fclose(f);
            return -1;
        }
        if (d->filled[fr][ci]) { fclose(f); return -1; }   /* duplicate row */
        d->residual[fr][ci] = e;
        d->length[fr][ci] = len;
        d->rest[fr][ci] = rest;
        d->pa[fr][ci] = a;
        d->pb[fr][ci] = b;
        d->filled[fr][ci] = 1;
        if (fr > maxf) maxf = fr;
        if (ci > maxc) maxc = ci;
    }
    fclose(f);
    if (maxf < 0 || maxc < 0) return -1;
    d->frames = maxf + 1;
    d->cols = maxc + 1;
    for (int fr = 0; fr < d->frames; fr++) {
        for (int ci = 0; ci < d->cols; ci++) {
            if (!d->filled[fr][ci]) return -1;
        }
    }
    return 0;
}

static int dbg_read_worst(const char* path, DbgDump* d) {
    FILE* f = fopen(path, "r");
    if (!f) return -1;
    char line[512];
    if (!fgets(line, sizeof(line), f)) { fclose(f); return -1; }
    int rows = 0;
    while (fgets(line, sizeof(line), f)) {
        int fr, ci;
        float e;
        if (sscanf(line, "%d,%d,%f", &fr, &ci, &e) != 3) { fclose(f); return -1; }
        if (fr < 0 || fr >= DBG_DUMP_FRAMES || ci < 0 || ci >= DBG_DUMP_COLS) {
            fclose(f);
            return -1;
        }
        if (d->worst_filled[fr]) { fclose(f); return -1; }
        d->worst_i[fr] = ci;
        d->worst_e[fr] = e;
        d->worst_filled[fr] = 1;
        rows++;
    }
    fclose(f);
    return rows == d->frames ? 0 : -1;
}

typedef struct {
    int w, h;
    unsigned char* data;
} Ppm;

static int ppm_read(const char* path, Ppm* out) {
    memset(out, 0, sizeof(*out));
    FILE* f = fopen(path, "rb");
    if (!f) return -1;
    char magic[4] = { 0 };
    int maxv = 0;
    if (fscanf(f, "%3s", magic) != 1 || strcmp(magic, "P6") != 0 ||
        fscanf(f, "%d %d %d", &out->w, &out->h, &maxv) != 3 ||
        maxv != 255 || out->w <= 0 || out->h <= 0) {
        fclose(f);
        return -1;
    }
    int c = fgetc(f);                      /* exactly one whitespace byte */
    if (c != ' ' && c != '\n' && c != '\r' && c != '\t') { fclose(f); return -1; }
    if (c == '\r' && fgetc(f) != '\n') { fclose(f); return -1; }
    size_t bytes = (size_t)out->w * (size_t)out->h * 3u;
    out->data = malloc(bytes);
    if (!out->data) { fclose(f); return -1; }
    if (fread(out->data, 1, bytes, f) != bytes) {
        free(out->data);
        out->data = NULL;
        fclose(f);
        return -1;
    }
    fclose(f);
    return 0;
}

/* Read the numeric export and the heat map and confirm that every colour
   decodes back to the residual recorded for that (frame, constraint).
   Returns the number of mismatches; prints one error line if any. */
static int dbg_verify(const char* dir) {
    char path[1024];
    DbgDump* d = calloc(1, sizeof(*d));
    if (!d) {
        fprintf(stderr, "verify: out of memory\n");
        return -1;
    }

    if (dbg_join(path, sizeof(path), dir, "residuals.csv") != 0 ||
        dbg_read_residuals(path, d) != 0) {
        fprintf(stderr, "verify: cannot parse %s/residuals.csv\n", dir);
        free(d);
        return 1;
    }
    if (dbg_join(path, sizeof(path), dir, "worst.csv") != 0 || dbg_read_worst(path, d) != 0) {
        fprintf(stderr, "verify: cannot parse %s/worst.csv\n", dir);
        free(d);
        return 1;
    }
    Ppm ppm;
    if (dbg_join(path, sizeof(path), dir, "heatmap.ppm") != 0 || ppm_read(path, &ppm) != 0) {
        fprintf(stderr, "verify: cannot parse %s/heatmap.ppm\n", dir);
        free(d);
        return 1;
    }
    if (ppm.w != d->frames * DBG_CELL || ppm.h != d->cols * DBG_CELL) {
        fprintf(stderr, "verify: heatmap is %dx%d, expected %dx%d (frames x constraints)\n",
                ppm.w, ppm.h, d->frames * DBG_CELL, d->cols * DBG_CELL);
        free(ppm.data);
        free(d);
        return 1;
    }

    int bad = 0;
    for (int fr = 0; fr < d->frames && bad < 8; fr++) {
        int wi = 0;
        for (int ci = 0; ci < d->cols; ci++) {
            float e = d->residual[fr][ci];
            float re = fabsf(d->length[fr][ci] - d->rest[fr][ci]);
            if (fabsf(re - e) > 1e-4f * (1.0f + fabsf(e))) {
                fprintf(stderr, "verify: frame %d constraint %d residual %.9g but |len-rest| = %.9g\n",
                        fr, ci, (double)e, (double)re);
                bad++;
            }
            if (d->residual[fr][ci] > d->residual[fr][wi]) wi = ci;

            /* every pixel of the cell must carry the ramp colour for this
               residual -- a single wrong pixel is a mismatch */
            int l = -1;
            int cell_ok = 1;
            for (int py = 0; py < DBG_CELL && cell_ok; py++) {
                for (int px = 0; px < DBG_CELL && cell_ok; px++) {
                    size_t off = (((size_t)ci * DBG_CELL + (size_t)py) * (size_t)ppm.w +
                                  ((size_t)fr * DBG_CELL + (size_t)px)) * 3u;
                    int li = ppm.data[off];
                    int g = (int)(180.0f * (float)(255 - li) / 255.0f + 0.5f);
                    if (ppm.data[off + 1] != g || ppm.data[off + 2] != 0) {
                        fprintf(stderr, "verify: frame %d constraint %d pixel (%d,%d) colour (%d,%d,%d) is not on the ramp\n",
                                fr, ci, px, py, ppm.data[off], ppm.data[off + 1], ppm.data[off + 2]);
                        cell_ok = 0;
                    } else if (l < 0) {
                        l = li;
                    } else if (l != li) {
                        fprintf(stderr, "verify: frame %d constraint %d cell is not a single colour\n", fr, ci);
                        cell_ok = 0;
                    }
                }
            }
            if (!cell_ok) { bad++; continue; }
            if (l != dbg_level(e)) {
                float decoded = (float)l / 255.0f * DBG_RAMP_MAX;
                fprintf(stderr, "verify: frame %d constraint %d residual %.9g but colour decodes to %.9g\n",
                        fr, ci, (double)e, (double)decoded);
                bad++;
            }
        }
        if (d->worst_i[fr] != wi || fabsf(d->worst_e[fr] - d->residual[fr][wi]) > 1e-6f * (1.0f + fabsf(d->worst_e[fr]))) {
            fprintf(stderr, "verify: frame %d worst recorded as c%d %.9g, per-constraint rows say c%d %.9g\n",
                    fr, d->worst_i[fr], (double)d->worst_e[fr], wi, (double)d->residual[fr][wi]);
            bad++;
        }
    }

    if (bad == 0 && d->frames > 0) {
        printf("verify_ok=1 frames=%d constraints=%d\n", d->frames, d->cols);
    }
    free(ppm.data);
    free(d);
    return bad;
}

/* A handful of FNV-1a hashes over the particle state after each frame.
   Same chain, same frames, export on or off -> same hashes. */
static void chain_state_hash(const Chain* c, char* out, size_t n) {
    uint64_t h = 1469598103934665603ull;
    const unsigned char* b = (const unsigned char*)c->p;
    for (size_t i = 0; i < sizeof(c->p); i++) {
        h ^= b[i];
        h *= 1099511628211ull;
    }
    snprintf(out, n, "%016llx", (unsigned long long)h);
}

static void run_scenario(Chain* c, int frames) {
    for (int i = 0; i < frames; i++) chain_step(c, 0.016f, 0.0f, -9.8f);
}

/* Independent recomputation: distances measured straight from the particle
   positions in double precision, rather than through chain_residual(). Used
   to confirm the export equals a hand-rolled per-constraint pass. */
static int manual_matches_export(const Chain* snaps, int frames, const DbgDump* d) {
    int bad = 0;
    for (int fr = 0; fr < frames && bad < 4; fr++) {
        double worst = -1.0;
        int worst_i = 0;
        for (int ci = 0; ci < d->cols; ci++) {
            double dx = (double)snaps[fr].p[ci + 1].x - (double)snaps[fr].p[ci].x;
            double dy = (double)snaps[fr].p[ci + 1].y - (double)snaps[fr].p[ci].y;
            double len = sqrt(dx * dx + dy * dy);
            double e = fabs(len - (double)snaps[fr].rest[ci]);
            if (fabs(len - (double)d->length[fr][ci]) > 1e-4 * (1.0 + len)) {
                fprintf(stderr, "manual: frame %d constraint %d length %.9g, export %.9g\n",
                        fr, ci, len, (double)d->length[fr][ci]);
                bad++;
            }
            if (fabs(e - (double)d->residual[fr][ci]) > 1e-4 * (1.0 + e)) {
                fprintf(stderr, "manual: frame %d constraint %d residual %.9g, export %.9g\n",
                        fr, ci, e, (double)d->residual[fr][ci]);
                bad++;
            }
            if (e > worst) { worst = e; worst_i = ci; }
        }
        if (fabs(worst - (double)d->worst_e[fr]) > 1e-4 * (1.0 + worst)) {
            fprintf(stderr, "manual: frame %d worst %.9g, export %.9g\n", fr, worst, (double)d->worst_e[fr]);
            bad++;
        }
        /* index may tie within a frame; only require a valid, maximal one */
        if (d->worst_i[fr] < 0 || d->worst_i[fr] >= d->cols ||
            fabs((double)d->residual[fr][d->worst_i[fr]] - worst) > 1e-4 * (1.0 + worst)) {
            fprintf(stderr, "manual: frame %d worst index %d is not maximal (c%d)\n",
                    fr, d->worst_i[fr], worst_i);
            bad++;
        }
    }
    return bad;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    Chain c;
    chain_init(&c, 8, 1.0f, 0.0f, 0.0f);
    check("init_len", c.n == 8 && fabsf(chain_seg_error(&c)) < 1e-5f);

    /* let it fall for a while */
    for (int i = 0; i < 200; i++) chain_step(&c, 0.016f, 0.0f, -9.8f);
    check("stays_in_stretch", chain_seg_error(&c) < 0.02f);
    check("pinned_anchor", fabsf(c.p[0].x) < 1e-6f && fabsf(c.p[0].y) < 1e-6f);
    check("falls_down", c.p[7].y < -1.0f);

    /* more iterations must not increase error */
    Chain a, b;
    chain_init(&a, 8, 1.0f, 0.0f, 0.0f);
    chain_init(&b, 8, 1.0f, 0.0f, 0.0f);
    a.iterations = 1; b.iterations = 16;
    for (int i = 0; i < 100; i++) { chain_step(&a, 0.016f, 0.0f, -9.8f); chain_step(&b, 0.016f, 0.0f, -9.8f); }
    check("more_iters_better", chain_seg_error(&b) <= chain_seg_error(&a) + 1e-6f);

    Chain cap;
    chain_init(&cap, 16, 0.5f, 0.0f, 0.0f);
    check("capacity_ok", cap.n == 16);

    /* ---- debug export ---- */

    /* reference run, export off */
    Chain off;
    chain_init(&off, 8, 1.0f, 0.0f, 0.0f);
    run_scenario(&off, 150);
    char h_off[32];
    chain_state_hash(&off, h_off, sizeof(h_off));

    /* same run, export on; keep quiet unless open fails */
    const char* dir = "out/debug";
    enum { FRAMES = 150 };
    Chain on;
    chain_init(&on, 8, 1.0f, 0.0f, 0.0f);
    DebugExport dbg;
    int opened = dbg_open(&dbg, dir) == 0;
    check("debug_opens", opened);
    if (opened) {
        for (int i = 0; i < FRAMES; i++) {
            chain_step(&on, 0.016f, 0.0f, -9.8f);
            dbg_record(&dbg, &on);
        }
    }
    int closed = !opened || dbg_close(&dbg) == 0;
    check("debug_writes", opened && closed);
    if (!opened) {
        /* the export checks below have nothing to read; fail the run */
        printf("all_ok=0\n");
        return 1;
    }

    char h_on[32];
    chain_state_hash(&on, h_on, sizeof(h_on));
    check("export_bit_identical", strcmp(h_off, h_on) == 0);

    /* the written export must parse into the shape we just recorded */
    DbgDump* dump = calloc(1, sizeof(*dump));
    char path[1024];
    snprintf(path, sizeof(path), "%s/residuals.csv", dir);
    int parsed = dump && dbg_read_residuals(path, dump) == 0;
    snprintf(path, sizeof(path), "%s/worst.csv", dir);
    parsed = parsed && dbg_read_worst(path, dump) == 0;
    check("export_parse", parsed);
    if (parsed) {
        check("export_shape", dump->frames == FRAMES && dump->cols == 7);
        int rows_ok = 1, worst_ok = 1;
        for (int fr = 0; fr < dump->frames; fr++) {
            int wi = 0;
            for (int ci = 0; ci < dump->cols; ci++) {
                if (dump->pa[fr][ci] != ci || dump->pb[fr][ci] != ci + 1) rows_ok = 0;
                if (dump->residual[fr][ci] > dump->residual[fr][wi]) wi = ci;
            }
            if (dump->worst_i[fr] != wi) worst_ok = 0;
        }
        check("export_row_labels", rows_ok);
        check("export_worst_matches_rows", worst_ok);

        /* the last exported frame must match the reference run, which ended
           on the same frame with export off */
        int ends_on_live = 1;
        for (int ci = 0; ci < dump->cols; ci++) {
            float len = 0.0f, e = 0.0f;
            chain_residual(&off, ci, &len, &e);
            if (fabsf(dump->length[dump->frames - 1][ci] - len) > 1e-4f * (1.0f + fabsf(len)) ||
                fabsf(dump->residual[dump->frames - 1][ci] - e) > 1e-4f * (1.0f + fabsf(e))) {
                ends_on_live = 0;
            }
        }
        check("export_matches_sim", ends_on_live);

        /* snapshot every frame up front so the manual pass sees the same state */
        Chain* snaps = malloc(sizeof(*snaps) * FRAMES);
        int snap_ok = snaps != NULL;
        if (snap_ok) {
            chain_init(&on, 8, 1.0f, 0.0f, 0.0f);
            for (int fr = 0; fr < FRAMES; fr++) {
                chain_step(&on, 0.016f, 0.0f, -9.8f);
                snaps[fr] = on;
            }
            char h_snap[32];
            chain_state_hash(&on, h_snap, sizeof(h_snap));
            snap_ok = strcmp(h_snap, h_on) == 0;
        }
        check("export_manual_consistent", snap_ok && manual_matches_export(snaps, FRAMES, dump) == 0);
        free(snaps);
    } else {
        check("export_shape", 0);
        check("export_row_labels", 0);
        check("export_worst_matches_rows", 0);
        check("export_matches_sim", 0);
        check("export_manual_consistent", 0);
    }

    /* heat map vs numbers must agree exactly; a tampered map must fail */
    check("heatmap_matches_numbers", dbg_verify(dir) == 0);

    if (parsed) {
        char ppm_path[1024];
        snprintf(ppm_path, sizeof(ppm_path), "%s/heatmap.ppm", dir);
        Ppm ppm;
        if (ppm_read(ppm_path, &ppm) == 0) {
            /* paint one cell a different colour; red/green stay a coherent
               ramp pair so only the value comparison can catch it */
            size_t off = ((size_t)(ppm.h / 2) * (size_t)ppm.w + (size_t)(ppm.w / 2)) * 3u;
            ppm.data[off] = (unsigned char)((ppm.data[off] + 40) % 256);
            ppm.data[off + 1] = (unsigned char)(180.0f * (float)(255 - ppm.data[off]) / 255.0f + 0.5f);
            FILE* f = fopen(ppm_path, "wb");
            if (f) {
                fprintf(f, "P6\n%d %d\n255\n", ppm.w, ppm.h);
                fwrite(ppm.data, 1, (size_t)ppm.w * (size_t)ppm.h * 3u, f);
                fclose(f);
            }
            check("heatmap_mismatch_detected", dbg_verify(dir) > 0);
            free(ppm.data);
        } else {
            check("heatmap_mismatch_detected", 0);
        }
    } else {
        check("heatmap_mismatch_detected", 0);
    }

    free(dump);
    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
