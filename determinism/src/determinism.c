
/* determinism: reproducible PBD rope solver plus exact run-to-run comparison.

   The solver is deliberately single-threaded, and both orders that matter for
   floating-point reproducibility are pinned down:

     - iteration order: a sweep is either forward (i = 0..n-2) or backward
       (i = n-2..0); DET_ORDER_ALT alternates on the iteration index.  The
       order never depends on scheduling.
     - summation order: diagnostics accumulate left-to-right.  There is no
       tree reduction, no atomics and no per-thread partial sum.
     - parallelism > 1 is rejected with DET_ERR_PARALLEL before touching any
       state, so changing the thread count cannot silently change rounding.

   det_trace_* captures the raw IEEE-754 bits of every coordinate per frame,
   so two runs -- including two separate processes -- can be compared exactly.
   det_sum_chunked/det_values_diff exist only to *quantify* what a parallel
   (chunked) reduction would change; the solver itself never uses them. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <unistd.h>
#include <sys/wait.h>

#define DET_MAX_P 32

typedef enum {
    DET_ORDER_FWD = 0,
    DET_ORDER_BWD = 1,
    DET_ORDER_ALT = 2,          /* fwd, bwd, fwd, ... on the iteration index */
    DET_ORDER__COUNT
} DetOrder;

typedef enum {
    DET_OK = 0,
    DET_ERR_ARG = -1,
    DET_ERR_PARALLEL = -2,      /* parallelism != 1 */
    DET_ERR_ORDER = -3,         /* scan plan was tampered with mid-run */
    DET_ERR_IO = -4
} DetStatus;

typedef struct {
    int n;                      /* particles, 2..DET_MAX_P */
    double seg;                 /* rest length, > 0 */
    int iterations;             /* sweeps per substep, >= 1 */
    int substeps;               /* integration substeps per frame, >= 1 */
    int parallelism;            /* must be exactly 1 */
    DetOrder order;
    double dt;                  /* frame duration, > 0 and finite */
} DetConfig;

typedef struct {
    DetConfig cfg;
    int n;
    double x[DET_MAX_P], y[DET_MAX_P];
    double px[DET_MAX_P], py[DET_MAX_P];
    double x0[DET_MAX_P], y0[DET_MAX_P];
    double rest[DET_MAX_P - 1];
    double w[DET_MAX_P];        /* inverse mass; 0 pins the particle */
    int frame;
    uint64_t declared_order_id; /* fingerprint of the plan fixed at init */
} DetRope;

/* ---------------- order fingerprint ---------------- */

static uint64_t fnv1a(uint64_t h, const void* p, size_t n) {
    const unsigned char* b = (const unsigned char*)p;
    for (size_t i = 0; i < n; i++) {
        h ^= (uint64_t)b[i];
        h *= 1099511628211ull;
    }
    return h;
}

static uint64_t order_id(const DetConfig* c) {
    uint64_t h = 1469598103934665603ull;
    int vals[4];
    vals[0] = (int)c->order;
    vals[1] = c->iterations;
    vals[2] = c->substeps;
    vals[3] = c->n;
    h = fnv1a(h, vals, sizeof(vals));
    h = fnv1a(h, &c->seg, sizeof(c->seg));
    h = fnv1a(h, &c->dt, sizeof(c->dt));
    return h;
}

/* ---------------- solver ---------------- */

static int order_id_valid(DetOrder o) { return (int)o >= 0 && o < DET_ORDER__COUNT; }

int det_rope_init(DetRope* r, const DetConfig* cfg) {
    if (!r || !cfg) return DET_ERR_ARG;
    if (cfg->n < 2 || cfg->n > DET_MAX_P) return DET_ERR_ARG;
    if (!isfinite(cfg->seg) || cfg->seg <= 0.0) return DET_ERR_ARG;
    if (!isfinite(cfg->dt) || cfg->dt <= 0.0) return DET_ERR_ARG;
    if (cfg->iterations < 1 || cfg->substeps < 1) return DET_ERR_ARG;
    if (!order_id_valid(cfg->order)) return DET_ERR_ARG;
    if (cfg->parallelism < 1) return DET_ERR_ARG;
    if (cfg->parallelism != 1) return DET_ERR_PARALLEL;  /* no parallel reduction */

    memset(r, 0, sizeof(*r));
    r->cfg = *cfg;
    r->n = cfg->n;
    for (int i = 0; i < r->n; i++) {
        r->x[i] = cfg->seg * (double)i;
        r->y[i] = 0.0;
        r->px[i] = r->x[i];
        r->py[i] = r->y[i];
        r->x0[i] = r->x[i];
        r->y0[i] = r->y[i];
        r->w[i] = 1.0;
    }
    r->w[0] = 0.0;                       /* pin the anchor */
    for (int i = 0; i < r->n - 1; i++) r->rest[i] = cfg->seg;
    r->frame = 0;
    r->declared_order_id = order_id(&r->cfg);
    return DET_OK;
}

void det_rope_reset(DetRope* r) {
    if (!r) return;
    for (int i = 0; i < r->n; i++) {
        r->x[i] = r->x0[i];
        r->y[i] = r->y0[i];
        r->px[i] = r->x0[i];
        r->py[i] = r->y0[i];
    }
    r->frame = 0;
}

static void integrate(DetRope* r, double h, double gx, double gy) {
    /* fixed index order; the pinned particle is skipped, not reordered */
    for (int i = 0; i < r->n; i++) {
        if (r->w[i] == 0.0) continue;
        double vx = r->x[i] - r->px[i];
        double vy = r->y[i] - r->py[i];
        r->px[i] = r->x[i];
        r->py[i] = r->y[i];
        r->x[i] += vx + gx * h * h;
        r->y[i] += vy + gy * h * h;
    }
}

static void solve_segment(DetRope* r, int i) {
    double dx = r->x[i + 1] - r->x[i];
    double dy = r->y[i + 1] - r->y[i];
    double d2 = dx * dx + dy * dy;
    if (!(d2 > 1e-24)) return;
    double d = sqrt(d2);
    double diff = (d - r->rest[i]) / d;
    double wa = r->w[i], wb = r->w[i + 1];
    double wsum = wa + wb;
    if (!(wsum > 0.0)) return;
    double ka = diff * (wa / wsum);
    double kb = diff * (wb / wsum);
    r->x[i] += dx * ka;
    r->y[i] += dy * ka;
    r->x[i + 1] -= dx * kb;
    r->y[i + 1] -= dy * kb;
}

static void sweep_forward(DetRope* r) {
    for (int i = 0; i < r->n - 1; i++) solve_segment(r, i);
}

static void sweep_backward(DetRope* r) {
    for (int i = r->n - 2; i >= 0; i--) solve_segment(r, i);
}

int det_rope_step(DetRope* r, double gx, double gy) {
    if (!r) return DET_ERR_ARG;
    if (!isfinite(gx) || !isfinite(gy)) return DET_ERR_ARG;
    if (r->cfg.parallelism != 1) return DET_ERR_PARALLEL;
    if (order_id(&r->cfg) != r->declared_order_id) return DET_ERR_ORDER;

    double h = r->cfg.dt / (double)r->cfg.substeps;
    for (int s = 0; s < r->cfg.substeps; s++) {
        integrate(r, h, gx, gy);
        for (int it = 0; it < r->cfg.iterations; it++) {
            int fwd;
            if (r->cfg.order == DET_ORDER_FWD) fwd = 1;
            else if (r->cfg.order == DET_ORDER_BWD) fwd = 0;
            else fwd = (it & 1) == 0;
            if (fwd) sweep_forward(r); else sweep_backward(r);
        }
    }
    r->frame++;
    return DET_OK;
}

/* Fault injection for the self-test: changes the live scan plan but leaves
   declared_order_id alone, so the next det_rope_step must reject it. */
int det_rope_force_order(DetRope* r, DetOrder order) {
    if (!r || !order_id_valid(order)) return DET_ERR_ARG;
    r->cfg.order = order;
    return DET_OK;
}

/* Diagnostics; every sum is strictly left-to-right. */
double det_rope_max_residual(const DetRope* r) {
    double worst = 0.0;
    for (int i = 0; i < r->n - 1; i++) {
        double dx = r->x[i + 1] - r->x[i];
        double dy = r->y[i + 1] - r->y[i];
        double e = fabs(sqrt(dx * dx + dy * dy) - r->rest[i]);
        if (e > worst) worst = e;
    }
    return worst;
}

void det_rope_seg_residuals(const DetRope* r, double* out) {
    if (!r || !out) return;
    for (int i = 0; i < r->n - 1; i++) {
        double dx = r->x[i + 1] - r->x[i];
        double dy = r->y[i + 1] - r->y[i];
        out[i] = fabs(sqrt(dx * dx + dy * dy) - r->rest[i]);
    }
}

/* ---------------- exact comparison ---------------- */

typedef struct {
    long long bit_differences;  /* coordinate words whose 64 bits differ */
    int coordinate_count;
    int first_coordinate;       /* 2*i+axis, or -1 */
    double max_abs;
    double l2;                  /* sqrt of left-to-right sum of squares */
    int nan_mismatch;
    int shape_mismatch;
    int frame_a, frame_b;
} DetDiff;

static uint64_t dbits(double v) {
    uint64_t b;
    memcpy(&b, &v, sizeof(b));
    return b;
}

void det_rope_diff(const DetRope* a, const DetRope* b, DetDiff* d) {
    memset(d, 0, sizeof(*d));
    d->first_coordinate = -1;
    if (!a || !b) return;
    int n = a->n < b->n ? a->n : b->n;
    d->coordinate_count = 2 * n;
    d->frame_a = a->frame;
    d->frame_b = b->frame;
    d->shape_mismatch = (a->n != b->n);
    double l2 = 0.0;
    for (int i = 0; i < n; i++) {
        double va[2] = { a->x[i], a->y[i] };
        double vb[2] = { b->x[i], b->y[i] };
        for (int axis = 0; axis < 2; axis++) {
            if (dbits(va[axis]) != dbits(vb[axis])) {
                d->bit_differences++;
                if (d->first_coordinate < 0) d->first_coordinate = 2 * i + axis;
            }
            double diff = va[axis] - vb[axis];
            if (isnan(diff)) { d->nan_mismatch = 1; continue; }
            double ad = fabs(diff);
            if (ad > d->max_abs) d->max_abs = ad;
            l2 += diff * diff;
        }
    }
    d->l2 = sqrt(l2);
}

/* ---------------- per-frame trace (raw bits) ---------------- */

typedef struct {
    int n;
    int nframes;
    int capacity;
    uint64_t* bits;             /* nframes * n * 2 words: x,y per particle */
} DetTrace;

int det_trace_init(DetTrace* t, int n, int capacity) {
    if (!t || n < 1 || n > DET_MAX_P || capacity < 1) return DET_ERR_ARG;
    memset(t, 0, sizeof(*t));
    t->bits = (uint64_t*)calloc((size_t)capacity * (size_t)n * 2, sizeof(uint64_t));
    if (!t->bits) return DET_ERR_ARG;
    t->n = n;
    t->capacity = capacity;
    return DET_OK;
}

void det_trace_free(DetTrace* t) {
    if (!t) return;
    free(t->bits);
    t->bits = NULL;
    t->n = t->nframes = t->capacity = 0;
}

int det_trace_capture(DetTrace* t, const DetRope* r) {
    if (!t || !r || t->n != r->n) return DET_ERR_ARG;
    if (t->nframes >= t->capacity) return DET_ERR_ARG;
    uint64_t* row = t->bits + (size_t)t->nframes * (size_t)t->n * 2;
    for (int i = 0; i < r->n; i++) {
        row[2 * i] = dbits(r->x[i]);
        row[2 * i + 1] = dbits(r->y[i]);
    }
    t->nframes++;
    return DET_OK;
}

typedef struct {
    int frames;                 /* compared frames = min(frames, ...) */
    long long bit_differences;
    int first_frame;            /* -1 if identical */
    int first_coordinate;       /* 2*i+axis, or -1 */
    double max_abs;
    double l2;
    int nan_mismatch;
    int shape_mismatch;
} DetTraceDiff;

void det_trace_diff(const DetTrace* a, const DetTrace* b, DetTraceDiff* d) {
    memset(d, 0, sizeof(*d));
    d->first_frame = -1;
    d->first_coordinate = -1;
    if (!a || !b) { d->shape_mismatch = 1; return; }
    d->shape_mismatch = (a->n != b->n);
    int n = a->n < b->n ? a->n : b->n;
    int frames = a->nframes < b->nframes ? a->nframes : b->nframes;
    if (a->nframes != b->nframes) d->shape_mismatch = 1;
    d->frames = frames;
    double l2 = 0.0;
    for (int f = 0; f < frames; f++) {
        const uint64_t* ra = a->bits + (size_t)f * (size_t)a->n * 2;
        const uint64_t* rb = b->bits + (size_t)f * (size_t)b->n * 2;
        for (int i = 0; i < n; i++) {
            for (int axis = 0; axis < 2; axis++) {
                int k = 2 * i + axis;
                if (ra[k] != rb[k]) {
                    d->bit_differences++;
                    if (d->first_frame < 0) {
                        d->first_frame = f;
                        d->first_coordinate = k;
                    }
                }
                double va, vb;
                memcpy(&va, &ra[k], sizeof(va));
                memcpy(&vb, &rb[k], sizeof(vb));
                double diff = va - vb;
                if (isnan(diff)) { d->nan_mismatch = 1; continue; }
                double ad = fabs(diff);
                if (ad > d->max_abs) d->max_abs = ad;
                l2 += diff * diff;
            }
        }
    }
    d->l2 = sqrt(l2);
}

#define DET_TRACE_MAGIC 0x44545231u  /* "DTR1" */
#define DET_TRACE_VERSION 1u

int det_trace_write(const DetTrace* t, const char* path) {
    if (!t || !path) return DET_ERR_ARG;
    FILE* f = fopen(path, "wb");
    if (!f) return DET_ERR_IO;
    uint32_t magic = DET_TRACE_MAGIC, version = DET_TRACE_VERSION;
    int32_t n = (int32_t)t->n, nframes = (int32_t)t->nframes;
    int ok = fwrite(&magic, sizeof(magic), 1, f) == 1 &&
             fwrite(&version, sizeof(version), 1, f) == 1 &&
             fwrite(&n, sizeof(n), 1, f) == 1 &&
             fwrite(&nframes, sizeof(nframes), 1, f) == 1 &&
             fwrite(t->bits, sizeof(uint64_t),
                    (size_t)t->nframes * (size_t)t->n * 2, f) ==
                 (size_t)t->nframes * (size_t)t->n * 2;
    if (fclose(f) != 0) ok = 0;
    return ok ? DET_OK : DET_ERR_IO;
}

int det_trace_read(DetTrace* t, const char* path) {
    if (!t || !path) return DET_ERR_ARG;
    FILE* f = fopen(path, "rb");
    if (!f) return DET_ERR_IO;
    uint32_t magic = 0, version = 0;
    int32_t n = 0, nframes = 0;
    int ok = fread(&magic, sizeof(magic), 1, f) == 1 &&
             fread(&version, sizeof(version), 1, f) == 1 &&
             fread(&n, sizeof(n), 1, f) == 1 &&
             fread(&nframes, sizeof(nframes), 1, f) == 1;
    ok = ok && magic == DET_TRACE_MAGIC && version == DET_TRACE_VERSION &&
         n >= 1 && n <= DET_MAX_P && nframes >= 1;
    if (ok) {
        if (det_trace_init(t, n, nframes) != DET_OK) ok = 0;
    }
    if (ok) {
        size_t words = (size_t)nframes * (size_t)n * 2;
        ok = fread(t->bits, sizeof(uint64_t), words, f) == words;
        t->nframes = nframes;
    }
    fclose(f);
    return ok ? DET_OK : DET_ERR_IO;
}

/* ---------------- reduction / parallelism quantification ----------------

   These do NOT feed the solver.  det_sum_chunked emulates the deterministic
   version of a segmented parallel sum: workers contiguous chunks, each summed
   left-to-right, then the partials summed left-to-right in chunk order.
   Comparing it with the serial sum quantifies what a change in parallelism
   costs; the value is stable for a fixed (n, workers). */

double det_sum_forward(const double* v, int n) {
    double s = 0.0;
    for (int i = 0; i < n; i++) s += v[i];
    return s;
}

double det_sum_chunked(const double* v, int n, int workers) {
    if (!v || n <= 0) return 0.0;
    if (workers <= 1) return det_sum_forward(v, n);
    if (workers > n) workers = n;
    double partial[DET_MAX_P];
    for (int w = 0; w < workers; w++) {
        int lo = (int)((long long)n * w / workers);
        int hi = (int)((long long)n * (w + 1) / workers);
        partial[w] = det_sum_forward(v + lo, hi - lo);
    }
    return det_sum_forward(partial, workers);
}

typedef struct {
    long long bit_differences;
    int first_index;
    double max_abs;
    double l2;
    int nan_mismatch;
} DetValueDiff;

DetValueDiff det_values_diff(const double* a, const double* b, int n) {
    DetValueDiff d;
    memset(&d, 0, sizeof(d));
    d.first_index = -1;
    double l2 = 0.0;
    for (int i = 0; i < n; i++) {
        if (dbits(a[i]) != dbits(b[i])) {
            d.bit_differences++;
            if (d.first_index < 0) d.first_index = i;
        }
        double diff = a[i] - b[i];
        if (isnan(diff)) { d.nan_mismatch = 1; continue; }
        double ad = fabs(diff);
        if (ad > d.max_abs) d.max_abs = ad;
        l2 += diff * diff;
    }
    d.l2 = sqrt(l2);
    return d;
}

/* ---------------- self-test ---------------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static DetConfig base_cfg(void) {
    DetConfig c;
    memset(&c, 0, sizeof(c));
    c.n = 8;
    c.seg = 1.0;
    c.iterations = 8;
    c.substeps = 1;
    c.parallelism = 1;
    c.order = DET_ORDER_ALT;
    c.dt = 0.016;
    return c;
}

#define TEST_FRAMES 240

static int run_frames(DetRope* r, int frames, double gx, double gy) {
    for (int f = 0; f < frames; f++)
        if (det_rope_step(r, gx, gy) != DET_OK) return -1;
    return 0;
}

static int build_trace(DetTrace* t, const DetConfig* cfg, int frames) {
    DetRope r;
    if (det_rope_init(&r, cfg) != DET_OK) return -1;
    if (det_trace_init(t, cfg->n, frames) != DET_OK) return -1;
    for (int f = 0; f < frames; f++) {
        if (det_rope_step(&r, 0.0, -9.8) != DET_OK) return -1;
        if (det_trace_capture(t, &r) != DET_OK) return -1;
    }
    return 0;
}

/* Simulation state only: config fields (incl. the tampered order) are
   deliberately excluded, and padding is never memcmp'd. */
static int same_positions(const DetRope* a, const DetRope* b) {
    if (a->n != b->n || a->frame != b->frame) return 0;
    return memcmp(a->x, b->x, sizeof(a->x)) == 0 &&
           memcmp(a->y, b->y, sizeof(a->y)) == 0 &&
           memcmp(a->px, b->px, sizeof(a->px)) == 0 &&
           memcmp(a->py, b->py, sizeof(a->py)) == 0;
}
static int trace_child(const char* path) {
    DetConfig cfg = base_cfg();
    DetTrace t;
    if (build_trace(&t, &cfg, TEST_FRAMES) != 0) return 2;
    int rc = det_trace_write(&t, path);
    det_trace_free(&t);
    return rc == DET_OK ? 0 : 3;
}

int main(int argc, char** argv) {
    if (argc == 3 && strcmp(argv[1], "--trace") == 0)
        return trace_child(argv[2]);

    DetConfig cfg = base_cfg();

    /* --- configuration is validated before any state is touched --- */
    DetRope probe;
    DetConfig bad = cfg;
    bad.parallelism = 4;
    check("reject_parallel_init", det_rope_init(&probe, &bad) == DET_ERR_PARALLEL);

    int bad_ok = 1;
    bad = cfg; bad.n = 1;
    bad_ok &= det_rope_init(&probe, &bad) == DET_ERR_ARG;
    bad = cfg; bad.seg = 0.0;
    bad_ok &= det_rope_init(&probe, &bad) == DET_ERR_ARG;
    bad = cfg; bad.seg = NAN;
    bad_ok &= det_rope_init(&probe, &bad) == DET_ERR_ARG;
    bad = cfg; bad.dt = 0.0;
    bad_ok &= det_rope_init(&probe, &bad) == DET_ERR_ARG;
    bad = cfg; bad.iterations = 0;
    bad_ok &= det_rope_init(&probe, &bad) == DET_ERR_ARG;
    bad = cfg; bad.parallelism = 0;
    bad_ok &= det_rope_init(&probe, &bad) == DET_ERR_ARG;
    check("reject_bad_config", bad_ok);

    /* --- two runs, compared bit-for-bit after every single frame --- */
    DetRope a, b;
    check("init_a", det_rope_init(&a, &cfg) == DET_OK);
    check("init_b", det_rope_init(&b, &cfg) == DET_OK);
    DetTrace ta, tb;
    check("trace_a", det_trace_init(&ta, cfg.n, TEST_FRAMES) == DET_OK);
    check("trace_b", det_trace_init(&tb, cfg.n, TEST_FRAMES) == DET_OK);

    long long total_bits = 0;
    double worst_abs = 0.0;
    int per_frame_ok = 1;
    for (int f = 0; f < TEST_FRAMES; f++) {
        int sa = det_rope_step(&a, 0.0, -9.8);
        int sb = det_rope_step(&b, 0.0, -9.8);
        if (sa != DET_OK || sb != DET_OK) { per_frame_ok = 0; break; }
        DetDiff d;
        det_rope_diff(&a, &b, &d);
        if (d.bit_differences != 0) per_frame_ok = 0;
        total_bits += d.bit_differences;
        if (d.max_abs > worst_abs) worst_abs = d.max_abs;
        det_trace_capture(&ta, &a);
        det_trace_capture(&tb, &b);
    }
    check("two_runs_per_frame_bitwise",
          per_frame_ok && total_bits == 0 && worst_abs == 0.0 &&
          a.frame == TEST_FRAMES && b.frame == TEST_FRAMES);

    DetTraceDiff td;
    det_trace_diff(&ta, &tb, &td);
    check("two_runs_trace_bitwise",
          td.bit_differences == 0 && td.max_abs == 0.0 && td.l2 == 0.0 &&
          td.first_frame == -1 && td.frames == TEST_FRAMES);

    /* bitwise comparison is sensitive: reset must reproduce exactly */
    det_rope_reset(&a);
    DetTrace again;
    check("trace_again", det_trace_init(&again, cfg.n, TEST_FRAMES) == DET_OK);
    for (int f = 0; f < TEST_FRAMES; f++) {
        det_rope_step(&a, 0.0, -9.8);
        det_trace_capture(&again, &a);
    }
    det_trace_diff(&ta, &again, &td);
    check("reset_reproduces_bits", td.bit_differences == 0 && a.frame == TEST_FRAMES);

    /* trace file round-trip */
    const char* path_a = "/tmp/mini_anim_det_roundtrip.bin";
    check("trace_write", det_trace_write(&ta, path_a) == DET_OK);
    DetTrace loaded;
    memset(&loaded, 0, sizeof(loaded));
    check("trace_read", det_trace_read(&loaded, path_a) == DET_OK);
    det_trace_diff(&ta, &loaded, &td);
    check("trace_roundtrip_bitwise", td.bit_differences == 0 && td.frames == TEST_FRAMES);
    det_trace_free(&loaded);
    remove(path_a);

    /* --- tampered iteration order must be rejected, state untouched --- */
    DetRope r;
    det_rope_init(&r, &cfg);
    run_frames(&r, 40, 0.0, -9.8);
    DetRope frozen = r;
    check("force_order", det_rope_force_order(&r, DET_ORDER_FWD) == DET_OK);
    check("tampered_order_errors", det_rope_step(&r, 0.0, -9.8) == DET_ERR_ORDER);
    check("tampered_order_no_state", same_positions(&frozen, &r));
    /* rejecting is repeatable, not a one-shot */
    check("tampered_order_errors_again", det_rope_step(&r, 0.0, -9.8) == DET_ERR_ORDER);
    /* putting the declared plan back makes it runnable again */
    det_rope_force_order(&r, DET_ORDER_ALT);
    check("restored_order_runs", det_rope_step(&r, 0.0, -9.8) == DET_OK);

    /* iteration order is part of the result: quantify a legitimate change */
    DetConfig fwd_cfg = base_cfg();
    fwd_cfg.order = DET_ORDER_FWD;
    DetTrace t_fwd, t_alt;
    int alt_ok = build_trace(&t_alt, &cfg, TEST_FRAMES) == 0;
    int fwd_ok = build_trace(&t_fwd, &fwd_cfg, TEST_FRAMES) == 0;
    DetTraceDiff order_d1, order_d2;
    det_trace_diff(&t_alt, &t_fwd, &order_d1);
    {
        DetTrace t_alt2, t_fwd2;
        int ok = build_trace(&t_alt2, &cfg, TEST_FRAMES) == 0 &&
                 build_trace(&t_fwd2, &fwd_cfg, TEST_FRAMES) == 0;
        det_trace_diff(&t_alt2, &t_fwd2, &order_d2);
        check("order_change_quantified",
              alt_ok && fwd_ok && ok && order_d1.bit_differences > 0 &&
              order_d1.max_abs > 0.0 &&
              memcmp(&order_d1, &order_d2, sizeof(order_d1)) == 0);
        det_trace_free(&t_alt2);
        det_trace_free(&t_fwd2);
    }

    /* --- parallelism changes are rejected and quantified, never random --- */
    DetRope p1, p4;
    det_rope_init(&p1, &cfg);
    run_frames(&p1, 20, 0.0, -9.8);
    bad = cfg;
    bad.parallelism = 4;
    check("parallel_run_rejected", det_rope_init(&p4, &bad) == DET_ERR_PARALLEL);
    /* even if the field is poked after init, step refuses and freezes state */
    p4 = p1;
    p4.cfg.parallelism = 4;
    DetRope p4_before = p4;
    check("parallel_step_rejected", det_rope_step(&p4, 0.0, -9.8) == DET_ERR_PARALLEL);
    check("parallel_step_no_state", same_positions(&p4_before, &p4));

    /* serial vs chunked reduction: the delta is exact and reproducible */
    double vals[4] = { 1e16, 1.0, -1e16, 1.0 };
    double fwd = det_sum_forward(vals, 4);
    double c2 = det_sum_chunked(vals, 4, 2);
    double c4 = det_sum_chunked(vals, 4, 4);
    check("reduction_order_detected", fwd == 1.0 && c2 == 0.0 && c4 == 1.0);
    DetValueDiff red;
    {
        double left[4], right[4];
        for (int i = 0; i < 4; i++) {
            left[i] = det_sum_chunked(vals, 4, 2);
            right[i] = det_sum_chunked(vals, 4, 4);
        }
        DetValueDiff red2;
        red = det_values_diff(left, right, 4);
        red2 = det_values_diff(left, right, 4);
        check("parallel_delta_quantified",
              red.bit_differences == 4 && red.max_abs == 1.0 && red.l2 == 2.0 &&
              memcmp(&red, &red2, sizeof(red)) == 0);
    }

    /* the rope's own residual sum: serial == workers=1, and the chunked
       variant is deterministic for each worker count (whatever its value) */
    {
        double res[DET_MAX_P];
        det_rope_seg_residuals(&p1, res);
        double s1 = det_sum_chunked(res, p1.n - 1, 1);
        double s1b = det_sum_forward(res, p1.n - 1);
        double r2 = det_sum_chunked(res, p1.n - 1, 2);
        double r2b = det_sum_chunked(res, p1.n - 1, 2);
        double r4 = det_sum_chunked(res, p1.n - 1, 4);
        double r4b = det_sum_chunked(res, p1.n - 1, 4);
        check("rope_serial_reduction_fixed", dbits(s1) == dbits(s1b) && dbits(r2) == dbits(r2b) &&
              dbits(r4) == dbits(r4b) && s1 == s1b);
    }

    /* --- one-ULP perturbation is detected and its size is exact --- */
    {
        DetRope q1, q2;
        det_rope_init(&q1, &cfg);
        det_rope_init(&q2, &cfg);
        run_frames(&q1, 15, 0.0, -9.8);
        run_frames(&q2, 15, 0.0, -9.8);
        q2.x[4] = nextafter(q2.x[4], INFINITY);
        double ulp = fabs(q2.x[4] - q1.x[4]);
        det_rope_step(&q1, 0.0, -9.8);
        det_rope_step(&q2, 0.0, -9.8);
        DetDiff d;
        det_rope_diff(&q1, &q2, &d);
        check("tiny_perturbation_detected",
              ulp > 0.0 && d.bit_differences > 0 && d.max_abs > 0.0);

        DetRope q3, q4;
        det_rope_init(&q3, &cfg);
        det_rope_init(&q4, &cfg);
        run_frames(&q3, 15, 0.0, -9.8);
        run_frames(&q4, 15, 0.0, -9.8);
        q4.x[4] = q4.x[4] + ulp;               /* same perturbation, again */
        det_rope_step(&q3, 0.0, -9.8);
        det_rope_step(&q4, 0.0, -9.8);
        DetDiff d2, d3;
        det_rope_diff(&q3, &q4, &d2);
        det_rope_diff(&q1, &q2, &d3);
        check("perturbation_delta_reproducible",
              memcmp(&d, &d2, sizeof(d)) == 0 && memcmp(&d, &d3, sizeof(d)) == 0);
    }

    /* --- the same input in two separate processes is bit-identical --- */
    {
        char path_b[256], path_c[256];
        snprintf(path_b, sizeof(path_b), "/tmp/mini_anim_det_%ld_b.bin", (long)getpid());
        snprintf(path_c, sizeof(path_c), "/tmp/mini_anim_det_%ld_c.bin", (long)getpid());
        char cmd[640];
        int rc_b, rc_c;
        snprintf(cmd, sizeof(cmd), "\"%s\" --trace \"%s\"", argv[0], path_b);
        rc_b = system(cmd);
        snprintf(cmd, sizeof(cmd), "\"%s\" --trace \"%s\"", argv[0], path_c);
        rc_c = system(cmd);
        int ok = rc_b != -1 && rc_c != -1 &&
                 WIFEXITED(rc_b) && WEXITSTATUS(rc_b) == 0 &&
                 WIFEXITED(rc_c) && WEXITSTATUS(rc_c) == 0;
        DetTrace proc_b, proc_c;
        memset(&proc_b, 0, sizeof(proc_b));
        memset(&proc_c, 0, sizeof(proc_c));
        ok = ok && det_trace_read(&proc_b, path_b) == DET_OK &&
             det_trace_read(&proc_c, path_c) == DET_OK;
        DetTraceDiff pd;
        det_trace_diff(&proc_b, &proc_c, &pd);
        DetTraceDiff pd2;
        det_trace_diff(&ta, &proc_b, &pd2);
        check("cross_process_bitwise",
              ok && pd.bit_differences == 0 && pd.max_abs == 0.0 &&
              pd2.bit_differences == 0 && pd2.frames == TEST_FRAMES);
        det_trace_free(&proc_b);
        det_trace_free(&proc_c);
        remove(path_b);
        remove(path_c);
    }

    det_trace_free(&ta);
    det_trace_free(&tb);
    det_trace_free(&t_alt);
    det_trace_free(&t_fwd);
    det_trace_free(&again);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
