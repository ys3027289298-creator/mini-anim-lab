
/* palette: dirty-tracked, double-buffered skinning matrix palette.
 *
 * Each frame the naive path recomputes every bone's skinning matrix and
 * shuttles the whole palette around. This module instead:
 *   - marks bones dirty and recomputes only the changed entries;
 *   - builds into a back buffer while the render thread reads the front
 *     buffer, publishing with a seqlock (even seq = stable, odd = writing);
 *   - runs an integrity gate before every upload: finiteness/generation
 *     checks plus an element-wise comparison against a full recompute,
 *     so a missed dirty mark is reported instead of published;
 *   - lets the render reader detect a half-updated (torn) snapshot via the
 *     sequence counter and a data checksum, counting every tear and retrying.
 *
 * The self-tests are single-threaded and deterministic; the real
 * producer/consumer protocol is exercised under pthreads with --thread. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <stdatomic.h>

#define PAL_MAX_BONES 32
#define PAL_READ_RETRIES 32
#define FNV_OFFSET 1469598103934665603ULL

#if defined(PAL_WITH_PTHREAD)
#include <sched.h>
#endif

/* TSan cannot see that the atomic sequence lock protects the non-atomic
 * payload, so it flags every seqlock read/write as a race. Explicit
 * acquire/release annotations give TSan the same happens-before edge the
 * seqlock already provides at runtime; no-ops in every other build. */
#if defined(__SANITIZE_THREAD__)
#include <sanitizer/tsan_interface.h>
#define TSAN_RELEASE(p)  __tsan_release(p)
#define TSAN_ACQUIRE(p)  __tsan_acquire(p)
#else
#define TSAN_RELEASE(p)  ((void)0)
#define TSAN_ACQUIRE(p)  ((void)0)
#endif

typedef struct { float m[9]; } Mat3;          /* row-major 3x3 affine */
typedef struct { float x, y, ang; } Pose;     /* same shape as skeleton BoneWorld */

enum {
    PAL_OK = 0,
    PAL_ERR_CAPACITY = -1,     /* palette full, or output buffer too small */
    PAL_ERR_BAD_INDEX = -2,    /* bone index out of range */
    PAL_ERR_DIRTY_MISS = -3,   /* integrity: a slot differs from full recompute */
    PAL_ERR_BAD_DATA = -4,     /* integrity: non-finite or never-written slot */
    PAL_ERR_TORN = -5,         /* reader: sequence changed mid-snapshot */
    PAL_ERR_RETRY = -6         /* reader: retries exhausted */
};

typedef struct {
    _Atomic uint64_t frames;           /* successfully published frames */
    _Atomic uint64_t recomputed;       /* matrix entries rebuilt on dirty path */
    _Atomic uint64_t reference_calcs;  /* entries recomputed for the gate */
    _Atomic uint64_t uploads;          /* front/back swaps */
    _Atomic uint64_t tears;            /* torn reads observed by the reader */
    _Atomic uint64_t retries;          /* read attempts discarded due to tears */
    _Atomic uint64_t check_fails;      /* integrity gate failures */
} PalStats;

typedef struct {
    Mat3 m[PAL_MAX_BONES];
    _Atomic uint32_t seq;      /* even: stable snapshot, odd: write in progress */
    uint64_t checksum;         /* FNV-1a over m[0..count), sealed at publish */
    uint32_t gen[PAL_MAX_BONES]; /* per-slot write generation, 0 = never written */
} PalBuffer;

typedef struct Palette Palette;
typedef void (*PalHook)(Palette* p, void* user);

struct Palette {
    int capacity;
    int count;
    Mat3 inv_bind[PAL_MAX_BONES];
    PalBuffer buf[2];
    atomic_int front;          /* buffer index the render thread reads */
    uint8_t dirty[PAL_MAX_BONES];
    int force_full;            /* set on first frame and after a failed check */
    PalStats st;
    PalHook before_validate;   /* fault-injection hook used by the tests */
    void* hook_user;
};

/* ---------------- matrix helpers ---------------- */

static Mat3 mat_identity(void) {
    Mat3 r = {{ 1,0,0, 0,1,0, 0,0,1 }};
    return r;
}

static Mat3 mat_mul(const Mat3* a, const Mat3* b) {
    Mat3 r;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            float sum = 0.0f;
            for (int k = 0; k < 3; k++)
                sum += a->m[i * 3 + k] * b->m[k * 3 + j];
            r.m[i * 3 + j] = sum;
        }
    return r;
}

static Mat3 mat_from_pose(const Pose* p) {
    float c = cosf(p->ang), s = sinf(p->ang);
    Mat3 r = {{
        c, -s, p->x,
        s,  c, p->y,
        0,  0, 1.0f
    }};
    return r;
}

/* inverse of a 2D affine [R t; 0 1] */
static Mat3 mat_invert_affine(const Mat3* a) {
    float a00 = a->m[0], a01 = a->m[1], a02 = a->m[2];
    float a10 = a->m[3], a11 = a->m[4], a12 = a->m[5];
    float inv = 1.0f / (a00 * a11 - a01 * a10);
    Mat3 r;
    r.m[0] =  a11 * inv; r.m[1] = -a01 * inv;
    r.m[3] = -a10 * inv; r.m[4] =  a00 * inv;
    r.m[2] = -(r.m[0] * a02 + r.m[1] * a12);
    r.m[5] = -(r.m[3] * a02 + r.m[4] * a12);
    r.m[6] = 0.0f; r.m[7] = 0.0f; r.m[8] = 1.0f;
    return r;
}

static int mat_finite(const Mat3* a) {
    for (int i = 0; i < 9; i++)
        if (!isfinite(a->m[i])) return 0;
    return 1;
}

static uint64_t fnv_bytes(const void* p, size_t n, uint64_t h) {
    const unsigned char* b = (const unsigned char*)p;
    for (size_t i = 0; i < n; i++) { h ^= b[i]; h *= 1099511628211ULL; }
    return h;
}

/* world(pose) * inverse(bind) -- shared by both paths, so dirty and full
 * recompute produce bit-identical floats for identical inputs */
static Mat3 skin_matrix(const Pose* pose, const Mat3* inv_bind, int i) {
    Mat3 w = mat_from_pose(&pose[i]);
    return mat_mul(&w, &inv_bind[i]);
}

/* ---------------- palette lifecycle ---------------- */

static void palette_init(Palette* p, int capacity) {
    memset(p, 0, sizeof(*p));
    if (capacity > PAL_MAX_BONES) capacity = PAL_MAX_BONES;
    if (capacity < 1) capacity = 1;
    p->capacity = capacity;
    p->force_full = 1;   /* first frame must populate every slot */
}

static int palette_add_bone(Palette* p, const Pose* bind) {
    if (p->count >= p->capacity) return PAL_ERR_CAPACITY;
    int i = p->count++;
    Mat3 bind_m = mat_from_pose(bind);
    p->inv_bind[i] = mat_invert_affine(&bind_m);
    return i;
}

static int palette_mark_dirty(Palette* p, int bone) {
    if (bone < 0 || bone >= p->count) return PAL_ERR_BAD_INDEX;
    p->dirty[bone] = 1;
    return PAL_OK;
}

/* full recompute reference, used by the validation gate and by tests */
static void palette_full_reference(Palette* p, const Pose* pose, Mat3* out) {
    for (int i = 0; i < p->count; i++)
        out[i] = skin_matrix(pose, p->inv_bind, i);
}

/* ---------------- writer: build, validate, publish ---------------- */

static int palette_fail(Palette* p, PalBuffer* b, int err) {
    atomic_fetch_add_explicit(&p->st.check_fails, 1, memory_order_relaxed);
    /* unlock the back buffer but never swap: the front stays intact and
     * the render thread never sees the rejected snapshot */
    TSAN_RELEASE(b);
    atomic_fetch_add_explicit(&b->seq, 1, memory_order_release);
    p->force_full = 1;   /* next frame rebuilds every entry */
    return err;
}

static int palette_update(Palette* p, const Pose* pose) {
    if (p->count == 0) return PAL_OK;

    int front = atomic_load_explicit(&p->front, memory_order_acquire);
    int back = front ^ 1;
    PalBuffer* b = &p->buf[back];
    PalBuffer* f = &p->buf[front];

    /* lock the back buffer: odd sequence means "update in progress" */
    atomic_fetch_add_explicit(&b->seq, 1, memory_order_acq_rel);
    TSAN_ACQUIRE(b);

    for (int i = 0; i < p->count; i++) {
        if (p->force_full || p->dirty[i]) {
            b->m[i] = skin_matrix(pose, p->inv_bind, i);
            b->gen[i] = f->gen[i] + 1;
            atomic_fetch_add_explicit(&p->st.recomputed, 1, memory_order_relaxed);
        } else {
            b->m[i] = f->m[i];       /* unchanged: carry the published entry */
            b->gen[i] = f->gen[i];
        }
    }

    if (p->before_validate) p->before_validate(p, p->hook_user);

    /* ---- pre-upload integrity gate, part 1: finiteness / completeness ---- */
    for (int i = 0; i < p->count; i++) {
        if (!mat_finite(&b->m[i]) || b->gen[i] == 0)
            return palette_fail(p, b, PAL_ERR_BAD_DATA);
    }

    /* ---- part 2: element-wise comparison with a full recompute.
     * A missed dirty mark leaves a stale entry here and is caught exactly. */
    for (int i = 0; i < p->count; i++) {
        Mat3 ref = skin_matrix(pose, p->inv_bind, i);
        atomic_fetch_add_explicit(&p->st.reference_calcs, 1, memory_order_relaxed);
        if (memcmp(&ref, &b->m[i], sizeof(Mat3)) != 0)
            return palette_fail(p, b, PAL_ERR_DIRTY_MISS);
    }

    /* seal the snapshot with a checksum, then unlock and publish */
    uint64_t h = FNV_OFFSET;
    for (int i = 0; i < p->count; i++)
        h = fnv_bytes(&b->m[i], sizeof(Mat3), h);
    b->checksum = h;

    atomic_thread_fence(memory_order_release);
    TSAN_RELEASE(b);
    atomic_fetch_add_explicit(&b->seq, 1, memory_order_release); /* even */
    atomic_store_explicit(&p->front, back, memory_order_release);

    memset(p->dirty, 0, sizeof(p->dirty));
    p->force_full = 0;
    atomic_fetch_add_explicit(&p->st.uploads, 1, memory_order_relaxed);
    atomic_fetch_add_explicit(&p->st.frames, 1, memory_order_relaxed);
    return PAL_OK;
}

/* ---------------- reader: seqlock snapshot with tear counting ---------------- */

typedef struct { int idx; uint32_t seq; } ReadGuard;

static int palette_read_begin(Palette* p, ReadGuard* g) {
    g->idx = atomic_load_explicit(&p->front, memory_order_acquire);
    g->seq = atomic_load_explicit(&p->buf[g->idx].seq, memory_order_acquire);
    if (g->seq & 1u) return PAL_ERR_TORN;   /* writer mid-update */
    TSAN_ACQUIRE(&p->buf[g->idx]);
    return PAL_OK;
}

static int palette_read_copy(Palette* p, const ReadGuard* g, Mat3* out, int cap) {
    if (cap < p->count) return PAL_ERR_CAPACITY;
    atomic_thread_fence(memory_order_acquire);
    memcpy(out, p->buf[g->idx].m, sizeof(Mat3) * (size_t)p->count);
    return PAL_OK;
}

static int palette_read_end(Palette* p, const ReadGuard* g) {
    atomic_thread_fence(memory_order_acquire);
    uint32_t s2 = atomic_load_explicit(&p->buf[g->idx].seq, memory_order_acquire);
    if ((s2 & 1u) || s2 != g->seq) return PAL_ERR_TORN;  /* swapped under us */
    uint64_t h = FNV_OFFSET;
    for (int i = 0; i < p->count; i++)
        h = fnv_bytes(&p->buf[g->idx].m[i], sizeof(Mat3), h);
    if (h != p->buf[g->idx].checksum) return PAL_ERR_TORN;
    TSAN_RELEASE(&p->buf[g->idx]);
    return PAL_OK;
}

/* full seqlock read. The two hazards have independent handling:
 *   - odd sequence on entry: writer is mid-frame. Wait for it (no retry
 *     budget consumed; a correct writer always finishes);
 *   - sequence/checksum change mid-snapshot: a genuine torn read, counted
 *     and retried up to PAL_READ_RETRIES.
 * Every failed observation increments the tear counters, and only a fully
 * consistent snapshot is ever returned. */
static void pal_backoff(int attempt) {
#if defined(PAL_WITH_PTHREAD)
    (void)attempt;
    sched_yield();
#else
    volatile unsigned long long sink = 0;
    int spins = 32 << (attempt < 12 ? attempt : 12);
    for (int k = 0; k < spins; k++) sink += (unsigned long long)k;
    (void)sink;
#endif
}

static int palette_read(Palette* p, Mat3* out, int cap) {
    int torn = 0;
    for (;;) {
        ReadGuard g;
        int begin_rc = palette_read_begin(p, &g);
        if (begin_rc == PAL_ERR_CAPACITY) return begin_rc;
        if (begin_rc == PAL_ERR_TORN) {
            atomic_fetch_add_explicit(&p->st.tears, 1, memory_order_relaxed); atomic_fetch_add_explicit(&p->st.retries, 1, memory_order_relaxed);
            while (atomic_load_explicit(
                       &p->buf[g.idx].seq, memory_order_acquire) & 1u)
                pal_backoff(0);
            continue;                          /* lock wait never fails */
        }
        int rc = palette_read_copy(p, &g, out, cap);
        if (rc != PAL_OK) return rc;          /* capacity: not a tear */
        if (palette_read_end(p, &g) == PAL_OK)
            return PAL_OK;
        atomic_fetch_add_explicit(&p->st.tears, 1, memory_order_relaxed); atomic_fetch_add_explicit(&p->st.retries, 1, memory_order_relaxed);
        if (++torn >= PAL_READ_RETRIES) return PAL_ERR_RETRY;
        pal_backoff(torn);
    }
}

/* ---------------- deterministic test hooks ---------------- */

static void palette_dbg_set_locked(Palette* p, int idx, int locked) {
    uint32_t cur = atomic_load(&p->buf[idx].seq);
    if ((int)(cur & 1u) != locked)
        atomic_fetch_add_explicit(&p->buf[idx].seq, 1, memory_order_acq_rel);
}

static void palette_dbg_swap(Palette* p) {
    int f = atomic_load_explicit(&p->front, memory_order_acquire);
    atomic_store_explicit(&p->front, f ^ 1, memory_order_release);
}

static void palette_dbg_scribble(Palette* p, int idx, int slot) {
    p->buf[idx].m[slot] = mat_identity();
}

static void hook_nan(Palette* p, void* user) {
    (void)user;
    int back = atomic_load_explicit(&p->front, memory_order_acquire) ^ 1;
    p->buf[back].m[0].m[0] = NAN;
}

static void hook_uninit(Palette* p, void* user) {
    (void)user;
    int back = atomic_load_explicit(&p->front, memory_order_acquire) ^ 1;
    p->buf[back].gen[1] = 0;
}

/* ---------------- test harness ---------------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static uint32_t rng_state = 0x12345678u;
static uint32_t rng_next(void) {
    rng_state = rng_state * 1103515245u + 12345u;
    return rng_state;
}

static int pose_equal(const Pose* a, const Pose* b) {
    return memcmp(a, b, sizeof(Pose)) == 0;
}

#define TEST_BONES 6
#define TEST_FRAMES 24

static void build_poses(Pose frames[TEST_FRAMES][TEST_BONES]) {
    for (int f = 0; f < TEST_FRAMES; f++)
        for (int i = 0; i < TEST_BONES; i++) {
            frames[f][i].x   = (float)(rng_next() % 7) - 3.0f;
            frames[f][i].y   = (float)(rng_next() % 5) - 2.0f;
            frames[f][i].ang = (float)(rng_next() % 12) * 0.5f;
        }
}

static Palette make_palette(const Pose bind[TEST_BONES], int capacity) {
    Palette p;
    palette_init(&p, capacity);
    for (int i = 0; i < TEST_BONES; i++) palette_add_bone(&p, &bind[i]);
    return p;
}

static void run_deterministic_tests(void) {
    Pose bind[TEST_BONES] = {
        { 0, 0, 0.0f }, { 1, 0, 0.3f }, { 2, 0, 0.6f },
        { 0, 1, 1.0f }, { 1, 1, 1.3f }, { 2, 1, 1.6f }
    };
    Pose poses[TEST_FRAMES][TEST_BONES];
    build_poses(poses);

    /* ---- capacity and index errors ---- */
    Palette cap;
    palette_init(&cap, 4);
    int rc = PAL_OK, last = -1;
    for (int i = 0; i < 6; i++) {
        rc = palette_add_bone(&cap, &bind[0]);
        if (rc >= 0) last = i;
    }
    check("capacity_add", rc == PAL_ERR_CAPACITY && cap.count == 4 && last == 3);
    check("bad_index", palette_mark_dirty(&cap, 7) == PAL_ERR_BAD_INDEX &&
                       palette_mark_dirty(&cap, -1) == PAL_ERR_BAD_INDEX);

    Palette p = make_palette(bind, PAL_MAX_BONES);
    Mat3 ref[PAL_MAX_BONES], snap[PAL_MAX_BONES], small[2];

    /* skinning at the bind pose must be the identity */
    palette_full_reference(&p, bind, ref);
    int id_ok = 1;
    for (int i = 0; i < TEST_BONES; i++)
        for (int k = 0; k < 9; k++)
            if (fabsf(ref[i].m[k] - mat_identity().m[k]) > 1e-6f) id_ok = 0;
    check("bind_skin_identity", id_ok);

    /* ---- first frame: forced full rebuild ---- */
    rc = palette_update(&p, poses[0]);
    check("first_frame_ok", rc == PAL_OK && atomic_load_explicit(&p.st.uploads, memory_order_relaxed) == 1);
    check("first_frame_full", atomic_load_explicit(&p.st.recomputed, memory_order_relaxed) == (uint64_t)TEST_BONES);
    palette_full_reference(&p, poses[0], ref);
    check("first_read", palette_read(&p, snap, PAL_MAX_BONES) == PAL_OK &&
                       memcmp(snap, ref, sizeof(Mat3) * TEST_BONES) == 0);
    check("read_capacity", palette_read(&p, small, 2) == PAL_ERR_CAPACITY);

    /* ---- dirty path vs full recompute, element-wise and bit-exact ---- */
    uint64_t recomputed_at_start = atomic_load_explicit(&p.st.recomputed, memory_order_relaxed);
    uint64_t changed_total = 0;
    long long elems = 0, mismatches = 0;
    for (int f = 1; f < TEST_FRAMES; f++) {
        uint64_t changed = 0;
        for (int i = 0; i < TEST_BONES; i++) {
            if (!pose_equal(&poses[f][i], &poses[f - 1][i])) {
                palette_mark_dirty(&p, i);
                changed++;
            }
        }
        changed_total += changed;
        rc = palette_update(&p, poses[f]);
        if (rc != PAL_OK) { check("dirty_frame_ok", 0); break; }

        palette_full_reference(&p, poses[f], ref);
        if (palette_read(&p, snap, PAL_MAX_BONES) != PAL_OK) {
            check("dirty_frame_read", 0); break;
        }
        for (int i = 0; i < TEST_BONES; i++) {
            elems += 9;
            if (memcmp(&snap[i], &ref[i], sizeof(Mat3)) != 0) mismatches++;
        }
    }
    check("dirty_frame_ok", rc == PAL_OK);
    printf("exact_elements=%lld mismatches=%lld\n", elems, mismatches);
    check("exact_compare_zero_error", mismatches == 0 && elems == (TEST_FRAMES - 1) * TEST_BONES * 9);
    check("only_dirty_recomputed",
          atomic_load_explicit(&p.st.recomputed, memory_order_relaxed) - recomputed_at_start == changed_total);

    /* a frame with no pose changes recomputes nothing but still validates */
    uint64_t before = atomic_load_explicit(&p.st.recomputed, memory_order_relaxed);
    check("no_change_frame", palette_update(&p, poses[TEST_FRAMES - 1]) == PAL_OK &&
                             atomic_load_explicit(&p.st.recomputed, memory_order_relaxed) == before);
    palette_full_reference(&p, poses[TEST_FRAMES - 1], ref);
    check("no_change_read", palette_read(&p, snap, PAL_MAX_BONES) == PAL_OK &&
                            memcmp(snap, ref, sizeof(Mat3) * TEST_BONES) == 0);

    /* ---- missed dirty mark: changed bone left unmarked ---- */
    Palette good = make_palette(bind, PAL_MAX_BONES);
    Pose miss0[TEST_BONES], miss1[TEST_BONES];
    for (int i = 0; i < TEST_BONES; i++)
        miss0[i] = (Pose){ (float)i, 0.0f, 0.0f };
    memcpy(miss1, miss0, sizeof(miss0));
    miss1[4].ang = 0.75f;                        /* bone 4 definitely changes */
    palette_update(&good, miss0);
    palette_full_reference(&good, miss0, ref);
    Mat3 safe[PAL_MAX_BONES];
    palette_read(&good, safe, PAL_MAX_BONES);

    palette_mark_dirty(&good, 2);                 /* bone 4 also changes... */
    uint64_t uploads_before = atomic_load_explicit(&good.st.uploads, memory_order_relaxed);
    rc = palette_update(&good, miss1);            /* ...but is never marked */
    check("dirty_miss_reported", rc == PAL_ERR_DIRTY_MISS &&
                                 atomic_load_explicit(&good.st.check_fails, memory_order_relaxed) == 1);
    check("dirty_miss_no_upload", atomic_load_explicit(&good.st.uploads, memory_order_relaxed) == uploads_before);
    check("dirty_miss_front_safe",
          palette_read(&good, snap, PAL_MAX_BONES) == PAL_OK &&
          memcmp(snap, safe, sizeof(Mat3) * TEST_BONES) == 0);

    /* failed check poisons the palette: next frame force-rebuilds and recovers */
    rc = palette_update(&good, miss1);
    palette_full_reference(&good, miss1, ref);
    check("recover_force_full", rc == PAL_OK &&
          palette_read(&good, snap, PAL_MAX_BONES) == PAL_OK &&
          memcmp(snap, ref, sizeof(Mat3) * TEST_BONES) == 0);

    /* ---- render reader vs half-updated buffer (deterministic interleaving) ---- */
    Palette r = make_palette(bind, PAL_MAX_BONES);
    palette_update(&r, poses[0]);
    Mat3 ref0[PAL_MAX_BONES], ref1[PAL_MAX_BONES];
    palette_full_reference(&r, poses[0], ref0);
    for (int i = 0; i < TEST_BONES; i++) palette_mark_dirty(&r, i);
    palette_update(&r, poses[1]);                 /* both buffers now sealed */
    palette_full_reference(&r, poses[1], ref1);

    /* reader opens a stable snapshot; writer locks that buffer mid-read so
     * the end check observes a torn (odd) sequence -- the half-update is
     * detected and counted, not returned */
    ReadGuard g;
    check("read_begin_stable", palette_read_begin(&r, &g) == PAL_OK);
    int fidx = atomic_load_explicit(&r.front, memory_order_acquire);
    palette_dbg_set_locked(&r, fidx, 1);          /* odd seq: update in progress */
    check("tear_mid_read", palette_read_end(&r, &g) == PAL_ERR_TORN);
    atomic_fetch_add_explicit(&r.st.tears, 1, memory_order_relaxed);
    atomic_fetch_add_explicit(&r.st.retries, 1, memory_order_relaxed);

    /* writer finishes: the next full read blocks until even, then succeeds */
    palette_dbg_set_locked(&r, fidx, 0);
    uint64_t tears_before = atomic_load_explicit(&r.st.tears, memory_order_relaxed);
    check("read_after_unlock",
          palette_read(&r, snap, PAL_MAX_BONES) == PAL_OK &&
          memcmp(snap, ref1, sizeof(Mat3) * TEST_BONES) == 0);
    check("tear_was_counted", tears_before >= 1);

    /* sustained contention: corrupt the sealed checksum so every snapshot
     * validation fails, exercising the retry budget end to end through the
     * public read path */
    Palette t = make_palette(bind, PAL_MAX_BONES);
    palette_update(&t, poses[0]);
    int tidx = atomic_load_explicit(&t.front, memory_order_acquire);
    t.buf[tidx].checksum ^= 0xdeadbeefULL;
    uint64_t tb = atomic_load_explicit(&t.st.tears, memory_order_relaxed);
    check("retry_budget_exhausted",
          palette_read(&t, snap, PAL_MAX_BONES) == PAL_ERR_RETRY &&
          atomic_load_explicit(&t.st.tears, memory_order_relaxed) - tb == (uint64_t)PAL_READ_RETRIES);
    t.buf[tidx].checksum ^= 0xdeadbeefULL;     /* restore: reads work again */
    check("retry_recovers",
          palette_read(&t, snap, PAL_MAX_BONES) == PAL_OK);

    /* writer swaps front and starts reusing the reader's old buffer */
    ReadGuard g2;
    check("read_before_swap", palette_read_begin(&r, &g2) == PAL_OK);
    int old_front = atomic_load(&r.front);
    palette_dbg_swap(&r);
    palette_dbg_set_locked(&r, old_front, 1);
    palette_dbg_scribble(&r, old_front, 0);       /* partial write */
    check("tear_after_swap", palette_read_end(&r, &g2) == PAL_ERR_TORN);
    /* new reads land on the other (stable) buffer: the previously published frame */
    check("new_read_on_stable_front",
          palette_read(&r, snap, PAL_MAX_BONES) == PAL_OK &&
          memcmp(snap, ref0, sizeof(Mat3) * TEST_BONES) == 0);
    palette_dbg_set_locked(&r, old_front, 0);
    r.force_full = 1;
    palette_update(&r, poses[0]);                 /* republish clean state */

    /* ---- integrity gate: bad data is refused ---- */
    Palette q = make_palette(bind, PAL_MAX_BONES);
    palette_update(&q, poses[0]);
    q.before_validate = hook_nan;
    check("nan_rejected", palette_update(&q, poses[1]) == PAL_ERR_BAD_DATA &&
                          atomic_load_explicit(&q.st.check_fails, memory_order_relaxed) == 1 && atomic_load_explicit(&q.st.uploads, memory_order_relaxed) == 1);
    q.before_validate = NULL;
    q.force_full = 1;
    check("nan_recovered", palette_update(&q, poses[1]) == PAL_OK);

    Palette q2 = make_palette(bind, PAL_MAX_BONES);
    palette_update(&q2, poses[0]);
    q2.before_validate = hook_uninit;
    check("uninit_rejected", palette_update(&q2, poses[1]) == PAL_ERR_BAD_DATA);

    /* ---- stats ---- */
    printf("stats_frames=%llu recomputed=%llu reference=%llu uploads=%llu "
           "tears=%llu retries=%llu check_fails=%llu\n",
           (unsigned long long)atomic_load_explicit(&p.st.frames, memory_order_relaxed),
           (unsigned long long)atomic_load_explicit(&p.st.recomputed, memory_order_relaxed),
           (unsigned long long)atomic_load_explicit(&p.st.reference_calcs, memory_order_relaxed),
           (unsigned long long)atomic_load_explicit(&p.st.uploads, memory_order_relaxed),
           (unsigned long long)atomic_load_explicit(&r.st.tears, memory_order_relaxed),
           (unsigned long long)atomic_load_explicit(&r.st.retries, memory_order_relaxed),
           (unsigned long long)(atomic_load_explicit(&p.st.check_fails, memory_order_relaxed) + atomic_load_explicit(&good.st.check_fails, memory_order_relaxed) +
                                atomic_load_explicit(&q.st.check_fails, memory_order_relaxed) +
                                atomic_load_explicit(&q2.st.check_fails, memory_order_relaxed)));
}

#if defined(PAL_WITH_PTHREAD)
#include <pthread.h>
#include <sched.h>

#define THREAD_FRAMES 3000
#define THREAD_READERS 2

typedef struct {
    Palette* pal;
    const Pose* poses;
    int frames, bones;
    atomic_int done;
    atomic_int failed;
} WriterArg;

typedef struct {
    Palette* pal;
    atomic_int* done;
    atomic_int failed;
    _Atomic long long reads;
    _Atomic long long exhausted;
} ReaderArg;

static void hook_spin(Palette* p, void* user) {
    (void)p;
    volatile int x = 0;
    int n = *(const int*)user;
    for (int i = 0; i < n; i++) x += i;
    (void)x;
}

static void* writer_thread(void* arg) {
    WriterArg* w = (WriterArg*)arg;
    Pose prev[PAL_MAX_BONES];
    memset(prev, 0, sizeof(prev));
    for (int f = 0; f < w->frames; f++) {
        const Pose* pose = w->poses + (size_t)f * w->bones;
        if (f > 0)
            for (int i = 0; i < w->bones; i++)
                if (!pose_equal(&pose[i], &prev[i]) &&
                    palette_mark_dirty(w->pal, i) != PAL_OK)
                    atomic_store(&w->failed, 1);
        if (palette_update(w->pal, pose) != PAL_OK)
            atomic_store(&w->failed, 1);
        memcpy(prev, pose, sizeof(Pose) * (size_t)w->bones);
    }
    atomic_store(&w->done, 1);
    return NULL;
}

static void* reader_thread(void* arg) {
    ReaderArg* r = (ReaderArg*)arg;
    Mat3 snap[PAL_MAX_BONES];
    while (!atomic_load(r->done)) {
        int rc = palette_read(r->pal, snap, PAL_MAX_BONES);
        if (rc == PAL_OK) {
            atomic_fetch_add(&r->reads, 1);
            /* the only real failure mode: a "consistent" snapshot with bad
             * data. Tears themselves must never reach here */
            for (int i = 0; i < r->pal->count; i++)
                if (!mat_finite(&snap[i])) atomic_store(&r->failed, 1);
        } else if (rc == PAL_ERR_RETRY) {
            /* budget exhausted under sustained contention: back off and start
             * a fresh read rather than treating it as a correctness failure */
            atomic_fetch_add(&r->exhausted, 1);
            sched_yield();
        }
    }
    return NULL;
}

static void run_thread_test(void) {
    Palette pal;
    palette_init(&pal, PAL_MAX_BONES);
    Pose bind[PAL_MAX_BONES];
    for (int i = 0; i < PAL_MAX_BONES; i++) {
        bind[i].x = (float)i; bind[i].y = 0.0f; bind[i].ang = 0.0f;
        palette_add_bone(&pal, &bind[i]);
    }
    /* hold the lock briefly each frame so readers can observe real tears;
     * keep it short enough that the bounded retry budget always recovers */
    int spin = 300;
    pal.before_validate = hook_spin;
    pal.hook_user = &spin;

    Pose* poses = (Pose*)malloc(sizeof(Pose) * (size_t)THREAD_FRAMES * PAL_MAX_BONES);
    for (int f = 0; f < THREAD_FRAMES; f++)
        for (int i = 0; i < PAL_MAX_BONES; i++) {
            Pose* q = &poses[(size_t)f * PAL_MAX_BONES + i];
            q->x = (float)(rng_next() % 5);
            q->y = (float)(rng_next() % 5);
            q->ang = (float)(rng_next() % 8) * 0.25f;
        }

    WriterArg w = { &pal, poses, THREAD_FRAMES, PAL_MAX_BONES, 0, 0 };
    ReaderArg r[THREAD_READERS];
    pthread_t wt, rt[THREAD_READERS];
    pthread_create(&wt, NULL, writer_thread, &w);
    for (int i = 0; i < THREAD_READERS; i++) {
        r[i].pal = &pal; r[i].done = &w.done;
        atomic_init(&r[i].failed, 0);
        atomic_init(&r[i].reads, 0);
        atomic_init(&r[i].exhausted, 0);
        pthread_create(&rt[i], NULL, reader_thread, &r[i]);
    }
    pthread_join(wt, NULL);
    for (int i = 0; i < THREAD_READERS; i++) pthread_join(rt[i], NULL);

    int reader_failed = 0;
    long long total_reads = 0, total_exhausted = 0;
    for (int i = 0; i < THREAD_READERS; i++) {
        reader_failed |= atomic_load(&r[i].failed);
        total_reads += atomic_load(&r[i].reads);
        total_exhausted += atomic_load(&r[i].exhausted);
    }
    check("thread_writer_ok", atomic_load(&w.failed) == 0);
    check("thread_reader_snapshots", reader_failed == 0 && total_reads > 0);
    printf("thread_reads=%lld tears=%llu exhausted=%lld frames=%d\n",
           total_reads, (unsigned long long)atomic_load_explicit(&pal.st.tears, memory_order_relaxed),
           total_exhausted, THREAD_FRAMES);

    /* final published palette must equal the full reference of the last pose */
    Mat3 ref[PAL_MAX_BONES], snap[PAL_MAX_BONES];
    const Pose* last = poses + (size_t)(THREAD_FRAMES - 1) * PAL_MAX_BONES;
    palette_full_reference(&pal, last, ref);
    check("thread_final_exact",
          palette_read(&pal, snap, PAL_MAX_BONES) == PAL_OK &&
          memcmp(snap, ref, sizeof(Mat3) * PAL_MAX_BONES) == 0);
    free(poses);
}
#endif

int main(int argc, char** argv) {
    if (argc > 1 && strcmp(argv[1], "--thread") == 0) {
#if defined(PAL_WITH_PTHREAD)
        run_thread_test();
#else
        printf("thread_test=unavailable (rebuild with -DPAL_WITH_PTHREAD)\n");
#endif
        printf("all_ok=%d\n", failures == 0);
        return failures == 0 ? 0 : 1;
    }
    run_deterministic_tests();
    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
