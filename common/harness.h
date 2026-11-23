/* harness.h: shared self-check scaffolding for the mini-anim-lab modules.
 *
 * Each module keeps its own freestanding main(). The unified runner compiles
 * the same source with extra flags that make main() print evidence *after*
 * the normal checks:
 *
 *   MAL_EVIDENCE_POSE    print mal_pose lines for a named pose
 *   MAL_EVIDENCE_TIMING  run a fixed workload, print mal_time fast/slow
 *   MAL_EVIDENCE_IMAGE   rasterize a pose, print mal_png with a CRC
 *   MAL_EVIDENCE_REPLAY  print the repro command for this build
 *
 * Evidence never calls mal_check(), so the check list is identical in every
 * build mode and the validator can compare it against the golden baseline.
 */
#ifndef MAL_HARNESS_H
#define MAL_HARNESS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#ifndef MAL_MODULE
#error "MAL_MODULE must name the module being built"
#endif

#if defined(__GNUC__)
#define MAL_UNUSED __attribute__((unused))
#else
#define MAL_UNUSED
#endif

#ifndef MAL_TIME_ITERS
#define MAL_TIME_ITERS 200000L
#endif

#ifndef MAL_TIME_TIER
#define MAL_TIME_TIER "fast"
#endif

static int mal_failures = 0;
static int mal_checks = 0;

static inline void mal_check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    mal_checks++;
    if (!ok) mal_failures++;
}

/* ---- pose evidence: named scalar values, order preserved --------------- */

#define MAL_POSE_MAX 64
static char mal_pose_names[MAL_POSE_MAX][24];
static double mal_pose_values[MAL_POSE_MAX];
static int mal_pose_count = 0;

static inline void mal_pose(const char* name, double value) {
    if (mal_pose_count >= MAL_POSE_MAX) {
        fprintf(stderr, "mal: pose evidence overflow in %s\n", MAL_MODULE);
        exit(2);
    }
    snprintf(mal_pose_names[mal_pose_count], sizeof mal_pose_names[mal_pose_count],
             "%s", name);
    mal_pose_values[mal_pose_count] = value;
    mal_pose_count++;
}

static inline void mal_emit_pose(void) {
    for (int i = 0; i < mal_pose_count; i++)
        printf("mal_pose k=%s v=%.9g\n", mal_pose_names[i], mal_pose_values[i]);
    unsigned long long h = 1469598103934665603ULL;
    for (int i = 0; i < mal_pose_count; i++) {
        unsigned long long bits;
        memcpy(&bits, &mal_pose_values[i], sizeof(bits));
        h ^= bits;
        h *= 1099511628211ULL;
    }
    printf("mal_pose_hash v=%016llx n=%d\n", h, mal_pose_count);
}

/* ---- timing evidence ---------------------------------------------------
 * Two builds of the same source, one with -O0 and one with -O2, run the same
 * iteration count so the fast/slow pair is a real measurement, not a guess.
 */

static inline double mal_now(void) {
    struct timespec ts;
    timespec_get(&ts, TIME_UTC);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

static inline void mal_emit_time(const char* tier, double seconds, long iterations) {
    printf("mal_time tier=%s ms=%.3f iterations=%ld\n", tier, seconds * 1000.0, iterations);
}

/* ---- image evidence: 128x128 line raster, 4x4 supersampled ------------ */

#define MAL_IMG 128
#define MAL_SS 4

typedef struct {
    unsigned char r[MAL_IMG * MAL_IMG];
    unsigned char g[MAL_IMG * MAL_IMG];
    unsigned char b[MAL_IMG * MAL_IMG];
} MalImage;

#define MAL_BG_R 10
#define MAL_BG_G 10
#define MAL_BG_B 16

static inline void mal_img_clear(MalImage* im) {
    for (int i = 0; i < MAL_IMG * MAL_IMG; i++) {
        im->r[i] = MAL_BG_R; im->g[i] = MAL_BG_G; im->b[i] = MAL_BG_B;
    }
}

static inline void mal_img_line(MalImage* im, double x0, double y0, double x1, double y1,
                                unsigned char cr, unsigned char cg, unsigned char cb) {
    const double half = MAL_SS / 2.0;
    const double inv = 1.0 / MAL_SS;
    const double vx = x1 - x0, vy = y1 - y0;
    const double vv = vx * vx + vy * vy;
    for (int j = 0; j < MAL_IMG; j++) {
        for (int i = 0; i < MAL_IMG; i++) {
            int hits = 0;
            for (int sy = 0; sy < MAL_SS; sy++) {
                for (int sx = 0; sx < MAL_SS; sx++) {
                    double px = i + (sx + 0.5) * inv;
                    double py = j + (sy + 0.5) * inv;
                    double wx = px - x0, wy = py - y0;
                    double t = vv > 0.0 ? (wx * vx + wy * vy) / vv : 0.0;
                    if (t < 0.0) t = 0.0;
                    if (t > 1.0) t = 1.0;
                    double dx = wx - vx * t, dy = wy - vy * t;
                    if (dx * dx + dy * dy <= half * half) hits++;
                }
            }
            if (hits == 0) continue;
            int k = j * MAL_IMG + i;
            int a = hits * 255 / (MAL_SS * MAL_SS);
            im->r[k] = (unsigned char)(((int)im->r[k] * (255 - a) + cr * a) / 255);
            im->g[k] = (unsigned char)(((int)im->g[k] * (255 - a) + cg * a) / 255);
            im->b[k] = (unsigned char)(((int)im->b[k] * (255 - a) + cb * a) / 255);
        }
    }
}

static inline unsigned long mal_crc32(const unsigned char* data, size_t n) {
    unsigned long crc = 0xFFFFFFFFUL;
    for (size_t i = 0; i < n; i++) {
        crc ^= data[i];
        for (int b = 0; b < 8; b++)
            crc = (crc >> 1) ^ (0xEDB88320UL & (unsigned long)(-(long)(crc & 1)));
    }
    return crc ^ 0xFFFFFFFFUL;
}

static inline void mal_emit_png(MalImage* im, const char* label) {
    unsigned long crc = mal_crc32(im->r, MAL_IMG * MAL_IMG);
    crc = mal_crc32(im->g, MAL_IMG * MAL_IMG) ^ crc;
    crc = mal_crc32(im->b, MAL_IMG * MAL_IMG) ^ crc;
    printf("mal_png label=%s w=%d h=%d crc=%08lx\n", label, MAL_IMG, MAL_IMG, crc & 0xFFFFFFFFUL);
}

#endif /* MAL_HARNESS_H */
