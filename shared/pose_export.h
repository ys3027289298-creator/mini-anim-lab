/*
 * pose_export.h - unified pose export format shared by skeleton/ik/curve/path.
 *
 * New format (opt-in via "<prog> --pose-export <outdir>"):
 *   <outdir>/pose.json   pose values, image hash, two timing tiers,
 *                        reproduction command and format metadata
 *   <outdir>/pose.ppm    P6 PPM rendering of the pose (deterministic)
 *
 * With no arguments every program keeps its historical stdout byte-for-byte.
 *
 * This header is fully self-contained (static inline, no extra build flags) so
 * each module directory still builds independently.
 */
#ifndef POSE_EXPORT_H
#define POSE_EXPORT_H

/* clock_gettime/mkdir under strict -std=c11 */
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 199309L
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>

#define POSE_FORMAT_VERSION 1
#define POSE_CANVAS_W 320
#define POSE_CANVAS_H 240

/* ------------------------------------------------------------------ timing */

static inline long pe_now_us(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (long)ts.tv_sec * 1000000L + ts.tv_nsec / 1000L;
}

/* ------------------------------------------------------------- mkdir -p-ish */

static inline int pe_mkdir_p(const char* path) {
    char buf[512];
    size_t n = strlen(path);
    if (n == 0 || n >= sizeof(buf)) return -1;
    memcpy(buf, path, n + 1);
    for (size_t i = 1; i < n; i++) {
        if (buf[i] == '/') {
            buf[i] = '\0';
            if (mkdir(buf, 0775) != 0 && errno != EEXIST) return -1;
            buf[i] = '/';
        }
    }
    if (mkdir(buf, 0775) != 0 && errno != EEXIST) return -1;
    return 0;
}

/* --------------------------------------------------------- json string write */

static inline void pe_write_json_string(FILE* f, const char* s) {
    fputc('"', f);
    for (const unsigned char* p = (const unsigned char*)s; *p; p++) {
        switch (*p) {
            case '"':  fputs("\\\"", f); break;
            case '\\': fputs("\\\\", f); break;
            case '\b': fputs("\\b", f); break;
            case '\f': fputs("\\f", f); break;
            case '\n': fputs("\\n", f); break;
            case '\r': fputs("\\r", f); break;
            case '\t': fputs("\\t", f); break;
            default:
                if (*p < 0x20) fprintf(f, "\\u%04x", *p);
                else fputc(*p, f);
        }
    }
    fputc('"', f);
}

/* ------------------------------------------------------------------- canvas */

typedef struct {
    int w, h;
    unsigned char px[POSE_CANVAS_H * POSE_CANVAS_W * 3];
    /* world -> pixel mapping (y flipped) */
    float xmin, ymin, scale;
} PeCanvas;

static inline void pe_canvas_init(PeCanvas* c, float xmin, float ymin,
                                  float xmax, float ymax) {
    c->w = POSE_CANVAS_W;
    c->h = POSE_CANVAS_H;
    memset(c->px, 255, sizeof(c->px));
    c->xmin = xmin;
    c->ymin = ymin;
    const float margin = 12.0f;
    float sx = ((float)c->w - 2 * margin) / (xmax - xmin);
    float sy = ((float)c->h - 2 * margin) / (ymax - ymin);
    c->scale = sx < sy ? sx : sy;
    /* store the un-flipped y origin separately via ymin; mapping below */
    (void)0;
}

static inline void pe_map(const PeCanvas* c, float wx, float wy,
                          int* px, int* py) {
    const float margin = 12.0f;
    *px = (int)(margin + (wx - c->xmin) * c->scale + 0.5f);
    float h = (float)c->h;
    *py = (int)(h - margin - (wy - c->ymin) * c->scale + 0.5f);
}

static inline void pe_put(PeCanvas* c, int x, int y,
                          unsigned char r, unsigned char g, unsigned char b) {
    if (x < 0 || y < 0 || x >= c->w || y >= c->h) return;
    unsigned char* p = &c->px[(y * c->w + x) * 3];
    p[0] = r; p[1] = g; p[2] = b;
}

static inline void pe_line(PeCanvas* c, float x0, float y0, float x1, float y1,
                           unsigned char r, unsigned char g, unsigned char b) {
    int ax, ay, bx, by;
    pe_map(c, x0, y0, &ax, &ay);
    pe_map(c, x1, y1, &bx, &by);
    int dx = abs(bx - ax), dy = -abs(by - ay);
    int sx = ax < bx ? 1 : -1, sy = ay < by ? 1 : -1;
    int err = dx + dy;
    for (;;) {
        pe_put(c, ax, ay, r, g, b);
        if (ax == bx && ay == by) break;
        int e2 = 2 * err;
        if (e2 >= dy) { err += dy; ax += sx; }
        if (e2 <= dx) { err += dx; ay += sy; }
    }
}

static inline void pe_circle(PeCanvas* c, float cx, float cy, float rad,
                             unsigned char r, unsigned char g, unsigned char b) {
    int pcx, pcy, prad;
    pe_map(c, cx, cy, &pcx, &pcy);
    prad = (int)(rad * c->scale + 0.5f);
    for (int y = -prad; y <= prad; y++) {
        for (int x = -prad; x <= prad; x++) {
            if (x * x + y * y <= prad * prad) pe_put(c, pcx + x, pcy + y, r, g, b);
        }
    }
}

static inline int pe_canvas_save(const PeCanvas* c, const char* path) {
    FILE* f = fopen(path, "wb");
    if (!f) return -1;
    fprintf(f, "P6\n%d %d\n255\n", c->w, c->h);
    if (fwrite(c->px, 1, sizeof(c->px), f) != sizeof(c->px)) {
        fclose(f);
        return -1;
    }
    fclose(f);
    return 0;
}

/* ------------------------------------------------------------------ sha-256 */

typedef struct {
    uint32_t state[8];
    uint64_t bitlen;
    uint8_t buf[64];
    size_t buflen;
} PeSha256;

static const uint32_t PE_SHA256_K[64] = {
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1,
    0x923f82a4, 0xab1c5ed5, 0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
    0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174, 0xe49b69c1, 0xefbe4786,
    0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147,
    0x06ca6351, 0x14292967, 0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
    0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85, 0xa2bfe8a1, 0xa81a664b,
    0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a,
    0x5b9cca4f, 0x682e6ff3, 0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
    0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
};

#define PE_ROTR(x, n) (((x) >> (n)) | ((x) << (32 - (n))))

static inline void pe_sha256_transform(PeSha256* c, const uint8_t* data) {
    uint32_t m[64];
    for (int i = 0, j = 0; i < 16; i++, j += 4)
        m[i] = ((uint32_t)data[j] << 24) | ((uint32_t)data[j + 1] << 16) |
               ((uint32_t)data[j + 2] << 8) | (uint32_t)data[j + 3];
    for (int i = 16; i < 64; i++) {
        uint32_t s0 = PE_ROTR(m[i - 15], 7) ^ PE_ROTR(m[i - 15], 18) ^ (m[i - 15] >> 3);
        uint32_t s1 = PE_ROTR(m[i - 2], 17) ^ PE_ROTR(m[i - 2], 19) ^ (m[i - 2] >> 10);
        m[i] = m[i - 16] + s0 + m[i - 7] + s1;
    }
    uint32_t a = c->state[0], b = c->state[1], cc = c->state[2], d = c->state[3];
    uint32_t e = c->state[4], f = c->state[5], g = c->state[6], h = c->state[7];
    for (int i = 0; i < 64; i++) {
        uint32_t S1 = PE_ROTR(e, 6) ^ PE_ROTR(e, 11) ^ PE_ROTR(e, 25);
        uint32_t ch = (e & f) ^ ((~e) & g);
        uint32_t t1 = h + S1 + ch + PE_SHA256_K[i] + m[i];
        uint32_t S0 = PE_ROTR(a, 2) ^ PE_ROTR(a, 13) ^ PE_ROTR(a, 22);
        uint32_t maj = (a & b) ^ (a & cc) ^ (b & cc);
        uint32_t t2 = S0 + maj;
        h = g; g = f; f = e; e = d + t1;
        d = cc; cc = b; b = a; a = t1 + t2;
    }
    c->state[0] += a; c->state[1] += b; c->state[2] += cc; c->state[3] += d;
    c->state[4] += e; c->state[5] += f; c->state[6] += g; c->state[7] += h;
}

static inline void pe_sha256_init(PeSha256* c) {
    c->state[0] = 0x6a09e667; c->state[1] = 0xbb67ae85;
    c->state[2] = 0x3c6ef372; c->state[3] = 0xa54ff53a;
    c->state[4] = 0x510e527f; c->state[5] = 0x9b05688c;
    c->state[6] = 0x1f83d9ab; c->state[7] = 0x5be0cd19;
    c->bitlen = 0;
    c->buflen = 0;
}

static inline void pe_sha256_update(PeSha256* c, const uint8_t* data, size_t len) {
    for (size_t i = 0; i < len; i++) {
        c->buf[c->buflen++] = data[i];
        if (c->buflen == 64) {
            pe_sha256_transform(c, c->buf);
            c->bitlen += 512;
            c->buflen = 0;
        }
    }
}

static inline void pe_sha256_hex(const PeSha256* src, char out[65]) {
    PeSha256 c = *src;
    uint8_t pad[64] = { 0x80 };
    uint64_t bitlen = c.bitlen + (uint64_t)c.buflen * 8;
    size_t padlen = c.buflen < 56 ? 56 - c.buflen : 120 - c.buflen;
    pe_sha256_update(&c, pad, padlen);
    uint8_t lenbe[8];
    for (int i = 7; i >= 0; i--) lenbe[i] = (uint8_t)(bitlen >> (7 - i) * 8);
    pe_sha256_update(&c, lenbe, 8);
    for (int i = 0; i < 8; i++)
        sprintf(out + i * 8, "%08x", c.state[i]);
    out[64] = '\0';
}

/* hash a whole file; returns 0 on success */
static inline int pe_sha256_file(const char* path, char hex[65]) {
    FILE* f = fopen(path, "rb");
    if (!f) return -1;
    PeSha256 ctx;
    pe_sha256_init(&ctx);
    uint8_t chunk[4096];
    size_t n;
    while ((n = fread(chunk, 1, sizeof(chunk), f)) > 0)
        pe_sha256_update(&ctx, chunk, n);
    fclose(f);
    pe_sha256_hex(&ctx, hex);
    return 0;
}

#endif /* POSE_EXPORT_H */
