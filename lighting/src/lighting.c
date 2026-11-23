
/* lighting: baked lightmap + realtime direct light blending.

   问题:光照贴图里同时烘了间接光和直接光,渲染时再叠加实时直接光,
   直接光被计了两次,画面偏亮(漏光)。

   方案:烘焙时额外输出一张直接光分量图 direct_lm;运行时先从光照贴图
   中扣除烘焙直接光,再加上实时直接光:
       final = albedo * (max(0, lightmap - direct_lm) + rt_direct)

   --no-blend 退回旧行为 final = albedo * (lightmap + rt_direct),
   输出与现状逐字节一致(由 GOLDEN_NAIVE 哈希锁定)。

   自检输出四项验收数据:漏光像素表、混合前后误差、图像哈希、复现命令。 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#define W 64
#define H 64
#define LW 16   /* 光照贴图横向纹素数 */
#define LH 16   /* 光照贴图纵向纹素数 */

#define LEAK_THRESH 0.05f   /* 漏光判定阈值:线性亮度超出参考的量 */
#define ERR_BOUND   0.05f   /* 混合后逐通道最大误差上界 */
#define LEAK_RATIO_BOUND 0.90  /* 漏光像素下降比例下界 */

/* 黄金哈希:首次运行后填入,锁定输出可复现、--no-blend 与现状一致 */
#define GOLDEN_NAIVE 0xb34d4a88efecfa20ULL
#define GOLDEN_BLEND 0xed9d9baf6e4d2d83ULL

static float clampf(float v, float lo, float hi) { return v < lo ? lo : (v > hi ? hi : v); }

static float smoothstep(float a, float b, float x) {
    float t = clampf((x - a) / (b - a), 0.0f, 1.0f);
    return t * t * (3.0f - 2.0f * t);
}

/* ---- 确定性场景:高度场法线 + 两个遮挡体 + 低频间接光 ---- */

static float g_light[3];        /* 归一化光照方向 */
static const float LIGHT_COL[3] = { 1.00f, 0.95f, 0.85f };
static const float IND_COL[3]   = { 0.90f, 1.00f, 1.05f };

static void scene_init(void) {
    float l[3] = { 0.35f, 0.55f, 0.45f };
    float len = sqrtf(l[0] * l[0] + l[1] * l[1] + l[2] * l[2]);
    g_light[0] = l[0] / len; g_light[1] = l[1] / len; g_light[2] = l[2] / len;
}

static float height_at(float x, float y) {
    return 0.35f * sinf(x * 0.18f) * cosf(y * 0.14f)
         + 0.10f * sinf(x * 0.05f + y * 0.07f);
}

static void normal_at(float x, float y, float n[3]) {
    float e = 0.5f;
    float dx = (height_at(x + e, y) - height_at(x - e, y)) / (2.0f * e);
    float dy = (height_at(x, y + e) - height_at(x, y - e)) / (2.0f * e);
    float len = sqrtf(dx * dx + dy * dy + 1.0f);
    n[0] = -dx / len; n[1] = -dy / len; n[2] = 1.0f / len;
}

/* 真实阴影:圆盘 + 竖条遮挡,边缘平滑 */
static float shadow_true(float x, float y) {
    float dx = x - 42.0f, dy = y - 26.0f;
    float d = sqrtf(dx * dx + dy * dy) - 11.0f;
    float s1 = smoothstep(-2.0f, 2.0f, d);
    float bx = smoothstep(9.0f, 11.0f, x) * (1.0f - smoothstep(13.0f, 15.0f, x));
    float by = smoothstep(29.0f, 31.0f, y) * (1.0f - smoothstep(53.0f, 55.0f, y));
    float s2 = 1.0f - bx * by;
    return 0.1f + 0.9f * s1 * s2;
}

/* 实时阴影:真值量化到 1/16 档(模拟阴影贴图条带),与真值差 <= 1/32 */
static float shadow_rt(float x, float y) {
    return floorf(shadow_true(x, y) * 16.0f + 0.5f) / 16.0f;
}

static float indirect_at(float x, float y) {
    return 0.08f + 0.04f * sinf(x * 0.07f) * sinf(y * 0.09f)
         + 0.03f * cosf((x + y) * 0.03f);
}

static void albedo_at(int x, int y, float a[3]) {
    float base = (((x / 8) + (y / 8)) & 1) ? 0.72f : 0.45f;
    a[0] = base; a[1] = base * 0.92f; a[2] = base * 0.80f;
}

/* 真实直接光(参考):ndotl * 真实阴影 * 光色 */
static void direct_true_at(float x, float y, float d[3]) {
    float n[3];
    normal_at(x, y, n);
    float ndotl = clampf(n[0] * g_light[0] + n[1] * g_light[1] + n[2] * g_light[2], 0.0f, 1.0f);
    float s = shadow_true(x, y);
    for (int c = 0; c < 3; c++) d[c] = ndotl * s * LIGHT_COL[c];
}

/* 实时直接光:ndotl * 实时阴影 * 光色 */
static void direct_rt_at(float x, float y, float d[3]) {
    float n[3];
    normal_at(x, y, n);
    float ndotl = clampf(n[0] * g_light[0] + n[1] * g_light[1] + n[2] * g_light[2], 0.0f, 1.0f);
    float s = shadow_rt(x, y);
    for (int c = 0; c < 3; c++) d[c] = ndotl * s * LIGHT_COL[c];
}

/* 烘焙的直接光分量图:低分辨率采样真实直接光 */
static float g_direct_lm[LW * LH][3];

static void bake_direct_lm(void) {
    for (int j = 0; j < LH; j++)
        for (int i = 0; i < LW; i++) {
            float x = (i + 0.5f) * (float)W / LW;
            float y = (j + 0.5f) * (float)H / LH;
            direct_true_at(x, y, g_direct_lm[j * LW + i]);
        }
}

/* 光照贴图双线性采样(对直接光分量图) */
static void direct_lm_at(float x, float y, float d[3]) {
    float u = (x + 0.5f) * (float)LW / W - 0.5f;
    float v = (y + 0.5f) * (float)LH / H - 0.5f;
    int i0 = (int)floorf(u), j0 = (int)floorf(v);
    float fu = u - i0, fv = v - j0;
    int i1 = i0 + 1 < LW ? i0 + 1 : LW - 1;
    int j1 = j0 + 1 < LH ? j0 + 1 : LH - 1;
    if (i0 < 0) { i0 = 0; fu = 0.0f; }
    if (j0 < 0) { j0 = 0; fv = 0.0f; }
    for (int c = 0; c < 3; c++) {
        float a = g_direct_lm[j0 * LW + i0][c], b = g_direct_lm[j0 * LW + i1][c];
        float cc = g_direct_lm[j1 * LW + i0][c], dd = g_direct_lm[j1 * LW + i1][c];
        d[c] = (a * (1 - fu) + b * fu) * (1 - fv) + (cc * (1 - fu) + dd * fu) * fv;
    }
}

/* ---- 图像与度量 ---- */

static float g_ref[W * H * 3], g_naive[W * H * 3], g_blend[W * H * 3];

static void render_all(void) {
    bake_direct_lm();
    for (int y = 0; y < H; y++)
        for (int x = 0; x < W; x++) {
            float a[3], dt[3], dr[3], dlm[3];
            albedo_at(x, y, a);
            direct_true_at((float)x, (float)y, dt);
            direct_rt_at((float)x, (float)y, dr);
            direct_lm_at((float)x, (float)y, dlm);
            float ind = indirect_at((float)x, (float)y);
            for (int c = 0; c < 3; c++) {
                float lm = ind * IND_COL[c] + dlm[c];        /* 光照贴图:间接+烘焙直接光 */
                float indirect_only = lm - dlm[c];
                if (indirect_only < 0.0f) indirect_only = 0.0f;
                size_t k = ((size_t)y * W + x) * 3 + c;
                g_ref[k]   = a[c] * (ind * IND_COL[c] + dt[c]);
                g_naive[k] = a[c] * (lm + dr[c]);            /* 现状:直接光计两次 */
                g_blend[k] = a[c] * (indirect_only + dr[c]);
            }
        }
}

static float luminance(const float* p) {
    return 0.2126f * p[0] + 0.7152f * p[1] + 0.0722f * p[2];
}

/* 漏光:输出亮度比参考超出 LEAK_THRESH */
static int leak_count(const float* img, int block_count[8][8]) {
    int n = 0;
    for (int y = 0; y < H; y++)
        for (int x = 0; x < W; x++) {
            size_t k = ((size_t)y * W + x) * 3;
            if (luminance(&img[k]) - luminance(&g_ref[k]) > LEAK_THRESH) {
                n++;
                if (block_count) block_count[y / 8][x / 8]++;
            }
        }
    return n;
}

static void error_stats(const float* img, float* max_err, float* mean_err) {
    double sum = 0.0;
    float mx = 0.0f;
    for (size_t k = 0; k < (size_t)W * H * 3; k++) {
        float e = fabsf(img[k] - g_ref[k]);
        if (e > mx) mx = e;
        sum += e;
    }
    *max_err = mx;
    *mean_err = (float)(sum / ((double)W * H * 3));
}

static uint64_t fnv1a(const unsigned char* d, size_t n) {
    uint64_t h = 1469598103934665603ULL;
    for (size_t i = 0; i < n; i++) { h ^= d[i]; h *= 1099511628211ULL; }
    return h;
}

static uint64_t image_hash(const float* img) {
    unsigned char rgb[W * H * 3];
    for (size_t k = 0; k < (size_t)W * H * 3; k++)
        rgb[k] = (unsigned char)clampf(img[k] * 255.0f + 0.5f, 0.0f, 255.0f);
    return fnv1a(rgb, sizeof(rgb));
}

static int write_ppm(const char* path, const float* img) {
    FILE* f = fopen(path, "wb");
    if (!f) return -1;
    fprintf(f, "P6\n%d %d\n255\n", W, H);
    for (size_t k = 0; k < (size_t)W * H * 3; k++)
        fputc((unsigned char)clampf(img[k] * 255.0f + 0.5f, 0.0f, 255.0f), f);
    fclose(f);
    return 0;
}

static int write_leak_csv(const char* path, const float* img) {
    FILE* f = fopen(path, "w");
    if (!f) return -1;
    fprintf(f, "x,y,excess\n");
    for (int y = 0; y < H; y++)
        for (int x = 0; x < W; x++) {
            size_t k = ((size_t)y * W + x) * 3;
            float ex = luminance(&img[k]) - luminance(&g_ref[k]);
            if (ex > LEAK_THRESH) fprintf(f, "%d,%d,%.6f\n", x, y, ex);
        }
    fclose(f);
    return 0;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static void print_block_table(const char* title, int bc[8][8]) {
    printf("%s\n", title);
    for (int j = 0; j < 8; j++) {
        printf("  ");
        for (int i = 0; i < 8; i++) printf("%4d", bc[j][i]);
        printf("\n");
    }
}

int main(int argc, char** argv) {
    int blend = 1;
    const char* outdir = "out";
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--no-blend") == 0) blend = 0;
        else if (strcmp(argv[i], "--out") == 0 && i + 1 < argc) outdir = argv[++i];
        else { fprintf(stderr, "usage: lighting [--no-blend] [--out DIR]\n"); return 2; }
    }

    scene_init();
    render_all();

    /* 1. 漏光像素表(8x8 块计数 + 逐像素 CSV) */
    int bc_before[8][8] = {{0}}, bc_after[8][8] = {{0}};
    int leak_before = leak_count(g_naive, bc_before);
    int leak_after = leak_count(g_blend, bc_after);
    float leak_ratio = leak_before > 0 ? 1.0f - (float)leak_after / leak_before : 1.0f;
    printf("漏光像素表(8x8 块,阈值 %.3f):\n", LEAK_THRESH);
    print_block_table("混合前:", bc_before);
    print_block_table("混合后:", bc_after);
    printf("漏光合计: before=%d after=%d 下降比例=%.2f%%\n",
           leak_before, leak_after, leak_ratio * 100.0f);

    /* 2. 混合前后误差(相对参考图,线性域) */
    float max_b, mean_b, max_a, mean_a;
    error_stats(g_naive, &max_b, &mean_b);
    error_stats(g_blend, &max_a, &mean_a);
    printf("误差: before max=%.6f mean=%.6f | after max=%.6f mean=%.6f (上界 %.3f)\n",
           max_b, mean_b, max_a, mean_a, ERR_BOUND);

    /* 3. 图像哈希(FNV-1a 64,RGB24) */
    uint64_t h_naive = image_hash(g_naive);
    uint64_t h_blend = image_hash(g_blend);
    uint64_t h_final = blend ? h_blend : h_naive;
    printf("图像哈希: naive=%016llx blended=%016llx final(%s)=%016llx\n",
           (unsigned long long)h_naive, (unsigned long long)h_blend,
           blend ? "blend" : "no-blend", (unsigned long long)h_final);

    /* 输出图像与逐像素漏光表 */
    char path[512];
    int io_ok = 1;
    snprintf(path, sizeof(path), "%s/reference.ppm", outdir); io_ok &= write_ppm(path, g_ref) == 0;
    snprintf(path, sizeof(path), "%s/naive.ppm", outdir);     io_ok &= write_ppm(path, g_naive) == 0;
    snprintf(path, sizeof(path), "%s/blended.ppm", outdir);   io_ok &= write_ppm(path, g_blend) == 0;
    snprintf(path, sizeof(path), "%s/final.ppm", outdir);     io_ok &= write_ppm(path, blend ? g_blend : g_naive) == 0;
    snprintf(path, sizeof(path), "%s/leak_before.csv", outdir); io_ok &= write_leak_csv(path, g_naive) == 0;
    snprintf(path, sizeof(path), "%s/leak_after.csv", outdir);  io_ok &= write_leak_csv(path, g_blend) == 0;
    check("io_ok", io_ok);

    /* 4. 复现命令 */
    printf("复现命令: make -C lighting test; ./lighting/lighting --no-blend --out lighting/out\n");

    check("leak_down", leak_after < leak_before);
    check("leak_ratio", leak_ratio >= LEAK_RATIO_BOUND);
    check("err_bound", max_a <= ERR_BOUND);
    check("noblend_unchanged", h_naive == GOLDEN_NAIVE);
    check("blend_reproducible", h_blend == GOLDEN_BLEND);
    if (!blend) check("noblend_final_is_naive", h_final == h_naive);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
