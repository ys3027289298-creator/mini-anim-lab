
/* quat: smallest-three quaternion compression (3 fixed-point values + 2-bit slot index).
 *
 * Encoding: find the largest-magnitude component, flip the whole quaternion so
 * that component is positive (q and -q are the same rotation), and store the
 * other three. Decode rebuilds the dropped component as
 * +sqrt(1 - x^2 - y^2 - z^2) and renormalizes, so the kept three need not be
 * exact. The dropped-component index (2 bits) has to travel with the data; the
 * size accounting below includes it.
 *
 * Error bound (see quat_angle_error_bound): each kept component lies in
 * [-1/sqrt2, 1/sqrt2] and is stored on a b-bit signed fixed-point grid, so it
 * can be off by h = 2^-b / sqrt2. The rotation-angle error is bounded by
 *   4 * asin( |p - q| / 2 ),   |p - q| <= sqrt(3) h + delta,
 *   delta = (3h + 3h^2) / ((1/4 - 3h - 3h^2) (1 + 1/sqrt(1/4 - 3h - 3h^2)))
 * with no dependence on the data; see the derivation at the function.
 */
#include <stdio.h>
#include <math.h>
#include <stdint.h>

#define QUAT_MIN_NORM   1e-6   /* below this norm the dropped component is unrecoverable */
#define QUAT_MIN_BITS   4      /* below this the analytic bound is vacuous */

typedef struct { double x, y, z, w; } Quat;

typedef struct {
    float q[3];   /* the three kept components, in encoded slot order */
    int   drop;   /* 0..3: index of the component that was left out */
} QuatCompressed;

static void q_normalize(Quat* q) {
    double l = sqrt(q->x * q->x + q->y * q->y + q->z * q->z + q->w * q->w);
    if (!(l > 0.0)) { q->x = q->y = q->z = 0.0; q->w = 1.0; return; }
    q->x /= l; q->y /= l; q->z /= l; q->w /= l;
}

static double q_norm(const Quat* q) {
    return sqrt(q->x * q->x + q->y * q->y + q->z * q->z + q->w * q->w);
}

/* angle of the relative rotation between two unit quaternions, in radians.
 * Uses 2*atan2(|v|,|w|) rather than 2*acos(w): acos has infinite slope at
 * w = 1, so a 1e-8 component error was read as ~5e-4 rad. */
static double q_angle_between(const Quat* a, const Quat* b) {
    double rw = a->w * b->w + a->x * b->x + a->y * b->y + a->z * b->z;
    double rx = a->w * b->x - a->x * b->w - a->y * b->z + a->z * b->y;
    double ry = a->w * b->y + a->x * b->z - a->y * b->w - a->z * b->x;
    double rz = a->w * b->z - a->x * b->y + a->y * b->x - a->z * b->w;
    if (rw < 0.0) { rw = -rw; rx = -rx; ry = -ry; rz = -rz; }
    return 2.0 * atan2(sqrt(rx * rx + ry * ry + rz * rz), rw);
}

static int q_compress(const Quat* src, QuatCompressed* out) {
    double n = q_norm(src);
    if (!(n >= QUAT_MIN_NORM)) return -1;   /* also rejects NaN */
    Quat q = *src;
    q_normalize(&q);
    double c[4] = { q.x, q.y, q.z, q.w };
    int drop = 0;
    for (int i = 1; i < 4; i++)
        if (fabs(c[i]) > fabs(c[drop])) drop = i;
    double sign = c[drop] < 0.0 ? -1.0 : 1.0;
    /* q and -q encode the same rotation: keep the copy whose largest
       component is positive so decode can take the +sqrt root */
    int k = 0;
    for (int i = 0; i < 4; i++) {
        if (i == drop) continue;
        out->q[k++] = (float)(sign * c[i]);
    }
    out->drop = drop;
    return 0;
}

static int q_decompress(const QuatCompressed* in, Quat* out) {
    if (in->drop < 0 || in->drop > 3) return -1;
    double c[4];
    int k = 0;
    double sum = 0.0;
    for (int i = 0; i < 4; i++) {
        if (i == in->drop) continue;
        c[i] = (double)in->q[k++];
        sum += c[i] * c[i];
    }
    if (sum > 1.0) return -1;   /* inconsistent payload: sqrt would be imaginary */
    c[in->drop] = sqrt(1.0 - sum);
    out->x = c[0]; out->y = c[1]; out->z = c[2]; out->w = c[3];
    q_normalize(out);
    return 0;
}

/* ---- rate / distortion model ------------------------------------------- */

#define QUAT_DROP_BITS 2

/* bits on the wire: 3 components * b bits + 2-bit slot index */
static double quat_compressed_bits(int bits_per_component) {
    return 3.0 * bits_per_component + QUAT_DROP_BITS;
}

/* max per-component error of a b-bit signed fixed grid over [-1/sqrt2, 1/sqrt2] */
static double quat_quant_step(int b) {
    return ldexp(1.0, -b) / sqrt(2.0);
}

/* Worst-case rotation angle (radians) for b-bit kept components, valid for
 * b >= 4. Set-up: q unit, d the dropped component, k the kept triple. Since d
 * is the largest, d >= 1/2 and |k|^2 = 1 - d^2 <= 3/4. Quantize k' = k + e,
 * |e| <= sqrt3 h. Decode forms p = (k', s) / sqrt(|k'|^2 + s^2) with
 * s = sqrt(1 - |k'|^2) (this is exactly renormalize of the rebuilt quaternion).
 *
 *   |k'|^2 = |k|^2 + 2 k.e + |e|^2 <= 3/4 + 3h + 3h^2 < 1, so s is real;
 *   |k'|^2 >= 1/4 - 3h - 3h^2, hence s >= 1/2.
 *   r = (k', s) has |r|^2 = |k'|^2 + s^2 = 1 + (2 k.e + |e|^2) = 1 + delta,
 *   |delta| <= 3h + 3h^2.
 *   p - q = (k', s)/|r| - (k, d), and (k, d) = q = (k, 1/|r|) + (0, d - 1/|r|):
 *      |p - q| <= |e|/|r| + |1 - 1/|r||   (the k-part)
 *               <= |e| + |1 - 1/|r||.
 *   |1 - 1/|r|| = |1 - |r|| / |r| = |delta| / (|r| (1 + |r|))
 *               <= (3h + 3h^2) / (r_min (1 + r_min))   with r_min = sqrt(1/4-3h-3h^2)
 *               (|r| >= s >= r_min, and t -> t(1+t) increasing for t > 0).
 *   |p - q| <= sqrt(3) h + delta with that denominator.
 * Two unit quaternions whose distance is |p - q| are rotations apart by
 * 4 * asin(|p - q| / 2). The dropped-side sign is forced positive on both
 * ends, so the quaternion dot product stays positive: no 2*pi wrap.
 * The bound is safe but not wildly loose: random + adversarial sweeps reach
 * ~0.6 of it. For b < 4 the radicand can go non-positive; only pi holds. */
static double quat_angle_error_bound(int b) {
    if (b < 4) return 3.14159265358979323846;
    double h = quat_quant_step(b);
    double r2_min = 0.25 - 3.0 * h - 3.0 * h * h;    /* |r|^2 lower bound */
    if (r2_min <= 0.0) return 3.14159265358979323846;
    double r_min = sqrt(r2_min);
    double num = 3.0 * h + 3.0 * h * h;               /* |delta| upper bound */
    double denom = r2_min * (1.0 + 1.0 / r_min);      /* |r| (1 + |r|) >= this */
    double mag = sqrt(3.0) * h + num / denom;
    if (mag >= 2.0) return 3.14159265358979323846;
    return 4.0 * asin(mag / 2.0);
}

/* b-bit signed fixed point over [-1/sqrt2, 1/sqrt2]; the offline encoder used
 * in the experiment. Kept values already lie in this range; clamping only
 * guards float rounding at the ends. */
static float quantize_fixed(float v, int b) {
    double lim = 1.0 / sqrt(2.0);
    double step = 2.0 * lim / ldexp(1.0, b);
    double qv = round((double)v / step) * step;
    if (qv > lim) qv = lim;
    if (qv < -lim) qv = -lim;
    return (float)qv;
}

/* ---- tests -------------------------------------------------------------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static void xorshift(uint32_t* s) {
    uint32_t x = *s;
    x ^= x << 13; x ^= x >> 17; x ^= x << 5;
    *s = x;
}

static double rand01(uint32_t* s) {
    xorshift(s);
    return (double)*s / 4294967295.0;
}

static Quat rand_quat(uint32_t* s) {
    Quat q;
    do {
        q.x = rand01(s) * 2 - 1; q.y = rand01(s) * 2 - 1;
        q.z = rand01(s) * 2 - 1; q.w = rand01(s) * 2 - 1;
    } while (q_norm(&q) < 0.1);
    q_normalize(&q);
    return q;
}

int main(void) {
    QuatCompressed c;
    Quat back;

    /* identity round trip */
    Quat id = { 0, 0, 0, 1 };
    check("identity_compress", q_compress(&id, &c) == 0);
    check("identity_drop_w", c.drop == 3);
    check("identity_kept_zero", c.q[0] == 0.0f && c.q[1] == 0.0f && c.q[2] == 0.0f);
    check("identity_decode", q_decompress(&c, &back) == 0);
    check("identity_error", q_angle_between(&id, &back) < 1e-6);

    /* 90 deg about z */
    Quat rz = { 0, 0, sin(0.7853981633974483), cos(0.7853981633974483) };
    check("rz90_compress", q_compress(&rz, &c) == 0);
    check("rz90_roundtrip", q_decompress(&c, &back) == 0 &&
          q_angle_between(&rz, &back) < 1e-6);

    /* q and -q encode the same rotation */
    Quat negid = { -id.x, -id.y, -id.z, -id.w };
    check("neg_identity_drop_w", q_compress(&negid, &c) == 0 && c.drop == 3);
    Quat negrz = { -rz.x, -rz.y, -rz.z, -rz.w };
    check("neg_rz_decode", q_compress(&negrz, &c) == 0 && q_decompress(&c, &back) == 0);
    check("neg_rz_same", q_angle_between(&rz, &back) < 1e-6);

    /* every component can be the dropped one */
    int all_slots = 1;
    for (int i = 0; i < 4; i++) {
        Quat q = { 0.1, 0.2, 0.3, 0.4 };
        if (i == 0) q.x = 0.9; else if (i == 1) q.y = 0.9;
        else if (i == 2) q.z = 0.9; else q.w = 0.9;
        q_normalize(&q);
        if (q_compress(&q, &c) != 0 || c.drop != i || q_decompress(&c, &back) != 0 ||
            q_angle_between(&q, &back) > 1e-6) all_slots = 0;
    }
    check("all_drop_slots", all_slots);

    /* degeneracies: too-small norm must be reported, unusable payloads rejected */
    Quat zero = { 0, 0, 0, 0 };
    check("zero_norm_compress_fails", q_compress(&zero, &c) == -1);
    Quat tiny = { 1e-9, 0, 0, 0 };
    check("tiny_norm_compress_fails", q_compress(&tiny, &c) == -1);
    Quat just_over = { 2e-6, 0, 0, 0 };
    check("just_over_threshold_ok", q_compress(&just_over, &c) == 0);
    QuatCompressed inconsistent = { { 0.8f, 0.8f, 0.8f }, 3 };
    check("decode_inconsistent", q_decompress(&inconsistent, &back) == -1);
    QuatCompressed badidx = { { 0.0f, 0.0f, 0.0f }, 7 };
    check("decode_bad_index", q_decompress(&badidx, &back) == -1);

    /* lossless (float) round trip over random unit quaternions */
    uint32_t seed = 12345u;
    double worst = 0.0;
    int rt_ok = 1;
    for (int i = 0; i < 100000; i++) {
        Quat q = rand_quat(&seed);
        if (q_compress(&q, &c) != 0 || q_decompress(&c, &back) != 0) { rt_ok = 0; break; }
        double e = q_angle_between(&q, &back);
        if (e > worst) worst = e;
    }
    check("random_roundtrip", rt_ok);
    check("random_roundtrip_exact", worst < 1e-6);

    /* distortion vs rate. For each b, measure the worst angle error under
     * random data and under adversarial +-h perturbations at every kept
     * component, and check it stays under the analytic bound. */
    printf("rate_distortion:\n");
    int bound_ok = 1;
    for (int b = QUAT_MIN_BITS; b <= 16; b += 2) {
        double h = quat_quant_step(b);
        double bound = quat_angle_error_bound(b);
        double measured = 0.0;
        uint32_t s = 987654321u;
        int decoded = 0, rejected = 0;
        for (int i = 0; i < 40000; i++) {
            Quat q = rand_quat(&s);
            if (q_compress(&q, &c) != 0) continue;
            /* real quantizer: snap each kept value to the grid, then decode */
            QuatCompressed snapped = c;
            for (int k = 0; k < 3; k++) snapped.q[k] = quantize_fixed(snapped.q[k], b);
            Quat dqs;
            if (q_decompress(&snapped, &dqs) == 0) {
                decoded++;
                double e = q_angle_between(&q, &dqs);
                if (e > measured) measured = e;
            } else rejected++;
            /* adversarial: each kept value may sit up to h on either side */
            for (int corner = 0; corner < 8; corner++) {
                QuatCompressed t = c;
                for (int k = 0; k < 3; k++)
                    t.q[k] += (corner & (1 << k)) ? (float)h : -(float)h;
                Quat dq;
                if (q_decompress(&t, &dq) != 0) { rejected++; continue; }
                decoded++;
                double e = q_angle_between(&q, &dq);
                if (e > measured) measured = e;
            }
        }
        double bits = quat_compressed_bits(b);
        printf("bits_per_component=%d total_bits=%.1f cr_vs_quat32=%.4f cr_vs_mat3=%.4f "
               "bound_rad=%.6g measured_rad=%.6g rejected=%d ok=%d\n",
               b, bits, bits / 128.0, bits / 288.0, bound, measured, rejected,
               (int)(measured <= bound * (1.0 + 1e-9)));
        if (measured > bound * (1.0 + 1e-9)) bound_ok = 0;
    }
    check("measured_within_bound", bound_ok);

    /* headline rate/error table for plausible animation budgets */
    printf("budget_table:\n");
    int table_ok = 1;
    for (int b = 4; b <= 16; b += 2) {
        double bits = quat_compressed_bits(b);
        double bound_deg = quat_angle_error_bound(b) * 180.0 / 3.14159265358979323846;
        printf("bits_per_component=%d total_bits=%.0f ratio_vs_mat3=%.3f max_error_deg=%.4f\n",
               b, bits, bits / 288.0, bound_deg);
        if (!(bits < 288.0)) table_ok = 0;   /* must beat whole-matrix storage */
    }
    check("beats_matrix_storage", table_ok);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
