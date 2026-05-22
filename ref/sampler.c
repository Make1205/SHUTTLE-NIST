/*
 * sampler.c - Discrete Gaussian sampler for SHUTTLE-NIST (ref version).
 *
 * Three modes selectable via SHUTTLE_SIGMA (see config.h / params.h):
 *
 *   sigma | small  |  k  | RCDT bits |  entries  | parameter set
 *   ------+--------+-----+-----------+-----------+---------------
 *    67   | 1.0469 |  64 |   93      |    11     | SHUTTLE-NIST-120
 *    79   | 1.2344 |  64 |   93      |    13     | SHUTTLE-NIST-180
 *    90   | 1.4062 |  64 |   93      |    15     | SHUTTLE-NIST-255
 *
 * All three modes share the same k=64, 93-bit RCDT (3 x 31-bit limbs)
 * layout. Output: signed int16_t samples with |r| <= 11*sigma:
 *   sigma=67 -> 737, sigma=79 -> 869, sigma=90 -> 990.
 *
 * Unified buffer design: one SHAKE-256 stream per sample_gauss_N call
 * produces ALL randomness, in repeating mini-batches.
 */

#include "sampler.h"
#include "approx_exp.h"
#include "approx_exp_constants.h"
#include <string.h>

/* ============================================================
 * RCDT table.
 *
 * Layout per entry: 3 limbs of 31 bits each (low, mid, upper). The high
 * bit of each limb is masked off in the comparison so the value lives in
 * a strict 31-bit range.
 *
 * Entry i encodes Pr[|X_small_sigma| >= i+1] * 2^93 for the mode's
 * small_sigma.
 * ============================================================ */

#if SHUTTLE_SIGMA == 67
/* SHUTTLE-NIST-120: small_sigma = 1.0469, 11 entries, 93-bit precision.
 * R_241 = 1 + 2^(-95.55), R_inf = 1 + 2^(-89.71). */
static const uint32_t RCDT_3x31[RCDT_ENTRIES][3] = {
    {0x4778B7BAU, 0x3AD4C2CBU, 0x395D0EFEU},
    {0x0587579EU, 0x46F607E2U, 0x0C9A2F23U},
    {0x27D78F61U, 0x09ADE3CAU, 0x01365D62U},
    {0x5BE1C32DU, 0x768A8E44U, 0x000C6D67U},
    {0x24EB781AU, 0x3D4D813CU, 0x000033E7U},
    {0x19E8EF98U, 0x48363C95U, 0x00000057U},
    {0x496FEE03U, 0x1DBADB68U, 0x00000000U},
    {0x7E997742U, 0x00081D40U, 0x00000000U},
    {0x5FEC7E81U, 0x000000E3U, 0x00000000U},
    {0x0503A6E6U, 0x00000000U, 0x00000000U},
    {0x000016B2U, 0x00000000U, 0x00000000U}
};
#elif SHUTTLE_SIGMA == 79
/* SHUTTLE-NIST-180: small_sigma = 1.2344, 13 entries, 93-bit precision.
 * R_361 = 1 + 2^(-93.76), R_inf = 1 + 2^(-89.06). */
static const uint32_t RCDT_3x31[RCDT_ENTRIES][3] = {
    {0x093D9B5BU, 0x43C2739EU, 0x4178E522U},
    {0x2B60EF0EU, 0x3E4D7A47U, 0x146F92FBU},
    {0x7B8A09CDU, 0x028ABCF5U, 0x039B90DAU},
    {0x5B345B31U, 0x7A31FDF8U, 0x00587AD4U},
    {0x218C0F73U, 0x3CDDB361U, 0x0004805AU},
    {0x2B78CE0CU, 0x5B98FA0DU, 0x00001EC8U},
    {0x6839942BU, 0x73ACA81DU, 0x0000006DU},
    {0x469FE738U, 0x661BE347U, 0x00000000U},
    {0x11ECAB98U, 0x0062989BU, 0x00000000U},
    {0x6FC7E575U, 0x0000316EU, 0x00000000U},
    {0x6E830A79U, 0x0000000CU, 0x00000000U},
    {0x00DE52E7U, 0x00000000U, 0x00000000U},
    {0x00000F92U, 0x00000000U, 0x00000000U}
};
#elif SHUTTLE_SIGMA == 90
/* SHUTTLE-NIST-255: small_sigma = 1.4062, 15 entries, 93-bit precision.
 * R_511 = 1 + 2^(-93.88), R_inf = 1 + 2^(-88.67). */
static const uint32_t RCDT_3x31[RCDT_ENTRIES][3] = {
    {0x7FA41C82U, 0x7E0897B9U, 0x476C5504U},
    {0x142B2B3CU, 0x1817F3ECU, 0x1B7CA058U},
    {0x2EC1E6CFU, 0x5C033B85U, 0x06E8F000U},
    {0x339216F2U, 0x0A968792U, 0x01192102U},
    {0x4BAA6424U, 0x1214A738U, 0x001BB4E8U},
    {0x3ADB1542U, 0x37B18CCEU, 0x0001ACA0U},
    {0x0EB480EFU, 0x63BC55E1U, 0x00000FC6U},
    {0x495FC076U, 0x17A34257U, 0x0000005AU},
    {0x41796E7CU, 0x1C0209E3U, 0x00000001U},
    {0x34D2DF92U, 0x01463540U, 0x00000000U},
    {0x2E7EB4B5U, 0x00019BE2U, 0x00000000U},
    {0x70C57416U, 0x00000139U, 0x00000000U},
    {0x482919B5U, 0x00000000U, 0x00000000U},
    {0x001403DEU, 0x00000000U, 0x00000000U},
    {0x00000358U, 0x00000000U, 0x00000000U}
};
#endif

/* ============================================================
 * sampler_sigma2: batched small-sigma CDT sampler.
 *
 * Two implementations behind an #if -- the 3x31 path masks each limb's
 * high bit so the comparison stays in a strict 31-bit range (matches
 * the AVX2 layout). The 3x32 path uses the full 32-bit limb (no mask).
 * ============================================================ */
int sampler_sigma2(int16_t *z_out, const uint8_t *rand) {
    for (int s = 0; s < GAUSS_BATCH; s++) {
        int group = s >> 3;
        int lane  = s & 7;
        int byte_base = group * 96 + lane * 4;

#if SHUTTLE_RCDT_LIMB_BITS == 31
        uint32_t v0 = load_le32(rand + byte_base)      & 0x7FFFFFFFU;
        uint32_t v1 = load_le32(rand + byte_base + 32) & 0x7FFFFFFFU;
        uint32_t v2 = load_le32(rand + byte_base + 64) & 0x7FFFFFFFU;

        int32_t z = 0;
        for (int i = 0; i < RCDT_ENTRIES; i++) {
            uint32_t cc;
            cc = ct_lt_u32(v0, RCDT_3x31[i][0]);
            cc = ct_lt_u32(v1 - cc, RCDT_3x31[i][1]);
            cc = ct_lt_u32(v2 - cc, RCDT_3x31[i][2]);
            z += (int32_t)cc;
        }
        z_out[s] = (int16_t)z;
#else
#  error "SHUTTLE-NIST requires SHUTTLE_RCDT_LIMB_BITS == 31"
#endif
    }
    return GAUSS_BATCH;
}

/* ============================================================
 * sampler_y: batched uniform Y_BITS-bit sampler.
 *
 * Produces GAUSS_BATCH = 16 uniform Y_BITS-bit values from exactly
 *   Y_RAND_BYTES = (GAUSS_BATCH * Y_BITS + 7) / 8 = ceil(16*Y_BITS/8)
 * input bytes, packed little-endian bit-wise. Implementation is a
 * generic accumulator loop -- portable for any Y_BITS up to ~24.
 * ============================================================ */
static void sampler_y(uint8_t y_out[GAUSS_BATCH], const uint8_t *rand) {
    const uint32_t y_mask = (1U << SHUTTLE_Y_BITS) - 1U;
    uint32_t acc = 0;
    unsigned int acc_bits = 0;
    unsigned int byte_idx = 0;
    for (int i = 0; i < GAUSS_BATCH; ++i) {
        while (acc_bits < (unsigned)SHUTTLE_Y_BITS) {
            acc |= (uint32_t)rand[byte_idx++] << acc_bits;
            acc_bits += 8;
        }
        y_out[i] = (uint8_t)(acc & y_mask);
        acc >>= SHUTTLE_Y_BITS;
        acc_bits -= SHUTTLE_Y_BITS;
    }
}

/* ============================================================
 * sample_gauss_attempt: one (x, y, sign_for_0, rej) attempt.
 *
 *   candidate = (x << K_BITS) | y                   in [0, 2^K_BITS * (RCDT_ENTRIES+1)]
 *   t         = y + (x << TWO_K_BITS)               = y + 2*k*x
 *   num       = y * t                                in uint32 (num < 2^20)
 *   a         = num / (2*sigma^2)
 *   accept iff rand_rej_63 < approx_exp(a_q60)       in Q63
 *
 * The input a_q60 = round(a * 2^60) to approx_exp is computed in
 * fixed-point with a Q80 reciprocal of 2*sigma^2:
 *
 *   I       = floor(2^80 / (2*sigma^2))                 ~ 64..66 bits
 *   M       = num * I    (full 128-bit integer multiply)
 *   a_q60   = M >> 20
 *
 * Error chain (relative to the ideal a):
 *   I underflow:           |2^80/N - I| < 1                   (Q80 ULP)
 *   M underflow:           |M_ideal - M| <= num <= 2^18       (Q80 ULP)
 *   shift truncation:      <= 1                               (Q60 ULP)
 *   total                  |a_q60 - a*2^60| <= 5/4 + 1 < 2   (Q60 ULP)
 *   propagation to exp(-a): <= 2 / 2^60 * exp(-a) = 2^-59 * exp(-a)
 *
 * For all three SHUTTLE-NIST modes (sigma in {67, 79, 90}), 2*sigma^2 is
 * not a power of two so the Q80 reciprocal path always applies. Max
 * num <= 63 * (63 + 64*RCDT_ENTRIES) <= 63 * (63 + 64*15) = 64449 << 2^20;
 * the 128-bit product M fits comfortably.
 * ============================================================ */
static int sample_gauss_attempt(int16_t *r,
                                uint32_t x,
                                uint32_t y,
                                const uint8_t *rej_rand,
                                const uint8_t *signs,
                                size_t idx) {
    int32_t candidate = ((int32_t)x << SHUTTLE_K_BITS) | (int32_t)y;

    uint64_t rand_tail = load_le64(rej_rand);
    uint8_t  sign_r0      = (uint8_t)(rand_tail & 1U);
    uint64_t rand_rej_63  = rand_tail >> 1;

    uint32_t t   = (uint32_t)y + ((uint32_t)x << SHUTTLE_TWO_K_BITS);
    uint64_t num = (uint64_t)y * (uint64_t)t;

    /* Q80 reciprocal path. Splitting I into (I_hi, I_lo) lets us cover
     * the 67-bit reciprocal of 2*sigma^2 with two 64-bit lanes.
     * The product is built as
     *
     *   M = num * I = num * I_lo  +  (num * I_hi) << 64
     *
     * `num` is bounded by 2^17 (see comment above) and `I_hi` < 2^3,
     * so `num * I_hi` always fits in a uint64 -- no nested 128-bit
     * multiplication is required. */
    wide_uint128 M_lo = (wide_uint128)num * APPROX_EXP_I_LO;
    uint64_t          M_hi = (uint64_t)num * APPROX_EXP_I_HI;
    wide_uint128 M = M_lo + ((wide_uint128)M_hi << 64);
    uint64_t a_q60 = (uint64_t)(M >> 20);

    uint64_t exp_val_q63 = approx_exp(a_q60);
    int accepted = (rand_rej_63 < exp_val_q63) ? 1 : 0;

    /* r=0 case: accept with prob 1/2 (folded-distribution correction). */
    uint64_t r_is_zero = (candidate == 0) ? 1ULL : 0ULL;
    accepted &= (int)(1 - (r_is_zero & (1 - (uint64_t)sign_r0)));

    uint8_t sign = (signs[idx / 8] >> (idx % 8)) & 1;
    int16_t cand16 = (int16_t)candidate;
    *r = sign ? (int16_t)(-cand16) : cand16;

    return accepted;
}

/* ============================================================
 * sample_gauss_N: unified SHAKE stream, mini-batch processing.
 * ============================================================ */
void sample_gauss_N(int16_t *r,
                    const uint8_t seed[SHUTTLE_SEEDBYTES],
                    uint64_t nonce, size_t len) {
    uint8_t buf[GAUSS_BUF_SIZE];
    stream256_state state;
    stream256_init(&state, seed, nonce);

    size_t sign_bytes = len / 8;
    size_t init_needed = sign_bytes + MINIBATCH_RAND_BYTES;
    size_t init_nblocks = (init_needed + STREAM256_BLOCKBYTES - 1)
                          / STREAM256_BLOCKBYTES;
    stream256_squeezeblocks(buf, init_nblocks, &state);

    /* len <= SHUTTLE_N = 256 across all NIST modes, so sign_bytes <= 32. */
    uint8_t signs[SHUTTLE_N / 8];
    memcpy(signs, buf, sign_bytes);
    size_t pos   = sign_bytes;
    size_t avail = init_nblocks * STREAM256_BLOCKBYTES - sign_bytes;

    size_t coefcnt = 0;
    int16_t z[GAUSS_BATCH];
    uint8_t y[GAUSS_BATCH];

    while (coefcnt < len) {
        if (avail < (size_t)MINIBATCH_RAND_BYTES) {
            if (avail > 0)
                memmove(buf, buf + pos, avail);
            pos = 0;
            stream256_squeezeblocks(buf + avail, SQUEEZE_NBLOCKS, &state);
            avail += SQUEEZE_NBLOCKS * STREAM256_BLOCKBYTES;
        }

        sampler_sigma2(z, buf + pos);
        pos   += SIGMA2_RAND_BYTES;
        avail -= SIGMA2_RAND_BYTES;

        sampler_y(y, buf + pos);
        pos   += Y_RAND_BYTES;
        avail -= Y_RAND_BYTES;

        for (int j = 0; j < GAUSS_BATCH && coefcnt < len; j++) {
            int accepted = sample_gauss_attempt(&r[coefcnt],
                                                (uint32_t)z[j],
                                                (uint32_t)y[j],
                                                buf + pos,
                                                signs, coefcnt);
            pos   += GAUSS_RAND_BYTES;
            avail -= GAUSS_RAND_BYTES;

            if (accepted)
                coefcnt++;
        }
    }
}
