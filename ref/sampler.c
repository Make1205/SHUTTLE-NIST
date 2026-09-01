/*
 * sampler.c - Discrete Gaussian sampler for SHUTTLE-NIST (ref version).
 *
 * Three modes selectable via SHUTTLE_SIGMA (see config.h / params.h):
 *
 *   sigma | small  |  k  | RCDT bits |  entries  | parameter set
 *   ------+--------+-----+-----------+-----------+---------------
 *    70   | 1.09375  |  64 |   93      |    12     | SHUTTLE-NIST-120
 *    82   | 1.28125  |  64 |   93      |    14     | SHUTTLE-NIST-180
 *    99   | 1.546875 |  64 |   93      |    17     | SHUTTLE-NIST-255
 *
 * All three modes share the same k=64, 93-bit RCDT (3 x 31-bit limbs)
 * layout. The finite table gives hard support |r| <= 831/959/1151; the
 * signature encoder separately enforces its configured 11*sigma vocabulary.
 *
 * Unified buffer design: one SHAKE-256 stream per sample_gauss_N call
 * produces ALL randomness, in repeating mini-batches.
 */

#include "sampler.h"
#include "approx_exp.h"
#include "approx_exp_constants.h"
#include "rcdt_tables.h"
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
            cc = ct_lt_u32(v1, RCDT_3x31[i][1] + cc);
            cc = ct_lt_u32(v2, RCDT_3x31[i][2] + cc);
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
 * For all three SHUTTLE-NIST modes (sigma in {70, 82, 99}), 2*sigma^2 is
 * not a power of two so the Q80 reciprocal path always applies. Max
 * num <= 63 * (63 + 64*RCDT_ENTRIES) <= 63 * (63 + 64*17) = 72513 << 2^20;
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
