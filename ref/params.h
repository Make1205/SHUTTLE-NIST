/*
 * params.h - Parameters for SHUTTLE-NIST signature scheme.
 *
 * Three parameter sets:
 *
 *   SHUTTLE_NIST_MODE=120 : n=256, q=32257, sigma=70, full NTT (base_deg=1)
 *   SHUTTLE_NIST_MODE=180 : n=256, q=15361, sigma=82, full NTT (base_deg=1)
 *   SHUTTLE_NIST_MODE=255 : n=256, q=7681,  sigma=99, full NTT (base_deg=1)
 *
 * HighBits(z) and hint coefficients are encoded by one merged rANS
 * context. Its frequency table is generated from the coefficient-count-
 * weighted mixture of the rounded-response and bucket-crossing PMFs.
 *
 * All three q satisfy q == 1 (mod 2n=512), so full NTT applies in every
 * mode.
 *
 * SHUTTLE_NIST_MODE is selected via config.h.
 */

#ifndef SHUTTLE_PARAMS_H
#define SHUTTLE_PARAMS_H

#include <stdint.h>
#include "config.h"

/* ============================================================
 * Shared constants (hash / seed sizes)
 * ============================================================ */
#define SHUTTLE_SEEDBYTES    32
#define SHUTTLE_SKSEEDBYTES  64
#define SHUTTLE_CRHBYTES     64
#define SHUTTLE_TRBYTES      64
#define SHUTTLE_RNDBYTES     32

/* Challenge hash seedC: 32/48/64 B for 120/180/255 (collision-resistance
 * matches classical security level). */
#if SHUTTLE_NIST_MODE == 120
#  define SHUTTLE_CTILDEBYTES  32
#elif SHUTTLE_NIST_MODE == 180
#  define SHUTTLE_CTILDEBYTES  48
#elif SHUTTLE_NIST_MODE == 255
#  define SHUTTLE_CTILDEBYTES  64
#endif

/* ============================================================
 * Mode-specific parameters
 * ============================================================ */
#if SHUTTLE_NIST_MODE == 120

#define SHUTTLE_N          256
#define SHUTTLE_Q          32257
#define SHUTTLE_QBITS      15            /* ceil(log2(32257)) = 15 */
#define SHUTTLE_BASE_DEG   1             /* full NTT, q == 1 mod 2n */
#define SHUTTLE_L          3             /* lenS */
#define SHUTTLE_M          2             /* lenE */
#define SHUTTLE_ETA        1
#define SHUTTLE_TAU        58
#define SHUTTLE_IRS_N      75
#define SHUTTLE_LOG_EXPONENT_REQUIRED 20940
#define SHUTTLE_LOG_SCALE_HI 6792
#define SHUTTLE_LOG_SCALE_LO UINT64_C(3884743587687722514)
#define SHUTTLE_ALPHA_H    128
#define SHUTTLE_ALPHA_1    8

/*   B_k = 25.79 -> 2579,   2579^2     = 6651241
 *   B_s = 3233  -> 323300, 323300^2   = 104522890000
 *   B_v = 4000  -> 400000, 400000^2   = 160000000000  */
#define SHUTTLE_BK_Q2          2579L
#define SHUTTLE_BS_Q2          323300L
#define SHUTTLE_BV_Q2          400000L
#define SHUTTLE_BK_LOW_Q2      0L

#elif SHUTTLE_NIST_MODE == 180

#define SHUTTLE_N          256
#define SHUTTLE_Q          15361
#define SHUTTLE_QBITS      14            /* ceil(log2(15361)) = 14 */
#define SHUTTLE_BASE_DEG   1             /* full NTT, q == 1 mod 2n */
#define SHUTTLE_L          4
#define SHUTTLE_M          3
#define SHUTTLE_ETA        1
#define SHUTTLE_TAU        80
#define SHUTTLE_IRS_N      89
#define SHUTTLE_LOG_EXPONENT_REQUIRED 28383
#define SHUTTLE_LOG_SCALE_HI 9321
#define SHUTTLE_LOG_SCALE_LO UINT64_C(2044287409645675025)
#define SHUTTLE_ALPHA_H    128
#define SHUTTLE_ALPHA_1    8

/*   B_k = 30.4  -> 3040,    3040^2    = 9241600
 *   B_s = 4312  -> 431200,  431200^2  = 185933440000
 *   B_v = 5250  -> 525000,  525000^2  = 275625000000 */
#define SHUTTLE_BK_Q2          3040L
#define SHUTTLE_BS_Q2          431200L
#define SHUTTLE_BV_Q2          525000L
#define SHUTTLE_BK_LOW_Q2      0L

#elif SHUTTLE_NIST_MODE == 255

#define SHUTTLE_N          256
#define SHUTTLE_Q          7681
#define SHUTTLE_QBITS      13            /* ceil(log2(7681)) = 13 */
#define SHUTTLE_BASE_DEG   1             /* full NTT, q == 1 mod 2n */
#define SHUTTLE_L          5
#define SHUTTLE_M          4
#define SHUTTLE_ETA        1
#define SHUTTLE_TAU        128
#define SHUTTLE_IRS_N      55
#define SHUTTLE_LOG_EXPONENT_REQUIRED 19690
#define SHUTTLE_LOG_SCALE_HI 13587
#define SHUTTLE_LOG_SCALE_LO UINT64_C(327583442694609146)
#define SHUTTLE_ALPHA_H    128
#define SHUTTLE_ALPHA_1    16

/*   B_k = 34.4  -> 3440,   3440^2    = 11833600
 *   B_s = 5556  -> 555600, 555600^2  = 308691360000
 *   B_v = 6638  -> 663800, 663800^2  = 440630440000 */
#define SHUTTLE_BK_Q2          3440L
#define SHUTTLE_BS_Q2          555600L
#define SHUTTLE_BV_Q2          663800L
#define SHUTTLE_BK_LOW_Q2      0L

#else
#  error "Unsupported SHUTTLE_NIST_MODE (expected 120, 180, or 255)"
#endif

#if SHUTTLE_LOG_EXPONENT_BITS < SHUTTLE_LOG_EXPONENT_REQUIRED
#  error "Logarithm exponent support is too short for this parameter set"
#endif

/* Pull in MONT / QINV / BARRETT_V / INVNTT_F / zetas table per mode. */
#include "ntt_constants.h"

/* ============================================================
 * Fixed-point bound scale.
 *
 * Bounds B_k, B_s, B_v are decimals (e.g. B_k = 25.79). Stored as
 * integers in Q2 (× 100); compare a regular int64 norm_sq against the
 * squared Q2 value by scaling norm_sq by SHUTTLE_BOUND_SCALE_SQ:
 *
 *   accept iff  SHUTTLE_BOUND_SCALE_SQ * norm_sq < SHUTTLE_*_Q2 * SHUTTLE_*_Q2
 *
 * Worst-case norm_sq <= B_v^2 ~ 4.4e8; * 10000 = 4.4e12 (fits in int64).
 * ============================================================ */
#define SHUTTLE_BOUND_SCALE    100L
#define SHUTTLE_BOUND_SCALE_SQ (SHUTTLE_BOUND_SCALE * SHUTTLE_BOUND_SCALE)

#define SHUTTLE_BK_SQ_FX     ((int64_t)SHUTTLE_BK_Q2     * (int64_t)SHUTTLE_BK_Q2)
#define SHUTTLE_BS_SQ_FX     ((int64_t)SHUTTLE_BS_Q2     * (int64_t)SHUTTLE_BS_Q2)
#define SHUTTLE_BV_SQ_FX     ((int64_t)SHUTTLE_BV_Q2     * (int64_t)SHUTTLE_BV_Q2)
#define SHUTTLE_BK_LOW_SQ_FX ((int64_t)SHUTTLE_BK_LOW_Q2 * (int64_t)SHUTTLE_BK_LOW_Q2)

#define SHUTTLE_NORM_LT_FX(norm_sq, bound_sq_fx) \
    ((int64_t)(norm_sq) * SHUTTLE_BOUND_SCALE_SQ < (bound_sq_fx))

#define SHUTTLE_NORM_GE_FX(norm_sq, bound_sq_fx) \
    ((int64_t)(norm_sq) * SHUTTLE_BOUND_SCALE_SQ >= (bound_sq_fx))

#define SHUTTLE_NORM_GT_FX(norm_sq, bound_sq_fx) \
    ((int64_t)(norm_sq) * SHUTTLE_BOUND_SCALE_SQ > (bound_sq_fx))

/* Full sk vector length = [alpha_1, s, e] */
#define SHUTTLE_VECLEN     (1 + SHUTTLE_L + SHUTTLE_M)

/* High-bits bucket count.
 *   W1_MAX = (q-1) / (2*alpha_h):
 *     120: 32256/256 = 126 -> 127 buckets [0,126], 7 bits
 *     180: 15360/256 =  60 ->  61 buckets [0, 60], 6 bits
 *     255:  7680/256 =  30 ->  31 buckets [0, 30], 5 bits
 *
 * The +1 rounding-edge bucket wraps modulo (W1_MAX+1) per rounding.c, so
 * W1_BITS just needs ceil(log2(W1_MAX+1)). */
#define SHUTTLE_W1_MAX     ((SHUTTLE_Q - 1) / (2 * SHUTTLE_ALPHA_H))

#if SHUTTLE_NIST_MODE == 120
#  define SHUTTLE_W1_BITS  7        /* W1_MAX = 126, 127 buckets fit in 7 bits */
#elif SHUTTLE_NIST_MODE == 180
#  define SHUTTLE_W1_BITS  6        /* W1_MAX = 60, 61 buckets fit in 6 bits */
#elif SHUTTLE_NIST_MODE == 255
#  define SHUTTLE_W1_BITS  5        /* W1_MAX = 30, 31 buckets fit in 5 bits */
#endif

/* ============================================================
 * mod 2q infrastructure
 * ============================================================ */
#define SHUTTLE_DQ          (2 * SHUTTLE_Q)
#define SHUTTLE_HALF_ALPHA_H (SHUTTLE_ALPHA_H / 2)
#define SHUTTLE_HINT_MOD    (2 * (SHUTTLE_Q - 1))
#define SHUTTLE_HINT_MAX    (SHUTTLE_HINT_MOD / SHUTTLE_ALPHA_H)

/* HighBits are hashed in the paper's scaled representation.  Every bucket
 * index is multiplied by alpha_h and serialized as a canonical little-endian
 * uint16_t, which is wide enough for all values modulo 2(q-1). */
#define SHUTTLE_WHASH_COEFF_BYTES 2
#define SHUTTLE_WHASH_PACKEDBYTES \
    (SHUTTLE_M * SHUTTLE_N * SHUTTLE_WHASH_COEFF_BYTES)

#if SHUTTLE_ALPHA_H == 128
#  define SHUTTLE_ALPHA_H_BITS 7
#else
#  error "Unsupported SHUTTLE_ALPHA_H (expected 128)"
#endif

#if SHUTTLE_ALPHA_1 == 8
#  define SHUTTLE_ALPHA_1_BITS 3
#elif SHUTTLE_ALPHA_1 == 16
#  define SHUTTLE_ALPHA_1_BITS 4
#else
#  error "Unsupported SHUTTLE_ALPHA_1 (expected 8 or 16)"
#endif

/* ============================================================
 * rANS hi/lo split parameters.
 *
 *   alpha_r  = 64 (fixed across modes): the rANS hi/lo split step for
 *              z^(1..lenS). HighBits live in the rANS vocabulary at
 *              scale r/alpha_r; LowBits are uniform in [-32, 32) and
 *              bit-packed at 6 bits/coef.
 *   alpha_0' = alpha_r / alpha_1: second-level split for z^(0). LowBits are
 *              bit-packed at log2(alpha_0') bits/coef.
 *     -> mode 120: alpha_0' = 64/8  = 8 (3 bits)
 *     -> mode 180: alpha_0' = 64/8  = 8 (3 bits)
 *     -> mode 255: alpha_0' = 64/16 = 4 (2 bits)
 *
 *   alpha_h = 2*alpha_r = 128 across all modes: the hint and z-hi scales
 *              are comparable. Their distributions are not identical, so
 *              the single-stream table uses their coefficient-count-
 *              weighted mixture.
 * ============================================================ */
#define SHUTTLE_ALPHA_R       64
#define SHUTTLE_ALPHA_R_BITS  6
#define SHUTTLE_HALF_ALPHA_R  (SHUTTLE_ALPHA_R / 2)

#define SHUTTLE_ALPHA_0P      (SHUTTLE_ALPHA_R / SHUTTLE_ALPHA_1)
#if   SHUTTLE_ALPHA_0P == 8
#  define SHUTTLE_ALPHA_0P_BITS 3
#elif SHUTTLE_ALPHA_0P == 4
#  define SHUTTLE_ALPHA_0P_BITS 2
#else
#  error "Unsupported SHUTTLE_ALPHA_0P (expected 4 or 8)"
#endif
#define SHUTTLE_HALF_ALPHA_0P (SHUTTLE_ALPHA_0P / 2)

/* ============================================================
 * Derived packing sizes (bytes)
 * ============================================================ */
/* eta=1: coefficients in {-1,0,1}, encode as 2 bits/coeff */
#define SHUTTLE_POLYETA_PACKEDBYTES   (SHUTTLE_N * 2 / 8)

/* Public key b: QBITS bits/coeff */
#define SHUTTLE_POLYPK_PACKEDBYTES    ((SHUTTLE_N * SHUTTLE_QBITS + 7) / 8)

/* z[0] after CompressY: |z[0]|_inf <= q/(2*alpha_1):
 *   mode-120: 32257/16 = 2016, 12 signed bits
 *   mode-180: 15361/16 =  960, 11 signed bits
 *   mode-255: 7681/32  =  240,  9 signed bits */
#if SHUTTLE_NIST_MODE == 120
#  define SHUTTLE_Z0_BITS  12
#elif SHUTTLE_NIST_MODE == 180
#  define SHUTTLE_Z0_BITS  11
#elif SHUTTLE_NIST_MODE == 255
#  define SHUTTLE_Z0_BITS  9
#endif
#define SHUTTLE_POLYZ0_PACKEDBYTES    ((SHUTTLE_N * SHUTTLE_Z0_BITS + 7) / 8)

/* z[1..L]: full-range signed coefficients. Z_BOUND <= 11*sigma + alpha_1*tau.
 *   mode-120: 11*70 + 8*58  = 1234, 12 signed bits
 *   mode-180: 11*82 + 8*80  = 1542, 12 signed bits
 *   mode-255: 11*99 + 16*128 = 3137, 13 signed bits
 * Use 13 bits uniformly. */
#define SHUTTLE_POLYZ_BITS  13
#define SHUTTLE_POLYZ_PACKEDBYTES     ((SHUTTLE_N * SHUTTLE_POLYZ_BITS + 7) / 8)

/* High-bits (w1) packing */
#define SHUTTLE_POLYW1_PACKEDBYTES    ((SHUTTLE_N * SHUTTLE_W1_BITS + 7) / 8)

/* Hint encoding (sparse-index list) -- placeholder; see polyveck_hint_pack_basic */
#define SHUTTLE_POLYVECH_PACKEDBYTES  (SHUTTLE_M * (SHUTTLE_N / 8 + 1))

/* ============================================================
 * rANS reservation budget for the merged single-stream design (design S).
 *
 * One rANS context encodes n*(lenS+1+lenE) symbols (HighBits of z^(0..L)
 * followed by hint). Its table and rate calculation use the analytic
 * mixture described above. The reservation is sized via a Gaussian
 * codelength approximation at p_rans^* = 2^-30. Numbers below are produced
 * by tools/SigSize.py. Regenerate with:
 *   python3 ../tools/SigSize.py
 * ============================================================ */
#if SHUTTLE_NIST_MODE == 120
#  define SHUTTLE_RANS_RESERVED_BYTES  469
#elif SHUTTLE_NIST_MODE == 180
#  define SHUTTLE_RANS_RESERVED_BYTES  675
#elif SHUTTLE_NIST_MODE == 255
#  define SHUTTLE_RANS_RESERVED_BYTES  923
#endif

#define SHUTTLE_RANS_BLOCK_BYTES     (2 + SHUTTLE_RANS_RESERVED_BYTES)

/* Packed size of the LowBits part of z^(0): ALPHA_0P_BITS per coef. */
#define SHUTTLE_POLYZ0_LO_PACKEDBYTES ((SHUTTLE_N * SHUTTLE_ALPHA_0P_BITS + 7) / 8)

/* Packed size of the LowBits part of z[1..lenS]: ALPHA_R_BITS per coef. */
#define SHUTTLE_POLYZ1_LO_PACKEDBYTES ((SHUTTLE_N * SHUTTLE_ALPHA_R_BITS + 7) / 8)

/* ============================================================
 * Public / secret key + signature sizes
 * ============================================================ */
#define SHUTTLE_PUBLICKEYBYTES  (SHUTTLE_SEEDBYTES \
                                   + SHUTTLE_M * SHUTTLE_POLYPK_PACKEDBYTES)

#define SHUTTLE_SECRETKEYBYTES  (SHUTTLE_SEEDBYTES \
                                   + SHUTTLE_TRBYTES \
                                   + SHUTTLE_SEEDBYTES \
                                   + SHUTTLE_L * SHUTTLE_POLYETA_PACKEDBYTES \
                                   + SHUTTLE_M * SHUTTLE_POLYETA_PACKEDBYTES)

/* Signature layout (single merged rANS stream):
 *   seedC || uint16 rans_len
 *   || rANS(z-hi  ⨁  hint) + pad to RANS_RESERVED
 *   || polyz0_lo_pack(lo(z^(0)))
 *   || L * polyz1_lo_pack(lo(z^(1..lenS))) */
#define SHUTTLE_BYTES  ( SHUTTLE_CTILDEBYTES \
                       + SHUTTLE_RANS_BLOCK_BYTES \
                       + SHUTTLE_POLYZ0_LO_PACKEDBYTES \
                       + SHUTTLE_L * SHUTTLE_POLYZ1_LO_PACKEDBYTES )

/* ============================================================
 * Gaussian sampler parameters (driven by SHUTTLE_SIGMA, set in config.h)
 *
 * All three SHUTTLE-NIST levels use BLISS-style sampling with k = 64 and
 * a 93-bit RCDT (3 x 31-bit limbs):
 *
 *   sigma | sigma_2 = sigma/k | RCDT entries
 *   ------+-------------------+---------------
 *    70   |   1.09375         |   12   (SHUTTLE-NIST-120)
 *    82   |   1.28125         |   14   (SHUTTLE-NIST-180)
 *    99   |   1.546875        |   17   (SHUTTLE-NIST-255)
 * ============================================================ */

#if SHUTTLE_SIGMA == 70
#  define RCDT_ENTRIES 12
#elif SHUTTLE_SIGMA == 82
#  define RCDT_ENTRIES 14
#elif SHUTTLE_SIGMA == 99
#  define RCDT_ENTRIES 17
#else
#  error "Unsupported SHUTTLE_SIGMA (expected 70, 82, or 99)"
#endif

#define SHUTTLE_RCDT_LIMB_BITS 31
#define SHUTTLE_RCDT_BITS      93

#define SHUTTLE_GAUSS_K        64
#define SHUTTLE_K_BITS         6
#define SHUTTLE_Y_BITS         6
#define SHUTTLE_TWO_K_BITS     7

/* Hard support of the convolution sampler before the signature encoder's
 * configured 11*sigma vocabulary check. */
#define SHUTTLE_GAUSS_MAX \
    (RCDT_ENTRIES * SHUTTLE_GAUSS_K + ((1 << SHUTTLE_Y_BITS) - 1))

#define SHUTTLE_TRUNC        11
#define SHUTTLE_BOUND        (SHUTTLE_TRUNC * SHUTTLE_SIGMA)

/* Packing bound for z coefficients (post alpha_1 stretch):
 *   |z[0]|_inf <= 11*sigma + alpha_1*tau   (after Stretch + IRS update)
 *   |z[i]|_inf <= 11*sigma + tau           for i >= 1
 * Use the larger of the two as a uniform bound. */
#define SHUTTLE_Z_BOUND      (SHUTTLE_BOUND + SHUTTLE_ALPHA_1 * SHUTTLE_TAU)

#endif /* SHUTTLE_PARAMS_H */
