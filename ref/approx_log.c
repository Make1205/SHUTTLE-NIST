/*
 * approx_log.c - Integer-only approximation of log2(U) for SHUTTLE IRS.
 *
 * The ordinary path consumes an 80-bit exponent block and a 57-bit mantissa
 * block (18 bytes in total). The exponent is geometric, obtained from the
 * number of leading zero bits. Additional 80-bit public-randomness blocks
 * are read lazily only after an all-zero block, up to the supported-range
 * cutoff. The mantissa is evaluated with the four-segment degree-13 Q62
 * polynomial specified in the paper. Keeping the exponent separate avoids
 * the dynamic-range overflow of a signed Q62 complete logarithm.
 */

#include "config.h"
#include "approx_log.h"

typedef __int128 wide_int128;

/* c[j][i] is the signed Q62 coefficient of x^i on segment j. */
static const int64_t LOG2_COEFFS[4][14] = {
    {
        INT64_C(0), INT64_C(1663314137230540005),
        INT64_C(-207914267153782834), INT64_C(34652377857794422),
        INT64_C(-6497320829889967), INT64_C(1299464001399016),
        INT64_C(-270720747985221), INT64_C(58008183734025),
        INT64_C(-12680634058615), INT64_C(2802591328231),
        INT64_C(-611697204761), INT64_C(123141961189),
        INT64_C(-19530076391), INT64_C(1717255301)
    },
    {
        INT64_C(1484631294131014398), INT64_C(1330651309784432187),
        INT64_C(-133065130978439107), INT64_C(17742017463685572),
        INT64_C(-2661302618120493), INT64_C(425808407352253),
        INT64_C(-70968007646487), INT64_C(12165731010280),
        INT64_C(-2128477788302), INT64_C(377485364884),
        INT64_C(-66839897685), INT64_C(11228353229),
        INT64_C(-1554229086), INT64_C(124985134)
    },
    {
        INT64_C(2697663385880076776), INT64_C(1108876091487026868),
        INT64_C(-92406340957251847), INT64_C(10267371217462318),
        INT64_C(-1283421402046772), INT64_C(171122852511038),
        INT64_C(-23767057140815), INT64_C(3395273701197),
        INT64_C(-495094418506), INT64_C(73261264270),
        INT64_C(-10884612082), INT64_C(1561954161),
        INT64_C(-191070343), INT64_C(14130118)
    },
    {
        INT64_C(3723267405961586381), INT64_C(950465221274594463),
        INT64_C(-67890372948185266), INT64_C(6465749804587722),
        INT64_C(-692758907616243), INT64_C(79172446438655),
        INT64_C(-9425290482070), INT64_C(1154114515161),
        INT64_C(-144257705268), INT64_C(18306992888),
        INT64_C(-2339854304), INT64_C(292286931),
        INT64_C(-31958889), INT64_C(2187745)
    }
};

/* Round-half-up division by 2^62. */
static inline int64_t round_shift_q62(wide_int128 value) {
    return (int64_t)((value + ((wide_int128)1 << 61)) >> 62);
}

static uint64_t load_be64(const uint8_t src[8]) {
    uint64_t value = 0;
    for (unsigned int i = 0; i < 8; ++i)
        value = (value << 8) | src[i];
    return value;
}

/* Return the leading-zero count of an 80-bit big-endian block. */
static unsigned int clz80(const uint8_t src[10]) {
    uint64_t hi = load_be64(src);
    uint32_t lo = ((uint32_t)src[8] << 8) | src[9];
    uint64_t hi_nz = (hi | (UINT64_C(0) - hi)) >> 63;
    uint32_t lo_nz = (lo | (UINT32_C(0) - lo)) >> 31;

    /* OR-with-one makes both CLZ calls defined on an all-zero limb. */
    unsigned int lz_hi = (unsigned int)__builtin_clzll(hi | UINT64_C(1))
                       + (unsigned int)(1U - hi_nz);
    unsigned int lz_lo = (unsigned int)__builtin_clz(lo | UINT32_C(1)) - 16U
                       + (unsigned int)(1U - lo_nz);
    uint32_t hi_mask = UINT32_C(0) - (uint32_t)hi_nz;

    return (lz_hi & hi_mask) | ((64U + lz_lo) & ~hi_mask);
}

/* Select one coefficient row with the input-independent masked scan used in
 * the paper's recommended implementation. */
static void select_coefficients(int64_t selected[14], uint32_t segment) {
    for (unsigned int i = 0; i < 14; ++i) {
        int64_t value = 0;
        for (uint32_t row = 0; row < 4; ++row) {
            uint64_t diff = (uint64_t)(row ^ segment);
            uint64_t nz = (diff | (UINT64_C(0) - diff)) >> 63;
            int64_t mask = (int64_t)(nz - 1U);
            value |= LOG2_COEFFS[row][i] & mask;
        }
        selected[i] = value;
    }
}

void approx_log2_sample(shuttle_log2_sample *out,
                        const uint8_t random[SHUTTLE_LOG_RANDOM_BYTES],
                        keccak_state *public_rng_state) {
    int64_t coeffs[14];
    unsigned int leading_zeros = clz80(random);
    unsigned int consumed = SHUTTLE_LOG_EXPONENT_BLOCK_BITS;

    out->exponent = leading_zeros + 1U;
    while (leading_zeros == SHUTTLE_LOG_EXPONENT_BLOCK_BITS
           && consumed < SHUTTLE_LOG_EXPONENT_BITS) {
        uint8_t extension[SHUTTLE_LOG_EXPONENT_BLOCK_BYTES];
        shake256_squeeze(extension, sizeof extension, public_rng_state);
        leading_zeros = clz80(extension);
        out->exponent = consumed + leading_zeros + 1U;
        consumed += SHUTTLE_LOG_EXPONENT_BLOCK_BITS;
    }

    /* The next eight bytes are read MSB-first; the low seven padding bits are
     * discarded, leaving a uniform 57-bit mantissa integer. */
    uint64_t mantissa =
        load_be64(random + SHUTTLE_LOG_EXPONENT_BLOCK_BYTES) >> 7;
    uint32_t segment = (uint32_t)(mantissa >> 55);
    uint64_t local = mantissa & ((UINT64_C(1) << 55) - 1U);
    int64_t x_q62 = (int64_t)(local << 7);
    select_coefficients(coeffs, segment);

    int64_t z = coeffs[13];
    for (int i = 12; i >= 0; --i)
        z = coeffs[i] + round_shift_q62((wide_int128)z * x_q62);

    out->fraction_q62 = z;
}
