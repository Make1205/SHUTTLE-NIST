/*
 * rejsample.c - Rejection-free iterative sampler for SHUTTLE.
 *
 * Each transition uses the paper's logarithmic interval test. A single
 * inner product and the (shared) squared shift norm determine all interval
 * boundaries. The logarithm is represented as an integer exponent plus a
 * Q62 mantissa, and every scaled comparison is performed in signed 128-bit
 * arithmetic. This avoids both reciprocal-rounding error and the wraparound
 * that occurs when a complete logarithm is forced into signed Q62.
 */

#include "params.h"
#include "rejsample.h"

typedef __int128 wide_int128;

static inline int64_t round_shift_q62(wide_int128 value) {
    return (int64_t)((value + ((wide_int128)1 << 61)) >> 62);
}

static void poly_cyclic_shift(poly *v, const poly *p, unsigned int j) {
    for (unsigned int k = 0; k < SHUTTLE_N; ++k) {
        if (k >= j)
            v->coeffs[k] = p->coeffs[k - j];
        else
            v->coeffs[k] = -p->coeffs[SHUTTLE_N + k - j];
    }
}

static void polyvec_cyclic_shift(polyvec *v,
                                 const polyvec *sk,
                                 unsigned int j) {
    for (unsigned int i = 0; i < SHUTTLE_VECLEN; ++i)
        poly_cyclic_shift(&v->vec[i], &sk->vec[i], j);
}

static int64_t polyvec_inner_product(const polyvec *z, const polyvec *v) {
    int64_t sum = 0;
    for (unsigned int i = 0; i < SHUTTLE_VECLEN; ++i)
        for (unsigned int k = 0; k < SHUTTLE_N; ++k)
            sum += (int64_t)z->vec[i].coeffs[k] * v->vec[i].coeffs[k];
    return sum;
}

/* Convert log2(U)=frac/2^62-a into 2r^2*ln(U), still in Q62.
 * The mode-specific scale is C_r=round(2^62*2r^2*ln(2)), split at bit 62
 * so that the fractional product stays within signed 128-bit arithmetic. */
static wide_int128 scaled_log_q62(const shuttle_log2_sample *sample) {
    wide_int128 scale = ((wide_int128)SHUTTLE_LOG_SCALE_HI << 62)
                      + SHUTTLE_LOG_SCALE_LO;
    wide_int128 fractional =
        (wide_int128)SHUTTLE_LOG_SCALE_HI * sample->fraction_q62
        + round_shift_q62((wide_int128)SHUTTLE_LOG_SCALE_LO
                          * sample->fraction_q62);
    return fractional - (wide_int128)sample->exponent * scale;
}

int irs_update_sign(int64_t inner_product,
                    int64_t shift_norm_sq,
                    const shuttle_log2_sample *log_sample) {
    /* Normalize v to make <y,v> nonnegative. The final orientation maps the
     * selected normalized shift back to the caller's original v. */
    int64_t negative_mask = inner_product >> 63;
    uint64_t t_abs = ((uint64_t)(inner_product ^ negative_mask)
                     - (uint64_t)negative_mask);
    int32_t orientation = (int32_t)(1 + 2 * negative_mask);
    wide_int128 lambda_q62 = scaled_log_q62(log_sample);
    uint32_t in_union = 0;
    wide_int128 even_magnitude = 0;
    wide_int128 increment = (wide_int128)2 * t_abs + shift_norm_sq;
    wide_int128 twice_norm = (wide_int128)2 * shift_norm_sq;

    /* N is odd for all parameter sets. Pair i represents the interval
     * [A_{2i+1}, A_{2i}), so the largest tested shift magnitude is N.
     * The squared-norm differences are advanced by the recurrence
     * m(k+1)-m(k)=2|<y,v>|+(2k+1)||v||^2. */
    for (unsigned int i = 0; i < (SHUTTLE_IRS_N + 1U) / 2U; ++i) {
        wide_int128 lower_magnitude = even_magnitude + increment;
        wide_int128 lower_q62 = -(lower_magnitude << 62);
        wide_int128 upper_q62 = -(even_magnitude << 62);

        uint32_t inside = (uint32_t)((lower_q62 <= lambda_q62)
                                   & (lambda_q62 < upper_q62));
        in_union |= inside;

        even_magnitude = lower_magnitude + increment + twice_norm;
        increment += twice_norm + twice_norm;
    }

    /* The paper's p_v is the probability of subtracting the normalized
     * shift: membership in the interval union therefore selects -v. */
    return orientation * (1 - 2 * (int32_t)in_union);
}

static void extract_challenge(unsigned int positions[SHUTTLE_TAU],
                              int32_t signs[SHUTTLE_TAU],
                              const poly *c) {
    unsigned int count = 0;
    for (unsigned int k = 0; k < SHUTTLE_N && count < SHUTTLE_TAU; ++k) {
        if (c->coeffs[k] != 0) {
            positions[count] = k;
            signs[count] = c->coeffs[k];
            ++count;
        }
    }
}

int irs_sign(polyvec *z,
             const poly *c,
             const polyvec *sk_stretched,
             keccak_state *rng_state) {
    unsigned int positions[SHUTTLE_TAU];
    int32_t challenge_signs[SHUTTLE_TAU];
    uint8_t random[SHUTTLE_TAU * SHUTTLE_LOG_RANDOM_BYTES];
    int64_t shift_norm_sq = polyvec_sq_norm(sk_stretched);
    polyvec v;

    extract_challenge(positions, challenge_signs, c);
    shake256_squeeze(random, sizeof random, rng_state);

    for (unsigned int step = 0; step < SHUTTLE_TAU; ++step) {
        polyvec_cyclic_shift(&v, sk_stretched, positions[step]);

        /* This path also supports signed challenges, although the current
         * challenge generator emits coefficients in {0,1}. */
        int32_t c_sign = challenge_signs[step];
        for (unsigned int i = 0; i < SHUTTLE_VECLEN; ++i)
            for (unsigned int k = 0; k < SHUTTLE_N; ++k)
                v.vec[i].coeffs[k] *= c_sign;

        shuttle_log2_sample log_sample;
        approx_log2_sample(&log_sample,
                           random + step * SHUTTLE_LOG_RANDOM_BYTES,
                           rng_state);
        int32_t update_sign = irs_update_sign(
            polyvec_inner_product(z, &v), shift_norm_sq, &log_sample);

        for (unsigned int i = 0; i < SHUTTLE_VECLEN; ++i)
            for (unsigned int k = 0; k < SHUTTLE_N; ++k)
                z->vec[i].coeffs[k] += update_sign * v.vec[i].coeffs[k];
    }

    return 1;
}
