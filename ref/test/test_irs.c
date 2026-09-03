/* Arithmetic regression tests for the logarithmic-domain IRS transition. */

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "../approx_log.h"
#include "../params.h"
#include "../rejsample.h"

#define PASS(msg) do { printf("  [PASS] " msg "\n"); } while (0)
#define FAIL(msg, ...) do { \
    printf("  [FAIL] " msg "\n", __VA_ARGS__); \
    return 1; \
} while (0)

static void set_exponent_bit(uint8_t random[SHUTTLE_LOG_RANDOM_BYTES],
                             unsigned int leading_zeros) {
    memset(random, 0, SHUTTLE_LOG_RANDOM_BYTES);
    if (leading_zeros < SHUTTLE_LOG_EXPONENT_BLOCK_BITS)
        random[leading_zeros / 8U] =
            (uint8_t)(UINT32_C(0x80) >> (leading_zeros % 8U));
}

static void set_mantissa(uint8_t random[SHUTTLE_LOG_RANDOM_BYTES],
                         uint64_t mantissa) {
    uint64_t encoded = mantissa << 7;
    for (unsigned int i = 0; i < 8; ++i)
        random[10U + i] = (uint8_t)(encoded >> (56U - 8U * i));
}

static int test_exponent_layout(void) {
    uint8_t random[SHUTTLE_LOG_RANDOM_BYTES];
    shuttle_log2_sample sample;

    for (unsigned int leading_zeros = 0;
         leading_zeros < SHUTTLE_LOG_EXPONENT_BLOCK_BITS;
         ++leading_zeros) {
        set_exponent_bit(random, leading_zeros);
        approx_log2_sample(&sample, random, NULL);
        if (sample.exponent != leading_zeros + 1U)
            FAIL("exponent at leading-zero count %u is %u, expected %u",
                 leading_zeros, sample.exponent, leading_zeros + 1U);
        if (sample.fraction_q62 != 0)
            FAIL("zero mantissa produced fraction %lld, expected 0",
                 (long long)sample.fraction_q62);
    }

    {
        static const uint8_t seed[] = "SHUTTLE lazy exponent test";
        uint8_t extension[SHUTTLE_LOG_EXPONENT_BLOCK_BYTES];
        keccak_state state;
        memset(random, 0, sizeof random);
        shake256(extension, sizeof extension, seed, sizeof seed - 1U);
        shake256_init(&state);
        shake256_absorb(&state, seed, sizeof seed - 1U);
        shake256_finalize(&state);
        approx_log2_sample(&sample, random, &state);
        unsigned int extension_zeros = 0;
        while (extension_zeros < SHUTTLE_LOG_EXPONENT_BLOCK_BITS
               && !(extension[extension_zeros / 8U]
                    & (UINT32_C(0x80) >> (extension_zeros % 8U))))
            ++extension_zeros;
        if (sample.exponent != SHUTTLE_LOG_EXPONENT_BLOCK_BITS
                               + extension_zeros + 1U)
            FAIL("lazy exponent extension is %u, expected %u",
                 sample.exponent, SHUTTLE_LOG_EXPONENT_BLOCK_BITS
                                  + extension_zeros + 1U);
    }

    PASS("geometric exponent and lazy 80-bit extension decode exactly");
    return 0;
}

static int test_log_polynomial(void) {
    static const int64_t segment_values[4] = {
        INT64_C(0),
        INT64_C(1484631294131014398),
        INT64_C(2697663385880076776),
        INT64_C(3723267405961586381)
    };
    uint8_t random[SHUTTLE_LOG_RANDOM_BYTES];
    shuttle_log2_sample sample;

    set_exponent_bit(random, 0);
    for (uint64_t segment = 0; segment < 4; ++segment) {
        set_mantissa(random, segment << 55);
        approx_log2_sample(&sample, random, NULL);
        if (sample.fraction_q62 != segment_values[segment])
            FAIL("segment %llu starts at %lld, expected %lld",
                 (unsigned long long)segment,
                 (long long)sample.fraction_q62,
                 (long long)segment_values[segment]);
    }

    /* The system libm reference is less precise than the implementation on
     * some platforms. This sweep is therefore a regression/monotonicity
     * test; the shipped coefficients carry the tighter analytical bound. */
    int64_t previous = -1;
    long double max_error = 0.0L;
    for (uint64_t i = 0; i <= UINT64_C(65536); ++i) {
        uint64_t mantissa = (i == UINT64_C(65536))
            ? ((UINT64_C(1) << 57) - 1U)
            : (i << 41);
        set_mantissa(random, mantissa);
        approx_log2_sample(&sample, random, NULL);
        if (sample.fraction_q62 < previous)
            FAIL("log polynomial is not monotone at sample %llu: %lld < %lld",
                 (unsigned long long)i,
                 (long long)sample.fraction_q62,
                 (long long)previous);
        previous = sample.fraction_q62;

        long double got = ldexpl((long double)sample.fraction_q62, -62);
        long double want = log2l(1.0L + ldexpl((long double)mantissa, -57));
        long double error = fabsl(got - want);
        if (error > max_error)
            max_error = error;
    }
    if (max_error > ldexpl(1.0L, -50))
        FAIL("sampled log2 error %.4Le exceeds the libm regression limit",
             max_error);

    printf("    maximum sampled |log2 error| = %.4Le\n", max_error);
    PASS("four-segment degree-13 log2 polynomial is monotone and accurate");
    return 0;
}

static int reference_update_sign(int64_t inner_product,
                                 int64_t shift_norm_sq,
                                 const shuttle_log2_sample *sample) {
    long double t = fabsl((long double)inner_product);
    long double log2_u = ldexpl((long double)sample->fraction_q62, -62)
                       - (long double)sample->exponent;
    long double lambda = 2.0L * SHUTTLE_SIGMA * SHUTTLE_SIGMA
                       * logl(2.0L) * log2_u;
    int in_union = 0;

    for (unsigned int i = 0; i < (SHUTTLE_IRS_N + 1U) / 2U; ++i) {
        long double even = 2.0L * i;
        long double odd = even + 1.0L;
        long double lower = -(2.0L * odd * t
                              + odd * odd * shift_norm_sq);
        long double upper = -(2.0L * even * t
                              + even * even * shift_norm_sq);
        in_union |= (lower <= lambda && lambda < upper);
    }

    int orientation = inner_product < 0 ? -1 : 1;
    return orientation * (in_union ? -1 : 1);
}

static int test_transition_intervals(void) {
    static const int64_t inner_products[] = {
        INT64_C(-500000), INT64_C(-1000), INT64_C(-1), INT64_C(0),
        INT64_C(1), INT64_C(1000), INT64_C(500000)
    };
    static const int64_t shift_norms[] = {
        INT64_C(1), INT64_C(64), INT64_C(1024), INT64_C(8192)
    };
    static const uint32_t exponents[] = {1U, 2U, 10U, 40U, 80U, 81U};
    static const int64_t fractions[] = {
        INT64_C(0), INT64_C(1) << 58, INT64_C(1) << 60,
        INT64_C(1) << 61, INT64_C(3) << 60
    };

#if SHUTTLE_NIST_MODE == 120
    const unsigned int expected_n = 75U;
#elif SHUTTLE_NIST_MODE == 180
    const unsigned int expected_n = 89U;
#else
    const unsigned int expected_n = 55U;
#endif
    if (SHUTTLE_IRS_N != expected_n)
        FAIL("IRS truncation index is %u, expected %u",
             SHUTTLE_IRS_N, expected_n);

    unsigned int checked = 0;
    for (unsigned int i = 0;
         i < sizeof inner_products / sizeof inner_products[0]; ++i) {
        for (unsigned int j = 0;
             j < sizeof shift_norms / sizeof shift_norms[0]; ++j) {
            for (unsigned int k = 0;
                 k < sizeof exponents / sizeof exponents[0]; ++k) {
                for (unsigned int f = 0;
                     f < sizeof fractions / sizeof fractions[0]; ++f) {
                    shuttle_log2_sample sample = {
                        fractions[f], exponents[k]
                    };
                    int got = irs_update_sign(inner_products[i],
                                              shift_norms[j], &sample);
                    int want = reference_update_sign(inner_products[i],
                                                     shift_norms[j], &sample);
                    if (got != want)
                        FAIL("transition mismatch: inner=%lld norm=%lld "
                             "exp=%u frac=%lld: got %d, expected %d",
                             (long long)inner_products[i],
                             (long long)shift_norms[j], exponents[k],
                             (long long)fractions[f], got, want);
                    ++checked;
                }
            }
        }
    }

    printf("    checked %u wide-range interval decisions with N=%u\n",
           checked, SHUTTLE_IRS_N);
    PASS("signed 128-bit interval decisions match an independent reference");
    return 0;
}

int main(void) {
    int ret = 0;

    printf("=== test_irs (MODE=%d, r=%d, N=%d) ===\n",
           SHUTTLE_NIST_MODE, SHUTTLE_SIGMA, SHUTTLE_IRS_N);
    ret |= test_exponent_layout();
    ret |= test_log_polynomial();
    ret |= test_transition_intervals();

    if (ret == 0)
        printf("\n=== All IRS tests PASSED ===\n");
    else
        printf("\n=== Some IRS tests FAILED ===\n");
    return ret;
}
