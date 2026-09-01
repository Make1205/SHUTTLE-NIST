/*
 * test_sampler.c - Test the SHUTTLE discrete Gaussian sampler.
 *
 * Tests:
 *   1) Basic statistical properties (mean ~0, variance ~sigma^2)
 *   2) Hard-support check: |r| <= SHUTTLE_GAUSS_MAX
 *   3) KAT reproducibility: same seed+nonce -> same output
 *   4) 4x KAT: sample_gauss_N_4x lane 0 matches sample_gauss_N
 *   5) Exhaustive ApproxExp accuracy over the configured (x,y) support
 *   6) Exact 3x31-bit RCDT comparisons at every threshold boundary
 */

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../sampler.h"
#include "../approx_exp.h"
#include "../approx_exp_constants.h"
#include "../rcdt_tables.h"
#include "../randombytes.h"

#define N 256

static void store_le32_test(uint8_t *out, uint32_t value) {
    out[0] = (uint8_t)value;
    out[1] = (uint8_t)(value >> 8);
    out[2] = (uint8_t)(value >> 16);
    out[3] = (uint8_t)(value >> 24);
}

static void set_rcdt_lane0(uint8_t rand[SIGMA2_RAND_BYTES],
                           uint32_t lo, uint32_t mid, uint32_t hi) {
    memset(rand, 0, SIGMA2_RAND_BYTES);
    store_le32_test(rand, lo);
    store_le32_test(rand + 32, mid);
    store_le32_test(rand + 64, hi);
}

static int test_rcdt_boundaries(void) {
    const uint32_t mask = UINT32_C(0x7FFFFFFF);
    uint8_t rand[SIGMA2_RAND_BYTES];
    int16_t samples[GAUSS_BATCH];

    for (int i = 0; i < RCDT_ENTRIES; ++i) {
        uint32_t lo = RCDT_3x31[i][0];
        uint32_t mid = RCDT_3x31[i][1];
        uint32_t hi = RCDT_3x31[i][2];

        set_rcdt_lane0(rand, lo, mid, hi);
        sampler_sigma2(samples, rand);
        if (samples[0] != i)
            return 0;

        if (lo != 0) {
            --lo;
        } else {
            lo = mask;
            if (mid != 0) {
                --mid;
            } else {
                mid = mask;
                --hi;
            }
        }
        set_rcdt_lane0(rand, lo, mid, hi);
        sampler_sigma2(samples, rand);
        if (samples[0] != i + 1)
            return 0;
    }
    return 1;
}

static int test_approx_exp_range(void) {
    long double max_a = 0.0L;
    long double max_rel = 0.0L;

    for (uint64_t x = 0; x <= RCDT_ENTRIES; ++x) {
        for (uint64_t y = 0; y < (UINT64_C(1) << SHUTTLE_Y_BITS); ++y) {
            uint64_t num = y * (y + 2 * SHUTTLE_GAUSS_K * x);
            wide_uint128 m_lo = (wide_uint128)num * APPROX_EXP_I_LO;
            uint64_t m_hi = num * APPROX_EXP_I_HI;
            wide_uint128 product = m_lo + ((wide_uint128)m_hi << 64);
            uint64_t a_q60 = (uint64_t)(product >> 20);
            long double a = (long double)num /
                            (2.0L * SHUTTLE_SIGMA * SHUTTLE_SIGMA);
            long double exact = expl(-a) * ldexpl(1.0L, 63);
            long double rel = fabsl((long double)approx_exp(a_q60) / exact
                                    - 1.0L);
            if (a > max_a) max_a = a;
            if (rel > max_rel) max_rel = rel;
        }
    }

    printf("ApproxExp: a_max=%.6Lf, max rel error=2^(-%.2Lf): ",
           max_a, -log2l(max_rel));
    if (max_a >= 11.0L || max_rel >= ldexpl(1.0L, -48)) {
        printf("FAIL\n");
        return 0;
    }
    printf("PASS\n");
    return 1;
}

int main(void) {
    uint8_t seed[SHUTTLE_SEEDBYTES];
    int16_t r[N];
    int ok = 1;

    /* Fixed seed for reproducibility */
    memset(seed, 0x42, SHUTTLE_SEEDBYTES);

    /* Test 1: Basic sampling and statistics */
    printf("=== SHUTTLE Gaussian Sampler Test (sigma=%d, N=%d) ===\n\n",
           SHUTTLE_SIGMA, N);

    ok &= test_approx_exp_range();

    printf("RCDT boundaries: ");
    if (!test_rcdt_boundaries()) {
        printf("FAIL\n");
        ok = 0;
    } else {
        printf("PASS\n");
    }

    sample_gauss_N(r, seed, 0, N);

    double mean = 0, var = 0;
    int max_abs = 0;
    for (int i = 0; i < N; i++) {
        mean += r[i];
        var += (double)r[i] * r[i];
        int a = abs(r[i]);
        if (a > max_abs) max_abs = a;
    }
    mean /= N;
    var = var / N - mean * mean;

    printf("Statistics:\n");
    printf("  Mean:     %.2f (expected ~0)\n", mean);
    printf("  Variance: %.2f (expected ~%.2f)\n", var, (double)SHUTTLE_SIGMA * SHUTTLE_SIGMA);
    printf("  Max|r|:   %d (hard support: %d)\n",
           max_abs, SHUTTLE_GAUSS_MAX);

    /* Test 2: Bound check */
    printf("\nBound check: ");
    for (int i = 0; i < N; i++) {
        if (r[i] > SHUTTLE_GAUSS_MAX || r[i] < -SHUTTLE_GAUSS_MAX) {
            printf("FAIL (r[%d] = %d out of bounds!)\n", i, r[i]);
            ok = 0;
        }
    }
    if (ok) printf("PASS\n");

    /* Test 3: KAT reproducibility */
    int16_t r2[N];
    sample_gauss_N(r2, seed, 0, N);
    printf("KAT test:   ");
    if (memcmp(r, r2, N * sizeof(int16_t)) != 0) {
        printf("FAIL (different output for same seed+nonce)\n");
        ok = 0;
    } else {
        printf("PASS\n");
    }

    /* Test 4: 4x KAT - lane 0 must match single call */
    int16_t r4x[4][N];
    sample_gauss_N_4x(r4x[0], r4x[1], r4x[2], r4x[3],
                       seed, 0, 1, 2, 3, N, N, N, N);

    printf("4x KAT L0:  ");
    if (memcmp(r4x[0], r, N * sizeof(int16_t)) != 0) {
        printf("FAIL (4x lane 0 != single)\n");
        ok = 0;
    } else {
        printf("PASS\n");
    }

    /* Check other lanes match individual calls */
    for (int lane = 1; lane < 4; lane++) {
        int16_t rs[N];
        sample_gauss_N(rs, seed, (uint64_t)lane, N);
        printf("4x KAT L%d:  ", lane);
        if (memcmp(r4x[lane], rs, N * sizeof(int16_t)) != 0) {
            printf("FAIL\n");
            ok = 0;
        } else {
            printf("PASS\n");
        }
    }

    /* Test 5: Larger statistical test with random seed */
    printf("\nLarge statistical test (10 x N=%d):\n", N);
    randombytes(seed, SHUTTLE_SEEDBYTES);
    double total_mean = 0, total_var = 0;
    int total_max = 0;
    int total_samples = 10 * N;

    for (int round = 0; round < 10; round++) {
        sample_gauss_N(r, seed, (uint64_t)round, N);
        for (int i = 0; i < N; i++) {
            total_mean += r[i];
            total_var += (double)r[i] * r[i];
            int a = abs(r[i]);
            if (a > total_max) total_max = a;
            if (r[i] > SHUTTLE_GAUSS_MAX || r[i] < -SHUTTLE_GAUSS_MAX) {
                printf("  ERROR: round %d, r[%d] = %d out of bounds!\n",
                       round, i, r[i]);
                ok = 0;
            }
        }
    }
    total_mean /= total_samples;
    total_var = total_var / total_samples - total_mean * total_mean;

    printf("  Mean:     %.2f\n", total_mean);
    printf("  Variance: %.2f (expected ~%.2f)\n", total_var,
           (double)SHUTTLE_SIGMA * SHUTTLE_SIGMA);
    printf("  Max|r|:   %d\n", total_max);

    printf("\n=== Overall: %s ===\n", ok ? "ALL PASSED" : "SOME FAILED");
    return ok ? 0 : 1;
}
