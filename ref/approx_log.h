/* Integer-only segmented approximation of log2(U) for the IRS transition. */

#ifndef SHUTTLE_APPROX_LOG_H
#define SHUTTLE_APPROX_LOG_H

#include <stdint.h>
#include "config.h"
#include "fips202.h"

#define SHUTTLE_LOG_MANTISSA_BITS 57
#define SHUTTLE_LOG_EXPONENT_BLOCK_BYTES \
    (SHUTTLE_LOG_EXPONENT_BLOCK_BITS / 8)
#define SHUTTLE_LOG_RANDOM_BYTES \
    (SHUTTLE_LOG_EXPONENT_BLOCK_BYTES + 8)

typedef struct {
    /* log2(U) = fraction_q62 / 2^62 - exponent. */
    int64_t fraction_q62;
    uint32_t exponent;
} shuttle_log2_sample;

#define approx_log2_sample SHUTTLE_NAMESPACE(approx_log2_sample)
void approx_log2_sample(shuttle_log2_sample *out,
                        const uint8_t random[SHUTTLE_LOG_RANDOM_BYTES],
                        keccak_state *public_rng_state);

#endif /* SHUTTLE_APPROX_LOG_H */
