/*
 * rejsample.h - Iterative Rejection Sampling for SHUTTLE.
 *
 * Implements the Simplified Log-Domain iterative rejection sampling.
 *
 * For each monomial x^{j_i} in the challenge c (tau terms):
 *   1. Compute v_i = x^{j_i} * sk (cyclic shift of secret key vector)
 *   2. Compute u = <z, v_i> / r^2 (inner product normalized by sigma^2)
 *   3. Sample r_rand uniform, compute ell = -ln(r_rand) via ApproxNegLn
 *   4. Accept/reject based on ell >= ln(M), determine sign via IntervalParity
 *   5. Update z = z +/- v_i
 *
 * All internal computations are constant-time (no secret-dependent branches).
 *
 * Returns 1 on success, 0 if any step rejects (caller must restart signing).
 */

#ifndef SHUTTLE_REJSAMPLE_H
#define SHUTTLE_REJSAMPLE_H

#include <stdint.h>
#include "params.h"
#include "poly.h"
#include "polyvec.h"
#include "fips202.h"

/*
 * irs_sign - Perform iterative rejection sampling over all tau challenge terms.
 *
 * Parameters:
 *   z             [in/out] Response vector, initially set to y (the commitment
 *                          randomness). On success, z is updated to z +/- v_i
 *                          for each accepted challenge monomial.
 *   c             [in]     Challenge polynomial with exactly tau nonzero
 *                          coefficients in {-1, +1}.
 *   sk_stretched  [in]     Stretched secret key vector [alpha_1, s, e] of
 *                          length VECLEN=6, each polynomial has N=256 int32
 *                          coefficients.
 *   gamma_q62     [in]     ||sk||^2 / (2 * sigma^2) in Q62 fixed-point.
 *   ln_M_q62      [in]     ln(M_alpha) in Q62 fixed-point.
 *                          For SHUTTLE with alpha >= 3, this is approximately 0.
 *   rng_state     [in/out] SHAKE-256 state for squeezing random bytes.
 *                          Must be initialized and finalized before calling.
 *
 * Returns:
 *   1 if all tau steps accepted (z is a valid response).
 *   0 if any step rejected (caller must restart with new y).
 */
#define irs_sign SHUTTLE_NAMESPACE(irs_sign)
int irs_sign(polyvec *z,
             const poly *c,
             const polyvec *sk_stretched,
             int64_t gamma_q62,
             int64_t ln_M_q62,
             keccak_state *rng_state);

#endif /* SHUTTLE_REJSAMPLE_H */
