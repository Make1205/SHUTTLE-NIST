/* Rejection-free logarithmic-domain iterative sampler for SHUTTLE. */

#ifndef SHUTTLE_REJSAMPLE_H
#define SHUTTLE_REJSAMPLE_H

#include <stdint.h>
#include "params.h"
#include "poly.h"
#include "polyvec.h"
#include "fips202.h"
#include "approx_log.h"

/* Internal transition helper exposed for arithmetic regression tests.
 * Returns +1 or -1, the coefficient multiplying the supplied shift v. */
#define irs_update_sign SHUTTLE_NAMESPACE(irs_update_sign)
int irs_update_sign(int64_t inner_product,
                    int64_t shift_norm_sq,
                    const shuttle_log2_sample *log_sample);

/* Apply exactly tau rejection-free transitions. */
#define irs_sign SHUTTLE_NAMESPACE(irs_sign)
int irs_sign(polyvec *z,
             const poly *c,
             const polyvec *sk_stretched,
             keccak_state *rng_state);

#endif /* SHUTTLE_REJSAMPLE_H */
