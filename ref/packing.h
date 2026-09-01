/*
 * packing.h - Serialization of keys and signatures for SHUTTLE-NIST
 * (design S, single merged rANS stream).
 *
 * Single rANS stream and frequency table: z-hi (HighBits of
 * z^(0..lenS)) and hint are encoded by one shuttle_rans_encode call.
 * The shared table is generated from their coefficient-count-weighted
 * mixture; the individual PMFs are not assumed to be identical.
 *
 * LowBits are bit-packed uniformly:
 *   - z^(0)        : ALPHA_0P_BITS/coef via polyz0_lo_pack
 *   - z^(1..lenS)  : ALPHA_R_BITS = 6/coef via polyz1_lo_pack
 *
 * Byte layout (see packing.c::OFF_* for the offsets):
 *   seedC || uint16 rans_len
 *   || rANS(z-hi || hint) + pad to RANS_RESERVED
 *   || polyz0_lo_pack(lo(z^(0)))
 *   || L * polyz1_lo_pack(lo(z^(1..lenS)))
 *
 * SHUTTLE_BYTES = CTILDEBYTES + RANS_BLOCK_BYTES
 *               + POLYZ0_LO_PACKEDBYTES + L * POLYZ1_LO_PACKEDBYTES.
 *
 * Signer-side reject conditions (caller restarts IRS):
 *   - An rANS symbol is outside the configured vocabulary (rc = -1).
 *   - The rANS stream exceeds its reservation (rc = -2).
 *
 * Verifier-side rejects:
 *   - Length field exceeds reservation.
 *   - rANS decode underflows, leaves trailing input, or fails the final-state
 *     check (x != L).
 *   - The fixed rANS slot is not the canonical encoding (including nonzero
 *     padding after the encoded stream).
 */

#ifndef SHUTTLE_PACKING_H
#define SHUTTLE_PACKING_H

#include <stddef.h>
#include <stdint.h>

#include "params.h"
#include "poly.h"
#include "polyvec.h"

#define pack_pk SHUTTLE_NAMESPACE(pack_pk)
void pack_pk(uint8_t pk[SHUTTLE_PUBLICKEYBYTES],
             const uint8_t rho[SHUTTLE_SEEDBYTES],
             const polyveck *b);

#define unpack_pk SHUTTLE_NAMESPACE(unpack_pk)
void unpack_pk(uint8_t rho[SHUTTLE_SEEDBYTES],
               polyveck *b,
               const uint8_t pk[SHUTTLE_PUBLICKEYBYTES]);

#define pack_sk SHUTTLE_NAMESPACE(pack_sk)
void pack_sk(uint8_t sk[SHUTTLE_SECRETKEYBYTES],
             const uint8_t rho[SHUTTLE_SEEDBYTES],
             const uint8_t tr[SHUTTLE_TRBYTES],
             const uint8_t key[SHUTTLE_SEEDBYTES],
             const polyvecl *s,
             const polyveck *e);

#define unpack_sk SHUTTLE_NAMESPACE(unpack_sk)
void unpack_sk(uint8_t rho[SHUTTLE_SEEDBYTES],
               uint8_t tr[SHUTTLE_TRBYTES],
               uint8_t key[SHUTTLE_SEEDBYTES],
               polyvecl *s,
               polyveck *e,
               const uint8_t sk[SHUTTLE_SECRETKEYBYTES]);

/* pack_sig
 *
 * Inputs:
 *   sig        : output buffer, SHUTTLE_BYTES bytes.
 *   c_tilde    : CTILDEBYTES challenge seed.
 *   z_1        : (1 + L) polys. z_1[0] is the CompressY'd Z_0, z_1[1..L] are
 *                the full-range z[1..L]. Passed as a pointer to the first
 *                polynomial of the array (caller supplies an array of length
 *                1 + L; we access z_1[0] ... z_1[L]).
 *   h          : polyveck of M hint polynomials (MakeHint output).
 *
 * Returns 0 on success; a negative code on rANS overflow / OOV. */
#define pack_sig SHUTTLE_NAMESPACE(pack_sig)
int pack_sig(uint8_t sig[SHUTTLE_BYTES],
             const uint8_t c_tilde[SHUTTLE_CTILDEBYTES],
             const poly *z_1,               /* array of 1 + L polys */
             const polyveck *h);

/* unpack_sig
 *
 * Outputs:
 *   c_tilde    : seedC.
 *   z_1        : array of 1 + L polys.
 *   h          : recovered hint.
 * Inputs:
 *   sig        : byte stream of length SHUTTLE_BYTES.
 *
 * Returns 0 on success; nonzero on malformed stream. */
#define unpack_sig SHUTTLE_NAMESPACE(unpack_sig)
int unpack_sig(uint8_t c_tilde[SHUTTLE_CTILDEBYTES],
               poly *z_1,                   /* array of 1 + L polys */
               polyveck *h,
               const uint8_t sig[SHUTTLE_BYTES]);

#endif /* SHUTTLE_PACKING_H */
