/*
 * shuttle_rans.h - rANS encoder/decoder for SHUTTLE-NIST (single merged stream).
 *
 * ============================================================
 * Design rationale
 * ============================================================
 *
 * SHUTTLE-NIST signatures carry a SINGLE rANS-encoded block that merges:
 *   - HighBits_{alpha_0'}(z^(0)) followed by HighBits_{alpha_r}(z^(1..lenS))
 *     (the z-hi sub-stream), and
 *   - hint h = MakeHint(...)
 *
 * The parameter set takes alpha_h = 2*alpha_r = 128 across all three
 * modes, so sigma_hint = 2r/alpha_h = r/alpha_r = sigma_z-hi. Hint and
 * z-hi share the same discrete-Gaussian PMF, and a single rANS context
 * with one frequency table encodes them end-to-end. The merge saves
 * ~16-21 B/sig of fixed overhead (1 flush + 1 length prefix) vs. a
 * two-stream design.
 *
 * Theoretical frequency table => OOV-free
 * ---------------------------------------
 * Frequency tables come from the theoretical discrete Gaussian PMF
 * (tools/gen_rans_tables.py), not from any empirical histogram. The
 * vocabulary radius is the max of the z-hi and hint l_infty tight bounds
 * (which coincide under alpha_h = 2*alpha_r), so every coefficient maps
 * to a real frequency bucket -- OOV failure is mathematically impossible.
 *
 * Failure modes
 * -------------
 * Encoder return values:
 *   0  on success
 *  -2  if the compressed stream would exceed max_bytes (overflow)
 *
 * Decoder return values:
 *   0  on success and x == L_ren at end-of-stream (state verified)
 *  -1  on underflow OR final-state mismatch (malformed signature)
 *
 * ============================================================
 * Engine conventions
 * ============================================================
 *
 * Byte-wise renormalization, 32-bit state, L = 2^23. Probability
 * precision t = 10 (table sums to 2^10 = 1024). The Python reference
 * lives in SHUTTLE-NIST/tools/rans.py and is bit-exact.
 *
 * Final-state flush is written little-endian.
 */

#ifndef SHUTTLE_RANS_H
#define SHUTTLE_RANS_H

#include <stddef.h>
#include <stdint.h>

#include "params.h"

/* Renorm threshold (shared across modes). */
#define SHUTTLE_RANS_L          (1u << 23)
#define SHUTTLE_RANS_BYTE_BITS  8
#define SHUTTLE_RANS_PROB_BITS  10
#define SHUTTLE_RANS_PROB_TOTAL (1u << SHUTTLE_RANS_PROB_BITS)

/* Encode/decode contracts for the single merged stream.
 *
 *   0  : success
 *  -2  : output would exceed max_bytes (encoder; signer restarts IRS)
 *  -1  : underflow or final-state mismatch (decoder)
 */
#define shuttle_rans_encode SHUTTLE_NAMESPACE(shuttle_rans_encode)
int shuttle_rans_encode(uint8_t *out, size_t *out_len, size_t max_bytes,
                        const int32_t *syms, size_t n);

#define shuttle_rans_decode SHUTTLE_NAMESPACE(shuttle_rans_decode)
int shuttle_rans_decode(int32_t *syms, size_t n,
                        const uint8_t *in, size_t in_len);

#endif /* SHUTTLE_RANS_H */
