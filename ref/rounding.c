#include <stdint.h>
#include "params.h"
#include "reduce.h"
#include "rounding.h"

/*
 * Rounding and hint helpers for SHUTTLE.
 *
 * SHUTTLE-NIST uses alpha_h = 128 uniformly, so 2*alpha_h = 256 across
 * all three modes. The natural range of HighBits is [0, W1_MAX] with
 * W1_MAX = (q-1)/256 = 126/60/30 for NIST-120/180/255 respectively. The
 * round-half-up edge case can produce W1_MAX+1; the highbits_mod_2q
 * path centers that back to 0 via a bit-mask (constant-time).
 *
 * All routines in this file operate without secret-dependent branches.
 */

/*************************************************
* Name:        highbits_mod_2q
*
* Description: Round-half-up HighBits for mod 2q commitment. Mirrors the
*              the spec's <x/alpha_h>.
*
*              Edge case: round(w/alpha_h) can reach HINT_MAX (at the top
*              of [0, 2q)); wrap it back to 0 to stay in [0, HINT_MAX).
*              Implemented with bit-mask for constant-time.
*
* Arguments:   - int32_t w: input, expected in [0, 2q).
*
* Returns bucket index in [0, HINT_MAX).
**************************************************/
int32_t highbits_mod_2q(int32_t w) {
  int32_t hb;
  int32_t edgecase;

  /* (w + alpha_h/2) >> log2(alpha_h) == (w + alpha_h/2) / alpha_h. */
  hb = (w + SHUTTLE_HALF_ALPHA_H) >> SHUTTLE_ALPHA_H_BITS;

  /* edgecase = -1 if hb >= HINT_MAX else 0. Using sign extension of the
   * difference: if (HINT_MAX - 1 - hb) < 0 then hb >= HINT_MAX. */
  edgecase = (SHUTTLE_HINT_MAX - 1 - hb) >> 31;
  hb -= SHUTTLE_HINT_MAX & edgecase;
  return hb;
}

/*************************************************
* Name:        lift_to_2q
*
* Description: Combine the "even part" (2 * u, with u in [0, q)) and the
*              "q-parity part" (q * bit) into a single [0, 2q) residue.
*              This is the per-coefficient operation behind compute_commitment
*              when moving the output of standard mod q NTT into the mod 2q
*              domain demanded by the derived public matrix hat_A.
*
* Arguments:   - int32_t u_mod_q: input in [0, q).
*              - int32_t parity:  0 or 1.
*
* Returns (2*u + q*parity) mod 2q in [0, 2q). Constant-time.
**************************************************/
int32_t lift_to_2q(int32_t u_mod_q, int32_t parity) {
  int32_t r;
  int32_t over;

  r = (u_mod_q << 1) + ((parity & 1) ? SHUTTLE_Q : 0);

  /* r is in [0, 2q) + [0, q) = [0, 3q). If r >= 2q, subtract 2q.
   * Constant-time using sign bit of (2q - 1 - r). */
  over = (SHUTTLE_DQ - 1 - r) >> 31;   /* -1 if r >= 2q else 0 */
  r -= SHUTTLE_DQ & over;
  return r;
}

/*************************************************
* Name:        make_hint_mod2q
*
* Description: Compute MakeHint difference per coefficient (spec Alg 9).
*              See rounding.h for the exact formula.
*
*              Returns the centered representative in (-HINT_MAX/2, HINT_MAX/2]:
*              the raw difference w_h - tilde_w_h lives in (-HINT_MAX, HINT_MAX),
*              and near bucket-wrap boundaries it would otherwise drift to
*              magnitude ~HINT_MAX. Centering unifies the two equivalent
*              representations into a single small-magnitude distribution
*              suitable for rANS encoding.
**************************************************/
int32_t make_hint_mod2q(int32_t w, int32_t z2_coef) {
  int32_t w_h = highbits_mod_2q(w);
  int32_t tilde_w = reduce_mod_2q(w - 2 * z2_coef);
  int32_t tilde_w_h = highbits_mod_2q(tilde_w);
  int32_t h = w_h - tilde_w_h;

  /* Center to (-HINT_MAX/2, HINT_MAX/2], constant-time. */
  int32_t half = SHUTTLE_HINT_MAX >> 1;
  int32_t upper_mask = (half - h) >> 31;      /* -1 if h > half */
  h -= SHUTTLE_HINT_MAX & upper_mask;
  int32_t lower_mask = (h + half) >> 31;      /* -1 if h < -half */
  h += SHUTTLE_HINT_MAX & lower_mask;
  return h;
}

/*************************************************
* Name:        use_hint_wh_mod2q
*
* Description: Recover signer's w_h from tilde_w and hint h (spec Alg 10
*              first half). Returns w_h in [0, HINT_MAX). Constant-time.
**************************************************/
int32_t use_hint_wh_mod2q(int32_t tilde_w, int32_t h) {
  int32_t tilde_w_h = highbits_mod_2q(tilde_w);
  int32_t w_h = tilde_w_h + h;
  int32_t over_mask, under_mask;

  /* Constant-time modular normalization to [0, HINT_MAX):
   *   if w_h >= HINT_MAX  ->  w_h -= HINT_MAX
   *   if w_h < 0          ->  w_h += HINT_MAX
   * We use sign-bit masking of (HINT_MAX - 1 - w_h) and (w_h) respectively. */
  over_mask  = (SHUTTLE_HINT_MAX - 1 - w_h) >> 31;   /* -1 if w_h >= HINT_MAX */
  w_h -= SHUTTLE_HINT_MAX & over_mask;

  under_mask = w_h >> 31;                             /* -1 if w_h < 0 */
  w_h += SHUTTLE_HINT_MAX & under_mask;
  return w_h;
}

/*************************************************
* Name:        recover_z2_coef_mod2q
*
* Description: Recover z_2 coefficient (spec Alg 10 second half).
*              See rounding.h for the exact formula and the +/-alpha_h/2
*              rounding-error bound.
**************************************************/
int32_t recover_z2_coef_mod2q(int32_t w_h, int32_t w_0, int32_t tilde_w) {
  int32_t w_app = SHUTTLE_ALPHA_H * w_h + w_0;
  int32_t diff  = w_app - tilde_w;
  int32_t over_mask, under_mask;

  /* diff currently in (-2q, 2q). Canonicalize to (-q, q]:
   *   if diff >  q  ->  diff -= 2q
   *   if diff <= -q ->  diff += 2q */
  over_mask  = (SHUTTLE_Q - diff) >> 31;          /* -1 if diff > q */
  diff -= SHUTTLE_DQ & over_mask;

  under_mask = (diff + SHUTTLE_Q) >> 31;          /* -1 if diff < -q */
  diff += SHUTTLE_DQ & under_mask;

  /* Integer division by 2 (arithmetic shift for two's complement). */
  return diff >> 1;
}
