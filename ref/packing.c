/*
 * packing.c - Serialization of keys and signatures for SHUTTLE-NIST.
 *
 * Public key format:
 *   rho (SEEDBYTES) || polypk_pack(b[0..M-1])
 *
 * Secret key format:
 *   rho || tr || key || polyeta_pack(s[0..L-1]) || polyeta_pack(e[0..M-1])
 *
 * Signature format (single merged rANS stream):
 *
 *   seedC || uint16 rans_len
 *     || rANS(zhi || hint) + pad to RANS_RESERVED
 *     || polyz0_lo_pack(lo(z^(0)))
 *     || L * polyz1_lo_pack(lo(z^(1..lenS)))
 *
 * The rANS block carries n*(lenS+1) + n*lenE coefficients in order
 *   [ HighBits(z^(0)), HighBits(z^(1)), ..., HighBits(z^(lenS)),
 *     hint(0), hint(1), ..., hint(lenE-1) ].
 * All share scale r/alpha_r so a single frequency table fits them all
 * (alpha_h = 2*alpha_r = 128 is what makes sigma_hint = sigma_z-hi).
 *
 * The rANS block is a fixed-size slot (length prefix + encoded stream +
 * zero padding). The fixed slot is what lets the signature have a single
 * on-the-wire length. Reservation budget SHUTTLE_RANS_RESERVED_BYTES is
 * sized via the codelength model for overflow probability p_rans^* = 2^-30 (see
 * tools/SigSize.py).
 *
 * Encoder-side OOV/overflow returns a negative code from pack_sig, and the
 * signer restarts with fresh masking randomness.
 */

#include <string.h>

#include "packing.h"
#include "params.h"
#include "poly.h"
#include "polyvec.h"
#include "shuttle_rans.h"

/* ============================================================
 * Public / secret key packing (unchanged).
 * ============================================================ */

void pack_pk(uint8_t pk[SHUTTLE_PUBLICKEYBYTES],
             const uint8_t rho[SHUTTLE_SEEDBYTES],
             const polyveck *b)
{
  unsigned int i;

  memcpy(pk, rho, SHUTTLE_SEEDBYTES);
  pk += SHUTTLE_SEEDBYTES;

  for(i = 0; i < SHUTTLE_M; ++i) {
    polypk_pack(pk, &b->vec[i]);
    pk += SHUTTLE_POLYPK_PACKEDBYTES;
  }
}

void unpack_pk(uint8_t rho[SHUTTLE_SEEDBYTES],
               polyveck *b,
               const uint8_t pk[SHUTTLE_PUBLICKEYBYTES])
{
  unsigned int i;

  memcpy(rho, pk, SHUTTLE_SEEDBYTES);
  pk += SHUTTLE_SEEDBYTES;

  for(i = 0; i < SHUTTLE_M; ++i) {
    polypk_unpack(&b->vec[i], pk);
    pk += SHUTTLE_POLYPK_PACKEDBYTES;
  }
}

void pack_sk(uint8_t sk[SHUTTLE_SECRETKEYBYTES],
             const uint8_t rho[SHUTTLE_SEEDBYTES],
             const uint8_t tr[SHUTTLE_TRBYTES],
             const uint8_t key[SHUTTLE_SEEDBYTES],
             const polyvecl *s,
             const polyveck *e)
{
  unsigned int i;

  memcpy(sk, rho, SHUTTLE_SEEDBYTES);
  sk += SHUTTLE_SEEDBYTES;

  memcpy(sk, tr, SHUTTLE_TRBYTES);
  sk += SHUTTLE_TRBYTES;

  memcpy(sk, key, SHUTTLE_SEEDBYTES);
  sk += SHUTTLE_SEEDBYTES;

  for(i = 0; i < SHUTTLE_L; ++i) {
    polyeta_pack(sk, &s->vec[i]);
    sk += SHUTTLE_POLYETA_PACKEDBYTES;
  }

  for(i = 0; i < SHUTTLE_M; ++i) {
    polyeta_pack(sk, &e->vec[i]);
    sk += SHUTTLE_POLYETA_PACKEDBYTES;
  }
}

void unpack_sk(uint8_t rho[SHUTTLE_SEEDBYTES],
               uint8_t tr[SHUTTLE_TRBYTES],
               uint8_t key[SHUTTLE_SEEDBYTES],
               polyvecl *s,
               polyveck *e,
               const uint8_t sk[SHUTTLE_SECRETKEYBYTES])
{
  unsigned int i;

  memcpy(rho, sk, SHUTTLE_SEEDBYTES);
  sk += SHUTTLE_SEEDBYTES;

  memcpy(tr, sk, SHUTTLE_TRBYTES);
  sk += SHUTTLE_TRBYTES;

  memcpy(key, sk, SHUTTLE_SEEDBYTES);
  sk += SHUTTLE_SEEDBYTES;

  for(i = 0; i < SHUTTLE_L; ++i) {
    polyeta_unpack(&s->vec[i], sk);
    sk += SHUTTLE_POLYETA_PACKEDBYTES;
  }

  for(i = 0; i < SHUTTLE_M; ++i) {
    polyeta_unpack(&e->vec[i], sk);
    sk += SHUTTLE_POLYETA_PACKEDBYTES;
  }
}

/* ============================================================
 * Signature packing (design S: single merged rANS stream + z bit-packs).
 * ============================================================ */

#define OFF_C_TILDE     0
#define OFF_RANS_LEN    (OFF_C_TILDE + SHUTTLE_CTILDEBYTES)
#define OFF_RANS_DATA   (OFF_RANS_LEN + 2)
#define OFF_Z0_LO       (OFF_RANS_DATA + SHUTTLE_RANS_RESERVED_BYTES)
#define OFF_Z1_LO       (OFF_Z0_LO + SHUTTLE_POLYZ0_LO_PACKEDBYTES)

#define MERGED_SYM_COUNT  ((SHUTTLE_L + 1 + SHUTTLE_M) * SHUTTLE_N)

int pack_sig(uint8_t sig[SHUTTLE_BYTES],
             const uint8_t c_tilde[SHUTTLE_CTILDEBYTES],
             const poly *z_1,
             const polyveck *h)
{
  unsigned int i, j;
  int rc;

  /* 1. seedC */
  memcpy(&sig[OFF_C_TILDE], c_tilde, SHUTTLE_CTILDEBYTES);

  /* 2. Split z^(0) and z^(1..lenS) into hi (rANS) + lo (bit-pack), and
   *    flatten the hint coefficients. The merged rANS input layout is:
   *      merged[0 .. n-1]                    = HighBits_{alpha_0'}(z^(0))
   *      merged[k*n .. (k+1)*n - 1]          = HighBits_{alpha_r}(z^(k)),
   *                                            for k in {1..lenS}
   *      merged[(lenS+1)*n + j*n .. ...]     = hint[j], j in {0..lenE-1}. */
  int32_t merged[MERGED_SYM_COUNT];
  int32_t lo_scratch[SHUTTLE_N];

  polyz0_split(&merged[0], lo_scratch, &z_1[0]);
  polyz0_lo_pack(&sig[OFF_Z0_LO], lo_scratch);

  for(i = 0; i < SHUTTLE_L; ++i) {
    polyz1_split(&merged[(i + 1) * SHUTTLE_N], lo_scratch, &z_1[1 + i]);
    polyz1_lo_pack(&sig[OFF_Z1_LO + i * SHUTTLE_POLYZ1_LO_PACKEDBYTES],
                   lo_scratch);
  }

  /* hint at the tail of the merged stream. */
  for(i = 0; i < SHUTTLE_M; ++i)
    for(j = 0; j < SHUTTLE_N; ++j)
      merged[(SHUTTLE_L + 1 + i) * SHUTTLE_N + j] = h->vec[i].coeffs[j];

  /* 3. rANS-encode the merged stream. */
  size_t rans_len = 0;
  rc = shuttle_rans_encode(&sig[OFF_RANS_DATA], &rans_len,
                           SHUTTLE_RANS_RESERVED_BYTES,
                           merged, MERGED_SYM_COUNT);
  if(rc != 0)
    return rc;

  sig[OFF_RANS_LEN + 0] = (uint8_t)(rans_len & 0xFF);
  sig[OFF_RANS_LEN + 1] = (uint8_t)((rans_len >> 8) & 0xFF);
  if(rans_len < SHUTTLE_RANS_RESERVED_BYTES)
    memset(&sig[OFF_RANS_DATA + rans_len], 0,
           SHUTTLE_RANS_RESERVED_BYTES - rans_len);

  return 0;
}

int unpack_sig(uint8_t c_tilde[SHUTTLE_CTILDEBYTES],
               poly *z_1,
               polyveck *h,
               const uint8_t sig[SHUTTLE_BYTES])
{
  unsigned int i, j;
  int rc;

  /* 1. seedC */
  memcpy(c_tilde, &sig[OFF_C_TILDE], SHUTTLE_CTILDEBYTES);

  /* 2. Merged rANS decode. */
  size_t rans_len = (size_t)sig[OFF_RANS_LEN + 0]
                  | ((size_t)sig[OFF_RANS_LEN + 1] << 8);
  if(rans_len == 0 || rans_len > SHUTTLE_RANS_RESERVED_BYTES)
    return -1;

  /* The fixed slot has a unique representation: every byte after the
   * declared rANS stream must be zero. */
  for(i = (unsigned int)rans_len; i < SHUTTLE_RANS_RESERVED_BYTES; ++i)
    if(sig[OFF_RANS_DATA + i] != 0)
      return -1;

  int32_t merged[MERGED_SYM_COUNT];
  rc = shuttle_rans_decode(merged, MERGED_SYM_COUNT,
                           &sig[OFF_RANS_DATA], rans_len);
  if(rc != 0)
    return rc;

  /* Re-encoding makes the rANS representation canonical, independently of
   * decoder implementation details. */
  {
    uint8_t canonical[SHUTTLE_RANS_RESERVED_BYTES];
    size_t canonical_len = 0;
    rc = shuttle_rans_encode(canonical, &canonical_len,
                             SHUTTLE_RANS_RESERVED_BYTES,
                             merged, MERGED_SYM_COUNT);
    if(rc != 0 || canonical_len != rans_len ||
       memcmp(canonical, &sig[OFF_RANS_DATA], rans_len) != 0)
      return -1;
  }

  /* 3. z^(0): combine hi (from merged stream) + lo (bit-packed). */
  int32_t lo_scratch[SHUTTLE_N];
  polyz0_lo_unpack(lo_scratch, &sig[OFF_Z0_LO]);
  polyz0_combine(&z_1[0], &merged[0], lo_scratch);

  /* 4. z^(1..lenS): combine hi + lo. */
  for(i = 0; i < SHUTTLE_L; ++i) {
    polyz1_lo_unpack(lo_scratch,
                     &sig[OFF_Z1_LO + i * SHUTTLE_POLYZ1_LO_PACKEDBYTES]);
    polyz1_combine(&z_1[1 + i],
                   &merged[(i + 1) * SHUTTLE_N],
                   lo_scratch);
  }

  /* 5. hint at the tail of merged. */
  for(i = 0; i < SHUTTLE_M; ++i)
    for(j = 0; j < SHUTTLE_N; ++j)
      h->vec[i].coeffs[j] = merged[(SHUTTLE_L + 1 + i) * SHUTTLE_N + j];

  return 0;
}
