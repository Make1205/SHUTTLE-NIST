/*
 * test_sig.c - round-trip tests for SHUTTLE-NIST pack_sig / unpack_sig
 *              under the single merged-stream rANS layout (design S).
 *
 * Tests:
 *   1. Minimal hand-crafted round-trip (all-zero z and h).
 *   2. 10 random rounds: z_1 is obtained by rounding response-Gaussian
 *      coefficients and h is drawn through the actual MakeHint map;
 *      reports the merged rANS length.
 *   3. Canonical-encoding rejection for nonzero padding and an extended
 *      rANS length prefix.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#include "../params.h"
#include "../poly.h"
#include "../polyvec.h"
#include "../packing.h"
#include "../rounding.h"
#include "../shuttle_rans.h"
#include "../rans_tables.h"

#define PASS(msg) do { printf("  [PASS] " msg "\n"); } while (0)
#define FAIL(msg, ...) do { printf("  [FAIL] " msg "\n", __VA_ARGS__); return 1; } while (0)

/* Response-Gaussian sampler on [-11*sigma,11*sigma] via cached inverse CDF. */
static int32_t sample_response_gaussian(void) {
  static int valid = 0;
  static double cdf[2 * SHUTTLE_BOUND + 1];
  if (!valid) {
    double Z = 0.0;
    for (int k = -SHUTTLE_BOUND; k <= SHUTTLE_BOUND; ++k)
      Z += exp(-(double)k * k /
               (2.0 * SHUTTLE_SIGMA * SHUTTLE_SIGMA));
    double running = 0.0;
    for (int k = -SHUTTLE_BOUND; k <= SHUTTLE_BOUND; ++k) {
      running += exp(-(double)k * k /
                     (2.0 * SHUTTLE_SIGMA * SHUTTLE_SIGMA)) / Z;
      cdf[k + SHUTTLE_BOUND] = running;
    }
    valid = 1;
  }
  double u = ((double)rand() + 1.0) / ((double)RAND_MAX + 2.0);
  for (int i = 0; i <= 2 * SHUTTLE_BOUND; ++i)
    if (u <= cdf[i]) return (int32_t)(i - SHUTTLE_BOUND);
  return SHUTTLE_BOUND;
}

/* Fill z_1[0..L] from the rounded response-Gaussian marginals. */
static void fill_z1(poly *z_1) {
  for (unsigned j = 0; j < SHUTTLE_N; ++j) {
    int32_t z = sample_response_gaussian();
    z_1[0].coeffs[j] =
        (z + SHUTTLE_ALPHA_1 / 2) >> SHUTTLE_ALPHA_1_BITS;
  }
  for (unsigned i = 1; i <= SHUTTLE_L; ++i)
    for (unsigned j = 0; j < SHUTTLE_N; ++j)
      z_1[i].coeffs[j] = sample_response_gaussian();
}

/* ------------------------------------------------------------ */
static int test_minimal(void) {
  printf("Test 1: minimal hand-crafted round-trip\n");

  uint8_t sig[SHUTTLE_BYTES];
  uint8_t c_tilde[SHUTTLE_CTILDEBYTES];
  poly    z_1[1 + SHUTTLE_L];
  polyveck h, h_rec;

  for (unsigned i = 0; i < SHUTTLE_CTILDEBYTES; ++i) c_tilde[i] = (uint8_t)(i * 17 + 3);
  for (unsigned i = 0; i < 1 + SHUTTLE_L; ++i) memset(&z_1[i], 0, sizeof(poly));
  memset(&h, 0, sizeof h);

  int rc = pack_sig(sig, c_tilde, z_1, &h);
  if (rc != 0) FAIL("pack_sig rc=%d", rc);

  uint8_t c_rec[SHUTTLE_CTILDEBYTES];
  poly    z_rec[1 + SHUTTLE_L];
  rc = unpack_sig(c_rec, z_rec, &h_rec, sig);
  if (rc != 0) FAIL("unpack_sig rc=%d", rc);

  if (memcmp(c_tilde, c_rec, SHUTTLE_CTILDEBYTES) != 0)
    FAIL("c_tilde mismatch%s", "");
  for (unsigned i = 0; i < 1 + SHUTTLE_L; ++i)
    for (unsigned j = 0; j < SHUTTLE_N; ++j)
      if (z_1[i].coeffs[j] != z_rec[i].coeffs[j])
        FAIL("z_1[%u][%u]: %d vs %d",
             i, j, z_1[i].coeffs[j], z_rec[i].coeffs[j]);
  for (unsigned i = 0; i < SHUTTLE_M; ++i)
    for (unsigned j = 0; j < SHUTTLE_N; ++j)
      if (h.vec[i].coeffs[j] != h_rec.vec[i].coeffs[j])
        FAIL("h[%u][%u]: %d vs %d",
             i, j, h.vec[i].coeffs[j], h_rec.vec[i].coeffs[j]);

  printf("    SHUTTLE_BYTES = %d\n", SHUTTLE_BYTES);
  PASS("minimal round-trip");
  return 0;
}

/* ------------------------------------------------------------ */
static int test_random_rounds(void) {
  printf("Test 2: 10 random rounds\n");

  size_t min_rans = (size_t)-1, max_rans = 0, sum_rans = 0;

  /* Offset of the merged rANS length-prefix in the signature. */
  const size_t off_rans_len = SHUTTLE_CTILDEBYTES;

  for (unsigned round = 0; round < 10; ++round) {
    uint8_t sig[SHUTTLE_BYTES];
    uint8_t c_tilde[SHUTTLE_CTILDEBYTES];
    poly    z_1[1 + SHUTTLE_L];
    polyveck h, h_rec;

    for (unsigned i = 0; i < SHUTTLE_CTILDEBYTES; ++i) c_tilde[i] = (uint8_t)rand();
    fill_z1(z_1);
    /* Draw hints through the implementation's MakeHint map. */
    for (unsigned i = 0; i < SHUTTLE_M; ++i)
      for (unsigned j = 0; j < SHUTTLE_N; ++j) {
        int32_t w = rand() % SHUTTLE_DQ;
        int32_t z2 = sample_response_gaussian();
        h.vec[i].coeffs[j] = make_hint_mod2q(w, z2);
      }

    int rc = pack_sig(sig, c_tilde, z_1, &h);
    if (rc != 0) FAIL("round %u: pack rc=%d", round, rc);

    uint8_t c_rec[SHUTTLE_CTILDEBYTES];
    poly    z_rec[1 + SHUTTLE_L];
    rc = unpack_sig(c_rec, z_rec, &h_rec, sig);
    if (rc != 0) FAIL("round %u: unpack rc=%d", round, rc);

    if (memcmp(c_tilde, c_rec, SHUTTLE_CTILDEBYTES) != 0) FAIL("round %u: c_tilde", round);
    for (unsigned i = 0; i < 1 + SHUTTLE_L; ++i)
      for (unsigned j = 0; j < SHUTTLE_N; ++j)
        if (z_1[i].coeffs[j] != z_rec[i].coeffs[j])
          FAIL("round %u: z_1[%u][%u] %d vs %d",
               round, i, j, z_1[i].coeffs[j], z_rec[i].coeffs[j]);
    for (unsigned i = 0; i < SHUTTLE_M; ++i)
      for (unsigned j = 0; j < SHUTTLE_N; ++j)
        if (h.vec[i].coeffs[j] != h_rec.vec[i].coeffs[j])
          FAIL("round %u: h[%u][%u]", round, i, j);

    size_t rans_len = (size_t)sig[off_rans_len]
                    | ((size_t)sig[off_rans_len + 1] << 8);

    if (rans_len < min_rans) min_rans = rans_len;
    if (rans_len > max_rans) max_rans = rans_len;
    sum_rans += rans_len;
  }

  printf("    merged rANS len: min=%zu max=%zu avg=%zu (reserved %d)\n",
         min_rans, max_rans, sum_rans / 10, SHUTTLE_RANS_RESERVED_BYTES);
  PASS("10 random rounds OK");
  return 0;
}

/* ------------------------------------------------------------ */
static int test_canonical_encoding(void) {
  printf("Test 3: canonical rANS encoding\n");

  const size_t off_rans_len = SHUTTLE_CTILDEBYTES;
  const size_t off_rans_data = off_rans_len + 2;
  uint8_t sig[SHUTTLE_BYTES];
  uint8_t bad[SHUTTLE_BYTES];
  uint8_t c_tilde[SHUTTLE_CTILDEBYTES] = {0};
  uint8_t c_rec[SHUTTLE_CTILDEBYTES];
  poly z_1[1 + SHUTTLE_L];
  poly z_rec[1 + SHUTTLE_L];
  polyveck h, h_rec;

  memset(z_1, 0, sizeof z_1);
  memset(&h, 0, sizeof h);
  if(pack_sig(sig, c_tilde, z_1, &h) != 0)
    FAIL("canonical setup pack failed%s", "");

  size_t rans_len = (size_t)sig[off_rans_len]
                  | ((size_t)sig[off_rans_len + 1] << 8);
  if(rans_len >= SHUTTLE_RANS_RESERVED_BYTES)
    FAIL("canonical setup has no padding: %zu%s", rans_len, "");

  memcpy(bad, sig, sizeof bad);
  bad[off_rans_data + rans_len] = 1;
  if(unpack_sig(c_rec, z_rec, &h_rec, bad) == 0)
    FAIL("nonzero rANS padding accepted%s", "");

  memcpy(bad, sig, sizeof bad);
  ++rans_len;
  bad[off_rans_len] = (uint8_t)rans_len;
  bad[off_rans_len + 1] = (uint8_t)(rans_len >> 8);
  if(unpack_sig(c_rec, z_rec, &h_rec, bad) == 0)
    FAIL("extended rANS input accepted%s", "");

  PASS("non-canonical encodings rejected");
  return 0;
}

/* ------------------------------------------------------------ */
int main(void) {
  int ret = 0;
  srand(0xC0DECAFEu ^ (unsigned)SHUTTLE_NIST_MODE);

  printf("=== test_sig (MODE=%d, SHUTTLE_BYTES=%d, RANS_RESERVED=%d) ===\n",
         SHUTTLE_NIST_MODE, SHUTTLE_BYTES, SHUTTLE_RANS_RESERVED_BYTES);

  ret |= test_minimal();
  ret |= test_random_rounds();
  ret |= test_canonical_encoding();

  if (ret == 0) printf("\n=== All sig tests PASSED ===\n");
  else          printf("\n=== Some sig tests FAILED ===\n");
  return ret;
}
