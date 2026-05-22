/*
 * test_sig.c - round-trip tests for SHUTTLE-NIST pack_sig / unpack_sig
 *              under the single merged-stream rANS layout (design S).
 *
 * Tests:
 *   1. Minimal hand-crafted round-trip (all-zero z and h).
 *   2. 10 random rounds: z_1 drawn so HighBits live in the merged
 *      vocabulary, h drawn from the same PMF (since sigma_hint =
 *      sigma_z-hi); reports merged rANS length.
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
#include "../shuttle_rans.h"
#include "../rans_tables.h"

#define PASS(msg) do { printf("  [PASS] " msg "\n"); } while (0)
#define FAIL(msg, ...) do { printf("  [FAIL] " msg "\n", __VA_ARGS__); return 1; } while (0)

#if SHUTTLE_NIST_MODE == 120
#  define R_NUM_SYMS  SHUTTLE_NIST120_RANS_NUM_SYMS
#  define R_SYM_MAX   SHUTTLE_NIST120_RANS_SYM_MAX
#  define r_syms      shuttle_nist120_rans_syms
#  define r_freqs     shuttle_nist120_rans_freqs
#elif SHUTTLE_NIST_MODE == 180
#  define R_NUM_SYMS  SHUTTLE_NIST180_RANS_NUM_SYMS
#  define R_SYM_MAX   SHUTTLE_NIST180_RANS_SYM_MAX
#  define r_syms      shuttle_nist180_rans_syms
#  define r_freqs     shuttle_nist180_rans_freqs
#else
#  define R_NUM_SYMS  SHUTTLE_NIST255_RANS_NUM_SYMS
#  define R_SYM_MAX   SHUTTLE_NIST255_RANS_SYM_MAX
#  define r_syms      shuttle_nist255_rans_syms
#  define r_freqs     shuttle_nist255_rans_freqs
#endif

/* Discrete-Gaussian sampler on [-M_voc, M_voc] via cached inverse CDF.
 * Models the real signer's distribution (truncated Gaussian, NOT the
 * quantized rANS table). The rANS rate is sized against this PMF. */
#define SG_CDF_SIZE  64

static int32_t sample_discrete_gaussian(double sigma, int M_voc) {
  static int valid = 0;
  static double cached_sigma = 0.0;
  static int cached_M = 0;
  static double cdf[SG_CDF_SIZE];
  if (2 * M_voc + 1 > SG_CDF_SIZE) {
    fprintf(stderr, "sample_discrete_gaussian: M_voc=%d exceeds SG_CDF_SIZE\n", M_voc);
    abort();
  }
  if (!valid || sigma != cached_sigma || M_voc != cached_M) {
    cached_sigma = sigma;
    cached_M = M_voc;
    double Z = 0.0;
    for (int k = -M_voc; k <= M_voc; ++k)
      Z += exp(-(double)k * k / (2.0 * sigma * sigma));
    double running = 0.0;
    for (int k = -M_voc; k <= M_voc; ++k) {
      running += exp(-(double)k * k / (2.0 * sigma * sigma)) / Z;
      cdf[k + M_voc] = running;
    }
    valid = 1;
  }
  double u = ((double)rand() + 1.0) / ((double)RAND_MAX + 2.0);
  for (int i = 0; i <= 2 * M_voc; ++i)
    if (u <= cdf[i]) return (int32_t)(i - M_voc);
  return (int32_t)M_voc;
}

static inline int32_t sample_gauss_zhi(void) {
  return sample_discrete_gaussian(
      (double)SHUTTLE_SIGMA / (double)SHUTTLE_ALPHA_R, R_SYM_MAX);
}

/* Fill z_1[0..L] so the merged stream's HighBits land in-vocabulary. */
static void fill_z1(poly *z_1) {
  for (unsigned j = 0; j < SHUTTLE_N; ++j) {
    int32_t hi = sample_gauss_zhi();
    int32_t lo = (rand() % SHUTTLE_ALPHA_0P) - SHUTTLE_HALF_ALPHA_0P;
    z_1[0].coeffs[j] = (hi << SHUTTLE_ALPHA_0P_BITS) + lo;
  }
  for (unsigned i = 1; i <= SHUTTLE_L; ++i)
    for (unsigned j = 0; j < SHUTTLE_N; ++j) {
      int32_t hi = sample_gauss_zhi();
      int32_t lo = (rand() % SHUTTLE_ALPHA_R) - SHUTTLE_HALF_ALPHA_R;
      z_1[i].coeffs[j] = (hi << SHUTTLE_ALPHA_R_BITS) + lo;
    }
}

/* ------------------------------------------------------------ */
static int test_minimal(void) {
  printf("Test 1: minimal hand-crafted round-trip\n");

  uint8_t sig[SHUTTLE_BYTES];
  uint8_t c_tilde[SHUTTLE_CTILDEBYTES];
  int8_t  irs[SHUTTLE_TAU];
  poly    z_1[1 + SHUTTLE_L];
  polyveck h, h_rec;

  for (unsigned i = 0; i < SHUTTLE_CTILDEBYTES; ++i) c_tilde[i] = (uint8_t)(i * 17 + 3);
  for (unsigned i = 0; i < SHUTTLE_TAU; ++i) irs[i] = (i & 1) ? (int8_t)1 : (int8_t)-1;
  for (unsigned i = 0; i < 1 + SHUTTLE_L; ++i) memset(&z_1[i], 0, sizeof(poly));
  memset(&h, 0, sizeof h);

  int rc = pack_sig(sig, c_tilde, irs, z_1, &h);
  if (rc != 0) FAIL("pack_sig rc=%d", rc);

  uint8_t c_rec[SHUTTLE_CTILDEBYTES];
  int8_t  irs_rec[SHUTTLE_TAU];
  poly    z_rec[1 + SHUTTLE_L];
  rc = unpack_sig(c_rec, irs_rec, z_rec, &h_rec, sig);
  if (rc != 0) FAIL("unpack_sig rc=%d", rc);

  if (memcmp(c_tilde, c_rec, SHUTTLE_CTILDEBYTES) != 0)
    FAIL("c_tilde mismatch%s", "");
  for (unsigned i = 0; i < SHUTTLE_TAU; ++i)
    if (irs[i] != irs_rec[i])
      FAIL("irs_signs[%u]: %d vs %d", i, irs[i], irs_rec[i]);
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
  const size_t off_rans_len = SHUTTLE_CTILDEBYTES + SHUTTLE_IRS_SIGNBYTES;

  for (unsigned round = 0; round < 10; ++round) {
    uint8_t sig[SHUTTLE_BYTES];
    uint8_t c_tilde[SHUTTLE_CTILDEBYTES];
    int8_t  irs[SHUTTLE_TAU];
    poly    z_1[1 + SHUTTLE_L];
    polyveck h, h_rec;

    for (unsigned i = 0; i < SHUTTLE_CTILDEBYTES; ++i) c_tilde[i] = (uint8_t)rand();
    for (unsigned i = 0; i < SHUTTLE_TAU; ++i)
      irs[i] = (rand() & 1) ? (int8_t)1 : (int8_t)-1;

    fill_z1(z_1);
    /* hint shares the same PMF as z-hi (sigma_hint = sigma_z-hi). */
    for (unsigned i = 0; i < SHUTTLE_M; ++i)
      for (unsigned j = 0; j < SHUTTLE_N; ++j)
        h.vec[i].coeffs[j] = sample_gauss_zhi();

    int rc = pack_sig(sig, c_tilde, irs, z_1, &h);
    if (rc != 0) FAIL("round %u: pack rc=%d", round, rc);

    uint8_t c_rec[SHUTTLE_CTILDEBYTES];
    int8_t  irs_rec[SHUTTLE_TAU];
    poly    z_rec[1 + SHUTTLE_L];
    rc = unpack_sig(c_rec, irs_rec, z_rec, &h_rec, sig);
    if (rc != 0) FAIL("round %u: unpack rc=%d", round, rc);

    if (memcmp(c_tilde, c_rec, SHUTTLE_CTILDEBYTES) != 0) FAIL("round %u: c_tilde", round);
    for (unsigned i = 0; i < SHUTTLE_TAU; ++i)
      if (irs[i] != irs_rec[i])
        FAIL("round %u: irs[%u]", round, i);
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
int main(void) {
  int ret = 0;
  srand(0xC0DECAFEu ^ (unsigned)SHUTTLE_NIST_MODE);

  printf("=== test_sig (MODE=%d, SHUTTLE_BYTES=%d, RANS_RESERVED=%d) ===\n",
         SHUTTLE_NIST_MODE, SHUTTLE_BYTES, SHUTTLE_RANS_RESERVED_BYTES);

  ret |= test_minimal();
  ret |= test_random_rounds();

  if (ret == 0) printf("\n=== All sig tests PASSED ===\n");
  else          printf("\n=== Some sig tests FAILED ===\n");
  return ret;
}
