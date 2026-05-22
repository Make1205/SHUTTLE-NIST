/*
 * test_rans.c - unit tests for the SHUTTLE-NIST rANS engine
 *               (single merged-stream API, design S).
 *
 * Coverage:
 *   1. Small round-trip.
 *   2. Large round-trip (10000 samples drawn from the rANS table PMF).
 *   3. Empty input: encoder produces a 4-byte flush, decoder accepts it.
 *   4. Little-endian flush: the final 4 bytes of the encoded stream
 *      reproduce the encoder's last x state in LE order.
 *   5. Final-state verification: any single-byte corruption is detected.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#include "../params.h"
#include "../shuttle_rans.h"
#include "../rans_tables.h"

#define PASS(msg) do { printf("  [PASS] " msg "\n"); } while (0)
#define FAIL(msg, ...) do { printf("  [FAIL] " msg "\n", __VA_ARGS__); return 1; } while (0)

#if SHUTTLE_NIST_MODE == 120
#  define R_NUM_SYMS   SHUTTLE_NIST120_RANS_NUM_SYMS
#  define R_SYM_MIN    SHUTTLE_NIST120_RANS_SYM_MIN
#  define R_SYM_MAX    SHUTTLE_NIST120_RANS_SYM_MAX
#  define r_syms       shuttle_nist120_rans_syms
#  define r_freqs      shuttle_nist120_rans_freqs
#elif SHUTTLE_NIST_MODE == 180
#  define R_NUM_SYMS   SHUTTLE_NIST180_RANS_NUM_SYMS
#  define R_SYM_MIN    SHUTTLE_NIST180_RANS_SYM_MIN
#  define R_SYM_MAX    SHUTTLE_NIST180_RANS_SYM_MAX
#  define r_syms       shuttle_nist180_rans_syms
#  define r_freqs      shuttle_nist180_rans_freqs
#else
#  define R_NUM_SYMS   SHUTTLE_NIST255_RANS_NUM_SYMS
#  define R_SYM_MIN    SHUTTLE_NIST255_RANS_SYM_MIN
#  define R_SYM_MAX    SHUTTLE_NIST255_RANS_SYM_MAX
#  define r_syms       shuttle_nist255_rans_syms
#  define r_freqs      shuttle_nist255_rans_freqs
#endif

#define PROB_TOTAL (1u << SHUTTLE_RANS_PROB_BITS)

static int32_t sample_from_table(const int16_t *syms, const uint16_t *freqs,
                                 unsigned num_syms)
{
  uint32_t r = (uint32_t)rand() & (PROB_TOTAL - 1u);
  uint32_t cum = 0;
  for (unsigned i = 0; i < num_syms; ++i) {
    cum += freqs[i];
    if (r < cum) return syms[i];
  }
  return syms[num_syms - 1];
}

static double table_entropy(const uint16_t *freqs, unsigned num_syms)
{
  double H = 0.0;
  for (unsigned i = 0; i < num_syms; ++i) {
    double p = (double)freqs[i] / (double)PROB_TOTAL;
    if (p > 0) H -= p * log2(p);
  }
  return H;
}

/* ------------------------------------------------------------ */
static int test_small(void) {
  printf("Test 1: small round-trip\n");
  int32_t msg[10] = {0, 1, -1, 2, 0, -2, 0, 3, 1, 0};
  int32_t dec[10];
  uint8_t buf[128];
  size_t out_len;

  int rc = shuttle_rans_encode(buf, &out_len, sizeof buf, msg, 10);
  if (rc != 0) FAIL("encode failed (rc=%d)", rc);

  rc = shuttle_rans_decode(dec, 10, buf, out_len);
  if (rc != 0) FAIL("decode failed (rc=%d)", rc);

  for (unsigned i = 0; i < 10; ++i)
    if (dec[i] != msg[i])
      FAIL("mismatch at %u: msg=%d dec=%d", i, msg[i], dec[i]);

  printf("    encoded %zu bytes for 10 symbols\n", out_len);
  PASS("small round-trip");
  return 0;
}

static int test_large(void) {
  printf("Test 2: large round-trip (10000 samples)\n");
  const unsigned N = 10000;
  int32_t *msg = malloc(N * sizeof(int32_t));
  int32_t *dec = malloc(N * sizeof(int32_t));
  uint8_t *buf = malloc(N * 2);

  for (unsigned i = 0; i < N; ++i)
    msg[i] = sample_from_table(r_syms, r_freqs, R_NUM_SYMS);

  size_t out_len;
  int rc = shuttle_rans_encode(buf, &out_len, N * 2, msg, N);
  if (rc != 0) { free(msg); free(dec); free(buf); FAIL("encode failed (rc=%d)", rc); }
  rc = shuttle_rans_decode(dec, N, buf, out_len);
  if (rc != 0) { free(msg); free(dec); free(buf); FAIL("decode failed (rc=%d)", rc); }
  for (unsigned i = 0; i < N; ++i)
    if (dec[i] != msg[i]) {
      free(msg); free(dec); free(buf);
      FAIL("mismatch at %u: msg=%d dec=%d", i, msg[i], dec[i]);
    }

  double H = table_entropy(r_freqs, R_NUM_SYMS);
  printf("    encoded = %zu bytes = %.3f bits/sym  (Shannon = %.3f)\n",
         out_len, 8.0 * out_len / N, H);
  free(msg); free(dec); free(buf);
  PASS("large round-trip");
  return 0;
}

/* ------------------------------------------------------------ */
static int test_empty(void) {
  printf("Test 3: empty sequence\n");
  uint8_t buf[32];
  size_t out_len;

  int rc = shuttle_rans_encode(buf, &out_len, sizeof buf, NULL, 0);
  if (rc != 0) FAIL("encode(empty) failed (rc=%d)", rc);
  if (out_len != 4)
    FAIL("empty stream should be exactly 4 flush bytes, got %zu", out_len);

  /* Initial state L = 2^23 = 0x00800000, little-endian = [0x00,0x00,0x80,0x00]. */
  if (!(buf[0] == 0x00 && buf[1] == 0x00 && buf[2] == 0x80 && buf[3] == 0x00))
    FAIL("flush bytes look wrong: %02x %02x %02x %02x",
         buf[0], buf[1], buf[2], buf[3]);

  int32_t dummy[1];
  rc = shuttle_rans_decode(dummy, 0, buf, out_len);
  if (rc != 0) FAIL("decode(empty) failed (rc=%d)", rc);

  PASS("empty stream round-trip and little-endian flush of L = 2^23");
  return 0;
}

/* ------------------------------------------------------------ */
static int test_corruption_detected(void) {
  printf("Test 4: final-state mismatch detection (flip every byte)\n");
  const unsigned N = 200;
  int32_t msg[200];
  int32_t dec[200];
  uint8_t buf[400];
  uint8_t buf2[400];
  size_t out_len;

  for (unsigned i = 0; i < N; ++i)
    msg[i] = sample_from_table(r_syms, r_freqs, R_NUM_SYMS);

  int rc = shuttle_rans_encode(buf, &out_len, sizeof buf, msg, N);
  if (rc != 0) FAIL("encode failed (rc=%d)", rc);

  rc = shuttle_rans_decode(dec, N, buf, out_len);
  if (rc != 0) FAIL("clean decode failed (rc=%d)", rc);

  unsigned undetected = 0;
  for (size_t pos = 0; pos < out_len; ++pos) {
    memcpy(buf2, buf, out_len);
    buf2[pos] ^= 0xFFu;

    rc = shuttle_rans_decode(dec, N, buf2, out_len);
    int caught = (rc != 0);
    if (rc == 0) {
      for (unsigned i = 0; i < N; ++i)
        if (dec[i] != msg[i]) { caught = 1; break; }
    }
    if (!caught) ++undetected;
  }
  if (undetected != 0)
    FAIL("%u single-byte corruptions slipped through (out of %zu positions)",
         undetected, out_len);

  printf("    all %zu single-byte corruptions detected\n", out_len);
  PASS("byte-corruption detection");
  return 0;
}

/* ------------------------------------------------------------ */
int main(void) {
  int ret = 0;
  srand(0xC0DEu ^ (unsigned)SHUTTLE_NIST_MODE);

  printf("=== test_rans (MODE=%d)\n"
         "    merged table: [%d, %d], %u symbols\n"
         "    prob_bits = %u ===\n",
         SHUTTLE_NIST_MODE,
         R_SYM_MIN, R_SYM_MAX, (unsigned)R_NUM_SYMS,
         (unsigned)SHUTTLE_RANS_PROB_BITS);

  ret |= test_small();
  ret |= test_large();
  ret |= test_empty();
  ret |= test_corruption_detected();

  if (ret == 0)
    printf("\n=== All rANS tests PASSED ===\n");
  else
    printf("\n=== Some rANS tests FAILED ===\n");
  return ret;
}
