#ifndef SHUTTLE_CONFIG_H
#define SHUTTLE_CONFIG_H

/* ============================================================
 * SHUTTLE-NIST parameter-set selector.
 *
 * Valid values (classical security level in bits):
 *   120 - SHUTTLE-NIST-120 (n=256, q=32257, sigma=67)
 *   180 - SHUTTLE-NIST-180 (n=256, q=15361, sigma=79)
 *   255 - SHUTTLE-NIST-255 (n=256, q=7681,  sigma=90)
 *
 * Override at build time with: -DSHUTTLE_NIST_MODE=120|180|255
 * ============================================================ */
#ifndef SHUTTLE_NIST_MODE
#define SHUTTLE_NIST_MODE 120
#endif

#if SHUTTLE_NIST_MODE == 120
#  define CRYPTO_ALGNAME       "SHUTTLE-NIST-120"
#  define SHUTTLE_NAMESPACETOP shuttle_nist120_ref
#  define SHUTTLE_NAMESPACE(s) shuttle_nist120_ref_##s
#elif SHUTTLE_NIST_MODE == 180
#  define CRYPTO_ALGNAME       "SHUTTLE-NIST-180"
#  define SHUTTLE_NAMESPACETOP shuttle_nist180_ref
#  define SHUTTLE_NAMESPACE(s) shuttle_nist180_ref_##s
#elif SHUTTLE_NIST_MODE == 255
#  define CRYPTO_ALGNAME       "SHUTTLE-NIST-255"
#  define SHUTTLE_NAMESPACETOP shuttle_nist255_ref
#  define SHUTTLE_NAMESPACE(s) shuttle_nist255_ref_##s
#else
#  error "Unsupported SHUTTLE_NIST_MODE (expected 120, 180, or 255)"
#endif

/* ============================================================
 * Discrete Gaussian sampler standard deviation (paper notation: r).
 *
 * SHUTTLE_SIGMA is derived from SHUTTLE_NIST_MODE.
 * ============================================================ */
#ifndef SHUTTLE_SIGMA
#  if SHUTTLE_NIST_MODE == 120
#    define SHUTTLE_SIGMA 67
#  elif SHUTTLE_NIST_MODE == 180
#    define SHUTTLE_SIGMA 79
#  elif SHUTTLE_NIST_MODE == 255
#    define SHUTTLE_SIGMA 90
#  endif
#endif

#endif
