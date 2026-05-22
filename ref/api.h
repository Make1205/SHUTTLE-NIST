/*
 * api.h - NIST-style API header for SHUTTLE signature scheme.
 *
 * This header exposes the standard SUPERCOP/NIST entry points
 * (crypto_sign_keypair, crypto_sign_signature, crypto_sign_verify,
 * crypto_sign, crypto_sign_open) with names mangled per SHUTTLE_NIST_MODE,
 * so that libshuttle_nist120_ref.so, libshuttle_nist180_ref.so and libshuttle_nist255_ref.so
 * can coexist in the same linker namespace.
 *
 * A caller compiled against mode 120 will have their `crypto_sign_keypair`
 * references rewritten by the preprocessor to `shuttle_nist120_ref_keypair`,
 * which matches the symbol exported from libshuttle_nist120_ref.so.
 */

#ifndef SHUTTLE_API_H
#define SHUTTLE_API_H

#include <stddef.h>
#include <stdint.h>
#include "params.h"

#define CRYPTO_PUBLICKEYBYTES SHUTTLE_PUBLICKEYBYTES
#define CRYPTO_SECRETKEYBYTES SHUTTLE_SECRETKEYBYTES
#define CRYPTO_BYTES          SHUTTLE_BYTES
/* CRYPTO_ALGNAME ("SHUTTLE-NIST-120" / "SHUTTLE-NIST-180" / "SHUTTLE-NIST-255") comes from
 * config.h (included via params.h). */

/* NIST-style entry points; mangled to match the library symbols. */
#define crypto_sign_keypair   SHUTTLE_NAMESPACE(keypair)
int crypto_sign_keypair(uint8_t *pk, uint8_t *sk);

#define crypto_sign_signature SHUTTLE_NAMESPACE(signature)
int crypto_sign_signature(uint8_t *sig, size_t *siglen,
                          const uint8_t *m, size_t mlen,
                          const uint8_t *sk);

#define crypto_sign_verify    SHUTTLE_NAMESPACE(verify)
int crypto_sign_verify(const uint8_t *sig, size_t siglen,
                       const uint8_t *m, size_t mlen,
                       const uint8_t *pk);

#define crypto_sign           SHUTTLE_NAMESPACE(sign)
int crypto_sign(uint8_t *sm, size_t *smlen,
                const uint8_t *m, size_t mlen,
                const uint8_t *sk);

#define crypto_sign_open      SHUTTLE_NAMESPACE(open)
int crypto_sign_open(uint8_t *m, size_t *mlen,
                     const uint8_t *sm, size_t smlen,
                     const uint8_t *pk);

#endif /* SHUTTLE_API_H */
