#ifndef BENCH_ADAPTER_H
#define BENCH_ADAPTER_H
#include "api.h"
#define B_N SHUTTLE_N
#define B_Q SHUTTLE_Q
#define B_PK CRYPTO_PUBLICKEYBYTES
#define B_SK CRYPTO_SECRETKEYBYTES
#define B_SIG CRYPTO_BYTES
#define b_keygen crypto_sign_keypair
#define b_sign crypto_sign_signature
#define b_verify crypto_sign_verify
static size_t b_unpadded(const uint8_t*p,size_t n){size_t a=SHUTTLE_CTILDEBYTES;size_t r=p[a]|((size_t)p[a+1]<<8);return n-SHUTTLE_RANS_RESERVED_BYTES+r;}
#endif
