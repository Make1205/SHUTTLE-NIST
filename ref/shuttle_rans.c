/*
 * shuttle_rans.c - Byte-wise rANS encoder/decoder for SHUTTLE-NIST.
 *
 * Single merged context: one rANS table encodes HighBits(z) followed by
 * hint h, because the parameter set takes alpha_h = 2*alpha_r so
 * sigma_hint = sigma_z-hi.
 *
 * State is uint32_t; bytes are emitted to the TAIL of the output buffer
 * during encode (LIFO), then compacted with memmove so the final layout
 * is a prefix of length *out_len. The 4-byte flush at the end is written
 * little-endian. Decoding reads the leading 4 little-endian bytes to
 * prime the state, then walks the stream.
 *
 * Canonical decode checks: after the last symbol is consumed, x must equal
 * L_ren (= 2^23) and the input must be consumed exactly.
 */

#include <string.h>

#include "params.h"
#include "shuttle_rans.h"
#include "rans_tables.h"

/* ============================================================
 * Per-mode table bindings.
 * ============================================================ */

#if SHUTTLE_NIST_MODE == 120
#  define SRANS_NUM_SYMS    SHUTTLE_NIST120_RANS_NUM_SYMS
#  define SRANS_SYM_MIN     SHUTTLE_NIST120_RANS_SYM_MIN
#  define SRANS_SYM_MAX     SHUTTLE_NIST120_RANS_SYM_MAX
#  define srans_syms        shuttle_nist120_rans_syms
#  define srans_freqs       shuttle_nist120_rans_freqs
#elif SHUTTLE_NIST_MODE == 180
#  define SRANS_NUM_SYMS    SHUTTLE_NIST180_RANS_NUM_SYMS
#  define SRANS_SYM_MIN     SHUTTLE_NIST180_RANS_SYM_MIN
#  define SRANS_SYM_MAX     SHUTTLE_NIST180_RANS_SYM_MAX
#  define srans_syms        shuttle_nist180_rans_syms
#  define srans_freqs       shuttle_nist180_rans_freqs
#elif SHUTTLE_NIST_MODE == 255
#  define SRANS_NUM_SYMS    SHUTTLE_NIST255_RANS_NUM_SYMS
#  define SRANS_SYM_MIN     SHUTTLE_NIST255_RANS_SYM_MIN
#  define SRANS_SYM_MAX     SHUTTLE_NIST255_RANS_SYM_MAX
#  define srans_syms        shuttle_nist255_rans_syms
#  define srans_freqs       shuttle_nist255_rans_freqs
#else
#  error "Unsupported SHUTTLE_NIST_MODE for rANS tables"
#endif

/* ============================================================
 * Per-mode context + lazy init.
 * ============================================================ */

typedef struct {
    const int16_t  *syms;       /* length num_syms, contiguous SYM_MIN..SYM_MAX */
    const uint16_t *freqs;      /* length num_syms, sum == 2^PROB_BITS */
    uint16_t       *cdf;        /* length num_syms+1, lazy-built */
    uint16_t       *sym_lookup; /* length 2^PROB_BITS, lazy-built */
    int16_t         sym_min;
    int16_t         sym_max;
    uint16_t        num_syms;
    uint8_t         initialized;
} rans_ctx_t;

static uint16_t g_cdf[SRANS_NUM_SYMS + 1];
static uint16_t g_sym_lookup[SHUTTLE_RANS_PROB_TOTAL];
static rans_ctx_t g_ctx = {
    .syms       = srans_syms,
    .freqs      = srans_freqs,
    .cdf        = g_cdf,
    .sym_lookup = g_sym_lookup,
    .sym_min    = SRANS_SYM_MIN,
    .sym_max    = SRANS_SYM_MAX,
    .num_syms   = SRANS_NUM_SYMS,
    .initialized = 0,
};

static void rans_ctx_init(rans_ctx_t *ctx) {
    uint32_t running;
    unsigned i, k;

    if (ctx->initialized) return;

    running = 0;
    ctx->cdf[0] = 0;
    for (i = 0; i < ctx->num_syms; ++i) {
        running += ctx->freqs[i];
        ctx->cdf[i + 1] = (uint16_t)running;
    }

    for (i = 0; i < ctx->num_syms; ++i) {
        uint32_t start = ctx->cdf[i];
        uint32_t end   = ctx->cdf[i + 1];
        for (k = start; k < end; ++k)
            ctx->sym_lookup[k] = (uint16_t)i;
    }

    ctx->initialized = 1;
}

/* Map a signed integer symbol to its index in ctx->syms. */
static inline int sym_to_index(const rans_ctx_t *ctx, int32_t sym) {
    int32_t idx = sym - ctx->sym_min;
    if (idx < 0 || idx >= (int32_t)ctx->num_syms || ctx->syms[idx] != sym)
        return -1;
    return (int)idx;
}

/* ============================================================
 * Core encoder / decoder
 * ============================================================ */

static int encode_core(rans_ctx_t *ctx,
                       uint8_t *out, size_t *out_len, size_t max_bytes,
                       const int32_t *syms, size_t n)
{
    uint32_t x = SHUTTLE_RANS_L;
    size_t idx_i;
    size_t write_pos;

    rans_ctx_init(ctx);

    write_pos = max_bytes;

    for (idx_i = n; idx_i > 0; --idx_i) {
        int si = sym_to_index(ctx, syms[idx_i - 1]);
        if (si < 0) return -1;

        uint32_t freq  = ctx->freqs[si];
        uint32_t start = ctx->cdf[si];

        uint32_t x_max = ((SHUTTLE_RANS_L >> SHUTTLE_RANS_PROB_BITS)
                          << SHUTTLE_RANS_BYTE_BITS) * freq;
        while (x >= x_max) {
            if (write_pos == 0) return -2;
            out[--write_pos] = (uint8_t)(x & 0xFFu);
            x >>= SHUTTLE_RANS_BYTE_BITS;
        }

        x = ((x / freq) << SHUTTLE_RANS_PROB_BITS) + start + (x % freq);
    }

    /* Final flush: write x as 4 little-endian bytes. */
    if (write_pos < 4) return -2;
    write_pos -= 4;
    out[write_pos + 0] = (uint8_t)(x >>  0);
    out[write_pos + 1] = (uint8_t)(x >>  8);
    out[write_pos + 2] = (uint8_t)(x >> 16);
    out[write_pos + 3] = (uint8_t)(x >> 24);

    {
        size_t bytes_used = max_bytes - write_pos;
        memmove(out, out + write_pos, bytes_used);
        *out_len = bytes_used;
    }
    return 0;
}

static int decode_core(rans_ctx_t *ctx,
                       int32_t *syms, size_t n,
                       const uint8_t *in, size_t in_len)
{
    rans_ctx_init(ctx);

    if (in_len < 4) return -1;

    uint32_t x = ((uint32_t)in[0] <<  0)
               | ((uint32_t)in[1] <<  8)
               | ((uint32_t)in[2] << 16)
               | ((uint32_t)in[3] << 24);
    size_t pos = 4;

    const uint32_t prob_mask = SHUTTLE_RANS_PROB_TOTAL - 1u;

    for (size_t k = 0; k < n; ++k) {
        uint32_t c = x & prob_mask;
        unsigned si = ctx->sym_lookup[c];

        uint32_t freq  = ctx->freqs[si];
        uint32_t start = ctx->cdf[si];

        syms[k] = ctx->syms[si];

        x = freq * (x >> SHUTTLE_RANS_PROB_BITS) + (c - start);

        while (x < SHUTTLE_RANS_L) {
            if (pos >= in_len) return -1;
            x = (x << SHUTTLE_RANS_BYTE_BITS) | in[pos];
            ++pos;
        }
    }

    if (x != SHUTTLE_RANS_L || pos != in_len) return -1;

    return 0;
}

/* ============================================================
 * Public API
 * ============================================================ */

int shuttle_rans_encode(uint8_t *out, size_t *out_len, size_t max_bytes,
                        const int32_t *syms, size_t n)
{
    return encode_core(&g_ctx, out, out_len, max_bytes, syms, n);
}

int shuttle_rans_decode(int32_t *syms, size_t n,
                        const uint8_t *in, size_t in_len)
{
    return decode_core(&g_ctx, syms, n, in, in_len);
}
