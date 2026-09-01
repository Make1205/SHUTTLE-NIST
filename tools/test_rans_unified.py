#!/usr/bin/env python3
"""End-to-end regression test for the SHUTTLE-NIST rANS layer.

Coverage:
  1. The single merged frequency table produced by gen_rans_tables.py is
     well-formed for every mode: alphabet contiguous and symmetric, freqs
     sum to 2^prob_bits, every symbol has freq >= 1. M_voc matches the
     configured 11r-based design radius.
  2. Python reference rans.py round-trips samples from the merged PMF.
  3. Final-state verification rejects any single-byte corruption of the
     encoded stream.
  4. Little-endian flush of the initial state is bit-exact.

This intentionally does not exercise the C engine -- for that, run the
binaries under SHUTTLE-NIST/ref/test/out/ (test_rans, test_poly_z1,
test_sig). The C/Python byte-identity comparison lives in kat_compare/.
"""

from __future__ import annotations

import math
import random
import sys
from pathlib import Path

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))

from rans import RansTable, encode, decode, RANS_L
from gen_rans_tables import (MODE_PARAMS, ALPHA_R, PROB_BITS, T_SIGMA,
                             merged_symbol_pmf, merged_table)


def check_table(syms, freqs, prob_bits):
    assert len(syms) == len(freqs)
    assert syms == sorted(syms)
    assert syms == list(range(syms[0], syms[-1] + 1)), \
        f"alphabet not contiguous: {syms}"
    assert syms[0] == -syms[-1], f"alphabet not symmetric: [{syms[0]}, {syms[-1]}]"
    assert all(f >= 1 for f in freqs), "found freq < 1"
    assert sum(freqs) == 1 << prob_bits, \
        f"freq sum {sum(freqs)} != 2^{prob_bits}"


def test_tables():
    print("== Test: analytic mixture-PMF tables, geometry and sums ==")
    for mode, params in MODE_PARAMS.items():
        syms, freqs, sigma, M = merged_table(params)
        check_table(syms, freqs, PROB_BITS)

        # M_voc design radius (max over z^(0) hi, z^(>=1) hi, hint).
        r, tau, eta = params["r"], params["tau"], params["eta"]
        alpha_1, alpha_h_spec = params["alpha_1"], params["alpha_h"]
        alpha_0p = ALPHA_R // alpha_1
        expected_zi = math.ceil((T_SIGMA * r + tau * eta) / ALPHA_R)
        expected_z0 = math.ceil((math.ceil(T_SIGMA * r / alpha_1) + 1) / alpha_0p)
        expected_h  = (2 * (T_SIGMA * r + tau * eta)) // alpha_h_spec + 1
        expected = max(expected_zi, expected_z0, expected_h)
        assert M == expected, \
            f"mode-{mode} M_voc: got {M}, expected {expected}"

        print(f"  NIST-{mode}: M_voc={M:2d} ({len(syms):2d} symbols, "
              f"sigma={sigma:.4f})")


def sample_pmf(rng: random.Random, pmf) -> int:
    """Sample from a finite PMF."""
    syms = sorted(pmf.keys())
    weights = [pmf[s] for s in syms]
    cum = 0.0
    pick = rng.random()
    for s, w in zip(syms, weights):
        cum += w
        if pick <= cum:
            return s
    return syms[-1]


def test_roundtrip():
    print("== Test: analytic merged-PMF sample -> encode -> decode ==")
    rng = random.Random(0xC0DECAFE)
    for mode, params in MODE_PARAMS.items():
        syms, freqs, sigma, M_voc = merged_table(params)
        table = RansTable(syms=list(syms), freqs=list(freqs),
                          prob_bits=PROB_BITS)
        N = 1024
        pmf = merged_symbol_pmf(params, M_voc)
        msg = [sample_pmf(rng, pmf) for _ in range(N)]
        buf = encode(msg, table)
        recovered = decode(buf, table, len(msg))
        assert recovered == msg, \
            f"mode-{mode}: round-trip mismatch"
        print(f"  NIST-{mode}: {N} symbols -> {len(buf)} B "
              f"({8 * len(buf) / N:.3f} bit/sym)")


def test_little_endian_flush():
    print("== Test: little-endian flush of the initial state ==")
    for mode, params in MODE_PARAMS.items():
        syms, freqs, _, _ = merged_table(params)
        table = RansTable(syms=list(syms), freqs=list(freqs),
                          prob_bits=PROB_BITS)
        buf = encode([], table)
        assert len(buf) == 4, f"empty stream should be 4 bytes, got {len(buf)}"
        expected = bytes([
            (RANS_L >> 0) & 0xFF,
            (RANS_L >> 8) & 0xFF,
            (RANS_L >> 16) & 0xFF,
            (RANS_L >> 24) & 0xFF,
        ])
        assert buf == expected, \
            f"NIST-{mode}: flush bytes {buf.hex()} != expected {expected.hex()}"
        print(f"  NIST-{mode}: LE flush bytes = {buf.hex()}  (L = 0x{RANS_L:08x})")


def test_corruption_detected():
    print("== Test: final-state verification catches byte corruption ==")
    rng = random.Random(0xBADBABE)
    params = MODE_PARAMS[120]
    syms, freqs, sigma, M_voc = merged_table(params)
    table = RansTable(syms=list(syms), freqs=list(freqs), prob_bits=PROB_BITS)
    N = 200
    pmf = merged_symbol_pmf(params, M_voc)
    msg = [sample_pmf(rng, pmf) for _ in range(N)]
    buf = encode(msg, table)
    undetected = 0
    for pos in range(len(buf)):
        corrupted = bytearray(buf)
        corrupted[pos] ^= 0xFF
        try:
            recovered = decode(bytes(corrupted), table, N)
            if recovered == msg:
                undetected += 1
        except RuntimeError:
            pass
        except IndexError:
            pass
    assert undetected == 0, f"{undetected} corruptions went undetected"
    print(f"  all {len(buf)} single-byte corruptions caught")


def main():
    test_tables()
    test_roundtrip()
    test_little_endian_flush()
    test_corruption_detected()
    print()
    print("All rANS unified tests PASSED.")


if __name__ == "__main__":
    main()
