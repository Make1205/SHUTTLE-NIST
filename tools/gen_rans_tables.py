#!/usr/bin/env python3
"""Generate static rANS frequency tables for SHUTTLE-NIST.

Background:

  A SINGLE rANS context encodes HighBits(z^(0..lenS)) and hint.  These
  blocks have comparable scales but not identical distributions: the
  response blocks are rounded discrete Gaussians, while MakeHint follows
  a bucket-crossing distribution.  The table is therefore generated from
  their coefficient-count-weighted mixture.

  Vocabulary radius:
    M_voc = max(
      ceil((11r + tau*eta)/alpha_r),               # HighBits(z^(i>=1))
      ceil((ceil(11r/alpha_1) + 1)/alpha_0'),      # HighBits(z^(0))
      floor(2*(11r + tau*eta)/alpha_h) + 1         # hint
    )

  In SHUTTLE-NIST with alpha_h = 2*alpha_r, the hint and z^(i>=1) bounds
  coincide (both = ceil((11r+tau*eta)/alpha_r)), so the max is taken
  over just two terms.

Quantization procedure:
  1. PMF eval: compute the rounded-response and bucket-crossing PMFs on
     [-M_voc, M_voc], then form their coefficient-count-weighted mixture.
  2. Round: g(s) = max(1, round(p(s) * 2^t)) with t = 10.
  3. Slack: adjust the largest bucket so sum(g) = 2^t exactly.

Usage:
  python3 gen_rans_tables.py --out ../ref/rans_tables.h
  python3 gen_rans_tables.py --mode 120 --format text   # human inspection
"""

from __future__ import annotations

import argparse
import math
import sys
from typing import Dict, List, Optional, Tuple


MODE_PARAMS: Dict[int, dict] = {
    120: {"n": 256, "lenS": 3, "lenE": 2, "eta": 1, "tau":  58,
          "r":  70, "alpha_1":  8, "alpha_h": 128},
    180: {"n": 256, "lenS": 4, "lenE": 3, "eta": 1, "tau":  80,
          "r":  82, "alpha_1":  8, "alpha_h": 128},
    255: {"n": 256, "lenS": 5, "lenE": 4, "eta": 1, "tau": 128,
          "r":  99, "alpha_1": 16, "alpha_h": 128},
}

ALPHA_R = 64         # cross-mode constant
T_SIGMA = 11         # SampleY truncation multiplier
PROB_BITS = 10       # quantization bits


def theoretical_pmf(sigma: float, M_voc: int) -> Dict[int, float]:
    """Discrete Gaussian PMF on the finite alphabet [-M_voc, M_voc]."""
    weights = {k: math.exp(-k * k / (2 * sigma * sigma))
               for k in range(-M_voc, M_voc + 1)}
    Z = sum(weights.values())
    return {k: w / Z for k, w in weights.items()}


def round_half_up_ratio(x: int, step: int) -> int:
    """Integer round-half-up of x/step, matching the C implementation."""
    return math.floor((x + step // 2) / step)


def rounded_response_pmf(r: int, step: int, M_voc: int,
                         inner_step: Optional[int] = None) -> Dict[int, float]:
    """PMF after the response coefficient is rounded into an rANS symbol.

    For z^(i), i >= 1, ``inner_step`` is absent and the map is
    RoundHU(x/alpha_r).  For z^(0), CompressY first applies
    RoundHU(x/alpha_1), after which the rANS split applies
    RoundHU(./alpha_0').
    """
    bound = T_SIGMA * r
    weights = {x: math.exp(-x * x / (2 * r * r))
               for x in range(-bound, bound + 1)}
    Z = sum(weights.values())
    pmf = {s: 0.0 for s in range(-M_voc, M_voc + 1)}
    for x, weight in weights.items():
        y = x if inner_step is None else round_half_up_ratio(x, inner_step)
        symbol = round_half_up_ratio(y, step)
        if symbol in pmf:
            pmf[symbol] += weight / Z
    total = sum(pmf.values())
    return {s: probability / total for s, probability in pmf.items()}


def hint_pmf(r: int, alpha_h: int, M_voc: int) -> Dict[int, float]:
    """Bucket-crossing PMF of one centered MakeHint coefficient."""
    bound = T_SIGMA * r
    weights = {z: math.exp(-z * z / (2 * r * r))
               for z in range(-bound, bound + 1)}
    Z = sum(weights.values())
    pmf = {h: 0.0 for h in range(-M_voc, M_voc + 1)}
    for z, weight in weights.items():
        delta = 2.0 * z / alpha_h
        k = math.floor(delta + 0.5)
        frac = delta - k
        if k in pmf:
            pmf[k] += (weight / Z) * (1.0 - abs(frac))
        if abs(frac) > 1e-15:
            adjacent = k + (1 if frac > 0 else -1)
            if adjacent in pmf:
                pmf[adjacent] += (weight / Z) * abs(frac)
    total = sum(pmf.values())
    return {h: probability / total for h, probability in pmf.items()}


def merged_symbol_pmf(params: dict, M_voc: int) -> Dict[int, float]:
    """Coefficient-count-weighted PMF of the merged rANS stream."""
    r = params["r"]
    lenS, lenE = params["lenS"], params["lenE"]
    alpha_1, alpha_h = params["alpha_1"], params["alpha_h"]
    alpha_0p = ALPHA_R // alpha_1
    p_z0 = rounded_response_pmf(r, alpha_0p, M_voc,
                                inner_step=alpha_1)
    p_zi = rounded_response_pmf(r, ALPHA_R, M_voc)
    p_h = hint_pmf(r, alpha_h, M_voc)
    block_count = 1 + lenS + lenE
    return {
        s: (p_z0[s] + lenS * p_zi[s] + lenE * p_h[s]) / block_count
        for s in p_z0
    }


def quantize(pmf: Dict[int, float], prob_bits: int) -> Tuple[List[int], List[int]]:
    """Quantize a PMF to integer frequencies summing to 2^prob_bits."""
    prob_total = 1 << prob_bits
    syms = sorted(pmf.keys())
    real_fs = [pmf[s] * prob_total for s in syms]
    freqs = [max(1, int(round(f))) for f in real_fs]
    diff = prob_total - sum(freqs)
    if diff != 0:
        max_idx = max(range(len(freqs)), key=lambda i: freqs[i])
        freqs[max_idx] += diff
        if freqs[max_idx] < 1:
            raise RuntimeError("quantization produced non-positive freq")
    assert sum(freqs) == prob_total
    return syms, freqs


def merged_table(params: dict) -> Tuple[List[int], List[int], float, int]:
    """Single shared rANS table for HighBits(z) || hint.

    The table probabilities are obtained from the actual rounded-response
    and bucket-crossing models, weighted by their coefficient counts.

    M_voc is the configured 11r-based design radius obtained from:
      - HighBits_{alpha_r}(z^(i)), i >= 1:  ceil((11r + tau*eta)/alpha_r)
      - HighBits_{alpha_0'}(z^(0)):         ceil((ceil(11r/alpha_1)+1)/alpha_0')
      - hint h: floor(2(11r + tau*eta)/alpha_h) + 1

    Symbols outside this design vocabulary are rejected by the encoder;
    their probability is accounted for separately from capacity overflow.
    """
    r = params["r"]
    tau, eta = params["tau"], params["eta"]
    alpha_1 = params["alpha_1"]
    alpha_0p = ALPHA_R // alpha_1
    sigma = r / ALPHA_R
    M_voc_zi = math.ceil((T_SIGMA * r + tau * eta) / ALPHA_R)
    M_voc_z0 = math.ceil((math.ceil(T_SIGMA * r / alpha_1) + 1) / alpha_0p)
    M_voc_h  = (2 * (T_SIGMA * r + tau * eta)) // params["alpha_h"] + 1
    M_voc = max(M_voc_zi, M_voc_z0, M_voc_h)
    syms, freqs = quantize(merged_symbol_pmf(params, M_voc), PROB_BITS)
    return syms, freqs, sigma, M_voc


def emit_table_block(mode: int,
                     syms: List[int], freqs: List[int],
                     prob_bits: int, sigma: float, file) -> None:
    """Write the single per-mode rANS frequency block."""
    prob_total = 1 << prob_bits
    macro_pref = f"SHUTTLE_NIST{mode}_RANS"
    var_pref = f"shuttle_nist{mode}_rans"
    print(f"/* ================================================================ */",
          file=file)
    print(f"/* SHUTTLE-NIST-{mode} rANS (merged z-hi || hint) frequency table.    */",
          file=file)
    print(f"/*   sigma = {sigma:.4f}, M_voc = {syms[-1]} (alphabet [-M_voc, M_voc]) */",
          file=file)
    print(f"/*   prob_bits = {prob_bits}, sum(freqs) = {prob_total}.                */",
          file=file)
    print(f"/* Generated by tools/gen_rans_tables.py from the analytic mixture PMF. */",
          file=file)
    print(f"/* ================================================================ */",
          file=file)
    print(f"#define {macro_pref}_PROB_BITS {prob_bits}", file=file)
    print(f"#define {macro_pref}_SYM_MIN  ({syms[0]})", file=file)
    print(f"#define {macro_pref}_SYM_MAX  ({syms[-1]})", file=file)
    print(f"#define {macro_pref}_NUM_SYMS {len(syms)}", file=file)
    print(f"", file=file)
    print(f"static const int16_t {var_pref}_syms[{len(syms)}] = {{", file=file)
    for i in range(0, len(syms), 12):
        line = ", ".join(f"{s:4d}" for s in syms[i:i+12])
        print(f"    {line},", file=file)
    print(f"}};", file=file)
    print(f"", file=file)
    print(f"static const uint16_t {var_pref}_freqs[{len(syms)}] = {{", file=file)
    for i in range(0, len(freqs), 12):
        line = ", ".join(f"{f:5d}" for f in freqs[i:i+12])
        print(f"    {line},", file=file)
    print(f"}};", file=file)
    print(f"", file=file)


def emit_header(modes: List[int], file) -> None:
    print("/* SHUTTLE-NIST rANS frequency tables.", file=file)
    print(" * Auto-generated by SHUTTLE-NIST/tools/gen_rans_tables.py.", file=file)
    print(" *", file=file)
    print(" * SHUTTLE-NIST uses design S (single merged rANS stream): HighBits", file=file)
    print(" * of z^(0..lenS) and hint h use one coefficient-count-weighted", file=file)
    print(" * mixture table. Response symbols use their exact rounded", file=file)
    print(" * discrete-Gaussian marginals and hints use the bucket-crossing PMF.", file=file)
    print(" *", file=file)
    print(" * Tables come from the analytic mixture PMF; the", file=file)
    print(" * alphabet uses the configured max(M_voc terms) design bound; the", file=file)
    print(" * encoder rejects any symbol outside this range.", file=file)
    print(" *", file=file)
    print(" * Do not edit by hand -- regenerate with:", file=file)
    print(" *   python3 SHUTTLE-NIST/tools/gen_rans_tables.py --out SHUTTLE-NIST/ref/rans_tables.h", file=file)
    print(" */", file=file)
    print("#ifndef SHUTTLE_RANS_TABLES_H", file=file)
    print("#define SHUTTLE_RANS_TABLES_H", file=file)
    print("#include <stdint.h>", file=file)
    print("", file=file)


def emit_footer(file) -> None:
    print("#endif /* SHUTTLE_RANS_TABLES_H */", file=file)


def text_summary(mode: int, params: dict) -> None:
    print(f"=== SHUTTLE-NIST-{mode} ===")
    syms, freqs, sigma, M = merged_table(params)
    print(f"merged: sigma = {sigma:.4f}, M_voc = {M}, alphabet size = {len(syms)}")
    print(f"top buckets:")
    for s, f in sorted(zip(syms, freqs), key=lambda t: -t[1])[:10]:
        print(f"  s={s:+4d} g={f:5d} p={f/1024:.4f}")


def main() -> None:
    ap = argparse.ArgumentParser(description="rANS theoretical-PMF table generator")
    ap.add_argument("--mode", type=int, action="append", choices=[120, 180, 255],
                    help="restrict to one mode; default emits all three")
    ap.add_argument("--format", choices=["text", "c"], default="c")
    ap.add_argument("--out", default=None,
                    help="output path (default stdout; only meaningful for --format c)")
    args = ap.parse_args()

    modes = sorted(set(args.mode)) if args.mode else [120, 180, 255]

    if args.format == "text":
        for m in modes:
            text_summary(m, MODE_PARAMS[m])
        return

    outf = sys.stdout if args.out is None else open(args.out, "w")
    emit_header(modes, outf)
    for m in modes:
        params = MODE_PARAMS[m]
        syms, freqs, sigma, _ = merged_table(params)
        emit_table_block(m, syms, freqs, PROB_BITS, sigma, outf)
    emit_footer(outf)
    if outf is not sys.stdout:
        outf.close()
        print(f"Wrote {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
