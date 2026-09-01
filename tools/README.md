# SHUTTLE Python Tooling

Offline helpers for the SHUTTLE rANS layer: analytic frequency-table
generation, a pure-Python reference encoder/decoder used to cross-check the
C implementation, and a KAT comparator for byte-identical regression.

## Files

| file | role |
|---|---|
| `rans.py` | Pure-Python rANS encoder/decoder (reference; cross-checks the C side). |
| `gen_rans_tables.py` | Emit `ref/rans_tables.h` from the analytic merged PMF. |
| `SigSize.py` | Recompute the signature-size decomposition and rANS reservations. |
| `kat_compare/` | KAT regression helpers. |

## Dependencies

Python 3.8+ (no external libraries).

## Quick start

```bash
cd SHUTTLE/tools

# Regenerate the C frequency tables (one-shot, deterministic).
python3 gen_rans_tables.py --out ../ref/rans_tables.h

# Inspect the generated tables in a human-readable format.
python3 gen_rans_tables.py --format text

# Run the Python reference self-test.
python3 rans.py
```

## Where the frequency tables come from

We do *not* sample an empirical histogram. The constituent PMFs are
evaluated analytically: response symbols use the exact rounded discrete-
Gaussian marginal, and hints use the bucket-crossing distribution of
`MakeHint`. The shared rANS table is their coefficient-count-weighted
mixture.

`M_voc` is derived from the configured `11*r` design bound. The C encoder
also checks every symbol and returns an error outside the table alphabet;
the signer treats that result as a restart condition.
All three modes use `alpha_h = 2*alpha_r`, so the two blocks have comparable
scales and can be encoded efficiently by one shared table; their PMFs are
not assumed to be identical.

Quantization is the standard "round each `p(s) * 2^t` to an integer,
nudge the largest bucket to absorb the slack" recipe, with `t = 10`.

## Reservation budgets

The merged-stream byte reservations in `ref/params.h` come from
`tools/SigSize.py::compute_merged_reservation`. Re-run it when changing a
parameter or regenerating the frequency table; the resulting value feeds
`SHUTTLE_RANS_RESERVED_BYTES`. The implementation target is `2^-30`.

## rANS reference

`rans.py` is a literal transliteration of the C engine. State is 32-bit,
renormalization is byte-wise, the final-state flush is little-endian
(matching ryg_rans / HAETAE), and the decoder verifies `x == RANS_L`
after the last symbol and requires exact input consumption. Signature
unpacking additionally checks zero padding and re-encodes the symbols to
enforce a unique byte representation.

```python
from rans import RansTable, encode, decode
table = RansTable(syms=[...], freqs=[...], prob_bits=10)
stream = encode([0, -1, 2, 0, 1, ...], table)
recovered = decode(stream, table, len(msg))   # raises on final-state mismatch
```
