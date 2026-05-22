# SHUTTLE-NIST: Reference Implementation

A portable C99 reference implementation of the **SHUTTLE-NIST** lattice-based signature scheme for three NIST classical-security levels: **120-bit (Level I)**, **180-bit (Level III)** and **255-bit (Level V)**.

This artifact is anonymized for review. It is a self-contained reference, intended for correctness and size validation, not for production deployment.


## 1. Algorithm Parameters

Symbols match the algorithm specification. All three modes use ring dimension `n = 256` and full NTT (`q ≡ 1 mod 2n`). All three modes share `alpha_h = 2 * alpha_r = 128`, which makes `sigma_hint = 2r/alpha_h = r/alpha_r = sigma_z-hi`; HighBits of `z` and `hint` therefore share a single rANS frequency table and are encoded in one merged stream.

| Symbol | Description | NIST-120 | NIST-180 | NIST-255 |
| --- | --- | ---: | ---: | ---: |
| `n` | Ring dimension | 256 | 256 | 256 |
| `q` | Modulus | 32257 | 15361 | 7681 |
| `lenS` | Polynomial components of `s` | 3 | 4 | 5 |
| `lenE` | Polynomial components of `e` | 2 | 3 | 4 |
| `eta` | CBD parameter of `s`, `e` | 1 | 1 | 1 |
| `tau` | Hamming weight of challenge `c` | 58 | 80 | 128 |
| `r` | Standard deviation of Gaussian sampling `y` | 67 | 79 | 90 |
| `alpha_1` | CompressY step on slot 0 | 8 | 8 | 16 |
| `alpha_h` | HighBits / MakeHint step | 128 | 128 | 128 |
| `B_k` | Secret-key rejection bound | 25.79 | 30.4 | 34.4 |
| `B_s` | Signature rejection bound | 3233 | 4312 | 5556 |
| `B_v` | Verification bound | 4000 | 5250 | 6638 |
| `\|seedC\|` | Challenge hash length (bytes) | 32 | 48 | 64 |
| | MLWE (Core-SVP, bits) | 121.7 | 185.4 | 258.1 |
| | MSIS (Core-SVP, bits) | 156.5 (124.6) | 225.4 (182.5) | 278.9 (227.8) |

rANS implementation parameters (single merged stream):

| Symbol | Description | NIST-120 | NIST-180 | NIST-255 |
| --- | --- | ---: | ---: | ---: |
| `alpha_r` | rANS hi/lo split step | 64 | 64 | 64 |
| `alpha_0'` | Slot-0 secondary hi/lo split, `= alpha_r / alpha_1` | 8 | 8 | 4 |
| `sigma_z-hi` | `= r / alpha_r` | 1.0469 | 1.2344 | 1.4062 |
| `sigma_hint` | `= 2r / alpha_h` (equals `sigma_z-hi`) | 1.0469 | 1.2344 | 1.4062 |
| `t` | rANS probability precision (sum of frequencies = `2^t`) | 10 | 10 | 10 |
| `L_ren` | rANS renormalization threshold | `2^23` | `2^23` | `2^23` |
| `M_voc` | Merged-stream vocabulary radius `[-M_voc, M_voc]` | 13 | 15 | 18 |
| `n_b` | Merged-stream symbol count `= n*(lenS+1+lenE)` | 1536 | 2048 | 2560 |

Sizes (measured from this implementation):

| | NIST-120 | NIST-180 | NIST-255 |
| --- | ---: | ---: | ---: |
| Public-key size (bytes) | 992 | 1376 | 1696 |
| Secret-key size (bytes) | 448 | 576 | 704 |
| **Signature size (bytes)** | **1153** | **1565** | **1965** |
| pk + sig (bytes) | 2145 | 2941 | 3661 |

Discrete-Gaussian sampler (BLISS-style): `r = 64 * sigma_2`, with `sigma_2 = 1.0469 / 1.2344 / 1.4062`. The reverse cumulative distribution table (`RCDT`) uses 3 limbs of 31 bits each (93-bit precision). Renyi divergence per entry:

| | order | divergence |
| --- | --- | --- |
| NIST-120 (`sigma_2 = 1.0469`, 11 entries) | 241 | `1 + 2^-95.55` |
| NIST-180 (`sigma_2 = 1.2344`, 13 entries) | 361 | `1 + 2^-93.76` |
| NIST-255 (`sigma_2 = 1.4062`, 15 entries) | 511 | `1 + 2^-93.88` |

Iterated Rejection Sampling (IRS): the IntervalParity walk uses 37 iterations per challenge monomial.


## 2. Repository Layout

```
SHUTTLE-NIST/
├── ref/                     # C reference implementation
│   ├── config.h             # mode selector (SHUTTLE_NIST_MODE = 120 / 180 / 255)
│   ├── params.h             # all per-mode parameters and derived constants
│   ├── sign.{c,h}           # high-level keypair / sign / verify
│   ├── packing.{c,h}        # pk / sk / signature serialization
│   ├── shuttle_rans.{c,h}   # rANS encoder / decoder (single merged stream)
│   ├── rans_tables.h        # auto-generated rANS frequency tables
│   ├── sampler.{c,h}        # BLISS-style discrete Gaussian sampler
│   ├── rejsample.{c,h}      # iterative rejection sampling (IRS)
│   ├── approx_exp.{c,h}     # constant-time fixed-point exp(-a)
│   ├── approx_log.{c,h}     # constant-time fixed-point -ln(r)
│   ├── poly.{c,h}           # single-polynomial operations
│   ├── polyvec.{c,h}        # polynomial-vector operations
│   ├── ntt.{c,h}            # NTT / inverse NTT
│   ├── reduce.{c,h}         # modular reduction helpers
│   ├── rounding.{c,h}       # HighBits / LowBits / MakeHint / UseHint
│   ├── fips202.{c,h}        # SHAKE-256 (FIPS 202)
│   ├── symmetric{.h,-shake.c}
│   ├── randombytes.{c,h}    # OS RNG wrapper
│   ├── gen_zetas.py         # regenerates ntt_constants.h
│   ├── gen_approx_exp_tables.py  # regenerates approx_exp_constants.h
│   ├── Makefile             # build driver (MODES = 120 180 255)
│   ├── test/                # functional + benchmarking tests
│   └── vector/              # KAT-vector generation
└── tools/                   # offline analysis utilities (Python)
    ├── gen_rans_tables.py   # regenerates ref/rans_tables.h
    ├── SigSize.py           # Shannon lower bound + rANS reservation sizing
    ├── rans.py              # pure-Python rANS reference
    ├── test_rans_unified.py # rANS table / round-trip / corruption tests
    └── kat_compare/         # cross-platform KAT comparison harness
```

Each mode is selected by the compile-time macro `SHUTTLE_NIST_MODE` (120 / 180 / 255). The Makefile builds all three side-by-side; symbols are namespace-prefixed as `shuttle_nist{120,180,255}_ref_*` so the three shared libraries can coexist.


## 3. Build Requirements

* A C99-capable compiler with support for `unsigned __int128` (tested with `gcc` and `clang` on x86-64 Linux).
* GNU `make`.
* Python 3.8+ (only needed if you want to regenerate the auto-generated headers `ntt_constants.h`, `approx_exp_constants.h`, `rans_tables.h`; the shipped headers are bit-exact and ready to compile).

No external C dependencies; everything required is in this repository.


## 4. Build

From `SHUTTLE-NIST/ref/`:

```bash
make MODES="120 180 255"        # builds all three modes
# or build a single mode:
make MODES="120"
```

This produces, under `ref/test/out/`:

* `libshuttle_nist{120,180,255}_ref.so` — shared libraries
* per-mode test binaries: `test_sampler<m>`, `test_ntt<m>`, `test_mod2q<m>`, `test_rans<m>`, `test_poly_z1<m>`, `test_sig<m>`, `test_sign<m>`, `speed_sampler<m>`, `speed_sign<m>`, `profile_sign<m>`, `diag_sign<m>`

Clean with `make clean`.


## 5. Correctness Testing

### 5.1 Full functional regression

From `SHUTTLE-NIST/ref/`:

```bash
make check MODES="120 180 255"
```

This runs, for every mode, the following test binaries in sequence and checks each for a `PASS` / `All * PASSED` line:

* `test_sampler<m>` — chi-squared / moment checks on the discrete Gaussian sampler against the target sigma.
* `test_ntt<m>` — forward / inverse NTT round-trip and convolution against a schoolbook reference.
* `test_mod2q<m>` — modular arithmetic, MakeHint / UseHint round-trip, `z_2` recovery.
* `test_rans<m>` — rANS encode / decode round-trip on small and large inputs, empty-stream little-endian flush, and single-byte corruption detection.
* `test_poly_z1<m>` — z-split / z-combine round-trip and full pack pipeline.
* `test_sig<m>` — `pack_sig` / `unpack_sig` round-trip with random inputs drawn from the theoretical distributions; reports observed rANS lengths against the static reservation.
* `test_sign<m>` — end-to-end keypair / sign / verify, a forgery detection test (signature and message tampering), and a 100-round stress test with random messages.

A successful run prints `All * tests PASSED` for each binary and ends with `All mode check runs completed`.

### 5.2 KAT (Known-Answer Test) generation

```bash
cd SHUTTLE-NIST/ref/vector
make            # builds the per-mode gen_vectors binaries
# then run, e.g.:
./gen_vectors_120 > shuttle_nist120.kat
```

`gen_vectors_<m>` uses a deterministic NIST-style DRBG (`rng.c`) so that the output is byte-identical across platforms; this is the intended cross-implementation check vector.

### 5.3 Python reference cross-check

```bash
cd SHUTTLE-NIST/tools
python3 test_rans_unified.py
```

This validates the auto-generated rANS frequency tables (alphabet geometry, frequency sums, vocabulary radii) and exercises the pure-Python rANS reference (`rans.py`) on a discrete-Gaussian sample.


## 6. Performance Testing

### 6.1 Signature throughput

From `SHUTTLE-NIST/ref/`:

```bash
./test/out/speed_sign120
./test/out/speed_sign180
./test/out/speed_sign255
```

Each binary measures cycles-per-call for `crypto_sign_keypair`, `crypto_sign_signature` and `crypto_sign_verify` using the platform's `rdtsc`-equivalent cycle counter (`test/cpucycles.c`). On x86-64 we pin to a single core and disable Turbo / SMT for reproducibility:

```bash
# Recommended:
sudo cpupower frequency-set --governor performance
taskset -c 0 ./test/out/speed_sign120
```

`speed_sign<m>` reports min, median and average over a large number of iterations. Use the median as the headline number.

### 6.2 Sampler microbenchmark

```bash
./test/out/speed_sampler<m>
```

Reports cycles per sample for the BLISS-style discrete Gaussian sampler in isolation.

### 6.3 Per-stage profiling

```bash
./test/out/profile_sign<m>
```

Breaks `crypto_sign_signature` down into:

1. preprocessing (key unpack, mu / rhoprime hashing);
2. `y` sampling;
3. commitment;
4. challenge;
5. iterative rejection sampling (IRS);
6. finalize (post-IRS arithmetic);
7. pack plumbing (header + low-bits packing);
8. merged-stream rANS encode.

Useful to identify hotspots before optimizing.

### 6.4 Signature-size reproducibility

```bash
cd SHUTTLE-NIST/tools
python3 SigSize.py
```

Prints the Shannon lower bound, the per-component byte breakdown, and the rANS reservation `R` at three target overflow probabilities (`2^-10`, `2^-20`, `2^-40`). The implementation uses `p_rans^* = 2^-20`.


## 7. Reproducibility Notes

* `randombytes.c` uses `getrandom(2)` on Linux. For deterministic fixed-seed runs, `vector/rng.c` provides a stand-alone NIST-style AES-CTR DRBG implementation; `vector/gen_vectors.c` uses it.
* All auto-generated headers are regenerable bit-for-bit from the Python scripts in `ref/` and `tools/`. Run:

```bash
python3 ref/gen_zetas.py
python3 ref/gen_approx_exp_tables.py
python3 tools/gen_rans_tables.py --out ref/rans_tables.h
```

These three scripts together regenerate every auto-generated header from first-principles math (no empirical samples).
* The IRS routine, the discrete Gaussian sampler and the inner approximations (`approx_exp`, `approx_log`) are constant-time with respect to all secret-dependent data; there are no secret-dependent branches and no secret-dependent table indexing.
* `__int128` is used in `approx_exp.c`, `approx_log.c` and `sampler.c` for 128-bit fixed-point arithmetic. Both `gcc` and `clang` lower this to a single `mulq` on x86-64.


## 8. Benchmark Results (to be filled in)

### 8.1 Platform

* CPU: _TBD_
* Microarchitecture: _TBD_
* Frequency policy: _TBD_ (e.g. fixed 2.x GHz, Turbo / SMT disabled)
* Compiler: _TBD_
* Compiler flags: _TBD_

### 8.2 Per-mode cycle counts

| Operation | NIST-120 | NIST-180 | NIST-255 |
| --- | ---: | ---: | ---: |
| `crypto_sign_keypair` (median) | _TBD_ | _TBD_ | _TBD_ |
| `crypto_sign_signature` (median) | _TBD_ | _TBD_ | _TBD_ |
| `crypto_sign_verify` (median) | _TBD_ | _TBD_ | _TBD_ |
| Sampler (cycles / sample) | _TBD_ | _TBD_ | _TBD_ |

### 8.3 Per-stage signing breakdown

| Stage | NIST-120 | NIST-180 | NIST-255 |
| --- | ---: | ---: | ---: |
| Preprocess | _TBD_ | _TBD_ | _TBD_ |
| Gaussian `y` sampling | _TBD_ | _TBD_ | _TBD_ |
| Commitment | _TBD_ | _TBD_ | _TBD_ |
| Challenge | _TBD_ | _TBD_ | _TBD_ |
| IRS | _TBD_ | _TBD_ | _TBD_ |
| Finalize | _TBD_ | _TBD_ | _TBD_ |
| Pack plumbing | _TBD_ | _TBD_ | _TBD_ |
| rANS encode | _TBD_ | _TBD_ | _TBD_ |


## 9. Comparison with Related Schemes (to be filled in)

### 9.1 Size

| Scheme | Level | pk (B) | sk (B) | sig (B) | pk + sig (B) |
| --- | --- | ---: | ---: | ---: | ---: |
| **SHUTTLE-NIST** | I (120) | 992 | 448 | 1153 | 2145 |
| **SHUTTLE-NIST** | III (180) | 1376 | 576 | 1565 | 2941 |
| **SHUTTLE-NIST** | V (255) | 1696 | 704 | 1965 | 3661 |
| _comparator A_ | I | _TBD_ | _TBD_ | _TBD_ | _TBD_ |
| _comparator A_ | III | _TBD_ | _TBD_ | _TBD_ | _TBD_ |
| _comparator A_ | V | _TBD_ | _TBD_ | _TBD_ | _TBD_ |
| _comparator B_ | I | _TBD_ | _TBD_ | _TBD_ | _TBD_ |
| _comparator B_ | III | _TBD_ | _TBD_ | _TBD_ | _TBD_ |
| _comparator B_ | V | _TBD_ | _TBD_ | _TBD_ | _TBD_ |

### 9.2 Performance

| Scheme | Level | keypair (cycles) | sign (cycles) | verify (cycles) |
| --- | --- | ---: | ---: | ---: |
| **SHUTTLE-NIST** | I (120) | _TBD_ | _TBD_ | _TBD_ |
| **SHUTTLE-NIST** | III (180) | _TBD_ | _TBD_ | _TBD_ |
| **SHUTTLE-NIST** | V (255) | _TBD_ | _TBD_ | _TBD_ |
| _comparator A_ | I | _TBD_ | _TBD_ | _TBD_ |
| _comparator A_ | III | _TBD_ | _TBD_ | _TBD_ |
| _comparator A_ | V | _TBD_ | _TBD_ | _TBD_ |
| _comparator B_ | I | _TBD_ | _TBD_ | _TBD_ |
| _comparator B_ | III | _TBD_ | _TBD_ | _TBD_ |
| _comparator B_ | V | _TBD_ | _TBD_ | _TBD_ |


## 10. License

This artifact is released for academic review purposes. License terms will be specified in the camera-ready version.
