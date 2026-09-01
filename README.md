# SHUTTLE-NIST: Reference Implementation

A portable C99 reference implementation of the **SHUTTLE-NIST** lattice-based signature scheme for three NIST classical-security levels: **120-bit (Level I)**, **180-bit (Level III)** and **255-bit (Level V)**.

This artifact is anonymized for review. It is a self-contained reference, intended for correctness and size validation, not for production deployment.


## 1. Algorithm Parameters

Symbols match the algorithm specification. All three modes use ring dimension `n = 256` and full NTT (`q ≡ 1 mod 2n`). HighBits of `z` and `hint` are encoded in one merged rANS stream. The shared frequency table is generated from the coefficient-count-weighted mixture of the rounded-response and bucket-crossing PMFs.

| Symbol | Description | NIST-120 | NIST-180 | NIST-255 |
| --- | --- | ---: | ---: | ---: |
| `n` | Ring dimension | 256 | 256 | 256 |
| `q` | Modulus | 32257 | 15361 | 7681 |
| `lenS` | Polynomial components of `s` | 3 | 4 | 5 |
| `lenE` | Polynomial components of `e` | 2 | 3 | 4 |
| `eta` | CBD parameter of `s`, `e` | 1 | 1 | 1 |
| `tau` | Hamming weight of challenge `c` | 58 | 80 | 128 |
| `r` | Standard deviation of Gaussian sampling `y` | 70 | 82 | 99 |
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
| `sigma_z-hi` | `= r / alpha_r` | 1.0938 | 1.2812 | 1.5469 |
| `sigma_hint` | `= 2r / alpha_h` (equals `sigma_z-hi`) | 1.0938 | 1.2812 | 1.5469 |
| `t` | rANS probability precision (sum of frequencies = `2^t`) | 10 | 10 | 10 |
| `L_ren` | rANS renormalization threshold | `2^23` | `2^23` | `2^23` |
| `M_voc` | Merged-stream vocabulary radius `[-M_voc, M_voc]` | 13 | 16 | 20 |
| `n_b` | Merged-stream symbol count `= n*(lenS+1+lenE)` | 1536 | 2048 | 2560 |

Sizes (measured from this implementation):

| | NIST-120 | NIST-180 | NIST-255 |
| --- | ---: | ---: | ---: |
| Public-key size (bytes) | 992 | 1376 | 1696 |
| Secret-key size (bytes) | 448 | 576 | 704 |
| **Signature size (bytes)** | **1175** | **1589** | **2013** |
| pk + sig (bytes) | 2167 | 2965 | 3709 |

The signature layout is `seedC || uint16(rans_len) || rANS payload/padding || low bits`. IRS update signs are not serialized: verification regenerates the challenge from `seedC` and uses its parity. The rANS payload must consume its declared length exactly, its fixed-slot padding must be zero, and decoding must reproduce the same canonical byte stream.

The BLISS-style discrete-Gaussian sampler uses `r = 64 * sigma_2`, with the exact values `sigma_2 = 70/64, 82/64, 99/64`. Its reverse cumulative distribution tables use 3 limbs of 31 bits (93-bit precision) and 12 / 14 / 17 entries, respectively. `ref/gen_rcdt_tables.py` deterministically regenerates the tables from the specified nonnegative discrete-Gaussian tails.

Direct finite-distribution evaluation, including the Q63 acceptance threshold and all `n*(1+lenS+lenE)` mask coefficients, gives `M1 - 1 < 2^-86.78, 2^-85.45, 2^-87.83` at Rényi orders 240, 360, and 510, respectively.

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
│   ├── rcdt_tables.h         # auto-generated 93-bit reverse-CDT tables
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
│   ├── gen_rcdt_tables.py    # regenerates rcdt_tables.h
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
* Python 3.8+ (only needed to regenerate the auto-generated headers; the shipped headers are ready to compile).

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
* `test_sig<m>` — `pack_sig` / `unpack_sig` round-trip plus rejection of non-canonical rANS lengths and padding.
* `test_sign<m>` — end-to-end keypair / sign / verify, forgery and canonical-encoding rejection tests, legacy-length rejection, and a 100-round stress test.

A successful run prints `All * tests PASSED` for each binary and ends with `All mode check runs completed`.

### 5.2 KAT (Known-Answer Test) generation

```bash
cd SHUTTLE-NIST/ref/vector
make NROUNDS=10  # builds generators and writes KAT files for all three modes
```

`gen_vectors_mode<m>` uses the deterministic SHAKE-256 DRBG in `rng.c`, so output is byte-identical across platforms.

### 5.3 Python reference cross-check

```bash
cd SHUTTLE-NIST/tools
python3 test_rans_unified.py
```

This validates the auto-generated rANS frequency tables (alphabet geometry, frequency sums, vocabulary radii) and exercises the pure-Python rANS reference (`rans.py`) on samples from the analytic merged PMF.


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

### 6.3 Signature-size reproducibility

```bash
cd SHUTTLE-NIST/tools
python3 SigSize.py
```

Prints the Shannon lower bound, the per-component byte breakdown, and the rANS reservation `R` at three target overflow probabilities (`2^-20`, `2^-30`, `2^-40`). The implementation uses `p_rans^* = 2^-30`.


## 7. Reproducibility Notes

* `randombytes.c` uses `getrandom(2)` on Linux. Deterministic vector generation substitutes the SHAKE-256 DRBG in `vector/rng.c`.
* All auto-generated headers are regenerable bit-for-bit from the Python scripts in `ref/` and `tools/`. Run:

```bash
python3 ref/gen_zetas.py
python3 ref/gen_approx_exp_tables.py
python3 ref/gen_rcdt_tables.py
python3 tools/gen_rans_tables.py --out ref/rans_tables.h
```

These four scripts regenerate every auto-generated header from first-principles math (no empirical samples).
* The IRS routine, the discrete Gaussian sampler and the inner approximations (`approx_exp`, `approx_log`) are constant-time with respect to all secret-dependent data; there are no secret-dependent branches and no secret-dependent table indexing.
* `__int128` is used in `approx_exp.c`, `approx_log.c` and `sampler.c` for 128-bit fixed-point arithmetic. Both `gcc` and `clang` lower this to a single `mulq` on x86-64.


## 8. Benchmark Results

> All numbers below are from **portable C reference implementations** (no AVX2 / SSE / NEON vectorization), measured by `rdtsc`-class cycle counter. Raw logs are in [`ref/speed_SHUTTLE.txt`](ref/speed_SHUTTLE.txt), [`ref/speed_HAETAE.txt`](ref/speed_HAETAE.txt), [`ref/speed_dilithium.txt`](ref/speed_dilithium.txt).

### 8.1 Platform

* OS: Ubuntu 24.04.2 LTS (kernel 6.11.0-26-generic, x86_64)
* CPU: 11th Gen Intel(R) Core(TM) i7-11700K @ 3.60 GHz (fixed)
* Turbo Boost: **disabled** (`intel_pstate/no_turbo = 1`)
* Hyper-Threading: **disabled** (one thread per core)
* Compiler: gcc 13.3.0 (Ubuntu 13.3.0-6ubuntu2~24.04)
* Build: each scheme's default Makefile (`CC ?= gcc`); no `-mavx2`, no hand-written assembly. Reference C only.

### 8.2 Per-mode cycle counts

Cycle counts have not been refreshed for the current parameter constants. Run `make speed-sign` on the platform above before publishing new figures; archived logs are not treated as current measurements.


## 9. Comparison with Related Schemes

> All schemes compared in this section are **reference (portable C) implementations without AVX2 / vector optimization**, measured on the platform of §8.1.

### 9.1 Size

| Scheme | Level | pk (B) | sig (B) | pk + sig (B) |
| --- | --- | ---: | ---: | ---: |
| **SHUTTLE-NIST** | I (120) | 992 | 1175 | **2 167** |
| HAETAE | I (120) | 992 | 1474 | 2 466 |
| Dilithium / ML-DSA | 2 | 1312 | 2420 | 3 732 |
| **SHUTTLE-NIST** | III (180) | 1376 | 1589 | **2 965** |
| HAETAE | III (180) | 1472 | 2349 | 3 821 |
| Dilithium / ML-DSA | 3 | 1952 | 3309 | 5 261 |
| **SHUTTLE-NIST** | V (255) | 1696 | 2013 | **3 709** |
| HAETAE | V (260) | 2080 | 2948 | 5 028 |
| Dilithium / ML-DSA | 5 | 2592 | 4627 | 7 219 |

### 9.2 Performance

No current SHUTTLE-NIST cycle comparison is reported until the revised constants are benchmarked on the platform in Section 8.1.


## 10. License

This artifact is released for academic review purposes. License terms will be specified in the camera-ready version.
