# Benchmark notes

Initial worktree: clean. `bench` was created from the recorded initial commit.
Only files under `bench/` are added. Algorithm source, parameters and original
build files remain unchanged. The runner extracts the selected commit with local
`git archive`, so an older upstream revision can be measured without resetting
the working branch or rewriting files. No network is used at run time.

Initial commit: `34bf327eee9459a9b29d271a43c305f4b0b8b374`.
Tested commit: `34bf327eee9459a9b29d271a43c305f4b0b8b374`.

## Timing source and measurement

The four unmodified timing files come from pq-crystals/dilithium
`bba9534a3c5cbef66481cd6518d786e23ae601b5`, `ref/test/` (there is no root
`test/` directory at that revision). They match the Antrag bench copies byte
for byte. RDTSC is used; USE_RDPMC is not defined. The unchanged print_results
function receives iterations+1 timestamps for exactly iterations operations,
subtracts its own counter overhead, and computes the median. Thousand cycles
means RDTSC ticks divided by 1000. Loop/call overhead is included consistently.

System GCC: -O3 -march=native -std=gnu11. Ref additionally disables automatic
loop/SLP vectorization; algorithm SIMD intrinsics are selected only in the
upstream SIMD implementation. Exact flags are recorded per variant in
config.json, metadata.json and build logs. Messages have 32 bytes. Signing and
verification receive serialized byte keys and include all decoding/expansion
performed by the upstream byte API. KeyGen returns serialized keys.
Expanded-key Sign is reported only when a supported API is available.

QUICK=1: 11 timing operations, upstream KAT/self-tests, 100 generated signatures.
Default FULL: 10001 timing operations, upstream KAT/self-tests, 10000 generated
signatures. Only QUICK is executed in cloud. Full-scale runs are for the user's
local machine. No test assertions are disabled. The correctness and length
statistics reuse the same signature pool, using a fixed-seed key per build/set.
All signatures self-verify, cross-verify when another build exists, and are
checked for rejection after changing the message. Equality expectations are
configured explicitly; floating-point implementations can use cross verification
when byte equality is inappropriate. No result is declared complete unless all
required commands finish with zero status. The CSV contains medians and sizes;
raw timer output and correctness evidence are retained.

A bench-only xorshift64* byte stream (seed 1234567, fixed update and byte order in
driver.c) supplies reproducible test entropy where the API accepts an external
provider. This is for correctness vectors, not production randomness. For NIST
APIs the timing provider is the original randombytes implementation, renamed only
in its provider translation unit and routed through a bench wrapper; the algorithm
call sites remain unchanged. Callback APIs use the repository SHAKE generator
with fixed test seeds and OS-seeded timing state as documented below.

Missing system-library headers/libraries are compile/link checked at startup,
before algorithm builds. The runner never installs packages. Bash, Python 3.8+,
Git and system GCC are required; additional apt dependencies are listed below.

## Scheme-specific configuration and format

Only the upstream portable reference exists. sampler_4x is scalar interleaving, not AVX2; AVX2 is unavailable (--). Contrary to the variable-length expectation, this revision always emits SHUTTLE_BYTES and verifies that exact size. The internal rANS payload is variable and padded; unpadded length reads its explicit two-byte little-endian length after seedC and excludes only unused rANS reservation bytes. KAT: original deterministic vector generator, ten rounds with original/tampered-signature/tampered-message checks. No expanded-key API. Timing uses native randombytes; fixed-seed checks use bench entropy.

## Native-source compilation changes

No algorithm source patches. All source copies come from the tested commit.
Additional source selection, preprocessor options and include paths appear in
config.json. Configuration-only exceptions, if any, are documented below.

## Timing file checksums

* cpucycles.c: `31cc4c16b1e896b2cf7cc154aa21430329c0047cf9898a1865f9bfe9904f25e4`
* cpucycles.h: `7f1e2b68cf3e3c76076364af3c5350de61442450c477aec809da8de10e131838`
* speed_print.c: `533fc7bd12204078a5390cfb58c05d550b214850876609f8e98f7377b9d4f1b3`
* speed_print.h: `3946b0226c63125664eb443a20d43563e6d1ba7fa3be5c1b952cd29d1477b2ef`

Additional apt dependencies: none beyond build-essential, python3 and git.

No precise upstream commit was supplied for SHUTTLE. The current fork commit
is measured as-is. README names Ji-Peng/SHUTTLE-NIST as the public artifact.
There is no avx2 directory, no intrinsic implementation, and README explicitly
states that no vectorized implementation is included. Comments describing
future AVX2-compatible sample layouts do not establish an AVX2 backend.

The premise that the public signature is variable-length does not hold at this
revision: sign.c writes *siglen=SHUTTLE_BYTES and verify rejects other sizes.
packing.c explicitly defines a fixed-size rANS slot with a uint16 length and
zero padding before packed polynomial tails. The bench reads this exact length;
it does not strip the end of the signature (which contains meaningful tails).
The hypothetical unpadded size is a statistic, not a changed accepted format.

The original vector/gen_vectors.c and vector/rng.c generate ten deterministic
keygen/sign/verify rounds, including altered payload and message rejection.
They execute unchanged, in the extracted source directory. Additional bench
checks cover the configured 100 or 10000 signatures. There are no system
library requirements beyond libm and the standard compiler/runtime.
