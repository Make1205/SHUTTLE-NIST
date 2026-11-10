# SHUTTLE — QUICK

QUICK is a workflow check only. Full-scale figures are measured locally.

CPU: AMD EPYC 9V74 80-Core Processor; affinity: 0. GCC: gcc (Debian 14.2.0-19) 14.2.0.

Tested commit: `34bf327eee9459a9b29d271a43c305f4b0b8b374`. Benchmark commit: `f76c032088c129b7139bb29f311a5c86fec5c9b1`.

Base flags: `-O3 -march=native -std=gnu11`. Each operation: 11 iterations; correctness and lengths: 100 signatures per build/set. Message: 32 bytes. Median units: thousand RDTSC ticks.

Timing uses fresh byte-key triples and discards 20 complete warmup triples.

Exact per-build options (in addition to base flags and local include paths):

* 120 ref: `-DSHUTTLE_NIST_MODE=120`; link `-lm`.
* 180 ref: `-DSHUTTLE_NIST_MODE=180`; link `-lm`.
* 255 ref: `-DSHUTTLE_NIST_MODE=255`; link `-lm`.

| 参数集 | n | q | pk | sig（仓库实际） | 论文 sig | ref KeyGen/Sign/Verify | AVX2 KeyGen/Sign/Verify |
|---|---:|---:|---:|---|---:|---|---|
| 120 | 256 | 32257 | 992 | 1175; unpadded avg 1144.6700, max 1157 | 1175 | 113.307/699.373/192.217 | -- |
| 180 | 256 | 15361 | 1376 | 1589; unpadded avg 1553.0700, max 1568 | 1589 | 193.387/1042.989/299.025 | -- |
| 255 | 256 | 7681 | 1696 | 2013; unpadded avg 1975.0800, max 1988 | 2013 | 291.277/1379.871/424.475 | -- |



Expanded-key Sign (thousand ticks; `--` means no supported expanded-key API):

| 参数集 | ref | AVX2 |
|---|---:|---:|
| 120 | -- | -- |
| 180 | -- | -- |
| 255 | -- | -- |

Correctness: all configured upstream KAT/self-tests passed; all generated signatures verified; changed messages rejected. Cross-build/equality details are in correctness.json.

Only the upstream portable reference exists. sampler_4x is scalar interleaving, not AVX2; AVX2 is unavailable (--). Contrary to the variable-length expectation, this revision always emits SHUTTLE_BYTES and verifies that exact size. The internal rANS payload is variable and padded; unpadded length reads its explicit two-byte little-endian length after seedC and excludes only unused rANS reservation bytes. KAT: original deterministic vector generator, ten rounds with original/tampered-signature/tampered-message checks. No expanded-key API. Timing uses native randombytes; fixed-seed checks use bench entropy. Timing update: fresh encoded keypair/signature per KeyGen/Sign/Verify triple, 20 complete warmup triples, 32-byte zero message. Ref has no bench-added auto-vectorization suppression. Algorithm source and native randomness are unchanged.

Build flags and commands per variant are in config.json, metadata.json and build logs. All source is extracted locally from the tested commit; algorithm files are unchanged. See bench/NOTES.md.
