# Timing protocol

System GCC -O3 -march=native plus required backend options; no extra bench
auto-vectorization suppression. Message 32 zero bytes. Every timing round
generates fresh encoded keys, signs from byte sk, verifies from byte pk. Discard
20 complete warmup triples. Native decoding, expansion, hashing and RNG remain
inside the public APIs. Use unchanged Dilithium bba9534 RDTSC/median files.
QUICK=1: 11 timings, 100 correctness/length signatures; FULL: 10001/10000.
Only QUICK runs in cloud. Native KAT/equality/cross checks remain mandatory.
