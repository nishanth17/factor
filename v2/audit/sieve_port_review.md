# C sieve ideas for Factor's Python 3 implementation

Reviewed 3 October 2026 against sibling `primesieve` commit
`5b4afb8f344ad5f6fbd20a186cd57b36182ea710`; the working tree was clean.
The C repository remains unchanged. This review does not reproduce its published
benchmarks or assume their speedups survive the Python runtime.

| Reference | Transfer candidate | Phase and decision |
| --- | --- | --- |
| C v2 `mark_erat_segment` | Odd-only candidates, private marking buffers, carry next strike across contiguous segments | Exact/private-state ideas used in Phase 1; rolling offsets and context reuse evaluated in P2.6 |
| C v2 Atkin path | Compact residue rows, signed recurrence seeds, shortened final segment | Preserve mod-60 recurrences with exact integer arithmetic and private rows; shorten final segment in Phase 1 |
| C v3 wheel-30 path | Eight residues per 30-number block; bounded ordered materialization | P2.7 experiment against bytearray slicing; Python per-strike bit updates may cost more than the memory saved |
| C v4 `sieve_context_create` / header | Reuse base primes, metadata, workspace; bound maximum endpoint; one owner per context | P2.6; add context limits and arbitrary restart tests before reuse |
| C v4 pre-sieve/strength-reduced wheel | Precomputed metadata, product-pattern initialization, phase-based strike updates | P2.7; compare pattern slicing/packed bitsets with Python bookkeeping cost |
| C v4 extraction and count/material policies | Separate marking from output; avoid materializing when only counts are needed | P2.7 diagnostics; factoring needs prime values, so include generation, decoding, traversal, and output memory |
| C v4 width-gated wheel-210 | Skip already-presieved sparse strikes without changing wheel-30 storage | Later P6.2 only if sparse-strike profiling identifies a bottleneck; C prime/width gates are not portable defaults |
| C v2/v3/v4 OpenMP scheduling | Independent segment ownership, read-only shared base primes, bounded ordered output | P2.8/P3.6/P6.3 evaluate Python processes and threads; native C thread scaling is not Python evidence |
| C v4 rejected bucket rewrite | Persistent sparse-prime scheduler changes architecture | Later optional P6.2; rejected there for scope, not proof that buckets never help Factor |
| C v4 compiler/SIMD/assembly/QoS choices | Native kernel and hardware tuning | Do not mechanically port to Python; revisit only if a measured native arithmetic boundary is introduced |

The C single-bound APIs use `p < n`, but their segmented API is inclusive
`[lo, hi]`. Factor v2 consistently uses `lo, hi)`: a C comparison must request
`segmented_sieve(lo, hi-1)`. Use exact `isqrt(hi-1)` and include the resulting
base-prime endpoint. Counts alone cannot prove prime-sequence correctness.

Relevant inspected sources:

- [C v2 sieve.c: rolling odd-index marking.
- C v3 sieve.c: wheel-30 and ordered bounded output.
- C v4 sieve.c and sieve.h: reusable contexts, wheel state, output ownership.
- C v4 README: retained/rejected experiments and the limits of its native benchmarks.

Parallelism needs two distinct measurements: matched candidate throughput and
first-valid-factor latency with cancellation. Start with independent ECM curves
or SIQS polynomial families rather than trying to parallelize individual Python
multiplications. Record the interpreter's actual GIL mode, use deterministic job
IDs/seeds, include cold spawn and warm worker reuse, and enforce aggregate
memory/CPU budgets. No default worker count is implied by this review.

## M9 follow-through: basic wheel-6 control

M9 reuses the original Python wheel-6 idea with bytearray slice marking, exact
`isqrt`, and a slot count that emits precisely `p < hi`. It reduces candidate
storage from roughly `hi/2` to `hi/3`, skipping multiples of both 2 and 3.
Every tiny endpoint residue and square boundary is checked against independent
trial division; larger lists are compared with an independent sieve.

This is a Python control improvement, not a mechanical port of the C wheel-30
or wheel-210 kernels. M9 measurements include full
prime materialization and complete factorization. P2.6 reusable contexts and
rolling strikes, and P2.7 pre-sieve/packed extraction remain open.
