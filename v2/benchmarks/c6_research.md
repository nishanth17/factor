# C6 / P4.5: bounded precomputed Lucas chains

The user activated C6 on 9 October 2026. The control is mainline
`bcf5f3d1e57304694b48ba6e7ef8b4ea2ffd0db0`; production remains unchanged.
The accepted A4 campaigns already reject checked PRAC promotion. This study
asks whether precomputation and compact execution change that decision.

## Sources and mathematical scope

| Source inspected | What it establishes | C6 use |
| --- | --- | --- |
| [GMP-ECM LucasChainGenerator, pinned 8ea5e214](https://github.com/sethtroisi/gmp-ecm/tree/8ea5e214fdf2f0ddf9415141b8dc039ed6f5874e/LucasChainGenerator) | McLaughlin/Zimmermann generator produces prime-chain codes described as optimal or near-optimal. Search enumerates increasing lengths with divisibility and reachability pruning; ties favor more doublings, then shorter codes. | Reproduce the upstream search at B1=2,000; preserve code and decoder inputs, independently verify all results. C6 certifies scalar action, not global minimality. |
| [GMP-ECM ecm.c at the same pin](https://github.com/sethtroisi/gmp-ecm/blob/8ea5e214fdf2f0ddf9415141b8dc039ed6f5874e/ecm.c) | A 64-bit code is decoded to offsets; stage one retains 16 projective points in a rolling buffer. The PRAC comments document a false-infinity example. | Extract decoder verbatim with assertions; compare the ring and last-use allocation on identical chains. Keep checks and recovery in the Python executor. |
| [CADO bytecode, pinned 692ecb7e](https://github.com/cado-nfs/cado-nfs/blob/692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b/sieve/ecm/bytecode.c) | PRAC rule bytecode separates generation, weighted costs, compression and integer checking. Its cache keys include the cost-model pointer. | Inspiration for separating chain choice from dispatch. No CADO code is copied. Our four-byte register instructions deliberately use one executor for PRAC and Lucas. |
| [Bernstein, Cottaar, Lange, *Searching for differential addition chains*](https://doi.org/10.1007/s40993-024-00604-8) | *Research in Number Theory* 11, article 45, published 27 March 2025; ePrint 2024/1044 is the earlier preprint. Minimum length means the continued-fraction subclass. Improved pruning and the left-length meet-in-the-middle variant are exact there; the left-interval variant can miss minima. Runtime is distinct from weighted M/S/constant costs. | Recognize existing CF records and evaluate Algorithm 1 with three live points. New offline search is conditional on the frozen compact-executor gate. No global-optimum claim. |
| [Kruppa, 2010 thesis](https://docnum.univ-lorraine.fr/public/SCD_T_2010_0054_KRUPPA.pdf), §4.5.1 | Precomputing PRAC rule sequences amortizes repeated stage-one schedules; interpreter compression is an implementation decision. | Explains why both fresh construction and reuse need measurement. |
| [GMP-ECM 6.2 maintainer announcement](https://gmplib.org/list-archives/gmp-discuss/2008-May/003158.html) | Technical discussion records a Lucas-chain generation bug affecting P+1/ECM. | Historical reminder to validate actions, not merely count instructions. |
| [Kushagra Singh's FLINT implementation blog](https://iamkush.me/implementing-the-ecm-stage-i-in-flint/) | Search located the author's account of ECM stage one and PRAC. Full page retrieval failed in this environment. | Discovery only; no correctness or PyPy performance claim relies on it. |

The paper's compressed CF state is `(a,b,c)`, initially `(1,2,3)` with
`c=a+b`. Each bit keeps either `a` or `b`, retains `c`, and computes its
sum with the retained entry using the other entry as known difference.
This allows three persistent points. General Lucas chains need not have
that shape; the GMP-ECM ring supports older differences and extra doublings.
PRAC can also use subtraction, so its record verifier checks both permitted
sum/difference directions. ECM x-coordinate chains are not a drop-in
implementation of Williams p+1 Lucas-sequence arithmetic.

## License and reproducibility

The exact GMP-ECM generator, decoder source and notices are versioned in
[inputs/upstream/c6_gmp_ecm](inputs/upstream/c6_gmp_ecm/NOTICE.md), under
LGPL-3.0-or-later with COPYING and COPYING.LIB. C6 invokes separate executables;
production does not link them. The only generator build change reduces the
six-million-entry allocation capacity to 4,096; search logic is unchanged.
One generator thread is used. CADO's LGPL-2.1 COPYING and source were inspected,
but no code was adapted. The Python allocator/verifier is independently written.

No licenses are inferred from a repository's popularity. This repository has
no top-level license granting blanket rights over all its files; upstream
notices apply to the vendored inputs and their extracted decoder.

## Independently verified contracts

* Every upstream offset record is decoded by upstream C, then translated to
  SSA, checked against its claimed integer values, and verified by A4's integer
  interpreter. A separate mutable-register interpreter checks the compact
  record, including overwritten and uninitialized slots and final scalar.
* Tests compare all 303 prime records with full-coordinate affine arithmetic,
  exercise small nonsingular curves, prime powers, CRT composite moduli,
  prime squares, proper nonunits, saturation, and the published regression
  `n=33554520197234177`, `sigma=2046841451`, B1=373. Projective validity is
  checked separately from equality; `(0,0)` cannot pass vacuously.
* Compact execution checks X and Z before discarding a point. A4 retains
  earlier X coordinates for exceptional recovery; a compact ring cannot
  discard those factors. The extra X GCD is charged, including during the
  one permitted checked-ladder retry. Infinity/order-two differences trigger
  recovery before ambiguous differential addition.
* A saturated prime power replays only from that power's starting point,
  at most its exponent times. A proper divisor is returned immediately;
  unresolved inputs remain explicit. B3's stage jobs/checkpoints are untouched.

## Frozen scope and finite allowances

`inputs/controls/c6_protocol.json` pins source hashes, certified inputs,
selection rules, bounds and seeds before accepted timing. Required inputs
are versioned; raw captures, failed runs and process transcripts stay in
`results/c6/` (unpacked upstream scratch is in local `audit/results/c6/`).

The scope is B1=2,000 / B2=147,396, eight curves per attempt and twenty seconds
wall/CPU per attempt, on ten existing certified inputs: balanced and 10-digit
small-factor cases at each input size 40/50/60/70/80 digits. Balanced target
factors are 20/25/30/35/40 digits. This is a finite feasibility comparison,
not adequate coverage for balanced 80-digit factoring. Seeds 41001/48920
are fixed before timing; the next two corpus seeds are reserved for any
qualifying fresh confirmation. Each timing sample repeats the complete
fixed cohort, so repeated samples are not new independent factoring trials.

Prime generation is offline, one thread, B1=2,000, with sixty-second CPU/wall
limits, a 512 MiB process-RSS watchdog and 16 MiB per-output-file limit.
The decoder is an isolated process consuming only those generated codes.
It uses a 64-element upstream array and assertions; Python rejects records
above 64 elements. This is not a public decoder for arbitrary hostile codes.
The parent also imposes the subprocess timeout and validates every output.

Python records have 32-bit scalars, at most 512 instructions and 16 retained
register points. Execution additionally retains one original point for retry
and one result temporary. CF execution has three persistent working points,
plus original/recovery and arithmetic temporaries. A bound-owned program
has at most 512 unique records and 1 MiB of bytecode. Catalog reads are capped
at 1 MiB. Construction temporarily owns the <=303-entry Lucas SSA catalog;
no curve points are cached, and the A4 compiler cache is cleared after each
program build. No full-lcm scalar is searched; the ladder/oracle may evaluate
its ordinary exact scalar as before.

Separate PyPy 3.11 workers isolate int and mpz traces. Each warmed group uses
at least three seconds of validated warmup and nine repeated cohort samples;
relative IQR above 15% extends to 5 seconds/18 samples, then 8/27. Instability
remaining at the cap cannot earn promotion. Cold processes and upstream C
construction are reported separately. The shared flock and process inventory
exclude B4/A6 heavy checks and timing overlap.

The frozen search gate requires compact Lucas to beat compact PRAC by 5%
with stable reused full-stage timings and be within 10% of the ladder.
Production consideration additionally needs a stable 10% full-stage and
complete-campaign gain plus fresh confirmation. A missed gate stops search
expansion; it does not imply that every possible differential chain is slow.
