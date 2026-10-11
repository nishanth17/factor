# QS/SIQS parameter audit — 10 October 2026

This is a source-to-current-control inventory for the continuing C10 roadmap
work, not a new parameter selection or performance result. Read it with the
[B1 calibrated controls](README.md#b1-joint-qsmpqssiqs-calibration--9-october-2026),
[C1 DLP research](c1_research.md) and [C1 complete-run report](c1_implementation_results.md).
The inspected external implementations use different arithmetic, memory
layouts, certainty rules and hardware. None supplies a universal optimum for
PyPy. No external code was copied or run for this audit.

`B` below means the **largest admitted factor-base prime**, unless a source
explicitly uses another definition. v2 `base_bound` is an exclusive prime
*value* bound; msieve's `fb_size` and Yamaquasi's factor-base size are prime
*counts*. v2 `half_width=M` sieves `[-M,M)`; an upstream interval may mean the
full `2M` positions. Hash slots, retained partials, owned atoms, cycle rows,
and matrix rows are different resources. Normalize these before any sweep.

## Current choices and source differences

| Parameter family | v2 current choice and evidence status | Inspected external choice or contrast | Experiment implication |
| --- | --- | --- | --- |
| Input-size selection | [`SIQSConfig`](../qs/siqs.py) defaults to base bound 1,000 and half-width 512; there is no automatic QS/DLP digit selector. B1 validated explicit balanced 30/40-digit bundles; C1 validated one explicit 40-digit DLP bundle. | [msieve 1.53 `mpqs.c`][msieve-release] interpolates a bit-size table; [Yamaquasi `siqs.rs`][yamaquasi-siqs] enables DLP by default above 256 bits, while msieve's `sieve.c` enables it from 282 bits. YAFU's current, unpinned [official configuration][yafu-config] lists SIQS overrides and a QS/NFS crossover option. | Treat external crossover points as candidate strata, including larger feasible inputs. Train complete-cost policies using only observable input/work metadata; preserve SLP and unresolved runs. |
| Factor base | v2 `base_bound=1,000` by default, with a 1,000,000 value cap; B1 uses 3,000/10,000 on selected balanced 30/40-digit controls. This is not a size-derived runtime default. | [msieve 1.53 `mpqs.c`][msieve-release] has `fb_size=3,000` at 200 bits, 55,000 at 283 bits and interpolates between entries. [Yamaquasi `siqs.rs`][yamaquasi-siqs] halves its chosen base *count* when DLP is active. [FLINT QS docs][flint-docs] expose factor-base count and restart with a larger base when necessary. | Compare actual prime counts and largest prime after construction, not the unlike numerical units. Pair base changes with interval, residual and relation targets; charge setup and matrix width. |
| Sieve interval and blocks | v2 default half-width 512; selected B1 balanced controls use 8,192 at 30 digits and 65,536 at 40; public half-width cap is 499,999. Default collector block width is 256, cap 4,096; selected B1 uses 4,096. | [msieve 1.53 `mpqs.c`][msieve-release] uses half-sieve 65,536 at 200 bits and 3×65,536 at 283 bits. [Yamaquasi `siqs.rs`][yamaquasi-siqs] uses 32,768-position blocks, one block through 180 bits and roughly 8–16 blocks in the 256–340-bit DLP band. | Measure full `2M` position counts, polynomial-switch cost, cache/RSS, and bucket/resieve costs. B1's wider-QS trial failed on useful rows; rerun only with a changed bottleneck or source-grounded joint bundle. |
| A factors, assignment and Gray quota | v2 default is 3 A factors, 16 families/pool entries, reference assignment, 4,096 coefficient trials for external-square MPQS and the full `2^(k-1)` Gray quota; streaming permits up to 32 factors and a finite family grant. B1 selected four/five flyer factors, 100,000-family allowance and eight/16 effective Gray polynomials at 30/40 digits. | [Yamaquasi `siqs.rs`][yamaquasi-siqs] chooses `k` by bit-size tiers (5 at 120–149 bits, 7 at 170–199, `bits/25` above 199) and scales its A count with size; it deliberately uses short intervals. [FLINT QS docs][flint-docs] generate small prime subsets for A. | Jointly test A quality, pool/tolerance, family count and Gray quota with width. Count duplicate or exhausted assignments, root preparation, checkpoint growth and useful relations per polynomial. |
| Multiplier selection | v2 normally uses multiplier 1; explicit selector tests 13 odd squarefree candidates through 31 with integer scoring against primes through 97. Its candidate set and score are not a source-derived optimum. | [msieve 1.53 `mpqs.c`][msieve-release] considers even and odd multipliers through 73 and scores up to 300 test primes, constrained by input/word capacity. [FLINT QS docs][flint-docs] expose multiplier-search effort. | Compare finite candidate sets and complete-run cost on training bands; independently verify roots, factor-base exceptions and proper GCD recovery before considering even multipliers. |
| Candidate scoring and tiny primes | v2 uses conservative exact integer score bounds at 32 units per bit, 10-bit mantissa lookup, `threshold_extra=0`, no skipped small primes and 64-entry metadata chunks by default. Optional powers/fixed scores, bytearray, cutoff, bucket and resieve policies have separate R2 evidence; selected B1 uses powers/bucket. | [SymPy QS][sympy-qs] uses scaled-log candidate scoring and `ERROR_TERM=25`; [FLINT intake][flint-intake] uses byte scores and blocked scanning; [msieve sieve][msieve-release] changes unsieved-small-prime treatment with input size. These numerical thresholds are not in v2 score units. | Normalize loss/false-candidate rates under exact division. Preserve conservative coverage unless explicitly testing a lossy arm; measure threshold, prime-power and A-exception contributions, not only candidate count. |
| Single-large-prime endpoint | v2 nested SIQS default residual limit is 10,000 (bare `SieveConfig`: 500); B1 balanced 30/40-digit controls use 9,000,000/100,000,000. C1's 40-digit DLP arm narrows SLP to 25,000,000. Residuals admitted as primes are proven. | [SymPy QS][sympy-qs] accepts a prime cofactor below `128B`; [FLINT intake][flint-intake] uses `<60B` plus a 30-bit small-cofactor gate, without v2's explicit intake proof. [msieve 1.53][msieve-release] uses bit-tiered LP multipliers (50 at 200 bits, 80 at 283). [Yamaquasi `siqs.rs`][yamaquasi-siqs] has a size-dependent multiplier. | Convert each rule to the actual `B` and certification domain. Track matched partials, unmatched occupancy, false candidates, proof cost, useful filtered rows and complete factors. |
| DLP endpoint/product and splitting | Opt-in C1 40-digit arm uses exclusive base bound 5,000, SLP 25M, endpoint 500,000 (`100×base_bound`), product 3.2B (`128×base_bound²`), candidate bound 0 and at most 131,072 splitting calls. These multipliers are **not** exact multiples of the largest admitted prime. v2's splitter uses at most two 2,048-evaluation rho attempts per call; this is a finite allowance, not a tuned splitter. | [msieve 1.53 `sieve.c`][msieve-release] enables DLP at 282 bits and screens products at about `SLP_limit^1.8`; [Yamaquasi `siqs.rs`][yamaquasi-siqs] uses about `100B²` near 200 bits with a separate endpoint limit. [JavaMath `TDiv_QS_2LP`][javamath-tdiv] chooses Hart below 46 cofactor bits, tiny ECM below 63, rho below 64, then nested QS, admitting at most 31-bit endpoints under its usual path. [YAFU intake][yafu-tdiv] separates product/endpoint bounds and uses micro-ECM. | Test endpoint, product, candidate and splitter allowances separately first, then source-supported combinations. Retain representative rejected residuals and charge every primality/split attempt. Never infer prime certainty from a native probable-prime check. |
| Partial/cycle retention | v2 SIQS defaults to 512 pending partials, 2,048 rows and 4,096 atoms; selected C1 40-digit DLP uses 8,192 unowned edges/rows and 32,768 atoms. `MAX_GRAPH_EDGES=65,536` and 256 atoms/cycle are safety caps; graph reservation is `32,768 + 4,096×edges` bytes. Fresh 50-digit C1 runs had median 21,504/28,672 evictions and no useful filtered matrix. | [YAFU filtering][yafu-filter] grows its relation list by 1.5× and cycle tables by 2×; [msieve 1.53 `sieve.c`/`relation.c`][msieve-release] starts its cycle table at 10,000 and doubles it, but its cycle path has a finite 100-edge-per-side heuristic. [Yamaquasi relations][yamaquasi-relations] uses dynamic partial maps and restricts DLP merges to the SLP-connected component as a yield heuristic. | C9 owns the 8,192/16,384/32,768 finite retention study with matched 256/512 MiB caps. Record peak live bytes, matching/eviction losses and post-filter rank; preserve complete all-component cycle correctness and ownership. An upstream hash-table size is not a retained-edge allowance. |
| Relation target, filtering and dependencies | v2 defaults to `row_excess=2`, `filter_row_growth=1`, batch width 256 and weight-two filtering; B1/C1 selected `row_excess=32`, growth 32 and batch width 4,096. Matrix hard caps are 65,536 rows and 100,001 columns, subject to a much tighter live-memory reservation. | [SymPy QS][sympy-qs] collects at least 105% of factor-base cardinality. [Yamaquasi relations][yamaquasi-relations] asks for 48 extra kernel relations beyond covered columns; [msieve 1.53][msieve-release] uses 64 extra and targets `fb_size+96` before filtering. Those targets reflect different solvers and row quality. | Vary solve cadence/oversampling only after measuring actual post-filter rank and proper-divisor yield. Charge filtering, provenance, matrix, extraction and repeated trivial dependencies; do not promote on raw relation count. |
| Recovery and stopping | v2 has `growth_steps=0` (so default `max_half_width=8,192` is dormant), `max_stalled=16` no-row windows and `max_trivial=128` trivial dependencies; selected B1/C1 bounds differ. These are finite stop/recovery policies, not estimated yield optima. | [Yamaquasi `siqs.rs`][yamaquasi-siqs] scales its planned A count by input bits; [FLINT QS docs][flint-docs] describe restarting with a larger base when relations are insufficient. Neither maps directly to v2's resumable work ledger. | Diagnose stalled windows versus low-quality dependencies before changing growth or stop rules. Charge repeated setup/old work and preserve resumable exhaustion semantics. |
| Capacity, resume and score safety | v2 has hard caps including 4,096 input bits, 1e12 proven residual endpoints, 1M factor-base prime bound, 1M-position collection windows, 4,096 sieve block width, 64 Hensel lift roots, 256 atoms/cycle, 65,536 graph/matrix rows, and serialized checkpoint allowances. The SIQS default is 32 MiB owned memory and a 1 MiB checkpoint; selected B1/C1 use 256 MiB. These protect exactness, bounded work or restoration; they are not claims of economic optimality. | Native growable lists, 32-bit/64-bit cofactor limits, GPU SLP witness-table slots and v2 owned-memory reservations do not measure the same object. The current, **not yet pinned** [CUDA-MPQS guide][cuda-guide] lists a 1,048,576-slot GPU SLP witness table; it does not implement a DLP graph. | Audit every cap by role: proof, representation, or policy. Change proof/storage limits only with independent correctness, peak simultaneous-memory, cancellation and charged-resume evidence. Measure process RSS separately from owned reservations. |

The [prior source/license review](qs_gnfs_research.md#external-implementations-and-licenses)
records inspected YAFU/msieve public-domain notices, Yamaquasi/SymPy BSD-3,
FLINT QS LGPL-3-or-later and JavaMath GPL-3-or-later. Dependencies have their
own terms; no source copying is proposed. The [official msieve 1.53 archive](c1_research.md#release-cross-check-and-current-research)
was checked against publisher SHA-256
`c5fcbaaff266a43aa8bca55239d5b087d3e3f138d1a95d75b776c04ce4d93bb4`.
The live CUDA guide must be pinned and license-checked before its settings
enter a reproducible comparison. Native timing is contextual, not a PyPy
prediction.

The first queue is evidence-led. C9 should separate 50-digit retention loss,
product/threshold loss and insufficient collection duration; its fresh 120 s
runs had evictions but zero surviving matrices, so a larger graph alone has no
established benefit. A feasible longer-window case can then compare linked
base/interval/A schedules against SLP under one resource cap. Revisit multiplier
selection or score policies only when their setup or candidate costs are
material in that band's profile. Defer solve-cadence tuning until actual useful
post-filter rows exist. These are ordered bounded gates, not an unlimited
search or a claim about 80–100-digit outcomes.

## Continuing experiment ledger (C10)

The first audit identifies the gaps; it does **not** choose larger defaults.
For each new feasible input band, arithmetic backend, collector or solver change:

1. Record its source pin, license, actual parameter units, hardware/backend,
   v2 default and selected control, and the local evidence status. Update this
   table when a source or implementation changes; keep accepted old controls
   versioned rather than rewriting their evidence.
2. Use stage attribution and representative rejected residuals to rank at
   most a few independent parameter families by expected complete-factor
   value. Freeze inputs/seeds, total work/wall/CPU, owned bytes, RSS reporting,
   checkpoint size, sample count and stop rule **before** training. C9 handles
   DLP storage/product economics; B1 is the 30/40-digit control; E1 combines
   selected challengers; C3/G1 own portfolio handoff. No blind Cartesian grid.
3. Tune on training inputs only; freeze each finite bundle before fresh
   confirmation. Compare SLP/DLP and QS/MPQS/SIQS under matched total
   resources. Count setup, sieve, splitting/certification, graph/eviction,
   filtering, matrix, provenance, extraction, classification and resume.
   Retain failures and unresolved cofactors. Use PyPy implementing Python
   3.11, at least three seconds of validated warmup, nine samples (more if
   unstable), separate cold/profile runs and the machine-wide timing lock.
4. Accept a band-specific opt-in or default only under the revised roadmap
   promotion policy with fresh complete-factor/completion evidence, exact
   reconstruction, proper divisors, certainty preservation and bounded
   checkpoint replay. Publish the chosen bundle and the rejected candidates;
   if pilots cannot reach a useful matrix within their declared caps, stop and
   record the limiting stage and the next concrete trigger.

The item stays open as a maintained decision ledger. A completed experiment
closes only its named size/workload tranche, not the whole parameter space.

[msieve-release]: https://sourceforge.net/projects/msieve/files/msieve/Msieve%20v1.53/msieve153_src.tar.gz/download
[yamaquasi-siqs]: https://github.com/remyoudompheng/yamaquasi/blob/3f95f43682ed15d8c1ed206a9a702dd655d7c8ad/src/siqs.rs
[yamaquasi-relations]: https://github.com/remyoudompheng/yamaquasi/blob/3f95f43682ed15d8c1ed206a9a702dd655d7c8ad/src/relations.rs
[yafu-config]: https://github.com/bbuhrow/yafu/blob/master/yafu.ini
[yafu-filter]: https://github.com/bbuhrow/yafu/blob/8110dfbd8c6f9486d93b1a02de6eb7b180e55a80/factor/qs/filter.c
[yafu-tdiv]: https://github.com/bbuhrow/yafu/blob/8110dfbd8c6f9486d93b1a02de6eb7b180e55a80/factor/qs/tdiv.c
[flint-docs]: https://github.com/flintlib/flint/blob/a4c9750d0d3d67bb01cf6d18c187591b313451c3/doc/source/qsieve.rst
[flint-intake]: https://github.com/flintlib/flint/blob/a4c9750d0d3d67bb01cf6d18c187591b313451c3/src/qsieve/collect_relations.c
[sympy-qs]: https://github.com/sympy/sympy/blob/2f22a5f81e2f4124380be3739a092e9ff20128de/sympy/ntheory/qs.py
[javamath-tdiv]: https://github.com/TilmanNeumann/java-math-library/blob/088d01fa97e7d0412c6cbbb8f463fa7e5f78ab97/src/main/java/de/tilman_neumann/jml/factor/siqs/tdiv/TDiv_QS_2LP.java
[cuda-guide]: https://github.com/drjanosch42/cuda-mpqs/blob/main/USER_GUIDE.md
