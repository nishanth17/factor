# Phases 3+ research and implementation reconciliation

Reviewed: 4 October 2026. Target: Python algorithms in `v2/` on PyPy
implementing Python 3.11. This is a research and planning change, with
source-frozen exploratory profiles and structural probes, but no algorithm
implementation, measured speedup comparison or default promotion.
Refreshed after the budget follow-up and subsequent SSS/checkpoint work;
the same report consolidates the initial review and the new findings.

The review compares the active worktree at HEAD
`f871a4a3bb05a2eab25490a4a0a1026f6c87c44f`, including existing uncommitted
work, with accepted milestone evidence and public research available on the
review date. HEAD alone cannot reproduce the reviewed worktree. Local content
hashes preserve that distinction; upstream references below are immutable
where source repositories provide commits.

Several implementation/benchmark files also changed elsewhere during this
pass. The local evidence records beginning and final hashes; affected paths
were reread and final tests/lint rerun. These concurrent changes were left
untouched and are distinct from this report's documentation additions.

The authoritative action list is [TODOS.md](TODOS.md). All fresh Phase 3
findings are folded into **P3.8**, as requested. P3.1–P3.7's existing text and
accepted records receive no edits from this research pass. Concurrent P3.5
work records its completed bounded challenger and upstream reproduction;
P3.8-R5 retains that evidence and owns the remaining fresh SIQS comparison.
**P3.4's implementation is complete**; its
documented meaningful large-number experiment remains open. **Phase 8 is the
M23 follow-up research/experiment backlog for Phase 2 optimization.** Its
ownership remains intact, including permission for useful independent spikes
before integrated portfolio evaluation.

Earlier reports remain useful controls:
[QS/SIQS research](quadratic_sieve_research.md),
[GF(2) matrix research](gf2_matrix_research.md), and
[Phase 2 optimization research](phase_two_optimization_research.md).
This pass reconciles and extends their decisions rather than resetting them.

## Findings and priority

The roadmap's algorithm families remain appropriate. The largest gaps are
the evidence needed to establish a useful SIQS workload, representation costs
before iterative matrix solving, and precise GNFS field/square-root contracts.
Advanced collectors and arithmetic kernels should remain independent measured
challengers. A bounded small GNFS reference need not wait for all of them.

| Priority | Finding | Action owner | Decision from this review |
| --- | --- | --- | --- |
| First | SIQS family/base/window/store limits jointly determine feasible coverage | P3.8-R1 | Establish a trained capacity control and fresh confirmation cohort |
| First | Candidate rejection, residual work and useful-row yield matter together | P3.8-R2 | Evaluate the present power-score and batch experiments end to end |
| First | Repeated preparation/filtering and dense provenance can limit capacity before solving | P3.8-R3 | Compare queues/batches, compact constraints and deferred lifting |
| First | Nonmonic norms, bad ideals and square roots need explicit supported-field contracts | P7.1–P7.4 | Specify them before enlarging GNFS coverage |
| Next | GMP availability and conversion overhead need actual PyPy evidence | P4.3 | Retain the Python-int reference; assess whole-stage boundaries |
| Next | Extra partials or kernels can still give only trivial congruences | P3.8-R4 / P5.4 | Measure verified independent dependencies and proper-divisor yield |
| Next | NFS polynomial quality and two-sided cofactoring require joint empirical selection | P7.6 | Use bounded trial sieving and staged cofactor experiments |
| Conditional | PRAC, p+1, pairing, Edwards, polynomial continuation and workers | P4–P6 | Strengthen exact contracts, then adopt/defer/reject separately |
| Conditional | Block Lanczos/Wiedemann, NumPy, SSS and richer collectors | P3.8 / P6.2 | Advance when workload evidence justifies their total costs |

These priorities are engineering inferences from the code and literature.
They are not universal bottleneck rankings or promises of speedup. The latest
user priority is the MPQS/SIQS/ECM optimization track below; GNFS contract
research remains available without delaying these targeted investigations.

## Below-100-digit optimization target

The final scope clarification targets **general factoring below 100 decimal
digits**, using MPQS/SIQS and ECM on this machine. **Do not use a 50-digit
cohort as the benchmark.** Previously collected 50-digit profiles remain
exploratory provenance; the planned repeated 50-digit timing baseline was
cancelled. Historical P3.4/P3.5 experiment records are unchanged. New workload,
capacity and collector actions are in P3.8; ECM changes are in P4/P5.

Read-only hardware inspection reports Apple M4 (`Mac16,13`), ARM64, 10
physical/logical CPU cores and 24 GiB installed RAM. This does not establish
free RAM, sustained clocks or the useful worker count. The project-local
PyPy 3.11 venv imports gmpy2 2.3.1/GMP 6.3.0. CUDA and x86 vector results are
algorithm/design evidence, not directly usable M4 kernels or speed forecasts.

### Capacity and polynomial quality come before a fair method comparison

There are stronger restrictions than the global timer. Current MPQS chooses
`A=q²` with q restricted to the factor base; its hard bound 100,000 implies
`A<10**10`. Classical MPQS can choose q outside the base near the square root
of the desired A. Because A is a known square, its contribution can be carried
as a modular square correction instead of demanding that it factor over the
base. This needs a revised collector/relation representation, independent
root and square-identity checks, and bounded coefficient construction. It is
not enough to change q selection while retaining the present "A must factor
completely over the factor base" constructor. Square corrections do not
require q itself to be prime, but root-building routines that assume primality
need their own validated domain and finite fallback; a probable-prime label
must not silently become a proof. [Zimmermann's MPQS formulas and dual batch
inversion](https://members.loria.fr/PZimmermann/talks/tiny-mpqs.pdf),
[primary implementation thesis, §§6.3–6.5](https://martinlauridsen.info/pub/bsc_thesis.pdf).

SIQS has a separate upper-range barrier. It permits at most eight distinct
factor-base primes in A, so `A<100000**8=10**40`. Even for the smallest
93-digit n and largest permitted half-width 499,999,
`floor(sqrt(2*n)/M)>2.8*10**40`; at 99 digits it exceeds `2.8*10**43`.
Thus **every 93–99-digit input has an unreachable intended A target under
the current caps**, independent of wall time. The computation uses integer
roots and h=1; positive multipliers only enlarge the target. This establishes
polynomial target mismatch, not mathematical impossibility of factoring with
an off-target A. The family/search and storage ceilings also remain finite.

The latest stable PARI/GP source, **2.19.0 released 28 September 2026**, is a
useful CPU scaling reference. Its function named MPQS implements a
self-initializing squarefree-A variant, so it belongs to the SIQS comparison
as well. `src/basemath/mpqs.h` requests the following native configurations
by the digits of **kN**, not necessarily n:

| kN digits | Factor-base prime count | A factor count | Half-width M |
| --- | ---: | ---: | ---: |
| 70 | 10,500 | 8 | 176,000 |
| 80 | 21,000 | 9 | 448,000 |
| 90 | 40,000 | 10 | 512,000 |

Factor's bound permits only 9,592 primes in total before residue selection,
and its eight-factor/499,999-width caps cannot express the latter regimes.
These are native design controls, not proposed PyPy defaults or measured
comparisons. PARI explicitly marks the table's 92-and-above entries as never
tested; the presence of a 99-digit row does not validate it. The reviewed
release archive matches its official SHA-256
`f317b9722eb5d9094a60303774f066f3a83e3ec1f170be8546c44d7583f30b6d`.
[Official release](https://pari.math.u-bordeaux.fr/download.html),
[reviewed source archive](https://pari.math.u-bordeaux.fr/pub/pari/unix/pari-2.19.0.tar.gz).

For A quality, PARI chooses all but one A prime then selects a final
compensating "flyer" prime to approach the actual product target. Factor
instead samples a pool ordered by distance from `p**factor_count` to target.
The PARI policy is a concrete challenger when product spread is poor; it is
not necessarily the first fix when existing A products are already close.
Stream an extendable assignment cursor with bounded duplicate state and
resident batches. Raising A factor count also increases the exponential Gray
family length, so preallocating all families/history defeats the capacity fix.

### Current QS bottleneck evidence and actionable changes

After the scope correction, two additional **instrumented diagnostics** used
one existing training input each at 40 and 80 digits, not a 50-digit benchmark.
Both ran the same frozen source and seed, B=50,000, M=32,768, 4,096-position
blocks, residual bound 10**8 and input-scaled A factor count (4 and 8).
Each had at least three seconds of validated bounded-job warmup and a
10-second active-run wall/CPU allowance, with ample work/storage allowances.
Every returned unresolved result reconstructed its input. These settings are
not claimed to be optimized, and the short warmup need not cover every later
path. cProfile changes PyPy/JIT behavior; generator call counts include
resumptions, cumulative timings overlap, and no timing share is a speedup.

| Diagnostic | 40 digits | 80 digits |
| --- | ---: | ---: |
| Scanned positions | 294,305 | 327,685 |
| Coarse candidates | 38,976 | 40,960 |
| Rejected at refined scoring | 36,556 | 40,960 |
| Admitted atoms | 1,401 | 0 |
| Usable full/combined rows | 70 | 0 |
| Unmatched partials | 1,328 | 0 |
| SLP matches | 3 | 0 |
| Calls to `Budget.consume` | 1,254,700 | 1,441,471 |

Neither produced a factor or established a production capability limit.
Collection occupied about 8.95/9.17 instrumented seconds; preparation and
filtering together were about 0.32 seconds at 40 digits and negligible at
80. Matrix elimination was negligible in both. The useful priority is
**collection and capacity**, with matrix representation addressed as a later
scaling barrier. In particular:

1. **Reduce accounting/polling frequency over bounded chunks.** `consume`
   repeatedly validates types and reads wall/CPU clocks; each trace records
   over a million calls. Batch exact work reservations and poll deadlines at
   declared finite chunk boundaries, preserving cancellation latency and
   refusal before any candidate/store commitment. `_divide` also charges a
   dense base-sized estimate before refined rejection: all 40,960 candidates
   at 80 digits were rejected, yet this configuration consumed about 15.8
   billion logical units. Align charges with actual refinement/sparse work
   rather than treating logical units as seconds or simply removing controls.
2. **Reuse polynomial prime-power plans.** `_sieve` restarts prime-power
   root lifting over the whole base for each working block and performs a
   separate base-root/hit pass. Cache or stream a bounded per-polynomial
   interval plan, combining first-level marking with hit metadata. Validate
   high valuations, p=2, p|A/N', singular roots, tails and cache identities.
3. **Use the sparse recovery already represented.** Bucket `_divide` still
   scans every base prime per surviving candidate despite recording a hit
   bitset. Iterate set bits plus sparse A support. `_resieve` also scans all
   candidates per prime to add mostly zero A exponents; initialize those
   contributions once, then append actual hit valuations. Keep exact full
   division as the independent oracle.
4. **Make skipped-prime scoring useful.** Omitting 2 in the present power
   policy grants a maximum-bit-length allowance, collapsing the coarse and
   refined thresholds to zero. Compare cheap exact tiny-prime corrections,
   bounded fixed-point logs and staged refinement. Native approximate scoring
   may intentionally miss candidates: measure that loss explicitly while
   retaining exact admission, residual checks and extraction.
5. **Tune useful partial yield, then add DLP.** Larger SLP bounds can fill the
   store without collisions; PARI's source explicitly warns about this.
   Measure matches, eviction, residual-splitting cost, graph cycle space and
   post-filter usable excess. DLP is a high-value scaling challenger with
   independently bounded prime/product limits and exact cycle corrections,
   not permission to count raw partials as progress. YAFU's implementations
   supply concrete resieving/cofactoring references. [Pinned collector](https://github.com/bbuhrow/yafu/blob/963dbe9c45283cc06e9a71d830b0676bc1b0d343/factor/qs/SIQS.c),
   [resieving](https://github.com/bbuhrow/yafu/blob/963dbe9c45283cc06e9a71d830b0676bc1b0d343/factor/qs/tdiv_resieve_32k.c).

Kleinjung/CADO family-wide hit scheduling remains conditional: its p>I
condition concerns the **whole polynomial interval**, not the 4,096-position
working block. With M=32,768 and B=50,000, no factor-base prime meets it.
Prime powers or larger future bases may qualify, but this is not the direct
fix for these traces' repeated root/polling cost.

Before matrix scaling, add stable mixed admission indices. The current
`_full + _combined` row ordering shifts old combined indices when new full
rows arrive, so naive append-count or dependency-mask caches are unsafe.
Atomic IDs omit exponent payloads; cache complete immutable row/base/store
identity. Then compare solve cadence, degree queues/disjoint merge batches,
live-column compaction, merge histories and sparse selected-bit extraction.
Structural probes found a 512-row sparse cycle taking 512 complete filter
rounds and 134,611,200 charged units; 32 one-bit rows with a single high column
label reserved 54.45 MB versus 50,944 bytes under compact labels. These are
scaling diagnostics of their captured source, not timing improvements or
present solver dominance. Concurrent work has since added touched-column
incidence updates and degree queues. The repeated-rebuild probe is therefore
historical; evaluate the new filter against its frozen control. Dense lifting
reservation and uncompacted labels remain separate concerns.

### ECM priorities across the whole range

ECM should be stratified by **factor size** as well as modulus size. An uneven
99-digit input with a 20-digit factor and a balanced 99-digit semiprime need
very different curve investment. Native curve tables describe heuristic
success distributions; even the expected number of curves leaves substantial
failure probability. Neither a few small tiers nor a short timeout fairly
represents serious ECM coverage. Preserve finite allowances, but permit
explicit workload-sized tiers and cumulative extension. The existing work
probe needs 4,081,645 logical units for one B1=50,000/B2=5,000,000 curve;
the 2,000,000-unit default stops during stage two. This is verified refusal
behavior, not an expected time-to-factor or a sub-100 benchmark.

Storage also needs a suitable input envelope. The current portfolio formula
reserves 9,304,064 bytes for an 11,000/1,900,000 ECM tier under its default
4,096-bit envelope, exceeding the default 8 MiB before execution. All inputs
below 100 digits fit 329 bits; the same tier reserves 2,665,104 bytes under
that campaign envelope. A 250,000/130,000,000 tier still reserves 16,216,200
bytes before sieve workspace/RSS. These are exact source-formula evaluations,
not memory measurements or proposed optimal tiers. Work, storage, target
factor size and curve counts must be configured together.

The highest-value ECM changes are a **compiled, reusable paired stage-two
coverage program and schedule amortization**, alongside a **coarse int/mpz
arithmetic boundary**. Current stage jobs regenerate small prime segments;
the eight-entry cache churns across multi-million B2 schedules, while its
power/gap interfaces are unused by production jobs. Keep immutable bound-owned
instructions reusable across curves and curve points private. CADO's planner
combines ± candidates and uses divisibility coverage to prune terms; copy
the coverage principle, not native absolute-B2-sized arrays.
[Stage-two planner](https://github.com/cado-nfs/cado-nfs/blob/692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b/sieve/ecm/stage2.c).

Point kernels currently allow roughly five-modulus-width intermediates before
final reduction. Compare selected early reductions, fused ladders and
unit-checked fixed-difference normalization; extra reductions may cost more.
The installed mpz doubling probe agrees with int, but public ladder validation
rejects an mpz modulus. Whole-loop specialization and canonical checkpoint
boundaries are required, rather than per-operation dispatch or a simple type
cast. Verified PRAC/prime-power chain programs must reach `stage_jobs.py`,
which presently calls the ladder directly; `multiply_prac` remains a wrapper.
The 2025 continued-fraction search is an offline/small-scalar challenger,
not an online search over a huge stage-one lcm. [Chain paper](https://doi.org/10.1007/s40993-024-00604-8),
[compact chain execution](https://github.com/cado-nfs/cado-nfs/blob/692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b/sieve/ecm/README.bytecode).

Pairing/common-Z, schedule reuse, coarse GMP and verified chains come before
polynomial/FFT continuation or another sparse solver. Assess p−1/p+1 by
marginal portfolio success, retaining structured controls separately. Compare
1/2/4 workers on the M4 with total CPU/RSS, startup and cancellation charged;
ten logical cores do not imply tenfold first-factor speedup.

Native allocators distinguish a cheap automatic pretest from an explicit
ECM-only/factor-target campaign. Pinned YAFU has target-factor states at
15/20/25/30/35/40/45 digits with rising B1, credits prior completed curves and
permits a custom pretest depth. yamaquasi likewise has separate automatic,
ECM-only and smooth-factor APIs. Use these policy distinctions, with the
M4's measured marginal costs, rather than repeatedly spending the same tiny
allowance or forcing ECM to seek half-size factors for every input.
[YAFU allocator](https://github.com/bbuhrow/yafu/blob/963dbe9c45283cc06e9a71d830b0676bc1b0d343/factor/autofactor.c),
[yamaquasi ECM](https://github.com/remyoudompheng/yamaquasi/blob/3f95f43682ed15d8c1ed206a9a702dd655d7c8ad/src/ecm.rs).

At larger factor targets the continuation algorithm must also scale.
yamaquasi switches direct stage-two products to polynomial root evaluation
at a native d1 threshold of 4,000. That is concrete evidence for bounded
product/remainder-tree and multipoint-evaluation experiments after paired
classical stage two; the numeric crossover depends on its native arithmetic.
It supplies no PyPy speed estimate. P6.2 already owns this optional extension;
P5.2 should record the trigger and preserve exact composite-modulus/nonunit
handling, node/coefficient/storage bounds and checkpoint reconstruction.
[Pinned continuation path](https://github.com/remyoudompheng/yamaquasi/blob/3f95f43682ed15d8c1ed206a9a702dd655d7c8ad/src/ecm.rs).

### Replacement evaluation design

Prespecify representative **30/40/60/70/80/90/99-digit** total-size classes,
balanced inputs and uneven classes with feasible smaller-factor sizes such
as 10/15/20/25/30/35/40 digits. Do not reuse inspected confirmation data as
untouched held-out evidence. Separate tuning, fixed-work collector/kernel
diagnostics and complete factorization. Freeze seeds/configurations/resource
allocations before fresh confirmation; record every time, work, schedule,
storage and yield stop. Work units differ by method, so equivalent wall/CPU/
memory resources need sufficient per-method logical allowances rather than
an arbitrary identical work integer that prevents one method starting.

Promotion still requires at least three seconds of validated PyPy warmup,
nine repeated samples with extensions for unstable measurements, exact
reconstruction/certainty, and the existing completion/regression policy. Upper
bands can use prespecified censored trials; report useful progress separately
and promise no successful-factor speed ratio when neither arm finishes.
Measure native PARI/YAFU/yamaquasi/GMP-ECM as explicitly different backend
reference arms when available, with builds/cores/memory recorded.

The practical crossover is machine and implementation dependent. Current
YAFU derives it in `tune()` rather than establishing its commented 97-digit
override as a universal rule. yamaquasi's author reports useful near-100-digit
native SIQS runs but labels its timing method nonrigorous. This motivates
aggressive MPQS/SIQS optimization throughout the requested range while
keeping the existing GNFS crossover experiment available at the upper end;
it does not establish that pure PyPy ECM/QS will quickly finish every balanced
99-digit input. [YAFU tuning](https://github.com/bbuhrow/yafu/blob/963dbe9c45283cc06e9a71d830b0676bc1b0d343/factor/tune.c),
[yamaquasi collector design and author timing qualifications](https://github.com/remyoudompheng/yamaquasi/blob/3f95f43682ed15d8c1ed206a9a702dd655d7c8ad/README_qsieve.md).

## Scope, source quality and reproducibility

The initial pass searched papers and author pages, read primary papers and
implementer reports, and fetched 52 selected source/documentation files from
eight pinned repositories. The refresh and below-100-digit extension expand
that total to **95 distinct repository files across ten pins**, plus four
selected files from the hash-verified PARI 2.19.0 release archive, additional
primary papers and an author-code archive. A separate recent RSA-oracle
repository was screened
at README level; it is not a factoring implementation comparison.
Inspection was targeted to relevant contracts and
algorithms; this was not a complete audit of all upstream code. No downloaded
implementation was executed or benchmarked. Repository trees and per-file
SHA-256 hashes are retained locally with the fetched excerpts.

| Source | Full commit pin | Commit date, UTC | Use and qualification |
| --- | --- | --- | --- |
| [CADO-NFS](https://github.com/cado-nfs/cado-nfs/tree/692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b) | `692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b` | 2026-10-02 | Current GNFS stage contracts, filtering/replay, characters, bad ideals, polynomial scores and square root |
| [YAFU](https://github.com/bbuhrow/yafu/tree/963dbe9c45283cc06e9a71d830b0676bc1b0d343) | `963dbe9c45283cc06e9a71d830b0676bc1b0d343` | 2026-09-24 | Native SIQS collection, resieving, batch cofactoring and relation filtering |
| [yamaquasi](https://github.com/remyoudompheng/yamaquasi/tree/3f95f43682ed15d8c1ed206a9a702dd655d7c8ad) | `3f95f43682ed15d8c1ed206a9a702dd655d7c8ad` | 2025-11-23 | Rust collector and GF(2) design; fixed-width and benchmark-method limits apply |
| [FLINT](https://github.com/flintlib/flint/tree/00cc19b350b5302876b4fa8c877df5a676c2061c) | `00cc19b350b5302876b4fa8c877df5a676c2061c` | 2026-10-02 | Native QS collector/filter/solver and lifecycle reference |
| [GMP-ECM mirror](https://github.com/sethtroisi/gmp-ecm/tree/8ea5e214fdf2f0ddf9415141b8dc039ed6f5874e) | `8ea5e214fdf2f0ddf9415141b8dc039ed6f5874e` | 2026-03-09 | ECM/p±1/continuation and Lucas chain code; README identifies an outdated mirror, not current official HEAD |
| [SSS author code](https://github.com/sbaresearch/smoothsubsumsearch/tree/8dbaf6d39ab88a40380965d25ec2c363d7f27358) | `8dbaf6d39ab88a40380965d25ec2c363d7f27358` | 2024-06-06 | Author implementation and experiment semantics, separate from Factor's adapter |
| [Historical msieve mirror](https://github.com/radii/msieve/tree/c8727d91305bdbe0972d160ef0ce61dd02ce9193) | `c8727d91305bdbe0972d160ef0ce61dd02ce9193` | 2011-02-21 | Historical QS/Lanczos/GNFS design; not a current msieve release |
| [Later msieve mirror](https://github.com/MersenneForum/msieve/tree/03e0ab9e13c4862e5bb2c51a4a39718aa2500e19) | `03e0ab9e13c4862e5bb2c51a4a39718aa2500e19` | 2016-10-26 | README/QS documentation only; still historical, not verified current official HEAD |
| [CUDA-MPQS](https://github.com/drjanosch42/cuda-mpqs/tree/b1a9c4500da15a2bea56a59ecad7ef61cf054756) | `b1a9c4500da15a2bea56a59ecad7ef61cf054756` | 2026-10-04 | Recent GPU design and failure diagnostics; author timings were not reproduced |
| [Standalone GPU block Wiedemann](https://github.com/drjanosch42/block-wiedemann/tree/ebfe814b759803f61771b36d1fd197b71b751d1d) | `ebfe814b759803f61771b36d1fd197b71b751d1d` | 2026-10-04 | README/usage refresh only; base-case generator and solver contracts do not justify another Phase 3 workstream |

Selected fetched paths, relative to these pins:

- CADO-NFS: `README.md`; `dev_docs/README.las`, `README.filter`,
  `README.renumber`; `linalg/bwc/README`, `linalg/characters.cpp`;
  `parameters/factor/params.c60`, `params.c80`; `sqrt/sqrt.cpp`,
  `sqrt/montgomery.md`; `sieve/ecm/facul.cpp`; `filter/merge.cpp`,
  `filter/replay.cpp`; `polyselect/murphyE.cpp`, `polyselect/E.sage`;
  `utils/badideals.cpp`.
- YAFU: `factor/qs/SIQS.c`, `filter.c`, `tdiv.c`,
  `tdiv_resieve_32k.c`; `factor/batch_factor.c`; `yafu.ini`; `CHANGES`.
- yamaquasi: `README_qsieve.md`; `src/relations.rs`, `src/matrix/gf2.rs`,
  `src/sieve.rs`.
- FLINT: `src/qsieve/collect_relations.c`, `factor.c`, `block_lanczos.c`,
  `linalg.c`.
- GMP-ECM mirror: `README`, `NEWS`, `ecm.c`, `pm1.c`, `pp1.c`, `stage2.c`,
  `LucasChainGenerator/README`.
- SSS: `sssif/sss.py`, `sssif/mstep.py`, `test.py`.
- Historical msieve: `Readme`, `Readme.qs`, `mpqs/relation.c`,
  `common/lanczos/lanczos.c`, `gnfs/sqrt/sqrt.c`.
- CUDA-MPQS: `README.md`, `docs/modules/matrix.md`, `docs/modules/sqrt.md`,
  `src/matrix/preprocess.cpp`, `src/matrix/merge_tree.cpp`,
  `src/sqrt/sqrt_step.cu`.

Additional refresh paths at the same CADO pin: `sieve/siqs-largesieve.hpp`,
`siqs-fill-in-buckets.inl`, `siqs-smallsieve-glue.hpp`, `siqs-smallsieve.cpp`,
`siqs-smallsieve.hpp`; `sieve/las.cpp`, `las-qlattice.cpp`, `las-duplicate.cpp`,
`las-todo-list.cpp`, `las-side-config.cpp`, `las-cofactor.cpp`;
`scripts/cadofactor/README.md`, `cadotask.py`; `filter/purge.cpp`;
`sqrt/crtalgsqrt.cpp`; `parameters/factor/params.c90`;
`sieve/ecm/README.bytecode`, `bytecode.c`, `ecm.cpp`, `ec_parameterization.hpp`,
`stage2.c`, `pp1.cpp`, `pp1_stage2.hpp`, `ec_arith_Montgomery.hpp`,
`ec_arith_Edwards.hpp`, `ec_arith_common.hpp`. The standalone GPU solver
adds `README.md` and `USAGE.md`. Duplicate `las.cpp` fetches are counted once.
Separate refresh manifests retain each URL, byte count and SHA-256, including
the chain-code archive and two new paper PDFs; no upstream code was executed.

The below-100-digit extension adds YAFU `factor/autofactor.c`, `factor/tune.c`
and `docfile.txt`; yamaquasi `src/lib.rs`, `siqs.rs`, `mpqs.rs`, `params.rs`,
`ecm.rs`, `ecm128.rs`, `README_ecm.md`, `benches/ecm.rs`,
`scripts/ecm_study.py` and `ecm_chains.py`; and the later msieve mirror's
`Readme`/`Readme.qs`. The official PARI archive contributes targeted review
of `src/basemath/mpqs.c`, `mpqs.h`, `ifactor1.c` and `doc/usersch3.tex`;
the large manual file was searched for relevant factoring contracts, not read
in full. Its public release hash and selected file hashes are retained. The
95-file total excludes PARI archive members and screened README-only oracle
research. This is targeted inspection, not exhaustive audits of these engines.

Source limitations matter. CADO's filtering README marks material as partly
obsolete; current merge/replay code is the more specific reference. The
five-line `README.las` is entirely an obsolete notice, so the refresh replaces
its lattice-contract citation with live sieve and task code. Official GMP-ECM
GitLab access returned an anti-bot page; the mirror remains explicitly
qualified, with current CADO ECM code providing additional evidence. The
polynomial-ranking article's HAL download was blocked, so its author listing,
CADO code and the RSA-240/250 paper support the limited ranking discussion
below. The 2023 alternative-sieving work was reviewed through the authors'
2024 seminar slides when the ePrint PDF could not be opened. Recent
deterministic high-order and odd-prime-power root papers were screened at
abstract level; this report does not claim a complete proof review of them.

Fetched snapshots, per-file hashes, source-fetch scripts, the pre-edit plan
and local preservation/verification records are kept in the ignored local
audit directory `phase_3_plus_20261004`. Public pins and review scope above
permit refetching without committing generated captures or upstream dumps.

## What the implementation already supplies

| Local observation | Implication for the plan |
| --- | --- |
| Exact normalized QS polynomial, A/sign/correction provenance, SLP combination, extraction and SIQS checkpoints exist | Preserve accepted P3.1–P3.4 controls; focus on workload and total cost |
| [Family selection](../qs/families.py) has finite pool/factor/family caps | Longer wall time alone cannot create additional assigned families |
| [Collector](../qs/sieve_collector.py) supports conservative and power score policies | Compare rejection/selectivity and setup costs; presence is not promotion |
| [Power sieve](../qs/power_sieve.py) caps root branching with a conservative fallback | Test high valuations and singular cases; do not remove its safety margin without proof |
| [Smooth batch](../qs/smooth_batch.py) uses bounded product/remainder trees | It finds non-base parts; full exponents still belong to relation recovery |
| [SSS adapter](../qs/sss.py) has finite work, an intentionally lossy SSSf option and serialized checkpoints | Treat it as an independent adaptation; validate assignment/solver replay and cumulative resources |
| [Pipeline](../qs/pipeline.py) prepares and filters each changed store when no solve is pending | Preparation caching and solve cadence are profile candidates |
| [Filter](../qs/linear_algebra.py) now has touched-column incidence updates/queues from concurrent work, while dense lift reservations remain | Evaluate the changed filter against the frozen rebuild control; improve representations before assuming sparse solving is the bottleneck |
| `matrix_workspace` reserves worst-case fill and quadratic original-relation masks | Sparse input alone does not imply sparse owned capacity |
| [ECM](../ecm.py) uses a ladder; `multiply_prac` delegates to it | P4.1 remains a genuine challenger, with valid fallback |
| Bounded stage jobs support rho, p−1 and ECM; no GNFS or p+1 engine exists | P5.1 and P7 remain implementation work, not documentation-only closure |
| Current [portfolio](../portfolio.py) and [CLI](../factor.py) explicitly select SSS/SSSf; portfolio checkpoints now use version 4 | Distinguish selectable methods from evidence for automatic dispatch promotion; reconcile public API documentation |
| The local PyPy venv imports gmpy2 2.3.1 / GMP 6.3.0; system PyPy has no gmpy2 | A supported-runtime comparison environment exists; backend/performance gates remain open |

The filter stores relation rows with constraint bits. P3.8's mathematical
interface uses constraint rows and relation columns, `M d = 0`. These are
transpose conventions for the same selection problem; freeze orientation and
conversion explicitly before reusing a sparse solver.

The worktree also contains larger-band runners and new corpora beyond M30's
accepted evidence. Their existence is not proof of coverage, stable timing or
publication readiness. Already opened/tuned inputs are reference/training
data for this review, not untouched future confirmation data.
Concurrent P3.5/README/changelog edits record complete 30-digit SSS/SSSf
splits and unchanged upstream reproduction. Those experiments are distinct
from this pass, which did not rerun their comparison; the declared challenger
still awaits a feasible trained SIQS comparison. Preserve their evidence
when executing P3.8-R5 rather than treating the adapter as a new proposal.

## Phase 3 carry-over owned by P3.8

### Establish useful workloads before algorithm rankings

The existing M30 evidence distinguishes a correct SIQS implementation from
its open large-number experiment. Small completed fixtures exercise identities
and checkpoints. Short unsuccessful probes do not establish the practical
limit, and extending a timer while leaving assignments finite is not a
large-number evaluation.

P3.8-R1 should train a feasible envelope jointly: factor-base cardinality and
bound, available A products near the integer target, window width, number of
polynomials/families, partial occupancy, atom retention and matrix reserve.
Report the precise stopping reason. Explicitly distinguish an unavailable A
target, exhausted family schedule, no useful excess, storage refusal and a
timer/work limit. Profile stage costs separately from timing evidence.

Compare complete SIQS and ECM-to-SIQS against the accepted portfolio under
identical resources. Use fresh held-out confirmation inputs after selecting
configurations, at least three seconds of validated PyPy warmup and nine
samples, with instability-driven extension. Keep owned-workspace reservations
and process/JIT RSS distinct. No proposed digit band is a capability claim.

### Candidate rejection and cofactoring form one pipeline

YAFU's collector divides the work into sieve bands, resieving/trial division
and bounded residual handling; its batch code supplies a product/remainder
design reference. These are useful experiments for the current collector,
with exact exponent reconstruction retained. Its native vector widths and
block sizes cannot select PyPy parameters. [Pinned trial division](https://github.com/bbuhrow/yafu/blob/963dbe9c45283cc06e9a71d830b0676bc1b0d343/factor/qs/tdiv.c),
[resieving](https://github.com/bbuhrow/yafu/blob/963dbe9c45283cc06e9a71d830b0676bc1b0d343/factor/qs/tdiv_resieve_32k.c)
and [batch factoring](https://github.com/bbuhrow/yafu/blob/963dbe9c45283cc06e9a71d830b0676bc1b0d343/factor/batch_factor.c).

Factor's conservative policy assigns an upper bound for all possible p-power
contribution at a root. That protects selection but may admit many false
candidates. The present power policy can sharpen this, with setup and root
branching costs. This is a hypothesis to measure, not a new measured finding.
Enumerated small windows must cover 2-adic valuations, singular roots, p|A,
p|N', tails and saturation before interpreting rejection counts.
The final worktree additionally prunes root classes that cannot hit the
current window before lifting further. The safe-coverage oracle must include
this window-local optimization and its work accounting.

Record the entire cascade: positions, score survivors, exact divisions,
residual primality/split work, admitted partials, unmatched occupancy, complete
relations, post-filter useful excess and verified proper factors. Optimize
these jointly. A faster score loop can still increase residual cost or worsen
matrix fill. Batch smooth-part detection cannot substitute for the exact
factor exponents needed by extraction.

yamaquasi describes a specialized 2-adic polynomial treatment and gives
useful sieve/layout hypotheses, while qualifying its own benchmark method.
Its fixed-width Rust arithmetic and polynomial normalization differ from
Factor's arbitrary-size integer contract. A related experiment must derive
all powers-of-two/A/sign corrections rather than copy formulas or thresholds.
[Pinned collector notes](https://github.com/remyoudompheng/yamaquasi/blob/3f95f43682ed15d8c1ed206a9a702dd655d7c8ad/README_qsieve.md).

### Filtering and provenance before bigger solvers

The frozen filter's repeated full-incidence rebuild was a concrete design
cost; concurrent implementation now adds degree queues and touched-column
updates. Evaluate that change before claiming a complete-pipeline improvement.
Compare the updated queue control with bounded
independent pivot batches, then compare capped higher-weight minimum-fill
merges. Independent pivots must have disjoint affected relations. Maintain
degree updates and lifting exactly; measure setup, fill and total pipeline
cost. The 2021 sparse Gaussian elimination study motivates these comparisons.
[Author paper](https://perso.lip6.fr/Charles.Bouillaguet/static/publis/merge.pdf).

Dense original-relation masks and full-index parity rows should be compared
with an immutable merge history and compact live constraint numbering.
Reconstruct only emitted dependencies through the retained inverse maps.
CADO's merge/replay separation is a design reference for this lifecycle.
[Pinned merge](https://github.com/cado-nfs/cado-nfs/blob/692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b/filter/merge.cpp)
and [replay](https://github.com/cado-nfs/cado-nfs/blob/692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b/filter/replay.cpp).

The current worst-case reservation is deliberately conservative. A new
representation needs a bound for all simultaneously live rows, provenance,
cache entries and reconstruction scratch. Simply shrinking that estimate is
not a capacity optimization. Likewise, a verified immutable preparation
cache must bind the factor base, polynomial and atom contents; it cannot
remove independent validation of untrusted checkpoints.

Compare re-filter cadence against the current changed-store control using
post-filter excess and useful dependency progress. Incremental elimination
state is another candidate, with its own checkpoint and invalidation costs.
Preserve parity-identical but arithmetically distinct relations; exact payload
duplicates and equal parity are different cases.

The existing matrix report already gives the key solver safeguards. Retain
original-operator checks for Gram-kernel outputs, rank-aware dense recovery,
bounded Four Russians tables and a genuine matrix-polynomial block Wiedemann
generator. Dense PLE/Four Russians remains a relevant intermediate challenger;
its native cache-tuned cutoffs do not transfer automatically. [Dense GF(2)
elimination paper](https://arxiv.org/pdf/1111.6549),
[CADO block Wiedemann stages](https://github.com/cado-nfs/cado-nfs/blob/692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b/linalg/bwc/README).

### Large primes need useful-dependency diagnostics

Classical multiple-large-prime work and CADO filtering both support measuring
the graph, excess and fill together. Implement P5.4 only with per-prime bounds,
a separate residual-product cap and finite splitting work. Include loops,
repeated primes and disconnected cycles; preserve exact exponent corrections
and referenced atoms through eviction. [Large-prime paper](https://ir.cwi.nl/pub/1367/1367D.pdf),
[Cavallar filtering paper](https://ir.cwi.nl/pub/4456/04456D.pdf),
[YAFU filter](https://github.com/bbuhrow/yafu/blob/963dbe9c45283cc06e9a71d830b0676bc1b0d343/factor/qs/filter.c).

Add proper-divisor yield per verified independent dependency to collector
diagnostics. A relation cycle is not yet a kernel dependency; a kernel
dependency is not yet a proper factor. Track repeated trivial congruences
when raising large-prime limits, while preserving both GCD signs.

CUDA-MPQS's recent notes describe merge histories and a source-specific
large-prime/trivial-square failure pattern. That motivates a diagnostic
fixture, not a general explanation or compulsory extra QS character columns.
Some GPU debug checks are compiled out on normal paths, so adopting a layout
must preserve Factor's independent verification. Timings remain author
reports on different hardware/arithmetic. [Pinned matrix notes](https://github.com/drjanosch42/cuda-mpqs/blob/b1a9c4500da15a2bea56a59ecad7ef61cf054756/docs/modules/matrix.md),
[square-root notes](https://github.com/drjanosch42/cuda-mpqs/blob/b1a9c4500da15a2bea56a59ecad7ef61cf054756/docs/modules/sqrt.md)
and [GPU root code](https://github.com/drjanosch42/cuda-mpqs/blob/b1a9c4500da15a2bea56a59ecad7ef61cf054756/src/sqrt/sqrt_step.cu).

### SSS, array and worker comparisons remain separate

The SSS paper reports complete factoring for smaller bands, while its
75–100-digit experiments compare one-hour smooth-relation collection on a
small sample. Its advertised ratios therefore do not establish full-factor
leadership at 100 digits. SSSf deliberately discards candidates, and factor-
base cardinality is not the same parameter as Factor's prime bound. Keep the
full-output comparison and useful-row accounting. [SSS paper, v2](https://arxiv.org/html/2301.10529v2).

The author implementation uses a different dependency/output environment and
is not internally constrained by Factor's resource contracts. Compare it as
a clearly labeled reference arm, with any external stopping/adaptation
disclosed. The worktree adapter is independent: recover forced-prime
quotients and all exponents, charge smooth-batch construction and residual
work, and validate both in-memory resume and the newly added serialized
assignment/store/solver reconstruction. Explicit CLI/portfolio selection
does not establish automatic dispatch superiority. The README's older
no-serialized-checkpoint description needs reconciliation with current code.
[Author SSS code](https://github.com/sbaresearch/smoothsubsumsearch/blob/8dbaf6d39ab88a40380965d25ec2c363d7f27358/sssif/sss.py),
[author test driver](https://github.com/sbaresearch/smoothsubsumsearch/blob/8dbaf6d39ab88a40380965d25ec2c363d7f27358/test.py).

For the optional NumPy arm, duplicated fancy indices do not accumulate all
increments with ordinary buffered `a[indices] += value`; use a verified
aggregation or `numpy.add.at`. Prove intermediate score/index bounds and
keep exact polynomial/factoring arithmetic in Python integers. PyPy's C-API
compatibility does not itself establish a speedup. [NumPy repeated-index
semantics](https://numpy.org/doc/stable/reference/generated/numpy.ufunc.at.html),
[PyPy NumPy guidance](https://doc.pypy.org/faq.html#should-i-install-numpy-or-numpypy).

For workers, reuse P6.3's common accounting. Stable assignment-derived seeds
can reproduce the assigned search across worker counts; arrival order and
first-factor stopping still need explicit semantics. Charge setup, IPC,
wasted assignments, cancellation and aggregate memory before promotion.

## Phase 4: arithmetic and chains

Production ECM uses the ladder; the prior PRAC counterexamples remain valid.
Near-optimal precomputed Lucas chains are worth a bounded challenger, but
abstract additions/squarings do not predict PyPy loop cost. Validate chain
integer invariants, terminal conditions and each complete prime-power action.
Test projective validity independently: cross-product equality alone accepts
`(0,0)` vacuously. [Pinned Lucas chain generator](https://github.com/sethtroisi/gmp-ecm/blob/8ea5e214fdf2f0ddf9415141b8dc039ed6f5874e/LucasChainGenerator/README),
[ECM survey](https://members.loria.fr/PZimmermann/papers/ecm.pdf).

P4.2 should record formula conventions, not just operation counts. Factor's
doubling uses `a24=(A+2)/4` and the squared difference. A formula using
`(A-2)/4` has the corresponding different term. Fixed-difference normalization
and fused loops must retain nonunit handling and scalar/chunk invariants.
Compare complete stage-one/two and factoring costs, including setup.

P4.3 belongs before GNFS scaling, while small correctness work can proceed
with Python integers. The system PyPy has no gmpy2 installation; the
project-local PyPy venv successfully imports gmpy2 2.3.1 with GMP 6.3.0.
Both use PyPy 7.3.23 implementing Python 3.11.15. This pass changed neither
environment. A future experiment should pin the supported build, keep
long-lived mpz values in specialized coarse kernels,
and measure conversions and bitset work separately. Division, roots, inversion
and certainty need exact contracts: ordinary mpz division can produce a
floating result. The documented GIL-release option is experimental and needs
an actual operation/thread test. [gmpy2 integer API](https://gmpy2.readthedocs.io/en/latest/mpz.html),
[context/GIL option](https://gmpy2.readthedocs.io/en/latest/contexts.html).

P4.4's existing reducer loss remains evidence. Extend exact reducer controls
to encoded identity, parameter/coordinate conversion, boundary widths and
canonical exits; prove GCD-preserving unit scaling where used. Retain native
`%` until a complete-loop/full-run experiment overturns the prior decision.
No new reducer benchmark was run here.

## Phase 5: complementary methods and continuation

P5.1's proposed module is corrected to `williams_pp1.py` for the active naming
convention. Specify binary Lucas evaluation, composition and doubling with
finite parameter trials. Discriminant GCD checks and a genuine stage two are
essential; a Jacobi symbol modulo composite n does not determine all hidden
prime Legendre symbols. Preserve element/group-order and prime-power conditions,
including saturation recovery. [Pinned p+1 reference](https://github.com/sethtroisi/gmp-ecm/blob/8ea5e214fdf2f0ddf9415141b8dc039ed6f5874e/pp1.c).

P5.2 should compare paired classical continuation before more complex
polynomial methods. Use an independent eligible-prime oracle, projective
cross-products and complete tails; tune D from actual table occupancy and
total baby/giant/product/replay memory. Positive initialization and exceptional
distances remain explicit. [Pinned stage-two reference](https://github.com/sethtroisi/gmp-ecm/blob/8ea5e214fdf2f0ddf9415141b8dc039ed6f5874e/stage2.c).

P5.3 needs an exact incremental B1 contract. Raising the bound includes
increased powers of existing primes, not just new primes. From 8 to 16, the
stage-one scalar ratio is `2*3*11*13`; the extra 2 and 3 matter. Preserve the
starting state and method-specific action: ordinary powering, Lucas
composition and elliptic scalar multiplication are different recurrences.
[Pinned ECM continuation](https://github.com/sethtroisi/gmp-ecm/blob/8ea5e214fdf2f0ddf9415141b8dc039ed6f5874e/ecm.c),
[p−1 reference](https://github.com/sethtroisi/gmp-ecm/blob/8ea5e214fdf2f0ddf9415141b8dc039ed6f5874e/pm1.c).

P5.4 retains the DLP implementation owner; P3.8-R4 owns current Phase 3
comparison/diagnostics. Do not silently expand the existing residual
prime-certification domain. Bound any rho/ECM/batch splitting and preserve
all atom references across graph eviction. Higher large-prime yield must
survive splitting, filtering, provenance and extraction costs.

## Phase 6: expensive options and credible comparisons

For Edwards/windowed/torsion ECM, success depends on the selected point's
order as well as curve order. State each family's congruence and formula
conditions, test small point orders independently, and compare success per
total CPU-second across many curves/seeds. The Edwards author resources
make both stage costs and success probability relevant; a best operation
count or favorable single curve is insufficient. [EECM author resources](https://eecm.cr.yp.to/),
[performance methodology](https://eecm.cr.yp.to/performance.html).

For polynomial continuation, compare paired classical stage two first.
Product/remainder trees and multipoint evaluation require coefficient,
node and storage limits; division over a composite modulus must specify
monic/unit assumptions and factor recovery. Bigint packing needs carry
bounds. FFT acceleration requires an exact reconstruction contract. Triple-
large-prime relations need general incidence/provenance, rather than a
two-endpoint cycle model. [Brent–Kruppa–Zimmermann continuation chapter](https://members.loria.fr/PZimmermann/papers/Chap8.pdf).

P6.3 should preserve assignment identity and committed/in-flight restart
state across 1/2/4 workers, while disclosing race-dependent stopping versus
a deterministic merge mode. Global work, CPU, deadline, retained data and
spill limits include every worker; restart cannot reset allowances. Compare
fixed-work throughput and first-valid-factor latency separately. Thread
scaling requires demonstrated backend GIL release on the actual hot path.

P6.4's older Python competitor set is useful for that class but insufficient
for a broad state-of-the-art claim. Add feasible native reference arms for
YAFU, FLINT, yamaquasi, GMP-ECM and CADO; label historical msieve explicitly.
Pin source/build, architecture, cores, arithmetic and output contract. Keep
factor-one, full decomposition and certainty labels distinct. [FLINT author
implementation blog](https://wbhart.blogspot.com/2017/02/integer-factorisation-in-flint.html),
[current CADO stage overview](https://cado-nfs.gitlabpages.inria.fr/).

FLINT's pinned QS lifecycle can reset retained relations when growing its
base. That is a design difference to account for, rather than a direct
checkpoint-preservation template. Factor needs an explicit preserve/remap
or bounded restart policy with all prior work charged. [Pinned lifecycle](https://github.com/flintlib/flint/blob/00cc19b350b5302876b4fa8c877df5a676c2061c/src/qsieve/linalg.c).

## Phase 7: GNFS contracts before scale

### P7.1: nonmonic norms and ideal identity

Let f have degree d, leading coefficient f_d and root alpha. The homogeneous
integer value is `F(a,b)=b**d*f(a/b)=f_d*Norm(a-b*alpha)`. Confusing these
objects loses a leading-coefficient correction. A scaled monic generator,
such as beta=f_d*alpha, changes basis and denominator bookkeeping. Record
that transformation and its modular evaluation explicitly. [RSA-240/250
paper](https://arxiv.org/pdf/2006.06197),
[CADO nonmonic square-root handling](https://github.com/cado-nfs/cado-nfs/blob/692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b/sqrt/sqrt.cpp).

Ordinary degree-one ideals need side, prime and affine/projective root
identity. Bad primes can need branch-specific valuation data; factoring
the norm alone does not supply this. A small initial reference may restrict
its supported fields/primes, provided it rejects or retries finitely and
documents the restriction. Leading coefficients, discriminants, order index
and nonunit denominators need explicit treatment. [CADO ideal numbering](https://github.com/cado-nfs/cado-nfs/blob/692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b/dev_docs/README.renumber),
[bad-ideal code](https://github.com/cado-nfs/cado-nfs/blob/692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b/utils/badideals.cpp).

### P7.2: bounded full relations first

Retain primitive pairs, signs, full homogeneous values, side labels and
known forced factors in the exact verifier. Full relations and exhaustive
small-region controls establish the first collector. Versioned spills need
polynomial/ideal-numbering identity, bounded disk/I/O and incomplete-write
recovery. Later special-q collection must retain its forced ideal exponent
and exact mapping from lattice coordinates back to `(a,b)`.
[CADO lattice construction](https://github.com/cado-nfs/cado-nfs/blob/692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b/sieve/las-qlattice.cpp),
[live sieve implementation](https://github.com/cado-nfs/cado-nfs/blob/692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b/sieve/las.cpp).

### P7.3: character placement and kernel correction

CADO can reimpose characters and omitted heavy constraints through a small
solve inside the provisional kernel span. This is a concrete alternative
to including every such constraint in the large operator. Both require
exact original-relation lifting, nonzero/independent output checks and a
declared character policy. A finite character screen does not prove the
algebraic product is a square. Factoring characters are distinct from
discrete-log Schirokauer maps. [Pinned character correction](https://github.com/cado-nfs/cado-nfs/blob/692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b/linalg/characters.cpp).

Keep the mandatory P7.4 field identity even after ideal and character parity
passes. Preserve original heavy constraints and every relation-map version
in checkpoints if the main solver omits them.

### P7.4: choose a finite square-root contract

Inert-prime p-adic lifting requires f to remain irreducible modulo an
auxiliary prime. Some irreducible fields, including quartic examples with
Klein-four Galois group, have no inert prime. Specify supported fields,
bounded search and refusal/retry or a separately validated alternative.
CRT approaches also need consistent root-sign reconciliation.
[Thomé's square-root paper](https://members.loria.fr/EThome/files/nfs-sqrt.pdf).

An explicit control is `x**4-10*x**2+1`, the minimal polynomial of
`sqrt(2)+sqrt(3)`: its splitting field has group V4, whose elements have no
4-cycle. Its irreducibility over the rationals therefore does not ensure
irreducibility modulo an auxiliary prime. This algebraic obstruction cannot
be repaired merely by extending a prime-search timer.

CADO's pinned `FindSuitableModP` searches primes below 1,000,000 and reports
failure when none is suitable; this is a finite reference policy, not a
recommended PyPy bound. Its root code also transforms nonmonic polynomials
and tracks denominator powers. Factor needs justified coefficient bounds,
documented square corrections and exact reconstruction verification before
mapping roots modulo n. [Pinned square-root implementation](https://github.com/cado-nfs/cado-nfs/blob/692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b/sqrt/sqrt.cpp).

Controls should include a no-inert-prime field, nonmonic scaling, failed
denominators, corrupted character data, nonsquare algebraic products and
trivial dependencies. Compute both GCD signs and validate proper divisors.
A square norm is not an algebraic square-root algorithm. This specification
should precede enlarging polynomial degrees or input bands.

### P7.5: retain the whole field and budget identity

Checkpoints must bind polynomials, basis/normalization, norms, ideal numbering,
characters, lifting history, arithmetic backend and consumed stage allowances.
Shared SIQS storage interfaces are useful engineering infrastructure; their
payload does not encode GNFS ideals or field elements. Finish the bounded
small general-composite pipeline through square roots before claiming scale.

### P7.6: polynomial quality and two-sided cofactoring

Size, skew and root quality all matter. Shortlist with calibrated estimates,
including Murphy E and optionally E', then spend a finite trial-sieve budget
on promising candidates. Recheck common roots and supported fields after
rotations/translations. Charge selection cost to complete factoring;
a better score alone does not establish a faster run. [Root-optimization
paper](https://arxiv.org/pdf/1212.1958),
[CADO Murphy E](https://github.com/cado-nfs/cado-nfs/blob/692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b/polyselect/murphyE.cpp),
[E/E' implementation](https://github.com/cado-nfs/cado-nfs/blob/692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b/polyselect/E.sage).

The 2023 alternative-sieving study motivates separating medium-prime sieving,
batched small-prime removal, staged tests and residual ECM, with a measured
first-side decision. Its reported local advantages include intentional
relation loss; complete-pipeline oversieving and matrix costs still determine
promotion. Compare cofactor strategy jointly with bounds and scores, rather
than copying a standalone threshold. [Authors' 2024 seminar slides](https://caramba.loria.fr/sem-slides/202402151030.pdf),
[CADO cofactor strategy](https://github.com/cado-nfs/cado-nfs/blob/692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b/sieve/ecm/facul.cpp).

Native parameter files are experiment references, not supported PyPy
presets. Special-q lattices, two-sided large-prime identities and sparse
operators all retain exact verification and resource gates. GNFS scaling
can reuse P3.8's measured kernels without waiting for every optional solver.

### P7.7: capability and crossover are measurements

General GNFS inputs and SNFS-friendly forms must be separate cohorts. Freeze
training/held-out data and report censored failures, completion, CPU/RSS/disk
and complete-factor time. Measure the SIQS/GNFS crossover under the actual
resource model; neither asymptotic L-notation nor record hardware timings
fixes a digit cutoff for this implementation.

## Refresh findings and action ownership

### P3.8: extend search and test family-wide hit scheduling

An extendable search should use finite resident batches with a cumulative
caller allowance, stable assignment cursor and bounded duplicate tracking.
Increasing A's factor count changes polynomial quality and the exponential
Gray-family size; it should not be used merely to manufacture more search
steps. YAFU's pinned SIQS generates further A-polynomials while collection is
needed and can raise its relation target after unsuccessful postprocessing.
Its optional time/relation stops are separate controls. These are design
references for P3.8-R1, not recommended Factor parameter values.
[YAFU collection](https://github.com/bbuhrow/yafu/blob/963dbe9c45283cc06e9a71d830b0676bc1b0d343/factor/qs/SIQS.c),
[documented stop/resume controls](https://github.com/bbuhrow/yafu/blob/963dbe9c45283cc06e9a71d830b0676bc1b0d343/docfile.txt).

CADO's experimental SIQS provides a distinct R2 candidate: for factor-base
primes/prime powers larger than the sieve interval, enumerate hits across
a Gray-code family using sorted CRT half-sums and buckets. Its coprime fast
path assumes `gcd(p,A)=1`; the local per-polynomial marker remains the oracle
for exceptions, labels, tails and prime powers. Half-sum tables are exponential
in the selected split dimension, so construction and queued hits need explicit
bounds. Test this only when profiles identify relevant family/root/scan cost.
[Pinned derivation](https://github.com/cado-nfs/cado-nfs/blob/692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b/sieve/siqs-largesieve.hpp),
[bucket implementation](https://github.com/cado-nfs/cado-nfs/blob/692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b/sieve/siqs-fill-in-buckets.inl).

The source follows Kleinjung's 2016 *Quadratic Sieving* work. The fresh
review inspected its code derivation, not the complete paper (AMS PDF access
was blocked). Inria reports experimental SIQS additions in 2025 and continued
2026 development, particularly for quadratic-field class groups; that is
not a completed PyPy factoring comparison.
[Publication DOI](https://doi.org/10.1090/mcom/3058),
[Inria activity report](https://radar.inria.fr/report/2025/caramba/index.html).

CUDA-MPQS documents two distinct provenance designs: V1 CSR/merge-tree replay
and V2 packed full exponents with incremental modular square-root payloads.
R3 can assess these separately if extraction/storage dominates. Incremental
residues still need independently checkable original parity, exponent packing,
square corrections and corrupt-state rejection; they cannot certify their
own history. The existing CADO replay reference remains useful.
[CUDA design distinction](https://github.com/drjanosch42/cuda-mpqs/blob/b1a9c4500da15a2bea56a59ecad7ef61cf054756/docs/modules/matrix.md).

### P4–P6: stronger bounded arithmetic challengers

Bernstein–Cottaar–Lange's 2025 *Searching for differential addition chains*
adds stronger low-memory pruning and meet-in-the-middle continued-fraction
chain search. Guaranteed and heuristic variants have different search
contracts. Their chain-length/node results are not full ECM timings; some
author prototype search bounds use floats. P4.1 should compare independently
verified bounded/offline integer chain records with guarded PRAC and the
ladder, charging search, loading, cache, dispatch and exceptional recovery.
CADO's compact bytecode is an execution representation reference; its
interpreter assumes the input bytecode is correct.
[Paper](https://doi.org/10.1007/s40993-024-00604-8),
[author manuscript](https://cr.yp.to/papers/dacmitm-20240627.pdf),
[author code archive](https://cr.yp.to/2024/dacbench-20240609.tar.gz),
[CADO bytecode specification](https://github.com/cado-nfs/cado-nfs/blob/692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b/sieve/ecm/README.bytecode).

P5.2 gains two separate candidates. CADO plans coprime baby distances and
can reuse a term when `v*w +/- u` is also divisible by another eligible
prime; coverage must be verified prime by prime, including wheel exceptions
and tails. Its arrays indexed by absolute B2 should become segmented/bounded
plans in Python. No-inversion common-Z tables instead trade setup/storage
for fewer cross-product multiplications. Verify
`X'_i-X'_j=(X_i*Z_j-X_j*Z_i)*product(Z_k for k != i,j)` and retain denominator
GCDs and mixed-factor saturation replay: the extra scaling may be a nonunit.
Neither optimization justifies dropping existing recovery.
[Planner](https://github.com/cado-nfs/cado-nfs/blob/692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b/sieve/ecm/stage2.c),
[common-Z/extraction implementation](https://github.com/cado-nfs/cado-nfs/blob/692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b/sieve/ecm/ecm.cpp).

After P5.1's binary Lucas correctness, rational starts such as 2/7 and 6/5
are finite parameter arms, with denominator/discriminant GCD recovery before
inversion. Their conditional order properties do not make repeated p−1
bases independent ECM-like trials.
[CADO p+1 implementation](https://github.com/cado-nfs/cado-nfs/blob/692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b/sieve/ecm/pp1.cpp).

For P6.1, CADO MISHMASH specifies a concrete mixed-coordinate package:
Edwards signed/double-base/precomputed blocks followed by differential
Montgomery blocks, with a Montgomery stage-two exit. Validate chains, maps,
coordinate tags and low point orders before comparing complete curve cost.
Farashahi–Fadavi–Sabbaghian's TCHES 2024 complete Montgomery laws are another
screen, but require finite-field congruence/character hypotheses and use
full coordinates. Composite-n Jacobi symbols do not establish all those
hypotheses; the paper does not supply a measured x-only ECM replacement.
[Mixed-chain specification](https://github.com/cado-nfs/cado-nfs/blob/692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b/sieve/ecm/README.bytecode),
[complete-law paper](https://doi.org/10.46586/tches.v2024.i4.737-762).

The current `parallel_candidates.py` gives each assignment 250,000 work
units with CPU/time caps disabled, then evaluates CPU acceptance afterward.
Its RSS observer samples and requests cooperative cancellation. This is
fixed-work feasibility evidence, not an enforced global production budget.
P6.3 needs parent-owned work leases, reconciliation of reservations, shared
deadlines, live/exited-worker CPU accounting and bounded queued/duplicated
storage. Copying `Budget` cannot aggregate worker CPU: `process_time()` is
per-process. Running futures also need an explicit cooperative stop protocol.
[Python clock contract](https://docs.python.org/3.11/library/time.html#time.process_time),
[executor cancellation semantics](https://docs.python.org/3.11/library/concurrent.futures.html#concurrent.futures.Executor.shutdown).

### P7: feedback, deduplication and supported root algorithms

CADO's `PurgeTask.request_more_relations` responds to useful post-dedup/filter
excess. A raw relation target triggers filtering; if the matrix is deficient,
collection resumes. Its workunit timeout cancels/resubmits an assignment,
rather than imposing a whole-job factoring deadline. P7.2–P7.5 should retain
this feedback under explicitly authorized finite assignments and distinguish
immutable field/store identity from extendable quotas and assigned ranges.
Changing bases/numbering/cofactor policies still needs preserve/remap/reverify
or charged restart. Native import checks sample relations; Factor retains
full independent verification. Numerical preset ratios are not PyPy defaults.
[Task feedback](https://github.com/cado-nfs/cado-nfs/blob/692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b/scripts/cadofactor/cadotask.py),
[continuation/import documentation](https://github.com/cado-nfs/cado-nfs/blob/692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b/scripts/cadofactor/README.md),
[workunit parameter semantics](https://github.com/cado-nfs/cado-nfs/blob/692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b/parameters/factor/params.c90).

Canonical `(a,b)` identity, sign and b=0 policy should precede exact dedup
of overlapping/retried discoveries. Optional online suppression is narrower:
a smaller special-q dividing a norm does not prove that an earlier assignment
covered the relation under its geometry and cofactor policy. CADO qualifies
its probabilistic rule and both-side limitations. Introduce prime special-q
before composite special-q, with exact determinant/congruence and inverse-
coordinate checks; native width/skew assumptions are not our proof.
[Duplicate checks](https://github.com/cado-nfs/cado-nfs/blob/692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b/sieve/las-duplicate.cpp),
[lattice construction](https://github.com/cado-nfs/cado-nfs/blob/692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b/sieve/las-qlattice.cpp).

GNFS relations can carry more than two retained large ideals across the two
sides. P7.6 therefore needs general side-labelled sparse incidence, unless
a restricted QS-edge reuse is proved. Separate each prime's lpb bound from
the whole residual mfb bound and effort on each side. Weak rejection tests
in upstream cofactor screening are not accepted-prime certification. Test
three-plus large ideals, powers, equal primes on different sides and distinct
roots over one prime against exact elimination and full norm reconstruction.
[Cofactor patterns](https://github.com/cado-nfs/cado-nfs/blob/692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b/sieve/las-cofactor.cpp),
[per-side configuration](https://github.com/cado-nfs/cado-nfs/blob/692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b/sieve/las-side-config.cpp).

CADO's separate CRT algebraic-root program has manual rational-root/GCD
integration and verification TODOs, plus an open-ended split-prime search.
It is not an already integrated general fallback. P7.4 must bound auxiliary
search, precision growth and sign reconciliation, then certify the exact
field identity. A heuristic coefficient estimate needs checked reconstruction
with finite growth, not promotion to a proved bound.
[CRT implementation and limitations](https://github.com/cado-nfs/cado-nfs/blob/692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b/sqrt/crtalgsqrt.cpp).

## Recent research screened for transfer

The September 2026 RSA-260 author report describes a completed GPU GNFS run
with approximately 4,900 GPU-days and stage engineering across the pipeline.
It does not claim a new factoring algorithm. The post is dated 9 September
and gives a run completion date of 3 September. Its resource/timing claims
were not independently reproduced here, and no public implementation link
was located in that post. The transferable lesson is complete-stage design
and accounting, not a PyPy performance forecast or project cost estimate.
[Author report](https://cognition.com/blog/factoring-rsa-260).

Guan's QS paper appeared on arXiv in September 2026 but carries a 2025 journal
DOI. Its four-hit marking unroll remains linear sieve work, and trivial
dependencies still exist. Treat it as a constant-factor loop experiment
only if profiling justifies it; retain exact both-sign GCD validation.
[Paper](https://arxiv.org/pdf/2609.06576v1).

Oznovich–Volk's 2025 result weakens a high-order target condition; Harvey–
Hittmeir's 2026 work removes it. These are relevant deterministic algorithmic
advances, but their abstracts do not establish a practical PyPy replacement
for heuristic ECM/SIQS/GNFS. Record them as screened theory and require a
relevant bounded implementation experiment before adding a core workstream.
[2025 paper](https://arxiv.org/abs/2506.07668v3),
[2026 paper](https://arxiv.org/abs/2601.11131v2).

Bernard–Fouque–Lesavourey's 2023 root paper concerns odd prime-power e and
generalizes methods originating in NFS square roots. It is not a drop-in
e=2 upgrade. CADO's `sqrt/montgomery.md` discusses related discrete-log root
work; it is not proof that a new general factoring square-root route is the
default. Keep this separate from the required P7.4 implementation.
[Root-paper abstract](https://arxiv.org/abs/2305.17425v2),
[CADO experimental root notes](https://github.com/cado-nfs/cado-nfs/blob/692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b/sqrt/montgomery.md).

GPU ECM and CUDA-MPQS are useful architecture references, with fixed-width,
native/GPU and verification differences. Neither supplies PyPy performance
evidence or an automatic new backend requirement. [GPU ECM paper](https://eprint.iacr.org/2020/1265.pdf),
[pinned CUDA-MPQS overview](https://github.com/drjanosch42/cuda-mpqs/blob/b1a9c4500da15a2bea56a59ecad7ef61cf054756/README.md).

## Suggested experiment order and evidence

1. Freeze a feasible SIQS control and untouched confirmation cohort under
   P3.8-R1. Identify dominant costs with separately instrumented profiles.
2. Compare candidate/power/resieve/batch variants and preparation/filter/
   provenance representations under R2–R3. Evaluate one causal change at a
   time, with total factorization and retained-resource accounting.
3. Specify P7.1–P7.4 field, ideal and finite square-root contracts while
   assessing P4.3 arithmetic availability. Build and integrate the bounded
   small GNFS reference; optional Phase 3 optimizations do not block it.
4. Promote DLP, dense/hybrid or iterative solving only when useful yield,
   solve cost or capacity justifies the challenger. Keep SSS/NumPy/workers
   separate, with explicit adopt/defer/reject outcomes.
5. Scale GNFS selection/collection/cofactoring and measure the crossover;
   then use Phase 8 for integrated Phase 2 portfolio tuning. Useful independent
   Phase 8 experiments may still proceed earlier.

Each promoted experiment needs exact reconstruction/certainty checks,
prespecified inputs/seeds/allowances, separate cold and warmed results,
at least three seconds validated warmup and nine samples, extensions for
unstable timing, and the roadmap's unchanged end-to-end promotion policy.
Include preparation, conversion, recovery, lifting, I/O, failed attempts and
finite refusals. Native results must disclose backend/build/cores separately.
Publish a concise summary with commands and limitations; retain useful raw
evidence locally or in durable external storage.

Before publication, every immutable corpus/baseline required by a loader must
be explicitly retained, and a committed-files-only checkout must run the
tests and import the benchmark runners. The present worktree contains
uncommitted implementation and generated evidence; this research pass does
not stage, commit or publish it.

## Budget follow-up: large-job capacity

The user's concern about larger jobs was checked on 4 October 2026 against
the current code, saved P3.5 captures and deterministic PyPy probes. These
checks establish restrictions and refusal behavior; they do not establish
time-to-factor after relaxing those restrictions. The action owner is P3.8-R1.

Finite controls are common in mature implementations: ECM curve/bound targets,
SIQS time/relation stops, memory limits and saved progress. There is no
mathematical requirement for Factor's particular 2,000,000-unit default or a
per-hit clock/type check. The concern is supported by numerical allowances,
hot-loop accounting and hard parameter ceilings, rather than by the existence
of cancellation/resume controls themselves. Preserve finite caller-selected
campaigns and amortized checks while fixing those restrictions.
[YAFU user controls](https://github.com/bbuhrow/yafu/blob/963dbe9c45283cc06e9a71d830b0676bc1b0d343/docfile.txt),
[native ECM allocation](https://github.com/bbuhrow/yafu/blob/963dbe9c45283cc06e9a71d830b0676bc1b0d343/factor/autofactor.c).

- `Budget()` defaults to 2,000,000 work units and 30 seconds each of wall and
  CPU time. Time caps can be disabled individually, and work can be increased,
  but this does not enlarge a method's finite assignment or storage limits.
- `SIQSConfig()` has 16 families and three A factors, hence at most 64
  polynomial steps and 65,600 sieve positions across their 1,025-position
  windows. It can stop earlier after 16 consecutive windows with no new full
  or combined rows. Default width growth is disabled. The default factor-base
  bound 1,000 also limits A to less than 1,000 cubed; for a 30-digit magnitude
  such as `10**29 + 1`, the configured width's A target is 873,464,053,710.
  Failure to approach the target is a parameter-quality issue, not proof
  that a polynomial cannot produce a factor.
- Separate implementation ceilings reject base bounds above 100,000,
  more than 64 families, more than 1 GiB SIQS workspace or more than
  16 MiB checkpoint capacity, even with a much larger global work/time
  allowance. Full-job resume rejects changed configuration, so increasing
  a global timer cannot extend an exhausted family schedule or enlarge its
  retained store in place under the current contract.
- A reproduced 2,000-row/500-column matrix with one set bit per row is
  refused by the default work budget before filtering begins. The first
  round reserves `2000 * (500 + 2000 + 1) = 5,002,000` units. With 10,000,000
  units, the identical matrix is accepted. This is a conservative accounting
  barrier, not a measured expensive round or a complete factoring experiment.
- The present worst-case matrix reservation for 32,768 rows and 5,000
  columns is 1,198.7 MiB, before collector and other job storage. A sparse
  one-bit-per-row probe is therefore refused even with 1 GiB available to
  filtering. Sparse input still carries dense lifting/fill reservations;
  a safe representation improvement is preferable to weakening the estimate.

Recorded evidence supports the concern independently of these probes.
`p35_30d_siqs_trained.json` contains 72 case/seed/sample outcomes: 54 wall
limits and 18 stalled-yield stops. In `p35_50d_diagnostic.json`, all eight
SIQS cases stop on stalled yield before their five-second ceilings. In the
60-digit diagnostic, each SSS arm has two memory refusals and six wall stops;
all eight SIQS cases reach the wall limit. These are historical source-pinned
captures, with the larger bands only one cohort; concurrent code changes mean
they must not be represented as fresh measurements of the current worktree.

The newer `phase_three_large` runner already provides much larger explicit
allowances (10**13 work, 1 GiB portfolio/768 MiB SIQS workspace, 30/60/300-second
bands and a separate cumulative continuation experiment). That is a useful
evaluation path, but it does not remove implementation ceilings or establish
successful larger-band coverage merely by existing. Retain finite controls,
make allowances appropriate to the declared workload, and distinguish time,
work, storage and assignment stops before interpreting algorithm capability.

Probe code, source/capture hashes and results are retained locally in
`phase_3_plus_20261004/budget_concern/`. No runtime default or algorithm was
changed by this follow-up.

Historical follow-up validation: the deterministic probes and lint passed.
That snapshot's `make -C v2 test` ran 205 tests with one
failure: `test_pause_and_budget_extension_match_uninterrupted_factor` in
`test_sss.py` expects replacing a job's budget to raise `ValueError`, but the
then-current SSS implementation did not raise. Concurrent SSS work has since
restored the original-budget guard. The refresh's initial checkpoint-refusal
test failure was also resolved by concurrent work; the final checks below
supersede both failures. Historical logs remain with the local probes.

## Verification of this pass

The original review progressed from 201 to 204 passing tests; its logs remain
historical. After this refresh and the concurrent SSS/filter work,
`make -C v2 test` passed all **219 tests** on PyPy 7.3.23 implementing
Python 3.11.15, and `make -C v2 lint` passed Ruff/format/pycodestyle checks
(72 formatted files). This research pass edited no runtime or test file.
These checks validate the inspected worktree, not new algorithm acceptance
or performance. The new cProfile runs are explicitly instrumented
diagnostics; no repeated speedup comparison or default promotion occurred.

The local refresh evidence retains source hashes, 40/80-digit profile captures,
exact resource/reachability probes, selected upstream manifests and the frozen
Python/corpus snapshot needed to reproduce diagnostics despite concurrent
changes. Previously collected 50-digit traces remain historical exploratory
evidence only. The new benchmark design covers the requested range and
excludes a standalone 50-digit cohort.

Only `TODOS.md`, this report and the report's explicit `.gitignore` retention
entry were changed among public project files by the research pass. Earlier
Phase 1/2 and P3.1–P3.7 sections, runtime/tests, `v1/`, historical research,
accepted behavior and default selection were left untouched by this pass.
