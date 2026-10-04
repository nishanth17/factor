# Factor: phased TODOs and decision gates

Generated run captures and detailed milestone journals are local records,
excluded from GitHub. Historical evidence references below appear as
plain labels; retained inputs and public research remain linked. See the
[benchmark guide](../benchmarks/README.md) for rerun commands and
[public changelog](../../CHANGELOG.md) for accepted behavior.

Date: 3 October 2026  
Source: codebase audit, reviewed commit `1272e033f889a105792cbb924bf8a12a46ac88ae`.

Research reconciliation: 4 October 2026. The
[Phases 3+ research pass](phases_three_plus_research.md) compares the active
worktree, accepted evidence and current primary literature/source code.
All new Phase 3 research findings live in P3.8; this pass edits none of
P3.1–P3.7's milestone records.
P3.4's bounded implementation and declared large-number experiment gate
are complete at M31; broader scaling/default promotion stays separate.
Phase 8 remains the M23 Phase 2 optimization follow-up. Research additions below close no implementation or experiment gate.

Start with correctness repairs in Phase 1. Establish bounded, reproducible
execution in Phase 2, add SIQS coverage in Phase 3, then build GNFS as a
committed workstream. Keep the factoring algorithms in Python on PyPy Python
3.11; compare built-in integers and optional gmpy2 arithmetic on that runtime
as separate performance tracks. Existing phase/task IDs stay stable; the
execution order below governs dependencies rather than numerical order.

This is an implementation backlog; checked items have acceptance evidence below. Findings and existing counterexamples come from the audit; task boundaries, numerical promotion thresholds, and experiment designs below are recommendations. Report statements are evidence to evaluate, not instructions to execute. This document alone does not authorize implementation or benchmarks; Phase 1 was subsequently implemented at your explicit request.


Implementation target: PyPy implementing Python 3.11 only in `v2/`. CPython
support and routine compatibility runs are retired; their historical results
remain evidence. Retain `v1/` and historical audit captures for provenance.
Large result captures are losslessly compressed; original paths and hashes
are recorded in the archive manifest.
Source filenames in the original report refer to v1;
the production modules now use snake_case. Maintain milestone evidence and
measured improvements/regressions in [public changelog](../../CHANGELOG.md).

Your C sieve versions were inspected at commit
`5b4afb8f344ad5f6fbd20a186cd57b36182ea710`. The
[C-to-Python transfer review](sieve_port_review.md) maps useful ideas to the
phases below. C thresholds and OpenMP speedups are reference evidence, not
Python defaults or measured Python improvements.

## How to use the gates

- **Acceptance gate (A):** the observable behavior required to close a TODO. Every returned split must satisfy `1 < g < n` and `n % g == 0`.
- **Experiment gate (E):** the comparison or adversarial experiment required before integration or changing a default. For correctness repairs, this means passing regression and fault-injection cases; speed is not a prerequisite for fixing incorrect behavior.
- **Phase exit:** the condition for moving the dependent production work into the next phase. Independent research spikes may proceed once their stated prerequisites exist.
- **Evidence:** store the implementation commit, command, corpus identifier, seeds, environment, raw results, and decision in a separate artifact for each task. Preserve the original audit JSON as historical evidence.

Performance promotion policy proposed for Phases 3–8: zero correctness failures; obey the same time, CPU, and memory limits; show either at least 10% lower end-to-end median time on a prespecified comparable cohort or at least 10 percentage points higher completion within budget; allow no more than a 5-percentage-point completion regression in another declared workload class. Check uncertainty with repeated seeds and a confidence interval for the relevant difference. If evidence is inconclusive, retain the baseline and expand the sample. Include setup, conversion, timeout, and recovery costs. These thresholds are project policy suggestions, not measured predictions.

## Execution order toward GNFS

**Decision (M18): GNFS is committed roadmap scope, not a conditional future
project.** Its implementation is not started. Sequence by these prerequisites:

1. P3.4's SIQS implementation, complete extraction/checkpoints and optional
   dispatch are tested. M31 completes the declared longer 30–80-digit and
   varied-input evaluation with one fresh 50-digit success and explicit
   negative outcomes; retain ECM defaults. Broader scaling remains open.
   P3.1–P3.3's exact relation/filter/dependency controls and their M28 repairs
   provide shared infrastructure and a comparison baseline.
2. Bring P4.3 forward: assess a coarse arithmetic boundary and GMP/`mpz` on
   PyPy Python 3.11 before GNFS scaling. Keep conversions outside hot loops
   and preserve canonical checkpoints. A small GNFS reference can start with
   Python integers while this proceeds; GMP availability or speed is not a
   substitute for algorithm correctness and does not require CPython support.
3. Build P7.1–P7.5: a bounded GNFS pipeline that completes small general
   composites through both rational and algebraic square roots. Start its
   contracts after P3.3; integrate after the SIQS baseline is working. Shared
   relation storage/filtering interfaces may be reused, but SIQS relations and
   GNFS ideal/character data have distinct mathematical contracts.
4. Scale P7.6–P7.7: improve polynomial selection, lattice sieving and sparse
   linear algebra, then measure a SIQS/GNFS dispatch crossover on held-out
   inputs. Define budgets and feasible bands before runs; no fixed digit
   cutoff or performance date is promised before measurement.

P3.6.1 is the immediate follow-up for the diagnosed P3.5/P3.6 costs; it can
start before P3.8 or P6.3. Further P3.5 SSS comparison, P3.6 parallel promotion,
P3.7 NumPy, and P4–P6 optional optimizations can proceed when useful; completing
every such experiment is not a prerequisite for GNFS. M23 consolidates
remaining P2 optimization in future Phase 8, after working SIQS and bounded
small GNFS, with individual
spikes eligible earlier when profiling justifies them. Small reference
correctness, scalable execution, and performance-based default promotion are
separate milestones. The eight-phase catalogue preserves existing
task/evidence IDs; Phase 8 is the user-requested future follow-up.

## Phase 1 — Immediate correctness and API repairs

**Status (2026-10-03): complete; Phase 2 now has separate evidence.** All nine
items passed the M10 reverification, then covered by 66 native tests. Those
tests passed
on PyPy 7.3.23 (Python 3.11.15), CPython 3.9.6, and CPython 3.14.7,
including the 1–20,000 reconstruction sweep and independent arithmetic/sieve
oracles. See the fresh PyPy evidence,
Python 3.9 evidence,
Python 3.14 evidence, and
[public changelog](../../CHANGELOG.md) for seeded benchmarks, regressions, limits,
and provenance. Subsequent Phase 2 implementation changed source hashes and
expanded the suite to 94 tests. All pass on the same three interpreters in
the M12 captures: PyPy,
Python 3.9, and
Python 3.14. Historical M9/M10
captures describe their measured versions. The later phase gates remain open;
preliminary parallel probes do not enable production concurrency.

**Prerequisite:** none. **Goal:** trust results before measuring speed. Source: correctness findings.

### P1.1 — Port the seven modules to explicit Python 3 arithmetic

- [x] Update `factor.py`, `ecm.py`, `pollard_rho.py`, `pollard_pm1.py`, `prime_sieve.py`, `utils.py`, and `constants.py`; replace obsolete syntax/APIs and document the supported Python versions. Classify each division as exact integer division, index arithmetic, or modular inversion. Keep the audit compatibility loader as an independent historical reference.
- **A:** all modules import natively; integer operations retain intended semantics; no factoring path converts a large integer through a float. Exact prime powers use multiplication and integer bounds, including `B1=243, p=3 → 243`.
- **E:** run native regression tests against the audit fixtures and the original 1–20,000 reconstruction sweep. Record intentional differences caused by repairs; do not blindly replace every `/` with `//`.

### P1.2 — Repair Suyama ECM initialization

- [x] Replace `ecm.py:287–295` rational setup with modular inversion and projective `P=(u³ mod n, v³ mod n)`. GCD-check `16*u³*v` and `A²−4`; return a proper factor or explicitly retry a singular/fully degenerate curve.
- **A:** `n=101, sigma=6` produces `a24=93` and affine `x=78`. Noninvertible denominators cannot raise an uncaught inversion error or produce an invalid accepted curve.
- **E:** compare setup and scalar results with an independent prime-field oracle; construct composite-modulus cases covering denominator GCDs `1`, a proper divisor, and `n`. Retain the independently checked ladder as the baseline.

### P1.3 — Repair and validate every sieve endpoint

- [x] Fix Atkin emission from `60*k+d` to `k+d`; adopt one convention, proposed `[lo, hi)`, and migrate all callers. Include base primes through `isqrt(hi-1)`. Replace estimated fixed output capacities with safe growth; handle empty/tiny ranges and remove shared scratch state.
- **A:** full prime sequences match an independent reference, including direct Atkin(100), the 3,500,000 dispatcher transition, and 59/60/61/67. No composite, duplicate, or out-of-range value appears. Tiny ranges no longer crash; `[3700,3722)` excludes `3721=61²` by composite marking, not merely by the upper endpoint.
- **E:** exhaust small intervals and sample ranges around prime squares and segment/dispatcher boundaries. Compare values, not just counts. Keep Atkin bypassed until it passes the same gate as Eratosthenes.

### P1.4 — Preserve unresolved cofactors and define input behavior

- [x] Replace ambiguous list/`−1` propagation in `factor.py:97–126` with a result carrying factor multiplicities, remaining cofactors, completion status, and certainty labels. Validate every recursive split. Reject zero; define one; either represent the negative sign explicitly or reject negatives consistently. Make formatting safe for failure/partial/empty results.
- **A:** `sign * product(factors) * product(remaining) == original_n` for every result. Completion means no unresolved composite remains; it does not silently certify probable primes. A failed rho call reaches the configured fallback, and a failed child retains its cofactor.
- **E:** force rho/ECM failure for `25013*25031`, `1000000000039*1000000000061`, and inputs with already-stripped factors. Inject invalid split values `1`, `n`, and a nondivisor. Exercise `0`, `1`, `−15`, and repeated-factor recursion.

### P1.5 — Repair p−1 batching and saturation handling

- [x] In `pollardPm1.py:75–102`, multiply every stage-2 term into a batch product, including the first prime. GCD-check each batch and the tail; retain enough state to split/replay saturated batches. On stage-1 GCD `n`, recover or retry rather than continue with the saturated residue. Accept explicit valid bounds instead of the small-input logarithmic formula.
- **A:** `n=607*1019, B1=10, B2=200` finds a proper divisor; `factorize_pm1(9)` returns a valid result/failure without a logarithm error. No eligible stage-2 relation is silently skipped.
- **E:** construct stage-1-only, stage-2-only, first-prime, tail-batch, and mixed-factor saturation cases. Verify each accumulated product against an independent per-prime calculation.

### P1.6 — Repair rho retries and cap recovery

- [x] Move failure return outside the offset loop in `pollard_rho.py`; give each attempt an explicit operation allowance. Use a fixed configurable GCD batch size and bounded recovery, with independent walk seeds/polynomials.
- **A:** the controlled `n=35, y=4, m=34` failure does not terminate the entire retry schedule. Exhaustion yields explicit failure after exactly the allowed work; no recovery scan can run indefinitely.
- **E:** force a failed first attempt followed by a successful attempt; test full exhaustion and saturated products. Log evaluations, GCDs, and retries. Leave batch-size optimization for Phase 2.

### P1.7 — Repair ECM curve limits and failure transitions

- [x] Replace the `<= MAX_CURVES_ECM` loop with an exact attempt count. Return a validated factor immediately, including on the final allowed curve. Skip stage 2 after fully saturated stage 1; implement local stage-2 replay or an explicit bounded retry outcome.
- **A:** budgets of zero, one, and two make exactly the permitted attempts. A factor found in stage 2 on the final curve survives. Saturation never becomes apparent success.
- **E:** inject stage-1 GCD `n`, stage-2 GCD `n`, setup failure, and final-curve success. Assert both outcomes and stage invocation counts.

### P1.8 — Clarify primality and helper contracts

- [x] Repair `utils.py` public entry points: defined large-input behavior, explicit composite/probable/proven classification, honored random-round counts, and documented deterministic ranges. Replace fragile binary search with `bisect`; reconcile `xgcd` documentation and all callers. Replace floating trial-division roots with exact arithmetic.
- **A:** `is_prime_fast(10**36+1)` does not crash; small results do not leak truthy numeric residues such as `6` for 29. Larger probable primes are labeled accordingly. `factorize_bf(10**400)`, empty searches, and duplicate boundaries work under documented contracts.
- **E:** use primes, known pseudoprimes, threshold neighbors, requested-round spies, and Bézout/inverse checks. Test the Atkin caller alongside any `xgcd` convention change. Treat the audit's unverified primality guarantees as unverified, not demonstrated counterexamples.

### P1.9 — Quarantine unsafe PRAC and make diagnostics reliable

- [x] Keep the ladder as default; guard unsupported PRAC scalars with a fallback or explicit error until Phase 4. Gate output behind verbosity and remove misleading performance claims from `README.md`. Convert the audit counterexamples into assertions for the repaired native modules.
- **A:** public scalar calls cannot hang on 2/4/6 or silently return the wrong multiple for 9. Quiet runs produce no unsolicited ECM output. Documentation explains partial results and primality certainty.
- **E:** time-box scalar-domain regressions and confirm the regression suite distinguishes `(0,0)` from a valid projective point. Preserve compatibility timings as historical diagnostics, not a native baseline.

**Phase 1 exit:** all nine acceptance/experiment gates pass; every complete factorization reconstructs its input with correct multiplicities and explicit certainty. Capture a corrected, native, seeded baseline. Do not compare later implementations with fast-but-empty failed results from the original code.

## Phase 2 — Bounded scheduling, reusable infrastructure, and benchmarks

**Prerequisite:** Phase 1. **Goal:** controlled time/memory and reproducible algorithm selection. Sources: ECM schedules, parameters, portfolio.

**Implementation status (M17, 2026-10-03): core exit passed; baseline frozen;
remaining tuning gates open.** The bounded API is implemented and
opt-in. `preprocessing.py`, `budget.py`, `portfolio.py`, `schedules.py`, and
`stage_jobs.py` cover exact powers, shared budgets, seeded dispatch, streaming,
chunk recovery, and JSON checkpoints. At M17 the PyPy suite passed 97 tests;
see acceptance evidence.
M20 later recorded 120 tests with P3.1 accepted; its
QS evidence is separate from Phase 2 tuning.

M13's arithmetic loops retain M12 candidate assignments, work accounting,
and checkpoint boundaries. M17 confirms them on a fresh, independently
certified corpus: 120 unique inputs, with no overlap with the old corpus,
twenty held-out inputs per selected band, five seeds, and nine repetitions.
Both complete and factor-one modes are measured, including cold startup.
See fresh confirmation and
remaining small-band controls.
The original 882-input corpus remains unchanged. Core acceptance and an
unchanged baseline do not imply that every optional tuning gate is complete.

- [x] Exact arithmetic, reconstruction, certainty, and independent oracles.
- [x] Shared budgets, explicit exhaustion, cancellation, and resume identity.
- [x] Bounded schedules, allocation limits, and corrupt-checkpoint rejection.
- [x] Fresh confirmation of the M13 loops without retuning parameters.
- [x] Freeze baseline configuration/source hashes and document larger limits.

The [frozen baseline](phase_two_m17_frozen_baseline.json) records production
and evaluation configurations, source hashes, exact commands, and decisions.
The large-band screen is one warm
repetition under five seeds: balanced 30-digit complete factoring succeeds
4/100 times and balanced 40–80-digit inputs have no completions under the
50 ms operation caps. Both 60/100-digit classes with five-digit small factors
complete 100/100 times. This is capped feasibility, not universal limits or
a performance ranking. No reconstruction or measured RSS gate failed across
the 32,470 final confirmation/control/screen samples. Broad repeated tuning
and external competitor evaluation remain open.

**Retained choices:** trial cutoff 25,000; rho batch 64; chunk size 16;
ECM B1/B2 2,000/147,396; cache and rolling off; serial execution; current
wheel-6/bytearray kernels. M17 evaluates two ECM curves under fixed limits;
the library's existing 32-curve default is unchanged. The bounded API stays
opt-in. No new parameter optimum, broad speed ranking, or competitor victory
is claimed. The original full items below remain unchecked wherever their
declared experiment gate still needs work.

**Research reconciliation (M23):** the
[Phase 2 optimization report](phase_two_optimization_research.md) reviews
pinned GMP-ECM, YAFU, FLINT, SymPy and primefac, Brent/Bernstein/ECM papers,
author errata and implementer blogs. The core remains accepted; M13 loops
and M14 parallel feasibility are not new TODOs. Remaining execution is
consolidated in [Phase 8](#phase-8--revisit-phase-2-portfolio-optimization),
without changing the historical P2 checkboxes or closing experiment gates:

| Historical gate | Future execution owner | Reconciliation |
| --- | --- | --- |
| P2.1 | P8.2 | Exact preprocessing exists; trial cutoff and proof-backed pruning remain experiments |
| P2.2 | P8.3–P8.4 | Bounded Brent/dispatcher exist; tune retries/batches and joint ECM allocation |
| P2.3 | P8.6 | Streaming/capped cache exist; setup and reuse promotion need changed workload evidence |
| P2.4 | P8.5 | Recovery/checkpoints exist; measure distinct GCD, polling, serialization and durable-write costs |
| P2.5 | P8.1; publication P6.4 | Corpus/runner exist; broad repeated and actual competitor comparisons remain open |
| P2.6–P2.7 | Conditional P8.6 | Keep rolling/cache off and wheel-6/bytearray; reuse validated arms only for a demonstrated bottleneck |
| P2.8 | Promotion P6.3; integration P8.7 | Feasibility passed; serial retained, first-factor benefit still required |

### P2.1 — Add exact preprocessing and classification reuse

M9 implements exact square splitting and per-factorization classification
reuse, with matched small-workload measurements in
the regression report. The bounded API now also
implements higher powers, power-of-two stripping, and bounded optional Fermat.
The retained cutoff still needs the declared training/held-out comparisons;
the entire P2.1 experiment gate is not yet closed.
M12 tried all four cutoffs on training data. M17 retains 25,000; it does not
use the fresh held-out confirmation to retune this cutoff.

- [ ] Strip powers of two efficiently; detect squares and higher perfect powers with exact root verification; reconstruct multiplicities after factoring the base. Cache classifications within one factorization. Add only a bounded optional Fermat close-factor path.
- **A:** primes, powers, mixed repeated factors, and values immediately beside powers remain distinct; output invariants hold for arbitrary-size integers. Fermat obeys its budget.
- **E:** compare preprocessing cost and saved work on power-rich and random corpora; try trial cutoffs 1,000/5,000/25,000/100,000. Promote a cutoff using full-run results rather than trial-division throughput alone.

### P2.2 — Build a budgeted portfolio dispatcher

The serial dispatcher and shared cancellation/work/wall/CPU limits are present,
including rho on large cofactors and explicit unresolved results. Candidate
assignments and consumed work survive resume. Tier/batch experiments remain
separate from production default selection.
M17 retains rho batch 64 and the existing ECM tier. Larger tier calibration
and matched held-out parameter selection remain open.

- [ ] Schedule trial division → exact powers/classification → short rho and p−1 → ECM tiers → later SIQS. Allow bounded rho on large cofactors. Inject an RNG/seed; carry one global deadline/work allowance across recursion and retries. Record stage outcomes and resume state. While SIQS is absent, terminate with the unresolved cofactor.
- **A:** every stage and child consumes the shared budget; cancellation returns reconstructible partial results. Remove the >80-digit jump to `B1=430,000,000`. Log bounds, seeds/sigma, and work.
- **E:** compare rho batches 32/64/128/256 and exploratory ECM `B1` tiers 2,000/11,000/50,000; fit `B2` and curve counts on training inputs. Test expired budgets, cancellation, and resumed results. These grids are candidates, not defaults to copy blindly.

### P2.3 — Stream and cache prime schedules within a memory cap

`SieveContext` and optional `ScheduleCache` provide bounded packed base primes,
private marking storage, half-open streams, and integer schedules. The bounded
stages retain one segment/batch and generate stage two on demand. Workspace
caps exclude consumer-retained values and interpreter/JIT RSS; cold/warm RSS
and disk/cache comparisons remain required experiment evidence.
M12 measured cold/warm RAM and packed-disk consumption, but a warm utility
gain did not establish a full factoring gain. M17 retains regeneration with
the optional cache disabled; broad cache promotion remains open.

- [ ] Build a segmented prime iterator with packed marking storage and exact prime-power schedules. Cache/reuse bounded prime-gap and distance schedules by bound and endpoint convention; avoid materializing the whole `B2` list. Generate stage 2 only when needed. Distinguish reusable integers from curve-specific points.
- **A:** peak schedule storage respects a configured cap; values above `2**32` cannot overflow packed storage; streamed and materialized reference schedules agree. Stage-1-only success avoids unnecessary stage-2 generation.
- **E:** compare bounded RAM cache, regeneration, and packed disk storage across repeated curves; measure total RSS, sieve CPU, startup, and complete prime consumption. Include cache cold/warm cases and interrupted generation.

### P2.4 — Add chunk checkpoints and a consistent recovery protocol

Stage-one chunks, fine saturation replay, stage-two term replay, and versioned
checkpoints are implemented. Native tests verify exponent action, final-chunk
recovery, and identical resumed work/events/results. Chunk tuning, replay cost,
and checkpoint-frequency measurements still require a recorded decision.
M12 swept chunk lengths 1/4/8/16/32/64; M17 retains 16 without claiming an
optimal replay/checkpoint frequency. Fine checkpoint-cost tuning stays open.

- [ ] Apply exact stage-1 prime powers in reusable moderate chunks for ECM and p−1. Save chunk-start state; periodically GCD-check; replay a saturated chunk at finer granularity. Checkpoint versioned modulus, schedule, seed/work position, and remaining budget.
- **A:** chunked execution matches the reference exponent action on nondegenerate cases; recovery yields a proper factor or bounded retry. Resume rejects incompatible/corrupt metadata and preserves result reconstruction.
- **E:** sweep chunk/checkpoint lengths; force saturation in early, middle, and final chunks. Compare a full run with pause/resume under the same work allocation. Measure replay and repeated ladder-start overhead before choosing chunk sizes.

### P2.5 — Freeze a benchmark corpus and runner

The frozen combined corpus has 20 training and at least 20 held-out inputs
per declared band, including forty added Carmichael controls. Oracles use
independent prime certificates and stay outside algorithm inputs. The isolated
runner records both modes, cold/warm costs, stage events, censored outcomes,
CPU, RSS, JIT/runtime/source metadata, and five seeds. Competitors are pinned
but have not been executed. Broad repeated evaluation remains an open gate.
M17 adds fresh confirmation on three bands and repeated current-source
controls on five more. Large-band screening remains feasibility evidence;
it cannot close broad repeated or competitor evaluation gates.

- [ ] Add proposed `benchmarks/` fixtures and runner with known factorizations hidden from algorithms. Split tuning from held-out evaluation. Pin competitors named by the audit and record whether factoring itself uses native code. Separate factor-one from complete factorization and cold from warm runs.
- **A:** raw results include runtime/backend versions, commit, CPU/core limits, time/memory caps, seed, certainty, stage timings, reconstruction, timeout, and peak RSS. Report completion, median, p90/p95, and censored timeouts; never average successful runs alone.
- **E:** start with at least 20 inputs per feasible workload band and five independent seeds for randomized methods; expand inconclusive cells. Cover balanced 20–80-digit bands under caps, unbalanced ~60/~100-digit inputs with 5–30-digit smaller factors, p±1 boundary cases, powers, primes/pseudoprimes, and close/random composites. Mark infeasible bands explicitly.

### P2.6 — Port reusable sieve contexts and rolling strikes selectively

Reusable contexts and an optional rolling-offset arm pass tiny/restarted/high
interval sequence tests. Rolling is disabled by default pending full factoring
and memory measurements. No native C crossover threshold was adopted.
M12's rolling arm regressed in training factoring. M17 retains plain strikes;
no rolling promotion is pending for the frozen baseline.

- [ ] Adapt the C v2 rolling next-strike state and v4 bounded reusable context to Python. Cache exact base primes and scratch storage by maximum bound; separate reusable integer schedules from modulus/residue-dependent values. Compare an odd-bytearray slicing baseline with carried offsets, including arbitrary interval restarts and short final segments.
- **A:** complete values match independent references on tiny segments, prime squares, random/high-offset intervals, and repeated context calls. Context limits are enforced; each concurrent caller owns its scratch state. Translate C inclusive `[lo, hi]` to Python `[lo, hi+1)` explicitly.
- **E:** measure fresh versus reused contexts, setup/mark/extraction separately, and total factor-stage consumption. Retain only wins after accounting for Python bookkeeping, regeneration, and cache memory. The C v4 context speedup is not a prediction for Python.

### P2.7 — Compare wheel/pre-sieve and packed-output ideas from C v3/v4

M9 retains a simple wheel-6 bytearray small sieve after exact sequence checks
and matched Python measurements. This becomes the small-sieve control.
Wheel-30, pre-sieve, integer-bitset, and packed-output experiment arms now exist
with independent sequence tests. They remain unpromoted pending the declared
complete-stage, consumption, and memory comparisons.
M12's alternatives lost prime-consumption comparisons; a small pre-sieve
factoring win was not broad evidence. M17 retains wheel-6/bytearray and defers
additional kernels until a measured bottleneck justifies them.

- [ ] Experiment one at a time with wheel-30 candidate mapping, pre-sieve pattern copies, packed prime gaps/64-bit arrays, and extraction strategies. Keep bytearray slice marking as the control. Count-only is a diagnostic workload; factorization needs actual primes. Defer width-gated wheel-210 and sparse buckets until profiling justifies their additional state.
- **A:** wheel-residue mapping, sign/phase resets, tail masks, and primes dividing the wheel validate independently; output is ordered, complete, and half-open. Packing above `2**32` cannot truncate values. Consumption/decoding is included in memory accounting.
- **E:** benchmark relevant B1/B2 ranges and repeated interval consumption, not only native C count rates. Measure Python-int bitsets versus bytearrays and complete factoring stages. Do not copy M4 C crossover/thread thresholds or rejected C controls without a Python A/B result.

### P2.8 — Establish whether parallel candidate search is worthwhile

`parallel_candidates.py` now implements fixed assignments and first-factor
cancellation for ECM/rho, threads, and spawned 2/4-process pools, including
cold startup, warm reuse, CPU, and aggregate RSS accounting. The M14 probe
completed 48 configurations and 432 samples, with no measured CPU/RSS gate
failures and complete worker RSS reporting. Larger fixed-work ECM throughput
improved with reused processes, but first-factor stopping and cold startup
did not justify a production switch. **Early feasibility gate passed;
decision: retain serial.** Broad held-out promotion remains a future gate.
See the final PyPy probe.

- [x] Run an early feasibility probe with identical independent ECM curve assignments in serial, threads, and 2/4 spawned processes. Record the actual GIL/build/backend mode. Test small and amortized job sizes, cold worker startup and warm reuse; later include seeded rho walks. Keep production serial unless a declared workload wins.
- **A:** every assignment is unique and reproducible; returned factors validate; exhausted jobs remain explicit failures. Compare the same candidate set without giving algorithms known factors. Aggregate RSS and total CPU budgets include all workers.
- **E:** measure wall time, CPU-seconds, startup, serialization/IPC, and first-valid-factor cancellation latency on PyPy Python 3.11. A fixed-work throughput probe does not establish early-stop benefit; both experiments are needed before adoption. PyPy threads are a comparison arm with a GIL, not an assumed bigint speedup. Historical CPython probes remain evidence, with no new CPython gate.

**Phase 2 core exit (M17): passed.** Cancellation, exhaustion, and resume
acceptance passed; owned workspace is capped and measured RSS stayed below
the declared evaluation limit. The corrected serial portfolio, corpus, and
baseline configurations are frozen and reproducible. Phase 3 relation-engine
work may use this baseline. The open P2.1–P2.7 experiment gates still prevent
claiming complete tuning or promoting unmeasured alternatives.

## Phase 3 — Add the missing balanced-composite engine

**Prerequisite:** Phase 2 schedules, contracts, and corpus. **Goal:** a correct Python QS/MPQS milestone followed by SIQS and a fair SSS challenger. Source: MPQS/SIQS design.

**Research gate (M19): complete; implementation gates remain open.** The
[quadratic sieve research report](quadratic_sieve_research.md) reviews eight
pinned implementations, distinguishes source findings from performance
claims, and maps concrete decisions to the tasks below. The
[source manifest](quadratic_sieve_research_sources.json) records retrieved
files and hashes; no competitor was installed or executed.

**Optimization follow-up (M21): complete as research.** The
[literature/blog follow-up](quadratic_sieve_research.md#optimization-follow-up-m21)
adds pinned FLINT source, Hart's implementation blog, polynomial-selection
research, sparse filtering and batch smoothness references. The
[supplement](m21_quadratic_sieve_optimization_sources.json) preserves fetched
hashes separately from M19. Fold the experiments into P3.2–P3.4, P5.4 and
P6.2; no implementation checkbox or performance gate is closed by research.

**P3.1 reference gate (M20): complete.** The exact `qs/` package, exhaustive
collector, independent root/factorization oracles, and provenance/cap fault
injection passed 120 PyPy tests and lint. Reference costs, raw measurements,
and source hashes are in the acceptance summary
and [public changelog](../../CHANGELOG.md). M20 did not supply complete extraction;
M26 accepts P3.3 below; M31 accepts P3.4's bounded integration and declared
large-band evaluation. Broader SIQS scaling and P3.7–P3.8 remain open;
P3.5/P3.6's bounded evaluations retain experimental arms.

P3.1–P3.3's exact pipeline is accepted. P3.4 implements SIQS families,
root reuse and optional bounded dispatch; M31 completes its declared
large-number evaluation with an explicit retain-ECM-default decision. P3.5's verified bounded challenger and declared comparisons
retain SSS as experimental. P3.6's bounded worker evaluation retains serial;
parallel promotion and P3.8 optimization experiments
remain open; P3.7 NumPy stays conditionally deferred. M26's weight-two
filtering and collector controls retain their measured scope below.
P5.4 double-large-prime work may move forward when useful partial yield
justifies it. SSS, double/triple large primes, sparse solvers, NumPy and
parallelism do not block the stated GNFS prerequisites. The M17 50 ms probes
are feasibility evidence; declare suitable Phase 3 resource budgets before
coverage experiments instead of inheriting those caps silently.

### P3.1 — Implement polynomial and relation identities

- [x] Add the `qs/` package and immutable relation/provenance structures.
  Build a factor base with explicit handling of 2, signs, and primes dividing
  `N′=h*n`; GCD-check `h`. Choose `A,B` with `B² ≡ N′ (mod A)` and exact
  `C=(B²−N′)//A`. Use normalized `F(x)=A*x*x+2*B*x+C`; record all exponents
  of `A*F(x)`, including A, and explicit square corrections for combined
  relations. Compute A targets with integer roots/comparisons.
- **A:** every atomic relation verifies `(A*x+B)²−N′ = A*F(x)` and its exact
  factorization; combined relations preserve their checked atomic provenance.
  Handle zero values before division, negative positions, half-open blocks,
  primes dividing A, and inversion failures explicitly. A normalized-F sieve
  cannot blindly copy full-square-difference skip-A behavior.
- **E:** compare a tiny exhaustive QS/MPQS collector with independent arithmetic
  and modular-root oracles, including negative values, repeated factors,
  multiplier factors, `p | A`, 2, empty blocks and tails. Gate all faster
  collectors and SSS adapters on the same exact verifier.

**P3.2 collector gate (M22): complete.** Reusable conservative score buffers,
exact root/bucket exponent recovery, bounded FIFO partial matching and pinned
combined provenance passed 134 PyPy tests and lint. Safe scoring missed no
admissible values in three training fixtures or two frozen held-out windows.
Matched collection costs regress 32–76% on these small fixtures; no production
or performance promotion is made. The acceptance summary
records raw measurements, source hashes, cold resources and separate profiles.
Filtering/extraction and serialized SIQS integration remain P3.3/P3.4 work.

**Carry-over (M24):** root-hit filtering and bucket recovery were evaluated;
a distinct candidate-only resieving pass was not implemented. That experiment
and controlled diagnosis of weak candidate rejection now belong to
[P3.3](#p33--implement-filtering-gf2-dependencies-and-factor-extraction).
P3.2 remains accepted as a correct bounded collector, without speed promotion.

### P3.2 — Build a bounded single-large-prime collector

- [x] Sieve reusable blocks with cached roots/logs; form big integers only for
  candidates. Compare list loops, bounded bytearray translation/slice updates,
  and array scores on PyPy. Handle prime powers, repeated hits and saturation.
  Store bounded one-residual-prime partials; combine matches with square
  corrections and provenance. Initially constrain residuals to the existing
  deterministic domain `r < 2**64` and a smaller configured collector bound.
  Separate the full window from working blocks and prime-metadata chunks;
  compare dense marking, sparse direct hits and bounded buckets. Compare
  full marking with omitting a small-prime prefix and recovering its exact
  contribution on candidates before expensive residual work.
- **A:** full/combined relations pass the verifier. Charge residual primality,
  GCD/inversion, combination and retained storage. Cap behavior and deterministic
  eviction cannot remove atoms referenced by accepted combinations. Overflow,
  threshold rounding and skipped-prime allowances have documented bounds;
  intentionally lossy scoring must be distinguished from exact hit coverage.
  Approximate or saturated scores alone cannot certify that division is
  complete. Any score-guided early exit needs a proved remaining-hit invariant
  and independent exponent-recovery checks, or retains full-division fallback.
- **E:** compare scores/candidates against exhaustive small-window enumeration.
  Sweep block widths and thresholds; count missed smooth values, false
  candidates, division cost, useful verified yield, slice allocation and RSS.
  When division dominates, compare root-hit filtering/resieving and bucket
  lookup with full factor-base trial division, including hit-storage costs.
  Unrolling must retain every hit and tail; no native cutoff is a PyPy default.
  Sweep small-prime cutoffs and metadata/block sizes on training inputs;
  include high prime powers, missed values and false candidates at each
  candidate stage. Freeze cutoffs before held-out completion measurements.

### P3.3 — Implement filtering, GF(2) dependencies, and factor extraction

**Core and M24 carry-over gates (M26): complete.** Exact provenance-aware
filtering, lifted bitset dependencies, modular square extraction and bounded
in-memory resume pass 151 PyPy tests and lint. An independent dense oracle
checks 264 filter/pivot kernel comparisons. The
acceptance summary records source hashes,
raw samples, limits and decisions. Complete QS on 16 fresh small balanced
inputs is 21.2% faster than M22 with identical frozen bucket/filter settings
(27.220 to 21.446 ms per cohort); a separate same-root cohort improves 21.3%.
These measurements include setup through verified, reconstructed factors.
Cold startup has no demonstrated gain. Candidate-only resieving and tighter
conservative scoring remain optional after showing no completion-time gain.
Tiny collector controls still trail exhaustive enumeration. P3.4 owns SIQS
families, serialized checkpoints and production crossover/dispatch promotion.

- [x] Remove exact duplicates and iterative singletons with provenance; use
  Python integer bitsets and dependency masks for elimination. Preserve direct
  dependencies between distinct relations of equal parity. Compare optional
  weight-two constraint elimination with singleton-only filtering before
  higher-way merges or a new solver. Retain full exponents
  and square corrections. Reconstruct X/Y modulo n without an unbounded
  integer product of all relation values. Try multiple dependencies and both
  `gcd(X−Y,n)` and `gcd(X+Y,n)`; resume after trivial dependencies.
- **A:** independently verify original-row parity, even exponent totals, and
  `X² ≡ Y² (mod n)`. Output only proper divisors. Charge matrix fill-in and
  provenance workspace as well as sparse input. Relation counts and completed
  elimination never signal factorization success by themselves.
  State matrix orientation; a prime column in two relation rows permits their
  XOR and constraint removal. Lift dependencies through every transformation;
  retain zero-row dependencies instead of dropping them as useless data.
- **E:** use independently known dependencies, corrupted provenance, duplicate
  rows, singleton cascades, signed values and trivial congruences. Complete
  small balanced fixtures before scaling. Compare pivot/filter strategies
  using post-filter dimensions, nonzeros, fill-in, memory and total completion.
  Sparse solvers wait for measured matrix cost or bitset infeasibility.
  Include equal-parity/different-exponent rows, weight-two cascades and empty
  reduced rows; compare lifted kernels with a tiny independent dense oracle.
  Measure collection/solve stopping rules using filtered usable-row excess,
  not just raw partial counts; include repeated filtering and trivial GCDs.

**P3.2 carry-over (M24): collector diagnosis and candidate-only resieving.**
M26 completes this follow-up after the exact extraction baseline above.
Its controlled comparisons retain M22's historical diagnosis:
M22's stable matched collection medians regress 32–76%, or 0.243–0.304 ms
per 257-position window. Safe scoring selects 256–257 positions, so weak
rejection warrants targeted diagnosis. These tiny fixtures do not establish
a production SIQS regression. Existing
profiles identify budget/type
validation and exact recovery as hypotheses, not proved causes; instrumentation
changes JIT behavior and cumulative times overlap.

- [x] Prototype a distinct bounded candidate-only resieving pass and compare
  it with full division, cached-root filtering and bucket recovery. Diagnose
  candidate selectivity, marking/hit storage, exact division, verification,
  budget checks and setup/reuse separately. Evaluate tighter conservative
  score bounds and small-prime omission without promoting lossy thresholds.
- **A:** recover every exponent, including omitted primes, repeated powers
  and A's factors; preserve signs, residual certainty and combined provenance.
  Prove score/overflow and remaining-hit bounds before any division early exit,
  otherwise retain full recovery. Preserve finite work/time/storage allowances
  and first-uncommitted-position refusal; test independent coverage, tails,
  singular roots, saturation, cancellation and resume.
- **E:** use unprofiled matched comparisons for causal tests, changing one
  component at a time. Keep identical inputs, assignments, output verification
  and budgets across exhaustive/full/root/bucket/resieving controls. Include
  setup, hit storage, conversion, output consumption, filtered usable-row
  excess, matrix/extraction costs, CPU and peak RSS. Measure cold and warmed
  runs separately on the supported PyPy, with validated warmup and stable
  repeated samples. Retain M22 windows as controls; tune separately and freeze
  before fresh held-out small balanced QS completion tests. Record an
  adopt/defer/reject decision under the end-to-end promotion policy; a utility
  win or higher raw relation count cannot close this follow-up. Broader SIQS
  crossover and dispatcher promotion remain P3.4 work.

**P3.1–P3.3 follow-up (M28): accepted.** Storage caps now attempt bounded
extraction from retained rows, with idempotent exhaustion and budget resume.
Live collector/preparation/matrix reservations are combined before allocation.
The audit passes 157 PyPy tests, lint, 214 root comparisons, 648 independent
collector windows and 320 lifted-kernel comparisons. The fixes cost 3.4–6.5%
on the fresh matched small QS cohort; no performance promotion is claimed.
Existing optional resieving/scoring and later matrix experiments remain
separate. At M28, P3.4 was the next checkpoint/dispatch/large-band milestone;
its bounded M31 acceptance is recorded below.

### P3.4 — Add SIQS self-initialization and bounded dispatch

**M31 bounded implementation/evaluation complete (2026-10-04):** the longer
trained/frozen comparison, varied-input study, cold/control measurements and
cumulative 50-digit continuations are complete and validated. The repaired
source is confirmed on the same inspected input/seed pairs without retuning.
An empty-store balanced 50-digit SIQS run succeeds in 1,183.355 seconds;
60–80-digit capped exploration remains unfinished. Accept this bounded
baseline, keep ECM automatic and SIQS opt-in, and retain general crossover,
arithmetic/matrix scaling and fresh promotion evidence as separate gates.
The checkout passes 230 PyPy tests and lint. Useful raw evidence and exact
measured source remain locally archived with hash verification; v1 is unchanged.
See [full results, commands and limits](../benchmarks/README.md#completed-larger-evaluation-and-filtering-repair-m31-4-october-2026).

The large-state audit repairs exact prime-power coverage, sparse partial-store
accounting, finite capacity, incidence filtering/cancellation and wide-mask
checkpoint fingerprints. All 2,010 old/new filter comparisons match exactly;
root/collector/dense-kernel oracles and pause/resume/corruption controls pass.
The measured filter control is 42.1% faster on its inspected representative;
this is not a global engine promotion. Measured v4 and repaired v5 source
identities stay distinct, including the final compatible digest-streaming
change. Original M30 evidence below retains its historical scope.

**Historical M30 review (2026-10-04, superseded by M31):** the earlier claim
that P3.4's experiment gates were complete
is withdrawn. The 0.2-second large-input runs are bounded diagnostic probes,
not adequate evidence of practical 50–60-digit factoring or the capped
70–80-digit exploration required below. Small 24–26-bit inputs establish
correctness and local controls. Longer wall/CPU budgets and workload-appropriate
base, A, interval, family and relation limits must be trained and frozen before
meaningful held-out comparisons. Record timer/work/storage/schedule exhaustion
separately; extending the timer alone leaves other finite limits in place.
Default promotion remains subject to the end-to-end policy.

Implemented: shared verified full/partial stores, exact CRT/Gray roots,
complete compact checkpoints with charged root/matrix/prefix reconstruction,
finite width/yield/trivial-dependency recovery, integer multiplier scoring,
and optional ECM-to-SIQS dispatch are implemented. 185 PyPy tests and lint
pass; the 214 root/648 collector/320 kernel controls remain valid.

The independently certified corpus freezes 18 training and 52 held-out inputs,
including 32 small 24–26-bit inputs and 20 balanced 30–80-digit inputs, at two
seeds. Matched small completion improves from 41/64 fresh-store input/seed
cases to 64/64 shared; the declared larger classes have no completion
regression against the fresh-store control. Frozen tuning gives 19.7 ms per
32-input completed SIQS cohort, versus 44.0 ms with the shared three-factor
control. QS is 18.6 ms and bounded ECM 4.0 ms on that small cohort. These
cohort figures are not per-input or v1 comparisons.

QS/MPQS/SIQS complete 0/40 large input/seed cases at the declared 0.2-second
wall/CPU caps; ECM completes 2/8 cases at 30 digits and none above. Retain ECM
as the default and add no digit cutoff: whole-portfolio promotion/scaling is
not declared complete. Scored multipliers and width recovery remain optional;
the extended recovery comparison uses five-second warmup and 31 stable
samples, with no meaningful benefit. All factoring outputs reconstruct,
including unresolved cofactors. See [commands, costs and limits](../benchmarks/README.md).

Snapshots retain checked-store, seed/family/Gray/block/resource identity,
compact pending elimination/extraction progress and integrity markers.
The matrix and caches are replayed under the resumed allowance. Width-only
growth retains the same base and needs no remapping; base growth/spill are
disabled and disk use is zero. Memory/checkpoint refusals preserve checked
in-memory state. The conditional remapping/spill clauses below therefore do
not require an additional backend for this accepted bounded implementation.
Earlier M29 evidence remains the family-only foundation, rather than the
current completion state.

**Foundation progress (M29):** exact squarefree CRT families, bounded seeded
A assignments, Gray/recentering updates and incrementally cached roots are
implemented. The collector independently certifies supplied complete roots.
Compact family-only checkpoints retain consumed resources and rebuild caches
under the same allowance. 167 PyPy tests and lint pass. A fresh small probe
completes 13/16 inputs equally in full/cached root arms, preserving three
cofactors at finite family exhaustion. Relation stores are still per
polynomial; shared store/full-job checkpoint/dispatch, multiplier tuning and
large-band experiment gates below remain open. No whole-SIQS promotion.

- [x] Generate reproducibly assigned families with cached CRT/Gray-code B
  updates and incremental roots, including recentering and `p | A` cases.
  Tune factor-base size, bounded multiplier selection, A target, block width
  and thresholds using training costs; retain `h=1` as a control. Integrate
  SIQS after bounded ECM under one setup-to-extraction allowance.
  Add a finite Knuth–Schroeppel-style score with multiplier-size penalties
  and exact residue/modulo-8 tests. Tune A factor count/range and reproducible
  A-set diversity; track symmetry/translation duplicates and distinct yield.
- **A:** roots agree with full recomputation and identities remain exact.
  Checkpoints retain family/Gray index, block position, seed state, verified
  relation-store identity and consumed resources. Bound full/partial relations,
  caches, matrix/provenance and checkpoint bytes. If spill is needed, use a
  versioned compact store with integrity markers and disk caps, not duplicate
  JSON matrix dumps. Exhaustion preserves cofactors and resumable state.
  Define a finite response to exhausted families, stalled useful yield or
  repeated trivial dependencies. Optional base/width growth charges setup,
  remaps verified relations and base identities, and never resets the budget.
- **E:** compare QS/MPQS, SIQS and bounded ECM on declared 30/40/50/60-digit
  bands, with capped 70/80-digit exploration. Freeze held-out inputs, including
  multiplier-sensitive/residue classes. Include setup, marking, candidate and
  residual work, filtering, matrix, extraction, cold/warm costs and memory.
  Validate pause/resume and corruption handling. Promote only under the
  end-to-end policy; measure the crossover instead of importing a digit table.
  Compare simple A schedules with diverse families, and fixed settings with
  bounded recovery. Include duplicate rates, useful-row excess, repeated
  solving, root-cache reconstruction and workspace needed during remapping.

### P3.5 — Evaluate Smooth Subsum Search on the same interface

**2026-10-04 bounded challenger evaluation:** the independent adapter and
unchanged upstream reproduction are implemented and verified under PyPy
Python 3.11. Repeated held-out small/30-digit comparisons and capped 40–60-digit
diagnostics are recorded in the [benchmark guide](../benchmarks/README.md).
Keep SSS/SSSf experimental: neither broad benefit nor a dispatcher policy is
established. Explicit `--method sss` / `sssf` and `PortfolioConfig(sss=...)`
provide opt-in use with shared allowances and checked full checkpoints;
automatic defaults stay unchanged. Immediate cost improvements and a fresh
comparison against a feasible, trained P3.4 control now belong to P3.6.1;
five-second larger probes do not close practical
50–60-digit performance gates.

**M19 historical research:** author-associated code was pinned and reviewed
without execution. Its paper's 75–100-digit results measure one-hour relation
yield, not completed factorizations. The new reproduction installs optional
SymPy/gmpy2 on PyPy and labels upstream settings/backend separately.

- [x] Reproduce compatible upstream settings as a labelled comparison arm.
  After P3.3, prototype an independently implemented collector adapter with
  our exact verifier, seeded assignments and common postprocessing. Include
  bounded product/remainder-tree work and full exponent recovery. Record all
  adaptations; do not call a changed collector an unchanged reproduction.
- **A:** challenger relations and final outputs pass independent verification;
  setup, residual handling, conversion and memory are included. Retain the
  supported PyPy runtime and disclose SymPy/gmpy2/backend availability. The
  upstream harness's NumPy statistics import is not a Factor dependency.
- **E:** after working SIQS, compare SSS/SSSf with single-large-prime SIQS
  under matched wall/CPU/memory/core budgets and held-out inputs. Report useful
  dependencies and complete factorization separately from collector-only
  diagnostics. Adopt only through the promotion policy; remain a challenger
  if upstream compatibility or end-to-end benefit is unestablished.

### P3.6 — Evaluate coarse parallel relation/candidate collection

**2026-10-04 bounded evaluation: accepted; retain serial.**
`ParallelSIQSJob` assigns the same finite SIQS families independently of worker
count and merges verified batches in assignment order. Parent-owned work leases,
aggregate cooperative CPU, bounded owned/result storage and quiescent checkpoints
retain consumed resources; incomplete private families replay with the same ID.
Eleven new tests cover serial/thread/spawned-process equivalence, changed-worker
restart, cap and budget refusals, pending admission, cleanup failures, and direct
residual splits without shortening fixed-work collection. The shared checkout
passes 226 PyPy tests and lint. Forty-two trained/held-out configurations provide
3,144 stored timed attempts with reconstruction and resource checks; source hashes
remain unchanged. See the [benchmark guide](../benchmarks/README.md#p36-coarse-siqs-workers-4-october-2026).
No worker arm passes the end-to-end gate. Small complete-factor cohorts are slower
with workers; the 13-digit challenger also has declared batch-cap refusals and
lower completion than native serial SIQS. This closes the bounded evaluator,
without promoting a dispatcher policy or practical large-band parallel factoring.
The API remains experimental. Immediate accounting, batch/capacity and stopping
improvements now belong to P3.6.1; broader ECM/portfolio parallel integration
remains P6.3. No held-out result retunes this configuration. RSS is a conservative
sum of process high-water marks rather than an OS-enforced or sampled simultaneous
ceiling. Worker startup before its first CPU publication and atomic operations
can overshoot cooperative time allowances; pool lifecycle costs are measured
separately in the cold arm.

- [x] After relation correctness exists, assign independent SIQS polynomial families or SSS search partitions to workers. Compare serial, thread, and process execution with bounded result batches; centralize or safely partition filtering/provenance. Reuse the parallel feasibility measurements from P2.8.
- **A:** no duplicated/lost candidate partitions; every returned relation passes the verifier; stop/resume preserves provenance and worker memory caps. Do not serialize every raw sieve position when only verified relations are needed.
- **E:** include worker startup, relation transfer, residual factoring, filtering, linear algebra, and early cancellation in complete time-to-factor. Keep serial collection unless the held-out promotion gate passes; candidate throughput alone is insufficient.

### P3.6.1 — Immediate diagnosis and improvement of P3.5/P3.6

**2026-10-04: diagnosis recorded; implementation and improvement gates open.**
This milestone owns the immediate SSS/SSSf and coarse SIQS worker cost work.
It starts from the accepted bounded implementations and does not depend on
finishing P3.8 or P6.3. P3.8 retains broader matrix/array reconciliation;
P6.3 retains general ECM/portfolio parallel execution and reuses this work.
P3.5/P3.6's historical acceptance and experimental/default decisions stand.

**Why now:** the [cost diagnosis](../benchmarks/README.md#p35p36-cost-diagnosis-4-october-2026)
finds concrete overhead in working code. Selected 30-digit SSS/SSSf stage
timers put about 79%/81% in collision generation and only 3%/5% in smoothness
trees. SSSf's accepted small cohort rejects 4,403 of 8,084 generated candidates
and needs 418 assignments versus SSS's 183. A small two-process profile makes
2,494 parent CPU polls; a medium worker family checks the shared cancellation
event 12,619 times. Complete-family publication delays extraction: an accepted
small first-factor cohort scans 2,048 positions in native serial versus
15,842 with four processes. Capacity also affects completion: P3.6's medium
first-factor cohort completes 6/8 versus native serial's 8/8 because of batch
refusals; fixed schedules complete 5/8. Historical larger SSS/SSSf probes show
memory refusals and a restrictive fixed residual bound relative to the base.
These are implementation, capacity and useful-yield issues to investigate
before expanding expensive parallel runs. Profiles identify targets, not
speedups or a general verdict on parallel SIQS.

- [x] Record the bounded initial diagnosis, separately from performance
  evidence: 16 configurations, validated warmups, 189 stage-timed attempts
  including noise extensions, and separate parent/worker profiles. Current
  loaded source hashes stay unchanged; every result reconstructs. Inputs
  are already inspected and cannot serve as fresh promotion evidence.
- [ ] **Freeze controls and isolate causes.** Retain an immutable current
  source/configuration control, then distinguish collision generation,
  budget/clock checks, worker collection, transfer, central verification,
  filtering/elimination, extraction and cancellation. Keep native serial
  and coarse serial controls, and compare 1/2/4 workers. Explain pure IPC
  and thread synchronization costs separately where measurable; their
  shares and any GIL effect were not isolated by the initial profiles.
- [ ] **Reduce accounting and cancellation overhead first.** Amortize
  shared CPU-array reads, event checks and clock reads over explicitly
  bounded chunks. Preserve exact work charges, a shared deadline, live
  aggregate CPU, final accounting and cancellation on a validated split.
  Declare maximum polling latency and cooperative overshoot, including
  worker startup and the largest atomic action. Exercise limits and resume
  at chunk boundaries; removing global checks is not an eligible speed arm.
- [ ] **Improve SSS collision generation.** Reduce repeated affine-root,
  signed-shift counter and invariant setup work within each assignment.
  Compare candidate sets against the current exact implementation, including
  distinct-prime collision counts, singular roots, dropped primes and the
  forced-divisor quotient. Preserve seeded assignment identity, full exponent
  recovery, finite scratch storage and refusal/replay behavior. Tree and
  matrix rewrites follow only if new profiles make them the dominant cost.
- [ ] **Address SSSf useful yield.** Train selection/filter settings as new,
  separately labeled arms; measure rejected candidates, useful verified
  rows/dependencies per assignment and complete factorization. Retain SSS
  as the control. The present filter loses yield, so fewer candidates at
  exponent recovery is not an improvement gate. Preserve the accepted
  configuration and never tune on held-out outcomes or known factors.
- [ ] **Investigate capacity refusals and insufficient useful yield.**
  Separate wall/CPU/work exhaustion, assignment exhaustion, candidate/tree
  caps, owned-memory refusals, partial/relation limits and missing useful
  dependencies. For SSS/SSSf, inspect the residual-bound/base-prime mismatch
  and the larger-probe tree/storage refusals; train bounded residual and
  capacity policies as separate arms. For P3.6, reproduce the medium batch
  refusals and compare finite chunking with separately trained larger caps.
  Measure completed factors and usable dependencies, not just candidate
  throughput; an unresolved result must retain its cofactor and spent budget.
  Investigate actual live/queued storage versus conservative reservations
  without presenting process-lifetime RSS sums as simultaneous peaks or
  enforced limits. Larger-band failures motivate investigation rather than
  establish practical completion at those sizes.
- [ ] **Publish smaller verified worker batches and extract earlier.**
  Compare block/polynomial chunks with complete-family tasks and train
  bounded batch capacities. Stable family/Gray/block IDs and committed versus
  pending cursors must prevent lost or repeated admission across worker
  counts and checkpoints. Account for all extra completed, cancelled and
  discarded work. Keep deterministic merge and any race-dependent stopping
  explicit; fixed-work mode must still finish the same mathematical schedule.
- [ ] **Reduce central and transport overhead where the diagnosis supports
  it.** Measure duplicate verification/setup, matrix growth and serialization
  independently. Reuse validation/preparation only under unchanged verified
  store/base identity; worker output and restored checkpoints still require
  exact verification. Bound queued and coexisting copies, preserve partial
  matching/provenance, and test reused pools separately from cold lifecycle.
- [ ] **Validate improvements and record adopt/retain/reject decisions.**
  Use independent training and fresh held-out inputs with repeated seeds,
  matched work/wall/aggregate-CPU/memory/core budgets, and a feasible trained
  SIQS control for SSS comparisons. Cover the diagnosed small cases and
  declared feasible larger bands within the sub-100-digit target; do not
  infer broad scaling from 8/13-digit worker tasks or five-second probes.
  Include setup through terminal classification, useful yield, wasted work,
  cancellation latency and memory. Separate cold startup, warmed execution
  and instrumented profiles; use at least three seconds of validated PyPy
  warmup and nine samples, extending unstable runs.
- **A:** zero arithmetic/output failures; proper splits and unresolved
  cofactors reconstruct; probable/proven certainty remains distinct.
  Exact relation/provenance, cross-worker resume, malformed checkpoints,
  cap/budget refusals and cancellation retain their contracts. Run
  `make -C v2 test` and `make -C v2 lint` for implementation changes.
- **E:** evaluate one change at a time. Apply the existing promotion policy:
  at least 10% lower complete-run median time or 10 percentage points higher
  completion, with no more than 5 points of regression in another declared
  class, under the same resource limits and with uncertainty reported.
  Keep serial defaults and SSS/SSSf experimental until their respective gates
  pass. Close this milestone only after the planned experiments, acceptance
  checks and explicit decisions are recorded; diagnosis alone does not
  complete the improvement work. Document API/accepted behavior and concise
  measurements, while keeping raw captures and profiles local.

### P3.7 — Optional NumPy spike for bounded SIQS sieve arrays

**Decision (M16): defer; no current dependency or Phase 2 exit gate.** The
current rho, p−1, ECM, roots, and primality paths require arbitrary-size exact
integers. NumPy's numeric integer types have fixed widths and can overflow;
`dtype=object` retains Python objects and is not an unboxed bigint backend.
Upstream NumPy runs through PyPy's C-API compatibility layer, so compatibility
does not establish a speedup. Current prime sieves already use bulk bytearray
slices; no NumPy improvement has been measured.

The plausible experiment is bulk score/root-hit storage for the future
P3.2 SIQS collector, once that work is correct and profiling identifies array
work as a bottleneck. Keep polynomial values, relation verification, modular
arithmetic, and factor extraction in exact Python integers. Sources inspected
2026-10-03: [PyPy's NumPy guidance](https://doc.pypy.org/faq.html#should-i-install-numpy-or-numpypy),
[NumPy integer overflow](https://numpy.org/doc/stable/user/basics.types.html#overflow-errors),
and [NumPy object dtypes](https://numpy.org/doc/stable/reference/arrays.dtypes.html).

- [ ] Only after P3.2 and profiling justify it, prototype bounded vectorized
  sieve-score updates and candidate extraction with an optional NumPy import.
  First verify availability on the supported PyPy Python 3.11 environment;
  record the exact NumPy build and any native dependencies.
- **A:** prove score/index bounds for every intermediate; test overflow,
  duplicate hits, tails, and dtype conversions against an independent integer
  oracle. Relations must pass the existing exact verifier. Preserve budgets,
  bounded memory, and the dependency-free path when NumPy is unavailable.
  Do not use floating-point factoring or ordinary numerical linear algebra
  for the P3.3 GF(2) dependency problem.
- **E:** compare PyPy plus NumPy with PyPy bytearray/array/list loops on the
  same blocks, seeds, and held-out corpus. Include import/cold startup,
  several seconds of validated JIT warmup, allocation, conversion, extraction,
  peak RSS, and complete time-to-factor. Tune separately; apply the existing
  end-to-end promotion policy. Reject or defer if only a microbenchmark wins.

### P3.8 — Evaluate GF(2) matrix and filtering optimizations

**2026-10-04 reconciliation and Phase 3 carry-over.** The
[new research report](phases_three_plus_research.md#phase-3-carry-over-owned-by-p38)
extends this section to own the broader Phase 3 reconciliation, including
collector and workload issues that determine whether matrix optimization is
useful. Immediate diagnosed P3.5/P3.6 cost improvements now belong to P3.6.1.
Keep P3.4's accepted implementation and P3.5–P3.7's existing ownership;
the follow-ups below add evaluation gates without reopening accepted work.
Power-aware sieving and batch smooth-part recovery are present in the reviewed
worktree. Concurrent P3.5 work records the completed bounded SSS adapter and
upstream reproduction; R5 retains that evidence and reuses P3.6.1's fresh
comparison. Broader performance and dispatcher promotion remain separate.
The refresh also finds serialized SSS checkpoints and explicit portfolio/CLI
selection in current code. These need current acceptance evidence and API
documentation; selectable execution is distinct from automatic promotion.

**Current performance target (user clarification, 4 October 2026): general
factoring below 100 decimal digits on the user's Apple M4, 24 GiB, 10-core
machine, with MPQS/SIQS and ECM.** Do not use a 50-digit cohort as the benchmark
or tune the portfolio to that size. Historical milestone records and the
already collected 50-digit exploratory traces remain provenance only. The
[sub-100-digit reconciliation](phases_three_plus_research.md#below-100-digit-optimization-target)
sets the new workload and implementation priorities. All fresh Phase 3 tasks
stay here; P4/P5 own ECM arithmetic and continuation changes.

**P3.8-R1 — Establish a feasible SIQS workload and capacity control**

**Implementation integrated (4 October 2026).** External-square MPQS,
streamed nearest/flyer assignments, larger A support, independent Gray/search
quotas, monotone checkpoint extension and sparse-aware work reservations now
pass the combined 258-test suite and lint, including the P2/P3.6.1 repairs.
Frozen confirmation covers one fresh balanced input per 30/40/60/70/80/90/99-digit
band at two seeds and nine validated warmed cohorts. The 30-digit control
completes in all four arms; upper one-second controls remain censored, with no
useful rows at 60–99 digits. The corpus also defines separate uneven/smooth/
close/power classes. Reference defaults remain unchanged. Checked items below
cover implemented capacity/reachability contracts; broader parameter selection,
Gray tuning and crossover gates remain open. Exact measured sources and
post-capture acceptance changes are recorded in the benchmark README.

- [ ] Freeze trained, feasible configurations before fresh held-out runs.
  Jointly choose factor-base bound and actual cardinality, reachable A targets,
  interval/family quotas, residual bounds, relation/partial stores and matrix
  workspace. Verify that finite family exhaustion does not masquerade as a
  wall-time limit. Keep already inspected M30/large-corpus data as tuning or
  reference data; create untouched confirmation inputs after selection.
- [ ] Extend corpus/runner coverage beyond its present 30–80-digit choices.
  Prespecify total-size bands such as 30/40/60/70/80/90/99 digits, with balanced
  semiprimes and separate uneven-factor classes. Stratify ECM by the smaller
  factor's size, not only n's size; include independently labeled p±1-smooth,
  close-factor and power controls. Keep inspected inputs in training. Report
  censored completion at upper bands, factor-one and complete factoring,
  verified useful yield and memory; no failed run establishes a time-to-factor
  ratio. Train a machine/backend-dependent crossover rather than importing a
  native engine's digit threshold. Upper bands may justify the existing GNFS
  roadmap, without making GNFS a prerequisite for MPQS/ECM optimization.
- [x] Repair polynomial reachability before ranking MPQS against SIQS.
  Current MPQS requires q in the factor base and A=q², hence A<10**10 at the
  maximum bound. Evaluate classical external-q square coefficients near
  sqrt(2N')/M, carrying their known square correction separately from the
  factor-base exponents; independently verify lifted roots and the exact
  relation identity. Keep bounded coefficient/root construction and explicit
  certainty semantics. Current SIQS permits at most eight A factors: even at
  bound 100,000 and maximum half-width 499,999, A<10**40 cannot approach the
  target for any 93–99-digit n. Separate an extendable, memory-bounded family
  generator from hard reference ceilings. This is a target-reachability
  barrier, not a proof that an off-target polynomial can never factor n.
- [ ] Compare A selection by actual product quality, including PARI's final
  compensating "flyer" prime, against the current nearest-p**s pool sampling.
  Jointly tune factor count, attainable A/target distribution, diversity,
  interval width and Gray-family reuse. Larger factor counts increase the
  exponential family length; stream/checkpoint families without preallocating
  the entire search history. For a classical MPQS arm, assess dual batch
  inversion only after its external-square representation is correct.
- [x] Address the verified large-job capacity restrictions (4 October 2026
  follow-up): default SIQS permits at most 64 polynomials, with no width
  growth, and terminates after 16 consecutive windows without new usable
  rows. Expose caller-selected finite search/storage allowances with a
  documented extension policy that retains checked relations and consumed
  resources. Separate implementation ceilings from job allowances: current
  maxima include base bound 100,000, 64 families, 1 GiB SIQS workspace and
  16 MiB checkpoints. Resume currently rejects a changed configuration.
  Evaluate sparse-aware work reservations: a 2,000-row/500-column one-bit
  matrix reserves 5,002,000 units for its first filtering round, exceeding
  the default 2,000,000 before any round executes. Reconcile matrix and
  collector capacities; 32,768 rows/5,000 columns reserve about 1,198.7 MiB
  for matrix work alone. Preserve proved bounds when changing representations.
  Details and the distinction from performance evidence are in the
  [budget follow-up](phases_three_plus_research.md#budget-follow-up-large-job-capacity).
- [x] Separate resident family-batch size from cumulative search allowance.
  Evaluate deterministic streaming/extension of assignments with a stable
  cursor, bounded duplicate tracking and checkpoint growth, preserving the
  verified store and consumed resources. Increasing A factor count changes
  polynomial quality and exponential Gray-family size; it is not a substitute
  for more search assignments. YAFU's continued collection and CADO's
  post-filter feedback are design references, not PyPy parameter presets.
- **A:** every outcome reconstructs n; distinguish timer, work, storage,
  candidate rejection and assignment exhaustion. Preserve complete extraction,
  shared allowances, checkpoints and certainty. A longer timer does not enlarge
  a finite search schedule. Required immutable loader inputs must be retained
  explicitly before publication, with clean-checkout imports/tests verified.
- **E:** measure complete SIQS and ECM-to-SIQS runs, setup, verified useful-row
  yield, filter excess, proper-factor yield and all stage costs. Use at least
  three seconds of validated PyPy warmup and nine samples; extend unstable
  measurements. Separate cold startup, instrumented profiles, owned workspace
  and process RSS. Do not infer 30–80-digit capability from tiny fixtures.

**P3.8-R2 — Compare the candidate cascade before expanding matrix work**

- [ ] Coarsen hot-loop work reservation and clock/cancellation polling over
  bounded metadata/root/hit chunks, retaining exact cumulative work and
  first-uncommitted-position semantics. Validate immutable metadata at its
  boundary rather than repeating public type checks per hit. Define maximum
  uninterrupted chunk/cancellation latency and reserve before execution;
  separate budget bookkeeping from actual arithmetic in causal comparisons.
  Current `_divide` charges a dense factor-base-sized amount before refined
  rejection, even for sparse/resieve recovery: align charges with performed
  evaluation, refinement and visited hit/division work under a documented
  model, without silently granting unaccounted work.
- [ ] Build bounded per-polynomial prime-power root/hit plans and reuse them
  across working blocks. Avoid restarting Hensel lifting for every prime in
  every block and avoid duplicate base-prime passes in power-score mode.
  Bind caches to polynomial, root, interval/valuation bounds and score policy;
  preserve singular/2-adic/exceptional-root fallbacks and charge replay/cache
  construction. Compare segmented streaming with bounded cached plans before
  attempting family-wide algorithms.
- [ ] Make exact recovery sparse in the already recorded hits. Bucket mode
  currently enumerates every factor-base prime per surviving candidate;
  iterate set bits plus cached sparse A support. Resieving currently loops
  over every candidate per prime to add mostly zero A contributions; seed
  A support once and append valuations only for actual candidate hits.
  Preserve repeated valuations, exponent ordering, complete coverage and the
  independent full-division oracle; include metadata allocation/reuse cost.
- [ ] Replace ineffective small-prime omission with a staged scoring
  experiment. In power-score mode, omitting 2 grants a maximum-bit-length
  allowance that collapses both coarse and refined thresholds to zero.
  Compare cheap exact tiny-prime corrections/refinement and bounded fixed-point
  log scores with explicit rounding/overflow bounds. If using an intentionally
  lossy tolerance arm like native sievers, label and measure missed admissible
  values separately; exact relation admission/extraction remains mandatory.
- [ ] Evaluate the present prime-power score policy against conservative
  root weights; then profile bounded resieving, small-prime omission and
  scalar versus batch smooth-part recovery. Record positions, score survivors,
  exact divisions, residual tests, admitted/matched partials and useful rows.
  A restricted 2-adic polynomial variant is a separate optional experiment.
- [ ] If family/root or large-factor-base-prime scanning costs justify it,
  compare per-polynomial marking with bounded family-wide hit scheduling,
  following Kleinjung's sorted CRT half-sums and CADO's experimental SIQS.
  Its p>I condition uses the complete polynomial interval, not the smaller
  working block; verify eligible prime/power counts before investing in it.
  Verify every polynomial/position hit, Gray label, tail and prime-power hit;
  handle p|A explicitly outside the coprime fast path. Bound half-sum tables
  and queued hits before allocation; include construction/replay and complete
  factorization. Experimental class-group/native results do not imply a
  PyPy factoring speedup. Sources are in the research report.
- **A:** independently enumerate small windows and high valuations, including
  p=2, p|A, p|N', singular lifts, tails and score saturation. A capped root
  lift needs a proven conservative fallback. Smooth-part recovery still needs
  complete exponent recovery; any changed polynomial includes its A, sign and
  powers-of-two corrections. Label intentional candidate loss explicitly.
- **E:** charge root/power setup, rejection, primality/cofactor work, tree
  construction and retained bytes. Compare complete factorization and the
  post-filter matrix under jointly trained thresholds; native cache sizes,
  score thresholds and batch sizes are hypotheses for PyPy, not defaults.

**P3.8-R3 — Reduce preparation/filtering/provenance cost first**

- [ ] Profile whole-store preparation on changed relation counts. Concurrent
  work now implements touched-column incidence updates and degree queues;
  evaluate these against the frozen repeated-rebuild control before claiming
  a full-pipeline gain. Compare bounded batches with disjoint incident
  relations, live-column
  compaction with an inverse map, and immutable merge histories with deferred
  lifting against dense original-relation masks. Evaluate verified immutable
  preparation caches and a trained re-filter/solve cadence separately.
- [ ] Give admitted full and combined rows stable mixed-sequence identities
  before introducing incremental preparation/dependency caches: the current
  `_full + _combined` ordering shifts old combined indices when a full row
  arrives. Atomic IDs alone exclude exponent payloads; bind caches to the
  complete immutable atom/base/store identity. Use bounded caches of tested
  dependencies keyed by stable selected row identities, not shifting masks.
  Maintain pivot/nonzero counters incrementally and iterate selected set bits
  for extraction; avoid repeated global sums/scans. Treat matrix scaling as
  the next potential barrier, not the present collection bottleneck.
- **A:** every dependency lifts and passes the original relation verifier.
  Cache identities include the base, polynomial and atom payloads; untrusted
  checkpoint loads retain complete reverification. Bound merge history, fill,
  cache, reconstruction and simultaneous live storage. Replace conservative
  quadratic reservations only with a demonstrated representation/workspace
  bound, rather than reducing a constant to admit larger inputs.
- **E:** include preparation, incidence rebuilds, fill-in, provenance, lifting
  and extraction in genuine matrix/full-factor comparisons. Independent pivot
  batches are conflict-free merges, not a claim of independent dependencies.
  Native/GPU sparse layouts support experiments; their timings do not choose
  Factor's solver or justify disabling its verification.

**Concurrent implementation note:** the initial sparse-cycle rebuild probe
describes the frozen pre-queue source. Its 512 full rounds/134,611,200 units
must not be presented as a fresh result for the touched-column implementation.
Dense lift-mask reservation and uncompacted labels remain separate capacity
concerns until changed representations and their gates are verified.

- [ ] If extraction/provenance storage dominates, separately assess retained
  merge histories and accumulated modular square-root payloads. CUDA-MPQS's
  V1 replay/tree and V2 packed-exponent/incremental-root designs are distinct
  references. Keep original parity/provenance independently checkable, prove
  exponent packing bounds, and retain square corrections and corrupt-state
  detection; an incremental residue cannot certify its own relation history.

**P3.8-R4 — Measure useful dependencies and large-prime economics**

- [ ] When R1–R3 identify insufficient useful yield, evaluate P5.4's bounded
  double-large-prime extension against optimized single-large-prime SIQS.
  Jointly tune each prime limit, residual-product cap, splitting allowance,
  unmatched occupancy, filtering and independent dependency output.
- **A:** preserve the admitted residual certainty domain or explicitly extend
  its proof contract. Recover exact square corrections for repeated primes
  and self-loops; support cycles in every connected component. Track verified
  nonzero/independent dependencies and both GCD signs. Extra QS character
  constraints require a proved purpose and a reproducer before adoption.
- **E:** report proper-divisor yield per verified dependency as well as raw
  relation/cycle counts. Diagnose repeated trivial congruences under raised
  large-prime bounds; recent GPU reports supply test hypotheses, not a general
  explanation or a substitute for a full-factor experiment.

**P3.8-R5 — Reconcile optional SSS, workers and array challengers**

Immediate diagnosis, SSS/SSSf cost improvements, worker accounting/batching,
and their fresh matched comparisons are owned by P3.6.1. This section consumes
that evidence for broader interface/array reconciliation.

- [ ] Retain P3.5's completed bounded adapter/upstream evidence and P3.6.1's
  eventual validated decisions, with unchanged author code as a separately
  labeled arm. Reconcile the common interface and record the
  forced-prime quotient, recovered exponents, SSSf candidate loss, tree caps
  and in-memory versus serialized resume scope. Reuse P3.6.1's checked stable
  assignments and aggregate-budget work for P3.6; derive seeds from assignment
  IDs, rather than worker IDs. Advance P3.7 only for a measured array bottleneck.
- **A:** recovered relations retain all forced factors; output validation
  cannot rely on stdout claims. Duplicate vector indices accumulate every hit
  (`numpy.add.at` or a verified aggregation), with bounded intermediate scores
  and indices. Parallel completion order and early stopping are explicitly
  distinguished from reproducible assignment identity and restart behavior.
- **E:** SSS's 75–100-digit paper runs measure one-hour collection, not complete
  factoring at those sizes. Compare complete outputs on matched inputs,
  cardinalities, budgets and supported PyPy dependencies. Include upstream
  adaptation, conversions, worker startup/IPC, wasted work and cancellation;
  retain the existing optional/defer decisions until their gates pass.
- [ ] Reconcile current SSS serialized checkpoints and explicit portfolio/CLI
  selection with README/API and current tests. Verify cumulative resources,
  interrupted assignment and solver reconstruction, terminal evidence and
  checkpoint-size refusal. The original-budget guard has been restored;
  preserve it for in-memory resume while separately validating serialized
  restoration. Neither API availability nor passing resume tests establishes
  automatic dispatch superiority or larger-band completion.

**Distinct milestone (M25 research; implementation open).** The
[matrix research report](gf2_matrix_research.md) reviews block Wiedemann,
block Lanczos, dense Four Russians/PLE, packed sparse products and stronger
filtering against papers, implementer blogs and seven pinned implementations.
[Source provenance](m25_gf2_matrix_sources.json) and
documentation verification preserve the
review evidence. Existing P3.1–P3.7, M24 carry-over, P5.4, P6.2 and P7.6 text
and acceptance gates remain unchanged; no earlier milestone is reopened.

**Sequence:** establish P3.3's exact filtering/extraction control first.
Matrix-only prototypes can then proceed; integrated promotion requires
P3.4's working SIQS. Start with measured matrix costs, shared kernels and
dense/hybrid methods; pursue iterative solvers when the solve cost or dense
memory requirement justifies them. Record adopt/defer/reject for each track.
This optional optimization milestone does not delay bounded GNFS work;
P7.6 may reuse its exact kernels and measured solver decisions later.
Keep Factor's algorithms in Python on supported PyPy Python 3.11.

**Control and shared kernels**

- [ ] Freeze an exact matrix interface and a representative QS/SIQS corpus
  with original relations and lifting maps. Declare orientation; in this
  milestone `M` has constraint rows and relation columns, and `M d = 0`
  selects relations. Profile the complete reference pipeline. Compare packed
  matrix-times-block kernels, CSR/column/tiled layouts and a dense heavy-row
  split, including transpose products and block widths.
- **A:** independent GF(2) products agree; repeated coordinates cancel modulo
  two, short tails are masked, and every transformation maps to original
  relations. A lane stores an independent binary vector, not an extension-
  field scalar. Bound input, layout, transpose and scratch workspace.
- **E:** compare entry-wise and bitset controls on genuine and synthetic
  fixtures; include empty, rectangular, skewed-density and duplicate-coordinate
  cases. Check `uᵀ(M v) = (Mᵀ u)ᵀv`. Charge conversion, tuning and both matrix
  orientations; native 64-bit/cache thresholds are not PyPy defaults.

**Dense/hybrid elimination and stronger filtering**

- [ ] Compare Python-bitset Gaussian elimination with bounded Four Russians
  tables, rank-aware echelon/PLE-style free-variable recovery and a sparse-to-
  dense trailing core. Separately compare higher-weight merges, minimum-fill
  or spanning-tree merge choices and independent pivot batches with P3.3's
  singleton/weight-two control. Evaluate exact disconnected components;
  surplus-clique pruning is a separate conditional, lossy candidate policy.
- **A:** cap `2**k` tables, merge degree, fill-in and provenance memory.
  Preserve zero/equal-parity dependencies and lift every emitted vector.
  Tiny exact dense results match independent rank/nullity/kernel oracles.
  Component decomposition is exact; pruning must explicitly retain sufficient
  useful excess and never claim to preserve all original dependencies.
- **E:** measure filtering, pivot search, allocation, conversion, elimination,
  recovery and lifting. Compare post-filter dimensions/nonzeros, observed
  fill and complete factor recovery; a smaller matrix or faster XOR alone
  does not pass the integration gate. Strassen recursion waits for a measured
  large dense-kernel bottleneck.

**Block Lanczos challenger**

- [ ] When sparse solving is justified, prototype finite-field block Lanczos
  with seeded finite retries, selected nonsingular Gram subblocks and original-
  operator kernel correction. After correctness, compare homogeneous recurrence
  and low-rank auxiliary updates, block width and packed-kernel variants.
  Otherwise record a profile/memory-based deferral with a trigger to revisit.
- **A:** handle self-orthogonal vectors and singular blocks without invalid
  inversion. Apply `Mᵀ(M V)` implicitly if required; do not form a dense Gram
  matrix. A Gram-kernel candidate must pass `M d = 0`: for
  `M = [[1], [1]]`, `Mᵀ M = [[0]]` has spurious solutions. Correct/filter
  candidate blocks, reject zero/dependent outputs and bound recovery.
- **E:** compare full recurrence, auxiliary operations, retained vectors and
  final verification with the dense/hybrid control. Use rank-deficient and
  singular-Gram fixtures, several fixed seeds and complete QS extraction.
  Measure restart rate and capacity separately from successful timing ratios.

**Block Wiedemann challenger**

- [ ] When justified, prototype seeded block Wiedemann with an explicit
  square-operator/kernel map, projected Krylov sequence, matrix-polynomial
  generator, reconstruction, lifting and original verification. Scalar
  Berlekamp-Massey applied independently to entries is not the block generator.
  Start with a correct base-case generator; compare Thomé-style divide-and-
  conquer only when generator cost dominates. Otherwise record a measured
  deferral and revisit trigger.
- **A:** check unused sequence terms and original `M d = 0`; square padding,
  permutations or symmetric embedding must return a nonzero relation vector.
  Bound projections, sequence extension, generator workspace, retries and
  checkpoint disk/RAM. Exhausted iterations do not prove rank or primality.
  State whether the API returns a verified batch or a complete basis.
- **E:** charge prep/projection, Krylov, generator, reconstruction, checking,
  lifting and I/O; sweep projection widths and inspect failure rates, repeated
  invariant factors and unlucky projections. Compare block Lanczos and dense/
  hybrid controls under matched resources. No universal solver ranking or
  matrix-dimension cutoff follows from native implementation timings.

**Resume, selection and promotion**

- [ ] Define stage checkpoints tied to matrix/atomic-relation hashes, field,
  orientation, permutations, lifting history, operator, seeds/projection
  blocks, widths, sequence extent and consumed allowances. Retain Krylov
  states needed by reconstruction. Freeze any layout/solver selection on
  training matrices, then record held-out adopt/defer/reject decisions.
- **A:** stop/resume and finite failure preserve the common budget and
  unresolved cofactor. Detect bit flips, stale matrices, corrupted provenance,
  incomplete sequences and incompatible checkpoints. Independently verify
  original parity, even exponents, `X² ≡ Y² (mod n)` and every proper divisor;
  handle multiple/trivial dependencies within bounded recovery. Research alone
  closes no gate: require a validated control, at least one justified bounded
  challenger and explicit evidence/decisions for the remaining tracks.
- **E:** use the same PyPy/input/seed/budgets; separate cold startup/JIT and
  stable warmed repeated runs. Record medians/spread/confidence intervals,
  completion/exhaustion, CPU, peak/aggregate RSS, disk and verification costs.
  Synthetic capacity tests do not establish SIQS speed. Integrated defaults
  require the unchanged end-to-end promotion policy above. Any independent
  sequence/process spike must also charge startup, matrix replication, IPC
  and cancellation against identical serial assignments; GPU/native libraries
  remain design references, with no new backend dependency here.

**Phase 3 exit:** at least one validated relation engine demonstrably improves balanced-composite completion under matched limits. The dispatcher can stop, resume, and preserve cofactors throughout the relation pipeline. SIQS and SSS need not both become defaults. The optional NumPy spike is not required for phase completion.

## Phase 4 — Optimize ECM chains and arithmetic after coverage exists

**Prerequisite:** Phases 1–2; end-to-end comparisons use the Phase 3 portfolio baseline. **Goal:** useful stage speedups that survive full-run costs. Sources: PRAC repair, coordinates, reduction, backends.

### P4.1 — Finish PRAC repair and precompute valid chains

- [ ] Use `audit/prac_reference.py` as a prototype, not a drop-in. Guard 0/1/2 and powers of two; select `k//2 < r < k` with `gcd(k,r)=1`, require terminal `d=e=1`, and validate integer/differential invariants. Separate chain generation from execution; cache exact-rational choices and use ladder fallback.
- **A:** zero nondegenerate mismatches against independent affine/ladder oracles; `(0,0)` is always detected and handled by factor extraction, recovery, or retry. The prototype's 797 exceptional cases are not counted as passing point equalities.
- **E:** extend the recorded 16,016 comparisons to larger fields, composite moduli, and actual prime-power schedules. Include chain construction, Python dispatch, exceptional recovery, and total stage-1 time. Enable PRAC only after the promotion gate.

**Cost-model follow-up (M27):** v1's `ADD_COST=6`, `DUP_COST=5` and
`lucas_cost` mirror GMP-ECM's abstract differential-add/double costs.
The review maps the Fibonacci/golden-ratio
split and existing prototype to this task; production still uses the ladder.

- [ ] Compare guarded PRAC with optional precomputed near-optimal Lucas
  chains, as in GMP-ECM 7.0.6. Model additions as 4M+2S and doubles as
  3M+2S, then calibrate costs on actual PyPy kernels and available backends.
  Keep integer/rational selection and charge generation, cache/storage,
  dispatch and exceptional recovery. Independently verify every chain;
  compare full stage-one and complete-factorization costs before promotion.

**2026-10-04 research refinement:** test projective validity separately from
cross-product equality, which can accept `(0,0)` vacuously. Exercise complete
prime-power chains over composite moduli, retaining nonunit/factor information
and charging finite exceptional recovery. GMP-ECM's chain generator is a
reference; abstract M/S counts need measured PyPy costs.

- [ ] Route verified chain execution into the actual bounded stage-one path:
  `multiply_prac` is currently a ladder wrapper and stage jobs call
  `scalar_multiply` directly. Compile small prime/prime-power chain records
  into bounded composed programs; retain chunk-start replay and factor/nonunit
  handling. Do not search near-optimal chains for an enormous full-lcm scalar
  online. Measure chain/schedule amortization across the sub-100 workload.

- [ ] Add bounded/offline continued-fraction chain search as a separate
  challenger, using Bernstein–Cottaar–Lange's 2025 pruning/meet-in-the-middle
  work and CADO's compact bytecode representation as references. Independently
  verify integer chain records and distinguish guaranteed from heuristic
  search termination/optimality. Charge generation, tables, dispatch and
  exceptional recovery; retain the ladder. Author prototype floating search
  bounds and shorter-chain/node-count results are not exact production
  contracts or measured full-ECM speedups.

### P4.2 — Compare fused and normalized Montgomery kernels

- [ ] Add explicit squares, fused ladder addition/doubling, and optional normalized fixed-difference kernels. Reuse setup inversion only where the algebra supports it. Maintain readable oracle formulas; specify intermediate-width and reduction bounds.
- **A:** kernel variants agree on valid points and correctly surface degeneracy/nonunits. Chunk transitions and normalized representations preserve the same scalar action.
- **E:** benchmark many moduli and actual curve states across 64–1024 bits; include normalization and chunk setup. Compare stage 1, stage 2, and full factorization, since fewer abstract operations may still cost more in Python.

- [ ] Record the exact a24/formula convention in every optimized kernel and
  checkpoint. The current doubling uses `(A+2)/4` with the squared difference;
  a reference using `(A-2)/4` needs the corresponding formula change. Prove
  fixed-difference normalization assumptions and retain failed-inversion GCDs.

- [ ] Measure selected intermediate reductions and fused arithmetic across
  the below-100-digit modulus range. Current point kernels reduce outputs,
  while intermediates can approach five times the modulus bit width. Compare
  late reduction, selected earlier reduction and whole-ladder fusion on int
  and mpz tracks; extra `%` operations may lose. Include normalization/replay
  and complete stage costs on the M4 rather than transplanting x86/GPU costs.

### P4.3 — Compare PyPy Python integers and optional GMP/`mpz`

**Sequence:** bring this task forward after the Phase 3 relation contracts,
before P7.6 GNFS scaling. Keep one supported runtime, PyPy Python 3.11.
Canonical scalar serialization and backend identity must survive checkpoints;
resume must reject an incompatible backend rather than silently converting it.

- [ ] Introduce a coarse arithmetic boundary with specialized hot loops, avoiding per-multiply virtual dispatch. Keep long-lived GMP values as `mpz`; small indices remain Python integers. Implement backend-consistent GCD, powering, inversion, roots, and exact division contracts.
- **A:** every available backend passes the same result/certainty suite; conversions and missing dependencies have explicit behavior. `divexact` follows a divisibility check; failed inversion retains factor information.
- **E:** compare PyPy Python 3.11 built-in integers and gmpy2 only where available on that runtime. Include import/startup, conversion, schedule, and whole-stage costs. Record unavailable environments without inventing speedup estimates; keep dependency-free and GMP results separate. No CPython support or comparison is required.

**2026-10-04 availability:** the project-local PyPy venv successfully imports
gmpy2 2.3.1 with GMP 6.3.0; the system PyPy has no gmpy2 installation. Both
implement Python 3.11.15 on PyPy 7.3.23. The available venv enables the GMP
comparison track; it does not close backend implementation or performance gates.

- [ ] Resolve the actual backend boundary: point formulas accept mpz, but
  the public ladder's integer validation rejects an mpz modulus. Keep typed
  whole-ladder/stage/product loops with long-lived coordinates/modulus inside
  the backend, and explicit canonical conversion/checkpoint boundaries.
  Validate these paths on the installed ARM64 PyPy/GMP build; a successful
  mpz doubling probe alone establishes neither engine compatibility nor speed.

- [ ] Pin any future PyPy/GMP build and compare specialized whole-stage kernels,
  preserving exact roots, division and certainty semantics. `mpz / mpz` is
  not an exact-integer division contract. Measure bitset XOR separately from
  modular arithmetic and charge all representation conversions. Test any
  experimental `allow_release_gil` use on actual operations before a thread
  arm; its existence alone does not establish useful parallelism.

### P4.4 — Keep Barrett/Montgomery reducers experimental until measured

- [ ] Add optional persistent modulus contexts with consistent encoded identities/parameters. For REDC enforce odd n and `0 ≤ t < n*R`; for the audit's Barrett reference enforce `0 ≤ t < 2**(2*k)`. Make lazy-reduction boundaries explicit.
- **A:** boundary/random arithmetic agrees with native `%`; no oversized intermediate violates correction bounds. Encoded subtraction and GCD checks behave correctly in rho, ECM, and p−1.
- **E:** repeat `audit/reduction_bench.py` with fused real algorithm loops and all available tracks, including setup/conversions. The audit's reducers lost near 166–200 bits: retain native `%` unless the full-run promotion gate overturns that result.

- [ ] Include encoded one (`R % n`), coordinate/parameter conversion and
  canonical exits in reducer oracles. Prove any GCD invariance using the unit
  scaling assumption; a reducer does not authorize dropping nonunit recovery.

**Phase 4 exit:** promote only variants with independently validated arithmetic and reproducible full-stage/portfolio benefit. A documented decision to retain the ladder or native `%` is a successful experiment outcome; a speedup is not guaranteed.

## Phase 5 — Complementary methods and stronger continuations

**Prerequisite:** Phase 2 recovery infrastructure; Phase 3 for relation extensions; Phase 4 chain verification before Lucas PRAC. Source: p+1, ECM continuation, SIQS.

### P5.1 — Implement Williams p+1 with binary Lucas evaluation

- [ ] Add proposed `williams_pp1.py`: GCD-check `A²−4`, use exact prime powers, binary Lucas identities, bounded parameter trials, and stage-1 checkpoints. Add a genuine Lucas stage 2 with accumulated `V_q(V_M(A))−2` terms or validated baby/giant steps; preserve recovery.
- **A:** Lucas identities agree with direct small-index evaluation; smooth p+1 and stage-2-only fixtures produce proper factors. Singular parameters and saturation produce bounded recovery/retry. A Jacobi symbol modulo n never certifies all unknown factors' Legendre symbols.
- **E:** use constructed p+1 and p−1 control cases immediately inside/outside bounds. Compare marginal portfolio completion per CPU-second; schedule p+1 only where its benefit survives its overhead.

- [ ] Specify Lucas composition/doubling and checkpoint parameter identities.
  Test prime-power multiplicities and actual element/group orders, rather than
  assuming that every parameter benefits from a smooth p+1. Keep bounded
  parameter trials and saturated-product replay in both stages.

- [ ] After binary correctness, compare fixed rational starts such as CADO's
  2/7 and 6/5 with bounded seeded starts. GCD-check denominators before
  modular inversion, preserve discriminant/saturation handling, and measure
  conditional order benefits. Repeated p−1 bases are not independent ECM-like
  smooth-order trials; distinguish recovery from additional useful coverage.

### P5.2 — Pair ECM stage-2 primes and tune table size

- [ ] Pair eligible `r−d`/`r+d` candidates using x-coordinate symmetry; cache the union of distances. Keep point tables curve-specific. Sweep D under a memory cap, retaining a positive initialization scalar for the existing recurrence or explicitly redesigning it.
- **A:** paired/unpaired schedules cover the same eligible primes; stage-2 terms and recovered factors validate. Product checkpoints handle early, late, and mixed-factor saturation; tail primes are included.
- **E:** measure baby/giant steps, relation products, replay, setup, and RSS across B1/B2/D grids. Choose D by total continuation cost and factor yield, not simply `isqrt(B2)`.

- [ ] Build an independent eligible-prime oracle for pairing, including
  projective cross-products, D exceptions, initialization and final buckets.
  Tune D using actual prime occupancy and simultaneous table/product/replay
  storage, not an asymptotic square-root estimate alone.

- [ ] Compare wheel/coprime-distance plans with explicit prime-to-term
  coverage certificates, including pruning when an existing `v*w +/- u`
  is divisible by another eligible prime. Cover wheel-divisor exceptions,
  initialization and tails; keep plan construction segmented/bounded rather
  than copying native arrays indexed by absolute B2. Separately evaluate
  no-inversion common-Z baby/giant tables. Verify the homogenized difference
  identity against original cross-products, retain denominator GCDs and
  mixed-factor saturation replay, and charge setup/table/product storage.
  Nonunit scaling over composite n is not projective equivalence.

- [ ] Compile reusable, immutable bound-owned prime-power and paired
  stage-two coverage programs, while keeping curve points private. Current
  production jobs repeatedly request 2,048-wide prime segments, churning the
  eight-entry schedule cache across multi-million B2 runs; cached powers/gaps
  are not consumed by those jobs. Compare bounded program blocks, packed
  schedules and regeneration, charging generation/reads/amortization across
  curves and interrupted resumes. Prioritize this alongside paired products:
  exploratory profiles show schedule generation can rival point arithmetic.
- [ ] Allocate ECM by target factor size and expected marginal success per
  total CPU-second across the sub-100 classes. Allow caller-selected finite
  curve/bound/storage tiers with sufficient work to finish them and cumulative
  extension; current 2,000,000-unit default cannot complete one 50,000/5,000,000
  curve in the existing deterministic work probe (4,081,645 units required).
  Native GMP-ECM curve tables supply hypotheses, not PyPy defaults, guaranteed
  success or evidence that ECM is economical for balanced 90–99-digit n.
- [ ] Size the portfolio's modulus envelope and storage together with tiers.
  All below-100-digit inputs fit 329 bits, while the default reservation uses
  4,096 bits. Exact current-formula probes reserve 9,304,064 bytes for an
  11,000/1,900,000 ECM tier with that default envelope, exceeding 8 MiB before
  execution; a 329-bit campaign reserves 2,665,104 bytes. Larger tiers still
  need explicit storage. These are conservative owned-workspace estimates,
  not RSS measurements; preserve validated bounds and checkpoint overhead.
- [ ] Distinguish cheap automatic ECM pretesting from explicit ECM-only or
  factor-target campaigns, crediting prior completed curves/bounds instead
  of restarting work. Current YAFU and yamaquasi provide concrete allocation
  examples; measure marginal success and sieve handoff on the M4 rather than
  copying their native thresholds. For larger target factors, assess P6.2's
  polynomial continuation when paired classical stage two becomes the
  algorithmic bottleneck, not merely by increasing B2 into a per-prime loop.

### P5.3 — Improve p±1 powering and continuation independently

- [ ] Compare p−1 per-prime powering with prime-power/chunk powering; tune B2/B1 and gap caching. Compare binary Lucas with validated cached Lucas PRAC for p+1. Share integer schedules, deadlines, and recovery tools while retaining distinct group recurrences.
- **A:** stage-1 actions and all stage-2 relations agree with each method's reference. Saturation can recover/retry without lost factors; no p−1 gap multiplication is copied into Lucas code as an ordinary-power update.
- **E:** measure stage-1/2 success gain separately on structured and random corpora. Include schedule/chain overhead and cache hit rate; require a full-run win before expanding default bounds.

- [ ] For a larger resumed B1, apply the ratio of the new and old exact
  prime-power schedules, including increased powers of old primes. For
  example, B1=8 to 16 needs extra factors 2 and 3 as well as new primes.
  Pin starting point/base and schedule extent; ordinary powering, Lucas
  composition and elliptic scalar action retain their distinct recurrences.

### P5.4 — Add double-large-prime SIQS and stronger filtering

**Sequence (M19):** eligible to move forward after working P3.4 SIQS when
profiling/yield justifies it; not a GNFS prerequisite. See the
[research report](quadratic_sieve_research.md#prioritized-experiments).

- [ ] Store residual pairs under explicit factorization/storage budgets;
  combine graph cycles with atomic provenance. Include repeated primes,
  self-loops, duplicate edges and disconnected cycles. Improve duplicate/
  singleton filtering and deterministic unmatched-partial eviction. Keep
  single-large-prime mode as the comparison baseline.
  Separate each large-prime limit from the residual product cap and splitting
  allowance. Document any residual-shape rejection proof; tune bounds jointly
  with thresholds and filtering using training data.
- **A:** cycle combinations reconstruct exact exponents and zero parity;
  residual work, graph size and RSS obey caps. Preserve square corrections
  and relation maps; all outputs pass the existing verifier. Do not assume
  every useful cycle touches the single-large-prime component. Cap/eviction
  cannot leave dangling provenance references.
- **E:** compare useful post-filter dependencies per CPU-second and total
  completion against optimized single-large-prime SIQS and SSS. Record raw
  partials, unmatched occupancy, cycle-space progress, cycle lengths and matrix
  weight; faster collection is insufficient if splitting, filtering or memory
  makes the pipeline slower. Graph cycle counts do not certify a useful
  GF(2) dependency or a proper divisor.

**2026-10-04 cross-reference:** P3.8-R4 owns the current Phase 3 comparison
and dependency-quality diagnostics; this section owns implementation of the
extension. Bound rho/ECM/batch residual splitting explicitly, preserve the
residual prime-certification domain, and retain referenced atoms during graph
eviction. Native triple-large-prime code is not a two-edge DLP template.

**Phase 5 exit:** complementary coverage or completion improves under the promotion policy without losing bounded execution. Adopt methods independently; a p+1 loss does not block a validated SIQS gain.

## Phase 6 — Research options, parallel execution, and claims

**Prerequisite:** a stable, measured portfolio through Phase 5. **Goal:** test higher-cost ideas with explicit stop/go decisions. Sources: coordinate families, continuations, benchmark design.

### P6.1 — Compare Edwards/windowed and torsion-aware ECM

- [ ] Prototype a complete stage-1 package: valid curve families, formula assumptions, signed windows, prime grouping, and table costs. Test a compatible conversion to Montgomery stage 2, including exceptional denominators and correct parameter scaling. Compare against corrected Suyama curves.
- **A:** independent point/map checks pass over prime and composite moduli; nonunits become factors/retries. A formula's exceptional cases are handled explicitly.
- **E:** compare empirical success within fixed work budgets and total time-to-factor across many curves/seeds, including setup/conversion. Promote only a whole-engine win; an operation-count advantage alone cannot pass.

- [ ] Distinguish curve-order torsion from the order of the selected point.
  State each family's congruence and formula conditions; compare success per
  total CPU-second over a distribution of curves and factors. Validate small
  point orders independently before extrapolating stage-one/two smoothness.

- [ ] Specify a mixed Edwards/Montgomery package, as in CADO MISHMASH:
  signed/double-base/precomputed blocks followed by differential Montgomery
  blocks and a compatible stage-two exit. Independently verify chain records,
  coordinate tags/maps and low-order exceptions. Screen the 2024 complete
  Montgomery laws only under their exact finite-field hypotheses; composite-n
  Jacobi symbols do not establish hidden-factor congruence/character conditions.
  Full-coordinate operation counts do not predict x-only ECM performance.

### P6.2 — Investigate polynomial continuations and richer relation collectors

**M19 QS refinement:** scalar exponent recovery remains the reference.
Batch smooth-part/product-remainder trees are eligible when candidate division
dominates, including within P3.5's bounded adapter. Charge tree construction,
node/bit/storage caps, candidate latency and full exponent recovery; a smooth
part or higher batch throughput alone does not close a factorization gate.
Triple-large-prime provenance is not a two-endpoint graph. Introduce a seeded
sparse GF(2) solver only after filtered-matrix/provenance cost justifies it,
with finite retries and original-matrix verification. These options remain
experiments and do not delay GNFS.

**M21 polynomial variant:** after ordinary P3.4 families work, optionally
evaluate Bradford–Monagan–Percival A=A0*q reuse across families. If q is
outside the factor base, explicitly record its exponent in A*F(x) and its
large-prime/provenance role; account for the changed partial and duplicate
policy. Compare against factor-base-smooth A with the same verifier, setup,
storage and extraction budgets. This is not a drop-in root-cache optimization.
For filtering, capped higher-way merges become a candidate only after the
weight-two control; include fill-in and retained provenance in the decision.

- [ ] Profile first, then choose one bottleneck: product/remainder trees, multipoint evaluation, Brent–Suyama extension, triple-large-prime collection, or sparse linear algebra. Prototype behind the validated interfaces. For bigint coefficient packing, prove carry/coefficient bounds and account for memory.
- **A:** continuation coverage or relation identities/dependencies independently validate; peak RSS and residual/matrix work remain bounded. Every prototype can fall back to the baseline.
- **E:** compare one change at a time against Phase 5. Advance sparse solvers only when filtered matrix cost dominates; advance a third large prime only when the full pipeline wins. Stop research that cannot pass the promotion policy.

- [ ] Compare paired classical stage 2 before polynomial/FFT continuations.
  Exact monic product/remainder operations over composite moduli need explicit
  nonunit handling and node/coefficient/storage bounds; floating-point FFT
  needs a proved exact reconstruction contract. Triple-large-prime relations
  require a general incidence/provenance model. Screen recent deterministic
  factoring/high-order work as theory; add a practical challenger only with
  relevant implementation and workload evidence.

### P6.3 — Add reproducible process-level parallelism

Immediate P3.5/P3.6 diagnosis and improvements are owned by P3.6.1 and can
proceed now. This milestone reuses those results for broader ECM/portfolio
parallel execution; it does not postpone their implementation.

- [ ] Assign distinct ECM curves/SIQS polynomial families to workers; share compact immutable schedules or bounded caches. Cancel promptly after a validated split, reconcile pending relations, and checkpoint worker assignments. Use threads only if measured backend operations release the GIL.
- **A:** no duplicated/lost assignments after restart; one validated split cancels remaining work safely; combined results reconstruct n. Aggregate memory and total CPU limits apply across workers.
- **E:** compare 1/2/4 workers under a fixed total workload, reporting wall time, CPU-seconds, setup/IPC, cancellation latency, and aggregate RSS. Compare serial/thread/process modes and both fixed-work throughput and first-valid-factor latency with the single-core baseline; avoid treating extra cores as an algorithmic speedup. Reuse P2.8/P3.6 findings instead of assuming parallelism helps.

- [ ] Make assignment IDs, seeds and committed/in-flight restart state stable
  across worker counts. Specify deterministic merge mode versus race-dependent
  early stopping. Charge all workers' consumed/wasted work, parent deadline,
  retained relations, spill files and cancellation; restarting a worker never
  resets a global allowance. Backend GIL release requires an observed test.

- [ ] Distinguish observed CPU/RSS acceptance gates from active global
  allowances: `parallel_candidates.py` uses per-assignment work and computes
  CPU pass flags afterward. Production workers need parent-owned work leases,
  committed/unspent/cancelled reservation reconciliation, a shared deadline,
  and live/exited-worker CPU accounting. `process_time()` is per-process;
  copying Budget cannot aggregate CPU. Bound duplicated/shared/queued memory
  and spill, disclose cooperative overshoot, and cap nested backend threads.

### P6.4 — Publish reproducible workload-specific results

- [ ] Update `README.md` with pinned configurations, corpus/runner links, certainty semantics, supported backends, measured limits, and timeout behavior. Publish per-workload results against the audit's competitor set: SymPy, primefac, labmath3, PyFactorise, numthy, and SSS, where feasible.
- **A:** every claim links to raw reproducible evidence; full factorization is distinguished from factor-one, and backend/core differences are disclosed. Balanced and unbalanced inputs remain separate; the old 56-digit README example does not stand in for a balanced semiprime.
- **E:** rerun final configurations on held-out inputs with repeated seeds and the declared promotion policy. Claim leadership only for the tested workload/resource class. Include the committed Phase 7 GNFS engine once its correctness and integration gates pass; distinguish a working reference from a performance-promoted configuration.

- [ ] Add separate native reference arms, where feasible: pinned YAFU, FLINT,
  yamaquasi, GMP-ECM and CADO-NFS; identify any historical msieve mirror.
  Disclose architecture, core/GPU count, backend, build and full-factor versus
  factor-one contracts. The September 2026 RSA-260 author report and recent
  CUDA-MPQS results are current research context, with unreproduced timings;
  they do not establish Factor's capability or a PyPy speedup.

**Phase 6 exit:** each experiment has a reproducible adopt/defer/reject decision; published claims match held-out evidence. Deferred experiments remain explicit TODOs with their failed/inconclusive gates recorded.

## Phase 7 — Implement and scale GNFS

**Status (M18): committed; not implemented.** **Prerequisites:** the Phase 2
bounded contracts and Phase 3 relation/filter/dependency infrastructure.
Complete SIQS is the integration/comparison baseline. Sequence P4.3 before
scaling; Phases 4–6 in their entirety are not prerequisites. The engine lives
in `v2/` on PyPy Python 3.11, with optional GMP arithmetic behind the measured
backend boundary. CADO-NFS supplies independently pinned reference results
and comparison data for our implementation.

Sources inspected 2026-10-03: [CADO-NFS's stage overview](https://cado-nfs.gitlabpages.inria.fr/),
[Zimmermann's implementation presentation](https://members.loria.fr/PZimmermann/talks/cado.pdf),
and [Bai, Brent and Thomé on polynomial root optimization](https://arxiv.org/abs/1212.1958).
The dependency order and acceptance gates below are Factor engineering
decisions; native reference timings are not predictions for PyPy.

### P7.1 — Establish exact polynomial and relation contracts

- [ ] Add GNFS polynomial/field and relation structures. Begin with bounded
  base-m polynomial selection for general inputs, a linear rational polynomial,
  and a common root modulo n. Verify irreducibility and nondegeneracy for the
  supported field representation. Store both homogeneous norms, signs,
  rational prime powers, algebraic ideal identities, and bad-prime metadata.
- **A:** independently verify common-root and homogeneous-evaluation identities
  using exact integers. Handle leading coefficients, ramified/bad primes and
  noninvertible denominators explicitly; return a proper factor or bounded
  retry. Equal norm primes cannot erase distinct algebraic ideals.
- **E:** use small general composites and independent polynomial/ideal oracles;
  include negative norms, primes dividing coefficients/discriminants and
  malformed relations. Special-form inputs alone do not establish GNFS.

- [ ] Specify nonmonic normalization exactly: for degree d and leading
  coefficient f_d, `F(a,b)=b**d*f(a/b)=f_d*Norm(a-b*alpha)`.
  Record any monic scaled generator/basis and denominator corrections. Ideal
  identities include side, prime, affine/projective root and bad-prime branch
  data; norm factorization alone does not identify every ideal valuation.
  A restricted first reference must explicitly bound/reject unsupported bad
  primes or fields and retain finite polynomial retry/factor recovery.

### P7.2 — Collect bounded rational/algebraic relations

- [ ] Start with a serial line-sieve reference over reproducible primitive
  pairs `(a, b)`. Factor both norms under bounded cofactor work and retain
  only verified full relations initially. Persist compact verified relations,
  collection position, and remaining budgets with versioned checkpoints.
- **A:** every stored relation reconstructs both norms and satisfies its
  ideal/root conditions. Duplicate detection, cancellation, storage caps and
  restart preserve provenance. Smoothness scores only select candidates;
  exact verification decides whether a relation is valid.
- **E:** compare sieve collection with exhaustive small-region enumeration.
  Measure useful relations, cofactor cost, CPU, RSS and disk consumption;
  validate pause/resume on identical assigned regions.

- [ ] Retain primitive-pair, sign, homogeneous-value and known special-q
  factors in exact relation checks. Tie spills to polynomial/ideal-numbering
  identities with bounded I/O and interrupted-write recovery. Establish the
  full-relation reference before tuning two-sided residual cofactoring.

- [ ] Canonicalize primitive `(a,b)` identity with a declared sign/b=0 policy
  and field/store identity; deduplicate retries and overlapping discoveries
  exactly before useful-yield measurements. Optional online suppression needs
  verified earlier assignment geometry, thresholds and cofactor policy;
  a smaller special-q dividing a norm does not prove prior discovery.
- [ ] Treat raw collection targets as filter triggers, then use deduplicated
  rows, live ideal columns, excess and verified dependencies to request
  further finite assignments when needed. Preserve the store/cursor between
  collection/filter rounds; fixtures with many duplicates/singletons must not
  confuse raw relation count with completion or finite-region exhaustion
  with algorithm failure. Final square-root/proper-divisor checks still apply.

### P7.3 — Filter and solve dependencies with GNFS constraints

- [ ] Reuse storage/provenance infrastructure, with distinct rational-prime
  and algebraic-ideal columns plus required sign/character constraints.
  Begin with exact Python-bitset elimination on small matrices and preserve
  the map from filtered rows back to original relations.
- **A:** independently recheck zero parity and all required character data
  against original relations. Ideal parity alone is not proof that the
  algebraic product is a square. Duplicate/singleton removal cannot lose
  exponent or dependency provenance.
- **E:** use known-dependency matrices, corrupted columns/characters and
  duplicate relations. Compare small results against an independent solver.
  Algebraic square-root checks in P7.4 remain mandatory after matrix checks.

- [ ] Choose explicit character placement: full-matrix constraints or a
  bounded correction solve inside the provisional kernel span, as in CADO.
  Reimpose any omitted heavy constraints with exact lifting, reject zero/
  dependent vectors, and verify originals. Finite character tests screen
  candidates; they do not prove the algebraic element is a square. Factoring
  characters and discrete-log Schirokauer maps have distinct contracts.

### P7.4 — Compute rational and algebraic square roots and split n

- [ ] Implement exact rational-root construction and an algebraic square-root
  algorithm with documented field representation, coefficient/denominator
  bounds, lifting or reconstruction, and correction factors. Map both roots
  through the common root modulo n; try both `gcd(X - Y, n)` and
  `gcd(X + Y, n)`. Retry trivial dependencies within finite allowances.
- **A:** certify the algebraic square-root identity in its field representation
  and independently check `X**2 % n == Y**2 % n`. Validate every returned
  split. Integer square roots of algebraic norms do not replace this step;
  zero parity, a large relation count, or a completed matrix is not success.
- **E:** complete small general composites through our entire GNFS pipeline.
  Test nonmonic polynomials, denominator failures and trivial congruences;
  compare root identities and factors with pinned independent references.

- [ ] Choose the reference square-root algorithm before field coverage grows.
  Inert-prime lifting needs irreducibility modulo an auxiliary prime, and
  some irreducible fields have no such prime. Bound that search and declare
  supported fields, a validated alternative or finite refusal/retry; CADO's
  pinned implementation also caps its search. Test a no-inert-prime quartic,
  such as `x**4-10*x**2+1`, nonmonic scaling and bad denominators. Use justified
  coefficient bounds and independently verify the reconstructed field identity
  and modular roots.
  A CRT alternative must reconcile root signs consistently. The 2023 odd-
  prime-power e-th-root paper is not an automatic e=2 implementation upgrade.

- [ ] A CRT alternative needs explicit finite split-prime search, precision
  growth and root-sign reconstruction, plus rational-root integration and
  exact field verification. CADO's separate CRT program is a reference with
  manual/integration limitations, not an already integrated general fallback.
  Heuristic coefficient estimates need checked reconstruction and bounded
  precision growth; test insufficient precision and failed sign recovery.

### P7.5 — Integrate the small GNFS engine with bounded dispatch

- [ ] Add GNFS as an explicit opt-in stage after the measured SIQS baseline,
  carrying one shared allowance across polynomial selection, collection,
  filtering, matrix work and square roots. Version checkpoints and immutable
  relation-store identities; retain unresolved cofactors on exhaustion.
- **A:** complete and factor-one outputs preserve certainty and multiplicity.
  Resume with the same assignments reproduces work/results; cancellation or
  corrupt/stale relation stores cannot silently discard factors or work.
- **E:** compare full pipeline runs with interrupted/resumed runs on small
  held-out general composites. Include setup, import, conversion, relation
  loading and output costs. Keep GNFS opt-in until P7.7 promotion passes.

- [ ] Checkpoint polynomials, field/basis, norm convention, ideal numbering,
  character policy, dependency lifting and backend identity together with all
  stage allowances. Reuse SIQS storage machinery through explicit interfaces;
  its relation payload cannot stand in for GNFS ideal/field data.

- [ ] Separate immutable field/ideal/store identities from extendable
  work/time/storage quotas and append-only assigned regions. Compare a single
  larger finite run with explicit quota/range extension through pause/resume,
  retaining prior consumption and verified relations. Changed bases, numbering
  or sieve/cofactor policies need versioned preserve/remap/reverify or charged
  restart; increasing a resource quota alone must not discard the store.

### P7.6 — Scale polynomial selection, sieving and matrix work

- [ ] Improve degree/skew/root-quality selection against the base-m control.
  Add bounded special-q lattice sieving and large-prime relation handling
  only after the reference verifier passes. Introduce sparse GF(2) solvers
  such as block Wiedemann or block Lanczos when matrix cost justifies them.
  Specify bounded disk spill, restart and total CPU/RSS accounting.
- **A:** optimized collectors/solvers preserve exact relation and dependency
  verification. Recheck sparse solutions against the original matrix; bound
  retries. Backend conversion, partial relations, buckets and matrix storage
  obey declared resource limits on the supported PyPy runtime.
- **E:** profile first, then compare one stage change at a time including
  polynomial-selection cost and final factor extraction. Compare Python-int
  and available GMP paths under identical inputs/budgets. Parallel scaling
  follows P6.3 accounting and is not a prerequisite for the small engine.

- [ ] Shortlist polynomials by measured size/root/skew estimates (including
  Murphy E and optional E') and bounded trial sieving; verify common roots
  after every rotation/translation and charge search cost. Prove special-q
  lattice mappings and preserve forced ideal exponents. Compare two-sided
  cofactor strategies jointly: medium-prime sieve, batch small-prime removal,
  staged tests, first-side choice and bounded ECM. The 2023 alternative-sieving
  study motivates this experiment; local collection gains with relation loss
  require full-pipeline confirmation. Keep native/GPU presets as references.

- [ ] Retain general side-labelled sparse ideal incidence for relations
  containing more than two large ideals; reuse a QS edge/cycle collector only
  under a proved restriction. Keep per-prime lpb and whole-residual mfb domains
  separate, with explicit accepted-cofactor certainty. Test repeated powers,
  three-plus large ideals, equal primes on distinct sides and different roots
  above one prime against exact elimination and full norm reconstruction.
- [ ] Introduce prime special-q before composite special-q and prove lattice
  determinant/congruence and inverse-coordinate mappings. Expand optional
  duplicate suppression only after overlapping/retried tasks, both-side
  ranges, projective roots and changed geometry pass independent coverage
  controls; native probabilistic suppression can lose useful relations.

### P7.7 — Measure coverage and the SIQS/GNFS crossover

- [ ] Freeze independent training/held-out general-composite bands. Grow from
  completed small fixtures into capped 60/70/80/90/100-digit exploration,
  with larger bands only under explicit budgets. Pin CADO-NFS for comparison
  and disclose its native backend and resource model separately from Factor.
- **A:** all successful outputs reconstruct inputs with correct certainty;
  failed/exhausted runs remain censored outcomes. A small correct GNFS engine
  does not imply practical coverage at every proposed larger band.
- **E:** report repeated seeds, completion, median/spread, CPU, aggregate RSS,
  disk and cold/warm costs for complete factoring. Apply the promotion policy
  before choosing a default crossover. Retain SIQS wherever GNFS does not win;
  record limits without abandoning the committed GNFS workstream.

- [ ] Separate general GNFS inputs from SNFS-friendly forms and report
  capability, completion and dispatch decisions independently. Neither an
  asymptotic L-notation comparison nor a published GPU record fixes a usable
  SIQS/GNFS digit crossover for bounded PyPy runs.

**Phase 7 milestones:** P7.1–P7.5 deliver a correct bounded small GNFS engine.
P7.6 delivers a validated scaling pipeline. P7.7 supplies a measured dispatch
decision and declared coverage. Report these exits separately; none is done
yet, and no calendar estimate or unmeasured digit threshold is claimed.

## Phase 8 — Revisit Phase 2 portfolio optimization

**Status (M23): future research/experiment backlog; not implemented.**
**Prerequisites:** the accepted M17 bounded contracts; working P3.4 SIQS and
P7.5 small GNFS for integrated portfolio decisions. Profiling or independent
preprocessing/rho experiments may proceed earlier when useful. This phase
does not postpone GNFS correctness or scaling and does not require every
optional P4–P6 optimization. **Goal:** resolve remaining Phase 2 tuning with
new evidence, retaining completed core work and earlier adopt/defer/reject
decisions. Sources and transfer limits are in the
[research report](phase_two_optimization_research.md).

Brent batches, prime-exponent power detection, streamed stages, chunk replay,
bounded caches, JSON checkpoints, M13 local loops and M14 worker feasibility
already exist. Preserve the [M17 baseline](phase_two_m17_frozen_baseline.json).
Keep historical P2 gates open until their declared experiments pass; P8 IDs
identify future execution and do not certify those gates by renumbering them.

### P8.1 — Freeze the new control and untouched evaluation data

- [ ] Reuse the runner/oracles and preserve M12–M17 captures. Freeze the
  current SIQS/small-GNFS portfolio, sources and configuration. Build fresh
  training and held-out inputs with independent certificates hidden from
  algorithms; M15/M17 confirmation data already informed this review.
  Profile stage/setup/recovery costs separately from timing. Execute the
  existing pinned competitor arms where feasible, with backend/core and
  output contracts disclosed. P6.4 owns publication of the final comparison.
- **A:** every successful output reconstructs its input with correct
  certainty/multiplicity. Record failed/exhausted/unavailable arms, source
  hashes, caps, seeds, JIT warmup/stability and cold startup. Do not replace
  the frozen M17 control or claim a ranking from a rho-only/native utility.
- **E:** repeat both complete and factor-one modes across balanced,
  unbalanced, powers, primes/pseudoprimes, p±1 boundaries and close/random
  cases. Start feasible cells with at least twenty inputs and five seeds,
  expand uncertainty, and report all-outcome medians/spread, completion,
  censored exhaustion, CPU and RSS. Keep tuning and final evaluation separate.

### P8.2 — Reduce preprocessing work with exact rejection proofs

- [ ] Tune the existing trial cutoff on training data. Compare a proven
  factor lower bound from completed trial division to reduce power exponents
  using `L**k <= n`; use exact comparisons and preserve progress proof on
  splits/resume. Compare small modular power-rejection filters before Newton
  roots, retaining final exact equality. For optional Fermat, compare one-time
  initialization and `D(a + 1) = D(a) + 2*a + 1` with square-residue rejection.
  Consider an extra Lucas/BPSW composite filter only if classification cost
  dominates, preserving configured Miller–Rabin rounds and certainty.
- **A:** no power or factor can be rejected incorrectly, including neighboring
  powers, mixed multiplicities, partial trial progress and resumed work.
  Filter passes are not proofs; arbitrary-size inputs never enter a float.
  Fermat remains finite; every split validates. A root-algorithm replacement
  needs a separate profiled justification and the same exact oracle.
- **E:** compare 1,000/5,000/25,000/100,000 trial cutoffs as candidates,
  with power-rich/random/adversarial-filter inputs and full-run held-out
  costs. Include filter setup/cache and power-loop work; leave Fermat off
  unless declared close-factor coverage justifies its portfolio cost.

### P8.3 — Calibrate bounded rho batches, walks and restarts

- [ ] Reuse the current Brent/local-loop control. Tune batches, walk allowance
  and restart allocation jointly after individual comparisons. Record GCD,
  modular-product, cycle-advance and saturation-recovery costs. Keep fixed
  assigned seeds and finite walk/replay allowances. Reducer implementation
  remains P4.4; process search remains P6.3.
- **A:** failure/recovery cannot hang, return `n` as a split, skip consumed
  work or lose resume identity. Force saturation and failed early attempts;
  include tail batches and cancellation. Native word-size overflow tricks
  cannot enter the arbitrary-size contract; exhausted walks are not primality
  proofs or successful factoring baselines.
- **E:** compare the existing 32/64/128/256 batch candidates and bounded
  restart policies under identical total budgets. Measure first-factor and
  complete runs, replay incidence and cancellation latency. Brent's paper,
  FLINT's batch 100 and Algorithmica's batch 1,024 do not establish a PyPy
  optimum. Retain 64 when differences are inconclusive.

### P8.4 — Tune ECM/p−1 allocation and relation-engine handoff

- [ ] Fit joint B1/B2/curve grids using measured PyPy stage costs and factor
  yield, with factor-size bands hidden from algorithms. Select a deterministic
  policy from observable input/configuration and completed-work metadata.
  Compare shorter preprocessing with earlier SIQS/GNFS handoff, including
  recursion/setup. Avoid duplicate p−1 base/bound work; any incremental
  prime-power or same-curve continuation implementation remains P5.3.
  P3.4/P7.7 retain dispatch/crossover ownership.
- **A:** one allowance spans all stages/children; policy identity and prior
  CPU/wall/work expenditure survive checkpoints. Model estimates do not
  certify hidden factor sizes or primality. Incremental B1 requires missing
  exponent ratios for old primes too; prove action before testing promotion.
- **E:** use P2.2's B1 2,000/11,000/50,000 as exploratory candidates and fit
  B2/counts jointly under memory caps. Compare marginal completion per
  CPU-second and full-run outcomes, conditioned on remaining cofactors.
  Include structured/random primes and relevant residue classes; neither
  GMP-ECM tables nor YAFU/FLINT native thresholds are production defaults.

### P8.5 — Measure recovery, polling and checkpoint costs separately

- [ ] Distinguish GCD/chunk polling, cooperative budget checks/atomic commits,
  caller durable checkpoint writes, and final JSON pack/unpack verification.
  Sweep existing chunks/batches and resumptions. Consider bounded product-tree
  saturation recovery only if linear replay dominates; include retained terms,
  reconstruction, node storage and fallback. The library does not write a
  checkpoint file every arithmetic iteration.
- **A:** preserve finite recovery, canonical/corrupt-state checks, exact
  reconstruction and consumed allowances. Document indivisible bigint calls
  and cancellation/deadline overshoot. A leaf GCD equal to `n` is not a proper
  factor. Version changed checkpoint or accounting semantics; do not remove
  validation to make resume look faster.
- **E:** compare chunk candidates 1/4/8/16/32/64 and stage-2 GCD batches on
  early/middle/tail and mixed-factor saturation. Include no-saturation controls,
  cold/warm pause/resume and JSON consumption. Measure CPU, RSS, replay and
  interruption latency. Changed work units need disclosure and actual wall/CPU
  comparisons, not numerical equivalence of different accounting schemes.

### P8.6 — Revisit context setup and schedules only for changed workloads

- [ ] Profile context construction before early classification/rho success.
  Compare lazy creation or staged bounds if wasted setup matters, preserving
  upfront configuration/cap checks and one-time charges. Reuse private buffers
  and capped integer prime/power/gap schedules; curve points and modular
  residues remain job-specific. Reuse the existing cache/rolling/wheel/packed
  arms only if new bounds, backend or reuse makes them relevant. Algorithmic
  stage-2 pairing/continuations remain P5.2–P5.3/P6.2.
- **A:** half-open complete ordered prime streams, arbitrary restarts and
  private ownership survive; no packing truncates values or cache reuses a
  different modulus's state. Memory reservations and resume/work identity
  remain valid. Record why M12's retained decision is being revisited.
- **E:** measure setup, marking, extraction, consumption/decoding and full
  factoring in cold and amortized regimes, with cache hit/miss and bounded
  RAM/disk costs. Include process/JIT RSS separately from owned workspace.
  Keep wheel-6/bytearray, cache off and rolling off without full-run promotion;
  use the [C transfer review](sieve_port_review.md), not its native cutoffs.

### P8.7 — Integrate accepted winners and close the reconciliation

- [ ] Integrate only independently accepted P8 variants and relevant P4–P6
  winners. Keep original implementation ownership for PRAC, kernels, GMP,
  reducers, p+1, pairing, curve families, polynomial continuations and workers.
  Record adopt/retain/defer/reject decisions against each historical P2 gate;
  update completion only when its acceptance/experiment evidence supports it.
- **A:** final sources/configuration are pinned and preserve result/certainty,
  caps, cancellation and canonical resume. Optional losses do not block other
  accepted work. M14 fixed-work throughput alone cannot enable workers;
  P6.3's first-factor, CPU and aggregate-RSS gate remains required.
- **E:** confirm final settings on P8.1's untouched inputs with repeated
  seeds and the shared promotion policy. Include cold import/setup, backend
  conversion, recovery, output consumption and complete factoring. Retain
  serial/dependency-free defaults where promoted evidence is absent, and
  disclose bounded opt-in versus any proposed API default change separately.

**Phase 8 exit:** every selected optimization has evidence or an explicit
retain/defer/reject disposition, with historical P2 tuning gates reconciled
honestly. No speed claim, new default, competitor execution or additional
implementation completion is established by this M23 research milestone.

## Phase checkpoints

| Phase | Concrete deliverable | Exit decision |
| --- | --- | --- |
| 1 | Native Python 3 repairs, result/certainty contracts, regression suite | All counterexamples and failure paths pass; freeze corrected baseline |
| 2 | Budgeted portfolio, bounded schedules, recovery/resume, corpus/runner | Time/memory limits and reproducibility demonstrated |
| 3 | Verified relation pipeline, SIQS, SSS comparison | Balanced-composite completion improves under equal limits |
| 4 | Validated PRAC/kernel/backend/reducer experiments | Promote measured winners; retain baseline for losses |
| 5 | p+1, paired/recoverable continuations, double-large-prime SIQS | Complementary completion gain without resource regression |
| 6 | Curve-family/research decisions, fixed-core workers, evidence-backed README | Held-out results support every adopted default and claim |
| 7 | Bounded GNFS reference, scaling pipeline, SIQS/GNFS crossover | Separate correctness, scaling and held-out dispatch gates |
| 8 | Reconciled Phase 2 calibration, preprocessing and overhead experiments | Fresh held-out promotion or explicit retain/defer/reject decisions |

Execution priority after Phase 3 is early P4.3 arithmetic work followed by
Phase 7. Existing phase numbers preserve references; optional Phases 4–6
experiments and future Phase 8 calibration do not postpone the GNFS pipeline.

## Existing evidence to reuse

Read the scripts before running them: the compatibility scripts load the original Python 2 modules and print diagnostic observations; they are not acceptance tests for a finished Python 3 port. Some depend on `lib2to3` and the original absolute source path. Build native assertions for Phase 1 and store new results separately.

- [validate.py](validate.py) and validation.json: setup, sieve, splitter, API, and injected-failure counterexamples.
- [extra_checks.py](extra_checks.py) and extra-validation.json: endpoints, independent affine oracle, PRAC failures/timeouts.
- [end_to_end_checks.py](end_to_end_checks.py) and end-to-end-validation.json: small reconstruction sweep, seed search, point timing diagnostics.
- [prac_reference.py](prac_reference.py), [test_prac_reference.py](test_prac_reference.py), and prac-validation.json: guarded chain prototype and separately counted exceptional projective pairs.
- [reduction_bench.py](reduction_bench.py) and reduction-validation.json: arithmetic checks and kernel measurements with disclosed exclusions.

The full audit report supplies research references and their limitations. Pin any external implementation before benchmarking it; re-check its contracts before adapting code.
