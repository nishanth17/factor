# C1 bounded DLP implementation and complete-factor evaluation

Status: the bounded C1 implementation and calibrated-SLP comparison are
complete. Adopt the explicit smaller-base DLP bundle for the tested balanced
40-digit class. Preserve defaults and the larger unresolved outcomes. Final
committed-files-only validation passes 526 tests, lint and required loaders.

## Why implementation proceeded

The [research review](c1_research.md) verifies all-component cycle mathematics
against primary literature and distinguishes pinned implementation policy from
mathematical guarantees. YAFU/msieve/Yamaquasi provide residual-bound and
smaller-base hypotheses, with native cutoffs explicitly excluded as PyPy
predictions. FLINT and JavaMath provide further certainty/ownership comparisons;
no external implementation code was adapted. Licenses and retrieved source
identities are recorded in that review.

The initial512-block screen could not decide DLP economics. Its final defer
was withdrawn. The longer bounded 40/50/60-digit study retained rejected
residuals and sampled positions as well as survivors. Both40- and both 50-digit
inputs produced independently verified DLP complete factors at prefixes where
SLP incidence had no dependency. Entire-census splitting/certification cost
was under 1 s per 40-digit input and under 3.3 s per 50-digit input, instrumented.
Those figures exclude collection and are not accepted performance evidence.
At 60 digits, both runs remained unresolved: many degree-one vertices survived,
and every LP-cancelled row was singleton-filtered. Raising a solver budget
would not repair those matrices.

The original literal follow-up gate charged both competing policy analyses,
SLP analysis and exhaustive extraction after a factor. It failed and remains
preserved. The explicitly post-observation attribution repair used retained
records only, stayed inside the original 4,500-second envelope, and charged
all processing for one candidate through its first complete factor. That
supports an investment go at 40/50, not a retrospective preregistered pass or
a factoring-speed claim. Full data and limitations are in
[the feasibility report](c1_followup_results.md).

## Implementation and acceptance

Production source `ebd049e` adds explicit `DoubleLargeSieveConfig`, separately
bounded proven endpoint pairs, bounded rho/square splitting, and a finite
all-component forest. The scalar residual contract remains unit-or-proven-prime.
Every combined row retains original atoms, exact factor-base exponents and
square corrections. Repeated primes, self-loops and parallel/disconnected
cycles are supported. A 256-atom cycle limit is an explicit yield loss.

Only unowned forest edges can be evicted. Planning is read-only, followed by
verification and a final resource check before publication. Owned paths survive
rebuilding. SIQS checkpoint version 4 reconstructs the forest and verifies every row,
closing edge, atom reference and polynomial ownership, preserving original
mixed row order. Prior work/wall/CPU, replay and started split attempts stay
charged. Ordinary SLP uses checkpoint version 3. The same explicit SIQS
configuration works through portfolio resume; no shared portfolio/stage
schema changed.

Independent tests enumerate the entire boolean incidence kernel of small
multigraphs and compare its span with emitted cycles, including eviction.
Exact N=91 fixtures expose proper 7/13 splits in a disconnected self-loop and
a disconnected parallel-edge cycle whose known square shares a large prime.
Further checks cover all four division backends, conservative and deliberately
narrowed candidate domains, corrupted corrections/stores, cap refusals,
cancellation, native/GMP arithmetic, and charged nested portfolio resume.

The final committed-only `40d8d71` archive passed **526 PyPy/GMP tests** in
80.333 seconds using `make -C v2 test` with the explicit PyPy 3.11 interpreter.
`make -C v2 lint` passed Ruff checks/formatting and pycodestyle for all 205
Python files. All **109 benchmark modules** imported, and required C1
controls, the owned SLP baseline, independent corpus certificates, selected
configurations and frozen source hashes loaded and verified. The receipt is
local `results/c1/final-committed-1.log`; the committed-only archive is retained
beside it. Later receipt edits are Markdown only, with tested runtime, test,
control and corpus bytes unchanged. Generated evidence is not committed.

## Frozen evaluation

The [implementation protocol](c1_implementation_protocol.md) fixes five arms:
committed calibrated SLP, graph-only SLP reuse, full DLP, narrowed DLP candidate
allowance, and DLP with half the numeric base bound. The [control notes](c1_implementation_controls.md)
explain the separate owned SLP source package, finite graph reservations and
charging. Original 40-digit input IDs collide across the two source corpora;
the `run-40-0-*` and `run-40-1-*` filenames and independently certified factors
identify the two distinct cases.

Training selects by completion, then total capped cost. Settings and source hashes froze in `a008760` before generation of the new
30/40/50-digit inputs. Confirmation includes at
least 3 s validated PyPy 3.11 warmup and 9 paired samples, extending to 18/27 on the
prespecified instability checks, within 10,800 wall/CPU seconds. All setup,
collection, splitting, certification, graph/eviction, provenance, filtering,
matrix, extraction, classification and final output validation are timed.
Cold processes and profiles are separate. RSS is the process high-water mark,
including prior arms in the same process; owned reservations are per-job.

## Cross-owner comparison boundary

A7's integrated reconciliation (`018236f`, read from mainline without changing
this branch) prepares calibrated 30-digit SSS comparison arms for E1. It
explicitly identifies the absence of a calibrated 40-digit SSS arm and assigns
that training/confirmation to E1 or an explicit challenger-class deferral.
C1's frozen five-arm study therefore does not invent larger SSS settings or
claim a DLP-versus-SSS crossover. P5.4's broader SSS comparison remains an E1
integration prerequisite. B3 released the performance window at 07:12:46 UTC;
C1's committed-only QA and subsequent study each hold the machine-wide lock.
No B3/A7 algorithm or shared portfolio/stage-job file was edited.

## Completed training and bounded profile

The fixed 20-run training matrix completed in 1,930.186 wall / 1,838.703 CPU
seconds. Every arm completed both 40-digit inputs. Total elapsed seconds
across those two single-run diagnostics were 20.160 for SLP, 20.055 for
graph-SLP, 15.492 for full DLP, 15.097 for narrowed DLP and 13.518 for the
smaller-base DLP arm. These exploratory values select a policy; they are
not accepted speed estimates. The half-base arm incurred no eviction on
either 40-digit input and retained less workspace than full DLP.

All five arms exhausted 120 seconds on both 50-digit training inputs. On the
first input, full/narrowed DLP and SLP filtered to zero. The half-base DLP
arm retained a nonempty 466-row / 711-column filtered matrix but no kernel.
It was closer to a useful matrix, not a complete-factor winner. Equal capped
failure costs leave the 50-digit DLP selection tied; the already implemented
arm-order tie break selects full DLP. No row-count surrogate changes that
predeclared ranking.

The one diagnostic profile used 120.019 wall / 119.439 CPU seconds, including
60 seconds of unprofiled collection followed by a profiled tail under the
same cumulative allowance. Sieve/root work dominates the profiled tail;
graph admission is small by comparison. The profile does not justify a
DLP-specific repair or a residual-certainty shortcut. Its function timings
are instrumentation, not accepted throughput estimates. No implementation
repair or repeated training matrix is undertaken. Broader sieve/root reuse
remains the existing R2/C8 workstream.

## Large-store charged resume

The separately frozen [resume control](inputs/controls/c1_resume_acceptance.json)
uses the already inspected first 50-digit input, half-base DLP, 512 MiB owned
allowance, 8 MiB checkpoint allowance, and 300 seconds cumulative wall/CPU.
These larger allowances are not substituted into timing arms. A pause at
60 seconds was serialized, released, read from JSON and restored with prior
work/wall/CPU charged. The checkpoint occupied 1,441,315 bytes; reconstruction
charged 32,155,909 additional work units.

The resumed run completed in 141.063 total wall / 140.888 CPU seconds, returning
6,234,745,488,095,579,115,919,703 × 9,192,770,129,255,987,774,468,263. Both
runtime labels remain probable-prime; independent corpus certificates prove
the expected factors. It retained 8,795 forest edges, including 889 owned edges,
after 20,480 evictions. Every returned divisor was proper and the outcome
reconstructed the original input. This establishes a concrete larger-input
production/resume witness, not a matched 50-digit speedup.

Training, profile and this acceptance check together used 2,191.268 active
wall seconds, within the 4,500-second envelope. Source/setting freeze and fresh
confirmation then followed. The [pre-generation seed-stability clarification](c1_seed_stability.md)
preserves pooled extension triggers and predeclares conditional fixed-seed
checks without additional runs or relaxed correctness/resource gates.


## Separating threshold, policy and duration losses

The retained census classifies every recoverable product-128 DLP residual by
its first exclusion from calibrated SLP. These are counts over the collected
prefix, not success-conditioned samples:

| Input | Recoverable DLP | Block-score loss | Refinement loss | Residual-policy loss |
| --- | ---: | ---: | ---: | ---: |
| 40/0 | 7,432 | 3,855 | 1,137 | 2,440 |
| 40/1 | 8,687 | 5,094 | 1,243 | 2,350 |
| 50/0 | 20,076 | 9,004 | 2,544 | 8,528 |
| 50/1 | 20,173 | 7,707 | 2,477 | 9,989 |
| 60/0 | 12,179 | 4,272 | 1,420 | 6,487 |
| 60/1 | 13,699 | 4,828 | 1,689 | 7,182 |

Thus changing only the residual policy would miss much of the opportunity.
The full DLP candidate allowance covers the declared product domain; the
narrowed training arm explicitly measures the alternative cost/loss tradeoff.
It did not supply a complete 50-digit training factor within 120 seconds.

Independent trial division cross-checks 32,768 randomly selected positions
per input, including rejected positions. Among SLP-threshold rejects, residual
bit-length quartiles were 63/70/76 and 63/69/75 at 40 digits, 76/84/90 and
77/84/91 at 50 digits, and 93/101/109 and 91/100/108 at 60 digits. Only
84/59, 33/17 and 4/7 sampled positions, respectively, passed both SLP scoring
stages. Every sampled residual within the widened product cap passed its
conservative DLP scoring checks. These finite samples do not prove the
unseen distribution, nor cover positions after the approximately 8.4-million
sampling cutoff. Large rejected residuals are not presumed semiprimes.

Duration and storage also matter. At the original 512-block prefix, both
40-digit DLP matrices vanished under singleton filtering; later prefixes
produced proper factors. The production 50-digit resume witness completes
only after the smaller-base collector's 120-second training window, with a
separately larger checkpoint/memory allowance. It demonstrates a viable
longer run, without isolating duration from the allowance change. Conversely,
the two 60-digit terminal graphs still had more than 90% degree-one vertices
and no post-filter matrix. More matrix-solving work alone cannot help those
particular captures. Neither observation establishes a universal digit cutoff.


## Fresh complete-factor confirmation

Runtime source remains `ebd049e`; selection, runner and fresh corpus were
committed in `a008760` before accepted timing. The machine was an Apple M4
with 24 GiB RAM, macOS 26.6.2 arm64, PyPy 7.3.23 implementing Python 3.11.15.
The frozen run completed in 5,536.298 wall seconds, below its 10,800-second
allowance. The parent used 5,028.868 CPU seconds; adding the separately timed
cold-child factoring work gives 5,281.032 CPU seconds. That subtotal excludes
cold-child interpreter/import CPU; cold process wall time is retained below.
All warmups were validated and lasted at least three seconds per case/arm.

The six new balanced inputs were disjoint from 2,048 previously versioned
inputs. Their Pocklington-conditioned generation population is disclosed;
this is not an RSA-distributed or broad general-composite population. The
algorithm receives n and its public configuration, never expected factors
or certificates. There are only two independent inputs per band, regardless
of the number of repeated runs.

| Fresh input | Samples per arm | SLP mean / median seconds | DLP mean / median seconds | Complete SLP / DLP |
| --- | ---: | ---: | ---: | --- |
| 30/0 | 27 | 0.4471 / 0.4452 | 0.3937 / 0.3919 | 27/27 / 27/27 |
| 30/1 | 9 | 0.6968 / 0.6826 | 0.4943 / 0.4857 | 9/9 / 9/9 |
| 40/0 | 9 | 10.9993 / 11.1321 | 7.7430 / 7.6862 | 9/9 / 9/9 |
| 40/1 | 9 | 12.4706 / 12.3636 | 8.2804 / 8.2175 | 9/9 / 9/9 |
| 50/0 | 9 | 120 capped / 120.0004 actual median | 120 capped / 120.0005 actual median | 0/9 / 0/9 |
| 50/1 | 9 | 120 capped / 120.0006 actual median | 120 capped / 120.0005 actual median | 0/9 / 0/9 |

DLP means the half-base bundle at 30/40 and the full-base bundle at 50.
The predefined capped-total-cost statistic reduces by 17.80% at 30 digits
(95% paired bootstrap interval 10.50–30.35%) and 31.73% at 40 digits
(29.21–33.99%). At 50 digits both capped-cost benefit and completion gain
are zero. The bootstrap resamples inputs and paired sample indices 10,000
times. Sample counts differ at 30 digits because the original stability rule
extended input 0 to 27; no sample was removed or reweighted after observation.
Completion is unchanged in every band. All 144 timed outcomes reconstruct.

Every cell passes the original pooled IQR/drift rule. The separately
predeclared fixed-seed analysis nevertheless exposes residual drift on
30/0 DLP: first/last-third shifts are 15.04%, 12.05% and 13.44% for seeds
7/29/47, with seed 7 also having 13.88% relative IQR. The original pooled
verdict is preserved; no conditional rescue rule is needed or used. This
extra diagnostic cautions against treating the 30-digit timing as a stable
per-seed latency estimate. Retain its positive aggregate observation and
zero completion regression, but make no new 30-digit preset promotion.
The finite 27-sample ceiling is respected. At 40 digits, the original pooled
checks pass at nine samples, and benefit is positive for every seed across
both inputs (31.54%, 33.23%, 30.38%).

The complete fresh factorizations, all runtime-labelled proven-prime, are:

| Input | Independently checked factors |
| --- | --- |
| 30/0 | 731970105903419 × 909211531732181 |
| 30/1 | 580764872393003 × 801667028477113 |
| 40/0 | 38889112373019322867 × 50447243852189729833 |
| 40/1 | 50345498199009279587 × 55469174876933793883 |

All 18 samples per arm at 50 digits end with `wall_limit`, no divisor, and
these complete unresolved cofactors retained:

- 50/0: `62991431486759765271261990605422558285400756698247`.
- 50/1: `78362163416325475516394492893621207568547247854013`.

## Useful dependencies, graph occupancy and total costs

Every algebraic mask is checked against original rows; every extracted square
independently reconstructs its exponents/correction, verifies x² = y² mod n,
and checks both GCD signs. The table counts recorded completed extraction
trials, stopped at the first proper factor. A proper-result trial counts once,
not once for each complementary divisor. These are repeated-run observations,
not an unbiased success probability over all possible kernel vectors.

| Band | SLP proper-result / verified trials | DLP proper-result / verified trials |
| --- | ---: | ---: |
| 30 | 36 / 102 (35.3%) | 36 / 192 (18.8%) |
| 40 | 18 / 30 (60.0%) | 18 / 66 (27.3%) |
| 50 | 0 / 0; no dependency | 0 / 0; no dependency |

Proper-result dependencies per total timed CPU-second are 0.08896 for SLP
and 0.13030 for DLP at 40 digits, including the entire pipeline and output
validation. The corresponding 30-digit observations are 1.96682/2.38908
(with the stability qualification above); both 50-digit values are zero.

Thus the factoring win does not require a higher success rate per dependency.
Thirty-digit input 1 alone produced 138 recorded trivial DLP trials across
nine runs, yet completed sooner. These verified trivial congruences are not
correctness failures or evidence for adding unproved character constraints.
`last_solver.dependencies` records only the last completed matrix solve;
unextracted masks are not represented as checked square congruences.

Median terminal matrix statistics show a substantial workload change:

| Input | SLP filtered rows × columns / nonzeros | DLP filtered rows × columns / nonzeros |
| --- | --- | --- |
| 40/0 | 504 × 501 / 9,662 | 321 × 314 / 10,674 |
| 40/1 | 482 × 478 / 9,064 | 313 × 314 / 10,266 |
| 50/0 | 0 × 0 / 0 | 0 × 0 / 0 |
| 50/1 | 0 × 0 / 0 | 0 × 0 / 0 |

DLP's smaller 40-digit matrices are denser. Refresh A8's workload evidence
before a later solver promotion; this study promotes no solver. At 40 digits,
DLP collected 264–287 cycles per run, including 213–242 involving a DLP atom.
Their aggregate mean length was 5.81 atoms. It retained 5,801–6,163 forest
edges, 5,183–5,595 unowned edges, and 249–291 connected components, with zero
evictions. All components participate in cycle detection. Largest reported
owned reservations were 119,903,744 bytes for DLP and 120,467,080 for SLP.

At 50 digits, DLP made 9,795–12,087 splitting attempts per run, formed
127–200 cycles (40–78 involving a DLP atom), and averaged 2.49 atoms per
cycle. All resulting final matrices were singleton-filtered to zero. Median
DLP position counts were 151.1/140.2 million versus SLP's 162.8/152.4 million;
median evictions were 21,504/28,672 versus 12,553/18,209. The 8,192-unowned-edge
policy caused eviction while the overall owned allowance still had room.
Largest reported owned reservations were 163,443,184 bytes for DLP and
129,062,256 for SLP, under the same 256 MiB cap. The larger graph competes
with SLP matches for retention and reduces positions processed within the
cap; raw extra cycles do not establish a factoring win.

No fresh run reported a failed split, split-quota rejection or overlong cycle.
Those finite failure paths remain explicitly tested. The process/JIT RSS
high-water mark reached 265.33 MiB across all prior arms in one process; it
is not a per-arm allocation estimate or a hard 256 MiB OS limit. Owned
reservations and RSS are distinct measurements. Every complete-run timer
includes splitting, certification, graph construction/rebuilding, provenance,
filtering, matrix work, extraction, final classification and output validation.
The result supports the jointly selected DLP/base bundle, not an isolated
causal percentage for the graph algorithm.

## Cold runs and decision

One separate cold process per band/arm used seed 7 and input 0. SLP/DLP
process times were 1.882/1.841 s at 30 digits, 12.679/9.740 s at 40 digits,
and 120.556/120.560 s at 50 digits. Both 50-digit processes remained unresolved.
These include interpreter, corpus/control loading and harness startup. One
process is a diagnostic, not a confidence-backed application cold-start win.
Instrumented residual censuses and the profile remain outside accepted timing.

**Adopt bounded opt-in DLP and the selected 40-digit balanced bundle.** Its
31.73% complete-run benefit, positive interval, unchanged completion, stable
pooled confirmation and independent arithmetic/provenance checks pass the
revised promotion policy for this narrow scope. Keep ordinary SLP as the
default and add no automatic digit cutoff. The 30-digit regression class has
no completion loss and a positive timing observation, with the seed-drift
qualification above. The selected 50-digit bundle earns no performance
promotion. The independently completing 50-digit resume witness establishes
useful larger-input capability under its stated allowances, not a matched win.

This finishes the bounded C1 implementation and calibrated-SLP evaluation.
P5.4's broader SSS comparison remains E1's explicitly recorded prerequisite;
A7 has no calibrated 40-digit SSS arm. Combined portfolio performance, new
50/60-digit calibration and 80–100-digit claims remain unpassed. Generic
sieve/root work is the profile's next actionable investigation; a particular
CRT optimization still needs its own eligibility and complete-run evidence.
Further occupancy/base tuning would require a new bounded training/fresh
study. No additional campaign, runtime default change, merge or push is made.

## Reproduction and retained evidence

Required controls, the owned SLP snapshot, and certified corpora live in
versioned `inputs/`. The original protocols and failed literal gate remain
unchanged. Local `results/c1/` retains `followup/`, `cost-attribution.json`,
`training-1/`, `dlp-profile.*`, `resume-acceptance/` and `confirmation-1/`.
Confirmation's initial captures occupy 711,049 bytes before the separate
quality-analysis file. Its manifest, `summary.json`, `quality-analysis.json`
and every complete/unresolved sample remain available locally.

From this branch, with the configured PyPy 3.11 environment, use a new ignored
output directory for any replay (existing captures refuse overwrites):

```sh
v2/.venv/bin/python -B -u -m v2.benchmarks.c1_implementation confirm \
  --output v2/benchmarks/results/c1/replay-confirmation
```

This replays the already inspected corpus; it is not a new held-out study.
The runner checks the frozen source hashes and acquires the shared performance
lock. Training and the separately bounded resume runner are documented by
`c1_implementation_protocol.md` and `inputs/controls/c1_resume_acceptance.json`.
