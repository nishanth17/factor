# C10: source-derived collection strategy, tranche 1

This is a continuing workstream. The first tranche evaluates bounded ordering
by actual A-product quality and the existing capped prime-power plan on a
settled C9 collector. It does not promote a new runtime default. Normalized
polynomials and per-job threshold adaptation are explicit investigations with
separate implementation gates below. Results and acceptance are recorded
separately from this premeasurement rationale.

## Control and prior evidence

The C10 branch starts at C9 commit
`129cb09e4ca684f3eb08cbbdbea3fb43848e5d95` (runtime frozen at `fbfd662`).
The primary checkout, C3 and C9 worktrees are read-only inputs. C9 owns residual
admission, splitting, matching, cycle construction and retention. Its final
decision and explicit CPU handoff are prerequisites for integrated timing.
The shared `/private/tmp/factor-performance.lock` also covers heavy QA.

The following evidence is not interchangeable:

| Record | Established | Limitation retained |
| --- | --- | --- |
| A1/R2 | Fixed scores and capped plans accepted as opt-ins; isolated 30-digit fixed+plan held-out benefit about 5.7%. | Defaults retained; previous tiny-prime skipping, batched recovery and charging losses remain losses. This is not evidence for every scoring/base/size combination. |
| B1 | Fresh balanced 30-digit SIQS bundle reduced time 36.6%. | Uneven/structured regressions prevent general promotion. The 40-digit MPQS 6.9% result failed its then-frozen gate. Larger censored probes remain censored. |
| C1 | Explicit half-base DLP 40-digit bundle reduced complete cost 31.73%, interval 29.21–33.99%. | 30-digit drift; fresh 50-digit 0/18 at 120 seconds per arm; retention evictions. Later changed-allowance/configuration completions are not isolated acceptance. |
| A7 | Relation, worker, checkpoint and loss-policy reconciliation accepted, 503 tests at that revision. | No collector/worker performance promotion. |
| C9.1 | Exact implementation and source/census work. | First confirmation stopped when retention direction changed; its revealed inputs are training thereafter. |
| C9.2 | User-directed growing storage implementation; selected aligned full-base 50-digit control. | 30/40-digit extended measurements remained unstable. The original 360-second 50-digit confirmation was stopped and superseded by a separately frozen 1,200-second study. No acceptance inferred from partial results. |
| C3 | Source-derived ECM presets promoted at explicit user direction after screening. | 324 calls/54 groups did not pass the fresh interval-width gate. Broader portfolio and schema-12 growth integration remain separate. |

References: [QS audit](../../qs_parameter_audit.md),
[parameter audit](../../parameter_audit.md), [A7](../a7/a7_r5_reconciliation.md),
[C1](../c1/c1_implementation_results.md),
[C9 growth](../c9/c9_growth_results.md). C3's `default_promotion.md` and
`round2_c9_interface.md` were read in its separate worktree. The primary
checkout's untracked `alpertron_gap_analysis.md` was inspected but not copied
or modified. AP01–AP40 are discovery inputs, not an exhaustive catalogue.

## Source provenance and executable paths

[Source manifest](../../inputs/controls/c10_source_pins.json) records exact
revisions, SHA-256, Git blob IDs and licenses for the inspected files. Git
blobs were checked against pinned repository trees; official msieve 1.53 was
checked against its release archive hash. Source snapshots stay in ignored
`results/c10/research/`. This implementation independently expresses integer
algorithms; no upstream source text is incorporated into runtime code.

| Source | Pin / license | Active path and branches |
| --- | --- | --- |
| [YAFU](https://github.com/bbuhrow/yafu/tree/8110dfbd8c6f9486d93b1a02de6eb7b180e55a80/factor/qs) | `8110dfbd8c6f9486d93b1a02de6eb7b180e55a80`; inspected QS files public domain; dependencies separate | `SIQS.c` setup and `siqs_sync`; `new_poly_a`, `computeBl`, `nextB`; `poly_roots_32k`, `med_sieve_32k`, `large_sieve`, `tdiv_scan`, `tdiv_small`, `tdiv_resieve`, `filter`. Native dispatch has architecture/build branches. DLP starts around 70 decimal digits in this source; these native thresholds are not v2 prescriptions. |
| [Official msieve](https://sourceforge.net/projects/msieve/files/msieve/Msieve%20v1.53/) | release 1.53, archive SHA-256 `c5fcbaaff266a43aa8bca55239d5b087d3e3f138d1a95d75b776c04ce4d93bb4`; public domain inspected code | `mpqs/mpqs.c` → `poly.c` → `sieve.c` / `sieve_core.c` → `relation.c` / `sqrt.c`. Its directory name does not mean every polynomial is independent MPQS: it has derived Gray polynomials. DLP activates at 282 bits. |
| [Yamaquasi](https://github.com/remyoudompheng/yamaquasi/tree/3f95f43682ed15d8c1ed206a9a702dd655d7c8ad/src) | `3f95f43682ed15d8c1ed206a9a702dd655d7c8ad`; BSD-3-Clause | `siqs` → `select_siqs_factors` / `select_a` → `prepare_a` / polynomial `next` → `sieve.rs` → `relations.rs`. Type 2 when discriminant is 1 mod 4; DLP default above 256 bits, with different base-count policy. |
| [FLINT](https://github.com/flintlib/flint/tree/a4c9750d0d3d67bb01cf6d18c187591b313451c3/src/qsieve) | `a4c9750d0d3d67bb01cf6d18c187591b313451c3`; LGPL-3.0-or-later | `qsieve_factor_with_tune` → multiplier / prime initialization → `qsieve_init_A` / `qsieve_next_A` / polynomial routines → collection → relation processing → linear algebra. Failed processing continues collection; failed overall round can enlarge the base 10%. This is a control loop, not a guarantee of success. |
| [JavaMath](https://github.com/TilmanNeumann/java-math-library/tree/088d01fa97e7d0412c6cbbb8f463fa7e5f78ab97/src/main/java/de/tilman_neumann/jml/factor/siqs) | `088d01fa97e7d0412c6cbbb8f463fa7e5f78ab97`; GPL-3.0-or-later | `SIQS` constructs `AParamGenerator02` and `SieveParamsFactory02`; generator and supplied sieve/trial-divider implementations form the executable strategy. `d=2` only for kN=1 mod 8. `Sieve03g`, `SingleBlockSieve`, `TDiv_QS_2LP` are selectable implementations, not simultaneously executed stages. |
| [Alpertron](https://github.com/alpertron/calculators/blob/93c5c8189cb2f149c8dae898eb11997c5f2e7980/siqs.c) | `93c5c8189cb2f149c8dae898eb11997c5f2e7980`; GPL-3.0-or-later | `FactoringSIQS` → `sieveThread` → `PerformSiqsSieveStage` → `SieveLocationHit` / `PerformTrialDivision` → full/partial relation routines → `LinearAlgebraPhase`. Packed paired lanes, small-prime periodic presieve and compact partial replay are active. `OldApplets/Siqs.java` is historical; disabled SQUFOF is not active DLP. |

Primary papers checked: Alford–Pomerance,
[Implementing the Self Initializing Quadratic Sieve](https://math.dartmouth.edu/~carlp/implementing.pdf),
§3–4 for polynomial scale/root identity; Carrier–Wagstaff,
[Implementing the Hypercube Quadratic Sieve](https://homes.cerias.purdue.edu/~ssw/qs4.pdf),
§2–3.1 for sieving/recovery and A construction. PDF hashes are respectively
`57c7f6e10a53b0cd9c61c26f7801284aea651b38d01848630a99c0b7399829ac`
and `ff7b5f1104a315cbec810456804c51dd7880ba1792803e4c1b5eb9655653dda7`.
They establish mechanisms, not Python timing expectations. The AMS Silverman
PDF could not be fetched (HTTP 403), so no claim is based on reading it.

Author notes were used as discovery material and checked against code:
[Yamaquasi QS notes](https://github.com/remyoudompheng/yamaquasi/blob/3f95f43682ed15d8c1ed206a9a702dd655d7c8ad/README_qsieve.md),
msieve `Readme.qs`, and
[YAFU documentation](https://sites.google.com/site/bbuhrow/yafudocumentation).
The latter's historical SLP/DLP description does not describe the pinned
current caller. Yamaquasi's displayed normalized-polynomial denominator omits
A in one explanatory equation; the executable polynomial and exact identity
below require `4A`. Informal timing/blog claims are not accepted measurements.

## End-to-end transfer map

| Stage | Executable source observation | v2 counterpart and transfer boundary |
| --- | --- | --- |
| Multiplier search / exceptions | msieve `choose_multiplier` tries a fixed set including even k; FLINT `qsieve_knuth_schroeppel` handles even n first; Alpertron `FactoringSIQS` tries 1, 2 and odd primes through 97; Yamaquasi `select_multiplier` searches an input-size-bounded integer range. | `select_multiplier` restricts scored candidates to distinct squarefree odd k<=255 and scores primes through 97. `build_factor_base` accepts bounded positive k, checks proper gcd first, and rejects gcd=n. Public QS rejects even n; the factoring preprocessor owns its factor 2. Extending scored candidates needs a defined even-D branch, not indexing the current odd mod-8 table. |
| Coupling multiplier to polynomial kind | Yamaquasi `expected_smooth_magnitude` uses a mod-8 weight of **3** for D=1, explicitly including its type-2 scale/valuation benefit; msieve's classic weight is 2. FLINT and Alpertron leave their initial half-log-two weight on even-D branches; msieve assigns no such contribution. | **C10-N08: missing representation-coupled scoring.** Current `select_multiplier` uses classic weights 2,1,1/2,1/2 for odd D. Do not copy Yamaquasi's 3 before implementing normalized polynomials, or normalize without considering its altered multiplier ranking. Defer together with N01; exact candidate/root/gcd and size-stratified complete-output checks apply. |
| Base size versus prime endpoint | Yamaquasi `FBase::new` targets an eligible-prime count rounded for SIMD; its `fb_size` adjusts the effective input for D=1 mod 8 and DLP. msieve/FLINT initialize count-oriented tables; Alpertron bounds a count formula with `MAX_PRIMES`. | `build_factor_base` includes eligible primes strictly below a prime-value bound. Native count/endpoint values cannot be compared numerically without constructing the actual base. Record both count and largest prime, singular multiplier primes and owned memory. |
| Powers and log units | Yamaquasi `smooths` reconstructs skipped small-prime contributions and adjusts near roots; JavaMath Factory02 scales logs around a byte threshold and has distinct sieve-hit, trial-division and smoothness bounds. msieve's main marking subtracts one log per prime hit, with allowances and later division. | `_weights`, `_power_marks`, `prime_power_roots`, `score.log_bounds` explicitly cover valuations, including bounded singular-lift fallback. `fixed` uses integer log intervals; bytearray saturation clips thresholds too. Native heuristic omitted powers are not equivalent conservative coverage. Powers' repeated-lift cost motivates the first policy arm. |
| Interval and distribution | Alpertron's half interval comes from a capped size formula; Yamaquasi uses 32,768-position blocks and `smooths` accounts for distance to a polynomial root; JavaMath derives a maximal Q/(dA) estimate from its sieve-array convention. | `half_width` covers `[-M,M+1)` (2M+1 integer positions), distinct from `block_width` and pipeline `batch_width`. `_bounds` checks endpoints/vertex, and `_candidate_passes` refines exact values. Width grows coefficient values and changes A target, cache reuse and useful yield jointly. |
| A pools and scheduling | msieve `poly_init` avoids very small A primes; `build_base_poly` corrects the last factor; FLINT uses separate core/correction index classes. Alpertron picks a target-root-centered index span and uses `getFactorsOfA` without complete-product ranking. | `family_assignments` is the bounded reference; `AssignmentStream` enumerates unique combinatorial cores with a deterministic permutation and disjoint corrective flyer. `PolynomialFamily` constructs CRT terms, fixes one sign, and uses Gray transitions. `diverse` sampling can spread Gray indices but does not constrain inter-A overlap. |
| Root initialization / updates | YAFU `computeBl` and polynomial root arrays, Yamaquasi `prepare_a` and `next`, Alpertron `sieveThread` prepare per-family update data. | `PolynomialFamily.__init__` already caches A inverses; `_make_step` applies centered-B deltas and exact branches for 2 and p|A. `_lift_inverses` caches derivative inverses within one polynomial. New large tables must justify their added construction/replay storage. |
| Presieve / density classes | Alpertron's `PerformSiqsSieveStage` repeats a periodic tiny-prime pattern and uses packed paired lanes. YAFU separates medium and bucketed large primes; msieve reuses block root offsets; Yamaquasi separates primes below a block, below an interval, and sparse larger hits. | `_sieve` / `_add` provide dense/sparse/bucket choices; general Python marking is not those specialized native kernels. A1a receives periodic presieve and precomputed transition-table work; SIMD transmutation/alignment details are native-only. |
| Candidate scan and exact recovery | YAFU scan → small-prime division → resieve/full residual work; msieve `scan_sieve_block` uses packed scan masks; JavaMath selectable sieves use a second filter before `TDiv_QS_2LP`; Alpertron hit handling invokes exact trial division and restores A factors. | `_candidate_passes`, `_resieve`, `_divide` use exact values, repeated division, signs, supported A and checked residuals. Different scan/recovery paths are available; choose by measured candidate density. Skipping a prime in marking never permits skipping it in reconstructed exponents. |
| Admission through useful matrix | Upstream SLP/DLP thresholds, opaque residual treatment and graph retention differ substantially. YAFU `siqs_sync` counts relations/cycles toward a base+extra target, then can raise the target on failed recovery. | C9's `divide_large_residual`/`admit_large_relation` certify bounds and endpoints; matching, cycles and retention stay fixed in C10. `QSJob._solve` prepares original provenance, filters after sufficient row growth, and checks usable row/column excess. Raw edges and cycles cannot replace this outcome. |
| Solve / extract / retry | FLINT invokes relation processing at its count target and can grow the base after failure; msieve and Alpertron use native linear algebra and square-root recovery after collection. | `filter_matrix`, `DependencySolver`, `DependencyExtractor` recheck original relations, parity and proper gcd; trivial dependencies can lead to more collection. Base-growth remapping is absent and would require a separate representation/resume design. C10 may investigate cadence; solver/storage ownership stays elsewhere. |

## Ranked gap ledger

Rank is a tranche decision, not a claimed percentage speedup. C9's instrumented
50-second diagnostic attributed 32.87 cumulative seconds to sieving and
23.56 to lifting in the aligned growth arm, versus 31.69/22.98 in fixed DLP.
Those nested times must not be added; profiler overhead prevents treating
23.56/50 as a reliable full-run speedup ceiling. Residual division, classification
and splitting were much smaller in that observation (1.33, 0.57, 0.53 seconds).
This favors testing existing root-plan reuse before building another residual
optimization. No v2 stage measurement yet establishes a product-ordering or
normalized-polynomial speedup.

| Rank / ID | Source evidence; exact v2 counterpart | Classification, benefit and disposition |
| --- | --- | --- |
| 1 / C10-N04 | Yamaquasi `prepare_a`/`next`, YAFU `computeBl`/root update, msieve `build_derived_poly`; v2 `SieveCollector._power_marks`, `_set_plan_interval`, `SIQSJob._activate` | **Uncalibrated policy + partial equivalence.** v2 already caches bounded per-polynomial power plans across blocks; no new kernel is required to try `power_plan_bytes`. Evaluate an isolated 2 MiB plan arm on C9's current `powers` scoring. Charge construction, fallback and resumed reconstruction. It preserves score coverage; plan capacity can waste memory or construction time. Family A inverses and live Gray updates already exist; larger precomputed transition tables and cross-polynomial lifting reuse remain A1a. |
| 2 / C10-N02 | Yamaquasi `select_a` (small-k product enumeration/ranking, larger random/corrective selection); FLINT `qsieve_init_A`/`qsieve_next_A`; JavaMath `AParamGenerator02.computeNextAParameter`; v2 `AssignmentStream._raw_item`, `family_assignments`, `PolynomialFamily` | **Missing complete-product priority; partial equivalent flyer.** Existing flyer already minimizes actual product error for a given core. Implement opt-in ordering of fixed 16-family batches by exact `abs(A-target)`, preserving the population and uniqueness. Measure setup overhead, A error, overlap and useful yield. No product-quality filter or diversity guarantee is claimed. Low implementation cost; C9 representation unchanged. |
| 3 / C10-N01 | Yamaquasi `PolyType`, `select_siqs_factors`, `prepare_a`, `_finish_polynomial`; JavaMath `SIQS` d choice and `SIQSPolyGenerator.computeBParameter`; v2 `Polynomial`, `polynomial_roots`, `a_target`, `PolynomialFamily`, `prime_power_roots`, `AtomicRelation`, checkpoints | **Missing mechanism.** Normalize discriminant polynomials to reduce value scale and exploit the power of two when permitted. High-confidence arithmetic; unknown useful-yield/total-cost gain. Independent exhaustive oracle now; gated representation tranche below. High implementation/replay cost, overlaps root kernel ownership and must coordinate A1a, A7 and C9. |
| 4 / C10-N03 | YAFU `siqs_sync` and trial-division cutoff; JavaMath `SieveParamsFactory02`; v2 `SieveCollector._bounds`, `_candidate_passes`, `SieveConfig.threshold_extra` | **Missing feedback policy; existing local value thresholds.** Adapt per job to measured candidate cost/useful yield, not just input digits. v2 already derives block/value thresholds exactly. Positive extra is intentionally lossy. Deferred controller pending rejected-position/useful-dependency audit and charged checkpoint state; see below. |
| 5 / AP10–14, N05 | Alpertron `FactoringSIQS`, `getFactorsOfA`; msieve `choose_multiplier`/`poly_init`; Yamaquasi `select_siqs_factors`; v2 `select_multiplier`, `build_factor_base`, `a_target`, config width/base/A count | **Equivalent mechanisms, uncalibrated joint policy.** Broader multiplier/base/window/factor-count search can improve smoothness but changes matrix size and per-family amortization. Preserve B1/C1 winners and losses. Even multiplier search is absent in the scoring API although native sources use even candidates. Needs explicit mod-8 accounting and proper-gcd tests; do not append native constants blindly. |
| 6 / A1a handoff, AP05–09 | Alpertron `PerformSiqsSieveStage` / root setup, YAFU dense/bucket SIMD dispatch, msieve `fill_sieve_block`; v2 `_sieve`, `_power_marks`, `_add`, family root transitions | **Missing kernels / native-only layouts.** Bounded periodic presieve, precomputed root deltas/inverses and density-aware marking are coherent kernel candidates. SIMD intrinsics and paired byte tricks are not direct PyPy transfers. Detailed handoff below; do not implement a competing kernel in C10. |
| 7 / N06 | msieve `fill_sieve_block` / `scan_sieve_block`, JavaMath skipped-prime compensation, Alpertron `PerformTrialDivision`; v2 `_sieve`, `_resieve`, `_divide`, `small_prime_cutoff`, score backends | **Equivalent options, previous losses and uncalibrated policy.** Scanning, root-guided exact recovery and resieving exist. Tiny-prime omission needs conservative compensation or a declared loss experiment; old losses remain. Saturation must clip thresholds consistently. A smaller candidate list alone is not a win. |
| 8 / N07 | FLINT relation-processing threshold and base-growth loop; msieve base+96 target; Alpertron max(33*base/32,base+50); v2 `QSJob._solve`, `row_excess`, `filter_row_growth`, `weight_two` | **Uncalibrated collection policy, separate solver ownership.** Cadence can avoid repeated filtering/solving or overcollection. Measure useful post-filter rows, original-row dependencies and extraction, not raw row counts. C10 may tune cadence later; matrix representation/solver changes belong to A8/B6/B7/C4/D1/D2. |
| separate / C9 deltas | `large_prime_collector`, `large_primes`, storage/checkpoint; native residual matching/graphs | **Owned by C9.** Keep exact certified endpoint/product domains and retention policies settled for C10. Native opaque residual matching is not equivalent to v2 certainty. No new split/match/graph changes here; C9a/C1a remain separate. |
| separate / C8 | Gray schedules versus CRT half-sum construction; v2 `PolynomialFamily` and `sss` schedulers | **Existing Gray mechanism; separate missing scheduler.** Product ordering neither supplies family-wide CRT half-sums nor establishes SSS superiority. C8/E1 remain separate. |

## Normalized / mod-8 polynomials: exact transfer contract

Let `D=k*n`, `A` be positive odd, `D=1 (mod 4)`, and choose odd `b` with
`b²=D (mod 4A)`. Then

```
P(x) = A*x² + b*x + (b²-D)/(4*A)
U(x) = 2*A*x + b
U(x)² - D = 4*A*P(x).
```

For D=1 mod 8, P is even for every integer x; for D=5 mod 8 it is odd.
The normalized derivative `2*A*x+b` is odd, so lifting roots at 2 is
nonsingular, unlike the always-even derivative of the standard F. This is a
second plausible kernel-cost benefit beyond smaller values, still unmeasured.
For a nonsingular odd factor-base prime p the roots are
`(±sqrt(D)-b)/(2A) mod p`. For p dividing A, solve the remaining linear or
constant polynomial; for p=2 enumerate residues. Singular prime-power lifting
must remain conservative. Multiplier gcds give a factor only if proper; even
input is handled before odd-discriminant construction. Even D and D=3 mod 4
cannot be forced into this representation by simply making b odd.

For the existing `F(y)=A*y²+2*B*y+(B²-D)/A`, choose delta in {0,1} such that
`b=B+delta*A` is odd. Substitution gives `F(2*x+delta)=4*P(x)`.
For example, D=1009, A=15 and B=2 give `F(y)=15y²+4y−67`,
`P(x)=15x²+17x−12`, and `F(2x+1)=4P(x)`. At x=0 the recovered identity is
`17²−1009 = 4*15*(−12)`: sign, A and the factor 4 all participate.
This is a parity subsequence with a different covered y interval. It does not
prove fourfold useful yield. Yamaquasi consequently targets approximately
`sqrt(D/2)/M` instead of `sqrt(2D)/M`, with M the half interval: at adjusted A
the typical size improvement is a factor of about two, not a timing claim.

JavaMath chooses d=2 only for D=1 mod 8, sets b odd, and evaluates
`Q=(d*A*x+b)²-D`; its trial divider restores fixed factors from `d*A`.
For d=2, `Q/(d*A)=2*P`. This differs from Yamaquasi's stored P, even though
the same square identity underlies both. Its comment says D=5 mod 8 has no
notable gain; that is an author observation, not evidence to exclude that
case from a v2 experiment.

JavaMath's A generator explicitly says its adjusted A preserves the maximal
size of Q/(dA), while Factory02 divides its logarithmic size estimate by d.
Those source-level heuristic choices must be audited with actual polynomial
values and recovered powers of two; the threshold estimate is not an algebraic
proof of smaller Q/(dA). This is another reason not to transplant its constants
independently of its sieve and recovery conventions.

`oracles.py` uses brute residue enumeration and exact integer division,
independent of v2's root/CRT implementations. Tests cover mod-8 parity,
both signs, p dividing A, p dividing D, p=2, square-modulus residues and full
signed reconstruction. A production tranche must additionally cover:

1. A tagged polynomial kind in identities, relation verification, extraction,
   checkpoints, legacy replay and worker assignments. Recover the factor 4
   exactly; it is a square but cannot disappear from reconstruction.
2. Exact target/rounding and interval bounds; Gray transitions preserving odd b
   and `b²=D mod 4A`; independent whole-family root/valuation checks.
3. Conservative score bounds and saturated lanes including unavoidable powers
   of two; candidate identity checks against exhaustive trial division.
4. All exceptional multiplier/A primes, zero values, negative values,
   interrupted setup and charged simultaneous checkpoint reconstruction.
5. Fresh complete-run comparison stratified by D mod 8, after source-informed
   size/width calibration. Do not mix this representation change into the
   first quality-order experiment.

## Product quality, spacing and diversity

These are distinct policies. Yamaquasi's small-factor-count path sorts actual
products by target distance, deduplicates, truncates and finally sorts A
numerically. Its larger path uses random products plus a corrective final
factor. FLINT splits core/correction index classes and accepts a broad target
window; msieve and JavaMath choose a corrective last factor. None establishes
that v2 should globally retain only the closest A values.

Alpertron `getFactorsOfA` excludes duplicate indices and one directional
adjacency condition within an A; it does not check distance to the previous
family. JavaMath's shared-q counting block is debug-only, so it is not an
active minimum-distance filter. Duplicate A prevention is not low overlap.
v2 `diverse=True` controls within-family scheduling; it is not a global
constraint on A-factor overlap.

The new `quality` stream permutes each fixed batch of 16 existing flyer sets.
Batch boundaries use the full combinatorial space, not the current grant,
so extensions preserve prefixes. There is no growing seen set. Actual products
and sort keys are charged before publishing a complete bounded cache; replay
reconstructs and charges it. The workspace reservation includes old/new
batches, sort keys and scratch. The metadata reservation also covers the live
cache while a checkpoint reconstructs assignments. Finite and growing C9
storage retain their existing None/count/memory semantics.

The independent census reports exact A error, unique products and both all-pair
and adjacent factor-overlap histograms. Full-batch population overlap is
unchanged; adjacency and incomplete-prefix overlap can change. If quality
improves only raw relations or worsens useful dependencies, defer it. An
explicit bounded diversity selector is a later arm only after this census
establishes a reproducible problem; arbitrary spacing is not assumed helpful.

## Per-job threshold adaptation

YAFU `siqs_sync` is active at threadpool synchronization (and the direct loop).
When optimization is enabled it compares relation rates at the initial cutoff,
initial+5 and initial−5, using 10 averaging units, or 5 above 300 bits. The
following `else if >320` branch is unreachable in this pin, so its apparent
2-unit setting is not active. A 1.02 SLP / 1.05 DLP bias favors the original
cutoff. Explicit small-cutoff override and restart-with-relations disable this
optimization. The feedback objective is relation rate, not useful dependencies
or complete-factor time; transplanting that metric would be unjustified.

v2's `_bounds` already responds to polynomial value distribution per block;
`_candidate_passes` refines with each exact value. `threshold_extra` is static
and intentionally lossy, while `small_prime_cutoff` can omit marking with
conservative compensation. Native skipped-prime expected-log adjustments are
not automatically conservative guarantees.

Proposed later controller: a finite pilot over a fixed set of extras including
zero, bounded epochs with deterministic tie-breaking, and a fixed held-out
portion of rejected positions recovered exactly. Record candidate/division
cost, full rows, admitted SLP/DLP, post-filter useful rows and kernel success.
Audit lost admissible atoms, matching opportunities and dependencies; a
residual-magnitude census alone is insufficient. `threshold_losses` supplies
only that smaller independent arithmetic diagnostic and says so explicitly.
Freeze the loss policy, pilot overhead and training/confirmation split before
testing. Persist/charge pilot position, chosen threshold, counters and
reconstruction; do not copy YAFU's restart disablement as silent state loss.
C9 must approve any change to its candidate/admission boundary. Defer this
controller until candidate recovery, rather than lifting, is measured as the
limiting stage or the audit demonstrates substantial useful surplus.

## Concrete A1a handoff

Start with the settled C9 50-digit full-base configuration and report exact
source/configuration identities. Consume C10's isolated capped-plan result
first: it may already remove repeated block lifting without new code.

Candidate kernel tranche: precompute a bounded table of CRT-term root shifts
(including the centered-B correction), beyond the existing cached A inverses
and per-step multiplication. Investigate reusing nonsingular Hensel derivative
inverses: at a root the derivative is ±2*sqrt(D) modulo p, independent of B.
The current collector caches these only within a polynomial and clears them
on switch. Then consider a bounded periodic presieve for tiny primes and
density-separated marking only if their affected-stage costs remain large.
Source anchors are Yamaquasi `prepare_a`/`next`, YAFU `computeBl`/root-update
arrays and Alpertron `PerformSiqsSieveStage`. Translate mathematical structure,
not licensed implementation text or native timing tables.

Required output: independently enumerated roots and valuations for exceptional
2/A/multiplier primes; unchanged candidate identities for conservative arms;
explicit loss audit otherwise; simultaneous live/temporary/replay storage;
charged fallback, interrupted setup and checkpoint rebuild. Preserve fixed
and growing collector controls. Profile setup, marking, root lifting, scanning
and exact recovery separately, then require matched complete outputs and
complete-factor timing. A1a owns the kernel implementation. C10 owns subsequent
size-band policy calibration. C8 owns CRT half-sums; A8 and solver owners own
matrix changes. No existing acceptance record is reopened.
