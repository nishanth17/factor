# Changelog

Start with [current status](v2/STATUS.md) for present defaults and work. Entries
below record decisions in their dated scope; later entries may supersede routes.

## B4a Phase 1 overhead continuation (9 November 2026)

- Amortize process-CPU clock reads while retaining per-action wall and
  cancellation checks plus root resource polls around foreign slices.
- Retain one fixed pthread team across drained generations; join on close
  and on partial startup failure. Add deterministic clock and lifecycle tests.
- Preserve the noisy same-cell re-baseline, matched lifetime comparison and
  failed harness attempts locally. No complete-ECM speedup is established.
  Stop before the full-ECM port on its canonical recovery/checkpoint design;
  setup, classical stage two and strict recovery remain Python.

## Experimental B4a pthread correctness (9 November 2026)

- Synchronize committed A7a into the B4a branch. Add C-only pthread dispatch
  with checked dynamic/static chunks, private GMP contexts, joined cleanup,
  deterministic faults/cancellation and saved-history resume coverage.
- Admit scalar/pthread correctness from committed source `0f2fee3`: native
  builds, 1,712 tests (20 optional skips), lint and 308 ordinary runner imports.
  Retain failed diagnostic receipts; no performance gain, default change,
  mainline integration or push is adopted.

## B4a native correctness preparation (8 November 2026)

- Correct the owned-group observer's Mach CPU-unit conversion without adding
  native-thread CPU twice. Add deterministic counter regressions.
- Add fixed 60/80/100-digit curve, nonunit recovery and native H3 lifecycle
  admission checks, and build header/linkage provenance. Runtime validation
  remains separate from source preparation.

## B4a H3 integration (8 November 2026)

- B4a branch: integrate H3 root resource accounting into optional native ECM
  ownership and resume; preserve kernels, historical protocols and SIMD loss.
  Native startup, engagement and worker economics require revalidation.


## A7a selected QS integration (9 November 2026)

- Adopt Y0 numerical policy and nearest-row process QS: serial below 40;
  C4/p4, C4/p8, C4/p16, C8/p16, C14/p64 and C8/p64 at 40/50/55/60/65/70.
  Ties choose the smaller size; above 70 is explicitly extrapolated.
- Add CLI `--qs-workers serial` and retain API `qs_workers=None` as S0.
  Preserve H3 resources, explicit and saved recipes, lazy pool startup, ECM
  defaults and existing portfolio SIQS eligibility. Restore all supported
  collector subtypes through the shared checkpoint decoder.
- Close the selected A7a scope with its G4/calibration evidence and retained
  censoring/drift limits. Defer 75/80 probes until after B5a and flag C3a/B12
  dispatch and cost tables for revalidation.

## Directional QS count extension (8 November 2026)

- Permit explicit current QS process counts through twenty for the bounded
  directional study. Automatic defaults and default search grids stay fixed.
- Extend exact fixed-work coverage, checkpoint restore and child cleanup
  validation to sixteen, eighteen and twenty spawned workers.

## Opt-in QS oversubscription experiments (8 November 2026)

- Permit explicit current QS process counts up to fourteen for the authorized
  worker-count race. Automatic choices and default search grids are unchanged.
- Extend fixed-work coverage, checkpoint restore and child cleanup checks to
  twelve and fourteen spawned workers; larger counts remain rejected.

## H3 observed resource limits (8 November 2026)

- Replace internal byte/count/work refusal ceilings and prepaid worker leases
  with one root wall/aggregate-CPU/current-tree-RSS owner. Omitted wall/CPU are
  unlimited; RSS defaults to floor(80% of physical memory). Preserve numerical
  recipes, actual cache/FIFO/packet capacities and exact validation.
- Use saved-history accounting without a durable ledger; retained live owners
  keep unsaved charges. Stage shares hand remaining root resources to fallback.
- Reject removed arguments by name and old resource checkpoint versions without
  migration. Make atomic checkpoint saves optional and uncapped by default;
  oversized saves record `checkpoint_size_skip` and preserve the result.
- Committed-only `5feff32` passes 1,404 tests (9 optional skips), lint and 267
  benchmark imports. Fixed-work QS/ECM confirmation on `508c490` matches base
  `40d2367`. No performance evidence or retuning is added; dependent costs
  require revalidation.

## Experimental QS integration with H3 (8 November 2026)

- Integrate local master into the isolated A7a branch. Replace current-worker
  reservations and synthetic charges with the observed root resource owner;
  retain numerical policies, packet/cache settings and exact validation.
- Version worker snapshots for saved-history accounting, retain arithmetic and
  corruption/cleanup coverage, and adapt historical recipes at recorder ingress.
  Preserve earlier diagnostic work separately. Handle failed process submission
  without broadcasting through a dead worker's stop lock. Pass 1,628 tests,
  lint and 295 benchmark imports. No new probe or timing claim.

## Experimental QS kernel specialization and restricted tuning (8 November 2026)

- Report both approved 75/80 feasibility probes on the isolated A7a branch.
  Retain incomplete internal memory outcomes, missing sub-limit attribution
  and the terminal path's discarded QS state. No automatic retry or adoption.
- Admit explicit QS bookkeeping grants through 32 GiB without changing the
  default. Probe cells can enforce a lower sampled live-RSS cancellation
  bound and retain their portfolio checkpoint, separately from recorder errors.
- Select default scalar marking, collection and native-stamp kernels at
  construction, including prime-power root lifting. Preserve explicitly
  selected alternatives, exact signatures, work accounting and resume behavior.
- Finish the approved four-row worker/packet calibration on the isolated A7a
  branch. Preserve noisy/inconclusive rankings, complete worker timelines and
  earlier failed or superseded measurements. Recipes remain experimental;
  automatic defaults and numerical policies are unchanged. Intermediate
  verification, fresh confirmation and feasibility probes remain gated.

## Experimental QS worker overhead checks (7 November 2026)

- Add explicit family-affinity placement, grouped exact marking leases and
  byte-sized recovery-word stamps on the isolated A7a branch. Preserve exact
  arithmetic, paid-prefix refusal, resume validation and conservative memory
  reservations. The experimental selectors remain disabled by default.
- Require worker active, queue-wait and publication timelines for G2H/G3
  process experiments; disable them for G4 confirmation. Reject incompatible
  observation settings before starting actors.
- Report the bounded overhead comparisons and full-factor remeasurement.
  Configuration tuning and fresh performance confirmation remain open.
- Close the final source gap with the requested interleaved fixed-work
  comparison. Preserve exact signatures/scanned positions and explicitly
  list default-path changes; this Tier 1 check is not a new speedup claim.

## QS reference allowance correction (7 November 2026)

- Apply the explicitly declared checkpoint allowance in the experimental
  mainline reference adapter, preserving its numerical policy. Retain thermal
  readings outside actor timing. G2 stopped before continuation; no production
  factoring defaults change.

## Bounded QS diagnostic review (7 November 2026)

- Review explicitly selected parent/worker CPU streams for bounded windows,
  retaining the unreviewed process inventory. Full profiles still require
  every owned process; missing targets and unflushed streams fail review.
- Finish the approved G2 breadth report on the isolated branch. Preserve
  incomplete mainline controls and historical recorder failures; architectural
  fixes, configuration tuning and fresh confirmation remain separately gated.

## Experimental QS campaign preparation (7 November 2026)

- Add bounded, explicitly gated stage proposals, checked coverage and durable
  per-stage accounting on the isolated A7a branch.
- Add independent certified confirmation inputs and a literal-mainline
  complete-QS reference adapter with bounded cold/warm ownership.
- Add optional worker service/queue spans and minute-level progress, including
  fixed/adaptive relation targets. Keep unresolved observer costs diagnostic.
- Pass the G1 preparation gate; coverage/search/confirmation and default
  adoption remain pending. No automatic factoring settings change.

## Experimental current-collector QS continuation (6 November 2026)

- On the isolated A7a branch, add bounded worker/cache configuration controls,
  checked resource accounting and a certified per-regime calibration corpus.
- Pair every clean benchmark repetition with parent/all-worker CPU tracing,
  owned RSS and parent/representative-worker heaps, analyzed before selection.
  Preserve recorder failures and unknown CPU charges separately.
- Add bounded joint designs with explicit geometry, conditional validity and
  pair-coverage checks. Record exploratory full 60/70-digit feasibility,
  cache/GC screens and the first joint 60-digit block. Worker
  selection, fresh confirmation and integration remain open; QS automatic
  settings retain their existing behavior.

## ECM package and bundled-data cleanup (6 November 2026)

- Move engine-specific execution into the flat `v2/ecm/` package, retaining
  old imports as compatibility aliases. Keep shared budgets/stages separate.
- Replace bundled JSON program tables and worker settings with lazy Python
  constants. Retain immutable benchmark JSON oracles and exact program checks.
- Replace runtime source scanning with explicit compatibility identifiers.
  Selector v3 and automatic metadata v2 reject older corresponding metadata;
  serial checkpoint behavior, numerical choices and arithmetic remain unchanged.

## B12/F4 closeout and cold serial admission (5 November 2026)

- Close and integrate the delivered process-ECM workstream into mainline
  with its selected numerical/cold/warm regimes. Keep ordinary QS fallback,
  explicit serial/resource settings and checked resume. The benchmark-only
  disabled-fallback configuration does not become a production default.
- Bypass policy loading for cold compact64 at 30 digits and fixed32 at 50/70.
  Keep warm/custom numerical schedules on ordinary selection. Admission-only
  checks reduce about 44 ms to 0.035 ms; 20-digit admission stays unchanged.
- Current committed-only QA passes 1,357 tests (9 optional skips), full lint
  and 265 runner/package imports. Retain ECM-only gains, quota failures and
  remaining cheap/warm regressions; make no new whole-portfolio or p90 claim.
- Preserve research, receipts and source identities in local ignored evidence;
  remove 21 clean F4 worktrees and 41 other temporary snapshots/caches after
  integration. Keep the shared lock file and other workstreams untouched.

## B12/F4 automatic default integration (5 November 2026)

- Include selector bookkeeping in parent-memory admission. Automatic and
  borrowed portfolio contexts decline workers safely when their combined parent
  reservation exceeds the worker ceiling; direct component calls refuse with
  `MemoryError`. Keep arithmetic, measured policy rows and serial fallback RNG.
- Add a source-isolated ECM-only comparison using selected cold/warm defaults
  and incumbent serial arithmetic on identical curve manifests. Preserve its
  seed adaptation explicitly, validate proper divisors independently and
  include declined-policy serial cost. Disable fallback in this endpoint.
- Enable automatic cost-gated process ECM in the normal auto CLI and bounded
  library API, with explicit serial and aggregate-memory controls. Keep exact
  arithmetic, caller limits, legacy checkpoints and stock fallback settings.
- Own automatic pool cleanup on every exit and before final resource snapshots;
  bind automatic policy, implementation, deterministic seed and memory on resume.
- Promote escalating140 at 60/80/100 omitted-tier inputs and package unchanged
  measured nominees/costs with integration source pins. Retained contexts remain
  available separately; cheap/unsupported/unfunded calls stay serial.
- Record the authorized adoption exception: balanced completion and precise
  tails remain unestablished, some cold medians regress and CPU increases.
  No new speedup claim, mainline merge, push, B12/F4 closure or E1.
- Committed-only `094ab1d` passes 1,337 tests (9 optional skips), full lint and
  264 benchmark imports; automatic snapshots use outer schema 13 so older
  readers reject their worker state. Computation and its lock are released.

## B12/F4 targeted policy and hard-case evidence (5 November 2026)

- Add bounded certified calibration/confirmation generation, independent tests,
  immutable corpora and a profile frozen before fresh input generation. Fit only
  an 81-digit cost anchor; holdouts validate within 4–17%.
- Record 1,248 checked cells / 1,838 timed calls against actual master. Fresh
  uniform completion is 22/36 versus 18/36; completion-sensitive mean latency is
  lower cold/warm, while some cold medians regress and CPU rises 4.37×.
- Keep a separate exploratory 415-curve hard pilot, all stock checkpoint refusals,
  censored tails, resource charges and the pending balanced guardrail explicit.
- Committed-only `9b7eaeb` passes 1,311 PyPy tests (9 optional skips), lint and
  263 benchmark imports. Computation is finished and the lock is released.
  No default promotion, merge, push, B12/F4 closure or E1 campaign.

## B12/F4 cofactor-aware selection and incumbent serial path (5 November 2026)

- Add bounded piecewise linear per-curve cost estimates between same-plan size
  anchors, partitioned at active PRAC eligibility boundaries. Retain bundle-bound
  crossovers, fixed engaged pool settings and explicit estimate provenance.
- Keep portfolio declines on incumbent serial arithmetic and RNG, with a sticky
  cofactor checkpoint decision. A declined call no longer pins a serial context;
  an engaged search never changes stream after a later refusal.
- Version selector profiles/checkpoints as v2, bind component/portfolio mode and
  add a source-pinned training catalogue with unchanged historical cost anchors.
  Native wire v7, mathematical SEARCHv1 and automatic defaults remain unchanged.
- Add independent cost, serial/fallback RNG, cancellation, ownership, resume and
  actual-cofactor spawn regressions. Committed-only `4a2f051` passes 1,308 PyPy
  tests (9 optional skips), lint, 262 imports and source-pin admission. No new
  performance result, default promotion, merge or B12/F4 closure.

## B12/F4 actual-master comparison (5 November 2026)

- Preserve benchmark work/CPU accounting through result/checkpoint construction
  failures; add three regression tests and finite censored refusal recording.
  Production arithmetic, resource limits and automatic defaults are unchanged.
- Version fresh independent certified/uniform inputs, the stock-master comparison
  policy and a separately frozen caller140 supplement/resource continuation.
- Record all 33 regimes, cold CLI and warm context outcomes: primary uniform
  completion is 18/36 versus master's 19/36 on the default schedule, or 24/36 with
  explicit caller140. Median uncertainty, the 80-digit losses, sparse actual
  worker engagement and pending balanced guardrail preclude default adoption.
- Final report/input committed checkout `20e6681` passes 1,287 PyPy tests
  (9 optional skips), full lint and 262 benchmark imports; the lock is released. Detailed receipts remain local. No B12/F4
  closure, merge or E1 campaign.

## B12/F4 measured opt-in catalogue (4 November 2026)

- Add lazy, source-checked selection for 33 exact mathematical regimes with
  separately calibrated cold/started/warm thresholds and explicit retained-pool
  intent. Preserve exact ECM arithmetic, SEARCHv1 identity, finite v7 ownership,
  actual incumbent control and automatic serial defaults.
- Version reproducible configuration snapshots, a disjoint certified/uniform
  corpus and censor-aware cluster analysis. Retain full failures and charges.
- Record repeated fixed-work scaling and 448 held-out complete calls, including
  median uncertainty, the 80-digit regression, explicit matched QS checkpoint
  funding and the pending balanced guardrail. No production adoption or closure;
  final committed-checkout tooling QA remains queued with the lock released.

## B4a optional ECM NEON screen closure (4 November 2026)

- Measure the fixed two-curve NEON prototype against GMP and the B4a scalar
  multiplication/reduction path. All six complete hybrid cases cost 2.27–3.14x
  the scalar path; close this prototype without integration or default changes.
- Pass exact targeted checks and lint before timing. Retain raw interleaved
  batches, uncertainty, setup failures and the original finite charge. Scalar
  C/GMP and OpenMP evaluation remains open.

## B4a optional ECM NEON preparation (4 November 2026)

- Prepare a separate two-curve product/square prototype at 200/266/333 bits,
  exact GMP controls and a bounded prospective Stage A runner. Include hybrid
  conversion/reduction costs and distinguish SIMD from layout effects.
- Preserve scalar C/GMP and OpenMP work. No executor integration or performance
  verdict; build/correctness admission remains separate from later timing.

## B12/F4 larger-quota follow-up (4 November 2026)

- Validate the 512-curve checkpoint/resource contract and matched spawn execution
  without arithmetic or default changes. Record 6.2% lower warm time from cached
  indices versus packets alone; cold improvement over the prior F4 incumbent is
  inconclusive. Preserve the recorder failure, late exit proof and finite charges.
- Keep complete-factor and engagement confirmation open; no default promotion.

## B12/F4 checked worker fast paths (4 November 2026)

- Add opt-in prospective work accounting, bounded immutable stage-two indices,
  compatible context binding and singleton-tail packet scheduling. Preserve exact
  arithmetic, mathematical v1 identities and accepted serial PRAC eligibility.
- Avoid per-group cached-index copies. Charge preparation and scratch, retain
  finite cache prefixes, and use exact existing paths on refusal/recovery.
  Resource/checkpoint v6 rejects older resource identities; defaults stay off.
- Add independent traversal/term/budget oracles, interrupted/refused preparation,
  schema/version checks and real spawn/drained restart tests. Full committed QA
  passes 1,046 tests (9 optional skips) and lint.
- Record warm component gains alongside inconclusive fresh CLI comparisons,
  incomplete factoring outcomes and the stopped larger-quota checkpoint refusal.
  No complete-factor speedup, engagement/default adoption or B12/F4 closure.

## B12/F4 bounded reuse and worker imports (4 November 2026)

- Add opt-in, separately funded decoded and compact metadata retention for
  streamed schedules, bounded parent manifest reuse and lightweight worker
  imports using the shared configuration implementation. Curve state stays private.
- Reuse existing reduced point kernels optionally in stage two. Preserve exact
  formulas, mathematical v1 streams, accepted PRAC guards and serial defaults;
  resource/checkpoint v5 rejects older resource contracts.
- Add independent cache/refusal/ownership tests and bounded common-seed parameter
  replay and parent-attribution tools. Full committed arithmetic QA passes 1,026
  tests (9 optional skips) and lint; diagnostic timing refusals retain elapsed
  durations. Record scoped cold/warm gains, losses and finite charges separately
  from unrun first-factor/full-factor/fresh-confirmation gates. B12/F4 remain open.

## B12/F4 wider-worker contracts and diagnostics (4 November 2026)

- Permit explicit, funded spawn configurations through 16 workers, retaining
  serial automatic defaults. Counts above ten are oversubscription experiments.
- Version the resource/checkpoint contract as v3 while preserving mathematical
  v1 curve identities. Normalize durable position keys before hashing and reject
  noncanonical aliases; retain old source pins rather than migrating checkpoints.
- Leave CPU for parent publication when issuing process leases, and expose the
  actual traversal termination reason even when an earlier factor was found.
- Add bounded CPU, GC, heap and RSS diagnostics, independent result oracles and
  cold/warm worker-count analysis. Full committed runtime QA passes on `a32b9f0`;
  wider scaling remains exploratory and B12/F4 remain open.

## B12/F4 larger-work headroom screen (4 November 2026)

- Add a pure receipt analyser with independent corruption/coverage tests.
- Record exploratory cold/warm scaling on 60/80-digit enriched inputs under
  larger fixed workloads. Four workers amortize startup in these observations;
  complete-factorization benefit and engagement thresholds remain unestablished.
- Retain checkpoint/CPU-lease admission refusals and finite charges. Executor,
  arithmetic and automatic defaults are unchanged; B12/F4 remain open.

## B12/F4 runtime QA and initial screen (3 November 2026)

- Pass 987 committed-checkout tests (9 optional skips), full lint and bounded
  startup/fixed-work/cold/warm screens; preserve failed QA receipts and charges.
- Correct two test fixtures, add chained A7 compatibility source identities,
  and include external result validation in the cold CLI measurement endpoint.
- Retain serial defaults: cold startup dominates the short screened prefix;
  warm scaling remains exploratory with engagement calibration pending.

## B12/F4 experimental preparation (3 November 2026)

- Add an opt-in spawn executor for independent ECM curves, versioned manifests,
  parent leases, checked publication and drained restart; retain serial defaults.
- Overlap funded worker startup with per-child deadlines and conservative
  failure accounting; add queued startup ownership tests. No speedup is claimed.
- Add opt-in bounded packets, deferred polling, private metadata/decoding reuse,
  pool-wide shutdown, checked stage-boundary queues and bounded engagement tables.
  Preserve mathematical v1 curve identities; version the extended wire/checkpoint
  contract as v2. Explicit wider-worker admission adds no measurement allowance.
- Add queued correctness/failure tests and prospective finite benchmark support.
  Arithmetic/process/full-suite QA and performance measurements are NOT RUN;
  B12/F4 remain open and this branch is not adopted or merged.

## C9a useful-yield continuation (4 November 2026)

- Remove the uncalibrated 15% canary veto from new C9a analyses, observer
  comparisons and feasibility admission. Retain legacy flags as diagnostics
  and provide explicit historical analysis replay. Keep instrumented runs and
  actual resource/recorder failures distinct; future performance judgments
  require measured variability and uncertainty rather than a universal cutoff.
  Post-correction committed-only mainline passes 1,114 tests (9 skipped), lint
  and all 243 benchmark imports; measured factoring arithmetic is unchanged.
- Integrate the tested C9a bundle into mainline and close the delivered stream
  by explicit direction. Keep the performance verdict exploratory and policy
  classes opt-in. Replace C9a in the ranked queue with C9b for current-mainline
  timing validation and a bounded remaining-headroom analysis.
- Mark bucket recovery metadata in bounded 60-bit integer words with block
  epochs; materialize exact wide masks only for selected candidates. Charge
  construction, periodic stamp copies and internal reads, keep mandatory
  polynomial-switch scratch funded, and retain the ordinary-list fallback.
  Checkpoint formats and mathematical collection order are unchanged.
- The extended 70-digit comparison completes both policies within 9,000s each:
  C9a/C0 take 1842.1/2742.1s on shared source `891eda5`. The final same-C9a
  hit-metadata full60 pair takes 159.5/132.2s on `ebb35fa`/`c8168dc`.
  Mixed/failed timing screens and one exposed modulus per band limit these
  raw differences to exploratory observations. No default adoption follows.
- Validate the final committed-only checkout with 1,113 tests (9 skipped),
  lint and all 243 benchmark imports. Calibrated whole-factor sampling and a
  complete mature heap inform the implementation; they are diagnostic data.

## C9a useful-yield preparation (3 November 2026)

- Prepare bounded complete prime-power root transfer across same-A
  polynomials, with independent uncached/modular-enumeration oracles. Preserve
  root traversal order and exceptional lifting; charge inverse construction,
  retain polynomial-switch scratch and discard derived caches on restore.
  Correctness and performance gates remain separate from adoption.

- Fuse exact dense hit/score updates, skip initialization overwritten by
  periodic scores, and replace repeated sparse-queue division with additive
  offsets and grouped reservations. Preserve event order, complete charges
  and cooperative polling. Avoid tiny root sorts and redundant lifting
  arithmetic while retaining independent polynomial-root verification.
  Benchmark warmup now retains every selected marking kernel.

- Remove repeated compact-factor decoding in private proof checks and exact
  canonical keys; retire obsolete raw proofs after packed publication. Reuse
  bounded Hensel inverse information and recovery masks across polynomials,
  hoist fixed-target queue eligibility, and batch scalar reject accounting
  with unchanged finite polling and exact refusal/resume. Keep the existing
  recomputation default and automatic routing pending matched measurements.

- Validate the final flat-queue bundle with 1,084 tests (9 skips), lint and 242
  committed benchmark imports. The matched exposed 60-digit incumbent/C9a pair
  finishes in 268.8/230.4s; both timing screens fail, so preserve the raw 14.3%
  difference without adoption. Keep fused scanning optional after mixed
  prefix evidence and retain every failed/rejected measurement.

- Replace slower packed-field queue replay with measured flat unsigned-word
  triples. Add the separate experimental `scan-report` nominee to activate
  exact fused scanning while preserving existing kernel and admission choices.

- Implement shared SIQS polling, checked same-A metadata reuse, packed sparse
  queues, redundant-refinement avoidance and compact residual coordinates.
  Preserve immediate work refusal, cancellation publication boundaries and
  independent checkpoint validation. Measurements remain pending.

- Complete the shared progress-path 60-digit pair: raw C9a time 320.9s versus
  incumbent 356.3s, with both timing gates failed. Add bounded CPU/heap evidence
  and identify resource-observer cost; retain exploratory scope without a
  production observation-policy change or adoption.

- Avoid copying cumulative atom/relation arrays and expanding partial IDs on
  every built-in QS/SIQS collection batch. Share the same traversal, exact
  publication and refusal/resume path with public collection; preserve full
  public snapshots, collector overrides and serialized checkpoint formats.

- Remove pivot sorting with an incrementally maintained presence bitset, avoid
  unnecessary panel-table invalidation, and bound
  cofactor lookahead checks in prepaid chunks. Add a rho-first split policy with
  retained SQUFOF fallback and reuse of paid splits for three-prime recovery.
  Prepare the separate `cost-report` nominee using scalar elimination and that
  split order; preserve earlier profiles and all exact validation contracts.

- Reject proven single-prime domain holes before hints/primality testing.
  Skip rho-first three-prime recovery below the proven residual-prime floor
  cubed while retaining the final useful two-prime rho walk, and
  validate narrower split memos on restore. Reuse batch scratch high-water
  bounds, avoid repeated ownership calculation and disabled recomputation
  scans, cache singular/nonsingular Hensel derivative results, and retain the
  dense-root view across sparse-queue blocks.

- Publish bounded opaque staging with only actual FIFO evictions instead of
  copying its map. Skip canonical-key construction for unmatched staged atoms
  after exact validation and the applicable nonunit check. Budget refusal
  preserves staging; interrupted allocation requires a prior checkpoint.

- Select serial recovery in the experimental `cost-report` nominee after an
  exact-prefix check and exposed complete-call diagnostic. Preserve shared
  batching as an option. The candidate still loses to the incumbent on this
  input and timing screens fail; no adoption or default change follows.

- Complete a larger bounded screen: four independent exposed inputs, eight
  successful calls at 50/55/60/62 digits. Candidate sieve coverage is lower,
  but raw complete time and RSS are greater in every pair. Preserve failed
  calibration gates and retain the incumbent pending targeted cost work.

- Complete the repaired 60-digit matched pair under the original 900-second
  and memory grants. Both factor successfully; C9a uses fewer sieve positions
  but greater observed complete time, RSS and charged work. Preserve failed
  timing gates and stop after the pair without adoption.

- Audit the actual 26,862-atom ownership stop. Charge shared slabs and polynomial
  contexts once, integer digit storage explicitly, and packing scratch separately.
  Rebuild compact in-memory records on restore and preserve cold-file identities.
  The same captured state replays at 367.9MB accounted workspace with matched
  helper allowances, versus 1,040.1MB previously. This is an ownership correction,
  not a measured RAM reduction, factoring speedup or default promotion.

- Extend one matched 60-digit comparison to 900 seconds per arm under unchanged
  memory grants. C0 factors the input; C1 stops at accounted memory before a
  matrix solve. Preserve failed timing screens and the remaining storage limit;
  stop after the pair without adoption or a default change.

- Add a separately versioned bounded-yield variant with capped opaque staging,
  exact promotion, shape-based incidence ownership and norm-bounded batch
  scratch. Preserve FIFO state in schema 14 and retain the incumbent reporting
  allowance in the new experimental profile. Fix sparse incidence ordering
  after light-row selection and preserve completed divisors across restore.
  The four-call rerun finishes matched prefixes under unchanged grants; early
  closed-row yield matches C0, charged work is greater and drift gates fail.
  Keep the variant experimental and stop before a larger campaign.

- Make the measurement warmup exercise selected coupled arithmetic. Complete
  the four-call matched feasibility screen; retain candidate ownership-cap
  stops and failed drift checks without a speedup claim. Pause larger runs
  for opaque-retention and memory-accounting review.

- Merge current C10a-01 mainline policies into the C9a branch. Preserve
  coupled arithmetic markers alongside yield schemas 11–13 and keep historical
  control decoding explicit; add combined checkpoint/resume regression cases.

- Fix expanded SIQS result statistics to report residual incidence without
  ordinary-forest attributes. Keep tuple-based legacy loaders independent of
  opt-in residual/storage imports. The PyPy unit suite passes 942 tests with
  9 skips after these repairs; performance acceptance remains pending.

- Prepare the separate expanded arm with lazy cycle arithmetic, integer
  incidence/provenance, three-prime and opaque residual contracts, bounded
  conversion/reduction, adaptive reporting and extraction-outcome feedback.
  Add recomputation/spill capabilities, finite asynchronous hints, sparse and
  structured filtering, a bounded panel solver, schema13 restore and a v4
  two-arm/parameter-selection protocol. New executable validation is pending;
  this records opt-in preparation, not acceptance or a performance result.

- Extend the opt-in combined challenger with occupied-column feedback, private
  checked-proof reuse, exact early rejection, bounded cofactor splitting and
  batch hints, direct/bidirectional cycle shortcuts, canonical congruences,
  base-index slabs and integer provenance handles. Add fill-aware filtering,
  reversible lighter-row attempts, joint extraction and schema-12 restore.
- Prepare a separately identified two-arm protocol with a four-call clean
  feasibility checkpoint. Expanded correctness checks are written but their
  execution is deferred until the shared quiet window; no new campaign ran.
- Add independent opt-ins for useful-row preparation targets, asymmetric
  proven large-prime endpoints, shorter exact cycles, immediate equal-parity
  dependencies and immutable packed exponent records. Preserve separate
  reporting/admission controls, original provenance and finite reservations.
- Add policy-bound schema-11 checkpoints, independent correctness oracles,
  an immutable mainline control, certified exposed inputs and bounded staged
  experiment tooling. Automatic defaults and routing are unchanged; no C9a
  performance campaign or adoption is included.

## C10a-01 mainline and automatic defaults (3 November 2026)

- Integrate exact parity polynomials, coupled odd/even multiplier ranking,
  reachable paired A selection, bounded modular root tables and guaranteed
  powers-of-two marking. Fresh automatic SIQS configurations at every size use
  the selected bundle, composing existing collector kernels. Exact standard
  fallbacks, explicit configurations and legacy checkpoints are preserved.
- Retain existing size parameters, finite grants and portfolio eligibility;
  this adds no calibrated 61+ preset. The all-size promotion is an explicit
  user-directed decision before fresh confirmation. Record exposed 60/62/65/70
  completions, the 65-digit host-state diagnosis and failed canary gate, and
  the user-censored 70-digit control in [results](v2/benchmarks/RESULTS.md#c10a-01-coupled-polynomial-defaults).
- Preserve source freezes, failed outcomes and research in the verified local
  archive before removing the task worktree and temporary state. No push or E1
  campaign is part of this integration.

## Study publication cleanup (2 November 2026)

- Remove the selective `v2/studies/` report/result collection from the current
  tree after exact local archival. Keep concise decisions in benchmark results
  and link detailed historical records through pinned Git revisions.
- Retain tested C10a-04 schedule/decision helpers under benchmarks. Retired
  campaign orchestration replays from its original revision; numerical inputs,
  frozen controls, production behavior and experimental verdicts are unchanged.

## A1b mainline evidence integration (2 November 2026)

- Integrate the completed fixed-score review, frozen runners/inputs and exact
  growing-DLP regression checks. All 30 prefix calls on nine exposed moduli
  reject the existing fixed/list bundle for advancement; retain current defaults
  and fixed opt-ins. No complete-factor speed claim or production change.
- Preserve original failures, missing clocks, source freezes and the explicit
  accounting repair. Update status, live backlog and results navigation without
  restoring the superseded documentation layout. See the
  [final decision](v2/benchmarks/RESULTS.md#a1b-existing-fixed-score-review).
- Verify committed-only tests/lint/imports and unchanged production source;
  archive local state before removing the A1b worktree and temporary QA state.

## README navigation cleanup (1 November 2026)

- Shorten the root and v2 overviews; put CLI, API and checkpoint instructions
  in focused guides, with package READMEs for algorithms and shared execution.
- Move benchmark measurements into a results guide, repair stale routing and
  module references, and retain prior detail through pinned Git history.
  Source, defaults, tests, corpora and experimental verdicts are unchanged.

## Roadmap documentation cleanup (1 November 2026)

- Split current status, live backlog and stable policy. Preserve historical
  quantitative narratives through a pinned Git revision and source/hash metadata
  in the [history index](https://github.com/nishanth17/factor/blob/23458a52e409add326c9fd63cf1f322c6e605569/v2/history/README.md). Remove the
  duplicate full snapshot and use Git history. No experimental verdict changed.

## C10a-07 combined defaults (1 November 2026)

- Adopt swapped adjacent sign coordinates for fresh
  automatic SIQS jobs and exact sparse prime/power queues in eligible collectors.
  Preserve A1a kernels, existing parameter/routing bands, total owned memory,
  explicit policies, legacy snapshots and charged exact fallback/resume.
- Record 23.1% observed wall / 23.0% CPU reduction on nine fresh 50/55/60 inputs as
  exploratory: the original drift gate remains unmet. The additional 18-call
  confirmation was not run; adoption is a separate policy decision.
- Close C10a-07 with findings and verified production integration. Preserve
  receipts/research/commits before removing temporary worktrees and scratch.
  See [adoption](v2/benchmarks/RESULTS.md#c10a-07-face-and-queue).

## A1a mainline adoption (31 October 2026)

- Select `density+presieve+scan_fused` for fresh native
  50–55-digit jobs, inclusive, with the default 1 GiB owned allowance. Preserve
  explicit configurations, other resource scopes and saved selectors.
- Keep selected-power periods and four-stride marking opt-in after their
  combined comparison lost at both sizes. Close the delivered integration
  and comparison roadmap item with its results; keep setup/root reuse and
  conditional two-root interleaving as future work.
- The 50-digit integration gate passed; 55 digits remains inconclusive due
  to drift, and 51–54 are policy extrapolations. Preserve the original
  evidence and source identities locally before removing task worktrees.
  See [adoption and replay](v2/benchmarks/RESULTS.md#a1a-integration-and-extension).

## A1a frozen comparison decisions (31 October 2026)

- Complete both separate 50/55-digit comparisons on fresh inputs. Accept the
  selected integration's scoped 50-digit complete-cost gate; retain control
  at 55 digits because drift persists after its prescribed full extension.
- Retain the prior integrated control at both sizes for the combined
  selected-power-period/four-stride extension: it increased complete cost.
  Preserve all 180 measured calls, original ledgers and frozen checkouts.
  [Report](v2/benchmarks/RESULTS.md#a1a-integration-and-extension) gives uncertainty, per-input
  results and protected accounts. Automatic defaults and E1 stay unchanged.

## A1a bounded collector extension preparation (30 October 2026)

- Add opt-in selected-power periodic initialization and four-stride per-root
  marking with exact recovery, bounded scratch and charged checkpoint rebuilds.
  Preserve the original A1a implementation and its frozen comparison.
- Prepare a separate combined comparison against the prior A1a source at
  50/55 digits. Performance acceptance and automatic adoption remain pending.
  Track quota-aware setup, reusable switch metadata and lifted-power reuse
  under A1a; defer same-modulus root interleaving to conditional A1c.

## A1a selected-kernel integration preparation (29 October 2026)

- Add explicit density, bounded periodic presieve and fused-scan selectors,
  preserving exact recovery, charged cancellation and saved-policy replay.
  Port the accepted combination without the other experimental kernels.
- Reconcile counted-base setup and historical source loaders. Prepare a frozen
  comparison against committed pre-integration mainline on fresh 50/55-digit
  inputs. Accepted timing and automatic adoption remain pending the manual
  execution handoff; existing defaults and E1 status are unchanged.

## Publication cleanup — 12 October 2026

- Keep detailed benchmark research and comparison reports in the local audit
  archive. Ignore external reference sources as a complete folder, preserving
  their source files, notices and licenses locally.
- Retain required numerical inputs and immutable baselines. Runtime C1 loaders
  validate numerical inputs without requiring local prose protocols; historical
  protocol hashes remain unchanged. Rename the optional comparison harness
  with an explicit source-migration record. Production arithmetic is unchanged.
- Validate the candidate checkout without ignored files: 743 PyPy tests pass
  with nine optional-backend skips; all 179 benchmark modules import. Lint
  passes, and the production AST is unchanged.

## C10a-04 default promotion (12 October 2026)

- Enable native growing-DLP SIQS defaults through 60 digits,
  preserving accepted 50/51 settings. Use charged exact 2500-entry compact bases
  at 52–55 and 5000 at 56–60, with finite budgets and saved-policy replay.
- Record preliminary and extrapolated coverage in the
  [promotion report](v2/benchmarks/RESULTS.md#c10a-04-larger-band-routing); close the delivered
  roadmap tranche without claiming the original campaign completed. E1 remains
  deferred, C10a-07 open, and 65/70 feasibility separate.

## Bounded C10a-04 follow-up (12 October 2026)

- Add a separate 60-digit compact/control study with fresh certified inputs,
  a ten-hour total budget and fixed preliminary early-stop rules. Preserve the
  original incomplete campaign, interrupted reservations and production
  defaults. This is experiment tooling, not a performance acceptance result.
- Record the [fresh preliminary result](v2/benchmarks/RESULTS.md#c10a-04-larger-band-routing):
  the fixed six-pair trend stop passed, with 6/6 compact completions against
  5/6 control and 36.16% lower capped time. Preserve the timeout and full
  reservations; keep production acceptance and promotion separate.

## C10a-04 calibration study freeze (11 October 2026)

- C10a-04 Stage 1 adds an isolated, source-pinned coupled SIQS calibration
  study, exact eligible-count resolution, durable cumulative accounting and
  interruption/resume validation. It freezes a Stage 2 experiment without
  promoting defaults, reopening completed studies or closing E1.

## Descriptive ECM preset names (11 October 2026)

- Replace upstream-named preset keys with `pretest8`, `escalating140`,
  `ladder115`, `ladder415`, `single74` and `single214`. Previous keys are
  rejected; numerical schedules, defaults and checkpoints are unchanged.
  Retain source provenance and historical evidence in benchmark documentation.

## C10a-05 larger-band default scope (11 October 2026)

- Extend the accepted native 1 GiB C9+C10+collection bundle to exactly
  51-digit inputs after a separate frozen study; preserve explicit settings,
  saved policies, 40/50-digit behavior and the previous 52+ default.
- Nine matched training groups on three certified 51-digit products found
  9/9 complete calls in each candidate arm versus 0/9 in the prior default.
  Collection pacing lowered complete training cost 28.19% versus its otherwise
  identical window-paced arm (95% input-cluster interval 27.43–29.07%).
  Nine held-out matched calls completed 9/9 with the selected bundle versus
  0/9 with the prior default. This confirms completion coverage, not a
  speedup against successful baseline calls.
- A revealed 55-digit product completed only under a longer bounded grant;
  a 60-digit product was censored at the existing 1,200-second automatic
  grant. Retain 52+ routing pending separate calibration. Preserve immutable
  C10a-05 controls and A7 history through a new chained source manifest.
  See the local study record.

## C10a-05 collection and solve pacing (11 October 2026)

- Separate shared SIQS polynomial-window completion from genuine collection
  exhaustion. Preserve fixed cadence/surplus, immediate checked square
  dependencies and a charged terminal attempt on retained originals. Keep
  pending solver/extractor state, cumulative allowances and finite stopping;
  schema 6 persists the new policy while legacy checkpoints keep window pacing.
- Adopt collection pacing only for the accepted native 50-digit bundle at
  exactly 1 GiB owned memory. Nine fresh matched pairs on three certified
  inputs reduced complete cost by 25.99% (95% interval 22.87–28.24%), with
  9/9 validated completions per arm. Retain 40-digit window pacing after its
  inconclusive eighteen-pair extension; small/30-digit controls also stay.
- Record six-source executable scheduling research, independent dependency
  witnesses, repeated serialized resume and finite benchmark leases. Preserve
  historical studies and original atoms. Close only C10a-05; C9a feedback,
  other C10a mechanisms and E1 combined acceptance remain open. No compounded
  portfolio speed claim. See the local study record.

## C10 bounded tranche integration (11 October 2026)

- Route omitted native 40/50-digit SIQS through the frozen growing-DLP bundles
  with C10's accepted 2 MiB power plan at the measured 1 GiB allowance.
  Preserve explicit configurations, smaller memory/work/time grants, legacy
  checkpoints and the reference `SIQSConfig()` default. Other size bands and
  backends retain their previous policy.
- Compose with C9's 50-digit selector without changing either frozen study.
  Verify the production configurations against C10's versioned plan controls;
  retain separate C1/C9/C10 timing comparisons and their balanced-input limits.
  Close this delivered C10 tranche and phase distinct untested collection
  mechanisms into conditional C10a gates with their original owners.

## C10 collection strategy tranche 1 (11 October 2026)

- On the settled C9 growing-DLP control, the isolated existing 2 MiB
  prime-power plan passed fresh complete-factor comparisons at 40 digits
  (21.83% lower wall cost; 95% interval 12.47–24.88%) and 50 digits (17.79%;
  interval 16.11–19.48%). Each band had 9/9 validated completions per arm on
  three independent certified inputs, with stable matched samples and about
  2 MiB additional peak owned workspace. The 30-digit control remains selected;
  no global default or combined portfolio policy changes.
- The fixed 16-family exact A-quality ordering did not advance from training.
  Record the scoped loss separately from open product-population and family-
  diversity mechanisms. Finish a nine-implementation source reconciliation,
  normalized/mod-8 arithmetic oracles and C10a conditional follow-ups while
  preserving A1a, C8 and C9 ownership. Keep C10 open after this bounded tranche.

## C10 collection strategy preparation (11 October 2026)

- Add Python-only experimental SIQS `assignment_policy="quality"`, ordering
  fixed batches of flyer families by exact A-product distance. Preserve full
  batch populations, extendable prefixes and distinct checkpoint identities;
  reserve cache coexistence and charge reconstruction for finite/growing jobs.
- Add independent normalized-polynomial/mod-8, root, reconstruction and
  declared-threshold-loss oracles, a pinned nine-implementation mechanism
  comparison, and a bounded isolated policy experiment harness. Existing
  defaults remain unchanged; the completed scoped decisions are above.

## Folder organization (10 October 2026)

- Group v2 algorithms under `ecm/`, `pm1/`, `rho/` and existing `qs/`;
  place shared helpers in `common/` and budgets/dispatch in `execution/`.
- Group benchmark runners, builders and research notes by subject and study;
  preserve immutable inputs and move the detailed study record behind an index.
- Record exact structural source migrations separately from frozen evidence;
  update module commands and loaders without changing factoring algorithms,
  defaults, work accounting or checkpoint formats.
- Consolidate ignore rules and place local journals under `v2/local/`.


## Current development

- Promote C9's confirmed native 50-digit growing-DLP bundle by default.
  Route omitted SIQS/automatic/CLI selections through the
  exact accepted configuration and finite 1 GiB / 1,200-second / 10^13-work
  grant; preserve explicit configurations, limits and concrete checkpoint
  settings. Measured benefit covers the balanced confirmation cohort;
  arbitrary 50-digit classes and broader combined E1 acceptance remain open.

- Integrate accepted C9.2 with C3 allocation/defaults. Preserve growing-DLP
  policy during implicit schema-12 resume and retain the finite resolved child
  allowance; reject changed policies. Cross off C9 and phase remaining
  literature transfers into C9a with separate attribution/acceptance gates.
  Integrated working/committed-only checks pass 690 tests, full lint and
  140 live benchmark imports plus required loaders.

- Add `GrowingDoubleLargeSieveConfig` for serial SIQS with allocator-grown
  partial history and no implicit FIFO, original-atom, row or cumulative
  split-count ceiling. Standalone growing memory/checkpoint quotas default to
  `None`; explicit caller quotas and finite portfolio ownership remain.
- Preserve proven endpoints, all-component cycles and exact original-atom
  square corrections. Schema 5 supports checked growing-store resume; legacy
  configurations retain schemas 3/4. Interrupted publication invalidates the
  live object while prior complete checkpoints remain recoverable.
- Complete the C9.2 source-derived yield/adaptive-storage tranche. Pinned
  literature and five leading implementations motivate history retention and
  separately tested endpoint alignment. Fresh balanced 50-digit confirmation
  reduces complete-call time by 38.10% versus calibrated SLP and 53.08% versus
  fixed DLP, with 9/9 validated completions in all arms. Adopt only the explicit
  measured bundle; its owned peak is 326 MiB versus 167/198 MiB for the controls.
  Retain 30/40 recommendations and defaults. The C9.2 report records costs,
  uncertainty, stopped studies and limitations; C10/E1 remain open.

- Integrate the delivered C3 allocation/handoff tranche and
  30/40-digit ECM defaults into mainline. Cross off the bounded C3 roadmap
  entry while retaining the unrun fresh/stability gates and broader joint
  calibration under the continuing follow-ups.

- Promote C3's size-scoped bounded defaults: compact64
  at 30 decimal digits and the 10/30/100 escalating prefix at
  40 digits. Preserve explicit tiers, numerical legacy resumes, all cumulative
  limits and SIQS settings; retain fixed32 for other sizes or tight workspace.
  This promotion precedes fresh confirmation and leaves calibration open.
  Committed-only checks pass 663 tests, lint and benchmark/input loaders.

- Complete the initial 54-group upstream ECM transfer screen: all 324 calls
  validate under unchanged SIQS settings and reserves. The escalating140 and ladder115 schedules show
  44.44%/47.39% lower revealed-training CPU cost
  than fixed32, with wide uncertainty. Stop this pass over
  runtime concerns; defer extension/fresh confirmation and retain defaults.
  Committed-only checks pass 653 tests, lint and all benchmark/input loaders.

- Reduce C3's next policy-screening scope after stopping the oversized
  513-group comparison. Preserve its 72 complete matched groups as
  inconclusive and retain the interrupted grant. Prepare a separate 54-group
  screen of four upstream-derived schedules against both current controls,
  with fixed SIQS bounds and finite uncertainty extensions; no policy gain is
  inferred from partial observations.

- Add explicit numerical ECM tier presets, retaining the caller's SIQS configuration, cumulative allocation,
  positive fallback floors and memory cap. Refuse incompatible workspace
  instead of reducing SIQS storage. Preserve numerical defaults; performance
  transfer and policy acceptance require matched measurements.

- Continue C3 after the negative quick8 result with bounded schedule
  and probability/cost research, 57 revealed training subjects and a validated
  seven-bundle attribution pilot. Reduce duplicate reserved-charge clock
  checks and reuse the allocation view, preserving the ordinary budget path
  and shared work ledger. Differential boundaries, 600 tests and lint pass;
  warmed policy selection and fresh confirmation remain open.
- Complete the frozen 2,052-call instrumented stopping calibration and commit
  the fitted compact64/wide64 table before comparison. Add complete matched
  group reservations, warm both size bands and SIQS paths, retain interrupted
  grants and preserve the existing comparison cap and decision gates.
  Prepare independently checked prime certificates for possible later fresh
  confirmation; no new policy or performance acceptance is claimed.

- Add opt-in cumulative ECM pretesting and finite campaign allocation, with
  protected work/wall/CPU admission floors for SIQS/SSS, lazy schedule setup,
  one-way handoff and partial-attempt evidence. Preserve shared child budgets,
  cancellation, exact reconstruction and probable/proven labels.
- Use schema 12 for explicit allocation identity and complete-config restore;
  charge snapshot verification/rebuilding and preserve admitted live fallback.
  Expose policy/tier/reserve CLI flags and retain schemas 2–11 and defaults.
- Complete pinned source/paper/blog research, a six-bundle frozen training
  study, commit-before-generation selection and 15 fresh certified controls.
  All 972 fresh calls complete; quick8 costs 51.56% more on the smaller cohort,
  with uneven regressions and remaining instability, so retain both numerical
  controls. Default-grant probes expose unfunded reserves rather than claiming
  early-refusal speed. Charged restore, 593 working/committed-only tests, full
  lint, 144 imports and 34 certified corpora pass. Keep broad C3/G1/E1
  calibration open; preserve v1 and C9/C10 collector ownership.

- Enable native reduced PRAC/batch 16 for supported fresh bounded ECM jobs.
  Keep B4 for unsupported jobs and insufficient explicit
  caps; select the separately approved finite 16 MiB automatic workspace only
  where supported. Preserve earlier construction/regression evidence.
- Add optional `--ecm-chain lucas|cf` and the `ecm_chain_family` API setting;
  `--ecm-chain off` retains the ladder. Keep eager certified recovery and finite
  run-local ownership. Default schema 10 and optional-family schema 11 pin
  execution identities; implicit resume restores saved caps/policy and legacy
  schemas remain compatible. No new search or arithmetic kernel is introduced.
- Freeze current-mainline 32-curve native regression controls and report charged
  portfolio and pause/resume comparisons separately from C6; retain the prior
  default-budget accounting diagnosis. Keep native claims separate from
  GMP and cold startup. Working and committed-only checks pass 567 tests, lint,
  benchmark imports and required catalogs. Close the bounded B3/P4.1 entry;
  broader allocation and population/default calibration remain C3/E1 work.

- Diagnose B3 recovery preparation and plan lifetime with matched PyPy
  ablations, fresh confirmation and exact original-harness replication.
  Retain accepted eager production/defaults: lazy setup has no supported
  full-run benefit. A charged finite cohort owner confirms a 2.48% native
  saving versus the ladder; public cross-call ownership remains a separate
  API gate. Identical Lucas execution shows native integers outperforming
  GMP here. Preserve all prior evidence and the unpromoted lazy fixture;
  working/committed-only checks pass 501 tests, full lint and 106 imports.

- Complete B3/P4.1 bounded PRAC/Lucas production integration with an explicit
  reuse option for B1=2,000, chunk 16, tiers of at least eight curves and
  40–80-digit cofactors. Native reduced PRAC and separate GMP tuple Lucas
  reuse C6 certificates; fresh/unsupported routes and all defaults retain B4.
- Bound run-owned plans, charged preparation/verification/eviction, simultaneous
  scratch and worst-case strict/saturated recovery. Schema 10 pins the opt-in
  identity and charges empty-cache rebuilding; schemas 2–9 remain compatible.
  Validate divisors, certified chunks and every unresolved reconstruction.
- Record production break-even and fresh equal-schedule portfolio/resume
  confirmation, including native whole-portfolio regressions and unequal-search
  rejection of initial capped gains. Retain the native default; keep GMP
  claims separate. See the B3 benchmark summary for intervals and limitations.
  Final and committed-only PyPy/GMP checks pass 499 tests and full lint.

- Integrate bounded two-large-prime SIQS as a configurable opt-in through the
  Python API and `--qs-dlp` CLI flag. Preserve single-large-prime defaults and
  charged DLP checkpoint/resume. C1's fresh balanced 40-digit complete-run
  comparison supports only the declared explicit bundle; wider tuning and
  portfolio/SSS dispatch remain open.

- Add experimental opt-in two-large-prime serial SIQS after the bounded C1
  residual study produced complete 40/50-digit witnesses. Preserve the failed
  literal gate and explicitly identify the corrected cost-attribution go as
  post-observation; no speedup or default promotion follows from that study.
- Preserve scalar residual certainty with separate proven endpoint pairs,
  finite splitting and all-component graph cycles. Retain exact original
  atoms through FIFO eviction and use charged DLP-only SIQS checkpoint4;
  default SLP remains checkpoint3. SSS/worker collection rejects the opt-in
  config. Independent complete graph-kernel tests cover disconnected cycles,
  repeated primes, loops, duplicates, exact corrections and retained ownership.
- Accept the bounded C1 40-digit balanced opt-in bundle after a 31.73%
  complete-run reduction against calibrated SLP (95% interval 29.21–33.99%),
  with unchanged completion. Retain the 30-digit positive observation with
  its fixed-seed drift limitation, all fresh 50-digit timeouts, and the
  separately completing 50-digit charged-resume witness. Keep SLP defaults;
  broader SSS/portfolio comparisons and larger-band promotion remain open.
- Reconcile A11 as complete using the existing three-mode QS/MPQS/SIQS CLI,
  matched library comparisons and nine passing CLI regressions. Correct stale
  roadmap status; automatic allocation and handoff remain C3/G1 work.

- Complete the bounded A7/P3.8-R5 interface/evidence reconciliation, preserving
  accepted repair decisions, immutable source/input pins and historical scope.
  Add both-mode SSS forced-valuation, interrupted/cumulative restore, authentic
  legacy migration and terminal/output validation coverage.
- Repair wide SSS solver fingerprints using the existing streamed hexadecimal
  encoding; keep explicit encoding identity and legacy decimal replay without
  changing the SSS schema version or old snapshot bytes.
- Fix B1 fresh-output validation to honor A10's wider deterministic certainty
  domains; reuse its complete-call adapter for separately labelled SSS arms.
  Prepare current calibrated controls and bounded worker arms for E1 without
  new timing or default promotion. Record worker schedule incompatibilities,
  missing 40-digit SSS calibration and fresh combined confirmation as deferrals.
  Flag A11's prior outstanding status as stale; automatic handoff stays C3/G1.

- Promote bounded p−1 recurrence and p−1-only chunk 64 to fresh defaults.
  Fresh complete stages save 9.73% CPU across
  nine matched cells; extended integer portfolio captures remain inconclusive.
  Record that limitation without claiming sustained portfolio improvement.
- Use portfolio checkpoint version 9 for the new executor; pin both settings,
  reserve finite table growth and charged resume reconstruction, and retain
  old execution automatically for old default snapshots and CLI resumes.
  Explicit legacy configurations preserve old work/state identity. Preserve
  ECM chunks, bounds, attempt/curve counts, exact integers and quiet calls.

- Integrate the completed C6 research modules, certified PRAC/Lucas/CF
  records, independent tests and pinned source controls into mainline for
  later experiments. Keep `ecm.py`, stage jobs, B4 arithmetic, checkpoints,
  production defaults and `v1/` unchanged by this merge. B3 still owns
  production routing, finite work reservations and resume behavior.
- Preserve both ignored C6 evidence sets in mainline before retiring the
  isolated checkout. The merged PyPy/GMP suite passes 479 tests and full
  lint. A committed-only archive of merge `ce7574e` passes the same 479 tests,
  full lint and 98 benchmark imports with required C6 catalogs. The C6
  worktree is retired; historical benchmark protocols retain their original
  source pins.

- Complete the reopened C6 executor optimization and isolated B4 comparison.
  Independently certify batched coordinate-factor checks with strict saturated
  replay; reproduce bounded CF-family minima without claiming global optimality.
- Confirm reused native reduced PRAC at 6.73% less complete-campaign time than
  pinned B4 [5.98%, 7.26%], or 3.98% less after charging first preparation.
  Separate GMP Lucas saves 8.42% reused and 5.97% with per-attempt construction.
  Preserve all 10/20 splits, unsuccessful attempts and unresolved cofactors.
- Retain the ladder for fresh native attempts and all production defaults;
  combined Lucas is inconclusive and larger CF search is not justified.
  Record generation, construction, bounded storage, 153 cold starts and
  separate profiles. B3 owns production work/resume/portfolio integration.
  A committed-only archive passes 431 PyPy/GMP tests, full lint and all 84
  benchmark imports; controls and certified inputs are versioned, raw evidence
  remains local. Preserve v1 and production kernels; C6 was still isolated
  at this study acceptance point.

- Record the initial conservative C6/P4.5 Lucas-chain comparison (historical).
  Preserve the pinned generator/decoder locally and 303 independently verified
  prime records; add reusable compact PRAC/Lucas, rolling and restricted
  three-point research executors with bounded nonunit/factor recovery.
- Retain the ladder: Lucas costs 2.03x native-int and 1.22x GMP complete ECM
  campaign time with identical completion on the frozen 40–80-digit cohort.
  Record generation, amortization, storage, cold costs and the conditional
  search stop. Apply the revised policy allowing confirmed sub-10%
  gains without rewriting the historical protocol. B3 production routing,
  budgets and checkpoint integration remain open; v1 and kernels are unchanged.
  A committed-only archive passes 411 PyPy/GMP tests, full lint and all 72
  benchmark imports, including frozen C6 inputs and independent certificates.

- Integrate accepted A6 bound-continuation and recurrence APIs into current
  B4 mainline at `62a57f2`. Preserve B4 arithmetic, direct p−1/portfolio/ECM
  defaults, bounds, attempts and checkpoint formats. Cross off the bounded
  A6 roadmap contribution; keep later ECM/Lucas migration and allocation open.
- Pass combined mainline test/lint and all 449 PyPy/GMP tests, lint, 82
  benchmark imports and A6/B4 frozen loaders from committed files. Verify all
  92 raw A6 evidence files before managed-worktree cleanup. Remove the A6
  checkout, fully merged branch and temporary test archives; preserve ignored
  mainline evidence, shared PyPy and unrelated C6 work. Historical A6 timing
  remains bound to `869ed42`; no combined B4/A6 speedup is claimed.

- Integrate B4's native selected-reduction ladder, confirmed at 6.79%
  complete-run saving (paired 95% interval 5.82–13.58%) in the frozen
  production bridge. Retain readable GMP arithmetic: whole-ladder fusion's
  bridge interval crosses zero. Select once per scalar action; preserve
  readable point formulas, exact `(A+2)/4` squared-difference convention,
  public APIs, native defaults, canonical coordinates and checkpoint formats.
- Extend independent field/composite/prime-power and bidirectional resume
  controls to production. Preserve finite work, cancellation and saturation
  recovery. Complete the bounded B4 roadmap entry; retain immutable research,
  source controls and certified corpora. Record final checks and the requested
  mainline integration/worktree cleanup in the benchmark acceptance record.
  Pass 414 tests, full lint, 76 benchmark imports and frozen loaders from
  committed files; merge to master and remove the B4 worktree/temp state and
  fully merged branch. Preserve v1 and unrelated C6/A6 work.

- Complete the requested fresh B4 bakeoff on a separately committed source,
  certified-input and sampling freeze. Reuse the same five candidates;
  counterbalance independent PyPy processes and extend only the unstable GMP
  whole-ladder comparison to 27 blocks. Fresh nine-block confirmation saves
  12.52% for native reductions (paired 95% interval 5.09–16.63%) and 6.05% for
  separate GMP whole-ladder fusion (1.52–17.61%), with 27/27 complete and
  matching work/results throughout. Both pass the revised scoped timing gate.
  Report noisy class losses and finite-workload limits. Keep production
  kernels/APIs/checkpoints unchanged; prepare integration review, no merge.
  Verify 414 tests, full lint, 75 benchmark imports and frozen loaders from
  committed files; release the shared performance window to C6/A6.

- Relax the roadmap's universal 10% time / 10-point completion promotion
  floor. Permit smaller repeatable complete-run gains supported by fresh
  confirmation and uncertainty beyond measurement noise; retain correctness,
  finite-resource, per-class regression and maintenance/risk checks. Preserve
  frozen experiment verdicts and make B4's positive native/GMP gains eligible
  for scoped integration review. No production change or automatic promotion.

- Complete the bounded B4/P4.2 arithmetic-kernel study on frozen mainline
  controls: explicit squares, fused step, whole ladder, selected reductions
  and unit-checked fixed-difference normalization. Pin primary research and
  implementation licenses; prove a24/normalization/width invariants and test
  independent affine/composite/prime-power, factor-recovery and resume controls.
- Retain the readable ECM kernel and native-integer default. Fresh full-run
  gains are 8.78% for native selected reductions and 4.73% for the separate
  GMP normalized arm, below the frozen 10% median gate; completion stays
  27/27 in every arm. Keep experimental code, certified corpora and immutable
  controls versioned, raw evidence local, and production APIs/checkpoints
  unchanged. Pass 409 committed-only PyPy/GMP tests, lint, 73 benchmark
  imports and all B4 frozen loaders. Stop after bounded comparison and fresh
  confirmation; no merge.

- Complete A6's reopened bounded p−1 bakeoff with bit-capped chunks,
  bounded even-gap recurrence and independently verified exact wheel/±
  relations. Accept recurrence with chunk 64 for scoped opt-in integration
  review: fresh 27-sample complete calls save 1.0–2.8%, aggregate 1.99%
  [1.85%,2.21%], with all medians and per-size-class intervals positive.
  Keep bit caps/wheels experimental; D=30/210 cost 76.2%/7.9% more.
- Version tuned work/state as `pm1-tuning-v1`, bind all configuration fields,
  charge table growth/recovery/reconstruction cumulatively, and preserve
  legacy arithmetic/work/checkpoint compatibility. The bounded legacy hook
  costs 0.78% CPU; existing production dispatch/APIs remain untouched.
  In-memory recurrence continuation saves 8.3–8.6%; verified resume costs
  18.1–21.8% more. The production-B2 portfolio bridge adds no completions
  over retained production and costs 3.6% more, so defaults/allocation stay.
- Broaden pinned literature/source/license coverage, including current
  ECM and special-form implementations, without copying code
  or claiming exhaustive coverage. Pass 14 additional adversarial tests;
  a committed-only checkout passes all 436 PyPy/GMP tests, lint, 73 benchmark
  imports and required certificate/input/selection loaders. Preserve local
  evidence and required versioned inputs; keep A6 isolated and unmerged.

- Add A6/P5.3's exact `inclusive-lcm-ratio-v1` integer schedule and opt-in
  `factorize_pm1_bounded` finite one-base campaigns. Increasing B1 includes
  higher powers of old primes and rebuilds stage two; a B2-only increase
  appends checked coverage. Stop after finite recovery, nonunit or exhaustion;
  never silently add bounds, bases or ECM curves.
- Bind the new schema-1 checkpoints to modulus/base/configuration, Python
  integers and `pm1-campaign-v1`. Reconstruct all saved actions and canonical
  numeric types under cumulative allowances before trusting arithmetic.
  Preserve proper divisors, unresolved reconstruction, certainty labels and
  quiet calls. Existing APIs, portfolio/ECM checkpoint versions, RNG and work
  accounting remain unchanged; the integer ratio unblocks later ECM scalar
  planning without implementing point/checkpoint migration.

- Apply the revised performance policy without changing frozen inputs or
  candidates. Accept chunk 64 for scoped p−1 integration review after all six
  27-repeat stage comparisons improve 6.1–9.2% with positive paired intervals.
  Preserve current defaults, streamed schedules, bounds and allocation.
  In-memory continuation saves 8.8–9.6% against fresh-each; verified resume
  costs 12.4–14.2% more. The declared portfolio cohort gains no completions.
  Pass 419 system-PyPy tests (three optional-GMP skips), lint, 422 tests from
  committed files with GMP enabled, all 70 benchmark imports and independent
  certificate/product loaders. Keep raw evidence local; leave ECM/Lucas
  migration and combined production promotion with their existing owners.

- Extend A10/P8.2.1 deterministic Miller–Rabin classification using the
  independently checked Sorenson–Webster bounds: first 12 prime bases for
  `n < 318665857834031151167461`, first 13 for
  `n < 3317044064679887385961981`. Preserve the faster smaller tests and
  reject certification at either strict endpoint by its preceding base set.
  The reported 23-digit factor is now proven; survivors beyond the final
  bound remain probable. Certificate generation stays B14.
- Share witness dispatch between direct and bounded classification. Fresh
  checkpoints pin `mr13-strict-v1`; old snapshots retain their original
  random witnesses, RNG/work accounting and conservative terminal labels
  across recursive children. Keep explicit probabilistic rounds and quiet
  calls. Document intentional fresh-run RNG/work changes and reject unknown
  policies or inconsistent retained evidence. Preserve v1; route its other
  capability candidates to existing roadmap owners.
- Retain the range-only implementation after frozen PyPy comparisons; no
  optional powering/early-exit/filter candidate met the complete-run gate.
  Wider-prime confirmation took 0.672 versus 1.946 ms per eight-input cohort.
  Complete-run medians improved in the declared prime/power-heavy workload,
  with remaining control-side tail spread disclosed. Version independent
  proof fixtures, immutable control sources, protocols and benchmark runners;
  keep raw evidence local. Pass PyPy tests and lint, including legacy/GMP
  resume checks. A committed-files-only archive passes 369 system-PyPy tests,
  372 GMP-enabled tests, all 59 benchmark imports and required input loaders.
  Certificate generation and combined portfolio promotion
  remain B14 and E1 respectively.

- Complete the bounded P4.3 int/GMP backend foundation across preprocessing,
  primality, rho, p−1, ECM, QS/MPQS/SIQS, SSS/SSSf, relation extraction and
  matrix bitsets. Add explicit backend selection, exact GMP helper contracts,
  canonical results and backend/build-bound checkpoints; preserve old integer
  resume and positional configuration APIs. Pending children retain the chosen
  backend after public splitters return canonical divisors. Missing GMP is an
  explicit error.
- Keep Python integers as the provisional default. The initial bounded ARM64
  PyPy study establishes no algorithm/size crossover: its QS fixture has only
  eight digits and its ECM bounds are small. In that declared scope,
  GMP uses approximately 1.8–2.0 times the portfolio time, with all
  factors, certainty, unresolved cofactors and work counts matching. Publish
  the reproducible runner, certified corpus, immutable control and limitations;
  keep captures local. The backend foundation passes 316 PyPy tests and lint;
  retain production-bound stage and larger QS measurements. Persistent mpz
  takes about 12 times native time in matched 50–100-digit bounded QS screens;
  a fivefold longer 50-digit run confirms the loss. Keep native loops as the
  provisional default. Further experiments stopped; disjoint
  confirmation and automatic selection remain open.
- Integrate the optional GMP boundary with verified PRAC and reusable ECM
  programs. Preserve native portfolio checkpoint schemas 4/5, use version 6
  for combined GMP identity, and accept pre-integration backend snapshots.
  Verify program resume on both backends and PRAC against an independent
  point oracle; retain the existing ladder and streamed defaults.
  The committed-files-only combined tree passes 358 PyPy/GMP tests, full lint
  and all 57 benchmark imports, with required corpus/control loaders verified.

- Complete the P4.1/A4 chain-correctness tranche and connect verified PRAC
  to `ecm.multiply_prac`. Separate bounded exact-rational generation from
  independent integer verification and checked point execution; retain a
  finite cache and one ladder retry. Expose proper factors or curve-retry
  outcomes through `prac.NonunitPointError`, never a degenerate point.
- Add independent affine, large-field, composite/prime-square and complete
  prime-power schedule controls, plus reproducible PyPy cost measurements.
  Compare full two-stage attempts on certified 40–80-digit composites with
  matched Python-int/gmpy2 ladder and checked-PRAC arms; validate every factor,
  unresolved cofactor and uncensored backend transition. Keep microbenchmarks
  separate from complete-campaign evidence.
  Preserve production ladder defaults. B3 integration remains open; defer
  reference ECM precomputed Lucas chains and compact execution to the separate
  P4.5/C6 experiment.

- Add P5.2 A3's opt-in immutable packed ECM prime/power blocks with bounded
  retention, charged regeneration and curve-private arithmetic. Preserve
  streamed defaults and version-4 snapshots; opt-in version-5 checkpoints pin
  schedule/configuration identity and retain cumulative resume allowances.
  Correct the resume cursor's odd-slot bound when prime 2 shares a tiny buffer.
- Add independent +/- prime-coverage certificates and affine cross-product
  oracles for the later paired executor. Document feasible explicit campaign
  storage/work and predeclared-curve continuation; keep paired execution,
  increased-B1 continuation, D tuning and allocation calibration open.

- Add independently certified 30–80-digit ECM workload fixtures and prepared
  finite deeper-campaign/large-bound runners with exclusive timing checks.
  Preserve exploratory completion records and verify a focused fixed-curve
  60/80-digit timing pass; defer broader timing and large-factor success claims
  after the expanded study ended early. Freeze the original
  control helpers separately without altering its immutable source snapshot.

- Expose QS/MPQS/SIQS through `--method qs|mpqs|siqs` after exact preprocessing,
  or `--siqs`
  as the rho/p−1/ECM portfolio fallback. Both reuse recursive factor validation,
  shared resource accounting and checked resume. Add finite `--qs-*` search
  flags (with `--siqs-*` aliases), help, examples
  and CLI regression coverage; retain existing automatic defaults pending
  measured dispatcher calibration. Track explicit CLI usability under P3.4,
  independently of automatic-dispatch performance gates.

- Condense benchmark history into dated README sections and add a QS guide.
  Organize required corpora, baselines and controls; preserve optional historical
  snapshots locally. Add explicit experiment-manifest selection and repair R3
  driver hash paths after the input relocation.

- Keep research notes and diagnostic tools in the ignored local audit tree.
  Move six required source/provenance inputs into benchmark inputs and retain
  the public acceptance plan as `v2/ROADMAP.md` (now reached through
  [current status](v2/STATUS.md)).

- Rewrite the repository and v2 guides around current methods, usage, result
  and resume contracts, measured comparisons and remaining work. Keep detailed
  experiments in the benchmark guide and label the historical emulated-v1 data.

- Organize required audit/benchmark JSON inputs under `inputs/` and ignore
  generated `results/` folders. Preserve historical raw evidence locally,
  remove disposable caches and update loaders and default output paths.

- Repair P2 preprocessing/cursor overhead, sparse QS exponent recovery,
  skipped-prime score refinement, repeated prime-power lifting inverses,
  sparse matrix label storage and dependency traversal.
- Amortize SSS collision accounting while preserving seeded candidates and
  exact completed-assignment work. Retain full-base SSSf smoothness testing
  after the disjoint-base candidate failed the measured improvement gate.
- Reuse verified QS preparation only for identical immutable rows, bases and
  provenance atoms under a bounded collector-owned cache. Keep public and
  checkpoint verification exact, and amortize native collector/solver polling
  with forced stage-boundary checks and the original shared work ledger.
- Add bounded family/Gray/block worker assignments and checkpoint version 2,
  retaining legacy complete-family resume. Poll external limits at explicit
  intervals of at most 64 atomic actions while checking every work reservation.
  Preserve central exact verification, aggregate CPU and charged replay.
- Complete the bounded P3.6.1 repair experiment gates. Fresh matched 30-digit
  complete cohorts improve by 70.7% for SIQS, 77.2% for SSS and 78.6% for
  filtered SSSf; all 72 timed attempts per arm complete. P2 small/medium
  complete cohorts improve by 18.7%/13.2%. Keep native serial and P2 batching
  defaults after worker/parameter controls, and retain upper-band limitations.
- Replace dense worker exponent-support reservations with an exact
  polynomial-norm bound. Optional 256-position chunks and a separately
  selected 1,024-atom cap resolve the tested historical medium batch refusals.
  Preserve assignment geometry across resume and release paused pool adapters.
  Final acceptance passes 259 tests and full lint; see the benchmark guide
  for uncertainty, rejected candidates, cold lifecycle and cooperative limits.

- Add an experimental coarse SIQS family executor for serial, thread and
  spawned PyPy process comparisons. Reverify bounded relation batches in a
  central store, preserve assignment identity across worker counts, charge
  parent-owned work leases and aggregate cooperative CPU, and drain cancelled
  work before returning. Add capped provenance checkpoints and charged replay
  of unfinished families. Fixed-work runs continue after direct residual splits.
- Retain serial collection after the P3.6 held-out experiment: worker arms do
  not pass the complete-factor promotion gate. Keep the standalone worker API
  experimental; medium-band batch-cap refusals remain explicit and resumable.

- Add an experimental bounded SSS/SSSf collector with seeded CRT/collision
  assignments, exact product/remainder-tree smoothness detection, complete
  exponent recovery and the shared QS relation/filter/extraction pipeline.
  Preserve candidate-prefix resume, finite storage and quiet calls.
- Add separate hash-checked upstream reproduction and trained held-out
  challenger runners. Optional comparison dependencies remain outside the
  factoring library; SSS has no automatic dispatcher promotion.
- Expose opt-in `--method sss` / `sssf` and `PortfolioConfig(sss=...)` under
  shared budgets and recursive result validation. Add capped SSS checkpoints
  with charged store/assignment/matrix reconstruction; retain version-2/3
  portfolio resume through version 4. Existing automatic defaults stay intact.
- Validate complete SSS/SSSf splits on the declared 30-digit held-out set.
  Keep larger capped runs diagnostic and retain the challenger pending a
  comparison with feasible trained SIQS; see the benchmark guide for limits.

- Preserve the original Python 2 implementation in v1; develop and test v2
  on PyPy implementing Python 3.11.
- Repair arithmetic boundaries, primality certainty, factor reconstruction,
  segmented sieve endpoints and saturated-batch recovery.
- Add bounded preprocessing and seeded rho/p−1/ECM portfolio execution with
  serializable checkpoints, shared allowances and explicit unfinished results.
- Add exact QS/MPQS polynomial and relation identities, bounded single-large-
  prime collection, provenance-aware filtering, GF(2) dependencies and modular
  factor extraction with in-memory resume.
- Remove redundant collector budget polling and proved-useless score marking.
  Keep experimental resieving and alternative scoring/filtering controls
  optional; see measurement scope.
- Keep code, tests, roadmap, research and required fixtures on GitHub.
  Exclude generated run captures, logs, profiles and local development records.

- Attempt factor extraction from retained QS relations before returning a
  storage-limit outcome; preserve resumable pending work and idempotent stops.
- Account for simultaneous QS collector, relation-preparation and matrix
  workspace before allocating private preparation storage.

- Add bounded experimental SIQS CRT/Gray families and incremental roots,
  with exact cached-root certification and compact family-only checkpoints.
- Complete bounded SIQS with shared verified relations, compact full-job
  checkpoints, charged root/matrix/prefix reconstruction and finite recovery.
- Add integer multiplier scoring, QS/MPQS/simple-family controls and optional
  ECM-to-SIQS dispatch under one allowance. Preserve version-2 portfolio resume.
- Validate disjoint certified inputs, including declared 30–80-digit capped
  comparisons. Keep ECM as the default: no broad SIQS crossover is established.

- Complete P3.4's bounded large-input evaluation with an independent certified
  30–80-digit and varied-shape corpus, trained frozen settings, cumulative
  checkpoint continuations, and declared negative outcomes.
- Add optional exact prime-power score marking, immutable base-column lookup,
  sparse matched-relation accounting and larger explicitly bounded capacities.
- Replace repeated global matrix rescans with exact deterministic incidence
  updates; the inspected 30-digit filtering control improves 5.920 → 3.425 s.
- Save wide binary dependency masks with compatible streamed hexadecimal
  digests and retain legacy checkpoint prefixes. A fresh balanced 50-digit
  SIQS run completes in 1,183.355 s; continuation time is recorded separately.

Historical P3.4/M31 decision (3 October 2026):
[P3.4](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p34--add-siqs-self-initialization-and-bounded-dispatch)'s
bounded implementation/evaluation is complete; ECM remains the default and
SIQS stays opt-in in that decision's scope. A general SIQS crossover and
practical 60–80-digit scaling remain unestablished. The earlier M30 completion
claim was corrected because 0.2-second probes and 24–26-bit fixtures were
insufficient; M31 supplies the declared longer experiment and its explicit
retain-baseline decision. See measurements.
At that date, P4.3 arithmetic assessment and bounded GNFS remained planned
work. PRAC optimization remained
[P4.1](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p41--finish-prac-repair-and-precompute-valid-chains);
production then used the safe Montgomery ladder. See [current status](v2/STATUS.md)
for later backend/PRAC decisions and scoped automatic SIQS routing.

## P3.8 R1 isolated capacity implementation — 4 October 2026

- Add opt-in streamed SIQS nearest/flyer A assignments with bounded resident
  state, up to 32 A factors, and separate cumulative search and Gray quotas.
- Add exact external-square MPQS coefficients and preserve square corrections
  through atomic/combined relations, extraction and checked checkpoints.
- Add explicit monotone direct-job resume extension without resetting checked
  stores or consumed resources; distinguish exhausted assignment space.
- Add actual-base capacity diagnostics and sparse initial-incidence work
  charging while retaining conservative matrix storage bounds.
- Add independently certified sub-100-digit training/confirmation generation,
  frozen-control checks and matched full-call/cold/profile runners. Experimental
  defaults and overall R1 completion remain gated on confirmation/integration.
- Frozen R1 confirmation completes the fresh 30-digit control in all four
  arms (18/18 attempts each), while all balanced upper bands remain censored
  under one-second caps. Nearest/flyer/reference controls, a 3.0% legacy
  regression, cold lifecycles and separate profiles are recorded without
  changing defaults or claiming a broad crossover.
- Preserve the exact measured source snapshot and a bounded restoration tool;
  account for verification-cache growth in explicit quota extensions.
- Integrate R1 with the P2/P3.6.1 repairs under a single writer. The combined
  candidate passes 258 tests and full lint; broad upper-band feasibility and
  crossover work remain open, with no automatic parameter promotion.
- Final R1 integration acceptance passes 259 tests, lint and all 40 benchmark
  imports from a clean committed candidate. The 58-input varied-shape screen
  validates 464 outcomes; a warmed final-runtime bridge completes 18/18.
  Structured controls and upper-band time limits remain separately labelled.


## P3.8 R3 isolated preparation and provenance work — 4 October 2026

- Preserve full and matched relations in one stable admission sequence; bind
  prepared identities to complete immutable base, polynomial and atom payloads.
- Retain mixed row order in SIQS/SSS version-2 and parallel version-3
  checkpoints. Read preceding supported formats with their original grouped
  order and fully reverify every relation before charged solver replay.
- Count pivot nonzeros incrementally. Preserve exact original-row dependency
  lifting and conservative fill/provenance workspace reservations.
- Add opt-in `filter_row_growth` and `tested_dependencies` SIQS/QS settings,
  retaining every-change filtering and disabled tested-dependency caching as
  defaults. Pending work resumes immediately; terminal windows/storage limits
  force the remaining extraction. Skipped checked trials still count against
  the existing trivial-dependency allowance.
- Add independent certified fixtures, immutable controls and bounded matrix,
  compaction, history, cache and complete-factor experiments. Decisions and
  confirmation scope are recorded in the benchmark guide; matrix scaling and
  conditional packed-exponent/root work remain open.

- Confirm the scoped opt-in cadence-32 setting on fresh 30-digit fixtures:
  38.7% lower complete-cohort median than R3 cadence 1, 54/54 complete per arm.
  Keep cadence 8 and tested-dependency cache promotions deferred; neither this
  experiment nor passing recovery tests changes automatic dispatch defaults.

- Integrate only the R3 delta with the repair/R1 work. The combined candidate
  passes 273 PyPy tests, lint and all 44 benchmark imports from committed inputs;
  a stable matched bridge validates 36 complete outcomes in every arm.

## P3.6.1 follow-up: restrictive allowances and repeated setup — 4 October 2026

- Treat worker assignment work as a ceiling, dispatch smaller remaining
  allowances, and wait for running lease refunds before refusing paused
  relation admission or matrix continuation. Keep exact charges and replay.
- Use local coordination state for serial/thread pools and reuse their
  parent's checked immutable base. Process ingress and returned relations
  retain full validation; aggregate limits and draining remain enforced.
- Cache one family identity per immutable base, with unchanged digest bytes,
  family order and checkpoint identities.
- Exclude disabled p−1/ECM stages from schedule and baby-table reservations,
  count a worker's shared family/collector base once, and reserve resieving
  exponent support from an exact norm bound before allocating block scratch.
- Preserve frozen post-R3 controls, measured prototypes, selected sources and
  a fresh independently certified confirmation corpus. Reject extra Hensel
  caches, ECM baby-step batches, single-hit valuation and CPU-read batching
  after their scoped comparisons missed the promotion gate.
- Confirm 12.3–12.7% lower complete serial-worker time and 32.2% lower time
  for the declared B=10000 chunked case on fresh inputs. Restrictive-capacity
  classes recover useful completion without raising their allowances; keep
  kernel/setup, fixed-work, cold-pool and complete-factor claims separate.
- Pass 282 PyPy tests, lint, exact cross-worker fixed-work comparisons,
  cooperative limit checks and retained-only test/benchmark imports.

## P3.8 R2 isolated collector evaluation — 4 October 2026

- Add opt-in integer fixed-point prime-power scoring with conservative rounding
  and saturating array scores, retaining exact relation verification.
- Add keyword-only `power_plan_bytes` with an upfront finite reservation,
  capped per-polynomial plan reuse, streamed fallback and charged checkpoint
  rebuilding. Preserve existing positional configuration arguments and
  first-uncommitted-position resume behavior.
- Evaluate sparse recovery, scoring, bounded batch/tiny-prime/grouped-charge
  challengers and whole-interval CRT eligibility against an immutable
  repair/R1/R3 control with independent certified training and held-out inputs.
  All complete-factor confirmation attempts finish; complete collector arms
  preserve atomic signatures and post-filter results.
- Retain defaults: fixed scoring/plans remain experimental, and no candidate
  passes the declared causal training promotion gate. Family-wide CRT and
  wide-block resieve capacity remain deferred. See the
  R2 measurements
  for commands, uncertainty, censored outcomes and cold/profile separation.
- Verify the isolated committed-files-only snapshot with 282 PyPy tests,
  full lint, all 45 benchmark imports and frozen-source/corpus checks.
  This records the original isolated acceptance; combined acceptance follows below.

## P3.8 R2 mainline integration — 4 October 2026

- Integrate R2 fixed-point scores and capped power plans as explicit opt-ins,
  preserving the latest repair/R3 contracts, the resieve support bound and
  existing defaults. Keep family-wide CRT and broader capacity changes
  conditional on new workload evidence.
- Move required R2 controls/corpora into versioned input folders without
  changing historical bytes; retain a separate repaired-mainline integration
  baseline and source/configuration freeze. Keep generated captures local.
- Validate the combined source with 292 PyPy tests, lint, all 48 benchmark
  imports and a committed-files-only checkout. Quiet matched comparisons
  complete every attempt; 1,188 collector outcomes retain exact atomic and
  post-filter results, including wide-block resieving under the same allowance.
  Previously inspected integration inputs supply no fresh promotion claim.
