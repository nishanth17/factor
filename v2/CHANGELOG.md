# v2 changelog

## A6 / P5.3 bounded p−1 correctness tranche — 9 October 2026

- Add the exact inclusive LCM ratio, including old-prime power increases,
  and an explicit finite one-base campaign. Reuse complete verified stage-one
  state; rebuild stage two after increased B1 and append only new B2 coverage
  when B1 is fixed. Saturation and nonunits stop this base after finite recovery.
- Introduce isolated schema-1 Python-integer checkpoints with complete
  assignment/configuration identity, charged deterministic reconstruction,
  canonical numeric types and cumulative work/wall/CPU allowances. No RNG is
  consumed. Preserve proper divisors, reconstruction and unresolved certainty.
  Existing direct, portfolio and ECM APIs/checkpoints remain unchanged.
- Inspect primary papers and pinned GMP-ECM, CADO-NFS, FLINT/PARI and SymPy
  implementations/licenses. Implement the integer identities independently;
  keep native polynomial/REDC techniques and ECM/Lucas migration separate.

- Reassess the original 10% gate under the revised roadmap without retuning:
  all six 27-repeat complete-stage comparisons confirm chunk 64 gains of
  6.1–9.2%, with positive paired intervals and declared stability. Accept this
  bounded configuration for scoped integration review; keep shared defaults,
  bounds, attempts, caches and allocation unchanged. Exact in-memory ratios
  save 8.8–9.6% versus fresh-each; charged checkpoint reconstruction loses
  12.4–14.2%. The declared portfolio cohort finds no extra completions.
- Pass 21 new independent/adversarial checks, full `make -C v2 test` with
  419 tests (three optional-GMP skips), and `make -C v2 lint`. A committed-file
  archive passes 422 PyPy/GMP tests, all 70 benchmark imports and independent
  certificate/product verification. Python sources and required inputs remain
  identical after the documentation-only acceptance update. Keep raw evidence
  in ignored results; this isolated tranche remains unmerged.

## A10 integration and research reconciliation — 9 October 2026

- Integrate source-verified strict 12/13-base Miller–Rabin ranges with settled
  B1/B2 mainline. Preserve B2 schedule/schema identities and legacy primality
  policy, random witnesses, cumulative work and conservative terminal labels;
  add int/GMP regression coverage for paired/wheel checkpoints 7/8.
- Cross off the completed bounded A10/P8.2.1 row. Publish an 11-project pinned
  QS/GNFS comparison and route prospective experiments to existing owners;
  no unmeasured candidate or combined portfolio/default gate is promoted.
- Pass 401 PyPy/GMP tests and lint; the committed-files-only archive passes
  398 system-PyPy tests (three optional-GMP skips), all 67 benchmark imports
  and required A10 proof/control/protocol/v1-adapter loaders. Preserve versioned
  inputs/citations/summaries; remove the isolated worktree and generated local
  captures after integration at the user's request.

## B2 aligned-wheel follow-up — 9 October 2026

- Add opt-in `ecm_pair_wheel`, with complete bounded nearest-center cells,
  coprime baby-point storage and shared certificate/saturation recovery.
  Preserve original paired execution, old checkpoint encodings and defaults.
- Add version-8 canonical wheel checkpoints and exact int/GMP continuation,
  including initialization at zero-predecessor giants and one-sided tails.
- Freeze fresh independent inputs, the original B2 control at `a670c4d` and
  three wheel choices per tier against accepted streamed/reusable controls.
  Complete 35 training and 20 held-out comparisons plus 180 separate cold
  starts. Extend unstable arms to 31/63 samples; every final capture meets
  the frozen stability criterion.
- Retain defaults. Selected W=30/210/840/1890 cuts campaign products 16.4%
  at exact coverage and reduces retained baby points from 1,024 to 216.
  Medium/campaign time improves 31.1%/4.8% against original pairing, but the
  uneven cohort regresses 44.0%. All selected wheels lose to reusable programs;
  their campaign cost is 23.6% higher. Completion is unchanged. Unpaired
  programs still save 10.8% versus streamed on this fresh nonsplitting campaign.
- Pass full test/lint commands, committed-only suites with 376 system-PyPy
  tests and 379 PyPy/GMP tests, and all 61 benchmark imports. Keep raw timing,
  profiles and work/storage evidence local; version required inputs and
  document both gains and losses with fixed-cohort limitations.
- Defer extended distance sets, relocation, graph matching and common-Z to
  C2; polynomial continuation to F3. Increased B1 still depends on A6.

## B2 paired ECM continuation — bounded opt-in tranche, 9 October 2026

- Add opt-in bounded reusable +/- stage-two execution through
  `PortfolioConfig.ecm_pair_distance`. Retain Python integers, the Montgomery
  ladder, unpaired streamed execution and existing production defaults.
- Keep immutable coverage records in the accepted finite program store and
  baby/giant points, term products and recovery in each curve's private job.
  Charge construction, decoding, arithmetic and recovery; preserve exact
  eligible primes, segmented tails and direct exceptions.
- Recover saturated paired terms by replaying both certified prime scalars.
  Validate every returned divisor. Add version-7 paired checkpoints, canonical
  int/GMP resume and predeclared finite campaign continuation with cumulative
  allowances. Preserve existing unpaired schemas and accounting.
- Freeze committed `94caf40` streamed/reusable controls, independent certified
  training/held-out inputs, D choices, budgets and selection criteria. Complete
  exclusive-window PyPy comparisons: 30 training arms, 15 held-out arms and
  135 separate cold starts. Extend four unstable training arms to 31 samples.
- Retain defaults after negative matched evidence. Selected paired execution
  is 81.3%, 154.2% and 45.0% slower than streamed complete factoring on the
  held-out small, medium and uneven cohorts; the nonsplitting campaign is
  15.0% slower. Completion is unchanged. Reusable unpaired programs also win.
  D choices 24/64/768/2048 are experimental selections, not recommendations.
  Segment boundaries eliminate pairing at the largest campaign choices;
  smaller D's product savings do not offset total overhead.
- Pass 367 system-PyPy tests (two optional skips), 370 PyPy/GMP tests with no
  skips and full `make -C v2 lint`, in the coordinated correctness window.
  Repeat both suites and import all 59 benchmark modules from committed files.
- Defer increased-B1 continuation because A6's exact schedule-ratio contract
  is absent. Adding curves or changing bounds on an exhausted checkpoint,
  wheel/common-Z, production PRAC, kernels and allocation remain separate.

## 9 October 2026 — QS snapshot ownership test

- Isolate the snapshot-release ownership assertion in a finite PyPy child
  with JIT disabled. Preserve JIT settings in the parent suite and normal
  production execution. Keep the deliberate retained-snapshot control and
  reconstruction of the unresolved prime fixture.
- Reproduce the original one-GC assertion on job 12 in three fresh processes
  using `--jit trace_eagerness=1`. A local heap capture identifies a compiled
  `JITFRAME` and the active tracer's `History → RefFrontendOp` as snapshot
  owners after Python's `del result`; three extra collections cannot release
  these live roots. No Python frame local owns the snapshot. This corrects
  the earlier assumption that more collections would suffice.
- With JIT disabled, all 120 original jobs pass; an in-memory control omitting
  `del result` fails immediately. The historical A10 failure has no heap
  capture, so its exact bridge/guard identity cannot be recovered. Raw
  reproduction logs and the heap remain in ignored local results.
- Production QS/SIQS code, arithmetic, budgets, checkpoints and defaults are
  unchanged. This fixes a test's ownership inference; it makes no factoring
  performance claim or roadmap promotion.

Run `make -C v2 test`, `make -C v2 lint`, or the focused regression with
`pypy3 -m unittest v2.tests.test_qs_pipeline.PipelineTests.test_batch_snapshots_released_before_more_collection`.

## 9 October 2026 — bounded B1 / P3.8-R1 calibration

- Add frozen joint QS/MPQS/SIQS calibration runners, controls, independent
  certified confirmation corpora and evidence-validation tests on integrated
  R2 mainline `94caf40`. Keep raw captures in ignored local results.
- Retain explicit balanced presets: 30-digit SIQS uses base 3,000, half-width
  8,192, four flyer-selected A factors and eight effective Gray polynomials;
  40-digit SIQS uses base 10,000, half-width 65,536, five A factors and sixteen.
  Fresh 30-digit balanced time falls 36.6%; 40-digit external-square MPQS's
  6.9% reduction against feasible SIQS misses the promotion threshold.
- Preserve runtime defaults, arithmetic, certainty labels and finite resume
  behavior. Uneven/structured regressions, failed wider-QS cohorts and censored
  60–99-digit probes limit promotion. Defer larger calibration, combined-source
  acceptance and portfolio handoff; no DLP, CRT, matrix or GNFS expansion.
- Verify 365 PyPy/GMP tests and selected early/deeper checkpoints. See the
  [benchmark record](benchmarks/README.md#b1-joint-qsmpqssiqs-calibration--9-october-2026)
  for configurations, timing uncertainty, commands and remaining gates.

This scoped changelog is under `v2/` to honor the B2 work boundary; the root
changelog and `v1/` are unchanged.
