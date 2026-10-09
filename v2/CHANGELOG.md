# v2 changelog

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
