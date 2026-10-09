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
