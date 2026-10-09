# v2 changelog

## 9 October 2026 — QS snapshot ownership test

- Harden the snapshot-release regression test with three bounded GC passes,
  following the installed Python 3.11 test-support cleanup convention.
  [PyPy documents](https://doc.pypy.org/cpython_differences.html) cases where
  several collections are needed before weak references clear.
- Add a deliberately retained-snapshot control that must still raise, and
  verify reconstruction of the normal unresolved prime fixture. Keep the
  assertion at the next collection boundary.
- Preserve production QS/SIQS code, arithmetic, budgets, checkpoints and
  defaults. The original intermittent A10 full-suite failure is retained
  locally; its exact transient cause remains unconfirmed. One hundred
  unchanged-case runs across default, disabled and accelerated JIT settings
  and forty observed jobs did not reproduce it. The three-pass variant passed
  forty accelerated runs and rejected all five in-memory retention controls.
  This is test hardening, not evidence of a repaired factoring leak or a
  performance improvement.
- Four further original A10 GMP suites (three default JIT and one accelerated)
  and a corrected failure-observer suite pass all 372 tests without recurrence.
  The repair passes all 358 PyPy/GMP tests, both 20-test pipeline suites with
  accelerated/disabled JIT, and full lint. The earlier observer-driver timeout
  is excluded from acceptance and retained locally.

Run `make -C v2 test`, `make -C v2 lint`, or the focused regression with
`pypy3 -m unittest v2.tests.test_qs_pipeline.PipelineTests.test_batch_snapshots_released_before_more_collection`.
