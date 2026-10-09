# v2 changelog

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
