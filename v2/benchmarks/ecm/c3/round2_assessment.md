# C3 round-two complete-call assessment v1

Frozen before fitting the sequence or collecting the three-arm comparison.
This assessor applies the gates in round2_training_protocol.md without ranking
instrumented calibration. It does not revise the calibration assignment,
weights, ladders, endpoints, minimum risk-set population or finite allowances.

Require all 57 subjects, all nine assigned seeds and all three arms. Reject
missing/duplicate rows, instrumented calls, seed mismatches and a changed
training source identity. Calls already validate proper divisors and exact
reconstruction before capture. An unresolved return pays at least its full
5/30-second grant in the time objective; returning early is not a speedup.

Weight modulus bands equally, kinds equally within each band, subjects equally
within kind and seeds equally within subject. Bootstrap 10,000 times using
Random(2026101109), sampling subjects within strata and matched seeds within
each sampled subject. Replay each bootstrap assignment for both arms. Report
CPU improvement intervals against historical control and optimized fixed32,
wall ratios and completion separately for every band/kind.

Only a positive lower 95% CPU-improvement endpoint against both controls,
wall ratios no greater than 1.05 and no class completion loss above five
percentage points can advance the fitted table to new confirmation. Otherwise
retain the fixed control and record a negative or inconclusive model result.
No replacement is selected from comparison results.

These nine-seed training comparisons are selection evidence. They do not
establish stable accepted performance. A candidate that passes needs a new
confirmation freeze before generating disjoint subjects, with new seeds,
at least nine matched observations and predefined finite instability
extensions. Report dispersion and any still-unstable cells there. No production
default or roadmap-completion change follows from this assessor alone.

Run after the committed fitted selection and complete comparison:

    pypy3 -B -m v2.benchmarks.ecm.c3.round2_report CAPTURES OUTPUT.json

The additive c3_round2_assessment_v1.json pins this assessor, its substantive
hand-constructed tests, this procedure and the immutable training control.
Generated assessment, rows and receipts remain in the ignored results tree.
