# C1 bounded DLP feasibility protocol — frozen before new probes

Control source: committed mainline `76f06b0110941017b926caa51d341db6221f470b`.
C1 works only in `codex/c1-dlp-siqs`. B3 owns production ECM; A7 owns optional
collector/checkpoint reconciliation. No shared API/checkpoint/portfolio changes
are authorized by a feasibility pass alone without reconciling those owners.

## Question and controls

Can rejected residuals produce useful extra constraints at an affordable total
cost? B1's 60-digit run retained 4,081 unmatched partials, one match and 47 rows
which all singleton-filtered away. This is not evidence of DLP economics.

Use the existing certified `p38_r1_training_corpus.json` balanced 30-, 40- and
60-digit fixtures, seeds 7 and 29. The first two bands use B1's selected SIQS
bundles from `b1_selected.json` and `b1_40d_selected.json`. The 60-digit control
is B1 `upper_config(60, 'siqs')`; it is a diagnostic configuration, not a
calibrated completion preset. All inspected fixtures remain training inputs.
Use python-int on PyPy implementing Python 3.11; refuse another interpreter.
No external/native timing predicts PyPy performance.

## Finite allowances

Six control starts, at most 5/30/30 wall and CPU seconds per 30/40/60-digit
start, 10^13 production work units, 256 MiB owned production workspace.
Six instrumented starts: at most 512 blocks (4,096 positions each) initially,
90 wall/CPU seconds and 10^13 production work units each. Record the 128-block
prefix as a duration comparison. Sampling/classification/splitting has a
separate 10^11-unit, 60-wall/CPU-second budget per start, included in the
instrumented process cost. A diagnostic refusal is censoring, never success.
No performance promotion is made from these instrumented runs.

For each base bound B, retain SLP's original prime bound B^2. DLP separately
limits *each* prime to 100*B and the product to (100*B)^2; do not conflate
these with the SLP prime cap. Allow at most 8,192 composite splitting attempts
per start, two deterministic Brent-rho attempts of at most 2,048 evaluations
each, batch 32 and recovery 64. Reserve their full work allowance before
splitting. Each recovered endpoint must independently classify PROVEN; an
integer square gets an exact-root check and one proven-prime check. Product,
proper split, prime bounds and original exponent identity are checked exactly.
No probable-prime residual, coverage shortcut, or opaque cofactor is admitted.

Retain at most 16,384 exact diagnostic records and 64 MiB estimated record
storage per start; at most 128 uniformly reservoir-sampled records in each
rejection stratum and 16 uniformly chosen positions per block, including
threshold rejects. Persist at most 16 MiB per cell and 256 MiB of study
captures. Stop on a cap; record unexamined/censored populations. Diagnostic
matrix workspace is capped at 128 MiB and 10^10 work units / 20 wall and CPU
seconds. These allowances are not production graph/storage guarantees.

## Sampling and measurements

The instrumented SLP run uses its unchanged score thresholds and admissions.
Audit every position passing the conservative DLP product threshold, including
SLP block/refinement rejects. Recover exact valuations with the existing root
hit map; independently full-divide the random 16 positions per block and
check agreement. Maintain separate counts for block threshold, refinement,
product cap, prime-too-large, split failure, non-semiprime and endpoint bounds.
Report residual bit lengths/shapes, cost of classification and splitting,
retained records, SLP unmatched occupancy/matches/evictions, post-filter rows,
independent dependencies and both-sign proper-divisor yield.

Analyze retained rows offline as full factor-base/large-prime GF(2) incidence;
this is a complete diagnostic matrix, not a graph collector. Compare full/SLP
rows with SLP+DLP rows from the same positions. Include all components and
self-loops. Validate every dependency against original incidence and every
square congruence using exact exponent sums, both GCD signs and reconstruction.
Do not call raw edge counts a factoring gain. Report ideal retained incidence
separately from actual bounded FIFO SLP matching. No-retention-loss diagnostics
are optimistic, not a prediction of a bounded production graph.

## Decision and stop rules

A cell may proceed to the **one** duration extension (4,096 total blocks,
180 wall/CPU seconds and 120 diagnostic wall/CPU seconds, same work/storage/
split caps) only if the first screen contains at least 32 certified DLP edges,
at least 16 repeated endpoint occurrences, no diagnostic cap, and splitting
plus certification CPU is below the complete SLP control CPU. This screen
permits measuring delayed graph formation; it does not establish feasibility.
At most six such extensions; no extra bands, bounds, seeds or splitter tuning.
Total experiment process allowance: 2,000 wall and CPU seconds. Stop the
campaign at that limit and defer remaining cells explicitly.

A **go** requires, for both seeds in at least one declared band: no correctness
or diagnostic censoring failure; at least max(16, 10% of SLP rows) additional
independent LP-cancelled constraints; a nonempty post-filter core or independently
verified dependency; and splitting/certification diagnostic CPU no more than
25% of the uninstrumented complete SLP control CPU. This is a deliberately
conservative investment screen for an invasive representation/resume change,
not a speed claim. At 60 digits with censored controls require an actual
verified useful dependency; extrapolated relation counts cannot pass.
If any requirement remains unsupported, finish with a bounded defer decision,
rejection economics and a specific future trigger. Do not implement DLP.

Only after a go: freeze a separate implementation/training/confirmation
protocol before further measurements. Preserve R3 identities, exact lifting,
independent complete-cycle oracle, residual certainty, charged resume and all
finite stores. Train changed settings on training inputs only; generate fresh
confirmation after selection. Matched total costs must include splitting,
certification, graph, filtering, provenance, matrix, extraction and final
classification. Use >=3 seconds validated warmup and >=9 samples, prespecified
stability/uncertainty rules and the revised roadmap promotion policy. A
feasibility pass does not close P5.4 implementation or R4 promotion gates.

Hold `/private/tmp/factor-performance.lock` for experiments and heavy checks.
Keep captures, profiles and transcripts under ignored `results/c1/`. Required
protocols/corpora/runners are versioned; preserve existing evidence. Do not
merge or push. Run tests/lint and committed-files-only tests/imports/loaders.
