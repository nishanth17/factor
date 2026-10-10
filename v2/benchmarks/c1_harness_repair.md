# C1 harness repair, before the unfinished 40/60-digit cells

The first campaign used frozen driver `1db45cf`. Both 30-digit cells were
saved and are retained. The next uninstrumented 40-digit control completed,
but B1's historical validator raised before serialization: it hard-codes
PROBABLE above 2^64. Mainline includes A10's accepted wider deterministic
range, so the result has PROVEN labels. The traceback remains in
`results/c1/screen-run.stdout`. No 40/60-digit residual probe had run.

C1 now validates reconstruction, proper divisors, multiplicities against the
certified corpus, and the strict accepted A10 endpoint locally. Historical
B1 code, inputs and captures are unchanged. A targeted regression checks both
sides of the certainty boundary. No collection/splitting/bound/selection rule
changes. A single continuation reuses the saved 30-digit cells and reruns the
failed-to-serialize 40-digit control once, then runs the four unfinished probes.
Thus at most seven control attempts, six diagnostic cells, and still at most
2,000 aggregate process seconds including the failed attempt (reserve 60
seconds for the initial interrupted campaign). No failed gate is rerun.
The continuation records both driver hashes and refuses changed runtime hashes.

Diagnostic `positions` counts scheduled block widths; a cap may interrupt the
last audited block. Population ratios must use the actual stratum counters,
not treat the unvisited tail (at most 4,096 positions) as sampled. Uniform
samples in that final censored block are a prefix of its random selection.
Matrix reservation refusals are recorded, not bypassed by lowering constants.
Further offline inspection of already-retained prefixes is diagnostic only and
cannot reverse a censored no-go or grant a collection extension.
