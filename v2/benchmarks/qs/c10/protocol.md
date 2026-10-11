# C10 tranche 1: premeasurement protocol

Hypotheses: (Q) ordering bounded batches by actual A-product error improves
useful yield enough to repay selection cost; (P) the existing 2 MiB capped
prime-power plan reduces repeated block lifting on the settled C9 collector.
Compare `control`, `quality`, `plan` separately. No combined arm, threshold
loss, new polynomial representation, base change, or residual-policy change.
Defaults remain unchanged. Quality's batch size 16 is frozen, not grid-searched.

The exact runtime, runners, protocol, source manifest, C9 decision and full
configuration dictionaries must be frozen in `c10_training.json` before any
measurement. C9 runtime is based on `129cb09e4ca684f3eb08cbbdbea3fb43848e5d95`;
freeze only after its final decision and CPU handoff. Record any later source
correction as an explicit new identity; never silently overwrite evidence.
30/40 controls use C1's accepted DLP bundles with C9's explicit growing storage;
50 uses C9's settled aligned full-base growing configuration. Growing storage
is user-directed implementation, not assumed performance acceptance.

Training has three previously revealed inputs per band. For 30 digits use
C1's two fresh inputs and C9 growth's one; 40 uses the three C9 growth inputs;
50 uses the two original C9 training inputs and the first C9 growth input.
All are training for C10. Seeds 7,29,47 cross each input: nine matched samples
per arm/band. Each new measurement process warms every active arm with at
least three seconds of validated wall AND CPU execution. Quiet PyPy 3.11 only.
Randomized deterministic arm order is fixed by band and group index.

Caps per call are 10/60/600 seconds for 30/40/50-digit training and
10/60/1,200 seconds for confirmation. These are observations, not feasibility
ceilings: prior 50-digit C9 progress justifies longer matched observations.
Each call also has 10^13 cumulative work units, 1 GiB owned memory and a
3 GiB process RSS watchdog. No work/time resets inside a call. Complete
factor cost includes setup, collection, solving, extraction, classification
and independent proper-divisor/reconstruction/certainty validation. Censored
or otherwise unresolved calls cost the full assigned cap in the selection
metric; report their actual spent work and reason as well.

The entire tranche has 14,400 seconds each of aggregate wall and CPU charged
work (waiting for others is excluded), including cold starts, warmups,
generation, profiling, QA, failed and interrupted calls. Prepay 900 seconds
for heavy QA; if insufficient, stop before new experiments and document the
need for an explicit budget amendment. Training can consume at most 5,400
additional seconds. Diagnostics can consume at most 600. Generation <=120.
Other work shares the remaining total. Captures are bounded to 256 MiB;
checkpoints to 16 MiB per diagnostic and total captures. An interrupted call
retains its full reservation. Reserve an entire matched group before starting
it; if it cannot fit, preserve previous groups and report resource censoring.
The lock must be held throughout heavy QA and experiment windows. Never
interrupt another owner's work or spin until the lock becomes free.

At nine samples extend to eighteen only if paired directions disagree within
an input, repeated-cell wall drift exceeds 15%, or paired effects differ by
more than 15 percentage points between the first and last three groups.
Extension is within the same aggregate caps. Eighteen is the final tranche
sample ceiling; remaining instability is reported and prevents a speed claim.
The 15% values detect instability; there is no minimum improvement percentage.

Select at most one challenger per band, first by number of independently
validated completions, then by sum of capped complete-run cost, with ties
retaining control and then preferring plan to quality. Advance only if it
beats control on that ordering. A training selection is not adoption.
Freeze selection, raw-evidence hashes and source/configuration identity
BEFORE generating a new corpus. Generate three independent balanced inputs
per selected band, disjoint from all existing corpus n values, using the
repository's independently verified prime certificates. The population is
conditioned on the certificate generator and is not claimed RSA-distributed.
No parameter changes after fresh values are revealed.

Confirmation compares control with that single selected challenger, nine
samples (three inputs x three seeds), with the same extension rule. Report
per-input paired effects, complete counts, capped work and elapsed cost;
hierarchical 95% bootstrap intervals resample inputs then matched samples.
Three independent inputs do not become nine independent inputs through seeds.
Scoped performance adoption requires no loss of completion, stable complete
outputs and an interval wholly below zero for challenger/control cost change.
Mixed/censored/unstable evidence is deferred or rejected, never a universal
speed claim. If no arm advances, record the negative training result and stop
confirmation; do not spend fresh inputs merely to fish for a winner.

Cold startup is one separately labelled subprocess per active arm/band, never
a warm sample. Profiles are separate instrumented calls, not timing evidence.
Record setup/collection/filter/solve/extraction timing, owned workspace, RSS,
residual-stage counters/profiles and checkpoint encode/restore costs. Measure
A-product error and overlap on fixed bounded prefixes. Raw relations/cycles,
kernel timing and residual magnitude alone cannot establish complete benefit.
Useful rows/dependencies and proper-factor outcomes must accompany yield claims.

Required validation: exact normalized-polynomial and root oracles, family
population/prefix identity, charged cache rebuild and failure atomicity,
finite and growing checkpoint replay, `make -C v2 test`, `make -C v2 lint`,
then the same tests and all benchmark imports from committed files only.
Publish concise decisions and limitations; keep raw output ignored. Do not
merge or push. C10 continues after this bounded tranche.
