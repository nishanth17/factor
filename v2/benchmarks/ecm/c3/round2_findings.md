# C3 continuing search: attribution findings

These are diagnostic findings, not a policy promotion. The 84-call pilot
frozen at `1864804` used six inputs, seeds 17/43 and seven bundles; every
call reconstructed exactly and completed. Cold/instrumented wall and CPU
observations cannot establish a speedup. All revealed subjects remain training.

## Coverage and economic witnesses

On `r2_40_balanced_0`, seed 43, the 32-curve control reached SIQS; expanding
its identical-bound prefix found a divisor on curve 60 and avoided SIQS.
On `r2_40_uneven_16_0`, seed 17, that happened on curve 64. Conversely,
seed 17 of the balanced input exhausted all 64 curves before SIQS. These
witnesses show both additional useful coverage and its cost on failures.
They justify measuring conditional yield; two hits do not establish a prior
or a general dispatch rule. Compact B2 and independent escalation remain
candidates. Economic25/50 are fixed cumulative clock ceilings, not fitted
probability sequences.

## Reservation accounting candidate

Matched cProfile runs used the first new balanced 30-digit subject, seeds
17/43/89, and at least three validated CPU seconds of warmup per arm. Both
profiles performed 48,654 portfolio advances and 48,351 arithmetic advances.
The protected path originally sampled each clock 402,614 times, versus 52,262
for the historical control. It constructed a policy wrapper repeatedly and
performed two base checks plus separate policy clock reads for each charge.

The candidate reuses one view per invocation and samples clocks once per
reserved charge on the shared base ledger. The matched protected profile
then sampled each clock 100,775 times. Remaining outer-loop checks are kept.
These counts explain the mechanism; cProfile times are not performance evidence.
The ordinary Budget.consume path, numerical schedule, work charges, seeds,
configuration fields and serialized formats remain unchanged. Cancellation
and deadlines remain cooperative; one reservation charge checks cancellation
once and uses one clock observation instead of several. A deadline reached
between former duplicate observations can therefore be seen on the next
atomic charge. No exact callback invocation count is promised.

The full pre-optimization source at `1864804` is an immutable versioned
baseline. A differential test compares refusals, reasons and unchanged charges
at fixed cancellation, actual-limit, policy-ceiling and protected-floor
boundaries. All 600 PyPy tests and full lint pass. The first full QA exposed
missing historical-path budget pins; an additive A7 manifest repairs those
pins while preserving every previous manifest and the frozen A7 arms.
Matched warmed complete calls remain required before a performance claim.

## Reproduction and local evidence

The pilot can be reproduced from commit `1864804`; its frozen hashes correctly
reject subsequent production edits. Earlier C3 evidence similarly belongs to
its recorded source commits. Do not rewrite historical controls to make a
new checkout pretend to be an old experiment.

Ignored local evidence is under `v2/benchmarks/results/c3/round2/`: `retry-v2/`
contains all pilot rows and its receipt; `profiles-before/` and `profiles-after/`
contain profiles and validated calls; `optimization-qa/` preserves the first
QA failure and `optimization-qa-retry/` the passing tests and lint. Including
the conservatively charged first attempt, C3 used about 228 active seconds
of the coordinated 600-second diagnostic lease and released it to C9.

Next, fit a small offline stopping sequence using measured conditional yield
and cost, then compare its compiled finite tiers against simple fixed
schedules on warmed complete calls. The model must earn its complexity before
any production dispatch change or roadmap completion.
