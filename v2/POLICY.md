# Factor policy

Start with [status](STATUS.md), then the [live backlog](ROADMAP.md).
Stable rules below retain the roadmap's promotion policy revised 9 October
2026 and its acceptance contracts. Historical protocols and verdicts are
preserved in the [history index](history/README.md); they are not rewritten.

## Scope and standing directions

Work in `v2/`; preserve `v1/` and historical captures as provenance. Use PyPy
implementing Python 3.11 only, with exact integer arithmetic, integer roots
and checked inversion. Optional native kernels/backends are allowed; retain
a working pure-Python implementation on PyPy. Python integers, gmpy2 and fused
native kernels remain distinct measured configurations. Keep canonical results
and backend-bound checked checkpoints. An explicit unavailable backend request
is an error, never a silent backend substitution or CPython fallback.

The objective is minimum complete-factorization wall time on the Apple M4 /
24 GiB / 10-core machine, with all cores available. Select useful concurrency
from measurements; using every core is permitted, not a requirement when
coordination, heterogeneous cores or memory pressure makes fewer faster.
Aggregate CPU time, peak memory and startup/cancellation waste are secondary
costs and finite resource constraints. Increased CPU consumption does not by
itself invalidate a confirmed wall-time gain.

The primary workload is uniformly sampled integers in each exact decimal-size
band: 60, 80 and 100 digits, with equal weight (one third each) across bands.
Do not condition this population on compositeness, parity, trial survival or
known factor size. Include prime and preprocessing exits with the same result
certainty contract in every arm. Freeze corpora and seeds before measurement;
publish per-band median and p90 time to complete factorization and completion
rates. Report equal-weight summaries without hiding a band behind easy inputs
from another. Repeated seeds do not supply independent moduli or reliable
tail estimates by themselves; size the fresh input cohort for p90 uncertainty.

Balanced-semiprime completion at 60 and 70 digits is a must-not-regress
constraint under matched declared grants. No completion-loss allowance from
the general promotion policy overrides this constraint. Prespecify its paired
completion/uncertainty decision and report balanced latency as well. Missing
or inconclusive guardrail evidence cannot be described as a passed guardrail;
retain the baseline or keep an improvement scoped and opt-in until resolved.
This governs adoption, not admission: a timed-out control does not bar a
correct, feasible candidate comparison under declared finite grants.
Near-term useful progress through 60–80 digits remains central; targeting
100 digits does not establish current feasibility or authorize a large run.

Timeouts and unresolved cofactors are incomplete outcomes, not fast completed
calls. Include them in completion denominators and retain their charges.
Use a prespecified censoring-aware quantile analysis; if a quantile is not
identified within the grant, report its bound or that it is unestablished,
not the median/p90 of successful calls alone. Report trial-sifted survivors,
uneven smaller-factor strata and recursive composites separately for diagnosis;
their oversampling must not change the primary random-integer weights.
The roadmap ranks expected benefit per effort under these choices, not measured
speedups. Existing historical verdicts and active freezes are unchanged.

For source-derived work, distinguish implemented controls, active worktree
coverage, retained losses and uncovered deltas. Trace the caller and build
configuration before treating an upstream option as production behavior.
Upstream overflow probabilities, digit cutoffs, torsion or success claims do
not supersede exact local invariants or demonstrated counterexamples. Preserve
read-depth/access limitations in the local source ledger; a downloaded tree
is not proof of an exhaustive implementation audit.

B4a/B5/B5a may select NumPy, fused native arithmetic or C/CFFI bulk kernels.
Choose one justified implementation with an exact Python fallback; do not
require parallel NumPy and CFFI projects. Label accelerated results,
dependencies and resources explicitly. Native code is authorized as an optional
backend, with production/default promotion subject to the same gates; a
feasibility result does not adopt it or authorize a wholesale rewrite. Kernel
and parallel speedups cannot be multiplied into an unmeasured whole-factor claim.

Keep a selected native prime-generation backend self-contained in Factor:
bundle its C sources and build adapter instead of depending on a sibling
checkout or an installed external prime-sieve library. Compiler/CFFI and any
optional OpenMP runtime remain explicit accelerated-backend requirements;
the Python/PyPy fallback must work without the native extension. Preserve the
source pin and adaptation record locally and validate the bundled source.

E1 combined confirmation is deferred. GNFS stays deferred
until H1 and a subsequent coverage review activate A9. Individual accepted
tranches may ship when their own gates pass; optional losses can settle a
predecessor. A task ID, research note or this backlog alone authorizes no run.
Default promotion requires the stated gates or an explicitly authorized,
separately recorded exception. Record every unmet gate and limitation.

## How to use the gates

- **Acceptance gate (A):** the observable behavior required to close a TODO. Every returned split must satisfy `1 < g < n` and `n % g == 0`.
- **Experiment gate (E):** the comparison or adversarial experiment required before integration or changing a default. For correctness repairs, this means passing regression and fault-injection cases; speed is not a prerequisite for fixing incorrect behavior.
- **Phase exit:** the condition for moving the dependent production work into the next phase. Independent research spikes may proceed once their stated prerequisites exist.
- **Evidence:** store the implementation commit, command, corpus identifier, seeds, environment, raw results, and decision in a separate artifact for each task. Preserve the original audit JSON as historical evidence.

**Performance promotion policy for Phases 3–8 — revised 9 October 2026;
wall-time objective and balanced guardrail clarified 2 November 2026:**
The 4 November timing correction removes universal percentage-based canary
vetoes from new performance decisions. The historical C9a 15% tolerance was
an engineering heuristic, not a threshold derived from host variability.
Preserve its raw observations and recorded flags for historical replay;
report the flags as diagnostic warnings in new analysis rather than discarding
otherwise validated timings. Actual correctness, resource and recorder
failures remain distinct. Calibrate variability with unchanged-code A/A and
matched references, quantify uncertainty, and justify any future rejection
rule prospectively from those measurements. A descriptive canary warning
alone neither proves nor disproves a speedup. Independent inputs and fresh
confirmation below still govern adoption.

Require zero correctness failures and the same finite work, wall/CPU time,
memory and storage ceilings within a matched comparison. For latency studies,
freeze a machine-wide aggregate CPU allowance appropriate to the permitted
concurrency and available cores; both arms receive that allowance even if a
serial arm cannot consume it before the wall deadline. CPU consumption need
not be equal. Compare equal-CPU efficiency separately when making that claim.
A repeatable reduction in complete-run wall-time median/p90 or increase in
completion within budget can qualify at any percentage.
Confirm the effect on fresh held-out inputs with repeated seeds, stable
samples and a prespecified uncertainty analysis, normally a 95% confidence
interval supporting a positive improvement. Account for run-to-run noise,
measurement order and drift; a positive point estimate alone is insufficient.
Include setup, conversion, timeout, recovery and output-validation costs.

Report absolute and relative benefits, per-class time/completion results and
uncertainty. Outside the protected 60/70-digit balanced-semiprime classes,
allow no more than a 5-percentage-point completion regression in another
declared class; explain any time, completion or resource regression
and narrow the promoted scope where warranted. Weigh the confirmed benefit
against implementation complexity, maintenance and arithmetic/resume risk.
Small, simple improvements can be worthwhile; more invasive changes need a
stronger benefit. Kernel-only wins cannot establish complete-factoring gains,
and claims remain limited to the confirmed workload/resource classes.

Freeze sources, inputs, seeds, bounds, budgets, metrics, sampling/stability
rules and regression limits before timing. A change-specific minimum useful
effect may be justified before measurement, but 10% time reduction and
10 percentage points of completion gain are no longer universal gates.
Inconclusive evidence retains the baseline; extend measurements only within
the declared study limits. Preserve previous frozen protocols and historical
verdicts. Reconsider their evidence through an explicit policy reassessment,
without retuning on held-out inputs or treating the policy change as automatic
production promotion.

## Dependencies and priority

Phases A–I are execution layers, distinct from the roadmap's numbered phases.
Rows in one layer have no dependency on one another's unfinished deliverables.
Each row names its predecessors; it may start as soon as those predecessors
are settled, without waiting for unrelated rows in an earlier layer. Deferred
GNFS rows retain their IDs in a separate later-tranche table; their original
letters no longer indicate immediate eligibility.

The plan includes two kinds of edge:

- **Required:** an interface, oracle, implementation or measurement needed by
  the successor. These are stated in the predecessor column.
- **Evidence/order gate:** a deliberate choice to evaluate a cheaper control
  before spending on a more expensive challenger. These are labeled in the
  reason column; they are not mathematical necessities.

E1's deferral does not bar justified isolated ECM, current-collector QS workers,
native boundaries, preprocessing/classification, rho, continuation or context work
under current committed controls. Those components require their selected
exactness, ownership and measurement contracts; later integrated portfolio
tuning retains D7 and combined acceptance retains E1.
Source bounds, expected curves, SIMD instruction support and core counts are
selection evidence, not proof of a local performance benefit.

An experimental predecessor is settled by an accepted implementation **or a
documented retain-baseline/defer/reject decision**. If a required mathematical
capability is missing, its dependent implementation stays deferred. A recorded
deferral cannot stand in for a correctness proof.

Task owners identify arithmetic, API and storage/provenance responsibilities;
they do not require separate chats or isolated successful performance studies.
A coherent workstream may select independently flagged mechanisms and compare
one combined candidate with the current control once every selected path is
exact and feasible under declared finite resources. Settle the relevant
interface/proof contract within that workstream; an unrelated global milestone
or standalone speedup is not an admission prerequisite.

Freeze the selected combination and a small targeted ablation plan before its
measurements. Use fixed-work prefixes or candidate/residual/relation replays
for concrete mechanism questions, preserving meaningful compatible flags.
Add whole-factor ablations when they resolve important interaction uncertainty;
do not automatically run every flag combination. Changed admission populations
need declared coverage checks rather than identical outputs. Adoption still
requires matched complete-cost/completion evidence and fresh independent-input
confirmation. Preserve old freezes, failures and scoped verdicts.

Priority means expected return per engineering/measurement effort, inferred
from current evidence, not a promised speedup. **High** gets the next available
slot; **medium** is useful independent work; **conditional** requires its stated
trigger; **low/conditional** is deliberately a later investment. Within a
phase, prefer higher-priority rows when slots are limited.

Before a suffix follow-up, reconcile its parent's final report and frozen
scope against AP01–AP40. Classify each candidate as covered, missing, rejected
or conditional. Implement only the uncovered delta. A suffix does not add a
prerequisite to the parent's closure. Changed integrations use a new source
identity; do not rerun every historical comparison unconditionally or add
component gains into a portfolio claim. Refresh representative matrices after
collection changes before promoting a solver. A timeout or zero useful rows
alone does not justify DLP or a more elaborate solver.

## Model and effort key

For **new** tasks, **Sol** means **GPT-6.1 Sol** (`gpt-6.1-sol`) and **Astra**
means **GPT-6 Astra** (`gpt-6-astra`). The completed R2 tranche used its
existing configuration; these recommendations apply to new tasks.

Use Sol for bounded implementation, measurement and integration; use Astra
for difficult algebraic/provenance invariants and new solver designs. `high`
fits controlled execution and reconciliation; `xhigh` fits coupled correctness
or selection decisions. Escalate to `max` only for a named unresolved proof,
failed invariant or derivation, not as the default for a long benchmark run.

These assignments are engineering recommendations, not results of a model
benchmark on Factor. They follow the general distinctions in
[OpenAI model-selection guidance](https://developers.openai.com/api/docs/guides/model-selection).
The [Sol](https://developers.openai.com/api/docs/models/gpt-6.1-sol) and
[Astra](https://developers.openai.com/api/docs/models/gpt-6-astra) model pages
confirm the model names and supported effort settings.

The roadmap's task rows are the source of truth for scope-specific assignments.
Use a fixed model/effort for the requested chat phase; these recommendations do
not assume autonomous model switching or authorize delegation. A later,
separately requested frozen-execution chat may use Sol / high when the row
allows it. C10a-01 retains Astra / xhigh throughout its current workstream,
including authorized execution. Expanded C9a admission/graph/storage research,
implementation and correctness use Astra / xhigh; Sol / xhigh remains suitable
for a separately scoped existing-policy cadence task. Reserve Astra / max for
a named unresolved invariant or derivation, not elapsed time or source volume.

## Parallel execution and ownership

| Lane | Scope and collision rule |
| --- | --- |
| QS collector | Use settled B1/C9/C10/A1a/C10a-07 decisions and the current committed control. A7a owns adapting existing workers to that collector; C10a/C9a/A1a/A1c retain their scopes. Preserve the active C10a-01 freeze. A1b is Done and A1b-F1 stays Deferred. B13/C8 retain their capacity/eligibility triggers. Coordinate shared collector, relation and checkpoint edits. |
| Matrix/provenance | A8/B6/B7/C4/D1/D2 share exact row/operator/lifting contracts. A selected C9a compact-record implementation may satisfy its B7 contract within the same workstream. Coordinate relation and checkpoint edits; ownership does not mandate a separate chat. |
| ECM arithmetic | P4.3/P4.1/P4.2/P4.4; separate oracle/record work from production `ecm/core.py` and `execution/stage_jobs.py` integration. |
| Schedules/continuations | P5.2/P5.3 and p+1 coordinate `execution/schedules.py`, stage jobs, backend boundaries and checkpoint versions. |
| CLI/dispatch | Completed A11 covers explicit QS/MPQS/SIQS access. C3a/G1 own remaining allocation/handoff deltas. Coordinate `factor.py`, configuration identity and usage guides. |
| Classification/proofs | A10a/B14 coordinate `common/utils.py`, result certainty, CLI and checkpoints. Other v1 transfers retain their arithmetic, schedule or sieve owner. |
| Evidence/integration | One integration writer reconciles shared configuration, portfolio, benchmark guide and roadmap decisions. |

Worktrees isolate edits, not CPU, RAM, disk or thermal conditions. Reserve one
machine-wide performance window. Pause competing tests, compression, profiling
and heavy correctness runs during accepted timings. Integrate sequentially and
compare the exact combined source. There is no benefit in launching every
eligible row at once. Current allocation belongs in STATUS.md.

## Acceptance, evidence and stop rules

1. **Separate workload classes.** Train on prespecified total-size bands such
   as 30/40/60/70/80/90/99 digits, with balanced inputs distinct from uneven
   smaller-factor bands and smooth/close/power controls. Feasibility probes
   choose finite allowances before fresh confirmation; they are not held-out
   performance evidence. Do not substitute a single 50-digit target for the
   primary 60/80/100-digit workload or its 60/70-digit balanced guardrail.
2. **Require exact outputs and finite resources.** Preserve proper divisors,
   complete reconstruction including unresolved cofactors, certainty labels,
   nonunit recovery, cumulative work/time/storage and checked resume. Add
   checkpoint identity/budget handling alongside each implementation, not
   after the optimization is finished.
3. **Measure the whole consequence.** Include setup, import/conversion where
   relevant, planning, failed attempts, splitting, filtering, lifting,
   serialization, replay and retained storage. Report owned bounds separately
   from observed RSS. Benchmarks must use matched inputs, seeds and budgets,
   at least three seconds of validated PyPy warmup and nine samples, extended
   when unstable. Cold startup and instrumented profiles stay separate.
4. **Use the revised promotion policy in “How to use the gates”.** Require zero correctness
   failures and a repeatable complete-run time or completion improvement
   supported by fresh confirmation and uncertainty beyond measurement noise.
   There is no universal minimum percentage gain. Keep the same finite
   resource ceilings and the declared parallel CPU allowance; report per-class
   regressions, preserve the 60/70 balanced completion guardrail and the
   5-percentage-point ceiling for other classes, and weigh benefit against
   maintenance/risk.
   Faster refusal, microkernel throughput, extra rows and synthetic matrix
   capacity alone do not establish a factoring win.
5. **Run relevant checks.** Implementation changes use `make -C v2 test` and
   `make -C v2 lint`; before publishing, verify tests and benchmark imports
   in a committed-files-only checkout. Update API docs, accepted behavior and
   measured summaries only after the corresponding gates pass.
6. **Use the new evidence layout.** Required immutable corpora, baselines and
   loader source snapshots belong in versioned `v2/benchmarks/inputs/`.
   New generated captures, profiles and scratch evidence go under ignored
   `results/` trees. The entire `v2/audit/` tree is local. Keep historical archives
   at their recorded restoration locations; do not recreate obsolete loose
   `_LOCAL` files or force-add generated outputs.
7. **Re-rank after each meaningful result.** If collection dominates, prefer
   calibrated collection/yield work; if sparse admission fails, require a
   storage proof; if stage-two schedules dominate, prefer reuse/pairing; if
   arithmetic dominates, use the backend/kernel evidence. Rejecting a
   challenger and retaining the baseline is a valid completed experiment.

## Experiment design and prompt review

The [repository working agreements](../AGENTS.md#experiment-design-and-prompt-review)
record the 1 November lessons from C10a-07. Keep their meaning in every design:

- Establish a small uninstrumented feasibility baseline on representative
  inputs using current committed mainline. Label historical comparisons.
- Measure diagnostic overhead on identical work with tracing enabled and
  disabled. Prefer short targeted traces or bounded counters when observation
  dominates. Instrumented completion is not clean performance evidence.
- Separate admission from adoption. Exact arithmetic, validated outputs,
  traversal/resume and feasible finite bounds admit a comparison. An unrelated
  reference need not factor every input or emit verified dependencies first.
  Missing explanatory attribution is uncertainty, not an automatic veto.
- Require matched complete-cost/completion evidence, independent inputs,
  uncertainty and disjoint fresh confirmation for adoption. Repeated seeds
  are not independent moduli; exploratory completion differences justify
  investigation without establishing a speedup.
- Give each stage a question, advancement/stopping decision, estimated total
  runtime and finite cumulative resources. Calibrate estimates from informative
  observations. Reassess recorder failures or unanswered timeouts before
  spending more hours; larger repeated campaigns are not a default remedy.
- Distinguish recorder storage/event exhaustion from production exhaustion.
  Preserve receipts, charges, failures and freezes. Revise a flawed protocol
  under a new identity; never selectively retry unfavorable observations,
  retune on confirmation or relax correctness.

## Evidence layout

Required immutable corpora, baselines and loader snapshots live in versioned
`v2/benchmarks/inputs/`. Generated captures, stdout, profiles, verification dumps
and scratch go in ignored `results/` trees. `v2/audit/` and `v2/local/LOG.md` are
local. Detailed research, comparison reports and result JSON stay in these
ignored trees; do not maintain a selective public `v2/studies/` collection.
Keep reusable runner code under `v2/benchmarks/` and concise decisions in
[benchmark results](benchmarks/RESULTS.md). Historical reports remain available
through [pinned Git history](history/README.md), with original source identities.
External sources stay in folder-ignored `v2/benchmarks/inputs/upstream/`;
keep licenses, inventories and comparison provenance locally. Do not force-add
generated evidence or create individual-file study ignore rules. Preserve
historical archives at their recorded restoration locations and lossless
compression manifests. Publish concise commands, findings and limitations.

The [history record](history/ROADMAP_CLEANUP.md) resolves the old removed-link
labels and records the historical Git revision, source hash and original link context.
The source audit used
`1272e033f889a105792cbb924bf8a12a46ac88ae`; the C-sieve review used
`5b4afb8f344ad5f6fbd20a186cd57b36182ea710`. Their native thresholds and OpenMP
results are reference evidence, not PyPy defaults. Audit compatibility scripts
load Python 2 sources, may depend on `lib2to3` and original absolute paths, and
are historical diagnostics rather than native acceptance tests. Independently
verify contracts and pin external implementations before adaptation/timing.
