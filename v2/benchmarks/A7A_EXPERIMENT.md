# A7a execution specification: overhead removal, then anchor tuning

Status: revised execution design, not a completed campaign or approval to
launch it. No successor stage is pre-approved. Follow the
[joint-search policy](../POLICY.md#joint-parameter-searches) and
[evidence tiers](../POLICY.md#evidence-tiers-and-sample-sizing).
Preserve historical freezes, failures, charges and captures unchanged.

## Deliverable, scope and source

Deliver cold CLI rows for **40, 50, 55, 60, 65 and 70 digits**. Add warm rows
only at **40 and 50 digits**. **62 digits is dropped** from new measurements;
retain its historical evidence without claiming an adopted interpolation.

Tune cold anchors at **40, 55 and 70** and warm anchors at **40 and 50**.
After architectural fixes, interpolate cold 50/60/65 from those anchors and
verify each of these three rows. They remain labeled "interpolated, verified",
not independently tuned optima. Warm/cold results are distinct lifetimes.
Report full recipes, input coverage, completion, elapsed time, CPU, aggregate
owned RSS, uncertainty, ties, exclusions and untested interactions. The claim
is "best found and confirmed within this search" at anchors, and the stated
verification claim at intermediate rows, not global optimality.

There is **no approved total budget for ordinary work**. All ordinary-stage
limits below are proposed maxima, not granted time. The duration column
contains provisional planning estimates; replace them with calibrated stage
estimates before approval.
The earlier 48-hour proposal is withdrawn; no stage may assume that allowance.

75/80 are separate single-candidate feasibility gates, never ordinary search
or confirmation rows. The 75-digit probe runs only after G3, with the G2H
fixes and frozen G3 recipes, and explicit approval at or after the G3 report.
It may run before or after G4; G4 does not depend on it. The 80-digit gate
requires the 75-digit report and a separate instruction. Neither probe starts
automatically.

Continue in `/Users/nishanthmohan/Documents/code/factor-a7a-regime-tuning`,
branch `codex/a7a-regime-tuning`. The implementation checkpoint is `276ebd8`,
based on mainline `e62ec18`; historical experimental source is `ae1cfd9`.
Later documentation commits do not establish fresh runtime evidence. Record
actual committed source, corpus hash and runtime versions in each new freeze.
Pin literal current mainline separately. Do not merge implementation changes,
push or change automatic defaults under this plan.

## Gates and finite stage contracts

**Stop, commit tracked changes, report and await explicit approval at every
exit gate. No stage starts without approval.** Plan revision is not approval
of preparation, QA, measurements or a successor. At G1, G2, G2H, any G2K, G3
and G4, deliver a progress table, tests/evidence, exact charges, remaining
resources and the priced proposed successor. Report the first informative
result during a stage, without treating that update as permission to advance.

| Stage / exit gate | Question and resulting decision | Expected duration before calibration | Finite new measurement ceiling |
| --- | --- | --- | --- |
| Preparation / G1 | Are the bounded driver, corpus, mainline adapter, counters and progress/thermal logs tested and ready? Report even if incomplete. | One working day of implementation; up to 1 h locked QA | 8 active implementation hours; 2 h QA wall, 20 core-hours, 2 GiB captures |
| Coverage / G2 | Do complete calls work across every admitted row, and what is the priced whole schedule and dominant cost? | 1–3 h, replace with the G1 frozen estimate | 5 h wall, 50 core-hours, 4 GiB captures |
| Architectural overhead / G2H | Which measured serial/coordination cost can be removed, and does the fixed bundle improve complete runs before tuning? | 1–2 h measurement; up to one working day engineering | 8 active engineering hours; 4 h measurement/QA wall, 40 core-hours, 8 GiB captures |
| Conditional cofactor engineering / G2K | If the measured cofactor share warrants a separate queue, is its exact contract and complete cost worthwhile? Otherwise defer. | Must be priced separately at G2H | At most 8 active engineering hours, 2 h measurement/QA, 20 core-hours, 2 GiB captures; requires separate approval and repricing |
| Anchor search and intermediate verification / G3 | Which (c) recipes survive, and does the interpolation rule pass? Freeze all recipes and references. | Price from post-fix complete costs before approval | 10 h wall, 100 core-hours, 16 GiB captures, including diagnostics/refinement |
| Independent confirmation / G4 | Do the frozen recipes beat literal mainline with adequate uncertainty and protected completion? | Price at G3, before approval | A separately priced wall/CPU/capture grant, approved at G3 before timing |

Ordinary ceilings include preparation QA, coverage, overhead, search,
confirmation, warmup, diagnostics, telemetry and cleanup, not just actor time.
Engineering time is recorded separately; no timed grant is an instruction to
spend it. The optional K stage is not an extra automatic allowance: fund it
by an explicitly approved redistribution or revised stage grant.
75/80 have separate limits below and cannot consume ordinary confirmation
reserves. Keep exact wall/CPU/output charges in each stage's own receipts,
and one cumulative resource total in the shared progress table. New grants
name only their stage and finite allowance; they do not require a recursive
historical charge chain or a sealed-ledger check. Preserve historical evidence.

At G2 publish expected cost and a conservative ceiling bill for **every**
remaining row/stage. At G2H update that bill after any source change, before
search approval. At G3 price confirmation again using frozen recipes and
literal mainline. Price the initial three-input confirmation block separately
from the optional three-to-six and six-to-nine extensions. If the initial
block cannot fit the budget actually approved by the user, stop and propose
cuts or a revised grant before spending the search allocation. If no total
is approved, publish the full forecast anyway and request only the named
next stage.
Do not discover missing coverage after a deep 60-digit campaign.

Stop earlier at the first finite limit, 1.5 times the frozen stage estimate,
or the first informative result contradicting its assumptions. Stop issuing
work, cancel/drain owned children, preserve checkpoints/charges, release the
owned lock after cleanup, report and await a revised instruction. Do not
silently extend a stage or redefine the estimate. All-arm uninformative
timeouts require a proposed redesign, not automatic larger limits.

## G1: preparation, bounded to a working day

Extend the existing benchmark runners; **do not build a general optimization
framework**. `joint_search.py` and `regime_search.py` are factories, not a
completed adaptive campaign. A small deterministic batch driver plus an
append-only coverage table and the existing ledger is sufficient. Manual
review at gates is deliberate. Cap driver work within the eight-hour G1
implementation timebox and report at G1 if unfinished; do not expand scope.

Required before coverage:

1. A fixed stage/row manifest, precomputed proposal list, gate-approval identity,
   source/corpus binding and per-stage recipe cap. A restart reads that stage's
   receipts, skips reported identities and preserves failures/charges. Reject a
   missing row, duplicate issue, unapproved stage or cap overrun. Reuse existing
   runner locking, accounting and cleanup rather than a new scheduler service.
2. A separately named immutable corpus with nine fresh certified balanced moduli
   per ordinary digit regime, generation seed `202611071442`. Validate prime
   certificates, digit/factor sizes and pairwise coprimality with existing
   corpora. Reserve it for confirmation; preserve `regime_corpus.json`.
   Primitives: `certified_prime`/`verify_certificates` in
   `v2.benchmarks.suites.build_phase_two_corpus`, `GenerationRandom` in
   `v2.benchmarks.qs.c9a.inputs`; prior local construction is
   `v2/audit/a7a/m_fresh_inputs_prepare.py`.
3. A literal committed-mainline adapter, pinned at execution. A7a `S0/Y0` is
   a same-policy explanatory serial control, not literal master. Use the same
   declared complete-QS CLI endpoint, cold ownership, seed and allowance; do
   not add unrelated portfolio engines to one arm. Measure decision/setup
   cost. Compare both serial controls during coverage; confirmation uses
   **literal mainline only**, not a third explanatory arm. This does not claim
   speedup against the strongest serial configuration without its own evidence.
4. Cheap timeline counters: worker active/idle/publication-blocked intervals;
   parent setup, scheduling, transport/decode, verification, graph admission,
   cofactor, filtering, solve and recovery; backlog and high-water marks;
   useful relations, bytes and work. Separate overlapping CPU spans from
   critical-path wall intervals. Existing Tier A's packet-wall-derived idle
   capacity and parent CPU/wall are proxies, not direct idle or Amdahl shares.
   Add only missing bounded outer-loop observations. Measure counters on/off
   on identical work; if their effect cannot be bounded within measured A/A
   noise, reduce them or confine them to diagnostics, keeping clean calls clean.
   Reuse these counters for the clean probes' 60-second progress and thermal
   log; preserve the last observation on timeout. No probe profiler, heap
   dumper, segment rotation or new capture-storage system is a prerequisite.
5. Exact-output, lease, crash, resume, telemetry/recorder failure, near-limit,
   no-leaked-child and pure planning checks. Code changes require full suite,
   lint and benchmark imports from a committed-files-only checkout under the
   shared lock. Documents alone do not require another factoring run.

G1 reports ready/missing items and tests, commits the bounded implementation,
and stops. No successor starts automatically. The 75-digit probe additionally
requires completion of G3 and explicit approval at or after its report.

## G2: complete coverage before deep tuning

Use certified balanced calibration inputs
`regime_calibration_{digits}_{0,1,2}`. Known factors validate outputs only.
Freeze seeds `104729`, `130363`, `155921`; match them between arms and do not
choose seeds from outcomes. For breadth use input 0 and seed `104729` on all
six cold rows and both warm rows. Arms: literal mainline, strongest valid
same-policy serial, retained process incumbent. Freeze exact serial/incumbent
recipes before results; deduplicate only genuinely identical executions.

| Digits | Actor wall ceiling | Anchor search or intermediate-check owner ceiling |
| --- | ---: | ---: |
| 40 | 30 s | 15 min cold anchor + 15 min warm anchor |
| 50 | 90 s | 15 min cold intermediate + 45 min warm anchor |
| 55 | 180 s | 75 min cold anchor |
| 60 | 300 s | 45 min cold intermediate |
| 65 | 900 s | 90 min cold intermediate |
| 70 | 1,800 s | 300 min cold anchor |

These owner ceilings total **10 hours** for G3, not ten hours per lifetime.
They include its controls, warmup, diagnostics, refinement and cleanup.
Subcaps/configuration maxima are upper bounds, not promises that every maximum
can be exercised simultaneously. Price the actual proposed cells first.

Ordinary calls use CPU at most ten times actor wall, aggregate RAM 22 GiB,
parent QS 4 GiB, worker/runtime 512/256 MiB, checkpoint 64 MiB, wire 4 MiB,
atom batch 2,048 and logical work 96 trillion. Enforce the enclosing owner cap
including setup/draining. Admit constructor/accounting bounds or report the
exclusion; do not silently shrink recipes. Freeze Y0 numerical policy and
resolved base/family/report/product/residual limits per regime. Mainline
keeps its literal production policy; different numerical work is not a pure
concurrency claim. Numerical/backend changes need a new identity and repricing.

Always set `endpoint="complete"`; `regime_search.cells()` defaults to
prefixes above 55 and must not be used unchanged. Each cold repetition uses a
fresh interpreter and pool with `samples=1` and zero executed warmup. Warm
rows use at least three seconds of validated warmup and a declared retained
lifetime; charge it. Estimate A/A timing noise on bounded identical work
expected to take under ten seconds, with three to five repetitions. Do not
add repeated long complete calls merely to estimate timing noise. Prefix
noise alone does not estimate full-call input effects: use paired complete
candidate/mainline ratios across inputs for that. Counterbalance arms and
retain near ties. Coverage input 0 establishes feasibility, not a selection.

Collect only representative family diagnostics, not a companion for every
candidate/repetition. Declare families by backend/execution path before runs:
serial/local packet, process with parent verification, process with worker
verification; add K only if separately admitted. Use at most one exploratory
CPU/memory diagnostic per family, normally a window up to 45 seconds. Full
CPU/heap diagnostics are reserved for finalists or a named anomaly, each with
a question and finite allowance. Measure observer overhead on identical work.
Record source/input/phase and limitations; one window does not cover all phases.

At G2 report coverage, the first per-arm timeline/cost table, measured cofactor
and graph-insertion shares, and the whole priced schedule; then stop. A
historical 65–80% parent share is a hypothesis to check on current committed
code, not a measured present serial fraction or guaranteed removable cost.

## G2H: remove architectural overhead before tuning

Use coverage timelines to rank the largest measured critical-path cost:
repeated setup/computation, allocation/GC, transport/copying, parent admission,
budget polling or scheduling stalls. Review heap owners and CPU attribution
together. Implement and validate one concrete dominant fix first. Review the
result before another fix; do not optimize a small matrix stage while the
parent stalls collection. Price and freeze a bounded next iteration; if its
need contradicts the stage assumptions or exceeds the cap, stop and report.

Measure the **N-independent-workers throughput upper bound** at N=1/4/8 using
matched deterministic disjoint packet work on three calibration inputs
per regime at 50 and 65 digits. Remove shared
parent admission/transport from that diagnostic, not arithmetic or validation.
Compare aggregate useful-work throughput and a local-packet/no-pool control
with the real timeline. This bounds a coordination opportunity on the actual
heterogeneous host; N times serial rate is not an empirical bound, and
independent full factorizations do not model cooperative collection. Include
startup, remaining serial stages and load imbalance in any completion model.

### Classify and screen before constructing a joint design

The following is the initial option inventory. Each option has a one-line
reason; alternative values for one mechanism share a row. Record every
available collector option, including inactive/unsupported ones, in the
frozen manifest.
Do not silently run the old 32-row all-binary design.

| Option(s) | Class | Reason / initial screen |
| --- | --- | --- |
| Per-packet budget leases | (a) | Replace repeated shared budget calls with bounded exact reservations/reconciliation for the same work. |
| Persistent worker state (`reuse_state`) | (a) | Avoid rebuilding identical immutable context between packets; charge retained memory. |
| Coefficient reuse (`coefficient_cache`) | (a) | Avoid recalculating identical coefficients without changing assignment coverage. |
| `density` | (a) | Skip redundant marking overhead while retaining exact hits. |
| `unroll4` | (a) | Change loop shape for identical marking work; no gain is presumed. |
| Packed versus object transport (`packed_transport`) | (b) | Alternative relation representation/publication with exact decode/admission. |
| Shared/private event/semaphore stopping | (b) | Alternative stop signalling under an identical admitted lease/check contract. |
| Parent/worker/batch64/batch256 verification | (b) | Alternative placement/batching of the same exact relation checks. |
| CPU snapshot versus repeated reads (`cpu_snapshot`) | (b) | Alternative charged CPU observation with a tested conservative staleness bound. |
| `presieve` versus `presieve_powers` versus neither | (b) | Alternative small-prime marking implementations under the same domain. |
| `scan_fused` | (b) | Fuse scan stages instead of separate passes over identical candidates. |
| `parity_two` | (b) | Alternative exact two-root marking implementation. |
| `packed_queue` | (b) | Packed versus object queue representation for identical entries. |
| `wide_hits` | (b) | Alternative hit computation at a fixed capacity and arithmetic domain. |
| `wide_inverses` | (b) | Alternative inverse computation/cache layout at fixed capacity. |
| `wide_transfers` | (b) | Alternative root-transfer representation at fixed capacity. |
| `pruned_transfers` | (b) | Omit provably redundant transfers while preserving exact roots. |
| Serial/process worker count | (c) | Startup, parent pressure and heterogeneous-core throughput vary with size. |
| Polynomials per packet, lookahead | (c) | Amortization and load balance interact with workers and family geometry. |
| Arrival versus ordered admission (`arrival_order`) | (c) | Relation trajectory/yield can interact with scheduling; label it trajectory-changing. |
| Nursery (`large_nursery`) | (c) | GC frequency and retained memory interact with packet size and concurrency. |
| State/root/cache capacities | (c) | Reuse and memory footprint vary with base size and concurrency. |
| `sparse_queue` byte ceiling | (c) | Queue capacity and rebuild cost depend on size and collector path. |
| `power_plan` byte ceiling | (c) | Preparation/reuse tradeoffs depend on base size and call lifetime. |
| Shared cofactor allocation K | Conditional (c) | A new exact queue requires measured cofactor headroom and its own G2K engineering gate. |
| Whole-family flag | Derived | Derive from packet width and actual family size; do not count a redundant independent dimension. |
| Numerical/base/large-prime policy | Frozen boundary | Changing mathematical effort needs a separately priced design, not a hidden concurrency comparison. |
| Native/NEON | Excluded | Deferred B5a implementation, not an available option in this process campaign. |

Screen (a)/(b) once at **50 and 65 digits**, on identical packet work and
fixed resources/seeds. Apply Tier 1: three independent calibration inputs per
representative regime, with one matched observation per input and comparison
against A/A noise. Maximum 28 distinct screen recipes per representative
regime, including baseline and alternatives. Repeat an input only to estimate
within-input noise for calls under ten seconds (three to five observations,
predeclared and priced); repetitions do not replace independent inputs.
Calibrate a fixed packet quota expected to take at
most five seconds per side, with a 45-second hard cap; price the full screen,
complete-run rechecks and QA inside G2H before approval. Freeze exact work
before results; a work prefix can diagnose identical-work overhead without
proving a full-factor ranking. Reject invalid signatures/bounds first.
Reuse earlier evidence only if committed runtime equivalence and these
work/input conditions are established; otherwise label it historical.

Enable (a) unless it clearly loses; inconclusive means the simpler valid path.
For (b), retain the supported winner; inconclusive means the simpler valid
implementation. Do not let this rule disable an already required correctness
or budget mechanism. A measured interaction (for example, reproducible
verification-placement reversal with worker pressure) may reclassify that
specific option as (c), after reporting evidence and approving a revised
finite design. Uncertainty alone is not interaction evidence. Do not promote
all inconclusive flags to joint dimensions.

Validate the fixed architectural bundle with matched complete calls on
three calibration inputs per regime at 50 and 65, then recheck the six cold
rows before freezing search source. At G2H
report attribution, screen decisions, fixed options, exclusions, full-run
results and the repriced search/confirmation bill; stop for approval. New
hotspots discovered later require a separately gated fix and retuning of
only affected rows, not an automatically repeated whole campaign.

### Conditional G2K: shared cofactor jobs

K is unimplemented and offers no guaranteed speedup. At G2H report measured
cofactor worker CPU share, critical-path/imbalance evidence and a realistic
recoverable-wall bound. Only if that supports a material opportunity, propose
a costed separate G2K stage and obtain approval before implementation.
Otherwise defer with the measured share. No nested pools. Specify exact
ownership, work/CPU leases, malformed/stale publication, cancellation, crash,
resume and overshoot; test alone/combined and run committed-checkout QA.
Keep K=0 evidence distinct. After G2K stop and report before anchor tuning.
Admitted allocations have total collection plus cofactor processes at most ten.

## G3: tune only (c), then verify interpolation

After G2H (and any approved G2K), freeze the implementation and (a)/(b) choices.
Initial (c) ranges:

| Tunable | Declared levels / conditional constraints |
| --- | --- |
| Collection workers | Serial, 1/2/4/6/8/10 initially; 3/5/7/9 as adjacent refinements. |
| Packet width | Valid family fractions g/16, g/8, g/4, g, plus adjacent valid powers of two; integer width <= actual gray limit g. |
| Lookahead | 0/2/4/8; inactive with no process pool. |
| Admission order | Ordered/arrival; arrival is a separately labeled mathematical trajectory. |
| Nursery | 4/16/64/128 MiB, selected before interpreter startup. |
| State cache | 16/64/128 MiB within charged per-worker memory; inactive without retained state. |
| Sparse queue | 0/8/32 MiB, only supported collector paths. |
| Power plan | 0/2/8 MiB, with preparation charged. |
| Other implemented root/inverse cache capacities | Only actual validated configuration fields; enumerate levels/dependencies at G2H or exclude with a reason. Never invent unsupported runner fields. |
| K, only after G2K | K=0/1/2 and collection allocation with total <=10; fixed quota, no oversubscription. |

Filter invalid/inactive dimensions before design construction. Larger caches
must coexist with worker/runtime/parent allowances. Serial has no process
scheduling dimensions. Reclassified options enter only with their measured
interaction and approved count/budget adjustment.

At each anchor/lifetime, form a small deterministic covering/space-filling
pool over (c), holding (a)/(b) fixed. Force serial plus one valid recipe at
1/2/4/6/8/10 workers. Choose remaining proposals by uncovered meaningful
parameter pairs, with frozen tie order; report actual conditional coverage.
Maximum **16 initial + 8 refinement** distinct recipes per anchor/lifetime;
for calibrated calls over five minutes, **8 + 4**. Controls count toward
these maxima; repetitions are not new configurations. Price the executable
subset inside the row cap, leaving room for independent-input survivors,
exploration and diagnostics before issuing it.

Use complete-factor input 0 to prioritize proposals, retaining uncertain/
near-tied candidates. Advance up to four to calibration inputs 1/2. A Tier 1
selection or rejection needs all three independent inputs and comparison
against measured noise; candidates not advanced remain deferred/untested,
not demonstrated losers. Inconclusive choices retain the simpler valid path.
Nine-input screening is not required. Reserve at least one proposal outside
the current neighborhood; refine one/two (c) levels and adjacent worker/
packet/cache levels.
An unrelated reference timeout does not block a valid candidate. Prefixes do
not rank complete-factor candidates without validated predictive agreement.
Finish a breadth round at every anchor before a second deep round at one.

Freeze this interpolation rule before intermediate outcomes: nearest cold
anchor by decimal digits, lower anchor on an exact tie. Inherit its worker,
lookahead, order and capacity levels; map packet fraction onto the target's
actual family using the declared valid integer levels. Apply only publicly
observable size/backend/allowance inputs. Verify cold 50/60/65 against literal
mainline on all three calibration inputs 0/1/2, at most one proposed recipe
plus the reference per row. If a check fails, stop at G3 and propose a separately
bounded new anchor; do not quietly optimize every intermediate row. Label
unmeasured boundaries and do not adopt a routing table before their checks.

At G3 freeze every recipe, reference commit, interpolation label, analysis,
confirmation cost tier and complete priced schedule. Include the exact resolved
75-digit probe recipe and its derivation before any probe result: the source/
identity of the cold 70-digit anchor, inherited worker count, lookahead,
admission order and capacity levels, fixed G2H implementation/(a)/(b) choices,
source packet fraction, actual 75-digit family and mapped valid integer packet
width, resolved Y0 numerical settings and admission under the probe allowances.
Configuration resolution does not authorize a factoring call. Report and stop;
the 75-digit call needs explicit approval at or after this report. G4 can
proceed independently under its own approval, before or after that probe.

## G4: independent confirmation against literal mainline

Use only the frozen candidate and literal mainline: two arms, not three.
Identical recipes need no duplicate timing and yield no speedup claim.
Retain same-policy serial coverage as explanatory evidence without treating
it as a separately confirmed comparison. No confirmation inputs may tune a
recipe, interpolation, work quota or seed.

Measure preliminary sigma during G2/G3: the across-input standard deviation
of log(candidate time / mainline time) on matched pairs. A/A estimates timing
noise, not candidate-specific input variability. Refine sigma using the first
fresh paired block; never assume a 5–10% value or transfer an incumbent's
variance to a changed recipe as if it were measured. Calibration inputs do not
count as fresh confirmation inputs.

Start each row with **three fresh independent inputs**, in committed order:

| Complete-call cost, frozen at G3 | Initial allocation per arm |
| --- | --- |
| Under about 10 s, where within-input noise needs estimation | 3 inputs x 3 interleaved repetitions; up to 5 only if chosen from calibration before confirmation |
| Otherwise | 3 inputs x 1 repetition |

One cold repetition means one fresh interpreter/pool, with no warmup. Warm
applies only at 40/50 with the separately frozen lifetime and at least three
seconds of validated, charged warmup. Freeze the repetition count before
confirmation; use one per input for calls at or above ten seconds. Treat
repetitions as a cluster and form one paired log-ratio observation per input.
Counterbalance arm order, preserve failures and reconstruct exact results.

Use prespecified **98 1/3% intervals** at looks after **3, 6 and 9 independent
inputs**. This is the approximately 98% per-look rule, with at least 95%
overall coverage for one comparison if its individual intervals are valid.
Freeze/test the interval procedure, distribution assumptions and decision
threshold before timing; adjust separately for joint claims across rows.
Extend from three to six, then nine **only** while the interval straddles the
decision and the approved stage grant funds the whole next block. Stop as soon
as the decision is clear. At nine, or the funded maximum, stop/report
inconclusive if ambiguity remains. Nine is a maximum, never an initial quota.
Do not retune, change repetitions or selectively retry from these results.

Keep timeouts censored: never use their cap as an exact completed time or
discard them to form favorable paired ratios. If censoring prevents the
predeclared ratio analysis, use the frozen completion/guardrail decision and
report its uncertainty; do not declare a latency win from successful pairs
alone. Protect 60/70 balanced completion. Published README/RESULTS performance
numbers and p90 claims need Tier 3 evidence and precision; this small default/
dispatch confirmation alone does not establish a population quantile.

At G2/G2H and again at G3, price each row as inputs times frozen repetitions
times the sum of **measured** candidate/mainline complete costs, plus warmup,
controls, diagnostics, telemetry and cleanup. Price the **three-input initial
block** at G3, and the **three-to-six and six-to-nine blocks separately**.
Report conservative actor-ceiling bills separately from expected runtime.
With no approved ordinary total, no extension is funded by an assumed campaign
allowance. Keep exact stage charges and update the one cumulative total in
the progress table. No automatic new grant or extension beyond nine.

G4 delivers the whole admitted table, including explicit infeasible/inconclusive
rows and untested interactions, and stops. Confirmation is not adoption or
permission to merge/change defaults.

## Separate 75/80 feasibility probes: one candidate, then stop

75 has a **two-hour actor budget**; 80 has **three times that, six hours**.
Run 75 **only after G3**, once the G2H fixes and frozen G3 recipes exist, and
only with explicit approval at or after the G3 report. It may run before or
after G4; **G4 does not depend on the probe**. After the 75 receipt, stop and
notify. The 80-digit probe still needs that report and a separate instruction.
Exactly one candidate, one input, one clean cold call per approved probe.
No second seed, companion, repetition, configuration, larger grant or
confirmation is implied by any outcome. Report admission refusal as an outcome.

| Field | 75 digits | 80 digits |
| --- | ---: | ---: |
| Input | `regime_calibration_75_0` | `regime_calibration_80_0` |
| Collection workers / packet width | Frozen G3 cold 70-digit recipe, mapped to 75 | Frozen G3 cold 70-digit recipe, mapped to 80 |
| Actor wall | 7,200 s | 21,600 s |
| Owner wall, including observation and cleanup | 7,500 s | 22,500 s |
| Actor / owner aggregate CPU | 72,000 / 75,000 s | 216,000 / 225,000 s |
| Aggregate owned RAM | 22 GiB | 22 GiB |
| Worker / runtime per worker | 768 / 256 MiB | 768 / 256 MiB |
| Parent QS / checkpoint | 4 GiB / 256 MiB | 4 GiB / 256 MiB |
| Logical work ceiling | 96,000,000,000,000 | 96,000,000,000,000 |

**75-digit candidate:** use the frozen G3 cold 70-digit anchor recipe,
carried to 75 by the declared interpolation rule. Inherit worker count,
lookahead, admission order and capacity levels; map its packet fraction onto
the actual 75-digit family using the declared valid integer levels. Keep the
G2H implementation bundle and (a)/(b) choices unchanged. Inactive settings
remain inactive; do not force an eight-worker or packet-16 fallback. Use seed
104729 and the Y0 numerical policy resolved for 75 digits, without retuning
base/report/product/residual limits. Keep the allowances above, including
actor 7,200 s, owner 7,500 s, worker/runtime 768/256 MiB and parent QS 4 GiB.
If the mapped recipe is refused, report the refusal and stop; do not substitute
another recipe. Record the exact resolved recipe and derivation in the G3
report before any probe result; do not choose it from 75-digit outcomes.

**80-digit candidate:** use the same derivation as the 75-digit probe:
the frozen G3 cold 70-digit anchor recipe, mapped to the actual 80-digit
family by the declared interpolation rule and valid integer packet levels.
Inherit worker count, lookahead, admission order and capacity levels; keep the
G2H implementation bundle and (a)/(b) choices unchanged. Use seed 104729 and
Y0 resolved for 80 digits without retuning base/report/product/residual limits.
Keep the 80-digit table allowances. Record the exact resolved recipe and its
derivation before the call; if refused, report the refusal and stop without
substituting another recipe.

If the 75-digit receipt motivates a different 80-digit candidate, propose it
in the 75-digit report for separate approval and freeze it before the 80-digit
call. A proposal does not authorize a switch. Never choose or tune the
candidate from 80-digit outcomes. The 75-digit report and separate instruction
remain prerequisites for the single clean 80-digit call.

### Clean observation and useful timeout receipts

Run **one clean call**, with the 60-second progress and thermal log only.
No vmprof, heap snapshot, segment rotation, JIT logging, per-packet dumps or
profile companion. CPU/heap profiling decisions come from approved 60–70-digit
work; they are not prerequisites for a two-/six-hour feasibility question.
The cheap progress/RSS observers must be checked on bounded ordinary work;
record their measured overhead without adding a second 75/80 call.

- Every 60 seconds, retain elapsed time, raw/accepted/usable relation counts,
  current target, usable/target fraction, graph/cycle progress, queue backlog,
  stage and owned CPU/RSS. Include cumulative and rolling 5/15-minute useful
  relations/s. Raw atoms are not matrix-ready independent relations; define
  numerator/target from actual collector state and record target revisions.
- Log bounded `pmset -g therm` before/after and every 60 seconds during the
  probe, with monotonic timestamps, command duration and missing/error status.
  Charge in-run telemetry and label its interference; before/after telemetry
  remains outside the actor measurement. This is partial evidence, not a
  guarantee of no throttling.
- On timeout, preserve the last validated progress/checkpoint, exact input/
  source, stop reason, reconstruction/unresolved cofactor, work/wall/CPU, RSS
  and matrix-stage status. Report collection-only ETA ranges from useful-
  relation-rate windows when meaningful. Late cycle yield and a changing
  target make linear extrapolation uncertain; filtering/solve/recovery are
  not included. Do not label that estimate a measured completion time.
- A missing telemetry sample is a recorded observation failure, not a reason
  to retry or duplicate the factor call. Preserve the production result and
  exact charges; hard production limits still stop production.

### Stage accounting and execution

Use the existing `combined.cell`, `configuration`, `settings`, `source()` and
`corpus_identity("regimes")` to build a single-cell **clean-mode** plan.
Required fields: complete endpoint, cold lifetime, samples=1, zero executed
warmup, the exact table allowances and candidate settings. Disable automatic
profile companions and profile analysis. Validate the 60-second log and
last-progress receipt on a short fixture first; unsupported fields must be
rejected rather than silently ignored.

Record one approved probe identity, regime and finite owner allowance in
`grant.json`. Keep exact charges in this stage's own receipts and update the
single cumulative total in the progress table. No M1/M2/M3 prior-charge chain,
`prior_ledger`, sealed-ledger admission check or new profile-storage work is
required for a new grant. Preserve historical records unchanged. Existing
finite runner output limits remain in force; no multi-GiB profile capture is
requested. Owner cleanup time is separate from the actor factoring allowance.

Keep generated plans/receipts and ordinary diagnostic captures in the main
repository's ignored `v2/benchmarks/results/a7a/`, using a new immutable identity. For ordinary
batches and a validated clean probe, the existing CLI entry point is:

```sh
/opt/homebrew/bin/pypy3 -m v2.benchmarks.qs.a7a.combined \
  --plan /absolute/path/to/frozen-plan.json \
  --output /absolute/path/to/new-receipts-directory
```

Run from the A7a worktree with PyPy implementing Python 3.11. Ordinary
profile analysis may need its optional reader; a clean probe does not. Static
freezing does not run factoring. The runner must acquire
`/private/tmp/factor-performance.lock` with nonblocking exclusive flock; an owner file is not admission. If occupied,
report without polling or claiming a reservation. Hold through cleanup/draining
and preserve the lock file. After the single probe, **stop and notify for every
outcome**. No automatic 80 call, search, confirmation or budget escalation.
