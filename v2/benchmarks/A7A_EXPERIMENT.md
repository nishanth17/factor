# A7a execution specification: overhead removal, then anchor tuning

## Directional extension amendment — approved after the reduced race

The reduced race is complete on `bce0ad0`: 36 factoring calls, 40 controls,
provisional C4/C4/C8/C14/C8 at 50/55/60/65/70. Its 10% freeze is historical.
The next, separately source-pinned extension uses a 5% input-0 race, followed
by inputs 1–2 only for survivors. Keep packets 8/16/16/64/64 and all other
options fixed. Validate explicit counts through 20 before execution.

At 65/70 test C16 against C14, continue C18 then C20 only while wall gains
exceed A/A drift, stopping at the first loss/tie. After a strict even bracket,
test best±1 only when both even neighbours are measurably worse. Skip odd tests
on a flat drift band and break ties toward fewer workers. Apply odd refinement
at lower sizes too, including C2 when needed to bracket C4 and endpoint C1/C3
if C2 strictly wins. A still-winning C20 is an unresolved boundary, not approval
for C21/C22. Report every comparison's wall/CPU/RSS and margin versus A/A drift.
Unimodality is a search heuristic, not a proved property of these inputs.

Resolve current-source references, the finite conditional catalogue, exact
noise-control placement, call/resource ceilings and expected runtime before
launch. Retain all outcomes and source labels; no reuse of G4 heldouts for
selection. No default adoption, merge, push or new 75/80 probe before B5a.

## Worker-count race amendment — approved after slim G4

Slim G4 is complete; its source-pinned evidence is separate. The later user
instruction authorizes a reduced worker-only race: 50/55 C4/C6/C8; 60
C4/C6/C8/C10; 65/70 C8/C10/C12/C14. Keep each target's frozen/interpolated
packet width and every other option fixed. Run all 18 input-0 cells; advance
only configurations with completed wall <=1.10 times the fastest completed
cell in that row to inputs 1–2, one repetition each, at most 54 factoring calls.
Retain every screen outcome and distinguish n=1 from n=3. Explicitly flag
C12/C14 as oversubscription. Report wall, whole-actor CPU and sampled peak RSS.
Use a separate current-source protocol and original calibration inputs, never
G4 heldouts. Validate opt-in fourteen-worker support and freeze aggregate
resources before launching. No automatic-default adoption or 75/80 probes.

## Slim G4 amendment — approved 8 November 2026

Defer all new 75/80 probes until after B5a. The earlier conditional probe
instructions below are historical and do not authorize another call.

After the H3 integration passes, confirm only cold 55, 65 and 70 against S0,
using three fresh independent balanced inputs per row, one repetition per
arm/input: 18 calls, no sample extensions. Frozen process recipes are C4/p16
at 55 and C10/p64 at 70. The unchanged nearest-anchor interpolation selects
70's C10/p64 for 65; label this row interpolated. Use target-resolved Y0 and
the frozen implementation/cache/queue/packet choices on both arms. S0 is the
same-policy serial control, replacing literal mainline for this scoped G4.
No tuning, default adoption, merge to master, push or successor follows.

Use committed confirmation inputs g4_{55,65,70}_{0,1,2}, in that order, matched
factoring/hash seeds 104729/130363/155921 and counterbalanced arm order.
Each run uses a fresh interpreter, no warmup, 4 MiB nursery, matched root
wall/CPU/RSS allowances and complete cold actor wall through validation and
shutdown. Root wall limits are 100/900/1800 seconds, CPU is ten times root
wall and RSS is 18 GiB. CPU/heap profiles, Tier-A hooks and worker timelines
are off; retain low-rate RSS and 60-second thermal observations.

Predeclare one 98 1/3% paired-log t interval (df 2) per row, threshold 1 for
candidate/S0 actor time. Compare the effect with measured A/A variability;
report ambiguous or censored results as inconclusive. Never replace a timeout
with a completed time, discard failures or select settings from these inputs.
Preserve the 70-digit completion guardrail; this reduced scope does not verify
fresh 60-digit completion or establish a population/p90 or master speedup.

Each row has opening/closing identical S0 fixed-prefix A/A pairs on its
calibration input, separate from the fresh confirmation sample. Freeze exact
prefix work/signatures before outcomes and retain the 10% drift stop, with no
automatic restart. Keep resource ceilings and exact resolved recipes in the
new source-pinned local protocol. Expected elapsed time is about 2–3 hours;
owner ceilings price all 18 calls, controls, telemetry and cleanup separately.
Stop and report after this slim G4.


H3 integration amendment (8 November): local master `755c312` replaces the
resource model used by the historical stages below. Internal logical-byte,
work/CPU lease and checkpoint refusal caps no longer admit or stop production
QS workers. Keep the frozen numerical recipes, packet/cache capacities and
historical receipts. New execution uses a fresh source/protocol identity,
root wall/aggregate CPU and the approved 18 GiB observed-RSS guard. Checkpoint
file limits apply only to optional atomic saves. Cost tables require
revalidation; the old ownership sums are historical, not current admission.


Status: restricted G3 reported on 8 November; the approved stage is finished
and stopped. All 100 calls complete on runtime `60bfd95`. Frozen experimental
process recipes: C4/p4 at 40 cold/warm, retained C4/p16 at 55 cold (ranking
inconclusive), C10/p64 at 70 cold. Pre-G3 cleanup and diagnostics are complete.
The subsequently approved generous-bounds 75/80 probes are also reported:
both incomplete under internal memory outcomes before solving. No further
probe is approved. Intermediate verification and G4 remain gated.
Follow the
[joint-search policy](../POLICY.md#joint-parameter-searches) and
[evidence tiers](../POLICY.md#evidence-tiers-and-sample-sizing).
Preserve historical freezes, failures, charges and captures unchanged.

## Deliverable, scope and source

The broader deliverable remains cold CLI rows for **40, 50, 55, 60, 65 and
70 digits**. The latest scope limits warm tuning to **40 digits**, superseding
the earlier warm-50 proposal. **62 digits is dropped** from new measurements;
retain its historical evidence without claiming an adopted interpolation.

Tune cold anchors at **40, 55 and 70** and the warm anchor at **40**.
Cold 50/60/65 proposals are currently **interpolated, unverified**; their
verification was excluded from the restricted G3 authorization. Once separately
authorized checks pass, label them "interpolated, verified",
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
Current G3 runtime and QA source is `60bfd95`; the earlier checkpoint above
is historical. Canonical local evidence is main's ignored
`v2/benchmarks/results/a7a/g3-restricted-20261107-r2/`. It contains the
source archive, executable controller, frozen plan, per-call commands,
receipts, progress/thermals, worker timelines, exact recipes, report and hashes.
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
   From the 7 November addendum, enable worker active, queue-wait and
   publication timelines in every subsequent G2H and G3 cell, on both arms
   of each comparison, and include their summaries in reports. The measured
   G1 overhead was approximately 1%. Keep worker timelines off for G4 clean
   confirmation. Preserve already launched historical cells' actual settings;
   do not relabel them as timeline-enabled observations.
   The stage driver, combined runner and direct cell entry point reject
   G2H/G3 process calls without a declared timeline directory and G4 calls
   with timelines enabled. Serial controls have no worker timeline.
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

Measure the **independent-packet throughput reference** at N=1/4/8/10 using
matched deterministic disjoint packet work on three calibration inputs
per regime at 50 and 65 digits. Remove shared
parent admission/transport from that diagnostic, not arithmetic or validation.
Compare aggregate useful-work throughput and a local-packet/no-pool control
with the real timeline. Report scanned positions per worker CPU-second for
every allocation, local/no-pool and S0; label whole-factor S0 separately when
it includes graph/solve work. Static independent streams can suffer tail
imbalance, so this reference is not a proven upper bound. N times serial
rate is not an empirical bound, and
independent full factorizations do not model cooperative collection. Include
startup, remaining serial stages and load imbalance in any completion model.

After the CPU-snapshot screen, the first worker-side fix is family-affinity
placement, then grouped exact lease charging. Validate any parent-budget fix
already in progress, but do not treat it as the main G2H fix: the G2 workers
were 90–97% busy while spending substantially more CPU on similar scanned
work. Isolate affinity from an explicitly wider credit window; preserve paid
packet identities, exact admission and aggregate memory/CPU reservations.
Do not issue the G2H gate report before implementing and measuring affinity,
or reporting a concrete reason it cannot be measured.

### Classify and screen before constructing a joint design

The following is the initial option inventory. Each option has a one-line
reason; alternative values for one mechanism share a row. Record every
available collector option, including inactive/unsupported ones, in the
frozen manifest.
Do not silently run the old 32-row all-binary design.

| Option(s) | Class | Reason / initial screen |
| --- | --- | --- |
| Per-packet budget leases | (a) | Replace repeated shared budget calls with bounded exact reservations/reconciliation for the same work. |
| Grouped exact lease charging (`grouped_lease_charging`) | (a) | Combine two already-polled affordable marking reservations; retain scalar paid prefixes at every boundary. |
| Persistent worker state (`reuse_state`) | (a) | Avoid rebuilding identical immutable context between packets; charge retained memory. |
| Coefficient reuse (`coefficient_cache`) | (a) | Avoid recalculating identical coefficients without changing assignment coverage. |
| Family-affinity placement (`family_affinity`) | (b) | Pin original packets to reuse contiguous family state; compare a matched credit window separately from the incumbent window. |
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
| `byte_stamps` | (b) | Alternative recovery-word epoch storage, with funded wraparound clears and the existing memory reservation. Screen its interaction with packed queues as one bounded memory-layout question. |
| `wide_hits` | (c), conditionally inactive | Raise the hit-helper admission ceiling from 32 to 128 MiB; this is a capacity choice, not a different hit algorithm. |
| `wide_inverses` | (c), conditionally inactive | Raise the inverse ceiling from 8 to 32 MiB and change coexistence reservations. |
| `wide_transfers` | (c) | Raise the transfer ceiling from 32 to 128 MiB and requested bytes per base entry from 4,096 to 8,192. |
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

The three `wide_*` classifications correct the initial inventory against the
implementation before their trials, rather than reclassifying an observed
interaction. With the frozen Y0 base-count bounds through 70 digits, hit
scratch fits the default 32 MiB ceiling and the inverse table fits 8 MiB.
Record exact capacity/coexistence eligibility before constructing G3; omit
inactive levels with that bound, rather than timing equivalent recipes.
`wide_transfers` can change admitted cache capacity within these regimes.

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

### Byte-stamp complete-call remeasurement — 7 November

This explicitly requested G2H block uses source `1b5ab81` and the committed
`regime_calibration_50_{0,1,2}` and `regime_calibration_65_{0,1,2}` inputs.
Keep it separate from the earlier layout screens. The only implementation
difference is the `byte_stamps` marking component; native stamps are the
control. Both arms retain Y0, packet 16, lookahead 4, ordinary placement,
scalar leases, worker verification, ordinary sparse queues, coefficient
reuse, nursery 4 MiB, worker 512 MiB/state 64 MiB, parent 4 GiB and total
22 GiB. CPU/heap captures are excluded from clean timing; Tier A, RSS and
all-worker timelines are enabled in both arms.

| Digits | Workers | Actor wall/CPU seconds | Owner wall/CPU seconds per call |
| --- | ---: | --- | --- |
| 50 | 4 | 45 / 450 | 65 / 510 |
| 65 | 8 | 400 / 4,000 | 420 / 4,060 |

Use one fresh cold call per arm/input: twelve calls, at most 2,910 owner
seconds, 27,420 CPU-seconds and 128 MiB new capture storage, charged within
the G2H grant. The 65-digit allowance funds outstanding private CPU leases;
the earlier underfunded collection screen is not a performance loss.

Visit input indices 0, 1, 2, running the 50-digit pair then the 65-digit pair
at each index. Native goes first for 50-digit indices 0/2 and 65-digit index
1; byte goes first in the other pairs. Match factoring seeds
104729/130363/155921 by index and set `PYTHONHASHSEED` to that same declared
value before each actor starts; spawned workers inherit it. This pins an
additional source of incidental container order without attributing earlier
variation to it. Earlier unset-hash observations are not pooled.

Validate complete reconstruction against the certified oracle. Retain
timeouts as censored. For each band report all three actor-time ratios,
worker CPU/work/scanned counts, RSS and worker timelines. Form a paired
log-ratio mean, across-input standard deviation and 98 1/3% paired-t interval
with two degrees of freedom, explicitly conditional on its distribution and
independence assumptions. An interval below one supports an exploratory gain,
above one supports a loss, and an interval spanning one is inconclusive.
There is no 5% practical-effect cutoff or all-three-faster requirement in this
block. Stop after the twelve calls; extensions and adoption are not automatic.
Use the existing combined runner with the frozen per-call JSON plans and
their recorded `PYTHONHASHSEED`; retain the enclosing grant and all charges.

### Final source-gap check — replacement instruction, 7 November

For the later default-off additions, the explicitly requested Tier 1 bridge
replaces remaining new full-factor gap calls. Preserve the completed and
interrupted calls from the replaced protocol; do not pool them with prefixes.
Compare intermediate `fabfa9f` with final `f589ea3`, whose runtime is `1b5ab81`
plus documentation. The preceding three-input matched full checks remain
labeled with their original sources.

Run 512 polynomials/p16 on each of the three 50- and three 65-digit calibration
inputs, three cold repetitions per source/input: 36 calls. Hold Y0, C4/C8,
lookahead 4, native stamps, scalar leases and all other mechanisms fixed.
Use packet CPU 5 and actor wall/CPU 30/300 seconds, with a 45-second owner
allowance per call. Enclosing grant: 1,640 wall seconds, 16,400 CPU-seconds,
128 MiB storage, charged within G2H. No CPU/heap capture; Tier A, RSS and
worker active/queue-wait/publication timelines are on both arms.

Visit repetitions 0/1/2, input indices 0/1/2, then sizes 50/65. Alternate
source order by repetition + input index + regime. Pin both factoring and
interpreter hash seeds to 104729/130363/155921 by input index. Require exactly
67,109,376 scanned positions per call and identical relation signatures
across all six calls for each input.

Geometrically average paired collector-time ratios within each input, then
equally across the three independent inputs. Report cold actor ratios and
within-input A/A ranges separately. Record this gap as closed when work and
signatures match and aggregate collector ratios at both sizes lie in
[0.97, 1.03]. Do not add an all-inputs or uncertainty-interval requirement to
this explicit Tier 1 point criterion, or interpret it as default adoption.
Audit the exact source diff and list default-path touches rather than
assuming default-off selectors leave the hot path unchanged. Stop and report;
G3 still requires its own gate.

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

## Pre-G3 cleanup and diagnostics — authorized 7 November

Retain the accepted G2H implementation bundle. Remove disabled mechanisms
from default per-block/per-root loops by construction-time specialization;
preserve explicitly selected alternatives. Validate committed sources with
full tests, lint and benchmark imports. Compare accepted `d24e6bb` against
the cleaned source at 50/65 digits: input 0, three repetitions/source,
512 polynomials/p16, equal signatures and 67,109,376 scanned positions.
Report timing ratios without adding a speedup or adoption claim.

Run two independent-packet diagnostics at 65 digits on calibration inputs
0/1/2, fixed 512 polynomials/p16:

- N=1 normal versus `taskpolicy -c background`; report positions per worker
  CPU-second and wall time. Background QoS influences placement; it does not
  establish exclusive E-core residency or a pure E/P efficiency ratio.
- N=1/4/8, native versus byte stamps. Report throughput/scaling and the
  per-worker buffer payload in bytes, accounted cache ownership and observed
  RSS separately. A byte-stamp scaling gain supports a footprint hypothesis;
  it does not alone establish the cause of every scaling loss.

The cumulative pre-G3 ceiling is 3,600 wall seconds, 36,000 CPU-seconds and
1 GiB output, including failed and successful QA and cleanup. Batch owner
reservations: prefix 750 seconds (16 calls), QoS 630 (10 calls), footprint
1,350 (22 calls); the remainder funds QA/preparation. Expected measurement
runtime is about 20 minutes. No revised attempts or enlarged limits are
implicit. Preserve all failed validation and stopped batches.

### Quiet window and drift gates — all remaining timed work

The user has cleared the other Factor workstreams. Hold the real nonblocking
performance flock through cleanup; other builds, suites and benchmarks may
not overlap cells. Record the process inventory. Capture `pmset -g therm`
before/after cells and every 60 seconds during a batch. Missing readings
remain explicit; telemetry does not establish absence of throttling.

Each batch begins and ends with an A/A pair using identical source, input,
seed, work, lifetime and settings. Define drift as `max(time)/min(time)-1`.
The later user amendment makes short pre-G3 and G3 40/55 checks report-only;
continue while retaining their uncertainty. For the long G3 70-digit batch,
stop/report above 10% within either pair, or between the geometric mean times
of the start/end pairs. The pre-G3 endpoint is full external cold actor wall;
G3 cold uses that endpoint and warm uses the timed complete library call.
No automatic retry/restart after a long-run drift stop. The long-run end gate
determines that batch's validity; retain observations if it fails, without
presenting them as a quiet-window-confirmed comparison.

Production calls keep Tier A, owned RSS and active/task-wait/publication
worker timelines on both arms. Independent packets have active private
streams, startup-barrier wait and final publication spans; they have no
production task queue. G4 timelines remain off. These diagnostics do not
add full CPU/heap captures or a 75/80 call.

## G3: restricted worker-count and packet-width search — then stop

The latest scope supersedes the earlier broader joint-search proposal. Tune
only workers **{4,6,8,10}** by two valid packet fractions **{g/4,g}** at
**40 cold, 40 warm, 55 cold and 70 cold**. Actual families are g=16 at 40
(widths 4/16) and g=64 at 55/70 (widths 16/64). These are (c) tunables:
process count changes CPU economics and packet width trades scheduling cost
against imbalance. Width must divide and not exceed the actual family.
K=0; total collection/cofactor allocation is at most ten. No other flags,
cache levels, lookahead, admission order or numerical work enter this grid.

Keep Y0, lookahead 4, ordered admission, 4 MiB nursery, 64 MiB state cache,
512/256 MiB worker/runtime, 4 GiB parent, 22 GiB aggregate, worker verification,
private semaphore, packed transport and coefficient reuse. Native stamps,
ordinary queues, scalar leases, no affinity and no CPU snapshot remain the
incumbent (a)/(b) choices unless the authorized byte-footprint diagnostic
establishes a clear choice; freeze that decision before G3 outcomes. Keep
numerical base/report/product/residual limits unchanged.

The completed pre-G3 diagnostic retained native stamps. Current default hot
loops use construction-time specialization through root lifting, validated
on `60bfd95` with exact matched 50/65 prefixes and committed-only full QA
(1,575 tests, nine optional skips, lint, 292 benchmark imports). An earlier
77-call cheap G3 batch on `d1ae1ae` was closed before any 70-digit cell to
finish the missing root specialization. Its outcomes and charges are retained,
not pooled with the current-source search. All cheap rows were rerun.
Represent whole-family packets canonically as explicit width g with
`whole_family=False`; setting that flag with the same width describes the
same allocation, not another tuning dimension.

Maximum eight distinct configurations per row, no extra refinement. At
40 cold/warm and 55 cold run all eight on three calibration inputs. At 70,
run all eight on input 0, then the incumbent C8/p16 and the best distinct
complete challenger on inputs 1/2 (four calls). Prioritize by complete cold
actor wall; if input 0 has no complete challenger, stop and report rather
than selecting from partial timings. Recipes not advanced have insufficient
coverage, not demonstrated losses. Require all three inputs for a Tier 1
choice, compare against A/A noise and retain incumbent when inconclusive.
Counterbalance order and match input/hash/factoring seeds. No fresh
confirmation inputs are used for tuning.

Each row adds four identical-incumbent A/A controls (start/end pairs).
Maximum calls: 28 each at 40 cold, 40 warm and 55 cold; 16 at 70 cold,
100 calls total. Price the executable manifest before starting:

| Row | Calls | Actor wall / CPU | Owner wall / CPU per call |
| --- | ---: | ---: | ---: |
| 40 cold | 28 | 30 s / 550 s | 45 s / 610 s |
| 40 warm | 28 | 30 s / 550 s | 45 s / 610 s |
| 55 cold | 28 | 100 s / 1,000 s | 120 s / 1,060 s |
| 70 cold | 16 | 1,800 s / 18,000 s | 1,820 s / 18,060 s |

The current calls reserve 35,000 wall seconds and 352,800 CPU-seconds.
Carry the earlier batch's exact 589.638380584 wall / 5,896.383805838 charged
CPU-seconds and reserve 100 wall / 1,000 CPU-seconds for static analysis.
The combined ceiling bill is 35,689.638380584 wall / 359,696.383805838
CPU-seconds, within 10 hours / 100 core-hours. There is no transfer of unused
pre-G3 allowance. Output ceiling 16 GiB. Actual expected
runtime is approximately 3–6 hours from accepted current costs, not a target
to spend. Funding must include outstanding leases; refusals are censored and
never automatically retried. Packet CPU stays at incumbent 30 seconds.
Warm includes validated >=3-second warmup, reported outside its timing;
all cold costs include spawn/import/JIT/preparation/validation/shutdown.

Report wall versus aggregate CPU for every worker/packet recipe, completion,
noise, timelines, RSS, exact charges, exclusions and all untested interactions.
No intermediate validation, broader tuning, G4, routing adoption, merge,
push or probes are authorized by this scope cut. Stop at the G3 report.
Interpolation remains a proposed routing rule: nearest cold anchor by size,
lower anchor on ties, inherit settings and map the selected fraction to the
actual family's valid integer packet widths; label rows interpolated. No
routing table is adopted without separately authorized intermediate checks.

Freeze the exact 75/80 probe derivations in the G3 report before outcomes,
using the frozen cold 70 anchor, its implementation bundle and Y0 resolved
at the target size. Configuration resolution is not a factoring call.
Keep the separate probe allowances and approval gates below unchanged.

### Reported G3 outcome and reproduction

The current 100-call search completed, reconstructed every input and drained
owned children. The long 70-digit A/A gate passed. Short-call drift is retained
as uncertainty, especially at 55 digits; it is not treated as proof of a new
worker-count gain. The report freezes full recipes, CPU/wall tradeoffs and
worker timelines, not merely worker counts. Intermediate proposals are
50: C4/p8, 60: C4/p16 and 65: C10/p64, all unverified.

The cumulative G3 charge, including the historical batch and analysis reserve,
is 12,083.295227043 wall / 120,832.952270429 conservative CPU-seconds.
Observed CPU is recorded separately. All historical failures remain retained.
The versioned driver/runner imports are covered by committed-only QA; this
reporting update changes documentation only.

Exact historical execution from the task worktree was:

```sh
PYTHONPATH=. /opt/homebrew/bin/pypy3 -B v2/audit/a7a/pre-g3/g3-r2.py
```

The archived copy is `controller.py` beside `plan.json`; `batch-helpers.py`
and `source.tar` preserve its dependencies and committed runtime. Paths in
these immutable scripts identify the original evidence directories. Do not
rerun into them. Any new run needs a new identity, fresh output directory,
approved grant and nonblocking performance lock. No successor is authorized
by this reproduction command.

The mapped 75/80 recipe is C10/p64 from the frozen 70 anchor (actual family
64 at both target sizes), with unchanged G2H choices and target-resolved Y0.
Prospective owned bounds are 23,965,302,784 / 24,477,270,016 bytes against
23,622,320,128 bytes allowed. These are accounting refusals, not observed
physical exhaustion. No probe ran; do not silently change worker count,
capacity or recipe. Exact derivations are saved before any probe outcome.

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
applies only at 40 under the latest scope, with the frozen lifetime and at
least three seconds of validated, charged warmup. Freeze the repetition count before
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

**Reported:** exactly one cold call per size on `5526118`; both ownership
admissions pass. Both return internal QS `memory_limit` events before any
filtering/matrix solve. Neither the 18 GiB RSS guard nor wall/CPU allowance
causes the stops. Feasibility is unestablished. The precise internal cap or
allocation error is unrecorded; do not substitute a guessed attribution.
The terminal path drops live QS state, so saved outer checkpoints retain
unresolved cofactors/accounting but cannot resume the collected rows.
Progress and thermals remain durable. No further run or retuning follows
automatically. Canonical local report:
`v2/benchmarks/results/a7a/probes-20261108/REPORT.md` in main.

The actors used 1,969.676906459 / 4,794.110493167 wall seconds, respectively.
Cumulative owner charges including prior validation are 6,973.427729834 wall
and 67,091.398186075 CPU-seconds, within 40,800 / 396,000. Runtime/QA source
`5526118` passes 1,579 tests (nine optional skips), lint and 292 imports.
The same-source 80 call reuses QA. Scripts, commands, exact recipes, source
archives, progress, receipts and hashes are in each ignored probe directory.
No test or benchmark overlapped a timed actor. The contracts below record
what was approved and executed; they are not a new grant.

Historical executed commands from the task worktree:

```sh
PYTHONPATH=. /opt/homebrew/bin/pypy3 -B v2/audit/a7a/probe75-20261108/generous-owner.py --digits 75
PYTHONPATH=. /opt/homebrew/bin/pypy3 -B v2/audit/a7a/probe75-20261108/generous-owner.py --digits 80
```

Each controller is archived as `owner.py` alongside `batches.py`, committed
`source.tar`, `plan.json` and `request.json` in its probe evidence directory.
Do not execute these historical identities again. A future run needs a new
approved identity, output directory and finite grant; retain old failures.

The latest 8 November generous-bounds amendment supersedes the earlier
24 GiB probe amendment. It raises only these probes' aggregate owned-memory
admission grant to **32 GiB**, worker scratch to **1 GiB per worker** and
parent QS to **6 GiB**. Runtime remains 256 MiB per worker and both
parent/worker checkpoint caps remain 256 MiB. The bookkeeping bound sums
concurrent caps; it is not expected physical use. The separate **18 GiB live
aggregate owned-RSS stop guard** remains unchanged. On exceeding it, cancel
cooperatively, preserve progress/checkpoint and report a memory outcome.

Recomputed with portfolio coexistence: C10/p64 needs 30,410,022,912 bytes
(28.322 GiB) at 75 digits and 30,921,990,144 bytes (28.798 GiB) at 80; both
fit 34,359,738,368 bytes. Earlier 22/24 GiB calculations, failed validation
and pre-lock attempts remain historical. No factoring call was issued by
those attempts. The latest instruction approves 75 now and 80 after its
report under the conditional continuation, with one fixed recipe per size.
No fresh G4 or intermediate-size work is authorized.

75 has a **three-hour actor budget**; 80 has **eight hours**. Owner wall is
actor plus ten minutes, including setup, observation and cleanup. CPU is ten
times actor wall, with exact preparation charges included in the cumulative
stage ledger. A predecessor's unused allowance cannot fund extra candidates.
Exactly one candidate, one input, one clean cold call per probe; no repetition,
profile companion, extra seed, retuning or larger grant is implicit.

| Field | 75 digits | 80 digits |
| --- | ---: | ---: |
| Input | `regime_calibration_75_0` | `regime_calibration_80_0` |
| Collection workers / packet width | Frozen G3 cold 70 C10/p64, mapped to 75 | Same anchor, mapped to 80 |
| Actor wall | 10,800 s | 28,800 s |
| Owner wall including observation/cleanup | 11,400 s | 29,400 s |
| Aggregate CPU ceiling | 108,000 s | 288,000 s |
| Aggregate owned RAM admission | 32 GiB | 32 GiB |
| Live aggregate owned RSS stop | Over 18 GiB | Over 18 GiB |
| Worker / runtime per worker | 1 GiB / 256 MiB | 1 GiB / 256 MiB |
| Parent QS / checkpoint | 6 GiB / 256 MiB | 6 GiB / 256 MiB |
| Logical work ceiling | 96,000,000,000,000 | 96,000,000,000,000 |

Use the frozen G3 cold-70 anchor, preserving its worker count, lookahead,
ordered admission, packet fraction, implementation bundle and (a)/(b)
choices. At both targets the actual family is 64, so fraction g maps to
packet width 64. Apply only the explicitly amended memory grants. Use seed
104729 and Y0 resolved at the target size without retuning base, reporting,
product or residual limits. Record exact derivations before outcomes.
A refusal ends that probe without another recipe. Never choose an 80-digit
candidate from 80-digit outcomes. Any proposal motivated by 75 must be
reported and approved before replacing the already frozen 80 candidate.

Poll aggregate owned RSS every two seconds at both sizes, covering the
factoring parent, workers, spawn helper and external owner. This fits the
existing bounded recorder even for eight hours; it is a sampled guard,
not an instantaneous physical-memory cap. Preserve the 60-second progress
and thermal log. Stop and notify after each receipt. The latest instruction
pre-approves the conditional 80 continuation after the 75 report; no broader
search, default adoption, merge or push is implied.

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
warmup, the exact table allowances and candidate settings. Set `rss_interval=60`
and a new exclusive `progress_path` for minute-level probe logging; omit
`timeline_path` so worker packet diagnostics remain off. Disable automatic
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
and preserve the lock file. After each single probe, **stop and notify for every
outcome**. No automatic 80 call, search, confirmation or budget escalation.
