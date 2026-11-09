# Factor status

Last updated: 8 November 2026 (America/Los_Angeles).

**H3 — Done, integrated locally on master on 8 November 2026** from
`codex/h3-resource-limits` at `e3c371e`, based on `40d2367`. Only root wall,
aggregate CPU, sampled current owned-process-tree RSS and an optional checkpoint
save cap may refuse resources. Wall/CPU default to unlimited; RSS defaults to
floor(80% of physical memory). Internal refusal ceilings and prepaid worker
leases are removed. Numerical recipes and retained capacities are unchanged.

Premerge committed-only source `e3c371e` passes 1,404 tests (9 optional skips),
full lint and all 267 benchmark imports. Fresh fixed-work QS/ECM equivalence
against `40d2367` also passes. Post-merge QA receipts, the retained incident
qualifications and branch adapters live in `v2/audit/h3_resource_limits/`.
No performance adoption or remote push is included.

| Dependent tables | Revalidation required after H3 |
| --- | --- |
| C3a | Portfolio allocation, pretest/fallback handoff and completion costs |
| B12 | Cold/pool/warm, startup/engagement and worker-selection cost tables |
| A7a | Packet/cache, worker economics and engagement/dispatch tables |
| B4a | Native startup/preparation, engagement and worker-economics tables |

Retain existing selected values until their owners validate the combined
committed source; historical measurements do not establish current costs.


The isolated A7a branch now integrates local master `755c312`, retaining its
G2H kernels and frozen G3 settings under H3 observed resources. The pre-H3
cap/diagnostic patch is preserved separately. Committed-only validation passes
1,628 tests (nine optional skips), lint and 295 benchmark imports. A final
format-only test change has identical ASTs to the full-suite source `5764f28`.
The integration also handles pool failure during submission without touching a
dead worker's stop lock. No 75/80 probe or performance adoption is added.
Historical timings require revalidation on the combined source; prior probe
reports retain their original resource model. Local integration receipts:
`v2/audit/a7a/master-integration-20261108/`.


**A7a slim G4 complete (8 November):** committed runtime `a51de80`, three
fresh balanced inputs at 55/65/70, frozen recipes versus matched strong serial
S0. At 55/65, confirmed wall speedups are 1.97x/3.10x with 2.34x/2.90x CPU.
The 65 recipe remains an interpolation. At 70, C10 completes 3/3 versus S0
2/3 within 1,800 s; no aggregate latency interval is assigned to this censored
row. Drift checks pass, all results reconstruct, and maximum sampled owned RSS
is 2.56 GiB. Process-arm inner CPU undercounts OS actor totals; reports use OS
CPU, and this experiment does not validate a binding CPU cap. Detailed local
receipts: `v2/benchmarks/results/a7a/g4-slim-20261108/` in the main repository.

**A7a reduced worker race complete:** runtime `bce0ad0`, 36/36 factoring
calls and 40 controls, with all results reconstructing. Input-0 pruning leaves
provisional process recipes C4/p8, C4/p16, C8/p16, C14/p64 and C8/p64 at
50/55/60/65/70. All four 70-digit worker counts are wall-time ties under the
measured 7.42% A/A drift; C8 wins the fewer-workers tie-break and uses less CPU
and RSS. Rejected counts generally have one input, survivors three; this is
internal calibration, not global optimality, routing adoption or fresh guardrails.
The campaign used 2.81 wall hours including its retained recorder-only launch
failure. Report and verified raw evidence remain local in the main repository:
`v2/benchmarks/results/a7a/worker-count-race-20261108-r2/`.

**Next authorized scope:** a separate directional worker extension, capped at
20. At 65/70 compare C16 with C14, continue C18/C20 only while gains exceed A/A
drift, then stop at the first loss/tie. Test best±1 only after both even
neighbours are measurably worse; skip odd tests on a flat drift band and favor
fewer workers. Lower sizes also get this refinement, with C2 when needed for
C4's lower bracket and C1/C3 if C2 strictly wins. Input-0 racing precedes
inputs 1–2 for the new 5% window. Explicit twenty-worker support and its
committed-source validation precede a separate frozen execution grant.
The running race's historical 10% rule is not rewritten. No automatic-default
adoption, merge or push; new 75/80 probes remain deferred until after B5a.


B12/F4's delivered process-ECM scope is integrated on mainline. Its separate
[ready-first startup follow-up](ROADMAP.md#b12-f4-in-progress) remains In progress:
real-spawn/arithmetic, failure/restart and full committed-checkout QA, then a
small startup comparison are pending. The normal bounded portfolio retains QS fallback.
Cold cheap defaults bypass the selector; the selected larger/warm worker
configurations remain cost-gated. ECM-specific execution now lives in the flat
`ecm/` package, with compatibility aliases for old imports, bundled Python
program/policy constants and explicit compatibility versions instead of runtime
source-file scanning. Numerical schedules and worker choices are unchanged.
Committed-only refactor source `08a9011` passes 1,360 tests (9 optional skips),
full lint and 265 benchmark imports. Local QA receipts: `v2/audit/ecm_layout/`.
That QA run released its owned performance lock; this is not a live
availability check.

Earlier admission-only verification on `c40298e` removes about 44 ms at cold
30/50/70 defaults, without importing the selector; 20-digit behavior is unchanged.
This refactor adds no new performance claim.

Matched ECM-only 60/80/100 fixtures with 18-digit smaller factors show
2.73/2.46/2.66x cold and 5.63/5.96/5.37x warm time-to-factor gains. These
exposed fixtures do not establish whole-portfolio, balanced or population-p90
gains. The [results](benchmarks/RESULTS.md#b12f4-ecm-only-selected-default-comparison-5-november-2026)
retain failures and cheap-path regressions; the [recorded exception](POLICY.md#recorded-automatic-ecm-adoption-exception--5-november-2026)
records component closure separately from those pending broader claims.
Original delivered-component evidence: `v2/audit/b12_f4/`; that task worktree
and its temporary checkouts are removed. The ready-first prototype remains
isolated on `codex/ecm-ready-start` with local `v2/audit/ecm_ready_start/` notes.

Read this first, then the [live backlog](ROADMAP.md) and [stable policy](POLICY.md).
[Results](benchmarks/RESULTS.md) hold numerical evidence; the
[completed register](ROADMAP.md#completed-register) holds delivered milestones.

## Target metric and workload

Minimum complete-factorization wall time on the Apple M4 / 24 GiB / 10-core
machine. All cores are available; report total CPU, peak memory and startup/
cancellation cost. Use equally weighted uniformly sampled 60/80/100-digit
integers, with median/p90 time and completion per band. Balanced-semiprime
completion at 60/70 is a must-not-regress constraint. Feasibility at 100 is
unestablished. Optional native backends retain a working pure-Python/PyPy
fallback. These declared objectives and censoring rules are in
[policy](POLICY.md#scope-and-standing-directions).

## Current automatic defaults

These are current committed behavior, not claims that every routed band passed
an experimental gate. Sources: [SIQS selector](qs/defaults.py),
[portfolio routing](portfolio.py), [ECM defaults](ecm/defaults.py) and
[chain eligibility](ecm/chain_options.py). Explicit settings and saved policies
win. Native means Python integers on PyPy implementing Python 3.11.

The bounded portfolio runs exact preprocessing, rho and p−1 serially, then
selects bounded process ECM or the incumbent serial path.
Omitted native configurations automatically enable SIQS at exactly 40 and
50–60 digits. Direct automatic SIQS calls use the same size selector. Explicit
`siqs=None` disables that fallback; reference `SIQSConfig()` keeps its own
fields. SSS remains opt-in. ECM workers use the selected automatic cost gate; QS
workers remain opt-in/experimental.

| Decimal digits / scope | Root resources / numerical recipe | Automatic SIQS selection | Adoption evidence and limits |
| --- | --- | --- | --- |
| Exactly 40, native | Root limits; calibrated recipe | Growing DLP; B=5,000, 5 A factors, pool 32; 2 MiB power plan; parity_two marking; window pacing; coupled scored multiplier | [C10 plan adoption](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#c10--source-derived-collection-strategies-closed-bounded-tranche). Growing versus fixed DLP was not separately confirmed at 40. Face/queue adoption below applies. |
| Exactly 50, native | Same root defaults; calibrated recipe | Growing DLP; B=30,000, 6 A factors, pool 32; 2 MiB plan; density+presieve+scan_fused+parity_two; collection pacing; coupled scored multiplier | [C9](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#c9--source-derived-large-prime-yield-and-dlp-economics-complete-and-integrated), [C10](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#c10--source-derived-collection-strategies-closed-bounded-tranche), [C10a-05](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#c10a--conditional-collection-policy-questions-open), [A1a adoption](benchmarks/RESULTS.md#a1a-integration-and-extension). A1a's 50-digit integration gate passed. |
| Exactly 51, native | Same root defaults; calibrated recipe | Same concrete bundle as 50 | [51-digit completion adoption](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#c10a--conditional-collection-policy-questions-open). A1a marking at 51 is an **unmeasured interpolation**. |
| 52–54, native | Same root defaults; calibrated recipe | Counted base 2,500 (endpoint ceiling 250,000), 7 A factors, pool 64; coupled odd/even scored multiplier (`multiplier=0`); 2 MiB plan; A1a marking plus parity_two; collection pacing | [C10a-04 adoption](benchmarks/RESULTS.md#c10a-04-larger-band-routing) and [A1a adoption](benchmarks/RESULTS.md#a1a-integration-and-extension): **unmeasured interpolated bands**. |
| Exactly 55, native | Same root defaults; calibrated recipe | Same counted bundle as 52–54 | Same [C10a-04 adoption](benchmarks/RESULTS.md#c10a-04-larger-band-routing) and [A1a adoption](benchmarks/RESULTS.md#a1a-integration-and-extension): **Adopted**; compact training unstable, A1a fresh integration Inconclusive due to drift. |
| 56–59, native | Same root defaults; calibrated recipe | Counted base 5,000, ceiling 250,000; 8 A factors, pool 64; coupled scored multiplier; 2 MiB plan; parity_two marking; collection pacing | [C10a-04 adoption](benchmarks/RESULTS.md#c10a-04-larger-band-routing): **unmeasured interpolated bands**. |
| Exactly 60, native | Same root defaults; calibrated recipe | Same counted bundle as 56–59 | [C10a-04 adoption](benchmarks/RESULTS.md#c10a-04-larger-band-routing): **Adopted**; preliminary fresh compact confirmation; original campaign incomplete. |
| Exactly 40 or 50–60, native, explicit reference recipe | Caller root limits; `recipe="reference"` | Same size/base/multiplier/growing policy; power plan off, parity_two marking, window pacing. The root guard still applies. | [C9/C10 routing scope](qs/CONFIGURATION.md#explicit-and-automatic-selection); resource variants are **unmeasured**, not additional accepted gains. |
| Other sizes, or other backends, when SIQS is selected | Root limits; reference size parameters | Reference size parameters with the coupled policies below; no automatic queue. No size-based portfolio SIQS enablement. | [Retained policy](qs/defaults.py) and [C10a-07 adoption](benchmarks/RESULTS.md#c10a-07-face-and-queue). No calibrated size preset at 41–49 or 61+; the coupled-policy promotion below still applies. |

Every fresh automatic SIQS configuration at every supported size/backend
selects `sign_policy="swap_high"`, `polynomial_policy="parity2"`,
`multiplier_policy="parity2_wide"`, `multiplier=0`,
`assignment_policy="windowed_pair"`, `root_update_policy="batch"` and composed
`parity_two` marking. The C10a-01 all-size promotion is an explicit policy
decision before fresh confirmation; the [results and limits](benchmarks/RESULTS.md#c10a-01-coupled-polynomial-defaults)
retain the 65-digit failed canary and censored 70-digit control. Other sizes
keep generic reference parameters with these selected policies;
there are no new calibrated 61+ presets or portfolio bands.
Native calibrated collectors above select a 32 MiB sparse-queue ceiling within
the existing owned allowance, with charged rebuild and exact fallback for
unsupported paths or capacity refusal. The combined C10a-07 policy is adopted.
Its 50/55/60 panel is exploratory; drift acceptance was unmet and the additional
18-call confirmation was not run. Other sizes,
backends and explicit recipes are unmeasured policy extrapolations. See the
[adoption and limitations](benchmarks/RESULTS.md#c10a-07-face-and-queue).

All calibrated bundles retain half-width 65,536, windowed paired assignments
and a 64-polynomial family prefix. Root wall and CPU are unlimited when omitted;
there is no work grant. Root RSS does not choose a different numerical recipe.
The [configuration guide](qs/CONFIGURATION.md#explicit-and-automatic-selection)
explains endpoint/product controls, explicit reference recipes and charged
base-count resolution. Historical experiment ledgers retain their own limits.

| ECM / shared portfolio scope | Numerical/cache eligibility | Current automatic choice | Adoption evidence and limits |
| --- | --- | --- | --- |
| Exactly 30 digits | Root observed resources; retained numerical eligibility | compact64: `(2000, 50000, 64)` | [C3 decision](benchmarks/RESULTS.md#c3-bounded-ecm-allocation-and-handoff-10-october-2026): **Adopted before fresh confirmation**; oversized comparison Inconclusive, source screen unstable. |
| Exactly 40 digits | Same root resources | escalating140: `(200, 7700, 10)`, `(2000, 81000, 30)`, `(10000, 554000, 100)` | Same [C3 adoption](benchmarks/RESULTS.md#c3-bounded-ecm-allocation-and-handoff-10-october-2026) and limitations. |
| Other sizes, or explicitly selected paired executors | Retained size/backend policy | Retained fixed32 numerical schedule | [C3 scope](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#calibrate-automatic-selection-and-handoff-c3-then-g1); not a measured optimum for every size/backend. |
| Supported native schedule with B1=2,000, at least 8 curves, chunk 16 and B2 below 2^64 | Compatible numerical schedule and retained program/chain caches | Reused PRAC/batch16; unsupported jobs use B4 ladder | [B3 adoption](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#b3-production-integration-acceptance-9-october-2026). Wider allocation/population confirmation remains open. |
| Bounded p−1 / portfolio p−1 | Root resources; finite numerical bounds | Direct gap recurrence, p−1 chunk 64 | [A6 adoption](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p53--improve-p1-powering-and-continuation-independently): stage evidence favorable; extended integer portfolio remains Inconclusive. |

## A7a 75/80 probes reported — 8 November 2026

Both approved, single-input cold probes are finished on `5526118`. Frozen
C10/p64 passes the amended 32 GiB ownership grant with 1 GiB worker / 256 MiB
runtime, 6 GiB parent and 256 MiB checkpoints. Both return internal QS
`memory_limit` events before filtering or a matrix solve; neither the 18 GiB
live-RSS guard nor wall/CPU grants stopped the calls. Feasibility remains
unestablished. The precise sub-limit is unrecorded; physical memory exhaustion
is not demonstrated.

Progress, thermals and outer checkpoints are saved. The existing terminal
internal-memory path clears the live QS job, so these checkpoints do not
retain collected rows. Precise refusal attribution and retaining that state
are proposed next work before another expensive call. The separate cooperative
RSS cancellation/checkpoint path passes focused tests. Committed-only QA:
1,579 tests (nine optional skips), lint and 292 imports. The same-source 80
probe reuses that passing QA. Both controllers exited after draining actors
and releasing their locks. No retry, extra candidate or adoption.

Canonical local reports in main's ignored archive:
`v2/benchmarks/results/a7a/probes-20261108/REPORT.md`.

## A7a restricted G3 reported — 8 November 2026

The approved worker-count/packet-width search completes all 100 calls on
committed runtime `60bfd95`: cold 40/55/70 and warm 40. Frozen experimental
process recipes are C4/p4 at 40 cold/warm, retained C4/p16 at 55 cold, and
C10/p64 at 70 cold. The 55-digit worker ranking is inconclusive under control
drift; unadvanced 70-digit recipes have only one input. These are calibration
choices, not a confirmed all-regime optimum or a mainline speedup claim.

Default block/root loops now use construction-time specialized kernels,
including transitive prime-power lifting. Explicit alternative mechanisms
remain available. Matched 50/65 prefixes preserve signatures and scanned
positions. Committed-only QA passes 1,575 tests (nine optional skips), lint
and 292 benchmark imports. Native stamps remain selected; the footprint and
background-QoS diagnostics do not establish a cache limit or a pure E/P ratio.

Owned actors drained and the controller exited, releasing its performance
lock. Detailed report and exact recipes:
`v2/benchmarks/results/a7a/g3-restricted-20261107-r2/REPORT.md` in the main
repository's ignored evidence archive. No G4, intermediate verification,
75/80 factoring, default adoption, merge or push was run within G3. Its static
75/80 derivations exceeded the original 22 GiB allowance. Subsequent amended
probe grants and outcomes are reported above; frozen G3 evidence is unchanged.

## A7a G2H complete — 7 November 2026

The approved overhead stage is complete on `codex/a7a-regime-tuning`.
Independent packet measurements include 1/4/8/10 workers and local execution.
Family affinity, grouped exact lease charging and byte-sized hit stamps are
implemented as explicit experimental choices. Matched full-factor checks and
the requested byte/native remeasurement do not establish a new wall-time gain;
their failures, uncertainty and original source identities remain preserved.
Current-source cold coverage completes at 40/50/55/60/65/70 digits.

The requested final source-gap check closes the remaining validation gap:
36 interleaved fixed-work calls, three inputs per representative size and
three repetitions per source/input, all match relation signatures and scanned
positions. Aggregate collector-time ratios meet the requested Tier 1 3%
criterion at both 50 and 65 digits. The diff audit lists default-path guards
and construction changes; it does not claim those paths are untouched.

Worker active, queue-wait and publication timelines are required in G2H/G3
process calls and disabled in G4 confirmation. The runner enforces this before
starting actors. Committed-only QA on `1b5ab81`: 1,572 tests, nine optional
skips, lint and benchmark imports pass. The measurement window is released.
Final closure, fixed options and diff audit:
`v2/audit/a7a/g2h-closure-20261107/REPORT.md`.
Earlier full evidence, timeline tables and repricing remain in
`v2/audit/a7a/g2h-report-20261107-r2/REPORT.md`.
Anchor tuning, fresh confirmation, probes and adoption remain pending;
G3 requires explicit approval.

## A7a G2 coverage report — 7 November 2026

The approved corrected continuation has finished all 24 coverage cells over
cold 40/50/55/60/65/70 and warm 40/50. All strong-serial/process calls complete;
literal mainline reaches its wall limit at 60 and exhausts its generic preset
at 65/70. One balanced input per size establishes feasibility and costs, not
selected recipes or fresh confirmation. Historical reference failures remain.

The representative CPU/heap capture was recovered through a scoped review
without repeating factoring. Fix `88430a5` distinguishes an intentionally
unsampled worker from a missing required stream; full captures still require
every owned process. Committed-only QA covers 1,548 tests (9 optional skips),
lint and 291 imports; two host-inspection permission errors passed when rerun
with host access on identical source. CPU accounting, sparse preparation and
budget checks supplied the overhead candidates reported above. Configuration
tuning, confirmation, probes and adoption remain pending. Local cost/resource
report: `v2/audit/a7a/g2-report-20261107/REPORT.md`. The owned lock is released.

## A7a preparation gate — 7 November 2026

G1 passes on isolated source `62854f9`: 1,546 committed-only tests (9 optional
skips), lint and 291 benchmark imports; fixed-work signatures survive observer
on/off checks, and literal-mainline cold/warm calls reconstruct correctly.
The driver binds approval to the exact stage/source/corpus and rejects missing
coverage, invalid grants and automatic retries of unfinished attempts. Fresh
certified confirmation inputs are reserved for the three/six/nine-input rule.
Observer overhead on short prefixes remains inconclusive; optional counters
stayed off in the separately approved G2 coverage. No per-regime winner,
large-probe completion or automatic-default change is established.
Local report: `v2/audit/a7a/g1-20261107/REPORT.md`.

## Work in flight

Inspected 6 November 2026; A7a results updated 8 November. These are isolated
experimental branches; their implementations are not mainline defaults.
Checkpoint summaries and explicit
remaining work live in [the roadmap's In progress section](ROADMAP.md#in-progress).
The snapshot does not assert that any process or performance window is active.

| ID | Implemented / remaining | Branch / inspected checkpoint |
| --- | --- | --- |
| [B12 / F4 ready-first](ROADMAP.md#b12-f4-in-progress) | Early dispatch/publication and partial-pool cleanup prototype; 15 mock-only tests pass. Real spawn/arithmetic, lifecycle/restart, full committed-only QA and cold-start comparison remain. | `codex/ecm-ready-start` / `636a755` |
| [F3](ROADMAP.md#f3-in-progress) | Exact polynomial continuation, P/G arithmetic, recovery/resume and opt-in profiles; full QA, selected screen, crossover evidence and integration remain. | `codex/f3-polynomial-stage2` / `66c07f3` |
| [A8 / D1](ROADMAP.md#a8-d1-in-progress) | Packed operators/lifting, block Lanczos, original correction/resume and ablations; genuine nonempty target-size fixtures, full QA, measurements and integration remain. | `codex/a8-d1-lanczos` / `9dd1968` |
| [F2](ROADMAP.md#f2-in-progress) | Suyama-11 and enhanced Edwards/setup/scalar packages; optimized edits are in flight, with full QA, feasibility/yield/cost evidence and integration pending. | `codex/f2-curve-packages` / `794e1d5` plus uncommitted edits |
| [A7a](ROADMAP.md#a7a-in-progress) | G3 and amended 75/80 probes reported; both probes incomplete under internal memory outcomes. Exact refusal attribution, state retention, intermediate verification and fresh confirmation remain. | `codex/a7a-regime-tuning`; probes/QA `5526118`, G3 `60bfd95`; ignored `probes-20261108/REPORT.md` |
| [B4a](ROADMAP.md#b4a-in-progress) | Optional scalar/native-team executor sources and rebased Python QA; native admission, complete-cost comparisons, recalibration and integration remain. | `codex/b4a-native-ecm` / `02e73d2` |

B12/F4's earlier delivered integration remains recorded, but its ready-first
startup follow-up is unfinished. The tuple-decoding microcheck is complete;
direct runtime tuple storage is the separate open B3b scope. C9a and C10a-01
remain delivered mainline scopes.

## Recent completed review

- **H3 — Done, 8 November:** locally integrated resource-limit simplification
  and 27 restored assertion-review cases. Saved-history option 2 and stage-share
  handoff are documented in the [resource API](docs/API.md). C3a, B12, A7a and
  B4a tables require revalidation on the combined committed source.
- **C9a — Done, 4 November:** delivered and integrated; policies remain opt-in.
  Complete60/70 evidence is exploratory; one exposed input per size and no
  disjoint fresh confirmation establish no population speedup.
  [Completion and evidence](ROADMAP.md#c9a-completion).
- **C10a-01 — Done, 3 November:** integrated coupled policies are adopted for
  fresh automatic SIQS before fresh confirmation. The failed 65 canary,
  censored 70 control and unmeasured all-size gains remain explicit.
  [Completion and evidence](ROADMAP.md#c10a-01-completion).
- **A1b — Done, 2 November integration:** retain defaults after the fixed/list
  screen lost at 50/55/60; no complete-factor claim.
  [Final decision](benchmarks/RESULTS.md#a1b-existing-fixed-score-review).

<a id="c10a-01-delivered-integration"></a>

## Next eligible items

The [ranked queue](ROADMAP.md#ranked-next-work--expected-benefit-and-cost) now
prioritizes current acceleration and shorter later experiment turnaround:

1. **B12/F4 ready-first:** validate real spawn/arithmetic, failure/restart and
   full committed-checkout QA, then compare cold startup — Sol / xhigh.
2. **A7a:** current-collector parallel QS preparation — Sol / xhigh.
3. **B4a:** coarse ECM C/CFFI, then selected native SIMD/threads — Sol / xhigh;
   separately selected new arithmetic proofs retain Astra / xhigh.
4. **B5a:** QS C/CFFI and ARM NEON, with A1c/B5 controls — Sol / xhigh.
5. **E6:** bundled native prime batches/reusable-context opportunity screen —
   Sol / high screen, Sol / xhigh ABI and implementation.
6. **B15:** native helper/SIMD opportunity selection across E4/A10a/A2a/C5 —
   Sol / high evidence review; Sol / xhigh for boundary/protocol design. Deliver
   at most one justified candidate or a no-go; implementation retains its owner.
7. **H2:** profile-driven overhead sweep after the B12/A7a/B4a integration
   checkpoint — Sol / xhigh. Rank time/heap costs, remove measured overhead,
   justify memory reservations and compare the combined result with the
   pre-sweep control; a no-conflict CLI/preprocessing/reconstruction subset
   may start earlier.
8. **C3a/C11:** policy/handoff screen on the settled worker cost identity —
   Sol / xhigh. Reuse accelerator receipts rather than repeating broad studies.

Then A1d, B3a, F3, broader C9b and the remaining ranked rows. F3 is the first
substantial Astra / xhigh algorithmic research choice; a concrete opportunity
can start before every conditional native port is settled.

A7a can prepare alongside B12 with pinned shared contracts; timed execution
uses the shared lock sequentially. Assess warm ECM→QS process reuse and parent
saturation. B4a/B5a need stable neighboring interfaces, not a prior worker win.
C5 is the separate low-level limb-buffer/modular-context follow-up. B4a already
owns reusable ordinary GMP objects/scratch; C5 research, specification,
implementation and handoff are outside that phase. Select C5 from remaining
allocation/dispatch/arithmetic cost against B4a's scalar baseline. Existing
Montgomery curves do not supply a residue backend, and earlier Python-reducer
losses do not settle the native question. B4a's native
thread screen centers on cold-CLI startup/heap savings against B12 processes,
using a C-side multi-curve team and a pthread alternative when OpenMP is absent.
E6/B15 at ranks 5–6 are conditional screens. B15 is not a C3a prerequisite;
its selected helper implementation does not inherit rank 6 automatically. H2
is rank 7, ahead of broader policy tuning; its removable overhead is unmeasured.
Reuse compatible H2 evidence in C11/C9b and retain their algorithmic owners.
Parallel latency gains do not imply lower CPU cost or higher batch throughput;
native-kernel gains and M4 SIMD remain unmeasured. x86-specific ports are outside
this priority increase. Small cost/A/A checks belong in each task; C9b is not
a gate on preparation. The local pressure test is
`v2/audit/research/acceleration_priority_20261104.md`.

This order is an expected-benefit judgment, not measured gains or campaign
authorization. Correctness/accounting repairs, C10a-06 and requested B14 proof
work take precedence. E1 remains Deferred and integrated tuning retains its
prerequisites. Closed C9a/C10a-01 evidence and active grants are unchanged.

## Selected ECM regimes

The bounded auto CLI/API use these defaults on supported macOS arm64/PyPy.
Explicit arithmetic/resource settings and checked resumes take precedence.
Worker counts are selected configurations, not a guarantee of optimality for
every input; the observable-cost guards and root resources still apply. Dependent C3a,
C3a, B12, A7a and B4a cost/dispatch tables require revalidation after H3.

| Total decimal digits | Numerical plan | Cold default | Retained warm preference |
| --- | --- | --- | --- |
| 20 | Existing serial plan | Serial; no selector | Serial; unmeasured worker scope |
| 30 | compact64 | Serial; no selector | Up to 12 workers, cost-gated |
| 40 | escalating140 | Up to 8 workers | Up to 14 workers |
| 50 | fixed32 | Serial; no selector | Up to 14 workers |
| 60 | escalating140 | Up to 10 workers | Up to 10 workers |
| 70 | fixed32 | Serial; no selector | Up to 10 workers |
| 80 | escalating140 | Up to 10 workers | Up to 10 workers |
| 100 | escalating140 | Up to 10 workers | Up to 10 workers |

compact64 is `(2000, 50000, 64)`; fixed32 is `(2000, 147396, 32)`.
escalating140 is `(200, 7700, 10)`, `(2000, 81000, 30)`,
`(10000, 554000, 100)`, as `(B1, B2, curves)`. Other sizes/cofactors use
the bounded existing cost model or serial fallback. Warm 30-digit easy calls
and some cold 40-digit calls regressed in the retained fixtures; hidden factor
sizes never select the runtime path. This table preserves QS numerical eligibility; H3 changes resource admission
and stage handoff independently.

## Known gaps and unresolved conflicts

- **Deferred:** E1 combined confirmation; GNFS until H1 and its coverage review.
- **Unmeasured/Inconclusive:** no general all-size speedup or broad crossover
  is established. Coupled SIQS policies apply at every size, but tuned size
  presets and automatic portfolio SIQS eligibility still end at 60. Other
  sizes/backends/grants retain the stated extrapolation and feasibility limits.
- **Open:** C9b timing variability and fresh confirmation; C10a-06 runner work;
  general R3/provenance and B14 proof scopes; current-worker/SSS comparison and
  upstream-license clarification. B12/F4's delivered ECM component remains
  recorded separately from its unfinished ready-first startup follow-up;
  A7a and B4a retain their own scopes.
- **Resolved:** C9a's historical preparation/campaign instructions are superseded
  by scoped closure. Its uncalibrated 15% canary flags remain historical warnings;
  future judgments require measured variability and uncertainty.
- **Unresolved conflicts:** none found in the checked current behavior/closure
  records. No conflicting experimental verdict was upgraded.
- **Local references:** earlier C9a preparation files are preserved under
  `v2/audit/c9a/archive-20261104/base/c9a/`, with former worktree roots recorded
  in `v2/audit/c9a/archive-20261104/relocations.json`. Detailed chronology and prior
  quantitative records remain in [pinned history](ROADMAP.md#documentation-history-and-local-evidence).
  No missing audit reference was discarded.
