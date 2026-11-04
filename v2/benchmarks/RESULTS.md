# Benchmark results and decisions

## B12/F4 complete capture and current-source headroom audit (4 November 2026)

The initial 62 CPU profiles and 83 heap dumps all had prior parsed summaries;
this audit reconciles their full histogram tails and adds 11 targeted current-
source diagnostic cells. All 189 retained raw captures (100 CPU / 89 heap) and
100 JIT-summary logs are now accounted for. No raw evidence or historical freeze
was overwritten. New cells validate identical declared curve results and drain
all their owned groups. Full capture coverage does not establish exhausted
headroom or an optimal configuration for every input.

Current-source CPU scopes cover L/H at 60/80/100 digits; P8 supplies four cases
and equivalent-manifest M0 supplies L80/H60 arithmetic. PRAC remains limited to
its accepted bounds/modulus guard. Current `consume` logical scopes occupy
roughly 24–37% of arithmetic samples, predominantly compiled JIT scopes; that
is not a removable-overhead estimate. All current trace-too-long aborts belong
to parallel parents, with none observed in arithmetic workers. A blanket global
JIT-limit increase is therefore not supported as an arithmetic optimization.

Two current M0 heap lanes show about 1.4/2.1 MB released across close and
0.1/0.3 MB live stage-two deltas at L60/H80. M0 parents retain portfolio imports;
these are cache/arithmetic snapshots, not the slim worker's exact import heap.
The full type tail does not identify a dominant missed allocation. Allocation
churn, inlined operation costs and memory bandwidth remain unresolved.

A deterministic budget prototype comparison fails its unchanged-clone timing
control. Its synthetic ratios are inconclusive for production cost; no budget
implementation or global default was changed. The retained receipts include an
initial adapter error and all unfavorable/null-control results. Remaining F4
candidates include an exact-source checked-budget comparison, immutable
stage-two index/group preparation, context lookup reuse and packet2. These
remain separate from B3a chains, B4a native/SIMD and parked E5a serde work.

New diagnostics use three measured calls after validated warmup (one at the
serial arithmetic probes), so no new speedup/p90/whole-factor claim follows.
The earlier nine-sample cold reuse result remains separately scoped. Maximum
cumulative booking is 1,610.21 wall / 11,947.36 CPU seconds within the unchanged
1,800 / 12,000 envelope; prior QA/failure charges remain. Local detailed coverage,
full per-file hashes and next decisions: `v2/audit/b12_f4/capture_coverage_v1/`.
B12/F4 remain open; no merge, defaults or E1 campaign.

## B12/F4 bounded metadata and import reuse (4 November 2026)

Opt-in resource/checkpoint v5 adds bounded decoded reuse for streamed schedules,
compact packed backing, parent manifest reuse, lightweight worker imports and
the existing reduced stage-two point kernels. Mathematical v1 curves, exact
arithmetic, stage-one acceptance guards and serial defaults remain unchanged.
Full committed arithmetic source `f9473f0` passes 1,026 tests (9 optional skips)
and full lint. Final recorder source `9065fca` passes 12 targeted checks and full lint;
25 arithmetic/execution/portfolio files are byte-identical to the full-suite source. Earlier failed QA/recorder/allocation receipts stay
charged and retained.

Nine fresh CLI samples per arm use identical declared curves and include
imports, spawn/JIT, preparation, validation and shutdown. B8 uses previous/off
settings on the same source; X8 combines slim imports, decoded 2 MiB L / 8 MiB H,
manifest 1 MiB, reduced point kernels and poll 64. M0 uses the same optimized
settings serially, a distinct equivalent-manifest control rather than incumbent
S0's legacy RNG stream.

| Case | Optimized M0 | B8 | X8 |
| --- | ---: | ---: | ---: |
|L60|1.881 s|1.130 s|0.924 s|
|L80|2.201 s|1.205 s|0.975 s|
|H100|3.118 s|1.451 s|1.168 s|
|C80|0.369 s|—|0.468 s|

X8 reduces the substantial fixed-curve CLI medians by 18–20% versus B8, using
about 2.4–2.9 times observed CPU versus optimized M0. Cheap C80 remains faster
serial. These are enriched, exposed fixtures with 20-digit smaller primes
(10 digits at C80), not uniformly sampled random integers. Only three independent
substantial moduli are covered; repeated timings are not independent inputs.
Conditional timing-block intervals support the scoped comparison, not population
p90, first-factor or complete-factorization gains.

At deeper `(50000,5000000)` bounds, one-curve common-seed warm replays improve
16–17% with compact metadata backing. Normal-bound segment/GCD changes are mostly
small; larger segments add little after backing exists. Packets and stage
splitting do not help the selected ablations. Optimized warm count candidates
still vary 10–16; close rankings and block variation remain unresolved. A 2/4/8
MiB nursery bracket yields no clear wall-time gain despite lower GC cost with
8 MiB; observed parent-plus-worker RSS is 1.26/1.55 GB in the 2/8 MiB counter
lanes, excluding controller/helper.

Eligibility clarification: PRAC requires both the accepted B1 schedule and
`10**39 <= n < 10**80`. L60/L80 use it; L100 uses the ladder. Residual profiles
put stage one near 43% and stage two near 55% of counter-phase CPU. Logical JIT
scope samples do not establish how much work-accounting cost is removable.
Warm prime-power recomputation is eliminated; scalar products and remaining
parent planning are small. Global serialization remains E5a.

The finite phase records 214 validated runtime cells and three fully charged
recorder/allocation failures; corrected designs use new identities. Including
retained QA/analysis maxima, booking is 1,310.21 wall / 11,147.36 aggregate CPU
seconds within 1,800 / 12,000; original QA and earlier phase ledgers are unchanged.
No calibrated engagement rule, defaults, merge, B12/F4 closure or E1 campaign.
Independent-input first-factor/complete-CLI comparisons and balanced guardrails
remain NOT RUN. Local report/review/handoff: `v2/audit/b12_f4/headroom_opt_v1/`.

## B12/F4 wider-worker CPU and heap screen (4 November 2026)

Runtime source `a32b9f0` uses explicit spawn processes and unchanged exact ECM.
The bounded screen covers 4/6/8/10/12/14/16 workers at 60/80/100 digits, with
20-digit smaller primes. L uses 128 curves at `(2000,200000)` and retains PRAC within the existing
40–80-digit modulus guard; L100 uses the ladder. H uses 32 curves at `(11000,1100000)` and the existing ladder. L/H share a
modulus at each size: three exposed substantial moduli, not six independent
inputs. All candidate arms preserve identical curves and results.

Cold component CLI medians below use nine fresh interpreters per arm, including
startup, factoring the declared curve set, validation and draining. M0 is the
serial equivalent-manifest arm, distinct from incumbent S0's RNG stream.

| Digits | L: M0 / P8 | H: M0 / P8 |
| --- | ---: | ---: |
| 60 | 2.138 / 1.284 s | 2.800 / 1.437 s |
| 80 | 2.549 / 1.279 s | 3.340 / 1.533 s |
| 100 | 3.025 / 1.364 s | 3.839 / 1.625 s |

P6 marginally wins L60's median, with unresolved P6/P8 timing uncertainty.
Other close count comparisons mostly remain unresolved; P8 is the leading cold
candidate, with roughly 2.3–3.1 times observed CPU versus M0. Cheap eight-curve
controls remain faster serial from a fresh CLI. Warm discovery winners vary
10–16 by size/bounds; close rankings and between-block variation prevent a
global optimum claim. H100's separate clean recheck favors P10 over P16
(0.650 versus 0.694 s mean), paying 5.26 versus 6.16 worker CPU seconds/call.

At H100, heap-observed P10/P16 parent-plus-worker RSS peaks are 0.95/1.46 GB,
excluding
controller and spawn helper. Stage two uses about 62% of measured compute CPU;
parent orchestration is about 0.5% of aggregate CPU samples. Heap-lane GC is
about 2.9%/5.1%, with roughly 26 MB live GC-managed heaps per worker. Capacity
pressure and parent JSON are not the observed dominant warm bottlenecks.
CPU/heap observation costs were measured separately; instrumented timings do
not rank candidates. JIT inlining limits operation-level attribution.

All 288 timing/diagnostic cells validate and leave no owned groups. Observed
cell totals: 942.59 wall / 5,362.38 aggregate CPU seconds; conservative cell CPU
charges 15,029.45 seconds. Observed CPU is a lower bound; unobserved helper
and teardown allowances remain charged. The separately funded ceiling is 3,600 wall / 33,000
CPU seconds with 2 GiB evidence, including previous profiles. Failed QA and a
wrong-ref validation attempt remain recorded and excluded from timing evidence.
Exact runtime-source full QA passes 1,004 tests (9 optional skips); later
analyser/adversarial tests and full lint pass in committed-only checkouts.

The opt-in contract now permits funded counts through 16 and versions resource
checkpoints as v3, repairing durable position-key hashing and parent CPU room.
Mathematical v1 identities and serial defaults remain unchanged. First-factor,
complete-portfolio, balanced-tail and fresh independent-input confirmation are
NOT RUN. No calibrated engagement policy or default adoption follows. Detailed
local report, all count tables and handoff: `v2/audit/b12_f4/counts16_v1/REPORT.md`;
preceding CPU/heap phase: `v2/audit/b12_f4/profiles-v1/REPORT.md`.
B12/F4 and E1 retain their open/deferred status.

## B12/F4 larger-work headroom screen (4 November 2026)

Runtime candidate `8a2bff7` and current incumbent `a2af71b` were pinned in
separate committed checkouts. These are exposed/enriched exploratory inputs,
with one observation per cell, explicit spawn and identical candidate manifests.
All 18 revised cells completed, including fresh CLI startup/output validation
and separate retained contexts with at least three validated warmup seconds.

| Input / fixed work | Cold M0 serial | Cold P4 | Warm M0 | Warm P4 |
| --- | ---: | ---: | ---: | ---: |
| 60 digits, smaller factor 18 digits; 128 curves, B1=2,000/B2=200,000 | 1.910 s | 1.106 s | 1.588 s | 0.508 s |
| 80 digits, smaller factor 20 digits; 32 curves, B1=11,000/B2=1,100,000 | 3.123 s | 1.414 s | 2.701 s | 0.811 s |

Cold incumbent S0 observations were 1.896/3.011 s, with different curve streams.
M0/S1/P2/P4 preserve all curve identities, splits, phases and chain counts.
The 60-digit manifest found a divisor on one curve; the 80-digit manifest found
none. Fixed-work differences do not establish first-factor or complete-factor latency.
Full-tier PRAC remains eligible at B1=2,000; every deeper-bound arm uses the
existing ladder path. The missing cheap warm M0 observation was 0.0591 s.

P4 cold observed CPU was 2.04/1.62 times M0. Warm child CPU counters indicate
about 3.9 active CPU equivalents; parent CPU was small. Helper CPU remains
incompletely observed and conservatively charged. Largest summed worker HWM
was 344 MB and parent HWM 89 MB, not simultaneous aggregate RSS. These results
justify a bounded first-factor/complete-call follow-up, without selecting a
production threshold or proving median/p90 or portfolio gains. Profiling and
heap dumps were not needed for this headroom decision and were not run.

Two initial refusals remain recorded: insufficient checkpoint ownership at
128 curves, then CPU leases filling the envelope without parent headroom.
The revised settings fund a 2 MiB checkpoint and cap individual curves at two
CPU seconds. These are scoped diagnostic settings, not new defaults. All
revisions and charges remain within 600 wall/2,400 aggregate CPU seconds.
Local freezes, detailed costs and limitations: `v2/audit/b12_f4/HEADROOM_REPORT.md`.
B12/F4 remain open; no merge, default adoption or E1 campaign accompanies this.

## B12/F4 runtime QA and initial process screen (3 November 2026)

Runtime source `8a2bff7` passes 987 tests (9 optional skips), full lint and
committed-only validation. The exposed 60/80-digit, ten-digit-smaller-factor
screen uses eight curves at `(2000,147396)`. All cells complete; the 60-digit M0/S1/P2/P4 arms
retain identical curves/results and full-tier PRAC/program eligibility.

One corrected cold CLI observation at 60 digits: S0 0.294 s, equivalent-manifest
M0 0.298 s, S1 0.507 s, P2 0.456 s, P4 0.457 s. S0 uses different serial seeds.
Separate warmed library calls were 0.0606/0.0374/0.0248 s for S1/P2/P4 after
validated retained-context warmup. Pure warm dispatch overhead and warm M0
remain unmeasured. This short prefix does not justify cold worker engagement;
warm scaling warrants targeted calibration, without adopting a count/threshold.
These single observations establish no speedup, complete-factorization benefit,
median/p90 or balanced-tail acceptance. Production defaults remain unchanged.

The first cold endpoint omitted external decoding/validation. Its freeze and
charges remain; a tested correction and separate five-arm cold freeze include
those costs. Earlier QA fixture/pin failures are retained. Historical A7 arms
remain intact through new chained source manifests. Local evidence and next
steps: `v2/audit/b12_f4/RUNTIME_REPORT.md` and `RUNTIME_HANDOFF.md`. No B–E or E1
campaign, wider-worker run, merge or push accompanies this screen.


These are scoped historical observations, not runtime guarantees or a list of
new experiments. [STATUS](../STATUS.md) records current automatic routing and
its evidence limits. Use the [runner guide](README.md) for reproduction and
[pinned history](../history/README.md) for frozen source/accounting detail.
Numbers from separate comparisons are not additive portfolio gains.

## Collection and larger-band decisions

### C10a-01 coupled polynomial defaults

**Integrated and enabled for fresh automatic SIQS at every size on 3 November
2026 by explicit user direction, before fresh confirmation.** The selected
bundle is parity2 polynomials, parity2_wide coupled ranking, windowed_pair A
selection, batch root updates and parity_two marking composed with existing
kernels. Standard/even-discriminant fallback, exact extraction, explicit
controls and legacy saved policies remain supported. Existing size parameters,
grants and portfolio bands remain unchanged; no calibrated 61+ route is added.
The [configuration guide](../qs/CONFIGURATION.md#coupled-polynomial-policies)
documents the API and checkpoint semantics.

One exposed balanced semiprime per size was used for complete-call feasibility,
with explicit larger-size configurations. These observations do not establish
a population speedup or performance of every automatic preset.

| Digits / comparison | Control wall seconds | Candidate wall seconds | Outcome and limitation |
| --- | ---: | ---: | --- |
| 60 | 432.656 | 204.734 | Both complete; exploratory. |
| 62 | 935.933 | 708.487 | Both complete; exploratory. |
| 65, original | 1,572.906 | 2,310.505 | Both complete; Low Power Mode enabled during candidate. |
| 65, matched rerun | 2,007.330 | 1,703.859 | Both complete; candidate raw cost 15.1% lower, control native canary drift fails. |
| 70, four-hour restart | User-stopped at 12,705.021 elapsed | 7,603.072 | Candidate complete; control censored. No successful-factor speed ratio. |

Returned factors reconstruct every completed input and match independent
corpus certificates. Runtime certainty remains probable-prime, not upgraded
by the corpus proofs. The 65-digit rerun reproduces every non-timing SIQS
statistic in both arms, including work, positions, families and matrix counts.
The original candidate's SIQS sentinel slows by 2.504 times; the system log
places the power-state change inside that run. No arithmetic regression was
identified. Rerun SIQS canaries pass, but control native post/pre is 0.695,
outside the frozen 15% drift gate. That failed check remains a limitation.
Native-only warmup diagnostics did not resolve the proxy's variability.

V4 short calibration comprises 38 cells, with 2,910 reserved wall/CPU seconds
and 414,862,762,264 settled work units. Earlier calibration, a failed revised
protocol, interrupted 70-digit launches and charged stops remain separate.
The completed 65 follow-up settles 8,930 reserved seconds and 841,500,280,597
work units; recorded numerical phases consume 4,339.326 wall and 4,337.977 CPU
seconds. An additional 2,144.422-second administrative GUI stall is recorded
separately; its CPU is unknown. No factor job was active during the stall.
Matched operational tiers were 1 GiB owned / 3 GiB RSS at 60, 2 / 5 GiB at
62/65 and 4 / 8 GiB at 70. They are grants, not digit-only memory theorems.

Fresh independent-input confirmation, the lower 50/55 supplement and general
all-size performance remain unmeasured. Censoring and the failed canary forbid
a confirmed speed claim. The all-size default decision is a policy promotion
with those evidence exceptions, not a relabeling of the experimental gates.
No further campaign or E1 execution is implied.

Candidate source `4b9a616`, literal control
`0076e256871722d52d849ab344152e0a7535c2ce`, and terminal handoff
`d37b809104c9bfd727c0ba0994121c9060a322ec` remain in Git history. Detailed
research, raw results, resource derivations and the verified archive chain are
local under `v2/audit/archives/c10a01-20261103/`; v12 SHA-256 is
`9d6aed9749b21dd6e177d626423862ab4de95a7897e3ea2240be308802ae2d89`.

### A1b existing fixed-score review

**Retain integrated SIQS defaults and fixed-score opt-ins.** The completed
[review and gap repair](https://github.com/nishanth17/factor/blob/23b71fa723b7982f1b1ceb18c426036e5e39290d/v2/studies/a1b_gap_closure.md) has 30 validated prefix
calls on nine exposed independent moduli. Fixed/list loses every pair;
input-level geometric wall/CPU ratios are 1.273/1.273 at 50 digits,
1.468/1.455 at 55 and 1.830/1.830 at 60. It bypasses sparse queues, so these
are bundled scoring/scheduling findings, not isolated arithmetic effects.
No band advances to a complete-factor pilot or fresh confirmation.

The [numerical report](https://github.com/nishanth17/factor/blob/23b71fa723b7982f1b1ceb18c426036e5e39290d/v2/studies/a1b_final_results.json) preserves the original
24 observations plus six previously unvisited 55-digit calls, source identities,
failures, missing historical clocks and explicitly revised phase accounts.
Stage A charges are 1465/1470 reserved wall/CPU seconds. Both evidence archives
are verified. Production arithmetic and defaults are unchanged. The study is
terminal; do not rerun historical timing to replace outcomes.

### C10a-07 face and queue

[Combined face/queue defaults are adopted](https://github.com/nishanth17/factor/blob/23b71fa723b7982f1b1ceb18c426036e5e39290d/v2/studies/c10a07_default_promotion.md) after production integration and committed
correctness/QA checks. The timing evidence below remains exploratory; the
optional 18-call confirmation was not run. Automatic SIQS jobs use the swapped
coordinate policy and exact queues where eligible, with bounded exact fallback.
No new 61–80 presets or portfolio stages are enabled.

The [completed C10a-07 panel](https://github.com/nishanth17/factor/blob/23b71fa723b7982f1b1ceb18c426036e5e39290d/v2/studies/c10a07_panel_results.md) retains 42 timed
calls at 45/50/55/60, preserving every historical freeze. Nine fresh active inputs
complete in every arm. Observed equal-band wall cost falls 22.4% for queue and
23.1% for face+queue, with only 0.86% incremental combined benefit. Missing parent
end drift keeps these exploratory; active continuation drift passes, while 45
is unresolved and its drift fails. The original report nominated queue alone;
the later combined adoption is a separate policy decision.
The [numerical summary](https://github.com/nishanth17/factor/blob/23b71fa723b7982f1b1ceb18c426036e5e39290d/v2/studies/c10a07_panel_results.json) records exact costs,
uncertainty, resource accounts and the verified 26.3 MB archive. Committed-only
857-test QA (nine skipped), lint and 215 benchmark imports pass.

The [C10a-07 recovery report](https://github.com/nishanth17/factor/blob/23b71fa723b7982f1b1ceb18c426036e5e39290d/v2/studies/c10a07_recovery_results.md) retains the
stopped historical quota/face diagnostics and the completed 24.7-minute cold
historical control. Its new current-mainline calibration uses short exact-work
prefixes, separately measured observer overhead and directly admitted face arms.
At preparation it nominated the existing A1a fused collector at 60; local norm
bounds and an 8 MiB power cache did not advance. Its original screen remains a
preserved preparation artifact, superseded for execution by the completed panel
above. It is not an instruction to run another comparison.

The separate [source-gap follow-up](https://github.com/nishanth17/factor/blob/23b71fa723b7982f1b1ceb18c426036e5e39290d/v2/studies/c10a07_source_gap_results.md)
adapts exact sparse block scheduling from the source comparison. Nine matched
prefix pairs on three exposed inputs show 40.4% lower median prefix cost against
the existing fused collector, with identical exact output/coverage counters and
an explicit separate helper reserve. This is not a complete-factoring gain.
Its original six-call complete-cost screen remains frozen and unstarted; the
consolidated panel superseded its execution. Do not add it to the completed
campaign. The later production adoption is recorded above; formal clean confirmation
remains unmet.

The [consolidated collection assessment](https://github.com/nishanth17/factor/blob/23b71fa723b7982f1b1ceb18c426036e5e39290d/v2/studies/c10a07_research_assessment.md)
records 42 source/code/evidence dispositions and the resulting Phase 3 order.
The panel's original recommendation is superseded by the combined adoption
above. The optional [18-call design](https://github.com/nishanth17/factor/blob/23b71fa723b7982f1b1ceb18c426036e5e39290d/v2/studies/c10a07_recovery_followup.md)
is retained for future validation; it was not run, is not required for closure
and schedules no new campaign. Keep normalized polynomial forms,
cross-B power roots and residual-yield policies distinct; avoid reimplementing
existing polling, fixed scores or batch arithmetic. Those ideas remain backlog
and were outside the completed panel's implementation scope. The source
assessment itself establishes no new performance claim.

### A1a integration and extension

The [A1a integration record](https://github.com/nishanth17/factor/blob/23b71fa723b7982f1b1ceb18c426036e5e39290d/v2/studies/a1a_frozen_results.md) compares the
selected density+presieve+scan_fused port with committed pre-integration
mainline at exactly 50 and 55 digits. Preparation includes an immutable owned
baseline, independent arithmetic/resume checks and a fresh-input freeze;
the [frozen execution report](https://github.com/nishanth17/factor/blob/23b71fa723b7982f1b1ceb18c426036e5e39290d/v2/studies/a1a_frozen_results.md) accepts a scoped
6.008% complete-cost reduction at 50 digits (99.5% interval: 4.373–7.795%
lower). The 55-digit result remains inconclusive after its one full stability
extension because repeated-cell drift persists. Subsequent
[adoption](https://github.com/nishanth17/factor/blob/23b71fa723b7982f1b1ceb18c426036e5e39290d/v2/studies/a1a_default_promotion.md) selects the
integration at native 50–55 digits, inclusive, at 1 GiB owned memory; 51–54
are unmeasured extrapolations and 55 retains its inconclusive timing label.
The completed isolated study's 5.945%/3.214% cost reductions do not
establish additive gains. Forty digits retains control; clean
60-digit results remain inconclusive and the 65-digit pilot showed no
completion within its allowance. Replay that original study from `d61df6c`,
preserving its controls and source hashes. E1 remains deferred.

The [bounded A1a extension](https://github.com/nishanth17/factor/blob/23b71fa723b7982f1b1ceb18c426036e5e39290d/v2/studies/a1a_frozen_results.md) separately compares
selected-power periods plus four-stride per-root marking with the exact prior
A1a implementation. Its original integration control, corpus and receipts
remain unchanged. Its separately excluded fresh panel completed all calls,
but the combined extension cost 5.285% more at 50 digits and 2.994% more at
55 digits (99.5% intervals: 4.008–6.245% and 2.432–3.501% higher). Both bands
retain the prior integrated control. This tests the combined extension;
individual-component effects and a direct extended-versus-mainline comparison
remain unmeasured. The extension stays opt-in; E1 remains deferred.

### C10a-04 larger-band routing

The separate [C10a-04 bounded follow-up](https://github.com/nishanth17/factor/blob/23b71fa723b7982f1b1ceb18c426036e5e39290d/v2/studies/c10a04_short.md)
compares compact with control at 60 digits, using fresh inputs and a ten-hour
total deadline. Fixed early trend stops are preliminary below nine samples;
they do not satisfy the ordinary accepted timing gate or promote defaults.
The original multi-band campaign remains incomplete and its evidence intact.
The [completed preliminary result](https://github.com/nishanth17/factor/blob/23b71fa723b7982f1b1ceb18c426036e5e39290d/v2/studies/c10a04_short_results.md) stopped
after six pairs: compact completed 6/6 versus 5/6 and used 36.16% less total
capped time, with 31.72–40.08% saving per input. The subsequent
[promotion](https://github.com/nishanth17/factor/blob/23b71fa723b7982f1b1ceb18c426036e5e39290d/v2/studies/c10a04_default_promotion.md) enables
50–60-digit routing, preserving 50/51 settings and documenting unstable 55
training and unmeasured intermediate sizes.
The [65-digit feasibility spike](https://github.com/nishanth17/factor/blob/23b71fa723b7982f1b1ceb18c426036e5e39290d/v2/studies/c10a04_spike65_results.md) timed
out in both arms at 1800 seconds; no routing above 60 is promoted.

## Collection pacing, DLP and integration

The C10a-05 report separates polynomial windows from
collection exhaustion. A single fixed challenger was compared against the
current C9/C10 plan: 63 training pairs and 27 fresh pairs, 180 validated
complete calls. Fresh 50-digit cost fell 25.99% (95% interval 22.87–28.24%)
over nine matched pairs on three independent certified inputs. Only that
native 1 GiB bundle adopts collection pacing; 40 digits remains on window
pacing after its unstable eighteen-pair extension. The report includes extra
collection, downstream stages, setup, work, storage and charged replay costs.
The 6244.182-second inclusive study charge includes prepaid feasibility/QA;
cold probes and censored profiles are separate. Reproduce the frozen runner
at `96e2aad`; current integration intentionally changes its strict source pin.
Other C10a items, C9a feedback and E1 remain open. These percentages are not
additive with earlier C1/C9/C10 studies.

The larger-band follow-up checks
50+ routing. Nine matched training groups on three 51-digit
products gave 9/9 complete candidates for C9/window, C10/window and
collection pacing, versus 0/9 for the actual 51-digit default. The direct
C10/window-to-collection complete-cost reduction was 28.19% in training
(95% input-cluster interval 27.43–29.07%). After committing selection,
nine fresh matched calls on three newly generated products were 9/9 complete
with the selected bundle versus 0/9 with the previous default. This is a
completion gain, not a speedup against successful baseline factoring.
Inclusive wall/CPU charge was 12,199.720/12,188.535 s, including 3,500 s
prepaid feasibility/QA; a six-group training lease resumed without resetting
charges. The 55-digit one-number pilot finished under a longer 1,200 s grant,
while the 60-digit pilot was censored at that grant. At this decision, only exact 51-digit native routing advanced; the later
C10a-04 extension above covers 52–60 with its stated limitations.

The C10 tranche 1 report records an isolated
2 MiB prime-power plan improvement on fresh complete factoring: 21.83% at
40 digits (95% interval 12.47–24.88%) and 17.79% at 50 digits (16.11–19.48%),
with nine matched validated calls per arm and three independent certified
inputs in each band. The 30-digit control stays selected. The tested 16-family
actual-A ordering did not advance, while normalized polynomials, independent
family diversity and per-job threshold feedback retain explicit C10a gates.
The 50-digit training extension stopped at 11/18 groups under its frozen
whole-group reservation; the fresh 9/9 result, not that partial training
interval, supports the scoped decision. Cold starts, profiles, work and
memory are reported separately. The C10 default decision
routes the exact 40/50-digit native bundles and plan automatically at the
confirmed memory allowance. The six conditional follow-ups were assigned to C10a; this
closes the delivered C10 tranche without claiming the whole size-policy space.
The subsequent C10a-05 decision above resolves only fixed collection pacing.

The subsequent C9 default promotion makes its
confirmed native 50-digit bundle automatic, preserving
explicit limits/settings and saved configurations. Its finite omitted grant
uses the confirmed resource window; this adds no new timing evidence.
The later C10 selector composes its accepted plan with that 50-digit bundle;
the C9 control and its default-promotion report remain historical evidence.

The C9 integration record documents C3/C9 resume
compatibility and preserved study freezes. Remaining literature transfers are
phased into roadmap C9a; E1 and C10a remain separate.

C9.2's adaptive-storage acceptance confirms an
explicit balanced 50-digit bundle: 38.10% less complete-call time than calibrated
SLP and 53.08% less than fixed DLP, with 9/9 valid completions in all arms.
Peak owned reservation rises to 326 MiB (SLP 167, fixed DLP 198). Three
independent inputs and three seeds per input support this narrow scope;
30/40 historical recommendations, C10a calibration and E1 combined acceptance
stay separate from the C9 result.

## C3 continuing optimization

Additional primary research and
diagnostic findings expand the search beyond quick8.
A frozen 84-call instrumented pilot validates all outcomes and finds two
extra-curve hits that avoid SIQS. Matched profiles identify redundant
reservation checks; their candidate reduction passes 600 tests and lint.
The offline stopping fit subsequently completes all
2,052 validated instrumented calibration calls in 2,893.780 active seconds.
It selects 64 curves at 2000/50000 for 30 digits and 64 at 2000/147396 for
40 digits, with nineteen independent root-ECM entrants per band. The table
and every calibration-capture hash are committed before comparison. These
are training projections; no warmed policy promotion follows from the fit.

The complete-call procedure froze
57 revealed subjects, nine seeds and three arms. It stopped after review of
the 513-group training scope: 72 complete matched groups (216 valid
calls) are preserved, with no comparison decision. Its interrupted group grant
remains charged; the protocol is permanently inconclusive and cannot resume.
See the scope decision.

The replacement upstream transfer screen
uses six revealed subjects, nine seeds and six arms: 54 groups, 324 calls.
It tests the protected source presets, retaining the
exact SIQS settings and admission bounds. Promising unstable results may extend
to 18/27 samples only within the separately fixed phase cap. This screening
cannot establish a broad default. A selected candidate still needs disjoint
fresh confirmation. See model limits and the
prepared certificate mechanism.
The initial source-screen results complete
all 324 calls. `escalating140` and `ladder115` show 44.44%/47.39% lower
training CPU cost than fixed32, but uncertainty remains wide. Time concerns
led to deferring the prescribed extension and fresh confirmation, selecting no
policy and retaining defaults. Final committed-only checks pass 653 tests, lint,
154 benchmark imports and 35 certified corpora.
The subsequent default promotion selects
compact64 at 30 digits and the escalating prefix at 40 digits,
subject to existing workspace admission. This separate decision overrides the
earlier defer decision; it does not alter the frozen assessment's null selection,
pass its stability gate, or claim a measured hybrid gain. No new timing search
or extension was initiated. Committed-only verification of the new defaults
passes 663 tests, lint, 154 imports, 35 certified corpora and 57 training inputs.
No round-two fresh inputs exist. Frozen earlier experiments use their source
commits; historical manifests are not repinned to newer production sources.

## C3 bounded ECM allocation and handoff (10 October 2026)

Research, protocol and
acceptance document the new opt-in cumulative pretest,
finite campaign and protected SIQS/SSS transition. Six bundles train on ten
historical inputs; the selected policy is committed before 15 fresh certified
inputs are generated. All 972 fresh calls validate and complete under matched
10^13-work, 5/30-second, 288-MiB native service allowances.

The initial bounded decision retained numerical defaults. Training-selected quick8 costs 51.56% more on the
smaller fresh cohort (cost-change 95% interval −5.36% to +202.05%); uneven
controls regress and seven cells remain unstable after the frozen extensions.
The larger band retains control at training. Explicit reservations do not make
the default 2M grant sufficient: four of five smaller service-policy probes
refuse unfunded fallback floors. Charged active/terminal restore checks pass.

Working/committed-only checks pass 593 tests, full lint, 144 benchmark imports
and 34 certified corpora. Required inputs are versioned; raw receipts remain
ignored and preserved locally. Reproduction commands and all memory, seed,
small-population and source-identity limits are in the acceptance report.
That tranche left broad C3 calibration, G1 and E1 open. The later
numerical promotion is recorded in the preceding C3 section.

## C10a-04 historical freeze

The original 55/60/65-digit campaign remains paused and incomplete. Its
frozen sources and protocols are retained in Git history at `225634b`, with
detailed research and every attempt preserved locally. The separate 60-digit
result informed the 50–60 default promotion; the 65-digit spike
established no completed-factor evidence at its 1800-second grant. See the
linked result and promotion summaries above. E1 remains deferred.

## Historical v1 comparisons


Historical M9 measurements on PyPy 7.3.23 / Python 3.11.15, macOS arm64:

| Complete-factorization batch | Emulated v1 | Repaired v2 (M9) | Elapsed-time reduction |
| --- | ---: | ---: | ---: |
| Original: 5 inputs × 5 seeds | 1.155 ms | 1.078 ms | 6.6% |
| Separate control: 56 inputs × 5 seeds | 21.091 ms | 16.994 ms | 19.4% |

Times cover entire batches, with at least three seconds of validated warmup and
15 samples. Every timed result has checked factors, multiplicities and exact
reconstruction. v1 uses a `lib2to3` syntax/integer-division adapter on the same
PyPy interpreter; this is not a native Python 2 comparison. These measurements
predate later v2 work and do not describe current performance on all inputs.
Invalid baseline outputs receive no speed ratio.

Later v2-to-v2 comparisons show where the newer work helps:

- A bounded 20-digit Phase 2 confirmation raised completion from **70.8% to
  91.1%** under matched 50 ms operation caps. Cohort cost includes unfinished
  outcomes and is not a successful-factorization speed ratio.
- P3.6.1 held-out complete 30-digit cohorts reduced elapsed time by **70.7%
  for SIQS, 77.2% for SSS and 78.6% for SSSf** against their frozen v2 controls;
  all 72 timed attempts per arm completed under the declared configurations.
- SIQS dependency cadence 32 reduced a later complete 30-digit cohort by
  **38.7%** against R3 cadence 1. The isolated R3 reference retained cadence 1;
  later size-scoped defaults have separate adoption evidence. Parallel workers have no demonstrated universal advantage.

See the [runner guide](README.md) and linked studies for protocols, commands,
uncertainty, rejected experiments and size/budget limitations. A recorded
balanced 50-digit success took **19 min 43 s**; one success does not establish
broad coverage. There is no supported “50–60 digits within a minute” guarantee.


The [earlier API guide](https://github.com/nishanth17/factor/blob/24ff153d26ada2ec4beeead07000d183e94fd355/v2/README.md) retains the complete stage/kernel
implementation narratives and scoped losses. Frozen studies and their source
identities remain unchanged.
