# Benchmark results and decisions

These are scoped historical observations, not runtime guarantees or a list of
new experiments. [STATUS](../STATUS.md) records current automatic routing and
its evidence limits. Use the [runner guide](README.md) for reproduction and
[pinned history](../history/README.md) for frozen source/accounting detail.
Numbers from separate comparisons are not additive portfolio gains.

## Collection and larger-band decisions

### C9a preparation

The subsequent current-source clean pair on `636484b` completes with C0/C1
raw times of 356.3/320.9s, about 9.9% less wall time for C9a on the exposed
balanced 60-digit modulus. Both timing gates fail; no adoption follows.
Collection is 237.5/234.7s, RSS 290.3/350.2MB, and each arm preserves its
historical exact work, coverage, rows and collector counts. Earlier pre-fix
times remain historical rather than matched causal estimates. Separate native
sampling at 50/60/70/80 and a 10,097-atom checked replay highlight polling,
scanning/marking, queue construction and polynomial switches. Two bounded
heap dumps show about 32.1MB incremental live GC memory, with no leak proof.
A nine-sample matched prefix diagnostic measures 0.385s with the RSS callback
versus 0.324s with a constant-false callback. That intervention preserves work
and outputs but is not a production resource-policy change or full-factor
attribution. Stop after bounded analysis; investigate observer cost before a
larger comparison. Details: `v2/audit/c9a/profiling-56.md`,
`verification-62.json` and `verification-64.json`.

The private progress-only collection path on `5b3a6cf` passes 1,070 tests
(9 skips), lint and committed-only benchmark imports. Nine alternating matched
samples after at least three seconds of warmup per mode reduce median
empty-return cost on the captured 26,862-atom / 24,196-partial store from
0.848ms for full snapshots to 0.000899ms for progress. Progress samples vary;
this is a return-cost observation, not whole-factor attribution. Matched
256-block exposed60 prefixes preserve exact serialized stores, work and stable
statistics but show no runtime benefit: C0 snapshot/progress medians are
0.324/0.346s, C1 0.394/0.412s. A separate JIT-off fixed-store retention check
returns to approximately 108MB post-GC live memory through 4,000 calls per
mode. Snapshot calls add allocation/GC work; retained store size does not grow.
This does not rule out growing-run retention issues or establish an RSS or
factoring speedup. Stop after these scoped checks; no adoption follows.
Details: `v2/audit/c9a/snapshot-50.md`, `verification-55.json`; raw replay and
retention evidence remain under `benchmarks/results/c9a/`.

The subsequent bounded 60-digit pair on `f0a5791` completes with C0/C1 raw
wall times of 411.0/487.1s. Both factors verify; both timing checks fail.
C1 uses 764.1M positions / 4,977 rows versus C0's 911.4M / 4,994, with
476.2/342.9MB process RSS. Collection remains the gap at 400.5/287.5s;
elimination is about one second in both. No adoption or independent-input
campaign follows. The diagnostic stage uses 83.6 seconds of its 180-second
allowance. On an actual mature 60-digit store, a cumulative snapshot takes
0.459ms and translates 24,196 partial IDs; the internal caller discards these
record arrays. This supports implementing a private progress-only collection
path while preserving public snapshots, not a measured full-factor gain.
Suppressing reporting counters showed no benefit. Short sentinel calibration
supports investigating rewarming under a separately identified protocol;
all original failures remain unchanged. Details: `v2/audit/c9a/followup-43.md`.

The cost follow-up prepares `cost-report` with scalar solving, rho-first
splitting and serial recovery, after removing repeated sorting, storage work
and impossible recovery. Source `e0acc17` passes 1,065 tests (9 skips), lint
and committed-only imports. Runtime/configuration equivalence is checked
against the serial diagnostic on `193b57c`. On the exposed 50-digit input,
C0/C1 take 24.2/36.2 raw seconds and both factor correctly. C1 uses 66.6M
positions / 1,678 rows versus C0's 71.0M / 1,704, with greater RSS. Serial C1
matches the shared candidate's terminal coverage/rows and charges less work;
all timing gates fail. A representative/late cofactor replay rejects shortened
rho effort: it loses eight useful decompositions and takes longer. No clean
overall speedup, larger-input acceptance or adoption follows. Stop after these
bounded checks; retain the incumbent. Local detailed evidence and source
dispositions: `v2/audit/c9a/cost-investigation-20-40.md` and `cost-review-30.md`.

The larger `7feeca4` screen completes all eight calls on four exposed inputs,
one seed each. At 50/55/60 uneven/62 digits, C0 versus C1 wall seconds are
24.6/49.7, 77.3/135.3, 236.3/375.3 and 777.9/1,059.6. Every factorization
matches independent certificates. Candidate terminal sieve coverage is lower
by about 6.1/8.2/15.3/20.9%, but raw time and process RSS are greater in every
pair. All eight timing gates fail; no clean speed ratio or population benefit
follows. Stage timers identify collection/elimination cost as the next questions.
Stop after this bounded tranche and retain the incumbent. Charges are 7,390
wall/CPU seconds; actual worker lifecycle totals are 2,769.7 wall / 2,787.6 CPU
seconds. Local details: `v2/audit/c9a/complete-screen-19-results.md`.

The repaired pair on `bb0764a` completes both 60-digit calls under the same
900-second and memory grants. C1 uses 764.0M positions and 618.4s; C0 uses
911.4M and 405.9s. Both factorizations match independent certificates. C1
needs 16.2% less terminal sieve coverage on this input but has greater raw time,
RSS and charged work. Accounted workspace is 586.6/629.4MB and RSS is
546.8/339.7MB for C1/C0. Both drift gates fail; no clean speed ratio or
population benefit follows. Stop after the pair, charging 2,010 wall/CPU
seconds; actual worker lifecycle totals are 1,032.7 wall / 1,039.4 CPU seconds.
No adoption follows. Local details: `v2/audit/c9a/feasibility-18.md`.

An actual-data audit reproduces the long-run memory stop exactly: 26,862 atoms,
2,666 rows and 579,211,587 positions. Packed factors occupy 158 slabs / 647KB,
while blanket exponent/context charges alone total 559.5MB. The `09a3a05`
repair charges shared slabs/contexts and integer digit storage explicitly.
Exact store replay with the original helper allowances reduces accounted
workspace from 1,040.1MB to 367.9MB, preserving all records and independently
verifying every closed row. The committed-only checkout passes 1,048 tests
(9 skips), lint and 241 imports. This is an accounting repair, not a measured
RAM reduction or factoring win. No new full comparison is run. Local evidence:
`v2/audit/c9a/ownership-audit-15-17.md`.

A subsequent 60-digit pair on `e4dec4a` (identical runtime to `44b1f94`)
grants 900 seconds each, a higher block ceiling and unchanged memory. C0
finds independently verified factors in 407.8s, with 4,994 closed rows. C1
stops at accounted memory after 362.4s with 2,666 rows and no matrix solve,
using 1.040GB accounted workspace / 349.5MB RSS. Both timing screens fail;
one modulus and unequal terminal coverage establish no yield or speed ratio.
The current candidate fails full-run feasibility under this grant where C0
completes. Stop after the pair. Charges are 2,010 wall/CPU seconds; actual
worker lifecycle totals are 778.5 wall / 783.6 CPU seconds. See local
`v2/audit/c9a/feasibility-14.md`; no adoption or default change follows.

The bounded-retention repair on `44b1f94` passes 1,040 unit tests (9 skips),
lint and committed-only checks. Its new four-call C10a-matched screen finishes
all 67,109,376-position prefixes under unchanged grants. At 60/65 digits the
candidate uses 210.8/192.4MB accounted workspace and 160.2/153.2MB process RSS,
without a memory stop. Its 125/111 closed relations match C0, while measured
charged work is 62.28/55.06 billion versus 7.88/7.92 billion for C0. No staged
opaque record is promoted or converted in these prefixes. All four timing
checks fail again; no complete factorization, speedup or adoption is claimed.
The run stops after four cells and analysis, charging 550 wall/CPU seconds;
actual worker lifecycle totals are 145.0 wall and 146.0 CPU seconds. Longer
yield remains uncertain. Local details: `v2/audit/c9a/feasibility-13.md`.

A follow-up used six short, equal-work prefixes (4,067,359 positions each)
at 60/65 digits, preserving the original memory grants. Full C1 retained
4,430/3,382 unique opaque cofactors: none repeated, converted or contributed
to a closed row. Disabling opaque retention preserved 6/7 closed rows and
reduced accounted ownership from 393.8/416.6MB to 185.6/252.3MB. Restoring only
the incumbent reporting bound also preserved those row counts, while reducing
charged factoring work by about 4.7/5.3 times. Cold diagnostic times are not
speed evidence; longer-horizon yield remains uncertain. The decision is to
revise opaque retention, reporting and conservative per-atom reservations
before raising limits or launching a larger campaign. The linked diagnostic
reserved another 130 wall/CPU seconds within Stage A's existing allocation;
details remain local in `v2/audit/c9a/opaque-investigation-12.md`.

The first matched C10a/C9a screen completed four calls on source `7fa31cb`
after 1,030 unit tests (9 skips), lint and committed-only checks passed.
Both arms used coupled arithmetic. At 60/65 digits C0 finished the fixed-work
prefixes with 125/111 closed relations; C1 stopped at its accounted memory caps
with 33/72 relations over shorter coverage. Opaque cofactors made up about
95–97% of C1's retained atoms. Process RSS stayed below those ownership grants.
All four canary drift checks failed; no clean timing ratio or speedup follows.
Stop before a larger campaign to investigate opaque retention and ownership
accounting. Detailed receipts and the finite decision remain local under
`v2/audit/c9a/feasibility-11.md` and `benchmarks/results/c9a/`.

The independently selectable useful-yield, admission, cycle, duplicate and
packed-record bundle is experimental. Retained complete-call evidence motivates
the comparison but does not measure this implementation. Current preparation
includes certified inputs, structural multiplier/base resolution, immutable
control and bounded runners. No new performance calls or complete-run captures
have executed; no speedup, larger-band feasibility or adoption is claimed.
See [runner commands](README.md#c9a-preparation-and-staged-execution).
The expanded combined bundle has source-level oracles and a revised two-arm
protocol. Before the C10a merge, its full unit suite passed 942 tests (9 skipped),
with lint and committed-only checks passing. Those checks establish no speedup.
The initial historical-control freeze was prepared but no feasibility worker
started before the request to merge current master. A comparison using coupled
policies in both arms requires a new matched source/configuration identity.

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
