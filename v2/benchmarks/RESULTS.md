# Benchmark results and decisions


## B12/F4 cold serial repair and component closure (5 November 2026)

The normal production admission path now skips catalogue/source validation for
its measured cheap cold plans: compact64 at 30 digits and fixed32 at 50/70.
No ECM or QS arithmetic is changed. Custom schedules, warm preference and
checked resumes retain ordinary selection. The usual portfolio keeps QS
fallback; the disabled-fallback ECM-only study remains a benchmark endpoint.

Source `c40298e` integrates current master `0efaba1e` and passes committed-only
QA: 1,357 tests (9 optional skips), full lint and 265 benchmark imports. The
initial import inventory incorrectly included a historical source snapshot
loaded as data; that harness failure and its full reservation remain retained.
The passing full suite is reused; only the corrected import check is rerun.

An admission-only check compares `ec4ee79` with `c40298e`, nine fresh processes
per arm/size after separately validated warmup. It never invokes ECM arithmetic
or fallback. Median cold production admission falls from 44.35/43.93/43.95 ms
to 0.033/0.035/0.035 ms at 30/50/70 digits, respectively. New calls construct
no selector and do not import its module. Whole admission-probe processes fall
from about 192–195 ms to 148–149 ms. Twenty-digit admission remains about
0.03 ms. These are admission costs, not new whole-factorization speed ratios.

Explicit direction closes the delivered B12/F4 component and authorizes its
mainline integration with the selected cold/warm regimes. The ECM-only results
below retain their source identities, quota failures and remaining regressions;
none is upgraded into a general portfolio, balanced-completion or p90 claim.
Native threads, QS workers, serialization and combined E1 work keep their own
scopes. Detailed evidence remains local under `v2/audit/b12_f4/`, including
`closeout_serial_bypass_v1/` and `closeout_serial_bypass_v2/`.

## B12/F4 ECM-only selected-default comparison (5 November 2026)

Source `ec4ee79` compares the unchanged selected cold/warm policy with actual
master `0efaba1e` serial ECM, using identical selected bounds and SEARCHv1
curve manifests. The control adapts curve seeds; it preserves incumbent
arithmetic and full-schedule program/PRAC eligibility. It is not stock-master
CLI behavior. No preprocessing, QS, SSS or balanced-semiprime run is included.

Fourteen previously exposed unbalanced moduli cover 30/40/50/60/70/80/100
digits (two per band). Nine fresh-interpreter cold samples and nine warm
calls after at least three validated seconds per retained context give 504
measured calls. All independent divisor, manifest and process-group cleanup
checks pass. The 18-digit-smaller-factor inputs at the three primary sizes
find factors in all samples:

| Total digits | Cold master / selected, seconds | Cold ratio | Warm master / selected, seconds | Warm ratio |
| --- | --- | --- | --- | --- |
| 60 | 3.0890 / 1.1322 | 2.73x | 2.7191 / 0.4832 | 5.63x |
| 80 | 2.6887 / 1.0917 | 2.46x | 2.3385 / 0.3925 | 5.96x |
| 100 | 2.9498 / 1.1099 | 2.66x | 2.6209 / 0.4877 | 5.37x |

All outcomes remain in the local tables. Both arms exhaust the same quotas
without a factor on the larger-factor 50/60/80/100 inputs and both 70-digit
inputs. Their faster search completion is not a time-to-factor result.
Cold 30/50/70 calls stay serial and incur selector overhead; both 30-digit
fixtures regress, as do the 50/70 fixtures. The cheap 40/15 fixture regresses
cold (0.2450 to 0.4569 seconds), while 40/19 improves. Warm 30/10 also
regresses (0.0041 to 0.0089 seconds). Thus neither universal improvement nor
optimality for every modulus follows from the prior grid selection.

Candidate engagement is 72/126 cold and 126/126 warm measured calls. At the
three successful primary-size fixtures, observed cold aggregate CPU rises
from 2.69-3.09 to 7.88-8.48 seconds; published warm per-call CPU rises from
2.33-2.72 to 3.59-4.47 seconds. Their cold pool handshake costs 0.185-0.189
seconds; observed drain tails are 3-12 ms and final close about 5 ms. Total
cold latency includes interpreter/import/JIT, preparation, validation and close.
Warm startup/warmup and whole-context CPU remain separate. Cancellation CPU
is included in total CPU but is not independently attributed after first factor.
Largest summed worker RSS high-water marks are 1.06 GiB, not a simultaneous
RSS measurement; maximum selected owned reservation is 2.26 GiB under 8 GiB.

The completed revision takes 846.147 wall / 2081.735 observed measurement-cell
CPU seconds. The original one-hour family, including the stopped different-
configuration diagnostic and one failed harness launch, totals 2245.864 wall /
15,858.490 conservative CPU charges under 3600/72,000. Final single-threaded lint passes in a further 3.407 wall /
3.297 CPU seconds, still within those ceilings. The harness launch
failed before ECM imports due to an inherited CPU limit; its full loan remains.
Twenty-three focused tests pass from committed source. Factorer arithmetic is
unchanged by this benchmark extension; its prior full validation remains
separately identified. The full suite was not repeated for this extension.

These are conditional ECM component measurements on exposed inputs, with one
warm pool lifetime per arm/input. Two independent moduli per band cannot
establish a population p90 or general adoption confirmation. No global tuning,
new default change, merge, B12/F4 closure or E1 execution follows. Detailed
protocol, all positive/negative cells, manifests, costs and receipts remain at
`v2/audit/b12_f4/selected_defaults_master_ecm_v2/` locally.

## B12/F4 automatic default integration (5 November 2026)

Automatic cost-gated ECM is promoted on the isolated worker branch by explicit
direction after the frozen comparison below. This is an adoption decision,
not a new performance result or a claim that the balanced completion guardrail
passed. It supersedes the earlier retain-opt-in recommendation; historical
verdicts below retain their original scope. Cold 60/100 medians regress, p90 is censored and observed CPU rises;
these limitations and the pending balanced guard are explicitly retained.

The normal auto CLI and bounded API now select the measured configuration and
own its pool lifecycle. Easy exits do not load the policy or multiprocessing;
unmeasured sizes and unfunded/unsupported calls remain serial. Omitted tiers
at 60/80/100 select the measured escalating140 policy; explicit and saved plans
win. Explicit caller grants and stock QS fallback/checkpoint settings remain.
New schema-13 automatic metadata binds source implementation, policy, seed
and memory; older readers reject it before ignoring worker state.
legacy serial checkpoints do not opt into workers on resume. Integer seeds
retain the measured SEARCHv1 mapping. The packaged production policy copies all
34 measured rows unchanged, with new integration pins; historical controls
and receipts remain immutable. Later code changes invalidate the old profiles.

Focused validation includes real spawn/drain, exact manifests, unchanged serial
RNG/curve receipts, memory/startup refusal, cancellation, error cleanup and no
orphan children. New automatic schema tests reject legacy-reader/downgrade
misinterpretation. Legacy numerical/checkpoint oracles explicitly request
serial; default-path tests check schema 13 and its mathematical subschema.

Implementation `0c4cd18` with fixture corrections `094ab1d` passes committed-only
PyPy 3.11 full validation: 1,337 tests (9 optional skips), lint and 264 benchmark
imports. Final validation costs 75.474 wall / 93.613 observed aggregate CPU
seconds. All eight QA attempts cost 316.460 wall / 395.016 CPU; full conservative
booking is 890 wall / 8,900 CPU under the finite 900 / 9,000 ceiling. The initial
wrong-regime spawn fixture, seven legacy event assumptions, and later schema
fixture failures stay retained with their full loans. No arithmetic or pool
lifecycle defect was found in those failures; frozen baselines remain unchanged.
No failed measurement was rerun because no new performance campaign ran.

Detailed local decision, attempts, resource ledger and handoff are under
`v2/audit/b12_f4/default_adoption_v1/`. The computation lock is released. This
branch changes defaults through the recorded adoption exception; it remains
unmerged/unpushed, B12/F4 remain Open, and E1 stays deferred.

## B12/F4 fresh cofactor-policy comparison (5 November 2026)

The corrected selector is measured from committed-only `9b7eaeb` against actual
master `0efaba1`, with unchanged native runtime `10f64d5` / `4a2f051`. One policy
was frozen before fresh input generation: explicit escalating140 at original
60/80/100-digit inputs, cost-gated process selection at actual cofactor sizes,
and literal stock SIQS fallback/checkpoint settings. This bundles quota,
SEARCHv1 curves and worker choices; it is not a pure parallel-scaling comparison.

All 1,248 process cells validate, including 1,838 timed factorization calls and
50 fixed-work observations. The primary panel has 12 independent unconditioned
uniform integers per size, with nine fresh cold CLI samples and nine warm-call
samples per input/arm. Warm calls use three independently prepared homogeneous
contexts with at least three validated warmup seconds each. Cold includes
interpreter/import/spawn/JIT, factoring, output and shutdown. Both arms receive
the same finite grants: 8 s factoring, 16 s cold censor/20 s cell guard, or 8 s
warm-call censor/60 s pool guard. Missing, refused and incomplete calls remain
censored; repetitions do not add independent moduli.

| Digits | Cold master → policy median (s) | Warm master → policy median (s) | Completion master → policy |
| --- | ---: | ---: | ---: |
| 60 | 0.3181 → 0.4508 | 0.0662 → 0.0244 | 7/12 → 8/12 |
| 80 | ≥16 → 0.7332 | ≥8 → 0.1480 | 4/12 → 7/12 |
| 100 | 0.5138 → 0.6163 | 0.2144 → 0.0738 | 7/12 → 7/12 |

Equal-weight primary completion rises from 18/36 to 22/36, a paired difference
of +11.1 percentage points (95% modulus-bootstrap interval +2.8 to +22.2).
Finite-horizon mean latency, charging each incomplete call its full horizon,
is 8.161 → 6.524 s cold, ratio 0.799 [0.617, 0.957], and 4.046 → 3.159 s warm,
ratio 0.781 [0.591, 0.944]. These are completion-sensitive restricted means,
not raw elapsed time summed across calls or a universal speedup multiple.
Raw primary cold cell wall is 528.0 s for master versus 560.9 s for the policy.
Overall median point estimates are 0.837 → 0.621 s cold and 0.505 → 0.086 s
warm, but their bootstrap intervals are unbounded with censoring. Cold medians
regress at 60 and 100 digits. Every primary p90 remains beyond 16/8 s; N=12 per
band cannot establish a precise tail. No universal median/p90 gain is claimed.

The separate certified ECM-enriched panel completes 11/12 versus master's 4/12:
4/4 versus 0/4 at 60, 3/4 versus 2/4 at 80, and 4/4 versus 2/4 at 100 digits.
It contains direct smaller factors of 15/18 digits and recursive cofactor factors
of 17/20 digits. These conditioned inputs do not replace the uniform population.
Its finite-horizon ratios are 0.197 [0.060, 0.443] cold and 0.168 [0.028, 0.428]
warm. Small per-band cohorts and the 80-digit median loss remain visible.
Auxiliary 30/40/50/70 uniform, prime and power cases all complete; their single
cold sample per case is exploratory and shows roughly 12–50 ms extra CLI cost.

Native receipts show workers engaged on 29/36 uniform moduli in both endpoints
(261/324 timings). Most use ten workers; one smaller-cofactor case uses eight
cold or fourteen warm. Eighteen warm rows lose coarse stage attribution when
stock SIQS checkpoint construction refuses, but their native receipts explicitly
show engagement. All refusals are retained: 26 cold/27 warm stock rows and 18
warm policy rows hit `full SIQS checkpoint exceeds checkpoint_bytes`. This is
the preserved 1 MiB fallback checkpoint limit, not host/worker RAM exhaustion.
No result, arithmetic, worker, watchdog or owned-process cleanup failure occurs.

Primary cold observed aggregate CPU is 511.0 → 2,235.1 s (4.37×); conservative
CPU charges are 1,159.0 → 16,898.7 s. Warm whole-context CPU is 864.5 → 3,128.1 s,
including variable validated warmup work, so it is not an equal-work efficiency
comparison. Peak declared portfolio ownership is 2.264 GiB under an 8 GiB grant,
plus a separate 256 MiB controller. One 81-digit diagnostic observes 827.4 MiB
simultaneous known-parent/worker RSS; it is not a whole-campaign peak or an RSS
theorem. Native ten-worker cold readiness is approximately 0.18–0.20 s in the
clean fixed-work screen; component close is about 5 ms. CLI import/JIT,
metadata, cancellation and selected-wrapper close remain in whole-call costs;
complete attribution is unavailable. Last-ECM selection decision medians are
0.744 ms cold and 0.153 ms warm, excluding separately charged lazy profile setup.

A separately frozen hard pilot compares identical existing 415-curve plans under
serial equivalent-manifest M0 and forced P10. Both complete an 80-digit input
with a 22-digit smaller factor (30.207 → 6.246 s) and a 100-digit input with a
25-digit factor (12.105 → 2.397 s). Stock master exhausts without completion.
These are one timing per input/arm: useful hard-case feasibility, not confirmed
fivefold speedups or coverage of 30–35-digit factors. Balanced 60/70-digit probes
remain incomplete in both arms; the protected guardrail is pending.

Calibration fits only a new 81-digit anchor; 55/75/90-digit model holdouts
predict serial body cost within 4–17%. The selected profile
`b12_f4_targeted_selected_v1.json` is committed before the fresh corpus and keeps
all earlier 33 nominees unchanged. Committed-only `9b7eaeb` passes 1,311 PyPy
3.11 tests (9 optional skips), lint and 263 benchmark imports: 74.135 wall /
88.670 aggregate CPU seconds. Later reporting changes are documentation only.

Measurement cells consume 4,362.476 wall / 10,259.762 observed CPU seconds;
conservative measurement CPU charges are 40,536.575 s. Final booking, including
preparation, QA and all analysis attempts, is 4,862.476 wall / 41,616.575 charged
CPU seconds within the frozen 9,900 / 104,000 ceiling. A preparation path error
and a resource-only auxiliary-grouping correction stay retained, without timing
reruns or primary-inference changes. Detailed report, tables, exact receipts,
ledgers and remaining limits: `v2/audit/b12_f4/targeted_policy_v2/REPORT.md`,
`TABLES.md`, `resources-r2.json` and `HANDOFF.md`. Computation is finished and
its lock is released. Keep this scoped policy opt-in; no defaults, merge, push,
B12/F4 closure or E1 campaign follows.

## B12/F4 cofactor and serial integration (5 November 2026)

The v2 selection adapter fills bounded intermediate-size cost estimates for the
same mathematical plan and active PRAC/ladder path, without extrapolation. It
retains pool-specific crossover estimates and publishes model/anchor provenance.
Portfolio declines now use incumbent serial curves, programs and RNG consumption;
serial resume stays on that branch. A declined call does not pin the context to
one worker. Later worker refusal cannot silently switch to a serial search.

`b12_f4_selection_v2.json` is a new source-pinned training input. It preserves all
33 literal cost/settings anchors from the immutable v1 catalogue. This is an
implementation change. Committed-only `4a2f051` passes 1,308 PyPy tests (9 optional
skips), full lint, 262 benchmark imports and v2 source-pin admission. Final QA
cost is 76.231 wall / 91.446 aggregate CPU seconds. All seven QA attempts cost
86.549 wall / 103.223 CPU seconds; failed attempts and full loans remain retained.
Booking is 645 wall / 1,190 CPU seconds under a separately recorded QA-only
675/1,200 ceiling; the initial 600-wall ceiling was revised after a strict test
style failure, without changing a factoring timeout or resetting any charge.
At that QA checkpoint, new complete-cost performance and fresh confirmation
were NOT RUN; the subsequent targeted comparison above supersedes that pending
status. The protected balanced guardrail remains pending. Tables below retain
earlier source `538462a`; they are not measurements of the changed adapter. Detailed local decision/QA ledger:
`v2/audit/b12_f4/selection_integration_v2/`.

## B12/F4 actual-master complete-cost comparison (5 November 2026)

The completed comparison uses actual committed master `0efaba1` and candidate
runtime/benchmark source `538462a`; its 99 factorer runtime files remain those
of `fc6fea3`. The control resolves stock numerical ECM, SIQS fallback and 1 MiB
QS checkpoint settings. Candidate arms replace only their declared ECM options
and use the frozen source-checked selection catalogue. Both arms receive the
same explicit machine-wide time/CPU/work/memory ceilings. This supersedes the
earlier R64 result for claims about stock settings; historical receipts remain.

All 4,384 process cells and 6,568 timed calls pass reconstruction, source/request
hash and owned-process cleanup checks. There are nine fresh cold CLI timings per
ordinary modulus/arm; warm execution uses three independent homogeneous context
lifetimes, each with at least three validated warmup seconds and three timings.
Cold latency includes interpreter/import/spawn, result construction and shutdown;
warm call latency is secondary, with preparation/warmup/close separately charged.
The cold factoring deadline starts at CLI entry (8 s), with a 16 s whole-call
censor/cleanup horizon; warm calls have 8 s and a 60 s whole-pool guard. Refused,
failed and incomplete outcomes remain in denominators and censored quantiles.

The unconditioned uniform panel has 12 independent moduli at each primary size
60/80/100, and four at auxiliary sizes 30/40/50/70. Completion counts below agree
across cold/warm endpoints. Medians include noncompletion; ≥16/≥8 are unidentified
cold/warm medians within the respective horizon. Repetitions do not add inputs.

| Digits | Cold master → policy (s) | Warm master → policy (s) | Completion master → policy |
| --- | ---: | ---: | ---: |
| 30 | 0.1611 → 0.1773 | 0.0004 → 0.0004 | 4/4→4/4 |
| 40 | 0.1731 → 0.1799 | 0.0012 → 0.0010 | 4/4→4/4 |
| 50 | 0.1862 → 0.1979 | 0.0052 → 0.0015 | 4/4→4/4 |
| 60 | 0.6348 → 0.9162 | 0.1449 → 0.2339 | 6/12→6/12 |
| 70 | ≥16 → ≥16 | ≥8 → ≥8 | 0/4→1/4 |
| 80 | 0.5589 → ≥16 | 0.1832 → ≥8 | 7/12→5/12 |
| 100 | 0.6112 → 0.4754 | 0.1193 → 0.1049 | 6/12→7/12 |

On the equal-weight primary population, the default-schedule policy completes
18/36 versus master's 19/36. Cold median is 0.916 versus 0.586 s; warm median is
0.256 versus 0.183 s. Its restricted-mean latency ratios are 1.059 cold
(95% modulus-cluster interval [0.889, 1.296]) and 1.058 warm [0.884, 1.301].
There is no broad improvement, and the 80-digit completion/latency loss matters.

A separately frozen supplementary panel applies the previously selected explicit
140-curve escalating schedule at all 36 primary uniform inputs. It was added
after some core outcomes were visible, without retuning its settings; it is a
supplementary comparison, not a pristine independent confirmation.

| Digits | Cold master → policy (s) | Warm master → policy (s) | Completion master → policy |
| --- | ---: | ---: | ---: |
| 60 | 0.6348 → 0.3555 | 0.1449 → 0.0463 | 6/12→9/12 |
| 80 | 0.5589 → 0.8683 | 0.1832 → 0.4477 | 7/12→8/12 |
| 100 | 0.6112 → 0.4344 | 0.1193 → 0.1341 | 6/12→7/12 |

Together, caller140 completes 24/36 versus 19/36. Cold median is 0.497 versus
0.586 s; warm median is 0.171 versus 0.183 s. Those point improvements are
uncertain: median-ratio intervals are [0, 3.315] cold and [0, 7.693] warm.
Finite-horizon mean latency is lower: 5.872 versus 7.720 s cold, ratio 0.761
[0.554, 0.946]; 2.996 versus 3.807 s warm, ratio 0.787 [0.586, 0.976]. The paired
completion difference is +5/36, interval [1/36, 9/36]. These are scoped,
censoring-aware comparisons, not an exact median speedup or a p90 claim.
Primary p90 remains beyond the 16/8 s horizons in both arms. The 80-digit
caller140 median regresses despite one additional completion.

The separate certified ECM-enriched panel (four independent inputs per size,
15/17/18/20-digit smaller factors) gives caller140 completion 4/4 versus 1/4 at
60, 4/4 versus 0/4 at 80, and 3/4 versus 0/4 at 100. Candidate cold medians are
1.004/0.808/1.321 s and warm medians 0.220/0.083/0.715 s. Control medians are
unidentified within these grants, so no exact speedup multiple is established.
All 33 recommended regimes were also compared on four other certified enriched
inputs per regime, including every loss. Larger schedules improve several
conditional 70/80/100-digit cells; they do not supply universal winners.

An integration limitation remains: exact-digit/options catalogue lookup often
misses the cofactor produced by preprocessing/recursion and falls back to SEARCHv1
serial. Caller140 engages on only 2/36 cold uniform moduli (18/324 timings),
including none at 100 digits; its warm panel has 9 engaged, 288 serial and 27
refusal rows without a decision receipt. The default panel has zero engaged
cold timings and nine warm. Thus quota/stream changes, serial choices and
parallel execution are bundled; these completion gains cannot be attributed to
workers alone. The opt-in serial SEARCHv1 stream differs from incumbent S0.

Matched primary cold observed aggregate CPU is 578.4 s for S0, 591.2 s for the
default-schedule policy and 912.5 s for caller140 (324 calls each); conservative
charged CPU is 1,226.4/1,239.2/2,532.1 s. Warm whole-context observed CPU is
941.2/1,000.5/1,272.4 s for 108 contexts per arm, including variable validated
warmup work, so it is not equal-work CPU efficiency. All measurement cells total
15,765.3 wall, 33,369.9 observed CPU and 112,869.3 conservatively charged CPU
seconds. Peak declared owned memory is 2.819 GiB. Parent/one-child high-water
marks (296/89.8 MiB) are not simultaneous aggregate RSS. Earlier separate
14–16-worker probes observed about 1.19–1.50 GiB parent-plus-worker RSS. Full
startup/shutdown is in endpoint/cell wall time; selected-wrapper close attribution
is missing, and native counters alone cannot establish total cancellation cost.

One cumulative supplementary reservation refused before dispatch; all 418
acknowledged cells stayed charged. A prospectively committed resource-only
revision raised that lane from 2,000 to 3,000 cumulative wall seconds within the
unchanged 18,000-wall/480,000-CPU campaign ceiling. It completed the remaining
14 cells without reruns or treatment changes. Stock checkpoint refusals remain
censored outcomes; no R64 substitution or successful-only analysis was used.
The 60/70 balanced probes remain incomplete: the protected guardrail is pending.
No default change, global optimum, B12/F4 closure or E1 campaign follows.

Final committed-only report/input checkout `20e6681` passes 1,287 PyPy tests
(9 optional skips), full lint and 262 benchmark imports: 75.649 wall / 90.367 CPU
seconds, with the full 90/220-second reservation retained. Final campaign booking
is 16,345.315 wall / 114,769.283 charged CPU seconds within 18,000 / 480,000. The
shared computation lock is released. Detailed all-regime tables, corpus proofs, receipts, uncertainty and
resource ledger are local under `v2/audit/b12_f4/master_comparison_v1/REPORT.md`,
`ALL_REGIMES.md`, `summary.json`, `cost-summary.json` and `HANDOFF.md`.

## B12/F4 configuration selection and held-out completion (4 November 2026)

Runtime `fc6fea3` includes committed production fallback `01d00f2`. It passes
1,270 PyPy tests (9 optional skips), lint and 261 benchmark imports from a
committed-only checkout. Subsequent benchmark analysis/input changes pass their
focused checks; final full committed-checkout QA is queued. Exact arithmetic,
SEARCHv1 identities, parent leases and serial defaults remain unchanged.

A finite campaign screens 0/1/2/4/6/8/10/12/14/16 workers and bounded reuse,
polling, packet, stage and budget options across 33 explicit size/schedule
regimes at 30/40/50/60/70/80/100 digits. Repeated finalists comprise 1,188
invocations; a targeted decoded-cache ablation adds 108, and missing-state cost
calibration adds 95. Cold timings use fresh interpreters/pools; warm timings use
three independent pool lifetimes with validated warmup and three samples each.
Timing-block labels do not establish temporal pairs: uncertainty resamples
independent processes/pool clusters. These are conditional work measurements,
with near ties, rather than a proof of a global optimum.

Cold cheap calls stay serial under the selected policy. Heavier cold nominees
usually use 4–12 workers; many retained-pool nominees use 12–16. Exactly 140
escalating curves at 60/80/100 select ten workers in both intents. For example,
80/100-digit caller140 fixed-work cold medians are 1.790/1.961 seconds versus
5.985/7.083 serial; warm medians are 1.120/1.382 versus 5.501/6.538. Mathematical
options and curves match within these work comparisons. Counts above ten
oversubscribe the M4; worker IDs imply no physical core class. Dropping the
near-tied 32 MiB stage-two index cache for the 80-digit deep warm nominee lowers
retained state; it does not establish that indexing is generally harmful.

Before fresh generation, the source-pinned catalogue, engagement rule, seeds,
exclusions and 12-modulus scope were frozen. Expected remaining serial ECM work
must exceed fitted cold/started/warm thresholds and a dispatch-wave floor;
cold intent has a one-second minimum. Retained intent assumes 64 compatible
ECM-reaching calls while paying actual first-call startup. Unknown regimes stay
SEARCHv1 serial, distinct from incumbent S0. Full configuration snapshots and the
frozen catalogue are versioned loader inputs.

All 448 fresh calls return checked results with source/seed/request hashes and
owned-group exit proofs. The primary selected scope has 12 independent moduli,
60/80/100 digits with 15/17/18/20-digit smaller factors, and nine fresh cold CLI
timings per modulus/arm. It uses explicit caller140 and matched **64 MiB QS
checkpoint funding**, labeled S0[R64]/POLICY[R64]. A prior unchanged stock
1 MiB QS checkpoint serialization failure remains charged and visible. The
larger grant changes neither numerical fallback settings nor the codec, and
these results do not describe stock resource defaults.

| Cold CLI endpoint | S0[R64] | POLICY[R64] |
| --- | ---: | ---: |
| Complete independent moduli | 10/12 | 11/12 |
| Median latency | 2.920 s | 0.937 s |
| Mean latency restricted at 20 s, preserving noncompletion | 5.343 s | 2.496 s |

The median ratio is 0.321 with modulus-cluster 95% interval **[0.179, 1.110]**:
the 3.12× point improvement does not establish a general median gain. The
prespecified finite-horizon mean ratio is 0.467 [0.140, 0.867]; the completion
difference is +1/12 [0, 3/12]. Four moduli per band give candidate/control
medians 0.688/0.848 seconds at 60, **1.272/1.048 at 80**, and 0.672/5.161 at 100.
The 80-digit regression remains part of the result. The 100-digit subgroup
cannot be selected retrospectively as a newly confirmed scope. Twelve moduli
support only descriptive, uncertain p90 reporting; repeated timings are not
independent inputs. The extra algorithm seed is exploratory and inconclusive.

The same-manifest M0 comparison uses only three timings per modulus and remains
diagnostic. A 52-input unconditioned uniform panel and 28 other enriched inputs
have one timing each: exploratory completion is 35/52 versus 31/52 and 15/28
versus 12/28 respectively. Uniform bands remain separate, with equal-weight
60/80/100 summaries; the pooled auxiliary-size panel is not portfolio acceptance.

Across the same 200 primary/exploratory calls per arm, observed aggregate CPU
is 896.52 seconds for POLICY versus 465.57 for S0. Conservative charges are
7,957.39/865.57 seconds, including retained startup/helper/shutdown allowances;
charges are not measured CPU consumption. Four separately priced RSS probes on
14–16-worker warm configurations observe simultaneous parent-plus-worker peaks
of 1.19–1.50 GiB, versus 2.24–2.83 GiB owned reservations. These omit separately
reserved helpers/controller; sampled RSS and payload bounds are not an exact
whole-process RSS theorem. Instrumented timings do not rank candidates.

Both balanced 60/70-digit pilots remain unresolved: 60 reaches the common
12-second body limit in QS, while 70 exhausts its numerical schedule without
a QS fallback. Both policy calls stay serial. This cannot pass the balanced
completion guardrail. No default promotion, B12/F4 closure, broad portfolio
claim, merge or E1 campaign follows. The current campaign books 14,499.59 wall
and 165,825.32 aggregate CPU seconds under finite 18,000/216,000 ceilings;
bookings retain failures and unused reservations rather than claiming actual
runtime. Earlier phases have separate preserved ledgers. The computation lock
is released; final tooling QA is queued for a later window. Detailed research,
all 33 regime choices/near ties, receipts, budgets and handoff remain local at
`v2/audit/b12_f4/regime_selection_v1/`.

## B12/F4 512-curve cache amortization (4 November 2026)

A resource-valid follow-up uses the unchanged committed runtime `b62ac7a` and
prior optimized F4 control `344f8f1`, with 5 MiB checkpoint ownership in every
arm. Pure checks in both sources confirm spawn and the 4,423,680-byte checkpoint
bound, including simultaneous buffers. Full committed QA remains 1,046 tests
(9 optional skips) and lint; no arithmetic/executor/default changes accompany
this measurement phase.

One exposed 80-digit modulus with a 20-digit smaller prime uses 512 identical
curves at `(2000,200000)`, seed 7919 and eight spawn workers. N0 is the prior F4
component control; K0 is current options-off; T adds tail-aware packet2; MT adds
2 MiB index preparation to T. N0 is not production S0. Every call checks all
512 curve identities/results. Two declared curves return proper divisors; this
fixed-work search continues after factor discovery.

Five retained-pool blocks include a repeated N0 drift anchor. Each validates
at least three warmup seconds before nine samples. The primary component CLI
comparison has nine interleaved fresh-process samples per arm, including all
startup/JIT/preparation, validation and shutdown costs. A blocking wait replaces
the earlier 50 ms controller polling grid.

| Arm | Warm median | Fresh component CLI median | Mean observed cold CPU |
| --- | ---: | ---: | ---: |
| N0 | 1.48508 s | 2.40574 s | 16.706 s |
| K0 | 1.55182 s | 2.48399 s | 16.945 s |
| T | 1.52758 s | 2.33625 s | 16.761 s |
| MT | 1.43326 s | 2.35565 s | 16.628 s |

Index preparation lowers warm time 6.2% versus T and 7.6% versus K0, with zero
measured index rebuilding after warmup. The warm gain versus N0 is only 3.5%,
while the end anchor drifts +3.3%; source/JIT/block effects remain visible.

Cold MT/K0 has matched timing-block ratio **0.96070**, conditional 95% interval
**[0.93870, 0.98646]**. MT/T is 0.98691 [0.95413, 1.01096]; MT/N0 is 0.98582
[0.94927, 1.01978]. The latter two include parity. Thus the new options-off
comparison supports a scoped cold improvement, but improvement over the actual
prior F4 incumbent remains inconclusive. K0's slower control must not manufacture
an incumbent gain. One exposed modulus supplies no independent-input, population
p90, first-factor or complete-factor claim, and the count is not a universal optimum.

Owned reservations are 1.799–1.818 GB within the 4 GiB executor ceiling, plus
128 MiB controller. Separate cold parent/worker high-water observations are
about 84–86 MB and 638–649 MB summed worker HWMs, not simultaneous aggregate RSS.
Worker readiness is about 0.17 s; whole interpreter/JIT costs remain in the CLI
endpoint. Normal close is about 4 ms; early-factor cancellation is unmeasured.

All 45 invocations return matching arithmetic results. One feasibility recorder
breaks on a transient group-proof permission error; an own-PGID-only late check
proves exit. Its full 20-wall/160-CPU grant remains charged and its failed timing
is excluded. The revised bounded proof passes normal, deadline and deterministic
transient-error checks; only unissued feasibility cells continue. All groups
exit and the lock is released. All earlier refused/failure receipts stay intact.

New-phase booking is 269.18 wall / 3,425.99 aggregate CPU seconds within 600 / 4,800;
observed runtime is 181.24 wall / 1,327.87 CPU seconds, with helper CPU partly
unobserved and conservatively reserved. Prior ledgers remain unchanged. Evidence
is about 3.4 MB under 128 MiB. No default promotion, merge, B12/F4 closure or E1.
The next gate is a small root-lineage-aware first/complete-CLI feasibility and
observable serial/pool engagement calibration. Existing component factor hits
do not guarantee identical portfolio search or completion. Detailed current
report, freezes and handoff: `v2/audit/b12_f4/long_quota_v2/`.

## B12/F4 checked fast paths and complete-cost screen (4 November 2026)

Committed executor source `b62ac7a` passes 1,046 PyPy tests (9 optional skips)
and full lint in a committed-only checkout. Opt-in resource/checkpoint v6 adds
prospective checked work accounting, immutable stage-two index/group preparation,
context binding and tail-aware two-curve packets. Cached indices are read without
per-group tuple copies. Exact arithmetic, mathematical v1 streams, accepted PRAC
eligibility and global serial defaults remain unchanged. All new flags default off.

N0 is prior F4 component source `344f8f1` with previously selected X8 settings;
K0 is the new source with the new flags off. Neither is production S0. T uses
packet2 with singleton tails; M prepares bounded stage-two indices; B uses the
prospective budget. All comparisons here use eight spawn processes. The earlier
4–16-worker results remain separate; this phase does not establish a new count
optimum.

Nine warm calls per cell follow at least three validated warmup seconds. L uses
128 curves at `(2000,200000)`; H uses 32 at `(11000,1100000)`. L/H share one exposed
modulus per size, each with a 20-digit smaller prime. Only L60/L80 use accepted
PRAC; L100 and H retain the ladder.

| Case | N0 warm median | K0 warm median | T warm median |
| --- | ---: | ---: | ---: |
| L60 | 0.33935 s | 0.33745 s | 0.30594 s |
| H60 | 0.43721 s | 0.43203 s | 0.40565 s |
| L80 | 0.41476 s | 0.41555 s | 0.40790 s |
| H80 | 0.51893 s | 0.52522 s | 0.49291 s |
| L100 | 0.52538 s | 0.47147 s | 0.47406 s |
| H100 | 0.61679 s | 0.61982 s | 0.60268 s |

T's equal-weight geometric ratio is 0.93979 versus N0 and 0.95669 versus K0:
about 6.0% and 4.3% lower warm fixed-work time. B at L60 and M at L80 lower
new-source medians about 11.8% and 15.3%, respectively. Combining all options
adds little: 2.2% below N0 and 0.45% below K0; binding alone loses 2.1% versus K0.
L100's source/control difference must not be credited to the packet policy.
Later metadata/packet combination blocks are exploratory. Nine timing repeats
on three moduli do not establish independent-input improvement or population p90.

The primary endpoint uses 81 fresh CLI invocations, nine samples per arm/case,
including interpreter/imports, spawn/JIT, preparation, validation and shutdown:

| Case | N0 cold median | T cold median | Additional challenger |
| --- | ---: | ---: | ---: |
| L60 | 0.92732 s | 0.92407 s | M+T: 0.95633 s |
| L80 | 1.03327 s | 1.03378 s | M: 1.04356 s |
| H100 | 1.20523 s | 1.20822 s | M+T: 1.22050 s |

All T conditional paired timing-block intervals include parity. The controller's
50 ms exit polling also limits small differences. **No convincing new cold-CLI
gain is established.** Extra preparation loses on these short calls. Observed
cold CPU is roughly 5.6/6.2/7.6 seconds; conservative charges retain unpublished
startup/helper/shutdown CPU. Separate parent/worker high-water observations are
about 80 MB and 613–623 MB summed worker HWMs; these are not simultaneous aggregate
RSS. Owned reservations span 1.80–1.93 GB within the 4 GiB executor ceiling, plus
128 MiB controller ownership. Warm reuse eliminates measured index rebuilding.

An eight-call ECM-only first/complete pilot validates certified outputs and
unresolved reconstruction. Both arms complete the 10-/70-digit 80-digit fixture
in about 0.49 s. Both exhaust the exposed 18-/42-digit 60-digit search without a
split. These single observations establish no complete-factor gain; SIQS/SSS
are explicitly disabled, so this is not unchanged-production-fallback evidence.

All 161 runtime cells validate and drain. A separately frozen 512-curve probe is
refused before worker launch: its checkpoint needs 4,423,680 bytes against 4 MiB
funded. This preparation error retains its full 35-wall/240-CPU grant; five later
cells remain NOT RUN. A future resource-valid freeze needs at least 5 MiB and
must reprice simultaneous checkpoint ownership before a new finite allocation.
The earlier failed focused schema check also stays retained/full charged.

New-phase maximum booking is 1,075.20 wall / 11,635.99 CPU seconds within the
1,800 / 12,000 envelope, separately from all previous phases. Observed QA/runtime
CPU lower bounds are about 166/3,861 seconds; helper publication remains incomplete.
No calibrated engagement rule, default adoption, merge, B12/F4 closure or E1.
Larger-quota amortization, matched production complete CLI, balanced 60/70 controls,
uniform random panels and fresh independent-input confirmation remain pending.
Detailed freezes, receipts, report and handoff: `v2/audit/b12_f4/headroom_impl_v2/`.

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
[pinned history](https://github.com/nishanth17/factor/blob/23458a52e409add326c9fd63cf1f322c6e605569/v2/history/README.md) for frozen source/accounting detail.
Numbers from separate comparisons are not additive portfolio gains.

## Collection and larger-band decisions

### C9a preparation

The bounded 4 November continuation has three separate matched observations:

| Comparison | Previous/control seconds | New/C9a seconds | Timing screens |
| --- | ---: | ---: | --- |
| Same-C9a full60, root transfer (`891eda5` / `ebb35fa`) | 168.974 | 157.717 | false / true |
| Same-C9a full60, hit words (`ebb35fa` / `c8168dc`) | 159.519 | 132.247 | false / false |
| Full70 policies, shared optimized `891eda5` | 2742.095 | 1842.083 | true / false |

Both same-C9a full60 pairs retain 764,077,765 positions, 4,977 rows, 5,830
polynomials and exact factors; additional helper construction/read work stays
charged. The final hit pair has 17.1% lower raw time and 1.5% higher peak RSS
(347,127,808/352,403,456 bytes). Nine warmed 2,048-block collector samples have
12.8%/17.4% lower median cost at 60/70 digits, favorable 9/9 at both sizes,
with identical stores/progress and explained extra work. This is component
evidence; no full70 word-helper result is available.

The separate extended full70 freeze grants 9,000 production seconds per arm,
4 GiB owned storage and 5 GiB process RSS. Both complete with proven factors.
C9a uses 34.3% fewer positions and 10.1% lower peak RSS, while its reserved
workspace is higher. The policy comparison uses the optimized shared runtime,
not a fresh master call. Earlier 900s outcomes and all charges remain preserved.
Each size has one exposed modulus/seed; repeated component samples are not
independent inputs. Historical timing flags are preserved, but the uncalibrated
15% heuristic is no longer a veto in new analyses: warnings accompany retained
validated observations. Independent-input variability and fresh confirmation
remain unmeasured. Paired gains are not multiplied into a population speedup.
Final committed-only
checks pass 1,113 tests (9 skipped), lint and 243 benchmark imports.

One calibrated whole-factor sample on `ebb35fa` attributes substantial cost
to marking, queue construction and GC (37.102s, 23.4% of traced wall). Its
complete mature heap census describes shallow categories, not retained
dominators or allocation rates. Diagnostic timing is separate from clean
evidence. Local detail: `v2/audit/c9a/optimization-83-89.md` and
`verification-89.json`; generated data: `results/c9a/hit-words-89/`, with
earlier worktree locations recorded in the local report. The bounded phase
stops here. The delivered C9a stream is integrated and closed by explicit
direction on 4 November; failed timing/fresh-confirmation gates remain unmet.
C9b owns subsequent validation and headroom analysis. Existing C9a policy
classes remain opt-in.

The final shared-source flat-queue pair on `180576c` factors the same exposed
60-digit input in 268.8/230.4s for incumbent/C9a: 14.3% lower raw candidate time.
Collection is 215.3/202.8s, filtering 9.07/0.90s and preparation 16.49/0.67s.
Peak process RSS is 289.0/331.1MB; accounted workspace is 615.0/567.7MB. Both use
the same 900s production allowance, work grant, input/seed and observer.
Exact reconstruction and independent proven-corpus agreement pass; terminal
work, positions, rows and semantic counters match each arm's earlier packed
queue run. Both revised timing screens fail: control samples are unstable;
C9a's stable post-run ratio is 0.775. Preserve these failures: the raw difference
does not establish a clean population speedup or adoption. Committed-only
1,084 tests (9 skipped), lint and 242 imports pass. Packed-field replay is
replaced by faster flat words; fused scanning remains optional after mixed
integration evidence. Local verification: `v2/audit/c9a/verification-77.json`;
full scope/retained and rejected changes: `v2/audit/c9a/optimization-65.md`.
No further run or quota/default change follows this finite phase.

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
