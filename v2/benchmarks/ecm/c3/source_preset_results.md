# Upstream preset transfer screen: initial results and stop decision

**Decision, 11 October 2026: defer promotion and retain numerical defaults.**
All 54 matched groups and 324 calls complete and validate, with a proper factor
and exact reconstruction in every call. Yamaquasi's escalating prefix and
Alpertron's first two ladder tiers are promising on these six revealed inputs.
The frozen assessment requests an 18-sample extension because uncertainty is
wide; it selects no policy. The user said the work was taking too long. This pass stops after the initial
assessment: no 18/27-sample batch and
no fresh confirmation are run. Preserve the recommendation and raw evidence;
do not reinterpret this as passing the stability gate.

## Complete initial assignment

Six revealed round-two subjects cover balanced/uneven_12/uneven_14 at 30 digits
and balanced/uneven_12/uneven_16 at 40 digits. Each has nine fixed seeds across
six arms, in rotated/reversed order. Every acquired arm/band has at least three
validated CPU and wall seconds of warmup, with forced SIQS warmed separately.
Cold startup and instrumented calibration are excluded from these call costs.
All arms receive the identical 10^13-work, 5/30-second CPU/wall, 288 MiB envelope,
frozen C1 relation configuration and positive fallback floors. The complete
report checks full SIQS/allocation identity within every matched group.

| Arm | Mean CPU seconds | Mean wall seconds | SIQS calls / 54 | Mean actual curves |
| --- | ---: | ---: | ---: | ---: |
| optimized fixed32 | 2.1081 | 2.1334 | 15 | 15.33 |
| fitted compact64/wide64 | 1.8103 | 1.8382 | 10 | 23.65 |
| Yamaquasi automatic eight | 3.9507 | 3.9805 | 51 | 7.85 |
| Yamaquasi escalating 140 | 1.1713 | 1.1797 | 1 | 35.85 |
| Alpertron 115 | 1.1091 | 1.1257 | 1 | 18.09 |
| GMP-ECM 74 | 1.8228 | 1.8316 | 0 | 10.46 |

Tier counts are ceilings, not curves forced after a factor is found. GMP-ECM's
row avoids SIQS entirely here but costs more than the two escalating schedules;
proper-factor yield alone does not establish better complete-call economics.

The declared equal-band/class/input weighting gives the following CPU cost
reductions. Negative percentages mean a regression. Intervals are the frozen
10,000 paired bootstrap draws; they are conditional on this small revealed
population and are not corrected for selection among the four source arms.

| Source arm | Reduction vs fixed32, 95% interval | Reduction vs fitted64, 95% interval |
| --- | --- | --- |
| Yamaquasi escalating 140 | 44.44% [14.37%, 66.47%] | 35.30% [1.52%, 59.49%] |
| Alpertron 115 | 47.39% [13.23%, 69.66%] | 38.73% [-2.37%, 64.47%] |
| GMP-ECM 74 | 13.54% [-16.82%, 37.92%] | -0.69% [-33.82%, 26.07%] |
| Yamaquasi automatic eight | -87.40% [-136.52%, -51.99%] | -118.23% [-194.16%, -69.27%] |

Both promising source arms pass the initial wall and completion gates, but
neither meets the prespecified <=20-percentage-point interval-width condition.
Alpertron's comparison against fitted64 also crosses zero. These observations
are hypotheses for confirmation, not accepted speedups or an optimum.

## Scope, resources and limitations

Fitted compact64 remains fastest at 30 digits: mean CPU 0.1367 seconds versus
0.1678 for Yamaquasi escalation and 0.1964 for Alpertron. At 40 digits the
corresponding means are 3.4840, 2.1747 and 2.0218 seconds. A size-dependent
combination is a plausible next hypothesis; it was not selected, constructed or
fresh-confirmed in this study. Expected factors never configure dispatch.

Exact SIQS parameters, cumulative reservations and the outer memory cap remain
unchanged. Curve workspace reservations range from 12,209,152 to 18,819,072
bytes across these arms. Maximum reported fallback-owned workspace is
123,567,712 bytes. No relation recovery is reported in these calls. The process
RSS high-water mark is 227,229,696 bytes in every arm's captures; it is cumulative
within the shared process and cannot attribute memory independently to an arm.
Existing recursive/resume/cancellation tests cover the allocation contract;
this six-subject screen is semiprime timing evidence, not a fresh resume study.

The six subjects are previously revealed nonuniform Pocklington-generated
training inputs, one per stratum. Repeating seeds does not create independent
inputs. Easy eight/ten-digit factors, broader structures and other input sizes
are absent. Native source curve/stage-two assumptions differ from PyPy. The
larger service grant is matched across arms and does not establish useful SIQS
reachability under the library's default 2M-work grant. Numerical defaults,
C3's calibration gate and broader G1/E1 remain unchanged/open.

## Reproduction and verification

Source and controls are committed before measurement at `e73dbac`; the screen
manifest is `inputs/controls/c3_source_presets_v1.json`, SHA-256
`9c4174e116d3750c36fd3b8623a5d7a371375739e236a06f9407fca7139e78c0`.
The full assignment and numerical tables are in the
[protocol](source_preset_protocol.md) and [source table](source_presets.md).
The superseded 513-group run remains separately
[inconclusive](round2_comparison_scope.md), with no selection from partial data.

    pypy3 -B -m v2.benchmarks.ecm.c3.source_preset_study measure CAPTURES --lease-seconds 1200
    pypy3 -B -m v2.benchmarks.ecm.c3.source_preset_study report CAPTURES

Use the frozen source commit for reproduction, not an automatic continuation of
the stopped local study. Raw group manifests, rows, assessment and stop receipt
remain under ignored `results/c3/source-presets-v1/`. Measurement including
warmup/capture takes 760.002 active seconds and assessment 12.931, for a
772.933-second cumulative ledger. All C3 heavy work is stopped and the machine
is explicitly released to C9.

Working checks and the final committed-files-only archive pass 653 PyPy 3.11
tests and full lint. The archive imports 154 benchmark modules and verifies
35 certified corpora, 57 training subjects, the 86-source training chain,
assessment/comparison controls and all 12 source-screen pins. No collector,
arithmetic kernel, production numerical default or v1 source is changed by the
preset addition.
