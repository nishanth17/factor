# C3 allocation and handoff acceptance

Decision date: 10 October 2026, America/Los_Angeles.

Adopt the opt-in allocation/resume contract; retain numerical defaults in both
bands. Fresh confirmation rejects promotion of the training-selected
quick8 bundle. C3's broad calibration gate remains open, with G1/E1 explicit.

## Implemented contract and research contribution

The opt-in `ECMAllocation` interface distinguishes cumulative automatic
pretesting from explicit finite campaigns. Optional work cannot consume the
selected relation engine's work/wall/CPU admission floor. Mandatory exact
classification/power checks, all engines, recursive children, failed work and
resume share one ledger. Existing simultaneous owned-memory validation also
covers the relation engine. A one-way transition records any abandoned curve's
seed/bounds/frontier and preserves its charges. Insufficient admission,
pretest exhaustion, schedule exhaustion and real resource/cancellation stops
remain distinct.

Schema 12 binds the complete configuration and policy identity. Python implicit
restore reconstructs that configuration; the CLI restores allocation, tiers,
fallback and memory while retaining its documented requirements for other
nondefault flags and total grants. Snapshot verification bytes and context
rebuilding are charged; serialization/verification elapsed time carries
forward. Live SIQS/SSS admission is retained across resume. Legacy schemas
2–11 retain their prior executors. `v1/`, numerical defaults and relation
collectors/splitters/graphs are unchanged.

[Research](research.md) traces finite effort and fallback paths through pinned
YAFU, Yamaquasi, FLINT, PARI and GMP-ECM sources, primary ECM probability and
stage-two papers, official parameter notes, an author erratum, technical blogs
and implementation discussions. Its practical contributions are separate
pretest/campaign modes, credit for completed work, residual-aware economics,
protected handoff and explicit setup amortization. Native curve tables,
input-digit ratios and timing claims remain hypotheses. The first three v1
bound pairs are tested in finite bundles; its much larger curve cap and
floating/Python-2 setup are not transferred.

The reachable-control diagnosis agrees with the earlier evidence: 32 unresolved
native curves charge 5,571,485 ECM units in the feasibility pilot, already above
the 2M library total. C1's two completing fresh 40-digit DLP inputs have median
work 11,460,521,604 and 11,877,411,389. Those units are conservative accounting,
not comparable CPU operations across engines. B3's final native portfolio
interval crossed zero; C1's balanced benefit does not establish broad dispatch.

## Frozen study and selection

Final implementation source is `a809391`, frozen at `abe740c`; control is
integrated `b3b3cfbea6105f08db0a8484ec09c8260d310e28`, loaded from committed
immutable production bytes. [Protocol](protocol.md),
[revision history](controls.md) and
[selection](../../inputs/controls/c3_selected.json) preserve the freeze chain.
The original pilot is cold feasibility evidence only. Revision-4 training was
stopped for the SSS nested-config restore repair and charged 300 seconds;
its partial rows do not enter selection. Revision-5 training was interrupted
by the user after 337 persisted calls. A separately committed recovery freeze
pins those receipts, charges 600 seconds against the remaining 2100-second
training grant, rewarms the unfinished input and executes only 23 missing
assignments. Its aggregate charge is 698.602 wall/697.617 CPU seconds, including
the conservative prior charge. Combined training stays inside the original
2400-second cap; no confirmation inputs existed during either repair/recovery.

Ten historical inputs (balanced, uneven, smooth p−1, powers and close factors
at each size), seeds 7/29 and three training samples per arm give 360 calls.
Every training call completed and independently validated. Training ranks by
completion, summed capped cost, then fixed arm order. The following sums are
selection evidence, not accepted performance estimates; each cell has 30 calls.

| Arm | ECM tiers | 30-band capped seconds | 40-band capped seconds |
| --- | --- | ---: | ---: |
| control | 32 × 2000/147396 | 4.916 | 43.743 |
| no_ecm | no curves | 3.812 | 103.100 |
| quick8 | 8 × 200/7700 | **3.497** | 46.815 |
| short4 | 4 × 2000/147396 | 4.057 | 47.060 |
| tiered | 8 × 2000/147396 + 2 × 11000/1873422 | 4.108 | 49.798 |
| wide1 | 1 × 50000/12746592 | 8.821 | 54.855 |

Selection commit `faabdea` chooses quick8 for the <=35-digit bundle and retains
control above it. Fresh generation follows in `5feb17a`, with 15 independently
certified inputs; no alternative winner may be selected from confirmation.
The generator's large-q Pocklington construction biases p−1. Actual smaller
factors in the two nominal 12-digit cases have 13/12 digits, and the nominal
40-digit 16-factor case has 17 digits. Fresh labels/certificates enter only
validation/reporting, never dispatch. There are nine <=35-digit subjects
(including a nine-digit close control) and six larger subjects, including a
40-digit prime. Repeated timings do not create more independent inputs.

All service arms receive 10^13 work, 5 seconds wall/CPU in the smaller bundle
or 30 seconds in the larger, 288 MiB owned memory, native integer arithmetic
and one process/core. They consume frozen C1 calibrated SLP/DLP bundles via
the API. The CLI tests validate the allocation/resume contract; the short CLI
example is not a reproduction of every C1 collector setting. Explicit policies
protect 500M work/1 second at the smaller size or 15B/10 seconds at the larger.
Pretests cap cumulative optional work at 2M and wall/CPU at 0.5 seconds.
No production allowance is raised.

Host: Apple M4, 10 logical CPUs, 24 GiB physical memory, macOS 26.6.2 (25G83).
Runtime: PyPy 7.3.23 implementing Python 3.11.15. Accepted calls are serialized
under `/private/tmp/factor-performance.lock`, separate from C9 and heavy QA.

## Fresh confirmation

Confirmation completed in **528.586 wall / 506.249 CPU seconds**, inside the
7200-second envelope. All **972/972 calls completed**, reconstructed and passed
factor/exponent/certainty validation. All 15 independent inputs received the
planned seeds and matched samples. Three inputs (nominal uneven-12 at 30
digits, the square and close control) extended both arms from 9 to 18 to 27
samples; the other twelve stopped at 9. Seven of 72 input/arm/seed cells still
exceeded 15% relative IQR at the maximum. They are reported unstable, with no
further measurements or reselection.

For the nine smaller-band subjects, quick8's capped complete-call cost is
**51.56% higher** than control. The conditional input-subject 95% interval for
cost change is **−5.36% to +202.05%** (equivalently saving −202.05% to +5.36%).
The completion-regression gate passes, but positive-saving and stability gates
fail. **Retain the control; do not promote quick8.** The two balanced subjects
alone show 14.78% higher cost, with a saving interval −58.59% to +16.51%.
The larger-band selection already retained control, so its fresh measurements
are validated outcomes without a challenger-speed comparison.

The table gives equally weighted means of per-input/seed medians. All calls
completed, so capped cost equals complete-call wall time. Raw call counts
include the finite stability extensions and are not independent subjects.

| Bundle | Complete / calls | Proper-factor yield / calls | Wall seconds | CPU seconds | Work units | Fallback executed / calls |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Smaller control | 405/405 | 405/405 | 0.054543 | 0.054304 | 2,550,170 | 9/405 |
| Smaller quick8 | 405/405 | 405/405 | 0.082667 | 0.082465 | 24,442,273 | 126/405 |
| Larger control | 162/162 | 135/162 | 1.790687 | 1.789924 | 3,423,102,797 | 45/162 |

The larger prime correctly contributes no proper factor despite completing.
For smaller subjects, the mean work count is 9.58 times control; that is an
accounting observation, not a cross-engine CPU conversion. Fewer cheap ECM
attempts do not imply an economical complete factorization: control records
153 successful ECM events versus quick8's 36 over their matched 405 calls.
Quick8 executes fallback fourteen times as often. This explains the direction
of the regression without asserting a fitted per-curve success probability.

| Smaller-band class | Independent inputs | Control milliseconds | Quick8 milliseconds | Cost increase |
| --- | ---: | ---: | ---: | ---: |
| Balanced | 2 | 187.323 | 215.000 | 14.78% |
| Uneven 8 | 1 | 4.465 | 5.158 | 15.53% |
| Nominal uneven 12 (actual 13) | 1 | 41.035 | 168.966 | 311.76% |
| Uneven 14 | 1 | 59.469 | 126.788 | 113.20% |
| Recursive | 1 | 8.347 | 9.874 | 18.30% |
| Square | 1 | 0.577 | 0.680 | 17.81% |
| Close | 1 | 0.356 | 0.422 | 18.58% |
| Smooth p−1 | 1 | 1.997 | 2.115 | 5.90% |

Every listed class retains 100% completion and proper-factor yield. Single-input
rows are regression controls; their degenerate subject-bootstrap intervals
would not quantify population uncertainty. The larger control's class costs
are 5.276850 seconds balanced, 0.006314 uneven-8, 0.118349 uneven-12,
0.063870 nominal uneven-16 (actual 17), and 0.001889 prime. Completion is 100%
in each class; proper-factor yield is 100% for composites and 0% for the prime.

Actual finished curves range 0–32 at 2000/147396 for control and 0–8 at
200/7700 for quick8; early resolved/structured inputs use zero curves.
Quick8's 126 handoffs all report schedule exhaustion and all execute SIQS.
There are no active/abandoned ECM jobs in these complete timing results;
partial-frontier handoff and cancellation are exercised independently in tests.
Both arms record two proper splitting events on every recursive call and
690400–693301 total work, without a fresh child allowance. Published SIQS
recovery counts are zero in these captures. Actual curve bounds/counts, seeds,
engine outcomes and published recovery statistics remain in full per-call
receipts, rather than inferred from configured tiers.

Largest published fallback workspace estimates are 43,722,520 bytes for the
smaller control, 44,211,480 for quick8, and 118,665,224 for the larger control.
Configured parent reserves are at most 13,724,672 bytes for control and
3,296,256 for quick8, within the shared 288-MiB owned cap. Maximum reported
process RSS is 309,870,592 bytes (295.52 MiB); this includes JIT and earlier
arms and must not be represented as owned workspace or a per-arm saving.
Maximum checkpoint JSON sizes are 12912/8637/14054 bytes respectively.

## Default allowance and resume diagnostics

Fifteen single-run 2M-work/2-second resource probes and three restore checks
finish in 7.006 seconds under the separate 600-second cap. These are censored
functional diagnostics, not accepted speed estimates; the two-second probe
clock does not redefine the production 30-second default.

At 30 digits, control completes all five historical cases under 2M work.
The selected quick8 service configuration completes only the power: four
others return `insufficient_fallback_work` before fallback starts, because its
500M admission floor cannot fit the caller's 2M grant. This includes a smooth
p−1 case that control resolves cheaply. The reserve deliberately prevents
optional work even when an unfunded floor makes that choice unhelpful; it
must not be recommended at the unchanged default grant. At 40 digits, the
retained control completes three of five probes; balanced and close inputs
stop inside ECM at 1,997,866 work, without executing fallback.

| Restore case | Initial reason / stage | Initial work | Cumulative work | Cumulative wall / CPU seconds | Result |
| --- | --- | ---: | ---: | ---: | --- |
| Smaller control | complete / terminal | 1,975,045 | 1,975,045 | 0.309839 / 0.306684 | complete |
| Smaller quick8 | insufficient fallback work / SIQS | 2,814 | 65,405,054 | 0.695847 / 0.684461 | complete |
| Larger control | work limit / ECM | 1,997,866 | 11,375,161,380 | 4.845417 / 4.818323 | complete |

The explicit enlarged **total** grant is 10^13 work / 30 wall and CPU seconds;
it preserves earlier expenditure and assignment. The terminal control restore
adds elapsed validation only, with no new search. Quick8 resumes its pending
handoff without restarting pretests. The larger legacy control continues its
saved ECM frontier before SIQS. New-policy live SIQS, both SSS modes, cancelled
ECM and recursive cumulative allowances have additional independent tests.

## Verification, limits and reproducibility

Working-tree checks pass **593 PyPy tests**, including 22 focused C3 tests,
and full Ruff/format/pycodestyle lint. A `git archive` of committed source and
inputs at `5feb17ab9c88eb78fd54e1f8624f898146a79b30` passes the same **593
tests** (46.986 seconds) and full lint; **144 benchmark modules** import,
**34 certified corpora** validate with their own formats, all six C3 bundles
construct in both bands, and the 71 C3 / 149 A7 source pins and required
catalogs load. The external PyPy/GMP/lint toolchain is installed separately;
no ignored source, corpus, control or evidence is needed by those checks.
The first extra loader sweep used the phase-two 16-bit/Pocklington validator
on A10's Lucas/32-bit-leaf corpus; the corrected sweep uses A10's unchanged
independent verifier. No source/corpus bytes or proof conditions were weakened.

All proper-divisor boundaries, multiplicities, exact reconstruction and
probable/proven labels retain their contracts. Expected certificates prove
validation truth; they never enter dispatch or upgrade returned labels. Tests
cover atomic reservations and cancellation precedence, no unused context
setup, partial curve frontiers, insufficient admission, real live fallback,
child budgets, deterministic resumed assignments, policy/config corruption,
CLI overrides and genuine pre-C3 checkpoint compatibility. `v1/` is unchanged.
The final documentation commit changes no validated executable source.

| Evidence | SHA-256 |
| --- | --- |
| `training-v5.json` | `bed2e3f317a5d8d233d6da5cc26b03f762f88af5a19df3f66ddf81726671f876` |
| `confirmation-v5.json` | `c94f7bd59d4a3d31a9fd561eb63db365a21814bbe5f86d1f270e76870c24bf2c` |
| `analysis-v5.json` | `e7370b2e41567420dc7888916e5ee7d508a2c6b64d66841deb28e81eda47f145` |
| `diagnostics-v5.json` | `ccf25f14f98df0899f389e0983e6e917d5ab2f8760e76b82e95f42cb74b74611` |

The input-subject interval is conditional on this small, deliberately mixed
and generation-biased corpus. Single-input classes are regression controls,
not population estimates. Matching the user seed does not force the same
SIQS assignment seed after policies consume different ECM RNG draws. The
experiment includes engine preparation, unsuccessful searches, recursion,
packing and independent output validation; configuration/certificate decoding
is outside both arms. Cold startup and instrumented profiles are not accepted
performance evidence here. No C3 profiler was run.

Published `fallback_owned_peak` values are maxima of reported workspace
estimates at job exit, not measurements of the true internal peak. The
288-MiB cap bounds configured simultaneous owned workspace; process RSS includes
both loaded implementations, PyPy/JIT and previous calls. RSS is a process
high-water, not a per-arm allocation or an OS-enforced cap. SIQS publishes
recovery counts; ECM's completed-event telemetry does not expose its complete
replay counters, so missing ECM data are unreported rather than zero recovery.

The pretest/campaign API is usable under explicit caller resources. Admission
floors never promise completion or create work/time. Keep the library's 2M
work default and the CLI's existing grants. Broad factor-size priors, economic
input-dependent deeper campaigns, all 31–35-digit inputs, general 50–99-digit
dispatch, dynamic stage-two geometry and same-curve extension remain unproved.
C1's 50-digit 120-second censoring and separate 141-second resume witness have
different configurations/allowances and do not define a universal cutoff.
G1 owns broader allocation; E1/H1 own combined fresh integration. C9/C10 retain
collector/configuration economics; this branch changes none of their engines.

The original implementation/recovery freezes and historical inputs remain
immutable. Required controls, baseline bytes and independently certified
confirmation corpus are committed. Raw per-call captures, research downloads,
logs, analysis and verification receipts stay ignored under `results/c3/`.
Recovery of the original interrupted training run additionally needs those
optional historical raw receipts; clean reproduction of training/confirmation
and all tests/imports uses committed inputs only. The managed worktree is kept
attached, and raw evidence is also copied into the main checkout's ignored
`v2/benchmarks/results/c3/` tree.

From the repository root, using PyPy 3.11 with the optional GMP test dependency:

```sh
make -C v2 test PYTHON=.venv/bin/python
make -C v2 lint
v2/.venv/bin/python -B -u -m v2.benchmarks.ecm.c3.c3_study confirm \
  --output v2/benchmarks/results/c3/repeat-confirmation.json
v2/.venv/bin/python -B -m v2.benchmarks.ecm.c3.analyze \
  v2/benchmarks/results/c3/repeat-confirmation.json \
  v2/benchmarks/inputs/controls/c3_selected.json \
  v2/benchmarks/results/c3/repeat-analysis.json
v2/.venv/bin/python -B -m v2.benchmarks.ecm.c3.diagnostics \
  v2/benchmarks/results/c3/repeat-diagnostics.json
```

Coordinate the machine-wide window for every benchmark command, including
short metadata/analysis commands. Measurement/diagnostic runners acquire the
lock themselves. Use a new output name; captures and fresh corpus cannot be
replaced. Retain the committed selected policy when reproducing confirmation;
never reselect from revealed fresh inputs. Production changes invalidate the
source pins and require a separately named protocol.
