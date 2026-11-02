# C10a-07 face and queue terminal results

Subsequent decision: [combined production defaults were adopted](c10a07_default_promotion.md).
The experimental findings, original recommendation and limitations below
remain unchanged. Evidence has moved to
`v2/audit/archives/c10a07-20261101/`; numerical metadata records its preserved path.

The bounded comparison is complete: **42 timed calls**, with all nine fresh
50/55/60-digit inputs factoring correctly in every arm. Queue reduced observed
complete wall cost by **22.4%** and CPU by **22.4%**; face+queue reduced them by
**23.1%** and **23.0%**, respectively. These are exploratory estimates, not
accepted default improvements: the first process lacks final drift evidence.
Nominate **queue alone** for production integration and clean fresh confirmation;
the combination adds only0.86% observed wall benefit beyond queue, without an
established incremental benefit. Defaults remain unchanged; E1 stays deferred.

[Exact numerical summary](c10a07_panel_results.json) contains input-level costs,
certainty/unresolved factors, counter and evidence identities, uncertainty,
resource accounts and archive hashes. No old quota grid, calibration or queue
prefix campaign was repeated. The two unstarted60-only screens were superseded
for execution; their freezes and accounts remain intact.

## Inputs and stages

Nomination used one exposed60-digit input, three calls: control1370.551s,
complement1323.801s, swap1305.662s. All completed exactly. Selection of
`swap_high` was committed as33bea9f before fresh generation. Those three calls
are one independent input, not three, and supply no acceptance inference.

Twelve fresh moduli were generated and independently proof-validated before
measurement: three each45/50/55/60, two balanced and one approximately35:65
factor-bit split per size,24 distinct ordinary prime factors. Digits, achieved
bit lengths, residues, exact products and exclusion against2237 prior moduli /
3344 factor entries were checked. This uses the existing Pocklington-conditioned
generator, not uniform RSA sampling. SIQS receives only N and seed; certified
input factors do not upgrade its literal probable-prime output labels.

The completed schedule contains39 fresh observations on10 independent inputs:
36 successful active-band calls on nine inputs, plus three unresolved45 calls
on one input. All three45 arms stop with native `stalled_yield` in under0.2s;
this is the configured finite stop, not time or memory exhaustion. The frozen
futility rule leaves the other two45 inputs / six cells unvisited. Queue is
bypassed at45; combined is an exact face alias there and was never duplicated.
There are no lost/censored timed cells, overshoots or selective timed reruns.

Integration checks, nomination, committed fresh generation, matched observations,
terminal analysis and archival are complete. The parent's native confirmation
report was not invoked without its prerequisites. A separately identified
continuation completed only unvisited cells. Production integration and adoption
confirmation remain conditional, separately frozen work; no additional campaign
was launched.

## Complete cost and uncertainty

Each primary queue/combined comparison has nine fresh independent inputs,
three per active size. Ratios are input-cluster bootstrap estimates within size
bands, with equal band weighting,50000 replicates and frozen seeds. All active
controls already complete9/9: these are cost observations, not completion gains.

| Arm versus current control | Wall ratio | Adjusted two-sided interval | CPU ratio | Completion |
| --- | ---: | --- | ---: | --- |
| Queue | 0.7761 | 0.7430–0.8184 | 0.7756 | 9/9 versus9/9 |
| Face + queue | 0.7694 | 0.7392–0.8099 | 0.7696 | 9/9 versus9/9 |

These99.17% intervals apply the frozen alpha0.05/6 allocation. They describe
input sampling conditional on the recorded clocks; they do not cover missing
temporal-drift evidence. Both arms meet the numerical wall/CPU bounds, with no
control-only completion. **The original adoption gate does not pass.**

Face's prescribed primary population has only10/12 observed inputs, including
the unresolved45 input. Its incomplete finite-budget service-cost estimate is
inconclusive and is not an overall factoring speedup. Active-band summaries
below are descriptive, three inputs per row; a negative reduction means slower.

| Digits | Face wall reduction | Queue wall reduction | Combined wall reduction |
| --- | ---: | ---: | ---: |
| 50 | -3.23% | 15.28% | 15.65% |
| 55 | 0.09% | 26.47% | 26.97% |
| 60 | 10.24% | 24.97% | 26.06% |

The first60 face-only call's29.2% timing difference did not replicate: the other
two60 face differences were1.8% lower and4.0% higher cost. The combination's
incremental wall ratio to queue is0.99139 (0.86% lower), CPU0.99234. This does not
establish that the face should accompany queue. No component percentages or
historical-study percentages are added.

Input-level wall seconds, ordinal within each size:

| Digits / input | Control | Face | Queue | Combined |
| --- | ---: | ---: | ---: | ---: |
| 50 / 0 | 23.394 | 23.692 | 20.598 | 20.371 |
| 50 / 1 | 38.402 | 40.248 | 32.829 | 32.951 |
| 50 / 2 | 47.872 | 49.607 | 38.678 | 38.458 |
| 55 / 0 | 138.561 | 137.484 | 99.227 | 99.824 |
| 55 / 1 | 100.256 | 99.260 | 74.947 | 72.721 |
| 55 / 2 | 112.332 | 114.038 | 83.428 | 83.722 |
| 60 / 0 | 947.671 | 671.343 | 810.273 | 795.553 |
| 60 / 1 | 610.513 | 599.613 | 442.456 | 432.514 |
| 60 / 2 | 626.836 | 651.614 | 427.291 | 425.999 |

## Repair, drift and correctness

After12 fresh calls, the first process exited during45 warmup: its whitelist
expected `stalled`, whereas native reference SIQS returns `stalled_yield`.
This was an avoidable harness preparation error, not a factoring correctness or
resource failure. A finite two-prefix audit confirmed the native stop after80
blocks in0.120s. Its40-second reservation and the failed warm's20s /1e13 work
reservation remain charged. The failed warm's actual clocks were not retained.

Repair source/protocol/tests a46e482 and control32e73eb preserve the old freeze.
Only warmup recognizes the valid native outcome; timed results remain literal.
Real45 readiness then passed all three arms, >=3 validated wall AND CPU seconds
with the original last-three stability rule. Every new measurement process
rewarmed. No saved timed cell was repeated and no cap or arithmetic gate changed.
The continuation makes sentinel starts durable and attempts ends on exit.

| Continuation band | End/start wall | End/start CPU | Drift gate |
| --- | ---: | ---: | --- |
| 45 | 1.2195 | 1.2170 | Fails1.20 limit |
| 50 | 0.9502 | 0.9584 | Pass |
| 55 | 0.9860 | 0.9852 | Pass |
| 60 | 1.1379 | 1.1193 | Pass |

The stopped parent has no end sentinels. They cannot be reconstructed from a
new process. Continuation-only active evidence has six inputs, not the required
nine. Neither omitting the old observations nor silently treating missing drift
as passed supports adoption. Retain the exploratory evidence and defaults.

Every reported divisor is proper and every result reconstructs exactly,
including unresolved cofactors. Independent final auditing checks42 receipt
hashes, labels, ownership and unique measured cells. All18 completed
control/queue and face/combined pairs match exact retained-state signatures and
candidate/coverage counters. Some original automatic integration comparisons
addressed top-level instead of nested counters; the independent audit checks
the actual nested values and finds no discrepancy. No executable was patched
inside either frozen identity.

Observation is disabled in timed calls. Queue preserves positions, root/power
hits, candidate refinement, admitted atoms and retained relations while avoiding
repeated empty sparse-prime visits. This supports the throughput mechanism, not
a better useful-relation yield. No face-level causal useful-dependency attribution
is available. Refined rejections precede ordinary prime division; high counts
justify no threshold relaxation. Mixed dependency touches overlap. Summed raw
plan/queue byte gauges are not resident peaks; the checked helper high-water
field and process RSS are reported separately.

## Identities, resources and archive

Runtime source remains mainline **e520ec9**, integration/protocol9757105 and
control8abbf09; selection33bea9f and corpus/order f53a294 precede fresh timing.
Mainline observed at reporting is090f8e5, differing only in public roadmap/study
documents; factoring runtime is unchanged. Neither historical92543797 controls
nor the later A1a adoption were substituted into old measurements.

Exact interpreter:
`/opt/homebrew/Cellar/pypy3.11/7.3.23_2/libexec/pypybin/pypy3.11`, SHA-256
`3c45fc49fa01571d19bac1594495d7300c30589c79c0cad2c5bee8d59f9e5ef2`.

- Parent freeze: `2d2564a0948411b7a3d144ff5178ed0d2d266d7109c00d1f5ab30c3679a1c5a2`.
- Continuation freeze: `3ed32cf987d2bb91f533ffc8beacbe83227a6bf5095b9c97a242ea4ffbf21ae8`.
- Fresh corpus: `c4e844fff20506fa852baa3a0de6302b040cf3ad069978ba3cba073fbcab40f3`.
- Schedule: `095d7909053ab85a8510e30b8f1e35b44cb3cad083b9cfe47ae0990b17310ee2`.

Recorded completed-call consumption, including checks/warm/generation, is
**13312.083 wall /13279.867 CPU seconds** (about3h42m wall). The lost warm actual
adds between0 and20 seconds to each. Conservative cumulative charges are
45430 wall/CPU seconds and14308059871447 work units;3950 wall/CPU seconds and
1705691950128553 work remain in the original accounts. Remaining phase allowances
are80s checks,3120s warm and750s unvisited confirmation; none is borrowed or spent
to fill a cap. All measured confirmation resources stay separate from exploration.

Committed QA recorded145.860s across the retained old, new and final source
verification runs, plus under5s
of local focused checks, within600s. Native analysis took1.242s within120s;
archive creation/verification0.716s. Fresh process RSS high-water379.1 MiB overlaps
arms; maximum reported native workspace580.6 MiB and checked queue helper peak
16.0 MiB remain below1 GiB total logical ownership and3 GiB RSS. Raw outputs use
1085416 bytes within512 MiB. No production memory exhaustion occurred.

Durable local archive:
`factor-c10a07-panel-continuation/v2/audit/results/c10a07_panel_continuation/terminal-evidence-v1.tar`,
26265600 bytes,916 source/evidence files verified individually. SHA-256:
`8635604d3d85ea8d1ea38c97499d7eae8e0994d71e1bd335cc7f017b72c4f607`;
manifest SHA-256:
`e96a6b41ffca2db51ad5f129595ea1dbe802a2b291640521b3e4408c9af97ca4`.
Receipts, ledgers, failed reservation, consoles, audit, source/control/corpus and
QA remain ignored locally and archived; none was force-added.

Committed-only32e73eb source passed857 PyPy tests (nine skipped), full lint,
215 benchmark imports and strict loader verification with the explicitly
retained parent evidence. One existing read-only `ps` test initially hit sandbox
permission denial and passed on its isolated rerun. No suite was selectively
repeated to improve performance data. Source-only tests/imports require committed
files; resumption provenance additionally requires the archived receipts.

To reproduce analysis, restore an archive copy and the pinned source/control
files at the documented worktree paths, then use the exact interpreter with
`-m v2.benchmarks.qs.c10a07_panel_continuation.study verify` and `report`.
The committed protocol supplies full execution/recovery commands. Preserve the
original archive; analysis on a copy may have a different elapsed-time field.
Do not rerun any measured stage to reconstruct an unfavorable observation.

## Terminal decision and next priority

Historical quota and face studies remain scoped inconclusive with no accepted
clean speedup; their shared three inputs cannot be pooled. The [retained recovery
report](c10a07_recovery_results.md) records their recorder failures, limits and
charges. This new panel supplies direct clean comparisons under its own identity.

The bounded C10a-07 panel and reporting are complete. **Nominate queue alone;
retain defaults and conclude adoption inconclusive under the original gates.**
The added face is not justified for adoption by its small incremental result.
A checked production port and separately frozen fresh bridge are the next useful
step if adoption is pursued, not another diagnostic or repeat of old grids.

The [reconciled follow-up design](c10a07_recovery_followup.md) proposes18 calls on
nine fresh inputs total across50/55/60, not nine per band. It excludes every
current input/factor and requires intact drift evidence. The24-input60-digit
balance/residue panel with a second sampling family remains a later conditional
generalization design. Explicit61/62 feasibility precedes any justified move
toward65; no larger sweep is scheduled. Special forms and portfolio behavior
remain separate. The [roadmap](../ROADMAP.md) records A1b/A1c and the other existing
owners with model/effort recommendations; A1a remains delivered. No new default,
merge, push, automatic follow-up campaign or E1 promotion occurred.
