# C3 complete-call comparison procedure v2

This additive procedure is frozen after the fitted table at ed38852 and before
any round-two comparison observations. It supersedes only the unrun comparison
measurement procedure in round2_train.py. The completed 2052-call calibration,
fit, policies, inputs, seeds, objective and round2_report.py assessor stay frozen.
The original v1 control and raw calibration hashes remain unchanged.

Session CPU costs varied substantially during the serialized C3/C9 work. A
matched three-arm group must therefore stay in one acquired process/window.
Reserve all three full service grants plus one second per arm before starting:
18 seconds for the 30-digit band and 93 for 40 digits. Rotate arm order by
sample/subject and reverse alternate samples, as in v1. A remaining 50 seconds
cannot admit even the first arm of a 40-digit group.

The assignment remains 57 revealed subjects x nine assigned seeds x three arms
(control, optimized fixed32, fitted), exactly 1539 measurement calls. No new
candidate or parameter adjustment is permitted. Every call uses the unchanged
training implementation, validates proper divisors and reconstruction, and
retains the 5/30-second CPU/wall and 10^13-work ceilings, frozen C1 settings and
288 MiB workspace. These observations compare complete recursive factoring;
root-curve calibration projections are not substituted for actual outcomes.

Every acquired process validates at least three seconds of both CPU and wall
warmup for each arm in both modulus bands. It additionally warms forced no-ECM
SIQS for the historical and current implementations in each band, so an ECM
hit cannot leave relation collection cold. Use the fixed r2_*_balanced_0
subjects and seed 17. Each of these ten warmup tasks prepays 60 seconds. Check
the 25-second continuation threshold between full service calls; a final
30-second call plus capture fits this grant. Incomplete warmup output fails the
phase. Warmups, setup, validation, captures, interruptions and final analysis
all count toward the unchanged cumulative 3600-active-second comparison cap.
Cold process startup precedes the lease and is excluded from warmed timings.

Measurement leases are 60..1200 seconds. Leave two seconds for receipt and lock
cleanup. Preserve 183 seconds of the phase cap for the final assessment lease
(180 seconds analysis plus three seconds overhead). Each admitted group/warmup
has a hard wall timer covering validation and capture as well as arithmetic.
Per-call arithmetic retains its cooperative CPU/wall limits. Exceeding a hard
grant is a failed/inconclusive study, never a fast successful observation.
Keep the shared machine lock and explicitly coordinate each C9/C3 handoff.

Before work begins, atomically replace the cumulative ledger with the whole
next grant prepaid. Normal completion replaces that grant with elapsed wall
expense. A killed process retains full prepayment. Fully completed groups may
resume after their assignment, source identity, same-lease identity and all
row hashes are checked. A partial/interrupted matched group or assessment
ends this protocol inconclusively: no partial-arm reuse, retry or favorable
observation deletion. An interrupted warmup can be reacquired with its prior
grant retained and a new full grant charged. Correctness/capture exceptions
permanently fail the phase. No extra measurement calls or phase extension are
allowed. A result file is valid only with a nonfailed completed ledger.

The final wrapper checks all 513 complete group manifests, rejects extra groups,
and invokes the unchanged paired_groups/assess functions under a prepaid
180-second wall grant. The same 10,000 clustered paired bootstrap draws, seed,
weights, positive CPU interval against both controls, wall ratio <=1.05 and
class completion-loss <=5 percentage-point gates apply. If the phase cannot
complete or fund assessment, record an inconclusive result and retain control.
No inference follows from a partial assignment. Eligible training results
still require a separate fresh confirmation freeze and acceptance experiment.

The additive c3_round2_comparison_v2.json pins this runner, its substantive
tests, this procedure, the fitted table and both earlier control manifests.
The runner requires those exact files to be committed before measurement.
Generated group directories, receipts and reports stay ignored.

    pypy3 -B -m v2.benchmarks.ecm.c3.round2_compare measure CAPTURES --lease-seconds 1200
    pypy3 -B -m v2.benchmarks.ecm.c3.round2_compare report CAPTURES --report OUTPUT.json
