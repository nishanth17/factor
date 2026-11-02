# A1b existing fixed-score decision

**Final status:** [gap closure and all-band decision](a1b_gap_closure.md)
now supersedes the pending status below. All 30 planned calls were measured;
retain defaults in every screened band. The original report is preserved.

**Subsequent prompt audit:** [requirement-by-requirement review](a1b_prompt_audit.md)
confirms actual 24-call prefix experiments but withdraws the broad completion
claim. The 55-digit screen, exact-composition evidence mapping and replacement
analysis accounting have unresolved gaps. The observations and historical
report below are preserved; they do not establish full prompt compliance.

**Retain current integrated SIQS defaults.** Existing fixed/list scoring loses
the completed 50- and 60-digit prefix screens; defer 55 digits because readiness
stopped before measurement. Fixed scoring remains opt-in. A1b's existing-option
review is terminal; no pilot, fresh campaign, default change or E1 work follows.
These are exploratory prefix findings, not population or complete-factor
slowdown claims. [Numerical results](a1b_results.json) retain input ratios,
counters, caps, failures and archive identities.

## Control, reconciliation and eligibility

Base **4a946b8d6294d51ac5d9c36f61880cc63d5b59fd** includes committed face/queue
integration 0946824 and historical-control fixes 308a4f1. Task branch is
`codex/a1b-fixed-score-review`, isolated in `factor-a1b`; C10a-01 changes are
excluded. Selection 268ca7a precedes candidate observations; runner/source freeze
ac70718 precedes timing. Reporting-only repair 106bf04 preserves the old freeze.
No uncommitted source or historical pre-queue control was substituted.

R2 isolated/integration baselines match their frozen hashes, including all 35
embedded source files each. Archived training, held-out and integration captures
were checked against their archival manifest. R2's fixed+plan nominal 30-digit
held-out benefit was about 5.68% versus its then-current arm; training missed the
then-frozen 10% gate. Preserve that retain-defaults verdict. The revised policy
has no universal 10% floor, and old conditional gains establish no current benefit.
C10's accepted 2 MiB power plan is held fixed, not retested or credited to A1b.
A1a 50 integration passed;55 drift remained inconclusive and the power-period /
unroll bundle lost 50/55. C10a-07 combined adoption is the current control;
its missing historical drift evidence and nonadditive percentages stay scoped
to that study. See [A1a](a1a_frozen_results.md), [C10a-07 results](
c10a07_panel_results.md) and [adoption](c10a07_default_promotion.md).

Actual caller paths, rather than configuration labels, govern eligibility:

| Existing option | Effective path and disposition |
| --- | --- |
| powers/list control | Queue eligible for aligned bounded windows, bucket recovery, zero cutoff/extra and supported roots. Dense marking plus queued sparse prime/power hits; exact native fallback on refusal. |
| fixed/list + current plan | Conservative scale 32 integer log bounds and tighter value thresholds; queue immediately falls back to direct marking. Keeps 50/55 density+presieve+scan_fused,60 sparse/scalar scanning. **Sole measured challenger.** |
| fixed/array | Matching score/threshold saturation at 2^32−1; direct fallback, per-hit saturation checks and no list density/presieve. No credible distinct advantage over fixed/list selected for measurement. |
| fixed/bytearray | Matching saturation at 255 preserves admissible coverage while allowing more false positives. Scale 32 makes filtering weak; direct fallback, translated slices and no list presieve/density. Rejected as nominee on mechanism, not a measured current slowdown. |
| historical fixed/no-plan, fixed+plan and marking labels | fixed+plan maps to the candidate with today's accepted plan. No-plan repeats the distinct reuse question; ordinary dense/sparse/bucket labels largely share the powers branch with bucket recovery. Do not promote aliases or add the losing A1a extension. |

Both policies support capped exact power plans. Presieve reserves periods and
scratch; fused scanning supports wider fixed-score lanes. Backends retain
conservative owned reservations, not an automatic packed-buffer discount.
Snapshotv 7 saves score/backend/kernel/plan/queue selectors and charges derived
rebuilds; incompatible score-policy substitution is rejected. Independent R2
integer-power/division oracles and current tests cover clipping with matching
threshold clipping, high valuations, p=2, A/multiplier divisors, signed windows,
tails, zeros, refinement/recovery, refusal/cancellation and charged restore.
New tests add fixed scoring with today's queue selector, presieve and scan.

Detailed matrix, provenance and focused literature/source reconciliation remain
local in `v2/audit/a1b/`. Pinned msieve 1.53 and Yamaquasi 270a6617 distinguish
integer log marking from sparse scheduling; native compact loops do not imply
a PyPy fixed-score win. [Bradford–Monagan–Percival §2.4](
https://www.daemonology.net/papers/maplent.pdf) describes statistical small-prime
allowances followed by per-hit refinement; this is not Factor's exact coverage
contract. Existing cutoff/refinement, powers, presieving and fused scans already
cover the relevant implemented mechanisms. CUDA-MPQS's audited sampled
saturation predicate supplies no global norm-bound theorem. No broad source
review or upstream benchmark was repeated; no source was copied into runtime.

## Remaining opportunity and diagnostic limits

One exposed 60-digit control, seed 7,2146 blocks: clean 2.4121 wall /2.4105 CPU
seconds, observed 2.7141/2.7122, identical work/counters. Observation adds 12.5%.
Exclusive bounded-wrapper wall clocks are:

| Observed category | Seconds |
| --- | ---: |
| Direct marking, prime powers/plans and block threshold | 0.7566 |
| Sparse queue construction | 0.7528 |
| Queue replay and dispatch | 0.2275 |
| Scanning, cursor/polling and collector glue | 0.2879 |
| Exact candidate refinement | 0.0295 |
| Recovery, including norm/residual processing | 0.1907 |
| Relation admission | 0.0107 |
| Outside wrappers: setup/switching, downstream and unattributed observer cost | 0.4586 |

Do not add overlapping native stage clocks to this table. Marking/lifting/plan
work and block thresholds are inseparable here; scanning includes result
construction. Native setup/root stage clocks are inclusive context, and this
prefix spends only~0.000046 s in preparation/filtering, exposing little downstream
solve work. Observer overhead cannot be assigned cleanly among categories.
For a hypothetical workload-preserving isolated elimination, affected observed
fractions give approximate ideal ceilings:27.9% for direct marking,11.7% for
scan/glue plus refinement, on this prefix only. The candidate changes both
candidate workload and scheduling, so these are **not ceilings on its bundle**
or whole-job/population theorems. No historical queue percentage was subtracted.

The hypothesis was that tighter scores might save enough refinement/recovery
to repay direct fallback. It did reduce division work, but coarse candidate
counts barely changed: mean 547036→539038 at 50 and 536312→534051 at 60.
Division-prime visits fell 116407→55008 and 36033→15325, respectively. Meanwhile
control avoided~2.91/10.12million sparse prime visits per prefix at 50/60;
fixed bypassed queues and its unchanged 2 MiB plan encountered mean 7818/230457
plan refusals, with exact streamed fallback. These are bounded cache refusals,
not production memory exhaustion or recorder failures. This identifies a
bundled integration tradeoff, not an isolated arithmetic penalty.

## Screen and advancement

The [frozen protocol](../benchmarks/qs/a1b_protocol.md), [literal configurations /
schedule](../benchmarks/inputs/controls/a1b_screen.json) and [source manifest](
../benchmarks/inputs/controls/a1b_source.json) precede timing. Every call uses
1 GiB total owned memory:8 MiB runner/observation reserve plus 1016 MiB native,
carved out **after** original 1 GiB default selection;3 GiB RSS and 512 MiB generated
output/evidence. Both arms retain 2 MiB plans,32 MiB queue selector, multiplier/A
policy, quota, swap_high, interval, LP domains and pacing. The only literal
change is `score_policy=powers`→`fixed`, both list/bucket. Queue fallback and
changed plan population are unavoidable bundled differences.

Each matched call covers 2146 blocks /8523841 positions /66 polynomials,
including family changes, with setup and exact output validation included.
Scanned positions, boundaries, final polynomial/family/cursor, root/power-hit
counts and admitted-atom counts match all 12 pairs; retained-atom hashes also
match. Candidate streams need not match. All 24 measured outcomes are native
`paused` with exact unresolved N; none is a complete factorization. Native
clocks excluding extra validation are also retained.

The frozen rule requires 2/3 independent wall wins, input-level geometric wall
ratio<1 and CPU<=1, exact coverage and readiness/stability. At 60 each input
uses its median ratio over seeds 7/29/47; lower bands use seed 7. No minimum gain
or significance threshold, and no favorable retries.

| Digits | Measured inputs / pairs | Input wall ratios | Geometric wall / CPU | Advancement |
| --- | --- | --- | --- | --- |
| 50 | 3 /3 | 1.2576, 1.2862, 1.2764 | 1.27334 /1.27340 | Retain:0/3 wall wins |
| 55 | 0 /0 | Unmeasured | — | Defer:readiness incomplete |
| 60 | 3 /9 | 1.8185, 1.8313, 1.8409 | 1.83022 /1.83029 | Retain:0/3 wall wins |

Complete-prefix wall sums control→fixed are 5.160→6.568 s at 50 and 21.856→40.119 s
at 60; CPU 5.158→6.566 s and 21.847→40.080 s.60's raw pair wall ratios range
1.8095–1.8727. Repeated seeds remain three independent moduli, not nine.
These descriptive ranges and clear observed losses justify shortlisting
rejection; there is no fresh confirmation, reliable population slowdown interval,
separate reliable per-size adoption claim or complete-factor cost inference.
No band supports a pilot; Stage B/C were not run.

## Stability, failures and accounts

Every measured arm/process accumulated>=3 validated wall AND CPU warm seconds
and passed last-three 20% range/median readiness.50 end/start wall ratios are
0.977/0.955;60 ratios 1.060/1.018, with CPU likewise within the frozen[0.8, 1.2].
No measured call was lost, timed out, refused, censored or selectively retried.

The warmup guard kept six end-observation slots reserved even after earlier
ends had settled. It stopped 55 fixed after two warm calls (<3 s and <3 calls),
despite 50 seconds remaining in the phase; this is an overconservative harness
allocation error, not candidate infeasibility or a native resource failure.
55 control has a durable start/end; fixed never became ready, so its missing
start/end remains missing. Rewarming another process cannot fill that evidence.
No 55 comparison was run and no cap was enlarged.

The diagnostic's first 2.78-second warm failed readiness; its receipt and 30-second
charge remain. The separately identified diagnostic continuation accumulated
warmup correctly and charged 60 seconds. The original terminal analyzer then
failed `len(generator)`, preserving its full 60-second reservation. Separate
reporting identity `a1b-terminal-audit-v1` materializes the same finite ratios
and applies unchanged formulas/rules as an ordinary runner-regression/reporting
QA repair (0.024 s); it schedules no new experimental call or reservation.
Original executable freezes, all failures and their charges remain untouched.

| Stage A account | Charged wall AND CPU seconds | Remaining |
| --- | ---: | ---: |
| Measured prefixes | 750 | 150 |
| Process-local warmup/end observations | 250 | 50 |
| Current-control diagnostics | 90 | 30 |
| Analysis, including failed reservation | 60 | 0 |
| Orchestration reservation | 90 | 0 |
| Total | **1240** | **230** |

Observed factoring/warm/diagnostic consumption is 119.019 wall /118.936 CPU
seconds, excluding separately bounded development/QA/reporting and unmeasured
orchestration overhead. Conservative work charges total 93955256878 units;
per-call/cumulative work/time caps remain finite.23,170-second whole-proposal
ceiling has 21,930 unspent; unused phase funds authorize no automatic follow-up.
Maximum observed native ownership 75,715,392 bytes, queue helper 17,140,480 bytes
inside native ownership, and shared process RSS136,970,240 bytes. Summed raw
plan/queue gauges are not resident peaks. Output and verified archive remain
well below 512 MiB. RSS is sampled process high-water, not a fresh per-arm peak.

Committed-only 106bf04 QA passes **834 PyPy 3.11 tests (9 skipped)**, full
`make -C v2 lint`,196 benchmark imports and source-freeze verification;47.127 s
wall under the separate 600 s QA allowance. Earlier fixture/preparation errors
and logs remain local. Exact interpreter SHA256:
`3c45fc49fa01571d19bac1594495d7300c30589c79c0cad2c5bee8d59f9e5ef2`.

Durable local archive:
`v2/audit/archives/a1b-20261102/a1b-terminal-evidence-v1.tar.gz` in the main
checkout,6,276,192 bytes;108 files individually verified, including committed
source, receipts, ledgers, research, primary references and QA. SHA256:
`72159c7ff09809027bf03945849b04f20a55e81a4ecb5c8466b6c82f567602e0`.
The numerical report includes the manifest hash and commands. Verify frozen
source with `pypy3 -B -m v2.benchmarks.qs.a1b verify`; this is a terminal study,
not an instruction to rerun its timing or failed analyzer. For recomputation on an archive copy, call the read-only
`v2.benchmarks.qs.a1b_terminal.audit()` function and save its returned JSON
to a new destination; do not replace the original terminal receipt.

**A1b-F1 (new, deferred): policy-bound fixed-score sparse queues**, with
independent weight/power/threshold coverage, finite scratch and charged rebuild
proofs. This task does not extend queue eligibility, alter quantization, rescale
bytes or redesign polynomial/recovery representations. Existing fixed/list is
not a useful nominee here;55 remains unmeasured. Retain defaults, preserve
fixed opt-ins, close this bounded review and keep new mechanisms with their
separate owners. No merge, push, E1 reopening or default change is authorized
by these findings.
