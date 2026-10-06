# Factor benchmarks

Validated factoring experiments on **PyPy implementing Python 3.11**. This guide
keeps the stage history, accepted changes and rejected experiments concise.
The [v2 guide](../README.md) covers usage; the [roadmap](../ROADMAP.md) records
remaining acceptance gates.

## SIQS CLI access (4 October 2026)

The initial SIQS comparison below is retained. The P3.4 usability follow-up
now also exposes `--method qs` and `--method mpqs`, with shared finite
`--qs-*` controls and checked resume for all three modes. Nine CLI test methods
cover extraction, recursive composite children, sign/multiplicity, fallback
order, mode constraints, help, interruption/resume and unchanged defaults.

Each selector was compared with its matching existing library mode on
`10002200057`, seed 7 and base bound 400, using the same allowances below.
Every mode received at least three seconds of validated PyPy warmup and nine
samples per arm; the unstable QS cohort was extended with five more warmup
seconds and 31 samples per arm. CLI/library factors, stage seeds, outcomes
and consumed work matched: QS 601,427, MPQS 944,027 and SIQS 1,524,681 units.
Both arms included formatting and checkpoint writes; startup was excluded
from these warmed samples. Local records are in `p34_cli_20261004` under
audit results. This verifies entry-point exposure, with no method ranking or
automatic-dispatch/default promotion claim.

`pypy3 -m v2.factor N --method siqs` selects SIQS after preprocessing;
`--siqs` enables the existing recursive rho/p−1/ECM → SIQS portfolio. CLI
regressions cover actual extraction, composite children, sign/multiplicity,
fallback order and checked resume: `pypy3 -m unittest v2.tests.test_siqs_cli`.

A plumbing comparison on `10002200057`, seed 7 and base bound 400 used the
same 200-million-unit, 30-second wall/CPU and 80-MiB portfolio allowances.
After 3.009 seconds of validated warmup, all nine warmed samples per arm
returned the exact factors with identical 1,524,681 consumed work units and
SIQS seeds. Both arms included output formatting and checkpoint writes.
Nine CLI cold-start checks were separate. Local runner/raw records are in
`siqs_cli_20261004` under audit results. This small library/CLI comparison
validates exposure; it establishes no method ranking, larger-input speedup
or automatic-default promotion.

## Run and interpret

From the repository root:

```sh
make -C v2 benchmark WARMUP_SECONDS=3 REPETITIONS=9
make -C v2 benchmark-phase-two WARMUP_SECONDS=3 REPETITIONS=9
make -C v2 benchmark-phase-three-siqs WARMUP_SECONDS=3 REPETITIONS=9
make -C v2 benchmark-phase-three-sss WARMUP_SECONDS=3 REPETITIONS=9
make -C v2 benchmark-phase-three-parallel WARMUP_SECONDS=3 REPETITIONS=9
```

Make creates unique timestamped captures in `results/`; override
`BENCHMARK_OUTPUT` with a path relative to `v2/`. Direct runners require an
existing output folder and a new output name. Use `--help` for their parameters.

Recorded results below use PyPy 7.3.23 / Python 3.11.15 on macOS arm64.
Warm evidence includes at least three seconds of validated warmup and nine
samples; unstable series extend warmup and sample count. M9 uses 15 samples.
Times cover whole cohorts unless stated otherwise. Compare matched inputs,
seeds, budgets and validation; report completion and stop reasons with timing.
Cold startup and instrumented profiles are separate evidence. Owned workspace
is distinct from process RSS. Faster refusal is not a time-to-factor gain.

## 3 October 2026 — correctness and early performance

### Phase 1 / M9: repaired v2 versus v1

| Complete batch | Emulated v1 | Repaired v2 | Less time |
| --- | ---: | ---: | ---: |
| Five inputs × five seeds | 1.155 ms | 1.078 ms | 6.6% |
| Separate 56-input control × five seeds | 21.091 ms | 16.994 ms | 19.4% |

Every timed answer has expected factors, multiplicities and reconstruction;
native results also require proven terminal factors. v1 uses `lib2to3` and
integer-division emulation on the same PyPy, not native Python 2. Invalid
baseline outputs receive no ratio. These are historical small-workload results.

Keep classification reuse, exact trial roots/square splitting, wheel-6 sliced
marking and reduced rho bookkeeping. Reject slower alternative GCD/root loops.
A longer segmented-sieve repeat found essentially unchanged performance;
search remains slower than v1's distinct-only shortcut while fixing its contract.

```sh
mkdir -p v2/benchmarks/results
pypy3 -m v2.benchmarks.regressions --legacy \
  --warmup-seconds 3 --repetitions 15 \
  --output v2/benchmarks/results/v1_comparison_UNIQUE.json
```

### Phase 2: bounded portfolio and schedules

Shared work/time/storage, streamed primes, saturation recovery and portable
resume replace uncontrolled retries/allocations. Independent 20-digit
confirmation under matched 50 ms operation caps improved completion from
**70.8% to 91.1%**; total attempt cost fell **39.283 → 16.362 ms**.
That cost includes unfinished outcomes and is not a successful-run speed ratio.
Broader schedule, parameter and competitor promotion remains separate.

### P3.1–P3.3: exact relations and ownership

Reference root/collector/kernel checks establish arithmetic correctness. Shared
postprocessing and ownership refinement reduce a complete 16-input 23–26-bit
cohort **27.220 → 21.446 ms (21.2%)**; ratio interval **0.726–0.863**.
Keep conservative collection and provenance-aware filtering. Resieving/tighter
scoring did not improve that complete-run comparison and remain optional.

## 4 October 2026 — relation engines and repair experiments

### M28: correction audit

214 root sets, 648 collector windows and 320 lifted kernels pass independent
checks. Storage-cap extraction and combined live reservations are repaired.
Matched root/bucket-recovery cohorts cost **6.5% / 3.4%** more than M26;
retain correctness repairs without claiming a speedup. Cold costs are inconclusive.

### M29–M30: families and shared SIQS stores

Root reuse lowers a setup-inclusive utility **1.287 → 0.818 ms (36.4%)**;
whole-family attempt cost is unchanged, with 13/16 completions in both arms.
Later shared stores complete **64/64** small input/seed cases versus **41/64**
for fresh stores. Tuned SIQS/QS/MPQS/ECM cohort times are
**19.7 / 18.6 / 20.2 / 4.0 ms**, all 64/64 complete.

Keep shared relations and validated full checkpoints. Scored multipliers,
recovery growth and simpler families do not establish default wins.
Initial 0.2-second large probes mostly fail and do not close practical large-size
gates; the earlier complete-gate claim was withdrawn.

<a id="completed-larger-evaluation-and-filtering-repair-m31-4-october-2026"></a>
### M31: larger budgets and continuation

Longer declared runs and trained finite schedules replace tiny diagnostic caps.
Repeated matrix-incidence work is repaired; one representative's matched cost
falls **5.920 → 3.425 s**. A fresh balanced 50-digit SIQS run completes in
**1,183.355 s (19 min 43 s)**. Checkpoints support longer continuation.

Keep stop reasons, failed attempts and success subsets separate. One success
is not broad 50-digit coverage; 60–80-digit capped failures remain. The old
“50–60 digits within a minute” promise is unsupported.

### P3.5: SSS and filtered SSSf

CRT collision enumeration, full prime-power smoothness and capped product trees
feed the common verified relation pipeline. Training and held-out comparisons
include QS/MPQS/SIQS, SSS/SSSf, multiplier and recovery controls.
Keep SSS/SSSf explicit; filtered SSSf can miss candidates. Initial overhead and
capacity problems motivate P3.6.1. Upstream comparisons remain capped diagnostics,
not a large-number ranking or automatic-dispatch justification.

<a id="p36-coarse-siqs-workers-4-october-2026"></a>
### P3.6: serial/thread/process workers

Coarse workers preserve central verification, aggregate CPU, cancellation and
assignment identity. Early fixed-work small medians are **630.8 ms serial**,
**1,289.7 ms with two threads**, and **1,359.4 ms with two processes**.
Medium cap refusals and poorer completion prevent worker promotion.

Keep serial defaults. Executor startup, IPC, larger retained batches and delayed
first extraction can dominate. Fixed-work throughput and first-factor latency
are separate experiments; more cores do not establish an algorithmic gain.

<a id="p35p36-cost-diagnosis-4-october-2026"></a>
### P3.5/P3.6 diagnosis

Profiles identify SSS collision/counter work, repeated verification and worker
polling/lock costs. Coarse publication scans more positions before discovering
a factor: a small cohort scans **2,048 native / 4,104 coarse serial / 15,842
four-process positions**. Try bounded polling and smaller result chunks first;
profile time is not performance evidence. Matrix/tree changes need separate gates.

### P3.8 R1: capacity and assignment policies

Frozen training/confirmation corpora each contain 58 disjoint inputs, but the
repeated balanced comparison has **one independent input per size band**.
Fresh 30-digit two-seed medians are **1.484 s flyer SIQS / 1.198 s external-square
MPQS / 0.279 s ECM / 0.277 s ECM then SIQS**, all 18/18 complete.

Flyer improves over nearest by 18.3% on that input; neither this nor
ECM-first completion establishes a population crossover. The fixed legacy
training control costs **1.701 → 1.752 s (3.0% slower)**. Keep defaults and
measured source snapshots; broader populations, Gray quotas and upper-band
capacity remain open. Cold two-seed lifecycles are 2.898 s SIQS / 0.822 s ECM.

<a id="fresh-confirmation-and-decisions"></a>
### P3.6.1: preprocessing, collection and worker repairs

Fresh four-input/two-seed complete cohorts compare frozen control to combined
repairs, including setup and terminal classification:

| Cohort | Before → after | Less time |
| --- | ---: | ---: |
| Small / medium P2 | 1.809 → 1.470 / 6.888 → 5.980 ms | 18.7% / 13.2% |
| Medium native SIQS | 86.587 → 35.436 ms | 59.1% |
| Medium coarse serial | 318.386 → 63.434 ms | 80.1% |
| Medium SSS / SSSf | 2.617 → 0.760 / 6.211 → 1.561 s | 70.9% / 74.9% |
| 30-digit SIQS / SSS / SSSf | 23.324 → 6.840 / 6.181 → 1.412 / 9.460 → 2.020 s | 70.7% / 77.2% / 78.6% |

All 72 timed 30-digit attempts per arm complete. Their 95% timing-ratio
intervals are **0.292–0.298 / 0.224–0.238 / 0.212–0.215**.
Keep sparse recovery, reused inversions, bounded preparation/polling and
compact labels. Reject the disjoint-base SSSf candidate; retain full-base
smoothness and the historical filtered control. Unfiltered SSSf does not
beat SSS by the promotion threshold. A censored 30-digit P2 comparison
completes 7/8 in both versions; its timing is not a full-factor gain.

### P3.8 R3: stable relations, matrix and cadence

Keep stable mixed-row admission, complete payload identities, checked replay
and incremental pivot counts. Reject live-column compaction, pivot batching and
merge histories as defaults; no consistent end-to-end gain appears. The mixed
supplement finds **83/92 verified kernels and 44/32 proper divisors**, unchanged
across arms. Preparation controls show no new cache win; dependency caching
adds overhead and remains off.

| Held-out cohort | Frozen control | R3 cadence 1 | Candidate |
| --- | ---: | ---: | ---: |
| Small | 42.450 ms | 42.205 ms | cadence 8: 41.560 ms |
| Medium | 81.620 ms | 71.483 ms | cadence 8: 64.570 ms |
| Nominal 20-digit | 810.340 ms | 739.495 ms | cadence 8: 575.072 ms |
| Nominal 30-digit | 7.745 s | 7.526 s | cadence 32: 4.611 s |

20-digit completion is **36/47/54 of 54**; other medium/30-digit arms complete
54/54. Cadence 32 passes its 30-digit gate: **26.9% training / 38.7% confirmation**
less complete-cohort time versus R3 cadence 1. Defer standalone cadence 8;
its isolated training gains miss the gate. Production cadence stays 1.

<a id="follow-up-restrictive-allowances-and-repeated-setup"></a>
### Follow-up: restrictive allowances and repeated setup

Fresh independent confirmation preserves exact work, factors and atom signatures.
Serial worker time drops **12.3% small / 12.7% medium / 32.2% with B=10000 chunks**;
serial-pool lifecycle drops 16.7%. Four-thread fixed work improves 79.7% against
its control but still trails serial. No extra complete-SIQS or fresh-process-pool
speedup is established; their timing-ratio intervals include 1.

Clip leases, reuse bounded local coordination/identities and exclude disabled
stage storage. Under matched restrictive allowances, completion changes
**0/8 → 8/8** small and **0/8 → 6/8** medium serial; thread/process lease and
10 MiB worker-cap controls reach 8/8. Resieve capacity reaches 8/8 collected
windows, not eight completed factorizations. Hensel caching, direct single-hit
valuation and other prototypes do not justify defaults. Unstable first-factor
thread results are excluded. Native serial and automatic dispatch stay unchanged.

## Files and retained evidence

| Folder | Contents and retention |
| --- | --- |
| `inputs/corpora/` | 13 versioned independent training/confirmation/comparison corpora |
| `inputs/baselines/` | 10 versioned source controls needed by tests, default runners and replay tools |
| `inputs/controls/` | Six versioned frozen policies, competitor metadata and provenance inputs |
| `history/` | Nine ignored optional candidate/experiment and duplicate measured-source snapshots, with a hash manifest |
| `results/` | Ignored raw captures, generated freezes, checkpoints, stdout and profiles |
| Repository `.local-evidence/` | Preserved historical archives and restoration metadata |

No benchmark evidence was deleted. Optional history remains local; published
inputs shrink from about **9.0 to 5.8 MiB**. The full former guide remains at
`v2/audit/benchmark-history.md`. Keep that guide, `history/`, and the raw archives
when moving machines. Immutable snapshots retain their exact bytes; current-source
runs are new experiments, not relabeled historical measurements.

`performance_followup --variant NAME` uses the local experiment manifest;
`--experiments PATH` accepts an external copy on a fresh checkout. Standard
runner defaults and tests use only versioned inputs. Replay tools may require
recorded runtime/source hashes; use each runner's explicit corpus/baseline
arguments when reproducing a historical candidate.

## P3.8-R2 bounded collector evaluation — 4 October 2026

R2 is isolated on top of the committed repair/R1/R3 integration control
`fba5a34a884fd129dfad5a07bba9f65f5657646f`, avoiding the concurrent P2/P3
repair checkout. The hash-checked `p38_r2_baseline.json` retains that control's
35 runtime modules. `p38_r2_frozen.json` binds the final 36 runtime modules,
three comparison drivers, training corpus and decisions; the fresh held-out
corpus binds that freeze. Required controls and certified corpora are committed;
raw captures, profiles and stdout remain local.

**Decision: retain defaults.** Integer fixed-point power scores and capped
per-polynomial power plans remain explicit experimental options. The public
adaptive/root settings and the configured streamed powers/bucket control stay
unchanged. Sparse hit recovery, cached A support and performed-work charging
come from the preceding repair control. R2 tests their coverage and measures
the additional options independently; it does not claim those repairs as R2
speedups or expand matrix capacity.

Each training band contains two independent balanced inputs; held-out bands
contain three, each run with seeds 7 and 29. Actual sizes are 8, 13, 21 and
30 digits in training, and 8, 13, 20–21 and 30 digits held out. Pocklington
certificates verify both factors independently. The sampler favors primes
with a large known factor of p-1: these are finite comparison cohorts, with
no RSA-distribution or upper-band claim. Complete-call measurements include
setup, recovery, filtering, extraction, primality labels and validation.
Every returned split is proper and every result reconstructs its input.

All arms have work limit 10^10 and 128 MiB of owned workspace. Wall and CPU
limits are each 0.2 seconds per smaller attempt and 2 seconds per 30-digit
attempt; collector-only calls use 30 seconds. Base bounds are 200/1000/3000/
10000 and half-widths 256/512/2048/8192. Blocks are 256 except 4096 for 30 digits.
The 30-digit configuration uses the already accepted explicit R3 cadence 32
in every arm, with row excess 32 and residual bound 10^7; it is no automatic
size policy. Plans reserve a full additional 1 MiB before setup. Every capture
records the complete per-arm configuration, including relation/store caps.

Each arm receives at least three seconds of validated PyPy 3.11 warmup and
nine interleaved samples; noisy comparisons extend to five seconds and fifteen
samples, up to three attempts. Final held-out comparisons pass stability:
90/90 complete outcomes per arm in each smaller band and 54/54 in 30 digits.
Training has 60/60 medium and 36/36 larger outcomes per arm. Small training
remains unstable after extensions and supports no promotion claim.

| Cohort | Training current / fixed + plans | Held-out current / fixed + plans | Held-out causal reduction, conditional 95% interval |
| --- | ---: | ---: | --- |
| Small | unstable | 6.429 / 6.264 ms | 2.6%, 1.3–8.7% |
| Medium | 17.162 / 16.766 ms | 28.496 / 26.886 ms | 5.7%, 0.8–7.1% |
| Nominal 20-digit | 276.733 / 263.307 ms | 272.581 / 233.373 ms | 14.4%, -9.4–17.3% |
| 30-digit | 1.024680 / 0.958610 s | 1.349218 / 1.272612 s | 5.7%, 1.4–8.2% |

Times are entire four-attempt training or six-attempt held-out cohort medians,
not per-input times. Intervals concern paired repeats conditional on these
fixed cohorts, not a population of future inputs. Training causal reductions
are 2.3% medium, 4.9% nominal 20-digit and 6.4% 30-digit, all below the frozen
10% time / 10-point completion gate. The larger held-out nominal 20-digit
median gain has an interval crossing zero and cannot overturn the training
decision. No held-out retuning occurred. Training predates the final
keyword-only compatibility declaration for `power_plan_bytes`; arithmetic,
configuration values and measured kernels were unchanged. Final confirmation,
cold and eligibility captures match the frozen source exactly.

The conservative root-weight, cutoff-5, exact tiny-prime singleton, batch
smooth-part and grouped-charge challengers all complete the matched 30-digit
training cohort, but increase its cost versus current by approximately 93.5%,
15.7%, 17.6%, 21.4% and 9.8%, respectively. Their bounded prototypes stay in
the benchmark package. Batch recovery includes tree construction and scalar
exponent recovery, with an independent residual check. Grouped charges retain
identical cumulative work and reserve before execution; groups contain at most
four root ranges (16,384 score updates), with cooperative cancellation between
groups. This finite operation bound is not a real-time latency guarantee.
Neither smaller bookkeeping cost nor smooth-part throughput earns promotion.

Fresh medium and 30-digit collector comparisons validate 594 complete outcomes
across repeats. Atomic signatures, admitted rows, partial occupancy and complete
post-filter statistics match the frozen control for all complete arms,
including medium resieving. Fixed scoring reduces 30-digit visited division
primes about 43–46% while preserving the same 12/28/16 rows and filtered excess
zero on the three first-polynomial fixtures. This reduction in division work
does not imply a comparable full-factor gain or more useful dependencies.
Wide-block 30-digit resieving refuses setup under the original dense capacity
reservation (0/54 collector outcomes); its short refusal times are censored
costs, with no time-to-factor ratio. Capacity repair was left with the parallel repair owner. The integrated
bridge below retains that owner's subsequently accepted sparse support bound.

The whole-interval CRT probe covers the first polynomial of each seeded
training workload. No factor-base prime exceeds the complete interval in
any band. Eligible higher powers do occur: in 30 digits their 27–28 hits are
at most 0.072% of all power-mark hits. This scoped diagnostic does not justify
a family-wide half-sum implementation, which remains deferred pending a
workload demonstrating the cost. Eligibility uses the whole interval, never
the smaller block width.

Peak held-out owned workspace for current is approximately 20.46/23.45/36.82/
79.03 MiB; the plan arms add the reserved 1 MiB. Observed process/JIT RSS peaks
at 479.08 MiB and is separate from the owned-workspace allowance. Nine fresh
processes per arm give cold first-split medians of 114.544/114.706/115.347 ms
(control/current/joint) on one small training input and 802.951/795.655/800.459
ms on one 30-digit input. These include process startup and imports, and
establish no cold promotion. Separate instrumented 30-digit profiles reach
their time cap; they locate costs but provide no performance evidence.

The original study is pinned to isolated commit `171e69c`; its control,
corpora and frozen policy now live under `inputs/baselines/`, `inputs/corpora/`
and `inputs/controls/`, with unchanged bytes. Historical source/configuration
pins remain separate from `p38_r2_integration_frozen.json` and its newer
pre-R2 mainline baseline. The integrated runner uses `--integration` for a
combined-source bridge on these previously inspected fixtures; this is not
a fresh held-out promotion experiment. New and partial capture paths refuse
overwrites. From the committed integration tree:

```sh
mkdir -p v2/benchmarks/results/p38_r2
pypy3.11 -m v2.benchmarks.p38_r2 --integration --phase factor --variants control,current,fixed,plans,fixed_plans,resieve --output v2/benchmarks/results/p38_r2/factor-new.json
pypy3.11 -m v2.benchmarks.p38_r2 --integration --phase collector --split held_out --bands medium,30d --variants control,current,fixed,plans,fixed_plans,resieve --output v2/benchmarks/results/p38_r2/collector-new.json
pypy3.11 -m v2.benchmarks.p38_r2_eligibility --output v2/benchmarks/results/p38_r2/eligibility-new.json
```

The acceptance suite covers exhaustive signed windows, p=2, p|A, ramified
p|N', high valuations, capped lifts, tails, integer rounding, score saturation,
cap refusal, reserve-before-norm evaluation, cancellation/resume, polynomial
switches, checkpoint cache rebuilding and unchanged positional construction.

The isolated committed-files-only acceptance snapshot passes **282 PyPy
tests**, full lint and imports of **all 45 benchmark modules**, and verifies
the frozen runtime/driver hashes, immutable control and both certified corpora.
It uses the existing development tooling environment without local captures
or scratch inputs. Integration with the parallel repair branch remains a
separate step; no unrelated cleanup or repair edits are included in R2.

## R2 combined mainline acceptance — 4 October 2026

The accepted opt-in R2 delta is integrated with the committed repaired/R3
control `328c823`, preserving its resieve support bound and all other repairs.
The four original R2 input files remain byte-identical in the versioned input
folders. A separate integration baseline/freeze pins the combined source and
runner hashes; generated captures stay under ignored `results/p38_r2/`.

A serial quiet-window bridge validates all outcomes under the original shared
work/wall/CPU/storage limits. Each arm receives at least three seconds of
validated PyPy warmup and nine samples, extended to fifteen when needed.
The six arms complete 60/60 medium and 36/36 nominal 20-digit and 30-digit
training outcomes apiece. Small integration timing remains unstable and earns
no timing claim. An earlier capture overlapped another session's verification
and is diagnostic only; the repeated medium/20-digit/30-digit capture supplies
the following performance evidence.

| Inspected cohort | Pre-R2 repaired control | Integrated streamed powers | Integrated fixed + 1 MiB plans |
| --- | ---: | ---: | ---: |
| Medium, four attempts | 18.320 ms | 17.496 ms | 16.658 ms |
| Nominal 20-digit, four attempts | 288.380 ms | 262.263 ms | 248.609 ms |
| 30-digit, four attempts | 1.031281 s | 1.014611 s | 0.904477 s |

These are integration results on previously inspected certified inputs, not
fresh held-out selection evidence. The 30-digit joint reduction versus the
integrated streamed arm is 10.9% (conditional repeat interval 8.5–11.3%).
Retain the original defaults and frozen promotion decision; this bridge does
not select a new policy or imply a population-level speedup.

The medium and 30-digit collector bridges complete 54/54 and 90/90 outcomes
per final arm, respectively. Across all repeat attempts, 1,188 outcomes match
the repaired control's atomic signatures, row/partial counts and full
post-filter statistics, including resieving. The repaired 30-digit resieve
arm no longer refuses setup in this scope; peak collector-owned workspace
is 28.30 MiB under the same 128 MiB allowance. Its whole-factor cohort is
1.017612 s versus 1.014611 s for integrated bucket recovery, so no resieve
default promotion is claimed. B13 can remain deferred unless a calibrated
workload encounters a new capacity refusal; family-wide CRT remains C8's
conditional follow-up.

The combined tree passes 292 PyPy tests, full lint and all 48 benchmark-module
imports. Committed-files-only validation checks all required loaders and
immutable controls, without local captures or other sessions' readability
changes. The API and checkpoint contracts are retained, and only the R2
delta is included in the integration commit.


## P4.3 arithmetic backends — 5 October 2026

P4.3 adds an explicit `python-int` / `gmpy2-mpz` boundary across the factoring
engines, preprocessing/primality, relation recovery/extraction, smoothness
products and matrix bitsets. Keep Python integers as the provisional default.
GMP is an available optional track, but its individual-operation wins do not pass the
whole-engine promotion gate on this ARM64 PyPy build.

The measured environment is PyPy 7.3.23 implementing Python 3.11.15,
gmpy2 2.3.1 and GMP 6.3.0. The project-local PyPy environment supplies GMP;
the system PyPy remains a dependency-free track. The optional package pin
is in `v2/requirements-gmp.txt`. No CPython, GIL-release or thread-performance
comparison is claimed.

Reproduce the comparison from the repository root:

```sh
make -C v2 benchmark-backends PYTHON="$PWD/v2/.venv/bin/python"   REPETITIONS=9 WARMUP_SECONDS=3   BENCHMARK_OUTPUT=benchmarks/results/p43_backend.json
```

The runner loads the hash-checked pre-P4.3 control at `9b2d380` from
[the immutable baseline](inputs/baselines/p43_before_sources.json).
[The declared corpus](inputs/corpora/p43_backend_corpus.json) selects the
first previously published held-out input in each of eleven bands and retains
independent recursive Pocklington/trial proofs. Three seeds are 104729,
130363 and 155921. This is a bounded backend study on inspected inputs,
not fresh population-level calibration or proof of practical large balanced
completion. Algorithms receive inputs/configuration/seeds, never the oracle
factorizations.

Every arm has a separate PyPy interpreter/JIT. Each receives at least three
seconds of validated workload warmup and nine samples, extended to 18 or 27
when relative standard deviation exceeds 15%. Arm order alternates by case.
Each timed answer matches the same validated reference, including unresolved
cofactors, certainty and exact work ledgers. Setup, schedules and representation
conversions stay inside the measured calls. Matrix oracle validation is outside
timing. Cold import/startup is measured separately, including shutdown.
Some small samples remain variable at the 27-sample cap; intervals are
conditional repeat-timing uncertainty, not input-population intervals.

Representative medians from the final full capture:

| Workload | Python integers | GMP | Observation |
| --- | ---: | ---: | --- |
| 128 modular inverses, 256-bit operands | 1.051 ms | 0.261 ms | GMP 75.1% lower time; conditional interval 74.5–76.2% |
| 128 modular powers, exponent 65537 | 0.434 ms | 0.277 ms | GMP 36.2% lower time; interval 34.3–39.2% |
| 128 seventh roots of 256-bit inputs | 0.422 ms | 0.221 ms | GMP 47.7% lower median time; some micro samples remain variable |
| 128 GCDs, 256-bit operands | 0.151 ms | 0.231 ms | GMP slower with conversions included |
| 320-row, 256-column filtering/elimination | 5.664 ms | 48.995 ms | GMP approximately 8.7 times the time |
| Complete rho job, 265-bit certified composite | 0.557 ms | 2.069 ms | GMP approximately 3.7 times the time |
| Complete p−1 job, same composite | 0.220 ms | 0.407 ms | GMP approximately 1.9 times the time |
| Complete ECM job, same composite | 0.936 ms | 2.648 ms | GMP approximately 2.8 times the time |
| Small SIQS setup through extraction | 1.572 ms | 20.006 ms | Optional GMP loses on this small fixture |
| Eleven-input, three-seed portfolio | 48.844 ms | 98.531 ms | GMP 101.7% more time; interval 84.5–121.0% |

A final result-boundary confirmation, which canonicalizes only exact `mpz`
instead of coercing other caller-supplied container values, measures the same
portfolio at **47.982 ms int versus 84.628 ms GMP**: GMP uses 76.4% more time
(interval 64.7–93.4%). The frozen integer control is 45.522 ms; its difference
from current integers is inconclusive (interval crosses zero). Both final
captures support retaining the integer default. Do not interpret the small
QS/SSS fixtures as scaling evidence or promote a native micro-optimization
from this backend experiment.

All 33 portfolio answers per repetition match across arms: **15 complete,
18 explicitly unresolved** under the same 200,000-work-unit cap. Time and CPU
caps are disabled for this fixed-work comparison; therefore it does not
establish completion within a deadline. Engine-owned workspace caps match,
and each isolated arm records process peak RSS separately. No owned-workspace
or reconstruction failure occurs. The separate cold lifecycle medians are
87.381 ms int and 110.814 ms GMP; neither is a warmed execution measurement.

Generated evidence remains local under ignored `results/`:
`p43_backend_final_20261005.json` and
`p43_final_result_boundary_20261005.json`. Earlier in-process captures are
diagnostic: mixed int/mpz traces contaminated the integer QS timings. They
are excluded from the accepted performance evidence. The versioned runner,
corpus and immutable control permit reruns without those captures.

Acceptance includes the same exact arithmetic/result/certainty cases on each
available backend, nonunit inverses with retained GCDs, checked exact division,
large root boundaries, saturated stage replay, canonical JSON, backend/build
mismatch rejection, old integer checkpoint compatibility, streamed/external
polynomials, matrix identities and serial/spawned-worker consistency.
The implementation passes 315 PyPy tests and full lint. Committed-files-only
verification remains pending; this worktree is retained and uncommitted.
Generated captures are never required inputs.

### Size-dependent follow-up: protocol

The first study cannot settle algorithm/size selection: complete QS uses an
eight-digit fixture, and ECM uses B1/B2=200/2000. A new proof-backed balanced
corpus has independent screen and confirmation inputs at exactly 3, 10, 20,
30, 40, 50, 60, 70, 80, 90 and 100 decimal digits. It uses recursive
Pocklington/trial certificates and is not an RSA-distribution sample. The
versioned runner is `p43_sizes.py`; process-local representation experiments
live in `p43_experiments.py`, with no production global switch.

The arms are current native integers, persistent `mpz`, native loops with
GMP powering/inversion/roots above 64 bits, and (for QS) persistent large
`mpz` with native small roots/offsets/matrix masks. A frozen `before-int` arm
checks default-path regressions. Full stage studies use B1/B2=2000/147396
and 11000/1000000, two-curve/base trials, and 8192-evaluation rho attempts.
QS studies include complete smaller splits, streamed large coefficients,
larger factor bases/windows and longer fixed-work collection. Unresolved
large runs compare time to identical work, not successful factoring time.

Each arm has its own PyPy interpreter/JIT, validated warmup and at least nine
samples. Variable samples extend to 18/27. Outputs include identical divisors,
cofactors, stage-state digests, collection counts and logical work. Setup,
conversion and result construction are measured. Timing fields are excluded
only from output equality. Cold startup remains separate.

Two early expanded captures overlapped the p5.2 benchmark. Their JSON now
explicitly labels them diagnostic, and they are excluded from performance
decisions. Cross-chat exclusive timing windows are coordinated. The runner
checks foreign benchmark/test workers before, during and after each arm and
aborts its own worker if overlap recurs. It saves completed arms incrementally;
an interrupted worker supplies no timing samples. Rows explicitly distinguish
complete comparisons from unfinished sets of arms.

### First exclusive-window batch: 6 October 2026

The 02:17–02:31 UTC window completed 15 screen comparisons, with five seconds
of validated warmup and 9–27 samples per arm. Every divisor, unresolved
cofactor, canonical stage-state digest and work count matched. The higher
ECM/p−1 tier uses B1/B2=11000/1000000; rho uses two 8192-evaluation attempts.
Medians below include setup, shared schedule reuse within each campaign,
arithmetic conversions and result construction.

| Stage | Input digits | Native int | Persistent mpz | Selective GMP helpers |
| --- | ---: | ---: | ---: | ---: |
| ECM, higher tier | 20 | 14.05 ms | 95.91 ms | 14.84 ms |
| ECM, higher tier | 50 | 116.57 ms | 448.26 ms | 119.83 ms |
| ECM, higher tier | 100 | 195.73 ms | 488.02 ms | 202.20 ms |
| p−1, higher tier | 20 | 4.89 ms | 16.12 ms | 3.99 ms |
| p−1, higher tier | 50 | 69.19 ms | 231.63 ms | 65.42 ms |
| p−1, higher tier | 100 | 112.26 ms | 248.33 ms | 96.81 ms |
| rho | 20 | 1.82 ms | 14.82 ms | 1.88 ms |
| rho | 50 | 2.96 ms | 16.11 ms | 2.96 ms |
| rho | 100 | 5.51 ms | 15.82 ms | 5.47 ms |

The lower ECM tier also rejects persistent mpz: 6.35, 4.00 and 2.46 times
native time at 20, 50 and 100 digits. The 20-digit ECM/p−1 trials find proper
factors; the larger balanced inputs remain unresolved after their finite
campaigns. These are matched campaign costs, not estimated times to factor
50- or 100-digit balanced composites. Selective p−1 helpers improve all six
screen cases, with higher-tier median reductions of 18.3%, 5.4% and 13.8%.
Disjoint-input confirmation and frozen native regression controls remain
required before promoting a helper policy. The 64-bit cutoff is experimental.
Some persistent-mpz samples still vary after extension to 27; the large
losses do not establish a precise universal slowdown ratio.

Five-arm QS correctness probes at 500 million work units compare native,
persistent mpz, selective helpers, native-small/mpz-large and frozen native
code. The 20-digit input completes. The 50-digit input collects 2,428,946
positions and 66 relations; the 100-digit input collects 1,703,949 positions
and no relations. Both larger inputs stop at their finite work allowance.
All five arms match the complete logical output and remain within owned
workspace limits. With the factor-base bound of 100000, sampled polynomial
values are 99 bits for the 50-digit modulus (166 bits) and 182 bits for the
100-digit modulus (331 bits); small factor-base primes are at most 17 bits.
Input digits alone therefore do not specify the arithmetic widths in QS.
These probes include cold execution and are correctness evidence only.
An earlier long warmed QS comparison was interrupted at the agreed window
boundary before its complete set of arms and is excluded from the timing
conclusions.

These arms use gmpy2 inside the Python algorithms. The
[gmpy2 overview](https://gmpy2.readthedocs.io/en/latest/overview.html) describes
a variable integer-performance crossover; it does not establish a threshold
for this PyPy implementation. The
[PyPy C-extension FAQ](https://doc.pypy.org/faq.html#do-c-extension-modules-work-with-pypy)
documents compatibility-layer/refcount costs. Binding and boxing overhead
are a plausible explanation for cheap repeated mpz operators losing while
larger single-call helpers win; that attribution remains an inference rather
than a measured profile. Native GMP-ECM/C QS programs and fused/CFFI kernels
are separate implementations and are not ranked by this substitution study.

Reproduce the completed stage screen with:

```sh
v2/.venv/bin/python -u -m v2.benchmarks.p43_sizes \
  --suites ecm pm1 rho --digits 20 50 100 \
  --warmup-seconds 5 --repetitions 9 \
  --output v2/benchmarks/results/p43_quiet_stages_screen.json
```

The proof-backed QS probe uses `--probe --suites qs --digits 20 50 100
--large-qs --qs-work 500000000 --backends python-int gmpy2-mpz helpers-gmp
gmp-small-native before-int`. Captures are ignored local evidence; the
versioned corpus and runners are sufficient to rerun the comparisons.

### QS size trend: 100-million-work screen

The next exclusive window completed four five-arm QS comparisons. Each arm
receives five seconds of validated warmup and 9 samples, extended to 18 for
the variable 40-digit native arm. All outputs, collection statistics and work
counts match the frozen native control. These runs stop at the work cap;
they measure the cost of identical bounded collection, not time to factor.

| Digits | Native int | Persistent mpz | Native loops, GMP helpers | Mpz with native small values | Frozen native | Mpz/native |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 40 | 0.685 s | 6.759 s | 0.654 s | 4.844 s | 0.638 s | 9.87x |
| 50 | 0.566 s | 6.686 s | 0.546 s | 4.104 s | 0.557 s | 11.81x |
| 80 | 0.414 s | 5.229 s | 0.421 s | 3.199 s | 0.416 s | 12.64x |
| 100 | 0.355 s | 4.428 s | 0.351 s | 2.750 s | 0.360 s | 12.46x |

The 50/80/100-digit cases share a factor-base bound of 100000, width 65536 and the
same 100-million-work allowance. Persistent mpz's penalty is roughly flat
across those sizes, unlike the shrinking relative gap in ECM/rho. Keeping
small roots, collector offsets and matrix masks native reduces the penalty
to 7.25/7.73/7.74x but does not reverse it. The 40-digit case uses a 50000
bound and width 32768, so it is not a controlled continuation of that trend.

Actual sampled F values are 81/99/149/182 bits at 40/50/80/100 digits; the
moduli are 132/166/266/331 bits. Collection visits 524296/454659/375958/307202
positions and finds 142/12/0/0 relations. The decrease in absolute runtime
does not mean larger integers factor faster: the fixed work cap permits
fewer positions and changes the mix of setup and collection costs.

Frozen/current native differences at 50 and 80 digits have conditional
repeat intervals crossing zero. The 40-digit native arm remains variable
(14.7% relative standard deviation), and its apparent 6.8% increase needs
a dedicated repeat before it can establish a default-path regression.
The helper arm supplies no consistent size trend or selected QS policy.

Reproduce this screen with:

```sh
v2/.venv/bin/python -u -m v2.benchmarks.p43_sizes \
  --suites qs --digits 40 50 80 100 --large-qs --qs-work 100000000 \
  --backends python-int gmpy2-mpz helpers-gmp gmp-small-native before-int \
  --warmup-seconds 5 --repetitions 9 --window-seconds 3600 \
  --output v2/benchmarks/results/p43_quiet_qs_trends_screen.json
```

### QS longer allowance: completed 50-digit comparison

The fivefold longer screen completes all five 50-digit arms with nine validated
samples each. Every arm performs 499999846 work units, visits 2428946 positions
in 612 blocks across 19 polynomials, and retains 66 relations. All stop at the
work allowance with the same unresolved cofactor and 70.51 MiB owned workspace.

| Arm | Median time | Relative to current native |
| --- | ---: | ---: |
| Current native | 2.837 s | 1.00x |
| Persistent mpz | 35.241 s | 12.42x |
| Native loops, selective GMP helpers | 2.962 s | 1.04x |
| Mpz with native small roots/offsets/masks | 21.359 s | 7.53x |
| Frozen native | 2.787 s | 0.98x |

Persistent mpz's conditional repeat interval is 12.01–12.81x native time.
Relative sample deviations are 4.6% native, 0.9% mpz, 3.5% helper, 0.5%
native-small/mpz-large and 1.2% frozen native. Longer collection does not
amortize away the GMP penalty: it changes from 11.81x at 100 million work
units to 12.42x at 500 million. The difference between those ratios is not
a separately randomized duration effect, but neither run supports a crossover.
Helper/current and frozen/current repeat intervals include zero difference.

The user ended the 100-digit extension and all further experiments. Its parent
and child were stopped and their absence verified. The incremental capture
retains the completed 50-digit comparison and explicitly labels the unfinished
100-digit set; the latter supplies no complete backend comparison. No new
30-digit or disjoint-input confirmation was performed after this stop request.

Reproduce the longer workload with the screen command above, replacing the
digits with `50`, work with `500000000`, and output with
`v2/benchmarks/results/p43_quiet_long_qs_screen.json`.

The available size evidence supports retaining native integer loops in this
PyPy implementation. The QS penalty is roughly flat across the common
50/80/100-digit short-run configuration, and the longer 50-digit run confirms
the large loss. These measurements do not select an automatic backend policy
or close the disjoint confirmation gate. The final implementation passes
316 PyPy tests and full lint, including canonical public divisors reentering
their selected backend. The worktree is retained and uncommitted; this task
does not merge it into mainline.
