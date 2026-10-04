# Benchmarks

Run from the repository root on PyPy implementing Python 3.11:

```sh
make -C v2 benchmark WARMUP_SECONDS=3 REPETITIONS=9
make -C v2 benchmark-phase-two WARMUP_SECONDS=3 REPETITIONS=9
make -C v2 benchmark-phase-three-reference WARMUP_SECONDS=3 REPETITIONS=9
make -C v2 benchmark-phase-three-collector WARMUP_SECONDS=3 REPETITIONS=9
make -C v2 benchmark-phase-three-pipeline WARMUP_SECONDS=3 REPETITIONS=9
make -C v2 benchmark-phase-three-siqs WARMUP_SECONDS=3 REPETITIONS=9
```

Use unique output names when supplying BENCHMARK_OUTPUT. Outputs remain
local and Git-ignored: raw JSON/gzip samples, stdout transcripts, profiles,
generated freeze files and validation/verification dumps. They are useful
working evidence, rather than source files to publish for every run.

Committed inputs include the independent Phase 2 corpora, competitor metadata,
exact old-source baselines in audit/, and the frozen P3.3 configuration used
by `phase_three_ownership.py`. Those snapshots are hash-checked before use.
The preserved v1 comparison is emulated Python 2 through `lib2to3`, not a
native Python 2 measurement.

## Recorded P3.3 result

A matched warmed comparison on 16 small held-out balanced inputs reduced
complete QS cohort median time from **27.220 to 21.446 ms**, about **21.2%**.
Both arms used the same frozen bucket recovery and filtering settings,
including setup through verified factor reconstruction.

The measurement used PyPy 7.3.23 / Python 3.11.15 on macOS arm64, default JIT,
at least three seconds of validated warmup, and nine samples. Training seed
329 preceded a configuration freeze; the final held-out seed was 335.
The fixed-cohort bootstrap median-ratio interval was 0.726–0.863.

Reproduce the final comparison with:

```sh
pypy3 -m v2.benchmarks.phase_three_ownership \
  --warmup-seconds 3 --repetitions 9 \
  --output v2/benchmarks/qs_comparison_LOCAL.json
```

These are 23–26-bit inputs, not evidence of large SIQS scalability. Cold
startup has no demonstrated gain. Tiny collector windows still trail
exhaustive enumeration. Resieving and tighter candidate scoring did not
improve complete-run time and remain optional.

Archive full raw results locally when needed. Publish selected summaries
with environment, commands, correctness checks and scope; avoid committing
a transcript or profiler dump for every development attempt.

## P3.1–P3.3 correction audit (M28)

157 PyPy tests and lint pass. Additional independent checks compare 214 root
sets, 648 collector windows and 320 lifted kernels. Fix storage-cap extraction
and combined live memory reservations. On a fresh 16-input 23–26-bit cohort,
matched M26/current warm medians are 12.773/13.598 ms with root recovery and
9.723/10.050 ms with bucket recovery: costs of 6.5% and 3.4%, respectively.
Nine samples follow at least three seconds of validated warmup per arm;
all outputs reconstruct and both factors are proven. The matrix utility is
0.113/0.117 ms with unchanged kernel code; no causal algorithm regression is
claimed. Cold lifecycle measurements are inconclusive. This is a correctness
repair, with its measured cost retained; no speedup or larger-size claim.

```sh
PYTHONDONTWRITEBYTECODE=1 pypy3 -m v2.benchmarks.phase_three_audit \
  --warmup-seconds 3 --repetitions 9 \
  --output v2/benchmarks/qs_audit_LOCAL.json
```

The hash-checked `qs_m26_baseline.json` is a required immutable source input.
It contains exactly the pre-fix QS module bytes, extracted from the local
preflight snapshot. Generated captures stay local. The recorded cold harness
used the larger preflight snapshot before this reader relocation; loading
that local diagnostic data is not a production-factorization startup cost.
P3.4 family/root reuse and P3.8 matrix experiments retain their own gates.


## P3.4 family/root foundation (M29)

167 PyPy tests and lint pass. Sixteen independent CRT families produce 60
polynomials whose roots agree with full recomputation and exhaustive modular
enumeration, including multipliers 1/2/3/9, recentering and primes dividing A.
Family refusal/resume, corrupt checkpoints and cached collector admission are
checked. P3.1–P3.3's 214 root, 648 collector and 320 kernel controls still pass.

A root utility including base/setup and all eight polynomials gives warm
medians **1.287 ms full / 0.818 ms cached**, 36.4% lower for the cache. This
utility uses n=4001*5003, bound 2000 and four A factors. It does not measure a
complete factorization speedup.

On 16 fresh 24–26-bit inputs at seed 2936, both arms complete **13/16**;
three keep their cofactors at finite family exhaustion. Whole-cohort attempt
medians, including unfinished attempts, are **48.693 / 48.621 ms**: unchanged
within sample spread. Both use bound 200, h=1, three A factors, eight families
(assignment seed 2934), width 513, bucket recovery, residual 1 and weight-two
filtering. Each input shares 200 million work units, ten wall/CPU seconds and
32 MiB owned workspace across setup/extraction. Stores are fresh per
polynomial. All completed splits prove both factors and reconstruct n;
unfinished records also reconstruct n.

Each arm has three seconds of validated warmup and nine samples of ten cohort
calls. Nine separate cold child lifecycles give median 304.6/303.7 ms,
CPU 301.7/300.3 ms and RSS 72.8/71.9 MiB. Cold startup shows no established
gain. The runtime is PyPy 7.3.23 / Python 3.11.15, macOS arm64, default JIT
and serial execution. This one-seed probe is experimental and has no dispatch
or large-band promotion claim.

```sh
make -C v2 benchmark-phase-three-families WARMUP_SECONDS=3 REPETITIONS=9 \
  BENCHMARK_OUTPUT=benchmarks/qs_families_LOCAL.json
```

Next P3.4 gates are shared relations, full job checkpoints, multiplier and
parameter experiments, bounded ECM-to-SIQS dispatch and larger comparisons.
The family-only checkpoint is not a resumable full relation job.


## Bounded P3.4 implementation and diagnostic probes (M30, 4 October 2026)

185 PyPy tests and lint pass. Full SIQS shares checked relations across
polynomials, preserves complete resumable job/store identity and rebuilds
charged root/matrix/elimination/extraction state from compact checkpoints.
Finite width/yield/trivial-dependency controls preserve cofactors. The
portfolio reaches optional SIQS after ECM under one allowance and retains
recursive reconstruction. Earlier portfolio checkpoints without SIQS remain
readable. The broad independent audit still passes 214 root sets, 648 collector
windows and 320 lifted kernels.

Freeze 70 independently certified, disjoint inputs before tuning: 18 training
and 52 held-out, including 32 held-out 24–26-bit inputs and 20 balanced
30/40/50/60/70/80-digit inputs. There are 830 recursive Pocklington/trial proof
nodes. Main large bands include residues 1/3/5/7 modulo 8; capped 70/80-digit
exploration includes 1/5. Algorithms receive n and seeded configuration;
known factors/certificates are used only for validation. Runtime probable
prime labels are preserved separately from the independent proofs.

Training controls vary base size, A count/pool, width/block, thresholds,
residuals, scored h, diversity and recovery. Each trained/measured arm uses
at least three seconds of validated PyPy warmup and nine samples; unstable
runs extend warmup and sample count. Freeze before opening held-out outcomes.
Two held-out seeds are 7 and 29. Common limits are 200 million work units and
64 MiB owned workspace per input, with ten wall/CPU seconds for small inputs
and 0.2 seconds for large inputs (large training uses 0.08 seconds).

Small frozen settings: bound 200, half-width 256, one A factor, eight families,
pool 16, h=1, bucket recovery, residual 1, block/batch 256, zero extra threshold,
weight-two filtering and fixed width. The shared/fresh-store control uses the
same settings with three A factors. Cohort times below average the two seed
medians; every time covers all outcomes, including unfinished attempts.

| Small held-out arm | Complete input/seed cases | 32-input cohort ms |
| --- | ---: | ---: |
| Fresh stores, three A factors | 41/64 | 194.1 |
| Shared stores, three A factors | 64/64 | 44.0 |
| Frozen tuned SIQS | 64/64 | 19.7 |
| QS | 64/64 | 18.6 |
| MPQS | 64/64 | 20.2 |
| Bounded ECM portfolio | 64/64 | 4.0 |

Shared stores add 35.9 percentage points of completion; the paired fixed-corpus
bootstrap interval is +25.0 to +48.4 points. Their lower attempt cost includes
unfinished fresh-store runs and is not a successful-factorization speed ratio.
Tuned/shared complete-time ratio has a paired interval 0.363–0.549. Tuned SIQS
versus QS is inconclusive (0.981–1.134); no winner is declared. Scored h costs
about 22.0 ms versus 19.4 ms for h=1 with identical completion, so h=1 remains
the control. Simple families complete 62/64 versus 64/64 diverse; keep the
finite diverse schedule without claiming a universal optimum.

One recovery arm remained noisy after its first extension. A separate matched
repeat uses five seconds of validated warmup and 31 samples per arm: fixed
17.124 ms, recovery 17.290 ms, both stable. The paired ratio interval
0.994–1.005 shows no meaningful recovery gain. Fixed width remains the default;
recovery is available when explicitly configured.

Large frozen settings use bound 2000, half-width 512 and four A factors,
with the same finite family/pool/bucket controls. Distinct completion counts
under the declared caps are:

| Balanced digits | QS | MPQS | SIQS | Bounded ECM |
| --- | ---: | ---: | ---: | ---: |
| 30 | 0/8 | 0/8 | 0/8 | 2/8 |
| 40 | 0/8 | 0/8 | 0/8 | 0/8 |
| 50 | 0/8 | 0/8 | 0/8 | 0/8 |
| 60 | 0/8 | 0/8 | 0/8 | 0/8 |
| 70, capped exploration | 0/4 | 0/4 | 0/4 | 0/4 |
| 80, capped exploration | 0/4 | 0/4 | 0/4 | 0/4 |

A matched fresh-store large control also has zero completions, so the shared
store causes no completion regression in these declared classes. Scored h
raises total large relation yield from 6 to 17 in a separate capped probe;
recovery gives 9, simple families 2. None complete, and yield alone does not
promote a default. The large cohort attempt medians are 93.2/687.8/1114.9/432.2
ms for QS/MPQS/SIQS/ECM. These include natural local exhaustion and capped
failures; their times cannot rank successful large factorizations.

Checkpoint roundtrip has three-second validated warmup and nine samples,
median 1.94 ms for the declared small representative, with a 6283-byte maximum
snapshot estimate. Root/store/prefix reconstruction charges additional work.
Fault tests separately retain pending matrix/extraction progress and reject
rehashed invalid arithmetic. Maximum measured owned workspace is about
10.9 MiB; this is distinct from process RSS. Base growth and disk spill are
disabled; width recovery retains the exact base and has no remapping cost.

Nine cold child lifecycles per arm cover one small representative and one
input in each larger band. Median lifecycle seconds are QS 0.291, MPQS 0.537,
SIQS 0.691 and ECM 0.411, with respective peak RSS about 122.4/128.6/132.2/131.8
MiB. This diagnostic harness also loads its frozen training capture; the
figures are not isolated production startup costs. Warm RSS is shared across
arms. Separate instrumented profiles expose marking, candidate/residual
recovery, verification, filtering and extraction; profile times are not
performance evidence.

Historical M30 decision, superseded by M31 below: the SIQS implementation
is tested; P3.4's large-number performance gate remained open. The earlier complete-gate claim is withdrawn.
The 0.2-second runs are diagnostic probes and do not assess practical
50–60-digit factoring or adequate capped 70–80-digit exploration. The small
successful cohort is 24–26 bits, roughly 7–8 decimal digits. Longer per-input
wall/CPU budgets must be paired with trained finite base, A, sieve interval,
family and relation limits; the current eight-family/512-half-width control
can exhaust its local schedule before a longer timer. Separate each stop
reason and preserve unfinished outputs in the next held-out evaluation.
SIQS remains an optional fallback; ECM stays the default, with no measured
SIQS/ECM crossover, imported digit cutoff or universal parameter optimum.

```sh
make -C v2 test
make -C v2 lint
make -C v2 benchmark-phase-three-siqs WARMUP_SECONDS=3 REPETITIONS=9 \
  BENCHMARK_OUTPUT=benchmarks/siqs_comparison_LOCAL.json
```

The retained `phase_three_p34_corpus.json` contains independent proof inputs.
The runner creates a unique freeze/capture; raw samples, extended controls,
profiles, journal and detailed acceptance summaries stay local. Large local
captures are losslessly archived with original hashes. The local control
provenance record corrects a copied config template in the extra fresh-store
capture; its executed override was shared_relations=False and all original
samples remain unchanged.


## P3.5 bounded SSS/SSSf challenger (4 October 2026)

The independent adapter uses exact CRT collisions, capped product/remainder
trees and complete prime-power recovery. Every admitted original polynomial
value passes the existing verifier; partial matching, provenance, filtering,
GF(2) solving and extraction use the common pipeline. Ten new tests cover an
independent scalar smoothness oracle, brute CRT/collision enumeration, signed
positions, powers above the upstream shortcut, storage refusal, candidate
prefix interruption and retained-budget resume. The shared worktree passes
202 PyPy tests and lint at that capture. The later explicit opt-in CLI/API
integration retains the same engine; automatic dispatch stays unchanged.

The retained certified `phase_three_p34_corpus.json` supplies 12 small training
inputs and 32 small held-out inputs, plus separate training and held-out
balanced digit bands. Its SHA-256 is
`47b4c57eeddcbe3a8a0a7b2bbd5c32877504a81d5e08af323c7d23f2a4d75972`.
Seeds are 7 and 29. Algorithms see n/configuration, with factor certificates
used separately for validation. Completed and unfinished outputs reconstruct
n; runtime probable-prime labels are preserved even when the corpus proves
the factors independently.

Shared per-input limits are two billion work units, five wall/CPU seconds,
64 MiB owned workspace, one core and an observed 512 MiB process-RSS ceiling.
Timing includes setup through extraction and terminal-factor classification.
Both arms use residual bound 10000, atom/full/partial caps 8192/4096/1024,
h=1 and weight-two filtering. RSS is cumulative within a process and is not a
per-job workspace estimate. Warm evidence has at least three seconds of
validated warmup and nine samples, with extensions when noisy. Cold startup
is not measured here and no cold-performance claim is made. Source hashes
identify a frozen local QS snapshot so concurrent P3.4 edits cannot change an
arm mid-run; captures report no source changes during execution.

SSS training tests bounds 400/1000 and selection sizes 3/6 using only the
small training inputs. It freezes bound 400, selection size 6 and 128 search
assignments. SIQS independently tests bounds 200/400 with one/three A factors,
half-width 256 and 64 families/pool entries, with power scoring and bucket
recovery; it freezes bound 200 and one A factor. All training measurements
are validated/warmed/repeated. SSSf is a predeclared filtered challenger,
not a separately optimized arm: seven draws, first-half factor base and
cutoff `10**max(1, decimal_digits//2-1)`. The cutoff deliberately loses yield.

| Small held-out arm | Completed input/seed cases | Warm 64-case cohort median |
| --- | ---: | ---: |
| Training-selected SIQS | 64/64 | 184.1 ms |
| Training-selected SSS | 64/64 | 509.2 ms |
| Predeclared SSSf | 63/64 | 1089.2 ms |

This cohort is 24–26 bits, about 7–8 decimal digits. SSS is slower here; SSSf
also retains one cofactor at assignment exhaustion. Times include every
attempt, so SSSf's figure is not a successful-factorization speed ratio.
An earlier roots/three-A-factor SIQS control took 351.6 ms; retain that local
capture as historical evidence rather than the final comparison.

At 30 digits, predeclared bound 5000 and 4096 assignments give SSS and SSSf
**8/8 complete cases in all nine samples**, with respective cohort medians
**7.066 / 11.175 seconds**. The original three-A-factor SIQS control completes
0/8 and stops on work/stalled-yield limits. Six training-only SIQS probes
vary A factors 3/4/5 and half-width 512/4096, using power scoring, bucket
recovery and 64 families/pool entries. None complete; highest checked-row
yield freezes four A factors and width 512. That stronger control also
completes **0/8 in all nine samples**: six wall limits and two stalled-yield
stops in the first cohort. Its 33.538-second cohort median includes failures
and cannot establish a successful-factorization speed ratio. Source hashes
for both 30-digit experiments agree on all arithmetic/QS modules.

SSS's successful 30-digit runs include verified relations, nontrivial kernels
and proper extraction, not just relation-yield counts. They support this
bounded challenger, but do not establish an advantage over feasible, fully
tuned SIQS. No held-out outcome is used to retune an accepted configuration.
P3.4's later feasible-control study and fresh confirmation inputs retain
their separate gates.

Larger runs are **one-cohort diagnostics**, with four held-out inputs and two
seeds in each band, using the same limits. Training-only SIQS probes freeze
one control per band; predeclared SSS bounds are 8000/22000/82000 for
40/50/60 digits. Each arm completes 0/8 in every band. At 40 digits all stop
on wall time; at 50 digits SIQS stalls while both SSS arms hit wall limits.
At 60 digits SIQS hits wall limits; each SSS arm has two memory refusals and
six wall stops. These are not practical 50–60-digit performance evidence.
The fixed residual bound is especially restrictive once it falls below the
largest base prime. The 60-digit memory refusal also exposes the tree/storage
capacity limit. Checked candidate/row yield without usable dependencies is
reported in raw diagnostics and does not count as a factoring success.

The separate unchanged-upstream arm hash-checks `sss.py`, `sssf.py` and
`mstep.py` at commit `8dbaf6d39ab88a40380965d25ec2c363d7f27358` of
[SmoothSubsumSearch](https://github.com/sbaresearch/smoothsubsumsearch/tree/8dbaf6d39ab88a40380965d25ec2c363d7f27358).
PyPy 7.3.23 / Python 3.11.15 on macOS arm64 runs SymPy 1.14.0, mpmath 1.3.0,
gmpy2 2.3.1 / GMP 6.3.0. Both functions complete 8/8 cases in all nine samples
on the first four small held-out inputs; warmed eight-case medians are about
45.7 / 63.6 ms. SSSf arguments are explicitly `(prop=2, digred=0)`; its
demonstration's `(10, 5)` is unsuitable for these small inputs. Global seeds,
floating-point setup, raw-prime-count parameter table, finite prime-power
shortcut and matrix helper using `gmpy2.isqrt` remain unchanged. Progress is
consumed by a bounded output validator. A child watchdog imposes 180 wall/CPU seconds
and 512 MiB RSS on the whole reproduction. This is a labelled compatibility
arm, with distinct backend/settings/cohort and no adapter speed ratio.
The upstream harness's NumPy statistics import is not required.

Adapter adaptations are local per-assignment seeds, sorted distinct indices,
distinct-prime collision counts, exact setup, filtered-base cardinality,
complete smooth-part powers, explicit finite capacities and common checked
postprocessing. They are independently implemented and are not described as
an unchanged reproduction. The adapter imports neither SymPy nor gmpy2.

Run from the repository root, retaining generated outputs locally:

```sh
pypy3 -m v2.benchmarks.phase_three_sss --train-only \
  --output v2/benchmarks/sss_training_LOCAL.json
pypy3 -m v2.benchmarks.phase_three_sss --train-siqs-only --band small \
  --output v2/benchmarks/siqs_small_training_LOCAL.json
pypy3 -m v2.benchmarks.phase_three_sss --band small \
  --training v2/benchmarks/sss_training_LOCAL.json \
  --siqs-training v2/benchmarks/siqs_small_training_LOCAL.json \
  --output v2/benchmarks/sss_small_LOCAL.json
pypy3 -m v2.benchmarks.phase_three_sss --train-siqs-only --band balanced_30d \
  --output v2/benchmarks/siqs_30d_training_LOCAL.json
pypy3 -m v2.benchmarks.phase_three_sss --band balanced_30d \
  --training v2/benchmarks/sss_training_LOCAL.json \
  --siqs-training v2/benchmarks/siqs_30d_training_LOCAL.json \
  --output v2/benchmarks/sss_30d_LOCAL.json
```

For larger diagnostic bands, train that band, then pass `--diagnostic-only`
with its `--siqs-training`. Use a new output for every run. The Make target
`benchmark-phase-three-sss` reproduces the initial default control; the
explicit training commands above reproduce the strengthened comparison.
To reproduce upstream, retrieve the three pinned files into a local source
directory, then use the existing PyPy development environment:

```sh
v2/.venv/bin/python -m pip install sympy==1.14.0 mpmath==1.3.0 gmpy2==2.3.1
v2/.venv/bin/python -m v2.benchmarks.phase_three_sss_upstream \
  --source-dir v2/audit/p35_upstream --mode sss \
  --output v2/benchmarks/upstream_sss_LOCAL.json
```

Repeat with `--mode sssf` and another output. The runner refuses changed
source bytes and requires the watchdog's own-child RSS monitor to work.
Generated captures, source freeze, upstream retrievals, transcripts and
detailed acceptance journal stay local. Required proof corpus and runners
remain repository inputs. Decision: retain the bounded challenger; defer
dispatch/promotion until the fresh feasible SIQS comparison and P3.6.1 gates
pass. P3.8-R5 reuses that evidence for broader reconciliation.

### Explicit opt-in dispatch and checkpoint acceptance

`--method sss` / `sssf` now selects bounded SSS after normal exact
preprocessing, with full recursive classification and reconstruction.
`PortfolioConfig(sss=...)` instead adds a fallback after the configured
rho/p−1/ECM schedule. Automatic defaults do not select SSS. Full SSS
checkpoints preserve checked relations, an unfinished assignment cursor and
solver/extraction progress, with charged reconstruction and bounded encoding
workspace. Portfolio version 4 retains compatible version-2/3 resume.
Six new tests cover both CLI methods, signs/multiplicity, resource refusal,
post-ECM dispatch, serialized resume during collection/elimination/extraction,
rehashed arithmetic/progress corruption and checkpoint capacity.

A matched integration probe uses eight of the already inspected small
held-out inputs and two seeds, 16 cases per cohort. Direct SSS plus terminal
classification and selected portfolio SSS both complete every case. Child
seeds, bound 400, 256 assignments, two-billion-work/five-wall-and-CPU-second
limits, 64 MiB SSS workspace within 80 MiB total, and one core match.
All outputs agree with the independent factor certificates. Both arms have
three seconds of validated warmup and nine stable samples. Cohort medians
are **69.79 / 69.40 ms**; the portfolio/direct timing ratio interval is
**0.936–1.081**, so no timing difference is established. This is an integration
overhead probe on inspected tiny inputs, not fresh engine-promotion evidence.
Cold startup remains unmeasured. The local dispatch source freeze, raw
comparison and full validation remain under the ignored audit/capture paths.


## Larger P3.4 evaluation and resumable runs (M31)

The independent `phase_three_p34_large_v3_corpus.json` freezes 80 unique
inputs and 4,875 recursive Pocklington/trial proof nodes: 18 training inputs
and 62 held out. Balanced semiprimes cover 30/40/50/60/70/80 decimal digits
and residues 1/3/5/7 modulo 8. They have distinct factors of comparable size
and exclude exceptionally close pairs. Additional cases include uneven
60/80-digit composites with exact 5/10/20-digit smaller factors, close pairs,
squares, primes, three-prime composites and small smooth p−1/p+1 factors.
Factors/proofs validate results and are never passed to a factoring method.
Pocklington construction is a curated certified corpus, not uniform sampling
from all primes; the close-pair construction is explicitly structured.

The `phase_three_large` runner trains finite base/A/width/residual/multiplier
choices and freezes configurations before held-out evaluation. It uses
10^13 work units and 1 GiB owned workspace per call (768 MiB reserved for
SIQS), with matched wall/CPU ceilings of 30 seconds at 30/40 digits,
60 seconds at 50/60 digits and 300 seconds for 70/80-digit exploration.
The optional prime-power sieve and enlarged finite relation/checkpoint caps
are explicit configurations; the library defaults are unchanged.

QS/MPQS/SIQS/ECM comparisons at 30–60 digits use four residue cases at seeds
7/29 plus one case at seed 47, after at least three seconds of validated
warmup. These are nine case/seed outcomes, not timing repetitions of one
input. The capped 70/80-digit exploration uses one fixed representative at
seed 7 for SIQS and ECM, with no timing promotion claim. Separate 30-digit
power-sieve controls use nine repeated calls and extend unstable runs.
Cold child lifecycles and instrumented profiles are separate measurements.

Unfinished primary balanced, varied and continuation calls save their full
checkpoints to sidecars with byte count, digest and write cost. Original
control/cold calls retain results and statistics without writing sidecars.
A fixed 50-digit SIQS representative can
continue to cumulative 600 and then 1,800 active wall/CPU seconds, stopping on
completion or an explicit non-time limit. Consumed work/time and rebuild
costs survive each resume; paused time is excluded. The continuation is a
separate time-to-factor experiment and cannot be mixed into fixed-limit
completion rates. A timeout remains unfinished; relation yield is not a
successful factorization. Memory, store, schedule and timer stops are recorded
separately, and all factors plus unresolved cofactors reconstruct the input.

```sh
make -C v2 benchmark-phase-three-large LARGE_PHASE=all \
  BENCHMARK_OUTPUT=benchmarks/large_siqs_UNIQUE_pypy.json
```

Individual phases are `train`, `balanced`, `varied`, `continuation`,
`controls` and `cold`; they preserve prior captures and require the same
frozen source and corpus. `--training-capture` permits inherited diagnostic
choices only with a retained capture hash and two fresh training calibration
calls per band; enlarged memory/checkpoint allowances are disclosed. The
main experiment freezes one source snapshot to exclude concurrent edits.
Raw journals, sidecar checkpoints, profiles and exact measured-source archives
remain local. The completed results and bounded acceptance decision follow.

### Completed larger evaluation and filtering repair (M31, 4 October 2026)

The initial frozen comparison validates 362 measured outputs and 145 saved
primary checkpoint sidecars. It includes 148 balanced, 152 varied, 24 control
and 36 cold attempts, plus two cumulative 50-digit continuation attempts.
These captures precede the incidence-filter repair and retain their original
source hashes. The later repaired SIQS/ECM confirmation uses the same frozen
settings, input/seed pairs and caps on already inspected cases; it is not a
fresh held-out promotion study. Successful and unfinished times remain
separate throughout.

On the original source, QS/MPQS complete none of the nine outcomes in each
30/40/50/60-digit band. They stop at their finite reference schedules; these
limited QS windows and MPQS q²-A schedules do not support a general method
ranking. SIQS completes 9/9 at 30 digits (descriptive median 5.947 s) and none
at the larger fixed limits. ECM completes 9/9 at 30 and 40 digits, 4/9 at
50 digits and none at 60/70/80 digits. The successful 50-digit ECM subset
ranges from 5.772 to 13.307 s; its median is not representative of failures
or a matched complete cohort. Stage, work, resource and terminal-event details
remain in the raw captures.

The varied study runs 38 independent held-out inputs at two seeds in both
arms. Preprocessing plus ECM and the same preprocessing plus optional SIQS
fallback each complete **63/76**. Both finish every case with a 5- or
10-digit smaller factor, close pairs, squares, primes and small smooth
p−1/p+1 factors. Each finishes 3/8 three-prime cases and 0/8 composites with
a 20-digit smaller factor. Fallback adds no completion here and its failed
calls consume the 60/300-second limit; the ECM-only arm exhausts its declared
six-curve schedule sooner. This varied schedule is weaker than the balanced
ECM schedule. Tiny preprocessing times are not evidence of an engine speedup.

The old-source prime-power score control completes 9/9 repeated 30-digit
calls in each arm: **5.808 s** with powers versus **20.080 s** with the
conservative upper-bound score. Both are stable after validated warmup.
Simple A schedules at 50/60/80 digits exhaust their finite families without
a factor; optional width recovery also adds no completion at the matched
60/60/300-second ceilings and is not promoted. Nine independent cold-child
lifecycles per arm include import/startup: median **7.188 s** for completed
SIQS and **0.803 s** for completed ECM at the 30-digit representative.
The separately instrumented 80-digit profile is diagnostic, not a timing arm.

A 50-digit continuation exposed repeated global matrix-incidence rescans:
the original attempt spent about 1,443 of its 1,802 active seconds in
filtering. Filtering now maintains column incidence bitsets and updates only
affected columns, with finite workspace, cooperative cancellation and the
same deterministic output/provenance. All 2,010 old/new oracle comparisons
match rows, masks, dependencies, statistics and reservations exactly. The
matched filtering-only control changes the filter function while retaining
all other modules/settings: **5.920 s → 3.425 s** median over nine stable
complete calls, a **42.1% elapsed reduction on this inspected 30-digit
representative**. The owned prior filter has a 4,096-row cap; this control
stays below it. This is not a whole-corpus speed ratio.

The retained 50-digit state then completes after **299.306 additional
seconds**, preserving consumed work, verified rows and duplicate identity.
Its cumulative active wall time is **2,100.709 s**; 299 seconds is a
continuation time. A separate empty-store run on the repaired source completes
in **1,183.355 s (19 min 43 s)** after 3.061 s of validated training warmup:

```text
62735374527037454574430698258919844832703035589529
= 6797018013217317517573283 × 9229837909071854200546963
```

It uses 4,497 verified relations, 548 distinct polynomials, zero duplicates,
3,766,141,943,212 work units, 649,573,392 bytes of reported owned peak
workspace and 312,770,560 bytes of process lifetime peak RSS. Engine stage
timers include 533.951 s collection, 261.600 s relation preparation,
339.775 s filtering, 33.985 s elimination and 0.030 s extraction. Stage timers
are engine-inclusive; outer classification/setup and snapshot costs are
included in elapsed time. Both returned factors retain runtime
`probable_prime` certainty; independent corpus certificates prove the inputs.
This is one inspected, explicitly configured success, not a general median
or the default CLI's timing.

Other large-state repairs enlarge explicit finite caps and account for
sparse matched-relation storage and simultaneous dense combination work.
The factor-base lookup is immutable and cached. Prime-power marking uses
exact bounded lifts and retains a proved conservative fallback for singular
branches. Full SIQS checkpoint capacity can be configured up to 16 MiB;
the real 6.2 MB state resumes without resetting resources. New solver-prefix
digests encode wide binary masks in hexadecimal, bypassing decimal digit
limits while accepting legacy decimal prefixes. The final checkout streams
those exact same hexadecimal digest bytes; a compatibility test compares
its hash to the as-measured recursive encoding. The final codec change is
storage handling and is not represented as a measured timing improvement.

Raw captures, journals, full measured-source snapshots, harnesses, profiles,
resource/provenance checks and checkpoint sidecars are retained locally in
lossless hash-verified archives. Original control/cold calls record statistics
without writing checkpoint sidecars; primary balanced/varied/continuation
calls and every unfinished repaired-confirmation call persist full states.
The library checkpoint always retains the full bounded state. Sidecar write
cost is reported separately and added to caller elapsed time, while paused
time remains excluded from cumulative active budgets.

The repaired-source confirmation contains 76 validated measured outcomes and
12 validated warmups. The table reports completion and, when there is a
success, the descriptive median of that successful subset:

| Decimal digits | Matched wall/CPU allowance | SIQS completed; success median | ECM completed; success median |
| --- | --- | --- | --- |
| 30 | 30 s | 9/9; 3.666 s | 9/9; 0.136 s |
| 40 | 30 s | 0/9 | 9/9; 5.806 s |
| 50 | 60 s | 0/9 | 4/9; 17.600 s |
| 60 | 60 s | 0/9 | 0/9 |
| 70 | 300 s | 0/1 | 0/1 |
| 80 | 300 s | 0/1 | 0/1 |

SIQS unfinished calls reach a cooperative wall-time stop; ECM unfinished
calls exhaust their finite curve schedule. Checkpoint/finalization work can
extend caller elapsed time beyond the cooperative timer boundary. Each
unfinished confirmation call retains a hash-checked full checkpoint and
all outputs reconstruct. The single-attempt 70/80 bands establish bounded
exploration only. These caps are evaluation choices; increased cumulative
allowances can resume a saved job, as demonstrated by the longer 50-digit
experiment. Complete-factor scaling at 60–80 digits remains open.

**Decision: accept the bounded P3.4 implementation/evaluation and retain ECM
as the automatic default.** SIQS stays opt-in; no broad crossover, digit
cutoff or general 50–80-digit runtime guarantee is established. The current
checkout passes 230 PyPy tests and lint. Further matrix/arithmetic
scaling and fresh default-promotion comparisons remain separately gated.
The measured runtime is PyPy 7.3.23 implementing Python 3.11.15 on macOS arm64,
with one core per factoring call and no JIT override.

The ignored local evidence directory is
`v2/audit/phase_three_m31_complete_evidence_20261004/`, with a lossless `.tar.gz`
and external hash manifest alongside it. It preserves the v4/v5 source
snapshots, original relative capture/checkpoint paths, all runner harnesses
and independent validation. The two snapshots use the unchanged frozen
corpus SHA-256
`28ae182a8c935dd30cf51c5b7f586c63172b7bafd84e3f625f86d444422e19f3`.
The archive records 557 losslessly verified files. Existing training and
historical archives are retained. Generated evidence remains local under
the repository's ignore policy.

The filtering-only repair and saved-state continuation can be reproduced with
the following command in the retained measured snapshot:

```sh
pypy3 -m v2.benchmarks.phase_three_filter --frozen FROZEN_JSON \
  --checkpoint CHECKPOINT_JSON --output benchmarks/UNIQUE.json
```

The fresh 50-digit harness uses stdin to preserve its
recorded runtime command: `pypy3 -u < harnesses/factor_p34_m31_fresh_50.py`.
Recreate a Git checkout at the recorded metadata head, overlay its measured
source snapshot and use the recorded PyPy runtime. Git metadata and optional
lint dependencies are not bundled into the source archive. These harnesses
record their original absolute snapshot paths; preserve those paths or adapt
them explicitly in a new capture. Choose fresh output paths
for every replay. Raw timing captures identify uncommitted source by hashes;
the Git HEAD alone is not the measured implementation identity.

## P3.6 coarse SIQS workers (4 October 2026)

**Decision: retain serial collection.** The standalone experimental
`v2.qs.parallel.ParallelSIQSJob` evaluates independent SIQS families in serial,
2/4 threads and 1/2/4 spawned processes. The parent re-verifies every exported
atom, matches partials centrally and runs the common filtering/GF(2)/extraction
pipeline. Batches commit in assignment order. No dispatcher default changes.
P2.8's reused-process ECM throughput findings do not establish a SIQS benefit.

The retained independent `phase_three_p36_corpus.json` has four training and
four held-out inputs in each of two bands: 25–26-bit/eight-digit small inputs
and 42–44-bit/13-digit amortized inputs. The builder certifies its small prime
factors by exact trial division independently of runtime classification.
Its SHA-256 is
`b3e14c14e674f62c55e4dcb93ebb1162687e762c245ca4e586025045322a7c2f`.
Seeds are 7 and 29; each cohort contains eight input/seed attempts. Algorithms
receive n/configuration, with certificates used separately for validation.

Training-only complete-factor probes select bounds 200 and 400 for small and
medium inputs respectively, from 200/400 and 400/1000 candidates. Other declared
settings stay fixed: half-widths 256/512, one/three A factors, four families,
16 pool primes, 512 exported atoms per family, and 50 million work units per
assignment lease. The score policy is `powers`, recovery is `bucket`, residual
bound is 10000, and central atom/full/partial caps are 8192/4096/1024. Native
SIQS uses the same frozen mathematical schedule and collector settings; its
existing block-level extraction remains the complete-factor baseline.
No held-out result changes the configuration.

Per-input allowances are two billion total work units, five wall seconds,
five aggregate CPU seconds, 512 MiB total owned workspace and an observed
1024 MiB RSS gate. Parent/worker local caps are 64/32 MiB. Conservative owned
reservations include checkpoint encoding, duplicated bases and three result/IPC
copies per slot. Work leases reserve before submission and refund only reported
unspent work; cancelled work and all worker CPU remain charged. CPU checks are
cooperative, including an initial worker-startup publication and a final
post-transfer snapshot. Atomic operations and startup before publication can
overshoot a small allowance. Raw captures report cancellation drain, worker
work/CPU, retained batches and complete stage costs. PyPy threads retain the
GIL; no backend release is assumed.

PyPy 7.3.23 / Python 3.11.15 on macOS 26.6.2 arm64 runs a frozen local source
snapshot; all measured source hashes agree with the resulting checkout.
Forty-two configurations (four training and 38 held-out) retain 3,144 timed
attempts. Each configuration has at least three seconds of validated warmup
and nine samples. The noisy small native baseline extends to a further
five-second warmup and 15 samples. All final attempts pass the declared drift
and relative-IQR stability checks. The runner validates every warmup and timed
output; an independent retained audit rechecks reconstruction, certainty,
unique cohort assignments and resource accounting. No arithmetic failure occurs.

Warmed complete-factor attempt-cohort medians, including setup, transfer,
central re-verification, residual handling, solving, extraction and terminal
classification:

| Arm | Small median (ms), complete cases | Medium median (ms), complete cases |
| --- | ---: | ---: |
| Native serial SIQS | 19.4, 8/8 | 128.3, 8/8 |
| Coarse serial | 114.1, 8/8 | 421.6, 6/8 |
| 2 threads | 384.6, 8/8 | 1487.0, 6/8 |
| 4 threads | 825.7, 8/8 | 3452.3, 6/8 |
| 1 process | 239.3, 8/8 | 578.2, 6/8 |
| 2 processes | 414.8, 8/8 | 1302.4, 6/8 |
| 4 processes | 789.0, 8/8 | 2900.4, 6/8 |

Medium failures are the same declared family-batch capacity refusal for both
seeds of one held-out input. Their original cofactors reconstruct; no private
prefix is reported as committed. Times include these failures and therefore
do not establish successful-factorization speed ratios. On the fully complete
small cohort, descriptive bootstrap 95% time-ratio intervals versus native
serial are 5.49–6.04 for coarse serial, 20.27–21.80 for two processes and
39.06–40.92 for four. These are small-cohort timing intervals, not evidence of
a general large-integer crossover. No arm meets the 10% median improvement or
10-percentage-point completion improvement gate; medium completion regresses
25 percentage points, beyond the permitted five-point cross-class regression.

The separate fixed-work experiment defers extraction and keeps scanning after
direct residual splits. Every successful schedule covers exactly 2052/16400
positions for small/medium inputs respectively. Complete schedules have
identical work counts and scanned positions across modes and all timed samples.

| Fixed-work arm | Small cohort median (ms) | Medium attempt-cohort median (ms) |
| --- | ---: | ---: |
| Coarse serial | 630.8 | 1100.4 |
| 2 threads | 1289.7 | 3061.4 |
| 4 threads | 1310.6 | 4047.6 |
| 1 process | 1052.2 | 1422.4 |
| 2 processes | 1359.4 | 2646.3 |
| 4 processes | 1792.2 | 3738.0 |

Small schedules complete 8/8 in every sample; medium schedules complete 5/8,
with later family-batch cap refusals in addition to the first-factor refusal.
The medium mixed cohort is diagnostic rather than complete fixed-work
throughput evidence. Candidate throughput alone cannot promote a worker mode.

Cold worker arms create and close a fresh executor for every input/seed job;
the parent is warmed while each worker interpreter/JIT is cold. Startup and
shutdown are included. The warmed arms reuse executors and retain distinct
timing scope.

| Cold coarse arm | Small cohort median (ms) | Medium attempt-cohort median (ms) |
| --- | ---: | ---: |
| Serial | 125.6 | 425.6 |
| 2 threads | 396.1 | 1496.4 |
| 4 threads | 832.7 | 3376.2 |
| 1 process | 2539.3 | 3840.6 |
| 2 processes | 2952.4 | 4414.4 |
| 4 processes | 4640.2 | 8383.5 |

Maximum per-input accounted CPU is 3.737 seconds; the separate cold lifecycle
CPU observation stays below 3.751 seconds. Maximum conservative owned reservation
is 377.4 MiB. RSS observations sum parent/worker process-lifetime high-water
marks, reaching 609.4 MiB; these are conservative bounds rather than sampled
simultaneous peaks or per-job workspace estimates. The sandbox disallows `ps`,
so no live aggregate-RSS sampling claim is made. Preliminary captures using
that unavailable sampler or earlier evaluator bugs are excluded. Profiles and
instrumented RSS cohorts are not used as performance evidence.

Eleven new acceptance tests cover fixed assignments, serial/thread/spawned
equivalence, every returned atom's verifier, changed-worker checkpoint restart,
retained resources, incomplete-family replay, pending admission, malformed
provenance, capacity refusal, CPU/wall/work limits, safe failure cleanup and
fixed-work residual-GCD continuation. The shared checkout passes 226 tests
and lint. Completed pending batches retain their cursor; unfinished private
families replay under the same ID with newly charged work. Serialized resume
rebuilds and verifies stores and solver caches; checksums do not authenticate
edited history. Larger feasible SIQS controls, finer batches, separately trained
capacities and earlier stopping now belong to immediate P3.6.1 work;
broader production portfolio integration remains P6.3. This study closes the
bounded evaluator and does not establish practical 30–60-digit parallel factoring.

Run from the repository root; captures and frozen source archives stay local:

```sh
pypy3 -m v2.benchmarks.build_phase_three_parallel_corpus \
  --output /tmp/p36_corpus.json
cmp v2/benchmarks/phase_three_p36_corpus.json /tmp/p36_corpus.json
pypy3 -m v2.benchmarks.phase_three_parallel --phase all \
  --output v2/benchmarks/p36_LOCAL.json
make -C v2 test
make -C v2 lint
```

For separate runs, first use `--phase train --output TRAINING_LOCAL.json`, then
`--phase warm` or `--phase cold` with `--training TRAINING_LOCAL.json` and a
fresh `--output`. The Makefile also provides `benchmark-phase-three-parallel`.
The accepted local capture hash is
`e538111a7468c2e988bb15e3231ef78be14949a4f6a9a6deb8e2f0864a99f840`.

## P3.5/P3.6 cost diagnosis (4 October 2026)

A separate instrumented diagnosis uses the current checkout, the frozen
P3.5/P3.6 configurations and one already inspected input at seed 7 per band.
SSS/SSSf cover small and 30-digit inputs; native/coarse serial, one/two processes
and two threads cover the P3.6 small/medium inputs. Each of 16 configurations
has at least three seconds of validated warmup and nine stage-timed repetitions.
Two noisy configurations receive five-second warmups and 15-sample extensions
until stable. Parent and worker cProfile captures are separate from these stage
timers and from the accepted performance cohorts. Every result reconstructs,
including refusals. Loaded source hashes remain unchanged during diagnosis.
These inspected inputs are for attribution, not fresh promotion or retuning.
The current sources differ from the historical P3.5 capture, so these profiles
do not replace its recorded timing/completion evidence.

The SSS collision generator is the first optimization target. On the selected
30-digit input, coarse stage timers attribute about 79%/81% of SSS/SSSf elapsed
time to collision generation, versus about 3%/5% to smoothness trees and
12%/8% to common solve work. On the small input, collision generation accounts
for about 62%/72%; trees account for about 6%/8%. Inclusive profiler ranks also
show frequent budget validation and clock reads: the selected 30-digit SSS
split uses 581,327 `Budget.consume` calls, with collision generation rebuilding
affine roots and signed-shift counters for each selected/dropped-prime choice.
Profiler timings change PyPy execution and can change which budget stops a
run; they are not speed ratios or predicted savings.

SSSf also has a useful-yield issue. In one accepted small 64-case cohort,
its cutoff rejects 4,403 of 8,084 generated candidates and requires 418 search
assignments, versus 183 for SSS. The accepted 30-digit cohort uses 202 versus
156 assignments. Two-stage smoothness work and discarded candidates explain
why fewer candidates reaching recovery does not imply faster factoring.
P3.5 is not uniformly slower: its recorded 30-digit splits complete while
the original SIQS controls fail; comparison with a feasible trained control
now belongs to the immediate P3.6.1 follow-up.

P3.6 has both synchronization cost and extra work before extraction:

- `_ParentBudget.consume` calls `CollectionPool.poll` at every parent work
  charge. Each process-mode poll reads synchronized CPU arrays and publishes
  the parent CPU value. A selected small two-process profile has 2,494 polls
  and about 12,500 synchronized array-item access calls. Polling is a major
  inclusive parent cost, nested inside verification and matrix work.
- Workers test the multiprocessing cancellation event on every budget charge,
  before the one-millisecond CPU-publication throttle. The selected medium
  first family makes 12,619 such checks. The event uses synchronization even
  in serial/thread modes; thread profiles show its lock cost too. These
  profiles do not isolate GIL cost from lock contention.
- Native serial extracts after a 256-position block; coarse workers publish
  complete families. In one accepted small eight-case first-factor cohort,
  native serial scans 2,048 positions, coarse serial 4,104, and four processes
  15,842. Extra worker families can finish before the ordered parent discovers
  a split. Their work remains charged even when retained as pending batches.
- Verified worker atoms are checked and admitted centrally, then relations
  are verified again for matrix preparation. Larger batches also enlarge the
  matrix: on the selected small input, the accepted native/coarse first solve
  has 51/121 input rows. The serial coarse path therefore already loses before
  process startup or IPC is added. Cold startup remains an additional large
  measured cost; this diagnosis does not quantify pure serialization overhead.

The next experiments should amortize shared-state/clock checks over explicitly
bounded chunks, then try smaller verified-result batches and earlier extraction.
Preserve exact work charges, aggregate CPU, cancellation latency, refusal/replay
semantics and bounded cooperative overshoot. Tune capacities and batch sizes on
training inputs, and evaluate each change independently on fresh held-out inputs.
SSS collision arithmetic/counter construction, SSSf useful yield, and worker
accounting, granularity and stopping are owned by
[P3.6.1](../audit/TODOS.md#p361--immediate-diagnosis-and-improvement-of-p35p36).
That milestone can start now, before P3.8/P6.3, and records reasons, concrete
TODOs and acceptance/performance gates. Broader matrix/array reconciliation
and ECM/portfolio parallelism reuse its results later. Matrix or tree rewrites
are lower priority on these profiles. No runtime default changes.

The local diagnostic driver, stage samples, parent/worker profiles and source
manifest are retained under `v2/audit/p35_p36_diagnosis_LOCAL/`; generated
captures and scratch code remain Git-ignored. The accepted benchmark runners
above remain the commands for uninstrumented performance comparisons.

## P3.8 R1 capacity controls (4 October 2026)

R1 was developed in an isolated snapshot alongside the P2/P3/P3.6.1 repair
pass. The combined candidate passes 258 tests and lint. Its new runner pins source hashes and the independent corpus hash.
Measurements from this snapshot cannot establish the performance of a later
integration. The training population contains balanced 30/40/60/70/80/90/99-digit
semiprimes, uneven classes identified by smaller-factor sizes 5/10/20/30 where
applicable, and separately labelled p−1/p+1-smooth, close and square controls.
Primes have independently checked trial/Pocklington certificates; the
Pocklington construction's large p−1 factor is an explicit sampling bias.
Known factors enter validation only, never configuration selection or methods.

`build_p38_r1_corpus` generates training and confirmation separately.
Confirmation requires an already written frozen-control file, whose hash is
embedded in its corpus. The runner refuses mismatched source or freeze hashes.
Use a new corpus seed after control selection; repeated timing samples do not
create additional independent inputs. Keep the corpus files and selected
configuration file as immutable inputs. Generated captures, stdout, profiles
and scratch files stay local.

```sh
pypy3.11 -m v2.benchmarks.build_p38_r1_corpus \
  --seed 38120261004 --count 1 \
  --output v2/benchmarks/p38_r1_training_corpus.json
pypy3.11 -m v2.benchmarks.p38_r1_capacity --phase train \
  --corpus v2/benchmarks/p38_r1_training_corpus.json --seconds 3 \
  --confirmation-small-seconds 10 --confirmation-large-seconds 1 \
  --frozen v2/benchmarks/p38_r1_frozen.json \
  --output v2/benchmarks/p38_r1_training_LOCAL.json
pypy3.11 -m v2.benchmarks.build_p38_r1_corpus \
  --seed 38120261005 --count 1 --frozen v2/benchmarks/p38_r1_frozen.json \
  --output v2/benchmarks/p38_r1_confirmation_corpus.json
pypy3.11 -m v2.benchmarks.p38_r1_capacity --phase confirm \
  --corpus v2/benchmarks/p38_r1_confirmation_corpus.json \
  --frozen v2/benchmarks/p38_r1_frozen.json \
  --output v2/benchmarks/p38_r1_confirmation_LOCAL.json
```

The training grid jointly varies the factor base, interval, integer A target,
factor count and nearest-product/flyer policy, with fixed finite matrix and
store caps. Actual base cardinality, theoretical product envelope, selected A
bounds, assignment limits, residual bounds, retained rows, filter statistics
and owned workspace are reported. `factor_base_product_envelope_applies=False`
identifies an external MPQS coefficient, for which factor-base product bounds
do not constrain A. Search quotas, per-family Gray quotas and resident storage
are separate controls; a longer timer cannot expand a schedule.

Confirmation supports SIQS, external-square MPQS, ECM and ECM-to-SIQS under
matched input/seed/work/wall/CPU/owned-memory limits. Preprocessing is common
across arms. Every result reconstructs, including unresolved cofactors, and
retains probable/proven labels. Full-call timing includes setup, classification
and normal checkpoint serialization. Warmups execute validated full calls for
at least three seconds. Nine samples are extended to fifteen with a five-second
warmup when unstable, up to three attempts. Summaries retain instability and
censored outcomes; a failed run supplies no successful-factorization ratio.

Use `--kinds all` for the separately labelled varied classes. `--phase screen`
is diagnostic coverage, not repeated timing evidence. `--phase cold` records
nine independent interpreter lifecycles including import and corpus-proof
verification. `--phase profile` records separately instrumented cProfile data;
its durations are excluded from performance evidence. Process RSS is a
lifetime peak, not owned workspace or a simultaneous multi-process peak.
All current R1 arms use one process; no machine-dependent portfolio crossover
or new default is inferred from a small completed band or upper-band censoring.

Confirmation caps and seeds are selected in the frozen control before its
corpus is generated. The documented exploratory protocol uses 10 seconds at
30 digits and one second at upper bands. Those short upper-band runs measure
bounded coverage/costs and censoring only, not practical completion. Longer
studies require a separately frozen protocol and new confirmation inputs.
`--seconds` applies to training/diagnostic screens; confirmation, cold and
profile modes use the frozen per-band allowances.

A separate direct driver compares the unchanged legacy SIQS schedule with
an immutable pre-R1 runtime, keeping inputs, seeds, resource caps and full
output validation matched:

```sh
pypy3.11 v2/benchmarks/p38_r1_regression.py \
  --corpus v2/benchmarks/p38_r1_training_corpus.json \
  --baseline-json v2/benchmarks/p38_r1_baseline.json \
  --output v2/benchmarks/p38_r1_regression_before_LOCAL.json
pypy3.11 v2/benchmarks/p38_r1_regression.py \
  --corpus v2/benchmarks/p38_r1_training_corpus.json \
  --output v2/benchmarks/p38_r1_regression_after_LOCAL.json
```

This runner includes setup, splitting and child classification under a shared
30-second cap, with at least three seconds of validated warmup and nine samples
plus noise extension. It measures a single inspected 30-digit control at two
seeds; it is not an independent performance-promotion population.

### Recorded R1 confirmation and its limits

The frozen snapshot used PyPy 7.3.23 implementing Python 3.11.15 on arm64.
Training and confirmation each contain 58 disjoint inputs; the repeated
balanced comparison uses **one independent input per size band**, two seeds
(7 and 29), and nine two-seed cohorts. All 28 groups met the declared stability
rule on their first nine samples after at least three seconds of validated
warmup. CPU-heavy work in the parallel repair task was paused for these runs.

| Fresh 30-digit arm | Median two-seed cohort | Completed / attempts |
| --- | ---: | ---: |
| Streamed SIQS, flyer | 1.484 s | 18 / 18 |
| External-square MPQS | 1.198 s | 18 / 18 |
| ECM | 0.279 s | 18 / 18 |
| ECM then SIQS | 0.277 s | 18 / 18 |

These are full calls under the frozen ten-second per-call allowance. Repeated
seeds on one input do not establish population-level superiority. ECM finished
before SIQS was needed, so this case does not measure the cost of handing off
from an unsuccessful ECM search. At every balanced 40/60/70/80/90/99-digit band,
all four arms completed **0/18** and reported `wall_limit` under the prespecified
one-second allowance. Post-limit result/checkpoint construction is included in
full-call time; for example the 40-digit SIQS cohort was 2.159 seconds. No failed
arm establishes a successful time-to-factor ratio or a practical upper-size
capability. Upper-band choices with zero useful yield are provisional controls,
not evidence of optimized configurations.

At 30 digits, SIQS retained a median 218.5 verified rows per attempt and MPQS
207; their maximum owned reservations were 170,091,456 and 166,338,760 bytes.
The corresponding process lifetime RSS observations were 114,884,608 and
120,782,848 bytes, measured separately from conservative owned reservations.
At 40 digits, the median retained counts were 24.5 and 34, but the first sample's
singleton filter removed every row. Balanced 60–99-digit attempts retained no
useful rows. The first 99-digit SIQS sample did reach A=8.012e44 against the
integer target approximately 8.173e44; target reachability alone did not produce
a useful relation. Exact integers, filter dimensions/excess, coarse stage costs,
result certainty, work and observed A ranges remain in the local raw capture.

The prespecified A-policy comparison at the same trained 30-digit settings
completed all 18 attempts per arm: nearest **1.767 s**, flyer **1.444 s**, and
reference **1.902 s** median two-seed cohorts. Nearest and flyer share stream and
Gray quotas; reference retains its 64-family/full-Gray schedule. Flyer reduced
the nearest cohort time by 18.3% on this one input; defaults remain unchanged.
Factor count and Gray quota were not independently tuned by this experiment.

```sh
pypy3.11 -m v2.benchmarks.p38_r1_policies \
  --corpus v2/benchmarks/p38_r1_confirmation_corpus.json \
  --frozen v2/benchmarks/p38_r1_frozen.json \
  --output v2/benchmarks/p38_r1_policies_LOCAL.json
```

The fixed legacy-schedule control used the repair task's immutable handoff as
its before snapshot. Both arms completed 18/18 attempts on the inspected
training input: **1.701 s before, 1.752 s after** (3.0% slower). This is a small
observed regression, not a speedup claim; no default change follows from it.
This direct control excludes serialized checkpoint construction. Nine cold
interpreter lifecycles per arm measured **2.898 s** SIQS and **0.822 s** ECM median
two-seed calls, including import, corpus-proof checking and output. Separately
instrumented 30/99-digit profiles identify collection, budget polling and sieve
work as substantial costs; their overhead changes progress and is not timing
evidence.

The exact measured sources are retained in `p38_r1_measured_sources.json`,
including the imported benchmark helpers. The later extension-validation fix
also reserves verification-cache growth when raising relation capacity. It
changes `capacity.py` after the frozen capture; the measured arms already use
the maximum cache and never extend a checkpoint. The original source hashes
and captures are preserved rather than relabelled as measurements of a changed
runtime. Recreate that snapshot before replaying the recorded frozen protocol:

```sh
pypy3.11 -m v2.benchmarks.restore_p38_r1_capture /tmp/p38-r1-replay
cd /tmp/p38-r1-replay
# Run the confirm/policies/cold/profile commands above with new output paths.
```

All five R1 JSON inputs (two corpora, controls, pre-R1 baseline and measured
sources) are explicitly retained by `.gitignore`. Raw captures, stdout and
profiles remain local. The R1 roadmap stays open for broader independent
populations, Gray/reuse tuning, feasible upper-band yield and an actual
machine-dependent ECM-to-SIQS crossover.

The recorded ECM summaries do not instrument owned workspace: their zero
placeholder means unavailable, not zero allocation. The reported process RSS
still includes ECM and JIT memory. Integration added the repair task's worker
memory/resume guards and the R1 cache-extension validation correction; those
changes pass combined acceptance and are distinct from the immutable timed
snapshot. No runtime or benchmark-driver writes overlap frozen captures.
