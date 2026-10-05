# Benchmarks

## File layout and retention

- Python runners and documentation remain in this directory.
- `inputs/` contains versioned corpora, immutable baselines, frozen controls
  and provenance manifests required by tests and benchmark loaders.
- `results/` holds ignored captures, generated freezes, stdout, profiles,
  checkpoints and scratch subdirectories. Use a unique output name per run.

The entire `../audit/` tree is local and Git-ignored. Research notes, citation
manifests and historical diagnostic scripts stay there; required loader inputs
have moved here without changing their bytes. The public roadmap is
[`../ROADMAP.md`](../ROADMAP.md). Historical raw captures remain in the ignored
`.local-evidence/` archive. Publish concise findings here and accepted behavior
in the changelog; retain exact baseline bytes for reproducible comparisons.

Run from the repository root on PyPy implementing Python 3.11. Make creates
output folders automatically; for direct runner commands, create them first:

```sh
mkdir -p v2/benchmarks/results v2/audit/results
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
exact old-source baselines in inputs/, and the frozen P3.3 configuration used
by `phase_three_ownership.py`. Those snapshots are hash-checked before use.
The preserved v1 comparison is emulated Python 2 through `lib2to3`, not a
native Python 2 measurement.

## Historical v1 comparison (M9, 3 October 2026)

These are historical repaired-v2 measurements, not a fresh run of the current
implementation. PyPy 7.3.23 / Python 3.11.15 on macOS arm64, one serial process,
at least three seconds of validated warmup per arm and 15 rotated-order samples:

| Complete-factorization batch | Emulated v1 median | Repaired v2 median | Less elapsed time |
| --- | ---: | ---: | ---: |
| Original: five inputs, five seeds | 1.155 ms | 1.078 ms | 6.6% |
| Independent control: 56 inputs, five seeds | 21.091 ms | 16.994 ms | 19.4% |

The control contains 30 balanced small semiprimes, ten unbalanced composites,
thirteen primes and three squares. Generation seed is 20261005; factorization
seeds are 0–4. All timed answers have exact expected factors, multiplicities
and reconstruction; native answers also require proven terminal factors.
v1 runs through the syntax/integer-division compatibility adapter with
`math.gcd` on the same PyPy. This is not a native Python 2 timing comparison,
and changed defaults mean the full-run ratio does not isolate a single repair.
Invalid outputs receive no timing ratio. No large-number claim follows.

The archived M9 report and full captures remain in local historical evidence;
the M8 sources are retained in `inputs/m8_source_snapshot.json`.
Run a new comparison of the current source, without expecting the old timings:

```sh
mkdir -p v2/benchmarks/results
pypy3 -m v2.benchmarks.regressions --legacy \
  --warmup-seconds 3 --repetitions 15 \
  --output v2/benchmarks/results/v1_comparison_UNIQUE.json
```

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
  --output v2/benchmarks/results/qs_comparison_LOCAL.json
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
  --output v2/benchmarks/results/qs_audit_LOCAL.json
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
  BENCHMARK_OUTPUT=benchmarks/results/qs_families_LOCAL.json
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
  BENCHMARK_OUTPUT=benchmarks/results/siqs_comparison_LOCAL.json
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
  --output v2/benchmarks/results/sss_training_LOCAL.json
pypy3 -m v2.benchmarks.phase_three_sss --train-siqs-only --band small \
  --output v2/benchmarks/results/siqs_small_training_LOCAL.json
pypy3 -m v2.benchmarks.phase_three_sss --band small \
  --training v2/benchmarks/results/sss_training_LOCAL.json \
  --siqs-training v2/benchmarks/results/siqs_small_training_LOCAL.json \
  --output v2/benchmarks/results/sss_small_LOCAL.json
pypy3 -m v2.benchmarks.phase_three_sss --train-siqs-only --band balanced_30d \
  --output v2/benchmarks/results/siqs_30d_training_LOCAL.json
pypy3 -m v2.benchmarks.phase_three_sss --band balanced_30d \
  --training v2/benchmarks/results/sss_training_LOCAL.json \
  --siqs-training v2/benchmarks/results/siqs_30d_training_LOCAL.json \
  --output v2/benchmarks/results/sss_30d_LOCAL.json
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
  --output v2/benchmarks/results/upstream_sss_LOCAL.json
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
  BENCHMARK_OUTPUT=benchmarks/results/large_siqs_UNIQUE_pypy.json
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
  --checkpoint CHECKPOINT_JSON --output benchmarks/results/UNIQUE.json
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
cmp v2/benchmarks/inputs/phase_three_p36_corpus.json /tmp/p36_corpus.json
pypy3 -m v2.benchmarks.phase_three_parallel --phase all \
  --output v2/benchmarks/results/p36_LOCAL.json
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
[P3.6.1](../ROADMAP.md#p361--immediate-diagnosis-and-improvement-of-p35p36).
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
  --output v2/benchmarks/inputs/p38_r1_training_corpus.json
pypy3.11 -m v2.benchmarks.p38_r1_capacity --phase train \
  --corpus v2/benchmarks/inputs/p38_r1_training_corpus.json --seconds 3 \
  --confirmation-small-seconds 10 --confirmation-large-seconds 1 \
  --frozen v2/benchmarks/inputs/p38_r1_frozen.json \
  --output v2/benchmarks/results/p38_r1_training_LOCAL.json
pypy3.11 -m v2.benchmarks.build_p38_r1_corpus \
  --seed 38120261005 --count 1 --frozen v2/benchmarks/inputs/p38_r1_frozen.json \
  --output v2/benchmarks/inputs/p38_r1_confirmation_corpus.json
pypy3.11 -m v2.benchmarks.p38_r1_capacity --phase confirm \
  --corpus v2/benchmarks/inputs/p38_r1_confirmation_corpus.json \
  --frozen v2/benchmarks/inputs/p38_r1_frozen.json \
  --output v2/benchmarks/results/p38_r1_confirmation_LOCAL.json
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
  --corpus v2/benchmarks/inputs/p38_r1_training_corpus.json \
  --baseline-json v2/benchmarks/inputs/p38_r1_baseline.json \
  --output v2/benchmarks/results/p38_r1_regression_before_LOCAL.json
pypy3.11 v2/benchmarks/p38_r1_regression.py \
  --corpus v2/benchmarks/inputs/p38_r1_training_corpus.json \
  --output v2/benchmarks/results/p38_r1_regression_after_LOCAL.json
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
  --corpus v2/benchmarks/inputs/p38_r1_confirmation_corpus.json \
  --frozen v2/benchmarks/inputs/p38_r1_frozen.json \
  --output v2/benchmarks/results/p38_r1_policies_LOCAL.json
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

## P2/P3.1–P3.4 and P3.6.1 performance repair pass (4 October 2026)

The repair pass starts from `performance_audit_baseline.json`, a hash-checked
snapshot of the active uncommitted implementation, not the older Git HEAD.
`performance_audit_candidate.json` preserves the combined repair/R1 runtime.
`performance_audit_frozen.json` records the selected controls before generation
of `performance_audit_confirmation_corpus.json`. The initial inspected
population is `performance_audit_corpus.json`; the additional P2 training
classes are in `performance_audit_training_shapes.json`. All are explicitly
retained by `.gitignore`. Raw captures, stdout, cProfile files and the detailed
local journal remain ignored. No v1 runtime was changed.

These are independently certified constructed populations, not uniform RSA
samples. Trial/Pocklington proofs and every reconstruction are checked outside
and inside the relevant run paths respectively. Pocklington generation biases
p−1 toward a large certified factor. Training has three inputs per small,
medium and 30-digit class; fresh confirmation has four, plus four inputs in
each of five P2 regression classes (powers, uneven factors, three factors,
100-bit primes and a deliberately shared p−1-smooth small factor). One input
per 40/60/80-digit band supplies explicitly censored diagnostics. Seeds 7 and
29 are repeated within each cohort; nine timing samples do not create nine
independent input populations. All 35 fresh inputs are disjoint from the
inspected and training populations.

The direct runner separates setup-through-classification complete attempts,
fixed-work collector/matrix controls, and worker lifecycle cases. It requires
PyPy implementing Python 3.11, at least three seconds of validated workload
warmup and nine samples. Relative IQR above 20% or first/last three-sample
drift above 15% triggers a five-second warmup and fifteen samples, for at most
three attempts. An unresolved cofactor remains in every reconstruction;
probable and proven factor labels remain distinct. Short diagnostic limits
and failed schedules never supply successful-factorization speed ratios.

### Diagnoses and selected repairs

The initial profiles identified repeated prime-segment generator setup,
perfect-power root attempts, collision Counter/set construction, repeated
clock/cancellation calls, dense exponent recovery, repeated Hensel inversions,
full relation re-verification and dense matrix-label reservations. Profiles
are instrumented cost evidence only. Exact work accounting still reserves
before mutation; external polling is amortized in the QS/worker paths with
forced stage/publication checks. P2 retains strict polling: an isolated
prototype in `performance_audit_p2_poll_candidate.json` failed the 10%
improvement gate on all eight training classes (ratios 0.944–1.041).

Single-component reversion controls use the original source AST in the
current runtime; other repairs and R1 remain active. The medium training
cohort (three inputs, two seeds) measured:

| Component control | Reverted | Current | Interpretation |
| --- | ---: | ---: | --- |
| Collision assignment including setup | 0.053723 s | 0.014757 s | 72.5% lower |
| P2 complete calls, old prime cursor | 0.004739 s | 0.004116 s | 13.1% lower |
| Sparse-label matrix/kernel control | 0.006642 s | 0.004486 s | 32.5% lower |
| Bucket recovery | 0.017927 s | 0.016600 s | 7.4% lower in this control |
| Resieve recovery | 0.020878 s | 0.020651 s | No material time gain here |
| Bucket score/refinement | 0.018276 s | 0.016600 s | 9.2% lower |
| Skipped-small-prime score refinement | 0.019414 s | 0.018716 s | 3.6% lower |
| Preprocessing-only necessary-power screens | 0.002964 s | 0.002976 s | No gain at this size |

The sparse support changes also remove dense-base recovery and reservation
assumptions; their acceptance is not a claim that every collector microcase
passes a 10% timing threshold. R1 separately replaces dense-bit work charges
with source/destination-word and nonzero charges. Tests preserve every lifted
kernel, exact relation, work refusal and caller-owned input.

Training did not justify changing P2 chunk/GCD/rho settings. Chunk sizes 8/32,
GCD batches 64/256 and rho batch 128 stayed within roughly 6% of the default
across the eight classes; the 1 MiB schedule cache made the small/medium
classes 40%/27% slower. The 30-digit cohort completed 5/6 attempts under every
arm, so its elapsed time is a bounded cohort cost, not a full-factor speedup.

SSSf selection/filter training retains the accepted filtered arm as a control
and selects six primes with filtering disabled as a separate challenger.
On the inspected 30-digit training cohort, all six attempts completed:
SSS 1.981 s, filtered seven-prime SSSf 3.239 s, six-prime unfiltered SSSf
2.129 s. Raising SSS's residual bound from 10,000 to 25 million increased the
same cohort to 2.437 s without improving completion, so the fixed default is
retained. These measurements precede the final native polling/R1 integration
and are selection evidence, not the final combined performance claim.

Worker training retains complete-family tasks by default. Width 128 improved
the small cohort but lost on medium; width 256 did not win across both.
Optional chunks preserve family/Gray/block identity and permit a finite
result cap to be met without dropping an unpublished prefix. Reservations
now bound the possible exponent support by the polynomial norm, not the
entire factor base. Resume rejects assignment/storage geometry changes,
permits work-lease/poll changes and smaller checkpoint caps, and releases
paused solver/extractor references to drained pools.

### Preparation, polling, transport and cooperative limits

The pre-tree-decision runtime is retained separately in
`performance_audit_pre_tree_decision.json`. On its feasible 30-digit training
cohort, **all 54 timed attempts per arm completed** (three independent inputs,
two seeds, nine cohorts). Current preparation/cache/polling took 7.570 s per
cohort; reverting only preparation reuse took 13.313 s, and reverting only
native polling took 15.187 s. Relative IQRs were 2.52%, 0.37% and 1.22%;
first/last drift was 2.00%, 0.20% and 1.26%. These isolate 43.1% and 50.2%
reductions against those respective controls. They do not add together.
The later full-base SSSf restoration does not affect this SIQS comparison.

The disjoint second SSSf base was rejected. A repeat with **ten seconds of
validated warmup and fifteen samples** measured the inspected small filtered
cohort at 0.139307 s with the full base versus 0.149328 s disjoint (+7.2%);
the six-prime unfiltered cohort was 0.024169 versus 0.023685 s (−2.0%).
The 30-digit filtered control was 2.481 versus 2.495 s. No arm cleared the
10% improvement gate. Both filtered small arms completed 5/6 attempts;
the six-prime unfiltered arms completed 6/6. Final source restores the full
second base. The initially generated, unused confirmation population was
retired locally; seed 361072026 generated the final fresh population after
this last decision.

A separate worker profile used the same medium input, first family and work
lease. Poll intervals 1 and 64 both scanned 4,100 positions, produced 248
verified atoms, charged 421,605 work units and reported 10,261,200 owned bytes.
Cancellation callback counts dropped from 10,675 to 168. Pickled output was
36,756 bytes in both cases. A 256-position chunk returned 13 atoms/2,318 bytes
with 33,685 charged units and 8,741,872 owned bytes. This is a smaller task,
not a claim that it does equal mathematical work faster.

`performance_costs.py` also echoes four identical 13-atom payloads per sample
with full payload equality checks. Serial encode/decode was 0.260 ms; warmed
1/2/4-process round trips were 0.955/0.944/1.268 ms. Thread-1/2 measurements
remained unstable after all three allowed extensions; their medians are not
accepted speed evidence. Thread-4 was 0.156 ms. Thread transport shares Python
objects while process transport serializes them. These controls include
queue/dispatch/wait/equality work, exclude factoring and cold startup, and do
not isolate the GIL or establish a share of full-factor latency.

`performance_limits.py` exercises warmed serial, 1/2/4-thread and
1/2/4-process pools under a 20 ms wall allowance and, separately, a 20 ms
aggregate CPU allowance. Every attempt reconstructs, reports the expected
limit, drains the pool and stays within its exact work ledger. All fourteen
cases became stable (two extended to fifteen samples). Across timed wall-cap
attempts, the largest return latency was 28.301 ms; across CPU-cap attempts,
the largest aggregate CPU was 31.179 ms. These are observations, not hard
overshoot guarantees. Startup is excluded here; atom/chunk costs, publication
throttling, parent waits and drain/barrier work remain cooperative. Reported
RSS values elsewhere are process-lifetime high-water marks, not simultaneous
resident peaks or enforced aggregate RSS caps.

### Fresh confirmation and decisions

Final captures use PyPy 7.3.23 implementing Python 3.11.15 on arm64. All
35 runtime hashes match `performance_audit_candidate.json` before and after
each capture. The table compares the original active-source snapshot with
the final combined repairs/R1 snapshot, using identical inputs, seeds and
resource limits. Values are seconds per eight-attempt cohort (four inputs,
two seeds), including setup and terminal classification. Every attempt in
these rows completed; all timing series passed the stated stability check.

| Complete cohort | Before | After | Reduction | 95% timing ratio interval | 95% input ratio interval |
| --- | ---: | ---: | ---: | --- | --- |
| Small P2 portfolio | 0.001809 | 0.001470 | 18.7% | 0.798–0.831 | 0.812–0.837 |
| Medium P2 portfolio | 0.006888 | 0.005980 | 13.2% | 0.852–0.889 | 0.846–0.863 |
| Medium native SIQS | 0.086587 | 0.035436 | 59.1% | 0.398–0.430 | 0.412–0.441 |
| Medium coarse serial executor | 0.318386 | 0.063434 | 80.1% | 0.184–0.210 | 0.190–0.202 |
| Medium SSS | 2.617038 | 0.760354 | 70.9% | 0.274–0.299 | 0.266–0.284 |
| Medium filtered SSSf | 6.211471 | 1.560663 | 74.9% | 0.245–0.256 | 0.231–0.247 |
| 30-digit feasible SIQS | 23.324077 | 6.839935 | 70.7% | 0.292–0.298 | 0.290–0.297 |
| 30-digit SSS | 6.180587 | 1.411689 | 77.2% | 0.224–0.238 | 0.221–0.229 |
| 30-digit filtered SSSf | 9.460341 | 2.019916 | 78.6% | 0.212–0.215 | 0.208–0.221 |

`performance_summary.py` separately resamples cohort timings 10,000 times and
enumerates paired input-cluster draws, keeping both seeds together. The
second interval uses sums of per-input/seed median times, so its center can
differ from the ratio of cohort medians. Neither interval turns four
constructed inputs into a broad population result. Completion differences
are zero in every class, including the censored P2 cases; empirical zero-width
completion intervals do not establish a zero population failure rate.

The five fresh P2 shape classes all complete 72/72 timed attempts per arm.
Ratios are 0.733 for powers, 0.847 for uneven factors, 0.913 for three factors,
0.934 for primes and 0.870 for the p−1-smooth class. Their respective timing
intervals are 0.704–0.747, 0.841–0.855, 0.906–0.920, 0.922–0.941 and
0.845–0.877. The 30-digit P2 portfolio completes 7/8 input/seed pairs in both
versions; 1.568→1.564 s is a censored cohort cost, not a speedup claim.

The separately trained six-prime, unfiltered SSSf arm takes 0.023086,
0.765811 and 1.445306 s on small/medium/30-digit cohorts. Against filtered
SSSf its ratios are 0.721/0.491/0.716, with timing intervals
0.648–0.750/0.444–0.504/0.702–0.721. It completes every timed attempt and
passes the improvement gate against that filtered control. It shows no
10% improvement over SSS (ratios 1.066/1.007/1.024). Retain this explicitly
labelled challenger and the historical filtered control; retain SSS as the
simpler control. On the declared feasible 30-digit comparison, SSS/SIQS is
0.206 (timing interval 0.202–0.215; paired-input interval 0.197–0.213), with
72/72 complete attempts each under equal 30-second wall/CPU, 10^13-work and
128 MiB outer allowances. The methods use their separately trained internal
capacities. This supports the optional SSS arm on this cohort; automatic
dispatch and a general sub-100-digit crossover remain unestablished.

### Worker throughput, startup and capacity

The separate final worker capture includes native serial, coarse serial,
1/2/4 threads and 1/2/4 processes. All first-factor attempts complete.
Warmed cohort medians are:

| Execution | Small | Medium |
| --- | ---: | ---: |
| Native serial | 0.008360 | 0.038709 |
| Coarse serial executor | 0.019361 | 0.059659 |
| Coarse threads 1 / 2 / 4 | 0.022259 / 0.036320 / unstable | 0.061008 / 0.164774 / 0.692656 |
| Coarse processes 1 / 2 / 4 | 0.069284 / 0.071724 / 0.115838 | 0.103412 / 0.161021 / 0.199135 |
| Chunk-256 serial executor | 0.012642 | 0.055042 |
| Chunk-256 processes 1 / 2 / 4 | 0.025138 / 0.060029 / 0.063828 | 0.141802 / 0.149872 / 0.123388 |

Small four-thread coarse and chunk runs remained unstable after all three
extensions and are excluded from speed decisions. Keep native serial as
the default and workers experimental. Chunks reduce excess collection but
trade it for repeated setup, verification and scheduling. In the first
medium timed sample, native serial scans 19,472 positions, coarse four-process
execution scans 130,382 (including 1,137,221 cancelled work units), and
chunked four-process execution scans 22,034. These illustrative counts are
not timing evidence or a promise of race-independent first-factor work.

Fixed schedules with a 1,024-atom batch cap finish all assignments across
serial, 1/2/4 threads and 1/2/4 processes. Every timed sample has identical
admitted-atom signatures and charged work for each input/seed pair, with no
pending assignments. Median cohort times are 0.323374 s serial,
0.319210/0.443979/0.896044 s threaded and 0.571900/0.441136/0.412573 s
process-based. More processes help relative to one process but do not beat
serial on this fixed workload.

Cold worker lifecycle samples create a fresh pool per eight-attempt cohort
inside the timed region, including shutdown; the parent is warmed and the
new children are cold. Coarse 1/4-process medians are 0.552951/0.708752 s
small and 0.828436/1.142542 s medium. Chunked medians are
0.420948/0.526257 s small and 0.830081/0.737185 s medium. These are not
per-input cold interpreter measurements. Per-job CPU counters exclude the
caller's final pool teardown; wall timings include it.

On the already-inspected historical medium P3.6 corpus, the repaired
512-atom coarse control completes 3/8 attempts; the other five explicitly
return `batch_limit`. This differs from the historical 6/8 because collector
yield has changed. Increasing the finite cap to 1,024 or publishing
256-position chunks completes 8/8 in every sample; native serial also
completes 8/8. Cohort costs are 0.118466/0.157470/0.137491/0.141313 s for
512-cap/1,024-cap/chunks/native respectively. The incomplete control is not
a successful-factorization speed baseline. Retain larger caps and chunks
as explicit capacity choices, with their own reserved-memory costs.

The single-change reservation control uses the same four-process pool,
10,000 base bound, 1,024-atom cap and 256-position assignment. The frozen
dense-support reservation refuses memory before dispatch in all nine
samples. The polynomial-norm support bound admits the assignment in all
nine, scans 256 positions, then pauses as requested, with 293,643,584 bytes
of accounted peak owned storage. This is a capacity repair, not a claim that
refusing less work is slower. Worker and parent quotas still include queued
copies and central verification storage.

### Larger diagnostics and acceptance scope

The fresh 40/60/80-digit 0.2-second probes complete no factors. Legacy SIQS
exhausts useful yield after 65,600 scanned positions with zero admitted atoms;
its fixed A-product envelope misses the target. SSS/SSSf exhaust wall time.
At 60/80 digits, the larger 82,000 base often consumes the allowance before
the first assignment finishes. Raising the outer cap from 64 to 128 MiB
does not remove this setup/time constraint. A 25-million residual cap admits
more atoms at 40 digits (640–888 in the diagnostic attempts), yet still
produces no factor. These are useful-yield/capacity diagnostics, not evidence
for a default change. R1's streamed/external-square polynomial choices and
its longer-cap experiments are recorded separately above.

A subsequent single-seed, two-second diagnostic (not warmed timing evidence)
reaches beyond setup. All 15 outputs reconstruct and retained atoms reverify.
The 40-digit SSS/128 MiB arm finds a proper split; the other 14 runs reach
wall limits, with no memory/tree-cap refusals. At 60 digits the four default
residual arms admit 17–26 atoms, whereas the larger residual admits 296 and
one matched pair. At 80 digits they admit 0–1 atoms versus two with the
larger residual, with no matches. At 40 digits the larger residual admits
6,037 atoms but evicts 4,420 partials and finds no split. More raw yield
therefore does not demonstrate more useful dependencies. A single
order-dependent diagnostic split cannot establish a memory-cap benefit.

Adopt the exact accounting, collision, sparse recovery, inverse reuse,
preparation and reservation repairs. Retain P2 batching, strict P2 clocks,
serial dispatch and coarse worker defaults; reject the disjoint SSSf tree
and automatic residual-bound increase. All 259 tests and full lint pass.
A temporary committed candidate containing only the 159 retained files also
passes the full suite, lint and all 40 benchmark-module imports; the main
Git index/history is unchanged. Arithmetic, reconstruction, checkpoint,
provenance, exact work, owned-storage and cancellation gates pass.

This completes the bounded P3.6.1 repair experiments, including explicit
negative decisions. Repeated filtering is now a substantial remaining
profiled cost; R3 owns its separately frozen identity/cadence/provenance
experiments. Practical larger-band scaling remains open. None of these
measurements claims that every remaining algorithmic bottleneck is removed.

Reproduce controls with new capture paths; the runner refuses overwrites:

```sh
pypy3.11 v2/benchmarks/performance_audit.py \
  --baseline-json v2/benchmarks/inputs/performance_audit_candidate.json \
  --corpus v2/benchmarks/inputs/performance_audit_confirmation_corpus.json \
  --bands small,medium --cases p2_full,sss,sssf,native,parallel_serial \
  --output v2/benchmarks/results/p361_confirmation_NEW_LOCAL.json
pypy3.11 v2/benchmarks/performance_audit.py \
  --baseline-json v2/benchmarks/inputs/performance_audit_pre_tree_decision.json \
  --split training --bands 30d --cases siqs_feasible --revert preparation \
  --output v2/benchmarks/results/p361_preparation_control_NEW_LOCAL.json
pypy3.11 -m v2.benchmarks.performance_costs \
  --output v2/benchmarks/results/p361_costs_NEW_LOCAL.json
pypy3.11 -m v2.benchmarks.performance_limits \
  --output v2/benchmarks/results/p361_limits_NEW_LOCAL.json
pypy3.11 -m v2.benchmarks.performance_capacity \
  --output v2/benchmarks/results/p361_capacity_NEW_LOCAL.json
pypy3.11 v2/benchmarks/performance_audit.py \
  --baseline-json v2/benchmarks/inputs/performance_audit_candidate.json \
  --corpus v2/benchmarks/inputs/performance_audit_confirmation_corpus.json \
  --bands 30d --cases siqs_feasible,sss_feasible,sssf_feasible \
  --output v2/benchmarks/results/p361_30d_NEW_LOCAL.json
pypy3.11 -m v2.benchmarks.performance_summary \
  --before v2/benchmarks/results/p361_confirm_baseline_30d_LOCAL.json \
  --after v2/benchmarks/results/p361_confirm_30d_LOCAL.json \
  --output v2/benchmarks/results/p361_summary_NEW_LOCAL.json
```

For fresh raw captures, run each `performance_audit.py` command again with
`performance_audit_baseline.json` and a distinct output path, then pass those
two paths to `performance_summary`. The illustrative summary command uses
local capture names, which are intentionally absent from a clean checkout.

`--baseline-json performance_audit_baseline.json` restores the original control;
`performance_audit_p2_poll_candidate.json` restores the rejected P2 prototype.
`--revert` isolates one frozen component without writing runtime files.
The SSS constructor control explicitly binds the original `super` owner when
loading its AST; a harness regression test checks that superclass setup runs.
`--profile` is a separately instrumented whole cohort after the timed samples;
its often much longer elapsed time is never included in their medians.

### R1 varied-shape screen and final integration acceptance

A separate one-second diagnostic screen replayed the exact measured sources
across all 58 confirmation inputs, both seeds and four arms: **464 validated
outcomes**, comprising 277 complete factorizations and 187 explicit wall-limit
results. It has no repeated warmed timing claim; its one-second cap also differs
from the ten-second 30-digit confirmation cap. Complete counts are stratified
below; factor-one counts equal complete counts in this semiprime/power corpus.

| Diagnostic class | SIQS | External MPQS | ECM | ECM then SIQS |
| --- | ---: | ---: | ---: | ---: |
| Balanced, seven size bands | 0/14 | 2/14 | 2/14 | 2/14 |
| Uneven, 5-digit smaller factor | 14/14 | 14/14 | 14/14 | 14/14 |
| Uneven, 10-digit smaller factor | 0/14 | 2/14 | 14/14 | 14/14 |
| Uneven, 20-digit smaller factor | 0/10 | 0/10 | 0/10 | 0/10 |
| Uneven, 30-digit smaller factor | 0/8 | 0/8 | 0/8 | 0/8 |
| p−1-smooth control | 14/14 | 14/14 | 14/14 | 14/14 |
| p+1-smooth control | 14/14 | 14/14 | 14/14 | 14/14 |
| Close factors | 3/14 | 4/14 | 5/14 | 5/14 |
| Square | 14/14 | 14/14 | 14/14 | 14/14 |

The p+1 control's smaller prime is also p−1-smooth, and its recorded completions
occurred in common p−1 preprocessing. These labels are not exclusive and do not
evaluate a p+1 algorithm. Close-factor controls are not guaranteed to fit the
fixed 1,000-step Fermat allowance. Structured completions at large total sizes
must not be read as balanced-semiprime capability or a crossover estimate.

```sh
# From the restored measured snapshot:
pypy3.11 -m v2.benchmarks.p38_r1_capacity --phase screen --kinds all \
  --seconds 1 --corpus v2/benchmarks/inputs/p38_r1_confirmation_corpus.json \
  --frozen v2/benchmarks/inputs/p38_r1_frozen.json \
  --output v2/benchmarks/results/p38_r1_shapes_LOCAL.json
```

The final combined runtime includes the repair task's measured rejection of
its optional SSSf disjoint-base optimization. A matched legacy-schedule bridge
on that exact final runtime completed **18/18**, with **1.718 seconds** median
two-seed cohort time after 4.902 seconds of validated warmup and nine stable
cohorts. This checks integration on an already inspected control; it is not new
held-out performance evidence. Its 35 runtime hashes match the repair task's
final freeze.

A temporary committed candidate containing only the 159 explicitly retained
project files passed **259 tests**, full lint and imports of **all 40 benchmark
modules**. Required immutable snapshots/corpora therefore work without ignored
raw evidence. Tests ran under PyPy implementing Python 3.11; the lint tools used
the existing development environment. The clean candidate, source manifest,
raw shape/bridge captures and check logs remain local. Main runtime and drivers
were unchanged during these captures; CPU-heavy runs were serialized with the
repair and R3 tasks.


### P3.8-R3 stable rows, preparation and provenance experiments

R3 froze the already integrated R1/repair source before changing identities or
adding caches. `p38_r3_baseline.json` contains all 35 hash-checked runtime
modules, rather than a moving import of the control. The subsequent repair
rollback affects only SSSf, which these SIQS comparisons do not exercise.
`p38_r3_measured_sources.json` preserves the first matrix/preparation capture;
`p38_r3_final_sources.json` preserves the native opt-in policy capture.
The frozen policy file also retains the final drivers with shared-pipeline
budget and coexistence checks. The profile driver and raw captures remain
local. Snapshots are data; restore
only their hash-checked source bytes into a temporary directory when replaying
that revision. Required inputs and runners are retained by explicit ignore
exceptions.

Training has two independent certified semiprimes in each of four classes,
with seeds 7 and 29. Fresh confirmation has three per class, generated only
after the policy freeze. The `20d` and `30d` labels denote nominal 20–21 and
30–31-digit ranges; actual lengths are stored in the fixtures' integers.
Pocklington generation biases p−1 toward a large certified factor. These tiny
cohorts establish scoped experiments, not a population ranking or a general
SIQS crossover. Known factors/certificates validate outputs and never enter
algorithm configuration. Full-call measurements include independent split
validation and budgeted certainty classification; every censored cofactor also
reconstructs the input.

Use PyPy implementing Python 3.11. All arms share 10^10 work units and 128 MiB
owned storage; wall and CPU caps are 0.2 seconds per smaller attempt and two
seconds per nominal 30-digit attempt. Matrix/preparation calls allow 30 seconds.
The 30-digit configuration uses base bound 10,000, width 8,192, four A factors,
32 candidates, 64 families and 4,096-position blocks. Smaller configurations
use bounds 200/1,000/3,000 and widths 256/512/2,048. The runner prints the sample
stability gate, validates at least three seconds of warmup for every arm and
interleaves nine samples; unstable sets extend to five seconds and 15 samples,
up to three attempts. The revised short comparison repeats the same two-input,
two-seed cohort five times within each sample, revalidating all 20 results.
Bootstrap intervals describe repeat uncertainty conditional on these fixed
inputs; they are not population confidence intervals. Censored cohort elapsed
time is never reported as time-to-factor improvement. Cold process launches
include startup/import/setup and first split, and are labelled separately from
warm complete-factor results. CPU-heavy phases ran serially with the other
sessions; profiles do not supply performance evidence.

The initial separately instrumented feasible 30-digit profile recorded 22
collection calls (1.450 seconds inclusive), 21 preparations (0.095 seconds) and
21 filters (0.205 seconds). Small extraction calls totalled about 0.8 ms.
Changed-prefix preparation covered 8/16/32/64/128 rows (or all available rows),
with 120 checked-cache hits on the first three stores. The sparse large stores
had one or zero admitted rows. This profile does not establish extraction as
the dominant cost or practical large-matrix capacity.

| Matrix/complete-solver fixture | Repeated rebuild | Frozen queues | R3 counters | Decision |
| --- | ---: | ---: | ---: | --- |
| 512-row cycle | 51.979 ms | 2.998 ms | 2.937 ms | Retain existing touched queues |
| 512-row cascade | 32.850 ms | 1.673 ms | 1.690 ms | Retain existing touched queues |
| 512-row dense | 34.489 ms | 35.849 ms | 34.555 ms | No broad queue win |
| 1,536 rows / 512 pivots | 3.518 ms | 10.434 ms | 4.344 ms | Adopt incremental pivot count |
| 512 sparse-label cycle | 776.203 ms | 8.563 ms | 8.531 ms | Retain R1 initial gap mapping |
| 128 verified full rows, all kernels extracted | 1.575 ms | 1.629 ms | 1.493 ms | No representation promotion |

The rebuild arm loads only the frozen M26 filtering function into otherwise
matched current limits; it is a function control, not a whole old-source
factoring claim. Synthetic rows use an independent set-based rank oracle; every
lifted mask is checked against original rows. Actual relation matrices also
extract every dependency through the original modular verifier. The separate
mixed runner interleaves 64 full and 64 matched rows and reserves preparation,
matrix and compaction/history storage simultaneously. Conversion, incidence,
fill, solver and deferred lifting all count in timing. The initial prefix comparison refreshed the budget at each stage; its
2.046/2.187 ms uncached and 2.040/2.219 ms R3 medians are diagnostic only. The
subsequent single-ledger preparation/filter/solve/extract repeat measured
1.821/2.076 ms uncached, 1.954/2.137 ms for the frozen cache and
1.863/2.085 ms for R3. Its peak owned reservations were 8.23/4.84 MiB including
the full reserved cache. It showed no cache win; the existing bounded preparation cache is retained on the repair pass's
separately documented larger repeated-call evidence.

- **Adopt:** stable mixed admission indices, complete payload identities,
  fully verified versioned checkpoint replay and the incremental pivot counter.
- **Reject as defaults:** live-column compaction, disjoint pivot batching and
  immutable merge histories. No consistent end-to-end benefit appeared;
  batching is conflict-free merging and supplies no independent-dependency
  claim. Dense masks remain the production provenance representation.
- **Defer:** replacing dense quadratic reservations, larger-matrix scaling,
  accumulated modular roots and packed exponents. The history prototype caps
  rows/columns at 4,096, nodes at twice its row count, and reserves the existing
  dense workspace plus history/scratch. No constant was reduced to admit an
  otherwise refused workload. Original parity/provenance remain independently
  checkable; incremental residues cannot attest their own history.

```sh
pypy3 -u -m v2.benchmarks.p38_r3 --phase matrix \
  --output v2/benchmarks/results/p38_r3_matrix_NEW_LOCAL.json
pypy3 -u -m v2.benchmarks.p38_r3 --phase prepare \
  --output v2/benchmarks/results/p38_r3_prepare_NEW_LOCAL.json
pypy3 -u -m v2.benchmarks.p38_r3_mixed \
  --output v2/benchmarks/results/p38_r3_mixed_NEW_LOCAL.json
pypy3 -u -m v2.benchmarks.p38_r3 --phase training \
  --output v2/benchmarks/results/p38_r3_training_NEW_LOCAL.json
pypy3 -u -m v2.benchmarks.p38_r3 --phase confirmation --split held_out \
  --bands small,medium,20d --variants control,current,cadence8 \
  --output v2/benchmarks/results/p38_r3_confirmation_NEW_LOCAL.json
pypy3 -u -m v2.benchmarks.p38_r3 --phase confirmation --split held_out \
  --bands 30d --variants control,current,cadence32 \
  --output v2/benchmarks/results/p38_r3_confirmation30_NEW_LOCAL.json
pypy3 -u -m v2.benchmarks.p38_r3 --phase cold \
  --output v2/benchmarks/results/p38_r3_cold_NEW_LOCAL.json
```

The committed training/held-out corpora are immutable. `--phase build` creates
new fixtures at a new path; it never overwrites retained corpora. Held-out
creation requires the frozen policy file and embeds its SHA-256. Captures
refuse overwrite and record runtime, control, driver, corpus and source hashes.


The mixed supplement found 83/92 verified kernels and 44/32 proper divisors
respectively, with identical results in every arm. Current medians were
1.472/1.787 ms, live compaction 1.592/1.876 ms, dense batch-32
1.435/1.752 ms and history-32 1.501/1.835 ms. Reservations including the
retained full provenance were 6.26/2.87 MiB for dense masks and
6.57/3.18 MiB for histories. The small batching changes do not meet a consistent
10% gate and history adds storage. Equal parity alone never deduplicates rows.

Native training froze cadence 8 for medium and nominal 20-digit experiments,
and cadence 32 for nominal 30-digit experiments. The production default stays
1. Training medium complete-cohort time was 36.229 → 32.042 ms (11.6%; paired
95% repeat interval 7.7–14.8%). The 20-digit completion rate changed from
47/60 to 56/60 (+15 percentage points; conditional interval +8.3–23.3).
Its elapsed-time reduction is censored and is not a successful time-to-factor
claim. The 30-digit complete-cohort comparison was 1.990 → 1.460 seconds
(26.6%; interval 26.4–32.2%), 36/36 complete in both arms. The small cohort
showed no qualifying gain; no class lost completion. Cadence-32 20-digit timing
remained unstable after extensions, so it was not selected for that class.
Pair stability requires both compared arms to pass; an unstable rejected arm
does not close or invalidate a different candidate's gate.

Tested-dependency caching stays off. Against R3 with caching off, the cache arm
adds roughly 1.7% on small, 5.8% on medium and 4.2% on the complete 30-digit
training cohort; it does not pass the 10% time or 10-point completion gate.
Full payload hashing/lookup costs are included. Its capped implementation is
available explicitly for future reuse experiments, and it shares the original
2 MiB verification-cache reservation. Cache refusal falls back to full trials;
checked trivial skips still count against the same finite logical retry limit.
Public and checkpoint preparation always reverify untrusted data.


Fresh held-out confirmation passed the stability and zero-error gates in every
compared arm. All short outcomes were complete (450/450 per arm, 15 samples of
30 attempts); medium and 30-digit outcomes were 54/54 per arm (nine six-attempt
cohorts). Every warmup was validated, including the extended short runs.

| Held-out class | Frozen control | R3 cadence 1 | Frozen candidate | Complete outcomes: control / R3 / candidate |
| --- | ---: | ---: | ---: | --- |
| Small, five repeated cohorts | 42.450 ms | 42.205 ms | cadence 8: 41.560 ms | 450 / 450 / 450 |
| Medium | 81.620 ms | 71.483 ms | cadence 8: 64.570 ms | 54 / 54 / 54 |
| Nominal 20-digit | 810.340 ms | 739.495 ms | cadence 8: 575.072 ms | 36 / 47 / 54 of 54 |
| Nominal 30-digit | 7.745 s | 7.526 s | cadence 32: 4.611 s | 54 / 54 / 54 |

Against the frozen control, the joint R3/cadence-8 medium arm reduces complete
cohort time 20.9% (conditional 95% interval 16.1–23.0%); the nominal 20-digit
completion gain is 33.3 points. Against R3 with cadence 1, cadence 8 changes
medium time by 9.7% and held-out 20-digit completion by 13.0 points. Its isolated
training gains were only 5.6% medium time and 8.3 points 20-digit completion.
Thus **defer a standalone cadence-8 promotion** despite the joint R3 option's
training/confirmation gains; do not attribute identity/order gains to cadence.
The small class has no qualifying gain and retains cadence 1.

Cadence 32 independently passes the declared nominal 30-digit gate: against
R3 cadence 1 it reduces complete-cohort median time **26.9% in training and
38.7% held out** (held-out interval 38.0–39.0%), with no completion loss. Against
the original frozen control the held-out reduction is 40.5% (39.7–40.8%).
**Adopt it only as an opt-in setting for this explicit feasible configuration.**
All other classes retain the production setting unless the caller explicitly
requests a candidate. This is no automatic size rule or parameter optimum.
Censored 20-digit elapsed times support no successful time-to-factor claim.

Cold first-split process medians for the two short training fixtures were
179.270/167.900 ms control and 179.490/167.281 ms R3, nine fresh PyPy processes
per cell. They reconstruct all outputs and validate proper splits against
certificates; they exclude complete-factor classification. No cold gain is
claimed. The acceptance suite covers corrupted/rehashed payloads, missing or
invalid mixed permutations, frozen legacy pending checkpoints, skipped-trial
allowances, cache saturation/refusal, exact brute-force kernels, inverse maps,
solver refusal and terminal/resumed opt-in behavior.


R3 integration acceptance preserves the repair/R1 shared budget, capacity,
known-square, portfolio and full-base SSSf interfaces. The final combined tree
passes **273 PyPy tests**, full lint and imports of **all 44 benchmark modules**
from a temporary committed-files-only snapshot. The isolated R3 commit also
passes 272 tests/lint and 42 clean imports. Raw evidence and temporary acceptance
snapshots remain local; required loaders use only explicitly retained inputs.

A serial matched integration bridge on the already inspected two-input,
two-seed 30-digit training cohort completes 36/36 in each arm after at least
three seconds of validated warmup and nine stable samples. Cohort medians are
1.995 seconds frozen control, 2.008 seconds R3 cadence 1 and 1.428 seconds
cadence 32. All 35 final runtime hashes match the bridge capture; the only
runtime difference from the isolated R3 source is the preserved repair SSSf
rollback, which SIQS does not exercise. This is integration verification on
inspected inputs, not a new held-out result or a change to the frozen policy.

## Follow-up: restrictive allowances and repeated setup

This second pass starts from the integrated post-R3 runtime, not the original
P3.6.1 control above. `performance_followup_baseline.json` retains its 35
modules; `performance_followup_experiments.json` retains hash-checked single
experiments, rejected prototypes and the exact `selected` implementation.
`performance_followup_frozen.json` fixes the seven changes, configurations,
measurement rules and new seed before confirmation-input generation. No
configuration was tuned on those fresh outcomes.

The audit covers lease fragmentation/refunds, paused parent solvers, disabled
stage allocations, duplicate live storage, resieving scratch, local versus
process synchronization, repeated base verification and identity serialization.
It also tests ECM baby-step batching, bounded Hensel-value caching, direct
valuation for single-hit primes and coherent CPU-counter reads. Existing P2
batch/GCD/polling controls, polynomial selection, sparse recovery, SSS trees,
preparation and matrix work remain documented in the first pass and R1/R3
sections. Library kernels have no automatic disk/network I/O; explicit CLI
output and portable checkpoint serialization are separate caller operations.

Accepted implementation behavior:

- Clip an assignment lease to available work. Wait for running leases to
  refund before refusing pending relation admission or solver continuation.
  Incomplete private work remains charged and replays under the same ID.
- Use local events/counters in serial/thread pools. Reuse their parent's
  checked immutable base, while serialized process ingress and returned
  relations retain full verification. No cross-job base cache is introduced.
- Retain one unchanged family digest per immutable base. This is a bounded
  repeated-setup optimization, with a separate kernel measurement below.
- Exclude disabled p−1 attempts and zero-curve ECM tiers from unused schedule
  and baby-table reservations; count shared family/collector base storage once.
- Bound resieving exponent support by the smallest-prime product fitting the
  actual `supported_A * abs(F(x))` norm bound. Reserve scratch before allocation,
  release it after the block, and preserve refusal/retry without admission.

### Fresh confirmation

The new seed is **361082026**, disjoint by integer value from 14 earlier
corpora. The retained corpus contains 35 certified inputs; these 22 paired
cells use **12 independent inputs** (four small, four medium, four nominal
30-digit) and seeds 7/29. The remaining generated shapes/larger fixtures are
not additional confirmation evidence. Each timing invocation isolates one
band/case/runtime in a fresh PyPy 3.11 process, alternating arm order. Worker,
native and fresh-pool cells use ten seconds of validated warmup and 15 samples;
capacity/setup/30-digit cells use at least three seconds and nine samples.
Unstable cells extend up to three attempts. Parent/worker profiles are not
used as timing evidence.

Times below are seconds per **eight-attempt cohort**, including required
output validation. Every timed attempt in these rows completes, 120/120 per
arm except nominal 30-digit SIQS at 72/72.

| Matched case | Post-R3 control | Selected | Interpretation |
| --- | ---: | ---: | --- |
| Small serial worker | 0.016483 | 0.014448 | 12.3% lower time |
| Medium serial worker | 0.024138 | 0.021066 | 12.7% lower time |
| Medium B=10000, chunked serial worker | 0.282389 | 0.191570 | 32.2% lower time |
| Medium fixed schedule, serial | 0.091003 | 0.092165 | No speedup claim |
| Medium fixed schedule, four threads | 0.470945 | 0.095560 | 79.7% lower fixed-work time |
| Medium fixed schedule, four processes | 0.135285 | 0.134768 | No speedup claim |
| Medium fresh serial pool lifecycle | 0.023474 | 0.019551 | 16.7% lower time |
| Medium fresh two-process pool lifecycle | 0.395658 | 0.399599 | No speedup claim |
| Nominal 30-digit feasible SIQS | 4.659316 | 4.595209 | No additional whole-factor speedup claim |

Serial/thread/process fixed schedules have identical atom signatures, exact
work totals and factors in both arms. The four-thread fixed-work improvement
does **not** make it faster than serial (0.095560 versus 0.092165 s). First-factor
two/four-thread controls remain unstable after extensions; their apparent
improvements are excluded from acceptance. Native serial negative controls
differ by −2.1% small and +1.1% medium. Automatic worker selection and native
serial defaults remain unchanged.

The 95% warmed timing ratio intervals (selected/control) are 0.867–0.894 for
small serial, 0.847–0.898 for medium serial, 0.656–0.712 for B=10000 chunks,
0.179–0.231 for the four-thread fixed schedule and 0.812–0.847 for fresh serial
pools. Paired input-cluster intervals for the first three are respectively
0.871–0.894, 0.859–0.889 and 0.663–0.681. These measure separate uncertainty
sources, not independent additional inputs. The SIQS timing interval
0.968–1.007 and fresh process-pool interval 0.987–1.040 support no speedup.
Fresh-pool timing includes create/run/close for one pool per cohort on a
warmed parent; it is not cold parent-interpreter startup. Earlier multi-case
medium comparisons showed large negative-control drift and are excluded.
Even isolated warmed intervals do not capture all between-process JIT or host
variation, and these small cohorts do not establish large-input scaling.

Capacity outcomes below are complete factors per eight-attempt cohort, except
the explicitly labeled collection row. Refusing faster is not a speedup.

| Matched finite allowance | Control | Selected |
| --- | ---: | ---: |
| Small serial, 1M total work, 10M assignment ceiling | 0/8 | 8/8 |
| Medium serial, same work limits | 0/8 | 6/8 |
| Medium four threads, 15M total work, 10M ceiling | 0/8 | 8/8 |
| Medium four processes, same work limits | 0/8 | 8/8 |
| Disabled huge p−1/zero-curve ECM bounds, small | 0/8 | 8/8 |
| Disabled huge bounds, medium | 0/8 | 8/8 |
| B=10000 chunked worker, 10 MiB worker cap | 0/8 | 8/8 |
| B=10000 resieve, 4096-position block, 64 MiB cap | 0/8 collected | 8/8 collected |

The selected resieve scans all 4096 positions each time, verifies recovered
exponents, and peaks at **23,159,008 owned bytes**; no factor-time claim follows
from this collection result. All spent work remains within the same limits.
The medium 1M-work misses remain explicit unresolved results. All other
declared fresh completion classes retain their completion rate.

### Cache scope, rejected experiments and reproduction

Repeated complete family setup improves 0.255556→0.138015 s on fresh medium
inputs (46.0%); every CRT term/inverse output and identity is checked. This
case initializes 64 families per attempt and is **not factorization**. An
additional training-only ablation removes just the digest cache from the
selected implementation: B=10000 chunked full-factor time is
0.124531→0.112137 s, 90/90 complete in each arm. Its 9.95% reduction does not
independently pass the 10% whole-factor promotion threshold. Retain the scoped
setup improvement without claiming a separate general factorization gain.
The whole-worker bundle above passes its declared complete-factor gate.

Single-change training rejects extra Hensel-value memoization (less than 1%
whole-factor benefit and slower collector cells), ECM baby batching (less than
4% across tested P2 shapes), and direct valuation for single-hit primes (5.4%
on the nominal 30-digit case). Fresh-process CPU-read batching improves the
tested medium cells by 3.7–9.9%, below the gate, and remains a prototype.
Keep the existing bounded derivative-inverse cache, exact P2 scheduling,
process synchronization, full-base SSSf pass and finite arithmetic/storage
allowances. No polynomial/search-policy, residual-cap or backend default is
promoted by this follow-up. Larger useful-yield and matrix-provenance scaling
remain the separate R1/R3/R4 work described above.

Run each case/arm separately; `none` selects the post-R3 control and `selected`
selects the exact accepted source overlay. Lease case names explicitly override
the generic runner's work header with the 1M/15M allowances recorded in the
`followup` metadata. Other configurations remain matched.

```sh
pypy3.11 v2/benchmarks/performance_followup.py --variant none \
  --corpus v2/benchmarks/inputs/performance_followup_confirmation_corpus.json \
  --split held_out --bands medium --cases parallel_base10k_chunk_serial \
  --warmup-seconds 10 --repetitions 15 \
  --output v2/benchmarks/results/followup_control_LOCAL.json
pypy3.11 v2/benchmarks/performance_followup.py --variant selected \
  --corpus v2/benchmarks/inputs/performance_followup_confirmation_corpus.json \
  --split held_out --bands medium --cases parallel_base10k_chunk_serial \
  --warmup-seconds 10 --repetitions 15 \
  --output v2/benchmarks/results/followup_selected_LOCAL.json
pypy3.11 -m v2.benchmarks.performance_limits \
  --output v2/benchmarks/results/followup_limits_LOCAL.json
make -C v2 test
make -C v2 lint
```

All 14 warmed wall/CPU limit cells pass stability and reconstruction checks,
without work overdraw or a locked/undrained pool. With a 20 ms allowance, the
largest final-sample wall-limit return is 24.580 ms; the largest aggregate CPU
use in CPU-limit cells is 21.206 ms. These are observations under cooperative
polling, not hard deadlines. The final sources pass **282 PyPy tests**, full
lint, and a separate committed-input-only checkout's tests and **44 benchmark
runner imports** (45 modules including the package). All 35 runtime hashes
match the frozen confirmation source; `v1/` is untouched. Raw captures,
profiles, local summaries and temporary checkouts remain excluded from Git.

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
