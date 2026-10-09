# Factor benchmarks

Validated factoring experiments on **PyPy implementing Python 3.11**. This guide
keeps the stage history, accepted changes and rejected experiments concise.
The [v2 guide](../README.md) covers usage; the [roadmap](../ROADMAP.md) records
remaining acceptance gates.

## A10 verified primality transfers (9 October 2026)

The frozen control is `94caf40`; required inputs are
[a10_before_sources.json](inputs/baselines/a10_before_sources.json),
[a10_primality.json](inputs/corpora/a10_primality.json) and the
[frozen protocol](inputs/corpora/a10_protocol.json). Expected truth comes from
independent trial, Pocklington and full n−1 Lucas proofs, or exact composite
divisors. These are test oracles; production certificate generation stays B14.

The supported sets are deliberately small:

| Strict upper bound | Bases | Guarantee/source |
| --- | --- | --- |
| `9080191` | 31, 73 | [Jaeschke 1993, p. 926](https://cr.yp.to/bib/1993/jaeschke.pdf), exhaustive computational limit |
| `4759123141` | 2, 7, 61 | Same primary computation |
| `2**64` | 2, 325, 9375, 28178, 450775, 9780504, 1795265022 | [Sinclair's 2011 record](https://miller-rabin.appspot.com/), also stated in [Forišek–Jančina 2015](https://ceur-ws.org/Vol-1326/020-Forisek.pdf); finite computational guarantee |
| `318665857834031151167461` | First 12 primes, through 37 | [Sorenson–Webster, Theorem 1.1](https://arxiv.org/abs/1509.00864), exhaustive computational result |
| `3317044064679887385961981` | First 13 primes, through 41 | Same theorem |

The final two endpoints equal `399165290221 * 798330580441` and
`1287836182261 * 2575672364521`. Tests independently check their strong
congruences; endpoint equality never uses the preceding set. The paper's
conjecture concerns search complexity, not these finite results. This task
checks the primary reasoning and independent counterexamples/proof fixtures;
it does not repeat the original enormous exhaustive searches. v1's higher
entries and log-based rules have no accepted unconditional guarantee here.

Relevant implementation comparison (versions pinned, no external code copied):

| Reference | Certainty and practical choices | PyPy transfer / license review |
| --- | --- | --- |
| [GMP 6.3.0 release](https://gmplib.org/download/gmp/gmp-6.3.0.tar.xz), [API](https://gmplib.org/manual/Number-Theoretic-Functions) (`pprime_p.c`, `millerrabin.c`) | Trial division, BPSW, then `max(0, reps-24)` random MR tests; return 0/1/2 separates composite/probable/definite. Default definite BPSW survivors have `n < 35*2**46`; a build macro extends this to `n < 2**64`. Neither certifies arbitrary integers. | Preserve our explicit random-round API. Native modular powers are a separate backend. Source is LGPL-3-or-later / GPL-2-or-later dual licensed. |
| [FLINT 3.6.0](https://github.com/flintlib/flint/blob/v3.6.0/src/ulong_extras/ll_is_prime.c), [proof path](https://github.com/flintlib/flint/blob/v3.6.0/src/fmpz/is_prime.c) | Exact Sorenson–Webster cutoff in double-limb MR; shared decomposition, small filters and native multi-base powering. Larger exact proofs use n±1 methods and APRCL; probable entry points stay separate. | Range dispatch transfers directly. Limb reducers, instruction parallelism and hashed tables need their own Python evidence. LGPL-3-or-later. |
| [PARI/GP 2.19.0 source release](https://pari.math.u-bordeaux.fr/pub/pari/unix/pari-2.19.0.tar.gz), [API](https://pari.math.u-bordeaux.fr/dochtml/html-stable/Arithmetic_functions.html#isprime) | Machine-word MR/BPSW dispatch (`n < 2**64` on a 64-bit build); `ispseudoprime` remains probable above that checked domain. `isprime` proves via n−1, APRCL or ECPP; `primecert` supplies checkable evidence. | Proof work belongs to B14. Keep Boolean helpers separate from certainty-bearing results. GPL-2-or-later. |
| [SymPy 1.14.0](https://github.com/sympy/sympy/blob/sympy-1.14.0/sympy/ntheory/primetest.py) | Small filters, hashed/fixed witnesses and the same 12/13-base bounds; GMP availability can bypass MR for strong BPSW. `_test` exits on a nontrivial square root of 1. | Compare `mr` with identical bases, not `isprime`'s backend shortcut. Test early exits and powering independently. BSD-3-Clause; copying would require notices. |
| [zmwangx/miller-rabin, `34cd694`](https://github.com/zmwangx/miller-rabin/tree/34cd694b34fe916dc58705f42c4dbc4001b7990a) | CPython C extension: 16-bit lookup, hashed 64-bit witnesses, preliminary division/Fermat tests and GMP powering; larger results are probabilistic. | Its CPython/GMP timings do not predict PyPy integer performance. Own code MIT, with separate bundled/GMP notices. |

[Mishra's implementation blog](https://neelmishra.github.io/blog/cp/number-theory-2/miller-rabin.html)
and [GMP implementation discussion](https://gmplib.org/list-archives/gmp-devel/2018-November/005073.html)
helped identify filters, powering and BPSW candidates; correctness rests on
primary papers and pinned source. BPSW is useful screening, not an arbitrary-size
proof. The blog's seven consecutive prime-base claim through the 13-base bound
is refuted by `341550071728321`, checked independently in the tests; no such
range is adopted. [AKS](https://annals.math.princeton.edu/2004/160-2/p12) gives an
unconditional general deterministic algorithm, but has no measured advantage
for this bounded Python workload. APRCL/ECPP and n−1 certificates add proof
capability rather than replacing a requested random test. No external
implementation was adapted; license notices remain with ignored local source
captures, and future source copying requires the stated obligations.

Shared `n-1 = d*2**s` decomposition, first-witness rejection, small division
filters and three-argument modular powering already transfer well to Python.
The measured challengers isolate direct built-in powering, Python binary
powering, early rejection when a square reaches 1, extra prime/square/GCD
filters and an always-13-base wider test. Hashed witness tables require a
separately verified table and memory/setup comparison; native limb reducers,
assembly and triple-base instruction parallelism are not assumed to transfer
to PyPy. Broader preprocessing/filter promotion remains E3. Screening by
BPSW followed by all required MR bases could preserve certainty, but adds a
Lucas implementation and its validation cost; BPSW alone cannot replace those
fixed bases or requested random rounds here.

Commands from the worktree root (PyPy 3.11 only):

```sh
pypy3 -m v2.benchmarks.a10_inputs
pypy3 -m unittest v2.tests.test_a10_primality -v
pypy3 -m v2.benchmarks.a10_primality --split training --arms control accepted native_pow binary_pow early_one extra_filters square_filter primorial_filter thirteen_bases --output v2/benchmarks/results/a10/training.json
pypy3 -m v2.benchmarks.a10_primality --split confirmation --arms control accepted --output v2/benchmarks/results/a10/confirmation.json
pypy3 -m v2.benchmarks.a10_primality --split training --factoring --arms control accepted native_pow binary_pow early_one --output v2/benchmarks/results/a10/factoring-training.json
pypy3 -m v2.benchmarks.a10_primality --split confirmation --factoring --output v2/benchmarks/results/a10/factoring.json
pypy3 -m v2.benchmarks.a10_primality --split confirmation --factoring --arms control accepted --warmup-seconds 8 --repetitions 45 --batch-seconds 0.5 --output v2/benchmarks/results/a10/factoring-stability.json
pypy3 -m v2.benchmarks.a10_primality --split confirmation --within-range --arms accepted gmp_same_bases sympy_same_bases --output v2/benchmarks/results/a10/references.json
pypy3 -m v2.benchmarks.a10_primality --cold --output v2/benchmarks/results/a10/cold.json
pypy3 -m v2.benchmarks.a10_primality --analyze v2/benchmarks/results/a10/factoring.json --output v2/benchmarks/results/a10/factoring-analysis.json
```

The runner acquires the shared nonblocking machine lock. Coordinate all heavy
checks too. Each arm receives at least three seconds of validated warmup;
seeded interleaving pairs at least nine approximately 100-ms samples. It
extends unstable groups by nine through 45 and reports any unresolved spread.
Cold process startup/import/testing and instrumented profiles are separate.
Matched factoring uses seeds 7/104729/130363, 1,000,000 work units, trial 30,000,
no rho/p−1/ECM attempts, and finite input/storage caps. Training/confirmation
numbers are disjoint except the declared reported regression. The control
labels wider primes probable; upgraded certainty is a capability change,
reported separately from speed. Native same-base timings include conversion
and are contextual evidence, not proof of a faster Python implementation.
Conditional 95% intervals use 10,000 paired-round bootstrap draws with seed
20261009; fixed inputs and starts are not resampled. The optional speed gate is
at least 10% lower complete-run median with an interval excluding zero, then
disjoint-input confirmation and no completion loss. All capture paths must be
fresh; raw evidence and detailed research stay local. Select the fastest stable
qualifying complete-training candidate before confirmation, or retain the
range-only baseline if none qualifies. Other confirmation arms are contextual.

**Measured decision:** retain the range-only implementation, existing small
filters, shared decomposition, modular-power helper and first-failing-witness
exit. No optional challenger qualified. PyPy 7.3.23 / Python 3.11.15 on macOS
26.6.2 arm64 ran in the coordinated A10 machine window, after B1 released and
before B2 started. Every warmup and sample validated results. The corpus has
365 independently checked classification cases, 283 proof nodes and 26
factoring fixtures; timing uses its frozen training/confirmation subsets.

| Warm cohort (whole cohort per call) | Frozen control median | Accepted median | Observed saving; conditional 95% interval | Samples / stability |
| --- | --- | --- | --- | --- |
| Confirmation: eight wider primes | 1.946 ms | 0.672 ms | 65.45%; 65.19–65.65% | 9; both stable |
| Complete training: 13 inputs × three seeds | 10.275 ms | 6.975 ms | 32.12%; 31.86–32.42% | 45; tails remain |
| Complete confirmation: 14 inputs × three seeds | 11.262 ms | 7.243 ms | 35.68%; 35.47–35.97% | 45; tails remain |
| Prespecified longer complete confirmation | 11.105 ms | 7.431 ms | 33.08%; 32.09–34.06% | 45; control remains unstable |

The longer follow-up was frozen in
[a10_stability_protocol.json](inputs/corpora/a10_stability_protocol.json)
before its outcome: eight seconds validated warmup and 500-ms batches. Its
range/median is 20.44% for control and 10.80% for accepted (limit 15%); initial
complete confirmation was 25.10% / 84.37%. Every sample is retained. These
conditional median intervals do not establish stable complete-run latency or
universal optimality. This prime/power-heavy workload disables rho, p−1 and
ECM; E1 still owns combined portfolio confirmation with B1/B2 changes.
Control wider survivors were probable after 40 random witnesses; accepted
ones receive deterministic 12/13-base proof. Their improved guarantee and
reduced witness/RNG work are intentional, not an equal-certainty shortcut.

Against accepted complete training, direct built-in powering saved 1.58%
(1.04–1.96%), early-square-one rejection saved 1.11% (0.54–1.55%), and Python
binary powering cost 6.34% (5.93–6.80%). None was stable or met the 10% gate.
The null selection was recorded before confirmation; contextual confirmation
then made all three slower than accepted. Always using 13 wider bases cost
4.24% on training wider primes. Extra 41–97 trial filters reduced their
43-divisible microcohort from 0.170 to 0.0149 ms (tails remain); a square check
reduced squares from 0.454 to 0.0337 ms (stable). Those specialized wins lack
complete-portfolio promotion evidence and belong to E3. A primorial GCD cost
38.86% on already cheap small-divisor inputs. No filter was promoted here.

Same-base native references used gmpy2 2.3.1 / GMP 6.3.0 and SymPy 1.14.0
with `ground_types=gmpy` / `gmpy2.mpz`. All inputs were inside the supported
range, with identical fixed witnesses and independently expected truth.
For the eight wider primes, medians were accepted 0.679 ms, GMP 0.466 ms and
SymPy 0.513 ms; all reference groups extended to 45 samples and retained tail
spread. The smaller-prime subgroup was 0.193 / 0.279 / 0.186 ms, showing
conversion/bridge costs can reverse the ordering. These are contextual native
backend observations, not a Python-backend adoption decision or a comparison
of unequal BPSW/random guarantees.

Nine cold processes (startup, import, one reported-prime classification)
had medians 33.59 / 33.72 ms for control/accepted. Separate instrumented
profiles used three seconds validated warmup per arm, then ten full-corpus
calls: 12,900 versus 6,570 MR witnesses across 780 bounded runs. Profiles
explain reduced classification work; their timings are not warmed evidence.
Raw captures, source hashes, selection records and profiles remain in the
ignored local `results/a10/` and `audit/a10/` directories.

Validation covers strict endpoints/neighbors, prime powers, Carmichael and
strong pseudoprimes, exact reconstruction/labels/unresolved cofactors,
explicit round/RNG spies, finite work/time/cancellation and checked resume
under legacy schemas 4/5/6 on int/GMP. `make -C v2 test` passed 369 tests
(two optional-GMP skips), the GMP-enabled repeat passed 372, and lint passed.
One earlier GMP suite failed the unchanged QS snapshot-lifetime assertion;
its isolated test and full repeat passed. Its cause is undetermined and
referred to B1; no QS implementation change is included in A10.

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
verification was pending at this first-study report boundary.
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
their selected backend. These timing captures precede the P4.1/P5.2 integration
and identify their measured sources; integration adds no new timing claim.
The user subsequently authorized merging the results and removing the worktree.

### Combined mainline acceptance — 5 October 2026

The committed-files-only integration candidate combines the backend foundation
with P4.1 PRAC and P5.2 reusable programs. It passes 358 PyPy/GMP tests, full
lint and all 57 benchmark imports. Both P4.3 proof-backed corpora and the
hash-verified P4.3/P5.2 loaders work without local captures. Independent PRAC
point checks, program resume on each backend and old backend checkpoint
compatibility are included. Native portfolio schemas 4/5 remain unchanged;
combined GMP snapshots use 6. Existing ladder, streamed and integer defaults
remain in place. Integration supplies no additional performance measurement.

All 12 local captures, including explicitly diagnostic and unfinished sets,
are copied and checksum-verified into ignored
`results/p43/worktree_results/` before managed-worktree removal. The preservation
manifest and combined-check logs stay in `results/p43/`. Required corpora,
baselines and runners are committed; generated evidence is never force-added.

## P4.1/A4 verified PRAC — 5 October 2026

`ecm.multiply_prac` now executes bounded verified records, including checked
exceptional recovery; it is no longer an alias for the binary ladder. This
accepts A4 correctness. Production standalone and bounded factoring retain
the ladder; B3 owns portfolio work charges, chunk replay and checkpoints.
Near-optimal/offline search is still a separate comparison, not an accepted
speedup.

The independent affine corpus has 16,016 comparisons and zero mismatches.
391 selected chain paths require checked recovery; these are not vacuous
`(0,0)` equalities or a recount of the original prototype's 797 failures.
Additional tests cover exhaustive small curves, fields through 521 bits,
composite and prime-square moduli, Suyama prime-power schedules, retained
nonunit factors, corrupt records and bounded storage/termination.
The exact false-infinity example documented in GMP-ECM's `ecm.c`
(`n=33554520197234177`, `sigma=2046841451`, B1=373) also passes the
independent affine oracle and finishes at a non-infinite point.

### Complete ECM attempts on 40–80-digit inputs

`p41_campaign.py --gmp` compares four arms: the actual `factorize_ecm` ladder,
checked PRAC with Python integers, the same ladder arithmetic with gmpy2,
and checked PRAC with gmpy2. These gmpy2 arms run Python ECM, not the
C GMP-ECM executable or its optional near-optimal chain interpreter.
The gmpy2 adapter uses private function bindings with identical
function code, replacing integer validation, GCD and inversion; it retains
`mpz` coordinates through setup, both stages and recovery. There is no global
backend patch or per-operation conversion. This benchmark adapter is not the
P4.3 production backend contract. It requires the optional dependency in the
same PyPy 3.11 interpreter and fails explicitly if unavailable.

The versioned `p41_ecm_40_80_corpus.json` contains 14 certified semiprimes:
exact 40/50/60/70/80-digit inputs, balanced cases, and controlled 10-/20-digit
small factors. Prime certificates are recursively checked with Pocklington
and trial-division leaves, independently of the ECM implementation. Fixtures
were selected before examining ECM outcomes. Repeated small factors across
sizes deliberately control factor difficulty; these are inspected fixtures,
not independent population samples. Known factors enter validators only.
The timer covers the known-composite ECM attempt; certificate and result
verification run outside it. Portfolio preprocessing, recursive dispatch and
checkpointing are not included.

All arms receive the same nine seeds, Suyama parameters, B1/B2 bounds, batch
size 128, maximum 32 curves and 60-second wall/CPU watchdogs. Every warmup and
sample is validated, preserving unresolved composites and reconstructing all
proper-factor results. gmpy2/Python pairs must also match every uncensored
factor and curve/stage transition. Each case/arm receives at least three
seconds of validated warmup. Groups with more than 25% timing spread repeat
the same nine paired seeds in reversed order, retaining the original samples.
Tables report the final nine-seed block for extended groups; the original
block remains in the evidence. Fresh-process cold attempts are separate.
Each cold sample repeats seed 41,001 on the balanced 40-digit case. Its
curve count can differ from the warmed table's median across nine seeds;
cold and warmed medians must not be subtracted to estimate startup overhead.

The `current` tier uses B1=2,000 and B2=147,396 on all 14 inputs. `factor20`
uses B1=11,000 and B2=1,873,422 on the five cases with 20-digit small factors.
The latter is a published GMP-ECM reference tier, not a calibrated optimum
for this Python engine. [Zimmermann's parameter table](
https://members.loria.fr/PZimmermann/records/ecm/params.html) distinguishes
factor size from input size and gives both newer and historical estimates.
32 curves are a finite comparison budget, not a claim of sufficient coverage
for every 20-digit factor. Balanced 50–80-digit cases are unresolved-work
controls at the current tier; they do not demonstrate practical extraction
of 25–40-digit factors with these bounds.

The experimental PRAC campaign builds its program once per attempt and
reuses it across curves. Construction, caller-record re-verification, dispatch,
all intermediate checks, stage-two prime generation and saturation recovery
are timed. Its cap is B1 <= 11,000 and 4,096 record references, in addition to
the scalar compiler's 512-entry cache and 512-instruction limit per record.
It avoids cache thrashing by retaining immutable bound-owned records and
replays a collapsed prime power only from that power's saved starting point.
This experimental composition does not modify `stage_jobs` or its accounting.

Reproduce using the same PyPy environment for all arms (the accepted gmpy2
configuration is gmpy2 2.3.1 / GMP 6.3.0):

```sh
v2/.venv/bin/python -B -u -m v2.benchmarks.p41_campaign --gmp --tier current --output v2/benchmarks/results/p41-current-new.json
v2/.venv/bin/python -B -u -m v2.benchmarks.p41_campaign --gmp --tier factor20 --output v2/benchmarks/results/p41-factor20-new.json
```

### Accepted current-bound results

All four arms find proper factors in the same **49/126 distinct seeded
attempts**. There are zero timeouts, invalid factors or mismatched backend
transitions; 432 GMP/Python pairs agree across primary and extended blocks.
All five 10-digit-factor fixtures complete for 9/9 seeds; the four unbalanced
20-digit-factor fixtures complete for 1/9 each. None of the balanced inputs
splits in nine 32-curve attempts at these bounds.

Median seconds per complete attempt, using the final nine-seed block where
extension was triggered (success counts are identical across the four arms):

| Input case | Factor digits | Successes / 9 | Int ladder | Int PRAC | gmpy2 ladder | gmpy2 PRAC |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Balanced 40 | 20 + 20 | 0 | 0.2536 | 0.3602 | 1.2179 | 1.3099 |
| Small-factor 40 | 10 + 30 | 9 | 0.0111 | 0.0226 | 0.0577 | 0.0576 |
| Balanced 50 | 25 + 25 | 0 | 0.2762 | 0.3981 | 1.2264 | 1.3240 |
| Small-factor 50 | 10 + 40 | 9 | 0.0120 | 0.0243 | 0.0558 | 0.0585 |
| Target-factor 50 | 20 + 30 | 1 | 0.2398 | 0.3645 | 1.2230 | 1.3306 |
| Balanced 60 | 30 + 30 | 0 | 0.3275 | 0.4728 | 1.2925 | 1.4235 |
| Small-factor 60 | 10 + 50 | 9 | 0.0150 | 0.0265 | 0.0579 | 0.0592 |
| Target-factor 60 | 20 + 40 | 1 | 0.3043 | 0.4525 | 1.2379 | 1.3697 |
| Balanced 70 | 35 + 35 | 0 | 0.3275 | 0.4955 | 1.2439 | 1.3765 |
| Small-factor 70 | 10 + 60 | 9 | 0.0171 | 0.0362 | 0.0513 | 0.0607 |
| Target-factor 70 | 20 + 50 | 1 | 0.3547 | 0.5241 | 1.2713 | 1.4057 |
| Balanced 80 | 40 + 40 | 0 | 0.3864 | 0.5799 | 1.2683 | 1.4285 |
| Small-factor 80 | 10 + 70 | 9 | 0.0197 | 0.0326 | 0.0593 | 0.0611 |
| Target-factor 80 | 20 + 60 | 1 | 0.4142 | 0.6073 | 1.2822 | 1.4416 |

On this finite corpus, checked PRAC/int takes 1.42–2.12 times the int ladder's
median attempt time. The all-mpz ladder takes 3.00–5.22 times the int ladder;
all-mpz checked PRAC takes 3.10–5.55 times. These are complete ECM attempts,
including unsuccessful searches. They establish no PRAC or all-mpz promotion
for this implementation. They do not calibrate selective GMP helpers or
coarser compiled curve kernels.

The first balanced-40 primary block showed roughly twofold timing drift
shared across arms. Its repeated block supplies the reported values; both
blocks are retained. Seed-dependent early factors also cause legitimate
variation in attempt duration. No small timing differences are claimed as
wins. Nine fresh-process balanced-40 attempts per arm give cold medians of
0.3590/0.5469/1.3394/1.5164 seconds in table order, including startup, imports,
certificate validation, construction and the first attempt.

### Larger-bound extension: completed subset

At B1=11,000 / B2=1,873,422, all four arms split **9/27** distinct seeded
attempts across the three completed cases (3/9 per case). There are no
invalid factors, timeouts or backend-transition mismatches in these captures;
108 backend pairs agree across primary and repeated blocks.

| Input case | Factor digits | Successes / 9 | Int ladder | Int PRAC | gmpy2 ladder | gmpy2 PRAC |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Balanced 40 | 20 + 20 | 3 | 1.9302 | 2.5034 | 9.4742 | 10.0334 |
| Target-factor 50 | 20 + 30 | 3 | 1.8026 | 2.5088 | 9.5594 | 10.0408 |
| Target-factor 60 | 20 + 40 | 3 | 4.5589 | 6.1432 | 19.4741 | 20.7905 |

Values are median seconds in the final nine-seed block. Checked PRAC/int
costs 1.30–1.39 times the paired int ladder; the gmpy2 ladder costs
4.27–5.30 times, and gmpy2 PRAC costs 4.56–5.57 times. The first two cases
and the third were captured in separate, explicitly coordinated timing
windows. Absolute times across those windows are not a hardware-scaling
comparison. The arithmetic source hashes match; the only driver change was
an optional cold-control case selector, with all timed function ASTs unchanged.

To close the existing work, the optional extension stopped at the completed
60-digit case boundary. **The larger-bound 70-/80-digit cases and its cold
controls are unfinished.** No result is inferred for them. This limitation
does not remove any case from the completed 14-case current-bound study.
The evidence retains every completed sample, original/repeated blocks,
window provenance and the explicit closeout record.

The checked implementation is a correctness foundation, not the preferred
performance path. GMP-ECM precomputed Lucas codes and compact execution are
explicitly deferred to **P4.5/C6** in the roadmap; no such executor or chain
corpus is included in this change.

### Kernel and standalone stage-one diagnostics

`p41_prac.py` separately measures 28-, 57- and 96-digit moduli. Those sizes
are **not** the 40–80-digit ECM campaign corpus and support no claim about
whole factoring performance. They isolate construction, fixed Suyama stages,
kernels and exceptional recovery. At B1=2,000 the records reduce the abstract
6/5-weighted operation total from 31,369 to 25,273 (19.4%). This is an abstract
count, not a 19.4% runtime gain.

With a 512-record cache (303 prime powers at B1=2,000), the diagnostic medians
for complete three-curve stage-one runs are:

| Modulus digits | Ladder, 16-power chunks | Ladder, individual powers | Checked PRAC, individual powers |
| --- | ---: | ---: | ---: |
| 28 | 6.533 ms | 6.787 ms | 11.834 ms |
| 57 | 12.938 ms | 13.382 ms | 24.491 ms |
| 96 | 25.209 ms | 26.172 ms | 44.785 ms |

B1=2,000 chain construction is 9.206 ms; 16 exceptional recovered operations
are 18.521 µs as a separate tiny-field diagnostic. Fresh-process startup plus
one checked multiplication is 37.280 ms median. Actual point addition/doubling
kernel ratios are about 1.15–1.18 in this diagnostic, rather than assuming the
abstract 6/5 model predicts total Python execution cost. Full construction,
recovery and sample arrays remain in ignored local evidence.

```sh
pypy3 -B -u -m v2.benchmarks.p41_prac --output v2/benchmarks/results/p41-diagnostic-new.json
```

Research references and the explicit B3/C6 follow-ups are recorded in
[ROADMAP P4.1](../ROADMAP.md#p41--finish-prac-repair-and-precompute-valid-chains).
GMP-ECM's near-optimal Lucas generator and newer continued-fraction searches
are useful follow-ups; neither provides evidence that its chain-selection
savings outweigh checked Python execution here.

## P5.2-A3 reusable programs and campaign feasibility — 5 October 2026

The bounded A3 tranche adds opt-in immutable packed prime/power blocks, a
run-local finite retention cap and independent +/- coverage fixtures. Production
continuation remains unpaired; D tuning, wheel/common-Z experiments and automatic
allocation belong to B2/C2/C3. `ecm_program_bytes=0` retains defaults, streamed
work accounting and version-4 checkpoints. An enabled store reserves its full
cap, charges construction and reads, and rebuilds after a version-5 resume;
completed curves, buffered actions and RNG progress remain credited. Increasing
an existing campaign's bounds or curve count is not supported by this tranche.

The driver pins the pre-change `9b2d380` portfolio/stage sources in
`inputs/baselines/p52_a3_baseline.json` and a separate immutable helper snapshot
in `p52_a3_dependencies.json`. The original baseline bytes remain unchanged.
The private control freezes arithmetic, budgets, preprocessing, sieves and result
classes; inactive QS fallback types only support its dataclass annotations.
The independent corpus in
`inputs/corpora/p52_a3_corpus.json` contains ten freshly generated inputs with
Pocklington certificates, seed 20261005052 and algorithm seeds 7/19/41. Complete
factor results reconstruct and match the certified factors. Schedule-only rows
match independent integer-index Eratosthenes/prime-power counts and digests;
finite campaign rows validate every divisor and reconstruct the input.

Four arms use matched inputs, seeds, bounds and finite work/time/storage grants:
the frozen baseline, current default, an 8 MiB program cap, and packing with only
the scratch reserve (no retained blocks). Confirmation on Apple M4 / 24 GiB with
PyPy 7.3.23 implementing Python 3.11.15 supplies at least three seconds of validated
warmup and nine samples per arm. All confirmation arms meet the relative-IQR
threshold of 0.15. The following medians compare retained programs with the frozen
baseline; intervals are unpaired bootstrap 95% intervals conditional on these
fixed cohorts, not estimates for future inputs.

| Cohort, per sample | Streamed baseline | Retained programs | Time reduction, 95% interval |
| --- | ---: | ---: | --- |
| Four 10-digit inputs, three seeds, 12/12 complete | 1.767 ms | 1.859 ms | -5.3%, -7.2 to -2.1% |
| Four 16-digit inputs, three seeds, 12/12 complete | 11.036 ms | 11.407 ms | -3.4%, -11.4 to 13.8% |
| One 80-digit input, two starts, three curves each; no splits | 812.181 ms | 744.879 ms | 8.3%, 7.7 to 9.1% |
| 11,000/1,900,000 schedule only, three passes | 160.669 ms | 121.047 ms | 24.7%, 22.4 to 26.8% |
| 50,000/5,000,000 schedule only, three passes | 455.596 ms | 324.290 ms | 28.8%, 27.1 to 32.8% |

The large campaign performs 54 curve attempts per arm across the nine samples,
with no factors. Its saving is a finite-curve execution result, not faster
factoring or improved success. Small complete factoring regresses; the medium
interval crosses zero. Packing without retention is about 4% slower on the large
campaign. Default streamed outcomes and work are unchanged; tiny default timing
differences earn no speed claim. The earlier 2,000/147,396 schedule comparison
remains unstable after eight seconds of warmup and 63 samples and is inconclusive.
Retain the default and bounds; the opt-in provides an explicit reuse control for
B2, not a general promotion.

A deterministic nonsplitting probe uses one 266-bit input, a 329-bit admission
envelope and three curves under 50,000,000 work units and 120-second wall/CPU caps.
The table includes full point/recovery work. Construction increases first-curve
work; later curves amortize it. Counts exclude one-time context setup.

| B1/B2 | Streamed units per curve | Program first / subsequent curve | Program accounting after three curves | Owned reserve, streamed / 8 MiB cap |
| --- | ---: | ---: | ---: | ---: |
| 2,000/147,396 | 110,517 | 124,445 / 45,014 | 414,528 bytes | 1,319,864 / 9,708,472 bytes |
| 11,000/1,900,000 | 1,458,993 | 1,602,357 / 446,275 | 1,888,800 bytes | 2,665,104 / 11,053,712 bytes |
| 50,000/5,000,000 | 4,080,627 | 4,434,273 / 1,122,672 | 4,345,712 bytes | 3,825,120 / 12,213,728 bytes |

Program accounting includes the fixed scratch reserve and per-block allowance;
the admission reserve includes the entire requested cap, not just populated
blocks. A 16 MiB workspace and 329-bit envelope admit these tiers with an 8 MiB
program cap. The 2,000,000-unit default cannot finish the tested full nonsplitting
50,000/5,000,000 curve. These counts are not universal minimums for finding a
factor. A predeclared 10,000-curve campaign is finite and admissible; a smaller
work grant pauses it and can be extended cumulatively under the identical config.
It does not guarantee completion or permit adding curves to an exhausted config.

Observed large-schedule worker RSS is 102.02 MiB baseline and 117.16 MiB with
programs. RSS includes the interpreter/JIT, independent oracle and warmup and is
separate from owned workspace. Nine separate cold worker starts per arm give
108.550/102.611/104.374/103.487 ms for baseline/default/programs/regenerated;
these include imports and corpus validation and establish no cold promotion.

Local ignored evidence lives in `results/p52/confirmation.json` and the earlier
`final-comparison.json`, with source hashes and prior noisy attempts preserved.
The final smallest-buffer checkpoint repair follows those measurements: an AST
comparison isolates it to `_verify_progress`; fresh execution and measurement
functions are unchanged. The repaired source passes the full acceptance
suite. The later private
control loader freezes the byte-identical helper sources to preserve the old
control after mainline arithmetic changes; future captures pin that newer
loader. Only the later focused 60/80-digit pass was collected in an exclusive
window; historical timing/cold source pins have not been rewritten.

```sh
mkdir -p v2/benchmarks/results/p52
pypy3 -B -m v2.benchmarks.p52_a3 --cases small medium large_campaign middle_schedule large_schedule --cold --probes --output v2/benchmarks/results/p52/confirmation-new.json
pypy3 -B -m v2.benchmarks.p52_a3 --cold --probes --output v2/benchmarks/results/p52/all-cases-new.json
```

Output paths refuse overwrites. Unstable arms extend to five seconds/31 samples,
then eight seconds/63 samples; unresolved instability prevents a timing claim.
Acceptance covers exact actions, independent point cross-products, coverage
exceptions/tails, atomic refusal, capped regeneration, mixed-factor replay,
quiet factoring, every resume phase, legacy schema, resealed bounds/endpoints,
wide packed words and cancellation with the smallest prime buffer. All 323
PyPy tests and full lint pass, including private-control isolation, independent
workload proofs and owned-worker overlap handling. Committed-files-only
acceptance is recorded with the branch.

## P5.2 realistic workload exploration and prepared protocol — 5 October 2026

The user ended the expanded study early, then requested a focused 60/80-digit
verification. The 256-curve campaigns, extra-seed sweep, larger-bound probes and
full twelve-input exclusive-window confirmation were prepared but **not run**.
Their commands below describe future experiments, not accepted results. The
focused fixed-curve comparison is accepted below; no population-wide factor-size
allocation or large-factor success claim is established by this tranche.

The initial middle-policy capture did complete nine validated samples per
arm after more than three seconds of warmup, on the twelve certified inputs
with seeds 7/19. Other benchmark jobs overlapped it, so all its timings remain
diagnostic. Its deterministic completion records remain useful: no input was
time-censored, every result reconstructed, and all three arms had identical
curve assignments and outcomes under the eight-curve 11,000/1,900,000 tier
and 50,000,000-unit grant. Each row counts unique input/seed starts; the nine
timing repetitions are not additional independent factoring trials.

| Smaller factor | Completed starts | Scope |
| --- | ---: | --- |
| 10 digits | 8/8 | One input in each 30/40/60/80-digit band, two seeds |
| 15 digits | 4/4 | Two 30-digit semiprimes, two seeds |
| 20 digits | 1/8 | Two 40-digit, one 60-digit and one 80-digit input |
| 30 digits | 0/2 | One balanced 60-digit semiprime, two seeds |
| 40 digits | 0/2 | One balanced 80-digit semiprime, two seeds |

These fixed-case counts do not estimate a population success probability.
In particular, eight unsuccessful curves do not show that a factor cannot
be found. The original A3 timing trend concerns repeated schedule reuse; its differing
bounds and campaign lengths do not isolate input bit length.

`p52_realistic` measures the full ECM-only portfolio on certified 30-, 40-, 60-
and 80-digit inputs, including preprocessing, classification, curve execution
and recursive split validation. Each band has a 10-digit-factor case, a 15/20-
digit target-factor case and a balanced semiprime. The target and balanced cells
at 30/40 digits use different inputs. The new frozen corpus uses generation seed
20261006053 and independent Pocklington proofs. Known factors are used only for
validation; terminal certainty labels and every unresolved cofactor are retained.

The two original algorithm seeds are 7/19. Warmed comparisons use at least three
seconds of validated warmup and nine samples, extending unstable cohort timing
to five seconds/31 samples and then eight seconds/63 samples. Cell timing is
separately labeled stable or inconclusive. These fixed inputs support scoped
completion/cost evidence, not a population-wide allocation policy.

| Policy | B1/B2 | Predeclared curves | Total work grant per input/seed |
| --- | --- | ---: | ---: |
| Pretest | 2,000/147,396 | 32 | 2,000,000 |
| Middle | 11,000/1,900,000 | 8 | 50,000,000 |
| Large | 50,000/5,000,000 | 4 | 50,000,000 |
| Deep 20-digit-factor comparison | 11,000/1,900,000 | 256 | 400,000,000 |

All use 120-second wall/CPU caps, a 329-bit envelope and 16 MiB owned workspace;
programs reserve 8 MiB. Trial bound is 5, rho/p−1 are disabled and there is no
sieve fallback. The pretest compares identical total allowances; changing the
work ledger can fund more curves and produce different terminal outcomes. The
campaigns compare identical declared curves, assignments and arithmetic outcomes.
Elapsed unresolved searches are search costs, not time to factor.

The deep comparison selects the four 20-digit-factor inputs (two 40-digit
semiprimes and the 60-/80-digit unbalanced cases). A separate completion-only
sweep adds algorithm seeds 41/73/101/137/179/223/269 without timing claims; combined
with the original two seeds it examines nine starts on each fixed input. The
native [GMP-ECM parameter guidance](https://github.com/sethtroisi/gmp-ecm/blob/main/README)
motivates deeper allowances, including roughly 74 expected curves for a 20-digit
factor at 11,000/1,900,000. This is a hypothesis for the PyPy engine, not a measured
success guarantee or copied default.

`p52_wider` separately probes balanced 60-/80-digit cases at
250,000/130,000,000 (four curves, one billion work units) and
3,000,000/5,700,000,000 (one curve, 50 billion units). Both reserve 512 MiB owned
workspace, with a 128 MiB program cap in the opt-in arm, and stop at 120 seconds
of wall or CPU use. Each is one cold seeded feasibility run: stage progress,
completed curves and the actual exhaustion reason are recorded, without a
comparative speed claim. These candidate bounds also come from the native table.

The initial realistic capture overlapped p4.3's study and was stopped. Its six
completed arm captures and source hashes remain under ignored `results/p52/`,
explicitly diagnostic for timing. `--check-quiet` checks for other benchmark/test
interpreter processes before, during and after each worker, terminating only the
owned worker on overlap. It supplements coordinated timing windows; it does not
measure every background operating-system activity. New output paths, including
per-arm partials, refuse overwrites.


### Focused 60/80-digit verification

The later quiet pass uses the middle tier on the balanced 60-digit (198-bit)
and 80-digit (263-bit) fixtures, with the same 11,000/1,900,000 bounds, two seeded
starts and eight predeclared curves per start. All three arms complete the same
32 curves per sample, with no factors, no time censoring and identical curve
assignments/outcomes. Default streamed work also matches the frozen control.
Each row below sums both starts, or 16 curves, for that input. Nine repetitions
per arm reuse these fixed starts; they are not independent factor-success trials.

Validated warmup is 8.106/8.052/6.973 seconds for baseline/default/programs,
respectively. All cohort and cell relative IQRs are below 0.05, with no extended
attempt needed. The monitor finds no competing benchmark/test interpreter
before, during or after workers. These are uninstrumented full-portfolio costs,
including preprocessing/classification and result validation, not profile data.

| Input | Frozen streamed | Current default | Retained programs | Program saving, conditional 95% interval |
| --- | ---: | ---: | ---: | --- |
| Balanced 60 digits, 16 curves | 3.567 s | 3.652 s | 3.007 s | 15.7%, 13.0 to 17.4% |
| Balanced 80 digits, 16 curves | 4.089 s | 4.205 s | 3.554 s | 13.1%, 11.5 to 14.3% |

Across both inputs, programs save 14.1% (12.3–16.2%). Current default medians
are 2.4%/2.8% slower by size; both conditional intervals cross zero, as does the
combined default interval (-5.4 to 1.4% reduction). No default speed claim is
made. Baseline uses 11,683,254/11,689,655 work units per start; programs use
4,737,592/4,743,993. Different work ledgers do not by themselves establish speed.
The same 16 MiB owned allowance and 8 MiB program cap apply throughout.

The observed relative saving falls by 2.63 percentage points at 80 digits.
A 10,000-draw bootstrap with seed 52080 resamples cohort indices jointly within
each arm, preserving correlation between its 60/80 timings, and independently
across arms. Its conditional difference interval is -4.52 to -0.09 points
(80 minus 60). The earlier independent-band calculation is preserved locally,
but the paired-cohort calculation is used for this difference. The per-size
intervals use the driver's 3,000-draw unpaired-arm bootstrap. Absolute savings
are similar, 0.560/0.534 seconds per 16 curves: the modest relative dilution is
consistent with arithmetic taking a larger fraction at bigger integers. Two
fixed inputs and two starts do not establish a general monotonic size trend.

Keep programs opt-in. This verifies useful finite-curve execution at 60/80 digits;
it does not show faster successful factoring of balanced 30/40-digit factors.
Deeper success calibration and huge-bound feasibility remain deferred. Raw local
captures are `results/p52/focused-60-80-quiet.json` with source/corpus/control
hashes, and `focused-60-80-paired-trend.json`. The source passes 323 tests, full
lint, both proof loaders, short validation-only smoke cases and all 50 benchmark
module imports in a committed-files-only archive, without generated evidence.

```sh
pypy3 -B -m v2.benchmarks.p52_realistic --check-quiet --policies middle --fixtures 60d_balanced 80d_balanced --output v2/benchmarks/results/p52/focused-60-80-new.json
```

The following commands describe the deferred broader protocol:

```sh
mkdir -p v2/benchmarks/results/p52
pypy3 -B -m v2.benchmarks.p52_realistic --check-quiet --output v2/benchmarks/results/p52/realistic-quiet-new.json
pypy3 -B -m v2.benchmarks.p52_realistic --check-quiet --policies deep --fixtures 40d_target 40d_balanced 60d_target 80d_target --output v2/benchmarks/results/p52/deep-quiet-new.json
pypy3 -B -m v2.benchmarks.p52_realistic --check-quiet --coverage-only --policies deep --fixtures 40d_target 40d_balanced 60d_target 80d_target --seeds 41 73 101 137 179 223 269 --output v2/benchmarks/results/p52/deep-seeds-new.json
pypy3 -B -m v2.benchmarks.p52_wider --check-quiet --output v2/benchmarks/results/p52/wider-quiet-new.json
```

### A3/A4 integration acceptance — 5 October 2026

The A3 merge preserves the accepted P4.1/A4 implementation and all six
struck-through columns in both completed execution-plan rows. A committed-files-
only merge candidate passes 337 tests under system PyPy, full lint, all 54
benchmark-module imports, both A3 proof/control loaders and short validation-only
smoke cases. System PyPy lacks optional gmpy2; all three otherwise-skipped GMP
campaign tests pass in the existing PyPy 3.11 venv where it is installed.

A3 source files are byte-identical to the measured branch. Compared with that
branch, mainline's common ECM functions have identical ASTs except for A4's
intentional `multiply_prac` implementation; the production ladder is unchanged.
Integration adds no timing claim and leaves historical source pins intact.
All 47 local A3 evidence files are copied and checksum-verified under ignored
`results/p52/` before managed-worktree removal. Combined check logs and the
preservation manifest stay local there as well.

```sh
make -C v2 test
make -C v2 lint
v2/.venv/bin/python -B -m unittest v2.tests.test_prac_campaign.GmpCampaignTests -v
```
