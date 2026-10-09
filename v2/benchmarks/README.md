# Factor benchmarks

Validated factoring experiments on **PyPy implementing Python 3.11**. This guide
keeps the stage history, accepted changes and rejected experiments concise.
The [v2 guide](../README.md) covers usage; the [roadmap](../ROADMAP.md) records
remaining acceptance gates.

## QS/GNFS research reconciliation (9 October 2026)

The [source-linked comparison](qs_gnfs_research.md) pins 11 repositories and
records license/attribution requirements, primary papers and independently
checked counterexamples. It routes residual handling to coordinated C1/E3,
DLP to C1/R4/P5.4, Four Russians to A8/B6, changed-workload sieve allocation
to E6, refreshed SSS to A7/R5 and polynomial/batch/general-incidence variants
to F5/F7/F6. GNFS remains A9 onward after H1's coverage review; F1 owns a
measured crossover. These are prospective experiments, not measured PyPy wins;
no new benchmark or implementation gate is closed by source research.

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
implementation was adapted; license obligations are identified above, and
future source copying requires checking and retaining the full source notices.

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
complete confirmation was 25.10% / 84.37%. Every sample was retained during
acceptance. These
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
Raw captures, source hashes, selection records and profiles were kept in the
isolated worktree's ignored `results/a10/` and `audit/a10/` directories during
acceptance. They were removed with that worktree at the user's request after
integration. Required controls, certified inputs, protocols, runners, source
citations and the acceptance summary remain versioned.

Validation covers strict endpoints/neighbors, prime powers, Carmichael and
strong pseudoprimes, exact reconstruction/labels/unresolved cofactors,
explicit round/RNG spies, finite work/time/cancellation and checked resume
under legacy schemas 4/5/6 on int/GMP. `make -C v2 test` passed 369 tests
(two optional-GMP skips), the GMP-enabled repeat passed 372, and lint passed.
One earlier GMP suite failed the unchanged QS snapshot-lifetime assertion.
B1 subsequently traced it to live PyPy JIT roots and isolated the ownership
assertion in a finite JIT-off child; mainline retains that accepted test fix.
No QS implementation change is included in A10.
A committed-files-only archive of `0ee86ff` also passed the 369/372-test
suites, all 59 benchmark imports and the corpus/control/protocol/v1-adapter
loaders. The final acceptance update changes documentation only.

**Mainline composition:** candidate `6bf7ca1` integrates settled B1/B2
`6ec8a01` with A10 and the roadmap research reconciliation. The combined
PyPy/GMP suite passes 401 tests and lint; a committed-files-only archive passes
398 system-PyPy tests (three optional-GMP skips), all 67 benchmark imports and
the A10 proof/control/protocol/v1-adapter loaders. A new integration regression
preserves fresh/legacy primality policies and work/labels in B2 schemas 7/8
on int/GMP. No new timings were run; historical A10/B1/B2 performance pins
remain unchanged and E1/H1 confirmation stays open.

## B1 joint QS/MPQS/SIQS calibration — 9 October 2026

The B1 study uses committed mainline `94caf40`, including integrated R2,
in an isolated worktree. Runtime arithmetic and production defaults are
unchanged. `inputs/controls/b1_frozen.json` pins runtime/driver hashes, both
previously inspected training corpora, seeds 7/29, 23 joint bundles, finite
allowances and selection criteria. `b1_40d_frozen.json` separately freezes a
bounded eight-bundle follow-up after the successful 40-digit feasibility
probe; it does not rewrite the original protocol. Selection files freeze
each mode before its new certified corpus is generated.

The 30-digit sweep uses three balanced inputs and both seeds; the 40-digit
follow-up uses one balanced training input and both seeds. Bound/interval,
nearest/flyer selection, A-factor count, Gray quotas, residual bounds and
row/partial/atom/matrix allowances vary together. This is a finite bundle
comparison, not an exhaustive parameter search. QS keeps its fixed A=1;
MPQS uses the accepted external-square representation. Powers/bucket
collection, exact verification, certainty labels and checked resume stay
intact. Known factors enter validation only. Pocklington generation biases
the corpus toward primes with a large known p-1 factor; no RSA-distribution
or population success claim is made.

Every call has a 10^13 work grant and finite wall/CPU limits. The common
owned-memory envelope is 256 MiB; individual bundles can reserve less.
Feasibility caps are 5 seconds at 30 digits, 30 at 40/60, and 15 at 70–99.
Warmed 30-digit timing admits bundles completing every screened start in
at most two seconds; the separately frozen 40-digit threshold is ten.
Each timed arm receives at least three seconds of validated PyPy warmup
and nine samples, extending unstable series to fifteen, up to three blocks.
Fresh comparisons interleave arms. Process inventory checks surround calls
and poll long calls; the same exclusive flock and owner record are shared
with B2/A10. These checks supplement coordination, not an OS-wide isolation
guarantee. Cold startup and instrumented profiles supply no performance claim.

The machine is Apple M4 / 24 GiB, using PyPy 7.3.23 implementing Python
3.11.15. The following fresh 30-digit times sum both seeded complete calls,
including setup, collection, filtering, extraction and classification.
Process-inventory checks outside those calls are excluded; their separate
outer timers remain in the captures. Each row is one independent input, not
18 independent factoring trials. Close-factor timing needed a second block
with five seconds of warmup and fifteen samples; other primary rows use nine.

| Fresh class | R2 control, seconds | Selected SIQS | Selected MPQS |
| --- | ---: | ---: | ---: |
| Balanced | 0.447480 | 0.283728 | 0.415146 |
| Uneven, five-digit smaller factor | 0.082175 | 0.429751 | 0.517652 |
| Uneven, ten-digit smaller factor | 0.474256 | 0.363277 | 0.426907 |
| p-1 smooth control | 0.504901 | 0.455407 | 0.576002 |
| p+1 smooth control | 0.694064 | 0.771634 | 1.115285 |
| Close factors | 0.462278 | 0.358044 | 0.418211 |
| Prime square | 0.051126 | 0.024225 | 0.092133 |

The balanced SIQS reduction is 36.6%, with a paired-repeat 95% interval of
35.4–37.5%; MPQS reduces it by 7.2% (1.6–8.8%). Intervals are conditional on
the fixed input/seeds. All R2/SIQS/MPQS classes complete both unique starts.
The five-digit-factor control's factoring timer initially failed stability
although its outer timer passed. Its reported cell comes from a separate
frozen-parameter repeat: over three seconds of actual factoring warmup and
nine stable samples per arm. The original unstable cell is excluded.

The selected 30-digit SIQS bundle uses base bound 3,000, half-width 8,192,
four A factors, flyer selection, eight effective Gray polynomials, residual
bound 9,000,000, 2,048 rows/partials and 8,192 atoms. MPQS uses base 3,000,
half-width 32,768 and the same stores/residual bound, with external square A
and one polynomial per coefficient. Both have a finite 100,000-family grant,
1 MiB checkpoint cap and 256 MiB owned cap. Actual assignment-space exhaustion
remains distinct from the configured family grant.

At 40 digits, six of eight training bundles complete both starts and pass
nine-sample stability; both fixed-QS bundles exhaust their windows. Selection
retains SIQS base 10,000, half-width 65,536, five flyer-selected A factors and
16 effective Gray polynomials. External-square MPQS selects the same base
and width. Both reserve 8,192 rows/partials, 32,768 atoms, a 100,000,000
residual bound, 256 MiB owned memory and 1 MiB checkpoints. On the fresh
balanced input, both complete 2/2 unique starts and all 18 timed calls.
The two-seed cohort medians are 15.086435 seconds for SIQS and 14.043059
for MPQS; MPQS's paired reduction is 6.9% (95% interval 5.4–9.5%). This
40-digit comparison uses the frozen feasible SIQS control, not the smaller
30-digit R2 configuration. Neither the point estimate nor interval reaches
the 10% timing gate, so it establishes no MPQS crossover or promotion.

A separately frozen wider-QS follow-up tests bases 3,000/10,000/30,000 and
half-widths 131,072/499,999, using 8,192 rows/partials and at most 32,768 atoms.
At the widest interval, training completion is 0/6, 2/6 and 4/6, respectively.
The frozen 30,000-base challenger fails both fresh starts: 1,072 verified rows
leave zero rows after singleton filtering. This is finite-window/useful-yield
exhaustion, not a matrix-memory failure or a time-to-factor result. The same
second fresh balanced input confirms the unchanged SIQS challenger at
0.496728 seconds versus 0.660333 for R2; no retuning used either fresh input.

The upper probes give SIQS/MPQS a 30-second cap at 60 digits and 15 seconds
at 70–99. At 60 digits each produces 47 verified rows, all removed by
singleton filtering, with owned peaks 81.49/66.62 MiB. SIQS uses about 90%
of its time in collection. At 70 digits SIQS has one row and MPQS zero;
80–99 have zero. A products approach the exact targets, the necessary matrix
reservation fits, and neither store nor matrix capacity stops these attempts.
Use post-filter output rows/columns to assess useful yield, not admitted-row
counts alone. These are censored feasibility probes; they establish no
upper-band time-to-factor ratio or impossibility claim.

**Decisions:** adopt the frozen protocols, independent confirmation inputs
and explicit balanced-input presets as the bounded B1 result. Retain runtime
defaults and exact powers/bucket collection. The 30-digit balanced SIQS win
does not establish a general portfolio policy: five-digit-factor and p+1
controls regress, populations are small, and no ECM-to-SIQS handoff is timed.
Defer larger calibration and combined E1/C3/G1/H1 confirmation. DLP remains
deferred pending an affordable recoverable-residual population, beyond these
low-yield probes. No calibrated matrix-capacity failure or eligible-base-prime
CRT cost triggers matrix redesign or C8; GNFS stays with its existing roadmap.
The small-store bundle fits 32 MiB but loses partials through FIFO eviction;
larger stores help its wide counterpart without beating the narrow selection.
These joint comparisons do not isolate a causal effect for every parameter.
The incomplete-QS tie-break uses admitted verified rows, although the freeze
describes them as useful rows; post-filter counts expose that limitation.
No incomplete QS selection supplies an accepted speed or promotion claim.

All three new 58-fixture certified corpora are mutually disjoint and disjoint
from training. Early checked resumes reject changed configuration and reduced
work grants. Eight distinct selected/control configurations also pass deeper
pause/restore and a second checkpoint, preserving cumulative resources and
reconstruction within the frozen 1 MiB checkpoint caps. Local raw JSON,
stdout and verification records remain in this worktree's ignored
`v2/benchmarks/results/b1/`. No cold-start or profile timing is pooled with
the warmed evidence. New integrated A10/B2 sources require a new confirmation;
these captures remain pinned to `94caf40` runtime arithmetic.
Required checks pass: `make -C v2 test PYTHON=.venv/bin/python` runs 365 tests
on PyPy/GMP, and `make -C v2 lint` passes. The v2-local changelog preserves
this task's directory boundary.

Native implementations inform the candidate design, not Python defaults:
[PARI's flyer selection](https://pari.math.u-bordeaux.fr/lcov-report/basemath/mpqs.c.gcov.html)
(development `31092-e6893b0017`),
[msieve's A-factor/reuse tradeoff](https://github.com/radii/msieve/blob/master/mpqs/poly.c),
and [Zimmermann's MPQS notes](https://members.loria.fr/PZimmermann/talks/tiny-mpqs.pdf)
relate coefficient quality, interval and reuse. Sources were read on
9 October 2026; no external engine was executed or timed.

Reproduce in a quiet machine window from this committed tree. Freeze and
selection outputs refuse replacement; use new paths for a new study.

```sh
pypy3 -m v2.benchmarks.b1_calibration freeze --protocol NEW_PROTOCOL.json
pypy3 -m v2.benchmarks.b1_calibration probe --protocol NEW_PROTOCOL.json \
  --output v2/benchmarks/results/b1/new-probes.json
pypy3 -m v2.benchmarks.b1_calibration train --protocol NEW_PROTOCOL.json \
  --output v2/benchmarks/results/b1/new-training.json
pypy3 -m v2.benchmarks.b1_calibration select --protocol NEW_PROTOCOL.json \
  --training v2/benchmarks/results/b1/new-training.json \
  --probes v2/benchmarks/results/b1/new-probes.json --output NEW_SELECTED.json
pypy3 -m v2.benchmarks.b1_calibration confirm --protocol NEW_PROTOCOL.json \
  --selected NEW_SELECTED.json --corpus NEW_CORPUS.json \
  --output v2/benchmarks/results/b1/new-confirmation.json
pypy3 -m v2.benchmarks.b1_40d freeze --protocol NEW_40D_PROTOCOL.json \
  --parent NEW_PROTOCOL.json
pypy3 -m v2.benchmarks.b1_40d run --protocol NEW_40D_PROTOCOL.json \
  --selected NEW_40D_SELECTED.json --corpus NEW_40D_CORPUS.json \
  --output v2/benchmarks/results/b1/new-40d.json
pypy3 -m v2.benchmarks.b1_qs_width freeze --protocol NEW_QS_PROTOCOL.json \
  --parent NEW_PROTOCOL.json
pypy3 -m v2.benchmarks.b1_qs_width run --protocol NEW_QS_PROTOCOL.json \
  --selected NEW_QS_SELECTED.json --corpus NEW_QS_CORPUS.json \
  --siqs-selected NEW_SELECTED.json \
  --output v2/benchmarks/results/b1/new-qs-width.json
pypy3 -m v2.benchmarks.b1_timer_confirmation --selected NEW_SELECTED.json \
  --corpus NEW_CORPUS.json --kind uneven_5 \
  --output v2/benchmarks/results/b1/new-uneven-timers.json
```

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

## P5.2-B2 aligned-wheel follow-up — 9 October 2026

`p52_b2_wheel.py` freezes the original paired engine at `a670c4d`, retains
the accepted `94caf40` streamed/reusable controls, and compares an opt-in
aligned-cell wheel. The existing B2 corpus has already been inspected, so
`build_p52_b2_wheel_inputs.py` generates fresh certified training and held-out
inputs with seed 2026100953. Required snapshots, corpus and protocol live in
versioned `inputs/`; raw captures remain local in `results/p52-b2-wheel/`.

The frozen workload sizes, tiers, algorithm seeds 7/19/41, 50-million-unit
work limit, 120-second wall/CPU limits, 16 MiB workspace, 8 MiB program cap and
1,024-slot segments match the original study. Wheel grids are 6/30/210 for
small, 30/210/420 for medium, 210/420/840 for uneven and 210/840/1890 for the
nonsplitting campaign. Structured inputs inherit the small selection. Old
paired controls use the previously selected D=24/64/768/2048. A regeneration
arm uses the middle wheel with only the minimum program scratch reservation.
Select the fastest stable training wheel with no completion regression against
either accepted control; freeze that choice before held-out timing. This is a
small fixed-cohort comparison, not factor-size or curve-allocation calibration.

At least three seconds of validated warmup and nine samples are required;
unstable arms extend to 5 seconds/31 samples, then 8 seconds/63 samples.
The relative-IQR threshold is 0.15, and wall/CPU-censored or nondeterministic
captures cannot pass. The same whole-factoring promotion rule applies: at
least 10% median time reduction with an interval above zero, or a 10-point
completion gain with an interval above zero, without a >5-point completion
regression elsewhere. Fixed nonsplitting campaigns cannot promote defaults.
All setup, tables, decoding, recovery, recursion and checkpoint costs are in
the end-to-end measurement. Cold starts and instrumented diagnostics are
separate; private-control reconstruction adds to cold startup and cannot be
interpreted as an algorithmic speed advantage.

Execution acquires the machine-wide flock and checks for competing benchmark/
test processes. Coordinate the window with B1 and other heavy work first.
The runner saves every attempt, pins source and input hashes, and refuses to
overwrite evidence. The original B2 controls/protocol/results are unchanged.

```sh
pypy3 -B -m v2.benchmarks.p52_b2_wheel --phase training --diagnostics --output v2/benchmarks/results/p52-b2-wheel/training.json
pypy3 -B -m v2.benchmarks.p52_b2_wheel --phase held_out --selection v2/benchmarks/results/p52-b2-wheel/training.json --cold --output v2/benchmarks/results/p52-b2-wheel/held-out.json
```

The implementation uses complete nearest-center cells, not Prime95's extended
distance/relocation matcher. Its sparse odd-multiple generator still pays for
discarded intermediate points. Those distinctions prevent attributing native
pairing percentages to this PyPy implementation. Substantial matching,
relocation, compact-map and common-Z extensions are deferred in C2; polynomial
continuation remains F3.

**Completed result: retain every production default.** Candidate `dba3cd2`
ran under the coordinated flock after B1's QA finished; A10 confirmed no heavy
work. No competing benchmark/test process was detected. Ordinary desktop
activity remained; the captures do not assert an otherwise idle operating
system. PyPy 7.3.23 implements Python 3.11.15. All 35 final training captures
and 20 final held-out captures satisfy the frozen stability rule. Training
needed 42 attempts: five arms finished at 31 samples and one at 63. Held-out
needed 25 attempts: medium wheel/programs/legacy finished at 31 and uneven
wheel at 63. All earlier attempts are retained. The uneven wheel's relative
IQR fell from 0.274 (9 samples) through 0.195 (31) to 0.141 (63); its final
median is 0.360 seconds. This extended fixed-cohort estimate is not evidence
of stability across all populations or machines.

Training selected W=30/210/840/1890 for small/medium/uneven/campaign, with
structured inputs inheriting 30. Held-out medians below measure complete
cohorts, including all declared seeds and all factoring/setup/serialization
costs. Each small, medium and uneven arm completes 9/9 starts; structured
completes 3/3. The campaign completes all four declared curves per start but
finds no factors (0/3 starts complete). No arm changes completion.

| Held-out class | Streamed seconds | Programs seconds | Original paired seconds | Wheel seconds | Wheel vs programs | Wheel vs original paired |
| --- | ---: | ---: | ---: | ---: | --- | --- |
| Small | 0.001245 | 0.001324 | 0.002392 | 0.002319 | 75.1% slower | 3.1% faster |
| Medium | 0.005141 | 0.005301 | 0.016067 | 0.011064 | 108.7% slower | 31.1% faster |
| Uneven | 0.166597 | 0.165668 | 0.249830 | 0.359836 | 117.2% slower | 44.0% slower |
| Fixed nonsplitting campaign | 1.588467 | 1.417380 | 1.840164 | 1.752086 | 23.6% slower | 4.8% faster |
| Structured | 0.000147 | 0.000146 | 0.000153 | 0.000151 | 3.6% slower | 1.0% faster, inconclusive |

All selected wheel arms also lose to streamed execution. Conditional 95%
bootstrap intervals for the wheel's time reduction versus programs are
[-79.9,-72.6]%, [-110.1,-107.6]%, [-133.7,-110.9]%, [-24.0,-22.4]% and
[-4.3,-1.7]% in table order. Versus original pairing, the medium interval is
[30.6,31.4]%, campaign [4.5,5.4]% and uneven [-54.9,-33.3]%. These resample
repeated timings of fixed cohorts; they are not independent input-population
confidence intervals. No speed or completion promotion gate passes. Reusable
unpaired programs retain value: this fresh fixed campaign saves 10.8% versus
streamed, consistent with the earlier 10.5% observation. Neither campaign
establishes a successful-factoring or global allocation advantage.

Separate three-curve campaign diagnostics explain the bounded layout change:

| Arm | Certified eligible primes | Product terms | Giant advances | Retained baby point slots | First/later curve work |
| --- | ---: | ---: | ---: | ---: | --- |
| Programs | 422,082 by unpaired schedule | 422,082 | 2,055 | 1,379 table slots | 1,602,356 / 446,274 |
| Original D=2048 | 422,082 | 422,082 | 1,383 | 1,024 plus sentinel | 2,024,552 / 587,082 |
| W=210 | 422,082 | 353,442 | 26,988 | 24 | 2,065,311 / 532,130 |
| W=840 | 422,082 | 353,064 | 6,747 | 96 | 2,205,886 / 518,573 |
| W=1890 | 422,082 | 352,884 | 2,997 | 216 | 2,049,600 / 516,419 |

Selected W=1890 eliminates 16.4% of products, representing 32.8% of primes
covered in pairs. Exact prime coverage is unchanged. It retains 216 rather
than 1,024 baby points, but still computes discarded intermediate multiples
and incurs more giant advances. Its final program store is 6,204,320 bytes
versus original pairing's 6,863,584 and programs' 1,888,800. Conservative owned
reserves are 12,845,944 / 12,948,696 / 11,053,712 bytes respectively. Held-out
campaign process peaks are 88.3 / 97.5 / 87.6 MiB; RSS also includes runtime,
JIT and private control loading, so it is not table storage. Regeneration-only
W=840 costs 2.746 seconds in training versus 1.843 with retention.

All 180 separate cold starts validate identical rows/work to their warm arm.
Campaign process-wall medians are 1.923 / 1.772 / 2.233 / 2.080 seconds for
streamed / programs / original paired / wheel. Small cold apparent wins cannot
promote the wheel: immutable controls pay private source reconstruction, while
the candidate imports normally. Separate validated profiles cover three
complete starts (12 campaign curves), not accepted timing ratios. The wheel
profile records 12,000 coverage reads, 3,000 coverage compilations and 1,411,536
sparse-slot binary searches. Packing/decoding and slot lookup remain concrete
costs despite fewer products. Compact indexed maps and richer matching require
their own C2 representation, recovery, memory and timing gates; this tranche
does not bundle that tuning.

Full `make -C v2 test` and `make -C v2 lint` pass. A committed-files-only archive
passes 376 system-PyPy tests (two optional skips), 379 PyPy/GMP tests (no skips),
and all 61 benchmark imports plus frozen proof/control loading. New cases
cover exact independent coverage across segment boundaries, affine points,
zero-predecessor initialization, mixed-factor and one-sided saturation, every
action's resumed execution, cancelled/refused work, regeneration, corrupted
metadata, canonical int/GMP campaigns and cumulative budget extension. No
v1 or SIQS implementation file changes.

Raw evidence is local under `results/p52-b2-wheel/`, including every attempt,
cold sample, work/storage diagnostic, profile and QA log, with `manifest.json`.
Training SHA-256 is
`63becb5e078b112ee44e958a4f50bd8b6ec4ee5bf42c26048555135043ac93cf`;
held-out SHA-256 is
`0c4485ab36648626d3abd5b1634590c0cf09a1255a3dd6125fbc03aa190c1ff0`.
Frozen legacy source / protocol / corpus hashes are respectively
`59b44b7c13ca8be758f062337a76c37b2fcecddce7da0b23bf9a4389a38b054f`,
`004c7e0c04c61b5562b57c413f81ac8f44d9efd1fedc304f71d24eae49a8337a`,
`b876defa8ad640219bdeb8d97f8d66a581fa3c110d1fc683309d80a576c9e0ab`.

## P5.2-B2 paired continuation — 9 October 2026

The bounded implementation and matched experiments are complete. Pairing stays
opt-in: all selected held-out configurations lose to both accepted controls,
with unchanged completion. No production promotion gate passes. Checks pass
367 system-PyPy tests (two optional skips), 370 PyPy/GMP tests without skips
and full `make -C v2 lint`.
Native integers, streamed execution, the ladder and production
parameters remain defaults. Required inputs are versioned; raw evidence stays
in ignored `results/` directories.

`p52_b2.py` compares immutable committed-mainline (`94caf40`) streamed and
reusable-program controls with three paired D choices and a regeneration-only
control. The private snapshot `inputs/baselines/p52_b2_mainline.json` freezes
all active ECM arithmetic, schedules, budgets and portfolio dependencies.
The older A3 snapshots are unchanged. Inactive QS types provide annotations
only; these runs cannot dispatch to SIQS/SSS. Both control and corpus hashes are
pinned in `inputs/controls/p52_b2_protocol.json`.

`inputs/corpora/p52_b2_corpus.json` supplies independent Pocklington/trial
certificates, generation seed 2026100952, disjoint training/held-out inputs and
algorithm seeds 7/19/41. Three inputs per small, medium and uneven cohort and
one fixed balanced campaign/structured input per split form a deliberately
bounded study. Repeated timing samples are not new independent inputs.

| Case | B1 / B2 / curves | Candidate D | Scope |
| --- | --- | --- | --- |
| Small balanced | 50 / 2,000 / 8 | 8, 16, 24 | Whole factoring |
| Medium balanced | 200 / 20,000 / 16 | 32, 64, 96 | Whole factoring |
| Uneven, 34-bit smaller prime | 2,000 / 147,396 / 8 | 128, 384, 768 | Whole factoring |
| Balanced 266-bit campaign | 11,000 / 1,900,000 / 4 | 512, 1,024, 2,048 | Finite campaign; not a success-rate claim |
| Prime-cube structure | 50 / 2,000 / 8 | Small-cohort choice | Whole factoring regression control |

Every row grants 50,000,000 work units and 120-second wall/CPU caps per input
and seed, a 329-bit envelope and 16 MiB owned workspace. Retained programs use
8 MiB. The paired dense table has D/2 points plus its recurrence state; its
coexisting certificates, products, recovery and checkpoint storage are reserved.
Wall/CPU censoring prevents an accepted timing claim. Partial results still
reconstruct and are validated against the certified factors.

Each arm receives at least three seconds of validated warmup and nine samples
on PyPy implementing Python 3.11. Relative IQR above 0.15 or changed deterministic
outcomes triggers 5 seconds/31 samples and then 8 seconds/63 samples; unresolved
instability stays inconclusive. Training chooses the fastest stable D without
a completion regression against either control, per bound tier; the structured
control inherits the small tier. The saved training output freezes that choice
before held-out execution and rejects subsequent source changes. The unchanged
roadmap promotion gate applies to fresh whole factoring; finite nonsplitting
campaign savings cannot promote defaults. Timing intervals are conditional on
the fixed cohort, and completion intervals cluster algorithm seeds by fixture.
This small study cannot calibrate population-wide factor-size allocation.

Costs include config, classification, setup, schedule/table construction,
products, recovery, recursion and final checkpoint serialization. Optional cold
rows use nine separate starts and are stored apart. Instrumented profiling,
work/phase counters and stage throughput must be reported separately. RSS
includes interpreter/JIT/warmup; it is not the owned-workspace cap.

The accepted runs followed explicit B1/A10 handoffs and held the machine-wide
exclusive flock at `/private/tmp/factor-performance.lock`. Reproductions must
coordinate the same window. The runner also checks for
competing benchmark/test interpreters and fails closed if inventory is denied.
It never kills unrelated work. Output paths refuse overwrites and preserve
unstable attempts as well as accepted captures.

```sh
mkdir -p v2/benchmarks/results/p52-b2
pypy3 -B -m v2.benchmarks.p52_b2 --phase training --diagnostics --output v2/benchmarks/results/p52-b2/training.json
pypy3 -B -m v2.benchmarks.p52_b2 --phase held_out --selection v2/benchmarks/results/p52-b2/training.json --cold --output v2/benchmarks/results/p52-b2/held-out.json
```

The optional diagnostic rows instrument at most three direct curves per arm
and cohort using seeds 7/8/9. They report setup/table/giant/recovery actions,
term products, work, peak baby-table slots and retained program storage without
time ratios. These standalone curve counts are separate from the whole
factoring cohort's seeds and recursive dispatch. They cannot establish a gain.

### Matched results and retained defaults

Production and runner source are frozen at `d7f873e`; documentation-only
acceptance follows that commit. PyPy 7.3.23 implements Python 3.11.15 on
macOS 26.6.2 arm64. All 30 training and 15 held-out final captures satisfy the
predeclared relative-IQR limit, deterministic output and uncensored-budget
requirements. Small `paired_1`, medium `paired_0`/`paired_1` and structured
`regenerated` extended from nine to 31 samples; all attempts remain local.
Every accepted capture has at least three seconds of validated warmup.

Training selects D=24/64/768/2048 for the four bound tiers before held-out
execution; the structured case inherits D=24. These select the least costly
tested paired option and do not recommend enabling it. Cohort medians include
all inputs and algorithm seeds, not one curve or one successful split:

| Held-out case | D | Streamed ms | Programs ms | Paired ms | Complete starts, each arm | Paired cost increase vs streamed / programs |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Small balanced | 24 | 1.207 | 1.238 | 2.189 | 9/9 | 81.3% / 76.8% |
| Medium balanced | 64 | 3.647 | 4.062 | 9.270 | 9/9 | 154.2% / 128.2% |
| Uneven | 768 | 154.148 | 155.629 | 223.469 | 9/9 | 45.0% / 43.6% |
| Balanced campaign | 2048 | 1685.998 | 1509.368 | 1938.743 | 0/3 | 15.0% / 28.4% |
| Prime cube | 24 | 0.146 | 0.148 | 0.153 | 3/3 | 4.8% / 3.2% |

The conditional 95% bootstrap intervals for paired cost increases against
streamed are 76.2–117.3%, 149.7–160.0%, 41.5–46.0%, 14.2–16.0% and 0.3–6.7%,
respectively. Intervals against programs also exclude a benefit. These describe
timing on fixed cohorts, not population factoring performance. All completion
differences are zero; fixture-cluster intervals are [0, 0] on this small set
and do not establish population equivalence. Training completion is also
matched: 9/9 small/medium, 8/9 uneven, 0/3 campaign and 3/3 structured. Campaign
rows are finite search costs; all four declared curves finish per start, and
no factor is found. No row is wall/CPU-censored.

The three-curve campaign diagnostics explain a limitation of this bounded
layout. Streamed and unpaired programs produce 422,082 terms. Paired D=512
produces 406,632 (3.7% fewer); D=1024 and D=2048 still produce 422,082 because
the 1,024-slot prime segments separate opposite signs. All paired choices
certify exactly 422,082 eligible primes. Their giant advances are 5,535/2,766/
1,383 and baby-table slots including the sentinel are 257/513/1,025. Lower D's
limited occupancy benefit does not pay for its extra recurrence/program work.
Cross-segment coalescing or a different segment grid was not evaluated here.

For D=2048, curve work is 2,024,553 initially and 587,083 on each subsequent
curve, versus programs' 1,602,357/446,275 and streamed's 1,458,993 per curve.
The paired store retains 6,863,584 bytes and records 923 coverage misses and
1,846 hits over those three curves. Its owned reserve is 12,948,696 bytes,
below the common 16 MiB cap; unpaired programs reserve 11,053,712 bytes.
Held-out campaign process RSS is 91.5 MiB paired, 85.0 MiB programs and 82.8
MiB streamed, including interpreter/JIT overhead. Regeneration-only training
costs 2.588 seconds per campaign cohort, versus 1.878 seconds for retained
programs at the same D=1024. The cap remains effective when retention is off.

Separate instrumented full-campaign profiles validate all three arms and
record 11,076 coverage calls / 2,769 certificate compilations in paired
execution, along with packing, decoding and replay bookkeeping. This locates
additional work but supplies no accepted timing ratio. Allocation tuning,
kernel changes and advanced wheel/common-Z plans stay in their own workstreams.

The 135 cold records are nine fresh starts for each held-out arm. Median
startup-plus-cohort seconds (streamed/programs/paired) are
0.163/0.166/0.152 small, 0.195/0.205/0.184 medium, 0.415/0.438/0.485 uneven,
2.028/1.849/2.312 campaign and 0.152/0.152/0.128 structured. These include
imports, proof checks and the controls' private snapshot reconstruction, so
small cold differences are harness costs, not evidence for promotion.

Raw captures, profiles and check logs remain local under ignored
`results/p52-b2/`. The training/held-out SHA-256 digests are
`569e06e22faf5210c03d9407de479151c99cc73559f81afab754d9ec019ea225` and
`4257fc45af79f0e8131aca86bedbd65118bcb6081e263f34b93c1d8589901034`.
The frozen protocol digest is
`bd31511643617b3a4df124454573335edb49c6cb3f149b700a6ce2cf50a24ab8`.
Required corpus/control/protocol inputs are versioned. The committed-files-only
verification includes the unchanged `v1/` hash fixtures; an initial archive
omitted them and its missing-file failure is retained alongside the repaired
archive checks. Both full test suites and all 59 benchmark imports pass there;
the proof/control loaders work without generated evidence. All 135 cold
outcome/work records agree with the corresponding warmed records. No factoring
source changed after training selection.

The design retains Montgomery's x-coordinate symmetry and projective cross
differences ([1987 paper](https://wstein.org/edu/124/misc/montgomery.pdf)).
GMP-ECM's [stage-two implementation](https://github.com/sethtroisi/gmp-ecm/blob/main/stage2.c)
illustrates separate memory/cost modeling and polynomial continuations; its
native thresholds and advanced pruning are not imported into this PyPy tranche.
Its [library contract](https://github.com/sethtroisi/gmp-ecm/blob/main/README.lib)
and [stage-one code](https://github.com/sethtroisi/gmp-ecm/blob/main/ecm.c)
explicitly account for higher powers of old primes when B1 increases. A6 has
not supplied that exact schedule-ratio contract here, so increased-B1 extensions
are deferred. No external code was copied into the implementation.

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

## B4 / P4.2 bounded arithmetic kernels — 9 October 2026

This separate study uses the immutable `bcf5f3d` mainline package in
[the versioned control](inputs/baselines/b4_mainline.json). The candidate
sources and required inputs are frozen at `08ccac2` in
[the source manifest](inputs/controls/b4_freeze.json),
[protocol](inputs/controls/b4_protocol.json) and
[certified corpus](inputs/corpora/b4_corpus.json). The
[source-linked research and proofs](b4_research.md) document GMP-ECM,
CADO-NFS, AVX-ECM and Yamaquasi revisions/licenses, native versus PyPy
applicability, a24 conventions, exceptional points, normalization and widths.
No upstream code is copied. C6 chains, A6 p−1, C5 reducers, new backends and
curve families remain separate; production source and `v1/` are unchanged.

### Frozen design and rerun

Five candidates compare with the readable ladder: integer `**2` squares;
a fused addition/doubling helper; an inlined whole ladder; that loop with
selected AA/BB/U/V reductions; and that loop with unit-checked normalized
fixed difference. Initial doubling in the fusion arms uses the explicit
square helper. Whole-run comparisons include its cost. Normalization pays
its checked inverse at every scalar entry; failed inversions retain the GCD
and return proper factors or retry through the private bounded adapter.

Full factoring uses nine independent fixtures per split and seeds 7/19/41:
two inputs in each 64/128/256/329 target-bit class, with certified small
factors of 22/26/28/30 bits, plus a prime-cube control. Actual products span
51–330 bits; target labels are generation classes, not exact product widths.
Training and confirmation each contain 27 input/seed trials per arm.
The stage controls are four independently certified balanced inputs and
three curve seeds. Stage timings separately pay setup/conversion and the
B1 lcm in stage one, and prime generation, table initialization, products
and recovery in stage two. Their fixed bounds are 1,000/50,000.
These direct stages complement actual chunked portfolio measurements.
Kernel diagnostics use 30 Suyama setups at ten widths from 64 through
1024 bits and three sigmas, with a fixed 128-bit scalar. Output action,
GCD/degeneracy status and complete/partial reconstruction are validated.
Kernel setup and conversion costs are included; they are not full-run gains.

All full arms use the same five-million-work, 30-second wall/CPU, 16-MiB
owned-workspace and 331-bit limits. ECM-only tiers are 200/20,000/16 in the
small class; 1,000/50,000/16 in the three wider classes; 50/2,000/8 for powers.
Native integers remain default; persistent mpz is a separate supported
track, not a backend crossover study. The available environment is PyPy
7.3.23 / Python 3.11.15, gmpy2 2.3.1 / GMP 6.3.0, macOS 26.6.2 ARM64 on M4.

Each worker receives at least three seconds of validated warmup and nine
samples. Relative IQR over 0.15 triggers 5 seconds/31 samples, then
8 seconds/63 samples; instability or censoring prevents promotion. Training
selects the fastest stable full arm per backend without a >5-point class
completion regression. Selection is saved before fresh inputs are measured.
The roadmap gate requires zero failures and >=10% fresh pooled median
reduction with a conditional 95% timing interval above zero, or >=10-point
completion gain with a fixture-cluster interval above zero. Fixed cohorts
and seeds limit population inference. Kernel/stage wins cannot promote.

The machine-wide lock is `/private/tmp/factor-performance.lock`, with
`factor-performance-owner.json` alongside it. C6 and A6 explicitly queued
heavy checks/timings under this same lock. The runner polls for competing
benchmark/test interpreters and fails closed if process inventory is denied.
Do not run full checks during its accepted timing window. Raw JSON, cold
captures, research downloads/manifests, logs and profiles stay in the ignored
`v2/benchmarks/results/b4/` tree. The controls/corpora/runner remain versioned.

```sh
mkdir -p v2/benchmarks/results/b4
# Use the project's GMP-enabled PyPy 3.11 venv for both declared tracks.
v2/.venv/bin/python -B -m v2.benchmarks.b4_profile --output v2/benchmarks/results/b4/new-profile.json
v2/.venv/bin/python -B -m v2.benchmarks.b4_study --phase training --output v2/benchmarks/results/b4/new-training.json
v2/.venv/bin/python -B -m v2.benchmarks.b4_study --phase held_out --selection v2/benchmarks/results/b4/new-training.json --scopes full --cold --output v2/benchmarks/results/b4/new-confirmation.json
v2/.venv/bin/python -B -m v2.benchmarks.b4_report --training v2/benchmarks/results/b4/new-training.json --confirmation v2/benchmarks/results/b4/new-confirmation.json --output v2/benchmarks/results/b4/new-report.json
```

The generation and freeze commands are one-shot creation tools, not required
rerun steps; they refuse overwriting the immutable inputs/manifest. Source
changes after freeze require a separately declared study rather than editing
the existing controls. Cold starts comprise nine separate processes per
confirmed arm/backend and include harness import/source/corpus verification;
they are reported apart from warmed evidence. RSS includes the interpreter,
JIT and warmup and is not the owned-workspace allowance.

### Separately instrumented baseline breakdown

Before candidate timing, private mainline profiles measured chunked portfolios
on the four balanced controls with algorithm seed 7 and the declared
class-specific tier caps. Instrumented timers attribute costs by job
phase; cProfile supplies independent function totals. These are diagnostic,
include instrumentation overhead, and are not comparable with warmed stage
or full-run times. The smaller 64-bit control ends earlier, so its cost shares
are not representative of a long campaign.

| Track / bits | Stage one % of instrumented total | Stage two incl. baby/giant setup % | Point add/double % | Ladder self % |
| --- | ---: | ---: | ---: | ---: |
| int / 64 | 21.6 | 53.9 | 11.0 | 10.4 |
| int / 128 | 43.3 | 45.7 | 17.8 | 20.4 |
| int / 256 | 42.8 | 45.9 | 23.9 | 16.2 |
| int / 330 | 45.0 | 44.5 | 28.1 | 14.7 |
| mpz / 64 | 38.0 | 54.2 | 43.6 | 1.9 |
| mpz / 128 | 49.7 | 49.1 | 53.5 | 1.9 |
| mpz / 256 | 49.8 | 48.8 | 53.5 | 2.2 |
| mpz / 330 | 49.8 | 48.9 | 53.4 | 2.2 |

On native 128/256/330-bit controls, stage one's point formulas take
36.1/49.2/55.5% of that phase and ladder self time 46.5/37.2/32.2%.
Prime scheduling is 0.8% and GCD 0.5–1.1%. In the stage-two term loop,
point adds take only 2.0–2.7%, prime scheduling 15.7–18.5% and batch/GCD
checks 4.9–7.8%; its residual includes cross products, modular product
accumulation, reservations, cursor handling and instrumentation overhead.
Baby-table point arithmetic is reported in its separate initialization phase.
On mpz, stage-one point formulas take 93.9–94.2%; stage-two recurrence point
adds only 4.6–4.7%, with most time remaining in cross/product arithmetic and
extension calls. Source-level operation counts therefore miss large costs.

For a component share s and fractional component saving r, Amdahl's
whole-run time saving is at most s*r and its ideal speedup is 1/(1−s).
The native point-only shares above cap total savings at 17.8–28.1%
(1.22–1.39x), even if their time vanished. Including all measured ladder
self time raises the diagnostic ceiling to 38.2–42.9% (1.62–1.75x).
The mpz point-only ceiling is about 53.4% (2.15x). A first cProfile-only
pass attributed native point shares of 13.9–26.9%, illustrating measurement
perturbation; neither ceiling is an uninstrumented forecast. Actual stages
and full runs below govern decisions.

### Matched candidate comparison and fresh confirmation

All 36 final training captures are stable: 27 use nine samples, eight extend
to 31, and one to 63. All four fresh full-run captures use nine stable
samples after validated warmup. The 46 training attempts are retained locally;
there are no censored or correctness-failing accepted samples. All six training arms and both confirmed selections/baselines complete their
27 trials with matching factors, certainty, logical work and
unresolved-cofactor reconstruction.

The table reports percentage **time reduction** against the corresponding
backend's baseline; a negative number is a loss. Kernel diagnostics aggregate
the ten declared widths, so their savings are not a below-100-digit prediction.

| Candidate | int kernel | int direct stages | int full training | mpz kernel | mpz direct stages | mpz full training |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Explicit squares | −48.61% | −79.92% | −45.77% | −5.02% | −2.17% | −3.73% |
| Fused step helper | −0.31% | +2.17% | +0.41% | +6.51% | +5.34% | +2.95% |
| Whole ladder | +0.73% | +2.44% | +1.22% | +6.58% | +7.17% | +5.54% |
| Selected reductions | +30.56% | +8.51% | +8.79% | +3.69% | −1.06% | −1.73% |
| Unit normalization | −7.78% | +2.15% | +0.85% | +3.73% | +9.93% | +5.67% |

For the native selected-reduction arm, summed stage-one medians fall from
31.55 to 26.77 ms (15.2%), while stage two is 26.87 versus 26.63 ms (0.9%).
For the selected mpz normalized arm, stage one is 106.05 versus 89.70 ms
(15.4%) and stage two 90.24 versus 86.48 ms (4.2%). These sums cover the
same twelve balanced input/curve trials. Phase medians are computed
separately and need not add exactly to the median combined cohort time.
The unchanged stage-two recurrence and term-product loop explain why a
large ladder diagnostic saving shrinks in actual factoring.

Training froze `reductions` for native int and `normalized` for mpz. The mpz
normalization/whole-ladder training medians differ by only about 0.12
percentage points; this is not evidence that normalization is universally
best. The deterministic predeclared selector still confirms the least
observed training cost without further tuning.

| Fresh full cohort | Baseline median | Selected median | Reduction | Conditional 95% interval | Completion |
| --- | ---: | ---: | ---: | ---: | ---: |
| int / selected reductions | 282.25 ms | 257.47 ms | 8.78% | 5.62% to 10.92% | 27/27 both |
| mpz / unit normalization | 627.79 ms | 598.10 ms | 4.73% | 3.02% to 6.51% | 27/27 both |

Fresh native class reductions are 7.74/1.13/6.75/9.46% for the
64/128/256/329 target classes and 1.73% for powers. The mpz reductions are
0.65/3.77/6.78/4.18%, with a 1.94% loss on the tiny power control.
Completion change is zero in every class; fixture-cluster intervals are
[0,0] on these finite cohorts. This does not prove zero population risk.
The larger classes dominate pooled time. Timing intervals resample fixed
cohort samples, not independently chosen machines or general composites.

Separate nine-process cold medians are 627.20/596.42 ms for native
baseline/reductions and 1,033.57/1,036.54 ms for mpz baseline/normalization.
They include immutable harness loading and certificate verification; they do
not establish application cold-start gains. They are excluded from promotion.

**Decision: retain the production baseline and native-int default.** Neither
selected arm reaches the >=10% fresh median threshold, and completion cannot
improve beyond the observed 27/27. The native timing interval crosses 10%,
so it does not establish a >=10% effect; the bounded study ends here as
requested. Explicit squares lose, fused/whole-ladder changes do not meet the
complete-run gate, and earlier reductions or normalization remain research
controls. No reducers, new backends, extra normalization policies or curve
families are added to force a win. B4's bounded comparison is settled;
broader workload/production-bound calibration and combined portfolio
confirmation remain separate work.

The eight new arithmetic tests cover independent affine/CRT/prime-power
oracles, exact readable-formula agreement, unit scaling, degeneracy,
failed-inversion factor recovery, cancellation and canonical resume on int
and available mpz. Final test/lint and committed-only checkout results are
recorded in the changelog. This branch is prepared for integration of the
study, fixtures and retain-baseline decision; no production kernel is
proposed for promotion and no merge is performed.

### Acceptance and retained local evidence

At source/test/document commit `93e8f69`, `make -C v2 test` with the
GMP-enabled PyPy interpreter passes **409 tests** from a committed-files-only
`git archive` checkout. All **73 benchmark modules** import there, required
B4 input/control loaders validate, and the frozen source manifest matches.
`make -C v2 lint` passes Ruff checking/formatting and pycodestyle in the
worktree. The later acceptance note changes documentation only. Production
ECM, stage jobs, schedules, arithmetic, portfolio and `v1/` have no diff
against mainline. The branch remains separate and unmerged.

Raw evidence is retained locally, including all 46 training attempts, four
fresh captures, 36 cold starts, two separately instrumented profile captures
and the committed-only verification logs. The terminal capture hashes are:

- `training.json`: `9c7977513b7238ab888b8282da672c23d764d371407241fd211e08f3e3867195`
- `confirmation.json`: `332c223612cca3a75dad4d283cf400dda7d7dcfd3716858e074a8db28c8aa406`
- `report.json`: `09e05266448be940b2f667dfd11f77e3835002a28ba2447ba8684f9d139d20eb`

These ignored captures are not prerequisites for importing or rerunning the
committed study. The frozen protocol, source control and certified inputs are
required and versioned. B4's exclusive window is released to C6/A6; no further
candidate experiments or heavy checks are planned in this tranche.
