# Factor v2

Integer factorization using **PyPy implementing Python 3.11** and exact
standard-library integers. v2 repairs the [original v1](../v1/) implementation
and adds bounded, resumable execution and experimental relation-based engines.

## Methods

| Method | Role and status |
| --- | --- |
| Exact preprocessing | Trial division, primality classification, powers of two and perfect powers; optional bounded Fermat search in the portfolio |
| Brent rho | Seeded walks with finite attempts, batched GCDs and bounded recovery |
| Pollard p−1 | Two-stage smoothness search with saturation replay; integrated into the bounded portfolio |
| Montgomery ECM | Suyama curves, binary ladder and recoverable stage-two batches; the default portfolio's final search stage |
| QS / MPQS | Reference polynomials and verified relation pipeline, useful as comparison controls |
| SIQS | Shared relations across CRT/Gray polynomial families, incremental roots, filtering, GF(2) dependencies and exact extraction; selectable CLI engine or optional portfolio fallback |
| SSS / SSSf | Experimental Smooth Subsum Search; SSSf adds an optional lossy candidate filter; explicitly selectable through the CLI |

The bounded default runs preprocessing, rho, p−1 and ECM serially. SIQS and SSS
require explicit selection; thread/process SIQS workers are experimental.
The compatibility `factorize()` interface remains available but has no shared
global deadline or resumable schedule.

## Improvements over v1

v2 fixes arithmetic and sieve boundaries, validates proper divisors, preserves
failed recursive cofactors, and distinguishes probable from proven primes.
Rho retries and ECM curve counts have explicit limits; saturated batches recover
or report failure. ECM uses the binary ladder; the unsafe original PRAC
chain is not the default.

`ecm.multiply_prac(k, x, z, n, a24)` now executes experimental verified PRAC
records. It returns a valid projective pair or raises
`prac.NonunitPointError`: its `factor` attribute is a proper divisor, or
`None` requests a curve retry. This is a deliberate change from the former
ladder wrapper; callers must handle that exception. Coordinates can differ
from the ladder by projective scaling. Supply a point on a nonsingular
Montgomery curve over an odd modulus, with `a24=(A+2)/4`, as returned by
`ecm.setup_curve`.

Generation uses exact rational splits for the odd part of scalars of at most
32 bits, with at most 30 candidates and 512 instructions each. This bounds
each prime/prime-power multiplier, not the integer being factored: the
modulus can have 40–80 digits or more. A 512-record
LRU cache contains immutable integer records, never curve points. Each
record receives a separate integer/differential verification. Larger
scalars use a checked ladder without chain search. Zero, one, powers of two,
infinity and the order-two point have explicit handling.

Execution checks intermediate projective states and reports nonunit factors.
An exceptional difference triggers at most 513 additional GCDs on retained
coordinates and one ladder retry from the original point; `(0,0)` is never
a successful point. `prac.get_chain`, `verify_chain`, `clear_cache` and
`cache_info` expose the bounded record interface. `prac.multiply(...,
chain=record)` re-verifies caller-supplied records before execution. Cost
weights are positive integers of at most 32 bits; defaults model 4M+2S for
addition and 3M+2S for doubling. They are not runtime speed estimates.

This completes the P4.1/A4 correctness tranche. Both `factorize_ecm` and the
bounded stage jobs still use the binary ladder. The [benchmark guide](
benchmarks/README.md#p41a4-verified-prac-5-october-2026) separates kernel
diagnostics from complete two-stage attempts on certified 40–80-digit inputs,
including matched optional gmpy2 arms. These are experimental comparisons;
B3 owns program composition, shared work accounting, checkpoint/replay
integration and complete-factorization comparisons. The standalone PRAC
helper has no portfolio budget or checkpoint contract.

The bounded portfolio adds one allowance across preprocessing, retries and
recursive children, with streamed prime schedules, controlled workspace and
validated checkpoints. Later work adds exact relation provenance, shared SIQS
stores, compact matrix storage and recovery. Performance repairs reduce repeated
primality/root work, sparse exponent recovery, prime-power inversions and
relation verification. New or untrusted data still receives exact verification.

## Run

Run all commands below from the repository root. `pypy3` must implement Python
3.11; `make -C v2 runtime` reports the interpreter. Use `PYTHON=/path/to/pypy3`
when selecting a different installation of the same supported runtime.

```sh
make -C v2 run N=626100403 SEED=7
make -C v2 run-bounded N=626100403 SEED=7
pypy3 -m v2.factor -15
pypy3 -m v2.factor 626100403 --bounded --seed 7 \
  --work-limit 2000000 --seconds 30 --cpu-seconds 30 --memory-mib 8
```

Omit the integer for an interactive prompt. `--verbose` shows diagnostics;
library calls remain quiet. Exit status is **0** for complete results, **1**
for unresolved composites and **2** for invalid input. Zero is rejected, one
has an empty factorization, and negative inputs retain their sign.

Select QS, MPQS or SIQS directly, or enable SIQS as the portfolio's final fallback:

```sh
pypy3 -m v2.factor 10002200057 --method siqs --seed 7
pypy3 -m v2.factor 10002200057 --method qs --qs-base-bound 400 --seed 7
pypy3 -m v2.factor 10002200057 --method mpqs --qs-base-bound 400 --seed 7
pypy3 -m v2.factor 10002200057 --siqs --ecm-curves 2 --seed 7
```

`--method qs|mpqs|siqs` runs exact preprocessing followed by the selected
engine, disabling rho, p−1 and ECM. `--siqs` runs preprocessing → rho →
p−1 → ECM → SIQS. These selections
use the existing recursive portfolio: a validated split sends both children
back through classification and factoring until only terminal factors or
explicit unresolved cofactors remain. Omit the integer for `Enter number:`.
Neither selection requires an additional `--bounded` flag.

QS/MPQS/SIQS selections default to 200 million work units, 30 wall/CPU seconds
and 80 MiB owned workspace; 16 MiB is reserved outside the sieve job for parent
and schedule/checkpoint coexistence. These are finite starting allowances,
not a calibrated general-number policy. For a longer attempt:

```sh
pypy3 -m v2.factor YOUR_INTEGER --siqs --ecm-curves 2 \
  --seconds 300 --cpu-seconds 300 --work-limit 1000000000
```

SIQS runs only if the earlier stages finish with work and time remaining.
Increasing the deadline does not extend an exhausted polynomial schedule.
Search controls include `--qs-base-bound`, `--qs-half-width`,
`--qs-max-half-width`, `--qs-factor-count`, `--qs-family-count`,
`--qs-pool-size`, `--qs-assignment-policy`, `--qs-polynomials-per-family`
and `--qs-residual-bound`; each also accepts its `--siqs-*` alias. They retain
the corresponding `SIQSConfig` defaults when omitted. QS uses one polynomial;
MPQS/SIQS use finite family schedules. Wider intervals need a compatible
maximum width. Streaming assignment (`nearest` or `flyer`) and Gray quotas
apply only to SIQS; extended quotas require a streaming policy. Invalid
mode/configuration combinations are rejected. These settings remain experimental.
Use `--verbose` to inspect method outcomes and work use. Measured automatic
handoff/default selection is open in the [roadmap](ROADMAP.md); explicit CLI
usability is tracked separately under P3.4.

Select an experimental SSS engine explicitly:

```sh
pypy3 -m v2.factor 626100403 --method sss --seed 7
pypy3 -m v2.factor 626100403 --method sssf --seed 7
```

These selections disable rho/p−1/ECM and use bounded execution after exact
preprocessing. Defaults are 200 million work units, 30 wall/CPU seconds and
80 MiB owned workspace. `--sss-base-bound` and `--sss-rounds` set finite search
parameters; success is not guaranteed.

## Library and result contracts

```python
from v2.factor import factorize, print_factorization

result = factorize(626100403, seed=7)
assert result.reconstruct() == result.original
print(print_factorization(result.original, result))

for factor in result.factors:
    print(factor.value, factor.exponent, factor.certainty.value)
```

`FactorizationResult` carries `original`, `sign`, prime factors and `remaining`
cofactors. `complete` means no unresolved composite remains; `proven` additionally
requires proven primality for every terminal factor. The deterministic
Miller–Rabin domain is strictly below `2**64`; larger surviving candidates are
labeled probable primes. An unresolved result still reconstructs its input.

Splitter functions return a proper divisor or `None`. Prime sieve bounds are
half-open; p−1/ECM B1 and B2 bounds are inclusive. Import through the `v2` package
rather than mixing both version directories on `sys.path`.

Use the bounded API for shared limits and resumable execution:

```python
from v2.budget import Budget
from v2.portfolio import PortfolioConfig, factorize_bounded

config = PortfolioConfig()
run = factorize_bounded(
    626100403, seed=7, config=config,
    budget=Budget(work_limit=1000, seconds=30, cpu_seconds=30),
)
assert run.result.reconstruct() == 626100403

resumed = factorize_bounded(
    run.result.original, config=config, checkpoint=run.checkpoint,
    budget=Budget(work_limit=2_000_000, seconds=60, cpu_seconds=60),
)
```

`run` exposes the result, stop reason, consumed work, events and checkpoint.
`PortfolioConfig(siqs=SIQSConfig(...))` enables optional SIQS after ECM;
`PortfolioConfig(sss=SSSConfig(...))` enables SSS. Their classes live in
`v2.qs` and `v2.qs.sss`, respectively. Advanced settings and exact relation
contracts are documented in the module docstrings and covered by the tests.

## Optional ECM programs and explicit campaigns

`PortfolioConfig(ecm_program_bytes=...)` opts into P5.2 A3's immutable packed
prime/power blocks. Each block owns half-open endpoints and, for stage one,
the inclusive B1 identity. Completed blocks are reused across curves and
recursive cofactors; points, residues, products and recovery remain private to
each job. `0` retains streamed execution and the existing version-4 checkpoint
schema. This is an experimental storage option, not a promoted default.

The program cap is part of `memory_bytes`, with a scratch reserve of
`4096 + 256 * segment_size` bytes and a 512-byte allowance per retained block
plus packed payloads. A full cap causes regeneration rather than eviction or
an unbounded allocation. Packed programs require B2 strictly below `2**64`;
the streamed integer API retains its existing endpoint domain. Owned reserves
are conservative estimates, not process RSS limits.

Generation reserves `segment_size + len(base_primes)` units per block,
plus one unit per prime for packing/reading and one per compiled stage-one
power. A retained block reserves one unit per decoded prime. Stage-one copying,
point arithmetic and recovery keep their existing charges. Program and streamed
work counts therefore differ; reduced work counts alone establish no speedup.

Programs are run-local and omitted from checkpoints. Opt-in snapshots use
version 5 and the `ecm-packed-blocks-v1` identity; old version-2/3/4 snapshots
remain readable with programs disabled. Resume preserves the prime buffer,
curve assignment, recovery and consumed allowances, but charges for regenerating
missing future blocks. Powers for an already-buffered resumed segment can be
recomputed under the existing copy reservation. A refusal during compilation
publishes no partial block or advanced cursor, although completed generation
work stays consumed.

For an explicit finite campaign, declare all curve tiers up front and choose
work, time and storage together. A 329-bit envelope admits every integer below
100 decimal digits and avoids the default 4096-bit coordinate reserve:

```python
config = PortfolioConfig(
    rho_attempts=0, pm1_attempts=0,
    ecm_tiers=((11_000, 1_900_000, 10),),
    max_input_bits=329, memory_bytes=16 * 2**20,
    ecm_program_bytes=8 * 2**20,
)
run = factorize_bounded(
    n, seed=7, config=config,
    budget=Budget(work_limit=50_000_000, seconds=300, cpu_seconds=300),
)
assert run.result.reconstruct() == n
```

These are caller-selected allowances, not calibrated factor-size tiers or a
success guarantee. Extend a paused campaign by increasing **total** allowances
under the identical configuration; completed curves and their RNG progress are
credited. An exhausted schedule stays exhausted. Adding curves/bounds to a
checkpoint, or extending B1 on the same curve, remains unsupported pending
B2/A6: increasing B1 needs missing powers of old primes as well as new primes.

`v2.ecm_programs.pair_coverage()` supplies bounded immutable +/- coverage
certificates for the later B2 implementation. It includes direct-scalar
exceptions, positive recurrence initialization and block tails. This bounded
compiler currently accepts D=0 (direct scalars), or even D>=2 with
`2*D < B1` for odd B1 and `2*D < B1-1` for even B1. Production
stage two still executes the existing unpaired terms; D tuning, paired recovery,
wheel pruning and common-Z tables retain their roadmap gates.

## Checkpoints and limits

```sh
mkdir -p v2/audit/results
pypy3 -m v2.factor 626100403 --bounded --seed 7 --work-limit 1000 \
  --checkpoint v2/audit/results/example-checkpoint.json
pypy3 -m v2.factor --resume v2/audit/results/example-checkpoint.json \
  --work-limit 2000000 --seconds 60 --cpu-seconds 60
```

Resume with the same configuration, including ECM curve count, memory, Fermat
and selected-method options. Increased budgets are **total allowances**, including
previously consumed resources. Paused time is excluded. Resuming an exhausted
local candidate schedule does not reset it; use different allowances in a new
run when that schedule has no remaining work. Checkpoints verify configuration,
checksums, reconstruction and arithmetic; checksums are not authentication.

For QS/MPQS/SIQS, repeat the same `--method` (or `--siqs` fallback) and search
options when resuming. For example, replace `siqs` with `qs` or `mpqs` in both
commands to resume that engine:

```sh
pypy3 -m v2.factor 10002200057 --method siqs --qs-base-bound 400 \
  --seed 7 --work-limit 130000 \
  --checkpoint v2/audit/results/siqs-checkpoint.json
pypy3 -m v2.factor --method siqs --qs-base-bound 400 \
  --resume v2/audit/results/siqs-checkpoint.json --work-limit 200000000
```

Time and cancellation checks are cooperative: an in-progress bigint operation
finishes before the next check. Workspace caps cover conservative owned storage,
not total process RSS, interpreter/JIT memory or caller-retained output.
Work units are algorithmic reservations, not seconds or comparable effort across
unrelated methods. The default bounded portfolio caps input size at 4096 bits;
that admission limit is not a practical factorization-size claim.

## Benchmarks

Historical M9 measurements on PyPy 7.3.23 / Python 3.11.15, macOS arm64:

| Complete-factorization batch | Emulated v1 | Repaired v2 (M9) | Elapsed-time reduction |
| --- | ---: | ---: | ---: |
| Original: 5 inputs × 5 seeds | 1.155 ms | 1.078 ms | 6.6% |
| Separate control: 56 inputs × 5 seeds | 21.091 ms | 16.994 ms | 19.4% |

Times cover entire batches, with at least three seconds of validated warmup and
15 samples. Every timed result has checked factors, multiplicities and exact
reconstruction. v1 uses a `lib2to3` syntax/integer-division adapter on the same
PyPy interpreter; this is not a native Python 2 comparison. These measurements
predate later v2 work and do not describe current performance on all inputs.
Invalid baseline outputs receive no speed ratio.

Later v2-to-v2 comparisons show where the newer work helps:

- A bounded 20-digit Phase 2 confirmation raised completion from **70.8% to
  91.1%** under matched 50 ms operation caps. Cohort cost includes unfinished
  outcomes and is not a successful-factorization speed ratio.
- P3.6.1 held-out complete 30-digit cohorts reduced elapsed time by **70.7%
  for SIQS, 77.2% for SSS and 78.6% for SSSf** against their frozen v2 controls;
  all 72 timed attempts per arm completed under the declared configurations.
- SIQS dependency cadence 32 reduced a later complete 30-digit cohort by
  **38.7%** against R3 cadence 1. It remains an explicit option; the default
  cadence stays 1. Parallel workers have no demonstrated universal advantage.

See the [benchmark guide](benchmarks/README.md) for protocols, commands,
uncertainty, rejected experiments and size/budget limitations. A recorded
balanced 50-digit success took **19 min 43 s**; one success does not establish
broad coverage. There is no supported “50–60 digits within a minute” guarantee.

## Development and benchmark runs

Factoring needs no third-party runtime packages. Set up the local lint tools:

```sh
pypy3 -m venv v2/.venv
v2/.venv/bin/python -m pip install ruff==0.14.14 pycodestyle==2.14.0
make -C v2 test
make -C v2 lint
make -C v2 validate
make -C v2 benchmark WARMUP_SECONDS=3 REPETITIONS=9
make -C v2 benchmark-phase-two WARMUP_SECONDS=3 REPETITIONS=9
```

Use `benchmark-phase-three-siqs`, `benchmark-phase-three-sss` or
`benchmark-phase-three-parallel` for focused comparisons. Make creates output
folders and unique timestamped names; direct runners need an existing output
folder. Use matched corpora, seeds and allowances, validate every outcome and
separate cold startup from warmed execution.

## Files and next work

`tests/` contains acceptance and arithmetic regressions; the
[QS guide](qs/README.md) covers relation engines. [Benchmark documentation](benchmarks/README.md) holds detailed results;
the [roadmap](ROADMAP.md) records acceptance gates. Benchmark `inputs/`
retain independent corpora, immutable sources, provenance and frozen controls.
The entire `audit/` tree is Git-ignored local research and diagnostic material;
its `results/` and benchmark `results/` hold generated captures, checkpoints,
profiles and scratch work. Accepted changes go in the
[root changelog](../CHANGELOG.md); the detailed development journal stays local.

Next priorities are measured arithmetic/backend improvements and a bounded
GNFS reference pipeline, followed by scaling and held-out dispatch experiments.
GNFS, double-large-prime SIQS and broader parameter/parallel promotion remain
roadmap work; existing small or configured-cohort wins do not close those gates.

## P3.8-R2 collector experiments

`SieveConfig(score_policy="fixed")` uses conservative integer bounds in
units of 1/32 bit. Prime weights round upward; norm thresholds round downward
and residual allowances upward. Eleven leading bits select a small integer
mantissa table, with adjacent bins bounding the exact logarithm. This avoids
floating-point arithmetic and large per-candidate powers. Exact exponent
recovery, residual checks and relation verification still decide admission.
`threshold_extra` uses these scaled units for this policy; positive values
remain intentionally lossy. Byte and unsigned-array scores saturate together
with their thresholds, admitting extra candidates without losing coverage.

The keyword-only `power_plan_bytes=0` preserves streamed lifting and existing
positional configuration calls. A positive allowance, up to
16 MiB, reserves that entire amount alongside the existing collector, store
and solver workspace before setup. Power/fixed scoring can retain bounded
modulus/root/weight plans for one polynomial. Fixed-polynomial and SIQS jobs
bind them to the entire polynomial interval so collection batches share
plans. Ordinary collector calls share plans only inside their checked span.
Changing the polynomial or extending the span clears the plans. Each plan's
construction and replay are charged; a cache-cap refusal falls back to
streamed lifting. Singular capped-lift fallbacks remain conservative.

Plans are disposable acceleration state. Checkpoints retain their configuration
and fully reverify restored relations, then rebuild charged plans; cached roots
are never trusted checkpoint evidence. These options preserve first-uncommitted
position semantics, exact cumulative work and finite cancellation polling.

The R2 benchmark harness compares these options with conservative scoring,
bucket/resieve recovery, tiny-prime corrections, scalar/batch smooth-part
recovery and grouped hit reservations. Its independent certified inputs,
frozen source control and separate cold/profile modes are described in
[the benchmark guide](benchmarks/README.md). Dispatcher defaults are unchanged.
