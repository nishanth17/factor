# Factor v2

Integer factorization using **PyPy implementing Python 3.11** and exact
Python integers or optional GMP integers. v2 repairs the
[original v1](../v1/) implementation
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

C6 adds reusable experimental machinery in `benchmarks.c6_chains`:
`compact(Chain)` produces an immutable four-byte register record;
`verify(Record)` independently proves its integer action, and
`Executor(record, backend)` verifies once before repeated point execution.
`build_program(B1, method, backend)` supports B1 <= 2,000 with checked PRAC,
compact PRAC, precomputed Lucas and rolling Lucas layouts. The catalog uses
303 independently verified prime records from a pinned GMP-ECM generator.
`ThreePointExecutor` accepts only recognized continued-fraction records.
These experiment APIs require the same nonsingular Montgomery curve, odd
modulus and `(A+2)/4` convention as `multiply_prac`. Results are valid X:Z
pairs or `NonunitPointError`; callers retain proper factors and handle retry.
Records have at most 512 steps and 16 retained point slots; programs own at
most 512 records, with no global point cache. See the [C6 study](
benchmarks/c6_research.md) for recovery, storage and reproducibility details.
No C6 candidate is routed into production stage jobs or checkpoint formats;
B3 integration remains open. The initial conservative C6 study retained the
ladder after full-stage and campaign losses; those results remain historical
controls. The user reopened executor optimization. See the [follow-up](
benchmarks/c6_optimization.md) for its independently proved factor coverage,
frozen comparisons and final decision.

`benchmarks.c6_fast.load_catalog()` verifies immutable scalar and coordinate
coverage certificates once. `build_program(B1, family, mode, backend, batch)`
constructs a bound-owned experimental program for binary, PRAC or upstream
Lucas records. `program(point, n, a24, extra)` returns `(point, factor)`;
a missing point with no factor denotes finite saturation/retry exhaustion.
A unit aggregate certifies the block's intermediate coordinates. A nonunit
aggregate replays the saved block through the unchanged strict recovery path,
including saturated products containing different proper factors.

The same curve/modulus preconditions apply. Programs have <=512 records,
<=512 operations per record, <=16 point registers and <=64 records per check
batch; generated source is capped at 128 KiB per record and 8 MiB per program.
Caller-owned programs can be reused; there is no global point or code cache.
`benchmarks.c6_cf` supplies the separately verified CF catalog and common-tuple
or three-point execution, with the same stage result and recovery contract.
These APIs do not provide a portfolio work ledger, cancellation checkpoint or
serialized resume format. B3 must reserve whole-block execution and bounded
replay work before integrating them.

`benchmarks.c6_b4` additionally pins the committed B4 ECM control in an
isolated module and supplies exact-residue early-reduction kernels and a
bounded independent D/A fusion pass. `build_program(family, mode, backend,
batch)` fixes B1=2,000 and accepts PRAC, Lucas or CF with `late`, `reduced` or
`fused` execution. It preserves the same guard/replay contract; the native
combined experiment does not modify production modules or checkpoint state.
Cold and reused programs have different measured costs and adoption scope.
The completed bounded study recommends reduced PRAC/batch 16 for B3's native
B4 integration experiment and tuple Lucas/batch 16 for its separate GMP work.
Native fresh construction retains the ladder; no production default changes.
The PRAC, Lucas and CF catalogs, executors and research runners are available
in mainline under `benchmarks.c6_*` for later experiments. Production ECM
continues to use the integrated B4 ladder until B3 completes stage-job,
budget and checkpoint integration.

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

The calibrated B1 balanced 30-digit SIQS preset is available explicitly.
Load its complete frozen configuration so the store, matrix and family
allowances match the measured bundle:

```python
import json
from pathlib import Path

import v2
from v2.budget import Budget
from v2.qs import SIQSConfig, SIQSJob, SieveConfig

selected_path = Path(v2.__file__).parent / (
    "benchmarks/inputs/controls/b1_selected.json"
)
selected = json.loads(selected_path.read_text())
settings = selected["configurations"][selected["selected"]["siqs"]]
settings["collector"] = SieveConfig(**settings["collector"])
config = SIQSConfig(**settings)
budget = Budget(work_limit=10**13, seconds=5, cpu_seconds=5)
result = SIQSJob(n, seed=7, config=config, budget=budget).run()
assert (result.divisor or 1) * result.cofactor == n
```

Here `n` is the integer to split. The result preserves an explicit unresolved
cofactor when its finite allowance ends. This preset reduced the tested
balanced cohort's time by 36.6%; uneven-factor regressions keep it an explicit
choice. See the [B1 measurements and limitations](benchmarks/README.md#b1-joint-qsmpqssiqs-calibration--9-october-2026).

Select an experimental SSS engine explicitly:

```sh
pypy3 -m v2.factor 626100403 --method sss --seed 7
pypy3 -m v2.factor 626100403 --method sssf --seed 7
```

These selections disable rho/p−1/ECM and use bounded execution after exact
preprocessing. Defaults are 200 million work units, 30 wall/CPU seconds and
80 MiB owned workspace. `--sss-base-bound` and `--sss-rounds` set finite search
parameters; success is not guaranteed.

## Arithmetic backends

`python-int` remains the dependency-free default. Explicit `gmpy2-mpz` selection
applies across preprocessing/primality, rho, p−1, ECM, QS/MPQS/SIQS, SSS/SSSf,
relation verification/extraction, smoothness trees and GF(2) matrix masks.
Moduli, residues, coordinates and polynomial coefficients retain `mpz` inside
their arithmetic loops. Loop indices, small-prime schedules, budgets and seeds
remain Python integers. Multiplication, reduction and XOR use those concrete
operand types directly; there is no backend callback per multiplication.

The tested optional build is **PyPy 7.3.23 / Python 3.11.15, ARM64,
gmpy2 2.3.1 / GMP 6.3.0**. Install the optional dependency into that PyPy
environment, then run:

```sh
v2/.venv/bin/python -m pip install -r v2/requirements-gmp.txt
v2/.venv/bin/python -m v2.factor 626100403 --bounded --backend gmpy2-mpz
v2/.venv/bin/python -m v2.factor 10002200057 --method siqs --backend gmpy2-mpz
```

An unavailable dependency raises an explicit error; GMP selection never falls
back to another backend or CPython. Importing and using the integer default
does not import gmpy2. GMP support is optional: faster individual operations
do not establish a faster factoring engine. See the
[matched backend study](benchmarks/README.md#p43-arithmetic-backends--5-october-2026).
The current selectors are explicit; there is no automatic digit threshold.
Algorithm/stage/size selection requires the separate production-bound and
larger-QS study. Aggregate portfolio timings do not establish that policy.

Library selection uses `factorize(..., backend="gmpy2-mpz")`, or
`PortfolioConfig(backend="gmpy2-mpz")` for `factorize_bounded`. A configured
SIQS/SSS fallback must select the same backend as its portfolio. Standalone
`SIQSConfig`, `SSSConfig` and `ParallelConfig` accept the same keyword.
Direct `factorize_rho`, `factorize_pm1`, `factorize_ecm` and `factorize_bf`
accept `backend=`; omission preserves the representation of an exact input.
`build_factor_base(..., backend=...)` and `stage_one_scalar(..., backend=...)`
also expose the boundary. Existing positional configuration arguments retain
their meanings.

Low-level arithmetic/point/polynomial objects may contain `mpz`. For helper
calls, convert once with `arithmetic.get_backend("gmpy2-mpz").integer(value)`.
`utils.gcd`, `utils.extended_gcd`, `utils.modular_inverse`, `utils.isqrt`, and
`preprocessing.integer_root` recognize exact GMP inputs. The shared
`arithmetic.pow` handles exact powers and modular powers; `arithmetic.divexact`
checks divisibility before invoking GMP exact division. `mpz / mpz` is never
used for factoring arithmetic. Failed inversion raises `NonInvertibleError`,
a `ValueError` subclass whose `divisor` retains the GCD, including saturation.
The [gmpy2 integer API](https://gmpy2.readthedocs.io/en/latest/mpz.html)
documents the underlying integer operations.

High-level factorization results, splitters' returned divisors, QS results
and serialized checkpoints contain canonical Python integers. Certainty,
witness selection, seeds, bounds and logical work reservations are shared
across the two tracks; GMP primality shortcuts do not upgrade classifications.

Unpaired GMP portfolio checkpoints use version **6**, including reusable ECM programs.
Native portfolio formats remain **4** for streamed execution and **5** for
programs. SIQS/SSS use **3**, parallel SIQS **4**, and polynomial families **2**.
They bind progress to the selected backend;
GMP identity includes gmpy2 and GMP versions. Resume rebuilds only arithmetic
values as `mpz`, retaining native counters/cursors and cumulative resources.
Backend/build mismatches are rejected. Older supported integer checkpoints
remain readable on `python-int`; they cannot silently become GMP jobs.
Pre-integration P4.3 version-5 backend snapshots also remain readable.
GIL tuning and thread promotion remain separate experiments.

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
Miller–Rabin domains have strict upper bounds: the existing two/three/seven
base tests cover smaller inputs through `n < 2**64`; the first 12 prime bases
(2 through 37) cover `n < 318665857834031151167461`, and the first 13 (through
41) cover `n < 3317044064679887385961981`. These wider guarantees rely on
[Sorenson–Webster's exhaustive computational results](https://arxiv.org/abs/1509.00864),
not a conjectured extension of v1's table. Both upper endpoints are composite
strong pseudoprimes to their preceding witness sets. Above the last bound,
survivors remain probable primes; certificate generation belongs to B14.
An unresolved result still reconstructs its input.

`utils.deterministic_bases(n)` exposes the strict witness dispatch (`None`
outside its supported domain, including `n < 2`); nonintegers are rejected,
and callers must handle small inputs/divisibility first. `classify_prime`
supplies certainty, while `is_prime` and
`is_prime_fast` remain Boolean convenience functions. `DETERMINISTIC_BASES`
retains the seven word-domain witnesses; use the selector with the expanded
`DETERMINISTIC_LIMIT`. Explicit
`use_probabilistic=True` still draws exactly `tolerance` witnesses for every
nontrivial survivor, including within the fixed ranges, and returns probable
status. A composite witness exits early; small-prime membership/divisibility
still decides trivial cases exactly.

Fresh bounded runs use the checkpoint identity `primality: mr13-strict-v1`.
Existing numeric schemas 4–6 and B2's paired/wheel schemas 7/8 retain their
backend and schedule identities; the primality policy composes with each.
A missing primality identity means `mr64-strict-v1`; it stays attached to the
whole resumed run, including pending children. Old random witness progress,
RNG consumption and conservative terminal labels remain unchanged. Unknown
identities or witness prefixes inconsistent with their policy are rejected.
No completed probable factor is silently upgraded on resume; start a fresh
run to request the wider policy. Stronger revalidation can reject a legacy
probable factor now proved composite. Historic source snapshots remain immutable.

The bounded classifier still reserves the small filter pass and then
`n.bit_length() + s` work units before each witness, where `n-1 = d*2**s`.
The new ranges use 12 or 13 fixed witnesses instead of `primality_rounds`
random witnesses, consume no RNG draws, and may therefore change later
random candidate seeds and work totals in fresh recursive runs.
`PortfolioConfig.primality_rounds` controls random tests outside that run's
policy; it does not request probabilistic mode inside fixed ranges. Decomposition
is retained once per cofactor; refused witness reservations occur before RNG
draws or progress mutation. Revalidation on resume retains existing wall/CPU
and work-accounting conventions. Cooperative cancellation and finite shared
allowances still apply; a single modular power remains atomic.

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

### Experimental SSS resume and loss policy

`SSSJob(n, seed=7, config=SSSConfig(...), budget=...)` returns the common
split/cofactor result; recursive terminal classification belongs to
`factorize_bounded`. `SSSConfig(mode="sssf", filter_bound=0)` retains two-stage
smoothness without the candidate cutoff. A positive `filter_bound` deliberately
loses candidates; it is an explicit policy, never an automatic default.
Collision trees inspect the forced-divisor quotient, while admission recovers
all exponents from the original polynomial value.

In-memory resume must extend the original `Budget` object. For serialized
resume, use `SSSJob.from_checkpoint(checkpoint, budget=total_allowance,
config=same_config)`. Native SSS versions 1/2 migrate to 3; backend/build and
configuration identities are checked. New solver fingerprints carry
`digest_encoding="hex-v1"` and avoid decimal-mask conversion limits; absent
tags retain legacy decimal replay, and unknown tags are rejected.
Restoration retains prior work/wall/CPU
and additionally charges setup, relation verification, assignment regeneration
and solver/extraction replay. A total grant equal to previously consumed work
can therefore refuse reconstruction. Checkpoint byte-cap refusal leaves the
in-memory job available; exhausted schedules do not restart on budget extension.
See the [A7 acceptance matrix](benchmarks/a7_r5_reconciliation.md) for worker
resource limitations and the prepared E1 comparison arms.

## Optional ECM programs and explicit campaigns

`PortfolioConfig(ecm_program_bytes=...)` opts into P5.2 A3's immutable packed
prime/power blocks. Each block owns half-open endpoints and, for stage one,
the inclusive B1 identity. Completed blocks are reused across curves and
recursive cofactors; points, residues, products and recovery remain private to
each job. `0` retains streamed execution and the native version-4 checkpoint
schema. GMP snapshots use version 6. This is an experimental storage option,
not a promoted default.

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
version 5 for native arithmetic or 6 for GMP and the `ecm-packed-blocks-v1`
identity; old version-2/3/4 snapshots
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
a separate ECM migration contract. A6 supplies the exact integer schedule
ratio, including missing powers of old primes as well as new primes. B2
supports continuation of the predeclared
finite campaign under cumulative work/wall/CPU allowances, including pauses
inside table construction, a paired block or scalar recovery. It does not
reinterpret an exhausted campaign as a new allowance of curves.

`v2.ecm_programs.pair_coverage()` supplies bounded immutable +/- coverage
certificates consumed by the opt-in B2 executor. It includes direct-scalar
exceptions, positive recurrence initialization and block tails. This bounded
compiler currently accepts D=0 (direct scalars), or even D>=2 with
`2*D < B1` for odd B1 and `2*D < B1-1` for even B1.

Set `ecm_pair_distance=D` together with a nonzero `ecm_program_bytes` to
execute those programs. `None` (default) preserves unpaired execution;
`0` is an explicit direct-scalar control. For example, add
`ecm_pair_distance=1024` to the 11,000/1,900,000 campaign above. This is a
caller-selected experimental setting, not a measured recommendation.
The frozen B2 study retains all production defaults: its selected paired
settings lose to streamed and reusable unpaired execution on held-out complete
factoring and a finite nonsplitting campaign. See `benchmarks/README.md` for
the D choices, costs and coverage/segmentation limitations.
The dense curve-private table holds even multiples through D, with a separate
2D recurrence step. Each certified pair contributes one cross-product;
singletons, block boundaries and inclusive tails remain covered. A saturated
product replays individual terms, then both certified primes for any saturated
term. A saturation never counts as a proper divisor.

Coverage records share the program retention cap and regenerate when it fills.
The additional conservative workspace reservation is
`8192 + 2048 * (segment_size + gcd_batch + 1)` bytes, alongside point tables,
the whole program cap, product/replay state and serialized output. Compilation
charges prime packing, certificate construction and record decoding; cached
reads charge decoding. Table additions/doublings and giant advances reserve
two units each; paired term/product actions reserve two units per term. Direct
scalars and recovery reserve their bit-length cost plus GCD/product cost.
These are algorithmic accounting units, not measured bigint operation counts.

Paired checkpoints use version 7 with `ecm-packed-pairs-v1`, exact D/config
and backend/build identity. They retain the active decoded coverage block,
curve-private table, products and recovery position. Resume verifies the prime
buffer and coverage certificates, record/table bounds and product consistency;
future missing programs rebuild under the cumulative allowance. Arithmetic
state is protected by the checksum, as for stage one. Old unpaired native and
GMP checkpoint formats remain unchanged. Performance promotion, wheel pruning,
common-Z, PRAC routing and allocation remain separate roadmap gates.

For the bounded wheel alternative, set `ecm_pair_wheel=W` instead of
`ecm_pair_distance`. W must be even, at least 2 and at most
`2 * segment_size`; a nonzero program cap is required. For example,
`ecm_pair_wheel=210` uses centers at multiples of 210 and retains only
coprime distances through 105. Program blocks end halfway between centers,
so their boundaries cannot split a pair. Initial/final partial cells and
wheel-divisor primes remain covered; primes assigned to center zero use
direct scalars. Initialization at the first positive center permits small B1;
the first W-to-2W giant transition uses doubling when its predecessor is zero.

The sparse table is generated by the existing odd-multiple recurrence with
two private scratch points. Skipped distances still incur construction work.
Workspace reserves conservatively allow a dense half-wheel plus scratch and
integer indices; fewer retained points do not imply the same reduction in RSS
or reserved bytes. Storage/decoding and mixed-factor replay use the original
program and paired contracts. Wheel checkpoints use version 8 and
`ecm-aligned-wheel-pairs-v1`; schemas 4–7 and their disabled-field encoding
remain readable and unchanged. Same-config finite campaigns can extend their
cumulative budgets. A6 supplies the exact increased-B1 integer ratio; ECM
checkpoint migration and point-specific recovery remain separate work.

This optional mode covers one nearest-center distance set. Extended sets,
relocation and overlapping-window graph matching belong to C2, as does
common-Z; polynomial continuation belongs to F3. No production default changes.
The fresh wheel comparison protocol is documented in `benchmarks/README.md`.
Its selected held-out wheel is 31.1% faster than original pairing on the medium
cohort and 4.8% faster on the fixed nonsplitting campaign, but 44.0% slower on
uneven inputs. Every selected wheel loses to reusable unpaired programs:
75.1%/108.7%/117.2% slower on small/medium/uneven complete factoring and 23.6%
slower on the campaign. Completion is unchanged. No promotion gate passes.

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

## B4 arithmetic kernels

`ecm.scalar_multiply` uses the confirmed selected-reduction ladder for native
integers and the readable baseline for an mpz modulus. Representation
selection happens once per scalar action. Native integers remain the default;
there is no new kernel selector or public factoring option.

Both paths preserve `a24=(A+2)/4` with the squared difference, exactly the
same canonical X:Z coordinates, scalar validation and infinity shortcuts.
The readable `point_add` and `point_double` formulas remain available as
oracles and for stage-two recurrences. No normalization/inversion is added to
production. Work charges, cancellation boundaries, saturation recovery and
checkpoint schema/backend identities are unchanged. Old and new checkpoints
resume in either engine without migration.

The final production bridge confirms a native complete-run saving of 6.79%
(paired 95% interval 5.82–13.58%) on the bounded cohort. GMP whole-ladder
fusion is retained as an experiment: its production bridge interval crosses
zero, so the readable GMP path remains in production.

The independent affine/composite/prime-power and exact readable-formula
controls exercise production as well as the frozen candidates. See the
[research/proofs](benchmarks/b4_research.md),
[source/license audit](benchmarks/b4_research_audit.md),
[fresh comparison](benchmarks/README.md#fresh-bakeoff-results-and-decision) and
[production integration](benchmarks/README.md#b4-production-integration-protocol).
The original frozen retain-baseline verdict and subsequent revised-policy
confirmation are historical records; their source controls and certified
inputs remain immutable. Experimental square, helper and normalized arms
remain reproducibility controls under `v2.benchmarks.b4_kernels`.

### A6: opt-in finite p−1 bound campaigns

```python
from v2.budget import Budget
from v2.pm1_bounded import PM1Config, factorize_pm1_bounded

config = PM1Config(bounds=((15, 15), (16, 16)))
run = factorize_pm1_bounded(17 * 1019, base=3, config=config,
                            budget=Budget(work_limit=100000))
assert run.divisor == 17
assert run.result.reconstruct() == 17 * 1019
```

This separate Python-integer API predeclares 1–64 monotone inclusive (B1,B2)
rungs for one explicit base. It stops on a valid factor, saturation, nonunit,
finite allowance or final campaign exhaustion. Existing `factorize_pm1`,
portfolio configuration, ECM campaigns and their checkpoints retain their
behavior. Repeated p−1 bases have the same p−1 smoothness structure and are
not independent ECM-like smooth-order trials.

A B1 increase applies the exact ratio M(new B1)/M(old B1), where
M(B)=lcm(1,…,B). This includes higher powers of old primes; it restarts stage
two at new B1 using the updated residue. A B2-only increase reuses checked
coverage and appends the new interval. Saturated chunks replay prime units,
and saturated stage-two batches replay terms, both under `recovery_limit`.
Saturation ends this base; increasing its bound cannot undo an identity
residue. Input and base sizes, prime workspace, chunk/batch sizes and
checkpoint output all have finite caps. The memory cap estimates owned
workspace, not process RSS.

`max_actions=N` pauses after at most N committed actions. Resume by passing
`run.checkpoint`, the same n/base/config, and an unused `Budget` containing
*total cumulative* allowances. Schema 1 binds `pm1-campaign-v1`,
`inclusive-lcm-ratio-v1` and `python-int`. No RNG is consumed. Checksums detect
accidental corruption; deterministic reconstruction verifies all saved
arithmetic and canonical numeric types before reuse, and consumes the same
action reservations in the cumulative budget. Reconstruction, context rebuilding and serialization time
are charged on resume. A small grant may be spent entirely on verification;
repeated pauses do not reset work or active-run wall/CPU usage. Paused time is
excluded. An incompatible identity is rejected. Exhaustion does not grant
new rungs, bases or ECM curves.

`PM1Run` exposes divisor, reason, cumulative work/time, verification work and
a reconstructible `FactorizationResult`. Its split pieces remain unresolved;
this API makes no primality assertion or certainty upgrade. Work units are
versioned for this API: prime segments cost segment_size+base-prime count;
compiled candidates cost one each; chunk powering/GCD costs one plus the sum
of exponent bit lengths; replay costs prime.bit_length()+1; stage-two terms
cost gap.bit_length()+1 even on a cache hit; transitions/GCDs cost one. Context
construction costs ceil(isqrt(max B2)/2). Every reservation precedes mutation.
Deadlines/cancellation are cooperative between bounded actions.

The group-independent ratio helpers and precise reuse rules are documented
in the [A6 research and contract](benchmarks/a6_pm1_research.md). They unblock
the scalar contract for later ECM work; ECM bound migration remains separate.

The separate `v2.pm1_tuning.PM1TuningConfig` is an opt-in execution
configuration for this same entry point. These A6 APIs are integrated into
mainline; their bounded correctness gate passes.
`chunk_size=64, gap_mode="recurrence"` is accepted for scoped integration
review after fresh 27-sample complete-call confirmation: 1.0–2.8% less CPU
than the frozen chunk-64 control, with positive aggregate intervals in every
size class. Bit caps and wheels remain experimental, unpromoted alternatives;
production defaults and allocation remain unchanged. See the
[complete comparison and limitations](benchmarks/README.md#a6-follow-up--bit-caps-recurrence-and-exact-wheel-execution).

```python
from v2.pm1_tuning import PM1TuningConfig

config = PM1TuningConfig(
    bounds=((2000, 20000),), chunk_size=64, gap_mode="recurrence"
)
```

It adds `chunk_bits` (zero, or 32–4096),
`gap_mode="cached"|"recurrence"`, `gap_entries` (1–256), and `wheel` (0, 30,
or 210). The prime-count chunk cap still applies; a bit cap bounds the sum of
factor bit lengths and therefore the product exponent's length. A nonzero
wheel replaces the gap executor with ordinary p−1 ± relations; it does not
change bounds or introduce Williams p+1/Lucas or ECM execution.

Paired tables retain coprime offsets and the small prime divisors of the
wheel, so exceptional primes remain covered. Records contain only eligible
primes. A two-prime trace term is the
product of their ordinary relations times a unit modulo n; singleton terms
are direct relations. Saturated batches replay each original prime under the
finite recovery limit. Table setup/inversion, plan construction, multiplications
and gap-table growth are charged before state mutation. The finite workspace
reserve includes baby/inverse tables, pending center records, replay metadata
and serialization copies; it is an owned-storage bound, not a process RSS cap.

Tuned checkpoints use `execution="pm1-tuning-v1"` and bind every configuration
field. They cannot resume as legacy campaigns or under different tuning.
Legacy `PM1Config` serialization, `pm1-campaign-v1` identity and action work
remain unchanged. Both modes reconstruct retained arithmetic and charge
verification to cumulative allowances. Increased B1 invalidates all tuned
residue tables; equal-B1 B2 extensions retain checked coverage and append only
the new interval. Defaults, RNG assignments and portfolio APIs are unchanged.

Tuned work units are separately identified by `pm1-tuning-v1`. Bit-capped
filling reserves `max(1, scanned_candidates)` before publishing a chunk;
powering and prime-unit recovery keep the legacy charges. Even-gap recurrence
keeps each gap's `gap.bit_length()+1` charge and adds one per newly retained
even power. Wheel setup reserves `2*n.bit_length()+2*D+4`, including unit
checking/inversion and both bounded power tables. A center plan with k primes
reserves `2*center.bit_length()+k+4` for its first giant, or
`2*max(1, distance_in_wheels.bit_length())+k+4` for a later giant. Collecting k
primes costs k+1; evaluating k paired/singleton records costs 2*k+1. Saturated
records replay original prime q with `q.bit_length()+1` per attempt. Finite
transitions, GCDs, context building and full checkpoint verification retain
explicit reservations. These are deterministic allowance units, not measured
CPU instructions; cache hits never reset or refund cumulative allowances.

## A6 production default promotion — 9 October 2026

Fresh bounded portfolio calls use `pm1_chunk_size=64` and
`pm1_gap_mode="recurrence"`. The shared `chunk_size=16` still governs ECM;
`pm1_chunk_size=None` inherits it. The recurrence retains at most 64 even
powers of the fixed stage-one residue. Oversized/odd gaps use the existing
64-entry exponent cache. Reserve each growth before mutation, in addition to
the legacy gap reservation. On resume, reserve/recompute even powers and
exceptional cached powers before use; this verification consumes cumulative
allowances. Conservative owned workspace grows by the table/prime-chunk
reserve, with all storage still inside the configured cap.

New portfolio snapshots with the optimized p−1 execution use version 9;
configuration identity pins both settings. Existing v2–v8 snapshots retain
legacy chunks/cache. Omitted-config library resume and CLI resume select the
saved p−1 settings automatically. An explicit custom config must match; for
an old checkpoint use `pm1_gap_mode="cached", pm1_chunk_size=None` with its
original remaining settings. Increasing a work allowance permits progress;
changing bounds, seeds, attempts or curves requires a new declared campaign.

`factorize_pm1_bounded(n)` now defaults to
`PM1TuningConfig(chunk_size=64, gap_mode="recurrence")`. Explicit `PM1Config`
continues to denote the legacy campaign executor and preserves its checkpoint
identity/work rules. Omitted-config resume of a legacy default campaign
retains that executor; custom campaigns still require the original config.
Its full deterministic resume reconstruction remains charged, and can cost
more than running a fresh final bound. No new continuation rung is allocated.

The user explicitly requested this default promotion. Fresh complete-stage
confirmation saves 9.73% CPU across nine cells; integer portfolio captures
remain inconclusive after extension. The [benchmark receipt](
benchmarks/README.md#a6-production-default-promotion--9-october-2026)
records both results. This is a user-directed default change, with the
legacy executor available for reproducibility and existing resumes.
