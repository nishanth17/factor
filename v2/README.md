# Factor v2 — PyPy Python 3.11

Phase 1 repairs the original implementation's arithmetic and failure behavior.
The algorithms remain Python implementations using standard-library integers.
PyPy implementing Python 3.11 is the only supported project runtime.
CPython is unsupported; earlier CPython measurements remain historical evidence.
The original Python 2 code is preserved in
[`../v1/`](../v1/) and its historical audit is in [`audit/`](audit/).

## Run

From the repository root:

```sh
make -C v2 run N=626100403 SEED=7
./v2/factor.py -15
```

Install PyPy with `brew install pypy3` on this Homebrew-based Mac.
`make -C v2 runtime` shows the selected interpreter. The Makefile defaults to
`pypy3`, and the executable `factor.py` uses a PyPy shebang. Use `PYTHON` only
to select another PyPy Python 3.11 executable. Library callers should use the
same supported runtime.

Omit the number to enter it interactively. Use `--verbose` for diagnostics.
The CLI exits with status 1 when a composite remains unresolved, and status 2
for invalid input. Zero is rejected; one has an empty factorization; negative
inputs retain their sign.

## Library API

```python
from v2.factor import factorize, print_factorization

result = factorize(626100403, seed=7)
print(print_factorization(result.original, result))
assert result.reconstruct() == result.original

for factor in result.factors:
    print(factor.value, factor.exponent, factor.certainty.value)
```

`FactorizationResult` always preserves the original integer through its sign,
factor multiplicities, and `remaining` cofactors. `complete` means no unresolved
composite remains. `proven` additionally requires proven primality for every
terminal factor. Larger probable primes are labeled explicitly in formatted
output. `utils.classify_prime()` supplies `COMPOSITE`, `PROBABLE`, or `PROVEN`;
the Boolean `is_prime()` convenience function does not certify large inputs.

The deterministic Miller–Rabin domain is strictly below `2**64`: two bases
`(31, 73)` below `9,080,191`, three `(2, 7, 61)` below `4,759,123,141`,
and the seven bounded bases above those thresholds. These strict domains
are documented in [SymPy's primality implementation](https://github.com/sympy/sympy/blob/master/sympy/ntheory/primetest.py).
Above that domain, the default runs 30 random rounds. Explicit probabilistic
mode honors the requested rounds for surviving candidates; trivial small-prime
proofs and compositeness witnesses can finish sooner.

## Modules and changed contracts

| Module | Role |
| --- | --- |
| `factor.py` | Complete/partial result API, exact trial division, CLI |
| `utils.py` | Primality evidence, exact GCD/inversion/powers, batch recovery |
| `prime_sieve.py` | Wheel-6 small sieve, odd-only segmented sieve, corrected Atkin |
| `pollard_rho.py` | Brent rho with bounded independent attempts and recovery |
| `pollard_pm1.py` | Two-stage p−1 with every relation included and saturation replay |
| `ecm.py` | Modular Suyama setup, Montgomery ladder, recoverable stage-2 batches |
| `constants.py` | Conservative configurable defaults |
| `budget.py` | Shared work, wall-time, CPU-time, and cancellation limits |
| `preprocessing.py` | Exact roots, power-of-two stripping, bounded Fermat |
| `schedules.py` | Capped prime streams, reusable contexts, integer caches |
| `stage_jobs.py` | Serializable rho, p−1, and ECM candidate actions |
| `portfolio.py` | Seeded bounded dispatch, recovery, checkpoint validation |
| `qs/` | Exact reference QS/MPQS polynomials, roots, and verified relations |

The original camel-case filenames have been replaced with snake_case. Import
these modules through `v2`; do not add both version directories to `sys.path`.
Splitter functions return a **proper divisor or `None`**, including for prime
and unit inputs. `factorize()` returns a result object, replacing lists/`−1`.
`xgcd(a, b)` now matches its stated coefficient-of-a contract; use
`extended_gcd()` for the complete Bézout identity. Atkin uses an explicit
modular inverse rather than relying on the old coefficient convention.

All sieve APIs use half-open bounds: `prime_sieve(hi)` lists `p < hi` and
`segmented_sieve(lo, hi)` lists `lo <= p < hi`. Factoring B1/B2 bounds remain
inclusive and callers translate them explicitly with `bound + 1`.

## Phase 3.1 reference relation API

`v2.qs` provides exact reference collection and verification. P3.3 adds the
filtering, GF(2) and factor-extraction pipeline described below. SIQS families
and dispatcher integration follow in P3.4.

```python
from v2.budget import Budget
from v2.qs import build_factor_base, collect_block, mpqs_polynomial

budget = Budget(work_limit=2_000_000, seconds=30, cpu_seconds=30)
setup = build_factor_base(101 * 103, bound=40, budget=budget)
if setup.divisor is None:
    base = setup.factor_base
    polynomial = mpqs_polynomial(base, half_width=16, budget=budget)
    block = collect_block(polynomial, base, -32, 33, budget=budget)
    for atom in block.relations:
        print(atom.position, atom.sign, atom.exponents)
```

`build_factor_base()` GCD-checks the multiplier and base primes, returning
either a base or a proper divisor. Prime bounds are half-open (`p < bound`).
It caches complete modular square roots, including 2 and primes dividing
the multiplier. `Polynomial(n, h, A, B)` requires positive A and exact
`B*B == h*n (mod A)` before deriving C. `qs_polynomial()` uses A=1;
`mpqs_polynomial()` reproducibly chooses a nonsingular factor-base prime q,
sets A=q squared nearest the integer target, and lifts B modulo A. This
single-polynomial selection is a reference control rather than SIQS families.
`polynomial_roots()` handles linear/constant cases for primes dividing A;
an `all_positions` marker avoids allocating every residue for a zero polynomial.

Atomic relations store the sign and full sparse exponents of **A*F(x)**,
including A, with immutable polynomial/position IDs. `verify_atomic()` checks
the exact integer identity and factorization. `parity_bits()` reserves bit
zero for sign and keeps full exponents available for future extraction.
`combine_relations()` accepts explicitly supplied atoms, checks each one,
and retains square corrections for even residual multiplicities. Its result
contains a verified full relation or a validated residual GCD divisor.
`verify_combined()` requires the original atom store and rechecks its IDs,
factorizations, combined fields, and modular congruence. Missing/corrupt
provenance raises `ValueError`; IDs identify positions, not authentication.

`collect_block()` exhausts signed positions in `lo, hi)`, dividing exact
values without score thresholds. Zero values are checked before division.
Full-only collection is the default; `residual_bound` permits small proven
prime residuals for reference fixtures. There is no partial matching store.
Rejected candidates, zeros, and accepted atoms count toward `scanned`.
Caps return a checked prefix and `next_position`, the first unprocessed x;
continue a refused position under an explicitly extended allowance. This
boundary is not yet a serialized Phase 3 checkpoint.

Reference limits are 4096 bits for n, A, B, and signed positions, multiplier
at most 1,000,000, factor-base bound at most 100,000, block width and retained
atoms at most 4096, residual at most 1,000,000 (within the deterministic `r < 2**64`
domain), and at most 256 distinct atoms per combination. The default owned
workspace reservation is 8 MiB, including referenced combination provenance;
it is not an RSS guarantee. Setup/combination memory failure raises
`MemoryError`; setup/verifier work exhaustion raises `BudgetExhaustedError`.
The collector returns an explicit work/time/cancellation/storage stop reason.
Pass one shared `Budget` through setup, polynomial selection, roots,
collection, verification, and combination. Calls remain quiet.

Run the focused tests and diagnostic costs with the supported PyPy runtime:

```sh
PYTHONDONTWRITEBYTECODE=1 pypy3 -m unittest v2.tests.test_qs -v
make -C v2 benchmark-phase-three-reference WARMUP_SECONDS=3 REPETITIONS=9 \
  BENCHMARK_OUTPUT=benchmarks/phase_three_reference_UNIQUE_pypy.json
```

This runner validates independent small-window/root oracles, records cold
startup separately, and measures an unchanged complete-factorization control.
It establishes reference costs; complete SIQS promotion still needs P3.4.
The [M20 acceptance summary and
[public changelog](../CHANGELOG.md) record the passed P3.1 gates and measured limits.

## Phase 3.2 bounded single-large-prime collection

`SieveCollector` caches normalized polynomial roots and integer log bounds,
reuses working buffers, and retains checked full relations and matching
single-large-prime partials. It collects a supplied QS/MPQS polynomial;
P3.3 supplies postprocessing. SIQS families and dispatcher integration remain
P3.4 work.

```python
from v2.budget import Budget
from v2.qs import SieveCollector, SieveConfig, build_factor_base, qs_polynomial

budget = Budget(work_limit=20_000_000, seconds=5, cpu_seconds=5)
setup = build_factor_base(1009 * 1013, bound=100, budget=budget)
if setup.divisor is None:
    base = setup.factor_base
    collector = SieveCollector(
        qs_polynomial(base), base,
        config=SieveConfig(residual_bound=500), budget=budget,
    )
    run = collector.collect(-128, 129)
    print(run.reason, len(run.full_relations), len(run.combined_relations))
```

Each `collect(lo, hi)` handles a signed, half-open window of at most
1,000,000 positions using working blocks of 1–4096 positions. Metadata
chunks contain 1–1024 cached prime records; the entire factor base remains
bounded and cached. Available score backends are `list`, `bytearray` and
`array`; marking policies are `dense`, `sparse` and `bucket`. Dense marking
walks root progressions; sparse marking uses direct hits when a prime's step
is at least the block width. Bounded buckets retain a prime-hit bitset per
position. Byte marking compares bounded translation/slice updates; other
backends use loops. No cutoff from a native implementation is assumed.

Scoring is deliberately conservative. For each block, endpoint/vertex
comparisons bound `|F|`; a sign crossing sets its lower bound to zero. For
prime p, exact powers give `e_max=floor(log_p(max |F|))`. Each root hit adds
`e_max*ceil(log2(p))`, an upper bound on all of p's possible contribution,
including arbitrarily high valuations within the finite input bounds. Logs
use integer bit lengths. The candidate threshold is a lower bound on
`log2(|F|)-log2(residual_bound)`. Omitting `p < small_prime_cutoff` subtracts
its complete upper allowance from that threshold. Thus `threshold_extra=0`
cannot miss admissible values through scoring. Large integers are evaluated
at a constant number of block extrema and at selected candidate positions.

Byte scores and thresholds both saturate at 255, preserving safe coverage
while potentially increasing false candidates. Array scores use unsigned
32-bit integers; the input, base and position limits bound their maximum
below `2**32`. Positive `threshold_extra` is explicitly lossy and can miss
smooth values. Small-prime omission and conservative scores may select most
positions; they promise coverage, not efficient rejection.

M26 adds `score_policy="adaptive"`, the default. When the conservative lower
threshold is already zero and `threshold_extra=0`, every nonnegative score
passes. Root/full/resieve division therefore bypasses useless weight/score
marking on that block; bucket division still builds its required hit metadata.
`score_policy="conservative"` retains the original marking control. This
bypass preserves candidate coverage; it does not tighten the threshold.

`score_policy="candidate"` instead refines selected candidates with
`floor(log2(abs(F(x))))`, computed from the exact norm already needed for
division. It retains the complete omitted-prime allowance and clips byte
thresholds consistently. The scores still upper-bound all valuations, so
zero-extra refinement cannot discard an admissible norm. It adds marking
and per-candidate work; use complete-run measurements to choose a policy.

`division="full"` tries every factor-base prime. `"roots"` excludes primes
only using exact cached congruences; `"bucket"` uses complete hit bitsets.
`"resieve"` makes a separate prime/root pass over selected candidate values,
retaining bounded exact exponents and residuals before admission. It includes
extra scratch storage and may refuse with `memory_limit` before any position
in that block is committed. All paths recover repeated factors, including
omitted small primes, and add A's cached exponents. Scores never justify a
division early exit. A must factor entirely over the supplied base; otherwise
construction raises `ValueError`. Accepted atoms pass `verify_atomic()`;
nonunit residuals must be proven primes at most the configured bound, which
is at most 1,000,000 and strictly within `r < 2**64`. Residual GCDs return
only checked proper divisors. Combinations retain square corrections and
pass the existing provenance/congruence verifier; no inversion is needed.

An unmatched residual retains one atom. A match consumes the pending entry
and pins both source atoms with its combined relation. At the partial cap,
FIFO eviction removes only the oldest unmatched atom. Accepted combinations
and their source atoms are never evicted. Full/combined relations and all
retained atoms have separate caps, each at most 4096. A zero partial cap
drops unmatched atoms and counts the loss; eviction affects later match
yield, independently of score coverage. Duplicate retained positions cannot
be admitted twice.

Work, wall/CPU and cancellation are shared through one `Budget`, including
root inversions, residual primality, GCDs, verification and matching. Store
caps refuse before publication and return `next_position`, the first
uncommitted position. Resume on the same collector after extending the
relevant allowance. When assigning a new budget, preserve `used`,
`prior_wall` and `prior_cpu` from the previous run while extending its limits.
The returned store is cumulative and diagnostics are per call. An empty
window completes without scoring. Setup failures raise `MemoryError` or
`BudgetExhaustedError`; collection returns a stop reason. Invalid arithmetic
or provenance propagates an exception without admitting that position.

The default owned-memory reservation is 8 MiB and includes cached metadata,
scores, buckets, bounded slices, bigint/combination scratch and retained
provenance. It is conservative workspace accounting, not a process RSS cap.
Release earlier result snapshots before resuming to exclude their retained
references from this accounting. The collector is serial and in-memory;
serialized family/checkpoint state belongs to P3.4. Do not change buffer or
marking configuration after construction.

```sh
make -C v2 benchmark-phase-three-collector WARMUP_SECONDS=3 REPETITIONS=9 \
  BENCHMARK_OUTPUT=benchmarks/phase_three_collector_UNIQUE_pypy.json
```

The runner separates tuning from frozen held-out collection, compares
candidate stages with independent exhaustive norm factorization, records
missed/false candidates, divisions, matches, storage and slice allocation,
and measures cold startup plus the unchanged complete-factorization control.
These are collector measurements; they do not establish a complete QS
factorization speedup. Evidence and the acceptance decision are recorded in
the M22 summary and [public changelog](../CHANGELOG.md).

## Phase 3.3 filtering and exact factor extraction

`QSJob` runs one supplied polynomial through collection, filtering, GF(2)
elimination and checked extraction under one shared allowance. It returns a
proper split or an explicit unfinished cofactor. It is an experimental
fixed-polynomial relation engine; the production dispatcher is unchanged.

```python
from v2.budget import Budget
from v2.qs import QSJob, SieveConfig, build_factor_base, qs_polynomial

budget = Budget(work_limit=200_000_000, seconds=10, cpu_seconds=10)
setup = build_factor_base(101 * 137, bound=100, budget=budget)
if setup.divisor is None:
    base = setup.factor_base
    job = QSJob(
        qs_polynomial(base), base, -256, 513,
        config=SieveConfig(residual_bound=1), budget=budget,
    )
    result = job.run()
    print(result.reason, result.divisor, result.cofactor)
```

`prepare_relations()` checks full atomic/combined payloads and provenance
before removing exact duplicates. Equal parity between distinct relations
is retained. Prepared rows pin their checked source atoms.
`filter_matrix()` uses relation rows and sign/prime columns;
iterative singleton removal preserves the kernel. Optional `weight_two=True`
XORs the two rows constrained by a weight-two column and XORs their original
row masks. Zero rows retain their dependencies. `DependencySolver` supports
highest/lowest pivots and lifts every dependency back to prepared input rows.

`extract_dependency()` independently checks original-row parity, even signed
exponent totals and `X*X % n == Y*Y % n`. It accumulates X, Y and all residual
square corrections modulo n, then tries both GCD signs. A returned divisor
does not certify either child prime. Complete balanced benchmark fixtures
separately prove both factors and reconstruct the input. Invalid provenance
or arithmetic raises `ValueError`; it never signals successful factoring.

Matrices, provenance and dependency trials have caps of 4096 rows/masks;
column indices are below 100001. Conservative reservations include incidence,
dense fill-in, lift masks, pivots and extraction temporaries, not just input
nonzeros. They share the configured memory allowance with the collector.
The default allowance is 8 MiB; measurements explicitly use 32 MiB. These
are owned-workspace estimates, not process RSS caps.

`QSJob` collects batches of 1–4096 positions inside its finite half-open
window. Filtered usable-row excess controls solves; zero-row dependencies and
the final window receive solve attempts independently of raw relation count.
Trivial GCDs permit further collection. Budget refusal retains elimination or
extraction progress and the first uncommitted collection position. Extend
the allowance through `job.budget`, preserving used work and prior active
wall/CPU; the next `run()` resumes in memory. Storage/window exhaustion keeps
the unresolved n explicit. A successful result satisfies
`result.divisor * result.cofactor == n`. Serialized checkpoints, new families
and bounded parameter growth remain P3.4 work.

```sh
make -C v2 benchmark-phase-three-pipeline WARMUP_SECONDS=3 REPETITIONS=9 \
  BENCHMARK_OUTPUT=benchmarks/phase_three_pipeline_UNIQUE_pypy.json
```

The runner imports the owned M22 collector snapshot in memory, shares P3.3
postprocessing across comparison arms, diagnoses changes on training data,
freezes recovery/filter settings, and then evaluates fresh balanced inputs.
It preserves exact recovery utilities, complete QS samples, cold CPU/RSS and
separate diagnostic profiles. The
M26 acceptance records 151 passing
PyPy tests and lint. On 16 fresh 23–26-bit balanced inputs,
the final collector fixes reduce complete QS cohort time from 27.220 to
21.446 ms with the same frozen bucket/filter settings, a 21.2% improvement.
This does not establish a cold-start gain or larger SIQS scalability. Tiny
collector controls still trail exhaustive enumeration. Resieving and tighter
candidate scoring remain optional; root division and singleton filtering
remain API defaults. See [public changelog](../CHANGELOG.md) for decisions and limits.

## Measured regression fixes

The M9 dispatcher classifies the input before trial division and reuses
classifications within each call. Trial division refreshes an exact square-root
bound only after division shrinks the remainder. Proven paths avoid allocating
random state, and composite squares split exactly before randomized searches.
Rho and ECM reuse the dispatcher's composite classification; standalone calls
retain their primality checks. Rho accounts for evaluations per completed
chunk while retaining finite recovery limits.

The small sieve marks wheel-6 candidates with bytearray slices, preserving
half-open endpoints. PyPy uses a one-coefficient Euclidean modular inverse and
a duplicate-safe search loop selected once at import. Historical CPython
fallbacks remain in the code without a support or testing commitment.
`binary_search(value, array)` returns the first
index with an item greater than `value`; `include_equal=True` instead returns
the first item greater than or equal to it. Sorted duplicate and empty
sequences are supported. GCD and integer roots remain standard-library exact
operations. ECM point formulas reuse coordinate sums/differences.

See the M9 report for matched before/after
and emulated-v1 results, rejected experiments, and residual costs. The original
source is unchanged. These small-workload gains do not close later phase gates.

## Limits and next work

### Bounded Phase 2 API

`factorize_bounded()` adds a serial rho/p−1/ECM portfolio with one shared work,
wall-time, and CPU-time allowance across preprocessing, retries, and children.
Exact higher powers preserve multiplicities; powers of two strip in one step.
Primality witnesses, trial chunks, rho batches, exponent chunks, and stage-two
prime batches have resumable boundaries. Optional Fermat steps are bounded.
Library calls remain quiet; inspect the returned events and stage timings.

```python
from v2.budget import Budget
from v2.portfolio import PortfolioConfig, factorize_bounded

config = PortfolioConfig()
run = factorize_bounded(
    25013 * 25031, seed=7, config=config,
    budget=Budget(work_limit=1000, seconds=30, cpu_seconds=30),
)
assert run.result.reconstruct() == 25013 * 25031

resumed = factorize_bounded(
    run.result.original, config=config, checkpoint=run.checkpoint,
    budget=Budget(work_limit=2_000_000, seconds=60, cpu_seconds=60),
)
```

Work units are conservative reservations for finite actions, including a
whole trial chunk even when it finds a factor early. They are not CPU seconds
or comparable across unrelated algorithms. Resume preserves these reservations,
candidate assignments, and RNG position; wall/CPU consumption includes active
earlier runs and serialization overhead. Paused time is excluded. Deadlines
and cancellation are cooperative: a running native bigint operation finishes
before the next check. Defaults cap input size at 4096 bits and arithmetic
chunks at 256 items. Primality testing checks the budget between witnesses.

Checkpoints contain versioned configuration/schedule identity, exact modulus,
seed/sigma or base, work position, remaining allowances, pending cofactors,
certainty, and checksums. Save with JSON; different configurations, corrupt
checksums, inconsistent reconstruction, and invalid terminal evidence reject.
Checksums detect corruption and are not authentication. An exhausted global
budget can resume with explicitly increased total allowances. Once every local
candidate is exhausted, an unresolved cofactor needs a new run with different
candidate allowances; resuming the same exhausted schedule does not reset it.
`stop_after_split=True` supports factor-one work; its checkpoint can resume a
complete factorization under the same configuration.

The M13 batch loops keep rho and stage-two arithmetic in local variables, then
commit the same state at the reserved boundary. They preserve M12 work counts,
candidate assignments, and checkpoint version 2; existing checkpoints remain
compatible. This reduces measured overhead without changing portfolio bounds.

M17 confirmed the loop change on an independent corpus and froze the
[Phase 2 core baseline](audit/phase_two_m17_frozen_baseline.json). With five
seeds, nine repetitions, and matched 50 ms operation caps, fresh balanced
20-digit complete factoring improved from 70.8% to 91.1% completion. Median
total cost including unfinished outcomes fell from 39.283 to 16.362 ms;
see confirmation evidence.
These are PyPy Python 3.11 results for the declared evaluation configuration,
with two ECM curves; the library's 32-curve default and opt-in API are unchanged.
The core exit passed; broader parameter and competitor gates remain open.

The single-pass large-band screen
completed 4% of balanced 30-digit inputs and none of the balanced 40–80-digit
inputs under those operation caps. This does not establish universal size
limits. Larger balanced coverage remains Phase 3 SIQS/SSS work.

The default 8 MiB memory allowance covers conservative owned scheduling,
cursor, point/batch, result, and checkpoint reserves. `SieveContext` reuses
packed base primes and private bytearray scratch; every prime API is half-open.
Stage two retains one segment and one GCD batch, rather than a complete B2
prime list. High output values stay Python integers. `ScheduleCache` optionally
caches packed integer primes, gaps, or prime powers; curve-specific residues
and points are never shared. Enable it with `schedule_cache_bytes`; zero is
the default. Each context/cache permits one active iterator and must be closed
when abandoned. Allocation checks reject bounds that cannot fit the workspace.
Consumer-retained output, interpreter/JIT memory, and total process RSS are
outside this workspace cap. Benchmark RSS caps are measured acceptance gates.

Run or save/resume from the repository root:

```sh
make -C v2 run-bounded N=626100403 SEED=7
./v2/factor.py 626100403 --bounded --work-limit 1000 \
  --checkpoint v2/audit/my_checkpoint.json
./v2/factor.py --resume v2/audit/my_checkpoint.json \
  --work-limit 2000000 --seconds 60 --cpu-seconds 60
```

CLI exit codes retain the complete/partial/invalid conventions. `--verbose`
shows the stop reason and bounded event history. Resume requires the same
ECM curve count, memory, and Fermat options used when saving. Existing
`factorize()` calls retain the Phase 1 compatibility interface and behavior;
use the bounded API or `--bounded` for the shared scheduler. Parallel workers,
alternate wheel/bitset kernels, and disk schedules remain experiment arms.

The binary ladder is the default. `multiply_prac()` currently delegates to
that ladder; the original PRAC chain is unsafe and will be revisited in Phase 4.
Degenerate projective `(0, 0)` outputs are not accepted as valid equal points;
ECM treats saturation as recovery/retry information.

Phase 1 uses modest, untuned ECM defaults (`B1=2000`, `B2=147396`, 32 curves)
and configurable rho operation/recovery limits. This avoids the original huge
digit-derived allocations while the budgeted portfolio is being developed.
The compatibility interface has no shared deadline or resumable schedule.
Its original segmented sieve APIs still materialize full output. The bounded
API above supplies streamed scheduling, shared limits, and p−1 integration.

The old “50–60 digits within a minute” claim is withdrawn. Balanced-composite
coverage needs the later SIQS/SSS work, and no Phase 1 microbenchmark establishes
such a performance guarantee. See [public changelog](../CHANGELOG.md) for measured changes,
including regressions, and [`audit/TODOS.md`](audit/TODOS.md) for the gated plan.

GNFS is committed roadmap work. Sequence: finish the Phase 3 SIQS relation
pipeline, bring P4.3's PyPy/GMP arithmetic assessment forward, then build the
Phase 7 GNFS reference through polynomial selection, rational/algebraic
relations, filtering, dependencies and both square roots. Scale and measure
its SIQS crossover after that. Phase numbers preserve existing task IDs;
optional ECM, NumPy and parallelism experiments do not block GNFS. Neither
GNFS nor GMP is implemented yet. See
[the execution order](audit/TODOS.md#execution-order-toward-gnfs) and
[GNFS milestones](audit/TODOS.md#phase-7--implement-and-scale-gnfs).

Your C sieve versions in `../../primesieve/v2`, `v3`, and `v4` are design references
(a sibling repository, not a runtime dependency). Private buffers, exact
endpoints, shortened final segments, and workload-specific output policies
inform this version. The basic wheel-6 sieve is now measured; reusable
contexts, rolling strikes, wheel-30/pre-sieve packing, and parallel candidate
searches are separately measured candidates with their own experiment gates.

## Validate and benchmark

From the repository root:

```sh
make -C v2 test
make -C v2 validate
make -C v2 benchmark
```

Validation and benchmark targets save timestamped artifacts in v2. Benchmarks
use nine repetitions and at least three seconds of validated workload warmup
per case, including JIT compilation. Inspect sample spread before claiming
steady-state gains. `WARMUP_SECONDS`, `REPETITIONS`, `VALIDATION_OUTPUT`, and
`BENCHMARK_OUTPUT` are configurable Make variables. Compare implementations on
the same PyPy/corpus/protocol; record cold startup separately. Use focused
tests while iterating and full validation at acceptance milestones.

The Phase 2 runner uses the independently certified
`benchmarks/phase_two_complete_corpus.json`, with separate training and held-out
inputs and five fixed seeds. Complete factorization and factor-one runs use
matching wall/CPU/RSS caps. Warmups verify answers against the hidden oracle;
comparison order rotates between repetitions. Cold samples cover one input
per band and mode under every requested seed, including worker startup,
serialization, and shutdown. Timeouts and exhausted runs stay in summaries.

```sh
make -C v2 benchmark-phase-two \
  BENCHMARK_OUTPUT=benchmarks/phase_two_NEW_pypy.json
PYTHONDONTWRITEBYTECODE=1 pypy3 -u -m v2.benchmarks.phase_two \
  --bands powers,random_small,primes,close_small \
  --seeds 5 --repetitions 9 --warmup-seconds 3 \
  --output v2/benchmarks/phase_two_NEW_small_pypy.json
```

Use `--corpus PATH` for a separate confirmation corpus. For broad comparisons,
`--warmup-scope corpus` warms every selected input and output mode rather than
one representative per band. It still validates hidden oracle results during
warmup. Inspect repetition medians for residual JIT drift before accepting a
timing claim; warmup can exceed its requested duration to finish a full pass.

Use fresh output names and run benchmarks sequentially. The
default sample caps are 50 ms wall/CPU and 256 MiB process RSS. RSS is a
measured acceptance gate; the library's workspace cap excludes interpreter
and JIT storage. Balanced large inputs can be infeasible under these caps.
Pinned competitors have not been executed. Experimental cutoffs, tiers,
cache/rolling/wheel kernels, and parallel workers require further evidence
before changing production defaults.

To compare the retained batch implementation with the verified M12 source
snapshot on identical assignments, use `--engines m12,bounded`. The control
loads the historical stage-job module in memory, verifies its dependency
hashes, and runs the same current portfolio interface. Cold control startup
includes snapshot verification/compilation; interpret cold costs separately.
The matched M13 comparison uses ten seconds of representative warmup, five
seeds, and nine repetitions. Large exploratory cohorts with many distinct
paths may need longer warmup than three seconds.

Use PyPy Python 3.11 with `lib2to3` for an emulated original-code comparison.
Match the interpreter for old/new comparisons:

```sh
pypy3 -m v2.benchmarks.regressions --legacy --warmup-seconds 3 \
  --repetitions 15 --output v2/benchmarks/regressions_NEW_pypy.json
```

Replace `NEW` with a fresh milestone identifier; preserve previous captures.
The regression runner verifies the frozen M8 source hashes and loads it in
memory. It validates all timed answers and rotates candidate order.

The historical audit scripts intentionally continue to load `v1/`. Their
captured JSON remains unchanged. New acceptance tests directly import `v2/`.

Benchmark runners, corpora and required immutable baselines are committed.
Run captures, gzip archives, profiler output and detailed validation
records stay local and are Git-ignored. See the
[benchmark guide](benchmarks/README.md) for commands and selected results.

PEP 8 verification tools are development-only dependencies:

```sh
pypy3 -m venv v2/.venv
v2/.venv/bin/python -m pip install ruff==0.14.14 pycodestyle==2.14.0
make -C v2 lint
```

The lint configuration excludes historical audit sources so provenance
artifacts retain their original contents. Existing copied production code,
new tests, and benchmark code are all formatted, spaced, and commented.
