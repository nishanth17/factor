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
filtering, GF(2) and factor-extraction pipeline described below. P3.4 provides
complete bounded SIQS jobs, checkpoints and optional dispatcher integration.

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
atoms at most 4096, residual at most 1,000,000,000,000 (within `r < 2**64`
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
It establishes reference costs. The bounded P3.4 comparisons are below;
broader portfolio promotion remains separate.
The M20 acceptance summary and
[public changelog](../CHANGELOG.md) record the passed P3.1 gates and measured limits.

## Phase 3.2 bounded single-large-prime collection

`SieveCollector` caches normalized polynomial roots and integer log bounds,
reuses working buffers, and retains checked full relations and matching
single-large-prime partials. It collects a supplied QS/MPQS polynomial;
P3.3 supplies postprocessing. P3.4 adds shared
family jobs and optional dispatcher integration below.

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

`score_policy="powers"` is an optional tighter sieve. Exact Hensel lifting
marks each prime power separately; a root class with no position in the
current block needs no higher lifts. Singular branches are capped at 64
roots, then receive a proved conservative allowance. This can admit extra
candidates but preserves every admissible norm at zero extra threshold.
Candidate refinement and byte saturation obey the same coverage contract.
Lifting work and both bounded root lists are reserved before use; configure
sufficient memory (64 MiB is suitable for the standalone test controls).
Independent checks cover signed windows, repeated prime powers, singular
roots, all division paths and score buffers. The default remains adaptive.

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
is at most 1,000,000,000,000 and strictly within `r < 2**64`. Residual GCDs return
only checked proper divisors. Combinations retain square corrections and
pass the existing provenance/congruence verifier; no inversion is needed.

An unmatched residual retains one atom. A match consumes the pending entry
and pins both source atoms with its combined relation. At the partial cap,
FIFO eviction removes only the oldest unmatched atom. Accepted combinations
and their source atoms are never evicted. Full/combined relations and all
retained atoms have separate caps, each at most 65,536. Defaults remain
small; increasing a count cap also requires sufficient owned memory. A zero partial cap
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
serialized full-job state is provided by SIQSJob. Do not change buffer or
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
row masks. Column incidence uses bounded row-index bitsets and updates only
affected columns after removals/merges. Singleton queues and a bounded heap
preserve the previous row/mask order, including serialized solver prefixes.
Zero rows retain their dependencies. `DependencySolver` supports
highest/lowest pivots and lifts every dependency back to prepared input rows.

`extract_dependency()` independently checks original-row parity, even signed
exponent totals and `X*X % n == Y*Y % n`. It accumulates X, Y and all residual
square corrections modulo n, then tries both GCD signs. A returned divisor
does not certify either child prime. Complete balanced benchmark fixtures
separately prove both factors and reconstruct the input. Invalid provenance
or arithmetic raises `ValueError`; it never signals successful factoring.

Matrices, provenance and dependency trials have caps of 65,536 rows/masks;
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
`result.divisor * result.cofactor == n`. Full relation-job checkpoints and
bounded width recovery are supplied by SIQSJob below.

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

P3.3 follow-up audit (M28) fixes storage-cap completion: relation, atom or
collector-memory exhaustion triggers one final extraction attempt from the
retained checked rows. Budget refusal preserves that attempt for resume;
unchanged unsuccessful storage stops return without repeating work. Combined
memory reservations include live collector buffers, preparation scratch and
matrix workspace, counting shared atoms/base once. Matched partials reserve their sparse exponent
lists as retained storage; dense combination scratch is temporary and
shares the live-store allowance. Checkpoint reconstruction applies the same
accounting and refuses before exceeding that allowance. A matrix memory refusal
keeps the unresolved cofactor explicit. The audit passes 157 PyPy tests and
lint, plus 214 root, 648 collector and 320 matrix-oracle comparisons.
See [benchmark scope and correction costs](benchmarks/README.md).

## Phase 3.4 bounded SIQS and complete checkpoints

`SIQSJob` shares verified full relations and single-prime partial matches
across exact CRT/Gray polynomial families. One budget covers multiplier/base
setup, roots, marking, exact recovery, filtering, elimination and extraction.
`mode="qs"` retains the fixed-polynomial control; `mode="mpqs"` visits a finite
schedule of lifted q-squared polynomials. The production portfolio uses SIQS
when explicitly configured after its ECM schedule.

```python
from v2.budget import Budget
from v2.qs import SIQSConfig, SIQSJob

config = SIQSConfig(base_bound=200, half_width=256)
job = SIQSJob(
    4001 * 5003, seed=7, config=config,
    budget=Budget(work_limit=200_000_000, seconds=10, cpu_seconds=10),
)
paused = job.run(max_blocks=1)
checkpoint = job.checkpoint()
resumed = SIQSJob.from_checkpoint(
    checkpoint, config=config,
    budget=Budget(work_limit=200_000_000, seconds=20, cpu_seconds=20),
)
result = resumed.run()
assert result.divisor * result.cofactor == 4001 * 5003
```

`run()` returns a proper split or the explicit unresolved n. `max_blocks`
pauses at collection boundaries. Work/time/cancellation refusal keeps the
first unpublished position and pending elimination/extraction state. Successful
and terminal local outcomes are idempotent. Limits include 1–8 A factors,
1–64 families, at most 128 polynomials per family, bounded block/row/atom
counts and an explicit total owned-memory allowance (32 MiB by default).
Reservations describe owned workspace, rather than process RSS.

`max_stalled` and `max_trivial` bound unproductive polynomials and repeated
trivial dependencies. `growth_steps` permits at most four finite width epochs,
bounded by `max_half_width`. Growth preserves the exact factor base and
verified store, charges assignments/roots again and retains consumed resources.
Base enlargement and disk spill are disabled; the disk allowance is zero.
Exhausted schedules, store limits or memory refusals retain an unfinished
cofactor. `diverse=False` and `shared_relations=False` expose simple/fresh-store
controls. Polynomial reflection/translation keys and checked atomic IDs track
duplicates; parity equality between distinct relations remains valid data.

`multiplier=1` is the fixed control; zero selects a finite integer
Knuth–Schroeppel-style score. Exact small-prime residues and the class modulo
8 contribute to the score, with a half-log multiplier-size penalty.
Integer-scaled logarithms apply only to bounded small integers. A multiplier
GCD can return a proper split. This independently implemented heuristic has
its own matched controls; it imports no native digit table or optimality claim.

Full checkpoints retain configuration, seed/assignment identity, width epoch,
family/Gray index, block position, checked-store identity, pending polynomial,
compact solver/extraction progress and consumed work/wall/CPU. Polynomial
provenance is stored once per polynomial; atoms and matches reference it.
New solver-prefix digests stream binary rows/masks as hexadecimal text,
with an explicit encoding tag. Earlier decimal digests remain readable; wide
masks avoid Python's decimal integer-string limit.
The matrix is reconstructed and the verified elimination prefix replayed;
completed trivial extraction trials are replayed too. Reconstruction costs
are charged to the resumed allowance. This bounds checkpoint size without
serializing duplicate matrices. A fresh resume Budget automatically receives
prior consumption; a used Budget must already retain those resources.

Checkpoint envelopes have a version, integrity markers and a default 1 MiB
cap, configurable from 4 KiB to 16 MiB. Larger envelopes reserve their
encoding workspace before a job begins; the default remains 1 MiB.
Retained roots/store/matrix and encoding scratch share the overall memory allowance. Incompatible configurations,
corrupt provenance, invalid prefixes and resource resets are rejected.
If a deliberately small checkpoint cap is insufficient, `checkpoint()` raises
`MemoryError` and leaves the checked in-memory state resumable. Integrity
markers detect corruption; they are not authentication signatures.

`PolynomialFamily`, `family_assignments()` and independently certified
`SieveCollector(..., precomputed_roots=...)` remain available as lower-level
interfaces. Gray transitions use the actual recentered B difference; 2 and
primes dividing A use exact complete branches. Family-only checkpoints remain
separate from full-job checkpoints. Caller-retained family state must share
memory with collection/postprocessing; `SIQSJob` performs that partition.

```python
from v2.portfolio import PortfolioConfig, factorize_bounded

portfolio = PortfolioConfig(
    memory_bytes=64 * 1024 * 1024,
    siqs=SIQSConfig(base_bound=200, half_width=256),
)
run = factorize_bounded(
    4001 * 5003, seed=7, config=portfolio,
    budget=Budget(work_limit=200_000_000, seconds=20, cpu_seconds=20),
)
assert run.result.reconstruct() == 4001 * 5003
```

The portfolio reserves its own retained state/context alongside SIQS and
keeps one allowance across methods and recursive children. It includes full
SIQS progress in its normal checkpoint. Version 3 adds SIQS state; version 4
adds opt-in SSS state. Version-2/3 checkpoints remain readable and upgrade on
save when their configuration matches. Runtime probable
prime labels remain distinct from proven labels and external corpus proofs.

The bounded P3.4 implementation and its declared large-number evaluation are
complete at M31. A repaired, explicitly configured SIQS job factors one
balanced 50-digit semiprime from an empty store in 1,183.355 seconds
(19 min 43 s); a separately resumed job preserves cumulative consumption.
The 30–80-digit comparison and varied-input study retain censored outcomes
and finite schedule/storage stops. The larger 60–80-digit bands establish
bounded exploration, with practical complete-factor scaling still open.
ECM remains the automatic default and SIQS is opt-in; no broad crossover or
general large-number timing guarantee is established. The current checkout
passes 230 PyPy tests and lint, including exact roots/store/matrix oracles,
cooperative filter cancellation, wide-mask and legacy checkpoint tests,
charged resume and recursive dispatch. See the
[completed results and limitations](benchmarks/README.md#completed-larger-evaluation-and-filtering-repair-m31-4-october-2026).
Earlier M29/M30 family/small-input measurements retain their historical scope.

## Phase 3.5 experimental Smooth Subsum Search

`SSSJob` is an optional serial SSS/SSSf challenger. It constructs signed
CRT/collision candidates for the exact A=1 polynomial, detects smooth parts
with capped product/remainder trees, recovers every prime exponent and uses
the same checked single-prime store, filtering, dependency solver and modular
extraction as QS/SIQS. It is available through an explicit CLI method or an
optional portfolio configuration. `auto` dispatch keeps its existing default.

```sh
pypy3 -m v2.factor 10002200057 --method sss --sss-base-bound 400 --seed 7
pypy3 -m v2.factor 10002200057 --method sssf --sss-base-bound 400 --seed 7
```

`--method sss` / `sssf` implies bounded execution. After sign/twos, primality,
trial division, exact powers and optional Fermat preprocessing, the selected
collector replaces rho/p−1/ECM. It returns full recursive factorization with
normal certainty labels and unresolved cofactors. `--sss-base-bound` defaults
to 1000 and `--sss-rounds` to 256. Selected methods default to 200 million work
units, 30 wall/CPU seconds and 80 MiB total owned memory; existing `auto`
limits stay unchanged. Explicit budget arguments always take precedence.
The CLI reserves 16 MiB for portfolio/context coexistence and gives the rest
to SSS. SSSf's two-stage processing has no lossy cutoff by default; a positive
cutoff remains an API configuration choice.

Save and resume with the same method, base, rounds, memory and Fermat options:

```sh
pypy3 -m v2.factor 10002200057 --method sss --sss-base-bound 400 \
  --work-limit 100000 --checkpoint v2/audit/sss_checkpoint_LOCAL.json
pypy3 -m v2.factor --method sss --sss-base-bound 400 \
  --resume v2/audit/sss_checkpoint_LOCAL.json --work-limit 200000000
```

For an SSS fallback after the existing bounded rho/p−1/ECM schedule, use:

```python
from v2.budget import Budget
from v2.portfolio import PortfolioConfig, factorize_bounded
from v2.qs.sss import SSSConfig

config = PortfolioConfig(memory_bytes=80 * 1024 * 1024,
                         sss=SSSConfig(base_bound=400))
run = factorize_bounded(10002200057, seed=7, config=config,
                       budget=Budget(work_limit=200_000_000))
```

`sss=None` is the default. Choose either `sss` or `siqs` as the relation
fallback. Set rho/p−1 attempt counts to zero and `ecm_tiers=()` to reach SSS
directly after preprocessing, as the selected CLI methods do. All stages and
recursive children share the same allowance. Opt-in availability is separate
from automatic performance promotion; larger crossover evidence remains open.

```python
from v2.budget import Budget
from v2.qs.sss import SSSConfig, SSSJob

budget = Budget(work_limit=200_000_000, seconds=10, cpu_seconds=10)
job = SSSJob(4001 * 4003, seed=7,
             config=SSSConfig(base_bound=400), budget=budget)
result = job.run(batch_limit=1)
if result.reason == "paused":
    result = job.run()
assert result.divisor is None or result.divisor * result.cofactor == job.n
```

The input is an odd positive integer with at most 4096 bits; callers retain
normal preprocessing and primality classification. `divisor` is either a
validated proper split or `None`, never a primality claim. Unfinished results
retain the complete input as `cofactor`. `next_position` counts completed
search assignments, rather than polynomial positions. The unfinished
assignment and its first uncommitted candidate remain in memory. To resume
a budget refusal, extend that same `Budget` in place; replacing it is rejected.
Consumption never resets. `checkpoint()` and
`SSSJob.from_checkpoint(snapshot, budget=...)` serialize the checked store,
completed assignment index, pending candidate cursor and solver/extraction
prefix. Resume regenerates the seeded pending assignment, verifies every
retained relation and rebuilds the matrix/prefix with charged work. It retains
consumed work/wall/CPU and requires matching configuration. `checkpoint_bytes`
defaults to 256 KiB and is capped at 1 MiB; checkpoint overflow raises
`MemoryError` while leaving the live job intact. Encoding/decoding coexistence
is reserved alongside the collector. An unchanged setup-memory or accepted-
store refusal is idempotent; other setup refusals
may repeat charged private work. Recreate the job to change its parameters.

`SSSConfig` bounds base size, assignments, selection size, collision count,
candidate batches, tree nodes/bits and owned memory. Its `collector` controls
residual bounds, full/partial/atom storage and deterministic eviction. A
candidate or tree overflow refuses the entire unpublished assignment; a
full store triggers the shared final extraction attempt. Setup, collision
counters, candidate/tree coexistence and stored provenance are reserved
together. The workspace reservation is distinct from observed process RSS.

`mode="sssf"` adds two-stage smoothness processing. `filter_divisor` controls
the first subset; positive `filter_bound` intentionally rejects candidates
whose remaining part is at least that bound, except values already smooth
over the subset. Zero disables this yield-losing cutoff. Both modes detect
all base-prime powers; the upstream finite prime-power shortcut is not used.
Exact division and the independent relation verifier still decide admission.
Assignments use local seeds and sorted indices; library calls stay quiet.

The adapter needs no third-party arithmetic packages. The separate,
hash-checked upstream reproduction runner uses optional SymPy 1.14.0,
mpmath 1.3.0 and gmpy2 2.3.1 in the project-local PyPy development environment.
It preserves the upstream source and separately records its floating-point
setup, raw-prime-count parameter table, global random seed and matrix helper.
Those settings are not the adapter's defaults. See the
[comparison commands and decisions](benchmarks/README.md).

## Phase 3.6 experimental SIQS workers

`v2.qs.parallel.ParallelSIQSJob` evaluates independent SIQS polynomial
families in serial, threads or spawned PyPy processes. It is a standalone
experimental API. The factoring CLI and automatic portfolio retain serial
collection. Use PyPy implementing Python 3.11 and a guarded entry point when
spawning processes:

```python
from v2.budget import Budget
from v2.qs.parallel import CollectionPool, ParallelConfig, ParallelSIQSJob


def main():
    allowance = Budget(work_limit=200_000_000, seconds=10, cpu_seconds=10)
    config = ParallelConfig(base_bound=200, family_count=16)
    with CollectionPool("process", 2) as pool:
        job = ParallelSIQSJob(4001 * 5003, seed=7,
                              config=config, budget=allowance)
        result = job.run(pool=pool)
    assert result.divisor is None or result.divisor * result.cofactor == job.n


if __name__ == "__main__":
    main()
```

`CollectionPool` accepts `serial` with one worker or `thread`/`process`
with one to four workers. It can be reused sequentially across jobs and must
be closed. Threads run under PyPy's GIL; no native bigint backend is assumed
to release it. The fixed family schedule depends on the input, configuration
and seed, independently of worker count. Completed batches merge in assignment
order, with exact re-verification and central single-prime matching, filtering,
GF(2) solving and extraction. Only verified atomic relations cross the worker
boundary. The result is a proper split or the original unresolved cofactor;
this API does not classify or recursively factor the split's children.

The parent charges each finite `assignment_work` lease before submission,
then refunds only reported unspent work. Cancelled private work is charged.
Aggregate CPU includes parent CPU, live worker publications and a final
post-transfer snapshot; threads are counted through the parent process.
Workers use the parent's remaining wall deadline. Limits are cooperative:
CPU publications are throttled to one millisecond and parent waits poll every
ten milliseconds. Running atomic integer operations and worker startup before
the first CPU publication can overshoot an allowance.
`ParallelConfig.poll_interval` defaults to 64 and accepts 1–64. Work
reservations are checked on every action; the first action, every polling
boundary, and publication of a completed batch check clocks/cancellation.
At most 63 additional bounded atomic actions can run between external polls.
This is an action bound, not a hard millisecond deadline: native integer work,
the one-millisecond CPU publication throttle, the ten-millisecond parent
wait, and process startup also contribute to cooperative overshoot. Setting
the interval to one provides the strict-polling comparison control.
Pool creation/teardown outside `run` is a caller lifecycle cost; the benchmark
reports cold startup and shutdown separately.

`memory_bytes` caps conservative aggregate owned reservations: parent storage,
worker workspace, duplicated bases and three bounded result/IPC copies per
slot. `parent_memory_bytes` and `worker_memory_bytes` also impose local caps.
This is separate from process RSS, whose observed high-water marks include
the interpreter/JIT and allocator. Result batches are capped by
`max_batch_atoms`; cap refusal publishes no incomplete worker prefix. There
are no spill files or nested backend threads.

`ParallelConfig.batch_width=0` preserves complete-family assignments.
Widths 1–4096 instead publish contiguous position blocks within each Gray
polynomial. Assignment IDs encode family, Gray index and block in that order;
the final block can be shorter. Smaller batches allow earlier extraction and
avoid overflowing a whole-family atom cap, at the cost of repeated setup and
transport. Reservations use the smaller of the atom cap and the maximum
number of positions in an assignment. Every returned atom is still checked
against its exact polynomial/window and independently verified centrally.
The per-atom exponent reservation uses an integer upper bound on the
polynomial norm and the product of the smallest factor-base primes; a relation
cannot contain the entire base when that product exceeds its norm.

`run(max_assignments=1)` pauses after one committed assignment and drains workers
before returning. In-memory resume extends the original `Budget` in place.
Completed pending batches and their admission cursor survive interruption;
incomplete private assignments replay from their beginning with the same ID and
newly charged work. A worker's `work_limit` refusal can require a larger
`assignment_work` lease in addition to a larger aggregate allowance.
`job.checkpoint()` returns byte-capped JSON-compatible state.
`ParallelSIQSJob.from_checkpoint(checkpoint, budget=new_allowance)` restores
cumulative work/wall/CPU, verifies retained and pending provenance, and charges
setup/store reconstruction. Worker count may change after restoration.
Solver caches rebuild from verified rows; they are not serialized. Checkpoint
creation requires quiescent workers. Checksums detect corruption rather than
authenticate history supplied by another party.
Worker checkpoints now use version 2. Version-1 complete-family checkpoints
remain readable with `batch_width=0`; new checkpoints bind chunk settings and
validate admitted prefixes against their saved pending atoms. Changing batch
width changes assignment identity and requires a new job.
In-memory resume rejects changes to assignment geometry or storage settings.
It permits changes to the work lease and polling interval, and a reduction of
the checkpoint byte cap. Drained jobs release pool adapters held by paused
solvers and extractors before returning.

`run(fixed_work=True)` defers extraction until all scheduled families commit
for the fixed-work experiment, recording direct residual splits while still
scanning the complete schedule. `result.stats["schedule_complete"]` distinguishes
complete scheduled work from an interrupted/refused attempt. Preserve this
flag on resume. Normal runs try
extraction after each assignment and cancel remaining work after a checked split.
Finite family/window/store limits may exhaust without a split. Measurements
and the adoption decision are in the [benchmark guide](benchmarks/README.md).

## Measured regression fixes

The P3.6.1 audit adds exact small-prime screens before higher-power roots and
an `isqrt` square fast path. `SieveContext.prime_segment(lo, hi)` materializes
one half-open range of at most twice its segment size, sharing the existing
single-consumer buffer. Budgeted prime cursors use this bounded path; prime
streams and cached schedules keep their existing interfaces.

Bucket division visits hit primes and the support of A; resieving recovers
exponents only at actual hits. Prime-power marking shares its first pass with
hit collection and caches bounded derivative inverses per polynomial.
Skipped small primes contribute their exact valuations to candidate refinement,
preserving the conservative acceptance test. Sparse matrix column labels are
compacted when fewer than half the labeled positions are used, with original
rows retained for dependency verification and included in memory reservations.
`filter_matrix().stats` adds `working_columns`; `input_columns` keeps the
original highest-label width. Extraction visits selected dependency bits.

The common QS pipeline reserves at most 2 MiB for verified preparation
records, bounded further by the relation cap. Reuse requires the identical
immutable base and relation objects and, for combined rows, the identical
referenced atoms still present in the current store. New or changed inputs
receive full exact verification. The collector owns the cache reservation
alongside matrix scratch and clears it before releasing pinned relations;
checkpoints rebuild it. Public `prepare_relations()` continues to verify all
inputs. Cache saturation falls back to full verification.

Native QS/SIQS/SSS collection and solving use the same maximum 64-action
polling interval as worker execution. Every work reservation remains exact;
stage transitions and publication force external-limit checks. Calls restore
the caller's original budget object on both success and failure. These
cooperative checks retain the atomic-operation and startup caveats above.

SSS collision counters count distinct primes directly and reserve unchanged
logical work in chunks of at most 64. Candidate order and seeded assignments
match the frozen control. SSSf removes the first subset completely before
testing the disjoint remaining base, avoiding repeated smoothness work. Its
optional lossy filter and experimental status remain unchanged. See the
[performance audit](benchmarks/README.md) for matched controls and decisions.

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
candidate assignments, and the then-current version-2 checkpoints. M30 adds
version 3 with compatible reading of those earlier checkpoints. This reduces
measured overhead without changing portfolio bounds.

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
limits. M30 also reports capped SIQS results above; broader balanced
coverage remains unestablished.

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

The old “50–60 digits within a minute” claim is withdrawn. M31 records one
balanced 50-digit success at 19 min 43 s; it does not
establish a general one-minute promise, and no Phase 1 microbenchmark
establishes that guarantee. See [public changelog](../CHANGELOG.md) for measured changes,
including regressions, and [`audit/TODOS.md`](audit/TODOS.md) for the gated plan.

GNFS is committed roadmap work. With the bounded P3.4 SIQS pipeline complete,
bring P4.3's PyPy/GMP arithmetic assessment forward, then build the
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

### P3.8 R1: bounded capacity and extendable SIQS jobs

The existing finite reference schedules remain the default. Opt-in
`SIQSConfig(assignment_policy="nearest")` streams unique A assignments from
an integer cursor; `"flyer"` chooses a final prime by the actual product's
distance from the integer target. Flyer primes come from outside the core
pool, so each assignment has a unique core without retaining a search history.
Both policies support 1–32 A factors and a pool of at most 4,096 primes.
`family_count` is a cumulative allowance (at most `2**31 - 1`), while only one
family and its root cache are resident. `polynomials_per_family=0` visits the
whole Gray family; a positive value limits Gray reuse independently of the
number of assignments. An s-prime family has `2**(s-1)` possible polynomials.
Streaming jobs require a fixed width, `growth_steps=0`, and `diverse=True`.
Their width and factor base cannot change on resume. The actual assignment
count and total combination space are reported separately from the requested
quota. Exhausting the pool's combination space reports
`assignment_space_exhausted`; a larger quota cannot create new subsets.

The factor-base implementation ceiling is 1,000,000. A job still has an
explicit finite bound, work/time allowances and memory checks before allocation.
Streaming jobs permit up to 4 GiB of owned workspace and 64 MiB of encoded
checkpoints. Checkpoint coexistence reserves eight times the configured byte
cap; raising that cap reduces the workspace available for collection. These
are implementation ceilings, not recommended parameters or process RSS limits.
The reference assignment API retains its eight-factor/64-family restrictions.

`SIQSConfig(mode="mpqs", external_coefficients=True)` searches coefficients
q near the square root of the integer A target, independently of the factor
base. The bounded search (`coefficient_trials`, default 4,096) checks q congruent
to 3 modulo 4, verifies a root, performs exact Hensel lifting and checks the
resulting polynomial identity. Its cursor advances only with a completed
polynomial/root preparation. A failed search reports `coefficient_limit`; extending `coefficient_trials`
replays that bounded search from the same cursor and retains prior charges.
Small coefficient primes are proven by the deterministic classifier; larger
ones remain labelled probable primes. Accepted root/inverse identities and
all relations are verified independently of that label. Coefficient certainty
never changes the certainty of an output factor.

`Polynomial(..., square_coefficient=q)` represents an explicit known square
factor of A. The constructor requires `q*q` to divide A. For an atomic relation,
`U² - h*n = sign * product(p**e) * residual * square_coefficient**2` exactly.
`Polynomial.supported_a` is A after removing that square. Collection recovers
all factor-base exponents of `supported_a * F(x)`; parity ignores the known
square and extraction multiplies its root. Combined relations and checkpoints
preserve the correction. The default correction is 1 and legacy polynomial
identities/checkpoints remain readable. No primality assumption is needed to
validate this representation.

Resume with the same configuration works as before. An explicit monotone
extension is available through `SIQSJob.from_checkpoint(...,
allow_extension=True, config=extended_config)`. It may increase streamed
`family_count`, external `coefficient_trials`, `max_stalled`, `max_trivial`,
collector atom/partial/relation
caps, memory and checkpoint space. Every other configuration field must match;
limits cannot decrease, and checkpoint/cache growth cannot reduce live
workspace. Increasing a relation cap below 1,024 also grows the reserved
verification cache (2,048 bytes per relation, up to 2 MiB); increase the memory
allowance to cover that growth.
The checked relation store, assignment/Gray/block position and consumed work,
wall and CPU resources are retained. Rebuilding roots, verifying restored
relations and replaying solver state consume the extended allowance. Raising a
timer alone does not enlarge a search schedule. The containing portfolio's
configuration-match contract is unchanged; this extension API is for a direct
SIQS job.

```python
from dataclasses import replace
from v2.budget import Budget
from v2.qs import SIQSConfig, SIQSJob

config = SIQSConfig(assignment_policy="nearest", family_count=128)
job = SIQSJob(n, config=config, budget=Budget(work_limit=10**10))
result = job.run()
checkpoint = job.checkpoint()
extended = replace(config, family_count=1024, max_stalled=256)
resumed = SIQSJob.from_checkpoint(
    checkpoint,
    config=extended,
    allow_extension=True,
    budget=Budget(work_limit=2 * 10**10, seconds=60, cpu_seconds=60),
)
result = resumed.run()
```

`v2.qs.capacity.capacity_report(base, config)` reports actual base cardinality,
exact attainable A-product bounds, search quotas and separate base, metadata
cache and matrix reservations. A target inside the product envelope does not prove
that a particular selected pool reaches it; observed A bounds appear in job
statistics. `matrix_alone_fits` is only a necessary capacity check: collector,
provenance and matrix objects coexist. The filter's initial incidence work now
scales with nonzeros and bounded 64-bit-word operations, while its conservative
fill-in/provenance storage reservation is unchanged. Reported work units are
algorithmic allowances, not measured CPU instructions.
