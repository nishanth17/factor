# Quadratic-sieve and relation engines

Exact, bounded factoring components for **PyPy implementing Python 3.11**.
This package supplies QS/MPQS controls, self-initializing quadratic sieve
(SIQS), experimental Smooth Subsum Search (SSS/SSSf), and shared relation,
filtering, GF(2) and extraction machinery.

## Engines

| Engine | Role |
| --- | --- |
| QS | Fixed-polynomial reference and bounded sieve collector |
| MPQS | Multiple independently chosen polynomials, useful as a comparison control |
| SIQS | CRT/Gray polynomial families, reused roots and shared verified relations |
| SSS | CRT collision candidates with exact smoothness testing |
| SSSf | SSS with an optional lossy candidate filter; missed candidates can affect completion |

These engines are experimental, with finite parameters and explicit unfinished
outcomes. The default portfolio uses rho, p−1 and ECM. Enable SIQS/SSS through
`PortfolioConfig`; the CLI also exposes `--method sss` and `--method sssf`.
See the [v2 guide](../README.md) for complete factoring and CLI examples.

## Pipeline

1. Build a factor base and exact modular roots, checking for setup divisors.
2. Construct polynomials with `B² ≡ h*n (mod A)` using integer arithmetic.
3. Collect candidates and verify sparse exponents of `A*F(x)`, including sign.
4. Retain full relations and bounded single-large-prime partials; combine
   matching partials with checked provenance and square corrections.
5. Filter singleton/weight-two rows and solve parity dependencies over GF(2).
6. Reconstruct the exact congruence of squares, take GCDs and validate a proper
   divisor. Trivial dependencies permit further work within the same allowance.

Every split must satisfy `1 < divisor < n` and divide n exactly. Unfinished jobs
preserve their cofactor. Relation parity alone is insufficient for verification:
full exponents, referenced atoms and original row identities remain available.

## Run a SIQS job

From the repository root, with the supported PyPy runtime:

```python
from v2.budget import Budget
from v2.qs import SIQSConfig, SIQSJob

n = 4001 * 5003
config = SIQSConfig(base_bound=200, half_width=256)
job = SIQSJob(
    n, seed=7, config=config,
    budget=Budget(work_limit=200_000_000, seconds=30, cpu_seconds=30),
)
result = job.run()
if result.divisor is not None:
    assert 1 < result.divisor < n
    assert result.divisor * result.cofactor == n
else:
    assert result.cofactor == n
print(result.reason, result.stats)
```

`QSJob` works on a supplied polynomial and signed half-open window.
`SIQSJob` manages families and shared stores; `SIQSConfig(mode="qs")` or
`mode="mpqs"` selects its reference controls. `SSSJob` and `SSSConfig` live in
`v2.qs.sss` and use the same verified postprocessing pipeline.

## Resume and limits

Use `job.run(max_blocks=...)` to pause SIQS at bounded block boundaries.
`job.checkpoint()` returns a portable checkpoint;
`SIQSJob.from_checkpoint(checkpoint, budget=..., config=...)` restores it with
charged replay and exact validation. New total allowances include prior
consumption. Explicit configuration extensions use `allow_extension=True`;
they must preserve compatible job identity. SSS exposes corresponding job
checkpoint APIs; `QSJob` supports in-memory continuation.

Base size, polynomial/family count, interval width, relation/partial storage,
trivial dependencies and workspace all have finite limits. Jobs report why they
stopped; increasing a timer alone cannot repair an exhausted local schedule.
Time/cancellation checks are cooperative and workspace caps are not process RSS
caps. Serial, thread and spawned-process SIQS workers live in `parallel.py`;
parallel execution is optional and has no universal speedup claim.

## Modules

| Files | Responsibility |
| --- | --- |
| `factor_base.py`, `polynomial.py`, `multiplier.py` | Exact setup, roots, polynomials and multiplier selection |
| `families.py`, `assignment_stream.py`, `external_square.py` | Family reuse and bounded streamed coefficient assignments |
| `reference_collector.py`, `sieve_collector.py`, `power_sieve.py` | Reference enumeration and scored candidate collection |
| `relations.py`, `smooth_batch.py`, `sss.py` | Verified relation storage and SSS candidate/smoothness work |
| `linear_algebra.py`, `extraction.py`, `pipeline.py` | Filtering, dependencies, congruences and divisor recovery |
| `siqs.py`, `capacity.py`, `checkpoint.py`, `sss_checkpoint.py` | Job orchestration, admission estimates and portable resume |
| `parallel.py` | Bounded worker coordination and central exact verification |

Public verification rechecks arithmetic and provenance. Internal reuse is
bounded and requires identical immutable inputs; caches are rebuilt on restore.
Dependency cadence and alternative assignment/scoring policies are explicit
experiments. Double-large-prime SIQS and GNFS remain roadmap work.

## Checks and evidence

```sh
make -C v2 test
make -C v2 lint
make -C v2 benchmark-phase-three-siqs WARMUP_SECONDS=3 REPETITIONS=9
make -C v2 benchmark-phase-three-sss WARMUP_SECONDS=3 REPETITIONS=9
```

The [benchmark guide](../benchmarks/README.md) summarizes measured improvements
and their limits. Keep small successful cohorts, capped failures and configured
large continuations distinct; no universal digit cutoff follows from them.
See the [roadmap](../ROADMAP.md) for remaining acceptance gates.
