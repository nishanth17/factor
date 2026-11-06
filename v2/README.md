# Factor v2

Factor is an integer-factorization library and command-line tool for **PyPy
implementing Python 3.11**. It combines small-factor searches with relation
engines, validates every split and recursively factors the remaining pieces.
Results preserve multiplicities, signs and unresolved cofactors, and distinguish
probable primes from proven primes.

The default Python-integer backend uses only the standard library. The active
implementation lives here; [v1](https://github.com/nishanth17/factor/blob/23458a52e409add326c9fd63cf1f322c6e605569/v2/history/V1_BASELINE.md) is the preserved
Python 2 baseline.

## Setup and first run

With PyPy implementing Python 3.11 installed, run from the repository root:

```sh
make -C v2 runtime
pypy3 -m v2.factor 626100403 --bounded --seed 7 \
  --work-limit 2000000 --seconds 30 --cpu-seconds 30
pypy3 -m v2.factor --help
```

Omit the number for an interactive prompt. Select a relation engine explicitly
with `--method qs|mpqs|siqs|sss|sssf`, or enable portfolio SIQS fallback with
`--siqs`. Explicit budgets can return a partial factorization with unfinished
cofactors; a digit count alone does not guarantee completion.

The [CLI guide](docs/CLI.md) describes methods, output and resource flags.
[Current defaults](STATUS.md#current-automatic-defaults) document automatic
routing; explicit settings and saved configurations take precedence.

## Algorithms available

| Area | Implementation |
| --- | --- |
| Preprocessing and certainty | Trial division, exact perfect-power detection, integer roots, prime sieves and primality classification. |
| Pollard rho | Seeded Brent cycle search with batched GCDs, retries and saturation recovery. |
| Pollard p−1 | Two-stage smoothness search with finite bounds and reusable prime-gap arithmetic. |
| Elliptic-curve method | Two-stage Montgomery ECM with Suyama curves, ladder/verified-chain scalar multiplication and optional paired continuation. |
| Quadratic-sieve family | QS, multiple-polynomial QS and self-initializing QS; CRT/Gray root reuse, prime-power scores, sparse hit queues and single/double-large-prime relations. |
| Experimental relation search | SSS/SSSf use CRT collision candidates and exact smooth-part recovery; optional collector policies extend intake, compact storage and graph processing. |
| Shared relation processing | Checked provenance, graph cycles, duplicate/singleton/weight-two filtering, GF(2) dependencies and congruence-of-squares extraction. |

These methods have different selection and resource contracts. Experimental
collectors and family workers are explicit options; GNFS is roadmap work.
The package guides below describe the available techniques and configuration.

## Python usage

Use the bounded portfolio when stages and recursive children must share one
work/time allowance:

```python
from v2.execution.budget import Budget
from v2.portfolio import factorize_bounded

n = 626100403
run = factorize_bounded(
    n,
    seed=7,
    budget=Budget(work_limit=2_000_000, seconds=30, cpu_seconds=30),
)
assert run.result.reconstruct() == n
print(run.reason, run.result.complete)
for factor in run.result.factors:
    print(factor.value, factor.exponent, factor.certainty.value)
print("Unresolved:", run.result.remaining)
```

`run.result.complete` means no unresolved cofactors remain; `proven` additionally
requires proven certainty for every factor. The compatibility API
`v2.factor.factorize()` remains available. Individual rho, p−1, ECM and QS jobs
return splits; recursive classification belongs to the portfolio.

The [API guide](docs/API.md) covers result types and low-level imports.
[Checkpoint guidance](docs/CHECKPOINTS.md) explains portable resume, total
allowances and charged validation. Saved jobs retain their arithmetic backend,
configuration and consumed work.

## Automatic ECM workers

The bounded auto CLI and `factorize_bounded()` use a finite, cost-gated
process executor for independent ECM curves on supported macOS/arm64 PyPy.
The normal portfolio retains preprocessing, recursive completion and its
existing SIQS/SSS fallback. The ECM-only benchmark runner is separate.

Use `--ecm-workers serial` or `ecm_executor=None` to keep incumbent serial
execution. `--ecm-memory-mib` / `ecm_memory_bytes` set the aggregate owned ECM
allowance, capped at 8 GiB; caller work, CPU, wall and parent-memory limits
still apply. Unsupported, unmeasured and unfunded calls stay serial.

For repeated library calls, a caller-owned context keeps its pool and immutable
programs warm:

```python
from v2.ecm.auto import selected_context

with selected_context(seed=7, prefer_warm=True) as context:
    run = factorize_bounded(n, seed=7, ecm_executor=context)
```

The context's seed, resource grants and selected worker configuration remain
bound for its lifetime. Do not share one active context across concurrent
calls. Use the same declared configuration and adequate finite budgets.
Cold CLI startup and warm retained-pool performance are separate endpoints.
The compact64 default at 30 digits and fixed32 defaults at 50/70 digits stay
serial cold before loading the cost catalogue. Twenty-digit automatic calls
also remain serial. Larger explicit numerical schedules retain selection. See the
[selected regime table](STATUS.md#selected-ecm-regimes) for cold/warm worker
counts and retained performance limitations.

Worker arithmetic remains exact. Stable curve identities, checked divisor
publication, bounded parent-owned leases, explicit spawn, cancellation/drain
and quiescent checkpoints preserve finite execution and checked resume.
Automatic checkpoint schema 13 binds policy, semantic implementation version,
seed and grants;
older serial checkpoints remain serial. Worker-count changes on resume are
rejected. The worker parent ceiling is 4 GiB including selector bookkeeping;
portfolio admission failures retain serial execution, while direct selected
component calls can raise `MemoryError`.

## Optional GMP backend

Install into the selected PyPy environment:

```sh
pypy3 -m venv v2/.venv
v2/.venv/bin/python -m pip install -r v2/requirements-gmp.txt
v2/.venv/bin/python -m v2.factor 626100403 --bounded --backend gmpy2-mpz
```

The default remains `python-int`; explicit `gmpy2-mpz` selection requires the
optional dependency. See [arithmetic backends](common/README.md#arithmetic-backends)
for compatible configurations and checkpoints.

## Methods and packages

| Path | Contents |
| --- | --- |
| [factor.py](factor.py), [portfolio.py](portfolio.py) | Public API/CLI, recursive bounded factoring and stage/fallback coordination. |
| [common/](common/README.md) | Exact arithmetic, preprocessing, prime generation and primality certainty. |
| [rho/](rho/README.md) | Brent rho splitting. |
| [pm1/](pm1/README.md) | Pollard p−1 and bounded campaigns. |
| [ecm/](ecm/README.md) | Montgomery ECM, verified chains and continuation programs. |
| [qs/](qs/README.md) | QS/MPQS/SIQS, SSS/SSSf, relation stores, filtering and extraction. |
| [execution/](execution/README.md) | Shared budgets, stage jobs, schedules and resource allocation. |
| [docs/](docs/) | CLI, API and checkpoint guides. |
| [tests/](tests/README.md) | Arithmetic, ownership, continuation and integration checks. |
| [benchmarks/](benchmarks/README.md) | Reusable runners, versioned inputs and measurement summaries. |

Collector configuration and storage contracts live in
[QS configuration](qs/CONFIGURATION.md) and [collector policies](qs/COLLECTORS.md).
[STATUS](STATUS.md), [ROADMAP](ROADMAP.md) and [POLICY](POLICY.md) track current
behavior, future work and experiment rules. [Results](benchmarks/RESULTS.md)
and the [changelog](../CHANGELOG.md) retain scoped decisions. Generated evidence,
research and scratch remain local in ignored results/audit directories.

## Development

```sh
pypy3 -m venv v2/.venv
v2/.venv/bin/python -m pip install ruff==0.14.14 pycodestyle==2.14.0
make -C v2 test
make -C v2 lint
```

Use the shared performance lock for full suites and measurements. See
[tests](tests/README.md) and [contributor guidance](../AGENTS.md) for the workflow.

## License

Active v2 code is licensed under the [MIT License](LICENSE). It permits reuse,
modification and redistribution, including commercial use, with the copyright
and permission notices retained. It is supplied without warranty.

This grant covers v2 project code and documentation. External dependencies
and linked papers retain their own terms. The preserved v1 baseline is outside
this new grant; its original source notices and provenance remain unchanged.
