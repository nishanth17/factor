# Python API and results

Import through the `v2` package from the repository root on PyPy 3.11.
The [CLI](CLI.md) and [package index](../README.md#methods-and-packages)
provide other entry points.

## Factorization results

```python
from v2.factor import factorize, print_factorization

result = factorize(626100403, seed=7)
assert result.reconstruct() == result.original
print(print_factorization(result.original, result))

for factor in result.factors:
    print(factor.value, factor.exponent, factor.certainty.value)
```

`FactorizationResult` exposes `original`, `sign`, `factors` and `remaining`.
Each terminal factor has a value, exponent and certainty label. Every result
reconstructs its original input, including unresolved cofactors. `complete`
means no unresolved composite remains; `proven` also requires proven certainty
for every terminal factor. A probable prime stays labelled probable.

`factorize()` is the compatibility API. Its per-method controls do not create
a single shared deadline or resumable portfolio. Use the bounded API below
for those contracts. Zero is invalid; one and negative inputs are defined.

## Bounded recursive factoring

```python
from v2.execution.budget import Budget
from v2.portfolio import PortfolioConfig, factorize_bounded

n = 626100403
run = factorize_bounded(
    n,
    seed=7,
    config=PortfolioConfig(),
    budget=Budget(work_limit=2_000_000, seconds=30, cpu_seconds=30),
)
assert run.result.reconstruct() == n
print(run.reason, run.work_used)
```

The run carries `result`, `reason`, cumulative work/time, events and checkpoint
state. One allowance covers preprocessing, retries and recursive children.
`Budget` reservations occur before atomic mutation. Work units are algorithmic
charges and are not comparable time units across different methods.

Omitting `config` permits automatic routing. `PortfolioConfig().for_input(n)`
resolves the root's concrete plan. Explicit configurations, `siqs=None` and
saved numerical settings preserve their intent. `asdict`/`dataclasses.replace`
copies freeze resolved fields; they do not request a new automatic selection.
Current size/resource choices and evidence live in [STATUS](../STATUS.md).

`PortfolioConfig(siqs=SIQSConfig(...))` supplies an explicit SIQS fallback;
`PortfolioConfig(sss=SSSConfig(...))` supplies experimental SSS. Import those
classes from `v2.qs` and `v2.qs.sss`. A fallback's arithmetic backend must
match its portfolio, and simultaneous storage must fit the outer allowance.

## Low-level components

| Import | Contract / guide |
| --- | --- |
| `v2.rho.brent.factorize_rho` | Proper divisor or `None`; [rho](../rho/README.md) |
| `v2.pm1.core.factorize_pm1` | Two-stage split; [p−1](../pm1/README.md) |
| `v2.pm1.bounded.factorize_pm1_bounded` | Finite same-base campaign and checked resume |
| `v2.ecm.factorize_ecm` | Finite curves; [ECM](../ecm/README.md) |
| `v2.qs.SIQSJob`, `QSJob` | Split/cofactor job results; [QS](../qs/README.md) |
| `v2.common.arithmetic`, `utils` | Exact arithmetic and [certainty](../common/README.md) |

Standalone splitters return a proper divisor or `None`; callers retain the
cofactor and perform recursive classification themselves. They do not inherit
the portfolio's shared deadline. Direct QS/SSS jobs likewise split a modulus;
full recursive terminal classification belongs to `factorize_bounded`.

## Backends and resume

Use `PortfolioConfig(backend="gmpy2-mpz")` or the compatibility API's
`backend=` argument for explicit GMP arithmetic. High-level results and saved
checkpoints contain canonical Python integers. Low-level objects may retain
`mpz`; see [arithmetic backends](../common/README.md#arithmetic-backends).

For serialized resume, pass `run.checkpoint` and a fresh `Budget` containing
total allowances, including previous consumption. `config=None` restores the
saved portfolio configuration; an explicit conflicting configuration is
rejected. [Checkpoint guidance](CHECKPOINTS.md) describes replay costs and
algorithm-specific continuation limits.
