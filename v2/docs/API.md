# Python API and results

H3 replaces internal refusal quotas with observed root limits. See the
[resource API](../README.md#resource-api) for defaults, removed arguments,
sampled RSS, stage handoff, saved-history accounting and optional atomic saves.


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

`factorize()` accepts a shared `budget` for its recursive calls. Use the bounded
API below for resumable portfolio state. Zero is invalid; one and negative inputs are defined.

## Bounded recursive factoring

```python
from v2.execution.budget import Budget
from v2.portfolio import PortfolioConfig, factorize_bounded

n = 626100403
run = factorize_bounded(
    n,
    seed=7,
    config=PortfolioConfig(),
    budget=Budget(seconds=30, cpu_seconds=30),
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
match its portfolio. Storage grows under the root observed RSS guard.

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
a new independent owner when called from a portfolio; they borrow the active
root. Direct QS/SSS jobs likewise split a modulus;
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

## Automatic and retained ECM contexts

`factorize_bounded(n)` defaults to `ecm_executor="auto"`. The call lazily selects
a funded process pool and owns its cleanup, including exception exits. Explicit
`ecm_executor=None` selects serial execution. Other context objects remain
borrowed; the caller closes them. A serial override retains the declared current
numerical schedule, not a historical version of all defaults.

Automatic ECM borrows the same root wall, aggregate CPU and observed RSS guard.
There are no prepaid startup/packet grants or internal ownership caps. Cost and
numerical eligibility still select the existing worker recipes. Removed
`ecm_memory_bytes` raises `RemovedResourceLimitError`.

```python
from v2.ecm.auto import selected_context

with selected_context(seed=7, prefer_warm=True) as context:
    for n in inputs:
        run = factorize_bounded(
            n, seed=7, ecm_executor=context,
            budget=Budget(seconds=30, cpu_seconds=4096),
        )
        assert run.result.reconstruct() == n
```

Warm intent uses the recorded horizon of at least 64 compatible ECM-reaching
calls; the first call still pays cold startup. Use cold intent for short or
uncertain reuse. The context fixes its first engaged worker configuration and
needs a parent-memory reserve matching any explicit portfolio workspace.
Automatic schema-13 snapshots require the same packaged policy/implementation,
seed on resume. Pre-H3 resource schemas are refused without migration.

Concurrent independent calls own separate contexts and resource owners. Service callers
coordinate machine-wide concurrency; nested process calls default to serial.
A shared retained context admits one active call at a time.
