# Pollard p−1

Two-stage Pollard p−1 and bounded same-base continuation for PyPy 3.11.
Use the [portfolio API](../docs/API.md) for recursive factoring under one shared
allowance. [Current status](../STATUS.md) records selected defaults.

## Standalone splitter

```python
from v2.pm1.core import factorize_pm1

n = 607 * 1019
divisor = factorize_pm1(n, b1=10, b2=200, base=2)
if divisor is not None:
    assert 1 < divisor < n and n % divisor == 0
```

B1 and B2 are inclusive. Stage two targets a prime divisor whose p−1 has an
eligible remaining large prime. The result is a proper divisor or `None`;
recursive classification belongs to the caller. Saturated stage-one chunks
and stage-two batches use finite recovery or stop that base.

## Finite bound campaigns

```python
from v2.execution.budget import Budget
from v2.pm1.bounded import PM1Config, factorize_pm1_bounded

config = PM1Config(bounds=((15, 15), (16, 16)))
run = factorize_pm1_bounded(
    17 * 1019, base=3, config=config,
    budget=Budget(work_limit=100000),
)
assert run.divisor == 17
assert run.result.reconstruct() == 17 * 1019
```

This Python-integer API predeclares monotone inclusive bound rungs for one
base. A B1 increase applies the exact ratio of least-common-multiple
exponents, including higher powers of old primes, and restarts stage two.
A B2-only increase appends coverage. Saturation ends that base; larger bounds
cannot undo an identity residue. Repeated p−1 bases are not independent
ECM-like smooth-order trials.

`max_actions` pauses at committed boundaries. Restore with the same modulus,
base and configuration, `checkpoint=run.checkpoint`, and an unused `Budget`
containing total cumulative allowances. Saved arithmetic is reconstructed
and verified with charges before reuse. A small grant may be spent entirely
on replay. New rungs or bases require a new declared campaign.

`PM1Run` carries the divisor, reason, reconstructible result, cumulative
work/time, verification work and checkpoint. Split pieces remain unresolved;
this low-level API makes no primality assertion.

## Execution choices

An omitted campaign config selects the recurrence/chunk-64 executor for fresh
calls. Explicit `PM1Config` preserves the legacy executor; `PM1TuningConfig`
from `tuning.py` exposes versioned alternatives. Saved jobs retain their
executor and accounting identity. Portfolio p−1 uses its own shared stage-job
integration and the selected fields in `PortfolioConfig`.

[core.py](core.py) contains the standalone splitter, [bounded.py](bounded.py)
the campaign state machine, [gaps.py](gaps.py) gap tables and exact extension
ratios, and [tuning.py](tuning.py) explicit execution alternatives.
Measurements and the adoption limits are in
[benchmark results](../benchmarks/RESULTS.md) and the
[historical A6 documentation](https://github.com/nishanth17/factor/blob/24ff153d26ada2ec4beeead07000d183e94fd355/v2/README.md#a6-production-default-promotion--9-october-2026).
