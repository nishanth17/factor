# Brent rho

Seeded Pollard rho with Brent-style walks, batched GCDs and finite saturation
recovery. The implementation is [brent.py](brent.py), targeting PyPy 3.11.

## Standalone splitter

```python
from v2.rho.brent import factorize_rho

n = 4001 * 5003
divisor = factorize_rho(n, seed=7, max_attempts=4, max_evaluations=5000)
if divisor is not None:
    assert 1 < divisor < n and n % divisor == 0
```

The return value is a proper divisor or `None`. Failed attempts continue to
the next seeded walk until the explicit attempt/evaluation limits end.
Saturated products receive bounded recovery; a recovery cap cannot hang or
turn n itself into a successful split.

## Controls and integration

`seed` or `rng` determines assignments. `max_attempts`, `max_evaluations`,
`batch_size` and `recovery_limit` bound the work. An optional `RhoStats`
instance records evaluations, GCDs, attempts and recovery. `backend=` selects
exact arithmetic; omission follows the exact input representation.

The standalone splitter does not provide recursive terminal classification,
a shared deadline or portable checkpoints. Use
[`factorize_bounded`](../docs/API.md#bounded-recursive-factoring) for those
contracts. Portfolio rho shares allowances with preprocessing, p−1, ECM and
relation fallbacks; tuning its batches/restarts is separate from this API.

See [current status](../STATUS.md) for portfolio choices and
[benchmark guidance](../benchmarks/README.md) for matched measurements.
