# Montgomery ECM

This package implements Suyama Montgomery-curve setup, scalar arithmetic,
two-stage ECM, verified stage-one plans and optional paired continuations.
Use the [bounded portfolio](../docs/API.md#bounded-recursive-factoring) when
curves must share a deadline, work allowance and checkpoints.

## Standalone split

```python
from v2.ecm import factorize_ecm

n = 4001 * 5003
divisor = factorize_ecm(n, seed=7, b1=2000, b2=50000, max_curves=8)
if divisor is not None:
    assert 1 < divisor < n and n % divisor == 0
```

`factorize_ecm` returns a proper divisor or `None`. Bounds are inclusive,
an unusable curve counts as an attempt, and the final permitted curve can
still return a factor. Singular setup and saturated batches produce proper
factors, checked recovery or a finite retry. The standalone call has a finite
curve count but no shared portfolio deadline or portable checkpoint.

## Bounded stage-one plans

`PortfolioConfig.ecm_chain_mode` accepts `auto`, `off` and `reuse`.
`off` keeps the B4 ladder. Supported native automatic jobs use reduced
PRAC/batch16; explicit Lucas/CF alternatives use `ecm_chain_family`.
Automatic GMP jobs keep the ladder, while explicit GMP reuse has its own
validated route. The [status table](../STATUS.md#current-automatic-defaults)
records current routing and the adoption boundary.

A plan requires compatible bounds, chunks, curve count and owned storage;
execution is restricted to supported 40–80-digit cofactors. Unsupported jobs
keep the ladder. See [eligibility](chain_options.py) and [plan bounds](chains.py)
for exact admission tests. The CLI exposes `--ecm-chain auto|off|prac|lucas|cf`.

Each invocation owns a finite plan LRU shared across its curves and children.
Plans contain immutable verified records, never shared curve points. Explicit
reuse needs program and chain caps that fit with the rest of the job. A miss
charges reading, verification and construction; a hit still incurs its charge.
Fast chunks reserve execution and bounded strict recovery before arithmetic.
These work units differ from ladder units.

Checkpoints pin the catalog, executor, backend/build, chunk and recovery
identities. Restore checks the committed prefix, rebuilds disposable plans
and charges the work. Legacy snapshots retain their original executor.

## Arithmetic helpers

`setup_curve`, `point_add`, `point_double`, `scalar_multiply`, `multiply_prac`,
`stage_one_scalar` and `stage_two` are exported from `v2.ecm`. Coordinates
represent X:Z on a nonsingular Montgomery curve over an odd modulus, with
`a24=(A+2)/4` as returned by setup. Projectively equivalent outputs need not
have identical coordinate pairs.

`prac.get_chain`, `verify_chain` and `multiply` separate bounded integer-chain
compilation, independent verification and point execution. The compiler caps
scalars, candidates, instructions and its record cache; larger unsupported
scalars use the checked ladder. Caller-supplied records are reverified.

`multiply_prac` can raise `prac.NonunitPointError`. A proper `factor` must be
returned to factoring; `factor=None` requests a curve retry. Invalid `(0,0)`
states are never accepted. The low-level helper does not supply a portfolio
budget, cancellation policy or serialized resume.

## Programs, pairing and campaigns

`ecm_program_bytes` reserves immutable packed prime/power blocks; points,
products and recovery stay private. Programs regenerate when their bounded
retention fills and rebuild with charges after resume. Packed endpoints
require B2 below `2**64`.

`ecm_pair_distance` enables the experimental bounded ± executor; `None`
keeps unpaired execution and `0` is a direct-scalar control.
`ecm_pair_wheel` selects the separate aligned-wheel alternative. Both need a
program cap, exact prime-coverage certificates and private table storage.
Saturated products replay terms and primes under the same recovery contract.
The two pairing layouts remain opt-in after their scoped losses; see
[historical decisions](../benchmarks/RESULTS.md).

For a deeper campaign, declare finite `(B1, B2, curves)` tiers through
`PortfolioConfig`. Budget extensions preserve those tiers and cumulative
consumption. Adding curves to an exhausted saved schedule or increasing B1
on the same completed curve needs a separate migration contract. A new tier
starts independent curves. [Execution](../execution/README.md) describes
protected handoff to relation engines.

## Source map

| Module | Responsibility |
| --- | --- |
| [core.py](core.py) | Setup, point formulas, ladder, standalone stages/splitter |
| [prac.py](prac.py) | Verified differential-chain generation and execution |
| [chains.py](chains.py), [chain_options.py](chain_options.py), [chain_records.py](chain_records.py) | Bounded production plans, records, selection and recovery |
| [programs.py](programs.py) | Reusable prime/power blocks and coverage certificates |
| [paired.py](paired.py), [wheel.py](wheel.py) | Optional paired stage-two executors |

Research-only catalog/executor APIs live under
[ECM benchmarks](../benchmarks/ecm/README.md). Current production behavior is
separate from historical kernel and reused-program measurements.

## Execution and bundled data

ECM-specific execution lives in this flat package: `auto.py`, `selection.py`,
`workers.py`, `worker.py`, `worker_state.py`, `stage_two.py`, `defaults.py` and
`presets.py`. `execution/` retains shared budgets, stage coordination and
contracts. Old `v2.execution.ecm_*` imports are compatibility aliases to these
modules, including for existing callers and test patch points.

`chain_data.py` and `cf_data.py` contain immutable integer/bytes tuples. Plans
load them lazily and still verify their mathematical actions and recovery
frontiers. The retained JSON controls are independent test oracles, not runtime
inputs. Preparation keeps its existing versioned work charges.

`worker_policy.py` holds immutable measured settings and returns owned copies;
its numerical rows and worker choices are unchanged. Runtime selection checks
explicit compatibility identifiers, not source files or their paths. Optional
external profiles remain bounded JSON, using selector format v3. Automatic
metadata uses `b12-f4-auto-v2`; old automatic v1 checkpoints and selector v2
profiles are rejected rather than silently reinterpreted. Existing serial
checkpoint handling and worker arithmetic/search contracts are unchanged.
