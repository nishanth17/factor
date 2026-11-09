# Factor v2

Factor is an integer-factorization library and command-line tool for **PyPy
implementing Python 3.11**. It combines small-factor searches with relation
engines, validates every split and recursively factors the remaining pieces.
Results preserve multiplicities, signs and unresolved cofactors, and distinguish
probable primes from proven primes.

The default Python-integer backend uses only the standard library. The active
implementation lives here; [v1](https://github.com/nishanth17/factor/blob/23458a52e409add326c9fd63cf1f322c6e605569/v2/history/V1_BASELINE.md) is the preserved
Python 2 baseline.

An experimental [coarse scalar C/GMP ECM packet](ecm/native/README.md) is
available through `PortfolioConfig(ecm_kernel="c-gmp")` or
`--ecm-kernel c-gmp`. It requires an explicit local build and retains the Python
fallback. A separate experimental OpenMP owner supplies static/dynamic/guided
continuation scheduling. Both builds and executable QA are pending; no performance
or default adoption is claimed.

A separate optional [ECM NEON feasibility prototype](ecm/native/neon/README.md)
provides two-curve products/squares at 200/266/333 bits. It is not integrated into
the executor; only a separately admitted Stage A screen may measure it.

## Setup and first run

With PyPy implementing Python 3.11 installed, run from the repository root:

```sh
make -C v2 runtime
pypy3 -m v2.factor 626100403 --bounded --seed 7 \
  --seconds 30 --cpu-seconds 30
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
root resource owner:

```python
from v2.execution.budget import Budget
from v2.portfolio import factorize_bounded

n = 626100403
run = factorize_bounded(
    n,
    seed=7,
    budget=Budget(seconds=30, cpu_seconds=30),
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

Use `--ecm-workers serial` or `ecm_executor=None` for serial execution.
Worker counts, numerical schedules, packet widths, caches and cost/size guards
retain their selected values. Worker admission no longer requires prepaid work,
CPU or estimated storage grants. Their dependent cost tables need revalidation;
H3 makes no performance claim.

For repeated library calls, a caller-owned context keeps its pool and immutable
programs warm:

```python
from v2.ecm.auto import selected_context

with selected_context(seed=7, prefer_warm=True) as context:
    run = factorize_bounded(n, seed=7, ecm_executor=context)
```

The context's seed and selected worker configuration remain
bound for its lifetime. Do not share one active context across concurrent
calls. Use the same declared numerical configuration and shared root limits.
Cold CLI startup and warm retained-pool performance are separate endpoints.
The compact64 default at 30 digits and fixed32 defaults at 50/70 digits stay
serial cold before loading the cost catalogue. Twenty-digit automatic calls
also remain serial. Larger explicit numerical schedules retain selection. See the
[selected regime table](STATUS.md#selected-ecm-regimes) for cold/warm worker
counts and retained performance limitations.

Worker arithmetic remains exact. Tokens identify assignments; they carry no
resource grants. The root observes parent, workers, owned helpers and descendants.
Live samples, process-lifetime acknowledgments and final libproc self CPU
before wait4 reaping reconcile
by process birth identity, including crashes and shutdown. Retained context
startup, each call and close own active intervals; idle time is excluded.

## Resource API

`Budget(seconds=None, cpu_seconds=None)` leaves wall and aggregate CPU unlimited.
The default RSS guard is `floor(4 * physical_memory / 5)`. Set `rss_bytes` to an
integer byte cap or explicitly pass `rss_bytes=None` to disable it. CLI equivalents
are `--seconds`, `--cpu-seconds`, `--rss-bytes` and `--no-rss-limit`. Library calls
remain quiet. `used` counts algorithmic work for diagnostics; it cannot refuse work.

macOS uses libproc through ctypes to sum current RSS over the owned process tree.
The parent is included; shared pages can be counted multiple times. Samples run
at entry and publication boundaries. Ordinary cooperative polls throttle
sampling to once per 50 ms; a long atomic operation can delay the next sample. Allocations, native operations and active worker stages can
overshoot the guard. Provider failures are explicit errors.

Only root wall, aggregate CPU, observed RSS and an optional checkpoint save cap
refuse resources. Internal stores, graphs, matrices and scratch grow as needed.
Partial FIFO capacities, caches, packet widths, numerical search effort and exact
representation/provenance checks retain their policies. Stage allocations order
work: exhausting a stage share hands remaining root resources to the next stage;
an underfunded fallback reserve no longer vetoes that handoff.

Removed arguments raise `RemovedResourceLimitError`, naming the argument and its
replacement. This includes `work_limit`, internal `memory_bytes`, `max_atoms`,
`max_relations`, child `checkpoint_bytes`, worker grants and the former CLI flags
`--work-limit`, `--memory-mib`, `--ecm-memory-mib`. Cache byte capacities remain
supported where they control retention. `default_siqs_config(recipe="calibrated")`
retains the omitted-grant recipe; `recipe="reference"` selects the former explicit
non-calibrated recipe without making a resource cap select mathematics.

Resource stops expose the limit name, observed/requested value, limit value and
scope. Cancellation, finite numerical exhaustion, `allocation_failure`, I/O
failure and unsupported encoding retain distinct identities. Interrupted
publication poisons its owner; it cannot be reused or checkpointed.

Checkpoints use resource version `observed-resources/saved-history-v3`. Older
resource schemas raise `UnsupportedResourceVersionError`; there is no migration.
Resumes inherit recorded cumulative charges plus the current invocation. Work
and resource consumption after the last valid save are **not recovered** after
a crash or skipped save. A caller retaining the live `Budget` keeps those unsaved
charges across retries. No durable ledger is built; benchmark/study ledgers remain
the authority for experiment charges. Paused intervals are excluded.

Saving is optional. `atomic_save(path, snapshot, checkpoint_bytes=None)` has no
size cap by default; CLI `--checkpoint-bytes` supplies one. Oversize records
`checkpoint_size_skip` with requested/limit bytes, preserves the old file and
continues factoring. Saves use only published state, a temporary sibling,
flush/fsync, atomic `os.replace` and directory fsync. Before replacement, an interrupted write leaves the previous checkpoint valid.
Atomic replacement exposes a complete old or new file; a later directory-fsync
failure can leave the complete new file installed. Allocation or I/O failure retains its own cause.

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

On the B4a task branch, optional native ECM packets and teams follow the H3
root resource API. Their former `memory_bytes` arguments are rejected; scratch
estimates are diagnostics, not admission caps. Native team checkpoints use a
new H3 identity and retain saved wall/CPU history. Historical B4a recorder
work frontiers remain benchmark controls and do not restore production limits.
Native builds and performance admission remain separate pending checks.

B4a's owned-group CPU observer converts Mach absolute counters through the
host timebase. Process totals include native-thread CPU; thread totals are not
added again. Native build manifests retain GMP header/runtime identities and
actual linked libraries. Missing explicit native backends report named errors.
