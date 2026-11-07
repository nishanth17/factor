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

## Opt-in current QS workers

`v2.qs.current_parallel.CurrentParallelSIQSJob` accepts a concrete current
`SIQSConfig` or `CountedSIQSConfig`, an explicit parent memory allowance,
`Budget`, and `QSWorkerConfig`. Growing DLP and the bounded C9a cost-report
configuration are supported. This experimental interface preserves current
width, multiplier, parity/sign, root, queue and exact recovery policies; it
uses one parent-owned relation/partial/cycle store.

`QSWorkerConfig(workers=1)` through `QSWorkerConfig(workers=10)` select spawned
PyPy3.11 processes, subject to the explicit aggregate ownership allowance;
`QSWorkerConfig(mode="serial")` executes the packet stream locally. Pass an
explicit `CollectionPool` to `run(pool=...)` to retain processes across jobs.
`factorize_bounded(..., qs_workers=settings, qs_pool=pool)` applies the opt-in
path only at SIQS handoff, including recursive cofactors. With no pool argument
it owns one lazy QS pool for the full call. Retained ECM pool ownership is
charged alongside QS workers using its finite memory cap or executor ownership
bound; a borrowed ECM context with unknown ownership is refused. An unused
automatically owned ECM context with an empty curve schedule needs no pool
reserve. Warm B12 ECM processes cannot yet dispatch QS packets.

Configuration and worker count/mode remain fixed across resume. Extend the
original budget or restore `CurrentParallelSIQSJob.from_checkpoint(...)` with
cumulative allowances. Quiescent v2 checkpoints retain paid prefixes, pending
admission, original provenance and split-call leases. Derived solver/root/proof
state rebuilds with additional charges. Lost reports retain full work and finite
split grants; unknown uncapped split counts stay explicit. Positive splitter
leases stop before additional recovery and resume the same unpublished cursor.
Quota scarcity can change scheduling/yield and must be separated from pure
concurrency. Nested recovery, spill, width growth and native-thread execution
are refused. Ownership accounting is prospective; it is not an OS RSS sandbox.

Defaults remain serial. The experimental `lease_budget=True` option uses plain
packet-local work
counters and parent aggregate CPU guarding. It retains immediate work refusal,
private CPU/wall limits, forced publication checks and the 64-action stop
cadence. Lost reports retain full grants. Stop latency is cooperative, including
indivisible integer operations and the funded publication tail. Enabled lease
policy is part of execution identity and cannot change across resume. Disabled
options are omitted from legacy settings/checkpoint identities. No default or
performance adoption follows from this option.

The opt-in `packed_queue` marking component stores each sparse hit in one
64-bit word instead of three 32-bit words. For example, select
`marking="presieve+unroll4+parity_two+packed_queue"` with a funded sparse
queue. Offset/index/weight widths are checked before construction; unsupported
shapes use the exact direct path. Hit order, work charges, polling and the
older conservative memory reservation remain unchanged. Derived queues are
rebuilt after interruption/restore; the marking selection pins configuration
identity. This changes no automatic default or admission policy.

Explicit QS worker aggregate grants may reach 22 GiB; the default remains
12 GiB. Prospective ownership still includes parent state, every worker,
decoded/in-flight publication, coexistence and checkpoint reserves. This is
an ownership bound, separate from measured RSS and physical-memory pressure.

The opt-in `wide_hits` marking component permits up to 128 MiB of lazy
native-word hit scratch per collector instead of the usual 32 MiB ceiling.
The actual reservation and construction work must fit the existing collector
grant, leaving polynomial-switch coexistence funded; otherwise the ordinary
exact mask path remains available. This can keep larger bases off repeated
large-integer mask operations. It follows `packed_queue` in component order.

`QSWorkerConfig(coefficient_cache=True)` retains one bounded parent CRT
template within each packet merge. It uses the factor base's checked column
map, selected roots and cached modular terms while preserving every original
per-coordinate work reservation and poll. The 64 KiB reserve is included in
aggregate ownership. The cache is discarded on merge exit and rebuilt after
resume; the enabled setting pins execution identity. The default is false.

The opt-in `wide_inverses` marking component permits up to 32 MiB of shared
Hensel-derivative inverse scratch instead of 8 MiB. A fixed target has at
most two opposite nonsingular derivatives per odd prime, so the existing
checked cache reuses inverses across polynomials and A values. Its complete
reservation and the next polynomial switch must fit the collector grant;
scarce memory retains independent exact inversion. No default changes.

`private_stop=True` requires process packet leases and gives each worker its
own synchronized stop flag. The parent sets or clears every flag, while each
worker retains the same 64-action polling cadence and exact budget checks.
This avoids sharing one flag lock among all worker readers. Broadcast is
cooperative and sequential; it does not promise a wall-time cancellation bound.
The option pins execution/checkpoint identity and must be selected before a
pool spawns any processes. Reused pools cannot change the policy after spawn.
Aggregate ownership reserves an additional 64 KiB per worker for these handles.
The disabled option preserves the legacy identity and remains the default.

With private process flags, `stop_backend="semaphore"` selects one nonblocking
semaphore per worker instead of an Event condition. Each flag has exactly one
worker reader. The parent sets idempotently during a run and clears only after
all packets and snapshot barriers drain. Setup checks also select that worker's
own flag; parent observation reads parent-owned state. These ownership rules
permit removing the Event condition lock while retaining the polling cadence.
The backend is pinned before spawn and across checkpoints. Default `"event"`
preserves existing private-flag identities. No raw shared-memory atomics or
change to arithmetic, work reservations or private deadline checks is involved.

`cpu_snapshot=True` also requires process packet leases. Parent verification
and solving read the worker CPU array under one lock against a run-local
baseline, instead of locking each scalar read. Aggregate charging, immediate
budget refusal and the parent polling cadence remain unchanged. The baseline
is rebuilt after each quiescent pool `begin`; the enabled option pins execution
identity and defaults off.

`reuse_state=True` installs an identity-bound immutable context on every worker
before sending packet deltas. `state_cache_bytes` bounds that context; aggregate
admission also reserves setup-message coexistence. Private family roots and
collector buffers can persist, with current-lease budget rebinding. Proven
prime/unit caches remain reusable; polynomial queues and packet-dependent
recovery hints are invalidated. Checkpoints store paid progress and context
identity, rebuild disposable caches, and charge new installation. A missing
installation acknowledgement retains prepaid setup work and conservative CPU
charges, and retires only the owned pool.
A retained-state pool pins its worker settings for its lifetime; changing or
disabling them requires a fresh pool. Independent jobs with the same settings
may install a different target/configuration context after checked retirement.

`packed_transport=True` publishes bounded varint relation rows. The legacy JSON
admission limit and work charge remain in force, so encoding cannot change a
packet's accepted prefix. Decoding checks dimensions and expansion against the
norm/support bound before exact relation verification. Mutable transport buffers
are private and reusable during an active run; quiescent pool finish releases
them. Published bytes do not alias those buffers. Pending packed checkpoint
payloads use strict base64 with an explicit encoding marker; disabled transport
retains the legacy JSON format.

`verification="worker"` fully checks atoms before worker publication and uses
the parent's checked ordered insertion without a redundant standalone check.
The parent still independently rejects inconsistent arithmetic; a worker claim
is not a proof. `verification="batch"` prevalidates bounded groups before
ordered insertion, reusing only object-bound proofs in the parent. Its
`verification_batch_atoms` allowance is included in aggregate memory admission.
Unsupported unchecked collector paths retain full parent verification. Default
`verification="parent"` preserves the original checks and legacy identity.

`lookahead=N` adds at most N queued/pending packet credits beyond the worker
count (maximum four extra credits per worker), with prospective response-memory
and work/CPU reservations. Complete portfolio calls then keep the executor
active instead of draining each worker-sized wave. Explicit fixed-work quotas
still bound total issuance; admission remains ordinal. `whole_family=True`
groups the entire configured family prefix without adding polynomials. Recovery
hint boundaries retain the original chunk scope. Flagged checkpoints validate
contiguous paid receipts for fully admitted partial packets as well as pending
work; legacy off identities and formats remain unchanged.

`arrival_order=True` is a separate, trajectory-changing experiment. Complete
packets can enter the graph in observed arrival order while paid witnesses
retain the bounded ordinal frontier. Partial packets retain their admission
cursor; checkpoints validate each arrival against its acknowledgement. This
option rejects fixed-work execution and is excluded from the combined
fixed-order comparison. It defaults off.

The numerical runner is
`pypy3 -B -m v2.benchmarks.qs.a7a.runner`; its `plan`/`freeze` commands are
metadata-only. `qa` and `run` require `--execute-phase3`; full runtime checks and
exploratory performance screens have run. Population confirmation, broader
size coverage and automatic worker selection remain pending. Detailed local
contracts,
research, budgets and handoff are under `v2/audit/a7a/`.
