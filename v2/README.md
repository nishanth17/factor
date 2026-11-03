# Factor v2

The active Factor implementation runs on **PyPy implementing Python 3.11**.
It uses exact integer arithmetic, validates every split and retains unfinished
cofactors. The default Python-integer backend needs no third-party packages.

## Setup and first run

Run commands from the repository root:

```sh
make -C v2 runtime
pypy3 -m v2.factor 626100403 --bounded --seed 7 \
  --work-limit 2000000 --seconds 30 --cpu-seconds 30
pypy3 -m v2.factor --help
```

Use `PYTHON=/path/to/pypy3` with `make` to select another PyPy 3.11 installation.
An omitted input opens an interactive prompt. Calls can return a complete
factorization or explicit unresolved cofactors when an allowance ends.

The [CLI guide](docs/CLI.md) covers method selection and resource flags.
The [API guide](docs/API.md) covers library results and bounded calls;
[checkpoint guidance](docs/CHECKPOINTS.md) covers saving and resuming progress.

## Methods and packages

| Package | Responsibility |
| --- | --- |
| [common/](common/README.md) | Exact arithmetic, preprocessing, prime sieves and certainty |
| [rho/](rho/README.md) | Seeded Brent rho with bounded recovery |
| [pm1/](pm1/README.md) | Pollard p−1 and finite bound campaigns |
| [ecm/](ecm/README.md) | Montgomery ECM, verified chains and continuation programs |
| [qs/](qs/README.md) | QS, MPQS, SIQS, SSS/SSSf and verified relation processing |
| [execution/](execution/README.md) | Shared budgets, schedules, stage jobs and allocation |

`factor.py` is the public API/CLI; `portfolio.py` coordinates bounded recursive
factoring. The unbounded compatibility API remains available. Use bounded
execution when a call needs shared work/time limits or resumable progress.

See [STATUS.md](STATUS.md) for automatic size/resource routing and its adoption
limits. Fresh automatic SIQS uses the C10a-01 coupled polynomial/multiplier/A/root
bundle at every size; see the [configuration semantics](qs/CONFIGURATION.md#coupled-polynomial-policies).
Existing size parameters, finite grants and portfolio routing remain scoped.
Explicit settings and saved configurations take precedence. Experimental
methods and adopted defaults retain their recorded evidence limits;
there is no general digit-based completion guarantee.

## Optional GMP backend

Install GMP support into a PyPy 3.11 environment with
`pip install -r v2/requirements-gmp.txt`, then select `--backend gmpy2-mpz` or
`PortfolioConfig(backend="gmpy2-mpz")`. A missing dependency is an explicit error.
The [shared arithmetic guide](common/README.md#arithmetic-backends) explains
backend and checkpoint compatibility.

## Development

```sh
pypy3 -m venv v2/.venv
v2/.venv/bin/python -m pip install ruff==0.14.14 pycodestyle==2.14.0
make -C v2 test
make -C v2 lint
```

Use the machine-wide performance lock for full suites and accepted timings.
See [tests](tests/README.md), [benchmark runners](benchmarks/README.md) and
[contributor guidance](../AGENTS.md) for checks and experiment rules.

## Status, results and history

- [Current status](STATUS.md): defaults, work in flight and eligible work.
- [Roadmap](ROADMAP.md) and [policy](POLICY.md): open gates and stable rules.
- [Benchmark results](benchmarks/RESULTS.md): scoped measurements and decisions.
- [Changelog](../CHANGELOG.md): accepted behavior changes.
- [History](history/README.md): pinned earlier documentation and evidence.

Generated captures, profiles and scratch stay in ignored `results/` folders.
Research and detailed journals remain local under `audit/` and `local/`.

### Experimental bounded ECM workers (B12/F4 preparation)

`v2.execution.ecm_workers.ECMExecutor` is an opt-in context for independent
curves using explicit `spawn` and PyPy implementing Python 3.11. This preparation
has not passed runtime arithmetic, multiprocessing or complete-checkout QA and
has no performance result. Automatic routing and worker defaults are unchanged.

The component accepts a resolved `PortfolioConfig`, a finite `Budget`, an
original/cofactor lineage string and an optional drained checkpoint. It returns
a checked divisor, reason, checkpoint and resource receipt. `WorkerSettings`
exposes 1/2/4 workers and finite work, CPU, memory, frame, checkpoint, startup and
drain limits. Only python-int is admitted initially. Full schedule eligibility
is retained in each worker; private program/chain stores are reused and bounded.
The caller owns `close()` or the context manager lifetime. One call may be
active per context. Process entry points must be protected by a main guard.

`factorize_bounded(..., ecm_executor=context)` is the narrow portfolio seam.
Its worker parent-memory reserve must include the literal portfolio allowance.
The context is borrowed; the caller closes it. Omit it for incumbent behavior.
No stock CLI worker flag is added. The experimental benchmark driver provides
fresh CLI measurements after runtime admission.

`EngagementRule` requires a matching schedule cost key, finite per-tier serial
costs and cold/warm thresholds. Missing/invalid calibration or unfunded startup
selects equivalent-manifest serial execution. `force_workers=True` is a separate
diagnostic arm. Neither supplies calibrated automatic settings. Rule estimates
use observable plans and pool readiness, not factor labels or future outcomes.

The opt-in search has versioned per-curve seeds independent of slot/count/order;
it differs from the incumbent shared RNG stream. ECM draws no longer move later
fallback seeds within the opt-in semantics. Sigma collisions are retained and
reported; IDs do not prove statistical independence. First-factor winners may
race. Checkpoints bind executor/rule/settings (including worker count), full
manifest, attempts, acknowledged results and cumulative charges. Resume rejects
a missing/different executor. Only drained checkpoints are supported; unfinished
curves restart setup with the same ID and retained prior charges, up to the
finite attempt cap (at most three attempts/curve and 32 drained runs per
cofactor search). These bounds also fund retained lease history. A parent-memory acknowledgement becomes durable only when
the caller atomically saves the returned checkpoint.

Missing reports consume full leases. Work units, observed CPU and conservative
CPU reservations are distinct. Shutdown CPU is reserved in full before spawn;
process RSS is an observation, not an owned-payload theorem. Cancellation polls
at existing stage actions; bigint operations and OS waits are cooperative.
Private channels are drained, then owned processes are terminated/killed with
finite joins. Runtime behavior, response tails and resource estimates remain
queued for validation. A zero/unfunded output allowance may refuse checkpoint
construction explicitly. No arithmetic state or live process handles are saved.
