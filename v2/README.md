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

## Current QS workers and selected defaults

`v2.qs.current_parallel.CurrentParallelSIQSJob` accepts a concrete current
`SIQSConfig` or `CountedSIQSConfig`, a root `Budget`, and `QSWorkerConfig`.
Growing DLP and the C9a cost-report configuration are supported. The interface
preserves width, multiplier, parity/sign, roots, queues and exact recovery,
with one authoritative parent relation/partial/cycle store.

`QSWorkerConfig(workers=1)` through `workers=20` select spawned PyPy3.11
processes; counts above the host's ten cores are explicit oversubscription
experiments. `mode="serial"` executes the same packets locally. Pass an explicit
`CollectionPool` to retain processes across jobs. The portfolio accepts
`factorize_bounded(..., qs_workers=settings, qs_pool=pool)` at SIQS handoff,
including recursive cofactors. Omitted `qs_workers` selects the digit-based
process recipe below; explicit `qs_workers=None` forces serial S0. Retained ECM and QS
processes share the root's observed wall, aggregate CPU and current-tree RSS.
There are no parent/worker memory reservations, provenance/history quotas,
prepaid work/CPU leases or aggregate logical-byte admission checks. Former
resource keywords raise `RemovedResourceLimitError` naming the argument.

Native automatic QS uses the Y0 numerical policy measured by S0 and the
worker calibration. Routing uses the original call's decimal input size and
the nearest row; exact ties choose the smaller size. Below 40 remains serial.

| Anchor digits | Workers / polynomials per task |
| --- | --- |
| 40 | C4/p4 |
| 50 | C4/p8 |
| 55 | C4/p16 |
| 60 | C8/p16 |
| 65 | C14/p64 |
| 70 | C8/p64 |

Above 70 inherits C8/p64: **extrapolated**, not validated larger-size
performance. Numerical policy retains Y0's existing 62/75/80 definitions
without adding calibration rows or running probes. Intermediate 50+ numerical
recipes inherit the nearest defined Y0 geometry; above 80 inherits 80. Existing numerical settings below 50,
non-native backends, explicit constructors and `recipe="reference"` remain
available. Automatic portfolio SIQS eligibility is unchanged; explicit
`--siqs` or `--method siqs` selects the engine outside those bands.

The selected worker bundle retains 64 MiB state cache, lookahead 4, packet
transport, worker verification, coefficient caching and private semaphore
stops. It keeps max_batch_atoms 2048 and all frozen capacities. Pools start
only when SIQS is reached. Worker timelines and full profiling stay off.
Use `--qs-workers serial` at the CLI or
`factorize_bounded(..., qs_workers=None)` in Python to force serial QS with
the same numerical policy. Saved automatic calls retain their concrete
worker recipe; saved serial calls stay serial. Explicit worker checkpoints
continue to require their original settings. Changing execution mode mid-resume
is refused. The checkpoint's `qs_routing` marks above-70 extrapolation.

Root wall/CPU limits are optional; RSS defaults to 80% of physical RAM.
For an explicit guard, use `Budget(rss_bytes=18 * 2**30, seconds=...)`.
OS observations include process startup and final CPU at reaping. Reported
work remains diagnostic; missing private reports leave its amount unknown
rather than inventing a full lease. Root limits are cooperative and may
overshoot between observations or indivisible arithmetic operations.

Configuration and execution settings remain fixed across resume. Version-3
current-worker checkpoints carry the H3 resource version, observed saved
history, published prefixes, pending admission, provenance and numerical
split-call allocations. Old resource schemas are rejected. Checkpoint encoding
has no internal size cap; optional file persistence uses
`atomic_save(..., checkpoint_bytes=...)`, which records an oversized-save skip
without discarding the factoring result. Derived solver/root/proof caches
rebuild with their actual execution cost. Width growth, nested recovery and
native threads remain unsupported by this experimental interface.

Scheduling and representation options remain explicit:

- `polynomials_per_task` selects a packet's original Gray prefix.
  `whole_family=True` groups a configured family without adding polynomials.
  `lookahead=N` adds at most four queued/pending credits per worker. Atom and
  wire capacities continue the same unpublished cursor when a packet fills.
- `reuse_state=True` retains an identity-bound context and reusable family and
  collector scratch. `state_cache_bytes` controls retention; an oversized
  context uses uncached packets. Recovery hints and polynomial queues are
  invalidated at their original boundaries. Missing installation acknowledgments
  retire only the owned pool. Retained pools pin worker settings.
- `lease_budget=True` retains the experimental packet-local diagnostic counter
  and 64-action root-cancellation cadence. It no longer introduces resource
  leases. `grouped_lease_charging=True` combines two marking-counter updates
  only within that polling interval and remains disabled by default.
- `packed_transport=True` publishes checked varint rows with private reusable
  buffers, unchanged numerical publication boundaries and strict checkpoint
  base64. `verification="worker"` verifies before publication; the parent still
  checks ordered insertion. `verification="batch"` prevalidates bounded groups.
  A bare `QSWorkerConfig()` uses full parent verification; the selected
  automatic bundle uses worker verification.
- `private_stop=True` gives each process its own stop flag and requires the
  packet-local counter. `stop_backend="semaphore"` replaces each Event with
  an exclusive-reader nonblocking semaphore. Both preserve cooperative polling.
  `cpu_snapshot=True` selects a private coherent counter snapshot; these
  diagnostics do not add worker CPU a second time to the observed root ledger.
- `family_affinity=True` pins a family to one process and requires retained
  state; use the same option on `CollectionPool`. Small credit windows can leave
  workers idle. `arrival_order=True` is a separate trajectory-changing option,
  rejects fixed-work execution and validates its arrival witnesses on restore.
- `coefficient_cache=True` retains a bounded parent CRT template during merge.
  `packed_queue` stores sparse events in 64-bit words; unsupported shapes use
  exact direct marking. `byte_stamps` selects checked byte epochs with charged
  wraparound clears. These mechanisms remain disabled unless selected.
- `wide_hits`, `wide_inverses` and `wide_transfers` enlarge the corresponding
  real cache capacities. `pruned_transfers` replays a complete cached root pair
  after window pruning. Their exact arithmetic and fallback paths are unchanged;
  removed total-memory estimates no longer select or refuse these caches.

The recorder entry point remains `pypy3 -B -m v2.benchmarks.qs.a7a.runner`.
Frozen pre-H3 recipes are translated only at the benchmark boundary; their old
resource fields remain historical metadata. Production constructors and
checkpoint restoration do not accept those fields. Historical G2H/G3 timings
need revalidation on the combined H3 source. No automatic worker selection or
performance adoption follows from this integration. Detailed local evidence
and handoff remain under `v2/audit/a7a/`.

## Experimental pthread ECM owner on B4a

The B4a branch adds `pypy3 -B -m v2.ecm.native.build --gmp-prefix PREFIX
--pthread` alongside the separate ordinary-GMP scalar build. Both retain O2.
`NativeECMExecutor(..., settings=WorkerSettings(workers=1|2|4),
schedule="dynamic"|"static", chunk_size=1..32)` reuses contexts/programs and
retains one fixed C-only pthread team across drained continuation-slice
generations, joining on close. Cheap wall reads remain at each action; process
CPU sampling is amortized to slice entry and one-millisecond wall intervals.
Checked publication, root wall/process-CPU/RSS and saved-history accounting
remain authoritative; native-thread CPU is included in process CPU once.
`NativeTeam.cancel()` cooperatively cancels active native work; close is
idempotent, active close/use-after-close reject, and close/cancel must not race.
The new `b4a-pthread-owner-h3-v3` checkpoints bind source, settings and chunk
size; old OpenMP checkpoints do not migrate. Scalar/pthread correctness is
validated initially on `0f2fee3`; the retained `c-gmp` stage-one control
keeps setup, classical stage two and recovery in Python. The separately selected
whole-ECM C owner is described below. See [native ownership and validation](ecm/native/README.md). No automatic
routing, performance adoption or mainline integration follows.


The B4a whole-ECM candidate is selected with
`PortfolioConfig(ecm_kernel="c-gmp-full")` or `--ecm-kernel c-gmp-full`.
It preserves scalar durable checkpoints and the persistent team's assignment
resume contract while executing setup, both ECM stages, GCD and strict recovery
in C. Python prepares immutable metadata and validates/publishes results.
The existing `c-gmp` stage-one executor remains a control; automatic routing
is unchanged. See [whole-native ownership and eligibility](ecm/native/README.md).

Explicit native builds prepare hash-checked CFFI ABI descriptors; first use
retains source/binary/GMP validation while avoiding declaration parsing.
Childless native resource polls use the same macOS observation APIs, with
the full tree provider for children or custom resource views. Compiler options
and the disabled limb-division experiment are documented in the native guide.
The separate `--coalesce-unit-yields` build experiment combines post-unit-GCD
yields within the same bounded slice while retaining checks and publication
ordering. It is disabled by default pending continuation validation.
These component experiments do not change automatic factoring selection.
