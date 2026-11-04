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
passes runtime arithmetic, multiprocessing and committed-checkout QA on the
current PyPy/macOS build. Initial scaling observations remain exploratory. Automatic routing and worker defaults are unchanged.

The component accepts a resolved `PortfolioConfig`, a finite `Budget`, an
original/cofactor lineage string and an optional drained checkpoint. It returns
a checked divisor, reason, checkpoint and resource receipt. `WorkerSettings`
exposes 1/2/4 workers and finite work, CPU, memory, frame, checkpoint, startup and
drain limits. Only python-int is admitted initially. Full schedule eligibility
is retained in each worker; private program/chain stores are reused and bounded.
Funded workers start before the parent collects readiness reports, so interpreter
imports may overlap. Every child retains its own startup deadline and lease;
unreported startups keep their full charges if a peer fails. This is unmeasured.
The caller owns `close()` or the context manager lifetime. One call may be
active per context. Process entry points must be protected by a main guard.

`factorize_bounded(..., ecm_executor=context)` is the narrow portfolio seam.
Its worker parent-memory reserve must include the literal portfolio allowance.
The context is borrowed; the caller closes it. Omit it for incumbent behavior.
No stock CLI worker flag is added. The experimental benchmark driver provides
fresh CLI measurements after runtime admission.

`WorkerSettings` additionally exposes bounded optimization choices:

- `packet_curves=1..4` batches commands while every job retains its own lease
  and streamed result. Admission shares ready work/grants among idle workers.
- `poll_every=1..64` controls clock/socket checks; work is checked on every
  reservation. Stage changes, handoff and publication force real checks.
  `parent_poll_seconds` bounds result-wait polling to 1–50 ms (default 10 ms).
- `metadata_slots=1..2` retains modulus-independent contexts in a private LRU.
  Optional `decoded_cache_bytes` retains immutable prime/power tuples, capped
  at 16 MiB/context and 64 blocks; refusal uses ordinary packed programs.
- `stage_policy="staged"` enables a bounded clean-boundary residue queue;
  `residue_capacity` includes in-flight stage-one jobs. Available stage-one
  trials have priority within a tier; continuations run when capacity fills.
- Workers outside 1/2/4, up to sixteen, require `allow_wide_workers=True` and funded
  ownership. This adds no wider-worker experiment allocation or automatic count.

Defaults remain one job/packet, polling every action, one metadata context,
no decoded cache and whole curves. These remain opt-in candidates; only the reference worker configuration has
initial scaling observations. Full-tier
PRAC eligibility remains intact; extra caches are separately funded rather than
shrinking arithmetic ownership. Metadata replaces/evicts before new allocation;
curve points never enter its cache. Every new process begins with empty caches.

`EngagementRule` accepts one row or a bounded `cost_table` of `CostRow` values.
Rows bind `Manifest.schedule` and the matching executor's `cost_identity` in
`execution`; missing/mismatched execution keys select serial. Literal committed
source/configuration pins are also required by the experiment freeze. Each row
provides finite per-tier serial costs and cold/warm thresholds; optional
`pool_seconds` distinguishes an already spawned pool without matching metadata.
Staged remaining work additionally needs `seconds_per_stage_two`. Missing,
invalid or unfunded estimates select equivalent-manifest serial execution.
`force_workers=True` is a separate diagnostic arm. No calibration is supplied.
The decision receipt reports cold/pool/warm state, estimate and threshold.

The mathematical manifest/curve stream remains `b12-f4-v1`, independent of
slot/count/order; it differs from the incumbent shared RNG stream. The extended
wire/executor/checkpoint contract is `b12-f4-v7`; v1–v6 checkpoints are rejected.
Pinned older sources remain available for those checkpoints. Version 3 fixes
JSON position-key hashing and withholds CPU for parent progress. Counts above
ten explicitly oversubscribe this M4; support does not establish a useful count.
Optional `packed_cache_bytes` (0–16 MiB) adds lazy compact prime/power retention
independently of chain eligibility. It backs streamed bounds or adds to an
existing program store without replacing accepted chains. An unfunded scratch
cap refuses before dispatch. It defaults to zero and adds separate ownership
for every retained context; its cold construction cost must be measured.

Decoded and stage-two index capacities may be explicitly funded up to 64 MiB
per retained context. Packed and manifest caches retain their 16 MiB maxima.
Larger requests must still fit simultaneous parent/worker/helper ownership;
resource/checkpoint v7 rejects v6 checkpoints and cost rows. Mathematical v1
curve identities and ordinary zero-cache defaults are unchanged.

Optional `decoded_cache_bytes` now reuses immutable primes/powers for streamed
bounds as well as packed programs. A funded prefix remains cached when a full
schedule cannot fit; there is no repeated scan eviction. `manifest_cache_bytes`
(0–16 MiB) separately bounds up to 16 parent identity manifests, with fresh
search/results on every call. Both default to zero and refuse retention safely.
`worker_imports="ecm"` uses the same configuration/context implementation while
avoiding unrelated relation-engine imports. `stage_two_kernel="reduced"` selects
the already-verified early-reduction point kernels for baby/giant steps. Defaults
remain `"portfolio"` and `"legacy"`; numerical schedules, curve streams, stage-one
eligibility and modulus-dependent state ownership stay unchanged. These choices
require separate complete-cost calibration; they are not automatic defaults.

Optional worker settings `budget_mode="prospective"`, `bind_context=True` and
`stage_two_cache_bytes` (0–64 MiB) are headroom candidates, with defaults
`"legacy"`, `False` and zero. The prospective budget preserves every exact work
check and the existing clock/cancellation cadence. Binding retains compatible
schedule references once per context, never curve state; 4 KiB per context is
separately reserved. The stage-two cache owns a funded prefix of immutable baby
indices and giant-group ends, with `4096 + 256 * segment_size` scratch and a
4,096-block ceiling. Unfunded scratch refuses before dispatch; full retention
falls back to the original continuation. Compilation checks cancellation every
256 primes and publishes only complete blocks. Each giant action, term order,
GCD batch and recovery path is preserved. Direct-scalar/paired/wheel paths keep
their existing implementation. Global serial execution is unchanged.

`packet_policy="tail_one"` optionally uses singleton packets when remaining
ready quota is at most twice the worker count; `"fixed"` remains the default.
Packet sizes 1–4 retain per-job leases and bounded frames/results. New stage-two
hit/miss/refusal counters and `stage_two_build_bytes` are cumulative preparation
observations. Build bytes include repeated allocations after context replacement;
they are not live memory or aggregate RSS. All flags bind the v7 executor/cost
identity; mathematical v1 manifests and curve streams remain unchanged. These
options need matched complete-cost evidence and are not adopted defaults.

ECM draws do not move later fallback seeds within the opt-in semantics. Sigma
collisions are retained; distinct IDs do not prove independence. First-factor
winners and staged search order can differ legitimately. Fixed-work runs retain
the same declared curves and existing exact stage arithmetic.

Checkpoints bind settings/count/rule, manifest, stage attempts, acknowledged
residues/results and cumulative charges. Worker-count changes remain rejected.
Only drained checkpoints are supported: unfinished whole/stage-one jobs restart
setup; unfinished continuations restart from the checked pre-table point. The
attempt cap is per curve-stage (maximum three), with at most 32 drained runs per
cofactor. Canonical residues bind identity, sigma, bounds, point and curve
parameter with a digest and range/unitness checks. These checks detect malformed
transport/state; they do not independently recompute stage one. Arithmetic
oracles and selected real-process cases pass; broader performance gates remain. Acknowledgements become
durable only when the caller atomically saves the returned checkpoint.

Missing reports consume full leases, including queued packet jobs. Work units,
observed CPU and conservative CPU reservations are distinct. Shutdown CPU is
reserved before spawn. Healthy fully drained pools can survive ordinary factor
cancellation; idle workers block. Closing broadcasts each stop phase and uses a
single pool-wide deadline for its joins. Bigint/OS waits and deferred polling
remain cooperative: response is bounded in charged actions plus an atomic
operation/wait, not a hard elapsed-time theorem. Owned memory includes retained
contexts, decoded caches, transport, residues and checkpoints; it is not an RSS
theorem. A possible runtime-managed spawn helper gets one additional worker
memory reserve and a full startup+shutdown CPU/work allowance per context. Its
CPU is explicitly unobserved in worker receipts; library close does not inspect
or stop shared multiprocessing internals. `checkpoint_bound` may refuse a large schedule before dispatch; increase
its explicit funded allowance rather than silently dropping history. A zero
output-work allowance may refuse checkpoint construction. No live process
handles or interpreter internals are serialized.


### Measured ECM selection context (experimental preparation)

`v2.execution.ecm_selection.SelectedECMContext(profile_path, profile_sha256,
seed=..., prefer_warm=False)` adapts the existing opt-in portfolio executor
interface. Use it as a context manager and pass it as `ecm_executor` to
`factorize_bounded`. Resolve an omitted configuration with
`PortfolioConfig().for_input(n)` before constructing explicit numerical settings;
this retains the production SIQS fallback. Use a finite shared `Budget`.

The adapter reads no profile and imports no multiprocessing machinery until
actual ECM is reached. A bounded, source-checked measured profile selects one
configuration at that first cofactor; the count stays fixed for the context
lifetime. Decisions use decimal size, exact mathematical options, remaining
quota and actual pool readiness. `prefer_warm=True` chooses the retained-library
candidate, while its first call still pays and checks cold startup. Unknown or
invalid cost estimates remain equivalent-manifest serial. This explicitly
versioned serial stream differs from legacy production ECM; production serial
remains the incumbent comparison. No automatic route or default is changed.

The profile is limited to 256 KiB and 256 rows. The parent separately funds
4 MiB decoded ownership, bounded input/hash scratch and three 16 KiB checkpoint
headers alongside native executor ownership. Source hashing streams small
buffers. Header work is reserved before arithmetic. Drained checkpoints bind
the profile, selected settings, native checkpoint and prior charges; changing
policy or worker count is rejected. Fresh contexts repay profile preparation
and process startup. Retained contexts reject incompatible ownership rather
than resizing their pools. No measured selection profile or adoption is implied
by this API; integrated runtime QA and configuration confirmation are separate
gates in the active selection campaign.
