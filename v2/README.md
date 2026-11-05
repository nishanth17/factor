# Factor v2

The active Factor implementation runs on **PyPy implementing Python 3.11**.
It uses exact integer arithmetic, validates every split and retains unfinished
cofactors. The default Python-integer backend needs no third-party packages.

## Setup and first run

Run commands from the repository root:

```sh
make -C v2 runtime
pypy3 -m v2.factor 626100403 --seed 7
pypy3 -m v2.factor 626100403 --seed 7 --ecm-workers serial
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

### Automatic bounded ECM workers

Normal `python -m v2.factor N` and `factorize_bounded(N)` now select cost-gated
process ECM automatically on PyPy 3.11 / macOS arm64. The packaged policy uses
actual cofactor size, exact configuration, remaining work and pool readiness.
It usually selects ten workers for the 140-curve plan; other declared regimes
retain their measured nominees. Worker IDs do not identify physical core types.
Primes, preprocessing exits, unsupported regimes/backends, nested processes and
unfunded calls stay serial. Catalogue/source incompatibility declines before
starting a new search. No external factorer is used.

Omitted numerical tiers select escalating140 at exactly 60/80/100 original
input digits, alongside the existing 30/40 choices. Explicit tiers and saved
plans win. Other sizes retain their existing numerical schedule. The 415-curve
hard pilot is not a default. SIQS routing, checkpoint limits and fallback
configuration are unchanged.

Use `--ecm-workers serial` or `factorize_bounded(..., ecm_executor=None)` to
select serial execution. Explicit contexts remain borrowed and caller-owned;
automatically created pools are closed on return and on errors. The legacy
`factorize()` compatibility helper remains serial; the bounded portfolio is
the default CLI and the automatic library entry point. `--bounded` still works.

Fresh omitted native non-QS grants are finite: 30 wall seconds, 4,096 aggregate
CPU seconds and 10^13 work units. Existing calibrated QS grants remain. Fresh
serial/automatic CLI comparisons receive the same omitted grants; explicit
Budgets and CLI limits are never increased. Missing finite wall/CPU limits,
very small startup grants or insufficient owned memory select serial behavior.

Omitted memory permits up to 8 GiB aggregate ECM ownership, including the
parent, all workers, spawn helper, metadata, transport and checkpoints. This is
an allowance, not an allocation or RSS claim. An explicitly supplied parent
workspace cap also limits automatic aggregate ownership unless a separate
`ecm_memory_bytes` / `--ecm-memory-mib` allowance is given. Parent/fallback
workspace is never silently reduced to fit workers. `--memory-mib` keeps its
existing workspace meaning; a tight cap can therefore select serial. Explicit
aggregate requests above 8 GiB are capped at the packaged policy's 8 GiB limit.

Automatic checkpoints use outer schema 13 with a checked mathematical
`portfolio_version` subschema. Older readers reject them. They bind
`b12-f4-auto-v1`, the packaged policy, implementation,
method seed, aggregate cap and native drained state. Legacy serial checkpoints
stay serial. Resume cannot change the automatic policy/cap/count or switch to
an explicit serial override. Integer method seeds retain the measured mapping;
other supported seed types derive a stored stream without advancing parent RNG.
Close/drain occurs before the final automatic snapshot/resource clocks.

For retained library use, `execution.ecm_auto.selected_context(...)` returns a
caller-owned context using the packaged policy. Use a context manager and pass
it as `ecm_executor`; `prefer_warm=True` retains the explicit 64-compatible-call
intent, while the first call still pays cold startup. Supply finite shared
per-call budgets that fund the selected pool. A context keeps its first engaged
worker configuration; it does not resize for arbitrary mixed workloads.

Adoption is explicitly authorized on 5 November despite the pending balanced
60/70 completion guardrail, censored p90, some cold median regressions and higher
CPU. The earlier fresh bounded comparison supports completion-sensitive gains;
it does not measure the new default wiring or every production grant. See
[results and the adoption exception](benchmarks/RESULTS.md#b12f4-automatic-default-integration-5-november-2026).

### Bounded ECM executor details

`v2.execution.ecm_workers.ECMExecutor` is an opt-in context for independent
curves using explicit `spawn` and PyPy implementing Python 3.11. This preparation
passes runtime arithmetic, multiprocessing and committed-checkout QA on the
current PyPy/macOS build. Initial scaling observations remain exploratory. The automatic adapter above
selects explicit settings; low-level WorkerSettings defaults remain unchanged.

The component accepts a resolved `PortfolioConfig`, a finite `Budget`, an
original/cofactor lineage string and an optional drained checkpoint. It returns
a checked divisor, reason, checkpoint and resource receipt. `WorkerSettings`
exposes 1/2/4 workers and finite work, CPU, memory, frame, checkpoint, startup and
drain limits. Only python-int is admitted initially. Full schedule eligibility
is retained in each worker; private program/chain stores are reused and bounded.
Funded workers start before the parent collects readiness reports, so interpreter
imports may overlap. Every child retains its own startup deadline and lease;
unreported startups keep their full charges if a peer fails. Startup costs are
recorded separately from useful execution.
The caller owns `close()` or the context manager lifetime. One call may be
active per context. Process entry points must be protected by a main guard.

The low-level executor remains explicitly configurable; the automatic wrapper
above selects it by default. `factorize_bounded(..., ecm_executor=context)` is
the narrow caller-owned portfolio seam.
Its worker parent-memory reserve must include the literal portfolio allowance.
The context is borrowed; the caller closes it. Omit it for incumbent behavior.
The stock CLI exposes automatic/serial selection as described above. The experimental benchmark driver provides
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
no decoded cache and whole curves. These remain opt-in choices with
source-pinned count/parameter measurements. Full-tier PRAC eligibility remains
intact; extra caches are separately funded rather than
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
`force_workers=True` is a separate diagnostic arm. The optional measured
catalogue below supplies source-specific calibration.
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

Native SEARCHv1 execution does not move the portfolio generator with its ECM
draws. Selected portfolio calls that decline workers instead preserve incumbent
serial curves and the incumbent fallback RNG consumption. Sigma collisions are
retained; distinct IDs do not prove independence. First-factor
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
actual ECM is reached. A bounded, source-checked profile selects a configuration
at the first engaged cofactor; declined calls do not pin a serial context.
The worker count stays fixed after engagement. Decisions use decimal size,
exact mathematical options, remaining quota and actual pool readiness. `prefer_warm=True` chooses the retained-library
candidate, while its first call still pays and checks cold startup. Unknown or
invalid cost estimates return portfolio execution to incumbent serial ECM.
Engaged searches retain the explicit SEARCHv1 stream; production serial remains
the incumbent comparison. No automatic route or default is changed.

The profile is limited to 256 KiB and 256 rows. The parent separately funds
4 MiB decoded ownership, bounded input/hash scratch and three 16 KiB checkpoint
headers alongside native executor ownership. Source hashing streams small
buffers and binds shared runtime imports, including the QS fallback. Header work is reserved before arithmetic. Drained checkpoints bind
the profile, selected settings, native checkpoint and prior charges; changing
policy or worker count is rejected. Fresh contexts repay profile preparation
and process startup. Retained contexts reject incompatible ownership rather
than resizing their pools. The optional measured catalogue below is a separate explicit input. It implies
no automatic adoption; runtime QA and configuration confirmation retain their
separate scopes.


### Measured opt-in ECM selection

`v2.execution.ecm_selection.SelectedECMContext(profile, profile_sha256, seed=0,
prefer_warm=False)` loads a bounded, source-checked catalogue lazily when ECM is
reached. Pass it as `factorize_bounded(..., ecm_executor=context)` and close it
with a context manager under a protected process entry point. The historical
training catalogue is `benchmarks/inputs/controls/b12_f4_selection_v2.json`, SHA256
`2c87e95c2cd74bee564adb093a2df71ecb9c45bd64dfb480a86dbcd31b0e5e29`.
It binds 99 runtime files and retains the 33 measured cost/settings anchors from v1. Their costs are transferred estimates for this changed adapter,
not measurements of its current complete cost. Source incompatibility is rejected.
The frozen v1 input remains immutable for its original sources; current selector
profiles/checkpoints use `b12-f4-selection-v2` and reject v1. Native worker wire v7
and mathematical SEARCHv1 identities are unchanged. No automatic routing or stock
CLI change follows.

Cold intent requires expected remaining serial ECM work above both a fitted
complete-cost crossover and a one-second minimum, with a nominal dispatch-wave
floor. Separate thresholds cover a started pool and compatible retained metadata.
`prefer_warm=True` selects an explicitly declared horizon of at least 64
compatible ECM-reaching calls. Actual first-call startup remains cold and fully
charged. Short or uncertain reuse should use cold intent. A context retains its
first engaged worker settings; mixed-regime pool resizing is not implemented.
Missing exact sizes use piecewise linear per-curve cost interpolation only
between same-plan anchors at most 20 decimal digits apart. Nearest-anchor ties
choose the lower size. Pool crossover estimates remain bound to the selected
bundle, with a recomputed dispatch-wave floor. PRAC/ladder eligibility boundaries
partition cost anchors; no extrapolation or interpolation across that boundary
is admitted. Decision receipts identify measured versus interpolated rates and
the rate/crossover anchors. A separate targeted calibration validates
55/75/90-digit fixed-work body estimates within 4–17%; other intermediate sizes
remain unmeasured.

For `factorize_bounded`, missing/invalid costs, fresh worker-memory refusal and
unfunded startup select the existing serial stage with its original curves,
programs and RNG consumption. A per-cofactor checkpoint flag keeps that decision
sticky across pause/resume. Once workers are selected, later refusal pauses the
same SEARCHv1 search; it cannot switch streams. Unknown serial calls do not pin
pool settings for a later eligible cofactor. Omit the context for unmodified S0.

`SelectedECMContext.run()` remains the direct component/fixed-work API and retains
SEARCHv1 serial replay. `run_portfolio()` is the narrow caller seam, returning
`None` on an incumbent-serial decision. Component and portfolio modes cannot be
mixed in one context; drained selector checkpoints bind the mode. Old-source
performance observations do not establish gains for this revised adapter.

Full constructor snapshots for reproducing the catalogue's mathematical regimes
live in `benchmarks/inputs/controls/b12_f4_selection_resolved_v1.json`. These
snapshots define measured caller configurations rather than new portfolio
presets. In particular, caller140 at 60/80/100 digits is explicit.

The earlier comparison of source `538462a` uses actual master stock QS
checkpoint/fallback settings, separating default-schedule and explicit caller140
arms. It finds conditional
completion gains and regressions, with no universal cold/warm median winner.
That source matched only exact cofactor sizes/options, leaving many calls on
its distinct SEARCHv1 serial stream. The v2 adapter above changes those two
integration behaviors. Its subsequent frozen comparison is linked below. A
nominated worker count therefore does
not mean workers engaged on every call. The earlier matched 64 MiB QS checkpoint
study remains a historical resource variant. See
[actual-master results](benchmarks/RESULTS.md#b12f4-actual-master-complete-cost-comparison-5-november-2026)
for population sizes, uncertainty, CPU/memory costs and pending acceptance gates.

A separately frozen opt-in profile,
`benchmarks/inputs/controls/b12_f4_targeted_selected_v1.json`, has SHA256
`16d52a6e55270d3cfea467f2a5ca41dd00ad1d8d052f2263861656a99b4b00cf`.
It adds one calibrated 81-digit anchor while retaining the earlier 33 bundles.
Use it explicitly with `SelectedECMContext`; it does not change global defaults.
The fresh comparison requests escalating140 at original 60/80/100-digit inputs
and preserves stock fallback/checkpoint settings. Cost admission uses actual
cofactor size and pool readiness, with incumbent serial execution on decline.
It improves completion-sensitive mean latency on the declared fresh uniform
panel while increasing CPU and regressing some cold medians. Warm results use
retained homogeneous contexts; they do not establish a cold CLI gain by
themselves. See the [fresh policy results](benchmarks/RESULTS.md#b12f4-fresh-cofactor-policy-comparison-5-november-2026)
for grants, uncertainty, refusals and the pending balanced completion guardrail.

Current automatic use loads the packaged `execution/ecm_worker_policy.json`
through `execution.ecm_auto`; no benchmark directory is needed at runtime. It
retains the measured 34 rows and new integration source pins. Earlier benchmark
profiles remain immutable for their pinned sources and reject changed code.
For current explicit retained contexts, use `selected_context()` above rather
than an old source-pinned benchmark profile.
