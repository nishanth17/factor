# QS and SIQS configuration

Use the [package guide](README.md) for jobs and relation invariants.
[STATUS](../STATUS.md#current-automatic-defaults) is the single current
size/resource routing table, including adopted and unmeasured bands.
The selectors are [defaults.py](defaults.py) and [siqs.py](siqs.py).

## Explicit and automatic selection

`SIQSJob(n)` can select a size-scoped bundle. Supplying `SIQSConfig(...)`
uses its explicit reference fields. `default_siqs_config` exposes automatic
selection; `CountedSIQSConfig` adds a positive eligible `base_count` below its
finite endpoint ceiling. Base resolution and restore are charged.

Explicit settings, CLI parameter overrides and saved concrete configurations
win over automatic choices. `asdict`/`replace` copies freeze the selected
values. Direct QS/MPQS controls and other backends retain their documented
scope. A finite default grant does not promise completion.

## Collector families

| Configuration | Role |
| --- | --- |
| `SieveConfig` | Single-large-prime reference and scored collection |
| `DoubleLargeSieveConfig` | Bounded DLP graph with fixed legacy quotas |
| `GrowingDoubleLargeSieveConfig` | DLP with allocator-grown stores and optional quotas |

A DLP configuration separately bounds `residual_bound` for SLP admission,
`large_prime_bound` for endpoints, `large_product_bound` for products,
`candidate_bound` for recovery intake and `split_call_limit` for splitting.
Endpoints must satisfy the exact certification contract; a failed split is
not a primality proof. Full original exponents and square corrections survive
self-loops, repeated primes, parallel edges and disconnected cycles.

Capping unowned partial edges does not cap pinned atoms, combined rows or
matrix storage. Lowering the SLP ceiling can deliberately lose otherwise
valid SLP matches; it is a retention policy, not an arithmetic shortcut.
DLP CLI use requires explicit endpoint/product bounds when manually selected.
SSS and coarse-worker exporters reject unsupported DLP configurations.

Growing stores have no implicit FIFO or count stop when their quotas are
`None`. They still use finite job allowances and per-residual effort. Explicit
outer portfolio memory constrains the child. Integer quota fields restore
finite application-specific stops; a larger store alone does not establish a
complete-factor benefit.

## Assignment and stopping

`mode="qs"`, `"mpqs"` and `"siqs"` select polynomial construction. SIQS
`assignment_policy` supports bounded reference, nearest and flyer streams;
`family_count` and `polynomials_per_family` govern finite search and consumed
Gray prefixes. `assignment_policy="quality"` is an explicit fixed-batch
actual-A ordering, not a proven diversity policy.

`half_width` and `max_half_width` must agree with the intended interval.
`factor_count` and `pool_size` affect reachable A products and family size.
Configuration/resume identity prevents changing a consumed schedule silently.

`filter_row_growth` and `row_excess` control repeated dependency work.
`solve_policy="collection"` separates a polynomial window from genuine shared
collection exhaustion, while preserving immediate checked square dependencies
and a charged terminal attempt on retained originals. `"window"` retains
per-window finalization. Pending extraction stays pending across resource
refusal; a budget extension cannot replace its saved policy.

## Scores and marking kernels

`score_policy="fixed"` uses conservative integer bounds in scaled units;
exact exponent recovery and relation verification decide admission.
Positive `threshold_extra` can deliberately lose candidates. Byte/array
saturation stays conservative with its threshold.

`power_plan_bytes` reserves bounded per-polynomial modulus/root/weight plans.
Construction and replay are charged; cap refusal falls back to streamed
lifting. Polynomial/span changes invalidate derived plans. Saved caches are
never accepted as verification evidence.

`marking` supports `density`, `presieve`, `presieve_powers`, `unroll4` and
`scan_fused` in valid duplicate-free combinations. The two presieve selectors
are mutually exclusive. Kernels specialize supported list-backed score/recovery
modes and retain exact scalar paths for unsupported or exceptional cases.
Scratch coexistence, rotation, cancellation and refused work stay bounded.
The selected-power/unroll extension remains opt-in after its combined loss;
see [A1a evidence](../benchmarks/RESULTS.md#a1a-integration-and-extension).

## Adjacent signs and sparse queues

Automatic SIQS uses `sign_policy="swap_high"`, an adjacent Gray-sign ordering;
explicit reference configs retain `"canonical"`. Saved schedules preserve their
identity. Nonadjacent jumps cannot reuse a one-sign update.

`sparse_queue_bytes` is a ceiling for exact per-polynomial prime/power hit
queues inside the existing owned allowance. Unsupported scores/recovery,
singular sparse roots, unaligned blocks or capacity refusal use the ordinary
exact path. All hits and candidate checks remain covered. Polynomial changes
release queues; interruption/restore rebuilds them with charges.

This selection was [adopted](../benchmarks/RESULTS.md#c10a-07-face-and-queue).
Its exploratory performance label and unmet drift/confirmation gates remain
attached to the evidence. Other scope extensions are not confirmed gains.

## Full parameter and historical reference

Configuration dataclasses and validation live in [siqs.py](siqs.py) and
[sieve_collector.py](sieve_collector.py). Versioned selected study bundles live
under [benchmark controls](../benchmarks/inputs/README.md). Load a whole frozen
bundle when reproducing its result; partial parameter copies change its scope.
[Earlier API details](https://github.com/nishanth17/factor/blob/24ff153d26ada2ec4beeead07000d183e94fd355/v2/README.md)
retain historical constants, schema narratives and example study loaders.
