# QS and SIQS configuration

Use the [package guide](README.md) for jobs and relation invariants.
[STATUS](../STATUS.md#current-automatic-defaults) is the single current
size/resource routing table, including adopted and unmeasured bands.
The selectors are [defaults.py](defaults.py) and [siqs.py](siqs.py).

## Explicit and automatic selection

`SIQSJob(n)` selects coupled policies with the existing size-scoped parameters. Supplying `SIQSConfig(...)`
uses its explicit reference fields. `default_siqs_config` exposes automatic
selection; `CountedSIQSConfig` adds a positive eligible `base_count` below its
finite endpoint ceiling. Base resolution and restore are charged.

Explicit settings, CLI parameter overrides and saved concrete configurations
win over automatic choices. `asdict`/`replace` copies freeze the selected
values. Direct QS/MPQS controls and other backends retain their documented
scope. A finite default grant does not promise completion.

## Coupled polynomial policies

Fresh `default_siqs_config(n)` and unconfigured SIQS jobs select
`polynomial_policy="parity2"`, `multiplier_policy="parity2_wide"`,
`multiplier=0`, `assignment_policy="windowed_pair"` and
`root_update_policy="batch"`. `parity_two` marking is composed with the
previous collector kernels. This applies at every supported size and backend;
performance outside the measured scope is unconfirmed. Size parameters,
owned allowances and finite grants are retained. The automatic portfolio still
introduces SIQS only at native 40 and 50–60 digits. Direct CLI `--method siqs`
or `--siqs` without QS parameter overrides uses these policies at other sizes
with the existing generic grants, not a newly tuned larger-size preset.

For `kN` congruent to 1 modulo 4, the parity family uses `A = 2a` for an odd
squarefree A product `a`, and odd centered `B`. The exact existing identity
`U(x)^2 - kN = A*F(x)` remains authoritative: `F(x) = A*x*x + 2*B*x + C`,
`C = (B*B-kN)/A` and `U(x) = A*x+B`. All integer x in the declared inclusive
window are admissible. Factor-base recovery retains the A contribution and
exact powers of two. Other discriminants, including even multipliers, use the
standard family. Multiplier ranking includes finite odd and even squarefree
candidates with representation-aware weights and size penalties; it may return
an exceptional proper GCD before collection. Supplying a positive multiplier
fixes k instead of invoking selection. Standard odd-only ranking remains
available with `multiplier_policy="standard"`.

`windowed_pair` searches a finite reachable product population with paired
completion; it can change the A population. `root_update_policy="cached"`
and `"batch"` preserve the selected physical trials while preparing bounded
modular Gray updates; batch construction is limited by the consumed family
quota. Refusal, cancellation and derived-state rebuild remain charged.
`parity_two` initializes guaranteed two-adic contributions only on supported
score paths and otherwise uses exact existing scoring. The narrower
`"parity2"` and norm-preview `"parity2_norm"` multiplier policies, and
`"paired_flyer"` assignments, remain explicit alternatives; norm preview is
not the selected default.

`SIQSConfig()` still constructs the literal standard reference policies.
Restoring a saved configuration never invokes today's automatic selector.
Checkpoint schemas 8/9/10 explicitly version the additional semantics; the
selected bundle uses schema 10. Incompatible representation, multiplier,
assignment or root-policy substitutions are rejected even when granting a
budget extension. Tables, queues and plans rebuild with their original policy
identities and charges. Legacy schemas retain their original arithmetic.
See the [promotion evidence and limits](../benchmarks/RESULTS.md#c10a-01-coupled-polynomial-defaults).

## Collector families

See [collector policies and storage](COLLECTORS.md) for the optional
useful-yield, expanded-intake and compact-record policy reference.

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
`assignment_policy` supports bounded reference, nearest, flyer, paired_flyer and windowed_pair streams;
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
`scan_fused`, `parity_two`, `packed_queue` and `wide_hits` in valid duplicate-free combinations. The two presieve selectors
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

The opt-in `packed_queue` component uses one checked 64-bit word per hit,
replacing three 32-bit words and the temporary append tuple. It preserves the
older conservative reservation, work/event order and exact direct fallback.
It follows `parity_two` in component order and changes no default.

`wide_hits` follows `packed_queue` and raises the optional native-word scratch
ceiling from 32 to 128 MiB. Its actual reservation, construction work and the
next polynomial's metadata coexistence must still fit the caller's grant.
Scarce resources retain exact large-integer masks. This selection is opt-in.

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
