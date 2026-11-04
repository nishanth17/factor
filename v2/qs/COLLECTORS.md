# Collector policies and storage

Detailed configuration and storage contracts for the optional QS/SIQS
collectors. Start with the [engine overview](README.md) and
[general configuration](CONFIGURATION.md) for usage and selection.

## Useful-yield policies


The built-in QS/SIQS pipeline consumes private collection progress without
copying the cumulative relation store after each batch. Public
`SieveCollector.collect()` still returns complete owned snapshots; callers
must release old snapshots to exclude those references from collector-owned
memory. Injected collector classes and overridden `collect` methods retain
their public dispatch. Checkpoints continue to encode the full checked store.

SIQS setup and polynomial switching share the existing 64-action cooperative
polling adapter. Every work reservation remains immediate; publication checks
poll explicitly, and callback changes reach the original budget. Switching
reuses checked factorization metadata only for the same supported A and
immutable factor-base primes, after validating all supplied roots. Sparse hit
queues use flat unsigned-word triples and residual incidence uses compact
private coordinates; checkpoint
formats and independent restore validation remain unchanged. These changes
are under measurement and do not establish a factoring speedup.

`qs.YieldSieveConfig` extends `GrowingDoubleLargeSieveConfig` with independent
experimental options. Pass it as `SIQSConfig(collector=...)`; automatic
configurations retain their existing policies.

| Field | Opt-in behavior |
| --- | --- |
| `preparation_policy="target-v1"` | Prepare at base-column count plus row excess; use bounded 5% top-ups after unsuccessful attempts. Immediate squares, equal-parity checks and final/storage-cap extraction remain available. |
| `large_product_bound=P`, `candidate_bound=C` | Admit recovered products up to P while holding the reporting allowance C fixed. `candidate_bound=0` instead derives reporting from the configured admission domain. |
| `small_prime_bound=S` | Require proven ordered endpoints p≤q with p≤S, q≤`large_prime_bound` and pq≤P, including equality and squares. SLP admission retains its separate `residual_bound`. |
| `cycle_policy="short-v1"` | Search all retained components for shorter exact paths, including earlier closing edges. `cycle_search_edges` bounds visits (default 4096, maximum 65536); the original forest path remains the fallback. |
| `duplicate_policy="pair-v1"` | Check equal-parity two-row dependencies before ordinary preparation. Keep both rows; bound representatives by `duplicate_bytes` (default 16 MiB, maximum 64 MiB). |
| `record_policy="packed-v1"` | Retain immutable packed exponent bytes with exact original identities, provenance and checked interchange. Existing conservative storage reservations remain in force. |

These `YieldSieveConfig` fields default to disabled policies. Its checkpoints use
schema 11 and bind the complete configuration; legacy checkpoints remain
readable. Restoring cannot silently change these policies. Packed records can
reduce payload storage but add decoding work; no performance or capacity gain
is established. The [benchmark runners](../benchmarks/qs/README.md) prepare
matched explicit experiments without changing production routing.

`qs.CombinedYieldSieveConfig` is the separate native-integer challenger. Selecting
this class enables the following combined policies; ordinary configurations
and the schema-11 policies above retain their defaults.

| Field / selected value | Behavior and reference alternative |
| --- | --- |
| `preparation_policy="active-v2"` | Count occupied columns, share checked row metadata and derive bounded top-ups from the filtered deficit. `fixed` retains ordinary cadence. |
| `recovery_policy="checked-v2"` | Exact impossible-domain rejection and bounded private atom/endpoint proofs. Full payload and provenance comparisons govern reuse; public checks and checkpoint admission remain independent. `reference` repeats arithmetic checks. |
| `split_policy="portfolio-v2"` | Bounded SQUFOF followed, only on failure, by two bounded Brent walks. Alternatives are `reference`, `rho-v2`, `ecm-v2` (one finite existing ECM fallback), and `rho-first-v3` (the two walks precede bounded SQUFOF fallback). All successful endpoints are proven and reconstructed. |
| `batch_policy="shared-v2"` | Batch GCD hints and a bounded certified-endpoint product. At most `batch_size=64` selected positions are looked ahead within a block. Alternatives are `known-v2` and `off`. Admissions remain in physical order. |
| `cycle_policy="short-v2"` | Direct-edge and minimal-path early exits, bounded two-frontier search and deterministic payload tie-breaking. `forest` retains the original path. |
| `duplicate_policy="congruence-v2"` | Exact unit-normalized congruence deduplication plus checked equal-parity pairs. Distinct square roots are preserved. `identity` disables both auxiliary indices. |
| `record_policy="indexed-v2"` | Base-index delta encoding, implicit exponent one, append-only slabs and integer provenance handles. `objects` retains ordinary records. |
| `filter_policy="light-v2"` | Fill-aware weight-two ordering and a reversible lighter-row view. Unsuccessful subset extraction retries the complete retained set. `fill-v2` uses all rows; `reference` keeps the original ordering. |
| `extraction_policy="joint-v2"` | Joint exponent products with bounded modular accumulation. `reference` uses independent prime powers. |

Checked compact records use exact byte/representation witnesses, including
the recomputation descriptor's own provenance, while public verification and
restore remain independent. Canonical duplicate keys retain the full exact
factor vector in base-index encoding. Expanded collectors retire raw proofs
only after successful packed publication. Serialized checkpoint formats are
unchanged. `recompute_prime_limit=0` retains every factor in the packed record;
the existing default remains 13 pending complete-cost comparison.

Power scoring may reserve a bounded cross-polynomial derivative-inverse cache
and an immutable recovery-mask prefix within the caller's memory allowance;
unfunded helpers keep exact uncached arithmetic. Scalar reject spans retain
unit work accounting, committed cursor/reporting state and the cooperative
polling bound. These implementation changes do not select optional collector
policies or alter automatic routing.

Dense bucket marking combines hit and score updates when the list backend
supports it. Periodic initializers own their complete active buffer slices.
Sparse queue construction uses exact additive offsets and grouped reservations
that retain each event's work and the original polling gap. Native root pairs
avoid sorting, but every polynomial root still receives independent arithmetic
verification. Benchmark warmup preserves the selected marking kernels.

Power scoring may also reserve up to 32 MiB for complete nonsingular odd-prime
root pairs across same-A polynomials. The exact B-shift preserves the original
root-class order; characteristic two, singular primes and unfunded caches use
ordinary lifting. Additional inverse construction remains charged. This helper
keeps polynomial-switch scratch funded, clears on target/A provenance changes
and rebuilds cold after restore. Its reservation is not an RSS estimate.

List-backed power scoring with bucket recovery may reserve another bounded
helper, up to 32 MiB after existing helpers and mandatory switch scratch.
For bases exceeding 60 entries, it marks recovery bits in 60-bit integer
words with block epochs and assembles exact wide masks only for selected
candidates. The word width keeps marking values in PyPy's native integer-list
range; it never truncates mathematical masks or changes admission bounds.
Construction, periodic stamp copies and internal materialization have explicit
work charges. Unfunded memory/work keeps ordinary lists; checkpoint restore
reconstructs derived scratch cold. No configuration field or schema is added.

The proof cache defaults to `checked_bytes=8 MiB` (maximum 64 MiB); the two
duplicate indices share `duplicate_bytes`. Batch sizes range from 2 through
256. These are finite experimental settings, not measured optima. Every helper
and its construction scratch is reserved within the existing owned allowance.
Compact `exponents` expose an immutable sequence with cheap `len` and iteration;
use `tuple(row.exponents)` when a concrete tuple is needed. External relation IDs
remain stable. Export through the checked checkpoint interface instead of
serializing implementation-specific dataclass fields directly.

Combined checkpoints use schema 12, bind all policy fields, and independently
rebuild records, graph ownership, endpoint proofs and solver progress. Batch
lookahead is disposable scratch: interruption preserves the first uncommitted
position and charged work; resume may repeat discarded lookahead. Combined
arithmetic is still materialized eagerly; the representation boundary does not
claim a lazy provenance solver. Performance scope and validation evidence are
recorded in [results](../benchmarks/RESULTS.md).

## Expanded collectors

`qs.ExtendedYieldSieveConfig` adds the separately versioned expanded
experimental arm. It preserves `CombinedYieldSieveConfig` and automatic
defaults. Its integer-handle incidence collector supports up to three proven
endpoints and explicitly typed opaque cofactor pairing. Opaque values are not
prime certificates; exact paired products and later proven endpoint conversion
retain original atoms. `verify_atomic` requires explicit `residual_arity=3`
and/or `opaque_pairs=True` to validate those extended contracts.

The expanded configuration includes deferred cycle arithmetic, direct indexed
recovery, exact piecewise reporting refinement, bounded reporting adaptation,
outcome-aware top-ups, component-local singleton/weight-two reuse, conditioned cofactor
splitting and memoization, dense-candidate base recovery, bounded cycle-basis
reduction, cheap-factor recomputation, sparse incidence, weight3–8 tree merges
and a bounded Four Russians pivot-panel solver. Original rows and exact lifting
remain available. A lighter view retries the full store when necessary.

`recovery_workers=1` enables one finite background hint task with a separately
reserved work allowance and a join on interruption; the default is serial.
Persistent cold factor records require an absolute existing `spill_directory`
and positive `spill_bytes`. The default uses no files. The spill cap is separate
from diagnostic output storage, survives relation resets and resume, and counts
orphaned partial writes. Files persist for external checkpoints; explicit
`SpillStore.cleanup()` invalidates those checkpoints. Missing, truncated or
corrupt referenced records are errors. Spilling reduces retained atom payload
ownership; construction/decode scratch and graph/solver ownership remain capped.

Expanded checkpoints use schema13. They serialize structural provenance
without forcing lazy modular arithmetic during emergency encoding. Restore
checks original arithmetic, incidence/conversion history, all policy fields,
reporting/outcome progress and pending solver/view continuation. Derived proof,
incremental and panel caches rebuild with charged work. Successful split proofs
are restored independently; volatile failed-split memo entries are discarded,
while earlier work and split-attempt charges remain consumed.
Pending solver checkpoints bind whether incremental contractions were used;
restore reproduces that choice. Panel solver `xors` counts logical reduction
actions, including one action for a whole panel; table construction and
primitive work are charged separately to the work budget.
Yield schemas 11–13 compose with the coupled polynomial policies by binding
their arithmetic marker and policy fields separately. Earlier yield checkpoints
retain standard arithmetic. Decoding into a historical runtime permits absent
standard policy fields and rejects unsupported coupled policies explicitly.

The expanded implementation remains experimental. Correctness and performance
evidence are tracked separately in [results](../benchmarks/RESULTS.md);
optional policies retain their declared resource, arithmetic and resume contracts.

Expanded `large_prime_graph` result statistics identify `incidence-v1` and
report retained atoms, residual rank, unowned atoms, closed relations, opaque
conversions and split calls. Ordinary collector forest statistics retain
their existing fields.

`qs.BoundedYieldSieveConfig` (`extended-v2`) is a separate experimental variant.
It stages unmatched opaque atoms in a FIFO capped by `opaque_cache_entries`
(default 256) and `opaque_cache_bytes` (default 4 MiB), outside retained incidence.
Exact repeats or proven factorizations promote the older atom before processing
the incoming one. Refusal preserves any completed promotion; closed rows keep
their original atoms. FIFO eviction can lose future pairing opportunities.

This variant charges graph vectors and lifting by their retained dimensions.
Packed slabs are charged once when allocated, with conservative object
overhead; polynomial contexts are shared only after exact field comparison.
Integer ownership follows digit storage, with separate bounds for atom/row
metadata and temporary packing. Allocations remain charged after an interrupted
packing operation, and resetting retained history releases its arena charges.
Its batch scratch uses exact norm/prime-product bounds, with a high-water
reservation established before each block. Work, memory and spill grants stay
binding. Schema 14 checkpoints save FIFO order and independently validate staged
arithmetic, caps and ownership on restore. Earlier variants/schemas remain
available. A completed checkpoint divisor survives factor-base reconstruction.
In-memory restore rebuilds packed factors and shared contexts. Existing cold
file references are preserved; older checkpoints do not encode cold references
for derived rows, so that restore path reserves resident row exponents instead.
Historical schema 14 payloads remain readable; internal ownership is recomputed
under the current implementation rather than treated as a serialized byte count.

`rho-first-v3` also retains a proper split found while classifying a two-prime
residual, so a three-prime retry can reuse it. Failure still proves no primality
claim; endpoint certification, work/time limits and split memo identities remain
binding. Shared lookahead prepays at most 64 cheap score comparisons per chunk;
an early full batch may leave up to 63 unused comparisons charged. Derived
lookahead never advances the committed position on refusal.

Panel elimination maintains a pivot-presence bitset, selects ordered keys
without sorting, and reuses its table until the selected immutable keys change.
It preserves elimination action order, provenance and checkpoint digests;
presence-mask work and scratch are reserved.

Expanded recovery rejects an oversized single-prime residual below the proven
residual-prime floor squared before hints or primality testing. The rho-first
policy avoids impossible three-prime recovery below the floor cubed while
retaining the final bounded rho walk when earlier splitting fails.
Successful two-prime memo entries remain valid under
a three-prime configuration; restore checks each entry's narrower arity and
exact endpoint proofs. Batch scratch bounds reuse their owned high-water mark.
Hensel lifting caches both invertible and singular derivative classes within
one polynomial; sparse queues retain their dense-root view until invalidation.
Unmatched opaque staging avoids an unused canonical duplicate key and copies
no FIFO map. Refused planning preserves every entry; an allocation failure
during publication poisons the collector and requires a prior checkpoint,
as retained relation publication does.

