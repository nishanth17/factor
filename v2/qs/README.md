# Quadratic-sieve and relation engines

QS/MPQS, self-initializing quadratic sieve (SIQS), experimental Smooth Subsum
Search (SSS/SSSf), and shared relation/filtering/extraction machinery for
PyPy implementing Python 3.11.

## Choose a route

| Engine | Role |
| --- | --- |
| QS | Fixed-polynomial reference collector |
| MPQS | Multiple independently chosen polynomials, including external-square coefficients |
| SIQS | CRT/Gray families, reused roots and shared verified relations |
| SSS | CRT collision candidates and exact smoothness recovery |
| SSSf | SSS with an optional lossy candidate filter |

The CLI exposes `--method qs|mpqs|siqs|sss|sssf` and `--siqs` portfolio fallback.
Native omitted configurations also enable automatic SIQS in the scoped bands
shown in [STATUS](../STATUS.md#current-automatic-defaults). Double-large-prime
collection is implemented; worker export and SSS support retain their narrower
contracts. GNFS remains deferred roadmap work.

## Run a job

```python
from v2.execution.budget import Budget
from v2.qs import SIQSConfig, SIQSJob

n = 4001 * 5003
job = SIQSJob(
    n, seed=7,
    config=SIQSConfig(base_bound=200, half_width=256),
    budget=Budget(work_limit=200_000_000, seconds=30, cpu_seconds=30),
)
result = job.run()
assert (result.divisor or 1) * result.cofactor == n
print(result.reason)
```

A job returns a split/cofactor result; recursive terminal classification belongs
to [`factorize_bounded`](../docs/API.md). `QSJob` works on an explicit polynomial
and signed half-open window. `SIQSConfig(mode="qs")` or `mode="mpqs"` selects reference
controls through the SIQS job interface.

`default_siqs_config(n, backend="python-int", memory_bytes=None)` exposes fresh
automatic selection. An explicit `SIQSConfig()` uses reference fields, while
omitting it allows size selection. The [configuration guide](CONFIGURATION.md)
explains collectors, kernels, quotas and exact fallback.

## C10a polynomial, collection and scheduling techniques

Fresh automatic SIQS selects the coupled C10a-01 policies below. Explicit
`SIQSConfig()` constructors keep reference choices, and saved configurations
keep their original identities. Automatic portfolio eligibility remains native
40 and 50–60 digits; selecting coupled policies at other sizes does not supply
a calibrated larger-size preset. See [configuration semantics](CONFIGURATION.md)
and [current defaults](../STATUS.md#current-automatic-defaults).

| Technique | Implementation and selection |
| --- | --- |
| Parity-aware polynomials | `polynomial_policy="parity2"` uses `A=2a` and odd centered `B` when the target `kN` is 1 modulo 4. The exact identity `(A*x+B)^2-kN=A*F(x)` retains the supported-A contribution and powers of two in recovery; unsupported discriminants use the standard family. |
| Coupled multiplier ranking | `multiplier_policy="parity2_wide"`, with `multiplier=0`, scores a finite odd/even squarefree population using representation-aware root weights and size costs. Exceptional proper GCDs can finish setup early. A positive explicit multiplier fixes the choice. |
| Reachable A products | `assignment_policy="windowed_pair"` uses bounded product windows and paired completion to search reachable factor-base products. It changes the A population; it is distinct from the explicit `quality` ordering and `paired_flyer` alternatives. |
| Bounded Gray-root tables | `root_update_policy="batch"` prepares modular coordinate deltas only for the consumed family prefix, reusing existing A inverses. `cached` is an explicit alternative. Physical trial order, exceptional roots and independent root verification are preserved. |
| Guaranteed two-adic scores | `parity_two` fuses certified universal two-power contributions into supported list-backed score initialization, with exact ordinary scoring elsewhere. It composes with existing marking kernels. |
| Adjacent sign faces | C10a-07's `sign_policy="swap_high"` changes the adjacent Gray-coordinate ordering while retaining one-sign root transitions. Reference configurations keep `canonical`; incompatible saved schedules are rejected. |
| Sparse prime/power queues | Bounded per-polynomial queues reuse exact hit schedules across aligned blocks. Calibrated automatic collectors select a 32 MiB ceiling; unsupported paths, singular sparse roots or capacity refusal fall back to ordinary marking. Rebuilds and replay remain charged. |
| Counted bases and collection pacing | C10a-04 resolves exact eligible-prime counts under an endpoint ceiling for larger calibrated bands. C10a-05's `solve_policy="collection"` separates polynomial-window exhaustion from genuine collection exhaustion, preserving immediate square checks and charged terminal extraction. Their size/resource selection is documented in STATUS. |

`parity2_norm` multiplier previews, alternative assignment policies and marking
kernels remain explicit choices. C10a-02/-03 and broader scheduling questions
remain scoped roadmap work; this list does not mark every C10a proposal
implemented or adopted. The [C10a-01](../benchmarks/RESULTS.md#c10a-01-coupled-polynomial-defaults)
and [C10a-07](../benchmarks/RESULTS.md#c10a-07-face-and-queue) records retain
confirmation limits and unmeasured extrapolations.

## Verified relation pipeline

1. Build the factor base and exact modular roots, checking setup divisors.
2. Construct polynomials with `B² ≡ h*n (mod A)` and check the norm identity.
3. Recover every candidate's sign and exponents, then verify residual admission.
4. Match single-large-prime partials or combine checked DLP graph cycles while
   retaining original atoms and square corrections.
5. Filter duplicates/singletons/weight-two rows and solve parity dependencies.
6. Lift to original relations, reconstruct a congruence of squares and validate
   proper GCD splits. Trivial dependencies permit finite further collection.

Parity alone does not verify a relation. Original exponents, atoms, row
identities and square corrections remain checkable through filtering and
extraction. Every split must satisfy `1 < g < n` and divide n exactly.

## C9a useful-yield and storage techniques

C9a is integrated and its delivered stream is closed. Its policy classes
remain opt-in through `SIQSConfig(collector=...)`: `YieldSieveConfig` exposes
independent mechanisms, `CombinedYieldSieveConfig` combines the checked
recovery/storage policies, and `ExtendedYieldSieveConfig` and
`BoundedYieldSieveConfig` add explicit expanded intake and ownership controls.
The [policy reference](../README.md#explicit-useful-yield-policies) lists fields,
alternatives and checkpoint schemas. Selecting a class does not imply that
every implemented mechanism is the measured nominee or an automatic default.

| Technique | Exact behavior |
| --- | --- |
| Useful-row preparation | Targets can use occupied columns and the filtered deficit, with bounded top-ups after unsuccessful extraction. Preparation, filtering and elimination remain separate costs; immediate square checks and final/storage-cap attempts are retained. |
| Separate reporting and admission | Candidate reporting, recovered-product admission, ordered asymmetric endpoints and SLP bounds have separate controls. Exact impossible-domain checks reject hopeless residuals before splitter work; successful splits reconstruct the residual and certify endpoints under the selected contract. |
| Bounded cofactor recovery | Selectable finite splitter portfolios, known-endpoint products, batch-GCD hints, split memoization and candidate-only smooth-part recovery reduce repeated work. Hints undergo exact checks, and batched admissions retain physical candidate order. |
| Shorter cycles and duplicates | Direct-edge/minimal-path exits and bounded searches cover retained graph components. Full congruence keys and equal-parity pair checks preserve original identities and distinct square roots; equal parity alone does not justify deleting a useful dependency. |
| Compact exact records | Packed/base-index exponent encodings, implicit exponent one, append-only slabs, integer provenance and lazy combined rows reduce repeated object construction. Representation witnesses permit checked private proof reuse; public and restore validation remain independent. Obsolete raw proofs retire only after packed publication. |
| Expanded residual ownership | Explicit policies support three proven endpoints and separately typed opaque residuals, incidence elimination, bounded opaque staging/conversion, cheap-factor recomputation and disk-backed cold records. Opaque residuals are not relabeled as proven primes; byte limits include simultaneous construction/replay scratch. |
| Filtering and extraction | Reversible lighter-row/component views, incremental singleton work and optional panel kernels preserve lifting to original rows. Joint exponent-product extraction keeps exact modular accumulation, both GCD signs and nonunit handling. |

Shared collector optimizations also apply to supported existing configurations:

- **Avoid repeated snapshots and scans.** Internal QS/SIQS collection consumes
  private progress while retaining the live store; public `collect()` still
  returns complete owned snapshots. Scalar reject spans and grouped queue
  reservations retain exact work, the cooperative polling bound and the first
  uncommitted position on refusal.
- **Fuse and reuse marking work.** List-backed first-level score/hit writes
  share a traversal, periodic initializers replace redundant clears, and
  sparse queues use flat unsigned-word triples with additive offsets. Small
  root pairs avoid sorting while independent polynomial checks remain.
- **Reuse exact roots.** Funded derivative inverses and recovery-mask prefixes
  avoid recomputation. Complete nonsingular odd-prime root pairs transfer
  across same-A polynomials by `-(B_new-B_old)/A` modulo each prime power.
  Root-class order is preserved; pruned subsets, characteristic two, singular
  cases and unfunded caches retain ordinary lifting.
- **Keep marking metadata native-width.** Supported list/power/bucket paths
  with more than 60 base entries use bounded 60-bit word lists and block
  epochs, assembling full exact masks only for selected candidate recovery.
  The width avoids wide-integer allocation at every hit; it does not truncate
  masks or impose an admission limit. Construction, stamp copies and internal
  materialization are charged. Root-transfer and hit-word helpers each have
  a 32 MiB ceiling, preserve switch scratch and use exact unfunded fallbacks.

These are implemented techniques, not additive performance guarantees.
[Results](../benchmarks/RESULTS.md#c9a-preparation) distinguish same-work
optimizations from changed intake/yield and record the exposed-input limits.
C9b owns further headroom and variability analysis. The historical 15% canary
heuristic is descriptive in new analysis; calibrated uncertainty and fresh
independent inputs govern broader performance conclusions.

## Resume and storage

`job.run(max_blocks=...)` pauses at SIQS block boundaries; `job.checkpoint()`
produces portable state. Restore with `SIQSJob.from_checkpoint(checkpoint,
budget=total_allowance, config=...)`. Exact relation, provenance and policy
verification plus reconstruction consume the cumulative allowance. Disposable
plans, presieve buffers, sparse queues and root/hit helpers are rebuilt with
charges. Optional helpers leave mandatory switch coexistence funded and keep
exact uncached paths when their memory/work reservation cannot be met.

Direct compatible quota extensions require `allow_extension=True`; they
preserve the modulus, family/collector identity and consumed prefix. An
exhausted local schedule does not restart merely because a timer increases.
`QSJob` supports in-memory continuation; SSS has its own checked serializer.

Fixed-store collectors have explicit row/atom/partial limits. Growing DLP
uses allocator-grown containers and optional quotas, with finite work/time
and per-residual split bounds. The outer portfolio resolves uncapped child
memory against its finite remaining allowance. Unsupported SSS/worker adapters
reject the growing collector rather than weakening ownership checks.

A quota refusal before publication remains resumable. Actual allocator failure
during publication raises `RelationPublicationError` and invalidates that live
object; restore the previous completed checkpoint. An OS kill requires an
already saved checkpoint. Expanded policies can explicitly select bounded
disk-backed cold factor records; this does not enable automatic spilling of
the growing relation graph.

## SSS and workers

`SSSJob` and `SSSConfig` live in `sss.py`. Positive SSSf `filter_bound`
deliberately loses candidates; zero keeps two-stage smoothness without that
cutoff. Collision trees inspect forced-divisor quotients, while admission
recovers exponents from the original norm.

In-memory SSS continuation extends its original budget; serialized restore
uses an unused total allowance and preserves backend/configuration and saved
solver encoding. Assignment regeneration and solver replay are charged.

`parallel.py` supplies bounded serial/thread/spawned-process family workers.
They centrally verify returned relations, use stable assignments and account
for aggregate resources and cancellation. Workers remain experimental; use
existing measured evidence rather than assuming extra workers improve factoring.

## Source map

| Modules | Responsibility |
| --- | --- |
| [__init__.py](__init__.py), [siqs.py](siqs.py), [defaults.py](defaults.py) | Public exports, QS/MPQS/SIQS orchestration and fresh automatic selection |
| [factor_base.py](factor_base.py), [polynomial.py](polynomial.py), [score.py](score.py) | Exact factor-base roots, norm identities and integer score bounds |
| [multiplier.py](multiplier.py), [norm_multiplier.py](norm_multiplier.py) | Reference/coupled multiplier scoring and explicit bounded norm previews |
| [count_policy.py](count_policy.py), [capacity.py](capacity.py) | Charged eligible-base counts, capacity diagnostics and monotone resume allowances |
| [families.py](families.py), [assignment_stream.py](assignment_stream.py), [product_window.py](product_window.py) | CRT/Gray families, bounded assignment streams and reachable product windows |
| [root_updates.py](root_updates.py), [coordinate_family.py](coordinate_family.py), [sign_schedule.py](sign_schedule.py) | Bounded modular Gray deltas, coordinate views and adjacent sign schedules |
| [external_square.py](external_square.py), [reference_collector.py](reference_collector.py) | External-square MPQS coefficients and exhaustive division control |
| [sieve_collector.py](sieve_collector.py), [power_sieve.py](power_sieve.py), [collector_kernels.py](collector_kernels.py) | Collector policies, exact Hensel lifting, conservative scores and marking/scan kernels |
| [sparse_queue.py](sparse_queue.py), [power_reuse.py](power_reuse.py), [hit_words.py](hit_words.py) | Per-polynomial flat queues, complete same-A higher-root transfers and native-width recovery metadata |
| [relations.py](relations.py), [storage.py](storage.py) | Exact relation/provenance verification and transactional publication contracts |
| [large_primes.py](large_primes.py), [large_prime_collector.py](large_prime_collector.py), [large_prime_checkpoint.py](large_prime_checkpoint.py) | Residual certification, DLP graph admission and checked DLP restore |
| [yield_support.py](yield_support.py), [yield_collector.py](yield_collector.py), [yield_feedback.py](yield_feedback.py) | Checked proof/congruence reuse, useful intake and bounded reporting/target feedback |
| [extended_policy.py](extended_policy.py), [extended_collector.py](extended_collector.py), [residual_relations.py](residual_relations.py) | Expanded/bounded policy classes, generalized transactional admission and residual type contracts |
| [cofactor_policy.py](cofactor_policy.py), [recovery_batch.py](recovery_batch.py) | Finite splitter strategies, batch hints, memoization and candidate smooth-part recovery |
| [opaque_staging.py](opaque_staging.py), [opaque_conversion.py](opaque_conversion.py), [yield_memory.py](yield_memory.py) | Bounded opaque FIFO/conversion and shape-based ownership/scratch reservations |
| [residual_kernel.py](residual_kernel.py), [short_cycles.py](short_cycles.py), [equal_parity.py](equal_parity.py) | Integer-handle residual incidence, bounded shorter paths and checked two-row dependencies |
| [packed_relations.py](packed_relations.py), [indexed_relations.py](indexed_relations.py), [factor_witness.py](factor_witness.py) | Packed/indexed factors, append-only provenance slabs and exact representation witnesses |
| [lazy_relations.py](lazy_relations.py), [cold_factors.py](cold_factors.py) | Deferred combined-row arithmetic, cheap-factor recomputation and bounded persistent cold records |
| [linear_algebra.py](linear_algebra.py), [structured_matrix.py](structured_matrix.py), [incremental_structure.py](incremental_structure.py) | Provenance-preserving GF(2) filtering, sparse/panel kernels and incremental singleton structure |
| [extraction.py](extraction.py), [pipeline.py](pipeline.py) | Verified preparation, dependency/extraction orchestration and validated congruence splits |
| [checkpoint.py](checkpoint.py), [extended_checkpoint.py](extended_checkpoint.py) | Versioned checked SIQS stores and expanded-policy restore |
| [sss.py](sss.py), [smooth_batch.py](smooth_batch.py), [sss_checkpoint.py](sss_checkpoint.py) | SSS/SSSf collision search, exact smooth-part trees and checked SSS continuation |
| [parallel.py](parallel.py) | Experimental bounded family workers and central result verification |

[QS benchmarks](../benchmarks/qs/README.md) map studies to runners.
[Results](../benchmarks/RESULTS.md) retain quantitative scope and limitations;
[STATUS](../STATUS.md) and [ROADMAP](../ROADMAP.md) govern current work.
