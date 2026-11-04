# Quadratic-sieve and relation engines

QS/MPQS, self-initializing quadratic sieve (SIQS), experimental Smooth Subsum
Search (SSS/SSSf), and shared relation/filtering/extraction machinery for
PyPy implementing Python 3.11. The Python-integer backend needs no third-party
packages.

## Setup and run

Run from the repository root with a PyPy installation implementing Python 3.11:

```sh
make -C v2 runtime
pypy3 -m v2.factor 626100403 --method siqs --bounded --seed 7 \
  --work-limit 2000000 --seconds 30 --cpu-seconds 30
pypy3 -m v2.factor --help
```

| Engine | Role |
| --- | --- |
| QS | Fixed-polynomial reference collector |
| MPQS | Independent polynomials, including external-square coefficients; [Silverman, *The Multiple Polynomial Quadratic Sieve* (1987)](https://doi.org/10.1090/S0025-5718-1987-0866119-8). |
| SIQS | CRT/Gray families, reused roots and shared verified relations; [Alford–Pomerance, *Implementing the Self Initializing Quadratic Sieve on a Distributed Network*](https://math.dartmouth.edu/~carlp/implementing.pdf). |
| SSS | CRT collision candidates and exact smoothness recovery; [Hittmeir, *Smooth Subsum Search* (2023)](https://arxiv.org/abs/2301.10529v2). |
| SSSf | SSS with an optional lossy candidate filter |

The CLI accepts `--method qs|mpqs|siqs|sss|sssf`; `--siqs` enables portfolio
fallback. See the [CLI guide](../docs/CLI.md) for resource and checkpoint flags.

For a bounded library job:

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

Jobs return a split and remaining cofactor. Recursive terminal classification
belongs to [`factorize_bounded`](../docs/API.md). Omitting the SIQS configuration
uses automatic selection; explicit `SIQSConfig()` keeps reference fields.
The [configuration guide](CONFIGURATION.md) describes defaults, optional
policies and supported combinations; [collector policies](COLLECTORS.md)
document the detailed intake, storage and checkpoint contracts.

## Techniques

The engines combine exact arithmetic with bounded collection and reuse.
The papers linked here describe algorithmic foundations; buffer layouts,
resource accounting and policy combinations are specific to this implementation.

| Area | Techniques available |
| --- | --- |
| Polynomial setup | Exact factor-base roots and eligible-prime counts; parity-aware polynomial forms; coupled multiplier scoring and optional norm previews; streamed A assignments with reachable product windows and paired completion; CRT/Gray families and adjacent sign schedules. |
| Root reuse | Bounded modular Gray-update tables, derivative-inverse reuse, cached prime-power plans and complete same-A higher-root transfers, with ordinary exact lifting for exceptional or unfunded cases. |
| Sieving | Conservative integer/fixed-point scores, prime-power marking, periodic presieves, guaranteed two-adic initialization, density/unrolled kernels, fused score/hit updates and candidate scanning. Flat sparse queues use additive offsets; native-width hit words and epochs defer wide-mask construction to selected candidates. |
| Recovery and admission | Separate reporting/product/endpoint bounds; exact impossible-domain rejection; bounded splitter portfolios, batch-GCD and known-endpoint hints, split memoization and candidate-only smooth-part recovery using product/remainder trees ([Bernstein, *How to Find Smooth Parts of Integers*, 2004](https://cr.yp.to/factorization/smoothparts-20040510.pdf)). Admissions reconstruct residuals and verify their declared types. |
| Relation graphs | Single/double-large-prime matching ([Lenstra–Manasse, *Factoring with Two Large Primes*, 1994](https://doi.org/10.1090/S0025-5718-1994-1250773-9)), optional three-endpoint and opaque-residual policies, incidence elimination, shorter-cycle searches, congruence deduplication and checked equal-parity dependencies. Original atoms, repeated primes and square corrections remain recoverable. |
| Records and ownership | Packed/base-index exponents, append-only slabs, integer provenance, lazy combined rows, exact representation witnesses and checked proof reuse. Optional bounded opaque staging/conversion, cheap-factor recomputation and disk-backed cold records control retained ownership. |
| Filtering and extraction | Useful-row/deficit targets and collection-aware solve cadence; duplicate, singleton and weight-two filtering; reversible sparse/component views, fill-aware merges and optional panel kernels; joint exponent products and exact congruence/GCD extraction. |
| Execution | Reused buffers, progress-only internal collection, grouped work reservations and bounded cooperative polling. SSS uses CRT collisions and product/remainder trees; experimental serial/thread/process family workers verify returned relations centrally. |

Some techniques are explicit policy choices. Extra reporting thresholds and
SSSf filters can deliberately lose candidates; growing-DLP and expanded
collectors have narrower worker/SSS compatibility. Public `collect()` returns
owned snapshots while internal jobs can retain the live store.

Every admitted relation retains checked sign, exponents and provenance. Parity
alone is insufficient: dependencies lift to original relations, both GCD signs
remain available, and every returned split must satisfy `1 < g < n` and divide
n exactly. Trivial dependencies permit finite further collection.

## Resume and resource bounds

`job.run(max_blocks=...)` pauses at block boundaries. Save with
`job.checkpoint()` and restore with `SIQSJob.from_checkpoint(...)`, supplying
the total allowance. Restore verifies arithmetic, provenance and saved policy
identities; compatible extensions use `allow_extension=True` and retain consumed
work and traversal progress. `QSJob` supports in-memory continuation; SSS has
its own checked serializer.

Work, time and owned storage remain finite. Derived root tables, presieves,
queues and hit helpers rebuild with charges and exact fallbacks. A refusal
before publication is resumable. `RelationPublicationError` invalidates the
live store; restore its previous completed checkpoint. Explicit cold-record
storage does not automatically spill the growing relation graph.

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

[QS benchmarks](../benchmarks/qs/README.md) describe the runners.
[Results](../benchmarks/RESULTS.md) retain measurement scope;
[STATUS](../STATUS.md) and [ROADMAP](../ROADMAP.md) track current work.
