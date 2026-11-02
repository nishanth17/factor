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

## Resume and storage

`job.run(max_blocks=...)` pauses at SIQS block boundaries; `job.checkpoint()`
produces portable state. Restore with `SIQSJob.from_checkpoint(checkpoint,
budget=total_allowance, config=...)`. Exact relation, provenance and policy
verification plus reconstruction consume the cumulative allowance. Disposable
plans, presieve buffers and sparse queues are rebuilt with charges.

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
already saved checkpoint. No disk-spill implementation is implied.

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
| [factor_base.py](factor_base.py), [polynomial.py](polynomial.py), [multiplier.py](multiplier.py) | Exact setup, roots, coefficients and multiplier scoring |
| [families.py](families.py), [assignment_stream.py](assignment_stream.py), [external_square.py](external_square.py) | CRT/Gray reuse and streamed assignments |
| [sieve_collector.py](sieve_collector.py), [power_sieve.py](power_sieve.py), [collector_kernels.py](collector_kernels.py) | Scoring, exact recovery and marking kernels |
| [large_primes.py](large_primes.py), [large_prime_collector.py](large_prime_collector.py), [relations.py](relations.py), [storage.py](storage.py) | Residual certification, graph cycles, provenance and ownership |
| [linear_algebra.py](linear_algebra.py), [extraction.py](extraction.py), [pipeline.py](pipeline.py) | Filtering, dependencies, congruences and splits |
| [siqs.py](siqs.py), [checkpoint.py](checkpoint.py), [count_policy.py](count_policy.py) | Jobs, concrete config and checked restore |
| [sparse_queue.py](sparse_queue.py), [coordinate_family.py](coordinate_family.py), [sign_schedule.py](sign_schedule.py) | Exact sparse scheduling and adjacent sign policies |
| [sss.py](sss.py), [smooth_batch.py](smooth_batch.py), [parallel.py](parallel.py) | Experimental candidates, smoothness arithmetic and workers |

[QS benchmarks](../benchmarks/qs/README.md) map studies to runners.
[Results](../benchmarks/RESULTS.md) retain quantitative scope and limitations;
[STATUS](../STATUS.md) and [ROADMAP](../ROADMAP.md) govern current work.
