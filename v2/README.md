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
limits. Explicit settings and saved configurations take precedence. Experimental
methods and adopted defaults retain their recorded evidence limits;
there is no general digit-based completion guarantee.

## Explicit useful-yield policies

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
is established. The [C9a runner](benchmarks/qs/c9a/) prepares matched explicit
experiments without changing production routing.

`qs.CombinedYieldSieveConfig` is the separate native-integer challenger. Selecting
this class enables the following combined policies; ordinary configurations
and the schema-11 policies above retain their defaults.

| Field / selected value | Behavior and reference alternative |
| --- | --- |
| `preparation_policy="active-v2"` | Count occupied columns, share checked row metadata and derive bounded top-ups from the filtered deficit. `fixed` retains ordinary cadence. |
| `recovery_policy="checked-v2"` | Exact impossible-domain rejection and bounded private atom/endpoint proofs. Full payload and provenance comparisons govern reuse; public checks and checkpoint admission remain independent. `reference` repeats arithmetic checks. |
| `split_policy="portfolio-v2"` | Bounded SQUFOF followed, only on failure, by two bounded Brent walks. Alternatives are `reference`, `rho-v2`, and `ecm-v2` (one finite existing ECM fallback). All successful endpoints are proven and reconstructed. |
| `batch_policy="shared-v2"` | Batch GCD hints and a bounded certified-endpoint product. At most `batch_size=64` selected positions are looked ahead within a block. Alternatives are `known-v2` and `off`. Admissions remain in physical order. |
| `cycle_policy="short-v2"` | Direct-edge and minimal-path early exits, bounded two-frontier search and deterministic payload tie-breaking. `forest` retains the original path. |
| `duplicate_policy="congruence-v2"` | Exact unit-normalized congruence deduplication plus checked equal-parity pairs. Distinct square roots are preserved. `identity` disables both auxiliary indices. |
| `record_policy="indexed-v2"` | Base-index delta encoding, implicit exponent one, append-only slabs and integer provenance handles. `objects` retains ordinary records. |
| `filter_policy="light-v2"` | Fill-aware weight-two ordering and a reversible lighter-row view. Unsuccessful subset extraction retries the complete retained set. `fill-v2` uses all rows; `reference` keeps the original ordering. |
| `extraction_policy="joint-v2"` | Joint exponent products with bounded modular accumulation. `reference` uses independent prime powers. |

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
claim a lazy provenance solver. Correctness execution for this expanded bundle
is queued for a quiet machine window; no performance acceptance is implied.

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
