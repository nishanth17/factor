# Factor

Integer factorization in Python, targeting **PyPy implementing Python 3.11**.
[Factor v2](v2/README.md) is the active implementation;
[v1](v1/) preserves the original Python 2 code for baseline comparisons.

v2 combines exact preprocessing, Brent rho, Pollard p−1 and Montgomery ECM
in a seeded, bounded, resumable portfolio. It also includes experimental
QS/MPQS, self-initializing quadratic sieve (SIQS), Smooth Subsum Search (SSS)
and filtered SSSf engines. Results preserve signs, multiplicities and unfinished
cofactors, with separate probable-prime and proven-prime labels.

## Run

Install a PyPy interpreter implementing Python 3.11, then run from this directory:

```sh
make -C v2 runtime
make -C v2 run N=626100403 SEED=7
pypy3 -m v2.factor 626100403 --bounded --seed 7 \
  --work-limit 2000000 --seconds 30 --cpu-seconds 30
pypy3 -m v2.factor 10002200057 --method siqs --seed 7
pypy3 -m v2.factor 10002200057 --method qs --qs-base-bound 400 --seed 7
pypy3 -m v2.factor 10002200057 --method mpqs --qs-base-bound 400 --seed 7
pypy3 -m v2.factor 10002200057 --siqs --ecm-curves 2 --seed 7
make -C v2 test
```

Omit the number for an interactive prompt. `--method qs|mpqs|siqs` selects
that engine after exact preprocessing; `--siqs` adds SIQS after rho/p−1/ECM.
These selections use shared limits and recursive factoring of the returned
children.

The factoring library uses standard-library integers and needs no third-party
runtime packages. Development lint tools have a separate local environment.
See the [v2 guide](v2/README.md) for setup, library calls, checkpoints and
method selection. The supported runtime is PyPy Python 3.11.

## What changed from v1

- Corrected arithmetic, sieve boundaries, retry and saturated-batch recovery;
  every accepted split is a proper divisor and every result reconstructs.
- Replaced ambiguous failure values with explicit complete or partial results.
- Added shared work/time/storage allowances, streamed schedules and checkpoints.
- Added verified relation collection, GF(2) solving and exact factor extraction
  for SIQS and related experimental engines.
- Reduced repeated preprocessing, root, relation and matrix work, with
  reproducible comparisons and independent validation.

A historical PyPy comparison measured **6.6% less elapsed time** on the original
25-factorization batch and **19.4% less** on a separate 280-factorization control.
These small-workload results compare repaired v2 with emulated v1 on the same
interpreter; they are not current large-number performance guarantees.
See [benchmark scope and results](v2/README.md#benchmarks).

## Repository

| Path | Contents |
| --- | --- |
| [v2/](v2/README.md) | Current methods, API, CLI and tests |
| [v1/](v1/) | Original baseline, preserved unchanged |
| [Benchmarks](v2/benchmarks/README.md) | Commands, measurements and limitations |
| [Roadmap](v2/ROADMAP.md) | Acceptance gates and planned work |
| [Changelog](CHANGELOG.md) | Accepted changes |

Benchmark `inputs/` retain reproducible corpora, immutable baselines and
required provenance. Generated captures, profiles, checkpoints and scratch work
belong in ignored `results/` folders. The entire `v2/audit/` tree and historical
raw evidence stay local. See
[contributor guidance](AGENTS.md).

The next major direction is a bounded GNFS pipeline, alongside arithmetic and
portfolio experiments. Large balanced inputs remain demanding; there is no
universal digit cutoff or guaranteed time to factor.
