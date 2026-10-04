# Factor

Integer factorization in Python, targeting **PyPy with Python 3.11**.

The active implementation is [v2](v2/). It combines exact preprocessing,
Brent rho, Pollard p−1 and ECM in a bounded, resumable portfolio. Results
preserve signs, multiplicities and unresolved cofactors, and distinguish
probable primes from proven primes.

An experimental QS/MPQS package supplies verified relation collection,
filtering, GF(2) dependencies and factor extraction. SIQS self-initialization
and production dispatch are the next milestone.

## Quick start

Install PyPy with Python 3.11; on macOS, `brew install pypy3`. From the
repository root:

```sh
make -C v2 run N=626100403 SEED=7
make -C v2 run-bounded N=626100403 SEED=7
make -C v2 test
```

See [usage and API contracts](v2/README.md) for library calls, resource limits,
checkpoints and development setup.

## Repository

| Path | Purpose |
| --- | --- |
| [v2/](v2/) | Active implementation and tests |
| [v1/](v1/) | Preserved original Python 2 baseline |
| [Roadmap](v2/audit/TODOS.md) | Accepted phases and remaining work |
| [Benchmarks](v2/benchmarks/README.md) | Runners and retained inputs |
| [Research](v2/audit/README.md) | Algorithm notes and reference prototypes |
| [CHANGELOG.md](CHANGELOG.md) | Public development summary |

Generated run captures, profiler output, transcripts and detailed local
journals are excluded from Git. Benchmark code, independent corpora and
required source baselines remain versioned. See
[contributor guidance](AGENTS.md).

CPython is unsupported. The original baseline is retained for comparisons,
not as a supported runtime or production implementation.
