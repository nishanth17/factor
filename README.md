# Factor

Integer factorization in Python, running on **PyPy implementing Python 3.11**.
The active implementation is [v2](v2/README.md); [v1 baseline](v2/history/V1_BASELINE.md) preserves
the original Python 2 baseline.

Factor combines exact preprocessing, Brent rho, Pollard p−1 and Montgomery
ECM with quadratic-sieve and experimental relation engines. Results preserve
factor multiplicities and unfinished cofactors, and distinguish probable
primes from proven primes.

## Quick start

With PyPy 3.11 installed, run from the repository root:

```sh
make -C v2 runtime
pypy3 -m v2.factor 626100403 --bounded --seed 7 \
  --work-limit 2000000 --seconds 30 --cpu-seconds 30
```

Omit the number for an interactive prompt. The Python-integer backend uses
only the standard library; GMP is optional.

## Documentation

| Guide | Use it for |
| --- | --- |
| [Getting started](v2/README.md) | Setup, methods and development commands |
| [Command line](v2/docs/CLI.md) | Engine selection and resource options |
| [Python API](v2/docs/API.md) | Results, bounded calls and library imports |
| [Checkpoints](v2/docs/CHECKPOINTS.md) | Saving progress and cumulative resume limits |
| [Current status](v2/STATUS.md) | Default routing, evidence limits and next work |
| [Benchmarks](v2/benchmarks/README.md) | Runners, measured results and reproduction |
| [Changelog](CHANGELOG.md) | Accepted changes |
| [Contributor guidance](AGENTS.md) | Project rules and experiment design |

Start with `v2/` for current usage. Its package guides describe each engine
and shared execution layer.

## License

Active v2 project code and documentation use the [MIT License](v2/LICENSE).
External dependencies and linked papers keep their own terms. The preserved v1
baseline is outside this new grant; its original provenance is unchanged.
