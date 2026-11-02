# Tests

Acceptance and regression tests for exact arithmetic, validated results,
finite budgets, checkpoint compatibility and experimental adapters.
Use PyPy implementing Python 3.11 and run from the repository root:

```sh
make -C v2 test
make -C v2 lint
```

The native track needs no optional GMP package. Backend-specific tests use
the configured PyPy environment and report unavailable optional dependencies.
The lint environment is described in [getting started](../README.md#development).

Full suites share the machine-wide performance window with timed experiments.
Acquire `/private/tmp/factor-performance.lock` exclusively and avoid competing
heavy jobs. Some acceptance tests inspect processes to enforce this boundary;
process-inventory access must be available. The suite is correctness evidence,
not a factoring benchmark or speed claim.

Fixtures, immutable baselines and numerical controls belong in versioned
[benchmark inputs](../benchmarks/inputs/README.md). Generated captures and
scratch stay in ignored `results/` trees. Before publishing, check a tree
containing only committed inputs; local audit material must not be required.
See [contributor guidance](../../AGENTS.md) for change-specific checks.
