# Contributor guidance

Explicit user instructions take precedence. Work on the active implementation
in `v2/`; preserve `v1/` as the original baseline.

## Code and arithmetic

- Target PyPy implementing Python 3.11. Do not silently fall back to CPython.
- Use PEP 8, descriptive names, relative imports and concise docstrings.
- Keep factoring arithmetic exact with integers, integer roots and inversion.
- Validate proper divisors and reconstruct every result, including unresolved
  cofactors. Keep probable-prime and proven-prime labels distinct.
- Preserve finite work/time/storage allowances and documented resume behavior.
- Keep library calls quiet unless explicit verbosity is requested.

## Checks and measurements

Use `make -C v2 test` and `make -C v2 lint` for relevant changes. Benchmark
behavior/performance changes with matched inputs, seeds, budgets and output
validation. Use at least three seconds of validated PyPy warmup and nine
samples; extend unstable measurements. Separate cold startup from warmed
execution, and instrumented profiles from performance evidence.

Document API changes in `v2/README.md`. Update roadmap completion only after
its acceptance and experiment gates pass. Summarize accepted behavior in
`CHANGELOG.md` and useful measurements in `v2/benchmarks/README.md`.

## GitHub contents

Commit code, tests, documentation, independent corpora and immutable baselines
required by test/benchmark loaders. Required inputs live in versioned
`v2/benchmarks/inputs/`; generated evidence belongs in folder-ignored
`results/` trees. The `v2/audit/` directory is local and Git-ignored; maintain
the public acceptance plan in `v2/ROADMAP.md`.

Keep generated captures, stdout transcripts, profiles, verification dumps,
scratch output and the detailed `v2/LOG.md` journal local. Do not add those
files by force. Preserve useful raw evidence locally or in durable external
storage; publish a concise summary with commands and limitations.

Before publishing, verify a checkout containing only committed files can run
the tests and import benchmark runners. Never include unrelated work in a
cleanup commit.
