# Contributor guidance

Explicit user instructions take precedence. Work on the active implementation
in `v2/`; preserve `v1/` as the original baseline.

## Code and arithmetic

- Target PyPy implementing Python 3.11. Do not silently fall back to CPython.
- Keep factoring arithmetic exact with integers, integer roots and inversion.
- Validate proper divisors and reconstruct every result, including unresolved
  cofactors. Keep probable-prime and proven-prime labels distinct.
- Preserve finite work/time/storage allowances and documented resume behavior.
- Keep library calls quiet unless explicit verbosity is requested.

## Readability and style

- Use PEP 8, descriptive names, relative imports and concise docstrings.
- Separate logical stages with blank lines: setup, validation, computation,
  recovery and result construction. Within functions, space distinct branches
  and loops when they represent separate steps; avoid blank lines after every
  statement or inside a tightly related operation.
- Keep expressions readable through sensible wrapping and clear intermediate
  names. Prefer small, focused changes over restructuring working algorithms.
- Comment nontrivial reasoning: arithmetic invariants, budget reservations,
  checkpoint/resume rules, recovery paths and deliberate performance choices.
  Explain why a step is necessary or what must remain true; do not narrate
  obvious assignments, conditions or loops.
- Use conventional mathematical names when they make formulas clearer, and
  explain their meaning near the formula when it is not already evident.
- Make tests and benchmark runners easy to scan by separating setup, execution
  and validation. Keep comments focused on the case or experiment being tested.
- Keep readability-only passes separate from behavior changes. Preserve APIs,
  algorithms, defaults, seeds, work accounting and serialized formats; leave
  immutable source snapshots and historical evidence unchanged.
- Verify readability-only changes with tests, lint and an AST comparison that
  ignores source locations. Review any intentional AST differences explicitly;
  a passing test suite alone does not establish unchanged behavior.

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
scratch output and the detailed `v2/local/LOG.md` journal local. Do not add those
files by force. Preserve useful raw evidence locally or in durable external
storage; publish a concise summary with commands and limitations.

Before publishing, verify a checkout containing only committed files can run
the tests and import benchmark runners. Never include unrelated work in a
cleanup commit.
