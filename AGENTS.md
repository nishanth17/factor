# Contributor guidance

Explicit user instructions take precedence. Work on the active implementation
in `v2/`; preserve `v1/` as the original baseline. Read [current status](v2/STATUS.md)
first, then the [live backlog](v2/ROADMAP.md) and [stable policy](v2/POLICY.md).

## Code and arithmetic

- Target PyPy implementing Python 3.11. Do not silently fall back to CPython.
- Keep factoring arithmetic exact with integers, integer roots and inversion.
- Validate proper divisors and reconstruct every result, including unresolved
  cofactors. Keep probable-prime and proven-prime labels distinct.
- Preserve finite work/time/storage allowances and documented resume behavior.
- Keep library calls quiet unless explicit verbosity is requested.
- Implement factoring searches in this repository. GMP/gmpy2 arithmetic is
  allowed; external factorers such as GMP-ECM, FLINT, YAFU, PARI and msieve are
  reference benchmarks, never production factoring delegates. Keep their
  wrappers and invocations confined to reference runners.

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
- Write project documentation in neutral terms. State decisions and evidence
  limits without personal attribution.
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
samples for performance claims; extend unstable measurements only within a
declared finite plan. Smaller feasibility screens are exploratory and cannot
establish a speedup. Separate cold startup from warmed execution, and
instrumented profiles from performance evidence.

Document API changes in `v2/README.md`. Update roadmap completion only after
its acceptance and experiment gates pass. Summarize accepted behavior in
`CHANGELOG.md` and useful measurements in `v2/benchmarks/README.md`.

## Experiment design and prompt review

- Measure the mechanism on inputs it can affect. ECM studies use smaller-factor
  strata appropriate to the declared bounds and curve effort, including failed
  searches; bounds do not guarantee success for a given factor size. QS studies
  use QS-bound inputs. Portfolio adoption additionally needs the declared
  uniform random-integer panel unless an explicit exception is recorded.
  Unaffected classes are do-no-harm guardrails, not mechanism speedup claims.
- Keep the user's optimization goal central, including useful progress toward
  60–80-digit factoring. Before recommending a task or execution prompt, review
  whether its admission gates and expected cost can answer that goal. Following
  a frozen protocol does not excuse recommending an uninformative campaign.
- Push back before using the performance window if a requested experiment
  cannot measure the requested benefit. Explain what would answer it; do not
  silently substitute a different endpoint or workload.
- Establish a small, uninstrumented feasibility baseline on representative
  inputs before launching a large diagnostic or matched campaign. Use current
  committed mainline for current improvement decisions; label historical
  comparisons explicitly and create new identities for changed integrations.
- Measure diagnostic overhead on identical work with tracing enabled and
  disabled before committing to long instrumented runs. Prefer short targeted
  traces or bounded counters when inner-loop observation dominates. Instrumented
  completion times cannot stand in for clean factoring performance.
- Separate admission from adoption. Exact arithmetic, valid outputs, checked
  traversal/resume and feasible resource bounds govern starting a comparison.
  An unrelated reference arm need not factor every input or produce verified
  dependencies before a correct candidate can be measured. Missing explanatory
  attribution is uncertainty; it is not automatically a reason to forbid a
  direct end-to-end comparison.
- Require matched complete-cost/completion evidence, adequate independent
  inputs, uncertainty analysis and disjoint fresh confirmation for adoption.
  Do not impose a universal minimum percentage gain. Treat promising exploratory
  results as reasons to investigate, without converting them into accepted
  improvements or counting repeated seeds as independent moduli.
- Give each stage a concrete question, an advancement/stopping decision, an
  estimated total runtime and finite cumulative resources. Calibrate those
  estimates from the first informative observations. Reassess before spending
  further hours when a stage only repeats recorder failures, times out without
  answering its question, or blocks a useful comparison on an unrelated result.
  Do not repeat a whole campaign with larger limits as the automatic response.
- Before any run longer than a few minutes, state its question and which
  outcome would change a decision in one sentence. Report the first informative
  result before continuing a long campaign. If both arms time out without
  answering that question, redesign rather than automatically extending them.
  Censored outcomes can still establish feasibility or resource limits; lack of
  a demonstrated gain is not evidence of no benefit.
- Keep identical-work throughput, time to first validated factor and complete
  factorization separate. Distinguish fresh-interpreter cold calls from calls
  using a retained warm context. Warmup belongs outside cold samples; neither a
  worker-throughput gain nor a warm-pool gain proves a cold CLI completion gain.
- Distinguish diagnostic storage/event exhaustion from production memory/work
  exhaustion. Recorder failures must be identifiable as harness failures and
  cannot support claims about the factorer's resource needs. Validate the
  recorder's useful observation window before scaling the campaign.
- Preserve receipts, charges, failures and historical freezes. If a protocol is
  flawed, document the reason and prepare a separately identified revision
  within the user's authorized scope. Keep confirmation inference separate
  from exploratory redesign; never selectively retry unfavorable observations,
  relax correctness or silently rewrite a control to manufacture acceptance.
  Carry earlier expenditure into revised protocols under the same campaign
  allowance; a new run identity does not reset the budget.

## Calibration and cheap paths

- Declare numerical work first: bounds, curve counts and quotas. Hold that work
  and its search identities fixed while comparing worker counts or scheduling.
  Worker economics and useful effort can interact; explicitly designed, bounded
  joint tuning is allowed, not an unnoticed change during a scaling comparison.
- When numerical work, backend, scheduling, source or cold/warm ownership
  changes, identify every dependent cost table, engagement threshold and worker
  selection needing recalibration. Do not reuse a grid winner as a universal
  optimum or route using known factors or future outcomes.
- Keep serial/easy paths cheap. Reject ineligible calls before expensive imports,
  catalogue loading, source hashing or pool preparation. Measure decision cost
  and the not-engaged path against committed production mainline; falling back
  to serial does not by itself establish unchanged serial cost.

## Merge and cleanup

- Before merging implementation changes, run the full suite, lint and checks
  from a committed-files-only checkout, including benchmark runner imports.
  Targeted tests alone are insufficient. Reuse passing checks only when their
  source remains applicable; documentation-only changes need not rerun factoring.
  Import data-only historical snapshots through their loaders, not as ordinary
  benchmark modules.
- Before removing a worktree, preserve ignored audit/results, receipts and useful
  raw evidence in the main repository's ignored archive with hashes and source
  identities. Retain historical failures and charges. Verify the archive before
  deleting reproducible snapshots or caches; never follow an environment symlink
  when deleting a worktree.
- Check active branches and execution prompts for dependencies on removed paths.
  Update owned active references and preserve relocation maps for immutable
  receipts; coordinate shared consumers without overwriting another task's work.
  Inspect matching worktrees before resuming or removing them. Keep the shared
  performance lock file and unrelated worktrees intact.

## Known host and backend constraints

Reuse these facts rather than researching them for every task; still record
actual software versions and recheck capabilities when a changed host/toolchain
or a new backend depends on them.

- Host: Apple M4 MacBook Air, model `Mac16,13`, 24 GiB, 4 performance and 6
  efficiency cores, fanless. All ten cores are available, but useful concurrency
  is measured; worker IDs do not identify physical core types. See Apple's
  [model identification](https://support.apple.com/en-us/102869) and
  [M4 Air announcement](https://www.apple.com/newsroom/2025/03/apple-introduces-the-new-macbook-air-with-the-m4-chip-and-a-sky-blue-color/).
- Record bounded `pmset -g therm` output before and after each timed cell,
  outside the measured interval; charge collection to the campaign ledger.
  Preserve failures or missing readings. This limited telemetry does not prove
  the absence of throttling on the fanless host. Apply new collection rules
  prospectively; do not rewrite earlier frozen measurement definitions.
- NEON is the current portable SIMD target, not the host's only capability:
  `hw.optional.arm.FEAT_SME` and `hw.optional.arm.FEAT_SME2` both report 1.
  Do not assume SVE or x86 IFMA code is available. A supported instruction set
  is not an implemented or adopted Factor backend.
- The retained B4a two-curve multiprecision ECM NEON prototype lost its scoped
  comparison. Keep that loss scoped to its kernel, widths and complete costs;
  reopening ECM SIMD needs a concrete changed hypothesis. Near-term SIMD work
  prioritizes small-integer QS kernels, not a blanket claim that all ECM SIMD
  must lose.
- Do not assume Apple Clang supplies a usable OpenMP runtime. Use pthreads for
  native threads unless a compatible project-local OpenMP runtime is explicitly
  provided and validated. Do not install system dependencies as a silent fix.

## GitHub contents

Commit code, tests, documentation, independent corpora and immutable baselines
required by test/benchmark loaders. Required inputs live in versioned
`v2/benchmarks/inputs/`; generated evidence belongs in folder-ignored
`results/` trees. The `v2/audit/` directory is local and Git-ignored; maintain
the current entry point in `v2/STATUS.md`, the live acceptance plan in
`v2/ROADMAP.md` and stable rules in `v2/POLICY.md`. Keep detailed research and
comparison documents under `v2/audit/`; do not publish selective report or
result collections in `v2/studies/` or another documentation directory.
Keep concise decisions and limitations in `v2/benchmarks/RESULTS.md`, reusable
runner code under `v2/benchmarks/`, and external reference sources under
folder-ignored `v2/benchmarks/inputs/upstream/`. Do not add individual-file
ignore rules for study names or external projects. Public docs describe this
project's behavior; preserve source inventories and comparison provenance locally.

Keep generated captures, stdout transcripts, profiles, verification dumps,
scratch output and the detailed `v2/local/LOG.md` journal local. Do not add those
files by force. Preserve useful raw evidence locally or in durable external
storage; publish a concise summary with commands and limitations.

Before publishing, verify a checkout containing only committed files can run
the tests and import benchmark runners. Never include unrelated work in a
cleanup commit.
