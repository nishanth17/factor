# C10 implementation validation

Validated source: `a4886ee`, following implementation/research commit
`f4583cd`, on `codex/c10-collection-strategies`. C9 runtime lineage is
`129cb09e4ca684f3eb08cbbdbea3fb43848e5d95` / frozen runtime `fbfd662`.
No C10 integrated timing has run at this preparation stage. C9's final
1,200-second-control decision is still pending.

On 11 October 2026 UTC, C9 explicitly released a bounded 900-second QA slot.
C10 acquired the shared performance lock for each heavy command and returned
the machine after 354.984 seconds of QA. The experiment ledger conservatively
prepays the full 900-second QA allowance, including these failures/retries.

Final committed-files-only checkout:

- PyPy 3.11.15 / PyPy 7.3.23 with gmpy2 2.3.1.
- `make -C v2 test PYTHON=<project PyPy environment>/bin/python`: **613 tests,
  no failures or skips**, 93.500 seconds reported by unittest.
- `make -C v2 lint`: ruff checks/format and pycodestyle pass.
- All **129 benchmark runner modules** import; the A7 current-source loader
  passes. Immutable baseline source files are loaded through their existing
  adapters/tests rather than imported as standalone package modules.
- Checkout was produced by `git archive HEAD`, without ignored captures,
  profiles, source caches or local journals. Only the external development
  environment was linked for lint/dependencies; required data were committed.

Working-tree verification included the full system-PyPy suite (608 tests,
8 optional-GMP skips), then the project PyPy/GMP suite (612 tests, no skips).
The subsequent confirmation-manifest test and GMP-ordering test passed in the
targeted suite before the final 613-test committed-only run. Initial lint
found eight overlong lines and one import-order issue; those were corrected.
The initial broad import script mistakenly included immutable baseline source
snapshots and failed on their adapter-relative imports. Restricting the check
to actual runners, while retaining baseline-loader tests, resolved that
verification-script error. Neither event is factoring performance evidence.

Coverage includes actual-A batch permutations/uniqueness, random access and
extension prefixes, one-factor and final partial batches, exhausted work
before and during private cache construction, old/new cache coexistence,
charged reconstruction, finite/growing checkpoint replay and policy mismatch
rejection. Exact normalized-polynomial checks enumerate exceptional roots,
mod-8 parity, nonsingular powers of two, signed reconstruction, multiplier
gcd exceptions and intentionally lost threshold reports. GMP and Python-int
produce identical product ordering.

Local raw logs are under ignored `v2/benchmarks/results/c10/qa/`. No merge or
push was performed. The primary checkout changed independently during this
work (C3 was integrated there); C10 did not edit it or either source worktree.
The final C3 worktree was no longer present during the read-only status check.

This accepts the bounded implementation's correctness checks. It does not
establish a speedup, optimal batch/cache size, a new default, fresh complete
factor acceptance, or completion of the continuing C10 workstream.
