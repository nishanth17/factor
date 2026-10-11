# C10 implementation validation

Validated source: `00a9d9e`, following implementation/research commit
`f4583cd`, on `codex/c10-collection-strategies`. C9 runtime lineage is
`129cb09e4ca684f3eb08cbbdbea3fb43848e5d95` / frozen runtime `fbfd662`.
No C10 integrated timing had run at this preparation stage. The completed
[tranche 1 results](tranche1_results.md) are recorded separately. C9 completed
its 1,200-second-control study, accepted the measured 50-digit growing-DLP
bundle and explicitly released the CPU. Its runtime and all configuration
dictionaries match this checkout's controls. The three exclusion-source
file hashes were reverified at handoff.

On 11 October 2026 UTC, C9 explicitly released a bounded 900-second QA slot.
C10 acquired the shared performance lock for each heavy command and returned
the machine after 354.984 seconds of QA. After C9's final CPU release, another
100.526 seconds covered the latest committed-only validation below, bringing
actual QA to 455.511 seconds. The experiment ledger conservatively prepays
the full 900-second QA allowance, including these failures/retries.

Preparation committed-files-only checkout:

- PyPy 3.11.15 / PyPy 7.3.23 with gmpy2 2.3.1.
- `make -C v2 test PYTHON=<project PyPy environment>/bin/python`: **615 tests,
  no failures or skips**, 93.924 seconds reported by unittest.
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

The final runner checks also reject a partial 12-sample extension and an
unbalanced nine-sample allocation as performance acceptance evidence.

Local raw logs are under ignored `v2/benchmarks/results/c10/qa/`. No merge or
push was performed. The primary checkout changed independently during this
work (C3 was integrated there); C10 did not edit it or either source worktree.
The final C3 worktree was no longer present during the read-only status check.

This accepts the bounded implementation's correctness checks. It does not
establish a speedup, optimal batch/cache size, a new default, fresh complete
factor acceptance, or completion of the continuing C10 workstream.

## Final committed checkout and acceptance checks

After committing the frozen controls, selected arms, certified fresh corpus,
source reconciliation and result documents at `a13dd3b`, a new `git archive
HEAD` checkout ran the project PyPy 3.11.15 / gmpy2 2.3.1 suite: **615 tests,
no failures or skips**, in 94.201 seconds reported by unittest. `make -C v2
lint` passed; all **129 benchmark runner modules** imported and the A7 current
loader passed. The ignored local QA ledger showed 556.360 seconds of actual QA
through that check. A further archive of documentation commit `ef3e6c1`
ran **615 tests without skips or failures** in 93.490 seconds, passed lint and
imported the same 129 modules plus A7. The ledger then totaled 656.366
seconds, below its conservatively prepaid 900-second aggregate charge.

The exact frozen source/configuration control and training selection still
load, the new corpus's independent prime certificates verify, and each fresh
40/50-digit comparison has nine balanced samples and passes its scoped complete
factor gate. This final validation note and the roadmap checkbox are
documentation-only changes after that committed-checkout run; runtime,
benchmark Python, controls and corpora remain byte-identical. No merge or push
was performed.

## Subsequent C9/C10 default integration verification

The later user-directed integration commits the scoped selector and new live
A7 manifest at `9935b3c`, on top of mainline C9 default routing `7f1c38f`.
It does not modify C10's frozen controls, corpus, timing runner or original
result. The selected production configurations exactly match both frozen C10
`plan` records; setting only `power_plan_bytes` back to zero recovers their
`control` records. C9's frozen 50-digit challenger is likewise unchanged
except for C10's plan allowance. The live A7 manifest chains to C9's manifest
and validates 25 current source hashes; the historical C10 timing loader
correctly refuses this changed production runtime.

The working-tree `make -C v2 test` passed 718 system-PyPy tests with nine
optional backend skips; `make -C v2 lint` passed. A `git archive HEAD` checkout
of `9935b3c`, with only the external project PyPy environment linked, passed
**721 PyPy 3.11/gmpy2 tests without failures or skips**, full lint and imports
of all 129 benchmark runner modules plus the A7 loader. A first sandboxed
committed-only run had one environmental `ps` permission error in the existing
C1 quiet-process harness; rerunning with process-list access passed. That
error was not a factoring failure. One certified fresh 40-digit and one
50-digit input also completed through omitted-config `SIQSJob` defaults with
proper divisors and exact reconstruction. Their single-call elapsed times are
functional smoke data, not a new timing comparison.

This receipt closes the delivered bounded C10 tranche after its original
nine-sample fresh gates and the subsequent default-route checks. The distinct
mechanisms and their revisit triggers remain under C10a or existing owners;
there is no claim that every implementation or size policy is adapted.
