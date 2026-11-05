# Factor benchmarks

Runners and required numerical inputs for PyPy implementing Python 3.11.
[Results and historical decisions](RESULTS.md) hold measurements;
[STATUS](../STATUS.md) owns current routing and next work.

## Before a measurement

Follow [policy](../POLICY.md#experiment-design-and-prompt-review): establish a
small uninstrumented feasibility baseline, measure observer overhead when
needed, and give every stage a question, stop decision and finite cumulative
resources. Freeze source/configuration, inputs, seeds and selection before
accepted timing. Use independent fresh confirmation and validate every output.

Accepted timings need at least three seconds of validated PyPy warmup and
nine observations under the machine-wide performance lock. Full test suites
share that lock; worktrees do not isolate CPU, RAM or thermals. Cold startup
and instrumented diagnostics remain separate from performance evidence.

## Commands

From the repository root:

```sh
make -C v2 benchmark WARMUP_SECONDS=3 REPETITIONS=9
make -C v2 benchmark-phase-two WARMUP_SECONDS=3 REPETITIONS=9
```

Focused Make targets include `benchmark-phase-three-siqs`,
`benchmark-phase-three-sss` and `benchmark-phase-three-parallel`.
Read the chosen runner before starting it; a command here does not authorize
a new campaign. Historical protocols retain their original source pins and
accounts. Replay from the pinned revision, or freeze a new current-source
comparison with its own identity.

## C10a-01 integrated handoff

The [coupled-default decision](RESULTS.md#c10a-01-coupled-polynomial-defaults)
records calibration, complete-call feasibility and limitations. Reusable runners
and protocols live in `qs/c10a01/`; required numerical inputs and literal
controls remain versioned. Historical measurements use candidate `4b9a616`
and control `0076e256871722d52d849ab344152e0a7535c2ce`; the 65-digit terminal
handoff is `d37b809104c9bfd727c0ba0994121c9060a322ec`. Replay a frozen
experiment from its pinned checkout. Its source guards intentionally reject
current automatic-default substitutions. The separate A7 compatibility
manifest preserves explicit E1 arms and grants; it launches no E1 execution.

Detailed research, source inventories, raw results, interrupted calls and
cleanup manifests are retained locally at
`v2/audit/archives/c10a01-20261103/`. This archive includes the verified evidence
chain, relocated upstream material and temporary-state captures. Historical
paths in sealed receipts are preserved; the relocation manifest maps them to
the archive. Worktree deletion does not discard that evidence.

## Find an experiment

| Folder | Contents |
| --- | --- |
| [suites/](suites/) | Phase-one/two suites, validation, regressions and shared corpus builders |
| [pm1/a6/](pm1/a6/) | Pollard p−1 tuning, reuse and production comparisons |
| [ecm/p41/](ecm/p41/) | Verified PRAC kernels and full ECM campaigns |
| [ecm/p52/](ecm/p52/) | Programs, pairing, wheels and realistic ECM portfolios |
| [ecm/b3/](ecm/b3/) | Production chain integration, coverage and recovery |
| [ecm/b4/](ecm/b4/) | Kernels, integration and candidate bakeoffs |
| [ecm/c3/](ecm/c3/) | ECM allocation research, protected handoff and frozen complete-call study |
| [ecm/c6/](ecm/c6/) | Chain catalogs, executors and cost studies |
| [qs/phase_three/](qs/phase_three/) | Relation engines and collector/pipeline studies |
| [qs/p38/](qs/p38/) | Capacity, collector and matrix experiments |
| [qs/b1/](qs/b1/) | QS/MPQS/SIQS calibration and width studies |
| [qs/c1/](qs/c1/) | Large-prime feasibility and implementation studies |
| [qs/c9/](qs/c9/) | Source-derived large-prime yield, adaptive storage and complete-factor economics |
| [qs/c10/](qs/c10/README.md) | Nine-source collection study, exact polynomial oracles, scoped 40/50-digit plan acceptance and C10a follow-ups |
| [qs/c10a05/](qs/c10a05/) | Fixed collection/solve pacing, six-source reconciliation, finite resume oracles and scoped 50-digit complete-cost acceptance |
| [qs/c10a05_plus/](qs/c10a05_plus/) | Frozen 51-digit default-scope extension, 55/60-digit feasibility limits and charged fresh completion evidence |
| [qs/a7/](qs/a7/) | Frozen comparison arms and reconciliation |
| [primality/a10/](primality/a10/) | Primality dispatch and validation |
| [infrastructure/](infrastructure/) | Arithmetic backends, parallelism, sieves and performance audits |
| [support/](support/) | Snapshot loaders and independent arithmetic oracles |

## Package guides and evidence

[ECM](ecm/README.md) and [relation-engine](qs/README.md) guides explain those
study families. [Inputs](inputs/README.md) describes versioned corpora,
controls and immutable baselines. Detailed research stays local under
`v2/audit/research/`; public [results](RESULTS.md) summarize decisions,
commands and adoption limitations. Detailed reports and generated result JSON
stay local. Reusable runners belong here; required inputs remain versioned.

Write new captures, stdout, profiles and verification dumps under ignored
`results/<subject>/<experiment>/<run>/`. Make targets create timestamped
outputs; direct runners may require an existing output folder. Do not
force-add generated results, local journals or external source checkouts.

## Completed A1b review

The [final decision](RESULTS.md#a1b-existing-fixed-score-review) retains current defaults.
[Original protocol](qs/a1b_protocol.md) and [55-digit continuation](
qs/a1b_continuation_protocol.md) preserve literal inputs, source freezes and
charges, including prior failures. These are terminal historical runners;
restore their archived evidence for audit, and do not rerun measurements.
Source-only verification requires the frozen continuation checkout at
`cffe8c9`, including its original source paths and pinned interpreter:

```sh
pypy3 -B -m v2.benchmarks.qs.a1b_continuation verify
```

## Structural migrations

The former `v2/studies/` collection is preserved in
[pinned Git history](../history/README.md). Its retired C10a-04 orchestration
replays from the original revision with unchanged source freezes. Only the
tested schedule and preliminary decision helpers remain in
[c10a04](qs/c10a04/README.md); this cleanup changes no experimental verdict.

[layout_migration.json](inputs/controls/layout_migration.json) records moved
paths, source hashes and exact relocated bytes. Accepted structural migrations
are distinct from a changed experiment. Frozen loaders retain private module
names and recorded compatibility chains; newer production defaults do not
repin an old measurement. Current modules use the grouped package paths shown
above. For syntax without a measurement, use a runner's `--help`.

The [earlier guide](https://github.com/nishanth17/factor/blob/24ff153d26ada2ec4beeead07000d183e94fd355/v2/benchmarks/README.md)
preserves older command/layout notes and numerical summaries.

B12/F4 preparation provides `ecm/b12_f4/runner.py`, `inputs.py` and the prospective
`inputs/controls/b12_f4_preparation.json`. The driver
requires a later runtime-admitted freeze and separates fresh CLI complete cost
from retained-context calls. S0 uses the pinned committed incumbent; S1/P2/P4
expose the new executor and M0 is its separate serial manifest equivalent.
Detailed source review, protocol, budgets and queued QA handoff are local under
`v2/audit/b12_f4/`. Runtime source `8a2bff7` passed QA and exploratory startup
and larger-work screens; B12/F4 and E1 remain open. See the scoped observations
in [RESULTS](RESULTS.md#b12f4-larger-work-headroom-screen-4-november-2026).

`python -m v2.benchmarks.ecm.b12_f4.headroom FREEZE RECEIPT_DIRECTORY` analyses
fixed-work receipts without importing factoring code. It checks declared curve
coverage, proper divisors and equivalent search results, and keeps warm-call
latency separate from whole-cell CPU including warmup. Failed/missing cells
remain visible. It does not turn a one-observation screen into a performance
claim. Use PyPy 3.11 for the runtime driver.

B12/F4 option preparation now includes bounded packet size, action polling,
private metadata/decoded caches and checked stage queues. The runner accepts
explicit per-cell `settings`/`rule` overrides from a later pinned freeze;
`P<count>` forced diagnostics must match that cell's declared count. Wider
counts still require explicit settings and a separately funded protocol.
`inputs/controls/b12_f4_options_preparation.json` is prospective and cannot be
executed by the runner. Runtime Q precedes any option comparison.

`python -m v2.benchmarks.ecm.b12_f4.worker_sweep RECEIPT_DIRECTORY --output FILE`
keeps cold CLI lifetime and warm per-call results separate when selecting
exploratory count nominees. It retains failed/censored rows and excludes
instrumented arms from selection. Repeated timings on one modulus do not create
independent inputs or establish a portfolio optimum. `profile_runner.py` requires
an inherited held performance-lock descriptor; observers and clean timings have
separate identities. Detailed local CPU/heap findings: `audit/b12_f4/profiles-v1/`.

`ecm/b12_f4/kernel_replay.py` compares selected GCD/segment/kernel settings on
common explicit curve seeds while binding replay IDs to both configurations.
It is a separate diagnostic contract, preserving the production RNG stream.
`parent_profile.py` records bounded coarse inclusive/exclusive times for parent
planning and repeated metadata work; compare its observer cost on identical
work before interpreting it. Instrumented times are not performance evidence.
Both tools require a held inherited performance lock and frozen finite requests.

B12/F4 reuse evidence separates genuinely fresh component CLI lifetime from
retained-pool calls. Its local entry point is
`v2/audit/b12_f4/headroom_opt_v1/HANDOFF.md`; concise measurements and limitations
are in [results](RESULTS.md#b12f4-bounded-metadata-and-import-reuse-4-november-2026).
`kernel_replay` preserves explicit common curve seeds across parameter changes
without rewriting production manifest identity. `parent_profile` bounds coarse
observation and must be paired with a clean arm. Censored diagnostic calls keep
elapsed durations and failure charges; they never rank as fast complete calls.

For capture coverage audits, `dump_analysis.heap_summary(..., limit=None)` and
`cpu_summary(..., limit=None)` retain every histogram entry within the existing
input-size bound. Default summaries still display 25 entries and now expose
`omitted_gc_bytes` / `omitted_leaf_samples`; a truncated display cannot silently
stand in for a full histogram. These counts do not resolve inlined JIT operations,
allocation lifetimes or heap-root ownership. The current local coverage record is
`v2/audit/b12_f4/capture_coverage_v1/REPORT.md`.


The current F4 fast-path evidence uses resource/checkpoint v6 and the existing
explicit per-cell settings interface: `budget_mode`, `stage_two_cache_bytes`,
`bind_context`, `packet_curves` and `packet_policy`. All default off/reference.
Current source `b62ac7a` passes full committed QA. The local current entry is
`v2/audit/b12_f4/headroom_impl_v2/HANDOFF.md`; [results](RESULTS.md#b12f4-checked-fast-paths-and-complete-cost-screen-4-november-2026)
separate warm component gains from inconclusive cold CLI and the small factoring
pilot. Before increasing curve quotas, preflight `settings.checkpoint_bound(count)`
and simultaneous owned buffers; a large time/work grant cannot fund a missing
checkpoint reservation. Historical freezes and refused cells are immutable.


The larger-quota follow-up is recorded at `v2/audit/b12_f4/long_quota_v2/HANDOFF.md`.
It uses the same committed runtime with a checked 5 MiB checkpoint grant and
blocking-wait cold timing. [Scoped results](RESULTS.md#b12f4-512-curve-cache-amortization-4-november-2026)
keep the prior F4 control, current options-off and packet/cache arms separate.
A recorder permission error remains fully charged, with explicit late exit proof;
no unfavorable timing is selectively rerun. Warm reuse benefits do not establish
first-factor or complete-factor latency gains. Root-lineage-aware complete-call
feasibility precedes engagement calibration and fresh confirmation.


Current selection-integration work uses source-pinned training input
`inputs/controls/b12_f4_selection_v2.json`; its cost/settings anchors are unchanged
from v1, with bounded same-path size interpolation and incumbent serial portfolio
declines. Frozen v1 sources/inputs/results remain separate. The subsequent targeted
comparison below measures the revised policy; old end-to-end timings remain
under their original source identities.
The local `v2/audit/b12_f4/selection_integration_v2/HANDOFF.md` records focused
validation and passing full committed-checkout QA (`4a2f051`: 1,308 tests,
9 optional skips, lint, 262 imports and source-pin admission). The computation
lock was released at that QA checkpoint; subsequent targeted validation is
recorded below.

Current complete-cost comparison entry point:
`v2/audit/b12_f4/master_comparison_v1/REPORT.md`, `ALL_REGIMES.md` and `HANDOFF.md`.
All 33 recommended regimes were compared with committed master `0efaba1` using
candidate runtime/benchmark `538462a`. Stock numerical fallback and QS checkpoint
settings remain literal. Separate cold CLI and retained-context endpoints cover
4,384 process cells / 6,568 timed calls, preserving refusals and uncertainty.
See [results](RESULTS.md#b12f4-actual-master-complete-cost-comparison-5-november-2026).

Reusable `ecm/b12_f4/selection_calibration.py` fits cold/started/warm costs;
`selection_analysis.py` resamples pool/modulus clusters and preserves censored
quantiles. `selection_inputs.py` generates bounded certified enriched inputs
separately from unconditioned uniform panels. `selection_child.py` runs the
source-pinned process request and preserves work/CPU charges even when result
construction raises. Optional benchmark refusal recording reconstructs the whole
input as unresolved rather than inventing partial progress. Runtime remains PyPy
3.11; controllers and text analysis are separate from factoring measurements.

Committed loader inputs include `inputs/controls/b12_f4_selection_v1.json`,
`b12_f4_selection_resolved_v1.json`, `b12_f4_master_comparison_v1.json`,
`b12_f4_master_uniform140_v1.json`, `b12_f4_master_uniform140_resource_r2.json`,
and `inputs/corpora/b12_f4_master_comparison_v1.json`. The resource revision
preserves prior charges and changes no per-call treatment. Preserve earlier
confirmation/R64 inputs and evidence under their original identities.

Final committed-only checkout `20e6681` passes 1,287 PyPy tests (9 optional
skips), full lint and 262 benchmark imports. Its receipt and exact commands are
in the local handoff. The computation lock is released. No additional performance
campaign or default promotion is queued by this report. B12/F4 remain open and
E1 stays deferred. That source leaves exact-cofactor coverage and the balanced
guardrail unresolved; the targeted comparison below addresses cofactor coverage
without establishing the balanced guardrail.

The targeted cofactor-policy comparison uses committed benchmark/helper source
`9b7eaeb` against actual master `0efaba1`, with the native implementation
unchanged from `10f64d5` / `4a2f051`. All 1,248 cells / 1,838 timed factor calls
validate. Nine fresh cold samples and three warm pool lifetimes × three samples
per primary modulus preserve separate endpoints and full failed-call charges.
See [fresh results](RESULTS.md#b12f4-fresh-cofactor-policy-comparison-5-november-2026)
and the local `v2/audit/b12_f4/targeted_policy_v2/REPORT.md`, `TABLES.md`,
`resources-r2.json` and `HANDOFF.md`. The lock is released; no further run or
default promotion is queued automatically.

Reusable `ecm/b12_f4/targeted_inputs.py` adds bounded certified calibration
primes, direct/recursive enriched inputs, unconditioned uniform integers and
cheap-exit controls. It requires PyPy 3.11, a frozen policy, a new output path
and an inherited held shared-lock descriptor. Its CLI is:

```text
pypy3 -m v2.benchmarks.ecm.b12_f4.targeted_inputs \
  --phase confirmation --seed SEED \
  --policy-freeze NEW_FROZEN_POLICY --output NEW_CORPUS
```

Generation/certification is runtime work: reserve it before execution and never
replace an exposed corpus. The committed loader inputs are
`inputs/controls/b12_f4_targeted_policy_v1.json`,
`b12_f4_targeted_policy_v2.json`, `b12_f4_targeted_selected_v1.json`,
`b12_f4_targeted_confirmation_v1.json`, and
`inputs/corpora/b12_f4_targeted_calibration_v1.json` /
`b12_f4_targeted_confirmation_v1.json`. The selected profile was committed
before fresh confirmation generation; only its 81-digit anchor was fitted.
Committed-only `9b7eaeb` passes 1,311 tests (9 optional skips), lint and 263
benchmark imports. Later report changes are documentation only. Historical
freezes, refused stock checkpoints, all charges and analysis corrections remain
local and immutable. B12/F4 stay open, including the balanced guardrail.

Automatic default integration now uses packaged production policy
`v2/execution/ecm_worker_policy.json`, copied from the frozen selected profile
without changing costs or nominees. Current direct library contexts can use
`v2.execution.ecm_auto.selected_context`; old source-pinned benchmark profiles
still belong to their historical commits. Use explicit `ecm_executor=None` for
serial numerical oracles and controls. The normal CLI's `--ecm-workers serial`
disables process selection while keeping declared default tiers/grants.

No new timing claim accompanies this wiring. The [adoption record](RESULTS.md#b12f4-automatic-default-integration-5-november-2026)
retains pending guardrails and cold/CPU losses. Local validation/charges:
`v2/audit/b12_f4/default_adoption_v1/`. Historical freezes stay immutable.

`ecm/b12_f4/adoption_profile.py` is a static maintenance helper. It transfers
unchanged measured rows with new integration source pins, keeps the reference
immutable and labels transferred costs as training estimates. Source repinning
is not new performance confirmation. After a deliberate code change, use a
reviewed measured reference and expected SHA, package into the production JSON
and update exactly the adapter policy constant, then perform locked correctness
QA. Existing automatic checkpoints require their original committed source.
