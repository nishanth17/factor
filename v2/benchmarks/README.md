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
`inputs/controls/b12_f4_preparation.json`. No measurements have run. The driver
requires a later runtime-admitted freeze and separates fresh CLI complete cost
from retained-context calls. S0 uses the pinned committed incumbent; S1/P2/P4
expose the new executor and M0 is its separate serial manifest equivalent.
Detailed source review, protocol, budgets and queued QA handoff are local under
`v2/audit/b12_f4/`. Runtime QA is the next phase; B12/F4 and E1 remain open.
