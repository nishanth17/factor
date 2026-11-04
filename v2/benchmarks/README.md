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
| [qs/c9a/](qs/c9a/) | Independent useful-yield/admission/graph/record policies, immutable control, bounded captures/replays and staged complete-cost preparation |
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

## C9a preparation and staged execution

The waste-removal bundle on `f28fd61` has matched, validated component
measurements against its parent `57b40ad`. After at least three seconds of
warmup per mode, nine samples of 128-block prefixes at 50/60/70/80 digits
have 22–32% lower median whole-prefix time and identical serialized stores,
positions and divisors. The 80-digit prefix has no retained atoms; this does
not establish useful large-input yield or complete-factor improvement.
An actual-score scan has 55% lower median cost; the 26,862-atom historical
graph replay has 23% lower median time and about 1.2 MB less live GC memory.
These component gains are not additive. Local receipts:
`results/c9a/waste-components-67/`; reproduction:
`pypy3 v2/audit/c9a/run_optimization65.py measurement components-67 680
pypy3 -B v2/audit/c9a/measure-waste-67.py` under the shared quiet lock.
The receipt name is immutable and cannot be reused.

The worker also supports the separately selected `timing_policy` value
`retire-rewarm-v1`. It preserves the immediate legacy post-run sentinel,
then releases the completed job, collects garbage and performs at least
three seconds of validated warmup before a second nine-sample sentinel.
Both screens, all warm work and cleanup costs are recorded. The original
15% drift/spread rule remains unchanged. Historical freezes retain the
legacy policy; a revised policy requires a new prospective freeze and a
matched harness for every source. The optional launcher `execution_root`
permits a verified archived library with that same shared harness.

The new expanded arm is prepared by `freeze --extended` under
`c9a-expanded-yield-v4`, after executable correctness and committed-only checks.
The v2 and v3 protocols retain their identities. The expanded primary remains
one C0/C1 comparison with the same first-four-prefix stop and cumulative stage
ceilings; no run has been launched. Eight input sizes remain 50–80 digits, and
80 still has no complete-factor cell.

`qs/c9a/extended.py` provides six separately named reporting/SLP/endpoint/product
profiles, base/interval alternatives, a bounded observed-population Pareto
ranking function and a matched future-integration builder. Its initial
three-prime envelope admits the geometric B,2B,4B example within 8B³; this is a
populated starting hypothesis, not a calibrated optimum. Wider unobserved
reporting remains unknown in ranking. Changed bases and intervals require
their own resolved literals and eligible populations. Optional parameter work
uses a new decision/identity within the existing diagnostic allocation, not an
automatic Cartesian sweep or extra multi-hour campaign.
`freeze --extended --profile NAME` selects one of the six joint hypotheses;
its literal configuration, identity and fresh confirmation preserve that choice.
Observed credits are a screening aid; changed admission can change graph
closures, and a wider population must be measured before selecting it.
The separately identified `--profile bounded-report` challenger retains the
incumbent reporting allowance, proven three-prime admission and the expanded
mechanisms, while selecting bounded opaque staging and shape-based ownership.
It requires a new integration freeze; it does not revise the historical full
reporting receipts or increase their grants.
The subsequent ownership audit uses one unchanged reproduction and an end-state
census, followed by exact store replay. Its accounted-byte comparison holds
helper reservations constant; replay RSS and encoded snapshot size are not
factoring memory or speed measurements. Current in-memory bounded storage
charges shared slabs/contexts once and keeps packing/restore scratch explicit.
The separately identified `--profile cost-report` nominee keeps that envelope
and ownership policy, selects scalar XOR elimination and serial recovery, and
uses `rho-first-v3` with SQUFOF as a bounded fallback. It follows real-matrix,
observed-cofactor and shared/serial recovery
cost checks; it does not revise earlier `bounded-report` configurations or
establish a factoring win. Shared batching and panel solving remain selectable.

Built-in QS/SIQS now consumes private progress after each collection batch;
the public collector still provides full owned snapshots. The scoped return
replay compares these entry points under the same committed runtime and checks
identical stores, work and progress. Its mature-store cost reduction does not
establish a whole-factor gain; short almost-empty prefixes showed no benefit.
Keep JIT-off retention checks and explicit GC observations separate from
uninstrumented timings. See [scoped observations](RESULTS.md#c9a-preparation).

Subsequent diagnostics use PyPy's built-in statistical sampler with Python/JIT
frames and a pinned local reader. Calibrate overhead on identical work and
retain the resulting signatures; sampling shares are not removable-cost bounds.
Keep heap dumps, explicit GC and resource-observer interventions separate from
clean performance calls. A growing heap needs object/lifetime attribution;
RSS high-water alone does not establish a leak. The measured callback cost and
failed timing gates remain attached to the current scoped results.

The primary configuration uses serial recovery and in-memory storage. Threaded
hints and persistent spill are implemented selectable capabilities; exercising
them needs a separately frozen configuration and finite production-disk grant
in addition to diagnostic output caps. A DLP incumbent cannot replay a 3LP or
opaque stream as though it admitted that population. Future C10a integration
binds the same supplied committed base to both arms under a new identity; it
does not alter the active C10a worktree or adopt its current nominee.
`freeze --extended --integration PATH` accepts a preparation containing the
exact current `source_commit`, `corpus_sha256`, per-fixture serialized
`incumbents` and independently resolved `resolved` metadata. Both arms then
execute that same committed runtime with their literal configurations; C0 does
not load the historical control archive. `preparation.structural` accepts
per-size configurations for resolving that geometry, and fresh confirmation
uses the same per-size incumbent and a distinct input identity. These helpers
have not been executed for a C10a integration.

Preparation delivers code and certified inputs without running A/B/M/R/F/C.
The literal configurations and actual input multiplier/base bounds live in
`qs/c9a/extended.py` and `inputs/c9a/structural.json`. The former
`configuration.py` plan retains its historical identity. Detailed research,
the decision note and a concrete execution handoff are local in `audit/c9a/`.

```sh
pypy3 -m v2.benchmarks.qs.c9a --help
pypy3 -m v2.benchmarks.qs.c9a freeze --extended \
  --output v2/benchmarks/results/c9a/expanded-v4-freeze.json
```

Run the queued correctness and committed-only checks before freezing. This
expanded bundle's executable QA is deferred until C10a finishes. Freezing
requires committed, clean preparation and never starts a campaign.
A separately authorized execution uses `run --stage A --freeze PATH
--output DIRECTORY --decision DECISION.json`, then a recorded usefulness/cost
decision for each later stage. There is no run-all command or full 80-digit
cell. In the combined protocol, the first A invocation stops after four clean
60/65-digit anchor prefixes (410 seconds maximum worker lifecycle total).
`--continue-early` requires their recorded feasibility review and reuses them
inside the 32-call 50–80-digit primary screen. B1 can follow the four initial
calls directly; it does not wait for every larger prefix or a successful
control. The complete 60/62-digit B1 anchors, broader B2, regression R,
conditional 75-digit F and disjoint one-band confirmation remain separately
decided stages. No component-by-component wins are prerequisites.
`observer` performs nine matched one-block overhead pairs; optional
`replay` uses a maximum of two validated complete streams and eight slots.
Captured fixed populations cannot establish additional I/E yield.

The runner reserves warmup/end work, lifecycle time and output storage before
each worker, validates factors against parent-held certificates, preserves
censored/lost cells and separates recorder exhaustion from native limits.
`select-confirmation`, `prepare-confirmation` and `freeze-confirmation` select
one band before fresh sampling and preserve earlier charges. Fresh required
inputs are committed before their freeze. These commands do not grant
permission to execute the preparation-only task's future campaigns.

The combined envelope preserves the incumbent domain and proposes
`S=max(floor(sqrt(P0)),4B)`, `L=max(L0,16B)`, and `P=max(P0,32B²)`, clipped only
to the exact ordered endpoint product. SLP has a separate allowance and
reporting retains the incumbent ceiling. Each cell records resolved literals.
These settings make the boundaries explicit; no measured optimum is claimed.
The current mainline/control runtime equivalence is recorded in
`inputs/c9a/control_equivalence.json`. A later runtime integration needs a new
matched identity. This QS study is not global portfolio acceptance for uniform
60/80/100-digit integers; balanced 60/70-digit completion must not regress.

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
