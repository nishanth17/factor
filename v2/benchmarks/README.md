# Factor benchmarks

A7a integrates local master `755c312`'s H3 resource model on its isolated
branch. Current-worker resources use observed root wall/CPU/RSS, while the
numerical recipe and packet/cache capacities stay fixed. Committed-only QA
covers 1,628 tests (nine optional skips), lint and 295 benchmark imports; the
final formatting-only test change preserves its AST. Earlier timing tables
remain historical and need revalidation. This integration validation itself establishes no speedup. The subsequent
[slim G4](RESULTS.md#a7a-slim-g4-8-november-2026) reports scoped S0 comparisons. The subsequent
[worker race](RESULTS.md#a7a-reduced-worker-count-race-8-november-2026)
reports provisional per-size process counts, including drift ties and CPU/RSS. Local receipts: `v2/audit/a7a/master-integration-20261108/`.


Runners and required numerical inputs for PyPy implementing Python 3.11.
[Results and historical decisions](RESULTS.md) hold measurements;
[STATUS](../STATUS.md) owns current routing and next work.

## H3 correctness and recorder boundary

H3's fixed-work oracle is `infrastructure/h3_correctness.py`. Under the shared
lock, run it on committed base `40d2367` and the candidate: three fixed QS inputs
each at 50/65 digits plus two nonempty witnesses, nine fixed ECM cases and four
manifest assignments compare exact JSON signatures and scanned positions. This
is a correctness comparison; no durations establish a speedup. Source-bound
receipts, including failures, remain local in `v2/audit/h3_resource_limits/`.

`infrastructure/resource_control.py` keeps historical recorder work frontiers
and source-pinned grant metadata at an explicit adapter boundary. Production
constructors reject removed resource arguments. Frozen controls, source pins,
experimental accounts and immutable baselines are retained; successful imports
do not authorize replaying an old frozen campaign against changed code.
Revalidate C3a allocation/handoff, B12 selection/startup, A7a packet/cache and
B4a native engagement costs on the final combined source before new decisions.

## Before a measurement

Follow [policy](../POLICY.md#experiment-design-and-prompt-review): establish a
small uninstrumented feasibility baseline, measure observer overhead when
needed, and give every stage a question, stop decision and finite cumulative
resources. Freeze source/configuration, inputs, seeds and selection before
accepted timing. Use independent fresh confirmation and validate every output.

Use the policy's evidence tiers: start with three independent paired inputs,
extending to six and nine only under a funded, prespecified uncertainty rule.
Warm timings need at least three seconds of validated warmup; cold timings use
a fresh interpreter with no warmup. Full test suites share the measurement
lock; worktrees do not isolate CPU, RAM or thermals. Instrumented diagnostics
remain separate from clean performance evidence.

## Current-collector QS regime studies

The current [gated A7a plan](A7A_EXPERIMENT.md) supersedes the older broad
designs described below. G2 covers cold 40/50/55/60/65/70 and warm 40/50;
62 is excluded, and 75/80 remain separately approved post-G3 feasibility
probes. Fix measured architectural overhead, classify options, then tune only
real tunables at anchors and verify interpolation. Fresh G4 inputs live in
`inputs/a7a/confirmation_g4.json` and remain outside tuning.

The restricted G3 is now reported: only workers 4/6/8/10 and packets g/4,g
were varied at cold 40/55/70 and warm 40. Experimental process recipes are
C4/p4 at 40 cold/warm, retained C4/p16 at 55 cold, and C10/p64 at 70 cold.
The 55 ranking is inconclusive under control drift. All 100 current-source
calls complete; these calibration choices do not establish a fresh mainline
speedup or a routing table. Runtime/committed-only QA source is `60bfd95`;
full tests, lint and benchmark imports pass. Intermediate verification and G4
remain gated. Subsequent probe outcomes are reported below. G3 evidence is
local in main's ignored `results/a7a/g3-restricted-20261107-r2/` directory.

The independent-packet diagnostic accepts an explicit bounded `--seconds`
override of 35 or 75 seconds. It changes the observation allowance, not the
fixed packet work. Background QoS is not core pinning. Buffer payload,
accounted cache ownership and observed RSS are distinct memory measurements;
the byte/native screen does not establish a working-set scaling limit.

Probe cells may declare `rss_stop_bytes` below their admission allowance.
The RSS observer then cooperatively cancels work at that sampled live bound;
`preserve_checkpoint=True` includes the validated portfolio checkpoint in the
sample result. Recorder failure and a live-memory stop remain distinct.
The amended 75/80 probes use a 32 GiB admission grant, 1 GiB worker /
256 MiB runtime, 6 GiB parent, an 18 GiB live guard and 256 MiB checkpoint
caps; this does not change G3's historical grants.

The approved 75/80 single-input probes are reported in main's ignored
`results/a7a/probes-20261108/REPORT.md`. Both pass admission and stop before
filtering/solving with generic internal `memory_limit` outcomes; neither
live RSS nor wall/CPU allowance causes the stop. Exact sub-limit attribution
is unrecorded. Outer checkpoints retain unresolved inputs but the terminal
QS path clears the collected state. Feasibility is not established; no retry
or default change follows. Committed-only probe QA passes on `5526118`.

G2's corrected coverage is reported locally in
`v2/audit/a7a/g2-report-20261107/REPORT.md`. All serial/process cells complete;
literal mainline is incomplete at 60/65/70. This is exposed-input feasibility,
not a confirmed speedup or routing table. Bounded capture review uses explicit
parent/representative-worker targets; full capture review requires all owned
streams. Representative diagnostics and finalist profiles follow the current
plan, rather than attaching a profile to every search loser.

Experimental cells can set `worker_observe` independently of `observe`:
`observe` retains owned RSS monitoring, while `worker_observe` controls the
per-relation worker timing/backlog counters. Omitted `worker_observe` preserves
historical behavior. Measure this observation cost before using the lighter
setting in clean comparisons; do not relabel historical observations. Worker
timeline capture composes with requested CPU/heap and recovery observers.
Short process screens may declare `features.assignment_cpu` to fund every
in-flight packet without assigning tiny residual CPU leases to lookahead
work. Hold it fixed across compared recipes and keep the aggregate CPU grant
unchanged. A local CPU-lease stop is not aggregate CPU exhaustion.

The following runners remain available for historical protocol reproduction;
their design sizes and per-candidate capture rule are not the current campaign.

The opt-in [A7a search](qs/a7a/regime_search.py) defines bounded worker/packet
neighborhoods and compatible feature removals. Its 40/50/55/60/62/65/70/75/80
anchors have separate configuration and resource identities; a prefix winner
does not establish complete-factor ranking or justify interpolation. The
versioned regime corpus partitions three independent calibration and three
confirmation moduli per anchor before timing. Confirmation outcomes must stay
outside tuning.

The [joint starting design](qs/a7a/joint_search.py) supplies 32 configurations
for declared 16/32/64-polynomial families: four worker, packet, verification and
stop-policy levels plus nineteen binary options. All raw factor pairs have
balanced coverage; higher interactions remain confounded. Conditional
inactivity, such as snapshots without leases or packed events without a queue,
must be reported separately from raw coverage. A partial block does not inherit
the full design's coverage. This is an exploratory starting design, not a
completed per-regime search or a configuration adoption.

A conditional eight-row refinement covers pairs of seven local factors around
an explicit anchor: workers, packets, verification, packed events, fused scan,
CPU snapshot and periodic prime-power initialization. It retains fast lease/
reuse/semaphore settings without declaring them globally dominant. Fund broader
exploration separately and preserve the aliases and incomplete coverage.

Numerical experiments may supply explicit collector endpoint, product and
reporting bounds, including `small_prime_bound`. The collector validates the
resulting configuration. Changing a factor-base size does not automatically
retune these literals: a coupled alternative needs a separate policy identity,
new collection observations and recalibration of dependent worker selections.

[Paired cells](qs/a7a/profile_pairs.py) reserve each clean timing and its
diagnostic companion together. Diagnostics trace the factoring execution in
the parent and every owned worker at 10 ms, pause sampling for parent and one
worker heap, then resume into a separate immutable stream. Two-Hz owned RSS
samples complement GC-object and logical-buffer data. Short calls retain a
labeled teardown heap if a live request cannot arrive in time. These captures
explain costs; their timings do not rank candidates. Source, flags, input,
pair identity and file hashes accompany the capture index. Worker finalizers
flush profiles before cleanup, and failed or unissued diagnostics remain
explicit. Every new optimization candidate receives its companion, including
rejected candidates.

On macOS the RSS recorder reads resident bytes directly from libproc for the
parent, owned workers and spawn helper; other platforms retain the `ps`
backend. Results identify the backend and missing PIDs. Recorder errors, empty
samples and diagnostic RSS-limit stops fail observation validation rather
than passing as factorer cancellation. Recorder sample exhaustion remains
separate from exceeding the actor memory allowance.

Clean and diagnostic receipts retain bounded host VM and process-CPU snapshots
before and after each call, alongside thermal readings. Their differences can
expose host compression, swapping or other CPU workloads. These include other
applications and are not the factorer's owned memory or exclusive CPU cost.
Missing or truncated telemetry remains explicit. These
bookends are outside the measured cell interval and charged to the campaign.
Prospective macOS bookends also query Foundation's thermal state and low-power
hint through a bounded, unprivileged command. Raw output and parse failures
remain explicit. This state is not CPU frequency or temperature; nominal can
also indicate unavailable thermal awareness. It cannot establish no throttling.
PyPy warmup and thermal conditioning are distinct. Declare one-off versus
sustained workloads and any finite preconditioning before selection. Use
frequent matched controls; a failed repeatability gate stops configuration
selection under that protocol. Retain the entire series instead of correcting
away drift or choosing a favorable subset retrospectively.

Agreement between endpoint controls does not bound intermediate machine state.
A serial arm can cool the host before a closing parallel control. Retain the
per-call state and ordering, and distinguish initial and sustained conditions.
Fresh-process complete actor costs include startup and teardown; they are
separate from warmed library intervals and production-CLI selector overhead.

Multi-sample cells retain the same repetition count in their diagnostic
companion, so later warmed execution and GC behavior remain visible. The
whole clean/diagnostic reservation must fit before either is issued.
Capture accounting tolerates disappearing temporary metadata names during
atomic replacement; missing immutable capture files still fail validation.
If a launched call aborts before worker CPU can be fully observed, the
unknown portion of its CPU grant remains charged separately from measured
CPU. Failed receipts retain detailed errors and traceback provenance.

Targeted diagnostic cells can add `packet_counters=true` to combine the
existing bounded recovery/split counters with CPU/heap/RSS capture. Clean
companions retain the configuration label without installing these timers.
Sample intervals exclude warmup from the counter summary; unmatched packets
remain counted separately. Recovery contains splitting. The combined profile
leaves stop callbacks unchanged; standalone packet diagnostics can still time
sampled stop calls. These instrumented spans include GC and diagnostic
interruptions. Calibrate their overhead on identical work before a long use.
The parent also counts reported residual shapes after schema decoding.
These counts describe incoming records; checked parent admission still
validates their arithmetic and determines actual retention and useful rows.

Warm complete-factor cells warm the same bounded factoring entry point with
a short partial allowance. They do not prime a bare-prefix pool under a
different coexistence contract. Warmup charges remain cumulative; neither
the control nor a candidate must solve during warmup to enter a comparison.
The default warmup target remains three seconds. A separately identified
selection protocol may request ten or twenty seconds through
`warmup_seconds`; these bounded partial calls remain charged, and the grant
must leave measured execution and cleanup funded. Record the warmup target
with each cost table; process warming does not establish thermal stability.

## B12/F4 delivered defaults

The normal bounded factorer uses the [selected regime table](../STATUS.md#selected-ecm-regimes)
with ordinary QS fallback intact. Cold 30/50/70 cheap plans bypass policy I/O;
20-digit automatic calls stay serial. A caller-owned warm context is a separate
endpoint. Latest committed-only QA and admission-only repair checks are in
[results](RESULTS.md#b12f4-cold-serial-repair-and-component-closure-5-november-2026).
Detailed evidence stays in local ignored `v2/audit/b12_f4/`. Component closure
does not establish whole-portfolio or population-tail gains.

## ECM-only matched execution

`ecm/b12_f4/selection_fixed_child.py` supports first-factor measurements with
`MASTER-MATCHED` and `SELECTED` request arms. The former loads the declared
incumbent source and runs its unchanged serial ECM with supplied manifest curve
seeds. The latter applies the packaged cold/warm policy, including its serial
decline. Numerical schedules and curve identities must match; this control is
explicitly distinct from the incumbent's legacy random stream and stock CLI.
First-factor requests require `siqs=None` and `sss=None`. The runner never
calls portfolio factorization or fallback for this endpoint.

Execute only under an inherited shared performance lock. Cold requests require
one fresh interpreter; warm requests retain their selected context and record
validated warmup, nine measured calls and shutdown separately. Independent
factor/certificate validation stays outside the factoring request. Detailed
protocols, source inventories and all receipts remain local under `v2/audit/`.

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

`qs.c9a.source_compare` supports matched committed-library comparisons through
the same external harness. `archive_library` records runtime file hashes and
`load_library` verifies bytes before private-package imports; `prefix_cost`
includes setup, collection, checked stores and helper work. `measure_rotating`
requires at least three CPU/wall seconds of warmup per mode and retains nine
rotated samples. Complete calls use `launch_pinned`; the worker entry point is
`pypy3 -B -m v2.benchmarks.qs.c9a.source_compare --worker REQUEST.json`.
The request contains the frozen cell and optional archived library; the driver
binds source identity and validates results against the fixture. Do not
substitute current code into a historical identity. Run performance work only
under the shared quiet lock and a declared finite cumulative ledger.

The final continuation uses the versioned C9a corpus, exposed_60_0/exposed_70_0,
seed 104729, list scores, `presieve+unroll4+parity_two` and the frozen C9a
configuration with recompute_prime_limit=0. Its new 2,048-block comparison
(`ebb35fa` / `c8168dc`) has nine favorable samples at both sizes; the full60
pair takes 159.519/132.247s with identical mathematical counters. Both timing
screens fail, so the raw difference is exploratory. The separate 9,000s-per-arm
full70 policy comparison on shared `891eda5` completes in 2742.095/1842.083s
(C0/C9a), with mixed timing validity. Neither is population confirmation.
Final checks: 1,113 tests (9 skipped), lint and 243 benchmark imports.
Prospective plans, exact invocation receipts, pinned settings and worktree
evidence locations are in local `v2/audit/c9a/optimization-83-89.md`;
final generated evidence is under `results/c9a/hit-words-89/`. Immutable
identities must not be reused. This continuation launches no further campaign.

New C9a analyses use `observational-v2`: the old `timing_valid` field retains
its historical 15% heuristic for provenance but is a diagnostic warning,
not a veto on otherwise validated timings. `resource_valid` separately records
real production resource/recorder compliance; instrumented data stays separate.
`summarize(..., legacy_timing=True)` explicitly reproduces the historical
analysis exclusion. New reports retain canary-warning counts, independent-modulus
uncertainty and all failures. Calibrate future timing decisions with unchanged
code and matched references; no universal percentage replaces that analysis.
Post-correction mainline publication checks pass 1,114 tests (9 skipped), lint
and all 243 benchmark imports from committed files only.

The final `180576c` matched pair uses flat queues in both arms and the literal
prior incumbent/`cost-report` configurations. On one exposed 60-digit input,
raw C0/C1 times are 268.807/230.364s. Both exact factorizations finish within the
same 900s/work/storage grants; unchanged terminal semantic counters and work
are verified against each arm's prior run. Both timing screens fail, so the
14.3% raw difference is qualified and is not adoption evidence. See
`results/c9a/final-pair-76/` and `v2/audit/c9a/verification-77.json`.
Historical invocation under the shared quiet lock:
`pypy3 v2/audit/c9a/run_optimization65.py execution final-pair-76 1900
pypy3 -B v2/audit/c9a/final-pair-76.py`.
The immutable identity cannot be reused. Final committed-only checks pass
1,084 tests (9 skips), lint and 242 imports. No additional campaign follows.

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

The subsequent whole 60-digit comparison on `359e935` finds correct factorizations
in 316.0/295.6/250.5s for the archived incumbent, optimized incumbent and C9a.
Both incumbent screens pass the revised timing policy; C9a's stable post-run
sentinel is 17.4% faster and fails the unchanged 15% drift gate. Its raw result
remains qualified. A targeted actual-event queue replay then rejects packed
64-bit fields: 77.4ms versus 42.2ms for tuples. Flat unsigned-word triples take
41.3ms with compact storage and replace that representation.
The distinct experimental `--profile scan-report` adds the existing exact
fused scanner to `cost-report` while preserving its other settings and all
pre-existing kernels. `cost-report` and automatic configurations retain their
policies. Nine warmed integration samples at each size preserve exact stores,
positions and divisors. At 60 digits flat queues improve 8/9 observations over
packed fields; adding fused scanning improves 7/9 over flat-only, while its 70-digit
median is worse. Its declared selection gate fails and is preserved. Keep
`scan-report` optional; the final current-source pair uses `cost-report` with
flat queues in both incumbent and candidate. This shared change requires a new
matched pair rather than comparing a flat candidate to an old packed control.
Local records: `results/c9a/whole-cost-70/`, `queue-replay-73/` and
`v2/audit/c9a/integration-74-plan.json`.

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
[pinned Git history](https://github.com/nishanth17/factor/blob/23458a52e409add326c9fd63cf1f322c6e605569/v2/history/README.md). Its retired C10a-04 orchestration
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

Final committed-only default-integration checkout `094ab1d` passes 1,337 tests
(9 optional skips), full lint and 264 benchmark imports. Its runtime files are
implementation `0c4cd18`; later changes fix legacy test-fixture intent only.
All source/default/API/checkpoint changes are validated. No new-source
performance measurement or balanced-completion claim is made. The shared
computation lock is released.

### A7a gated preparation

The isolated A7a branch provides `coverage_plan.prepare(identity, reference_source,
folder)` to construct the proposed G2 manifest without factoring. The reference
source is a full committed mainline hash, not a working-copy path.
`stage_driver.validate(plan, approval)` checks exact source/corpus bindings,
complete endpoints, coverage, configuration caps and whole-stage reservations.
The approval object names `stage`, `plan_sha256=stage_driver.digest(plan)` and
`approved=True`; it must reflect explicit approval of that exact stage.
`stage_driver.execute(plan, approval, output)` records issued attempts before
running and refuses automatic retries of unfinished work. It starts no successor.
The separately certified `inputs/a7a/confirmation_g4.json` corpus is reserved
for three/six/nine-input confirmation; generation is exclusive and bounded.
Benchmark observers record parent stages, worker service/queue spans and
minute-level progress. After the G1 calibration and explicit addendum, worker
active/queue-wait/publication spans are required on both arms of G2H/G3
process experiments and disabled for G4 confirmation. The stage driver,
combined runner and direct cell entry point reject conflicting settings.
Historical cells retain their actual observation settings.

G2H includes independent packet references through ten workers, matched
family-affinity comparisons, exact grouped marking charges and an explicit
byte-stamp collector selector. The requested byte/native remeasurement uses
complete cold factorizations at 50 and 65 digits, three independent paired
inputs per band and matched factoring/interpreter hash seeds. Clean timings
exclude CPU/heap captures; timeline tables and full evidence remain local.
See [the execution specification](A7A_EXPERIMENT.md) for the protocol and
gates. No automatic factoring defaults are changed.

G2H closes with a separate intermediate/final fixed-work bridge: three inputs
at each of 50/65 digits, three interleaved repetitions per source/input, exact
relation signatures and scanned positions. Its requested Tier 1 criterion
checks aggregate collector-time ratios within roughly 3%; cold actor times
and default-path diff touches are also reported. This is not a statistical
equivalence or full-factor speedup claim. G3 remains separately gated.
