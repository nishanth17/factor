# Factor live backlog

Start with [STATUS.md](STATUS.md) for current defaults, work in flight and next
eligible items. [POLICY.md](POLICY.md) holds stable rules. Updated 2 November
2026 from committed mainline `41997da` and C10a-01 roadmap review `4d65d33`.

Status vocabulary: **Done**, **Open**, **Conditional** (requires its trigger),
**Deferred**, **Rejected**, **Inconclusive**. Done describes the delivered scope;
it can retain an Inconclusive experimental verdict. Adoption, retention as
opt-in, and closure are recorded explicitly.

| ID scheme | Meaning and relationship |
| --- | --- |
| P1–P8, including P3.8-R1–R5 and P8.2.1/P8.2.2 | Capability/acceptance catalogue. Numbers preserve ownership and historical gates. |
| A–I execution layers, including suffixes | Dependency-based work packages across P phases. Each row maps to its P scope; layers are not a second list of capabilities. GNFS rows retain their letters while Deferred. |
| C10a-##, AP##, MG##, M## | Existing subtranches, source-transfer inventory, gap inventory and milestone evidence IDs. They retain their identities and linked historical meaning. |

## QS optimization priorities — coverage reconciliation, 2 November 2026

Source: C10a-01 worktree review `4d65d33`, including the priority update
`1d21bbf` and implementation handoff `f306359`. This integration brings over
roadmap guidance only. Implementation statements and the 859-test record below
refer to that worktree; its code is not integrated into mainline. Production
defaults and historical experimental verdicts remain unchanged. Detailed
research and handoff reports stay local.

This review updates the priorities and preparation dependencies of the
execution rows below. Their arithmetic, resource, comparison and adoption
gates still apply. Model assignments in this review govern the named new QS
work; other roadmap assignments retain their scope.

The objective is the fastest practical Python QS implementation under declared
runtime, correctness and resource conditions, with useful progress through
60–80 digits. This review adds concrete implementation nominees; it establishes
no new speedup or default. Detailed source pins, primary-paper coverage,
mechanism dispositions and access limits remain in local audit material.

Use the following additions with the existing task IDs. A prior isolated loss
retains its tested scope. An unrelated completed campaign is not a prerequisite
for an exact, feasible local comparison. E1 remains deferred. The older C10a-01
calibration is closed; its later 60/62/65/70 complete-factor stages, including
the funded 70-digit pair, remain unrun. Revised code needs a new source/work
freeze before measurement; this backlog allocates no additional experiments.

| Priority / owner | Concrete next deliverable | Acceptance and discriminating evidence |
| --- | --- | --- |
| High / C10a-01 with A1a | Precompute bounded modular Gray shifts while retaining centered B and its wrap correction. Consider small-factor batch-inverse table construction after the simple exact cache. | Same physical trials, direct CRT/brute-root checks, exceptional roots, quota-aware storage, cancellation and charged replay. Compare setup plus complete family cost. Existing inverse caching alone does not implement this. |
| High / C10a-01 with C10a-02 | Improve the reachable A-product population beyond paired-flyer minimization within a fixed core. Specify finite feasible-core/product-window selection; retain actual-A/norm-aware multiplier ranking as a separately selected extension. | Exact target reachability, unique families, finite exhaustion and norm/useful-relation distributions. Do not repeat retired quality-order, face or quota grids. Wider local-density ranking is already implemented. |
| High design / F7 with C9a | Integrate exact SIQS smooth-part screening with direct-versus-tree dispatch. Separate base-product capacity from candidate-batch capacity and retain scalar exponent recovery. | Current base endpoint and tree-node limits are both barriers. Prove simultaneous ownership and pending-candidate resume; compare identical candidates against current bucket recovery, then whole cost. Residual-prime batch cofactoring is a separate domain. |
| High / A1c with A1a/B5 | Select a Python bulk marking/scanning or dense two-root interleaving kernel; consider fused root-switch/marking without changing physical coverage. | Exact saturation, repeated roots, tails and mutation/reservation boundaries. Preserve prior unroll/presieve losses. Native SIMD and byte constants are not Python defaults. |
| High preparation / A8 | Freeze genuine matrix/operator/lifting fixtures and packed forward/transpose products. Include current relation/preparation ownership and potential dense fill. | Independent dense-oracle products and original-relation kernel/extraction checks. Refresh representative fixtures after collector changes; synthetic capacity tests alone establish no factoring win. |
| High / B6 and C4 | Add a bounded Python-int Four-Russians/PLE challenger and fill-aware higher-weight filtering with compact reverse substitution. | Measure preparation, tables, fill, solver, lifting and extraction together. Stronger filtering need not wait for Four Russians to win. Native panel size and digit cutoffs require local selection. |
| High preparation / D1 | Design and implement the selected genuine block-Lanczos recurrence against A8's contract, with packed products and bounded restarts. Consider reduced-vector storage during design. | Singular-block and original-operator correction tests, checked resume and all live storage. B6/C4 are comparison controls, not required successes. Capacity can motivate sparse work before solve CPU dominates. |
| High design / B7 with C9a | Compact relation/exponent/provenance storage and replay; preflight or stream large checkpoint serialization. | Prove aliasing, pinned ownership, replay and interrupted scratch; verify original identities and corrupt-state rejection. Lower reservation coefficients or encoded-size quotas alone are insufficient. |
| High / C9a with C10a | Useful-row/solve pacing, then recovered-product intake, unequal LP endpoints and shorter cycle representatives. | Count retained useful relations and full cost, including rejected positions, repeated solves and failures. Preserve all-component cycles and certainty. |
| Investigate / F6 and C9a/F7 | General higher-arity LP incidence; size-aware residual splitting and residual batch cofactoring. | Explicit residual-domain, certification, provenance and capacity designs. Do not represent 3LP relations with a two-endpoint graph or assume native splitting cutoffs transfer. |
| Retain alternatives / A1a, C8, C4, D2, F4 | Cross-B prime-power reuse, eligible hypercube scheduling, hybrid matrix cores, joint exponent products, Block Wiedemann and bounded process workers. | Each has a distinct mathematical or ownership delta; use affected workload evidence and finite designs. Do not duplicate existing per-polynomial plans/queues or require unrelated successes. |
| High design / C11 and P6.4 | Current Python competitor and complete-cost comparison protocol. | Separate PyPy native-int, optional GMP, native wrappers and multiworker tracks; include cold/warm setup, CPU, RSS, completion and censoring. Independent moduli and fresh confirmation are required for a fastest/speedup claim. |

C10a-01 retains its implemented normalized form, wider exact multiplier ranking,
per-call residue reuse, paired A completion and guaranteed-two-power scoring.
Additional arithmetic/root changes belong in a newly identified coupled arm.
Batch, marking and matrix components may join a correct, feasible bundle without
first proving an isolated speedup; use informative ablations rather than a
Cartesian sweep. Do not make all roadmap implementation a prerequisite for the
existing four-size feasibility handoff, or require 65-digit completion before
an independently valid 70-digit comparison.

### Model and effort for the new QS work

For these recommendations, **Astra** means `gpt-6-astra` and **Sol** means
`gpt-6.1-sol`. The source review checked official [model-selection guidance](https://developers.openai.com/api/docs/guides/model-selection)
and the [Astra](https://developers.openai.com/api/docs/models/gpt-6-astra) and
[Sol](https://developers.openai.com/api/docs/models/gpt-6.1-sol) pages
on 2 November 2026. Both support high, xhigh and max. The assignments below are
engineering judgments about these tasks, not results of a Factor model benchmark.

**The C10a-01 workstream uses Astra / xhigh throughout**, including code,
correctness QA, calibration design and any authorized execution. Other or later
workstream recommendations retain their own scope. These assignments authorize
no delegation or experiment.

| Workstream | Design, implementation and correctness | Frozen execution and analysis | Reason / escalation boundary |
| --- | --- | --- | --- |
| C10a-01: Gray cache, attainable A population, actual-norm ranking and table setup | **Astra / xhigh** throughout | **Astra / xhigh** for this task | Coupled coordinates, arithmetic, finite population, cache ownership and resume; retain the workstream assignment. |
| A1c / A1a: bulk marking, paired roots and fused switch/marking | **Sol / xhigh** | **Sol / high** | Bounded kernel changes with existing independent oracles. Use Astra / xhigh for a specifically unresolved coverage or mutation/reservation proof. |
| A1a: cross-B lifted roots; C8: hypercube/CRT scheduling | **Astra / xhigh** for new arithmetic and its implementation | **Sol / high** after the protocol is frozen | Characteristic two, singular lifts, coordinate eligibility and table/queue identity create new mathematical obligations. |
| F7 / C9a: SIQS direct/tree smoothness and residual batch domains | **Astra / xhigh** | **Sol / high** | Wider-base construction, complete valuations, batch ownership, pending transitions and exact LP admission are coupled. |
| A8: matrix/operator/lifting contract and independent fixtures | **Astra / xhigh** | **Sol / high** for frozen fixture measurements | Orientation, original-kernel correction and lifetime accounting establish the validity of every solver comparison. |
| B6 / C4: M4RI/PLE and structured filtering | **Astra / xhigh** | **Sol / high** | Rank-deficient panels, fill, compact substitution and original-relation lifting need coordinated proofs and implementation. |
| D1 / D2: Lanczos or Wiedemann | **Astra / xhigh** | **Sol / high** | Singular-block recurrence or block polynomial generators, original-operator recovery, finite restarts and checkpoints. Reserve max for a named unresolved proof, not routine long runs. |
| B7: compact records, provenance and checkpoint preflight | **Astra / xhigh** | **Sol / high** | Aliasing, pinned atoms, reconstruction, corruption handling and simultaneous scratch determine correctness and capacity. |
| C9a / C10a: useful-row pacing and existing-domain policy | **Sol / xhigh** | **Sol / high** | Existing exact relation machinery supports bounded policy changes. New residual-certification, graph or provenance semantics move that subtask to Astra / xhigh. |
| F6: higher-arity LP incidence; C9a: new splitter/certification domains | **Astra / xhigh** | **Sol / high** | General incidence is not the existing DLP edge graph; certainty, repeated factors and lifting require new contracts. |
| C4: joint exponent products and hybrid recovery | **Astra / xhigh** for changed arithmetic/lifting; **Sol / xhigh** for an already specified kernel | **Sol / high** | Preserve sign, nonunit handling and original square identities while changing computation or representation. |
| B12 / F4: bounded process workers | **Sol / xhigh** | **Sol / high** | Lease accounting, cancellation, in-flight ownership and aggregate CPU/RSS dominate the implementation risk. Escalate an unresolved ownership invariant to Astra / xhigh. |
| C11 / P6.4: competitor protocol, configuration and budget selection | **Sol / xhigh** | **Sol / high** | Freeze comparable inputs/backends, censoring, warmup, work accounts and uncertainty. Use Astra / xhigh when interpreting new arithmetic or conflicting causal evidence. |

For future tasks, escalation to Astra / max requires a concrete failed invariant,
counterexample or unresolved derivation. More elapsed time, more source files or
a longer experiment alone is not a reason to increase reasoning effort. Model
recommendations do not replace independent oracles or matched performance evidence.

### C10a-01 implementation remaining after the coverage review

Already implemented and independently checked: parity2 normalization with
standard fallback, wider exact multiplier ranking, per-call residue reuse,
paired-flyer A completion and guaranteed-two-power initialization. The last
full arithmetic suite passed: 859 tests, with nine optional-backend skips.
The revised combined bundle still has no timed calibration or complete-factor
result; the earlier prototype calibration is a different source identity.

The next concrete implementations are:

1. **Centered modular Gray-shift cache.** Preserve the existing B and U trials;
   include wrap correction, exceptional roots, bounded consumed-coordinate
   tables, policy identity, construction/refusal charges and checked resume.
2. **Better attainable A populations.** Specify and implement a finite
   feasible-core/product-window policy beyond the current paired completion.
   Preserve uniqueness/exhaustion and test actual reachable products against an
   independent Cartesian oracle. The current closest pair need not be close
   to the target for every core.

Two additional alternatives belong to this workstream but still need concrete
selection/specification: **actual-A/norm-aware multiplier ranking**, and
**small-factor/batched-inverse root-table construction**. These are not already
provided by wider local-density scoring or cached A inverses. Choose exact,
bounded rules before implementing them; no prior standalone speedup is required
for a correct, feasible opt-in experiment.

Each selected change also needs integration with independent CRT/root/valuation
checks, caches/queues/plans, compatibility/resume, ownership and work bounds.
Then refresh the literal challenger, short-calibration protocol and later-stage
source/work freezes. The 60/62/65/70 complete-factor pilots remain later work.
Bernstein SIQS integration and new matrix solvers have explicit separate roadmap
items; their completion is not a prerequisite for this arithmetic handoff.

## Open and conditional execution items

Priority and predecessor order below are retained from the existing plan.
A predecessor settles through an accepted result or recorded retain/defer/reject
decision, subject to the [policy](POLICY.md#dependencies-and-priority).
Each live row links its complete original gates, quantitative context and
research. The [history index](history/README.md) supplies local research paths.

### Execution phase A — independent foundations that can proceed now

| ID / scope | Deliverable and boundary | Predecessors | Status, trigger and existing priority | Model / effort | Evidence and research |
| --- | --- | --- | --- | --- | --- |
| A5 — P5.1 Williams p+1 binary baseline | Implement exact binary Lucas stages 1 and 2, bounded parameter trials, discriminant checks, saturation recovery and checkpoints. | Existing P2 bounded recovery | **Open.** Priority / trigger: Medium, complementary coverage. This is independent of SIQS and PRAC. Binary Lucas supplies the correctness control required before optimized Lucas chains or special starts are ranked. | Sol / **xhigh** | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-a--independent-foundations-that-can-proceed-now) |
| A8 — P3.8 matrix control + remaining R3 diagnosis | Freeze exact matrix/operator/lifting interfaces and genuine post-filter fixtures; measure solve/filter/extraction cost, fill, provenance peaks and capacity refusals; compare packed sparse products using reviewed implementations references (AP16). | Accepted R3 identities/lifting and P3.3 control | **Open.** High preparation under the [QS coverage review](#qs-optimization-priorities--coverage-reconciliation-2-november-2026). Freeze representative fixtures before comparisons. Refresh after C9/C10 change the matrices. Establish whether solving or memory blocks useful factoring; synthetic capacity cases cannot establish a current speedup. A verified dense-memory refusal can justify B7/D1/D2 preparation without dominant solve CPU. | Astra / **xhigh** | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-a--independent-foundations-that-can-proceed-now) |

### Execution phase B — exploit the settled controls

| ID / scope | Deliverable and boundary | Predecessors | Status, trigger and existing priority | Model / effort | Evidence and research |
| --- | --- | --- | --- | --- | --- |
| B5 — P3.7 optional NumPy spike | Test bounded vectorized score/root-hit updates and candidate extraction against PyPy bytearray/array/list controls. | A1 plus a measured remaining array bottleneck | **Conditional.** Priority / trigger: Conditional. Verify PyPy availability, overflow bounds, duplicate-hit accumulation and full-run conversion/import costs; otherwise defer. It is not required for Phase 3 exit. | Sol / **xhigh** | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-b--exploit-the-settled-controls) |
| B6 — P3.8 Four Russians | First bounded dense/hybrid challenger against the bitset control; account for tables, conversion, recovery, lifting and peak simultaneous storage. | A8 plus a representative solve-cost/memory case | **Open.** High challenger preparation under the [QS coverage review](#qs-optimization-priorities--coverage-reconciliation-2-november-2026). First dense-kernel comparison; use A8 fixtures. Existing identities and real matrix controls are the prerequisites. It can overlap independent collector and ECM work; it does not wait for R4, R5 or NumPy. | Astra / **xhigh** | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-b--exploit-the-settled-controls) |
| B7 — remaining P3.8-R3 capacity/provenance | Prove bounded compact partial/relation/provenance layouts, replayable records, merge histories or normalized square-root payloads against reviewed implementations (AP02/AP15); independently verify original identities and lifting. | A8 for matrix/provenance work, or C9a for a scoped partial-store implementation; demonstrated capacity/provenance bottleneck | **Open.** High design under the [QS coverage review](#qs-optimization-priorities--coverage-reconciliation-2-november-2026). A compact C9a partial-store challenger need not wait for a solver. Charge replay, pinned ownership and simultaneous scratch; reducing a reservation coefficient alone is not a storage implementation. Keep original-row verification and nonunit recovery. | Astra / **xhigh** | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-b--exploit-the-settled-controls) |
| B8 — Williams p+1 starts (P5.1) / Lucas optimization (P5.3) | Compare bounded rational/seeded p+1 starts and validated Lucas-chain execution against the binary control; validate p+1 bound extension separately from ordinary powering. | A5, A4; A6's schedule-ratio contract | **Conditional.** Priority / trigger: Medium/conditional. Binary correctness and verified chain machinery must exist first. This targets marginal p+1 coverage; the ECM recurrence is not a drop-in Lucas implementation. | Astra / **xhigh** | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-b--exploit-the-settled-controls) |
| B12 — P6.3 worker-contract preparation | Specify reusable assignment IDs, immutable schedules, parent-owned leases, aggregate CPU/RSS, cancellation and restart contracts; reconcile existing QS workers. Prepare bounded fixtures, not a new worker default. | A2, A3, A7 | **Open.** Priority / trigger: Medium, early preparation. Contracts can proceed once backend, program and repair interfaces settle. Production portfolio experiments wait for the stable serial tranche in F4; existing QS workers are reused. | Sol / **xhigh** | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-b--exploit-the-settled-controls) |
| B13 — P3.8-R2 wide-block resieve capacity | Reuse the accepted sparse resieve support bound. Only if a new calibrated workload still refuses setup, prove a bound for its actual scratch/storage representation, verify coverage and refusal/resume, and compare complete bucket/resieve runs under matched budgets. | A1; a new calibrated wide-block capacity refusal | **Conditional.** Priority / trigger: Deferred unless newly triggered. The combined 30-digit bridge resolves the isolated refusal with the accepted repair. A larger workload needs new evidence and a proof; unrelated feasible collection work proceeds independently. | Astra / **xhigh** | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-b--exploit-the-settled-controls) |
| B14 — P8.2.2 checked prime certificates | Add bounded n−1/Pocklington certificate generation and a separate exact verifier beyond deterministic ranges, with recursive work/storage caps and requested-proof CLI/resume. Record historical APRT-CLE (AP28) as a later coverage candidate only if this first tranche leaves a measured gap; distinguish its probabilistic fallback. | A10; existing shared-budget/checkpoint controls | **Open.** Priority / trigger: High for certified results. Proof is an explicit capability beyond v1 parity. Start after the range contract settles, without waiting for GNFS; certify supported cases rather than claiming arbitrary-size proof from extra probable-prime rounds. | Astra / **xhigh** | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-b--exploit-the-settled-controls) |

### Execution phase C — measured current-engine challengers

| ID / scope | Deliverable and boundary | Predecessors | Status, trigger and existing priority | Model / effort | Evidence and research |
| --- | --- | --- | --- | --- | --- |
| C2 — P5.2 advanced pairing/common-Z | Separately compare extended distances, relocation/matching and coverage pruning with reference ECM/other pinned references; compare projective terms, the calculator reference-style affine differences, derived batch inversion and no-inversion common-Z tables (AP29). | B2 aligned-cell decision plus remaining schedule/product cost | **Conditional.** Priority / trigger: Conditional. Coordinate geometry and coordinate choices without conflating them. Bound planning, tables, inversions and replay; independently certify prime coverage and nonunit recovery. Prior B2 losses remain the control; polynomial continuation stays F3. | Astra / **xhigh** | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-c--measured-current-engine-challengers) |
| C4 — P3.8 other dense/hybrid/filtering | Compare PLE/free-variable recovery, sparse-to-dense cores, components and stronger bounded filtering; integrate any accepted provenance representation. | A8; B7 decision for changed representations | **Open.** **High preparation.** Four Russians remains the first dense-kernel comparison; stronger filtering may proceed independently. Reuse its control before adding more interacting transformations; charge fill, retained history and recovery, not just matrix dimension. | Astra / **xhigh** | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-c--measured-current-engine-challengers) |
| C5 — P4.4 reducers | Revisit persistent Barrett/Montgomery contexts only in actual fused engine loops, with exact encoded identities, width bounds and canonical exits. | A2, B4 plus a remaining reduction bottleneck | **Conditional.** Priority / trigger: Low/conditional. Earlier reducers lost near 166–200 bits. Backend/kernel results must provide a reason to reopen them; native `%` remains the default if the whole-run gate fails. | Astra / **xhigh** | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-c--measured-current-engine-challengers) |
| C8 — P3.8-R2 family-wide CRT hit scheduling | Use calibrated whole-polynomial-interval prime/power eligibility and root costs to decide whether bounded CRT half-sum scheduling is warranted; if justified, verify every Gray/position hit, exceptional roots and table/queue limits, then compare complete factoring. | A1, B1; demonstrated family/root or eligible-prime scanning cost | **Conditional.** Priority / trigger: Low/conditional. The isolated first-polynomial probe found zero eligible base primes and sparse eligible power hits. Retain the deferral until a calibrated workload supports the investment; block width is not the eligibility interval. | Astra / **xhigh** | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-c--measured-current-engine-challengers) |
| C10a — conditional collection-policy follow-ups | Remaining collection-policy questions are C10a-01/-02/-03/-06 and explicit larger-size feasibility below. C10a-04/-05/-07 are Done for delivered scopes; retain their unmet gates and scoped decisions. | Closed C10 tranche evidence; coordinate A1a and C9a boundaries | **Open.** Each subtranche requires its own stage-cost/useful-yield/exact-loss trigger. Fixed 16-family quality ordering was Rejected for promotion; broader product selection remains Conditional. | Astra / **xhigh** for representation; Sol / **high** for frozen policies | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-c--measured-current-engine-challengers) |
| C11 — Cross-engine constant provenance and continuing calibration | Maintain the cross-engine parameter ledger against pinned reviewed implementations and other primary implementations. Route allocation deltas to C3a, classification deltas to A10a, inverse microkernels to A2a and GCD/replay cadence to E5. Record AP34/AP35 native fixed-limb/FFT/NTT dispositions under existing Python/PyPy scope and prior losses; C9a/C10a retain new QS intake/collection ownership. | Reuse accepted A2/A6/A10/B3/B4/C6 controls; coordinate E3/E4/C3/E6/E1/G1; C9a/C10a retain remaining QS/DLP ownership | **Open.** Priority / trigger: High, continuing. Inherited bounds, abstract chain weights, batch sizes, duplicated standalone/portfolio settings and conservative memory coefficients are not proven optima. Begin with observed cost/yield gaps and missing provenance; retain defaults on inconclusive evidence. No repeated completed chain search or blanket Cartesian sweep. | Sol / **high** for source reconciliation and frozen configuration studies; Astra / **xhigh** for changed arithmetic, certainty or recovery proofs | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-c--measured-current-engine-challengers) |

#### New source-derived deltas — separate follow-up items

| ID / scope | Deliverable and boundary | Predecessors | Status, trigger and existing priority | Model / effort | Evidence and research |
| --- | --- | --- | --- | --- | --- |
| C3a — source-derived ECM policy follow-up | After C3, reconcile reviewed implementations effort, curve/bound and pretest/campaign policies (AP04) with its final selection. Compare only missing source-supported alternatives, prior-work credit or handoff behavior; retain C3's frozen evidence. | C3 decision; a documented uncovered policy or changed-cost hypothesis | **Conditional.** Priority / trigger: Conditional follow-up. Existing C3 coverage can resolve candidates without another study; any new selection needs fresh complete-factor evidence. | Sol / **xhigh** | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#new-source-derived-deltas--separate-follow-up-items) |
| C9a — remaining large-prime yield, target feedback and compact storage | Attribute the literature-verification gaps: yield-dependent retry targets, wider recovered-product intake, asymmetric ordered endpoints, shorter cycles/algebraic duplicates and compact replayable partial records. Coordinate scheduling/scoring with C10a and storage/provenance proofs with B7; C1a owns anchor-only SLP. | Completed C9.2 control and transfer audit; reconcile the post-C10a-07 collector; measured affected-stage cost or useful-yield loss | **Open.** Priority / trigger: Open; ranked follow-up. Start with target feedback/pacing attribution, then residual domains and cycle weight. Splitter/batch/disk work remains profile-gated. Freeze each selected tranche and confirm fresh complete factors under matched total resources; preserve C9 acceptance. | Sol / **xhigh** policy design and implementation; Astra / **xhigh** for changed residual/provenance proofs; Sol / **high** frozen execution/reporting | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#new-source-derived-deltas--separate-follow-up-items) |
| A1a follow-ups — quota-aware setup and lifted-root reuse | Bound setup by consumed quota, reuse validated switch metadata, and separately investigate prime-power roots across B transitions. Existing base-prime Gray reuse and per-polynomial plans are controls. | Delivered A1a remains closed; current evidence of material setup/lifting cost after the C10a-07 queue decision | **Conditional.** Priority / trigger: Conditional backlog. Separate the two hypotheses; preserve larger-cache losses and existing C8 half-sums. No reopening of delivered A1a or current face/queue freeze. | Sol / **xhigh** design/implementation; Astra / **xhigh** only for an unresolved lifted-root/exceptional-case proof; Sol / **high** frozen execution/reporting | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#new-source-derived-deltas--separate-follow-up-items) |
| A1b-F1 — fixed-score sparse-queue integration | A separate mechanism: policy-bound fixed-score queue weights/eligibility, exact threshold/power coverage, scratch refusal and charged rebuild. Do not modify historical A1b controls. | Completed A1b existing-option decision; a new credible cost hypothesis and separately frozen design | **Deferred.** The implemented fixed/list bundle lost all 50/55/60 screens. No queue extension or additional campaign is admitted by that result. | Sol / **xhigh** integration/ownership; **high** frozen evaluation | [A1b decision and boundary](benchmarks/RESULTS.md#a1b-existing-fixed-score-review) |
| A1c — same-modulus two-root interleaving | Assess one shared stride loop for two distinct roots of a prime or power, including sorted translated offsets, tails and scalar exceptional-root fallback. Specify reservation/cancellation boundaries and charged interrupted-block replay before implementation. | Settled A1a kernels and demonstrated remaining dense-marking cost | **Conditional.** Priority / trigger: Conditional follow-up. Source-supported but deferred from the current two-kernel extension: interleaving changes the existing root-by-root mutation/reservation boundary. Preserve low-allowance behavior and compare complete factoring against the settled A1a control. | Sol / **xhigh** design and implementation; **high** frozen experiments; Astra / **xhigh** only for an unresolved atomicity proof | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#new-source-derived-deltas--separate-follow-up-items) |
| A2a — inverse-microkernel follow-up to completed backend study | Compare batched binary/Montgomery inversion and small-integer quotient shortcuts (AP32/AP33) only for newly profiled inversion-dominant sizes, against native Euclid/GMP including conversions and nonunits. | Completed A2/B4; C11 profile and a changed workload/cost hypothesis | **Conditional.** Priority / trigger: Low/conditional. Preserve backend/reducer losses. Native fixed-limb/SIMD/FFT/NTT references (AP34/AP35) remain deferred C11/C5 dispositions, not a new Python limb engine. | Sol / **xhigh**; targeted Astra proof review if required | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#new-source-derived-deltas--separate-follow-up-items) |
| A7a — SIQS worker-granularity follow-up | Compare historical SIQS family-sharing/partitioning with existing bounded QS workers (AP37), only after measured setup/coordination cost motivates the alternative. F4 integrates any selected worker result. | Completed A7; B12 contracts; stable serial control and F4 experiment gate | **Conditional.** Priority / trigger: Conditional. No generic worker reimplementation or changed A7 closure. Charge shared stores, barriers, aggregate work/RSS, cancellation and restart; Java thread behavior is not PyPy evidence. | Sol / **xhigh** | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#new-source-derived-deltas--separate-follow-up-items) |
| A10a — classification follow-up to completed primality transfer | Evaluate MR-residue proper-factor recovery, optional BPSW composite rejection and grouped small-prime remainders (AP25–AP27) only where E3/C11 profiles identify an uncovered benefit. | Completed A10; E3 frozen control or a justified isolated oracle spike; classification-cost evidence | **Conditional.** Priority / trigger: Conditional. Preserve deterministic ranges, requested MR rounds, RNG/work identity and certainty. A BPSW pass is not proof; prior primorial losses remain the control. B14 retains certificate ownership. | Sol / **xhigh** | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#new-source-derived-deltas--separate-follow-up-items) |
| A11a — optional input/browser follow-up | Separately select validated factor hints/opt-in lookup, bounded batch/expression input or browser progress/offline/save-continue (AP38–AP40). Give each selected subtranche its own API and acceptance boundary. | Completed A11; explicit product selection, outside the optimization critical path | **Conditional.** Priority / trigger: Optional; no factoring speed claim. Preserve quiet library APIs, exact reconstruction, certainty, finite input/output/network allowances and checked checkpoints. | Sol / **high**; **xhigh** for changed accounting/resume | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#new-source-derived-deltas--separate-follow-up-items) |
| C1a — cheaper SLP anchor reuse follow-up | Compare an anchor-only SLP implementation or changed retention regime (AP01) with ordinary SLP and the existing graph-SLP/DLP controls. Prove residual cancellation and bounded pinned ownership without adopting source dependency-discard policies. | Completed C1; C9 decision; a new cost/retention hypothesis | **Conditional.** Priority / trigger: Low/conditional. C1 graph-SLP training totaled 20.055 s versus 20.160 s for SLP, with no established win; DLP already retains anchors. Do not repeat that arm or promote on raw cancellation counts. | Astra / **xhigh** | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#new-source-derived-deltas--separate-follow-up-items) |

### Execution phase D — sparse-solver alternatives

| ID / scope | Deliverable and boundary | Predecessors | Status, trigger and existing priority | Model / effort | Evidence and research |
| --- | --- | --- | --- | --- | --- |
| D1 — P3.8 block Lanczos | Use pinned official reviewed implementations to implement bounded packed-product Block Lanczos (AP16/AP17): seeded recurrence, selected nonsingular blocks, original-operator kernel correction, lifting, finite retries and checked resume. | A8 operator/lifting contract; B6/C4 are comparison controls, not required successes | **Open.** **High preparation; first sparse candidate by scheduling preference.** Compare total solve/recovery/storage with feasible bitset/dense/hybrid controls. Capacity can justify the work before solve CPU dominates. Every candidate must satisfy original `M d = 0`; no solver fixes absent useful relations. | Astra / **xhigh**; **max** for a specific unresolved proof | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-d--sparse-solver-alternatives) |
| D2 — P3.8 block Wiedemann | Use pinned NFS reference for genuine Block Wiedemann: projected Krylov sequences, a block polynomial generator, reconstruction, original-kernel verification, finite I/O/checkpoints and independent-sequence accounting. Reuse A8 packed-product oracles; the calculator reference does not implement this solver. | A8 operator/lifting contract; B6/C4 are comparison controls, not required successes | **Conditional alternative to D1; no D1 dependency.** Compare generator/reconstruction/I/O and failure costs, retained state and complete factors. Prefer Lanczos first for the current single-machine tranche unless measured capacity/sequence/restart needs favor Wiedemann; start with an exact base-case generator. | Astra / **xhigh**; **max** for a specific unresolved proof | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-d--sparse-solver-alternatives) |

### Execution phase E — current-engine acceptance, control freeze and P8 experiments

| ID / scope | Deliverable and boundary | Predecessors | Status, trigger and existing priority | Model / effort | Evidence and research |
| --- | --- | --- | --- | --- | --- |
| E1 — final R1/R5 and P4/P5 acceptance | Integrate selected changes, including accepted A10/B14 capabilities and A11 CLI integration; reconcile v1 transfer decisions and recalibrate affected SIQS/ECM/p±1 parameters; freeze selections; run fresh complete-factor comparisons, R5 SSS/worker challengers, resume checks and clean-checkout validation; publish adopt/defer/reject decisions. | A1, A7, A10, A11 and chosen P3.7/P3.8/P4/P5 predecessor decisions; B14 if proof support is selected; B13/C8 and other unchosen branches have recorded deferrals | Deferred; resumption requires explicit authorization. Original combined confirmation and all unrun gates remain required for experimental acceptance. | Sol / **xhigh**; targeted Astra / **xhigh** review for new arithmetic/provenance | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-e--current-engine-acceptance-control-freeze-and-p8-experiments) |
| D7 — P8.1 fresh portfolio control | Freeze the accepted current SIQS/MPQS, rho, ECM and p±1 portfolio with certified training/untouched confirmation inputs; profile separately and retain source/configuration identity for every later challenger. | E1; existing SIQS/P2 controls | **Open.** Priority / trigger: High enabling value, bounded scope. P8 tuning uses the accepted current engines before GNFS. Keep this control immutable; future GNFS requires a separately versioned extension and crossover evaluation. | Sol / **high**; **xhigh** for protocol decisions | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-e--current-engine-acceptance-control-freeze-and-p8-experiments) |
| E3 — P8.2 preprocessing | Compare bounded multiplier Fermat/Lehman, shared power screens, inherited trial coverage and high-valuation extraction (AP18/AP22–AP24). Separately evaluate opt-in algebraic a^k±1, sums/differences of powers and historical Fibonacci/Lucas families (AP19–AP21); A10a owns AP25–AP27 classification implementation. Reuse accepted A10/B14 certainty controls. | D7 | **Conditional.** Priority / trigger: Medium; profile-gated. Broad cheap-path savings may help, but prior repair decisions are retained. Full integrated evaluation now has the right control; a justified isolated oracle spike could occur earlier. | Sol / **xhigh** | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-e--current-engine-acceptance-control-freeze-and-p8-experiments) |
| E4 — P8.3 rho calibration | Tune bounded batch/walk/restart policies under identical total budgets and assigned seeds; measure first-factor and complete runs including saturation recovery. | D7 | **Open.** Priority / trigger: Medium. Parameter tuning is relatively contained, but must show marginal portfolio value on the new workload. It reuses current Brent/local-loop controls rather than rebuilding rho. | Sol / **high**; **xhigh** for changed recovery semantics | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-e--current-engine-acceptance-control-freeze-and-p8-experiments) |
| E5 — P8.5 recovery/polling/checkpoint costs | Separate polling, atomic commits, saturation recovery, durable writes and JSON verification; compare source-derived ECM GCD cadence/replay (AP31) and bounded recovery trees only where profiles justify them. Preserve strict local factor recovery, deadlines and charged resume. | D7 | **Conditional.** Priority / trigger: Medium; conditional on changed cost. Existing repair gains remain the starting point. Revisit only unresolved or newly dominant overhead, preserving validation and disclosing work-unit changes. | Sol / **xhigh** | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-e--current-engine-acceptance-control-freeze-and-p8-experiments) |
| E6 — P8.6 contexts/schedules for changed workloads | Reuse C10-diagnosed/A1a-implemented periodic presieve, paired-polynomial, bulk-mask, density kernels and Gray-shift tables (AP05–AP09); avoid duplicate implementation. Test changed-workload shared contexts, lazy setup/private buffers and hit prefilters against accepted exact controls and prior cache/wheel losses. | D7 | **Conditional.** Priority / trigger: Conditional. This is workload-specific tuning after control freeze, not another implementation of P5.2 pairing. Retain prior cache/wheel/rolling and small-prime-skipping decisions unless new evidence overturns them. | Sol / **high**; **xhigh** for ownership/rejection changes | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-e--current-engine-acceptance-control-freeze-and-p8-experiments) |

### Execution phase F — independent P6 challengers

Use a current committed control for independent F6/F7 feasibility; deferred
E1 is not their prerequisite. Other P6 rows retain their stated predecessors.
Choose challengers with concrete cost/yield hypotheses and finite designs.

| ID / scope | Deliverable and boundary | Predecessors | Status, trigger and existing priority | Model / effort | Evidence and research |
| --- | --- | --- | --- | --- | --- |
| F2 — P6.1 Edwards/windowed/torsion-aware ECM | Compare complete curve-family/stage-one packages from reference ECM and other pinned references with accepted Suyama, including the calculator reference different rational Montgomery family (AP30), Edwards/windowed variants, validated stage-two conversion and setup costs. | E1; A2/B3/B4 accepted contracts or retained controls | **Conditional.** Priority / trigger: Conditional. First establish how far ordinary backend/chain/kernel improvements go. Curve-order torsion, chosen-point order and exceptional maps need whole-engine evidence, not operation counts. | Astra / **xhigh** | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-f--independent-p6-challengers) |
| F3 — P6.2 polynomial ECM continuation | Choose a justified reference ECM/primary-source Brent–Suyama, product/remainder-tree or multipoint continuation challenger. Assess required exact polynomial kernels, cached reciprocal reduction and coefficient/storage bounds using AP36 references; the calculator reference polynomial calculator is not an implemented polynomial ECM engine. | E1, B2; C2 decision and a remaining continuation bottleneck | **Conditional.** Priority / trigger: Conditional, potentially important at larger targets. Paired classical stage two is the baseline. Advance when per-prime continuation remains the limiting algorithm, not merely because larger bounds are available. | Astra / **xhigh** | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-f--independent-p6-challengers) |
| F4 — P6.3 production ECM/portfolio workers | Reuse bounded parent-owned leases, stable assignments, cancellation/restart and settled QS worker decisions; consume any selected A7a family-granularity result (AP37). Compare 1/2/4 production workers for fixed-work and first-factor cases with total CPU and aggregate RSS. | E1, B12; E2 if GNFS jobs are included | **Conditional.** Priority / trigger: Conditional measured parallelism. Stable serial engines and worker contracts precede credible comparisons. Reuse P3.6/P3.6.1 results; extra cores and earlier thread improvements do not establish a first-factor win. | Sol / **xhigh** | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-f--independent-p6-challengers) |
| F5 — P6.2 cross-family A=A0*q polynomial reuse | Compare bounded reuse across QS families with calibrated factor-base-smooth A, explicitly carrying any external q exponent, partial-relation role and duplicate policy. | E1; B1 calibrated QS control | **Conditional.** Priority / trigger: Conditional on polynomial/root setup cost. This changes relation semantics and is not a drop-in root cache. Ordinary families and the R1 external-square distinction provide the control, not a correctness shortcut. | Astra / **xhigh** | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-f--independent-p6-challengers) |
| F6 — P6.2 triple-large-prime/general QS relations | If DLP economics justify another extension, implement general sparse incidence/provenance and bounded residual splitting; compare useful dependencies and complete factors with the calibrated DLP control. | C1/current DLP control and an explicit higher-arity domain/capacity design; E1 is not required for isolated feasibility | **Low/conditional.** A third large prime increases splitting, storage and verification complexity. Do not reuse a two-endpoint edge-cycle algorithm as though it represented every higher-arity relation. | Astra / **xhigh** | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-f--independent-p6-challengers) |
| F7 — P6.2 batch smooth-part/remainder-tree revisit | Compare bounded SIQS direct/tree smooth-part screening on the changed larger-basis workload; retain scalar exact exponent recovery and include full pipeline latency/storage. | Current collector control, SIQS-domain ownership proof and a changed-workload cost hypothesis; E1 is not required | **Open.** **High design priority.** R2 already evaluated a narrower batch-recovery prototype. This is a changed-workload revisit, not duplicate work or an automatic reversal of a loss. Sparse solvers and higher merges remain owned by the existing matrix rows. | Astra / **xhigh** | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-f--independent-p6-challengers) |
| F8 — P6.2 specialized fast P−1/P+1 continuation | Investigate Montgomery–Kruppa reciprocal-Laurent-polynomial construction and geometric-progression evaluation through pinned reference ECM the implementation, separately from F3's ECM/Brent–Suyama continuation. | A6 control and measured remaining PM1 stage-two cost; E1 for integrated acceptance; A5 only for the P+1 arm | **Conditional.** Priority / trigger: Conditional, later algorithmic challenger. Potentially useful at larger continuation bounds, with exact polynomial/backend and storage costs established first. No C2 geometry dependency and no native speedup extrapolation. | Astra / **xhigh** | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-f--independent-p6-challengers) |

### Execution phase G — integrated portfolio allocation

| ID / scope | Deliverable and boundary | Predecessors | Status, trigger and existing priority | Model / effort | Evidence and research |
| --- | --- | --- | --- | --- | --- |
| G1 — P8.4 joint allocation and handoff | Fit bounded rho/p−1/p+1/ECM allocations and QS/MPQS/SIQS handoff using observable input and completed-work metadata. Consume accepted engine changes without reimplementing them. | D7, E1; E3–E6 and selected F2–F7 decisions | **Open.** Priority / trigger: High downstream value. Tune the settled current engines before GNFS. Known factor sizes stratify evaluation but stay hidden from the dispatcher. Deferred F1 later supplies any measured GNFS crossover update. | Sol / **xhigh** | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-g--integrated-portfolio-allocation) |

### Execution phase H — final Phase 8 reconciliation

| ID / scope | Deliverable and boundary | Predecessors | Status, trigger and existing priority | Model / effort | Evidence and research |
| --- | --- | --- | --- | --- | --- |
| H1 — P8.7 integrate and confirm selected winners | Integrate selected P8 and relevant P4–P6 winners, freeze the final configuration, confirm on untouched data and reconcile each historical P2 gate with adopt/retain/defer/reject evidence. | G1, E3, E4, E5, E6 decisions; selected P6 implementations | **Open.** Priority / trigger: Required for the selected final tranche. Component wins do not prove combined performance. Optional losses do not block other accepted changes, and completed P2 work is not reopened by renumbering. | Sol / **xhigh** | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-h--final-phase-8-reconciliation) |

### Execution phase I — reproducible publication

| ID / scope | Deliverable and boundary | Predecessors | Status, trigger and existing priority | Model / effort | Evidence and research |
| --- | --- | --- | --- | --- | --- |
| I1 — P6.4 final comparative publication | Publish source/configuration pins, commands, independent corpus/runner links, backend/core disclosures, supported limits and complete versus factor-one outcomes; run feasible pinned competitor arms. | H1; selected P6 decisions and clean-checkout validation; F1 only for later GNFS claims | **Open.** Priority / trigger: Required for comparative claims. Final publication follows the measured integrated result, although individual milestone evidence can be published earlier. P8.1 supplies the protocol; P6.4 owns presentation/comparison rather than duplicating it. | Sol / **high**; **xhigh** for disputed comparisons | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-i--reproducible-publication) |

### Later tranche — GNFS after current-engine optimization

**Deferred:** activate only after H1 and the coverage review.

| ID / scope | Deliverable and boundary | Predecessors | Status, trigger and existing priority | Model / effort | Evidence and research |
| --- | --- | --- | --- | --- | --- |
| A9 — P7.1 GNFS contracts + P7.4 root-method scope | Define polynomial/field, nonmonic norm, ideal/root and bad-prime identities; choose a finite supported algebraic-root strategy before expanding fields. Use exact small oracles. | H1; coverage review activates GNFS; existing P2 and P3.3/P3.4 controls | **Deferred.** Priority / trigger: Deferred, later coverage investment. Start only after current-engine optimization and fresh confirmation, if the coverage review justifies GNFS. Root-method constraints must inform field selection early. | Astra / **xhigh** | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#later-tranche--gnfs-after-current-engine-optimization) |
| B9 — P7.2 full-relation collector | Build bounded serial line sieving with exact two-norm/ideal verification, primitive-pair deduplication, collection cursors and bounded storage. Start with full relations. | A9 | **Deferred.** Priority / trigger: Deferred with GNFS. High for the GNFS reference. Stable field/ideal identities must precede stored relations. Full relations provide a tractable control before partials and special-q scaling. | Astra / **xhigh** | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#later-tranche--gnfs-after-current-engine-optimization) |
| B10 — P7.3 GNFS matrix/character control | Implement exact bitset constraints and original-relation lifting on independently generated fixtures, with explicit sign/character placement and kernel correction. | A9; existing P3.3 solver control | **Deferred.** Priority / trigger: Deferred with GNFS. High for correctness. Fixture-based solver development can run alongside collection. It does not require Four Russians, Lanczos or Wiedemann; ideal parity alone is insufficient. | Astra / **xhigh** | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#later-tranche--gnfs-after-current-engine-optimization) |
| B11 — P7.4 algebraic/rational root implementation | Implement the selected bounded root method against known-square field fixtures, including coefficient/precision bounds, signs, denominators and modular mapping. Defer full pipeline acceptance to C7. | A9 | **Deferred.** Priority / trigger: Deferred with GNFS. High and proof-sensitive. Root arithmetic can develop alongside B9/B10 using independent fixtures. An unsupported field or failed auxiliary-prime search must produce finite refusal, not an invalid root. | Astra / **xhigh**; **max** for a named unresolved proof | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#later-tranche--gnfs-after-current-engine-optimization) |
| C7 — P7.3–P7.5 small GNFS integration | Connect collected relations, constrained dependencies and both square roots; validate complete small general composites, corrupt-state checks, cumulative allowances and explicit opt-in dispatch. | B9, B10, B11; existing working SIQS baseline | **Deferred.** Priority / trigger: Deferred with GNFS. High correctness milestone once activated. This is the first whole-pipeline GNFS correctness gate; a matrix success or integer norm root cannot substitute for a proper divisor. | Astra / **xhigh** | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#later-tranche--gnfs-after-current-engine-optimization) |
| D3 — P7.6 polynomial selection | Compare bounded degree/skew/root-quality search and trial sieving with the base-m reference; verify every common root and charge selection cost. | C7, A2 decision | **Deferred.** Priority / trigger: Deferred with GNFS. High scaling candidate after correctness. Better norms can affect the entire collector. Compare on a fixed reference collector so improvements are attributable; GMP is optional. | Astra / **xhigh** | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#later-tranche--gnfs-after-current-engine-optimization) |
| D4 — P7.6 prime special-q lattice sieve | Implement bounded prime-special-q assignments, lattice mappings, forced ideal exponents, buckets/spills and exact duplicate coverage on fixed accepted polynomials. | C7, A2 decision | **Deferred.** Priority / trigger: Deferred with GNFS. High potential, substantial implementation. The full-relation verifier is the control. This branch can run independently of D3 using frozen polynomials; composite special-q follows prime-q correctness within this track. | Astra / **xhigh** | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#later-tranche--gnfs-after-current-engine-optimization) |
| D5 — P7.6 two-sided cofactoring and large ideals | Evaluate bounded residual strategies and general side-labelled ideal incidence against the full-relation control; preserve repeated powers, certainty and distinct roots above the same prime. | C7, A2 decision | **Deferred.** Priority / trigger: Deferred with GNFS. Conditional on collection/cofactor cost. Use fixed polynomial and sieve fixtures for independent work. QS DLP code supplies implementation experience, not a valid general GNFS edge-graph contract. | Astra / **xhigh** | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#later-tranche--gnfs-after-current-engine-optimization) |
| D6 — P7.6 matrix/storage scaling | Profile genuine GNFS matrices, retain ideal/character constraints and lifting, and validate bounded layouts/spill. Reuse accepted P3.8 kernels; coordinate any new sparse solver with D1/D2. | C7, A2 decision | **Deferred.** Priority / trigger: Deferred with GNFS. Conditional on actual matrix pressure. GNFS matrix preparation can run alongside D3–D5 against frozen fixtures. No duplicate Lanczos/Wiedemann implementation or compulsory sparse solver is introduced. | Astra / **xhigh** | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#later-tranche--gnfs-after-current-engine-optimization) |
| E2 — P7.6 scaled GNFS integration | Integrate accepted selection/sieve/cofactor/matrix branches, retain rejected controls, validate full reconstruction and bounded restart, and freeze the scaling configuration before crossover measurement. | D3, D4, D5, D6 decisions; D1/D2 only if their solvers are selected | **Deferred.** Priority / trigger: Deferred with GNFS. High downstream milestone. Independent stage wins can interact or lose useful relations. This is the serial scaling gate; broad process parallelism is optional and must use F4 accounting if later adopted. | Astra / **xhigh** | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#later-tranche--gnfs-after-current-engine-optimization) |
| F1 — P7.7 SIQS/GNFS crossover | Freeze scaled GNFS and calibrated SIQS arms; compare fresh general-composite bands, complete/factor-one outputs, exhaustion, CPU/RSS/disk and cold/warm costs before any default handoff. | E2, E1 | **Deferred.** Priority / trigger: Deferred with GNFS. High decision value. A correct small GNFS engine is not evidence of a useful digit crossover. This needs both measured engines and remains independent of optional P6 research. | Sol / **xhigh** | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#later-tranche--gnfs-after-current-engine-optimization) |

## Collection subtranches and retained follow-ups

| ID / scope | Status and trigger | Predecessors and boundary | Model / effort | Evidence and research |
| --- | --- | --- | --- | --- |
| C10a-01 — normalized odd/mod-8 forms with multiplier ranking | **Open, worktree implementation in progress.** Opt-in parity2/coupled scoring and the revised bundle are checked in the C10a-01 worktree. Next: bounded Gray shifts and attainable A populations; select actual-norm ranking and cheaper root tables explicitly. Mainline code integration, refreshed source/work freezes, revised calibration, complete factors and adoption remain open. | Closed C10/C10a-07 control; A1a kernels, C8 scheduling and C9a intake keep their owners. | Astra / **xhigh** throughout this workstream | [Gates](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#c10a--conditional-collection-policy-questions-open); [42-mechanism assessment](https://github.com/nishanth17/factor/blob/23b71fa723b7982f1b1ceb18c426036e5e39290d/v2/studies/c10a07_research_assessment.md); [current handoff and remaining implementation](#c10a-01-implementation-remaining-after-the-coverage-review) |
| C10a-02 — actual-A-product population and inter-family diversity | Conditional on product error/overlap predicting useful post-filter work. The 128-A census demotes an immediate rewrite. | C10's rejected fixed 16-family ordering is a scoped loss; C10a-07 owns its completed traversal decision. | Sol / high census; xhigh selector/resume; high frozen execution | [Gates and census](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#c10a--conditional-collection-policy-questions-open) |
| C10a-03 — bounded per-job threshold feedback | Conditional on recovery cost or exact rejected-position audit supporting a declared loss budget. Preserve exact refinement; high rejection counts alone justify no omissions. | Settled collection control; C9a residual admission/target feedback remains separate. | Sol / xhigh bounds/coverage/implementation; high execution; Astra / xhigh only for unresolved rejection proof | [Gates](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#c10a--conditional-collection-policy-questions-open) |
| C10a-06 — separate runner tranche | Open: repeated lease markers are fixed only in C10a-05's new runner. Reconcile the separate runner before new execution; keep the frozen C10 runner unchanged. | C10/C10a-05; scope is runner/accounting, not a collection-policy promotion. Model was not separately specified; parent recommendation applies. | Sol / high frozen policy/runner; xhigh for changed accounting | [Recorded open item](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#c10a--conditional-collection-policy-questions-open) |
| C10a-04/G1 — larger-size feasibility after delivered routing | Conditional: declare one finite explicit 61/62-digit feasibility stage after current 60 evidence; practical completion/memory may justify a matched stage toward 65. Separately approved longer 65 or bounded 70 pilot remains possible. | Routing ends at 60. The 65 spike timed out both arms at 1,800 s; 8,500-entry/nine-A candidate has no completed-factor evidence. 61–65 interpolation is unvalidated. No automatic 61–80 sweep. | Sol / xhigh design/config reconciliation; high frozen execution; Astra / xhigh only for new arithmetic | [65 spike](benchmarks/RESULTS.md#c10a-04-larger-band-routing); [retained scope](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#c10a--conditional-collection-policy-questions-open) |
| A1a — direct extended versus pre-integration comparison | Conditional on a useful distinct hypothesis and new freeze; independent studies cannot establish this comparison by chaining gains. | Delivered A1a stays Done; selected-power-period/unroll4 extension remains opt-in after its losses. | Sol / high frozen comparison; xhigh only for changed semantics | [Open gate](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#a1a--selected-collector-kernel-delta-delivered-31-october-2026) |

A1a's quota-aware setup and lifted-root reuse are separate hypotheses. Bind
switch metadata to checked base/A/backend/B transitions, and lifted roots to
polynomial/discriminant/power/interval eligibility. Test characteristic two,
A/multiplier divisors, singular roots, recentering, invalidation, simultaneous
old/new storage and failed construction. Preserve order, quotas, public
validation and charged rebuild/refusal; derived state cannot serialize trust.
The old larger-cache losses remain relevant. The original sparse-queue six-call
screen and C10a-07 eight-call recovery screen were replaced by the completed
panel; their unrun historical designs are not live A1a work. See the
[panel](benchmarks/RESULTS.md#c10a-07-face-and-queue) and
[adoption](benchmarks/RESULTS.md#c10a-07-face-and-queue).

C9a's ordering remains target feedback/pacing attribution, then residual domains
and cycle weight; compact storage, splitter, batch and disk work need affected
cost/yield evidence. The accepted cohort's 5,402 solver calls across nine
factorizations motivates attribution, without proving avoidable work. Compare
wider recovered-product intake separately from asymmetric smaller-endpoint
bounds. Sample rejected positions and late prefixes as well as survivors.
Preserve all components, repeated primes, self-loops, parallel edges, original
atoms and square corrections. Do not sum overlapping dependency touches.
Coordinate admission/scoring changes through a joint freeze and checkpoint
identity; C10a owns solve scheduling, B7 representation proofs and F6 higher
arity. F7 batch intake needs a new domain/ownership proof: existing SmoothBatch
arithmetic does not cover the 5,000-prime SIQS base. Full
[C9a gates](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#c9a--remaining-large-prime-yield-and-storage-transfers-open)
remain attached.

## Numbered phase catalogue and acceptance gates

The execution rows own scheduling. This catalogue preserves each P ID and the
remaining acceptance boundary; full original A/E gates, quantitative evidence,
proposed grids, limitations and research are linked per item. Historical P2
unchecked gates stay open until settled by their mapped P8/P6 work.

### Phase 1 — Immediate correctness and API repairs

[Original phase scope and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#phase-1--immediate-correctness-and-api-repairs)

Done for the bounded core; see the completed register below.

### Phase 2 — Bounded scheduling, reusable infrastructure, and benchmarks

[Original phase scope and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#phase-2--bounded-scheduling-reusable-infrastructure-and-benchmarks)

| ID | Status | Execution owner / predecessors | Remaining scope and trigger | Full gates and evidence |
| --- | --- | --- | --- | --- |
| P2.1 | Open | P8.2 / E3 | Exact preprocessing exists; trial cutoff, proof-backed pruning and bounded Fermat tuning remain. | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p21--add-exact-preprocessing-and-classification-reuse) |
| P2.2 | Open | P8.3–P8.4 / E4 / G1 | Shared-budget dispatcher exists; tune retries/batches and joint allocation without hidden factor labels. | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p22--build-a-budgeted-portfolio-dispatcher) |
| P2.3 | Conditional | P8.6 / E6 | Streaming/capped caching exists; changed workloads must justify setup/reuse experiments. | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p23--stream-and-cache-prime-schedules-within-a-memory-cap) |
| P2.4 | Open | P8.5 / E5 | Recovery/checkpoints exist; separate GCD, polling, serialization, durable-write and replay costs. | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p24--add-chunk-checkpoints-and-a-consistent-recovery-protocol) |
| P2.5 | Open | P8.1 / D7; P6.4 / I1 | Corpus/runner exists; broad repeated and actual competitor comparisons remain. | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p25--freeze-a-benchmark-corpus-and-runner) |
| P2.6 | Conditional | P8.6 / E6 | Measured setup/consumption bottleneck; retain cache/rolling off and exact half-open streams. | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p26--port-reusable-sieve-contexts-and-rolling-strikes-selectively) |
| P2.7 | Conditional | P8.6 / E6 | Changed profile must justify wheel/presieve/packing work; retain wheel-6/bytearray and prior losses. | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p27--compare-wheelpre-sieve-and-packed-output-ideas-from-c-v3v4) |
| P2.8 | Conditional | P6.3 / F4; P8.7 / H1 | Feasibility is Done; first-factor benefit, aggregate CPU/RSS and cancellation still govern worker promotion. | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p28--establish-whether-parallel-candidate-search-is-worthwhile) |

### Phase 3 — Add the missing balanced-composite engine

[Original phase scope and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#phase-3--add-the-missing-balanced-composite-engine)

| ID | Status | Execution owner / predecessors | Remaining scope and trigger | Full gates and evidence |
| --- | --- | --- | --- | --- |
| P3.7 | Deferred | B5 | Measured remaining array bottleneck and usable PyPy dependency required; exact overflow/duplicate-hit and conversion controls. | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p37--optional-numpy-spike-for-bounded-siqs-sieve-arrays) |
| P3.8 | Open | A8/B6/B7/C4/D1/D2; collection owners below | Freeze real post-filter matrices, exact products and original-relation lifting. Capacity refusal can justify sparse work without dominant solve CPU. | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p38--evaluate-gf2-matrix-and-filtering-optimizations) |

### Phase 4 — Optimize ECM chains and arithmetic after coverage exists

[Original phase scope and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#phase-4--optimize-ecm-chains-and-arithmetic-after-coverage-exists)

| ID | Status | Execution owner / predecessors | Remaining scope and trigger | Full gates and evidence |
| --- | --- | --- | --- | --- |
| P4.2 | Conditional | B4/C11 | Bounded B4 is Done. Wider 64–1024-bit formula/reduction/normalization comparisons need a remaining arithmetic bottleneck, fresh whole-engine evidence and exact a24/width/checkpoint/nonunit controls. | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p42--compare-fused-and-normalized-montgomery-kernels) |
| P4.3 | Deferred | A2/A2a/C11 | Foundation and stopped study are Done. Further experiments stopped; independent 3–100-digit stage/size confirmation, selective GMP-helper policy and automatic backend selection remain unsettled. Resume only as a separately selected tranche. | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p43--compare-pypy-python-integers-and-optional-gmpmpz) |
| P4.4 | Conditional | C5 | Remaining reduction bottleneck after A2/B4; exact encoded-domain/unit-scaling, intermediate-width and recovery proofs before whole-loop comparison. | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p44--keep-barrettmontgomery-reducers-experimental-until-measured) |

### Phase 5 — Complementary methods and stronger continuations

[Original phase scope and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#phase-5--complementary-methods-and-stronger-continuations)

| ID | Status | Execution owner / predecessors | Remaining scope and trigger | Full gates and evidence |
| --- | --- | --- | --- | --- |
| P5.1 | Open | A5 then B8 | Binary Lucas stages 1/2, discriminant/parameter checks, exact powers, genuine stage two, bounded trials/saturation/resume; prove marginal portfolio value. | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p51--implement-williams-p1-with-binary-lucas-evaluation) |
| P5.2 | Conditional | C2/C3a/F3; G1 allocation | A3/B2 delivered programs/pairing retain scoped decisions; increased-B1 ECM points, advanced pruning/common-Z, polynomial continuation and allocation remain separate. | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p52--pair-ecm-stage-2-primes-and-tune-table-size) |
| P5.3 | Open | B8/F8; G1 allocation | A6 p−1 tranche is Done. Production-bound p−1 B2/B1 and gap/cache calibration remain with C11/G1; Lucas recurrence, p+1 extension and ECM old-point migration need independent composition/ratio/checkpoint proofs. | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p53--improve-p1-powering-and-continuation-independently) |
| P5.4 | Conditional | C9a/C1a; A8/B7/C4 | C1/C9 delivered DLP tranches are Done. Yield/intake/retention and stronger filtering keep separate triggers, exact residual certainty and original atoms. | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p54--add-double-large-prime-siqs-and-stronger-filtering) |

### Phase 6 — Research options, parallel execution, and claims

[Original phase scope and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#phase-6--research-options-parallel-execution-and-claims)

| ID | Status | Execution owner / predecessors | Remaining scope and trigger | Full gates and evidence |
| --- | --- | --- | --- | --- |
| P6.1 | Conditional | F2 | E1 and whole curve-family/stage-one cost/yield justification; distinguish curve torsion from selected-point order and validate stage-two maps. | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p61--compare-edwardswindowed-and-torsion-aware-ecm) |
| P6.2 | Conditional | F3/F5/F6/F7/F8 | Separate exact polynomial ECM, cross-family A=A0*q, higher-arity incidence, batch and specialized fast P±1 continuations. Each needs its own cost/domain trigger. | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p62--investigate-polynomial-continuations-and-richer-relation-collectors) |
| P6.3 | Conditional | B12/A7a/F4 | Stable serial control and worker contracts; compare 1/2/4 workers, fixed work and first-factor latency under aggregate limits. | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p63--add-reproducible-process-level-parallelism) |
| P6.4 | Open | I1 | H1, selected P6 decisions and committed-only validation; pinned feasible competitor comparisons with complete/factor-one and backend/core disclosures. | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p64--publish-reproducible-workload-specific-results) |

### Phase 7 — Implement and scale GNFS

[Original phase scope and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#phase-7--implement-and-scale-gnfs)

| ID | Status | Execution owner / predecessors | Remaining scope and trigger | Full gates and evidence |
| --- | --- | --- | --- | --- |
| P7.1 | Deferred | A9 | H1 and coverage review first; exact field/nonmonic norm/ideal/bad-prime contracts and bounded root-method coverage. | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p71--establish-exact-polynomial-and-relation-contracts) |
| P7.2 | Deferred | B9 | A9; bounded full rational/algebraic relation collection, primitive-pair deduplication and checked two-norm/ideal resume. | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p72--collect-bounded-rationalalgebraic-relations) |
| P7.3 | Deferred | B10/C7 | A9; exact ideal/sign/character constraints and original-relation lifting. Ideal parity alone does not prove algebraic squareness. | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p73--filter-and-solve-dependencies-with-gnfs-constraints) |
| P7.4 | Deferred | B11/C7 | A9 root scope; rational and algebraic roots with bounded coefficient/precision/sign/denominator recovery and modular mapping. | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p74--compute-rational-and-algebraic-square-roots-and-split-n) |
| P7.5 | Deferred | C7 | B9/B10/B11; complete small general composites, both square roots, proper splits, corrupt-state checks and bounded opt-in dispatch. | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p75--integrate-the-small-gnfs-engine-with-bounded-dispatch) |
| P7.6 | Deferred | D3–D6/E2 | C7 and A2 decision; selection, prime special-q before composite-q, general side-labelled cofactoring and matrix/storage scaling. | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p76--scale-polynomial-selection-sieving-and-matrix-work) |
| P7.7 | Deferred | F1 | E2/E1; fresh SIQS/GNFS crossover with completion/censoring, CPU/RSS/disk and cold/warm costs; separate SNFS-friendly forms. | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p77--measure-coverage-and-the-siqsgnfs-crossover) |

### Phase 8 — Revisit Phase 2 portfolio optimization

[Original phase scope and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#phase-8--revisit-phase-2-portfolio-optimization)

| ID | Status | Execution owner / predecessors | Remaining scope and trigger | Full gates and evidence |
| --- | --- | --- | --- | --- |
| P8.1 | Open | D7 | Deferred E1 must settle before freezing accepted current engines and untouched evaluation; preserve M12–M17 evidence. | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p81--freeze-the-new-control-and-untouched-evaluation-data) |
| P8.2 | Open | E3/A10a; A10 Done; B14 open | D7 for integrated tuning; justified isolated oracle spikes may precede it. Preserve exact rejection/certainty and prior losses. | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p82--reduce-preprocessing-work-with-exact-rejection-proofs) |
| P8.2.2 | Open | B14 | A10 is settled; start n−1/Pocklington with independent verification, bounded recursive proof/search/storage and requested-proof resume. | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p822--generate-and-independently-verify-prime-certificates-b14) |
| P8.3 | Open | E4 | D7; identical total grants/seeds, full recovery and resume; demonstrate marginal first-factor/complete value. | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p83--calibrate-bounded-rho-batches-walks-and-restarts) |
| P8.4 | Open | G1/C3a; C3 Done | D7/E1 and selected decisions; charge prior work and fit only observable metadata, with matched fresh confirmation. | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p84--tune-ecmp1-allocation-and-relation-engine-handoff) |
| P8.5 | Conditional | E5 | D7 and changed cost; preserve deadlines, saturation recovery, atomic commits, charged replay and corrupt-state checks. | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p85--measure-recovery-polling-and-checkpoint-costs-separately) |
| P8.6 | Conditional | E6/C11/A2a | D7 for changed-workload tuning; continue parameter provenance without repeating accepted chains/kernels or copying native thresholds. | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p86--revisit-context-setup-and-schedules-only-for-changed-workloads) |
| P8.7 | Open | H1 | G1 and E3–E6/selected P6 decisions; integrate confirmed winners and reconcile each historical P2 gate on untouched data. | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p87--integrate-accepted-winners-and-close-the-reconciliation) |

### P3.8 carry-over and mathematical boundaries

| Existing ID | Status / next scope | Trigger, predecessors and owner | Evidence |
| --- | --- | --- | --- |
| P3.8-R1 | Open beyond completed B1 calibration: extend independent corpus/runner coverage and larger feasible calibration; actual-A product/interval/factor-count/reuse and classical MPQS dual-batch roots remain coupled questions. | Useful finite collection, stratified smaller-factor ECM classes and untouched confirmation. Coordinate C10a-02/-07 decisions, C3a/G1 and explicit larger-size feasibility; no native crossover threshold. Sol / high calibration; xhigh selection. | [Original carry-over](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p38--evaluate-gf2-matrix-and-filtering-optimizations) |
| P3.8-R2 | Conditional B13/C8 capacity and whole-interval CRT; bounded A1 integration is Done. | New calibrated refusal or eligible-prime/root-cost evidence. Astra / xhigh. | [R2 and deferrals](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p38--evaluate-gf2-matrix-and-filtering-optimizations) |
| P3.8-R3 | Conditional cadence-8/general filtering-policy reassessment and larger provenance/fill bounds. Bounded mixed-row/cadence-32 opt-in work is Done. | Fresh scoped evidence under revised policy without held-out retuning; A8/B7/C4 own storage/filtering. Preserve dense reservations until new representations are proved. Astra / xhigh representation; Sol / high frozen cadence study. | [R3 gates and rejected layouts](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p38--evaluate-gf2-matrix-and-filtering-optimizations) |
| P3.8-R4 | Conditional residual-classification/certification reuse after delivered C1; preserve certainty. | Profile must show repeated intake/atomic primality cost. A base-squared shortcut needs complete eligible-prime coverage including A/multiplier exceptions; opaque composite cofactors need a separate type/checkpoint proof. C9a/E3 coordinate; Astra / xhigh changed proof. | [R4 gates](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p38--evaluate-gf2-matrix-and-filtering-optimizations) |
| P3.8-R5 | Deferred E1 final SSS/worker/current-control comparison after Done A7 reconciliation. | Worker flyer/large-interval compatibility and 40-digit SSS calibration need explicit decisions. Clarify upstream root license before further adaptation; use assignment IDs for seeds, include transport/conversion/loss/replay costs. Sol / high reconciliation; xhigh accounting. | [R5/E1 scope](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#a1c--same-modulus-two-root-interleaving-conditional); [pinned arms](benchmarks/inputs/controls/a7_r5_e1_arms.json) |

For matrix work, `M` has constraint rows and relation columns; `M d = 0`
selects original relations. A8 freezes orientation, transpose products, lifting
and genuine fixtures, including the independent bilinear-product identity.
B6 precedes C4 by investment policy. D1/D2 need B6/C4 decisions or documented
dense-memory infeasibility; D2 has no D1 dependency. Four Russians must bound
`2**k` tables, fill, simultaneous provenance and recovery. Lanczos must correct
Gram-kernel candidates against original `M`, including self-orthogonal and
singular-block cases. Wiedemann needs a genuine block generator, finite
sequences/reconstruction/I/O and original verification, starting from an exact
base case. All solver checkpoints bind original relation/matrix identities,
orientation, seeds, widths, lifting and consumed allowances. Exact square
congruences and proper divisors remain required; fewer rows or faster XOR
alone cannot promote a solver. The
[full matrix A/E gates](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#a1c--same-modulus-two-root-interleaving-conditional)
retain the counterexamples, widths, references and finite failure contracts.

B14 starts with n−1/Pocklington rather than expanding automatically to APRT-CLE
(AP28). Independently verify recursive prime factors, exponent divisibility,
factored-part bounds and witnesses. Proof exhaustion leaves probable certainty;
compositeness witnesses return to factoring with work retained. Complete and
fully proven results differ. The exact example, corruption/resume cases and
coverage limits remain in its [full gates](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p822--generate-and-independently-verify-prime-certificates-b14).

### Optional product work

A11a is Conditional on explicit product selection, outside optimization priority.
AP38 covers validated hints and opt-in lookup with finite network/cache work;
AP39 covers bounded batch/expression parse/exponent/output and per-item failures;
AP40 covers browser progress/offline/save-continue using checked cumulative
checkpoints. Hidden corpus hints cannot enter timings. Sol / high; xhigh for
changed accounting/resume. [Full boundaries](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#optional-product-follow-ups--separate-from-factoring-optimization).

## Completed register

- **A1b existing fixed-score review** — 1 November 2026 study; 2 November integration; **Done, retain defaults**. All 30 prefix calls on nine exposed moduli; fixed/list loses at 50/55/60 and no band advances. Harness/accounting repairs and exact growing-DLP checks complete; no complete-factor claim or production change. [Final decision](benchmarks/RESULTS.md#a1b-existing-fixed-score-review), [numerical evidence](https://github.com/nishanth17/factor/blob/23b71fa723b7982f1b1ceb18c426036e5e39290d/v2/studies/a1b_final_results.json).

One line per delivered item. Dates below are the recorded acceptance/decision
dates; later integration dates are retained when distinct. Detailed results,
intervals, sample/input counts, timings, failures and rejected arms remain at
each linked location. Remaining scopes are listed above.

- **P1.1–P1.9** — 3 October 2026; **Done**; Exact arithmetic/API/sieve/recovery repairs; PyPy core accepted. [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#phase-1--immediate-correctness-and-api-repairs)
- **P2 core / M13 / M14 / M17** — 3 October 2026; **Done**; Bounded core and fresh loop confirmation accepted; optional tuning gates remain open. [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#phase-2--bounded-scheduling-reusable-infrastructure-and-benchmarks)
- **P3.1** — 3 October 2026; **Done**; Exact polynomial/relation contracts accepted. [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p31--implement-polynomial-and-relation-identities)
- **P3.2** — 3 October 2026; **Done**; Bounded SLP collection accepted. [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p32--build-a-bounded-single-large-prime-collector)
- **P3.3 / M28** — 3 October 2026; **Done**; Verified filtering, lifting and extraction; retained resieve/scoring decisions. [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p33--implement-filtering-gf2-dependencies-and-factor-extraction)
- **P3.4 / M31** — 3 October 2026; **Done**; Bounded SIQS/evaluation; fresh 50-digit split in 1,183.355 s; broad 60–80 scaling unestablished. Current scored-multiplier routing is in STATUS. [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p34--add-siqs-self-initialization-and-bounded-dispatch)
- **A11 / P3.4 CLI** — 9 October 2026 (reconciled); **Done**; Explicit QS/MPQS/SIQS usability accepted; allocation stays separately owned. [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p34-usability-follow-up--explicit-qsmpqssiqs-cli-access-a11)
- **P3.5** — 3 October 2026; **Done**; Bounded SSS/SSSf adapter/upstream comparison; experimental, no automatic promotion. [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p35--evaluate-smooth-subsum-search-on-the-same-interface)
- **P3.6** — 3 October 2026; **Done**; Bounded worker study retains serial; complete-factor promotion gate failed. [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p36--evaluate-coarse-parallel-relationcandidate-collection)
- **P3.6.1** — 4 October 2026; **Done**; Repair tranches accepted; fresh complete cohorts improve 70.7% SIQS, 77.2% SSS, 78.6% SSSf; native serial retained. [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p361--immediate-diagnosis-and-improvement-of-p35p36)
- **P3.8-R1 implementation** — 4 October 2026; **Done**; Streamed assignments/external-square MPQS/capacity/resume; upper one-second bands censored. [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p38--evaluate-gf2-matrix-and-filtering-optimizations)
- **P3.8-R3 bounded tranche** — 4 October 2026; **Done**; Mixed identities/recovery accepted; cadence-32 opt-in saves 38.7%, 54/54 per arm; general policy open. [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p38--evaluate-gf2-matrix-and-filtering-optimizations)
- **A1 / P3.8-R2** — 4 October; reverified 5 October 2026; **Done**; Combined integration 758d5b3; fixed-score/plans originally opt-in; retain-default verdict preserved. [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-a--independent-foundations-that-can-proceed-now)
- **A2 / P4.3** — 5 October 2026; **Done**; Backend foundation; wider study stopped. Native integers retained, disjoint calibration open. [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p43--compare-pypy-python-integers-and-optional-gmpmpz)
- **A3 / P5.2 programs** — 5 October 2026; **Done**; Reusable programs opt-in; small complete factoring loses, streamed default retained. [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-a--independent-foundations-that-can-proceed-now)
- **A4 / P4.1 correctness** — 5 October 2026; **Done**; Verified bounded PRAC records with independent point/exception controls; production routed separately under B3. [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p41--finish-prac-repair-and-precompute-valid-chains)
- **A6 / P5.3 p−1** — 9 October 2026; **Done**; Adopted recurrence/chunk-64 default; fresh stages save 9.73%; extended integer portfolio Inconclusive. [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p53--improve-p1-powering-and-continuation-independently)
- **A7 / P3.8-R5 reconciliation** — 9 October 2026; **Done**; Interface/resume reconciliation accepted; SSS/workers stay experimental and E1 comparison deferred. [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-a--independent-foundations-that-can-proceed-now)
- **A10 / P8.2.1** — 9 October 2026; **Done**; Verified deterministic-range/v1 transfers; range-only selection retained, B14 proofs remain open. [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p821--transfer-verified-v1-advantages-starting-with-primality-a10)
- **B1 / P3.8-R1 calibration** — 9 October 2026; **Done**; Balanced 30-digit saving 36.6%; regressions/small populations retain defaults; 40-digit 6.9% misses frozen gate. [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-b--exploit-the-settled-controls)
- **B2 / P5.2 pairing** — 9 October 2026; **Done**; Both bounded pairing layouts lose, remain opt-in; no completion gain. [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p52--pair-ecm-stage-2-primes-and-tune-table-size)
- **B3 / P4.1 production** — 9 October acceptance; 10 October 2026 integration; **Done**; Adopted native PRAC/batch16 default with B4 fallback; wider allocation/confirmation open. [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#b3-production-integration-acceptance-9-october-2026)
- **B4 / P4.2** — 9 October 2026; **Done**; Native bridge saves 6.79% (95% 5.82–13.58%); GMP whole-ladder bridge Inconclusive, readable control retained. [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p42--compare-fused-and-normalized-montgomery-kernels)
- **C6 / P4.5** — 9 October 2026; **Done**; Reused native PRAC gains 6.73%, GMP Lucas 8.42%; fresh native retains ladder in this isolated study. [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p45--experiment-with-reference-ecm-precomputed-lucas-chains-c6)
- **C1 / P3.8-R4 / P5.4** — 10 October 2026; **Done**; Explicit calibrated 40-digit DLP saves 31.73% (95% 29.21–33.99%); original SLP defaults retained for that decision. [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p54--add-double-large-prime-siqs-and-stronger-filtering)
- **C3 bounded allocation** — 10 October 2026; **Done**; Adopted compact64 at 30/escalating140 at 40 before fresh confirmation; oversized/source-screen gates Inconclusive. [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#calibrate-automatic-selection-and-handoff-c3-then-g1)
- **C9 / C9.2** — 11 October 2026; **Done**; Explicit 50-digit growing DLP saves 38.10% versus SLP / 53.08% versus fixed DLP, 9/9 per arm; higher memory. [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#c9--source-derived-large-prime-yield-and-dlp-economics-complete-and-integrated)
- **C10** — 11 October 2026; **Done**; 2 MiB plan saves 21.83%/17.79% at 40/50; scoped automatic integration. Fixed 16-family A-order promotion Rejected. [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#c10--source-derived-collection-strategies-closed-bounded-tranche)
- **C10a-05** — 11 October 2026; **Done**; 50-digit collection pacing saves 25.99% (95% 22.87–28.24%); 40-digit extension Inconclusive, window retained. [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#c10a--conditional-collection-policy-questions-open)
- **C10a-05 size-scope extension** — 11 October 2026; **Done**; 51-digit bundle completes 9/9 versus 0/9 prior default; completion gain, not speedup against successful controls. [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#c10a--conditional-collection-policy-questions-open)
- **C10a-04 delivered routing** — 12 October 2026; **Done**; Closed: compact defaults through 60 despite unstable 55 training, preliminary 60 confirmation and unmeasured interpolation; original campaign incomplete. [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#c10a--conditional-collection-policy-questions-open)
- **A1a delivered integration/extension** — 31 October 2026; **Done**; Native 50–55 adoption: 50 passes (6.008%); 55 Inconclusive, 51–54 unmeasured. Combined extension loses 5.285%/2.994%, opt-in. [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#a1a--selected-collector-kernel-delta-delivered-31-october-2026)
- **C10a-07** — 1 November 2026; **Done**; Closed after combined face/queue adoption; 23.1% wall / 23.0% CPU exploratory on nine 50/55/60 inputs; drift gate unmet, extra 18 calls unrun. [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#c10a--conditional-collection-policy-questions-open)
- **C11 first-pass audit** — by 1 November 2026; **Done**; Parameter/source reconciliation only; ongoing calibration remains Open. [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#c11--cross-engine-constant-provenance-and-calibration-open)

## Phase checkpoints

[Original phase exit gates](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#phase-checkpoints).

| Phase | Concrete deliverable | Exit decision |
| --- | --- | --- |
| 1 | Native Python 3 repairs, result/certainty contracts, regression suite | All counterexamples and failure paths pass; freeze corrected baseline |
| 2 | Budgeted portfolio, bounded schedules, recovery/resume, corpus/runner | Time/memory limits and reproducibility demonstrated |
| 3 | Verified relation pipeline, SIQS, SSS comparison | Balanced-composite completion improves under equal limits |
| 4 | Validated PRAC/kernel/backend/reducer experiments | Promote measured winners; retain baseline for losses |
| 5 | p+1, paired/recoverable continuations, double-large-prime SIQS | Complementary completion gain without resource regression |
| 6 | Curve-family/research decisions, fixed-core workers, evidence-backed README | Held-out results support every adopted default and claim |
| 7 | Bounded GNFS reference, scaling pipeline, SIQS/GNFS crossover | Separate correctness, scaling and held-out dispatch gates |
| 8 | Reconciled Phase 2 calibration, preprocessing and overhead experiments | Fresh held-out promotion or explicit retain/defer/reject decisions |

Execution priority after Phase 3 is current-engine calibration, selected
P4–P6 improvements and fresh Phase 8 portfolio confirmation. Existing phase
numbers preserve references. GNFS remains deferred until H1 and a coverage
review justify activating the later tranche.

## Historical heading anchors

These compatibility anchors preserve the old public fragments. Follow “Detail”
for the original narrative; STATUS and the live rows above govern current work.

<a id="factor-phased-todos-and-decision-gates"></a>

- Factor: phased TODOs and decision gates — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#factor-phased-todos-and-decision-gates).

<a id="how-to-use-the-gates"></a>

- How to use the gates — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#how-to-use-the-gates).

<a id="execution-order-optimize-current-engines-then-revisit-gnfs"></a>

- Execution order: optimize current engines, then revisit GNFS — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-order-optimize-current-engines-then-revisit-gnfs).

<a id="master-execution-sequence-and-model-recommendations"></a>

- Master execution sequence and model recommendations — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#master-execution-sequence-and-model-recommendations).

<a id="how-to-read-the-phases"></a>

- How to read the phases — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#how-to-read-the-phases).

<a id="current-foundation-already-available"></a>

- Current foundation: already available — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#current-foundation-already-available).

<a id="model-and-effort-key"></a>

- Model and effort key — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#model-and-effort-key).

<a id="dependency-overview"></a>

- Dependency overview — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#dependency-overview).

<a id="parallel-execution-and-ownership"></a>

- Parallel execution and ownership — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#parallel-execution-and-ownership).

<a id="acceptance-evidence-and-stop-rules"></a>

- Acceptance, evidence and stop rules — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#acceptance-evidence-and-stop-rules).

<a id="sources-and-verification-scope"></a>

- Sources and verification scope — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#sources-and-verification-scope).

<a id="p11--port-the-seven-modules-to-explicit-python-3-arithmetic"></a>

- P1.1 — Port the seven modules to explicit Python 3 arithmetic — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p11--port-the-seven-modules-to-explicit-python-3-arithmetic).

<a id="p12--repair-suyama-ecm-initialization"></a>

- P1.2 — Repair Suyama ECM initialization — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p12--repair-suyama-ecm-initialization).

<a id="p13--repair-and-validate-every-sieve-endpoint"></a>

- P1.3 — Repair and validate every sieve endpoint — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p13--repair-and-validate-every-sieve-endpoint).

<a id="p14--preserve-unresolved-cofactors-and-define-input-behavior"></a>

- P1.4 — Preserve unresolved cofactors and define input behavior — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p14--preserve-unresolved-cofactors-and-define-input-behavior).

<a id="p15--repair-p1-batching-and-saturation-handling"></a>

- P1.5 — Repair p−1 batching and saturation handling — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p15--repair-p1-batching-and-saturation-handling).

<a id="p16--repair-rho-retries-and-cap-recovery"></a>

- P1.6 — Repair rho retries and cap recovery — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p16--repair-rho-retries-and-cap-recovery).

<a id="p17--repair-ecm-curve-limits-and-failure-transitions"></a>

- P1.7 — Repair ECM curve limits and failure transitions — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p17--repair-ecm-curve-limits-and-failure-transitions).

<a id="p18--clarify-primality-and-helper-contracts"></a>

- P1.8 — Clarify primality and helper contracts — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p18--clarify-primality-and-helper-contracts).

<a id="p19--quarantine-unsafe-prac-and-make-diagnostics-reliable"></a>

- P1.9 — Quarantine unsafe PRAC and make diagnostics reliable — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p19--quarantine-unsafe-prac-and-make-diagnostics-reliable).

<a id="p21--add-exact-preprocessing-and-classification-reuse"></a>

- P2.1 — Add exact preprocessing and classification reuse — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p21--add-exact-preprocessing-and-classification-reuse).

<a id="p22--build-a-budgeted-portfolio-dispatcher"></a>

- P2.2 — Build a budgeted portfolio dispatcher — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p22--build-a-budgeted-portfolio-dispatcher).

<a id="p23--stream-and-cache-prime-schedules-within-a-memory-cap"></a>

- P2.3 — Stream and cache prime schedules within a memory cap — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p23--stream-and-cache-prime-schedules-within-a-memory-cap).

<a id="p24--add-chunk-checkpoints-and-a-consistent-recovery-protocol"></a>

- P2.4 — Add chunk checkpoints and a consistent recovery protocol — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p24--add-chunk-checkpoints-and-a-consistent-recovery-protocol).

<a id="p25--freeze-a-benchmark-corpus-and-runner"></a>

- P2.5 — Freeze a benchmark corpus and runner — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p25--freeze-a-benchmark-corpus-and-runner).

<a id="p26--port-reusable-sieve-contexts-and-rolling-strikes-selectively"></a>

- P2.6 — Port reusable sieve contexts and rolling strikes selectively — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p26--port-reusable-sieve-contexts-and-rolling-strikes-selectively).

<a id="p27--compare-wheelpre-sieve-and-packed-output-ideas-from-c-v3v4"></a>

- P2.7 — Compare wheel/pre-sieve and packed-output ideas from C v3/v4 — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p27--compare-wheelpre-sieve-and-packed-output-ideas-from-c-v3v4).

<a id="p28--establish-whether-parallel-candidate-search-is-worthwhile"></a>

- P2.8 — Establish whether parallel candidate search is worthwhile — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p28--establish-whether-parallel-candidate-search-is-worthwhile).

<a id="p31--implement-polynomial-and-relation-identities"></a>

- P3.1 — Implement polynomial and relation identities — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p31--implement-polynomial-and-relation-identities).

<a id="p32--build-a-bounded-single-large-prime-collector"></a>

- P3.2 — Build a bounded single-large-prime collector — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p32--build-a-bounded-single-large-prime-collector).

<a id="p33--implement-filtering-gf2-dependencies-and-factor-extraction"></a>

- P3.3 — Implement filtering, GF(2) dependencies, and factor extraction — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p33--implement-filtering-gf2-dependencies-and-factor-extraction).

<a id="p34--add-siqs-self-initialization-and-bounded-dispatch"></a>

- P3.4 — Add SIQS self-initialization and bounded dispatch — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p34--add-siqs-self-initialization-and-bounded-dispatch).

<a id="p34-usability-follow-up--explicit-qsmpqssiqs-cli-access-a11"></a>

- P3.4 usability follow-up — explicit QS/MPQS/SIQS CLI access (A11) — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p34-usability-follow-up--explicit-qsmpqssiqs-cli-access-a11).

<a id="p35--evaluate-smooth-subsum-search-on-the-same-interface"></a>

- P3.5 — Evaluate Smooth Subsum Search on the same interface — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p35--evaluate-smooth-subsum-search-on-the-same-interface).

<a id="p36--evaluate-coarse-parallel-relationcandidate-collection"></a>

- P3.6 — Evaluate coarse parallel relation/candidate collection — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p36--evaluate-coarse-parallel-relationcandidate-collection).

<a id="p361--immediate-diagnosis-and-improvement-of-p35p36"></a>

- P3.6.1 — Immediate diagnosis and improvement of P3.5/P3.6 — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p361--immediate-diagnosis-and-improvement-of-p35p36).

<a id="p37--optional-numpy-spike-for-bounded-siqs-sieve-arrays"></a>

- P3.7 — Optional NumPy spike for bounded SIQS sieve arrays — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p37--optional-numpy-spike-for-bounded-siqs-sieve-arrays).

<a id="p38--evaluate-gf2-matrix-and-filtering-optimizations"></a>

- P3.8 — Evaluate GF(2) matrix and filtering optimizations — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p38--evaluate-gf2-matrix-and-filtering-optimizations).

<a id="c10--source-derived-collection-strategies-closed-bounded-tranche"></a>

- C10 — Source-derived collection strategies (closed bounded tranche) — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#c10--source-derived-collection-strategies-closed-bounded-tranche).

<a id="c10a--conditional-collection-policy-questions-open"></a>

- C10a — Conditional collection-policy questions (open) — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#c10a--conditional-collection-policy-questions-open).

<a id="a1a--selected-collector-kernel-delta-delivered-31-october-2026"></a>

- A1a — Selected collector-kernel delta (delivered 31 October 2026) — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#a1a--selected-collector-kernel-delta-delivered-31-october-2026).

<a id="a1b--existing-fixed-score-promotion-follow-up-after-completed-a1-and-c10"></a>

- A1b — Completed existing fixed-score review — [Final decision](benchmarks/RESULTS.md#a1b-existing-fixed-score-review); [original gates](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#a1b--existing-fixed-score-promotion-follow-up-after-completed-a1-and-c10).

<a id="a1c--same-modulus-two-root-interleaving-conditional"></a>

- A1c — Same-modulus two-root interleaving (conditional) — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#a1c--same-modulus-two-root-interleaving-conditional).

<a id="p41--finish-prac-repair-and-precompute-valid-chains"></a>

- P4.1 — Finish PRAC repair and precompute valid chains — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p41--finish-prac-repair-and-precompute-valid-chains).

<a id="b3-production-integration-acceptance-9-october-2026"></a>

- B3 production integration acceptance (9 October 2026) — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#b3-production-integration-acceptance-9-october-2026).

<a id="p42--compare-fused-and-normalized-montgomery-kernels"></a>

- P4.2 — Compare fused and normalized Montgomery kernels — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p42--compare-fused-and-normalized-montgomery-kernels).

<a id="p43--compare-pypy-python-integers-and-optional-gmpmpz"></a>

- P4.3 — Compare PyPy Python integers and optional GMP/`mpz` — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p43--compare-pypy-python-integers-and-optional-gmpmpz).

<a id="p44--keep-barrettmontgomery-reducers-experimental-until-measured"></a>

- P4.4 — Keep Barrett/Montgomery reducers experimental until measured — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p44--keep-barrettmontgomery-reducers-experimental-until-measured).

<a id="p45--experiment-with-reference-ecm-precomputed-lucas-chains-c6"></a>

- P4.5 — Experiment with reference ECM precomputed Lucas chains (C6) — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p45--experiment-with-reference-ecm-precomputed-lucas-chains-c6).

<a id="p51--implement-williams-p1-with-binary-lucas-evaluation"></a>

- P5.1 — Implement Williams p+1 with binary Lucas evaluation — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p51--implement-williams-p1-with-binary-lucas-evaluation).

<a id="p52--pair-ecm-stage-2-primes-and-tune-table-size"></a>

- P5.2 — Pair ECM stage-2 primes and tune table size — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p52--pair-ecm-stage-2-primes-and-tune-table-size).

<a id="p53--improve-p1-powering-and-continuation-independently"></a>

- P5.3 — Improve p±1 powering and continuation independently — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p53--improve-p1-powering-and-continuation-independently).

<a id="p54--add-double-large-prime-siqs-and-stronger-filtering"></a>

- P5.4 — Add double-large-prime SIQS and stronger filtering — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p54--add-double-large-prime-siqs-and-stronger-filtering).

<a id="c9--source-derived-large-prime-yield-and-dlp-economics-complete-and-integrated"></a>

- C9 — Source-derived large-prime yield and DLP economics (complete and integrated) — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#c9--source-derived-large-prime-yield-and-dlp-economics-complete-and-integrated).

<a id="c9a--remaining-large-prime-yield-and-storage-transfers-open"></a>

- C9a — Remaining large-prime yield and storage transfers (open) — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#c9a--remaining-large-prime-yield-and-storage-transfers-open).

<a id="c1a--cheaper-slp-anchor-delta-after-c9"></a>

- C1a — Cheaper SLP anchor delta (after C9) — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#c1a--cheaper-slp-anchor-delta-after-c9).

<a id="p61--compare-edwardswindowed-and-torsion-aware-ecm"></a>

- P6.1 — Compare Edwards/windowed and torsion-aware ECM — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p61--compare-edwardswindowed-and-torsion-aware-ecm).

<a id="p62--investigate-polynomial-continuations-and-richer-relation-collectors"></a>

- P6.2 — Investigate polynomial continuations and richer relation collectors — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p62--investigate-polynomial-continuations-and-richer-relation-collectors).

<a id="p63--add-reproducible-process-level-parallelism"></a>

- P6.3 — Add reproducible process-level parallelism — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p63--add-reproducible-process-level-parallelism).

<a id="p64--publish-reproducible-workload-specific-results"></a>

- P6.4 — Publish reproducible workload-specific results — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p64--publish-reproducible-workload-specific-results).

<a id="p71--establish-exact-polynomial-and-relation-contracts"></a>

- P7.1 — Establish exact polynomial and relation contracts — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p71--establish-exact-polynomial-and-relation-contracts).

<a id="p72--collect-bounded-rationalalgebraic-relations"></a>

- P7.2 — Collect bounded rational/algebraic relations — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p72--collect-bounded-rationalalgebraic-relations).

<a id="p73--filter-and-solve-dependencies-with-gnfs-constraints"></a>

- P7.3 — Filter and solve dependencies with GNFS constraints — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p73--filter-and-solve-dependencies-with-gnfs-constraints).

<a id="p74--compute-rational-and-algebraic-square-roots-and-split-n"></a>

- P7.4 — Compute rational and algebraic square roots and split n — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p74--compute-rational-and-algebraic-square-roots-and-split-n).

<a id="p75--integrate-the-small-gnfs-engine-with-bounded-dispatch"></a>

- P7.5 — Integrate the small GNFS engine with bounded dispatch — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p75--integrate-the-small-gnfs-engine-with-bounded-dispatch).

<a id="p76--scale-polynomial-selection-sieving-and-matrix-work"></a>

- P7.6 — Scale polynomial selection, sieving and matrix work — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p76--scale-polynomial-selection-sieving-and-matrix-work).

<a id="p77--measure-coverage-and-the-siqsgnfs-crossover"></a>

- P7.7 — Measure coverage and the SIQS/GNFS crossover — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p77--measure-coverage-and-the-siqsgnfs-crossover).

<a id="p81--freeze-the-new-control-and-untouched-evaluation-data"></a>

- P8.1 — Freeze the new control and untouched evaluation data — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p81--freeze-the-new-control-and-untouched-evaluation-data).

<a id="p82--reduce-preprocessing-work-with-exact-rejection-proofs"></a>

- P8.2 — Reduce preprocessing work with exact rejection proofs — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p82--reduce-preprocessing-work-with-exact-rejection-proofs).

<a id="p821--transfer-verified-v1-advantages-starting-with-primality-a10"></a>

- P8.2.1 — Transfer verified v1 advantages, starting with primality (A10) — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p821--transfer-verified-v1-advantages-starting-with-primality-a10).

<a id="p822--generate-and-independently-verify-prime-certificates-b14"></a>

- P8.2.2 — Generate and independently verify prime certificates (B14) — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p822--generate-and-independently-verify-prime-certificates-b14).

<a id="later-integrated-preprocessing-tuning-e3"></a>

- Later integrated preprocessing tuning (E3) — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#later-integrated-preprocessing-tuning-e3).

<a id="p83--calibrate-bounded-rho-batches-walks-and-restarts"></a>

- P8.3 — Calibrate bounded rho batches, walks and restarts — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p83--calibrate-bounded-rho-batches-walks-and-restarts).

<a id="p84--tune-ecmp1-allocation-and-relation-engine-handoff"></a>

- P8.4 — Tune ECM/p−1 allocation and relation-engine handoff — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p84--tune-ecmp1-allocation-and-relation-engine-handoff).

<a id="calibrate-automatic-selection-and-handoff-c3-then-g1"></a>

- Calibrate automatic selection and handoff (C3, then G1) — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#calibrate-automatic-selection-and-handoff-c3-then-g1).

<a id="c3a--source-derived-policy-delta-after-c3"></a>

- C3a — Source-derived policy delta (after C3) — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#c3a--source-derived-policy-delta-after-c3).

<a id="p85--measure-recovery-polling-and-checkpoint-costs-separately"></a>

- P8.5 — Measure recovery, polling and checkpoint costs separately — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p85--measure-recovery-polling-and-checkpoint-costs-separately).

<a id="p86--revisit-context-setup-and-schedules-only-for-changed-workloads"></a>

- P8.6 — Revisit context setup and schedules only for changed workloads — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p86--revisit-context-setup-and-schedules-only-for-changed-workloads).

<a id="c11--cross-engine-constant-provenance-and-calibration-open"></a>

- C11 — Cross-engine constant provenance and calibration (open) — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#c11--cross-engine-constant-provenance-and-calibration-open).

<a id="p87--integrate-accepted-winners-and-close-the-reconciliation"></a>

- P8.7 — Integrate accepted winners and close the reconciliation — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p87--integrate-accepted-winners-and-close-the-reconciliation).

<a id="optional-product-follow-ups--separate-from-factoring-optimization"></a>

- Optional product follow-ups — separate from factoring optimization — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#optional-product-follow-ups--separate-from-factoring-optimization).

<a id="existing-evidence-to-reuse"></a>

- Existing evidence to reuse — [Detail](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#existing-evidence-to-reuse).
