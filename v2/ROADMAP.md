# Factor roadmap

Updated 9 November 2026.

This document ranks future work, records work in progress and lists completed scopes.
[STATUS.md](STATUS.md) owns current behavior and defaults; [POLICY.md](POLICY.md)
owns operating and acceptance rules; [RESULTS.md](benchmarks/RESULTS.md) owns measurements.
A completed component can have separate unfinished follow-ups.

Tables show primary model/effort assignments; linked contracts retain proof and
separately approved frozen-execution exceptions. No scope or acceptance gate is
removed by this navigation cleanup.

## Future work

### Ranked queue

Relative priority among remaining candidates is unchanged; numbering is compacted
after removing active and completed scopes. Dependencies and measured cost determine
when an item is eligible. Conditional items need their stated trigger.

Correctness/accounting repairs and explicitly requested proof work take precedence.

| Order | Item and next deliverable | Dependency or decision | Model / effort |
| ---: | --- | --- | --- |
| 1 | [E6: Native prime batches and reusable contexts](https://github.com/nishanth17/factor/blob/cb37da3ea5fbbeb44923faa87fd45d4e5d926bc4/v2/ROADMAP.md#execution-phase-e--current-engine-acceptance-control-freeze-and-p8-experiments) | Current context/cache contracts; D7 for integrated tuning. Screen generation/setup cost first; include conversion and unused preparation. | Sol / high screen; xhigh implementation |
| 2 | <a id="b15--native-helpersimd-opportunity-selection"></a>[B15: Select one native helper opportunity](https://github.com/nishanth17/factor/blob/cb37da3ea5fbbeb44923faa87fd45d4e5d926bc4/v2/ROADMAP.md#execution-phase-b--exploit-the-settled-controls) | Current receipts; C5 selection needs a useful B4a context. Reuse current cost evidence; choose at most one candidate or no-go. Implementation stays with its owner. | Sol / high; xhigh boundary design |
| 3 | <a id="h2--profile-driven-overhead-sweep-of-integrated-production"></a>[H2: Integrated production overhead sweep](https://github.com/nishanth17/factor/blob/cb37da3ea5fbbeb44923faa87fd45d4e5d926bc4/v2/ROADMAP.md#h2--profile-driven-overhead-sweep-of-integrated-production) | Full sweep follows the B12/A7a/B4a integration checkpoint; independent CLI/preprocessing work can start earlier. | Sol / xhigh |
| 4 | [C3a + C11: Portfolio opportunity screen and one policy delta](https://github.com/nishanth17/factor/blob/cb37da3ea5fbbeb44923faa87fd45d4e5d926bc4/v2/ROADMAP.md#new-source-derived-deltas--separate-follow-up-items) | C3 decisions and current cost evidence; no E1 gate for component work. Revalidate dispatch/allocation costs after H3/A7a and relevant worker/native changes; reuse compatible profiles. | Sol / xhigh; Astra for new proofs |
| 5 | <a id="a1d--prime-major-bucket-generation-across-polynomials"></a>[A1d: Prime-major polynomial buckets](https://github.com/nishanth17/factor/blob/cb37da3ea5fbbeb44923faa87fd45d4e5d926bc4/v2/ROADMAP.md#new-source-derived-deltas--separate-follow-up-items) | Current C10a-01 family/queue identity. Stable family/queue contract; charge construction and unused preparation. | Sol / xhigh; Astra for new arithmetic |
| 6 | [B3a: Larger-bound checked ECM chains](https://github.com/nishanth17/factor/blob/cb37da3ea5fbbeb44923faa87fd45d4e5d926bc4/v2/ROADMAP.md#execution-phase-b--exploit-the-settled-controls) | C3a/C11 cost evidence and B3/C6 oracles. Current stage-one cost must justify extending existing verified machinery. | Sol / xhigh; Astra for new chain mathematics |
| 7 | [C9b: Broader QS validation and remaining headroom](https://github.com/nishanth17/factor/blob/cb37da3ea5fbbeb44923faa87fd45d4e5d926bc4/v2/ROADMAP.md#new-source-derived-deltas--separate-follow-up-items) | Integrated C9a; coordinate shared H2 evidence. Use the integrated collector and fresh inputs; routine pilot checks remain inside each task. | Sol / xhigh |
| 8 | [C2: Common-Z and classical continuation changes](https://github.com/nishanth17/factor/blob/cb37da3ea5fbbeb44923faa87fd45d4e5d926bc4/v2/ROADMAP.md#execution-phase-c--measured-current-engine-challengers) | B2 decision and current continuation cost; coordinate C3a/F3. Keep coordinate and geometry effects attributable against the current backend. | Astra / xhigh |
| 9 | <a id="f7--siqs-directtree-screening-and-distinct-residual-batches"></a>[F7: SIQS direct/tree smooth-part screening](https://github.com/nishanth17/factor/blob/cb37da3ea5fbbeb44923faa87fd45d4e5d926bc4/v2/ROADMAP.md#execution-phase-f--independent-p6-challengers) | Current collector and SIQS-domain contract; no E1 gate for isolated work. Proceed when residual recovery/division is material; retain exact domain coverage. | Astra / xhigh |
| 10 | [A1a; C4/B7/B6/D2: Setup/transition reuse, then matrix work](https://github.com/nishanth17/factor/blob/cb37da3ea5fbbeb44923faa87fd45d4e5d926bc4/v2/ROADMAP.md#new-source-derived-deltas--separate-follow-up-items) | A1a/C10a-01 setup identity; A8 fixtures and B7 provenance contracts. Matrix work uses A8 contracts and genuine fixtures; capacity pressure can advance it independently. | Sol / xhigh metadata; Astra / xhigh arithmetic/proofs |
| 11 | [A5 → B8: Williams p+1, then starts and chains](https://github.com/nishanth17/factor/blob/cb37da3ea5fbbeb44923faa87fd45d4e5d926bc4/v2/ROADMAP.md#execution-phase-a--independent-foundations-that-can-proceed-now) | P2 recovery; B8 additionally needs A4 and A6 contracts. Establish the binary correctness control before optimized starts/chains. | Sol / xhigh binary; Astra / xhigh starts/chains |
| 12 | [F8: Specialized P−1/P+1 continuation](https://github.com/nishanth17/factor/blob/cb37da3ea5fbbeb44923faa87fd45d4e5d926bc4/v2/ROADMAP.md#execution-phase-f--independent-p6-challengers) | A6 control and polynomial contracts; E1 only for combined acceptance. Current PM1/PP1 costs must justify it; A5 is required only for PP1. | Astra / xhigh |


### Other open and conditional items

These retain their existing dependencies and triggers; they are not an extra ranked campaign.
Links preserve the detailed task contracts without duplicating them here.

| Item | Status / next deliverable | Dependencies | Model / effort |
| --- | --- | --- | --- |
| [B3b](https://github.com/nishanth17/factor/blob/cb37da3ea5fbbeb44923faa87fd45d4e5d926bc4/v2/ROADMAP.md#execution-phase-b--exploit-the-settled-controls) | Open: Use direct tuples for ECM instructions; compare preparation, execution, strict recovery and memory. | Current B3/C6 controls; coordinate the representation with B4a | Sol / xhigh |
| [B5](https://github.com/nishanth17/factor/blob/cb37da3ea5fbbeb44923faa87fd45d4e5d926bc4/v2/ROADMAP.md#execution-phase-b--exploit-the-settled-controls) | Conditional: Screen one NumPy bulk-kernel alternative against B5a/A1c on the actual PyPy toolchain. | Current collector cost share and exact A1c oracle; select B5 or B5a from boundary/toolchain evidence, not both by default | Sol / xhigh |
| [B13](https://github.com/nishanth17/factor/blob/cb37da3ea5fbbeb44923faa87fd45d4e5d926bc4/v2/ROADMAP.md#execution-phase-b--exploit-the-settled-controls) | Conditional: Revisit wide-block resieve capacity only after a new relevant refusal. | A1; a new calibrated wide-block capacity refusal | Astra / xhigh |
| [B14](https://github.com/nishanth17/factor/blob/cb37da3ea5fbbeb44923faa87fd45d4e5d926bc4/v2/ROADMAP.md#execution-phase-b--exploit-the-settled-controls) | Open: Add bounded n−1/Pocklington certificates and an independent exact verifier. | A10; existing shared-budget/checkpoint controls | Astra / xhigh |
| [C5](https://github.com/nishanth17/factor/blob/cb37da3ea5fbbeb44923faa87fd45d4e5d926bc4/v2/ROADMAP.md#execution-phase-c--measured-current-engine-challengers) | Conditional: Select a native limb-buffer or modular-context improvement beyond B4a reusable mpz. | A2/B4 oracles; B4a scalar native baseline and useful current allocation, dispatch or modular-arithmetic headroom | Astra / xhigh |
| [C8](https://github.com/nishanth17/factor/blob/cb37da3ea5fbbeb44923faa87fd45d4e5d926bc4/v2/ROADMAP.md#execution-phase-c--measured-current-engine-challengers) | Conditional: Assess family-wide CRT/hypercube scheduling after eligible-prime and occupancy checks. | A1/B1 and a concrete family/root or eligible-prime scanning hypothesis | Astra / xhigh |
| [C10a](https://github.com/nishanth17/factor/blob/cb37da3ea5fbbeb44923faa87fd45d4e5d926bc4/v2/ROADMAP.md#execution-phase-c--measured-current-engine-challengers) | Open: Remaining collection population, representation and scheduling deltas; C10a-01 is Done. | Closed C10 tranche; current committed control and explicit shared policy/source identities | Sol / xhigh existing policy; Astra for new proofs |
| [A1c](https://github.com/nishanth17/factor/blob/cb37da3ea5fbbeb44923faa87fd45d4e5d926bc4/v2/ROADMAP.md#new-source-derived-deltas--separate-follow-up-items) | Open: Select a bulk mark/scan or same-modulus two-root kernel; keep A1d scheduling separate. | Settled A1a/queue controls, exact kernel oracles and a current marking/scanning cost hypothesis | Sol / xhigh |
| [A2a](https://github.com/nishanth17/factor/blob/cb37da3ea5fbbeb44923faa87fd45d4e5d926bc4/v2/ROADMAP.md#new-source-derived-deltas--separate-follow-up-items) | Conditional: Select one measured inverse, GCD or quotient helper at a useful batch boundary. | Completed A2/B4; C11 profile and a changed workload/cost hypothesis | Sol / xhigh |
| [A10a](https://github.com/nishanth17/factor/blob/cb37da3ea5fbbeb44923faa87fd45d4e5d926bc4/v2/ROADMAP.md#new-source-derived-deltas--separate-follow-up-items) | Conditional: Evaluate measured classification or witness-batch improvements with exact factor recovery. | Completed A10; E3 frozen control or a justified isolated oracle spike; classification-cost evidence | Sol / xhigh |
| [A11a](https://github.com/nishanth17/factor/blob/cb37da3ea5fbbeb44923faa87fd45d4e5d926bc4/v2/ROADMAP.md#new-source-derived-deltas--separate-follow-up-items) | Conditional: Separately select factor hints, batch input, browser progress or save/continue features. | Completed A11; explicit product selection, outside the optimization critical path | Sol / high |
| [C1a](https://github.com/nishanth17/factor/blob/cb37da3ea5fbbeb44923faa87fd45d4e5d926bc4/v2/ROADMAP.md#new-source-derived-deltas--separate-follow-up-items) | Conditional: Compare cheaper SLP anchor reuse only under a changed cost hypothesis. | Completed C1; C9 decision; a new cost/retention hypothesis | Astra / xhigh |
| [D7](https://github.com/nishanth17/factor/blob/cb37da3ea5fbbeb44923faa87fd45d4e5d926bc4/v2/ROADMAP.md#execution-phase-e--current-engine-acceptance-control-freeze-and-p8-experiments) | Open: Freeze the current portfolio and independent training/confirmation controls. | E1; existing SIQS/P2 controls | Sol / high |
| [E3](https://github.com/nishanth17/factor/blob/cb37da3ea5fbbeb44923faa87fd45d4e5d926bc4/v2/ROADMAP.md#execution-phase-e--current-engine-acceptance-control-freeze-and-p8-experiments) | Conditional: Select a measured preprocessing or classification improvement. | Current committed preprocessing/classification oracle and observed cheap-path opportunity for isolated design/feasibility; D7 for integrated portfolio tuning | Sol / xhigh |
| [E4](https://github.com/nishanth17/factor/blob/cb37da3ea5fbbeb44923faa87fd45d4e5d926bc4/v2/ROADMAP.md#execution-phase-e--current-engine-acceptance-control-freeze-and-p8-experiments) | Open: Screen rho batching, restart effort or a justified whole-batch native kernel. | Current committed rho/recovery control for isolated screening; D7 for integrated portfolio tuning | Sol / high |
| [E5](https://github.com/nishanth17/factor/blob/cb37da3ea5fbbeb44923faa87fd45d4e5d926bc4/v2/ROADMAP.md#execution-phase-e--current-engine-acceptance-control-freeze-and-p8-experiments) | Conditional: Remove measured polling, recovery, checkpoint or durable-I/O overhead. | Current committed boundary/accounting contracts and a measured remaining cost for isolated feasibility; D7 for integrated portfolio tuning | Sol / xhigh |
| [E5a](https://github.com/nishanth17/factor/blob/cb37da3ea5fbbeb44923faa87fd45d4e5d926bc4/v2/ROADMAP.md#execution-phase-e--current-engine-acceptance-control-freeze-and-p8-experiments) | Open: Select a measured serialization, encoding or copying improvement across engine boundaries. | E5/C11 measured serialization cost; coordinate B12/F4, A7a, B7/C9a and E6 ownership | Sol / xhigh |
| [F5](https://github.com/nishanth17/factor/blob/cb37da3ea5fbbeb44923faa87fd45d4e5d926bc4/v2/ROADMAP.md#execution-phase-f--independent-p6-challengers) | Conditional: Evaluate cross-family A=A0*q reuse with exact external-factor handling. | E1; B1 calibrated QS control | Astra / xhigh |
| [F6](https://github.com/nishanth17/factor/blob/cb37da3ea5fbbeb44923faa87fd45d4e5d926bc4/v2/ROADMAP.md#execution-phase-f--independent-p6-challengers) | Conditional: Design and compare higher-arity large-prime relations with exact provenance and recovery. | C1/current DLP control and an explicit higher-arity domain/capacity design; E1 is unnecessary for isolated feasibility | Astra / xhigh |
| [F9](https://github.com/nishanth17/factor/blob/cb37da3ea5fbbeb44923faa87fd45d4e5d926bc4/v2/ROADMAP.md#execution-phase-f--independent-p6-challengers) | Conditional: Screen decomposable genus-two/Kummer factoring only if separately selected. | A separately selected larger-factor opportunity after cheaper B3a/F2/F3 controls; independent oracles and a plausible local cost model | Astra / xhigh |
| [G1](https://github.com/nishanth17/factor/blob/cb37da3ea5fbbeb44923faa87fd45d4e5d926bc4/v2/ROADMAP.md#execution-phase-g--integrated-portfolio-allocation) | Open: Tune joint effort allocation, recursive credit and relation-engine handoff. | D7, E1; E3–E6 and selected F2–F7 decisions | Sol / xhigh |
| [H1](https://github.com/nishanth17/factor/blob/cb37da3ea5fbbeb44923faa87fd45d4e5d926bc4/v2/ROADMAP.md#execution-phase-h--final-phase-8-reconciliation) | Open: Integrate selected winners and confirm the combined portfolio on untouched inputs. | G1, E3, E4, E5, E6 decisions; selected P6 implementations | Sol / xhigh |
| [I1](https://github.com/nishanth17/factor/blob/cb37da3ea5fbbeb44923faa87fd45d4e5d926bc4/v2/ROADMAP.md#execution-phase-i--reproducible-publication) | Open: Publish reproducible comparative results after integration and confirmation. | H1; selected P6 decisions and committed-only validation; F1 only for later GNFS claims | Sol / xhigh |
| [C10a-02](https://github.com/nishanth17/factor/blob/cb37da3ea5fbbeb44923faa87fd45d4e5d926bc4/v2/ROADMAP.md#collection-subtranches-and-retained-follow-ups) | Conditional: a distinct product/diversity or population hypothesis. | Current C10a-01 results and useful post-filter evidence. | Sol / xhigh |
| [C10a-03](https://github.com/nishanth17/factor/blob/cb37da3ea5fbbeb44923faa87fd45d4e5d926bc4/v2/ROADMAP.md#collection-subtranches-and-retained-follow-ups) | Conditional: exact pre-norm screening or bounded threshold feedback. | Current collector, rejected-position coverage and exact refinement. | Sol / xhigh |
| [C10a-06](https://github.com/nishanth17/factor/blob/cb37da3ea5fbbeb44923faa87fd45d4e5d926bc4/v2/ROADMAP.md#collection-subtranches-and-retained-follow-ups) | Open: reconcile the separate historical runner/accounting path. | C10/C10a-05; preserve frozen sources and failed receipts. | Sol / xhigh |
| [C10a-04/G1](https://github.com/nishanth17/factor/blob/cb37da3ea5fbbeb44923faa87fd45d4e5d926bc4/v2/ROADMAP.md#collection-subtranches-and-retained-follow-ups) | Conditional: larger-size portfolio feasibility and eligibility. | A7a worker routing is adopted; new 75/80 probes wait until after B5a. | Sol / xhigh |
| [A1a](https://github.com/nishanth17/factor/blob/cb37da3ea5fbbeb44923faa87fd45d4e5d926bc4/v2/ROADMAP.md#collection-subtranches-and-retained-follow-ups) | Conditional: a distinct direct comparison under a new hypothesis. | Delivered A1a/A7a remain closed; preserve scoped extension losses. | Sol / high |
| [F3a](https://github.com/nishanth17/factor/blob/cb37da3ea5fbbeb44923faa87fd45d4e5d926bc4/v2/ROADMAP.md#f3a--polynomial-representation-accounting-and-native-backend-headroom) | Conditional: polynomial representation, accounting and backend improvements | F3 contracts and measured complete-boundary costs | F3 assignment |
| [P3.8-R1/R2/R3/R4](https://github.com/nishanth17/factor/blob/cb37da3ea5fbbeb44923faa87fd45d4e5d926bc4/v2/ROADMAP.md#p38-carry-over-and-mathematical-boundaries) | Remaining base calibration, queue scheduling, provenance and classification reuse | C10a-02/G1; A1d/B13/C8; A8/B7/C4; C9b/E3 respectively. Delivered tranches stay closed. | Existing arithmetic and proof owners |

### Deferred work

| Item | Restart condition |
| --- | --- |
| E1 / P3.8-R5 | Explicit authorization for combined SSS/worker/current-control acceptance. |
| A1b-F1 | A new credible fixed-score queue hypothesis; the completed existing-option loss stays closed. |
| P4.3 broader reducer work | A concrete changed workload/domain hypothesis; retain scoped earlier losses. |
| 75/80-digit QS probes | After B5a, under a separately selected finite probe plan. |
| [A9/B9/B10/B11/C7/D3/D4/D5/D6/E2/F1 — GNFS](https://github.com/nishanth17/factor/blob/cb37da3ea5fbbeb44923faa87fd45d4e5d926bc4/v2/ROADMAP.md#later-tranche--gnfs-after-current-engine-optimization) | H1 and the coverage review must activate GNFS; original stage dependencies still apply. |

## In progress

B12/F4 and B5a entries refreshed from their session handoffs on 9 November
2026. Other workstreams retain their recorded inspection dates.

Implementation snapshot checked on 6 November 2026 against mainline
`bf243dba2edef25c74812887f45473036b4db845`. The remaining workstreams below are
unfinished tasks with isolated implementations; none of their task code is
integrated by this documentation update. A8 and D1 share one worktree but keep
their separate contracts. This is a checkpoint inventory, not a claim that a
process is currently running or that an experimental acceptance gate passed.
In progress remains acceptance-open; phase-catalogue gates below still apply.

A7a is Done. Current preparation positions retain B4a (3), F3 (11), F2 (14), and
A8/D1 within the matrix branch (16). Adding H2 at rank 7 shifts the previous
rank-7-and-later items by one; their relative order and task scopes are retained.
Detailed local entry points below are in each named task worktree, not mainline.

**B12 → F4 has an unfinished ready-first startup follow-up.** Its earlier
process-ECM delivery and selected defaults remain recorded in the
[completed register](#completed-register), but do not validate this new
prototype. Native ECM continues under B4a; the selected current-collector
QS worker integration is recorded under completed A7a. Broader portfolio, balanced completion and population-p90
evidence remain separate acceptance follow-ups.

| ID / scope | Deliverable and boundary | Predecessors | Status, trigger and existing priority | Model / effort | Evidence and research |
| --- | --- | --- | --- | --- | --- |
| B12/F4 — ready-first worker startup follow-up | Dispatch/publication while peers start, with checked readiness, H3 observed root accounting, original startup deadlines and partial-pool cleanup. Barrier remains default. | Delivered B12/F4 ECM core; current flat ECM package and unchanged curve manifests | **In progress; H3-adapted prototype.** Seventeen mock-only checks and static checks pass. Real-spawn/arithmetic, lifecycle/reuse tests, committed-only full QA, matched startup costs and engagement-table recalibration remain. | Sol / **xhigh** validation/design; separately frozen execution can use Sol / **high** | [implemented and remaining work](#b12-f4-in-progress) |
| B5a — native QS packet and ARM SIMD | Implement an opt-in persistent exact mark/queue/scan packet with shared buffers; compare scalar C, compiler vectorization and NEON on one checked ABI. | Current Y0 and A1a/A1c/A1d oracles; reviewed B4a native infrastructure and H3 root contracts | **In progress; Phase 1 review complete.** Source audit, research, opportunity ledger and staged protocol exist. Native implementation, exactness/lifecycle admission, complete packet costs, integration and confirmation remain; no native speedup is established. | Sol / **xhigh** boundary/design/implementation; separate arithmetic proofs retain Astra / **xhigh** | [implemented and remaining work](#b5a-in-progress) |
| F3 — P6.2 exact polynomial ECM continuation | Select bounded product/remainder-tree and multipoint continuation from GMP-ECM and Yamaquasi. Compare schoolbook/Karatsuba and exact Kronecker packing on native ints, with optional coarse GMP multiplication as a separate track. Include monic reciprocal reduction, blocked/tree storage, exact evaluation/recovery and a separately flagged Brent–Suyama power/Dickson extension; see the current gates below. | Current B2/classical control and C11/C3a bound-specific opportunity; selected polynomial oracles/contracts; E1 only for combined acceptance; C2 win unnecessary | **In progress; opt-in polynomial continuation.** Exact degree-one continuation, P/G arithmetic, recovery/resume and explicit ablations are committed. Full QA, the unexecuted prospective screen, whole-continuation crossover and integration remain pending. | Astra / **xhigh** research, implementation, exactness and experiment design; optional separate Sol / **high** frozen execution | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-f--independent-p6-challengers) ; [implementation and remaining work](#f3-in-progress) |
| A8 — P3.8 matrix/operator/lifting control | Freeze genuine original-relation and post-filter fixtures, orientation, lifting maps and packed forward/transpose products, including dense-prefix/sparse-tail and selected locality layouts. Independently check dense-oracle products, repeated-coordinate cancellation, short tails, bilinear identities and original-relation extraction; include preparation ownership and possible dense fill. | Accepted R3 identities/lifting and P3.3 control; fixtures refresh after collector changes | **In progress with D1.** Exact operator/lifting adapters, independent oracles and fixtures are implemented. Current nonempty target-size fixtures, heavy QA and capacity measurements remain open; no matrix-memory improvement is established. | Astra / **xhigh** design, implementation and correctness; optional separate Sol / **high** frozen fixture measurements | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-a--independent-foundations-that-can-proceed-now) ; [implementation and remaining work](#a8-d1-in-progress) |
| D1 — P3.8 packed block Lanczos | Implement one genuine selected-block recurrence with packed forward/transpose products, singular-block handling, original-operator kernel correction, lifting, bounded seeded restarts and checked resume. Consider Thome's reduced-vector variant at recurrence selection; independently verify every original M d = 0. | A8 operator/lifting contract; use B6/C4 as comparison controls, not required successes | **In progress after the A8 contract.** Genuine selected-block recurrence, original-kernel correction, bounded restarts/resume and default-off ablations are implemented. Complete QA, representative matrix measurements, confirmation and integration remain open. | Astra / **xhigh** design, implementation and proof; **max** only for a named unresolved invariant; optional separate Sol / **high** frozen execution | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-d--sparse-solver-alternatives) ; [implementation and remaining work](#a8-d1-in-progress) |
| F2 — P6.1 torsion-aware curve/formula packages | Compare complete packages against Suyama: GMP-ECM small-parameter/batch formulas (MG17), CADO Montgomery12/16 and mixed Montgomery/twisted-Edwards bytecode (MG18), Edwards/windowed or sliding-window alternatives, Suyama-11 and the calculator reference family (AP30). Screen JKL twisted-Hessian/small-parameter families with their field-dependent torsion conditions before selecting a package. Cost curve generation, chain setup, nonunit maps and continuation conversion together. | A2/B3/B4 retained arithmetic/oracle contracts; C11 current factor-size/cost hypothesis; E1 only for later combined acceptance | **In progress; opt-in package preparation.** Suyama-11, enhanced Edwards/setup/scalar options and scheduled execution/plan reuse are committed. Further optimization edits are uncommitted at this snapshot; heavy QA, feasibility, yield/cost evidence and integration remain pending. | Astra / **xhigh** research, implementation and correctness; optional separate Sol / **high** frozen execution | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-f--independent-p6-challengers); [JKL candidate](https://doi.org/10.1112/S1461157016000231) ; [implementation and remaining work](#f2-in-progress) |
| B4a — P4.2 fused native ECM chunk boundary | Select one unchanged bounded ECM ladder/verified-program native executor under the [coarse-boundary contract](#b4a--coarse-native-arithmetic-and-whole-portfolio-costs). Keep canonical state, private scratch and the Python fallback. | Existing B3/B4 exact oracles and current committed serial identity; a measured Python arithmetic-loop share; funded native memory/work and cancellation/restore contract | **Delivered within the explicit native scope.** Coarse observed C segments, native metadata/reuse and lazy imports pass committed-only QA and 1,512 fresh affected-ECM calls. All three 60-digit integration pairs and one conditional 70-digit pair pass canonical-work, quiescence and cleanup checks. Deep retained 60-digit performance and broader population/p90 promotion remain unresolved. Historical losses and charges are preserved; no automatic adoption. | Sol / **xhigh** unchanged arithmetic, boundary design and implementation; Astra / **xhigh** only for a separately selected new encoded-domain/formula proof; optional separate Sol / **high** frozen execution | [current native-boundary gates](#b4a--coarse-native-arithmetic-and-whole-portfolio-costs); [CFFI API mode](https://cffi.readthedocs.io/en/latest/overview.html#purely-for-performance-api-level-out-of-line) ; [implementation and remaining work](#b4a-in-progress) |

<a id="b12-f4-in-progress"></a>

### B12 / F4 — ready-first startup follow-up

- **Implemented on the task branch:** dispatch/publication while peers start,
  checked readiness, original startup deadlines and partial-pool cleanup.
  The prototype now uses H3 root process tracking and startup/result CPU
  acknowledgments, without prepaid worker/work/CPU/memory grants. Barrier
  remains the default; curve manifests and arithmetic are unchanged.
  Seventeen mock-only checks and static checks pass; real workers are unvalidated.
- **Remaining:** reconcile with current mainline; real-spawn early-publication
  and cleanup checks; identical-manifest serial/P1/P2/P4 arithmetic; persistent
  reuse and worker failure/restart regressions; committed-files-only validation,
  full suite, lint and benchmark imports. The prepared focused launcher does
  not replace the full suite. Then compare matched barrier/ready-first cold
  startup, first dispatch/validated factor and total return/drain/close costs.
  Recalibrate startup, engagement and worker cost tables on the H3-combined
  source before adoption. No startup speedup is established.
- **Checkpoint:** `codex/ecm-ready-start` at `dd37423`, with local master
  `755c312` integrated. Worktree:
  `/Users/nishanthmohan/Documents/code/factor-ecm-ready-start`.
  Latest local entry: `v2/audit/ecm_ready_start/H3_MERGE.md`; older
  `REBASE.md` and `HANDOFF.md` retain their historical source identities.
  The delivered B12/F4 ECM core stays Done. The completed tuple-decoding
  microcheck and open B3b direct-tuple implementation are separate scopes.

<a id="b5a-in-progress"></a>
<a id="bulk-kernels-pypy-numpy-and-optional-native-acceleration"></a>

### B5a — native QS packet and ARM SIMD

- **Completed preparation:** current Y0 execution/representation audit,
  implementation/literature review, opportunity ledger, ranked boundary design
  and staged protocol. Short current-source prefixes support prioritization;
  no native QS implementation or speedup is established. The msieve/YAFU
  reference-rate comparison remains unmeasured after build failures and a
  busy-lock repair gate.
- **Remaining implementation:** build the opt-in scalar persistent exact
  mark/queue/scan packet with reusable score/event/hit/scratch buffers and the
  existing sparse queue. Reuse the reviewed B4a build/provenance/ABI/lifetime
  infrastructure. Preserve Y0, exact score widths and hit masks, H3 root
  accounting, action-quantum polling, cancellation and dirty-block replay.
  Compare scalar C, compiler-vectorized C and NEON through the same ABI.
- **Remaining validation and decisions:** Q0 exactness/build/lifecycle admission;
  Q1 clean feasibility and baselines; Q2 prepared-kernel screen; Q3 complete
  packet costs including imports, conversions, output and polling. Advance to
  Q4 serial and Q5/Q5O process/ordered-parent integration only through their
  gates, then applicable full-factor feasibility and fresh confirmation.
  Candidate-only hit reconstruction versus a full native sidecar is an
  explicit comparison. Persistent root shifts and a whole native worker remain
  conditional on the first packet's measured costs. Reference header repair
  is a separate pending comparison task. No later campaign is authorized by
  this roadmap entry; new 75/80 probes remain deferred until after B5a.
- **Checkpoint:** `codex/b5a-native-qs`, Phase 1 documentation commit `f238b3f`,
  production base `cb37da3`. Worktree:
  `/Users/nishanthmohan/Documents/code/factor-b5a-native`.
  Local entry: `v2/audit/b5a/opportunity_review_20261109/PHASE1_HANDOFF.md`;
  `OPTIMIZATION_PLAN.md` and `EXPERIMENT_PLAN.md` hold the implementation and
  gated execution details. Production source is unchanged by Phase 1.

<a id="f3-in-progress"></a>

### F3 — exact polynomial ECM continuation

Target useful 20–35-digit smaller factors within larger inputs; production
remains classical continuation until an explicit integration/adoption decision.
Establish a complete-continuation crossover before enlarging B2 or claiming
that a packed multiplication backend pays off.

- **Implemented on the task branch:** exact nonnegative Kronecker packing,
  monic reciprocal reduction, product/remainder trees, degree-one wheel
  continuation, checked normalization, bounded saturation recovery and
  source/backend-bound block resume. Python-int P and optional coarse-gmpy2 G
  paths are available in the task environment. Default-off arithmetic,
  packing/reuse and redundant-work profiles have focused exactness checks;
  the latest handoff reports 63 focused tests, including P-only validation.
- **Remaining:** committed-only full QA and runner lifecycle checks; execute
  one selected source-pinned screen, not every historical panel; measure
  complete construction/evaluation/recovery cost and establish the crossover
  at useful deeper bounds. Freeze later independent comparisons from those
  results; reconcile B4a costs before production integration or promotion.
  Power/Dickson extension implementation and the F3a representation/accounting/
  scalar-backend candidates remain separate, unimplemented follow-ups.
- **Checkpoint:** `codex/f3-polynomial-stage2` at `66c07f3`.
  Local entry: `v2/audit/f3/CURRENT.md`, then
  `HEADROOM_OPTIONS_HANDOFF_20261106.md`. Its replacement 36-cell screen is
  unexecuted and replaces an earlier panel within the existing allowance;
  it is not an additional grant or performance result.

<a id="a8-d1-in-progress"></a>

### A8 / D1 — matrix controls and packed block Lanczos

Capacity can justify sparse work before solve CPU dominates. Synthetic
operator gains or an unrelated control's failure do not settle adoption.

- **Implemented on the task branch:** constraint-by-relation operators,
  packed forward/transpose products, independent oracles, exact fixture and
  filtering/lifting adapters, and a genuine Thomé selected-block recurrence.
  Original-kernel correction, finite seeded restarts, transactional resource
  handling, checked resume and exact extraction adapters are implemented.
  Default-off kernel/layout, joint terminal recovery and further recurrence/
  recovery ablations have focused correctness checks. The solver returns some
  verified independent dependencies, not a complete basis or proved nullity.
- **Remaining:** current genuine nonempty target-size fixtures, particularly
  60–80-digit matrices; committed-only full QA and runner/resource checks;
  matched complete preparation/solve/correction/lift/extraction measurements,
  including restart costs, useful dependencies and peak memory. The existing
  dense-fill/provenance reservation is still a pipeline constraint: no
  capacity improvement is established. Integration and adoption are pending.
  Compact provenance, stronger filtering, Four Russians and a second sparse
  solver retain B7/C4/B6/D2 ownership; they are not delivered by this prototype.
- **Checkpoint:** `codex/a8-d1-lanczos` at `9dd1968`.
  Local entries: `v2/audit/a8_d1/HANDOFF.md` and
  `implementation-i5.json`; retained historical fixture replays and small
  exact tests do not replace target-size performance evidence.

<a id="f2-in-progress"></a>

### F2 — torsion-aware curve and formula packages

Current Suyama already has favorable divisibility. Compare complete cost and
chosen-point behavior, distinguish point order from curve order, and preserve
finite family exhaustion and nonunit/close-factor recovery. Higher torsion
does not establish a free increase in B1.

- **Implemented on the task branch:** a Suyama-11 chosen-point package with
  full twist/point checks, bounded deterministic sampling, duplicate accounting,
  saturation recovery, attempt-boundary resume and checked Montgomery
  continuation. Enhanced setup/auxiliary-reuse, extended-Edwards stage one,
  signed-scalar options, coordinate scheduling, cached addition terms and
  bounded scalar-plan reuse are committed with independent tiny oracles and
  explicit experiment preparations.
- **Remaining:** finish and commit the in-flight optimized setup/quotient/
  sparse-execution edits before treating them as a validated checkpoint;
  perform full QA and runner lifecycle checks; establish an uninstrumented
  complete-call feasibility baseline, then measure setup, stage one, conversion
  and successful/unsuccessful searches under matched resources. Family-yield
  inference needs independent moduli and residue-class coverage, not repeated
  timings or tiny-field order checks. Fresh confirmation, backend recalibration
  and integration remain open. Alternative families/formulas and identified
  preparation/plan-capacity headroom remain conditional research candidates.
- **Checkpoint:** `codex/f2-curve-packages` at `794e1d5`, with additional
  uncommitted implementation edits observed. Local entries:
  `v2/audit/f2/focused/README.md` and `v2/audit/f2/headroom/README.md`.
  Existing prospective panels are unexecuted; no package/yield win is established.

<a id="b4a-in-progress"></a>

### B4a — optional coarse native ECM

- **Post-integration tuning delivered:** from `2b46e3f`, analyze raw and
  refreshed actual-production CPU/heap/JIT evidence and test thirteen B4a
  mechanisms. Keep signed submul, private numerator, supported GMP division
  entry and checked resource-query transport. The final 648-call cold/warm
  panel, 90 recovery and 96 cheap cells pass with explicit uncertainty.
  Typed ownership graphs and full allocation histories separate retention,
  churn and RSS; unresolved JIT/GC object churn remains qualified. No new
  representation/reducer, prime algorithm, routing or automatic default.
  Current delivery, costs, resources and the revised unexecuted 3F plan:
  `v2/audit/b4a/post_integration_tuning_20261110/HANDOFF.md`.

- **3P delivered:** explicit scalar/pthread/OpenMP production owner, portfolio,
  publication and checked checkpoint/replay paths. H3 repair `98bf836` and
  native-unit adapter `4f30d87` are included. OpenMP remains an explicit optional
  LLVM/libomp dependency; its compatibility benchmark import re-exports the
  production owner.
- **Current recipe:** generic O2/ordinary GMP/coalesced, scalar 256/team 128,
  pthread four for the shallow nine-curve anchor and eight for sustained
  deep 37-curve work, dynamic/chunk 1/capacity 32. Whole-owner support through
  sixteen is validated; the 216-call count screen supports no larger count.
  Counts and ownership remain explicit. No automatic numerical/routing change.
- **Validated:** implementation `0c7c327` passes 1,778 committed-only tests
  (14 optional skips) with all three native builds. Formatting-only `fb52b8c`
  has identical test ASTs and passes full lint/format and 311 Git-listed
  ordinary imports. Actual teams, failures, limits, cancellation, paid replay
  and cleanup pass. The exact observer remains enabled; a bounded inherited-
  lock preflight holds the performance window through all final stage cleanup.
- **Boundary implementation closeout:** from `410eba0`, establish current removable
  critical-path ceilings against matched A/A for preparation/registration,
  generation/refill/settlement and conversion/export. Implement material exact
  bulk alternatives or document applicable retain-path evidence. Coarse
  resource-observed C coordination, lazy relation imports and reused accepted
  catalogue bytes are implemented through `6a34003`. Later refill/import
  challengers are admitted and tested without useful improvement. Source-labelled
  individual and combined ceilings retain model/systematic limits. Committed-only
  `3bb5a69` passes 1,805 tests (14 optional skips), lint, 317 ordinary imports
  and all three reproducible native builds. Final fresh confirmation is delivered
  below. Existing whole-classical C/GMP arithmetic and previous wins/losses
  remain controls.
  Current request: `v2/audit/b4a/native_boundary_closeout_20261111/REQUEST.md`.
  Persisted cache data is project-local, checked and included in the finite
  capture allowance. An unresolved material path prevents closeout.
- **Final 3F delivered:** 1,512 fresh affected calls cover 6/9/6 independent
  inputs at 60/80/100 digits, both anchors and cold/retained ownership. Zero
  censors; pooled primary contrasts support scoped gains, while 60-digit deep
  retained performance remains inconclusive under its frozen stop rule.
  Recovery/cheap guards pass. All three 60-digit integration pairs precede one
  clean 70-digit pair; canonical QS work agrees, native workers stay parked
  with zero observed CPU advance, and final ownership reaches zero. The slower
  native 70-digit total is retained as descriptive QS-dominated evidence;
  no additional integration question justifies more calls. Population/p90 and
  broad balanced promotion remain separate unmet gates. Current tables,
  limitations, readiness and exact separate charges:
  `v2/audit/b4a/confirmation_3f_20261111/3F_HANDOFF.md`.
- **Preserved:** earlier source freezes, failures, charges, compiler/expression
  decisions and scoped SME/NEON losses. C5 arithmetic/SIMD, B3 layouts,
  C2 continuation, E6 prime generation and C3a/C11 routing retain their owners.

<a id="qs-optimization-priorities--coverage-reconciliation-2-november-2026"></a>
<a id="model-and-effort-for-the-new-qs-work"></a>

<a id="completed-register"></a>

## Completed

One entry per delivered scope. Detailed measurements, failures and adoption limits are in
[RESULTS.md](benchmarks/RESULTS.md); older records remain in the
[pinned register](https://github.com/nishanth17/factor/blob/cb37da3ea5fbbeb44923faa87fd45d4e5d926bc4/v2/ROADMAP.md#completed-register). Follow-ups retain separate IDs above.

| Scope | Delivered outcome |
| --- | --- |
| <a id="a7a-completion"></a><a id="a7a-in-progress"></a>A7a | Current-collector process QS, Y0 and selected routing integrated; serial override retained. |
| <a id="h3--simplify-resource-limits"></a>H3 | Observed root resources, allocator-grown internal state and checked saved-history resume. |
| B12 → F4 / P6.3 ECM | Reproducible process ECM and selected cost-gated routing; ready-first remains separate. |
| A1b existing fixed-score review | Retain defaults; scoped challenger lost. |
| P1.1–P1.9 | Core exact arithmetic, API, prime-sieve and recovery repairs. |
| P2 core / M13 / M14 / M17 | Bounded core and loop verification. |
| P3.1 | Polynomial and relation identities. |
| P3.2 | Single-large-prime collection. |
| P3.3 / M28 | Filtering, lifting and exact extraction. |
| P3.4 / M31 | SIQS and bounded polynomial evaluation. |
| A11 / P3.4 CLI | Explicit QS/MPQS/SIQS CLI access. |
| P3.5 | SSS/SSSf adapters; remain opt-in. |
| P3.6 | Historical bounded-worker study; retained serial decision. |
| P3.6.1 | SSS/worker diagnosis and repairs. |
| P3.8-R1 implementation | Streamed assignments and external-square implementation. |
| P3.8-R3 bounded tranche | Mixed identities/recovery; scoped cadence option remains opt-in. |
| A1 / P3.8-R2 | Exact sieve/resieve repair and scope decisions. |
| A2 / P4.3 | Scoped reducer comparison; retained baseline. |
| A3 / P5.2 programs | Reusable stage-two programs remain opt-in; streamed default retained. |
| A4 / P4.1 correctness | Checked scalar-chain correctness control. |
| A6 / P5.3 p−1 | Selected p−1 continuation and schedule control. |
| A7 / P3.8-R5 reconciliation | Historical worker/SSS scope reconciliation. |
| A10 / P8.2.1 | Classification and primality transfer. |
| B1 / P3.8-R1 calibration | Bounded SIQS base/polynomial calibration. |
| B2 / P5.2 pairing | Scoped pairing comparison and retained controls. |
| B3 / P4.1 production | Verified chain execution at selected bounds. |
| B4 / P4.2 | Scoped fused/native bridge comparison; readable control retained. |
| C6 / P4.5 | Reused PRAC/Lucas comparison; fresh native ladder retained in that study. |
| C1 / P3.8-R4 / P5.4 | Scoped DLP calibration; original SLP defaults retained for that decision. |
| C3 bounded allocation | Original bounded allocation control; C3a/G1 remain separate. |
| C9 / C9.2 | Growing-DLP comparison with its higher memory cost recorded. |
| C10 | Source-derived collection strategy tranche. |
| C10a-05 | Selected collection pacing. |
| C10a-05 size-scope extension | Recorded larger-size pacing extension. |
| C10a-04 delivered routing | Selected compact counted-base routing. |
| A1a delivered integration/extension | Selected marking integration; scoped extension losses retained. |
| C10a-07 | Combined face/queue integration. |
| C11 first-pass audit | Parameter/source reconciliation; ongoing atlas remains separate. |
| <a id="c10a-01-completion"></a>C10a-01 | Coupled polynomial, multiplier, population and root policies. |
| <a id="c9a-completion"></a>C9a | Combined useful-yield/storage implementation; opt-in policies retained. |

<!-- Compatibility fragments: detailed contracts and history are in the pinned source above. -->
<a id="a1a--selected-collector-kernel-delta-delivered-31-october-2026"></a>
<a id="a1a-follow-ups--checked-transitions-and-lifted-roots"></a>
<a id="a1b--existing-fixed-score-promotion-follow-up-after-completed-a1-and-c10"></a>
<a id="a1c--same-modulus-two-root-interleaving-conditional"></a>
<a id="a7a--current-collector-parallel-relation-collection"></a>
<a id="a7a--done-current-collector-parallel-qs"></a>
<a id="acceptance-evidence-and-stop-rules"></a>
<a id="b3-production-integration-acceptance-9-october-2026"></a>
<a id="b4a--coarse-native-arithmetic-and-whole-portfolio-costs"></a>
<a id="broader-portfolio-opportunities-and-relative-value"></a>
<a id="c10--source-derived-collection-strategies-closed-bounded-tranche"></a>
<a id="c10a--conditional-collection-policy-questions-open"></a>
<a id="c10a-01--coupled-polynomial-multiplier-a-population-and-root-setup"></a>
<a id="c10a-01-implementation-remaining-after-the-coverage-review"></a>
<a id="c11--cross-engine-constant-provenance-and-calibration-open"></a>
<a id="c1a--cheaper-slp-anchor-delta-after-c9"></a>
<a id="c3a--source-derived-policy-delta-after-c3"></a>
<a id="c5--low-level-native-limb-buffers-and-modular-contexts"></a>
<a id="c5--native-montgomeryredc-residue-contexts"></a>
<a id="c9--source-derived-large-prime-yield-and-dlp-economics-complete-and-integrated"></a>
<a id="c9a--combined-useful-yield-and-storage-workstream"></a>
<a id="c9a--remaining-large-prime-yield-and-storage-transfers-open"></a>
<a id="c9b--mainline-validation-and-remaining-headroom"></a>
<a id="calibrate-automatic-selection-and-handoff-c3-then-g1"></a>
<a id="collection-subtranches-and-retained-follow-ups"></a>
<a id="coverage-decisions-from-the-expanded-source-review"></a>
<a id="current-foundation-already-available"></a>
<a id="dependency-overview"></a>
<a id="documentation-history-and-local-evidence"></a>
<a id="e6--self-contained-native-prime-generation"></a>
<a id="ecm-opportunity-and-continuation-gates"></a>
<a id="execution-order-optimize-current-engines-then-revisit-gnfs"></a>
<a id="execution-phase-a--independent-foundations-that-can-proceed-now"></a>
<a id="execution-phase-b--exploit-the-settled-controls"></a>
<a id="execution-phase-c--measured-current-engine-challengers"></a>
<a id="execution-phase-d--sparse-solver-alternatives"></a>
<a id="execution-phase-e--current-engine-acceptance-control-freeze-and-p8-experiments"></a>
<a id="execution-phase-f--independent-p6-challengers"></a>
<a id="execution-phase-g--integrated-portfolio-allocation"></a>
<a id="execution-phase-h--final-phase-8-reconciliation"></a>
<a id="execution-phase-i--reproducible-publication"></a>
<a id="existing-evidence-to-reuse"></a>
<a id="f3a--polynomial-representation-accounting-and-native-backend-headroom"></a>
<a id="factor-live-backlog"></a>
<a id="factor-phased-todos-and-decision-gates"></a>
<a id="historical-heading-anchors"></a>
<a id="how-to-read-the-phases"></a>
<a id="how-to-use-the-gates"></a>
<a id="later-integrated-preprocessing-tuning-e3"></a>
<a id="later-tranche--gnfs-after-current-engine-optimization"></a>
<a id="master-execution-sequence-and-model-recommendations"></a>
<a id="model-and-effort-key"></a>
<a id="new-source-derived-deltas--separate-follow-up-items"></a>
<a id="numbered-phase-catalogue-and-acceptance-gates"></a>
<a id="open-and-conditional-execution-items"></a>
<a id="optional-product-follow-ups--separate-from-factoring-optimization"></a>
<a id="optional-product-work"></a>
<a id="original-h3-acceptance-checklist-retained-local-draft"></a>
<a id="p11--port-the-seven-modules-to-explicit-python-3-arithmetic"></a>
<a id="p12--repair-suyama-ecm-initialization"></a>
<a id="p13--repair-and-validate-every-sieve-endpoint"></a>
<a id="p14--preserve-unresolved-cofactors-and-define-input-behavior"></a>
<a id="p15--repair-p1-batching-and-saturation-handling"></a>
<a id="p16--repair-rho-retries-and-cap-recovery"></a>
<a id="p17--repair-ecm-curve-limits-and-failure-transitions"></a>
<a id="p18--clarify-primality-and-helper-contracts"></a>
<a id="p19--quarantine-unsafe-prac-and-make-diagnostics-reliable"></a>
<a id="p21--add-exact-preprocessing-and-classification-reuse"></a>
<a id="p22--build-a-budgeted-portfolio-dispatcher"></a>
<a id="p23--stream-and-cache-prime-schedules-within-a-memory-cap"></a>
<a id="p24--add-chunk-checkpoints-and-a-consistent-recovery-protocol"></a>
<a id="p25--freeze-a-benchmark-corpus-and-runner"></a>
<a id="p26--port-reusable-sieve-contexts-and-rolling-strikes-selectively"></a>
<a id="p27--compare-wheelpre-sieve-and-packed-output-ideas-from-c-v3v4"></a>
<a id="p28--establish-whether-parallel-candidate-search-is-worthwhile"></a>
<a id="p31--implement-polynomial-and-relation-identities"></a>
<a id="p32--build-a-bounded-single-large-prime-collector"></a>
<a id="p33--implement-filtering-gf2-dependencies-and-factor-extraction"></a>
<a id="p34--add-siqs-self-initialization-and-bounded-dispatch"></a>
<a id="p34-usability-follow-up--explicit-qsmpqssiqs-cli-access-a11"></a>
<a id="p35--evaluate-smooth-subsum-search-on-the-same-interface"></a>
<a id="p36--evaluate-coarse-parallel-relationcandidate-collection"></a>
<a id="p361--immediate-diagnosis-and-improvement-of-p35p36"></a>
<a id="p37--optional-numpy-spike-for-bounded-siqs-sieve-arrays"></a>
<a id="p38--evaluate-gf2-matrix-and-filtering-optimizations"></a>
<a id="p38-carry-over-and-mathematical-boundaries"></a>
<a id="p41--finish-prac-repair-and-precompute-valid-chains"></a>
<a id="p42--compare-fused-and-normalized-montgomery-kernels"></a>
<a id="p43--compare-pypy-python-integers-and-optional-gmpmpz"></a>
<a id="p44--keep-barrettmontgomery-reducers-experimental-until-measured"></a>
<a id="p45--experiment-with-reference-ecm-precomputed-lucas-chains-c6"></a>
<a id="p51--implement-williams-p1-with-binary-lucas-evaluation"></a>
<a id="p52--pair-ecm-stage-2-primes-and-tune-table-size"></a>
<a id="p53--improve-p1-powering-and-continuation-independently"></a>
<a id="p54--add-double-large-prime-siqs-and-stronger-filtering"></a>
<a id="p61--compare-edwardswindowed-and-torsion-aware-ecm"></a>
<a id="p62--investigate-polynomial-continuations-and-richer-relation-collectors"></a>
<a id="p63--add-reproducible-process-level-parallelism"></a>
<a id="p64--publish-reproducible-workload-specific-results"></a>
<a id="p71--establish-exact-polynomial-and-relation-contracts"></a>
<a id="p72--collect-bounded-rationalalgebraic-relations"></a>
<a id="p73--filter-and-solve-dependencies-with-gnfs-constraints"></a>
<a id="p74--compute-rational-and-algebraic-square-roots-and-split-n"></a>
<a id="p75--integrate-the-small-gnfs-engine-with-bounded-dispatch"></a>
<a id="p76--scale-polynomial-selection-sieving-and-matrix-work"></a>
<a id="p77--measure-coverage-and-the-siqsgnfs-crossover"></a>
<a id="p81--freeze-the-new-control-and-untouched-evaluation-data"></a>
<a id="p82--reduce-preprocessing-work-with-exact-rejection-proofs"></a>
<a id="p821--transfer-verified-v1-advantages-starting-with-primality-a10"></a>
<a id="p822--generate-and-independently-verify-prime-certificates-b14"></a>
<a id="p83--calibrate-bounded-rho-batches-walks-and-restarts"></a>
<a id="p84--tune-ecmp1-allocation-and-relation-engine-handoff"></a>
<a id="p85--measure-recovery-polling-and-checkpoint-costs-separately"></a>
<a id="p86--revisit-context-setup-and-schedules-only-for-changed-workloads"></a>
<a id="p87--integrate-accepted-winners-and-close-the-reconciliation"></a>
<a id="parallel-execution-and-ownership"></a>
<a id="phase-1--immediate-correctness-and-api-repairs"></a>
<a id="phase-2--bounded-scheduling-reusable-infrastructure-and-benchmarks"></a>
<a id="phase-3--add-the-missing-balanced-composite-engine"></a>
<a id="phase-4--optimize-ecm-chains-and-arithmetic-after-coverage-exists"></a>
<a id="phase-5--complementary-methods-and-stronger-continuations"></a>
<a id="phase-6--research-options-parallel-execution-and-claims"></a>
<a id="phase-7--implement-and-scale-gnfs"></a>
<a id="phase-8--revisit-phase-2-portfolio-optimization"></a>
<a id="phase-checkpoints"></a>
<a id="ranked-next-work--expected-benefit-and-cost"></a>
<a id="sources-and-verification-scope"></a>
<a id="why-acceleration-leads-and-where-the-argument-stops"></a>
