# Factor live backlog

Start with [STATUS.md](STATUS.md) for current defaults, work in flight and next
eligible items. [POLICY.md](POLICY.md) holds stable rules. Updated 4 November
2026 from mainline `7a3ca71`, preserving pending A7a/B5a review additions.
B12/F4 remains ongoing; C9a and C10a-01 are delivered scopes.

Status vocabulary: **Done**, **Open**, **Conditional** (requires its trigger),
**Deferred**, **Rejected**, **Inconclusive**. Done describes the delivered scope;
it can retain an Inconclusive experimental verdict. Adoption, retention as
opt-in, and closure are recorded explicitly.

| ID scheme | Meaning and relationship |
| --- | --- |
| P1–P8, including P3.8-R1–R5 and P8.2.1/P8.2.2 | Capability/acceptance catalogue. Numbers preserve ownership and historical gates. |
| A–I execution layers, including suffixes | Dependency-based work packages across P phases. Each row maps to its P scope; layers are not a second list of capabilities. GNFS rows retain their letters while Deferred. |
| C10a-##, AP##, MG##, M## | Existing subtranches, source-transfer inventory, gap inventory and milestone evidence IDs. They retain their identities and linked historical meaning. |

The [ranked queue](#ranked-next-work--expected-benefit-and-cost) owns preparation
order; execution rows own scope, triggers, dependencies and model/effort.
[STATUS](STATUS.md) owns current defaults, ongoing work and the declared metric.
[POLICY](POLICY.md) owns acceptance, resources, stop rules and ownership.
Completed scopes remain in the [register](#completed-register); their detailed
records are preserved in results and pinned Git history. No run is scheduled
by this document.

<a id="qs-optimization-priorities--coverage-reconciliation-2-november-2026"></a>
<a id="model-and-effort-for-the-new-qs-work"></a>

## Ranked next work — expected benefit and cost

This is the authoritative preparation queue as of 4 November under the
[declared workload and wall-time objective](STATUS.md). C10a-01 historical freezes remain immutable; new comparisons use current
committed code. The short C3a/C11 screen comes first for broad routing value
at low initial cost; ECM workers are the next major implementation, followed
by current-collector QS workers and the separate C9b validation/headroom item.
Both worker paths can prepare independently once their contracts are settled.
Ranking weighs reachable latency/completion benefit, workload coverage and
engineering/measurement cost. These are judgments, **not measured gains**.

| Rank | Existing owner / next deliverable | Perceived benefit and why this position | Confidence / cost and next decision |
| --- | --- | --- | --- |
| 1 | **C3a + C11: current portfolio/ECM opportunity screen, then one missing policy delta** | Broad leverage across ordinary, uneven and recursive inputs; excessive ECM also delays hard-tail fallback. Current 60–100-digit routing still uses fixed32, while source-derived deeper presets already exist. Check joint B1/B2/count, completed-work credit and observable-cost handoff before buying new arithmetic. C11 also attributes preprocessing, classification, rho/PM1, polling and setup so the broader queue follows actual complete costs. | High confidence in a coverage gap; medium confidence in a winning policy. Small initial screen, medium implementation. Sol / xhigh; stop or narrow after the cost/yield decision. |
| 2 | **B12 → F4: opt-in independent ECM curve workers** | Main broad-workload concurrency lever: independent curves can reduce the long ECM portion of random-integer calls. Rank for wall latency while reporting total CPU, startup and cancelled work. | High confidence in parallelizable work, medium in complete-factor benefit; medium ownership/resume cost. Sol / xhigh. Measure 1/2/4 workers, then wider M4 occupancy if funded and useful. |
| 3 | **A7a: current-collector parallel QS relation collection** | Main concurrency lever for the protected balanced tail. Adapt A7 processes to the integrated collector; native shared-memory threads are a conditional alternative with B5a's coarse GIL-free packets. Existing workers reject growing DLP and constrain windows, so current-policy compatibility comes first. | High confidence in parallelizable collection, unmeasured current complete-factor gain; medium/high integration cost. Sol / xhigh. Compare current serial, adapted one-worker and 2/4-worker runs at 60/70; price parent admission, startup, transport and aggregate memory. The process path does not wait for native kernels. |
| 4 | **C9b: mainline timing validation and remaining QS headroom** | C9a is integrated and closed. Establish an unchanged-code A/A timing baseline, then price remaining queue construction, cofactor-recovery allocation and family-wide reuse on the final mainline collector. Coordinate existing A1d/C8/B5a owners rather than duplicating their implementations. | Medium confidence in remaining opportunity; its size is unmeasured after the final hit-word change. Sol / xhigh for bounded attribution/protocol work; separately scoped arithmetic changes retain their owners. Stop after a small informative screen unless a concrete full-cost hypothesis advances. Fresh inputs and uncertainty govern performance confirmation. |
| 5 | **B4a: fused native ECM chunks through CFFI, optional backend** | Move one bounded existing ladder/verified-program chunk into native code so arithmetic stays behind one call. Current GMP representation still runs operator loops in Python; this is a distinct opportunity from a backend toggle or a per-multiply wrapper. It may also enable native worker threads. | High confidence in the boundary gap, unmeasured local gain; medium engineering/ownership cost. Sol / xhigh with unchanged formulas and existing oracles. First compare scalar native execution including conversion and recovery; SIMD/new reducers require their own evidence. |
| 6 | **A1d: prime-major multi-polynomial bucket generation** | Concrete new scheduling gap: current queues are polynomial-local; YAFU builds large-prime buckets across a bounded Gray-code polynomial batch. Amortize prime/root traversal while retaining the existing queue and C10a-01 tables. This differs from C8's CRT half-sum search. | High confidence in the code gap, medium in payoff; medium integration/storage cost. Sol / xhigh. Count actual root/prime visits and charge preparation plus consumed polynomials before complete-factor comparison. |
| 7 | **B3a: selected larger-bound stage-one chain coverage** | A concrete uncovered delta: accepted plans and automatic eligibility stop at B1 = 2,000. Deeper 10,000/11,000/50,000 presets use the ladder. Reusing checked chains may make those bounds economical before a new curve family is needed. | High confidence in the gap, medium in benefit; medium cost if existing verification suffices. Sol / xhigh. Advance only where C3a shows worthwhile larger-bound stage-one cost. |
| 8 | **F3: exact polynomial ECM stage-two feasibility** | Highest-upside ECM algorithmic research for larger continuation bounds and roughly 20–35-digit factors inside 60–99-digit composites. Product/remainder trees and packed multiplication may change the cost curve; they have no local crossover evidence yet. | Medium opportunity confidence; high implementation/proof cost. Astra / xhigh. Compare complete continuation against the actual classical control; research/prototype need not wait for E1 or a C2 win. |
| 9 | **A1c, with B5/B5a as separately labeled alternatives: bulk marking/scanning** | Attacks QS hot-loop cost. Test existing PyPy kernels first; a coarse NumPy or C/CFFI scan, root-update or queue/mark packet can quantify the acceleration tradeoff. Preserve accepted queue eligibility and price buffer conversion. Native speedups remain a separate implementation track. | Medium confidence; low initial scan screen, medium integration. Sol / xhigh. Rank the candidate by measured hot-loop share and complete boundary cost, not a speculative 5–20× kernel gain. |
| 10 | **C2: common-Z/coordinate and selected classical continuation changes** | Smaller experiment than F3; can reduce repeated cross-product work. Current stage two already uses baby/giant steps and inversion-free terms, and earlier pairing losses remain relevant. | Medium confidence in opportunity, lower in a win; medium cost. Astra / xhigh. Keep geometry and coordinate changes separately attributable. |
| 11 | **F2: complete torsion/curve/formula packages** | Potentially improves success per CPU second through both curve distribution and execution. Suyama already has favorable divisibility; merely selecting a curve with 12-torsion is not a new advantage. | Medium opportunity confidence; high proof and statistical cost. Astra / xhigh. Compare setup, stage one, stage-two maps and success, not multiplication counts alone. |
| 12 | **F7: SIQS direct/tree smooth-part screening** | Algorithmic leverage if candidate recovery/division materially costs the current larger-base collector. Exact product/remainder primitives may be shared with F3 after domain-specific proofs. | Medium/conditional benefit; high proof/integration cost. Astra / xhigh. The previous narrower-domain loss remains a control. |
| 13 | **A1a remaining transitions/setup; then A8 → C4/B7/B6/D1/D2** | Reusing trusted metadata can be inexpensive. Matrix/filtering work has strong capacity value but a smaller immediate time ceiling while collection dominates; a demonstrated memory/solve wall moves the matrix branch up. | Setup medium confidence; solver payoff conditional on real fixtures. Sol / xhigh for checked metadata; Astra / xhigh for lifting/filtering/solver proofs. |
| 14 | **A5 → B8: Williams p+1 and selected starts/chains** | Adds complementary coverage, particularly smooth p+1 factors. Expected broad benefit is less established than improving existing ECM/collection paths. | Medium coverage confidence, low quantified portfolio confidence; medium/high cost. Sol / xhigh binary control, Astra / xhigh new starts/chains. |
| 15 | **F8: specialized fast P−1/P+1 continuation** | Important distinct polynomial algorithm at large bounds, but improves a narrower smoothness class. Reuse verified kernels where useful; ECM continuation is not its correctness template. | Conditional on measured PM1/PP1 continuation cost; high cost. Astra / xhigh; A5 is required only for PP1. |
| 16 | **Other triggered deltas: E3/E4/E5/E6, A10a/A2a/C5, C10a-02/-03, C1a/C8, F5/F6/F9** | Select only an uncovered, observed cost/capacity/certainty gap. Rho/PM1 and cheap-path work can matter greatly in their own strata; their current complete-cost share or useful success must justify expansion. Existing losses and C10a-01 overlap argue against blanket sweeps. | Profile/capacity-dependent; task rows retain assignments. C11 reconciles source provenance without turning every reference mechanism into an implementation project. |


Required correctness, accounting and requested-proof work takes precedence over
this optimization queue: B14 remains High for certified results; C10a-06 and
any discovered invalid result/charge are repairs. A11a is optional product work.
E1 stays Deferred; D7/E3/E4/G1/H1/I1 retain their integrated acceptance
dependencies, and GNFS stays Deferred. Component C3a/A7a/F2/F3/F4 work does
not silently resume E1.

All M4 cores are available, with 1/2/4-worker screens preceding a justified
wider configuration. Optimize latency within declared aggregate grants and
report CPU/memory cost; neither ten workers nor a tenfold gain is assumed.
Select one justified Python, NumPy or native kernel candidate per task while
retaining the Python fallback. Refresh this order from current complete-cost
evidence and the balanced completion guardrail, not hypothetical workload changes.

## Broader portfolio opportunities and relative value

The [ranked queue](#ranked-next-work--expected-benefit-and-cost) and execution
rows retain the selected opportunities: C3a/C11 portfolio policy, B4a native
ECM, B5a QS packets, E3 preprocessing, E4 rho, A10a classification, A2a/C5
arithmetic helpers, E5/E6 polling/allocation/reuse and B12/F4/A7a workers.
Their cost/domain triggers remain binding. The earlier comparative rationale
and environment observations are in the [4 November source](https://github.com/nishanth17/factor/blob/7a3ca71c435bd5c5967da5737af728aa5b9bfb01/v2/ROADMAP.md#broader-portfolio-opportunities-and-relative-value).
Detailed inventories stay in local audit storage.

## Coverage decisions from the expanded source review

The [4 November coverage record](https://github.com/nishanth17/factor/blob/7a3ca71c435bd5c5967da5737af728aa5b9bfb01/v2/ROADMAP.md#coverage-decisions-from-the-expanded-source-review)
preserves implemented controls, scoped losses, source-access limits and
mechanism dispositions. A1d retains multi-polynomial buckets; F9 retains
conditional genus-two/Kummer research. Remaining mechanisms keep their existing
execution owners. C9a and C10a-01 are delivered; broader changes require a
distinct current-source hypothesis and their stated correctness/resource gates.
Primary-source pins and detailed comparisons remain local.

## Open and conditional execution items

The rows below state the current priorities and predecessor contracts.
A predecessor settles through an accepted result or recorded retain/defer/reject
decision, subject to the [policy](POLICY.md#dependencies-and-priority).
Each live row links its complete original gates, quantitative context and
research. The [history index](https://github.com/nishanth17/factor/blob/23458a52e409add326c9fd63cf1f322c6e605569/v2/history/README.md) supplies local research paths.

### Execution phase A — independent foundations that can proceed now

| ID / scope | Deliverable and boundary | Predecessors | Status, trigger and existing priority | Model / effort | Evidence and research |
| --- | --- | --- | --- | --- | --- |
| A5 — P5.1 Williams p+1 binary baseline | Implement exact binary Lucas stages 1 and 2, bounded parameter trials, discriminant checks, saturation recovery and checkpoints. | Existing P2 bounded recovery | **Open.** Priority / trigger: Medium, complementary coverage. This is independent of SIQS and PRAC. Binary Lucas supplies the correctness control required before optimized Lucas chains or special starts are ranked. | Sol / **xhigh** | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-a--independent-foundations-that-can-proceed-now) |
| A8 — P3.8 matrix/operator/lifting control | Freeze genuine original-relation and post-filter fixtures, orientation, lifting maps and packed forward/transpose products, including dense-prefix/sparse-tail and selected locality layouts. Independently check dense-oracle products, repeated-coordinate cancellation, short tails, bilinear identities and original-relation extraction; include preparation ownership and possible dense fill. | Accepted R3 identities/lifting and P3.3 control; fixtures refresh after collector changes | **Open; High preparation.** Shared foundation for C4/B6/D1/D2, with current solve/filter/extraction costs and capacity observations. Synthetic fixtures establish correctness/capacity only. A capacity case can justify sparse work before solve CPU dominates; unrelated comparison success is unnecessary. | Astra / **xhigh** design, implementation and correctness; optional separate Sol / **high** frozen fixture measurements | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-a--independent-foundations-that-can-proceed-now) |

### Execution phase B — exploit the settled controls

| ID / scope | Deliverable and boundary | Predecessors | Status, trigger and existing priority | Model / effort | Evidence and research |
| --- | --- | --- | --- | --- | --- |
| B3a — P4.1 larger-bound reuse of checked ECM chains | Extend existing verified plan/executor machinery only to selected useful bounds beyond B1 = 2,000. Compare streamed/ladder, fresh and reused plans with setup amortization, prime-power/frontier validation, strict saturation recovery, finite cache ownership and resume identities. This is a follow-up; completed B3/C6 decisions stay closed. | C3a/C11 establish material stage-one cost at a worthwhile larger bound; existing B3/C6 oracles | **Conditional; high-value opportunity check.** Larger bounds currently fall back to the ladder. Do not run a universal optimal-chain search or automatically enlarge catalogues; unchanged arithmetic and checked generation are the first candidate. | Sol / **xhigh** existing verified machinery; Astra / **xhigh** only for separately selected new chain mathematics; optional separate Sol / **high** frozen execution | [current ECM gates](#ecm-opportunity-and-continuation-gates) |
| B4a — P4.2 fused native ECM chunk boundary | Select one unchanged bounded ECM ladder/verified-program native executor under the [coarse-boundary contract](#b4a--coarse-native-arithmetic-and-whole-portfolio-costs). Keep canonical state, private scratch and the Python fallback. | Existing B3/B4 exact oracles and current committed serial identity; a measured Python arithmetic-loop share; funded native memory/work and cancellation/restore contract | **Open opportunity design; conditional prototype/adoption.** A narrow new native track, not resumption of the stopped broad backend sweep. CFFI is an interface, not a speedup guarantee. Charge setup, conversions, strict GCD/recovery and unused work. Coordinate B12/F4 native-thread feasibility without requiring their prior win. | Sol / **xhigh** unchanged arithmetic, boundary design and implementation; Astra / **xhigh** only for a separately selected new encoded-domain/formula proof; optional separate Sol / **high** frozen execution | [current native-boundary gates](#b4a--coarse-native-arithmetic-and-whole-portfolio-costs); [CFFI API mode](https://cffi.readthedocs.io/en/latest/overview.html#purely-for-performance-api-level-out-of-line) |
| B5 — P3.7 optional NumPy bulk-kernel spike | Compare one coarse threshold/offset scan first, then justified bulk mark/root-hit updates against A1c. Verify an installable PyPy 3.11/arm64 build, actual SIMD dispatch, exact dtype/overflow/duplicate semantics and zero-copy/copy costs. Keep reusable arrays and avoid one call per sparse prime. | Current collector cost share and exact A1c oracle; select B5 or B5a from boundary/toolchain evidence, not both by default | **Conditional; cheap screen eligible.** Native NumPy kernels can use ARM SIMD although PyPy calls cross cpyext. No assumed gain and no dependency-free default change. Full conversion, hit metadata and whole-factor cost govern advancement. | Sol / **xhigh** design/implementation; optional separate Sol / **high** frozen execution | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-b--exploit-the-settled-controls) |
| B5a — optional C/CFFI bulk kernels and ARM SIMD | Select one queue-compatible native QS scan/root/mark packet under the [bulk-kernel contract](#bulk-kernels-pypy-numpy-and-optional-native-acceleration), including conditional reentrant GIL-free packets for A7a. | Current measured kernel share and exact A1a/A1c/A1d oracle; declared buffer/ownership ABI and queue-compatible representation; actual PyPy 3.11/arm64 toolchain | **Open investigation; conditional prototype/adoption.** Include list-to-buffer copies or keep a funded persistent representation. The current sparse queue requires list scores; a new buffer must not silently disable it. No native speedup, SIMD multiplier or global default migration is assumed. | Sol / **xhigh** boundary/design/implementation; separately selected new arithmetic requires Astra / **xhigh**; optional separate Sol / **high** frozen execution | [current acceleration gates](#bulk-kernels-pypy-numpy-and-optional-native-acceleration); [ARM intrinsics](https://arm-software.github.io/acle/neon_intrinsics/advsimd.html) |
| B6 — P3.8 Python-int Four Russians/PLE | Implement one bounded panel-table/PLE challenger against Python-int bitset elimination. Check rank-deficient panels, free-variable recovery and original kernels; charge finite 2**k tables, conversion, pivot search, full provenance, lifting and all simultaneous scratch. | A8 operator/lifting contract and representative fixtures; C4 filtering can proceed independently | **Open; High challenger preparation.** First dense-kernel comparison by evaluation preference. Select panel sizes from local evidence; native widths and digit cutoffs are references. Compare complete preparation-to-extraction cost and factoring impact. No R4/R5, NumPy or prior filtering win is required. | Astra / **xhigh** design, implementation and correctness; optional separate Sol / **high** frozen execution | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-b--exploit-the-settled-controls) |
| B7 — P3.8-R3 compact records, provenance and checkpoint bounds | Implement selected compact atom/exponent/partial/provenance records, bounded replay and merge histories; preflight or stream checkpoint encoding before oversized intermediate allocation. Verify original identities, exact lifting, corruption handling, pinned aliases, eviction and interrupted restore with simultaneous scratch. | Existing exact stores/checkpoints; A8 for matrix representations; C9a's contract for a selected partial-store integration | **Open; High design.** A separately selected partial-store bundle may settle its scoped representation under this contract without completing all B7 or a separate speedup study. A global format change needs its own identity. Lower reservation coefficients or encoded-size quotas alone are insufficient; retain nonunit recovery. | Astra / **xhigh** research, design, implementation and correctness; optional separate Sol / **high** frozen execution | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-b--exploit-the-settled-controls) |
| B8 — Williams p+1 starts (P5.1) / Lucas optimization (P5.3) | Compare bounded rational/seeded p+1 starts and validated Lucas-chain execution against the binary control; validate p+1 bound extension separately from ordinary powering. | A5, A4; A6's schedule-ratio contract | **Conditional.** Priority / trigger: Medium/conditional. Binary correctness and verified chain machinery must exist first. This targets marginal p+1 coverage; the ECM recurrence is not a drop-in Lucas implementation. | Astra / **xhigh** | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-b--exploit-the-settled-controls) |
| B12 — P6.3 ECM/QS worker contracts | Specify stable method/assignment IDs, independent RNG streams, immutable schedules, parent-owned leases, aggregate CPU/RSS, bounded packets, cancellation waste and restart. Reuse A7 infrastructure for F4 ECM and A7a current-collector QS; settle only the contract needed by each workstream. | A2/A3/A7 retained contracts; current committed serial identity for the selected engine; no cross-engine worker dependency | **Open; ongoing in `codex/b12-f4-ecm-workers`. High for wall latency.** Isolated ECM and QS preparation can proceed independently without E1. Aggregate resources, parent bottlenecks and first/full-factor cancellation determine value; measure 1/2/4 workers before a funded wider-core configuration. | Sol / **xhigh** design, implementation and contracts; optional separate Sol / **high** frozen execution | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-b--exploit-the-settled-controls) |
| B13 — P3.8-R2 wide-block resieve capacity | Reuse the accepted sparse resieve support bound. Only if a new calibrated workload still refuses setup, prove a bound for its actual scratch/storage representation, verify coverage and refusal/resume, and compare complete bucket/resieve runs under matched budgets. | A1; a new calibrated wide-block capacity refusal | **Conditional.** Priority / trigger: Deferred unless newly triggered. The combined 30-digit bridge resolves the isolated refusal with the accepted repair. A larger workload needs new evidence and a proof; unrelated feasible collection work proceeds independently. | Astra / **xhigh** | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-b--exploit-the-settled-controls) |
| B14 — P8.2.2 checked prime certificates | Add bounded n−1/Pocklington certificate generation and a separate exact verifier beyond deterministic ranges, with recursive work/storage caps and requested-proof CLI/resume. Record historical APRT-CLE (AP28) as a later coverage candidate only if this first tranche leaves a measured gap; distinguish its probabilistic fallback. | A10; existing shared-budget/checkpoint controls | **Open.** Priority / trigger: High for certified results. Proof is an explicit capability beyond v1 parity. Start after the range contract settles, without waiting for GNFS; certify supported cases rather than claiming arbitrary-size proof from extra probable-prime rounds. | Astra / **xhigh** | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-b--exploit-the-settled-controls) |

### Execution phase C — measured current-engine challengers

| ID / scope | Deliverable and boundary | Predecessors | Status, trigger and existing priority | Model / effort | Evidence and research |
| --- | --- | --- | --- | --- | --- |
| C2 — P5.2 advanced pairing/common-Z | Compare selected distances/relocation/matching/coverage pruning and independent coordinate choices: current projective cross-products, affine/batch inversion and CADO-style inversion-free common-Z tables (AP29). Cost prefix/suffix normalization against saved term multiplications and bounded block/table reuse. | B2 aligned-cell retain/loss decision; current bound-specific continuation cost; coordinate C3a/F3 | **Conditional; medium ECM challenger.** Current stage two already has baby/giant steps; paired executors are retained opt-ins after scoped losses. Certify coverage, nonunit/saturation recovery, storage and replay. Neither a C2 win nor affine inversion is a prerequisite for F3 feasibility. | Astra / **xhigh** | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-c--measured-current-engine-challengers) |
| C4 — P3.8 structured filtering, hybrid cores and extraction | Implement fill-aware bounded higher-weight elimination, including a selected minimum-weight merge-tree/MST heuristic, with compact reverse substitution; compare components, ordering and sparse-to-dense cores. Preserve rank/nullity and lift to original relations. Joint exponent-product/Horner extraction is an optional exact subtranche, preserving both GCD signs and nonunit handling. | A8 contract/fixtures; settle the selected B7 representation contract within any changed-provenance implementation | **Open; High preparation.** Stronger filtering need not wait for B6 to win. B6 remains the first dense-kernel comparison. Measure dimension and weight together, including preparation, fill, retained substitution history, solver, lifting and extraction; fewer rows alone are insufficient. | Astra / **xhigh** changed filtering/provenance/extraction; Sol / **xhigh** only for a separately specified existing kernel; optional separate Sol / **high** frozen execution | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-c--measured-current-engine-challengers); [merge-tree reference](https://github.com/cado-nfs/cado-nfs/blob/692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b/filter/mst.c) |
| C5 — P4.4 reducers | Revisit persistent Barrett/Montgomery contexts only in actual fused engine loops, with exact encoded identities, width bounds and canonical exits. | A2, B4 plus a remaining reduction bottleneck | **Conditional.** Priority / trigger: Low/conditional. Earlier reducers lost near 166–200 bits. Backend/kernel results must provide a reason to reopen them; native `%` remains the default if the whole-run gate fails. A B4a fused native context is a different cost regime; settle its scalar baseline before selecting new Montgomery/Barrett arithmetic. | Astra / **xhigh** | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-c--measured-current-engine-challengers) |
| C8 — P3.8-R2 family-wide CRT/hypercube hit scheduling | Assess whole-polynomial-interval prime/power eligibility and occupancy before bounded CRT half-sum or hypercube scheduling. Verify every Gray/position hit, exceptional roots, recentering, finite tables/queues and charged replay; retain compact hit-layout alternatives with B7. | A1/B1 and a concrete family/root or eligible-prime scanning hypothesis | **Conditional; Low/conditional.** The old first-polynomial probe found zero eligible base primes and sparse eligible powers. Block width is not the eligibility interval. Existing per-polynomial queues/plans are controls, not this family-wide mechanism. | Astra / **xhigh** arithmetic, implementation and ownership; optional separate Sol / **high** frozen execution | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-c--measured-current-engine-challengers) |
| C10a — collection representation, population and scheduling | Own C10a-01/-02/-03/-06 and explicit larger-size feasibility below. Coordinate useful-row preparation/filter scheduling and score/recovery domains with the integrated C9a collector and C9b follow-up; A1a/A1c retain transition/kernel ownership; A1d owns multi-polynomial bucket scheduling. | Closed C10 tranche; current committed control and explicit shared policy/source identities | **Open.** C10a-04/-05/-07 remain Done for delivered scopes. C10a-01's selected bundle is integrated/default and its four-size exploratory feasibility is terminal; no duplicate A/root/multiplier project. Correct, feasible components may share a bundle without prior isolated wins. | Astra / **xhigh** throughout C10a-01; Sol / **xhigh** separately scoped existing-policy design; optional separate Sol / **high** frozen execution for other scopes | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-c--measured-current-engine-challengers) |
| C11 — workload/cost atlas, parameter provenance and competitor protocol | Maintain pinned parameter/caller provenance and a small current complete-cost/yield atlas stratified by modulus size, smaller-factor size and input distribution. Attribute preprocessing/classification, rho/PM1/ECM/QS, polling, allocation, setup and native-boundary costs separately. Use the equal-weight 60/80/100 random-integer workload and 60/70 balanced guardrail; report rough survivors and recursive multiprimes separately. Reconcile current chains and continuation before selecting C3a/B3a/F2/F3/F4. Prepare Python competitor complete-cost comparisons with P6.4/I1. | Current committed controls; existing C3 findings and source presets; coordinate C3a/B12/B3a/B4a and QS A7a/C9a/C10a | **Open; High targeted selection.** This atlas supports rank 1 and does not mandate a full Cartesian campaign. Native-int, GMP, NumPy/native kernels, wrappers and worker counts remain separate tracks; source tables are hypotheses. A specifically selected F3 backend spike does not restart the stopped broad P4.3 study. Final comparative claims retain I1 gates. | Sol / **xhigh** protocol/selection design; Sol / **high** source reconciliation or separately frozen execution; Astra / **xhigh** separately selected new arithmetic/proofs | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-c--measured-current-engine-challengers) |

#### New source-derived deltas — separate follow-up items

| ID / scope | Deliverable and boundary | Predecessors | Status, trigger and existing priority | Model / effort | Evidence and research |
| --- | --- | --- | --- | --- | --- |
| C3a — current ECM effort, recursive credit and handoff delta | Reconcile completed C3 first; profile the current control and select at most one uncovered joint B1/B2/count, completed-search credit or observable-cost handoff hypothesis. Separate ECM/SIQS/classification RNG streams with B12 so changed curve counts do not silently change fallback seeds. Specify valid descendant credit for completed work, finite fallback floors and exact checkpoint lineage. | C3 final decisions; C11 targeted current feasibility/cost evidence; a named uncovered delta; no E1 prerequisite for component work | **Open; High screen/design; implementation conditional on opportunity.** Existing C3 already covers empirical survival/cost scheduling and several deeper presets. Preserve its negative, unstable and early-adoption outcomes. A failed curve is evidence conditional on a workload prior, not a posterior factor-size certificate. Global budget carry-over already exists; descendant search credit and same-curve extension are different mechanisms. | Sol / **xhigh** existing-arithmetic policy/lineage design and implementation; Astra / **xhigh** only for a separately selected new probability/arithmetic proof; optional separate Sol / **high** frozen execution | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#new-source-derived-deltas--separate-follow-up-items) |
| C9b — mainline timing validation and remaining QS headroom | Diagnose unchanged-code A/A timing stability, establish an uninstrumented current-mainline baseline, then selectively attribute queue construction, recovery allocations and family-wide reuse. Compare any justified change at complete cost, including setup, unused preparation, memory and exact validation. | Integrated C9a; existing A1d/C8/B5a/F7 ownership; a newly declared finite protocol | **Open; ranked 4.** The earlier whole profile predates the final hit-word optimization. Measure current bottlenecks before selecting another implementation. Preserve historical observations and flags; use independent inputs, uncertainty and disjoint confirmation for adoption. No run is authorized by this roadmap item alone. | Sol / **xhigh** for attribution/protocol and existing-code work; arithmetic/proof extensions retain their existing task assignments | [current follow-up scope](#c9b--mainline-validation-and-remaining-headroom) |
| A1a follow-ups — checked transitions, quota-aware setup and lifted roots | Reuse validated switch metadata/internal certificates without removing public root verification; bound setup by consumed quota. Investigate cross-B lifted-power roots and fused/paired-polynomial switching as explicit alternatives. Coordinate existing C10a-01 Gray tables and A1c kernels to avoid duplicate implementation. | Delivered A1a stays closed; current setup/lifting hypothesis and current C10a-01 integration identity | **Open design.** Compare setup plus actual consumed traversal; preserve larger-cache and presieve/unroll losses. Bind certificates/caches to base, A/B, representation and backend; prove characteristic-two/singular cases, invalidation, old/new coexistence and charged refusal/rebuild. Derived trust is rebuilt on resume. | Sol / **xhigh** checked metadata/fused kernel work; Astra / **xhigh** cross-B lifting arithmetic and its implementation; optional separate Sol / **high** frozen execution | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#new-source-derived-deltas--separate-follow-up-items) |
| A1b-F1 — fixed-score sparse-queue integration | A separate mechanism: policy-bound fixed-score queue weights/eligibility, exact threshold/power coverage, scratch refusal and charged rebuild. Do not modify historical A1b controls. | Completed A1b existing-option decision; a new credible cost hypothesis and separately frozen design | **Deferred.** The implemented fixed/list bundle lost all 50/55/60 screens. No queue extension or additional campaign is admitted by that result. | Sol / **xhigh** integration/ownership; **high** frozen evaluation | [A1b decision and boundary](benchmarks/RESULTS.md#a1b-existing-fixed-score-review) |
| A1c — bulk marking/scanning and same-modulus two-root interleaving | Select a bounded PyPy bulk mark/scan or paired-root kernel; consider fused root-switch/marking with A1a. A1d owns prime-major multi-polynomial scheduling. Verify sorted translated offsets, tails, duplicate/singular roots, saturation, hit metadata and reservation/cancellation/interrupted-block replay. Two roots and two polynomials are distinct hypotheses. | Settled A1a/queue controls, exact kernel oracles and a current marking/scanning cost hypothesis | **Open; High design/implementation.** Keep queue eligibility and backend consequences explicit, preserving low-allowance behavior and physical coverage. Native SIMD/byte constants are references. Retain prior unroll/presieve losses and compare complete factoring after a bounded identical-work screen. | Sol / **xhigh** design, implementation and correctness; optional separate Sol / **high** frozen execution; a named unresolved proof requires a separately selected Astra / **xhigh** scope | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#new-source-derived-deltas--separate-follow-up-items) |
| A1d — prime-major multi-polynomial buckets | Transpose prime/polynomial queue construction across a funded Gray batch under the [batching contract](#a1d--prime-major-bucket-generation-across-polynomials). Preserve C8's separate CRT scope. | Current queue/control, current family/Gray identity and selected C10a-01 table interface; source/policy identity fixed before measurements | **Open design; medium/high opportunity.** Verify identical per-polynomial hit coverage, ordering, tails/powers/exceptions, quota/cancellation and charged rebuild. Reserve only funded batches, charge unused precomputation and compare build plus consumed traversal. No C8 win or new polynomial representation is required. | Sol / **xhigh** design/implementation/ownership; Astra / **xhigh** only for a separately selected new arithmetic proof; optional separate Sol / **high** frozen execution | [current batching gates](#a1d--prime-major-bucket-generation-across-polynomials); [pinned executor](https://github.com/bbuhrow/yafu/blob/c7d90283972820f99b3b4920bfc1eff05c1f6342/factor/qs/update_poly_roots_32k.c#L2747) |
| A2a — scoped inverse/GCD/quotient helper follow-up | Select one profiled helper: batched binary/Montgomery inversion or small-quotient shortcuts (AP32/AP33) against current PyPy Euclid/GMP, or a scoped native GCD helper against math.gcd/GMP. Include conversion, allocation, operand-size class and nonunit recovery; packed/tree operands differ from scalar moduli. | Completed A2/B4; C11 profile and a changed workload/cost hypothesis | **Conditional.** Priority / trigger: Low/conditional. Preserve backend/reducer losses. Native fixed-limb/SIMD/FFT/NTT references (AP34/AP35) remain deferred C11/C5 dispositions, not a new Python limb engine. Current GCD already uses math.gcd/GMP; optimized native internals are the control, and large packed/tree operands differ from 60–99-digit scalar moduli. | Sol / **xhigh**; targeted Astra proof review if required | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#new-source-derived-deltas--separate-follow-up-items) |
| A7a — current-collector parallel SIQS collection | Adapt A7 collection to integrated C9a/C10a policies; compare processes and conditional B5a native threads under the [current-worker contract](#a7a--current-collector-parallel-relation-collection). | Completed A7; selected B12 lease/RNG/resource contract; current committed serial feasibility; B5a ABI only for the native-thread alternative; no F4 win or E1 prerequisite | **Open; High, rank 3.** Reuse checked merge/resume and one common partial/cycle store. Growing DLP is currently rejected and worker windows are smaller than current serial defaults. Price parent admission, graph/filter/solve, transport, lost work, aggregate RSS and startup. Existing Python thread mode is not multicore throughput evidence. | Sol / **xhigh** implementation, ownership and experiment design; separately selected new graph/provenance mathematics uses Astra / **xhigh**; optional separate Sol / **high** frozen execution | [current parallel QS gates](#a7a--current-collector-parallel-relation-collection); [A7 framework](qs/parallel.py) |
| A10a — classification follow-up to completed primality transfer | Evaluate MR-residue proper-factor recovery, optional BPSW composite rejection and grouped small-prime remainders (AP25–AP27), or a scoped native witness batch after pricing early rejection only where E3/C11 profiles identify an uncovered benefit. | Completed A10; E3 frozen control or a justified isolated oracle spike; classification-cost evidence | **Conditional.** Priority / trigger: Conditional. Preserve deterministic ranges, requested MR rounds, RNG/work identity and certainty. A BPSW pass is not proof; native library repetitions are not automatically the requested MR-round contract; prior primorial losses remain the control. B14 retains certificate ownership. | Sol / **xhigh** | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#new-source-derived-deltas--separate-follow-up-items) |
| A11a — optional input/browser follow-up | Separately select validated factor hints/opt-in lookup, bounded batch/expression input or browser progress/offline/save-continue (AP38–AP40). Give each selected subtranche its own API and acceptance boundary. | Completed A11; explicit product selection, outside the optimization critical path | **Conditional.** Priority / trigger: Optional; no factoring speed claim. Preserve quiet library APIs, exact reconstruction, certainty, finite input/output/network allowances and checked checkpoints. | Sol / **high**; **xhigh** for changed accounting/resume | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#new-source-derived-deltas--separate-follow-up-items) |
| C1a — cheaper SLP anchor reuse follow-up | Compare an anchor-only SLP implementation or changed retention regime (AP01) with ordinary SLP and the existing graph-SLP/DLP controls. Prove residual cancellation and bounded pinned ownership without adopting source dependency-discard policies. | Completed C1; C9 decision; a new cost/retention hypothesis | **Conditional.** Priority / trigger: Low/conditional. C1 graph-SLP training totaled 20.055 s versus 20.160 s for SLP, with no established win; DLP already retains anchors. Do not repeat that arm or promote on raw cancellation counts. | Astra / **xhigh** | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#new-source-derived-deltas--separate-follow-up-items) |

#### A7a — current-collector parallel relation collection

The wall-time objective makes this an immediate QS workstream, ranked third,
above the separate C9b validation/headroom item. Completed P3.6/A7 supplied
bounded workers, stable assignments, checked parent admission, drained
checkpoints and cancellation accounting; their
historical losses and experimental status remain. The missing deliverable is
current-collector compatibility and a fresh complete-factor latency comparison.

Current `ParallelConfig` rejects `GrowingDoubleLargeSieveConfig`, caps
half-width at 8,192 and total owned memory at 1 GiB, and the pool permits at
most four workers. Current calibrated serial collectors use growing DLP and
half-width 65,536. These are integration gaps, not evidence that merging is
intrinsically expensive or that current workers already accelerate that control.
Do not weaken the serial control to fit historical worker limits.

**Execution alternatives:** adapt the pure-PyPy process path first unless an
already available native boundary makes threads cheaper to investigate.
Native shared-memory threads are an explicit conditional A7a alternative,
coordinated with B5a; they do not gate process-worker preparation. The existing
thread executor still runs Python sieve/recovery code under PyPy's GIL.
Thread scaling requires substantial work inside a verified GIL-free native
packet, with immutable shared tables, private mutable buffers and bounded batch
publication to the parent collector. Keep exact Python validation and measure
the remaining serialized recovery/admission cost. A CFFI call per prime or
Python callback per candidate is not an adequate coarse boundary.

Reconcile the selected current base/multiplier, assignment/flyer, root, score,
queue and graph identities before building. Reuse A7's task transport and
parent-owned leases; keep exact atom provenance, cross-worker partial matching,
duplicate/repeated-prime handling, bounded publication and restart. Fund
worker-private state, pending batches, parent graph/filter/solve and checkpoint
coexistence. A selected compact-record contract can be settled with B7 within
this workstream; a global storage rewrite is not a prerequisite.

Select task granularity from current setup, coordination and tail-waste costs:
different A families, private B-polynomial snapshots within a shared A family,
or disjoint interval tiles. Prefer bounded dynamic leases on the heterogeneous
M4; do not require every decomposition. Keep score/root/queue scratch private
rather than introducing shared per-hit writes by default. Price pending-batch
ownership and ordered-merge stalls before changing deterministic admission or
checkpoint order. All workers contribute partials to the common graph; useful
rows and filtered progress, rather than raw atom counts alone, guide stopping.

Start from a small uninstrumented serial feasibility baseline at 60 and, when
informative under a declared grant, 70 digits. Consume compatible C10a-01
and C9a evidence with their source identities and limitations; preserve their
historical freezes rather than duplicating the size studies.
Use current serial and an adapted one-worker arm to expose integration cost,
then 2/4 workers to measure scaling. Preserve a deterministic assignment stream
for fixed-work checks and declare complete-factor seeds/order separately.
If collection policy must change, label the bundle; its gain does not isolate
the effect of concurrency. A historical exporter cannot stand in for the
strongest current serial comparison.

For a justified native-thread screen, compare one native worker with 2/4 native
workers to separate concurrency from the native kernel's scalar gain; also
compare against the strongest valid process path and current serial control.
Declare a serial engagement rule and separate cold CLI/pool startup from warm
library calls. Reuse A7/B12 ownership and cancellation contracts, keep native
source self-contained and retain the Python fallback. SIMD is independent of
threading; a scalar native packet can be sufficient. Nested backend pools need
an explicit aggregate grant rather than multiplying worker counts silently.

Freeze a finite staged campaign after feasibility: each stage names its
question, admission/stop decision, runtime estimate and cumulative wall/CPU,
memory/work/storage allowances. Complete-factor latency includes parent work,
merge/filter/solve, process/JIT startup, copies, cancelled work and draining.
Check parent throughput with short diagnostics before buying long scaling
runs. Fixed-work throughput is exploratory, not a complete-factor speedup.
Above four workers requires explicit pool/resource-limit changes and new
ownership checks; all ten cores are available, with useful occupancy measured.

Promotion uses fresh independent inputs and the 60/70 balanced completion
guardrail. Historical family-sharing or alternate chunk sizes enter only if
current coordination/setup evidence motivates them. Literature and current
implementations establish feasible decomposition, not a measured M4 speedup.
The local source/decision ledger is
`v2/audit/research/qs_parallel_review_20261104.md`.
F4 later composes selected worker results; it does not gate isolated QS
feasibility or require an ECM win.
Keep current C9a/C10a source identities explicit and integrate shared edits
sequentially. No campaign is launched by this roadmap change.

#### C9a — combined useful-yield and storage workstream

Delivered scope: [C9a completion](#c9a-completion). The original combined
contract and chronological decisions remain in the [preserved source](https://github.com/nishanth17/factor/blob/7a3ca71c435bd5c5967da5737af728aa5b9bfb01/v2/ROADMAP.md#c9a--combined-useful-yield-and-storage-workstream).
C9b owns current validation/headroom; B7, C9b, F6 and F7 retain their distinct
representation, attribution, higher-arity and smoothness/batch boundaries.

#### C9b — mainline validation and remaining headroom

**Open.** Start from the final integrated C9a collector, retaining its exact
mathematical work and validation. The existing whole-factor profile predates
the last hit-metadata optimization, so its percentages cannot establish the
remaining performance ceiling.

First investigate the changing timing baseline with unchanged-code A/A checks
and representative uninstrumented work. Distinguish engine/JIT warmup from
host variability. The historical 15% engineering heuristic is not a new
acceptance gate: retain its observations/flags, report warnings, and derive
any future rejection rule from measured variability and uncertainty. Declare
the question, stopping decision, expected runtime and finite
cumulative CPU/wall/storage limits before execution; no campaign starts here.

Then screen queue construction, cofactor-recovery bookkeeping and allocation,
buffer reuse, and polynomial-family root/bucket preparation against current
complete costs. Price setup, unused precomputation, memory, cancellation and
checked resume. Advance only a measured concrete hypothesis; retain a finite
fallback when a cache or batch is unfunded. A1d owns multi-polynomial buckets,
C8 family-wide CRT/hypercube scheduling, B5a native kernels and F7 smoothness
screening; C9b coordinates attribution and selection rather than duplicating
those workstreams.

Further adoption needs matched complete-cost/completion evidence on independent
inputs, uncertainty and disjoint fresh confirmation, including the 60/70-digit
balanced guardrail. Keep component, diagnostic and whole-factor claims
separate; the 80-digit band first needs a feasible uninstrumented baseline.
An unrelated control failure must not block measuring a valid candidate.
Stop when a screen is uninformative rather than enlarging it automatically.

#### A1d — prime-major bucket generation across polynomials

The current sparse schedule key contains one polynomial identity and the
collector releases it on a polynomial switch. YAFU's explicit polybatch root
executor nests a polynomial batch inside a prime block; this is a concrete
alternative to repeatedly loading the same prime/root metadata. It differs
from a larger single-polynomial queue, two-root interleaving, C10a-01's modular
Gray tables and C8's sorted half-sum search for family-wide eligible primes.

Select a small funded batch and exact enumeration order. Compare prime/root
visits, root-update/setup work, bucket generation/application and consumed
positions against the current queue; include unused precomputation when
factoring stops early. Prototype native Python storage first or explicitly
label a B5a accelerated boundary. Prove all per-polynomial hits, signed windows,
powers, singular roots, A/multiplier exceptions and family boundaries. Pending
atoms must identify the right polynomial; quota refusal/cancellation/resume
must neither publish future rows early nor repeat completed work. Derived
buckets may rebuild under charged bounds. Do not inherit probabilistic bucket
capacity or dropped-overflow policies from native references.

#### A1a follow-ups — checked transitions and lifted roots

Internal certificates may reuse already validated metadata while public root
verification remains intact. Bind them to the checked base, A/B transition,
representation, backend and consumed Gray quota. Corrupt or mismatched
certificates must reject; resume rebuilds derived trust under charged bounds.
Keep cross-B prime-power reuse distinct from C10a-01's base-root Gray tables
and existing per-polynomial plans. Test characteristic two, A/multiplier
divisors, singular roots, recentering, invalidation, interrupted construction
and simultaneous old/new ownership. Fused switching/marking and paired
polynomials need explicit physical-trial ordering; pairing two roots is
different from pairing two polynomials.

The delivered A1a integration stays Done. Preserve larger-cache, local-bound
and combined presieve/unroll losses. Old sparse-queue/face preparation screens
were superseded by the completed C10a-07 panel; they are not new A1a work.
See the [retained decisions](benchmarks/RESULTS.md#c10a-07-face-and-queue).

### Execution phase D — sparse-solver alternatives

| ID / scope | Deliverable and boundary | Predecessors | Status, trigger and existing priority | Model / effort | Evidence and research |
| --- | --- | --- | --- | --- | --- |
| D1 — P3.8 packed block Lanczos | Implement one genuine selected-block recurrence with packed forward/transpose products, singular-block handling, original-operator kernel correction, lifting, bounded seeded restarts and checked resume. Consider Thome's reduced-vector variant at recurrence selection; independently verify every original M d = 0. | A8 operator/lifting contract; use B6/C4 as comparison controls, not required successes | **Open; High preparation.** First sparse candidate by scheduling preference. Compare all recurrence/recovery/restart and retained-vector costs with feasible dense/hybrid controls. Capacity may justify work before solve CPU dominates; native solver names or SpMV throughput are insufficient. | Astra / **xhigh** design, implementation and proof; **max** only for a named unresolved invariant; optional separate Sol / **high** frozen execution | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-d--sparse-solver-alternatives) |
| D2 — P3.8 block Wiedemann | Implement a genuine projected Krylov, block-polynomial-generator/lingen and reconstruction path. Verify original kernels, independent sequences and bounded I/O/checkpoints; cost retained state, sequence, generator, reconstruction and failures together. Select a subquadratic approximant/generator only when the simpler exact generator is material; mod-2 polynomial arithmetic has a separate domain from F3 mod-N kernels. | A8 operator/lifting contract; B6/C4 are comparison controls, not required successes | **Conditional alternative to D1; no D1 dependency.** Prefer Lanczos for the current single-machine tranche unless capacity, sequence parallelism or restart evidence favors Wiedemann. Begin with an exact base-case generator; the calculator reference does not supply this solver. | Astra / **xhigh** design, implementation and proof; **max** only for a named unresolved invariant; optional separate Sol / **high** frozen execution | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-d--sparse-solver-alternatives); [generator algorithm](https://members.loria.fr/EThome/files/jsc.pdf) |

### Execution phase E — current-engine acceptance, control freeze and P8 experiments

| ID / scope | Deliverable and boundary | Predecessors | Status, trigger and existing priority | Model / effort | Evidence and research |
| --- | --- | --- | --- | --- | --- |
| E1 — final R1/R5 and P4/P5 acceptance | Integrate selected changes, including accepted A10/B14 capabilities and A11 CLI integration; reconcile v1 transfer decisions and recalibrate affected SIQS/ECM/p±1 parameters; freeze selections; run fresh complete-factor comparisons, R5 SSS/worker challengers, resume checks and clean-checkout validation; publish adopt/defer/reject decisions. | A1, A7, A10, A11 and chosen P3.7/P3.8/P4/P5 predecessor decisions; B14 if proof support is selected; B13/C8 and other unchosen branches have recorded deferrals | Deferred; resumption requires explicit authorization. Original combined confirmation and all unrun gates remain required for experimental acceptance. | Sol / **xhigh**; targeted Astra / **xhigh** review for new arithmetic/provenance | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-e--current-engine-acceptance-control-freeze-and-p8-experiments) |
| D7 — P8.1 fresh portfolio control | Freeze the accepted current SIQS/MPQS, rho, ECM and p±1 portfolio with certified training/untouched confirmation inputs; profile separately and retain source/configuration identity for every later challenger. | E1; existing SIQS/P2 controls | **Open.** Priority / trigger: High enabling value, bounded scope. P8 tuning uses the accepted current engines before GNFS. Keep this control immutable; future GNFS requires a separately versioned extension and crossover evaluation. | Sol / **high**; **xhigh** for protocol decisions | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-e--current-engine-acceptance-control-freeze-and-p8-experiments) |
| E3 — P8.2 preprocessing | Compare inherited negative trial coverage, its exact bound on possible power exponents, high-valuation removal and same-integer valuation-gcd pruning before costlier Fermat/Lehman/Hart/SQUFOF or structured a^k±1, power-sum/difference and historical Fibonacci/Lucas candidates (AP18–AP24). A10a owns classification deltas. Parent non-power/compositeness results do not automatically hold for children; C3a coordinates lineage. | Current committed preprocessing/classification oracle and observed cheap-path opportunity for isolated design/feasibility; D7 for integrated portfolio tuning | **Conditional; medium/profile-gated.** Keep prior primorial/root/pretest losses. Use ordinary integers, rough survivors, powers, small residuals and close/special-form strata; prevent cheaper screening from delaying expensive-path completion. No restart of E1/D7 is implied. | Sol / **xhigh** | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-e--current-engine-acceptance-control-freeze-and-p8-experiments) |
| E4 — P8.3 rho calibration | Tune existing Brent batch/walk/restart and finite effort under identical total budgets/seeds; compare a selected native whole-batch kernel only if rho cost matters. Include local saturation replay, smaller-factor strata, failure expenditure and the next engine. Existing Brent/local loops are controls. | Current committed rho/recovery control for isolated screening; D7 for integrated portfolio tuning | **Open preparation; medium/conditional experiments.** Rho has useful small-factor and residual coverage, not a practical balanced 60–80-digit growth path. Independent-walk parallelism has different scaling from ECM curves. Bigger iteration counts and a new kernel each need marginal complete-cost evidence. | Sol / **high**; **xhigh** for changed recovery semantics | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-e--current-engine-acceptance-control-freeze-and-p8-experiments) |
| E5 — P8.5 recovery/polling/checkpoint costs | Attribute polling, atomic state commits, replay, result validation, serialization and durable I/O separately. Reuse existing C3 reservation/clock repairs and rho chunking; select only remaining overhead. Fund bounded chunks, preserve deadline/cancel precedence, refusal states, strict factor recovery and charged resume; native buffers need deterministic release. | Current committed boundary/accounting contracts and a measured remaining cost for isolated feasibility; D7 for integrated portfolio tuning | **Conditional; medium/profile-gated.** A smaller clock-call count or larger atomic unit is not a free win if it changes responsiveness, work coverage or validation. Price checkpoint policy and restoration honestly; no automatic rewrite of the shared budget API. | Sol / **xhigh** | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-e--current-engine-acceptance-control-freeze-and-p8-experiments) |
| E6 — P8.6 contexts/schedules and bundled prime generation | Select reusable immutable prime/power/product metadata, lazy setup, compact allocation, warm worker contexts or a source-bundled C11/CFFI prime-batch provider. Reuse existing wheel-30/segmented context and sorted-output machinery; adapt it to finite Factor packets. Separate modulus state, cache ownership, native buffers and semantic identity; preserve prior Python wheel/cache/rolling losses. A1/A1d/C10a retain QS kernels. | Current context/cache contracts and C11 prime-generation/consumption cost evidence; D7 for integrated tuning; coordinate B4a/B12/F4/A7a | **Conditional; cheap opportunity screen eligible.** Bundle source/build adapter with a Python fallback; no sibling-runtime or external prime-sieve-library dependency. Price setup, extraction/copies, memory and unused preparation. Native threads start at one inside factoring workers; wider prime generation needs its own complete-cost evidence and aggregate grant. | Sol / **high** cost screen; **xhigh** native ABI, ownership, rejection/resume and implementation; optional separate Sol / **high** frozen execution | [prime-provider gates](#e6--self-contained-native-prime-generation); [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-e--current-engine-acceptance-control-freeze-and-p8-experiments) |

### Execution phase F — independent P6 challengers

Use current committed controls for isolated F2/F3/F4/F6/F7 feasibility;
deferred E1 is not a blanket prerequisite. The rows retain their exactness and
engine-specific contracts; later combined acceptance keeps E1. Choose
challengers with concrete cost/yield hypotheses and finite designs.

| ID / scope | Deliverable and boundary | Predecessors | Status, trigger and existing priority | Model / effort | Evidence and research |
| --- | --- | --- | --- | --- | --- |
| F2 — P6.1 torsion-aware curve/formula packages | Compare complete packages against Suyama: GMP-ECM small-parameter/batch formulas (MG17), CADO Montgomery12/16 and mixed Montgomery/twisted-Edwards bytecode (MG18), Edwards/windowed or sliding-window alternatives, Suyama-11 and the calculator reference family (AP30). Screen JKL twisted-Hessian/small-parameter families with their field-dependent torsion conditions before selecting a package. Cost curve generation, chain setup, nonunit maps and continuation conversion together. | A2/B3/B4 retained arithmetic/oracle contracts; C11 current factor-size/cost hypothesis; E1 only for later combined acceptance | **Open opportunity design; conditional implementation.** Isolated feasibility may precede E1. Baseline Suyama already biases divisibility; distinguish rational torsion, reduction group order and chosen-point order. Require success per complete CPU/wall cost, finite family exhaustion and close-factor saturation recovery, not a promised free B1 increase. | Astra / **xhigh** research, implementation and correctness; optional separate Sol / **high** frozen execution | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-f--independent-p6-challengers); [JKL candidate](https://doi.org/10.1112/S1461157016000231) |
| F3 — P6.2 exact polynomial ECM continuation | Select bounded product/remainder-tree and multipoint continuation from GMP-ECM and Yamaquasi. Compare schoolbook/Karatsuba and exact Kronecker packing on native ints, with optional coarse GMP multiplication as a separate track. Include monic reciprocal reduction, blocked/tree storage, exact evaluation/recovery and a separately flagged Brent–Suyama power/Dickson extension; see the current gates below. | Current B2/classical control and C11/C3a bound-specific opportunity; selected polynomial oracles/contracts; E1 only for combined acceptance; C2 win unnecessary | **Open; High algorithmic feasibility, conditional full integration.** Target useful 20–35-digit smaller factors within 60–99-digit inputs. Current production is classical baby/giant continuation, not a promoted paired executor. Establish a whole-continuation crossover before enlarging B2 or claiming GMP FFT pays off. | Astra / **xhigh** research, implementation, exactness and experiment design; optional separate Sol / **high** frozen execution | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-f--independent-p6-challengers) |
| F4 — P6.3 opt-in ECM workers and integrated portfolio parallelism | Reuse B12 leases and stable disjoint curve IDs for opt-in process-level ECM. Compare 1/2/4 workers on fixed committed curves and first-factor latency, charging startup/JIT, cancellation, in-flight waste, total CPU and aggregate RSS; only then consider additional heterogeneous cores. After B4a, compare native-thread execution only when GIL release and thread-private state are verified. Later budgeted mixed-engine racing retains QS/A7a/G1 and engine-specific contracts. | B12 selected lease/RNG/checkpoint contract and current committed serial control for ECM; E1 for later combined portfolio acceptance; E2 only if GNFS included | **Open; ongoing in the B12/F4 worktree.** Isolated ECM feasibility need not wait for E1 or a new curve family. Fixed-work throughput does not establish first-factor or full-factor gains. Match total CPU for efficiency claims; report wall-latency experiments under their declared parallel CPU grant. | Sol / **xhigh** implementation/ownership and design; optional separate Sol / **high** frozen execution | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-f--independent-p6-challengers) |
| F5 — P6.2 cross-family A=A0*q polynomial reuse | Compare bounded reuse across QS families with calibrated factor-base-smooth A, explicitly carrying any external q exponent, partial-relation role and duplicate policy. | E1; B1 calibrated QS control | **Conditional.** Priority / trigger: Conditional on polynomial/root setup cost. This changes relation semantics and is not a drop-in root cache. Ordinary families and the R1 external-square distinction provide the control, not a correctness shortcut. | Astra / **xhigh** | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-f--independent-p6-challengers) |
| F6 — P6.2 higher-arity large-prime QS relations | Specify a genuine sparse LP-incidence/provenance domain, exact repeated-factor corrections and bounded residual splitting for a selected 3LP/general relation challenger. Compare useful dependencies and complete factors against calibrated DLP. | C1/current DLP control and an explicit higher-arity domain/capacity design; E1 is unnecessary for isolated feasibility | **Low/conditional investigation.** A residual census may nominate work; native 135-digit results do not establish a local crossover. Additional splitting/storage/certification must be charged. A two-endpoint cycle forest does not represent general higher arity. | Astra / **xhigh** research, design, implementation and correctness; optional separate Sol / **high** frozen execution | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-f--independent-p6-challengers) |
| F7 — P6.2 SIQS direct/tree smooth-part screening | Integrate bounded SIQS direct/tree smooth-part screening under its [domain and recovery contract](#f7--siqs-directtree-screening-and-distinct-residual-batches). Residual-prime batches require a separate domain proof. | Current committed bucket/queue collector, a larger-basis cost hypothesis and the selected SIQS-domain ownership contract; E1 is unnecessary | **Open; High design.** Compare identical candidate captures against current bucket recovery, then whole latency/storage/completion. Preserve the old R2 loss in its narrower domain. Residual-prime batch cofactoring has a distinct C9a/F7 domain and may join a justified combined candidate; new solvers and higher arity keep their owners. | Astra / **xhigh** research, design, implementation and correctness; optional separate Sol / **high** frozen execution | [current integration gates](#f7--siqs-directtree-screening-and-distinct-residual-batches); [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-f--independent-p6-challengers) |
| F8 — P6.2 specialized fast P−1/P+1 continuation | Investigate Montgomery–Kruppa reciprocal-Laurent-polynomial construction and geometric-progression evaluation through pinned GMP-ECM pm1fs2.c, separately from F3 ECM/Brent–Suyama. Reuse verified polynomial kernels only where their exact domains agree; do not transfer ECM extension or bound semantics blindly. | A6 control and measured remaining PM1 stage-two cost; A5 only for PP1; selected polynomial contracts; E1 only for combined acceptance | **Conditional; later algorithmic challenger.** Larger-bound smoothness coverage may justify it. Compare complete construction/evaluation, backend and storage cost; C2 geometry is not a dependency. Generic polynomial ECM continuation is not a substitute for this specialized family. | Astra / **xhigh** | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-f--independent-p6-challengers) |
| F9 — P6.2 decomposable genus-two/Kummer factoring | Screen Cosset/GMP-HECM: special decomposable Jacobians represent two elliptic searches through Kummer arithmetic. Require exact parameter generation, exceptional maps, both image curves, bounded recovery/checkpoints and probability per complete cost. | A separately selected larger-factor opportunity after cheaper B3a/F2/F3 controls; independent oracles and a plausible local cost model | **Low/conditional research; implementation not selected.** A historical native HECM win does not establish a PyPy crossover at 20–35-digit factors. Treat this as a distinct algorithm, not ordinary Edwards ECM or a reason to port Python limbs. | Astra / **xhigh** if explicitly selected; optional separate Sol / **high** frozen execution | [current coverage](#coverage-decisions-from-the-expanded-source-review); [Cosset](https://arxiv.org/abs/0905.2325) |

#### F7 — SIQS direct/tree screening and distinct residual batches

`SmoothBatch` already implements Bernstein smooth-part arithmetic. The SIQS
delta is placement, bounded dispatch and domain/ownership support, not a
duplicate generic algorithm. Its current 100,000-prime endpoint and 8,192
total-node cap are independent barriers: a 5,000-leaf base tree needs 10,005
nodes. Separate the base-product builder's capacity from candidate-batch
capacity while preserving legacy refusals and exact bounds.

Compare direct radical reduction and remainder-tree screening on identical
current candidates against bucket recovery, which already avoids full-base
division. Account for scalar exponent recovery on admitted values, signed/zero
norms, duplicates, all valuations, A contributions and normalization. Keep
queues and unrelated scoring choices unchanged where supported; disclose any
fallback as a bundled change. Prove pending-batch commit/refusal/resume and
simultaneous radical/tree/candidate/output ownership before scaling.

The old R2 small-domain batch loss remains evidence for that prototype;
it is not a universal rejection of larger-base dispatch. Short candidate
comparisons explain a mechanism; fresh complete-cost/completion evidence
governs adoption. Residual-prime-product batch cofactoring uses a different
prime domain, latency and certainty contract, owned with C9a. It may enter a
justified combined candidate under that explicit contract; it is not supplied
by the existing base-smoothness routine. Higher arity remains F6.

#### ECM opportunity and continuation gates

C3a/C11 start from the actual committed engine: Suyama Montgomery X:Z,
classical baby/giant continuation with projective cross-products, optional
paired/wheel executors after scoped losses, and reused PRAC only for the
supported B1 = 2,000 schedule. Other-size automatic ECM remains
`(2000, 147396, 32)`; deeper explicit presets already exist. Completed C3
already investigated cost/survival and joint depth/count policies at smaller
sizes. Its unstable results and promotion before fresh confirmation stay
recorded. None establishes a current optimal 60–99-digit policy.

**Opportunity before implementation.** Reuse valid retained costs first;
otherwise design a small uninstrumented current feasibility screen, then
targeted identical-work attribution with measured observer overhead. Separate
planning/setup, stage one, continuation, recovery, repeated child work and
fallback. A fraction f of serial cost bounds the saving from eliminating that
component to f; a hypothetical kernel speedup r gives total speedup
`1 / ((1 - f) + f / r)` only when coverage and the remaining cost are unchanged.
Policy, success-probability and worker changes require their own model. State
which gap can plausibly change completion before scheduling a large campaign.

**Workloads and models.** The primary population remains the policy-defined
60/80/100-digit random-integer bands. Diagnostic strata may include
60/80/99-digit ordinary integers, declared
trial-sifted survivors, recursive composites with several medium factors,
and planted smaller factors around 15/20/25/30/35 digits at compatible total
sizes. Include 60–80-digit balanced controls, prime/power/close-factor cases
and roughly 35:65 splits; these are strata, not a mandatory Cartesian grid.
Report preprocessing exits and conditional ECM entrants separately. Hidden
factor sizes label offline evidence, never runtime dispatch. Give median,
tail/mean cost, complete reconstruction and censored completion by class;
do not quote a random-integer asymptotic as the measured service distribution.
Author ECM parameter tables are version/backend/curve-family dependent and
guide a local grid; their expected curves are not native-PyPy promises.

**Policy and recursion.** Distinguish completed-search credit from cumulative
resource charges, which already survive a split. Record completed bounds,
curves/families and child lineage before deciding what credit is justified;
an interrupted or factor-producing curve need not have completed the same
search on both children. Do not inherit parent compositeness, reset budgets,
or silently skip necessary child classification. A prior plus an appropriate
failure likelihood is needed for a posterior; absent a defensible prior,
use a declared empirical observable-state policy and sensitivity analysis.
Current shared random draws couple curve counts to fallback seeds: a new
method-stream identity/checkpoint version must be explicit, with old restore
semantics preserved. Compare completed policies under fresh independent
inputs, not isolated curve timings or revealed training winners.

**Same-curve extension is a distinct conditional C3a/P5.2 delta.** Reusing a
stage-one point at a larger B1 needs the exact lcm-exponent ratio, checked
prime powers, nonunit/saturation recovery and funded migration. Extending B2
needs exact already-covered intervals and resumable continuation state.
Neither descendant effort credit nor existing p−1 extension proves this.
Implement only when saved preparation/point work justifies the state cost.

**F3 exactness and cost.** Begin with degree-one continuation, matched against
the current classical coverage oracle. Separately flag power/Dickson
Brent–Suyama extensions and label additional success outside literal prime
coverage. Compare blocked product accumulation modulo a baby polynomial,
monic reciprocal/remainder trees and direct or transposed multipoint
evaluation; select from actual degree/bound costs, not a full implementation
of every variant. For nonnegative canonical coefficients and lengths a,b,
each convolution coefficient is at most
`min(a, b) * (N - 1)**2`; choose a packing base strictly above that bound.
Prove any signed/truncated/reciprocal variant independently. Monic division
works over the composite-modulus ring; arbitrary coefficients need not be
units. Preserve denominator/GCD recovery, saturated products and exact replay.

Native Python-int packing is a valid control; optional coarse gmpy2 packing
may benefit from GMP's Toom/FFT dispatch only at suitable measured sizes.
Count packing/unpacking, conversion, multiplication, reduction, tree build,
evaluation, validation and recovery. A fast multiplication alone is no F3
acceptance. Compare retained/rebuilt/blocked trees with simultaneous scratch,
owned-memory and RSS bounds; keep this analysis proportionate to the selected
sizes. Disk-backed trees are optional only when a real memory wall warrants
their I/O and checkpoint contract. Direct NTT/CRT or floating convolution is
a later selected backend with an exact reconstruction proof, not a necessary
first implementation. F8's reciprocal-Laurent/geometric algorithm is distinct.

For every selected stage, follow POLICY.md for validity and write the
question, admission/advance/stop rule, expected runtime and finite cumulative
CPU/wall/work/storage grants before execution. Calibrate larger-bound work
charges and memory from the actual algorithm; do not make an accounting cap
silently bind before the intended time allowance. A valid candidate can be
compared even if an unrelated reference fails. This review launches no runs.

#### B4a — coarse native arithmetic and whole-portfolio costs

CFFI supplies a boundary, not an arithmetic algorithm or automatic SIMD.
Current optional GMP values remain in their backend through Python operator
loops; a useful native challenger moves a bounded ladder/program chunk behind
one call, keeping scratch and modular state native. Avoid a call per modular
multiply, a Python callback per point operation, or casting interpreter-private
bigint layouts. Do not downcast modular values or combined exponents into C
words. Serialize canonically at checked boundaries and compare the
same program/coverage, including strict frontier/GCD recovery. A ladder-only
microcontrol does not replace the current eligible PRAC production control
for whole-factor comparisons; label any executor change. Native wrappers
around another factorer remain C11 reference tracks rather than this kernel.

Start with scalar compiled arithmetic and existing exact algorithms. A native
Montgomery/Barrett context is a later C5 candidate under a changed cost regime;
previous Python reducer losses remain historical controls. ECM's many point
operations make it the first selected boundary candidate. PM1 stage one already
combines a prime-power chunk into one modular pow; rho already uses Brent
batches and local variables. Both can improve, but merely adding their existing
optimizations or wrapping one pow/GCD is not a new gain.

Preflight the actual PyPy/CFFI/compiler/GMP build before selecting a kernel.
The reviewed standalone PyPy 3.11 binary finds CFFI and GMP headers/library;
it does not find gmpy2 in that interpreter's environment. Other study
environments have separate identities. No package was installed here. Do not
claim documentation-only native APIs, thread release or zero-copy behavior
were measured locally. Persistent native state requires exact serialization,
owned allocation/RSS bounds, explicit cleanup, refused allocation handling,
finite atomic latency and cumulative work/CPU charges. Separate funded
preflight refusal from native allocator/process death: the latter recovers
from a last committed checkpoint, not an assumed native rollback. GC pressure
hints do not enforce a memory quota. B12/F4 may test native threads after verifying
thread-private state, release behavior and shared-ledger CPU semantics.

C11's initial cost atlas includes ordinary/easy exits, rough survivors,
prime-heavy children, smooth p±1, small-factor rho, close/power/special-form
controls and balanced hard inputs. PM1 is property-specific rather than
uniformly low potential: smooth p−1 can make even a large balanced factor
accessible. Repeated PM1 bases are correlated. Rho's work grows roughly with
the square root of the factor; tune marginal effort without promising it
will replace ECM/QS on large balanced composites.

E3/A10a/E4/E5/E6 may design and screen a justified isolated delta against
current committed contracts before deferred E1/D7; integrated acceptance
retains D7. Recoverable completed trial coverage can pass to child divisors,
whereas parent compositeness and non-power results cannot. If all primes up to
B were excluded for the tested n, a nontrivial power n = a**k requires
`(B + 1)**k <= n`; use exact comparisons to prune exponents. Valuation-gcd
pruning applies only to the same tested integer, not a cofactor after those
small factors were removed. Validate and charge restored coverage evidence. Trial-product
or perfect-power algorithms must beat current screens on their actual size/
reuse class, including setup. Keep native witness batching's early-rejection
waste, requested MR rounds, certainty and RNG prefix visible. Context/cache,
polling, validation and checkpoint changes need observed whole-cost headroom;
do not weaken exactness or quietly change cancellation responsiveness.

The broad review does not restart the stopped backend sweep, select every
cheap-path candidate or authorize a campaign. Follow the existing stage gates,
clean feasibility, observer calibration, finite cumulative grants and fresh
confirmation rules for whichever single hypothesis is selected next.

#### E6 — self-contained native prime generation

Use the existing C11 wheel-30 segmented implementation as the first native
prime-provider candidate. Bundle a pinned source adaptation, header and CFFI
adapter in Factor when selected; runtime use must not require a sibling
repository or external prime-sieve library. Retain the current Python stream.
This direct native boundary is distinct from earlier Python ports and their
recorded losses; its current complete-factor benefit remains unmeasured.

Return bounded packed batches with explicit status and capacity. Adapt the
inclusive C range endpoint to Factor's half-open contract, distinguish empty
output from allocation/range failure, preserve integer widths and fund native
workspace plus output/copy coexistence. Size reusable buffers for the selected
packet, not a broad materialization workload. Check cancellation and cumulative
work at bounded boundaries and resume from canonical checked positions.
Contexts are private to one caller. Use one inner thread inside ECM/QS workers
unless a separately funded concurrency plan shows a benefit; the standalone
CLI's adaptive OpenMP policy is not a library or Factor policy. C11 attributes
generation, repeated setup and consumption before raising this conditional
item above the current worker and factoring-kernel priorities.

#### Bulk kernels: PyPy, NumPy and optional native acceleration

B5/B5a/A1c share the exact current kernel/capture oracle and stay separate
implementation tracks. The reviewed host reports Apple M4/arm64 with NEON,
SME and SME2; x86 AVX/AVX-512 implementations are algorithm/layout references,
not runnable acceleration for this host. The inspected PyPy AArch64 JIT lacks
the vector backend found in its x86/zarch sources. Do not claim that changing
a JIT flag will vectorize these sieve loops, or that Python bigints become
packed SIMD lanes. ARM features being present does not establish useful SME
compiler/OS integration or a suitable algorithm for it.

Select a cheap contiguous threshold/offset scan first when its current cost
share warrants it. The existing packed-word scan and sparse queue are controls;
this is not an untouched Python element loop. Keep original score values,
conservative thresholds, sorted offsets, duplicate/singular roots, exact hit
metadata, tails, saturated/wide scores, cancellation and charged replay.
Bulk marking is a separate next candidate: dense small-prime strides and
sparse irregular updates have different SIMD/memory behavior. Do not silently
replace arbitrary-width hit masks or exact scores with wrapping fixed dtypes.

B5a is a QS native-kernel task as well as an instruction-set investigation.
Choose one useful coarse boundary from current cost evidence:

| Candidate | Plausible benefit | Required comparison and integration boundary |
| --- | --- | --- |
| Contiguous scan and offset compaction | Regular loads/comparisons are a natural SIMD target | Compare the actual configured scan and available packed-word control; include list packing/copy, output conversion, thresholds and tails. A tiny scan share cannot justify a large project. |
| Batched modular root shifts | Contiguous root/prime/delta arrays permit lane-wise add/subtract and conditional reduction | Reuse A1a/A1d/C10a tables; prove lane bounds, exceptional primes and exact roots. Repacking dataclasses every polynomial may erase the gain. |
| Sparse queue application or dense marking | Scalar compiled loops can remove Python event/update overhead even without SIMD | Keep exact repeated-offset accumulation, power weights and first-level factor hits. Irregular scatter is not automatically a vector win. |
| Fused marking/queue application/scan | Persistent buffers can amortize conversions across stages | A justified combined packet may proceed without each isolated kernel winning. Charge all sidecars, synchronization, unused preparation and native-to-Python recovery. |

The current sparse queue explicitly requires powers/list/bucket settings and
updates scores with ordinary list addition. Changing only the backend flag or
relaxing the eligibility guard is insufficient: fixed-width buffers need their
proved update/clipping rules, and exact factor-hit masks have arbitrary width.
Preserve the accepted queue or label and measure a combined scheduling change
against the actual production control. Use funded full-width bitsets or a
separately checked sparse/replayed hit sidecar; never truncate prime metadata.
A scan view may clip safely, while a persistent clipped score used by later
norm refinement needs its own threshold/coverage argument.

On M4, compare scalar C, compiler vectorization and NEON. Contiguous scans and
root shifts are different from irregular score writes. Duplicate destination
offsets make naive gather-add-scatter lose updates; serialize, combine or
prove conflict handling. Unsigned saturation, signed comparisons and wrapping
byte arithmetic are different semantics. Test threshold/width boundaries,
misalignment, tails and coincident roots against exact oracles.

AVX2/AVX-512 are optional x86 ports, not M4 accelerators. Keep a baseline
translation unit and verify the needed CPU/OS feature subsets before dispatch;
one AVX flag does not imply all AVX-512 byte/conflict operations. SIMDe or a
similar portability layer is optional source reuse, with generated code and
whole-cost behavior checked on each target. Do not extrapolate vector-width
ratios, run x86 assembly as an ARM benchmark or require an unavailable x86
machine before measuring the valid M4 candidate. SME remains later and
separately justified. Packed matrix/SpMV kernels retain A8/B6/D1/D2 ownership.

For B5, verify an actually usable NumPy build on PyPy 3.11/arm64. Its native
loops may use ARM SIMD despite cpyext call overhead; include import, reusable
array/buffer views, dtype width, copies, temporary masks and output conversion.
No one-call-per-prime design or assumption that every operation is vectorized.

For B5a, select a coarse CFFI buffer ABI and compare compiled scalar C,
compiler-vectorized C and justified explicit NEON. Inspect compiler reports
or generated instructions to establish vectorization. Price buffer pinning,
copies, FFI, build/startup, scratch, cleanup and bounded interruption; an unsafe
overflowing kernel is not a control. Reuse exact Python recovery and validation.
An optional kernel prototype does not authorize replacing the complete engine
with C or changing dependency-free defaults. Match backend/worker resources
and report whole-factor consequences before any adoption decision. GPU/SME
and hand-written multiprecision SIMD remain later profile-driven investigations.

A7a owns a conditional native-thread application of this ABI to QS collection.
Make selected native packets reentrant with private mutable scratch, immutable
shared tables, explicit buffer lifetimes and bounded interruption; verify GIL
release and avoid Python API access inside concurrent native work. Price the
Python recovery/parent bottleneck. B5a supplies the kernel contract; A7a owns
leases, shared relation admission and complete-factor scaling. Neither task
requires a second pool implementation or waits for the other's unrelated path.

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
| I1 — P6.4 final comparative publication | Publish source/configuration pins, independent corpora/commands, backend/core disclosures, supported limits and complete/factor-one outcomes. Consume C11's pinned feasible Python competitor protocol, with distinct native-int/GMP/native-wrapper/multiworker classes, cold/warm costs and censoring. | H1; selected P6 decisions and committed-only validation; F1 only for later GNFS claims | **Open; Required for comparative claims.** C11 may prepare the protocol and bounded competitor feasibility earlier. Final fastest/speedup claims need matched complete cost, fresh independent-input confirmation and uncertainty; no claim follows from feature counts or native reference timings. | Sol / **xhigh** comparison/protocol design; optional separate Sol / **high** frozen execution and publication | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#execution-phase-i--reproducible-publication) |

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
| C10a-02 — actual-A population and inter-family diversity | **Conditional uncovered delta.** C10a-01 now owns the implemented reachable-window/paired A population and norm-preview experiments. Reconcile those results before proposing broader product/diversity changes; the old 128-A census and rejected ordering retain their scopes. | Require a distinct product/norm/overlap hypothesis and useful post-filter evidence. Preserve unique families, deterministic finite exhaustion and resume prefixes; do not repeat retired quality-order, face or quota grids. | Sol / **xhigh** a separately scoped selector/census design; Astra / **xhigh** for a separately selected new coupled-arithmetic proof; optional separate Sol / **high** frozen execution | [Gates and census](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#c10a--conditional-collection-policy-questions-open) |
| C10a-03 — conservative pre-norm screening and bounded threshold feedback | **Conditional.** Investigate cheaper exact pre-norm rejection on captures including rejected positions, or explicitly selected threshold feedback with a declared coverage/loss contract. Preserve exact refinement and scoped local-bound losses. | Current collector; C9a retains residual admission/target policy. High refined-rejection counts do not prove division dominates or authorize omissions; saturation and small-prime heuristics need independent coverage checks. | Sol / **xhigh** existing-bound design/implementation; Astra / **xhigh** a separately selected new rejection proof; optional separate Sol / **high** frozen execution | [Gates](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#c10a--conditional-collection-policy-questions-open) |
| C10a-06 — separate runner/accounting tranche | **Open.** Reconcile repeated lease-marker handling in the separate runner; C10a-05 fixed its own new runner only. Preserve frozen C10 sources and historical failures. | C10/C10a-05; new protocol/source identity for a changed runner. Runner repair establishes no collection-policy promotion. | Sol / **xhigh** accounting changes; Sol / **high** separately frozen execution | [Recorded open item](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#c10a--conditional-collection-policy-questions-open) |
| C10a-04/G1 — larger-size feasibility after delivered routing | **Conditional follow-up after delivered C10a-01 feasibility.** Consume its literal 60/62/65/70 results before scheduling duplicate larger-size probes. One exposed modulus gives feasibility, not general calibration or promoted routing. | Automatic routing still ends at 60. The historical 65 spike timed out under an older collector; unvalidated 61–80 interpolation is not a default. A valid 70 comparison does not require 65 completion; subsequent studies need explicit finite grants. | Sol / **xhigh** a separately selected design/configuration task; optional separate Sol / **high** frozen execution; historical C10a-01 used Astra / **xhigh** throughout | [65 spike](benchmarks/RESULTS.md#c10a-04-larger-band-routing); [retained scope](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#c10a--conditional-collection-policy-questions-open) |
| A1a — direct extended versus pre-integration comparison | Conditional on a useful distinct hypothesis and new freeze; independent studies cannot establish this comparison by chaining gains. | Delivered A1a stays Done; selected-power-period/unroll4 extension remains opt-in after its losses. | Sol / high frozen comparison; xhigh only for changed semantics | [Open gate](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#a1a--selected-collector-kernel-delta-delivered-31-october-2026) |


<a id="c10a-01-implementation-remaining-after-the-coverage-review"></a>

### C10a-01 — coupled polynomial, multiplier, A-population and root setup

Delivered scope: [C10a-01 completion](#c10a-01-completion). The original
arithmetic/ownership contract and terminal evidence limits remain in the
[preserved source](https://github.com/nishanth17/factor/blob/7a3ca71c435bd5c5967da5737af728aa5b9bfb01/v2/ROADMAP.md#c10a-01--coupled-polynomial-multiplier-a-population-and-root-setup).
C10a-02 owns a distinct population/diversity hypothesis; A1a transitions,
A1c kernels, A1d batching, C8 family scheduling and F7 intake keep their scopes.

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
| P3.7 | Conditional | B5/B5a | Measured remaining array bottleneck and usable PyPy dependency required; exact overflow/duplicate-hit and conversion controls. | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p37--optional-numpy-spike-for-bounded-siqs-sieve-arrays) |
| P3.8 | Open | A8/B6/B7/C4/D1/D2; C9a/C10a collection interfaces | Real original-relation fixtures and exact forward/transpose products support fill-aware filtering, bounded Python-int Four Russians/PLE and genuine packed sparse solvers. Preserve compact substitution/provenance and original lifting. A8 is required; a B6/C4 win is not. | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p38--evaluate-gf2-matrix-and-filtering-optimizations) |

### Phase 4 — Optimize ECM chains and arithmetic after coverage exists

[Original phase scope and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#phase-4--optimize-ecm-chains-and-arithmetic-after-coverage-exists)

| ID | Status | Execution owner / predecessors | Remaining scope and trigger | Full gates and evidence |
| --- | --- | --- | --- | --- |
| P4.2 | Conditional | B4/B4a/C11 | Bounded B4 is Done; B4a is a separate optional fused native boundary. Wider 64–1024-bit formula/reduction/normalization comparisons need a remaining arithmetic bottleneck, fresh whole-engine evidence and exact a24/width/checkpoint/nonunit controls. | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p42--compare-fused-and-normalized-montgomery-kernels) |
| P4.3 | Deferred | A2/A2a/C11 | Foundation and stopped study are Done. Further experiments stopped; independent 3–100-digit stage/size confirmation, selective GMP-helper policy and automatic backend selection remain unsettled. Resume only as a separately selected tranche. | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p43--compare-pypy-python-integers-and-optional-gmpmpz) |
| P4.4 | Conditional | C5 | Remaining reduction bottleneck after A2/B4; exact encoded-domain/unit-scaling, intermediate-width and recovery proofs before whole-loop comparison. | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p44--keep-barrettmontgomery-reducers-experimental-until-measured) |

### Phase 5 — Complementary methods and stronger continuations

[Original phase scope and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#phase-5--complementary-methods-and-stronger-continuations)

| ID | Status | Execution owner / predecessors | Remaining scope and trigger | Full gates and evidence |
| --- | --- | --- | --- | --- |
| P5.1 | Open | A5 then B8 | Binary Lucas stages 1/2, discriminant/parameter checks, exact powers, genuine stage two, bounded trials/saturation/resume; prove marginal portfolio value. | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p51--implement-williams-p1-with-binary-lucas-evaluation) |
| P5.2 | Conditional | C2/C3a/B3a/F3; G1 allocation | Existing classical and optional paired decisions remain controls. Joint bounds/count, descendant search credit, selected higher-bound chains, old-point extension and polynomial continuation have distinct proofs; component feasibility does not wait for E1. | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p52--pair-ecm-stage-2-primes-and-tune-table-size) |
| P5.3 | Open | B8/F8; G1 allocation | A6 p−1 tranche is Done. Production-bound p−1 B2/B1 and gap/cache calibration remain with C11/G1; Lucas recurrence, p+1 extension and ECM old-point migration need independent composition/ratio/checkpoint proofs. | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p53--improve-p1-powering-and-continuation-independently) |
| P5.4 | Open follow-ups; C9a delivered stream Done | C9b with C10a/B7; C1a; A8/C4 | C1/C9 and the integrated C9a stream remain Done in their scopes. C9b owns current-mainline timing/headroom analysis and fresh performance confirmation; broader matrix/storage and higher arity retain their owners. C9a policies remain opt-in and its experimental gates remain unmet. | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p54--add-double-large-prime-siqs-and-stronger-filtering) |

### Phase 6 — Research options, parallel execution, and claims

[Original phase scope and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#phase-6--research-options-parallel-execution-and-claims)

| ID | Status | Execution owner / predecessors | Remaining scope and trigger | Full gates and evidence |
| --- | --- | --- | --- | --- |
| P6.1 | Conditional | F2 | Whole curve-family/stage-one cost/yield, nonunit maps and chosen-point order; isolated feasibility can precede E1, which remains for combined acceptance. | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p61--compare-edwardswindowed-and-torsion-aware-ecm) |
| P6.2 | Open selected F3/F7 design; other implementations Conditional | F3/F5/F6/F7/F8/F9 | F3 owns exact polynomial ECM kernels/trees; F8 owns distinct fast P±1; F9 is the conditional genus-two branch. F7 and C9a own their SIQS/batch domains; F6 higher arity has its own proof. Isolated F2/F3/F4/F6/F7 feasibility uses current controls without E1; later combined acceptance and F5 retain their gates. | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p62--investigate-polynomial-continuations-and-richer-relation-collectors) |
| P6.3 | Open preparation; conditional adoption | B12/A7a/F4 | Independent ECM and current-collector QS worker workstreams can prepare before E1. Reuse A7; settle growing-DLP/current-policy compatibility for QS. Compare 1/2/4 workers, then justified wider M4 occupancy, with total CPU, aggregate memory and first/full-factor wall time; later integrated acceptance retains E1. | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p63--add-reproducible-process-level-parallelism) |
| P6.4 | Open | C11 protocol preparation; I1 final publication | Prepare comparable pinned Python competitor feasibility now, separating runtime/backend/worker classes and complete versus factor-one outcomes. Final claims retain H1, selected P6 decisions, committed-only validation and fresh matched confirmation. | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p64--publish-reproducible-workload-specific-results) |

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
| P8.3 | Open | E4 | Current committed control for isolated screening; D7 for integrated tuning. Match total grants/seeds, full recovery/resume and marginal complete value. | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p83--calibrate-bounded-rho-batches-walks-and-restarts) |
| P8.4 | Open | G1/C3a; C3 Done | D7/E1 and selected decisions; charge prior work and fit only observable metadata, with matched fresh confirmation. | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p84--tune-ecmp1-allocation-and-relation-engine-handoff) |
| P8.5 | Conditional | E5 | Current committed contracts and observed cost for isolated screening; D7 for integrated tuning. Preserve deadlines, saturation recovery, atomic commits and charged replay. | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p85--measure-recovery-polling-and-checkpoint-costs-separately) |
| P8.6 | Conditional | E6/C11/A2a/B4a | Current contracts for isolated context/boundary screening; D7 for integrated tuning. Preserve prior cache/wheel losses and complete cold/warm costs. | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p86--revisit-context-setup-and-schedules-only-for-changed-workloads) |
| P8.7 | Open | H1 | G1 and E3–E6/selected P6 decisions; integrate confirmed winners and reconcile each historical P2 gate on untouched data. | [Gates and evidence](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p87--integrate-accepted-winners-and-close-the-reconciliation) |

### P3.8 carry-over and mathematical boundaries

| Existing ID | Status / next scope | Trigger, predecessors and owner | Evidence |
| --- | --- | --- | --- |
| P3.8-R1 | Open beyond completed B1 calibration: larger feasible calibration and independent corpus/runner coverage; reconcile C10a-01's implemented A/root/multiplier bundle before new population work. | C10a-01 retains its completed four-size feasibility; C10a-02 keeps uncovered product/diversity hypotheses. Coordinate C3a/G1 and stratified factor balance; no native crossover threshold or duplicate larger-size campaign. Historical C10a-01 used Astra / xhigh; Sol / xhigh separate selection. | [Original carry-over](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p38--evaluate-gf2-matrix-and-filtering-optimizations) |
| P3.8-R2 | Open A1d batching design; conditional B13/C8 capacity and whole-interval CRT. Bounded A1 integration stays Done. | A1d exact queue/family contract and current traversal-cost hypothesis: Sol / xhigh. B13/C8 require their distinct refusal/eligibility evidence: Astra / xhigh. | [R2 and deferrals](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p38--evaluate-gf2-matrix-and-filtering-optimizations) |
| P3.8-R3 | Delivered C9a cadence/storage Done; A8/B7/C4 filtering/provenance follow-ups remain Open. Old cadence-8/layout losses stay recorded. | C9b owns further attribution on the integrated collector. Use current preparation/filter/solve counters and useful progress; no automatic old-campaign repeat. Changed representations retain exact bounds, original-row lifting and simultaneous ownership. Existing owner assignments govern new mechanisms. | [R3 gates and rejected layouts](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p38--evaluate-gf2-matrix-and-filtering-optimizations) |
| P3.8-R4 | Conditional residual-classification/certification reuse after delivered C1; preserve certainty. | Profile must show repeated intake/atomic primality cost. A base-squared shortcut needs complete eligible-prime coverage including A/multiplier exceptions; opaque composite cofactors need a separate type/checkpoint proof. C9a/E3 coordinate; Astra / xhigh changed proof. | [R4 gates](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#p38--evaluate-gf2-matrix-and-filtering-optimizations) |
| P3.8-R5 | Deferred E1 final SSS/worker/current-control comparison after Done A7 reconciliation. | Worker flyer/large-interval compatibility and 40-digit SSS calibration need explicit decisions. Clarify upstream root license before further adaptation; use assignment IDs for seeds, include transport/conversion/loss/replay costs. Sol / high reconciliation; xhigh accounting. | [R5/E1 scope](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#a1c--same-modulus-two-root-interleaving-conditional); [pinned arms](benchmarks/inputs/controls/a7_r5_e1_arms.json) |

For matrix work, `M` has constraint rows and relation columns; `M d = 0`
selects original relations. A8 freezes orientation, transpose products, lifting
and genuine fixtures, including the independent bilinear-product identity.
B6 is the first dense-kernel comparison by investment preference; stronger
C4 filtering can proceed independently on A8's contract. D1/D2 use feasible
B6/C4 controls without requiring their success; capacity can justify sparse
work before solve CPU dominates. D2 has no D1 dependency. Four Russians must bound
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

<a id="c10a-01-completion"></a>

- **C10a-01** — 3 November 2026; **Done, integrated and adopted before fresh confirmation**; coupled policies at every fresh automatic SIQS size, with existing size presets/grants and portfolio bands. Exposed 60/62/65/70 completions retain the failed 65 canary and censored 70 control; general all-size performance remains unconfirmed. [Results](benchmarks/RESULTS.md#c10a-01-coupled-polynomial-defaults).

<a id="c9a-completion"></a>

- **C9a** — 4 November 2026; **Done, integrated; policies remain opt-in**; exact complete60/70 observations are exploratory, with one exposed modulus per size and no fresh confirmation. Historical timing flags remain recorded; future timing judgments use calibrated variability rather than the uncalibrated 15% veto. C9b owns further validation/headroom. [Results](benchmarks/RESULTS.md#c9a-preparation).

## Documentation history and local evidence

The [pre-cleanup roadmap](https://github.com/nishanth17/factor/blob/7a3ca71c435bd5c5967da5737af728aa5b9bfb01/v2/ROADMAP.md) and
[status](https://github.com/nishanth17/factor/blob/7a3ca71c435bd5c5967da5737af728aa5b9bfb01/v2/STATUS.md) preserve detailed results,
chronology and prior instructions. Git history holds the exact committed bytes;
no duplicate historical report tree is published. Pending A7a/B5a additions
are retained in the live contracts above.

Former C9a worktree roots are recorded in local
`v2/audit/c9a/archive-20261104/relocations.json`; original preparation files are
under `v2/audit/c9a/archive-20261104/base/c9a/`. Final local evidence starts at
`v2/audit/c9a/optimization-83-89.md`. C10a-01 evidence is retained under
`v2/audit/archives/c10a01-20261103/`. The [public results](benchmarks/RESULTS.md)
retain concise numerical findings and their limitations.

## Phase checkpoints

[Phase exit gates](POLICY.md#phase-checkpoints) live with the stable acceptance
rules. The [original gates](https://github.com/nishanth17/factor/blob/4a946b8d6294d51ac5d9c36f61880cc63d5b59fd/v2/ROADMAP.md#phase-checkpoints)
retain their historical context.

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

- [x] C9a — Useful-yield, admission and storage stream (Done; integrated) — [scope and evidence limits](#c9a--combined-useful-yield-and-storage-workstream).
- [ ] C9b — Mainline validation and remaining headroom — [scope](#c9b--mainline-validation-and-remaining-headroom).

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
