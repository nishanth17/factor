# Factor benchmarks

Run benchmarks with PyPy implementing Python 3.11, from the repository root.
Use `make -C v2 benchmark` for the introductory suite. Every accepted timing
requires validated warmup, at least nine observations and output validation;
see [contributor guidance](../../AGENTS.md).

At explicit user direction, [combined face/queue defaults are adopted](
../studies/c10a07_default_promotion.md) after production integration and committed
correctness/QA checks. The timing evidence below remains exploratory; the
optional18-call confirmation was not run. Automatic SIQS jobs use the swapped
coordinate policy and exact queues where eligible, with bounded exact fallback.
No new61–80 presets or portfolio stages are enabled.

The [completed C10a-07 panel](../studies/c10a07_panel_results.md) retains42 timed
calls at45/50/55/60, preserving every historical freeze. Nine fresh active inputs
complete in every arm. Observed equal-band wall cost falls22.4% for queue and
23.1% for face+queue, with only0.86% incremental combined benefit. Missing parent
end drift keeps these exploratory; active continuation drift passes, while45
is unresolved and its drift fails. The original report nominated queue alone;
the later combined adoption is a separate user-directed decision.
The [numerical summary](../studies/c10a07_panel_results.json) records exact costs,
uncertainty, resource accounts and the verified26.3 MB archive. Committed-only
857-test QA (nine skipped), lint and215 benchmark imports pass.

The [C10a-07 recovery report](../studies/c10a07_recovery_results.md) retains the
stopped historical quota/face diagnostics and the completed 24.7-minute cold
historical control. Its new current-mainline calibration uses short exact-work
prefixes, separately measured observer overhead and directly admitted face arms.
At preparation it nominated the existing A1a fused collector at60; local norm
bounds and an8 MiB power cache did not advance. Its original screen remains a
preserved preparation artifact, superseded for execution by the completed panel
above. It is not an instruction to run another comparison.

The separate [source-gap follow-up](../studies/c10a07_source_gap_results.md)
adapts exact sparse block scheduling from the source comparison. Nine matched
prefix pairs on three exposed inputs show 40.4% lower median prefix cost against
the existing fused collector, with identical exact output/coverage counters and
an explicit separate helper reserve. This is not a complete-factoring gain.
Its original six-call complete-cost screen remains frozen and unstarted; the
consolidated panel superseded its execution. Do not add it to the completed
campaign. Production integration and clean confirmation remain separate gates.

The [consolidated collection assessment](../studies/c10a07_research_assessment.md)
records 42 source/code/evidence dispositions and the resulting Phase 3 order.
The panel's original recommendation is superseded by the combined adoption
above. The optional [18-call design](../studies/c10a07_recovery_followup.md)
is retained for future validation; it was not run, is not required for closure
and schedules no new campaign. Keep normalized polynomial forms,
cross-B power roots and residual-yield policies distinct; avoid reimplementing
existing polling, fixed scores or batch arithmetic. Those ideas remain backlog
and cannot delay the current study's report. This closes the scoped
research assessment only. No new performance claim or campaign was added.

The [A1a integration protocol](qs/a1a_integration/protocol.md) compares the
selected density+presieve+scan_fused port with committed pre-integration
mainline at exactly 50 and 55 digits. Preparation includes an immutable owned
baseline, independent arithmetic/resume checks and a fresh-input freeze;
the [frozen execution report](../studies/a1a_frozen_results.md) accepts a scoped
6.008% complete-cost reduction at 50 digits (99.5% interval: 4.373–7.795%
lower). The 55-digit result remains inconclusive after its one full stability
extension because repeated-cell drift persists. Subsequent
[user-directed adoption](../studies/a1a_default_promotion.md) selects the
integration at native 50–55 digits, inclusive, at 1 GiB owned memory; 51–54
are unmeasured extrapolations and 55 retains its inconclusive timing label.
The completed isolated study's 5.945%/3.214% cost reductions do not
establish additive gains. Forty digits retains control; clean
60-digit results remain inconclusive and the 65-digit pilot showed no
completion within its allowance. Replay that original study from `d61df6c`,
preserving its controls and source hashes. E1 remains deferred.

The [bounded A1a extension](qs/a1a_extension/protocol.md) separately compares
selected-power periods plus four-stride per-root marking with the exact prior
A1a implementation. Its original integration control, corpus and receipts
remain unchanged. Its separately excluded fresh panel completed all calls,
but the combined extension cost 5.285% more at 50 digits and 2.994% more at
55 digits (99.5% intervals: 4.008–6.245% and 2.432–3.501% higher). Both bands
retain the prior integrated control. This tests the combined extension;
individual-component effects and a direct extended-versus-mainline comparison
remain unmeasured. The extension stays opt-in; E1 remains deferred.

The separate [C10a-04 bounded follow-up](../studies/c10a04_short.md)
compares compact with control at 60 digits, using fresh inputs and a ten-hour
total deadline. Fixed early trend stops are preliminary below nine samples;
they do not satisfy the ordinary accepted timing gate or promote defaults.
The original multi-band campaign remains incomplete and its evidence intact.
The [completed preliminary result](../studies/c10a04_short_results.md) stopped
after six pairs: compact completed 6/6 versus 5/6 and used 36.16% less total
capped time, with 31.72–40.08% saving per input. The subsequent
[user-directed promotion](../studies/c10a04_default_promotion.md) enables
50–60-digit routing, preserving 50/51 settings and documenting unstable 55
training and unmeasured intermediate sizes.
The [65-digit feasibility spike](../studies/c10a04_spike65_results.md) timed
out in both arms at 1800 seconds; no routing above 60 is promoted.

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
| qs/c10/ | Nine-source collection study, exact polynomial oracles, scoped 40/50-digit plan acceptance and C10a follow-ups |
| qs/c10a05/ | Fixed collection/solve pacing, six-source reconciliation, finite resume oracles and scoped 50-digit complete-cost acceptance |
| qs/c10a05_plus/ | Frozen 51-digit default-scope extension, 55/60-digit feasibility limits and charged fresh completion evidence |
| [qs/a7/](qs/a7/) | Frozen comparison arms and reconciliation |
| [primality/a10/](primality/a10/) | Primality dispatch and validation |
| [infrastructure/](infrastructure/) | Arithmetic backends, parallelism, sieves and performance audits |
| [support/](support/) | Snapshot loaders and independent arithmetic oracles |

Runners and builders live beside their experiments. Detailed protocols,
research, comparisons and acceptance reports are preserved locally under
`v2/audit/research/`. Public summaries below retain measured scope and limits.

The C10a-05 report separates polynomial windows from
collection exhaustion. A single fixed challenger was compared against the
current C9/C10 plan: 63 training pairs and 27 fresh pairs, 180 validated
complete calls. Fresh 50-digit cost fell 25.99% (95% interval 22.87–28.24%)
over nine matched pairs on three independent certified inputs. Only that
native 1 GiB bundle adopts collection pacing; 40 digits remains on window
pacing after its unstable eighteen-pair extension. The report includes extra
collection, downstream stages, setup, work, storage and charged replay costs.
The 6244.182-second inclusive study charge includes prepaid feasibility/QA;
cold probes and censored profiles are separate. Reproduce the frozen runner
at `96e2aad`; current integration intentionally changes its strict source pin.
Other C10a items, C9a feedback and E1 remain open. These percentages are not
additive with earlier C1/C9/C10 studies.

The larger-band follow-up checks the user's
requested 50+ routing. Nine matched training groups on three 51-digit
products gave 9/9 complete candidates for C9/window, C10/window and
collection pacing, versus 0/9 for the actual 51-digit default. The direct
C10/window-to-collection complete-cost reduction was 28.19% in training
(95% input-cluster interval 27.43–29.07%). After committing selection,
nine fresh matched calls on three newly generated products were 9/9 complete
with the selected bundle versus 0/9 with the previous default. This is a
completion gain, not a speedup against successful baseline factoring.
Inclusive wall/CPU charge was 12,199.720/12,188.535 s, including 3,500 s
prepaid feasibility/QA; a six-group training lease resumed without resetting
charges. The 55-digit one-number pilot finished under a longer 1,200 s grant,
while the 60-digit pilot was censored at that grant. Only exact 51-digit native
routing advances; 52+ remains unchanged pending size-specific calibration.

The C10 tranche 1 report records an isolated
2 MiB prime-power plan improvement on fresh complete factoring: 21.83% at
40 digits (95% interval 12.47–24.88%) and 17.79% at 50 digits (16.11–19.48%),
with nine matched validated calls per arm and three independent certified
inputs in each band. The 30-digit control stays selected. The tested 16-family
actual-A ordering did not advance, while normalized polynomials, independent
family diversity and per-job threshold feedback retain explicit C10a gates.
The 50-digit training extension stopped at 11/18 groups under its frozen
whole-group reservation; the fresh 9/9 result, not that partial training
interval, supports the scoped decision. Cold starts, profiles, work and
memory are reported separately. The C10 default decision
routes the exact 40/50-digit native bundles and plan automatically at the
confirmed memory allowance. The six conditional follow-ups were assigned to C10a; this
closes the delivered C10 tranche without claiming the whole size-policy space.
The subsequent C10a-05 decision above resolves only fixed collection pacing.

The subsequent C9 default promotion makes its
confirmed native 50-digit bundle automatic at the user's direction, preserving
explicit limits/settings and saved configurations. Its finite omitted grant
uses the confirmed resource window; this adds no new timing evidence.
The later C10 selector composes its accepted plan with that 50-digit bundle;
the C9 control and its default-promotion report remain historical evidence.

The C9 integration record documents C3/C9 resume
compatibility and preserved study freezes. Remaining literature transfers are
phased into roadmap C9a; E1 and C10a remain separate.

C9.2's adaptive-storage acceptance confirms an
explicit balanced 50-digit bundle: 38.10% less complete-call time than calibrated
SLP and 53.08% less than fixed DLP, with 9/9 valid completions in all arms.
Peak owned reservation rises to 326 MiB (SLP 167, fixed DLP 198). Three
independent inputs and three seeds per input support this narrow scope;
30/40 historical recommendations, C10a calibration and E1 combined acceptance
stay separate from the C9 result.

## C3 continuing optimization

Additional primary research and
diagnostic findings expand the search beyond quick8.
A frozen 84-call instrumented pilot validates all outcomes and finds two
extra-curve hits that avoid SIQS. Matched profiles identify redundant
reservation checks; their candidate reduction passes 600 tests and lint.
The offline stopping fit subsequently completes all
2,052 validated instrumented calibration calls in 2,893.780 active seconds.
It selects 64 curves at 2000/50000 for 30 digits and 64 at 2000/147396 for
40 digits, with nineteen independent root-ECM entrants per band. The table
and every calibration-capture hash are committed before comparison. These
are training projections; no warmed policy promotion follows from the fit.

The complete-call procedure froze
57 revealed subjects, nine seeds and three arms. It was stopped after the user
questioned the 513-group training scope: 72 complete matched groups (216 valid
calls) are preserved, with no comparison decision. Its interrupted group grant
remains charged; the protocol is permanently inconclusive and cannot resume.
See the scope decision.

The replacement upstream transfer screen
uses six revealed subjects, nine seeds and six arms: 54 groups, 324 calls.
It tests the protected source presets, retaining the
exact SIQS settings and admission bounds. Promising unstable results may extend
to 18/27 samples only within the separately fixed phase cap. This screening
cannot establish a broad default. A selected candidate still needs disjoint
fresh confirmation. See model limits and the
prepared certificate mechanism.
The initial source-screen results complete
all 324 calls. `escalating140` and `ladder115` show 44.44%/47.39% lower
training CPU cost than fixed32, but uncertainty remains wide. After the user raised
a time concern, defer the prescribed extension and fresh confirmation, select no
policy and retain defaults. Final committed-only checks pass 653 tests, lint,
154 benchmark imports and 35 certified corpora.
The user subsequently requested default promotion:
use compact64 at 30 digits and the escalating prefix at 40 digits,
subject to existing workspace admission. This separate decision overrides the
earlier defer decision; it does not alter the frozen assessment's null selection,
pass its stability gate, or claim a measured hybrid gain. No new timing search
or extension was initiated. Committed-only verification of the new defaults
passes 663 tests, lint, 154 imports, 35 certified corpora and 57 training inputs.
No round-two fresh inputs exist. Frozen earlier experiments use their source
commits; historical manifests are not repinned to newer production sources.

## C3 bounded ECM allocation and handoff (10 October 2026)

Research, protocol and
acceptance document the new opt-in cumulative pretest,
finite campaign and protected SIQS/SSS transition. Six bundles train on ten
historical inputs; the selected policy is committed before 15 fresh certified
inputs are generated. All 972 fresh calls validate and complete under matched
10^13-work, 5/30-second, 288-MiB native service allowances.

Retain numerical defaults. Training-selected quick8 costs 51.56% more on the
smaller fresh cohort (cost-change 95% interval −5.36% to +202.05%); uneven
controls regress and seven cells remain unstable after the frozen extensions.
The larger band retains control at training. Explicit reservations do not make
the default 2M grant sufficient: four of five smaller service-policy probes
refuse unfunded fallback floors. Charged active/terminal restore checks pass.

Working/committed-only checks pass 593 tests, full lint, 144 benchmark imports
and 34 certified corpora. Required inputs are versioned; raw receipts remain
ignored and preserved locally. Reproduction commands and all memory, seed,
small-population and source-identity limits are in the acceptance report.
C3 broad calibration, G1 and E1 remain open; no collector or v1 change.

## Inputs and generated evidence

`inputs/{baselines,controls,corpora}/` contains required numerical inputs,
controls and immutable baselines. Keep those historical bytes and hashes unchanged.
External source regeneration requires local files under folder-ignored
`inputs/upstream/`; those files, notices and licenses are not part of a fresh
checkout. Historical protocols still require their original pinned sources.

Write captures, profiles, stdout, verification dumps and temporary experiment
output under `results/<subject>/<experiment>/<run>/`. The entire `results/`
folder is ignored. Optional historical material belongs in ignored `history/`;
local journals and personal scratch work belong in `../local/`.

## Structure migration

Detailed migration verification is retained in the local research archive.

The folder migration changes import and command paths; it does not supply new
performance evidence. For example:

```sh
pypy3 -m v2.benchmarks.suites.phase_one --help
pypy3 -m v2.benchmarks.ecm.p41.p41_campaign --help
pypy3 -m v2.benchmarks.pm1.a6.a6_pm1 --help
```

[The migration manifest](inputs/controls/layout_migration.json) records old source
paths/hashes and exact relocated bytes. Source-pin validation in the migrated studies
accepts only explicitly recorded structural migrations. Other historical studies
retain their original freeze requirements: use their pinned checkout to repeat
an old measurement, or create a new protocol/freeze for current-source timings.
New captures always identify the actual current source bytes. Frozen source
loaders retain their private historical module names; mixed arms adapt only live
imports to that private topology.

### ECM preset naming maintenance (11 October 2026)

The active helper uses descriptive schedule names (see `../README.md`).
Frozen C3 runners and evidence keep their original keys and source hashes;
replay them from their pinned checkout. The new
`inputs/controls/a7_r5_preset_names_sources.json` chains the prior A7 manifest
for the unchanged numerical defaults and renamed preset API. No new performance
claim or measurement accompanies this naming change.

## C10a-04 historical freeze

The original 55/60/65-digit campaign remains paused and incomplete. Its
frozen sources and protocols are retained in Git history at `225634b`, with
detailed research and every attempt preserved locally. The separate 60-digit
result informed the user-directed 50–60 default promotion; the 65-digit spike
established no completed-factor evidence at its 1800-second grant. See the
linked result and promotion summaries above. E1 remains deferred.


### C10a-01 normalized polynomials: calibration and design handoff

Current complete-factor outcomes and the retained 65-digit host-confound
diagnosis are summarized in [benchmark decisions](RESULTS.md). The
[host-audited recheck protocol](qs/c10a01/regression65_recheck_protocol.md)
preserves the original observations, records power/runtime state and a native
canary, and admits a separate finite pair. The earlier unrun statements below
describe their original source freezes; they are historical, not current
campaign status. No defaults are promoted.

The [current v4 short-calibration result](../studies/c10a01_v4_calibration.md)
records 38 completed fixed-work cells, the frozen exploratory nominee, all
charges and a verified but unrun 60/62/65/70 complete-factor pilot runner.
Its 70-digit 30.340-second result is a prefix, not a factorization. The
50/55 lower-band supplement remains unrun pending Astra xhigh design review.

The [results and commands](../studies/c10a01_results.md) report 16 scored
fixed-position prefixes, four forced-policy diagnostics, four observer pairs
and two capacity fixtures. These exploratory observations establish no speedup.
The opt-in arithmetic uses current combined-policy control; all eight exposed
inputs, certificates, literal configurations and numerical allowances are
versioned. Detailed sources, research and raw receipts remain local.

The later `qs.c10a01.pilots` runner provides readiness checks and separately
funded complete-factor stages at 60/62, 65 and 70 digits. They were not run in
preparation. Stage C uses an explicitly forced eligible multiplier because the
balanced scored-selector cell otherwise repeats standard fallback. The
prospective v2 protocol isolates each cell in a fresh worker and explicitly
suballocates parent memory from the unchanged matched tiers. Native memory
censoring preserves the next cell. The original unrun freeze remains intact.
No automatic retry, memory escalation, default promotion or
untouched-confirmation claim.


The [combined v3 implementation](../studies/c10a01_coupled_implementation.md)
adds wider coupled ranking/residue reuse, paired A completion and guaranteed
powers-of-two initialization. `configuration_v3.py` constructs literal matched
control/challenger configurations. This bundle is correctness-validated but
unmeasured. Earlier calibration and source freezes remain historical; the old
pilot loader rejects changed sources. No new performance call is admitted by
the configuration helper or compatibility-manifest update.


`configuration_v4.py` adds centered root tables, batch table construction,
windowed-pair A selection and bounded actual-norm ranking. It retains v3 and
literal current controls and exposes component overrides for targeted
ablations. Structural QA covers all eight exposed 60/62/65/70 inputs; it is
not a timing result. The original calibration allocation is closed. No revised
calibration or later complete-factor pilot ran for this source; old frozen
pilot readiness continues to reject changed source. Freeze new source/work
accounts before measurements, including the separately funded 70-digit stage.

The v4 implementation passed 869 tests (nine optional-backend skips), lint and
206 benchmark/study imports in a committed-files-only export. Verification
included rejecting the old pilot source freeze; it executed no performance
calls. Detailed implementation/ownership proofs and retained QA failures are
local audit evidence.


The [v4 experiment design](qs/c10a01/protocol_v4.md) proposes a separate
3400-second wall/CPU calibration: clean full-arm/control prefixes first,
then targeted component diagnostics, followed by one frozen nominee. Its
machine-readable plan preserves the 60/62/65/70 complete-factor grants and
makes the revised whole ceiling 38600 seconds, including the original
calibration envelope. No new experiment was run. Validate the committed
numerical/source design with:

    pypy3 -m v2.benchmarks.qs.c10a01.plan_v4 check

The checker is not an execution runner. Bounded worker/ledger adapters,
readiness evidence and a post-calibration nominee release remain explicit
prerequisites. The owned mainline control imports without Git history or
archived worktrees; the old pilot source freeze remains unchanged.

Design commit `14caea4` passed 875 tests (nine optional-backend skips), lint,
208 benchmark/study imports and the numerical/source checker from a
committed-files-only export. These checks include all proposed literal
configurations and isolated-control loading; they ran no performance calls.

The active design now also includes the [lower-band supplement](qs/c10a01/protocol_lower_v4.md):
three exposed inputs each at 50/55 for complete-cost screening, plus short
51/52/55/56 boundary checks. Its separate 3700-second wall/CPU proposal raises
the whole envelope to 42300 seconds without reducing any 60–70 grant. It
preserves existing lower-band presieving and requires fresh per-band evidence
before claiming gains across 50–60. No experiments ran; pending boundary
inputs, nominee release and execution adapters remain explicit prerequisites.
