# Factor benchmarks

Run benchmarks with PyPy implementing Python 3.11, from the repository root.
Use `make -C v2 benchmark` for the introductory suite. Every accepted timing
requires validated warmup, at least nine observations and output validation;
see [contributor guidance](../../AGENTS.md).

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
| [qs/c10/](qs/c10/tranche1_results.md) | Nine-source collection study, exact polynomial oracles, scoped 40/50-digit plan acceptance and C10a follow-ups |
| [qs/c10a05/](qs/c10a05/results.md) | Fixed collection/solve pacing, six-source reconciliation, finite resume oracles and scoped 50-digit complete-cost acceptance |
| [qs/c10a05_plus/](qs/c10a05_plus/results.md) | Frozen 51-digit default-scope extension, 55/60-digit feasibility limits and charged fresh completion evidence |
| [qs/a7/](qs/a7/) | Frozen comparison arms and reconciliation |
| [primality/a10/](primality/a10/) | Primality dispatch and validation |
| [infrastructure/](infrastructure/) | Arithmetic backends, parallelism, sieves and performance audits |
| [support/](support/) | Snapshot loaders and independent arithmetic oracles |

Each experiment keeps its milestone identifier, runners, builders and research
notes together. The [study record](docs/studies.md) preserves detailed accepted
results, rejected experiments, commands and limitations.

The [cross-engine parameter audit](parameter_audit.md) reconciles constant
provenance and missing tuning evidence under roadmap C11. QS/SIQS and DLP keep
their [existing parameter audit](qs_parameter_audit.md) and C10a/C9a follow-up ownership.
These ledgers distinguish exact invariants, finite caps and measured choices
from uncalibrated defaults; adding an audit does not change production settings.

The [C10a-05 report](qs/c10a05/results.md) separates polynomial windows from
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

The [larger-band follow-up](qs/c10a05_plus/results.md) checks the user's
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

The [C10 tranche 1 report](qs/c10/tranche1_results.md) records an isolated
2 MiB prime-power plan improvement on fresh complete factoring: 21.83% at
40 digits (95% interval 12.47–24.88%) and 17.79% at 50 digits (16.11–19.48%),
with nine matched validated calls per arm and three independent certified
inputs in each band. The 30-digit control stays selected. The tested 16-family
actual-A ordering did not advance, while normalized polynomials, independent
family diversity and per-job threshold feedback retain explicit C10a gates.
The 50-digit training extension stopped at 11/18 groups under its frozen
whole-group reservation; the fresh 9/9 result, not that partial training
interval, supports the scoped decision. Cold starts, profiles, work and
memory are reported separately. The [C10 default decision](qs/c10/default_promotion.md)
routes the exact 40/50-digit native bundles and plan automatically at the
confirmed memory allowance. The six conditional follow-ups were assigned to C10a; this
closes the delivered C10 tranche without claiming the whole size-policy space.
The subsequent C10a-05 decision above resolves only fixed collection pacing.

The subsequent [C9 default promotion](qs/c9/c9_default_promotion.md) makes its
confirmed native 50-digit bundle automatic at the user's direction, preserving
explicit limits/settings and saved configurations. Its finite omitted grant
uses the confirmed resource window; this adds no new timing evidence.
The later C10 selector composes its accepted plan with that 50-digit bundle;
the C9 control and its default-promotion report remain historical evidence.

The [multi-source implementation gap audit](implementation_gap_audit.md)
extends the [Alpertron inventory](alpertron_gap_analysis.md), distinguishes
missing mechanisms from existing/planned work, and ranks follow-ups against
current evidence. Roadmap A1b separately reviews remaining R2 fixed-score combinations
against the accepted C10 plan control; it does not repeat that plan study.

The [C9 integration record](qs/c9/c9_integration.md) documents C3/C9 resume
compatibility and preserved study freezes. Remaining literature transfers are
phased into roadmap C9a; E1 and C10a remain separate.

C9.2's [adaptive-storage acceptance](qs/c9/c9_growth_results.md) confirms an
explicit balanced 50-digit bundle: 38.10% less complete-call time than calibrated
SLP and 53.08% less than fixed DLP, with 9/9 valid completions in all arms.
Peak owned reservation rises to 326 MiB (SLP 167, fixed DLP 198). Three
independent inputs and three seeds per input support this narrow scope;
30/40 historical recommendations, C10a calibration and E1 combined acceptance
stay separate from the C9 result.

## C3 continuing optimization

[Additional primary research](ecm/c3/round2_research.md) and
[diagnostic findings](ecm/c3/round2_findings.md) expand the search beyond quick8.
A frozen 84-call instrumented pilot validates all outcomes and finds two
extra-curve hits that avoid SIQS. Matched profiles identify redundant
reservation checks; their candidate reduction passes 600 tests and lint.
The [offline stopping fit](ecm/c3/round2_fitted.md) subsequently completes all
2,052 validated instrumented calibration calls in 2,893.780 active seconds.
It selects 64 curves at 2000/50000 for 30 digits and 64 at 2000/147396 for
40 digits, with nineteen independent root-ECM entrants per band. The table
and every calibration-capture hash are committed before comparison. These
are training projections; no warmed policy promotion follows from the fit.

The [complete-call procedure](ecm/c3/round2_comparison_protocol.md) froze
57 revealed subjects, nine seeds and three arms. It was stopped after the user
questioned the 513-group training scope: 72 complete matched groups (216 valid
calls) are preserved, with no comparison decision. Its interrupted group grant
remains charged; the protocol is permanently inconclusive and cannot resume.
See the [scope decision](ecm/c3/round2_comparison_scope.md).

The replacement [upstream transfer screen](ecm/c3/source_preset_protocol.md)
uses six revealed subjects, nine seeds and six arms: 54 groups, 324 calls.
It tests the protected [source presets](ecm/c3/source_presets.md), retaining the
exact SIQS settings and admission bounds. Promising unstable results may extend
to 18/27 samples only within the separately fixed phase cap. This screening
cannot establish a broad default. A selected candidate still needs disjoint
fresh confirmation. See [model limits](ecm/c3/round2_model_limits.md) and the
[prepared certificate mechanism](ecm/c3/round2_confirmation_input_design.md).
The [initial source-screen results](ecm/c3/source_preset_results.md) complete
all 324 calls. Yamaquasi escalation and Alpertron show 44.44%/47.39% lower
training CPU cost than fixed32, but uncertainty remains wide. After the user raised
a time concern, defer the prescribed extension and fresh confirmation, select no
policy and retain defaults. Final committed-only checks pass 653 tests, lint,
154 benchmark imports and 35 certified corpora.
The user subsequently requested [default promotion](ecm/c3/default_promotion.md):
use compact64 at 30 digits and the Yamaquasi escalating prefix at 40 digits,
subject to existing workspace admission. This separate decision overrides the
earlier defer decision; it does not alter the frozen assessment's null selection,
pass its stability gate, or claim a measured hybrid gain. No new timing search
or extension was initiated. Committed-only verification of the new defaults
passes 663 tests, lint, 154 imports, 35 certified corpora and 57 training inputs.
No round-two fresh inputs exist. Frozen earlier experiments use their source
commits; historical manifests are not repinned to newer production sources.

## C3 bounded ECM allocation and handoff (10 October 2026)

[Research](ecm/c3/research.md), [protocol](ecm/c3/protocol.md) and
[acceptance](ecm/c3/acceptance.md) document the new opt-in cumulative pretest,
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

`inputs/{baselines,controls,corpora,upstream}/` contains versioned inputs needed
to reproduce the experiments, including licenses and immutable source snapshots.
Keep those historical bytes and hashes unchanged.

Write captures, profiles, stdout, verification dumps and temporary experiment
output under `results/<subject>/<experiment>/<run>/`. The entire `results/`
folder is ignored. Optional historical material belongs in ignored `history/`;
local journals and personal scratch work belong in `../local/`.

## Structure migration

See the [verification record](docs/layout.md) for checks and reviewed differences.

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
