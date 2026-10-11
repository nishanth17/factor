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
| [qs/a7/](qs/a7/) | Frozen comparison arms and reconciliation |
| [primality/a10/](primality/a10/) | Primality dispatch and validation |
| [infrastructure/](infrastructure/) | Arithmetic backends, parallelism, sieves and performance audits |
| [support/](support/) | Snapshot loaders and independent arithmetic oracles |

Each experiment keeps its milestone identifier, runners, builders and research
notes together. The [study record](docs/studies.md) preserves detailed accepted
results, rejected experiments, commands and limitations.

The [cross-engine parameter audit](parameter_audit.md) reconciles constant
provenance and missing tuning evidence under roadmap C11. QS/SIQS and DLP keep
their [existing parameter audit](qs_parameter_audit.md) and C10/C9 ownership.
These ledgers distinguish exact invariants, finite caps and measured choices
from uncalibrated defaults; adding an audit does not change production settings.

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
