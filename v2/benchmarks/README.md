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
