# Folder migration verification — 10 October 2026

The source tree groups algorithms under `ecm/`, `pm1/`, `rho/` and `qs/`,
shared arithmetic under `common/`, and budgets/schedules/dispatch under
`execution/`. Four Python modules remain directly in `v2/`. The tracked
benchmark root now contains only its README and package initializer; research
notes, runners and builders live beside their related experiments.

Imports and benchmark module commands changed. `v2.factor` and `v2.portfolio`
remain the factoring entry points. The [API guide](../../README.md#source-layout)
and [benchmark index](../README.md) describe the new paths.

## Validation

- PyPy 7.3.23 implementing Python 3.11.15: 568 tests passed, with eight optional
  GMP-backend skips, in a temporary Git checkout containing only versionable
  files. The checkout contained no local audit, history, result or environment
  folders. The existing process-isolation test required access to `ps`.
- `make -C v2 lint` passed Ruff checks, Ruff formatting and pycodestyle.
- All 140 benchmark modules/packages imported. Immutable snapshot loaders,
  mixed historical/current QS arms, chain catalogs, A7 prepared arms and B4
  source/input validation passed in the clean checkout.
- Standalone audit, follow-up and R1 regression drivers still select their
  runtime before importing its implementation. Their help entry points and
  the relocated suite, ECM and p−1 commands were checked.
- All 141 existing benchmark input files remain byte-identical. `v1/` is
  unchanged. Folder ignore rules cover `results/`, `history/` and `v2/local/`.

Run the normal checks from the repository root:

```sh
make -C v2 test
make -C v2 lint
```

## Reviewed structural differences

An AST comparison ignoring source locations found unchanged production code
once imports and the relocated ECM catalog path were normalized. The remaining
intentional research/test differences are module/worker command paths, stable
root resolution, recursive production-source enumeration, private snapshot
import adaptation, explicit migration-pin validation, and tests for those
contracts. A legacy phase-two snapshot fingerprint now uses the shared benchmark
root instead of constructing an incorrect nested input path. Scratch schedule
output now goes beneath the ignored results folder.

[The migration manifest](../inputs/controls/layout_migration.json) preserves
original source identities and records exact relocated bytes. It supplements
existing pins for explicitly supported structural migrations; frozen source
snapshots, certificates, protocols and historical timings were not rewritten.
New measurement captures identify current bytes. Historical studies with other
freeze requirements still require their pinned checkout or a new freeze.

This is structural verification, not new performance evidence. No arithmetic
algorithm, seed, factoring default, work allowance or checkpoint format was
changed. Raw verification output and the detailed AST review remain local under
`benchmarks/results/layout/2026-10-10/`.
