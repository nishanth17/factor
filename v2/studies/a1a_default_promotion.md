# A1a mainline adoption — 31 October 2026

Fresh native Python-integer jobs at exactly 50 through 55 decimal digits,
inclusive, select `density+presieve+scan_fused` with the default 1 GiB owned
allowance. This applies through `default_siqs_config`, `SIQSJob`, automatic
portfolio selection and the CLI. Explicit configurations, other owned-memory
allowances/backends and saved selectors keep their prior behavior. Existing
count/base, Gray quotas, finite work/time limits and serialized formats remain
unchanged. Forty digits and 56–60 digits keep their previous marking.

This adoption follows the
[frozen comparisons](a1a_frozen_results.md):

| Comparison | 50 digits | 55 digits |
| --- | --- | --- |
| Integration vs owned pre-integration mainline | Pass: 6.008% lower complete cost; 99.5% interval 4.373–7.795% lower | Inconclusive: all-data drift persists after the one full extension; descriptive estimate 6.291% lower |
| Combined extension vs owned prior A1a | Retain control: 5.285% higher complete cost | Retain control: 2.994% higher complete cost |

All 180 measured calls completed and reconstructed exactly. The 51–54-digit
policy has no separate timing confirmation, and adoption does not convert
the inconclusive 55-digit result into an accepted speed claim. No new
performance measurement, component-level inference, chained percentage or
portfolio/E1 claim accompanies this change. E1 remains deferred.

## Extension configuration

Selected-power periodic scoring and per-root four-stride marking remain
available through explicit configuration. They are not automatically selected:

```python
from dataclasses import replace
from v2.qs import SIQSJob, default_siqs_config

config = default_siqs_config(n)
config = replace(config, collector=replace(
    config.collector, marking="density+presieve_powers+unroll4+scan_fused",
))
job = SIQSJob(n, config=config)
```

Pass that same config as `PortfolioConfig(siqs=config)` for a portfolio call.
Checkpoints retain the concrete selector and rebuild derived state under
remaining charged allowances; an old sparse checkpoint remains sparse.

## Historical replay and cleanup

The experiments retain their original source, control, corpus and decision
identities in the execution report. Reproduce source verification from their
exact committed versions, rather than treating current defaults as a frozen
historical arm:

```sh
# From the main repository; use a new path for read-only source verification.
git worktree add --detach /private/tmp/a1a-integration-replay fee6a41c9af2702c34f964536e71ad3011d4cb22
# cwd: /private/tmp/a1a-integration-replay
/opt/homebrew/bin/pypy3 -B -m v2.benchmarks.qs.a1a_integration.study verify

# From the main repository; separate replay path and source identity.
git worktree add --detach /private/tmp/a1a-extension-replay 80d8cda1b11f2ca2724d428078c206930999e1ac
# cwd: /private/tmp/a1a-extension-replay
/opt/homebrew/bin/pypy3 -B -m v2.benchmarks.qs.a1a_extension.study verify
```

Verification creates no accepted timings. The original ledgers are terminal;
do not rerun generation, freeze, measurements or final inference against them.
Current A1a runners continue to reject changed source closures. A new A7 live
manifest records this promotion, retaining the old manifests and explicit A7
numerical arms. The historical C10a-04 template explicitly retains sparse
marking when reading the newer default factory.

Before removing task-owned worktrees and temporary QA checkouts, preserve
their useful receipts, research, journals and finite-account ledgers under
the main repository's ignored `v2/audit/a1a-archive/20261031/`. Its manifest
records original paths, revisions, archive hashes and verification. The
original isolated study at `d61df6c`, integration at `fee6a41`, extension at
`80d8cda` and report at `9afb19e` remain in Git. Task worktrees can therefore
be removed without discarding either source history or experimental evidence.

The delivered selected-kernel and bounded-extension item is closed in the
roadmap. Quota-aware setup, reusable switch metadata, lifted-power root reuse
and conditional same-modulus two-root interleaving retain separate future
gates and model/effort recommendations.

## Verification

The promoted source passes 823 PyPy tests with nine optional-backend skips,
plus `make -C v2 lint`. Boundary tests cover both signs, both ends of every
size band, other backends and adjacent owned-memory grants. CLI/portfolio
checks preserve the promoted, explicit extension and legacy sparse selectors
across resume. A committed-files-only checkout repeats the full tests and
lint, imports all 194 benchmark modules and verifies the strict A7 chain and
historical private adapters. Local QA receipts retain earlier charges.
