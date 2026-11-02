# ECM experiments

Runners, builders, certified controls and historical executors for ECM.
Production code is in [v2/ecm](../../ecm/README.md); these research APIs do
not supply its shared-budget/checkpoint contract unless explicitly adapted.

| Directory | Study family |
| --- | --- |
| [p41/](p41/) | Verified PRAC arithmetic and complete ECM campaigns |
| [p52/](p52/) | Reusable programs, paired/wheel continuation and realistic portfolios |
| [b3/](b3/) | Production chain integration, coverage and strict recovery |
| [b4/](b4/) | Formula/reduction candidates and full-production bridges |
| [c3/](c3/) | Finite allocation, protected handoff and source-policy screens |
| [c6/](c6/) | Certified PRAC/Lucas/CF catalogs and reused-program executors |

Run modules through `pypy3 -m v2.benchmarks.ecm.<study>.<module>` from the
repository root. Read the selected runner and its frozen controls first.
For command syntax without a measurement:

```sh
pypy3 -m v2.benchmarks.ecm.p41.p41_campaign --help
```

C6's `c6_chains`, `c6_fast`, `c6_cf` and `c6_b4` modules expose experimental
record verification and bounded program/executor variants. Their current
imports include the `ecm.c6` package path. Caller-owned research programs lack
the portfolio's ledger and serialized resume; kernel or reuse results cannot
be treated as a complete-factor default claim.

Use the [benchmark workflow](../README.md) for finite plans, locks and capture
layout. [Results](../RESULTS.md) preserve decisions and scoped measurements;
[STATUS](../../STATUS.md) owns current automatic routing. Historical source
pins and prior losses remain unchanged by new experiments.
