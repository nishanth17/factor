# Budgets, stages and allocation

This layer coordinates finite work across factoring stages and recursive
children. The public entry point is
[`factorize_bounded`](../docs/API.md#bounded-recursive-factoring) in
[portfolio.py](../portfolio.py).

| Module | Responsibility |
| --- | --- |
| [budget.py](budget.py) | Shared work, wall/CPU time and cooperative cancellation |
| [work_budget.py](work_budget.py) | Bounded work adapters |
| [schedules.py](schedules.py) | Prime contexts, exact powers and bounded reuse |
| [stage_jobs.py](stage_jobs.py) | Resumable rho/p−1/ECM stage state and recovery |
| [allocation.py](allocation.py) | Cumulative pretest ceilings and protected fallback admission |
| [ecm_defaults.py](ecm_defaults.py), [ecm_presets.py](ecm_presets.py) | Size selection and explicit numerical ECM tiers |

## One shared allowance

`Budget` carries `work_limit`, `seconds`, `cpu_seconds` and an optional
`cancelled` callback. Operations reserve work before mutation; refused work
leaves checked state available for continuation. Work units are versioned
algorithmic charges, not elapsed time or comparable bigint-operation counts.
Time/cancellation checks are cooperative.

Portfolio limits include preprocessing, attempts, recovery, child calls and
restore verification. Memory bounds conservatively cover simultaneous owned
workspace; process/JIT RSS is separate. Cached schedules are shared only
where their immutable integer identities match; curve/residue state stays
private. See [checkpoint rules](../docs/CHECKPOINTS.md).

## Pretesting and fallback

`ECMAllocation(mode="pretest", pretest_work=...)` caps optional search against
the original cumulative ledger. Earlier preprocessing, children and resumes
count toward that ceiling. Optional wall/CPU ceilings use `pretest_seconds`
and `pretest_cpu_seconds`. Mandatory classification/power checks still use
the shared total allowance.

`mode="campaign"` executes predeclared finite ECM tiers without a separate
pretest ceiling. Either mode can protect an enabled SIQS/SSS fallback with
`fallback_work`, `fallback_seconds` and `fallback_cpu_seconds`.
A floor is an admission reservation, not extra resources or a completion
estimate. An unfunded floor returns an explicit insufficient-fallback reason;
without a fallback, reaching a pretest ceiling returns `pretest_exhausted`.

```sh
pypy3 -m v2.factor 1000036000099 --siqs --ecm-policy pretest \
  --pretest-work 1000000 --fallback-work 5000000 \
  --work-limit 20000000 --seconds 30 --cpu-seconds 30 --memory-mib 80
```

A fresh cofactor hands off once; resumed admitted relation work does not need
a new full reservation. Partial attempts keep their consumed charges.

## Numerical tiers and saved policy

Omitted ECM tiers resolve once from the original root input. Explicit tuples,
including `()` and a fixed legacy schedule, take precedence; children and
checkpoints keep that concrete plan. Paired executors and insufficient
coexistence space preserve the applicable retained schedule.
[STATUS](../STATUS.md#current-automatic-defaults) owns the current routing table.

`with_ecm_preset(config, name)` replaces only tiers on an explicit protected
SIQS configuration. It requires positive work/wall/CPU fallback floors and
unchanged child storage that can coexist with the larger ECM workspace.
Names are `pretest8`, `escalating140`, `ladder115`, `ladder415`, `single74`
and `single214`; old upstream-named keys are rejected. Frozen studies keep
their original names and source revisions.

Explicit allocation checkpoints bind configuration and policy identity.
Implicit API resume restores the saved configuration; conflicting explicit
fields fail validation. New defaults do not reinterpret older saved jobs.
Detailed method-specific continuation lives in the [ECM](../ecm/README.md),
[p−1](../pm1/README.md) and [QS](../qs/README.md) guides.
