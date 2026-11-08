# Command-line guide

H3 replaces internal refusal quotas with observed root limits. See the
[resource API](../README.md#resource-api) for defaults, removed arguments,
sampled RSS, stage handoff, saved-history accounting and optional atomic saves.


Run from the repository root with PyPy implementing Python 3.11.
[Getting started](../README.md) covers setup; `pypy3 -m v2.factor --help`
lists all current flags.

## Bounded factoring

```sh
pypy3 -m v2.factor 626100403 --bounded --seed 7 \
  --work-limit 2000000 --seconds 30 --cpu-seconds 30 --memory-mib 16
```

The bounded portfolio shares its allowances across preprocessing, attempts
and recursive children. `--seed` makes candidate assignments reproducible.
`--verbose` enables diagnostics; library calls are quiet by default.

Zero is rejected, one has an empty factorization, and negatives retain their
sign. Exit status is 0 for complete results, 1 for unresolved composites and
2 for invalid input. Completion does not imply every terminal prime is proven.

## Choose an engine

| Selection | Behavior |
| --- | --- |
| `--bounded` | Preprocessing, rho, p−1 and ECM, with configured relation fallback |
| `--method qs` | Preprocessing followed by the fixed-polynomial control |
| `--method mpqs` | Preprocessing followed by multiple-polynomial QS |
| `--method siqs` | Preprocessing followed by SIQS |
| `--siqs` | Add SIQS after rho/p−1/ECM |
| `--method sss` / `--method sssf` | Select experimental Smooth Subsum Search |

Explicit method selections also enable bounded execution. Direct relation
methods disable rho/p−1/ECM; all routes recursively validate returned splits.
Omitted native configurations enable SIQS in the calibrated bands documented
in [current status](../STATUS.md#current-automatic-defaults). That table also
records larger omitted resource grants and unmeasured routing scopes.

```sh
pypy3 -m v2.factor 10002200057 --method siqs --seed 7
pypy3 -m v2.factor 10002200057 --method qs --qs-base-bound 400 --seed 7
pypy3 -m v2.factor 10002200057 --method mpqs --qs-base-bound 400 --seed 7
pypy3 -m v2.factor 10002200057 --siqs --ecm-curves 2 --seed 7
```

## Resource and search options

`--work-limit`, `--seconds`, `--cpu-seconds` and `--memory-mib` set explicit
allowances. They take precedence over omitted limits. Time checks are
cooperative; the memory allowance bounds owned workspace, not total process
RSS. Increasing a deadline cannot restart an exhausted candidate schedule.

QS controls include `--qs-base-bound`, `--qs-half-width`, `--qs-max-half-width`,
`--qs-factor-count`, `--qs-family-count`, `--qs-pool-size`,
`--qs-assignment-policy`, `--qs-polynomials-per-family` and
`--qs-residual-bound`; corresponding `--siqs-*` aliases are accepted.
Streaming assignments and Gray quotas apply to SIQS. Wider intervals need a
compatible maximum width. See [QS configuration](../qs/CONFIGURATION.md).

Explicit DLP selection uses `--qs-dlp` with both `--qs-large-prime-bound` and
`--qs-large-product-bound`; it is restricted to SIQS. Optional candidate/split
limits narrow that finite search. Invalid combinations are rejected.

SSS controls include `--sss-base-bound` and `--sss-rounds`. SSSf's candidate
filter may lose candidates; it remains an explicit experimental policy.

## ECM policies and backends

`--ecm-curves` selects a curve count. Repeated `--ecm-tier B1,B2,curves`
flags declare finite independent tiers and cannot be combined with
`--ecm-curves`. `--ecm-policy pretest` needs `--pretest-work`; optional
`--fallback-work`, `--fallback-seconds` and `--fallback-cpu-seconds` protect
configured SIQS/SSS admission. `--ecm-policy campaign` uses the declared tiers.
These policies require the auto portfolio. See [execution](../execution/README.md).

`--ecm-chain auto|off|prac|lucas|cf` chooses bounded chain routing and implies
`--bounded`. Supported schedules, cofactor size and memory determine whether
a chain plan is used; [ECM](../ecm/README.md) describes the fallback rules.

`--backend python-int|gmpy2-mpz` explicitly selects arithmetic. GMP must be
installed in the chosen PyPy environment; no runtime or backend fallback occurs.

## Save and resume

Use `--checkpoint PATH` to save state and `--resume PATH` to restore it.
Resume grants are cumulative totals. Preserve the saved method and compatible
configuration; see [checkpoints](CHECKPOINTS.md) for examples and restore rules.

## Automatic ECM selection

The normal auto command now uses the bounded portfolio and automatic cost-gated
process ECM. `--bounded` remains accepted. Omitted tiers choose escalating140
at exactly 60/80/100 original digits; explicit tiers and saved configurations
win. Unsupported/unfunded regimes and cheap calls stay serial. Direct QS/SSS
methods keep their own execution path.

Use `--ecm-workers serial` for an explicit serial run, or leave the default
`--ecm-workers auto`. `--ecm-memory-mib` sets aggregate ECM ownership including
the parent and workers. Omitted ordinary memory permits up to 8 GiB; an explicit
`--memory-mib` workspace cap is also respected conservatively unless a separate
aggregate cap is granted. Worker count is selected by the packaged policy.

Fresh non-QS omitted grants are 30 wall seconds, 4,096 aggregate CPU seconds and
10^13 work units in both auto and serial mode. Explicit limits are never raised;
a small explicit CPU/work grant can select serial. Calibrated QS grants and
stock fallback/checkpoint settings remain. This is finite bounded execution,
not a guarantee of completion. The legacy `factorize()` API remains available.

Automatic snapshots use schema 13; resume retains the policy, implementation,
method seed and resource identity. Old serial snapshots stay serial. Switching
an automatic snapshot to `--ecm-workers serial` is rejected; use a fresh run
for a new execution policy. See [adoption limits](../POLICY.md#recorded-automatic-ecm-adoption-exception--5-november-2026).
