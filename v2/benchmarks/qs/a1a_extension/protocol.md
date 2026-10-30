# A1a bounded extension: frozen combined comparison

Prepared 30 October 2026 after explicit authorization to extend preparation.
Sol xhigh implements, verifies and freezes this comparison. Stop before its
accepted calibration or measurements; the user manually switches to Sol high
for execution, analysis and reporting. No automatic adoption is enabled.

## Sources and exact scope

Control is the exact selected A1a implementation at committed integration
`fee6a41c9af2702c34f964536e71ad3011d4cb22`, replayed from the hash-checked
owned `inputs/baselines/a1a_extension_mainline.json`. Its own private package
retains its configuration, budget, arithmetic and relation classes. Control
marking is `density+presieve+scan_fused`. The sole challenger uses
`density+presieve_powers+unroll4+scan_fused` in the frozen extension source.

Change marking alone: native exact 50/55 decimal digits, 1 GiB owned memory,
half-width 65536, block width 4096, list/powers/bucket collection, a 2 MiB
power plan and collection pacing. The 50-digit config retains bound 30000,
multiplier 1, six A factors, pool 32 and 32 actual Gray variants. The 55-digit
counted config retains exactly 2500 eligible entries under ceiling 250000,
selected multiplier, seven A factors, pool 64 and 64 Gray variants. Counted
setup is included in each job's budget and measured endpoint.

This tests the combined extension only. Individual selectors are implemented
and tested, but there are no component-effect looks, adaptive ablations or
post-outcome arm choices. A failed or inconclusive combined gate retains the
prior A1a control. No 40/51/52–54/56–60 policy, explicit caller config, lower
grant, other backend, checkpoint policy or E1/portfolio gate is promoted here.
Quota-aware setup, reusable switch metadata, lifted-power reuse and A1c
two-root interleaving remain roadmap work with separate evidence gates.

The original integration at `fee6a41` remains separately replayable in its
preserved checkout, including all original source/control/corpus hashes and
the existing generation charge. This extension does not repin, replace or
consume that comparison's ledger. Passing this incremental comparison alone
does not establish a complete extended-versus-pre-integration-mainline claim.

## Fresh inputs and sampling

Commit implementation, tests, runner, protocol and frozen control before new
generation. Use seed 20261030025055 for nine certified independent products
per band: four odd residues crossed with balanced and 40:60 factor-bit shapes,
plus one extra balanced residue-1 case. Exclude the original integration's
entire fresh panel and its constituent primes, in addition to its historical
exclusion inventory. Required proofs and the new corpus are committed before
timing. The generation cap is 115 seconds, one million draws, 10000 attempts
per stratum and 8192 certificate entries; its reservation is 120 seconds.
These conditioned p-1 prime populations do not represent all semiprimes.

Both arms run seed 7 on every input twice, in opposite orders. Alternate first
arm by input index and repeat; one process owns a whole band/stage. There are
18 pairs/36 calls per band and nine independent inference units. Never tune
on either fresh panel or replace an input based on an outcome.

## Warmup, validation and stability

In each measurement process, calibrate both privately loaded arms on committed
revealed warmup fixtures, not fresh inputs. Alternate order and use equal
1056-block prefixes, 30-second wall/CPU call caps, at most 20 calls per arm per
stage. Require at least three validated wall AND CPU seconds per arm, and
three consecutive semantically equal prefixes with wall and CPU range/median
each <=10%. Calibration failure stops the band and retains every charge.
Calibrate and measure in the same process; cold imports/startup and profiles
are outside this warmed endpoint.

Accepted calls include config/job/count setup, complete factoring, exact
reconstruction/proper-divisor validation, certainty classification and graph
teardown with explicit collection. A call has 10^13 work units, wall AND CPU
caps of 480 seconds at 50 digits or 1159 at 55, 1 GiB owned memory and a 3 GiB
cooperative RSS cap. Reservations add five seconds. Record unresolved factors
and all work on timeout, storage refusal, cancellation or interrupted calls.

After the complete base schedule, inspect correctness and repeated-cell
stability only. Require wall and CPU range/median <=15% for each identical
input/seed/arm cell, with unchanged deterministic paths. If unstable, allow
exactly one full opposite-order extension for every input and both arms in a
newly calibrated process. Retain the original adverse/censored/noisy calls;
all-data stability is required. No selective reruns, effect look before the
extension decision, retuning, or second extension is permitted.

## Acceptance and resource accounting

Use the same fixed two-arm analysis as the original integration: average
opposite orders within each independent input, 50000 deterministic input-level
bootstrap draws for a 99.5% interval, and a one-sided exact sign-flip test over
the nine paired input means. Accept a band only if every scheduled call
completes, deterministic paths agree, warmup and all-data stability pass,
the cost-ratio-change interval's upper endpoint is below zero and exact
p<=0.005. There is no arbitrary minimum percentage improvement. Report both
bands and every loss, failure or inconclusive result; no additive claim.

Separate output: `v2/benchmarks/results/a1a_extension/frozen/`. Durable ledger
identity is the new frozen control hash. Full reservations remain charged,
including unsuccessful preparation within the runner; lost measured receipts
become censored outcomes. Existing/unfinished stages cannot silently restart.
Admit a whole matched pair before starting it. Use the shared performance
lock and inspect coordination before taking ownership. Preserve prior owner
metadata on exit. Disk allowance is 128 MiB, including the owned runtime.

| Protected phase | Wall and CPU seconds each | Work units |
| --- | ---: | ---: |
| Generation | 120 | 0 |
| Warmup 50, both possible stages | 2800 | 800000000000000 |
| Warmup 55, both possible stages | 2800 | 800000000000000 |
| Base 50 | 17460 | 360000000000000 |
| Base 55 | 41904 | 360000000000000 |
| Stability extension 50 | 17460 | 360000000000000 |
| Stability extension 55 | 41904 | 360000000000000 |
| Orchestration | 600 | 0 |
| Final analyses | 600 | 0 |
| Total | 125648 | 3040000000000000 |

Orchestration prepays 60 seconds per stage and records/charges an overrun
before failing. Each final analysis prepays 300 seconds. These are separate
finite allowances for the explicitly authorized extension, with no transfer
from the original integration ledger. They are worst-case reservation limits,
not runtime predictions. There is one final inferential look per band.

## Commands after the manual model switch

From the committed extension checkout, use PyPy implementing Python 3.11:

```sh
pypy3 -B -m v2.benchmarks.qs.a1a_extension.study verify
pypy3 -B -u -m v2.benchmarks.qs.a1a_extension.study run --digits 50
pypy3 -B -u -m v2.benchmarks.qs.a1a_extension.study run --digits 55
```

Run `extend --digits BAND` only when that band's base-finished receipt declares
it required. Then run `analyze --digits BAND` once. Stop for audit on an
unfinished/failed stage; preserve charges and receipts. Do not rerun generation
or reset either ledger. Handoff records the final exact source/control/panel
identities, verification receipts and remaining allowances.
