# A1a integration: frozen two-arm confirmation

Prepared 29 October 2026. The user stops preparation before accepted timing,
then manually switches to Sol high for execution, analysis and reporting.
No accepted measurement or calibration has run during preparation. E1 remains
deferred. This protocol tests one completed A1a combination, with no search,
threshold tuning or additive-gain claim.

## Sources, scope and arithmetic bridge

Control is committed pre-integration master
`92543797f4ed3fae347ecf5a3e938361d7ec7047`, replayed from the hash-checked owned
`inputs/baselines/a1a_integration_mainline.json`. Its private package preserves
its own configuration, relation, budget and arithmetic classes. The selected
arm uses the frozen integration source with
`marking="density+presieve+scan_fused"`. All other actual fields agree.

Use exactly native 50 and 55 decimal digits, 1 GiB owned memory, half-width
65536, 4096-position list/powers/bucket blocks, 2 MiB power plans and collection
pacing. The 50-digit bundle uses bound 30000, multiplier 1, six A factors and
pool 32. The 55-digit bundle uses CountedSIQSConfig: exactly 2500 eligible base
entries under endpoint ceiling 250000, selected multiplier, seven A factors,
pool 64, 64 Gray variants. Setup and count resolution share each job's budget.

The completed study at `d61df6c` supports its own fixed 50/55 workloads (5.945%
and 3.214% lower complete cost). Its 55-digit resolver ran before SIQSJob;
mainline's counted setup and memory preflight run inside the job. Independent
count/root/configuration/replay checks bridge arithmetic; these new timings
must establish integrated economics. Detailed source reconciliation stays local
under `v2/audit/a1a_integration/`; no historical controls or sources are repinned.

Current automatic defaults remain unchanged during preparation. A passed gate
permits a subsequent reviewed selector only at exactly 50/55, native arithmetic,
the measured configuration and 1 GiB envelope. It must preserve explicit
configurations and smaller work/time/memory grants. Omitted-budget routing and
any future allowance guard must be reconciled in standalone and portfolio entry
points before adoption. Saved jobs replay saved marking, with charged disposable
cache rebuilds and no serialized trust. 40 and 51 digits, 52-54/56-60 policies,
other backends and all other routing retain their current behavior.

## Fresh population and fixed order

Commit this source/protocol/control before generating new inputs. Generate
nine independently certified products per band with seed 20261030015055, crossing
balanced and 40:60 factor-bit shapes with odd residues 1/3/5/7, plus one extra
balanced residue 1. Exclude all collected historical factors/products, including
the completed A1a worktree, via the committed exclusion set. Generation has
115 seconds, one million draws, 10000 attempts per stratum and 8192 certificate
entries. Proofs are verified independently. Commit the complete fresh corpus
before timing, with no outcome-based replacement. These conditioned p-1 prime
populations are not arbitrary semiprimes.

Each base stage runs both arms on every fresh input, seed 7, then repeats every
input in the opposite arm order. Alternate the first order by input index.
There are 18 matched pairs/36 calls per band, nine independent statistical units.
Run one band per process. Print progress only at complete matched-pair boundaries.

## Warmed execution and stopping rules

Each measurement process first warms both separately loaded runtimes using
committed revealed inputs, never the fresh panel. Both arms run 1056-block
prefixes with 30-second wall/CPU caps, alternating order. Require at least three
validated wall AND CPU seconds per arm and three consecutive semantically equal
prefixes whose wall and CPU range/median are each <=10%. At most 20 calls/arm
per stage. Failure stops that band, retains every charge and keeps its default.
Calibration and accepted measurements occur in the same process.

Accepted calls include configuration/job setup, count resolution, complete
factoring, proper-divisor and exact-output validation, certainty classification
and graph teardown with explicit gc.collect(). Preserve unresolved cofactors.
Work cap is 10^13 per call, wall/CPU caps 480 seconds at 50 and 1159 at 55. Owned
memory is 1 GiB and RSS is 3 GiB, checked cooperatively. Source loading/imports and
orchestration have separate receipts/prepayments; CLI cold startup and profiles
are outside the warmed endpoint and supply no speed claim.

After the whole base schedule, inspect only stability and correctness, never a
cost-effect interval. Repeated identical input/seed/arm cells require wall and
CPU range/median <=15% and unchanged deterministic paths. If unstable, exactly
one full opposite-order schedule extension is justified: rerun every input and
both arms twice, under the original caps, in a newly calibrated process.
All original adverse/noisy/censored rows remain. No selective rerun, parameter
change or second extension. Persistent all-data drift prevents acceptance.

One final inferential look per band follows its required schedule. Average
opposite orders within each independent input. Resample the nine input means
50000 times (seed 2026103055) for a 99.5% percentile interval of the total-cost
change. An exact one-sided input sign-flip test must have p <=0.005, and the
interval's upper endpoint must be <0. No 10% floor. Every call in both arms must
complete, with no per-input completion loss, valid reconstruction/certainty,
stable warmup and no final drift. Incomplete outcomes remain censored and
cannot establish complete-factor speed or coverage benefit. A failed gate
retains that band's current default; scope is never broadened.

## Protected finite resources and durable receipts

Wall and CPU ceilings are identical, separately protected; unused accounts
cannot be transferred. Successful calls retain their time prepayments; work
settles actual use. Failure/interruption retains full reserved work/time plus
any observed overshoot. Existing C10a-04 Ledger supplies fsynced atomic receipts,
hashes and reservations; no historical ledger enters this new inference.

| Account | Wall and CPU seconds each | Work ceiling |
| --- | ---: | ---: |
| Generation | 120 | 0 factoring units; separate finite draw cap |
| Warmup 50, including possible extension warmup | 2800 | 8e14 |
| Warmup 55, including possible extension warmup | 2800 | 8e14 |
| Base 50 | 17460 | 3.6e14 |
| Base 55 | 41904 | 3.6e14 |
| Extension 50, conditional | 17460 | 3.6e14 |
| Extension 55, conditional | 41904 | 3.6e14 |
| Orchestration | 600 | 0 |
| Final analysis | 600 | 0 |
| Total | 125648 | 3.04e15 |

Evidence plus the materialized owned control is capped at 128 MiB. Each group
requires room for both full call reservations before it starts. Measured rows
carry fixture/seed/arm/order, actual configuration, base identity, work, elapsed
wall/CPU, stats, output and memory. Stage receipts record whole-stage time and
outside-call overhead, covered by a 60-second prepayment per stage. An overhead
overrun, mathematical disagreement, failed accepted call or interrupted stage
requires audit and stops execution. Never restart an unfinished stage: a cold
process cannot reuse an earlier process's warmup. Recoverable receipt publication
windows retain their original outcome; lost measured outcomes become censored,
with full charges, and stop for audit. Do not reconstruct or reset ledgers.

Inspect shared coordination state and acquire
`/private/tmp/factor-performance.lock` before generation, heavy QA or timing.
The runner checks competing benchmark/test processes and owns the lease for
its entire stage. Output ownership uses its own nonblocking lock. No messages
to other sessions, parallel heavy checks, agents, merge or push.

## Commands after the manual model switch

Use PyPy implementing Python 3.11. In the integration worktree:

```sh
pypy3 -B -m v2.benchmarks.qs.a1a_integration.study verify
pypy3 -Bu -m v2.benchmarks.qs.a1a_integration.study run --digits 50
pypy3 -Bu -m v2.benchmarks.qs.a1a_integration.study run --digits 55
```

After each base stage, inspect its `base-finished-N.json`. Only if its frozen rule
requires extension, run `extend --digits N`. Then run `analyze --digits N`.
The default output directory is
`v2/benchmarks/results/a1a_integration/frozen/`: `ledger.json`, per-call JSON,
environment, stage-finished markers and terminal `decision-N.json`.

Finish with concise public results and limitations, local detailed receipts,
tests/lint and committed-only imports/loaders. Preserve the separate clean
60-digit inconclusive result and 65-digit noncompletion. No automatic adoption,
roadmap acceptance closure or E1 claim follows merely from preparation.
