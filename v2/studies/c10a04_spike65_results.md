# 65-digit feasibility result — 12 October 2026

Neither arm completed the one fresh 65-digit product within its matched
1800-second wall/CPU grant. Retain current routing above 60 digits. This
spike supports neither a complete-factor speedup nor default interpolation
through 65. A longer grant might finish, but no extension or retry is implied.
The 70-digit feasibility step remains open on the roadmap.

| Observation | Compact | Control |
| --- | ---: | ---: |
| Complete | no | no |
| Stop reason | wall limit | wall limit |
| Inclusive/capped wall seconds | 1800.014 | 1800.001 |
| CPU seconds | 1749.139 | 1709.104 |
| Charged factoring work | 432808207231 | 100377204722 |
| Base entries | 8500 | 1612 |
| Multiplier | 19 | 1 |
| Exclusive base endpoint | 187952 | 30000 |
| Recorded owned workspace, bytes | 749009760 | 55847160 |
| Additional helper reserve, bytes | 16008192 | 0 |
| Retained relations | 6946 | 18 |

Compact used nine A factors, pool 64, half-width 65536 and at most 64 Gray
variants per family. Control reused the accepted 50/51-digit bundle at 65,
with six A factors and pool 32; it was not the automatic 65-digit reference.
Compact's last filter retained 4339 rows and 5474 columns; no solver result
was recorded. Control's last filter retained zero rows/columns and recorded
one dependency. More raw rows do not establish a factor or a speedup.

Both unresolved outcomes reconstruct the complete input exactly. Receipt
hashes, configuration, certainty and finite work/memory limits pass. Owned
workspace plus helper reserve stayed below 1 GiB. Recorded process RSS
high-water observations were 374521856 and 397885440 bytes, below 3 GiB;
they are cumulative observations from one sequential process, not independent
arm memory measurements. Wall-cap polling overshoot was about 14 ms and 1 ms.

The plan/runner were committed as `4edc39b` before once-only generation;
the fresh certified input was committed as `7ef1237` before measurement.
Prime/product exclusions and exact parity with the production counted policy
were independently checked. One conditioned balanced product and seed 7,
compact first then control, provide feasibility evidence only.

The first launch stopped during warmup when another checkout ran tests;
neither timed call had started. A new charged warmup was authorized using the
same committed input after those tests ended. Continuation
`d9cf145` retained the failed row, full 5-second/10^13-work reservation and
original deadline. Both new warmups passed. Each timed call retains its full
1805-second reservation. Final ledger totals are 4350 wall/CPU seconds and
70534349373872 work, including preparation, generation and the failed attempt.
There are no pending reservations, new failures or automatic extensions.

Raw receipts, failed ledger snapshot, console and final seal remain local
under `v2/benchmarks/results/c10a04/spike65/` and `v2/local/`. The seal records
all output hashes and detailed stages/work/memory. Its SHA-256 is
`0a76b2c8491b0f9a5f774e95032232919668dcde327cfd0c501acba3868c4c78`.
The compact receipt is
`ca9136f72f115e90b3549c8afdd42f46867cda08a3ad671e44cdd7f9e11050f8`;
control is `f1aa5d1695223f4acf7c04ee1584dcc91fe6a2d02838d6be1cee0dcec5fd7710`.
Neither this result nor progress counters justify extrapolated 61–65 defaults.
The separate 50–60 promotion remains as documented, with E1
deferred and C10a-07 open.
