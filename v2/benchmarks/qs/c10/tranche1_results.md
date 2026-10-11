# C10 tranche 1: bounded collection-policy study

This is the first bounded C10 tranche. The continuing C10 collection-strategy
and size-band work remains open. The [source mechanism study](research.md),
[source expansion](source_supplement.md), [frozen protocol](protocol.md) and
[C10a follow-ups](c10a_followups.md) give the arithmetic and source context.
Raw rows, profiles, checkpoint captures, expense records and cold starts are
local in ignored `v2/benchmarks/results/c10/experiment/`.

## Frozen controls and implementation

The runtime derives from C9's `129cb09e4ca684f3eb08cbbdbea3fb43848e5d95`
and unchanged runtime `fbfd662`; C9's final accepted 50-digit decision is
`868c26622031cddec5da61ec18e54c8b3809700b`. The exact configurations and
source hashes were committed in `c10_training.json` at `06a663b`. The training
selection, with evidence hashes, was committed at `b8419d7`; the six independent
certified 40/50-digit confirmation inputs at `3871843`. The 50 later C3/C9 n
values in the versioned exclusion file were excluded from fresh generation.

The isolated `quality` arm reorders fixed 16-family batches of the same flyer
population by exact `abs(A-target)`; it does not filter A products or enforce
inter-family diversity. The isolated `plan` arm enables the existing 2 MiB
prime-power plan for the same `powers` score policy. Other C9 residual,
storage, base and interval settings are identical to control. The opt-in
quality stream has exact integer assignment and checkpoint tests; the plan
uses the previously accepted bounded fallback. Neither arm changes the global
default in this tranche.

PyPy 3.11.15 with gmpy2 2.3.1 ran validated warmups of at least three seconds
of both wall and CPU time for every active arm in each process. Cold process
starts and instrumented 30-second profiles were captured separately. Every
timed row validates complete factorization or its unresolved cofactor,
proper divisors, reconstruction and certainty labels. Caps are 10/60/600
seconds per training call and 10/60/1,200 for confirmation at 30/40/50 digits;
work, owned memory, process RSS and aggregate charges follow the protocol.

## Training on previously revealed inputs

Training uses three known inputs per band and seeds 7, 29 and 47. The table
reports complete-call wall sums, including setup, collection, filtering,
solving, extraction and classification. All reported calls completed.

| Digits | Matched groups | Control | Plan | Quality | Frozen selection |
| --- | ---: | ---: | ---: | ---: | --- |
| 30 | 18 | 7.404 s | 7.669 s (+3.58%) | 7.746 s (+4.63%) | Control |
| 40 | 18 | 164.528 s | 127.968 s (−22.22%) | 166.610 s (+1.27%) | Plan |
| 50 | 11 | 1138.928 s | 923.405 s (−18.92%) | 1138.749 s (−0.02%) | Plan |

At 30 digits, within-input paired directions and repeated-cell drift remain
mixed after the full 18-group extension; no challenger advanced. At 40 digits,
the plan is faster on each of the three known inputs, but the quality arm's
mixed directions and drift forced the full 18-group extension. At 50 digits,
the quality arm's mixed directions forced extension. Its 11/18 groups are
complete and valid, with a 4/4/3 input allocation; the next full three-arm
reservation of 1,806 seconds would exceed the 5,400-second training phase
limit. The 50-digit training interval is therefore descriptive and
resource-censored. The plan advanced by the frozen completion-first, capped-
cost ranking, not by a training speed-acceptance claim.

The final 50-digit training process exited with `FileExistsError` while saving
its second `train-yield-50.json` marker. The prior lease had already written
that ignored marker. This happened after all 33 rows and their expense records
were saved; each of the 11 matched groups has all three arms. No call was
rerun or discarded to erase this record. The declared phase reservation
prevents another group irrespective of that marker-writing error. The runner
source and frozen identities were not tuned after observing data.

The 30-digit plan's +3.58% is a numerical training difference, not an
established slowdown: its paired interval crosses zero. Each of the 756
polynomials in those matched calls covered five blocks; the corresponding
40- and 50-digit polynomials covered 33 blocks each. The plan built 91,582
per-prime power plans and replayed them 457,910 times at 30 digits. Control
and plan had the same 3,780 blocks, candidate counts, relations and solver
dependencies. Plan collection/preparation totaled 6.075/0.495 seconds,
versus 5.913/0.437 for control. Building a per-polynomial plan has much less
reuse opportunity in the short 30-digit interval; these accounts support an
amortization explanation, but cannot isolate a causal time penalty from run
noise. At 40 and 50 digits the corresponding plan replay/build ratios were
33, and the fresh comparisons establish the scoped complete-cost benefit.

The plan preserves the control's candidate, relation and final dependency
counts on every matched training set. At 40 digits its measured collection
time totals 113.943 versus 150.421 seconds; CPU sums 122.020 versus 156.524
seconds and charged work 194.4 versus 208.2 billion units. Peak owned
workspace is 122.66 versus 120.56 MB. At 50 digits, collection is 600.078
versus 819.157 seconds, CPU 877.220 versus 1079.951 seconds, work 879.8
versus 960.7 billion units, and peak owned workspace 343.45 versus 341.36
MB. These are observed stage accounts, not a theoretical speedup ceiling.
Process RSS is a high-water mark shared across arms in each PyPy process and
is not a per-arm memory comparison.

Quality ordering produces 2318 versus 2272 relations and 134 versus 98 final
solver dependencies across the 18 30-digit calls, but complete cost is worse.
It also sent 6,200,864 versus 5,789,002 positions through candidate recovery;
5,969,784 versus 5,560,464 were rejected by exact refinement. Thus the
extra candidates almost entirely failed the exact value check. Quality
collection totaled 6.168 versus 5.913 seconds. Lazy batch sorting occurs
inside collection, so its cost cannot be separated from the changed candidate
work by these stage timers. Attributing the excess to sorting alone would
ignore the measured 411,862 additional candidate checks. The ordering changes
which early polynomials are visited before a successful stop; it does not
prove that low product error causes useful yield.
At 40 digits it produces 5940 versus 5924 relations and 86 versus 90 final
dependencies, also with worse complete cost. At 50 digits its 11 groups have
18060 versus 18048 relations and 48 versus 46 final dependencies, with near-
equal complete cost. These totals are repeated measurements on three known
inputs, not independent yields on 18 or 11 different numbers.

## Mechanism and replay diagnostics

Across nine input fixtures, exact product error for the first four and eight
families under `quality` is 1.9% and 6.0% of the flyer-order error. At every
complete 16-family batch boundary the product population and total error are
identical. Across 128-family prefixes, all-pair factor overlap is identical;
adjacent pairs sharing at least one A prime number 606 for quality versus 552
for flyer. The tested ordering changes early prefixes and adjacency, but does
not guarantee family diversity. These arithmetic diagnostics do not replace
useful-dependency or complete-cost evidence.

Separate 50-digit checkpoint probes encoded 6.7–7.2 KB payloads and charged
about 3.0 million work units for each restore. The plan probe included cached
plan reconstruction/refusals; quality rebuilt assignment state under the same
bounded allowance. All three resumed with the same identity and cumulative
budget. Three 30-second instrumented profiles stopped incomplete at their
declared caps. Their nested cumulative times and separate cold starts are
attribution only and are excluded from the warmed performance comparison.

## Independent 40-digit confirmation

The plan and control each completed 9/9 validated factor calls on three new
certified balanced inputs, each crossed with seeds 7, 29 and 47. No extension
was required: the frozen drift and paired-direction checks returned no reasons.
Control totaled 78.406 seconds and plan 61.290 seconds, a **21.83% reduction**
in complete-call wall cost. A 4,000-draw hierarchical input-then-sample
bootstrap gives a 95% cost-change interval of −24.88% to −12.47%. Each input
favors the plan: −11.14%, −21.70% and −24.93%, respectively. Completion and
the entire planned sample allocation satisfy the scoped acceptance gate.

CPU totals are 74.776 versus 58.494 seconds. Charged work is 104.0 versus
97.6 billion units. Collection stages total 71.283 versus 53.879 seconds;
the plan preserves 12,142,217 candidates, 3,014 relations, 2,125 DLP cycles,
2,954 final filter rows and 64 final solver dependencies across the nine
calls. Peak owned workspace is 123.05 versus 125.15 MB. Process RSS is a
shared high-water mark, not a per-arm delta. The result accepts the isolated
2 MiB plan for this explicit 40-digit growing-DLP bundle and certified-input
population. It does not promote a global default or a combined policy.

## Independent 50-digit confirmation

The plan and control each completed 9/9 validated factor calls on three new
certified balanced inputs and the same three seeds. The frozen stability checks
returned no reasons, so no extension was needed. Control totaled 1422.122
seconds and plan 1169.093 seconds: a **17.79% reduction** in complete-call wall
cost. The hierarchical 95% interval is a 16.11–19.48% reduction. Each input
favors the plan: 17.00%, 18.08% and 18.22% reductions. The balanced complete
sample plan, completion parity, stability and interval pass the scoped gate.
This fresh acceptance is independent of the resource-censored training
extension and does not promote its partial training interval.

CPU totals are 1349.994 versus 1110.575 seconds. Charged work is 999.9 versus
893.3 billion units. Collection totals 1003.842 versus 751.799 seconds; other
stage times remain comparable. The plan preserves 76,767,640 candidates,
14,568 relations, 8,940 DLP cycles, 13,940 final filter rows and 31 final
solver dependencies. Its 923,465 plan refusals used the exact fallback, with
no loss of complete outputs. Peak owned workspace is 343.34 versus 345.44 MB,
about 2.10 MB more for the plan. Process RSS high-water marks are shared by
arms within each process and cannot be read as an isolated memory delta.

The isolated 2 MiB plan is accepted for the measured 40- and 50-digit growing-
DLP configurations on this certified balanced-input population. The 30-digit
control remains selected; the fixed 16-family quality ordering does not
advance. There is no automatic default change, no combined C9/C10/E1
portfolio claim and no claim about larger sizes, different backends or input
populations. C10a records conditional representation, diversity, threshold,
coupled-size and fixed-cadence follow-ups without taking over other owners.

## Resource, validation and interpretation limits

The complete tranche charged 7,537.840 seconds against the frozen 14,400-second
aggregate wall/CPU allowance, including the conservatively prepaid 900-second
QA block, 3,771.316 seconds of training, 91.822 of diagnostics, generation and
2,774.575 of confirmation. The training phase stayed below 5,400 seconds;
its 50-digit extension was stopped by the predeclared whole-group reservation,
not by a factorization failure. The 50-digit confirmation used two process
leases, each with new validated warmup. Cold starts and instrumented profiles
were kept separate. Captured local evidence uses about 1.49 MB, below the
256 MiB cap; diagnostic checkpoints are far below the 16 MiB cap.

The source and runner identities were frozen before timing. No arm, target,
seed, observation cap or acceptance threshold was tuned after seeing results.
Training fixtures were previously revealed; only the certified disjoint
confirmation corpus supports adoption. The generator conditions on a large
factor of p−1 and is not RSA-distributed. Three independent numbers per band
remain three independent numbers despite nine seeded calls; the hierarchical
interval and per-input effects make that limit explicit. C10 remains a
continuing workstream after this bounded tranche.

Committed-files-only checkouts at `a13dd3b` and the final documentation commit
`ef3e6c1` each ran **615 PyPy 3.11 tests without failures or skips**, passed
`make -C v2 lint`, and imported all 129 benchmark modules plus the A7 current
loader. The frozen control/selection hashes and certified fresh corpus were
checked again after measurement.
The detailed [validation record](validation.md) distinguishes the earlier
implementation check from this final committed-source check.

## Subsequent integration boundary

The historical retain-defaults verdict above describes this frozen first
tranche. The user subsequently directed a [scoped automatic selection](default_promotion.md)
of its accepted 2 MiB plan in the exact native 40/50-digit growing-DLP bundles.
That routing is a later integration decision. It changes no first-tranche
controls, inputs, samples, source pins, expense records or timing claims.
