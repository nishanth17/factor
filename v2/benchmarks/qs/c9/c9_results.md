# C9.1 useful large-prime yield: evidence and acceptance

Status: C9.1 stopped during confirmation at the user's direction; no promotion.
C9.2 reopens implementation around growing storage and upstream strategy gaps.
This report
separates source-supported hypotheses, instrumented attribution and fresh
complete-call evidence. The [research](c9_research.md) and immutable initial
[protocol](c9_protocol.md) explain the investment and finite limits.

## Historical control and ownership

The parent is integrated `b3b3cfb`; C9 works on `codex/c9-large-prime-yield` in
an independent checkout. C3's portfolio/handoff work has separate ownership.
No active C10 chat was found. The joint arms consume only C1's already studied
half-base configuration; broader collector calibration remains C10. No sieve,
polynomial, matrix, provenance, splitter or checkpoint implementation changes
are required by the selected policy hypothesis. Existing bounds encode the
policy and their values already participate in configuration identity.

C1's **31.73%** fresh 40-digit complete-cost improvement remains an accepted
opt-in result with its original inputs and 29.21–33.99% interval. Its separate
fresh 50-digit SLP/DLP comparisons completed **0/18 per arm at 120 seconds**.
This was full-base 30,000 DLP versus the 900M SLP policy on two fresh
balanced inputs (nine calls each); DLP median evictions were **21,504 and
28,672**, respectively. C1's half-base resume witness
completed in **141.063 cumulative seconds**, using 512 MiB, an 8-MiB checkpoint
and a 300-second allowance. It changed both configuration and allowances.
None of these is replaced by the new two-input training result below.

## Individual-arm attribution

Source/protocol were frozen in `22bd076` before valid training; exact
configurations and fixture certificates are versioned in
`inputs/controls/c9_training.json`. Two previously inspected C1 training
inputs use seeds 7/29, 240-second work/CPU/wall caps, 256-MiB owned memory,
8,192 rows and 32,768 atoms. Current full-base means numeric bound 30,000,
interval half-width 65,536 and six A factors. DLP uses endpoint 3M/product
115.2B; alignment changes only its SLP residual ceiling from 900M to 3M.
The low-bound SLP arm also narrows its scoring domain and is not residual-only
attribution. Both retain exact integer recovery and proven endpoints.

| Arm | Input 0 seconds | Input 1 seconds | Complete | Capped total seconds |
| --- | ---: | ---: | ---: | ---: |
| SLP, 900M | 110.839 | 126.447 | 2/2 | 237.285 |
| SLP, 3M | 115.938 | 86.878 | 2/2 | 202.816 |
| DLP, current | 207.483 | 107.456 | 2/2 | 314.939 |
| DLP, aligned SLP endpoint | 127.586 | 69.902 | 2/2 | 197.489 |
| DLP, 16,384 unowned edges | 31.881 (memory refusal) | 17.268 (memory refusal) | 0/2 | 480.000 |

These are instrumented, single-seed training calls, not accepted speedups.
The aligned DLP lead over calibrated SLP is only 2.63% in this small training
sum and reverses on input 0. The predeclared selection rule chooses 3M SLP;
it does not cherry-pick a different control for each input. All failures cost
the full 240-second cap when comparing outcomes, while actual resource use
is also retained.

| DLP arm | Final rows (0 / 1) | Filtered rows | Evictions | Splits | Failed splits | Owned peak MiB |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Current | 1587 / 1610 | 1237 / 1329 | 78848 / 75776 | 27272 / 26703 | 0 / 0 | 193.74 / 194.85 |
| Aligned | 1612 / 1633 | 1442 / 1457 | 12288 / 12288 | 15488 / 15578 | 0 / 0 | 219.30 / 222.41 |
| 16K retention | 293 / 276 | 0 / 0 | 0 / 0 | 5019 / 5063 | 0 / 0 | 255.39 / 254.93 |

The 16K refusal applies to the frozen 8,192-row/32,768-atom reservation
bundle. It does not prove that 16K retention is always uneconomic: a smaller
matrix-row reservation or different finite storage design could make it fit.
Upstream YAFU/msieve preserve history through growing memory and/or savefiles;
FLINT also grows relation arrays. Their initial capacities are not eviction
limits. Budget-aware growth, packed storage and bounded disk retention remain
unmeasured algorithm candidates; admission alignment does not settle them.
That is an untested coordinated resource bundle, not a reason to rescue this
arm with more memory. The two permitted follow-ups instead use the historical
half-base configuration, supported by C1's prior useful-prefix witness.

Alignment retained more useful combined rows and pinned cycle paths, so its
owned memory **increased** despite fewer evictions and lower full-call cost.
Current DLP still held 4,000/4,584 SLP atoms above its DLP endpoint at completion;
the capacity arm held 8,607/8,581 when refused. This supports a collision-history
and admission interaction, rather than a claim that eviction count alone
causes failure. No split failed, no duplicate was rejected and no DLP cycle was
lost to the length cap in these runs. Certification was not weakened.

Delayed graph structure matters. On input 0, current DLP had 472 rows at 60.2s
and 969 at 120.9s, with no surviving filtered matrix at either prefix. It
completed at 207.5s with 1,237 filtered rows. Alignment had 754 rows but no
filtered survivors at 60.3s; at 120.5s it had 1,535 rows and 1,328 survivors,
then completed at 127.6s. On input 1, aligned DLP already had 1,231 filtered
rows at 60.8s and completed at 69.9s. Both ended with a verified useful
dependency and proper factor. Empty short prefixes were not used as stop rules.

Graph diagnostics count original retained atoms, including closing edges,
not just the acyclic spanning forest. On the two completed DLP calls, current
retained cycle ranks were 1,587/1,610 and 2-core edge counts 2,254/2,258;
alignment produced 1,612/1,633 and 3,185/3,235 respectively. These counts do
not replace matrix kernel or proper-factor evidence. SLP retained incidence
can have extra combinations which its disjoint-pair matcher did not publish:
31/38 extra LP basis vectors for the 900M arm and 104/123 for the 3M arm.
This measures a matching-policy opportunity, not additional useful factor-base
dependencies or a measured speedup. Reusable SLP representatives remain a
separate bounded C9 candidate; this study does not claim globally optimal SLP.

The initial window used **1002.353 wall / 965.113 CPU seconds** including
prefix diagnostics. One earlier attempted launch was aborted by the process
guard when a C3 metadata-only benchmark command appeared; no sample was
accepted and a conservative **240 seconds each** is charged for that attempt.
Raw manifests and the abort receipt remain ignored under `results/c9/`.

Instrumentation limitation: the initial prefix inspector builds additional
bounded adjacency/core scratch outside the production owned-reservation
counter. The 256-MiB number is the job's reservation, not a measured total for
the instrumented harness. RSS remained below the separate 1-GiB abort bound.
This limitation prevents treating the prefix calls as accepted timings;
fresh timing calls omit the extra graph inspector. Census diagnostics have
separate explicit simultaneous-storage reservations.

## Sampling coverage and loss attribution

The retained C1 census visits every candidate in the conservative 128B²
report domain until its 65,536-record cap. Independent Bernoulli sampling
uses geometric skips with probability 1/256 per position, independent of
success and block length, and fully trial-divides sampled positions to
cross-check bucket recovery. Its 32,768-position validation prefixes end at
8,424,776 and 8,446,400 positions. The full censuses continue to 138,602,985
and 161,894,304 positions. The sample is representative of those **early
prefixes**, not of all later collection or arbitrary input bands.

Of the sampled positions, 32,735/32,751 fail the original SLP selection and
33/17 pass. Every sampled exact residual inside the widened product bound
passes its conservative widened score/refinement check. Across the larger
retained census, recoverable DLP losses relative to SLP divide into block
score **9004/7707**, refinement **2544/2477**, and residual policy
**8528/9989**. Endpoint rejections were **136/146** of **20212/20319** split
attempts. These retained observations establish missing SLP-domain opportunity;
they are not fresh estimates or proof that widening pays for its cost.
C9's isolated alignment arms hold DLP score and product domains fixed.

## Fixed-prefix retention counterfactual

A separate 7.079-second replay independently validates every retained atom
and its proven endpoints, replays the current transactional forest with an
8,192-unowned-edge cap, then reconstructs factor-base rows from original atoms
and verifies each extracted square/GCD. The offline analysis reserves raw
records, simultaneous forest/replay maps and matrix/provenance within a
separate 512-MiB diagnostic allowance. It is not a complete-call timing arm:
acquisition, score recovery and the stored input history were already paid by C1.

| Same 65,536-record prefix | Input 0 | Input 1 |
| --- | ---: | ---: |
| SLP records above DLP endpoint | 35017 | 34765 |
| Distinct such endpoints | 34923 | 34666 |
| Independent repeat combinations lost by alignment | 94 | 99 |
| Optimistic no-eviction LP rank, current / aligned | 6177 / 6083 | 6472 / 6373 |
| Finite retained LP cycles, current / aligned | 1136 / 2427 | 1193 / 2494 |
| Finite evictions, current / aligned | 56320 / 18432 | 56320 / 18432 |
| Filtered rows × columns, current | 401 × 767 | 578 × 938 |
| Filtered rows × columns, aligned | 2376 × 1577 | 2449 × 1605 |
| Algebraic nullity, current / aligned | 0 / 800 | 0 / 845 |
| Proper factors, current / aligned | 0 / 1 | 0 / 1 |

This directly measures the deliberate loss: alignment removes 94/99 possible
high-endpoint repeat combinations from an unlimited-history graph. Under the
finite FIFO policy it instead preserves enough more connected history to add
1,291/1,301 cycles and useful factor-base dependencies. The current finite
prefix has nonempty filtered matrices but full row rank, so counting survivors
alone would still overstate progress. One disconnected LP cycle occurs on
input 0 under both policies; the all-component implementation preserves it.
Aligned input 1 requires four extracted dependencies (three trivial, then a
proper factor). The other 799/841 algebraic dependencies are verified kernels
but are **not** reported as extracted square congruences or useful factors.
Every returned/unresolved output reconstructs and matches the independent
certified training fixture.

The source record SHA-256 values are
`5f48825c298e75a937a8d81ab8a22491cef6bab8283386587d693aab05e3b74f`
and
`edcecfec1fed9f76afe609d1b6927ae3fde0cd6c8834beda844c3b3b0c5e0a3b`.
The 65,536-record cutoffs, input order and original sampling scope remain C1's.
Unlimited-history rank is only an optimistic diagnostic, not a feasible
memory-matched challenger.

## Frozen joint collector study

The second freeze is committed in `9006eae`; exact settings and the completed
individual-arm gate are in `inputs/controls/c9_combinations.json`. It compares
C1's historical half-base DLP against that same bundle with endpoint alignment.
Both use base 15,000, endpoint 1.5M/product 28.8B, interval half-width 65,536,
six A factors and unchanged finite stores. The SLP ceilings are 225M and
1.5M respectively. These remain instrumented training calls, not fresh timings.

| Joint arm | Input 0 seconds | Input 1 seconds | Complete | Total seconds |
| --- | ---: | ---: | ---: | ---: |
| Historical half-base DLP | 77.510 | 86.588 | 2/2 | 164.098 |
| Half-base DLP plus endpoint alignment | 45.603 | 53.616 | 2/2 | 99.219 |

The prescribed ranking selects **half-base plus alignment** as the DLP
challenger and full-base **3M SLP** as its calibrated SLP comparison. The
candidate's training cost falls 39.54% against the same half-base collector,
51.08% against selected SLP and 68.50% against current full-base DLP. These
percentages select a candidate; they do not establish accepted benefits.
SLP calibration compares the two frozen endpoint policies at the existing
full base, not every possible SLP base/interval combination. Broader joint
SLP/DLP optimization remains C10.

The half-base comparison isolates alignment at that collector: both controls
evict 20,480 edges; both aligned calls finish without eviction. Filtered rows
rise from 806/874 to 840/901, while owned peaks rise from 179.81/185.49 to
190.08/193.55 MiB. Split attempts fall from 11,246/11,374 to 7,067/7,077;
all splits succeed or reject an out-of-bound endpoint, with no failed split.
The final solver produces 2/2 dependencies for the control and 2/3 for the
aligned bundle; every complete output has a proper factor and independently
certified reconstruction. The retained incidence cycle ranks are 858/911 and
861/914, with 2-core edges 1,747/1,850 and 2,352/2,475. The aligned calls retain
9,939/9,929 original atoms, distinct from their 861/914 combined rows.

Cross-base effects are not isolated C9 credit. On input 0, half-base alignment
scans **137.4M** positions versus full-base alignment's **106.6M**, yet uses
51.35B work units versus 117.97B and has a much smaller parity universe.
More reports or positions alone would rank those outcomes incorrectly.
The collector and residual-policy interaction is the useful finding.

## Separate profiles, resource checks and charged restore

The joint window took **438.160 wall / 436.818 CPU seconds**, including all
four calls, the repeated census, two profiles and repeated restore. A second
census enforces the 512-MiB simultaneous-storage reservation before graph
publication, checks RSS and repeats independent output validation. It takes
7.913 wall / 7.870 CPU seconds; its process maximum RSS is 707.45 MiB, below
the separate 1-GiB abort. Both census attempts are charged. PyPy heap/JIT RSS
is distinct from live owned reservations, and this diagnostic is not a
256-MiB timing competitor.

Both 60-second profiles are intentionally capped and uncompleted. Current
full-base DLP spends 38.44 cumulative profile seconds in `_sieve`, including
27.91 in `prime_power_roots`; exact residual handling takes 1.90, classification
0.63 and splitting 0.59. Half-base alignment has corresponding sieve/root
values 36.52/25.00 and residual/split values 0.75/0.34. These nested times
must not be added. The profiler also magnifies the process guard (about 5.1s
per profile), so these are qualitative attribution, not production speed
fractions. Combined with zero failed training splits, they do not justify
new cofactoring or a residual-certainty shortcut. Scoring/root optimization
belongs to C10.

The selected training input is paused twice and restored using the **same**
300-second wall/CPU and 10¹³-work budget. Old stores are released before
restore, which independently verifies atoms, graph paths and provenance.
The 10/20-second checkpoints contain 429,090/845,045 bytes; rebuilds charge
10,211,836/19,059,716 work units. Post-restore cumulative wall is 10.306/20.338s.
The run completes with the expected proper divisor at **46.173 wall / 45.066
CPU seconds**, charging 51,380,746,662 work units. Its 512-MiB owned and
8-MiB checkpoint allowances make this a correctness/resource diagnostic only.
The process maximum RSS is 714.38 MiB and includes preceding census/profiles;
it is not isolated per-job RSS or a speed comparison with the 256-MiB arms.

## Confirmation identity and finite accounting

Selection was committed as `7085ba7` **before** generating the seven fresh
certified inputs committed in `ff90aec`. The selected control records every
runtime/benchmark Python hash and all collector settings. Three independent
50-digit and three 40-digit inputs each use seeds 7/29/47; one independent
30-digit input supplies repeated regression samples. Nine matched calls do
not mean nine independent integers. Pocklington generation conditions primes
on a large factor of p−1; this population is not RSA-distributed. Expected
factors and certificates enter independent validation only.

`inputs/controls/c9_accounting.json` records prior receipts and conservative
charges: **2,310 wall / 2,270 CPU seconds** before confirmation, including the
entire 600-second heavy-QA allowance prepaid for working and committed-only
verification. Unused QA reservation is not refunded. The persistent ignored
`expenses.json` adds generation, separate cold calls, every warmup and every
confirmation/extension call; an interrupted call keeps its full reservation.
The total remains bounded by 14,400 seconds on each resource. Lock waiting
and ordinary source/documentation work are excluded.

## Reproduction

Use the committed C9 branch with PyPy 7.3.23 implementing Python 3.11.15.
The measured host is Apple M4 (10 logical CPUs), 24 GiB, macOS 26.6.2.
The runners acquire the machine-wide lock and publish ownership; coordinate
with other factoring studies even for metadata-only benchmark commands.
Do not regenerate or overwrite the frozen controls/corpus to repeat a run.

```sh
pypy3 -B -m v2.benchmarks.qs.c9.c9_study train \
  --output v2/benchmarks/results/c9/repeat/training
pypy3 -B -m v2.benchmarks.qs.c9.c9_combinations run \
  --output v2/benchmarks/results/c9/repeat/combinations
pypy3 -B -m v2.benchmarks.qs.c9.c9_confirmation confirm --digits 30 \
  --output v2/benchmarks/results/c9/repeat/confirmation
pypy3 -B -m v2.benchmarks.qs.c9.c9_confirmation confirm --digits 40 \
  --output v2/benchmarks/results/c9/repeat/confirmation
pypy3 -B -m v2.benchmarks.qs.c9.c9_confirmation confirm --digits 50 \
  --output v2/benchmarks/results/c9/repeat/confirmation
pypy3 -B -m v2.benchmarks.qs.c9.c9_confirmation summarize \
  --output v2/benchmarks/results/c9/repeat/confirmation
```

These are reproducibility runs on now-public inputs, not new held-out
confirmation. The original study's expense ceiling and selection remain
immutable. Start another study only with its own declared resources and
identity. The training runner accepts only its initial source hash; the
combination runner checks its additional freeze; confirmation checks every
runtime/benchmark Python hash. Documentation changes do not alter these
identities. `confirm` resumes saved calls, refuses a foreign manifest and
does not append samples to a finished band.

Optional fixed-prefix replay requires the ignored C1 files
`50-0-records.json` and `50-1-records.json` with the hashes above. Their absence
does not prevent versioned training, confirmation, tests or imports. To
repeat that diagnostic, supply their directory to `c9_census --retained`
and an ignored output path. Profiles/resume use `c9_diagnostics profile` or
`resume --config v2/benchmarks/inputs/controls/c9_combinations.json
--arm half_aligned --output <ignored-path>/resume.json`; choose an empty
output directory because checkpoint publication refuses replacement.
Generated JSON, profiles, transcripts and scratch remain ignored.

## Measurement limits

Accepted calls omit the extra adjacency/core inspector. They retain native
collector counters, final forest/unowned occupancy, filtering/solver/extraction
statistics, peak owned reservations and process maximum RSS. Detailed endpoint
degree/core/rank measurements come from instrumented training and the separate
fixed-prefix oracle. RSS is a **cumulative process maximum within a band**,
including warmup and earlier arms, not independently isolated arm allocation;
the cold diagnostics use separate processes.

Final retained atoms equal admitted atoms minus evictions. Exact peak atom
and unowned occupancy are not sampled inside accepted calls. Their finite
bounds remain explicit: peak unowned occupancy equals the configured cap
when any eviction occurs; otherwise it lies between the final count and
that cap. Since published rows and pinned atoms are never removed, peak
retained atoms lie between the final count and that count plus
`max_partials - final_unowned`. Report these as bounds, not observed peaks.
This avoids adding a second, uncharged graph snapshot to a timed call.

## User-directed implementation revision

The user redirected the study on 10 October 2026: implement adaptive storage,
remove arbitrary normal-path retention/memory ceilings where the upstream
implementations have none, and audit other consequential strategy gaps.
The 30-digit band finished 18 matched samples, the 40-digit band nine. The
50-digit band stopped after four complete matched samples (12 calls); the
next DLP call was interrupted and its full 242-second reservation remains
charged in the new ledger. All completed outputs were validated. These
partially revealed inputs are now training/regression evidence, not fresh
confirmation for any new implementation. No C9.1 policy is promoted.

The original source/configuration freezes, receipts and raw outcomes remain
unchanged. The revised protocol will account for all prior expenditure and
freeze a new implementation/selection before generating different fresh
inputs. The source-backed growth mechanism, reusable matching and arbitrary
stopping/storage restrictions are reopened; a failed 16K reservation bundle
is not a reason to preserve a fixed-cap production strategy.

