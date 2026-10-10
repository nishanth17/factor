# C1 initial screen: results and limitations

**Initial screen only; final decision withdrawn.** The user reopened the
assessment because the short prefixes, excessive product bound and dense
diagnostic reservation do not justify a general DLP deferral. The separately
frozen follow-up tests source-informed bounds and larger collection windows.
C1, P5.4 implementation and R4 performance promotion remain open. No production source, default, API, checkpoint version,
portfolio policy, R3 provenance contract or `v1/` file changed.

Worktree: `/Users/nishanthmohan/.codex/worktrees/c1-dlp-siqs/factor`.
Branch: `codex/c1-dlp-siqs`, based on committed mainline `76f06b0`.
Protocol frozen in `74709c6`, initial runner `1db45cf`, current-range validator
repair `15cd091`/`762539c`. The original checkout and B3/A7 worktrees were
preserved. No merge or push is part of this task.

## Scope and evidence

The [protocol](c1_protocol.md) fixes six training cells: one certified balanced
input per 30/40/60-digit band, seeds 7/29. B1's selected 30/40-digit SIQS
configurations are controls; its 60-digit probe configuration is explicitly
uncalibrated for completion. Inputs and certificates were already versioned.
The [research review](c1_research.md) records primary papers, pinned YAFU,
msieve, Yamaquasi, FLINT, JavaMath and SymPy sources, licenses, blogs and
technical discussion. No external implementation was adapted or benchmarked.

The exact residual census runs beside unchanged SLP. With base bound B,
SLP keeps its B² prime cap; DLP allows each prime at most 100B and product at
most (100B)². All widened-threshold survivors are examined, including SLP
rejects. Each block also samples up to 16 uniformly chosen positions for
independent full-division checking. Rejection-stratum reservoirs hold at most
128 values each; failures and unresolved outputs remain in local captures.
This is stratified per-block sampling, not an unweighted random sample of all
positions: one-position polynomial tails have higher inclusion probability.
The last block of a censored cell is incomplete. No population-rate confidence
interval or threshold-loss estimate is inferred from these samples.

The splitter permits two 2,048-evaluation Brent attempts per composite, with
separate product, endpoint, certification, call-count and work bounds. Every
endpoint is PROVEN and its product checked. Root-hit recovery is cross-checked
with full factor-base division on 34,466 sampled positions across the six cells.
The full incidence oracle includes every component; independent exhaustive
small-matrix tests cover disconnected triangles, parallel edges, self-loops,
square corrections, both GCD signs and corrupted identities. Diagnostics never
enter production stores or checkpoints.

Raw data live in `v2/benchmarks/results/c1/screen/`; source hashes, immutable
protocol hashes and both driver hashes are retained. The original launch
correctly refused A7's lock. One subsequent control failed B1's stale >64-bit
label assertion after A10 expanded the proven range; the
[harness repair](c1_harness_repair.md) preserves that failure, adds a local
validator and continues only unfinished cells. No bound or gate was retuned.

## Residual economics and stopping reasons

The following counts are diagnostic, not accepted performance measurements.
Full/SLP/DLP records are mathematically checked residual shapes, not numbers of
useful factor-base dependencies. Seeds 7 and 29 share the same 40/60-digit
polynomial prefix; they are not independent input observations.

| Band | Audited blocks, seeds 7 / 29 | SLP records, including full | DLP records | Composite split attempts | Endpoint-bound rejects | Decision |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| 30 | 102 / 98, final blocks censored | 2,293 / 2,201 | 3,247 / 3,216 | 8,192 / 8,192 | 4,941 / 4,972 | Split-call cap and cost screen fail; full diagnostic matrix also refuses storage. |
| 40 | 512 / 512 | 1,951 / 1,951 | 2,344 / 2,343 | 6,231 / 6,231 | 3,884 / 3,885 | Recoverable population exists; complete diagnostic matrix is censored by its storage bound. Feasibility remains inconclusive. |
| 60 | 512 / 512 | 144 / 144 | 114 / 119 | 345 / 345 | 199 / 197 | Only one repeated DLP endpoint occurrence; no extra LP-cancelled constraint or useful dependency. |

At 30 digits there are 15,295/15,198 primes above the SLP bound and
19,343/19,186 product-bound rejections among widened/refined candidates.
Only 39.6%/39.3% of attempted composites fit the two endpoint bounds; four
per seed have more than two prime factors. No rho split fails in this band.
Splitting plus classification costs 0.385/0.296 diagnostic CPU seconds,
already above the frozen 25%-of-control investment screen. Those costs are
instrumented and include classification that a future design might avoid;
they are not a DLP-versus-SLP timing ratio or a proof against every splitter.

At 40 digits 13,014 primes exceed the SLP prime cap and 17,322 candidates
exceed the product cap per seed. Three rho attempts fail per seed; about 62%
of the composites split into endpoints outside the allowed range. Most valid
DLPs are lost at the existing SLP block threshold (1,678/1,677), with 337 at
refinement and 329 at residual policy. These are expected consequences of
SLP's smaller admissible domain, not missed SLP relations or a threshold bug.
The 128-to-512-block prefix grows DLP records from 605 to 2,344/2,343 and
SLP records from 495 to 1,951. Raw growth alone does not establish useful yield.

At 60 digits 911 candidates are primes above the SLP cap and 1,549 exceed
the product cap; 32/29 splitting attempts fail. DLP grows from 35/38 to
114/119 records between 128 and 512 blocks. Only three original full rows
cancel LP parity in either the SLP or SLP+DLP incidence matrix, and all three
are singleton-filtered away. Additional DLP constraints, verified dependencies,
square congruences and proper divisors are all zero. The two seeded prefixes
contain 368/378 distinct LP labels across 258/263 total records. Adding edges
mostly adds vertices. Threshold expansion recovers candidates but does not
solve this occupancy problem at the measured duration.

Rejected random-position samples have typical bit lengths near 51, 70 and
101–102 in the three bands, compared with product caps of 37, 40 and 47 bits.
These medians describe the stratified samples only. They do not classify large
unsplit residuals as primes or as semiprimes. Split failures are unresolved,
not proven non-semiprimes. The rejection samples include those failures and
large residuals, rather than retaining only recovered DLPs.

The full offline dense-incidence reservation is 255.6–265.9 MiB at 30 digits
and about 179.2 MiB at 40, above its declared 128 MiB cap. These are **diagnostic
oracle reservation refusals**, not observed production memory exhaustion and
not evidence for changing R3's bound or activating B7. The constants were not
lowered. Consequently these full DLP matrices have unknown post-filter yield;
zero must not be substituted for unavailable results. The user subsequently authorized a separately bounded follow-up, including
complete analysis of these retained populations and longer collection. The
original frozen gate remains reproducible; its failed extension criterion
is not a scientific argument against the newly authorized experiment.

None of the cells qualifies for the one permitted duration extension. At 30
digits the split cap stops collection; at 40 the diagnostic capacity refusal
prevents a complete feasibility result; at 60 endpoint repetition is below
the frozen extension criterion. No extra band, prime cap, splitter or seed was
introduced. Retained record reservations peak at 34.6 MiB and process RSS at
215.6 MiB; interpreter/JIT RSS is not owned-workspace accounting. The six cell
captures occupy about 13.8 MiB, below the finite retained-byte envelope.

## Subsequent complete retained-record analysis

After the user reopened the assessment, the separately frozen follow-up's
complete spanning-forest analyzer removed the dense LP-matrix artifact.
All original captures remained unchanged. Local `results/c1/retained-full/`
pins original record hashes, runtime/source hashes and driver commit `5f376c3`.
Eighteen targeted PyPy 3.11 tests passed, including generic GF(2) comparisons,
disconnected cycles, a 2,001-edge cycle, parallel edges, self-loops, exact
square corrections, finite refusals and unchanged SLP behavior/work.

| Initial records, seeds 7 / 29 | SLP LP-cancelled rows | DLP LP-cancelled rows | DLP dependencies | Verified complete factors |
| --- | ---: | ---: | ---: | --- |
| 30 digits | 121 / 104 | 316 / 288 | 92 / 64 | 880909969535437 × 902418641416693, both captures |
| 40 digits | 55 / 55 | 57 / 57 | 0 / 0 | None at the short prefix |
| 60 digits | 3 / 3 | 3 / 3 | 0 / 0 | None at the short prefix |

Every enumerated dependency was checked through original incidence, exact
exponents, square corrections and both GCD signs. Dependency counts include
trivial outcomes; the listed factorizations establish at least one nontrivial
outcome in each 30-digit capture. Same-position SLP-only incidence has no
dependency in all six captures. All real cycles in this population happen to
touch the SLP component; disconnected completeness is covered independently.
The revised diagnostic reservations fit within about 25 MiB per analysis.
This is a different sparse representation, not a reduction of R3's production
matrix constants. Single instrumented processing costs are not timing claims.

This result materially corrects the first screen: useful DLP relations already
exist in the current collector's rejected population. It does not establish
an end-to-end speedup, nor does the short 40/60-digit prefix settle delayed
matching. The longer, tighter-product follow-up remains necessary.

## Complete-factor evidence and limits

Uninstrumented SLP controls complete both 30-digit starts and both 40-digit
starts, with exact reconstruction and current PROVEN labels on both factors.
The 60-digit controls both stop at the 30-second wall allowance with the entire
input unresolved. They retain 47 rows each, one SLP match, 4,096/4,075 unmatched
partials and zero post-filter rows, closely reproducing B1's diagnosis.

Actual SLP matching/eviction differs from the optimistic no-eviction offline
oracle. At the short probe horizon, 30-digit SLP has 56/47 matches, 68/0
unmatched evictions, and 2,048/2,038 pending atoms. At 40 digits it has ten
matches, no eviction and 1,886 pending atoms; at 60, zero matches, no eviction
and 141 pending atoms. In the complete controls, 30-digit SLP evicts
1,402/1,770 unmatched atoms and 40-digit SLP evicts 5,098/5,801, yet all four
runs complete. Longer collection helps these controls; the short diagnostic
prefix must not be mistaken for their complete useful yield.

There is no production DLP factorization result, no fresh held-out DLP
confirmation and no accepted timing comparison. The controls and instrumentation
are a bounded feasibility study, not a performance claim; their one sample per
seed does not meet the >=3-second warmup/nine-sample promotion standard. No
historical native speedup, raw cycle total or cold/warmed mixture is promoted.
The revised roadmap policy is preserved: any future promotion needs matched
total resources, complete-run benefit with uncertainty and fresh confirmation.

## Diagnosis motivating the user-authorized follow-up

1. At 30 digits, an alternative finite cofactor policy must reduce wasted
   endpoint-bound splits enough to justify its complexity against an already
   completing SLP control. Asymmetric caps from the primary paper are a future
   preregistered hypothesis, tested explicitly in the follow-up, not retroactively applied to this screen.
2. At 40 digits, make a bounded complete offline assessment of the retained
   population feasible under an explicitly proved representation bound before
   collecting more data or building a production graph. Do not confuse this
   diagnostic representation issue with production matrix capacity. Establish
   useful rows/dependencies and whole-pipeline cost, not just affordable splits.
3. At 60 digits, demonstrate sufficient repeated endpoints and useful
   dependencies within an explicitly funded collection duration. The present
   screen supports neither extrapolation to completion nor impossibility of
   later DLP gains. Simply enlarging stores cannot create missing collisions.
4. Only a new justified gate can authorize production graph implementation.
   It must preserve proven residuals, original atom ownership under eviction,
   exact repeated-prime corrections, all-component cycles, independent lifting,
   cancellation and charged resume. Coordinate any shared schema with A7/B3.

## Reproduction and validation

```sh
pypy3 -B -m v2.benchmarks.c1_feasibility --output v2/benchmarks/results/c1/reproduction
make -C v2 test PYTHON=/path/to/pypy3.11-with-gmpy2
make -C v2 lint
```

The runner requires PyPy 3.11 and acquires the machine-wide performance lock.
Use a new output directory; saved captures are exclusive-create. Its fixed
protocol and corpus loaders work from committed files. Run heavy validation
inside the same lock, as coordinated with A7 and B3. The scoped tests verify
sampling/recovery, finite split refusal, the full independent incidence oracle,
exact extraction and the A10 certainty boundary. Final full and committed-only
validation results are recorded in the benchmark guide after they pass.
