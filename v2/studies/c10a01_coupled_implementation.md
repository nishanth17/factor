# C10a-01 combined implementation handoff

The four nominated extensions are implemented as opt-ins on the existing
exact parity representation. Test one combined challenger first. Independent
component oracles establish correctness; selective ablations can explain a
promising or surprising bundle result. A factorial performance sweep is not
required for admission. No new timed calibration or complete-factor pilot
was run for this implementation.

## Implemented interfaces

| Component | Explicit selector | Contract |
| --- | --- | --- |
| Wider coupled ranking and residue reuse | `multiplier_policy="parity2_wide"`, `polynomial_policy="parity2"`, `multiplier=0` | All 61 squarefree multipliers through 97, including 20 even choices; 300 primes through 1987. Exact common-denominator scoring equals the direct rational oracle. Per-call symbol reuse reduces Euler tests from 18,239 to 7,774; this is an operation count, not a timing claim. Proper exceptional GCDs are recovered. |
| Actual A-product selection | `assignment_policy="paired_flyer"` | At least three A factors. Each deterministic core has two disjoint completion sets; an exact monotone scan chooses the closest product in their Cartesian product. Distinct cores give distinct A values. The finite population intentionally differs from `flyer`. |
| Guaranteed powers of two | `collector.marking="parity_two"` | Fuse the guaranteed one or two score levels into buffer initialization, retaining hit metadata, higher-power roots and conservative overmarks. Unsupported backend/policy/presieve cases use the existing exact route. |

The local symbol cache ends before factor-base setup. Reuse into base setup
was not added: it would need an additional ownership/identity contract and is
not necessary for the wider selector. Existing `standard`, `parity2`, `flyer`
and other explicit selectors keep their behavior. No defaults changed.

The combined literal configurations are constructed by
`v2.benchmarks.qs.c10a01.configuration_v3.configuration(digits, arm)`, with
`digits` in `(60, 62, 65, 70)` and `arm` equal to `"control"` or `"coupled"`.
The control equals the previous literal control exactly. Resource allowances,
base counts, intervals, LP domains, pacing, queues, plans and backends remain
matched. Only the disclosed coupled policies differ. Fixed multipliers still
bypass scoring; even or otherwise ineligible discriminants use standard-form
polynomials under the established arithmetic contract.

Jobs using any new selector serialize as SIQS schema 9 with
`coupled_semantics="c10a01-coupled-v3"`. Older schemas reject injection of each
new selector. Policy substitution on resume is rejected even when resource
extension is allowed. The A-stream identity is versioned separately; queue and
power-plan formats do not change. No new derived cache is serialized.

## Exact verification and structural results

The whole PyPy 3.11 suite passes: 859 tests, with nine optional-backend skips;
lint passes. Independent checks include direct target Euler/Fraction scores,
Cartesian pair minimization, original norm identities, complete score/hit
buffers versus the previous kernel, direct division, a 40-level two-adic
valuation, small proper-factor extraction, cancellation, pending activation,
charged resume and individual legacy-policy injection rejection. Existing
CRT, singular-root, queue, presieve, boundary and routing checks remain in use.
The committed-only verification receipt is retained with the local evidence.

Setup-only checks use all eight previously frozen exploratory inputs, without
sieving or passing their known factors to the algorithm. The resulting scored
choices and actual base endpoints are:

| Digits | Scored k, balanced / uneven | kN mod 8 | Base endpoint, balanced / uneven | New fixed-work / complete factoring |
| --- | --- | --- | --- | --- |
| 60 | 97 / 3 | 1 / 1 | 103788 / 102644 | Unmeasured / unrun |
| 62 | 13 / 7 | 1 / 1 | 127914 / 127104 | Unmeasured / unrun |
| 65 | 3 / 1 | 1 / 1 | 189020 / 188160 | Unmeasured / unrun |
| 70 | 1 / 29 | 7 / 1 | 237180 / 237234 | Unmeasured / unrun |

All bases attain the requested counts; all paired completion envelopes contain
the target. Three indexed assignments per input satisfy exact norm identities
at both interval boundaries and interior signed positions. This does not mean
every core admits a close completion. Sampled relative A errors range from
under 1 ppm to about 95.1%; the scan is optimal within its declared pair set,
not over all possible A products. Useful relation yield remains an experimental
question. Do not describe this change as a demonstrated A-quality improvement.

The scored balanced 70-digit cell still chooses standard fallback. A combined
comparison there tests the wider selection/A policy but does not exercise
normalization. An eligible forced diagnostic must be separately declared; the
old forced-k=31 result cannot be relabelled as the revised scored policy.
Even-k and residue-5 branches have exact correctness coverage; new performance
observations for those branches remain unmeasured.

## Resource and work amendment

The prior conditional relation/preparation/dense-fill bounds remain applicable:
there is no new relation payload, exponent convention, matrix representation
or solver. With at most ten primes below 250000, k at most 97 and the fixed
interval, the prior conservative 512-bit arithmetic bound still holds.
Additional persistent accounting relative to the parity2 prototype is
`32768 + (8192 + 64*F + 256*r) + 4096` bytes for score records, paired metadata
and the kernel. The separate 8 MiB selector scratch is released before base
construction; take the maximum of phases rather than adding it to the matrix
peak. These are ownership reservations, not exact physical PyPy RAM bounds.

| Digits | Extra persistent bytes | Matched owned / RSS tier | Whole-schedule pair-scan work bound |
| --- | ---: | --- | ---: |
| 60 | 367104 | 1 GiB / 3 GiB | 73,485,300,000 |
| 62 | 431360 | 2 GiB / 5 GiB | 98,983,500,000 |
| 65 | 591360 | 2 GiB / 5 GiB | 140,233,500,000 |
| 70 | 719616 | 4 GiB / 8 GiB | 192,131,700,000 |

The extra persistent amount is below 0.001 GiB at every size. The prior
central population scenarios retain more than 0.23, 1.08, 0.66 and 2.28 GiB
owned headroom respectively after the pilot parent's reserve. Those scenarios
are empirical assumptions, not bounds on relations needed to factor. Previously
infeasible simultaneous-capacity cases still refuse. No new RSS measurement
or admission is claimed. The 24 GiB host reported 68% available at this check;
headroom must be rechecked before any future stage, with cooperative-overshoot
limitations unchanged.

For N up to 4096 bits, the wider selector charges at most 14,008,178 units,
within its standalone 32,000,000-unit limit. The actual common denominator
has 2799 bits; the analytic allocation/charge proof uses 6601. The existing
10-billion-unit setup reserve covers this change. For a paired scan, at most
`F-1` visits each charge at most `18*r+3` units. The table bounds all 100000
families. A conservative future allowance can add that whole-schedule scan
bound to the earlier per-call limit; actual work is charged only once:

| Size | Per-call wall AND CPU grant | Proposed coupled/control common work guard |
| --- | ---: | ---: |
| 60 | 1,800 seconds | 1,769,476,365,809,652,160 |
| 62 | 2,400 seconds | 5,168,273,203,786,896,190 |
| 65 | 3,600 seconds | 25,835,282,003,227,202,268 |
| 70 | 7,200 seconds | 81,438,250,284,525,778,332 |

Including the existing warmup and analysis reservations gives prospective
stage work guards A = 13,875,501,469,193,096,700;
B = 51,670,565,176,454,404,536;
C = 162,876,501,739,051,556,664.
These are derived design values, not an admitted v3 execution freeze or
newly measured rates. The funded complete-factor time design remains A =
9,200 seconds, B = 7,600 seconds and C = 14,800 seconds, including the explicit
70-digit pair. With the original calibration ceiling the envelope remains
35,200 wall AND CPU seconds. A fresh v3 calibration allocation and source-pinned
protocol are still needed; no extra calls or time are silently granted here.

## Reproduction and remaining experiment work

From the task checkout, using PyPy implementing Python 3.11:

```sh
pypy3 -m unittest v2.tests.test_c10a01_coupled -v
make -C v2 test
make -C v2 lint
```

Serialize heavy QA and performance through `/private/tmp/factor-performance.lock`.
The original calibration ledger remains closed at 3,124 charged wall AND CPU
seconds and 138,632,524,390 charged work units. Its failures, receipts and
archives are unchanged. Development failures are also retained:
initial fixture assumptions, four stale compatibility-pin errors and one lint
issue were corrected before the passing run.

The old pilot runner intentionally rejects current source drift. Its freeze
remains usable in the historical source checkout, not as a v3 execution plan.
Next experimental work is a separately identified, bounded calibration of the
combined arm, then a new freeze for the funded four-size design. Include one
informative eligible 70-digit diagnostic if scored fallback persists. Choose
selective ablations from those results; do not automatically multiply arms or
repeat unsuccessful pairs. No defaults, merge, push, E1 or complete-factor
pilot is part of this implementation handoff.
