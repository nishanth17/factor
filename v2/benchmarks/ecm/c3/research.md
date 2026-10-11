# C3 research and mechanism comparison

Inspected 10 October 2026, America/Los_Angeles. This is a mechanism transfer
review, not a native-factorer speed ranking. No upstream code was copied or
executed. Downloaded source bytes remain local under results/c3/research;
[URL/hash/failure receipts](../../inputs/controls/c3_research_sources.json)
are versioned for reproducibility. Existing source pins are deliberately reused.

## Probability and economic model

[Lenstra 1987](https://pages.cs.wisc.edu/~cs812-1/Lenstra1987.pdf),
§§2.5, 2.8–2.11, separates finite independent attempts from conjectural
smoothness-based cost. Cost depends primarily on the smaller factor, which
is hidden, and modular arithmetic on the full input. Complete factoring adds
recursive costs; a failed campaign proves neither primality nor absence of a
factor size. For a fixed hypothetical factor size and per-curve success q,
k independent curves have model success 1−(1−q)^k. For example q=1/100 gives
about 7.7% after eight curves and 27.5% after 32. That illustration is not a
fitted q for v2; the marginal value depends on the prior distribution and the
cost of the alternative engine. Survivorship changes that prior after failures.

[Zimmermann–Dodson, 20 years of ECM](https://members.loria.fr/PZimmermann/papers/40760525.pdf),
§§2–3, explains prime-power stage one, one-large-prime continuation, polynomial
stage two, blocking and Brent–Suyama extensions. Stage-two success and cost
must be fitted jointly with B1, memory and curve count. The paper's curve
family experiment has an [author erratum](https://members.loria.fr/PZimmermann/papers/)
under 2006: many curves at one prime do not establish average torsion across
prime residue classes. C3 consequently varies inputs as well as seeds and
retains existing independently verified curve/coverage oracles.

The [GMP-ECM parameter note](https://members.loria.fr/PZimmermann/records/ecm/params.html)
explicitly scopes its 2019 optimization to GMP-ECM 7 and a 512-bit modulus.
Its expected-count interpretation assumes an unknown target factor size and
roughly independent trials. At one expected count, the model miss rate is
about exp(−1), not a guarantee that all smaller factors were removed.
The pinned README advises increasing a target-size tier or moving to MPQS/NFS
and continuing against a composite remaining cofactor. Its tables differ
across versions/continuations; native expected counts do not establish PyPy
curves or a cross-engine work conversion. In v2, work is conservative policy
accounting; CPU/wall and actual search coverage must decide the economics.

## Implementation paths, assumptions and bounded transfers

| Pinned source / license | Actual caller and mechanism | v2 gap, transfer and rejection experiment |
| --- | --- | --- |
| [YAFU 8110dfbd8c6f9486d93b1a02de6eb7b180e55a80](https://github.com/bbuhrow/yafu/tree/8110dfbd8c6f9486d93b1a02de6eb7b180e55a80), factor/autofactor.c, schedule_work and get_next_state; public-domain file notice, external dependencies separate | Automatic scheduling computes target depth, credits completed bound/count effort and switches to a sieve. Explicit pretesting is separate. Source uses light 2/9, normal 4/13, deep 1/3 of input digits; docfile says 2/7 for light. Small target depth can skip ECM. These ratios express a workload prior, not known factors. Native small-ECM/GMP-ECM and calibrated sieve paths differ. | v2 already has finite tiers/cursors but no protected handoff. Transfer explicit effort modes and prior-work credit, not stale ratios or t-level tables. Compare no-ECM/short4/32 under equal complete-call budgets. Reject automatic adoption if uneven completion loses >5 points or confirmation is inconclusive. |
| [Yamaquasi 3f95f43682ed15d8c1ed206a9a702dd655d7c8ad](https://github.com/remyoudompheng/yamaquasi/tree/3f95f43682ed15d8c1ed206a9a702dd655d7c8ad), src/lib.rs factor_impl → src/ecm.rs ecm_auto/ecm_only, src/params.rs stage2_params; BSD-3-Clause | Auto mode tries a cheap ECM pretest then SIQS; explicit ECM-only escalates deeper. At 65–160 bits ecm_auto uses eight curves B1=200, B2≈7700. Above 370 bits it invests far more because SIQS is expected to be uneconomical. Native fixed-width Montgomery/Edwards arithmetic and FFT polynomial continuation change its costs. P−1 is marked once in shared preferences. | v2 lacks this explicit policy split. quick8 tests the mechanism at its lowest tier with existing Suyama x/z arithmetic, not an imported millisecond claim. Compare training complete outcomes and marginal factor yield against short4/32. Keep high-size escalation deferred without v2 feasibility evidence. |
| [FLINT a4c9750d0d3d67bb01cf6d18c187591b313451c3](https://github.com/flintlib/flint/tree/a4c9750d0d3d67bb01cf6d18c187591b313451c3), src/fmpz_factor/factor_no_trial.c → factor_smooth.c → ecm.c; LGPL-3.0-or-later inspected headers | A finite factor-bit table targets roughly one-third success, B2=100B1, with rough tuning above 62 bits. Tiered smooth factoring recursively removes found factors before a relation fallback. Native GMP limbs, not Python objects, determine arithmetic and storage. | Treat finite tiers and residual-aware termination as candidates. v2 already classifies children and shares Budget; never reset pretest spending on recursion. tiered/wide1 test 11000/50000 with original v1 B2 pairs and finite counts. Reject any broad transfer justified only by an upstream bit cutoff. |
| [PARI development 31092-e6893b0017](https://pari.math.u-bordeaux.fr/lcov-report/basemath/ifactor1.c.gcov.html), src/basemath/ifactor1.c ellfacteur and ifac_crack; GPL-2.0-or-later header, coverage snapshot 6 October 2026, locally hash-pinned | ellfacteur(N,insist) distinguishes a finite pre-MPQS path from deeper ECM. Normal/insist phases use disjoint curve seeds, and repeated calls on factors try fresh curves. ifac_crack reaches normal ECM, MPQS and later forced ECM according to flags. Stack-managed GEN arithmetic/native kernels and its continuation differ. | Preserve deterministic new curve assignments and checkpoint identity; do not interpret a changed bound as extending a saved curve. Campaigns have finite v2 tiers even though insist may keep escalating upstream. Check interrupted campaigns against uninterrupted seed/attempt prefixes and forbid config drift on restore. |
| [GMP-ECM 8ea5e214fdf2f0ddf9415141b8dc039ed6f5874e](https://github.com/sethtroisi/gmp-ecm/tree/8ea5e214fdf2f0ddf9415141b8dc039ed6f5874e), README, rho.c prob/ecmprob; LGPL-3.0-or-later library files, program licensing separate | Smoothness estimates include continuation/torsion assumptions; block count trades polynomial-stage time for storage. It is an ECM executor with explicit caller-selected bounds, not a complete SIQS portfolio. | Reuse existing independent coverage and finite program/plan accounting. Fit bounds/count/setup jointly; do not copy expected curves, polynomial benefit or memory coefficients. C3 does not reopen B2 geometry/pairing or P5.3 continuation. |

FLINT's `factor_no_trial.c` computes the smooth-factor goal from the current
cofactor as max(bits/3−17, 2). Its remaining <=128-bit path tries rho; at
>=118 bits it tries four additional 1000/100000 ECM curves before qsieve.
That caller, rather than the factor_smooth table alone, establishes dispatch.
The thresholds refer to native arithmetic and do not become v2 defaults.
PARI is pinned by the displayed revision plus downloaded content hash; it is
not represented as a fetched Git checkout. Library dependencies and file-level
license notices remain separate. Factor has no top-level license in this
checkout; public source availability is not permission to copy code here.

## Blogs, author notes and measured pitfalls

The [2016 YAFU user experiment](https://inaz2.hatenablog.com/entry/2016/01/14/230032)
preserves a 77-digit trace: multiple bound tiers accrue t-level effort before
SIQS. It illustrates real caller transitions; its hardware, native engine,
old revision and one input cannot calibrate v2. The current pinned source,
not that trace, establishes the policy path and documentation discrepancy.

[FLINT developer discussion #2625](https://github.com/flintlib/flint/issues/2625)
reports that continuing ECM toward a goal based on the original 420-bit input
wastes work after smaller factors reduce it to a 232-bit cofactor. The author
reports substantial improvement from earlier termination and later corrects
the implemented heuristic to bits/3−17. This is a measured pitfall and a reason
to test recursive portfolios, not a transferred timing claim or theorem.
The actual pinned caller and remaining-cofactor path must be read together.

Zimmermann's parameter note and erratum are author evidence for backend and
sampling dependence. Search also found a recent LLM-assisted general blog
claiming marginal-curve stopping; it is not used as authority, because it does
not establish the pinned source implementation or a reproducible comparison.
A relevant Mersenne Forum page failed retrieval; no unsupported forum claim
is included. Important mechanisms above trace to papers and source bytes.

## Existing evidence and reachable v2 behavior

B1 provides explicit calibrated balanced bundles and documented uneven and
structured regressions. A7/R5 preserves SSS/worker interface limitations and
final E1 ownership. C1's completing 40-digit half-base DLP is scoped; 50-digit
120-second censoring and a separately enlarged 141-second resume witness do
not establish a universal time cutoff or dispatch rule.

Direct review of C1's retained confirmation captures confirms median work
11,460,521,604 and 11,877,411,389 for its two fresh 40-digit DLP inputs,
with seed-dependent range 11,214,965,862–11,940,277,559. These are conservative
units, not billions of comparable ECM operations. The B3 bridge uses 32M
shared units to allow its 32 native curves; roughly 5.6M units on unresolved
cases cannot fit the library default 2M. The C3 pilot records actual current
configuration charges as a fresh reachability check.

B3's native default confirmation timing interval crosses zero; short curves
and early successes can fail to amortize certified chain preparation. A3
retains streamed defaults for small calls; B2 retains unpaired defaults after
losses. B4's accepted kernels and A6's finite exact p−1 executors are reused.
A10 primality labels and A11 explicit relation CLI remain intact. Existing
captures are training/regression evidence, never fresh confirmation.

The old portfolio eagerly builds a context for the largest configured endpoint,
then runs trial/powers/Fermat/rho/p−1/ECM before optional SIQS. Work/time refusal
halts inside the active stage; it does not execute fallback. CLI --siqs's 200M
work grant is still far below the C1 40-digit witness. Simultaneous fallback,
ECM plan, context and output memory must fit the explicit cap. A reservation
cannot manufacture a usable service allowance, so C3 exposes insufficient
fallback admission separately from schedule exhaustion and retains defaults.

## Preserved v1 policy and transfer limits

At the immutable integrated base `b3b3cfb`, `v1/ecm.py:compute_bounds`
chooses 2000/147396 through 30 digits, 11000/1873422 through 40,
50000/12746592 through 50, then 250000/128992510, 1000000/1045563762
and 3000000/5706890290 through 60/70/80 digits. Its final bound cap is
430000000/20000000000. `v1/constants.py` permits 10000 curves; the
`curves <= MAX_CURVES_ECM` loop can actually enter a 10001st attempt.
`v1/factor.py:factorize` recursively tries trial/rho/ECM and returns failure
when ECM fails; it has no SIQS fallback or cumulative protected allowance.
The original implementation uses Python 2 and floating logarithms/roots in
stage setup, so it is neither run as a PyPy 3.11 control nor copied into v2.

The useful hypothesis is the joint bound pair, not a digit oracle or the large
curve cap. C3 tests only the first three pairs, with bounded 8+2 or one-curve
campaigns under the same total service resources. Complete factoring and
marginal proper-factor yield can reject those transfers even if individual
ECM attempts succeed. Higher tiers would require a new finite feasibility
protocol; preserved v1 remains unchanged.

## Chosen tranche and deferrals

Implement cumulative optional-stage ceilings, work/wall/CPU reservations,
one-way handoff with partial-attempt evidence, exact shared recursive budgets
and schema-12 policy identity. Existing legacy schemas retain their executors.
The bounded six-arm protocol compares native costs/yield and source-grounded
finite tiers, freezes selection before fresh generation and retains the
baseline when evidence is inconclusive. No expected-curve table enters dispatch.

Defer a fitted factor-size posterior, arbitrary input-dependent deeper tiers,
new curve families, dynamic stage-two geometry, same-curve extension and
cross-engine work normalization. Each needs new measured costs/yield or a
separate arithmetic contract; none is necessary to provide a usable protected
handoff. G1/E1 retain broad multi-engine/population calibration.
