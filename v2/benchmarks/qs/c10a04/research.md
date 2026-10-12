# C10a-04: joint SIQS calibration, Stage 1

Starting mainline is `e0dafa1ebce2307f1f2ac3d787af89dd54f5ecfc`.
The active implementation is v2. This study consumes completed B1, C1, C9,
C10 and C10a-05 results; it changes no production defaults or arithmetic.
E1 remains deferred. All sources/configurations are separately frozen for
this study rather than modifying historical controls.

## Reconciliation with completed work

| Evidence / mechanism | Classification for this stage |
| --- | --- |
| B1 nearest/flyer, three-to-five-factor/Gray reuse calibration | Already covered at balanced 30 digits; uneven/structured regressions remain. Its 40-digit MPQS result missed its gate; larger censored probes remain unresolved. No new MPQS arm. |
| C1 half-base DLP | Accepted only in its balanced 40-digit regime; fresh 50-digit 0/18 under 120-second grants and fixed-store losses remain. This is not evidence against every smaller-base configuration or larger allowance. |
| C9 full-base aligned growing DLP | Accepted 50-digit bundle supplies the collector/control. C9a owns residual/splitter/graph changes. |
| C10 powers plan | Accepted 2 MiB plan is fixed in every arm; its old quality-16 ordering loss does not reject a different population algorithm. That algorithm belongs to C10a-02. |
| C10a-05 / plus | Collection pacing accepted at exactly 50/51 digits and 1 GiB. Larger pilots reused bound 30000/M=65536/six A factors: 55 digits completed in 421.872 s; 60 digits remained unresolved at 1200.042 s, with 352 rows and collection dominating. Both are revealed training, never fresh confirmation. |
| Count-based bases, odd multiplier scoring, interval/A coupling | Existing mechanisms; missing joint calibration at larger sizes. These are this study's hypotheses. |
| Normalized polynomials, new population/diversity, threshold feedback | Conditional C10a-01/-02/-03 work, excluded here. A1a retains new collector kernels. |
| Even multipliers | Uncovered opportunity. Existing explicit target construction is insufficient validation for a new scored branch; no even-k candidate arm is included. |

References: [B1 results](../../docs/studies.md#b1-joint-qsmpqssiqs-calibration--9-october-2026),
[C1 results](../c1/c1_implementation_results.md),
[C9 results](../c9/c9_growth_results.md),
[C10 reconciliation](../c10/candidate_reconciliation.md),
[C10a-05](../c10a05/results.md), [larger pilots](../c10a05_plus/results.md),
[QS parameter ledger](../../qs_parameter_audit.md).
The roadmap remains open until fresh experiment gates pass.

## Executable sources and representation

The [source manifest](../../inputs/controls/c10a04_sources.json) reuses C10's
pinned release/revision/file digests and licenses. Ignored local snapshots
were re-read and checked against those digests. The only new file is YAFU
`siqs_aux.c`, fetched at the same revision and checked against its GitHub tree
blob. Native sources were inspected, not executed or copied into v2.

| Implementation / pin | Executable coupling and limits |
| --- | --- |
| [msieve 1.53](https://sourceforge.net/projects/msieve/files/msieve/Msieve%20v1.53/), archive SHA-256 `c5fcbaaff266a43aa8bca55239d5b087d3e3f138d1a95d75b776c04ce4d93bb4` | `mpqs.c::get_sieve_params` interpolates base cardinality, LP multiplier and sieve size; `build_factor_base` selects k then fills an eligible count. `choose_multiplier` includes even k, a size penalty and classic mod-8 contributions (even D gets none). `sieve.c::do_sieving` rounds total positions into native blocks, passes the resulting half interval to `poly_init`, and derives LP allowance from the actual endpoint. `poly.c::poly_init` uses sqrt(2D)/M, available bit bins and A-factor count; `build_base_poly` corrects the last factor and `build_derived_poly` reuses A. Base-plus-excess collection and downstream filtering belong to its own solver. Public domain inspected code; dependencies separate. |
| [YAFU](https://github.com/bbuhrow/yafu/tree/8110dfbd8c6f9486d93b1a02de6eb7b180e55a80/factor/qs), `8110dfbd8c6f9486d93b1a02de6eb7b180e55a80` | `SIQS.c::siqs_static_init` calls `siqs_aux.c::get_params`, rounds B for SIMD, selects k through `choose_multiplier_siqs`, builds the base and computes the A target. `get_params` has SMALL_SIQS_INTERVALS, AVX512 and generic tables plus overrides: the backup table and commented historical tables are not the active universal policy. Main interval is native block size times block count. Active D=1 mod 8 code halves A target and uses Q2; its multiplier contribution is 2.625 ln(2), explicitly ad-hoc in source. `new_poly_a` has small-input and general paths, non-small-prime pools and corrective factors; `computeBl`/`nextB` reuse A. Comments about two-entry distance do not prove an independent overlap controller in every branch. Public domain headers; no table or 2.625 score transferred. |
| [Yamaquasi](https://github.com/remyoudompheng/yamaquasi/tree/3f95f43682ed15d8c1ed206a9a702dd655d7c8ad/src), `3f95f43682ed15d8c1ed206a9a702dd655d7c8ad` | `fbase::select_multiplier` / `expected_smooth_magnitude` searches an input-bounded range and assigns 3 to D=1 mod 8 because of Type2 representation. `siqs::fb_size` adjusts effective size for that representation and DLP; `FBase::new` rounds eligible count to multiples of eight. `interval_size` returns total positions in 32768-position blocks. `select_siqs_factors` uses mm/2 as M and a representation-dependent target; `nfactors`, `a_value_count`, `select_a` and `prepare_a` couple available primes, products and 2^(r-1) reuse. `siqs_calibrate` is a separate diagnostic path, not the normal caller. LP bounds derive from actual endpoint with separate product policy; `relations` owns its targets. BSD-3-Clause. |
| [FLINT](https://github.com/flintlib/flint/tree/a4c9750d0d3d67bb01cf6d18c187591b313451c3/src/qsieve), `a4c9750d0d3d67bb01cf6d18c187591b313451c3` | `qsieve_factor_with_tune` → `qsieve_init_with_tune` → `qsieve_knuth_schroeppel` → prime initialization → `qsieve_init_A`/`qsieve_next_A`. Even input first yields 2; the even-k score retains initial ln(2)/2 rather than msieve's zero. `primes_init.c` forms A target from sqrt(2kn)/(sieve_size/2). Count-oriented prime storage and A bit bins constrain factors; the initializer admits products inside a target window and uses corrective factors. Processing retries can collect further rows; the outer failure path grows prime count by 10%. That growth changes the base and is not copied into the fixed-base v2 checkpoint contract. LGPL-3.0-or-later. |
| [Alpertron](https://github.com/alpertron/calculators/blob/93c5c8189cb2f149c8dae898eb11997c5f2e7980/siqs.c), `93c5c8189cb2f149c8dae898eb11997c5f2e7980` | Active C `FactoringSIQS` computes bounded count and half-interval formulas plus A-factor count, selects from 1/2/odd-prime multipliers, then builds eligible primes. `getFactorsOfA` samples a span around the target's rth root with within-A spacing; `sieveThread` reuses Gray variants. Relation target is max(33B/32,B+50), with finite native matrix/partial arrays. Packed lanes and root kernels are A1a opportunities, not calibration changes. `OldApplets/Siqs.java` and disabled SQUFOF are historical/optional rather than this active path. GPL-3.0-or-later. |

[Alford–Pomerance](https://math.dartmouth.edu/~carlp/implementing.pdf), §§3–4,
and [Carrier–Wagstaff](https://homes.cerias.purdue.edu/~ssw/qs4.pdf), §§2–3,
support A scaling and amortized family construction; their constants do not
predict PyPy optima. [Kechlibar](https://www.karlin.mff.cuni.cz/~krypto/Implementace_MPQS_SIQS_files/main_file.pdf),
§§2.3.8–10, distinguishes product quality and uniqueness. The concrete gap
addressed by [Percival's author note](https://www.daemonology.net/blog/2005-11-30-quadratic-sieve-constant.html)
is initialization cost versus interval/base size: it motivates comparing two
intervals on the same enlarged base, not copying its claimed native gains.
Existing C10 literature/supplemental sources answer the other gaps; no extra
implementation is needed merely to enlarge the bibliography.

## Units, invariants and hypotheses

v2 uses D=kN, U=Ax+B, F=(U²-D)/A, with exact reconstruction U²-D=A F.
This is the classic representation; Yamaquasi's normalized form has denominator
4A and YAFU's active Q2 branch changes its target. Their score bonuses cannot
be transplanted independently. v2's M is a half-width, covering 2M+1 positions;
block_width=batch_width=4096 is separate. Its A target is
`max(1,isqrt(2D)//M)`. Count means eligible prime entries including 2 and
singular multiplier primes, excluding the matrix sign column; endpoint means
the largest actual prime. `base_bound` is exclusive.

The study helper chooses the existing odd squarefree multiplier, scans a
finite prime sieve using exact Euler tests and proper gcd checks, then returns
endpoint+1 for an exact eligible count. It charges the scan and its 16,008,192
byte reservation inside inclusive cost. SIQS builds/validates the actual roots.
No production arithmetic, API, checkpoint format or default changes.
Exceptional p=2, p dividing k and p dividing A retain existing kernel branches;
independent exhaustive residue/count tests cover the helper. Even-k helper
counts are tested arithmetically, but this does **not** establish an even-k
scoring/root/relation performance arm. Full branch validation remains a gate.

| Arm | Falsifiable hypothesis and potential loss |
| --- | --- |
| control | Existing 50-digit base may be too small at larger N. Increasing time alone may still suffice; retain this arm at every matched grant. |
| compact | Odd k, 2500/5000/8500 eligible primes, 7/8/9 A factors, pool 64, M=65536: richer smoothness and longer A amortization should improve useful row/dependency formation enough to repay larger root/marking and matrix costs. Reject for a band if complete/capped cost and coverage fail; raw rows alone cannot win. |
| wide | Same count/k policy and pool/factor count, M=131072/131072/262144: fewer transitions per position should repay larger F values, marking, power-plan fallback and recovery. Relative to compact this tests interval/A-target coupling; relative to control it remains a whole joint bundle. |

The native size trends justify trying moderate increases; the numeric counts
are local hypotheses, not copied optima. Gray quota 64 retains full reuse at
seven factors and a bounded existing diverse subset at eight/nine. Actual
control quota is 32 at six factors despite configured 64. Each arm retains
100000 families and its exact combinatorial reachable space; count changes do
not imply that all theoretical A products are visited. The flyer pool algorithm
is unchanged; enlarging its existing pool is not C10a-02 population selection.

All arms retain absolute SLP/DLP endpoint 3,000,000, product 115,200,000,000,
scoring/admission/splitter, growing storage, 1 GiB owned memory, 1 MiB checkpoint,
32-row excess/growth, weight-two and collection pacing. Relative LP/base
ratios consequently shrink as bases enlarge; this is a deliberately frozen
C9 boundary and may limit higher bands. No residual-bound performance claim
follows. Growing counts avoid a new count eviction policy, but matrix and
owned-byte refusals remain possible. Peak RSS is process high-water, separate
from owned reservation. No capacity increase is silently added.

65 digits probes beyond the revealed pilots; it is a finite informative band,
not a universal upper size. Partial graph rows may all disappear during
singleton filtering, so Stage 1 screening records actual filtered matrix,
cycles, dependencies, collection cost, finite search and storage stops.
[Protocol](protocol.md) and [handoff](handoff.md) freeze the mechanical next
stage; no accepted training or fresh campaign is run in Stage 1.
