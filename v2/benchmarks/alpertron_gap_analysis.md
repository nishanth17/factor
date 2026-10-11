# Alpertron repository comparison and transfer phases

Date: 10 October 2026. Factor control: `b3b3cfb`.
Upstream: [alpertron/calculators][upstream], commit
`93c5c8189cb2f149c8dae898eb11997c5f2e7980` (13 September 2026).

This is a repository-wide inventory and factoring-path source comparison, not
a line-by-line correctness proof of every calculator. The downloaded archive
contains all 307 Git tree blobs: 65 C sources, 22 headers, 22 historical Java
sources, 39 JavaScript files, six assembly files and the remaining build,
documentation, web and media assets. The ordinary `rg --files` listing has 304
entries because it omits three hidden files. Build manifests, preprocessor
guards and callers distinguish active implementations from dormant experiments
and historical applets. No upstream program was executed or code copied.
There are no new runtime measurements or claimed percentage improvements.

The comparison uses v2's actual equivalents: exact native/GMP arithmetic,
accepted PRAC/Lucas-chain integration, Suyama ECM, paired stage-two challengers,
Gray-code SIQS, fixed scores/powers/buckets/resieving, certified DLP with an
all-component forest, singleton/weight-two filtering and Python-int elimination.
The [QS source reconciliation](qs/qs_gnfs_research.md),
[QS parameter ledger](qs_parameter_audit.md) and
[cross-engine ledger](parameter_audit.md) retain earlier sources and decisions.

Scheduling update: C3/C9 are already in progress and retain their original
scope. New source-derived work is separately tracked as C3a/C9a, with A1a,
A2a, A7a, A10a, A11a and C1a for deltas related to completed parents. Each
follow-up first reconciles the parent's final results and implements only
uncovered scope; already-covered/rejected/deferred decisions avoid duplicate
experiments. AP identifiers remain source-candidate IDs, not task IDs.

## Source coverage

| Area | Sources inspected or inventoried | Disposition |
| --- | --- | --- |
| Integer factoring | `factor.c`, `ecm.c`, `siqs.c`, `siqsLA.c`, associated headers and `commonstruc.h` | Trace setup, scheduling, collection, recovery, filtering, solving and factor construction; candidates below. |
| Shared arithmetic | `bignbr.c`, `bigint.c`, `modmult.c`, `MontgomeryMult.c`, `ClassicalMult.c`, `karatsuba.c`, `fft.c`, `fftUtil.c`, `division.c`, `sqroot.c`, `ntt.c`, assembly and WASM generator | Separate native kernels from transferable algorithms; ordinary integer GCD is binary GCD. |
| Other arithmetic/calculators | Polynomial/rational modules; `gcdrings.c`; discrete logs, quadratic equations, continued fractions, sums of squares/cubes, partitions and Ulam graphics | Polynomial kernels are conditional future references. Other solvers use different mathematical domains or consume the factoring engine; they are not missing generic integer-factor algorithms. |
| Product/input | `expression.c`, `parseexpr.c`, `batch.c`, `inputstr.c`, `ecmfront.c`, `fromBlockly.c`, linked-number/output modules, browser workers, storage and service worker | Batch input, hints, expressions and browser UX are optional capability work, separate from factoring speed. Existing v2 checked resume remains the control. |
| Historical Java | `OldApplets/ecm.java`, `ecmc.java`, `Siqs.java` plus expression/polynomial and other applets | Additional Fibonacci/Lucas special forms, APRT-CLE and real SIQS thread coordination; explicitly historical. |
| Build/validation/assets | Makefile, BAT/Visual Studio manifests, tests/coverage, GitHub workflows, license, HTM/manifests/media | Establish reachability and license; no arithmetic gain inferred from assets or CI tooling. |

## Phase 1: yield/allocation deltas after running C3/C9

Finish the existing C3/C9 scopes first. C3a/C9a consume their final controls
and investigate only source-derived gaps those sessions did not cover; they are
optional and do not become new parent acceptance requirements. C1a is a separate
anchor-only follow-up to completed C1. C10 can fit remaining collection bundles
on the settled C9 control without requiring these deltas to win.

| ID / owner | Source difference and v2 gap | Benefit hypothesis and acceptance boundary |
| --- | --- | --- |
| AP01 — C1a after C9 | [SLP anchor reuse][partials]. Alpertron keeps the first equal-cofactor partial after a match. v2's SLP `_admit` deletes the pending partner; the DLP forest already retains anchors, and C1 already evaluated graph-SLP reuse. | r distinct equal-residual atoms permit r−1 anchor pairs instead of floor(r/2) disjoint pairs. These span the LP cancellation space before factor-base filtering; they do not guarantee more useful dependencies. A cheaper anchor-only implementation is a follow-up candidate, not a previously untested yield mechanism. Reuse the C1 graph-SLP result; compare bounded anchor lifetime, row/atom caps, duplicates, eviction and checkpoint ownership only if a changed cost/retention hypothesis justifies it. |
| AP02 — C9a; B7 for representation proofs | [Compact replayable partials][partials]. Source stores residual/sign, hash link, U=abs(Ax+B) and A-generation seed, then reconstructs the old smooth factors when matching. v2 retains full atom/exponent payloads with conservative reservations. | More pending opportunities within the same real memory cap may reduce eviction loss. Prove reconstruction, identity, certainty and corruption checks; charge replay and simultaneous scratch. A smaller reservation constant alone is not an implementation. Preserve independently checkable original relations. |
| AP03 — C9 control; C9a only for uncovered deltas; C10 coordination | [Residual opportunity and retention policy][setup]. Source uses an LP bound 100 times the largest actual base prime, equal-cofactor hashing and a much larger fixed partial array. | Use as a normalized control for endpoint/product/admission/retention/duration studies, including opaque equal-cofactor matching only through a separate verified type. Larger arrays are not proof that capacity is the sole bottleneck. Source does not provide an active DLP graph or a stronger splitter. |
| AP04 — C3a after C3; C11 coordination | [ECM factor-size effort table and ECM-to-SIQS policy][ecm-policy]. Escalating curve counts/bounds contrast with v2's bounded fixed defaults. | Fit finite tiers and marginal success using current engine costs, credit completed work and reserve handoff resources. Source's total-input-digit switch and B2=100B1 are candidate policies, not optimal Python constants. Factor-size labels stratify evaluation, never feed the dispatcher. |

### Why AP01 is a real difference

For distinct atoms e0,...,e(r−1) sharing a residual, pairs e0+ei for i>0
form a basis of the even-weight subspace of the residual-incidence row over
GF(2). Each non-anchor atom occurs in exactly one basis pair, proving
independence at this stage. With four atoms this gives three pairs rather than
two. Factor-base parity can still make rows redundant, and the source later
discards some equal-parity/empty rows that v2 correctly preserves. Retain v2's
stronger dependency handling. This is a candidate improvement to SLP, not an
unimplemented feature of the accepted DLP forest.

The [completed C1 implementation study](qs/c1/c1_implementation_results.md)
already compared graph-only SLP reuse: its two 40-digit training runs totaled
20.055 seconds versus 20.160 for ordinary SLP. These single-run diagnostics
establish no repeatable speedup. AP01 therefore means a potentially cheaper
SLP-specific implementation or a changed retention regime, not repeating the
same graph arm or claiming that the algebra alone supersedes its result.

## Phase 2: collection quality and kernels

C10 owns diagnosis and AP10–AP14 policy studies; A1a owns selected AP05–AP09
collector-kernel implementations after completed A1. E6 reuses those results
for shared context integration. C8 retains its larger CRT half-sum scheduler:
A1a modular root-shift tables have a distinct root-cost/setup/storage gate and
do not require eligible sparse-hit primes.

| ID / owner | Source-derived candidate | v2 comparison and required boundary |
| --- | --- | --- |
| AP05 — A1a; C10 diagnosis/E6 integration | [Periodic tiny-prime presieve][sieving]: build a short score pattern and replicate it across the interval. | Generic score translation already exists; periodic multi-prime initialization is the gap. Verify all residues, powers, saturation, block offsets, setup charges and resumed positions. |
| AP06 — A1a; separate from C8 CRT scheduling | [Precomputed Gray-transition root shifts][roots]: Bainv2 for each CRT term and base prime. | v2 already caches A inverses and shifts roots, but computes the new delta/inverse product per prime. Bound the extra table and account for centered-B wrap corrections, p=2 and p dividing A. |
| AP07 — A1a; C10 diagnosis/E6 integration | [Two related polynomials in one sieve pass][sieving]: two byte lanes packed into a short. | v2 processes polynomial windows separately. Test whether shared traversal saves enough on PyPy; prove lane isolation/carry behavior and independent conservative scores, roots, work and pause cursors. Native short packing is not automatically the right Python representation. |
| AP08 — A1a; C10 diagnosis/E6 integration | [Bulk candidate extraction][scan]: inspect high-bit masks across groups before individual positions. | v2's per-position scan is the comparison; numthy already motivates a related candidate. Handle clipped thresholds/tails and preserve exact candidate identities, position charges and interruption behavior. |
| AP09 — A1a; C10 diagnosis/E6 integration | [Prime-density-specific marking][sieving]: unrolled short-prime loops, one/two-hit large-prime cases and grouped integer log weights. | v2 already has dense/sparse/bucket paths. Evaluate only the missing fused/specialized traversal and weight reuse, including table/build cost; do not present generic buckets as new. |
| AP10 — C10 | [Joint size-derived collection policy][setup]: base cardinality, interval, A-factor count and family size. | v2 has explicit calibrated bundles, not a proven size optimum. Normalize prime count versus prime-value bound, whole interval versus block, and all row/atom/memory caps. Compare useful filtered rank and complete factors. |
| AP11 — C10 | [Broader multiplier search][setup]: candidates include 2 and primes through 97; scoring examines many more test primes. | v2's optional selector tests 13 odd squarefree candidates through 31 against primes through 97. Compare search cost and resulting smoothness, with exceptional roots and even multipliers verified before use. |
| AP12 — C10 | [Scoring and omitted-hit policy][sieving]: heuristic byte thresholds and omission of A-prime marking. | v2 already has fixed scores, prime-power plans and explicit threshold controls. Source-specific tiny/A-prime omission and grouped-log bundles need a new declared accuracy/loss study; never describe heuristic missed hits as exact rejection. |
| AP13 — C10 | [A-factor pool/spacing policy][roots]: narrow randomized pools and separated sorted indices. | v2 already supports bounded A selection and unique assignment streams. Compare polynomial quality and family coverage; the new candidate is pool/spacing policy, not CRT/Gray generation itself. |
| AP14 — C10 | [Relation targets and solve cadence][setup]: max(33/32 times base cardinality, base+50) before LA. | v2 already exposes excess/growth controls. Train oversampling and retry cadence using post-filter rank, trivial dependencies and extraction cost. Do not copy a target tied to a different solver. |

## Phase 3: matrix/storage work when justified

These are genuine missing implementations, but the existence of native sparse
solvers does not establish immediate benefit. Refresh A8 when settled C9/C9a/C10 produce
different matrices or expose a storage refusal. Retain B6's first dense/hybrid
comparison; if a representation prevents that comparison, B7 may precede it.
Do not require implementing every dense challenger before making a measured
decision that sparse solving is needed.

| ID / owner | Source-derived candidate | Required boundary |
| --- | --- | --- |
| AP15 — B7 with C9a | [Normalized square-root payloads][insert]: divide out square contributions on admission and retain sparse parity plus a modular square-root value. | v2 keeps exact exponents/original atoms and square corrections. A compact representation needs independent original-identity verification, retained lifting evidence, nonunit factor recovery and simultaneous memory proofs. Do not discard empty/equal-parity dependencies as the source can. |
| AP16 — A8/D1/D2 | [Packed sparse matrix-times-block products][la]: one machine word holds 32 independent GF(2) vectors. | Compare Python-int bitset elimination with bounded incidence/CSR/tiled scatter/gather products, both orientations and useful block widths. Alpertron stores relation rows; translate orientation explicitly into v2's constraint-row/relation-column M. No dense Gram construction. |
| AP17 — D1, alternative D2 | [Block Lanczos][la]: selected nonsingular subblocks, recurrence history and final correction against the original operator. | A genuine solver gap. Seeded finite retries, singular-block handling, original M d=0, lifting and checked resume are mandatory. Source's unbounded retry paths are not transferable. Block Wiedemann is absent here; study CADO separately rather than attributing it to Alpertron. |

## Phase 4: preprocessing and proof coverage

E3's existing frozen-control gates remain; A10a owns new classification
implementation deltas after completed A10. Special-form detectors can be a
separate opt-in tranche with structured and generic populations; a very large
gain on a recognizable family does not imply broad random-semiprime benefit.

| ID / owner | Source-derived candidate | v2 comparison and required boundary |
| --- | --- | --- |
| AP18 — E3 | [Bounded multiplier Fermat/Lehman search][lehman] with congruence/square-residue screening and incremental differences. | v2 has finite plain Fermat, default off. Charge tested multipliers/positions, prove exact square witnesses and compare near-square, moderately unbalanced and generic controls. Do not run the source pretest before every curve without measuring its tax. |
| AP19 — E3 special forms | [a^k±1 detection, Cunningham/cyclotomic and Aurifeuillean splits][special]. | New structured-input capability. Exact identity, proper divisor and reconstruction checks; bounded recognition and coefficient work. `PowerPM1Check` is not Pollard p−1. |
| AP20 — E3 special forms | [Sums/differences of powers][special-sums], including small-exponent algebraic identities. | Additional structured families beyond perfect powers; retain exact root/equality certification and bounded recognition. |
| AP21 — E3 special forms | [Fibonacci/Lucas recognition and index-divisor splitting][java-special]. | Historical Java only; absent in current C factor.c and v2. Use exact recognition/fast-doubling or checked supplied family metadata; do not port floating-log rejection as an exact guarantee. These Lucas numbers are distinct from B3 scalar-chain work. |
| AP22 — E3 | [Shared/adaptive modular perfect-power rejection][power]. | v2 already has exact integer Newton roots and a fixed modulus screen for several exponents. Candidate is multi-exponent shared screening and extra certified moduli, with setup/crossover costs. |
| AP23 — E3/G1 | [Per-cofactor trial-coverage inheritance][trial]. | Source carries excluded-prime bounds through splits; v2 restarts the trial cursor after reclassification of a reduced cofactor. Prove exclusion inheritance and serialize checked coverage, including partially committed chunks. |
| AP24 — E3 | [High-valuation extraction by exponent doubling][trial]. | Compare p,p²,p⁴,... with v2's repeated exact division where multiplicities are high; account for existing perfect-power detection and failed product tests. |
| AP25 — A10a; E3/C11 diagnosis | [Recover factors from primality-test residues][carmichael]. | v2's boolean strong-MR helper discards intermediate witnesses. Reuse nontrivial square roots of 1 or checked root collisions; source's extra 20-base search is a separate cost policy. Preserve certainty/RNG/work/checkpoint contracts. |
| AP26 — A10a; E3/C11 diagnosis | [BPSW-style composite prefilter][bpsw]: base-2 strong MR plus strong Lucas-Selfridge. | v2 has supported deterministic MR ranges and requested random rounds above them. An optional composite rejection path may reduce expensive later rounds; a BPSW survivor remains probable and cannot replace a supported proof or requested rounds silently. |
| AP27 — A10a; E3/C11 diagnosis | [Small-prime groups per remainder][bpsw]: three-prime products before expensive tests. | A10 already evaluated broader primorial-filter ideas. This is a narrower changed-workload candidate, not a new generic primorial discovery or reason to reverse that verdict automatically. |
| AP28 — later B14 scope | [Historical APRT-CLE/Jacobi-sum primality][java-proof]. | Additional proof-method coverage beyond planned n−1/Pocklington, not a new requirement for B14's first tranche. Source falls back to probabilistic testing through the same return interface: separate proof/probable outcomes and require an independent exact verifier and finite budgets. |

## Phase 5: conditional ECM/arithmetic and product ideas

| ID / owner | Source-derived candidate | Disposition |
| --- | --- | --- |
| AP29 — C2 with C11 | [Affine stage-two difference terms][ecm-stage2] versus v2 projective cross-products. | Source performs individual inversions. Batch inversion is a derived challenger, not an upstream feature; no-inversion common-Z is another distinct existing C2 candidate. Charge inversions, detect nonunits and retain local mixed-factor recovery. |
| AP30 — F2 | [Alternative rational Montgomery curve/point family][ecm-curve]. | Source uses u=2t/(3t²−1), A=(1−6u²−3u⁴)/(4u³), point (3u²+1:4u), unlike v2 Suyama. Prove exclusions/maps and compare whole-engine success; no torsion or yield advantage inferred from the formula alone. |
| AP31 — E5/C11 | [Stage-one/stage-two GCD cadence and replay][ecm-stage1]. | Source's coarse product/final-GCD approach may reduce checks but worsens detection/replay granularity. Existing v2 chunked recovery is stronger. Revisit only from profiles, preserving finite polling and exact saturation recovery. |
| AP32 — A2a after a C11 profile | [Batched binary/Montgomery inverse][inverse], Kaliski/Savas-Koc style coefficient updates. | Different from v2's native Euclid or GMP inversion. Benchmark inversion-dominant supported sizes including nonunits and conversion; do not replace native bigint algorithms on operation counts. |
| AP33 — A2a; C10/C11 diagnosis | [Small-integer Euclid quotient shortcuts][inverse-small]. | Source special-cases tiny quotients. Low priority: compare actual prime-root/setup use against native/GMP inverse, including dispatch cost. |
| AP34 — C11/C5 deferred native reference | [Fixed-limb unrolled multiply/reduce][native], assembly/SIMD and shared scratch layout. | These are C/backend advantages. Python ints and GMP already supply native arithmetic. Prior reducer/backend losses stand; do not implement a Python limb engine merely to copy C loops. |
| AP35 — C11 far-later arithmetic reference | [Cached constant FFT transforms][fft] and [three-prime CRT NTT][ntt], specialized squaring and unbalanced multiplication. | Relevant to much larger arithmetic/native whole kernels. NTT is included in the Visual Studio project, absent from Makefile ECM sources, and has no caller in this snapshot. Floating FFT needs an exactness proof; integer/GMP multiplication remains the control. |
| AP36 — F3/future GNFS only if required | [Polynomial kernels][poly]: Karatsuba/Kronecker/FFT multiplication, cached reciprocal reduction, Newton division, modular factorization, Hensel lifting and integer LLL/van Hoeij factorization. | These serve a polynomial calculator, not an existing polynomial ECM or GNFS factoring engine. Reuse a specific required mathematical primitive only after its field/coefficient contract is established; no reason to build all this now. |
| AP37 — A7a; F4 integration | [Historical SIQS family partitioning and private root/sieve state][java-workers]. | v2 already has bounded coarse SIQS worker experiments. Source gives an alternative family-sharing granularity, not proof of an unimplemented generic parallelism capability. Charge synchronization, shared stores, aggregate work/RSS and cancelled work. Current C `numberThreads=1` is not a working multithreaded comparator. |
| AP38 — A11a optional product subtranche | [Caller-supplied known-factor hints][known] and source Cunningham's optional external known-factor lookup. | Validate every hint by gcd/division and reconstruct unresolved cofactors; external lookup would require an explicit opt-in with bounded network/cache behavior. Never use hidden corpus factors in measured comparisons. Keep outside completed A11 scope. |
| AP39 — A11a optional product subtranche | [Batch and expression input][batch]. | Convenience for structured input/reproducible campaigns; not an arithmetic speedup. Bound parser/exponent/output work and preserve per-item failures, certainty and budgets. |
| AP40 — A11a optional product subtranche | [Browser progress, save/continue, worker and offline UI][front]. | Useful only if a browser product is desired. Existing v2 checkpoint/resume remains the exact-state control; source curve/factor persistence is not a stronger checkpoint contract. |

## Ideas already covered, weaker, dormant or unrelated

- Ordinary integer GCD in `bignbr.c:740` is binary GCD, not Lehmer or half-GCD.
  v2 already calls native `math.gcd` or GMP. Gaussian/quaternion GCDs in
  `gcdrings.c` solve different domains; they are not faster integer GCDs.
- PRAC, prime-power stage one, wheel sieving and ± stage-two pairing already
  have v2 implementations/experiments. Source wheels 210/2310 and native
  cutoffs are calibration inputs, not new algorithms. Do not repeat completed
  chain search or reverse a prior complete-run loss without a changed trigger.
- CRT/Gray families, log sieving, root-guided division, sparse hits, buckets,
  resieving and prime-power handling are already present. The narrower missing
  kernels are AP05–AP09; source floating limb remainders require their own exact
  bounds if ever adapted.
- Source SIQS matches equal residuals; it does not provide the accepted v2
  all-component DLP graph. SQUFOF and its prospective intake call in current
  `siqs.c` are under `#if 0`; the Java SIQS defines SQUFOF without an intake
  caller, and Java ecm's call is a debugging benchmark. Not an active production
  splitter advantage. C9 or an uncovered C9a delta may compare a verified active SQUFOF source.
- Source singleton removal is less extensive than v2 singleton/weight-two and
  zero-dependency handling. Its parity-only duplicate/empty-row rejection can
  discard dependencies; extraction checks only one sign. Preserve v2's original
  identities, both GCD signs and stronger recovery.
- No Block Wiedemann or Four Russians implementation was found in this tree.
  No generic GNFS engine was found. Discrete-log, quadratic-equation, continued-
  fraction and sum-of-squares/cubes calculators are different capabilities, not
  missing stages in generic integer factorization.
- Coverage examples can supply independently validated regression inputs, not
  trusted arithmetic oracles or independent held-out benchmark populations.
  Media, localization, web packaging and static-analysis integration establish
  no factoring performance improvement.

## Wider implementation comparison and sparse-solver decision

The next source pass should extend existing pinned research, not restart it:

| Source | Highest-value comparison | Existing owner |
| --- | --- | --- |
| [GMP-ECM][gmp] | Whole curve/stage-one packages, stage-two geometry/memory selection, Brent–Suyama/polynomial continuation and expected-effort policies; separate native modular kernels from Python-transferable algorithms. | C3a/C11 after C3; C2/F2/F3 when their triggers pass. |
| [YAFU][yafu] and official [msieve 1.53 release][msieve-release] | SIQS polynomial/marking/resieve paths, compact relation stores, residual screens, tiny cofactor factoring, cycle/filtering economics and stopping criteria. Compare source behavior beyond constants already in the ledgers. | Settled C9 control, C9a uncovered deltas and C10; B7 for changed storage. |
| [msieve Block Lanczos][msieve-la], Yamaquasi and Alpertron | Packed sparse products, recurrence/singular blocks, filtering/recovery costs and actual matrix memory. The linked msieve mirror is historical; reconcile with the official release before copying. | A8 then selected B6/B7/C4/D1 work. |
| [CADO-NFS][cado] | Genuine Block Wiedemann pipeline, block generator, independent sequences, reconstruction, checked checkpoints and sparse layouts. | D2; future GNFS remains separately deferred. |

Block Lanczos is the more direct first sparse-solver implementation candidate for
this single-machine QS project. Block Wiedemann is a peer alternative, especially
if independent Krylov work, restart/I/O or future scale is decisive; it does not
depend on implementing Lanczos first. Neither is a guaranteed earlier win than
collection. For a matrix with k retained relations and s nonzeros, dense bitsets
plus full provenance can approach quadratic storage, whereas sparse block
methods keep sparse incidence and a bounded number of block vectors, plus their
method-specific history/reconstruction state. This can make a previously
infeasible matrix feasible, even when solver CPU is not dominant.

Refresh A8 with representative post-filter matrices, original relations/lifting
maps, solve/filter/extract times, observed fill, provenance/storage peaks and
refusals. A capped synthetic capacity stress test is useful for future feasibility
but cannot establish current full-factor speed. If measured solve fraction is f,
eliminating all solve time alone saves at most f of current total runtime; a faster
solver cannot create missing relations. Conversely, a proved matrix-memory
refusal justifies a storage/sparse-solver investigation without waiting for solve
time to dominate. These are the explicit triggers for advancing D1/D2.

The [C1 evidence](qs/c1/c1_implementation_results.md) supports the current
collection-first placement: all five arms exhausted 120 seconds on both
50-digit training inputs; full/narrowed DLP and SLP on the first input filtered
to zero, while half-base DLP retained 466 rows/711 columns without a kernel.
At 60 digits every LP-cancelled row in the two census captures was removed by
singleton filtering. The one instrumented 50-digit tail was dominated by
sieve/root work, not graph admission; it is diagnostic, not accepted timing.
The separate larger-store 50-digit resume witness completed in 141.063 seconds
with changed allowances, so it does not isolate duration, retention or solver
benefit. Do not extrapolate these captures into a universal size cutoff.

All accepted changes retain exact arithmetic, proper divisors, independent
original-relation/kernel verification, certainty labels, finite work/time/storage
and checked resume. Use matched validated PyPy 3.11 runs, at least three seconds
of warmup and nine stable samples, fresh inputs and the roadmap's uncertainty
policy. There is no universal 10% promotion floor. This source inventory closes
no acceptance gate. Alpertron is GPL-3.0-or-later; copying code requires compatible
distribution terms. Independently implemented mathematical ideas still need
source attribution.

[upstream]: https://github.com/alpertron/calculators/tree/93c5c8189cb2f149c8dae898eb11997c5f2e7980
[partials]: https://github.com/alpertron/calculators/blob/93c5c8189cb2f149c8dae898eb11997c5f2e7980/siqs.c#L1752
[setup]: https://github.com/alpertron/calculators/blob/93c5c8189cb2f149c8dae898eb11997c5f2e7980/siqs.c#L2181
[sieving]: https://github.com/alpertron/calculators/blob/93c5c8189cb2f149c8dae898eb11997c5f2e7980/siqs.c#L156
[roots]: https://github.com/alpertron/calculators/blob/93c5c8189cb2f149c8dae898eb11997c5f2e7980/siqs.c#L3000
[scan]: https://github.com/alpertron/calculators/blob/93c5c8189cb2f149c8dae898eb11997c5f2e7980/siqs.c#L3326
[insert]: https://github.com/alpertron/calculators/blob/93c5c8189cb2f149c8dae898eb11997c5f2e7980/siqs.c#L2643
[la]: https://github.com/alpertron/calculators/blob/93c5c8189cb2f149c8dae898eb11997c5f2e7980/siqsLA.c#L57
[ecm-policy]: https://github.com/alpertron/calculators/blob/93c5c8189cb2f149c8dae898eb11997c5f2e7980/ecm.c#L64
[ecm-stage1]: https://github.com/alpertron/calculators/blob/93c5c8189cb2f149c8dae898eb11997c5f2e7980/ecm.c#L440
[ecm-stage2]: https://github.com/alpertron/calculators/blob/93c5c8189cb2f149c8dae898eb11997c5f2e7980/ecm.c#L585
[ecm-curve]: https://github.com/alpertron/calculators/blob/93c5c8189cb2f149c8dae898eb11997c5f2e7980/ecm.c#L779
[lehman]: https://github.com/alpertron/calculators/blob/93c5c8189cb2f149c8dae898eb11997c5f2e7980/factor.c#L868
[special]: https://github.com/alpertron/calculators/blob/93c5c8189cb2f149c8dae898eb11997c5f2e7980/factor.c#L293
[special-sums]: https://github.com/alpertron/calculators/blob/93c5c8189cb2f149c8dae898eb11997c5f2e7980/factor.c#L809
[java-special]: https://github.com/alpertron/calculators/blob/93c5c8189cb2f149c8dae898eb11997c5f2e7980/OldApplets/ecm.java#L2300
[java-proof]: https://github.com/alpertron/calculators/blob/93c5c8189cb2f149c8dae898eb11997c5f2e7980/OldApplets/ecm.java#L2378
[java-workers]: https://github.com/alpertron/calculators/blob/93c5c8189cb2f149c8dae898eb11997c5f2e7980/OldApplets/Siqs.java#L3307
[power]: https://github.com/alpertron/calculators/blob/93c5c8189cb2f149c8dae898eb11997c5f2e7980/bignbr.c#L1456
[trial]: https://github.com/alpertron/calculators/blob/93c5c8189cb2f149c8dae898eb11997c5f2e7980/factor.c#L1969
[carmichael]: https://github.com/alpertron/calculators/blob/93c5c8189cb2f149c8dae898eb11997c5f2e7980/factor.c#L1768
[bpsw]: https://github.com/alpertron/calculators/blob/93c5c8189cb2f149c8dae898eb11997c5f2e7980/bignbr.c#L2265
[inverse]: https://github.com/alpertron/calculators/blob/93c5c8189cb2f149c8dae898eb11997c5f2e7980/modmult.c#L435
[inverse-small]: https://github.com/alpertron/calculators/blob/93c5c8189cb2f149c8dae898eb11997c5f2e7980/modmult.c#L346
[native]: https://github.com/alpertron/calculators/blob/93c5c8189cb2f149c8dae898eb11997c5f2e7980/MontgomeryMult.c#L52
[fft]: https://github.com/alpertron/calculators/blob/93c5c8189cb2f149c8dae898eb11997c5f2e7980/fft.c#L390
[ntt]: https://github.com/alpertron/calculators/blob/93c5c8189cb2f149c8dae898eb11997c5f2e7980/ntt.c#L498
[poly]: https://github.com/alpertron/calculators/blob/93c5c8189cb2f149c8dae898eb11997c5f2e7980/multpoly.c#L626
[known]: https://github.com/alpertron/calculators/blob/93c5c8189cb2f149c8dae898eb11997c5f2e7980/factor.c#L2095
[batch]: https://github.com/alpertron/calculators/blob/93c5c8189cb2f149c8dae898eb11997c5f2e7980/batch.c
[front]: https://github.com/alpertron/calculators/blob/93c5c8189cb2f149c8dae898eb11997c5f2e7980/ecmfront.js#L181
[gmp]: https://github.com/sethtroisi/gmp-ecm/tree/8ea5e214fdf2f0ddf9415141b8dc039ed6f5874e
[yafu]: https://github.com/bbuhrow/yafu/tree/8110dfbd8c6f9486d93b1a02de6eb7b180e55a80
[msieve-release]: https://sourceforge.net/projects/msieve/files/msieve/Msieve%20v1.53/msieve153_src.tar.gz/download
[msieve-la]: https://github.com/radii/msieve/blob/c8727d91305bdbe0972d160ef0ce61dd02ce9193/common/lanczos/lanczos.c
[cado]: https://github.com/cado-nfs/cado-nfs/tree/692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b/linalg/bwc
