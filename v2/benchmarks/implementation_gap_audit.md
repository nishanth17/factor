# Multi-source factoring implementation gaps and next priorities

Audit date: 10 October 2026. Initial mainline inspected: `cde025d`;
C3 integration/status reconciled against `d273cb8` before delivery.
This extends the [Alpertron inventory](alpertron_gap_analysis.md),
[QS source reconciliation](qs/qs_gnfs_research.md),
[QS parameter ledger](qs_parameter_audit.md) and
[cross-engine ledger](parameter_audit.md). It changes no algorithm, default
or historical acceptance decision. MG identifiers below identify research
candidates, not additional implementation tasks.

## Decision and limits

**Continue C10, with A1a implementing its selected collector kernels.**
Separately, A1b can review promotion of existing scoring options on stable
mainline; it is not bundled into or blocked on C10. Identify costly collection
functions on the settled C9 control, then choose implementation changes
using that evidence. Polynomial normalization, actual-product A
selection and threshold adaptation deserve explicit investigation. They
were less concretely developed in the earlier inventory than generic
base/interval tuning and the Alpertron kernels.

**Scope reconciliation with the supplied running C10 prompt:** it already
explicitly includes MG01–MG03 and traces the wider collection pipeline
covered by MG04–MG10. These findings support its existing research scope;
they are not additional tasks or requirements to implement every candidate.
C10 selects and evaluates a coherent tranche, with A1a retaining selected
kernel implementation. Reserve C10a for an explicitly deferred or uncovered
collection delta identified against C10's delivered report, rather than
duplicating its current investigation. A1b remains separate.

This is not an exhaustive proof that every useful idea has been found.
Recursive trees were inventoried for six repositories; 445 selected Git
files were downloaded and checked against their Git blob hashes. Another
21 QS/Lanczos/build/documentation files were extracted from the official
msieve 1.53 archive. Targeted arithmetic paths, callers, configuration
branches and v2 counterparts were inspected; downloading a file does not
mean every line was audited. The earlier Alpertron inventory covers its
complete 307-blob tree separately. No upstream program was executed and
no external implementation was copied into the factoring runtime.

There are no new performance measurements in this report. Rankings are
engineering judgments from existing evidence, affected stages, applicability
and implementation cost. A missing mechanism is not automatically a gain.
Native operation counts, SIMD performance and upstream presets do not
predict PyPy performance. GNFS field selection/lattice sieving and every
GMP integer kernel were not freshly audited in this pass.

## Controls and evidence used for the ranking

Mainline retains the accepted C1 result: 31.73% lower complete-run cost in
its fixed 40-digit confirmation cohort. Its unresolved historical 50-digit
120-second comparisons did not establish a universal feasibility limit.
See [C1's result and limitations](qs/c1/c1_implementation_results.md).

The separate C9 branch was inspected through `e065a69`, including
`v2/benchmarks/qs/c9/c9_growth_results.md`. C9.2's implementation and training
are complete; fresh confirmation is still in progress at this checkpoint.
Its full-base aligned growing-DLP arm was selected on revealed inputs.
The two instrumented calls took 47.475/54.697 seconds versus SLP's
72.088/84.721 seconds. These are selection diagnostics, not accepted speed
estimates or a population forecast. The half-base combination was worse
in aggregate; do not assume B1/C1's old half-base choice remains optimal.

In one growing unaligned call, collection took 41.40 of 54.825 seconds,
about 75.5%. The aligned call spent 34.30 of 47.475 seconds in collection.
Profiles identify sieve/root work as the next actionable investigation.
Splitter failures and cycle-length rejections were zero in the initial grid;
the resulting matrices had roughly 1,500 rows. This supports investing in
collection before replacing successful splitters or building large solvers.
It does not establish that every input has the same stage shares or that
matrix capacity will never matter.

For scale only: if a pure kernel change halved that 41.40-second collection
stage and changed nothing else, total cost would fall about 37.8%. That is
an Amdahl calculation, not an expected improvement. Yield changes can alter
all stages and require complete-factor comparisons.

C3's source-preset study and allocation/resume implementation were also
read in its separate worktree and reconciled after integration at `d273cb8`.
Its delivered implementation tranche is complete, with user-directed
compact64/escalating140 defaults in the 30/40-digit bands. That promotion
is distinct from completed fresh confirmation. C3a must consume this final
coverage rather than repeat the same source-derived allocation study.

## Sources and provenance

| Implementation | Pinned identity and inspected areas | License evidence / qualification |
| --- | --- | --- |
| [YAFU][yafu] | `8110dfbd8c6f9486d93b1a02de6eb7b180e55a80`; SIQS setup, polynomial generation/roots, small/medium/large sieve and recovery, filtering, batch factoring, micro-ECM, dispatcher, bundled Lanczos; `include/qs_impl.h` and build files | Inspected SIQS files have public-domain notices and upstream acknowledgements. Bundled libraries have separate notices. |
| [Official msieve 1.53][msieve-release] | `msieve153_src.tar.gz`, SHA-256 `c5fcbaaff266a43aa8bca55239d5b087d3e3f138d1a95d75b776c04ce4d93bb4`; `mpqs/`, `common/lanczos/`, `Readme.qs`, `Makefile` | Inspected source files have public-domain notices. This archive is distinguished from the historic radii and MersenneForum mirrors in earlier ledgers. |
| [Yamaquasi][yama] | `3f95f43682ed15d8c1ed206a9a702dd655d7c8ad`; SIQS/MPQS, sieve, relation handling, fbase, rho/ECM, sparse/dense algebra and dispatcher | BSD-3-Clause. |
| [FLINT][flint] | `a4c9750d0d3d67bb01cf6d18c187591b313451c3`; QS setup/collection/polynomials/storage/filter/Lanczos/extraction and `fmpz_factor` dispatch/pretests | Inspected files specify LGPL-3.0-or-later; COPYING and COPYING.LESSER retained. |
| [GMP-ECM mirror][gmp] | `8ea5e214fdf2f0ddf9415141b8dc039ed6f5874e`; ECM/PM1/PP1, batch curves, classical/fast stage two, polynomial/remainder and arithmetic helpers | Inspected library files specify LGPL-3.0-or-later. This is a pinned mirror, not a verified latest official tip. The [author's download page][gmp-official] points to Inria GitLab, whose web access returned an anti-bot page. |
| [CADO-NFS][cado] | `692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b`; cofactor strategies, batch cofactors, ECM/mixed coordinates, filter/merge and block-Wiedemann documentation | Root COPYING is LGPL-2.1; check selected component notices before any adaptation. The full Block Wiedemann implementation was not line-by-line audited here. |
| [JavaMath][java] | `088d01fa97e7d0412c6cbbb8f463fa7e5f78ab97`; SIQS A generation, polynomial generation, sieve/thresholds, trial division, partial collectors/solvers and matrix solvers | GPL-3.0-or-later in source notices. Variant classes are not all reached by the same constructor. |
| [Alpertron][alpertron] | `93c5c8189cb2f149c8dae898eb11997c5f2e7980`; prior complete repository inventory and AP01–AP40 | Preserve the earlier active/historical/build distinctions and component notices. |

The SIQS primary references remain [Alford–Pomerance][siqs-paper] and
[Carrier–Wagstaff][qs-book]. For specialized P±1 continuation, the source
names Montgomery–Kruppa; their [ANTS-VIII presentation][fast-stage-paper]
explains reciprocal Laurent polynomials, structured construction and
geometric-progression evaluation. CADO's batch implementation names
Bernstein's Algorithm 7.1 and explicitly separates smoothness screening
from subsequent factoring. Existing author notes and technical discussions
are discovery aids; implementation claims below use source and callers.

## Most useful remaining collection work

| ID / owner | Concrete difference, source and v2 comparison | Priority / implementation boundary |
| --- | --- | --- |
| MG01 — C10 | **Discriminant-normalized polynomials.** Yamaquasi `src/siqs.rs:polytype`, `Poly::eval`, `_finish_polynomial` use `G=Ax²+Bx+C`, `B²−4AC=kN` when `kN≡1 mod4`. YAFU `tdiv_small.c` has the `knmod8==1` recovery branch; JavaMath `SIQS.java` selects `d=2` for that case. v2 `qs/polynomial.py:Polynomial` fixes `F=Ax²+2Bx+C`, and `a_target` fixes its associated target. | **High-value research; medium economic confidence.** Smaller typical norms can increase admissible smooth values without enlarging the residual graph. The roadmap already mentions a restricted 2-adic variant; the new contribution is a concrete cross-source implementation hypothesis. This is not a wholly absent roadmap idea. Verify factor 4, A support, p=2, exceptional roots, multiplier coupling, Gray shifts, recovery and checkpoint identity. |
| MG02 — C10 | **Select by actual A-product quality.** Yamaquasi `siqs.rs:select_a` exhaustively ranks products in a small eligible regime, otherwise samples and narrows to suitable products. v2 `assignment_stream.py` ranks individual primes by `abs(p**factor_count−target)`, traverses distinct combinations, and optionally chooses a nearest flyer outside the core pool. | **Useful incremental challenger.** Individual-prime closeness does not rank whole products. v2 already has unique combinations and a flyer, so neither generic random A selection nor a nearest last factor is a new gap. Compare finite product search, distribution of actual norms, family diversity/overlap and setup amortization. Do not copy an upstream retry loop without finite refusal. |
| MG03 — C10 | **Per-job threshold adaptation.** YAFU `SIQS.c:siqs_sync` samples `tf_small_cutoff` at its original value and ±5, then selects from relation throughput when optimization is enabled. Overrides and resumed jobs disable it. v2 has configurable fixed scoring/cutoff policies and offline studies, not an equivalent online selection policy. | **Plausible, after cheaper fixed bundles.** Distinguish threshold units from v2's numeric small-prime cutoff. Use useful yield and total cost, charge the learning prefix, and specify persisted observations/decisions on resume. YAFU's raw-relations objective and 2%/5% switching biases are not v2 acceptance rules. |
| MG04 — A1a; C10 diagnosis | **Precomputed Gray root shifts and family setup reuse.** Yamaquasi stores `deltas_mod_p` and pool-prime inverses; FLINT `compute_poly_data.c` builds `A_inv2B` and updates roots by addition/subtraction. v2 `families.py:PolynomialFamily` already caches A inverses and updates roots, but multiplies each centered-B delta by an inverse on every step. | **Strong match to observed root cost.** AP06 already owns the core mechanism. Measure table/setup/storage costs and handle centered-B wrap, A primes and 2 independently. Pairwise pool inverses are an additional setup option, not a replacement for the accepted inverse backend. |
| MG05 — A1a; C10 diagnosis | **Avoid repeated derivation work on trusted internal switches.** v2 `sieve_collector.py:set_polynomial` revalidates all roots and refactors supported A for each switch. It correctly caches A support for candidate recovery; A support itself is not missing. A certified family transition might reuse more metadata. | **Derived candidate, profile first.** This is an inference from v2's callers, not a copied upstream algorithm. Preserve full validation of public/untrusted construction and checkpoints, and original relation verification. Do not simply delete checks or change charged work without documenting the revised contract. |
| MG06 — A1a; C10 diagnosis | **Cheaper staged candidate screening.** Yamaquasi `sieve.rs:smooths` combines bulk masks, tiny-prime corrections and an integer distance-to-root correction before factor recovery. Java's Factory02 distinguishes sieve-hit, trial-test and final smooth-residual domains. v2 already has a safe block minimum and exact per-candidate norm refinement. | **Potential saving only in the missing stage.** Investigate an inexpensive conservative piecewise/root-distance bound before bigint norm evaluation or division; otherwise declare and measure deliberate loss. Do not describe staged scoring itself as missing or transfer upstream heuristic thresholds as proofs. AP08 owns bulk-mask implementation. |
| MG07 — A1a/B7 | **Packed multi-level hit tables.** Yamaquasi `sieve.rs:SieveTable`/`SieveTableLarge` group hit records by prime bit length, with different bucket widths and compressed prime offsets. v2 already has sparse/dense/bucket marking and hit-guided resieving, but different Python object/bitset layouts. | **Conditional on hit-storage/traversal cost.** Native packing/locality may not win in PyPy. The small-table overflow buffer is finite and excess hits can be omitted; a transfer must preserve exact coverage with checked spill/refusal or explicitly measure the loss. Generic buckets are not a new mechanism. |
| MG08 — A1a; C10 diagnosis | **Periodic presieving, paired traversal, density loops and scan fusion.** Alpertron AP05/AP07–AP09, YAFU specialized sieve/recovery files and Java sieve variants provide concrete alternatives to per-position/per-prime Python traversal. | **Already planned; still important.** Choose the costly function first. A periodic multi-prime initializer differs from R2's losing exact tiny-prime singleton arm. Native paired lanes and SIMD are design references, not guaranteed Python representations. |
| MG09 — C10 | **Multiplier, base/interval and family quotas together.** Official msieve `mpqs.c:choose_multiplier` and size-interpolated parameters, YAFU setup, Java Factory02 and Yamaquasi presets couple norm, base and candidate domains. v2 has controls and measured bundles, but no demonstrated optimum across changed DLP regimes or larger bands. | **High-value calibration, not a new generic algorithm.** Normalize count/value and interval/block units; evaluate wider/even multipliers with exact exceptional roots. Include A-factor count, family reuse and optional MG01. C9's half-base training loss is a concrete reason to revisit combinations. |
| MG10 — C10 with A8/C4 | **Collection stopping and postprocessing cadence.** Upstreams request relation excess and increase it after unsuccessful postprocessing. v2 has row/excess/cadence controls and prepared caches; C9 notes forced polynomial-window solves can repeat preparation. | **Worth measuring; no solve-bottleneck claim.** Count repeated preparation, useful nullity and trivial dependencies before redesign. Original relations remain available for later cycles even if a current filter discards singletons. Solver/storage changes retain their existing owners. |

### Existing scoring option: separate A1b promotion review

R2 already implemented fixed-point scoring and capped per-polynomial power
plans. Its [historical record](docs/studies.md#p38-r2-bounded-collector-evaluation--4-october-2026)
reports approximately 5.7% lower complete-cohort time on both medium and
30-digit held-out cases, with conditional intervals 0.8–7.1% and 1.4–8.2%.
The nominal 20-digit interval crosses zero. The frozen study retained
defaults because training improvements missed its then-required 10% floor.
This is a preserved historical verdict, not proof the option is unhelpful.

C9's frozen arms use list scores, sparse marking, streamed `powers`, zero
threshold extra and zero power-plan bytes. A1b separately reviews promotion
of existing fixed-score and plan options on stable current mainline. It
does not require C10's research or C9's completion. A claim about C9's new
growing-DLP configuration additionally needs that settled control. Reuse
applicable retained evidence; a new current-default performance claim needs
current matched complete-factor evidence. Separate scoring from plan cost
where useful; the old report's combined result does not attribute all of
its gain to either component. C10 consumes the decision without duplicating
the comparison. No new scoring algorithm is assigned to A1b.

The old conservative root-weight, cutoff-5, exact tiny-prime singleton,
batch smooth-part and grouped-charge prototypes increased matched 30-digit
training cost by about 93.5%, 15.7%, 17.6%, 21.4% and 9.8%. Preserve those
losses. They do not reject a substantially different presieve or changed
larger workload automatically, but they do rule out presenting the same
prototype as an untested free improvement. Apply the revised promotion
policy without a universal 10% improvement floor.

## Storage, residuals and linear algebra

| ID / owner | Remaining difference and evidence | Placement |
| --- | --- | --- |
| MG11 — C9a/B7 | **Compact relation arenas and deferred replay.** FLINT `large_prime_variant.c:qsieve_add_relation` stores relation descriptors plus offsets into a separate growing factor arena. YAFU's enabled `SPARSE_STORE` stores selected factors and reconstructs cheap factors/A contributions during filtering (`tdiv_small.c`, `tdiv.c`, `filter.c`). Alpertron AP02 instead reconstructs from its seed/U payload. v2 keeps full atom/exponent objects. | **Potentially high when memory or replay cost limits useful yield.** C9 growth preserves history but does not implement these representations. Compare arena, compressed immutable records and recomputation as distinct designs. Charge pinned data, factor reconstruction, scratch, serialization and corruption detection. A lower conservative reservation alone is not compact storage. |
| MG12 — C9a | **Size-dependent residual splitting portfolios.** Yamaquasi uses rho then ECM fallback; Java's 2LP path selects Hart, tiny ECM, rho or nested SIQS by residual size; YAFU uses native micro-ECM. CADO facul chooses strategies by cofactor-size tuples and resumes later methods after partial splits. v2's certified square/rho splitter has fixed attempts and per-attempt work. | **Defer replacement on present evidence.** C9's initial grid has no failed splits. Broader domains or a measured splitting/certification share could change the ranking. Preserve proven endpoints; native probable-prime filters and fixed-width speed do not establish an equivalent certainty contract. |
| MG13 — C9a/F7 | **Batch cofactor smoothness screening.** CADO `sieve/ecm/batch.cpp:smoothness_test` constructs product/remainder trees and repeatedly strips GCDs; `find_smooth_groups` amortizes one tree across survivor lists. YAFU `factor/batch_factor.c` also implements batch cofactors. This is distinct from batching the factor-base division of already selected sieve candidates. | **Conditional future delta.** F7 already holds the losing candidate-recovery revisit. A residual-prime-product batch needs a separate workload, batch-latency/storage and exact-valuation hypothesis; do not label it the same rejected experiment or infer a current win. C9a owns intake semantics; F7 supplies the tree experiment. |
| MG14 — A8/B6/C4/D1/D2 | **Matrix representations and algorithms.** Official msieve has packed/cache-blocked sparse products and Lanczos with removed-row correction; Yamaquasi/Alpertron have sparse/dense kernels; CADO describes the Block Wiedemann sequence/generator/reconstruction pipeline. v2 uses Python-int elimination and exact filtering/lifting. | **Real gaps, already planned, conditional on genuine matrix cost/capacity.** A8 can prepare representative matrices now. Four Russians/PLE, structured merging and compact provenance can precede a full sparse solver when they fit. A verified memory refusal can justify Lanczos/Wiedemann without dominant solve CPU; current small matrices alone do not. |
| MG15 — A8/C4 | **Joint modular exponent-product extraction.** FLINT `square_root.c` forms bit-plane products and evaluates exponents by Horner, with word-sized product accumulation. v2 `extraction.py:extract_dependency` calls modular power for each base prime's half exponent. | **Low until extraction is costly; cheap isolated oracle.** Exact joint exponentiation works for nonunits as well. Apply existing residual square corrections after the Horner product, preserve full parity/identity checks and both GCD signs. Native word flush thresholds do not transfer unchanged. |

Opaque equal-composite residual matching, SLP anchor reuse, graph growth and
all-component cycle identities were already audited under C1/C9/C1a. They
are not newly discovered gaps. C9's reusable graph-SLP combination lost its
revealed aggregate comparison; a cheaper anchor-specific arm needs a changed
cost hypothesis. A path heuristic that ignores disconnected cycles or a
collector that drops all-even rows is not an unconditional improvement to
v2's stronger exact dependency contract.

## ECM, P±1, preprocessing and arithmetic

| ID / owner | Remaining difference and source | Placement |
| --- | --- | --- |
| MG16 — F8 | **Specialized fast P−1/P+1 stage two.** GMP-ECM `pm1.c`/`pp1.c` dispatch to `pm1fs2`/`pp1fs2` in `pm1fs2.c`, selecting convolution parameters and NTT/non-NTT paths. Montgomery–Kruppa exploit reciprocal Laurent polynomials and geometric-progression evaluation. v2 PM1 uses streamed/gap/wheel continuations; Williams P+1 remains A5. | **A separate algorithmic gap.** F3's generic polynomial ECM/Brent–Suyama wording does not identify this specialized method adequately. Add F8, first gated by actual PM1 stage-two cost and feasible exact polynomial kernels; P+1 coverage additionally needs A5. It does not need C2's ECM geometry decision. No native high-B2 speed claim transfers directly. |
| MG17 — F2 | **Batch-specific curve/formula packages.** GMP-ECM `batch.c` and `ecm.c` expose small-coordinate and small-parameter stage-one packages plus a reusable exponent. v2 already batches/compiles prime powers, but not those curve/formula packages. | **Conditional curve experiment already within F2.** Sharing an exponent is not the same as choosing a new curve/point/formula. Prove nonsingularity, chosen-point behavior and nonunit recovery; measure whole-engine success including setup and compatible stage two. |
| MG18 — F2 | **Mixed Edwards/Montgomery and torsion families.** CADO `ecm.cpp:bytecode_mishmash_interpret_mixed_repr` is reached by the `FULLMONTYTWED` stage-one branch; parameterization dispatch includes Brent–Suyama, Montgomery12/16 and Z6 constructions. | **Known roadmap scope, not a new discovery.** F2 already explicitly names MISHMASH. Existing PRAC/chain integration does not provide coordinate maps or those formulas. Point order and curve order must be distinguished. |
| MG19 — E3/C3a/C11 | **Factoring pretest/dispatch packages.** FLINT `factor_no_trial.c` combines proof/power checks, smooth-factor work, native two-limb rho/ECM/SQUFOF and QS fallback. `factor_smooth.c` uses size-derived trial/ECM effort. Alpertron adds the previously inventoried multiplier Fermat/Lehman and structured forms. | **Mostly existing capability or policy scope.** v2 already has exact power screens, trial division, rho, PM1, ECM and validated fallback. SQUFOF/Hart are missing methods, but there is no present evidence that adding them beats tuning existing PyPy pretests. E3/C3a own uncovered cases; source digit thresholds are not optimal defaults. |
| MG20 — A10a/B14 | **Classification/proof deltas.** MR-residue factor recovery, optional BPSW rejection, grouped trial remainders and prime certificates remain useful candidates from the earlier comparison. | **Known scope.** A10 already settled deterministic ranges; do not reopen it as a generic primality rewrite. B14 adds independently verified proof capability, not a demonstrated broad factoring speedup. Keep PRP labels distinct from proof. |
| MG21 — C11/A2a/C5 | **Native reciprocal division, limb reducers, fast multiplication and inversion.** Official msieve precomputes reciprocal division constants; YAFU/FLINT use native word/GMP operations. v2 already dispatches native math or persistent GMP arithmetic through a coarse backend boundary. | **Conditional kernel references, not a reason to build Python limbs.** Accepted GCD/root/backend losses and reducer losses remain the control. Revisit quotient/inverse work only for a newly costly fused path; include conversion and nonunit costs. Core GMP GCD algorithms were not independently re-audited here. |
| MG22 — A7a/F4/A11a | **Parallel campaigns, disk relation exchange and product adapters.** msieve documents external relation-file merging; YAFU coordinates thread-local polynomial work and savefiles; Alpertron has browser/input/progress features. | **Known optional capability or worker scope.** v2 already has bounded worker contracts and checked resume. Full disk spill/distributed campaigns require independent storage/identity/accounting work. Native threads and batch UI are not current serial speed evidence. |

YAFU's `USE_SS_SEARCH` and `USE_POLY_BUCKET_SS` definitions are commented
out at this pin, unlike enabled `SPARSE_STORE` and `USE_BATCH_FACTOR`.
Subset-sum polynomial bucket code is therefore an optional source path,
not evidence that the ordinary build gets that benefit. Keep C8/SSS
comparisons separate from ordinary Gray-shift tables.

## Order of work and what could change it

1. **Now, independently: A1b promotion review and C10 research/oracles.**
   A1b uses stable mainline and existing scoring/plan options. C10 pins its
   source control and assesses the base/interval/A/family choices that C9
   changed the economics of, using settled C9 for interacting comparisons.
   Neither task waits for broad completion of the other. Serialize heavy
   timing. Do not wait for every possible source gap to be catalogued before
   implementing an evidenced C10 change.
2. **Next implementation: C10's selected policy/algorithm or A1a's selected
   collector kernel.** Root-shift/setup reuse is strongly aligned with the
   observed cost; periodic presieve and scan/density kernels are other
   candidates. MG01/MG02 target useful yield directly but need more arithmetic
   and recovery integration. A function-level profile and source/oracle checks
   should choose the first tranche, not an invented speedup estimate.
3. **Conditional next: C9a compact storage or C10 stopping/cadence.** Advance
   whichever has a demonstrated limiting cost. Growth has reduced the urgency
   of simply increasing quotas. It has not solved compact ownership, replay
   or durable capacity at arbitrarily larger inputs.
4. **A8 preparation can overlap without heavy timings.** Preserve real
   matrices/provenance and establish cost/capacity. Select B6/C4/B7/D1/D2
   from that evidence; no current result makes a large solver the immediate
   highest-return implementation.
5. **Other engine work remains separate.** A5 is an independent complementary
   capability. F2/F3/F8 may have large benefits in their applicable regimes,
   but need their workload/backend gates and do not presently outrank
   collection for the measured SIQS problem. C3a first consumes C3's final
   scope. C11 calibration remains continuing rather than an indiscriminate
   sweep of every constant.

This ordering can change for a concrete reason: fresh C9 confirmation loses,
collection stops dominating, a stable fixed-score comparison wins, memory
refusals appear, or genuine matrix/continuation cost becomes dominant.
It should not change merely because another implementation has a prestigious
algorithm name. No universal 10% floor or 120/360-second ceiling is introduced.

## Validation and retained evidence

Independent checks ran with PyPy implementing Python 3.11.15:
16,218 normalized-polynomial integer identities, 7,632 modular-root sets
including p=2 and A-prime degeneracies, and 2,000 joint exponent products
including nonunits, zero exponents and square corrections. They validate
the small arithmetic examples, not an implemented family generator,
checkpoint contract or performance claim. A simple example is
`G(x)=13x²+5x−1`, for which `(26x+5)²−77=52G(x)` and every G(x) is odd.

Downloaded sources, trees, manifests, the archive and check script are local
ignored evidence under `v2/audit/multisource-gap-20261010/`. The manifest
records file sizes, Git blob IDs and SHA-256 hashes; the msieve archive
identity is above. Reproduce the arithmetic check with:

```sh
pypy3 v2/audit/multisource-gap-20261010/check_identities.py
```

Only documentation is changed. No implementation tests, benchmark timings,
production promotion, merge or push are claimed. Public roadmap checkboxes
and completed-task records stay unchanged.

[yafu]: https://github.com/bbuhrow/yafu/tree/8110dfbd8c6f9486d93b1a02de6eb7b180e55a80
[msieve-release]: https://sourceforge.net/projects/msieve/files/msieve/Msieve%20v1.53/
[yama]: https://github.com/remyoudompheng/yamaquasi/tree/3f95f43682ed15d8c1ed206a9a702dd655d7c8ad
[flint]: https://github.com/flintlib/flint/tree/a4c9750d0d3d67bb01cf6d18c187591b313451c3
[gmp]: https://github.com/sethtroisi/gmp-ecm/tree/8ea5e214fdf2f0ddf9415141b8dc039ed6f5874e
[gmp-official]: https://members.loria.fr/pzimmermann/records/ecm/devel.html
[cado]: https://github.com/cado-nfs/cado-nfs/tree/692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b
[java]: https://github.com/TilmanNeumann/java-math-library/tree/088d01fa97e7d0412c6cbbb8f463fa7e5f78ab97
[alpertron]: https://github.com/alpertron/calculators/tree/93c5c8189cb2f149c8dae898eb11997c5f2e7980
[siqs-paper]: https://math.dartmouth.edu/~carlp/implementing.pdf
[qs-book]: https://homes.cerias.purdue.edu/~ssw/qs4.pdf
[fast-stage-paper]: https://antsmath.org/ANTSVIII/files/kruppa.pdf
