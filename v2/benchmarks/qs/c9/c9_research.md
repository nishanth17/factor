# C9 source-derived useful relation yield

Research freeze: 10 October 2026. No upstream implementation was copied or
executed. Source archives and raw evidence are ignored; versioned
`inputs/controls/c9_source_pins.json` records revisions and inspected-file
hashes. This extends, rather than replaces, [C1 research](../c1/c1_research.md)
and the [parameter ledger](../../qs_parameter_audit.md).

## Arithmetic and economic model

A fully divided polynomial value supplies an exact residual. In the prime
endpoint contract, an SLP relation is an edge (1,p), a DLP relation an edge
(p,q). Every even-incidence set cancels all endpoint parities. A loop p=p
already supplies a square; parallel edges supply a two-edge cycle; disconnected
cycles are just as valid as cycles through 1. The independent cycle-space
dimension is E−V+C. This is only the first kernel: factor-base parity must
also cancel, the original atoms and repeated square factors must reconstruct,
and both GCD signs must be tested for a proper divisor.

[Lenstra–Manasse, 1994, §§1–3](https://www.researchgate.net/publication/266349860_Factoring_with_Two_Large_Primes)
describes that graph construction, k−1 independent combinations from k equal
partials, nonuniform endpoint occurrence and initially delayed combined yield.
Their sieve discussion charges increased false reports and splitting. A larger
admissible product increases opportunities and processing cost together.
Their hardware-specific crossover is not a Python digit threshold.

[Boender–te Riele, 1996, §§5–7](https://ir.cwi.nl/pub/1367/1367D.pdf)
separates the factor-base bound, one-prime limit and double-prime region.
Its SGI/Cray evidence shows that memory, reporting thresholds and implementation
change the crossover. Its DLP prediction applies to similar sizes and fixed
parameters. Section 6 explicitly rejects prime residuals beyond the common
endpoint bound because those partials contributed little in its experiments. This supports
measuring long prefixes, rather than extrapolating a 60-second empty matrix
or copying an 80-digit threshold.

[Cavallar, 2000, §2.1–2.2](https://ir.cwi.nl/pub/2078/2078D.pdf)
explains recursive singleton removal and the size/weight tradeoff of merging.
Removing a singleton is exact for the **current batch**, but not a proof its
partial can never become useful after more collection. Her NFS ideal and
character requirements do not transfer wholesale to QS. Filtering redesign
remains with the matrix owners.

[Gower–Wagstaff, Square Form Factorization, §§3–5](https://homes.cerias.purdue.edu/~ssw/squfof.pdf)
gives SQUFOF's arithmetic, queue and multiplier obligations. It is a plausible
bounded residual splitter, not evidence to replace successful v2 rho. Batch
smooth-part extraction would additionally require finite batch storage and
amortization; higher-arity relation admission remains F6.

A useful statistical approximation is that an endpoint seen with a small
Poisson mean λ contributes about λ records but only λ²/2 equal-endpoint
pairs. Keeping many very rare endpoints therefore has poor immediate return
per stored atom. DLP changes the calculation because two distinct endpoints
can become useful through a path; connected-component development and its
retained history matter, not just pair collisions. This approximation ranks
hypotheses; it does not prove an optimal cutoff or justify deleting degree-one
edges. The fixed-prefix replay measures the actual lost matches and gained
cycles instead of assuming uniform prime frequencies.

## Implementation comparison: follow the actual intake

B means the largest factor-base prime below, unless explicitly stated. v2's
`base_bound` is an exclusive numeric limit; upstream factor-base cardinalities
are not interchangeable with it. All native byte/log thresholds need conversion.

| Implementation / source path | Candidate, residual and effort policy | Matching, retention and useful dependencies | Transfer and cost |
| --- | --- | --- | --- |
| [**YAFU** `8110dfb`](https://github.com/bbuhrow/yafu/blob/8110dfbd8c6f9486d93b1a02de6eb7b180e55a80/factor/qs/tdiv.c); `factor/qs/SIQS.c`, `siqs_aux.c`, `tdiv_scan*.c`, `tdiv.c`, `filter.c` | SIQS setup chooses bit/digit-dependent base, reporting thresholds and LP multiplier. Default DLP starts at 70 digits; forced paths and higher arity are separate. SLP intake is residual < L. DLP requires >B², <L^e and <64 bits; e is 1.75/1.8/1.9/1.95 across bit tiers <270/<290/<330/otherwise. Base-2 PRP rejection precedes `getfactor_uecm`; unsuccessful splitting drops the candidate. Both endpoints must be <L. This path validates divisibility, not v2-style endpoint certificates. | Buffered/savefile relation records, growable vertex hash/cycle tables, duplicate/singleton passes and cycle construction. Table doubling is vertex storage, not a fixed retained-edge count. Relation-list growth and disk records are separate. Collection starts from base cardinality +64, with progress/extra-relation updates in callers. | Common SLP/DLP L motivates endpoint alignment. Do not copy PRP rejection, unsigned limits or native micro-ECM speed claims. Growing disk-backed stores would require a new finite storage contract. |
| [**msieve** `03e0ab9`](https://github.com/MersenneForum/msieve/blob/03e0ab9e13c4862e5bb2c51a4a39718aa2500e19/mpqs/sieve.c), cross-checked with 1.53; `mpqs/mpqs.c`, `sieve.c`, `relation.c`, `sqrt.c` | Size-interpolated base/cardinality/interval/LP multiplier; report threshold derived from residual allowance. DLP starts at 282 bits. SLP first, then product about L^1.8, base-2 rejection and SQUFOF; two endpoints <L. Full squares handled. Failed split discarded. Custom multiprecision and native word arithmetic, not PyPy int. | Text savefile records persist independently of the cycle table (initial allocation 10,000, growing). Filter removes singleton/duplicate relations, roots every component, constructs and sorts cycles. Traceback arrays impose 100 edges per side. Default target is fb_size+96 (64 extra solver relations). Verification reconstructs saved relations before square root. | Common L and excluding rare large-large pairs are useful policy hypotheses. v2 already has all-component cycles and exact verification, and must retain disconnected cases rather than adopt path-loss heuristics silently. Native slots/records are not 4-KiB Python forest reservations. |
| [**Yamaquasi** `3f95f43`](https://github.com/remyoudompheng/yamaquasi/blob/3f95f43682ed15d8c1ed206a9a702dd655d7c8ad/src/fbase.rs); `src/siqs.rs`, `fbase.rs`, `pollard_rho.rs`, `relations.rs` | Default DLP above 256 bits with half-size base preset. L=min(B*large_factor,2³²−1); multiplier near 166 bits is floor(166*66/100)=109. Nominal DLP score allowance B²*floor(bits/2) below 256 bits. **Caller distinction:** `cofactor` does not enforce that nominal product: it admits lucky reports up to L², splitting above B² with rho then ECM fallback after a composite screen. Endpoints individually ≤L. | Packed partial hash map and ordered double maps; greedy propagation from SLP-connected component, shortest available representatives, equal-endpoint squares. Duplicate endpoint-pair map keys overwrite. All-even combined-factor rows can be discarded as a heuristic. No finite FIFO cap. `gap` requests covered columns+48 cycles; final sparse solve/extraction follows. | A score target is not an admission bound. C9 keeps conservative v2 coverage and explicit product admission. Common endpoint ceiling and retention are candidates; its ignored disconnected cycles and all-even-row discard are incompatible with supported v2 arithmetic cases. |
| [**JavaMath** `088d01f`](https://github.com/TilmanNeumann/java-math-library/blob/088d01fa97e7d0412c6cbbb8f463fa7e5f78ab97/src/main/java/de/tilman_neumann/jml/factor/siqs/tdiv/TDiv_QS_2LP.java); `SIQS.java`, `SieveParamsFactory02`, `TDiv_QS_2LP`, `CongruenceCollector01`, `PartialSolver01` | Factory separates sieve-hit, trial-test and smooth residual bounds using floating powers of N; callers predivide small primes. After exact division: residual 1 is smooth, residual <B² or passing PRP is SLP, normal endpoints limited to 31 bits. Composite splitting routes Hart <46 bits, tiny ECM <63, rho <64, otherwise nested SIQS; too-large/failing outcomes rejected. Repeated factors become `Smooth1LargeSquare`. | SIQS constructor actually uses Collector01, not every available CycleFinder. Factor→partial lists gather connected partial sets for a partial solver; only independent partials retained after finding a smooth combination. Dynamic Java containers. The serial CombinedFactorAlgorithm big-argument caller supplies 10 extras and block Lanczos (small arguments use TDiv_QS_Small/Gaussian elimination); after failure another 10 are requested. | Confirms splitting effort should depend on residual size, and surplus needs extraction feedback. Native/JVM backend and probable-prime semantics differ. No evidence for new v2 splitter while current split failures are zero. Related-set solve can become expensive; do not substitute it for the bounded forest. |
| [**FLINT** `a4c9750`](https://github.com/flintlib/flint/blob/a4c9750d0d3d67bb01cf6d18c187591b313451c3/src/qsieve/large_prime_variant.c); `qsieve/collect_relations.c`, `factor.c`, `large_prime_variant.c`, `linalg.c` | Small-prime refinement precedes full division. Full residual 1 accepted; small residual (30-bit bound) <60B and coprime to multiplier becomes a partial. No explicit primality certification at that intake; equal opaque cofactors can cancel exactly. SLP implementation, not a DLP splitter. fmpz/native-limb paths and architecture-specific arithmetic. | Partial hash counts and growable in-memory relation arrays feed sorting, merging equal cofactors, deduplication and light filtering. Factor caller triggers processing at about 1.10*base cardinality + multiplier primes + extras, then matrix reduction/Lanczos. Insufficient processed relations can cause base growth/restart. | Matching opaque composite residuals is mathematically possible but requires a new relation type/compatibility policy in v2. Keep prime certainty now. Its target and storage design support measuring post-filter surplus, not copying 60B or growing without limits. |

Consequential source corrections: FLINT at this pin stores relation records in
growable in-memory arrays (`large_prime_variant.c`), not an on-disk partial
file. Its old documentation/storage terminology must not substitute for the
current code.

The older ledger describes Yamaquasi's
`double_large_factor` as a product limit. It is a **candidate-score target**
in the inspected caller; `fbase::cofactor` explicitly retains larger recovered
products if endpoints fit. C1's strict v2 product policy remains exactly what
its historical experiment tested. No historical result is rewritten.

Two caller examples help normalize these choices. JavaMath Factory02 uses
`sieveHitExponent = 0.21 + max(0, bits−150)*0.045/150`, while trial-test and
smooth-residual exponents are 0.12 below 250 bits and 0.16 thereafter. Near
50 digits this gives roughly a 10^6 exact residual ceiling, even though the
initial sieve-hit allowance is much larger. Skipped-small-prime compensation
and post-filtering explain why those domains differ. Java's `BigInteger`,
custom unsigned division and word Montgomery splitters supply a different
cost model; its `findSingleFactor` calls are not v2's explicit shared quotas.
For Yamaquasi near 166 bits, the nominal DLP scoring factor is 83B² and the
common endpoint multiplier is 109B; the recovered product can still exceed
the nominal score factor. Neither source's logarithm/byte constants map to
v2's 32-units-per-bit conservative integer score. Yamaquasi's cofactor
caller tries rho three times at 2,048 evaluations below 40 residual bits,
three at 4,096 below 48 bits, and once at 8,192 otherwise, before its ECM
fallback. That size-dependent effort schedule is meaningful; its numerical
limits and native ECM economics have not been calibrated for PyPy.

YAFU's setup selects DLP at 70 decimal digits, ordinarily postpones TLP to
an effectively disabled 221-digit cutoff unless forced, and then explicitly
assigns its bit-tiered `dlp_exp` in `SIQS.c`. Suggestive exponent columns in
architecture-specific parameter tables are not the operative assignment.
Higher-arity batch paths and GPU variants are distinct callers, not the DLP
path evaluated here. FLINT may continue collection when relation processing
returns insufficient rows; its later failed solve/extraction path can restart
with a larger base. These control-flow distinctions prevent mistaking a
processing trigger for a guaranteed stopping target.

## Storage counts that must not be equated

A v2 `max_partials=8192` DLP configuration caps **unowned spanning-forest
edges**, not all relations. Cycle paths pin forest edges, closing edges remain
in the original-atom store, smooth atoms need no nontrivial LP edge, and
published combined rows have a separate cap. A graph with 8,000 unowned and
1,600 pinned edges can therefore coexist with 1,600 rows and more original
atoms. The 4-KiB/edge reservation covers simultaneous old/new forest maps and
scratch; it is not a native packed-edge size or a measured Python object size.

YAFU/msieve hash buckets locate **vertices**; their allocated vertex arrays,
saved original observations, discovered cycle lists and final matrix columns
are different counts. A cycle can consume several original records. Yamaquasi
stores relation representatives in an SLP hash and DLP endpoint-pair maps;
replacement changes retained original information, not just table occupancy.
JavaMath's factor-to-partial adjacency contains references to shared partials;
summing list lengths double-counts them. FLINT's partial hash table and its
relation-array capacities are likewise separate. An LP-degree-one edge and
a factor-base matrix singleton occur in different incidence systems. Neither
raw record count nor LP cycle rank is a useful-dependency count.

## Source identities, terms and regime

The GitHub HEAD check returned the same YAFU/Yamaquasi/JavaMath/FLINT pins as
C1 on this date. `radii/msieve` still points at historical `c8727d9` (2011);
C9 also inspected `MersenneForum/msieve` at
`03e0ab9e13c4862e5bb2c51a4a39718aa2500e19` (2016). The retained official 1.53
archive has SHA-256
`c5fcbaaff266a43aa8bca55239d5b087d3e3f138d1a95d75b776c04ce4d93bb4`.
The relevant 1.53 rules above agree with the maintained mirror. A maintained
repository label does not imply a recent QS algorithm change.

YAFU/msieve inspected files have public-domain notices; dependencies have
separate terms. Yamaquasi BSD-3-Clause; JavaMath GPL-3-or-later; FLINT QS
LGPL-3-or-later. Factor has no top-level license. These are independently
implemented mathematical/policy ideas, not permission to paste source.
Upstream resource regimes include growable native/Java memory and disk files,
64-bit cofactoring and probable-prime/coverage assumptions. C9 uses PyPy 3.11
native integers, proven residual endpoints, finite owned reservations, shared
wall/CPU/work budgets, and Apple M4 / 24-GiB hardware. No native benchmark is
used to estimate its speed or select a digit cutoff.

## Technical notes and practical pitfalls

[Programming Praxis's large-prime exercise](https://programmingpraxis.com/2013/04/02/quadratic-sieve-with-large-primes/)
illustrates sorting equal partials into complete relations. It is a discovery
and teaching reference, not evidence for a production retention capacity.
The [FLINT development discussion](https://groups.google.com/g/flint-devel/c/qDqu17-J_O8)
points toward expensive cofactoring, relation processing and tuning. Follow-up
inspection of FLINT's actual intake and caller above distinguishes its SLP
implementation from proposed DLP work. Yamaquasi's author comments on density,
shorter paths and rare disconnected cycles are useful hypotheses; the exact
cycle oracle, not an empirical comment, decides supported arithmetic cases.

The maintained [msieve author's `Readme.qs`](https://github.com/MersenneForum/msieve/blob/03e0ab9e13c4862e5bb2c51a4a39718aa2500e19/Readme.qs)
adds an operational warning: distributed savefile progress is not additive,
and duplicate records can inflate the apparent relation target until combining
fails. This is why C9 distinguishes records, unique original atoms, graph
cycles, filtered rows and independent dependencies. Persisting more records
can help collision history, but also costs sorting, duplicate removal and
verification. YAFU's old `attempted_squfof` counter now measures micro-ECM calls;
its disabled (`&& 0`) DLP batch branch is not the active splitter. Source
variable names alone are insufficient attribution.

A base-2 probable-prime rejection before splitting is also a deliberate yield
heuristic: some composite products pass that test. For example,
2047 = 23×89 satisfies 2¹¹ ≡ 1 (mod 2047), hence the base-2 Fermat condition.
Such a screen can discard a valid two-prime residual when both endpoints fit.
v2 only skips splitting after its declared **proven** classification; the
native shortcut cannot be imported as a certainty-preserving optimization.
FLINT's equal-opaque-cofactor cancellation is a different algebraic idea:
exact equality can make the product square without proving each cofactor
prime, but v2 would need an explicit new relation contract to use it.

Concrete example: with v2 base_bound=30,000, a prime residual 15,485,863 is
admissible SLP under 900M but cannot touch a DLP edge bounded by 3M. It can
still pair with another 15,485,863 partial, so discarding it is **lossy policy**,
not a proof it is useless. A residual p*q with both ≤3M and product ≤115.2B
continues through the exact splitter under either endpoint-aligned or current
DLP policy. The candidate allowance remains 115.2B in both; score/refinement
coverage is unchanged. An SLP-only low-bound comparison does narrow scoring
and must not be attributed purely to residual admission.

## Ranked transfers before new training

| Rank / mechanism | Rationale and existing support | Missing algorithm or setting; expected cost | Falsifying/supporting experiment |
| --- | --- | --- | --- |
| 1. Align SLP admission to DLP endpoint | Excludes singletons that cannot join DLP paths; v2 already exposes separate exact bounds. | Uncalibrated setting. Fewer retained large SLPs/rebuilds, but loses rare SLP pairs. No new arithmetic or graph implementation. | Frozen 50-digit full-base DLP, residual_bound=3M vs 900M, equal product/candidate/8192-edge/memory/time limits; complete calls and filtered dependencies, not merely eviction. |
| 2. Retain 16,384 unowned edges | More collision history may permit graph growth; v2's existing capacity supports it. | Setting, not new graph. More objects, rebuild/traversal work, simultaneous reservations and possible 256-MiB refusal. | Isolate capacity at unchanged domains/base; compare with rank 1, then at most one selected combination. |
| 3. Endpoint/product effort shaping | Rare high endpoints tend to need more reports/splitting per useful collision. Existing explicit bounds and rho counters. | Settings initially; source product targets are not equivalent. Rejects may lose recoverable cycles. | Representative recovered/rejected census by size, current failures/endpoint rejects, then a separately frozen product arm only if diagnostic allowance and evidence justify it. |
| 4. Smaller base plus selected policy | C1 half-base 50-digit prefix survived filtering and separate resume completed. | Joint collector setting, C10 interface. More difficult smoothness but smaller parity universe; source and target change together. | Only after individual C9 controls; new explicit freeze, charge all resources and label bundle effect. |
| Defer: new splitter, opaque cofactors, leaf pruning, shortest-cycle replacement, disk spill | Potential savings depend on observed split cost, unsupported composites, eviction relevance or long paths. | New algorithms with certification/provenance/ownership risk. Existing splitter succeeds; graph cycles are short in C1. Leaf status in a prefix does not prove future uselessness. | Require failures/cost/lost-match attribution before implementation. F6 and matrix-owner gates remain intact. |

## Frozen v2 route and remaining algorithm gaps

`SieveCollector` forms a conservative block threshold and refines it against
an exactly evaluated norm. Powers/bucket recovery divides every hit valuation
to completion and adds the supported A exponents. For DLP, candidate allowance
zero resolves to `max(residual_bound, large_product_bound)`. Thus alignment
preserves this experiment's score domain because the product bound dominates
both SLP ceilings. `divide_large_residual` then applies the separate exact
product and SLP bounds; a prime residual above the SLP ceiling is rejected
before graph mutation.

`split_two_primes` tests a repeated square first, otherwise reserves worst-case
work for two Brent-rho walks of at most 2,048 evaluations each (including
recovery), batch 32, seed 7. Each endpoint must fit its bound and receive a
**proven** classification, with exact product reconstruction. A failed walk
never establishes primality. The 131,072-call allowance is cumulative and a
started attempt survives cancellation/restore. Upstream probable-prime
shortcuts are therefore not drop-in replacements.

There is also a useful arithmetic explanation for this frozen residual
population. Under the checked complete factor-base division invariant, each
remaining prime is at least the numeric base bound B. Here B is 15,000 or
30,000 and the product cap is 128B² < B³. Consequently an admissible composite
cannot contain three such prime factors, counting multiplicity. This explains
why `not_two_primes` rejection is not an expected dominant loss in these
bundles. It does **not** remove endpoint certification or verification:
other public configurations allow different products, and corrupted or
incompletely recovered atoms must still be rejected.

SLP-only `_admit` removes its pending representative after a match. For k
observations at one endpoint without intervening eviction it publishes
floor(k/2) disjoint pairs. The
DLP forest retains pinned paths and can publish k−1 independent combinations.
This is an existing **SLP matching-policy gap**, not a DLP correctness gap.
The completed training captures quantify 31/38 additional LP incidence basis
vectors among retained high-bound SLP atoms, and 104/123 among retained
low-bound SLP atoms, beyond the published row counts. Those vectors have not
been lifted to new factor-base dependencies or timed to complete factors.
A reusable SLP representative is a plausible later bounded algorithm transfer;
it requires explicit reference ownership, duplicate/restore identities and
arithmetic-oracle tests. The fixed C9.1 grid is not enlarged after selection.

DLP planning is transactional: reroot the smaller component, find its unique
forest path, or rebuild after FIFO eviction of one eighth of the unowned
edge allowance. Atoms used by published cycles are pinned. A closing edge
is retained as original provenance rather than added to connectivity. The
256-atom combination cap is an explicit yield loss; all-component paths,
loops and parallel edges are supported. Verification, scratch/cap checks and
a final budget poll precede publication. Full smooth atoms use the existing
full-row path. There is no opaque-cofactor or higher-arity relation type.

The pipeline prepares original-atom rows, removes recursive singletons,
merges weight-two columns, lifts Gaussian-elimination kernels and reconstructs
squares before both proper-divisor GCDs. Growth/excess controls are 32 in this
freeze, but **polynomial-window finalization also forces a solve** on changed
stores. They are not a global promise to wait for 32 surplus rows. On training
input 0 full-base alignment prepares/filters 686 times, whereas the half-base
bundle's smaller parity universe costs far less overall. Changing that cadence
could help some late prefixes, but it may delay useful early dependencies;
collection scheduling and matrix/provenance ownership gates remain separate.
A future test should compare cumulative complete costs and preserve final
storage-cap extraction, cancellation and charged resume, not merely reduce
solve counts.

## Existing v2 evidence and coverage

Read production `large_primes`, `large_prime_collector`, `sieve_collector`,
`siqs`, checkpoints and independent DLP tests, plus versioned C1 configurations
and local raw training/confirmation/resume captures. C1 50-digit full-base is
30,000/65,536 interval half-width/6 A factors; half-base changes the numeric
base to 15,000 and its residual bounds, not just duration. Both use 8,192
unowned partials, 32,768 atoms, 8,192 rows. Training used 256 MiB; the completing
resume used 512 MiB, 8-MiB checkpoint and a 300-second cumulative allowance.

C1's full residual census includes both selected and rejected domains. Its
32,768 uniformly chosen positions per input cross-check exact trial division,
including rejects, but stop near 8.4M sampled-prefix positions. The 50-digit
recoverable-DLP census losses were 9004/7707 at block score, 2544/2477 at
refinement and 8528/9989 at SLP residual policy. These are **old collector/census**
measurements, not estimates for new later prefixes. C9's first policy comparison
preserves the DLP score domain, so no unseen rejected-position factorization
assumption is required. Split failure, duplicate loss and long-cycle losses
were zero in C1 fresh confirmation; filter loss and eviction were nonzero.

The new [bounded protocol](c9_protocol.md) fixes longer observations and the
initial finite experiment. Findings and decisions will be appended separately;
this ranking is a research investment decision, not a performance claim.
