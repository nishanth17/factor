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
parameters. This supports measuring long prefixes, not extrapolating a 60-second
empty matrix or copying an 80-digit threshold.

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

## Implementation comparison: follow the actual intake

B means the largest factor-base prime below, unless explicitly stated. v2's
`base_bound` is an exclusive numeric limit; upstream factor-base cardinalities
are not interchangeable with it. All native byte/log thresholds need conversion.

| Implementation / source path | Candidate, residual and effort policy | Matching, retention and useful dependencies | Transfer and cost |
| --- | --- | --- | --- |
| **YAFU** `8110dfb`; `factor/qs/SIQS.c`, `siqs_aux.c`, `tdiv_scan*.c`, `tdiv.c`, `filter.c` | SIQS setup chooses bit/digit-dependent base, reporting thresholds and LP multiplier. Default DLP starts at 70 digits; forced paths and higher arity are separate. SLP intake is residual < L. DLP requires >B², <L^e and <64 bits; e is 1.75/1.8/1.9/1.95 across bit tiers <270/<290/<330/otherwise. Base-2 PRP rejection precedes `getfactor_uecm`; unsuccessful splitting drops the candidate. Both endpoints must be <L. This path validates divisibility, not v2-style endpoint certificates. | Buffered/savefile relation records, growable vertex hash/cycle tables, duplicate/singleton passes and cycle construction. Table doubling is vertex storage, not a fixed retained-edge count. Relation-list growth and disk records are separate. Collection starts from base cardinality +64, with progress/extra-relation updates in callers. | Common SLP/DLP L motivates endpoint alignment. Do not copy PRP rejection, unsigned limits or native micro-ECM speed claims. Growing disk-backed stores would require a new finite storage contract. |
| **msieve** official 1.53 and maintained mirror `03e0ab9`; `mpqs/mpqs.c`, `sieve.c`, `relation.c`, `sqrt.c` | Size-interpolated base/cardinality/interval/LP multiplier; report threshold derived from residual allowance. DLP starts at 282 bits. SLP first, then product about L^1.8, base-2 rejection and SQUFOF; two endpoints <L. Full squares handled. Failed split discarded. Custom multiprecision and native word arithmetic, not PyPy int. | Text savefile records persist independently of the cycle table (initial allocation 10,000, growing). Filter removes singleton/duplicate relations, roots every component, constructs and sorts cycles. Traceback arrays impose 100 edges per side. Default target is fb_size+96 (64 extra solver relations). Verification reconstructs saved relations before square root. | Common L and excluding rare large-large pairs are useful policy hypotheses. v2 already has all-component cycles and exact verification, and must retain disconnected cases rather than adopt path-loss heuristics silently. Native slots/records are not 4-KiB Python forest reservations. |
| **Yamaquasi** `3f95f43`; `src/siqs.rs`, `fbase.rs`, `pollard_rho.rs`, `relations.rs` | Default DLP above 256 bits with half-size base preset. L=min(B*large_factor,2³²−1); multiplier near 166 bits is floor(166*66/100)=109. Nominal DLP score allowance B²*floor(bits/2) below 256 bits. **Caller distinction:** `cofactor` does not enforce that nominal product: it admits lucky reports up to L², splitting above B² with rho then ECM fallback after a composite screen. Endpoints individually ≤L. | Packed partial hash map and ordered double maps; greedy propagation from SLP-connected component, shortest available representatives, equal-endpoint squares. Duplicate endpoint-pair map keys overwrite. All-even combined-factor rows can be discarded as a heuristic. No finite FIFO cap. `gap` requests covered columns+48 cycles; final sparse solve/extraction follows. | A score target is not an admission bound. C9 keeps conservative v2 coverage and explicit product admission. Common endpoint ceiling and retention are candidates; its ignored disconnected cycles and all-even-row discard are incompatible with supported v2 arithmetic cases. |
| **JavaMath** `088d01f`; `SIQS.java`, `SieveParamsFactory02`, `TDiv_QS_2LP`, `CongruenceCollector01`, `PartialSolver01` | Factory separates sieve-hit, trial-test and smooth residual bounds using floating powers of N; callers predivide small primes. After exact division: residual 1 is smooth, residual <B² or passing PRP is SLP, normal endpoints limited to 31 bits. Composite splitting routes Hart <46 bits, tiny ECM <63, rho <64, otherwise nested SIQS; too-large/failing outcomes rejected. Repeated factors become `Smooth1LargeSquare`. | SIQS constructor actually uses Collector01, not every available CycleFinder. Factor→partial lists gather connected partial sets for a partial solver; only independent partials retained after finding a smooth combination. Dynamic Java containers. Starts full solver at base cardinality+caller-specified extras, then adds extras after failure. | Confirms splitting effort should depend on residual size, and surplus needs extraction feedback. Native/JVM backend and probable-prime semantics differ. No evidence for new v2 splitter while current split failures are zero. Related-set solve can become expensive; do not substitute it for the bounded forest. |
| **FLINT** `a4c9750`; `qsieve/collect_relations.c`, `factor.c`, `large_prime_variant.c`, `linalg.c` | Small-prime refinement precedes full division. Full residual 1 accepted; small residual (30-bit bound) <60B and coprime to multiplier becomes a partial. No explicit primality certification at that intake; equal opaque cofactors can cancel exactly. SLP implementation, not a DLP splitter. fmpz/native-limb paths and architecture-specific arithmetic. | Partial hash counts and growable in-memory relation arrays feed sorting, merging equal cofactors, deduplication and light filtering. Factor caller triggers processing at about 1.10*base cardinality + multiplier primes + extras, then matrix reduction/Lanczos. Insufficient processed relations can cause base growth/restart. | Matching opaque composite residuals is mathematically possible but requires a new relation type/compatibility policy in v2. Keep prime certainty now. Its target and storage design support measuring post-filter surplus, not copying 60B or growing without limits. |

Consequential source corrections: FLINT at this pin stores relation records in
growable in-memory arrays (`large_prime_variant.c`), not an on-disk partial
file. Its old documentation/storage terminology must not substitute for the
current code.

 the older ledger describes Yamaquasi's
`double_large_factor` as a product limit. It is a **candidate-score target**
in the inspected caller; `fbase::cofactor` explicitly retains larger recovered
products if endpoints fit. C1's strict v2 product policy remains exactly what
its historical experiment tested. No historical result is rewritten.

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

Concrete example: with v2 base_bound=30,000, a prime residual 50,000,003 is
admissible SLP under 900M but cannot touch a DLP edge bounded by 3M. It can
still pair with another 50,000,003 partial, so discarding it is **lossy policy**,
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
