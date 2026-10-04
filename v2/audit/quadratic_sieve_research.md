# Quadratic sieve research for Factor Phase 3

Date: 3 October 2026. Milestones: M19, M21. Scope: source and paper review for
Factor's PyPy Python 3.11 implementation in `v2/`.

## Findings and implementation decision

Build an exact, bounded relation pipeline, then SIQS self-initialization.
The most useful early performance experiments are polynomial/root reuse,
bounded sieve-score storage, candidate selection, and avoiding unnecessary
factor-base division. Single-large-prime relations and Python-integer GF(2)
bitsets provide the first complete baseline. Design the relation interface
to accommodate double-large-prime SIQS and Smooth Subsum Search (SSS).

The [M21 follow-up](#optimization-follow-up-m21) adds concrete experiments
from FLINT, implementation blogs, polynomial-selection literature, sparse
filtering and batch smoothness research. Fold small-prime omission, duplicate
avoidance, weight-two elimination and bounded collection recovery into the
current plan. Keep more elaborate polynomial variants and sparse solvers
conditional on measured benefit.

This review identifies implementations to study and hypotheses to measure;
it does not establish a fastest implementation on our machine. No competitor
was installed or executed, and no Factor runtime code changed. Published
native timings, Python timings, and relation-throughput measurements have
different meanings. The [M17 baseline][baseline] remains the control.
The separately accepted [M20 reference collector][p31-acceptance] supplies
the exact small-window control for P3.2; it is not a complete SIQS benchmark.

The implementation work is folded into [TODOS.md][todos]. GNFS follows the
existing [dependency order][gnfs-order]: working SIQS, early arithmetic/GMP
assessment, small GNFS, then scaling. Optional QS optimizations do not become
blanket prerequisites for GNFS.

## Sources and reproducibility

GitHub files were fetched from immutable revisions using HTTPS; selected
routines, documentation, and dependency/license declarations were inspected.
This is a targeted design review, not a correctness audit of every file.
The compact [source manifest][manifest] records the URL, revision, byte size,
SHA-256, and unexecuted status for each fetched file. Temporary downloaded
copies are removed after verification; neither full repositories nor native
binaries are retained. Existing Phase 2 competitor pins remain unchanged.

| Implementation | Reviewed revision | Role in this pass |
| --- | --- | --- |
| [YAFU][yafu] | `963dbe9c45283cc06e9a71d830b0676bc1b0d343` | Current native SIQS engineering reference; commit dated 24 September 2026 |
| [Msieve mirror][msieve] | `c8727d91305bdbe0972d160ef0ce61dd02ce9193` | Historical polynomial, relation and matrix reference; commit dated 21 February 2011 |
| [Yamaquasi][yama] | `3f95f43682ed15d8c1ed206a9a702dd655d7c8ad` | Rust QS/SIQS reference with explicit storage and matrix choices; commit dated 23 November 2025 |
| [SymPy][sympy] | `06fdc356176f6ae82c865a365701f64c5ee67aa9` | Current Python single-large-prime and bitset reference; reused M12 pin |
| [primefac][primefac] | `2401341d8689ea2ab517ef311f4f91ecd33d14d0` | Python SIQS baseline and parameter reference; reused M12 pin |
| [PyFactorise][pyfactorise] | `e96079fea53bb279fe520b3e530c368333ac548f` | Related Python polynomial/matrix implementation; reused M12 pin |
| [numthy][numthy] | `1496df97faec010c909e54e14d6bf9cc4399e5fe` | Python bytearray collector and triple-large-prime reference; reused M12 pin |
| [SSS author-associated code][sss] | `8dbaf6d39ab88a40380965d25ec2c363d7f27358` | Alternative collector and original comparison harness; reused M12 pin |
| [FLINT][flint] | `00cc19b350b5302876b4fa8c877df5a676c2061c` | M21 addition: maintained native single-large-prime SIQS, polynomial selection and execution-policy comparison; commit dated 2 October 2026 |

These are complementary references, not a performance ranking. In
particular, the old Msieve mirror's README superlatives cannot establish
present-day leadership. YAFU's newer source is more useful for current native
engineering, while Python references are closer to our supported runtime.

### Native SIQS references

**YAFU.** Its [trial-division source][yafu-tdiv] separates omitted-small-prime
checks, ordinary division, resieving, bucket-hit lookup, and residual
cofactorization. It also contains three-large-prime splitting; its simpler
[function documentation][yafu-doc] mainly describes double-large-prime SIQS.
[Polynomial generation][yafu-poly] and [root setup][yafu-roots] expose cached
CRT components and root-update constants. [SIQS orchestration][yafu-siqs]
and [batch cofactorization][yafu-batch] show the additional queue and budget
surface of advanced collection.

For Factor, borrow the stage separation and measure each stage. The source's
SIMD kernels, prime cutoffs, block widths, and native residual splitters are
not PyPy defaults. Native build configuration must be disclosed in a future
reference comparison. A native factorizer remains an external reference;
Factor's collector and factor extraction stay in Python.

**Msieve.** [Readme.qs][msieve-qs] describes self-initializing MPQS with two
large primes and separate collection/postprocessing. The
[relation code][msieve-rel] tracks cycles and duplicate relations;
[GF(2) processing][msieve-gf2] calls block Lanczos. It is a useful
algorithmic reference, but
the pinned mirror is historical. Its old NFS crossover discussion and
distributed-computing assumptions do not set our thresholds.

**Yamaquasi.** Its [QS notes][yama-qs] disclose fixed-width Rust arithmetic,
in-memory relations, and a Gauss/block-Lanczos switch. The notes explicitly
qualify their benchmark methodology. [Sieve source][yama-sieve] partitions
small and sparse large-prime hits into different storage strategies.
[Relations][yama-rel] retain exponents and verify modular identities; their
large-prime combination strategy concentrates on the component connected
to single-large-prime relations. [Matrix code][yama-gf2] uses randomness for
block Lanczos.

Our graph oracle must include cycles disconnected from the single-prime
component. Choosing to skip such cycles may be a measured yield tradeoff;
it is not a general proof that they never occur. Likewise, Rust's matrix
size threshold and storage bounds cannot be transferred to Python objects.

### Python references

**SymPy.** [qs.py][sympy-qs] separates polynomial generation, integer-scaled
log scores, exact candidate division, matching residual primes, and bitset
elimination. The polynomial stores the full square difference, whereas our
planned normalized polynomial separates `A`. This distinction matters when
handling primes dividing `A`. Its `qs` wrapper returns a set without
multiplicity; `qs_factor` returns a mapping whose keys may still be composite.
Its optional arithmetic helper backend must be recorded in comparisons.
The [QS tests][sympy-tests] provide additional boundary/reference fixtures.

Use it as a design and comparison reference. A wrapper return value is not
independent primality or reconstruction evidence, and its public interface
does not supply Factor's shared budget/checkpoint contract.

**primefac and PyFactorise.** [primefac's SIQS][primefac-code] identifies its
PyFactorise ancestry. It uses parameter tables, polynomial families and
integer-bitset elimination; the reviewed collector accepts fully smooth
values rather than implementing the large-prime combination used by SymPy.
It computes some targets through `float(n)`. [PyFactorise][pyfactorise-code]
is useful for inspecting the related CRT/Gray-code transitions and matrix
layout, but its custom modular powering is not a measured improvement for
Factor.

Compare these related implementations without counting them as independent
algorithm families. Their parameter tables are initial research controls,
not calibrated settings. Factor's polynomial targets and identities must
remain exact over the declared input domain.

**numthy.** The reviewed [SIQS routines][numthy-code] use integer roots for
polynomial targeting, cached saturating byte-subtraction translation tables,
slice updates, and threshold-mask extraction. They handle up to three large
primes through streaming parity elimination and combination masks. Residual
factoring scans a precomputed allowed-prime list; polynomial generation uses
`SystemRandom`. The [README][numthy-readme] advertises dependency-free SIQS.

The bytearray translation approach is a high-priority PyPy experiment.
Its slices allocate and copy, so fewer Python loop iterations do not prove a
win. Seeded assignments, memory accounting, candidate verification, and
resume semantics still need our own implementation. Triple-large-prime
yield must be weighed against residual factoring and growing provenance.

### SSS evidence and qualification

The [SSS paper, version 2][sss-paper] compares selected Python SIQS baselines
using shared linear algebra. Its 30–70-digit experiments complete
factorizations; the 75–100-digit experiments compare relations after one
hour. Thus the reported approximately 5–7-fold gains in the latter range are
collection evidence, not measured complete-factorization speedups. The
primefac baseline lacks a large-prime variant, and no tested method uses
multiplier selection. The paper provides no PyPy-specific result. These
limitations make SSS a serious challenger rather than a demonstrated winner
over Factor or modern native SIQS.

The pinned [collector][sss-code] forms CRT/collision candidates and checks
smooth parts with product/remainder trees. Its limited prime-power product
can miss smooth candidates; exact verification is still required for accepted
relations. It depends on SymPy, prints progress, uses process-exit paths and
floating-point initialization, and does not expose our bounded API.
The [matrix helper][sss-matrix] imports `gmpy2.isqrt`; [test.py][sss-test]
imports NumPy for statistics and uses unseeded input generation. These are
properties of the upstream harness, not new Factor dependencies.

Use two clearly labelled future arms: reproduce compatible upstream settings,
then compare an independently implemented collector adapter with Factor's
exact verifier, seeded assignments, budgets and common postprocessing.
Record every adaptation. Do not quietly remove a dependency or replace the
matrix code and call the result an unchanged upstream reproduction.

### Foundational and recent optimization papers

[Boender and te Riele's large-prime study][large-primes] explains why
single/double-large-prime performance depends on implementation, parameters
and available memory. Its experiments concern different hardware and old
implementations. It supports measuring the tradeoff, not copying a digit
crossover or advertised multiplier into our dispatcher.

[Guan, Zhuang and Mastorakis][recent-paper] discuss four-hit loop unrolling
and binary search in parameter tables. Their pseudocode updates all four
hits and handles tails; merely stepping by `4*p` would omit valid hits.
For a fixed unroll factor, `O(n/4)` is still `O(n)`. The tables describe
bit lengths, and the implementation is compiled C. Our decision is to defer
unrolling until PyPy profiling identifies marking-loop overhead. Tiny
parameter tables are not an initial optimization target. This paper is an
implementation idea source, not evidence of a new factoring asymptotic or
of a speedup for this project.

## Optimization follow-up (M21)

**Decision: yes, fold the literature into the existing tasks.** The useful
additions concern useful relations per total CPU-second, duplicate work and
matrix cost. A new optimization phase is unnecessary. This follow-up adds a
ninth implementation reference, three research papers and course notes,
with public implementation write-ups as supporting evidence. The FLINT
revision and fetched-byte hashes are in the [M21 supplement][supplement];
the original M19 manifest and acceptance capture remain historical evidence.
No upstream implementation was executed and no local speedup is claimed.

### Additional sources and what they establish

**FLINT source and an implementer's blog.** The reviewed
[collector][flint-collect] omits the smallest marking primes, checks them on
candidates, tests sieve roots before division and uses score-based early
termination. It blocks the sieve and separates dense from sparse hits.
[Polynomial setup][flint-poly] selects distinct A coefficients near an
integer target; [multiplier selection][flint-multiplier] balances small-prime
residues, the class modulo 8 and the size penalty for multiplying n.
[Orchestration][flint-factor] counts full relations and graph cycles,
filters before solving, and enlarges the factor base after failure or family
exhaustion. These are inspectable engineering references, not PyPy defaults.

[William Hart's February 2017 post][hart-blog] explains the FLINT work he
implemented: small-prime omission, cache blocking, polynomial-selection
variants and robustness on exhausted polynomial schedules. It helps connect
the algorithms to implementation costs. Its performance discussion is
historical, and its OpenMP and machine-word techniques do not establish
Python thresholds. The current pinned source, rather than the blog's account
of future work, establishes which paths were inspected in this pass.

**Bradford, Monagan and Percival, sections 2.3–2.5.** Their
[Maple/C implementation paper][bmp-paper] discusses A's factor count,
duplicate relations, omitted-small-prime estimates, division early exit and
cycle-based collection progress. More factors in A amortize setup over more
polynomials but also change smoothness behavior. Their polynomial strategy
fixes a product A0 and varies a prime q above the factor-base bound, reusing
initialization across A=A0*q. This is a separate variant to evaluate; it is
not a drop-in version of Factor's planned factor-base-smooth A families.
Historical timings and tuned statistical allowances are not portable results.

**Cavallar, sections 2–3.** [Sparse filtering research][filter-paper] studies
the tradeoff between reducing matrix dimensions and increasing nonzeros.
Eliminating a constraint present in exactly two relations needs their XOR
and introduces no fill-in; higher-frequency merges can make rows denser.
The paper concerns NFS ideal relations. Applying its GF(2) reduction ideas
to QS is an engineering inference; NFS ideals and extraction rules are not
being copied into the QS relation contract.

**Bernstein, smooth-part algorithms.** The
[2004 author-hosted draft][smoothparts-paper] and
[2006 course notes, sections 14.2–14.5][bernstein-notes] give a primary basis
for the existing product/remainder-tree experiment. Computing smooth parts
in batches and recovering individual prime factors are distinct operations.
These asymptotic algorithms support P6.2's investigation, but do not supply
a batch size or a measured PyPy advantage for our candidate stream.

The author's [Programming Praxis QS exercise][praxis-blog] supplies a small
worked implementation and a pointer to Contini's thesis. It is useful for
understanding the initial pipeline, not a competitive SIQS baseline. The
Contini PDF could not be retrieved in this pass; no new technical finding is
attributed to unread thesis text. The online write-ups are selected for their
authors' implementation experience rather than general claims of speed.

### Changes worth making to the plan now

The following are Factor experiment designs inferred from those sources,
with correctness and resource requirements added for this repository.

1. **Separate the sieve window from its working blocks (P3.2).** Sweep both
   sieve-position blocks and chunks of factor-base metadata. Compare dense
   marking, sparse direct hits and bounded buckets; count every hit, cache
   construction, Python object and decode cost. A sparse bucket can cost more
   than direct marking in Python. Retain a simple loop control and choose
   cutoffs from PyPy measurements rather than native cache sizes.
2. **Make small-prime omission a measured candidate policy (P3.2).** Compare
   full marking with omitting a bounded prefix and recovering exact small
   factors before costly residual work. Measure candidate loss and rejection
   at both stages, including unusually high small-prime powers. An average
   contribution is a heuristic, not an upper bound. Keep score-guided division
   early exit conditional on a proved remaining-hit invariant; saturated or
   approximate scores cannot by themselves certify complete exponent recovery.
3. **Tune family diversity as well as Gray-code switching (P3.4).** Use
   reproducible, bounded A schedules and canonical family identities. Measure
   setup amortization, A distance from its integer target, distinct relations,
   repeated A sets and symmetry/translation duplicates. Vary the number and
   range of factors in A. Explicit Knuth–Schroeppel scoring gets an h=1
   control, a finite multiplier list and setup charges; exact residue tests
   and integer magnitudes never require converting arbitrary n to float.
4. **Add weight-two elimination before a new solver (P3.3).** Keep maps to
   original verified relations and an explicit matrix orientation. A prime
   column occurring in two relation rows can be removed by replacing those
   rows with their XOR. Capture dependencies when a row becomes zero.
   Distinct relations with identical parity already give a two-row dependency;
   parity equality alone is not permission to delete one as a duplicate.
   Compare dimensions, nonzeros, provenance bytes and extraction cost. Higher
   merges remain a capped experiment if the simple reduction is insufficient.
5. **Stop and recover using useful relations (P3.3–P3.4).** Track distinct
   full/combined relations, partial matches, filtered dimensions and attempted
   dependencies. For a double-large-prime multigraph, E−V+C measures its cycle
   space, not successful factors. Test periodic bounded filtering and a
   declared excess of usable rows before solving. Continue after trivial
   dependencies; cap every retry. If families are exhausted or yield stalls,
   compare a fixed configuration with a finite factor-base/width growth policy
   under the original budget. Remap relations and base identities explicitly
   on growth; never silently reset consumed resources.
6. **Tune large-prime policy through extraction (P5.4).** Separate each
   residual-prime limit from the composite product cap and splitting budget.
   Screen impossible residual shapes only with a documented mathematical
   justification. Compare raw partials, matches/cycles, cycle lengths, matrix
   weight and completed factorizations. More partials can merely create a
   larger unmatched store or denser combinations. Tune against the optimized
   single-large-prime control, then freeze policy before held-out evaluation.

### Conditional additions and limits

The Bradford–Monagan–Percival A0*q experiment belongs in P6.2 after ordinary
SIQS works. A q outside the factor base appears in A*F(x), so its exponent
and large-prime/provenance role must be explicit. It can consume a large-prime
entry and change matching and duplicate rules. Reusing CRT data without
recording that factor would violate the relation contract. Compare it with
factor-base-smooth A families using the same verifier and total budgets.

Batch smoothness, triple large primes, block Lanczos/Wiedemann, SIMD-inspired
rewrites and parallel collection retain their existing conditional gates.
For a first Python baseline, weight-two filtering and integer bitsets are
less implementation work than a new randomized sparse solver. Native
reciprocal-division and limb kernels are arithmetic-backend references;
reimplementing machine-word assumptions in arbitrary-size Python arithmetic
requires its own proof and measurement. These additions do not postpone GNFS.

## Proposed contracts and architecture

The following are Factor design decisions derived from exact identities,
our [budget contract][budget], and the reviewed stage boundaries. The current
[reference API][qs-reference] implements small polynomial/relation checks
and exhaustive collection. Complete sieving, filtering, extraction and the
optimization proposals still require their implementation/acceptance gates.

### Exact polynomial identity

Start with odd, composite, non-power `n` after existing preprocessing. For
an optional bounded positive multiplier `h`, inspect `gcd(h, n)` first and
validate any split. Set `N_prime = h*n` and use:

```text
B*B == N_prime (mod A)
C = (B*B - N_prime) // A
F(x) = A*x*x + 2*B*x + C
U = A*x + B
U*U - N_prime = A*F(x)
```

Check exact divisibility before constructing `C`. A relation records the sign
and full exponents of `A*F(x)`, including `A`; parity alone is insufficient.
Bound candidate bit sizes and treat `F(x) == 0` separately before division.
The collector's interval is explicitly half-open, with negative positions,
empty blocks, and tails mapped consistently.

For an odd factor-base prime `p` not dividing `A`, derive roots using
`(sqrt(N_prime) - B)*inverse(A, p)` and the negative square root. For odd
`p | A`, the normalized polynomial reduces to `2*B*x + C` modulo `p`;
handle its linear root or degenerate case separately. Primes dividing
`N_prime`, 2, and failed inversions have explicit branches. Never blindly
copy a full-square-difference sieve's skip-`A` behavior into a normalized
`F` sieve.

Compute the target near `isqrt(2*N_prime)//M` and compare candidate `A`
values with integer arithmetic. Cache modular square roots and inverses
per family. Gray-code changes to `B` must include recentering/translation in
both the root updates and persisted polynomial identity. Verify incremental
roots against full recomputation before timing their gain.

### Relations, partials and verification

Use immutable atomic relation IDs and polynomial/position provenance.
Store sparse exponent pairs, parity bits, and explicit square corrections
for combined partials. Keep root-hit data ephemeral and separate from the
mathematical relation. SIQS and SSS adapters submit the same independently
verifiable modular relation contract; GNFS ideal/character data remain a
distinct contract even where storage infrastructure is shared.

Initial residual primes stay in the existing deterministic primality domain
`r < 2**64` and under the collector's smaller configured limit. Charge
residual classification. Larger probable residuals require a separate
documented certainty policy; relation admission must not manufacture proven
prime labels.

For matching partial residual `r`, inspect `gcd(r, n)` before inversion.
When it is 1, either normalize `U1*U2` by `r` modulo `n` or retain the
equivalent square correction explicitly. Verify the constituent exact
factorizations and combined congruence without constructing an unbounded
integer product. A residual atom that has not been validated cannot become
a verified full relation merely because it matches another residual.

For later double-large-prime collection, graph edges retain original atomic
relations. Include repeated primes, self-loops, duplicate edges, isolated
cycles and corruption in the independent combination oracle. A third large
prime changes the combination problem; do not force it into a two-endpoint
graph representation.

### Filtering, dependencies and extraction

Start with exact-duplicate and iterative-singleton removal, with maps to
original relations. Preserve dependencies between distinct rows of equal
parity. Compare weight-two constraint elimination before higher-way merges
or a sparse solver. Use Python-integer parity rows and separately tracked
dependency masks. Reserve workspace for fill-in and provenance, not just
initial sparse rows. Keep complete exponents after parity filtering.

For each dependency, recheck original-row parity and even exponent totals.
Build `X` and `Y` modulo `n` from original relations and half-exponents,
including square corrections. This avoids taking an integer square root of
an enormous product of all relation values. Check `X*X % n == Y*Y % n`,
try both GCD signs, validate the divisor, and continue bounded collection
after trivial dependencies. Relation counts and matrix completion are
progress measures, not factorization success.

Sparse solvers become experiments when filtered matrix/provenance costs
dominate or bitset workspace is infeasible. Record seeds and finite retries;
verify every returned vector against the original matrix. Do not port a
native matrix-size switch as a PyPy threshold.

### Storage, budgets and resume

Apply one shared budget to factor-base setup, family construction, blocks,
candidate division, residual work, filtering, elimination and extraction.
`Budget` is cooperative, so each atomic block and bigint/tree operation also
needs bounded size. Track wall time, CPU and owned workspace; measure process
RSS separately. Existing workspace estimates do not equal an RSS guarantee.

Bound full relations, unmatched partials, caches, provenance masks and
checkpoint output. Persist polynomial family/Gray index, block position,
seed state, verified relation-store identity, configuration and consumed
allowances. Define cap behavior and deterministic unmatched-partial eviction;
never evict an atom still referenced by an accepted combination.

Keep checkpoint metadata compact. If relation spill is justified, use a
versioned compact store with explicit disk caps and integrity/commit markers,
not repeated pretty-printed JSON matrix dumps. Interrupted writes must not
silently become verified relations. No disk store is mandatory for the first
small in-memory baseline. The existing default 8 MiB workspace is a constraint
to evaluate, not permission for unlimited relation storage.

## Prioritized experiments

Every performance result must validate its outputs and include total
time-to-factor under the [promotion policy][promotion]. The table states
priority and decision criteria, not an expected speedup.

| Priority | Experiment | Main comparison and gate | Task |
| --- | --- | --- | --- |
| First | Polynomial families and root reuse | Full recomputation versus incremental updates; exact identities and roots agree, including `p | A` and recentering | P3.1, P3.4 |
| First | Sieve storage | PyPy list loops versus bounded bytearray translation and `array` scores; include slice allocation, saturation, tails, scanning and RSS | P3.2 |
| First | Candidate threshold and prime powers | Exhaustive small-window oracle; measure missed smooth values, false candidates, repeated factors and useful verified relations | P3.2 |
| First | Single-large-prime storage | Full-only versus matched partials; include primality, combination, cap/eviction and provenance costs | P3.2 |
| Next | Avoid full factor-base division | Root-hit filtering/resieving and buckets versus full trial division; include hit-storage construction and decoding | P3.2 |
| First | Small-prime omission and division early exit | Full marking/division control; quantify missed values, exact recovery and both candidate-stage costs; prove any score-based exit invariant | P3.2 |
| Next | Filtering and pivot order | Compare singleton-only and weight-two reduction; preserve equal-parity dependencies; measure nonzeros, provenance and final factors | P3.3 |
| Next | Multiplier and family diversity | `h=1` and simple A schedules versus bounded scored multipliers/diverse families; count setup, duplicates and larger norms | P3.4 |
| Next | Collection stopping and recovery | Filtered row excess and bounded growth versus fixed settings; charge repeated filtering/solving and retained-state remapping | P3.3, P3.4 |
| Challenger | SSS/SSSf adapter | Upstream reproduction and matched collector/postprocessing arms; useful dependencies and final completion, not only relation counts | P3.5 |
| Extension | Double large primes | Single-prime control versus complete graph-cycle handling; charge residual splitting, filtering, memory and final extraction | P5.4, eligible after P3.4 |
| Conditional | Batched smooth-part detection | Only when division dominates: scalar verification versus bounded product/remainder trees, including tree allocation and candidate latency | P6.2 |
| Conditional | A0*q polynomial reuse | Explicitly account for q outside the factor base; compare useful yield, partial policy, setup and extraction with ordinary SIQS | P6.2 |
| Conditional | GMP arithmetic | Persistent big operands on PyPy versus built-in integers; include conversions, imports, checkpoint serialization and full factoring | Early P4.3 |
| Conditional | Processes, NumPy, triple large primes, sparse solvers | Profile-backed one-change comparisons; preserve exact verifier and aggregate budgets | P3.6, P3.7, P6.2 |
| Low | Four-hit unrolling and table lookup | First prove identical hits/tails; prioritize only if these loops actually dominate | P3.2 |

The bytearray experiment must quantify candidate loss caused by scale,
rounding, skipped small primes and finite prime-power allowances. An exact
acceptance verifier prevents false successes but cannot reveal relations
missed by the score filter. Preserve an exhaustive small-window reference
and separate full score-equivalence experiments from intentionally lossy
heuristics. No unchecked byte wraparound is acceptable.

For batch smoothness, distinguish detecting the smooth part from recovering
its full prime exponents. Product/remainder trees can avoid many divisions
while moving cost into bigint multiplication, allocation and recovery.
Predeclare bit/node/storage caps, charge each batch, and verify exponent
recovery independently. Batch throughput alone does not justify default use.

[PyPy's FAQ][pypy] documents C-extension compatibility costs and the need to
warm representative Python work. Therefore NumPy remains the optional P3.7
experiment. Upstream NumPy use for statistics does not imply it speeds up
sieving. [gmpy2's overview][gmpy] identifies `mpz` as GMP-backed arithmetic;
backend availability and the actual PyPy crossover remain unmeasured here.
GMP is an arithmetic library, not a replacement for the SIQS algorithm.

## Measurement plan and phase sequence

1. **P3.1:** exact polynomial/relation types and a tiny exhaustive collector.
   Use independent arithmetic fixtures before making performance claims.
2. **P3.2–P3.3:** bounded full/single-large-prime collection, filtering,
   dependencies and complete factor extraction. Complete small balanced
   fixtures and validate interruption/cap behavior.
3. **P3.4:** SIQS families, dispatch/checkpoints and held-out coverage.
   Compare with QS/MPQS and the frozen bounded ECM portfolio. Then bring
   early P4.3 arithmetic assessment forward as already scheduled.
4. **P3.5:** a collector spike may begin after P3.3; complete comparison uses
   working SIQS. Double-large-prime P5.4 may also move forward when partial
   yield justifies it. Neither is required to start GNFS contracts or engine
   work after the stated prerequisites pass.
5. **Conditional optimizations:** batch trees, parallelism, NumPy, third
   large primes and sparse solvers follow measured bottlenecks.

Use the independent M15 corpus for compatible existing comparisons and add
a separately frozen Phase 3 corpus where needed. Keep known factors hidden
from collectors. Start evaluation with small balanced inputs, then declared
30/40/50/60-digit bands; 70/80-digit runs remain capped exploration.
Include difficult control classes, repeated factors and negative inputs at
the complete-result boundary. Include separately generated `n` residue
classes modulo 8 and multiplier-sensitive cases without tuning on held-out
results. Do not promise completion at a proposed digit band.

The M17 50 ms operation caps were Phase 2 feasibility probes; they are not
automatically appropriate SIQS performance budgets. Choose wall/CPU/work,
workspace/RSS and disk limits before measurements and disclose changes.
Operation counts of different algorithms are not equal CPU costs. Preserve
the old-budget feasibility arm if needed, and use comparable declared
resources for the coverage/promotion arm.

Record factor-base setup, polynomial/root initialization, marking, candidate
scan, exact division, residual factoring, combination, filtering, matrix and
extraction costs. Report verified useful relations/dependencies, matrix
dimensions/nonzeros, memory peaks, storage, completion and budget exhaustion.
Collector-only diagnostics cannot substitute for complete-factorization
promotion. Native library output and total factorization have their own
comparison arm; invoking a general dispatcher is not a SIQS-only comparison.

Keep all Python arms on PyPy Python 3.11; record optional backend/build versions
and availability. Separate cold startup from several seconds of representative
JIT warmup and check sample stability. Match seeds, input assignments,
output obligations and resources. Retain failures as censored outcomes;
compare timings only for validated compatible successes. Nine samples and
the existing independent-input uncertainty analysis are starting controls,
not substitutes for checking spread and representative coverage.

## Status and remaining uncertainty

Research passes M19 and M21 are complete. P3.1's separate M20 acceptance is
recorded in the backlog; P3.2–P3.7 and the later implementation/promotion
gates remain open. This documentation update closes no implementation gate.
The next implementation task is P3.2, using the accepted exact reference.

Unknown until implementation and matched measurement: the best score format,
factor-base/block sizes, multiplier benefit, single/double-large-prime
crossover, SSS advantage, GMP availability/performance, feasible digit bands,
matrix solver threshold and useful parallel grain size. No literature claim
closes these gates for Factor.

The source review did not execute upstream code or verify every external
claim. The scanned Alford/Pomerance and poorly extracted Carrier/Wagstaff
PDFs were discovered but are not relied on for findings in this report.
The accessible Boender/te Riele paper supplies the foundational large-prime
discussion. Any later code adaptation must retain the source's applicable
attribution; the source manifest preserves the inspected license/header
locations rather than assigning one license to YAFU's entire distribution.

[baseline]: phase_two_m17_frozen_baseline.json
[todos]: TODOS.md#phase-3--add-the-missing-balanced-composite-engine
[gnfs-order]: TODOS.md#execution-order-toward-gnfs
[promotion]: TODOS.md#how-to-use-the-gates
[manifest]: quadratic_sieve_research_sources.json
[supplement]: m21_quadratic_sieve_optimization_sources.json
[p31-acceptance]: phase_three_m20_p31_summary.json
[qs-reference]: ../README.md#phase-31-reference-relation-api
[flint]: https://github.com/flintlib/flint/tree/00cc19b350b5302876b4fa8c877df5a676c2061c
[flint-collect]: https://github.com/flintlib/flint/blob/00cc19b350b5302876b4fa8c877df5a676c2061c/src/qsieve/collect_relations.c
[flint-poly]: https://github.com/flintlib/flint/blob/00cc19b350b5302876b4fa8c877df5a676c2061c/src/qsieve/compute_poly_data.c
[flint-multiplier]: https://github.com/flintlib/flint/blob/00cc19b350b5302876b4fa8c877df5a676c2061c/src/qsieve/knuth_schroeppel.c
[flint-factor]: https://github.com/flintlib/flint/blob/00cc19b350b5302876b4fa8c877df5a676c2061c/src/qsieve/factor.c
[hart-blog]: https://wbhart.blogspot.com/2017/02/integer-factorisation-in-flint.html
[bmp-paper]: https://www.cecm.sfu.ca/~pborwein/MITACS/papers/percival.pdf
[filter-paper]: https://ir.cwi.nl/pub/4456/04456D.pdf
[smoothparts-paper]: https://cr.yp.to/factorization/smoothparts-20040510.pdf
[bernstein-notes]: https://cr.yp.to/2006-aws/notes-20060309.pdf
[praxis-blog]: https://programmingpraxis.com/2013/03/19/quadratic-sieve/
[budget]: ../budget.py
[yafu]: https://github.com/bbuhrow/yafu/tree/963dbe9c45283cc06e9a71d830b0676bc1b0d343
[yafu-doc]: https://github.com/bbuhrow/yafu/blob/963dbe9c45283cc06e9a71d830b0676bc1b0d343/docfile.txt
[yafu-tdiv]: https://github.com/bbuhrow/yafu/blob/963dbe9c45283cc06e9a71d830b0676bc1b0d343/factor/qs/tdiv.c
[yafu-poly]: https://github.com/bbuhrow/yafu/blob/963dbe9c45283cc06e9a71d830b0676bc1b0d343/factor/qs/new_poly.c
[yafu-roots]: https://github.com/bbuhrow/yafu/blob/963dbe9c45283cc06e9a71d830b0676bc1b0d343/factor/qs/poly_roots.c
[yafu-siqs]: https://github.com/bbuhrow/yafu/blob/963dbe9c45283cc06e9a71d830b0676bc1b0d343/factor/qs/SIQS.c
[yafu-batch]: https://github.com/bbuhrow/yafu/blob/963dbe9c45283cc06e9a71d830b0676bc1b0d343/factor/batch_factor.c
[msieve]: https://github.com/radii/msieve/tree/c8727d91305bdbe0972d160ef0ce61dd02ce9193
[msieve-qs]: https://github.com/radii/msieve/blob/c8727d91305bdbe0972d160ef0ce61dd02ce9193/Readme.qs
[msieve-rel]: https://github.com/radii/msieve/blob/c8727d91305bdbe0972d160ef0ce61dd02ce9193/mpqs/relation.c
[msieve-gf2]: https://github.com/radii/msieve/blob/c8727d91305bdbe0972d160ef0ce61dd02ce9193/mpqs/gf2.c
[yama]: https://github.com/remyoudompheng/yamaquasi/tree/3f95f43682ed15d8c1ed206a9a702dd655d7c8ad
[yama-qs]: https://github.com/remyoudompheng/yamaquasi/blob/3f95f43682ed15d8c1ed206a9a702dd655d7c8ad/README_qsieve.md
[yama-sieve]: https://github.com/remyoudompheng/yamaquasi/blob/3f95f43682ed15d8c1ed206a9a702dd655d7c8ad/src/sieve.rs
[yama-rel]: https://github.com/remyoudompheng/yamaquasi/blob/3f95f43682ed15d8c1ed206a9a702dd655d7c8ad/src/relations.rs
[yama-gf2]: https://github.com/remyoudompheng/yamaquasi/blob/3f95f43682ed15d8c1ed206a9a702dd655d7c8ad/src/matrix/gf2.rs
[sympy]: https://github.com/sympy/sympy/tree/06fdc356176f6ae82c865a365701f64c5ee67aa9
[sympy-qs]: https://github.com/sympy/sympy/blob/06fdc356176f6ae82c865a365701f64c5ee67aa9/sympy/ntheory/qs.py
[sympy-tests]: https://github.com/sympy/sympy/blob/06fdc356176f6ae82c865a365701f64c5ee67aa9/sympy/ntheory/tests/test_qs.py
[primefac]: https://github.com/lucasaugustus/primefac/tree/2401341d8689ea2ab517ef311f4f91ecd33d14d0
[primefac-code]: https://github.com/lucasaugustus/primefac/blob/2401341d8689ea2ab517ef311f4f91ecd33d14d0/primefac.py
[pyfactorise]: https://github.com/skollmann/PyFactorise/tree/e96079fea53bb279fe520b3e530c368333ac548f
[pyfactorise-code]: https://github.com/skollmann/PyFactorise/blob/e96079fea53bb279fe520b3e530c368333ac548f/factorise.py
[numthy]: https://github.com/ini/numthy/tree/1496df97faec010c909e54e14d6bf9cc4399e5fe
[numthy-code]: https://github.com/ini/numthy/blob/1496df97faec010c909e54e14d6bf9cc4399e5fe/numthy.py
[numthy-readme]: https://github.com/ini/numthy/blob/1496df97faec010c909e54e14d6bf9cc4399e5fe/README.md
[sss]: https://github.com/sbaresearch/smoothsubsumsearch/tree/8dbaf6d39ab88a40380965d25ec2c363d7f27358
[sss-code]: https://github.com/sbaresearch/smoothsubsumsearch/blob/8dbaf6d39ab88a40380965d25ec2c363d7f27358/sssif/sss.py
[sss-matrix]: https://github.com/sbaresearch/smoothsubsumsearch/blob/8dbaf6d39ab88a40380965d25ec2c363d7f27358/sssif/mstep.py
[sss-test]: https://github.com/sbaresearch/smoothsubsumsearch/blob/8dbaf6d39ab88a40380965d25ec2c363d7f27358/test.py
[sss-paper]: https://arxiv.org/html/2301.10529v2
[large-primes]: https://ir.cwi.nl/pub/1367/1367D.pdf
[recent-paper]: https://arxiv.org/pdf/2609.06576v1
[pypy]: https://doc.pypy.org/faq.html
[gmpy]: https://gmpy2.readthedocs.io/en/latest/overview.html

## Matrix optimization follow-up (M25)

The dedicated [GF(2) matrix report](gf2_matrix_research.md) adds an extensive
review of block Wiedemann, block Lanczos, Four Russians/PLE elimination,
packed sparse products, dense heavy-row kernels and stronger filtering.
It inspects eleven papers, comparison slides, three implementer blog posts
and seven pinned implementations;
[retrieval/source hashes](m25_gf2_matrix_sources.json)
preserve the review provenance.

These new experiments belong to the distinct
[P3.8 milestone](TODOS.md#p38--evaluate-gf2-matrix-and-filtering-optimizations).
Existing milestone text/gates and the earlier report remain unchanged.
At this follow-up, P3.1 and P3.2 have separate accepted reference/collector
evidence; P3.3 exact filtering/extraction and P3.4 working SIQS remain open.
P3.8 follows that exact control and measures complete filtering/solve/lifting
costs before promotion. Its optional optimizations can feed later GNFS scaling
without making every challenger a prerequisite for bounded GNFS.

The main caution is mathematical: GF(2) Gram-kernel candidates can fail the
original parity matrix, and random projection/iteration success alone does
not certify rank, a complete basis or a proper factor. The new report and
milestone require original verification, finite retries, bounded checkpoint
storage and held-out end-to-end evidence. M25 changes documentation only;
it closes no implementation or speed gate.
Verification
records preservation and content checks.
