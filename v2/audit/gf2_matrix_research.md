# GF(2) matrix optimization research (M25)

Date: 3 October 2026. Scope: exact dependency solving and filtering for QS,
SIQS and reusable future GNFS infrastructure on PyPy Python 3.11.

## Findings and decision

Block Wiedemann is worth evaluating, alongside block Lanczos and dense
Four Russians elimination. The useful choice depends on the filtered matrix,
available memory and execution topology. Create the distinct
[P3.8 milestone][p38]; preserve every existing milestone, including P3.3's
reference implementation and M24 collector carry-over. This review neither
implements a solver nor closes an implementation or performance gate.

The recommended experiment order is: measure the reference pipeline; improve
packed matrix products and dense/hybrid elimination; evaluate stronger
filtering; compare block Lanczos and block Wiedemann when sparse solving
is justified. These are project recommendations inferred from the sources,
not measured PyPy improvements. Native implementations demonstrate mature
algorithm choices; their word sizes, thresholds and timings are not Factor
defaults. There is no defensible universal "best solver" ranking here.

| Candidate | Why investigate it | Priority and limit |
| --- | --- | --- |
| Python-integer bitsets, echelon/free-variable recovery | Retain an exact control; whole-row XOR avoids per-entry Python arithmetic | First, after P3.3; compare recovery storage as well as elimination |
| Four Russians/M4RI-style tables; PLE/PLUQ | Tabulate combinations of pivot rows and reuse bulk XOR; rank-aware decomposition can expose free variables | First dense challenger; table creation and memory must pay for themselves |
| Packed sparse matrix times block; dense heavy-row split | Reuse each index traversal for several independent GF(2) vectors; adapt layout to uneven QS density | Shared prerequisite for iterative solvers; include transpose and layout construction |
| Bounded higher-weight merges, fill-aware pivot batches, components | Reduce dimensions without letting nonzeros and lifting history overwhelm the solve | Separate follow-up to P3.3's singleton/weight-two control |
| Block Lanczos | Established QS/NFS method with a short recurrence and comparatively little sequence storage | Sparse challenger for one machine; singular block handling and final correction are essential |
| Block Wiedemann | Split Krylov generation, polynomial generation and reconstruction; independent sequence blocks can suit separated resources | Sparse challenger; generator, reconstruction, checkpoints and I/O count |
| Faster matrix-polynomial generator | Thomé-style divide-and-conquer can reduce the generator stage's asymptotic cost | Only after a correct base-case generator and measured generator dominance |
| GPU, MPI and native exact libraries | Useful design references for layouts, checking and scheduling | Reference evidence here; no new native solver/backend or hardware requirement |

## Reconciliation with completed and planned work

| Existing work | Evidence/status at this review | P3.8 boundary |
| --- | --- | --- |
| P3.1 | [M20][m20]: exact polynomial/relation contracts and exhaustive collector accepted | Reuse exact atomic relations and the verifier |
| P3.2 | [M22][m22]: bounded single-large-prime collector accepted; its small-window controls do not establish a production speedup | No collector rewrite or changed acceptance status |
| P3.3 | Open in the backlog: filtering/extraction/pipeline source work is in progress in the shared workspace; [M24][m24] adds the separate collector diagnosis/resieving carry-over | Accept and measure this control first; P3.8 owns new matrix optimization experiments |
| P3.4 | Open: SIQS families, bounded dispatch and serialized integration | Integrated P3.8 promotion requires working SIQS; matrix-only prototypes can follow P3.3 |
| P5.4 | Open: double-large-prime SIQS and graph-cycle filtering | Preserve its collector scope; evaluate matrix methods on its fixtures later if available |
| P6.2 and P7.6 | Existing conditional sparse solver investigations/scaling work | Reuse P3.8 findings and kernels when eligible; preserve their text and gates |
| Phase 8 | Existing future Phase 2 optimization follow-up | No scope or status changes |

At M25's preflight, accepted evidence reaches the exact relation contracts
and bounded collector. Filtering, linear-algebra, extraction and pipeline
sources are being developed in the shared workspace; their P3.3 acceptance
and matrix performance evidence are not yet recorded in the backlog. M25
records observed source changes without editing or accepting that work.
Accordingly, these sources justify experiments, not a claim that matrices
currently dominate Factor. P3.8 is an optional
optimization milestone; completing all its challengers is not a prerequisite
for the committed bounded GNFS work. P3.3/P3.4 retain their existing ownership.

## Literature and implementer blogs

The [source manifest][manifest] records URLs, resolved URLs, fetched sizes,
SHA-256 hashes, revision dates and review scopes. Eleven papers were read for
the relevant algorithms/contracts, plus Kruppa's comparison slides and three
direct implementer blog posts. This is targeted source review, not a claim of
complete proof audits or independent reproduction of external benchmarks.

| Primary source | Relevant result and qualification |
| --- | --- |
| [Wiedemann, 1986][wiedemann] | Sparse black-box linear algebra via projected sequences and minimal polynomials; useful scalar reference, with randomized failure handling |
| [Coppersmith, 1994][coppersmith] | Block Wiedemann over GF(2), including the matrix-sequence generator and several dependency candidates; historical 32-bit packing is an implementation setting |
| [Montgomery, 1995][montgomery] | Block Lanczos for factoring dependencies; handles self-orthogonality through selected nonsingular subblocks, not the ordinary real-valued Lanczos recurrence |
| [Kaltofen and Saunders, 1991][kaltofen] | Black-box solving, rank and preconditioning analysis; hypotheses for these randomized tasks must be checked rather than inferred from a returned kernel vector |
| [Thomé, 2002][thome] | Subquadratic matrix-sequence generator; reported large-characteristic experiments do not establish a GF(2)/PyPy crossover |
| [Albrecht, Bard and Hart, arXiv:0811.1714][m4rm] | Dense M4RM/Strassen-Winograd hierarchy and data-locality costs; distinguishes multiplication from the elimination task |
| [Albrecht, Bard and Pernet, 2011][ple] | Rank-aware PLE and Four Russians elimination, Gray-code tables and cache effects; relatively sparse matrices still challenge dense storage |
| [Bertolazzi and Rimoldi, 2012][bertolazzi] | Dense rectangular decomposition avoiding column exchanges; a useful pivot/layout challenger, not proof that removing swaps explains every speed gain |
| [Bouvier, 2013, HAL v3][bouvier] | Purge/clique weighting, higher-way merges and minimum-weight spanning trees; reducing dimensions can increase matrix weight |
| [Bouillaguet and Zimmermann, 2021][sge] | Structured elimination using low-cost independent pivot batches and bounded passes; native OpenMP scaling is separate from a serial PyPy adaptation |
| [Valenta et al., 2016][valenta] | Complete factoring trade-offs between filtering, oversieving and linear algebra; historical EC2 comparison favored Msieve over CADO on their configuration |
| [Kruppa, WCNT 2011 slides][kruppa] | Block Lanczos versus Wiedemann, block-size overhead, homogeneous recurrence and low-rank auxiliary updates; workloads/hardware determine the useful choice |
| [Albrecht, M4RI 20121224][blog-tuning] | Maintainer describes retuning base cases and memory behavior after a competing dense decomposition; profile implementation costs before crediting an algorithm label |
| [Albrecht, Gröbner-basis GF(2) M4RI post][blog-hybrid] | Specialized sparse/dense structure can change the winner; Gröbner matrices have different structure from QS, so transfer requires a measured fixture |
| [Albrecht, Bertolazzi code post][blog-lu] | Implementer follow-up on the competing LU code and retuned M4RI; historical results are context, not a current library ranking |

The fetched M4RM manuscript includes a 2013 revision date; its arXiv identifier
originates in 2008. Bouvier's retrieved HAL v3 is 28 PDF pages including archive
material, despite the older author page's 22-page description. The published
SGE paper is cited as 2021, consistent with its authors' publication records;
the current CADO source header labels the reference 2020. These date/page
differences are preserved rather than treated as new performance evidence.

## Pinned implementation evidence

Revisions were fetched on this review date. CADO's GitHub revision was also
confirmed against its Inria GitLab API. The Msieve, Yamaquasi and FLINT revisions
reuse the prior QS research pins so the follow-up does not silently change
comparison sources. All downloaded upstream files were inspected as data;
none was imported, compiled, installed or executed.

| Implementation | Exact revision | Inspected matrix design |
| --- | --- | --- |
| [CADO-NFS][cado] | `692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b` | BWC staged/checkpointed pipeline, dense/sparse bucket and sliced products, purge/merge/replay; commit 2 October 2026 |
| [Msieve][msieve] | `c8727d91305bdbe0972d160ef0ce61dd02ce9193` | `common/lanczos/`: preprocessing, dense-row handling, reordering, sparse products, retry and original-matrix verification |
| [Yamaquasi][yama] | `3f95f43682ed15d8c1ed206a9a702dd655d7c8ad` | `gf2.rs`: dense Gaussian control, packed block Lanczos, dense heavy rows plus tiled coordinates, final kernel correction |
| [FLINT][flint] | `00cc19b350b5302876b4fa8c877df5a676c2061c` | QS block Lanczos, sparse reduction and final dependency checks; retained 2 October 2026 pin |
| [M4RI][m4ri] | `35f06e132363d12592cdd94c98473d36e736cef9` | Four Russians/Gray-code elimination, PLE/PLUQ and triangular kernel recovery; commit 25 August 2026 |
| [Wouters blanczos][blanczos] | `232a57090e07e95a436ae03f4e6536e30d6f727d` | Density-dependent cache blocks with relative coordinates in the active source; alternative CSR class; native 64-vector solver; commit 15 July 2025 |
| [Januszewski block-wiedemann][gpu] | `509fb79bdeb231bce4d0e22c36b814d8d3bbd390` | CUDA sparse-block autotuning, three-stage solver, CPU/Python references and saved sequences; commit 27 July 2026 |

The new GPU solver is linked from its author's [Paderborn software page][gpu-home].
Its README says the production generator is currently block Berlekamp-Massey;
the Python reference contains a divide-and-conquer structure but uses NumPy
and is described as slow. This is useful recent design evidence, not a
validated pure-Python accelerator or an independently verified SOTA claim.
Its README also loosely describes rectangular input while the Python class
documents square input: a Factor adapter must define its square operator and
kernel mapping explicitly. No CUDA path is added to the supported roadmap.

Wouters' README describes per-block CSR, but the pinned solver includes
`sparse.hpp`, whose active multiplication loops use blocked relative-coordinate
pairs. A separate `csr.hpp` implementation is present. This review follows
the inspected call path rather than assigning the README's layout to it.

Wouters' `Ncol >= Nrow + 64` and at-most-64 output vectors are wrapper
restrictions, not mathematical requirements of every block Lanczos solver.
Yamaquasi uses system randomness; Msieve's outer recovery can retry until a
dependency set is obtained. Factor needs explicit seeds and finite allowances
instead of copying those execution policies. Inspected headers/licenses
include Msieve's public-domain notice, Yamaquasi's BSD-style notice, Wouters'
BSD-3-Clause, M4RI GPL-2.0-or-later, the CADO files' LGPL-2.1-or-later, FLINT
LGPL-3.0-or-later and the GPU project's LGPL-3.0-only/CUDA exception. Review
the exact file notices before any future source adaptation.

## Exact matrix and lifting contracts

For this report define `M` with **rows = sign/prime parity constraints** and
**columns = verified relations**. A dependency is a nonzero binary vector `d`
with `M d = 0`. P3.3 can store the transpose, but each adapter must declare
orientation, dimensions, permutations and the meaning of output indices.
Repeated sparse coordinates cancel modulo two; blindly deduplicating them
with a set gives the wrong coefficient when multiplicity is even.

Lift every reduced dependency through filtering and permutations to original
atomic relations. Preserve full exponents and square corrections, including
dependencies created by empty reduced relations and distinct equal-parity
relations. Recheck original parity, even exponent totals, `X² ≡ Y² (mod n)`
and each proper GCD split. A valid kernel vector can give a trivial congruence;
retain bounded attempts at several distinct vectors/combinations and collect
more useful relations when the existing pipeline's recovery policy permits.

**Normal-equation caveat:** `ker(M)` is contained in `ker(Mᵀ M)`, but equality
fails over GF(2). For `M = [[1], [1]]`, `Mᵀ M = [[0]]`, although `M [1]` is
nonzero. This direct mathematical counterexample requires correction of an
iterative candidate block, not merely a Gram-kernel test. Yamaquasi explicitly
forms `M Y`, computes its small kernel `K`, and returns nonzero columns of
`Y K`; [Msieve][msieve-lanczos] and [FLINT][flint-lanczos] also correct and
check candidates against their original operator.

Apply `Mᵀ(M V)` implicitly when a method needs a symmetric operator; explicitly
forming the Gram matrix can destroy sparsity. Block Wiedemann can instead
use a documented square padding/completion of `M`, retaining permutations and
discarding artificial-coordinate solutions. The symmetric embedding
`[[0, M], [Mᵀ, 0]]` is another exact option: its relation-coordinate component
must be nonzero and satisfy `M d = 0`, and its larger dimension/work must be
charged. These embeddings are experiments with explicit lifting checks,
not interchangeable unverified adapters.

An iterative solver may return a verified independent batch without finding
the entire nullspace. State that contract, report bounded failures separately,
and do not certify exact rank/nullity or factorization success from the number
of vectors found. Structural matching gives a sparsity-based rank bound;
actual GF(2) cancellations still require arithmetic checks.

## Dense, hybrid and filtering experiments

### Dense Four Russians and kernel recovery

For a stripe with `k` independent pivot rows, build the `2**k` XOR combinations
using Gray-code updates, then clear the stripe in other rows by a table lookup
and bulk XOR. The [PLE paper][ple] and [M4RI elimination source][m4ri-echelon]
provide the contracts and rank-deficient cases. A Python-integer adaptation
is plausible because an XOR covers many bits, but table construction, bigint
allocation and pivot extraction could dominate. Measure them on supported
PyPy; native SSE/AVX and cache constants do not transfer automatically.

Compare existing dependency masks with echelon/PLE-style free-variable
back-substitution, emitting a bounded batch of solutions. The
[M4RI kernel source][m4ri-solve] recovers free-variable vectors through a
triangular solve and inverse permutation. For right-kernel solving, elementary
equation-row operations do not require carrying a full identity matrix just
to identify relation columns. Filtering still needs its separate lifting
history. Avoid allocating a complete basis when the contract asks for a
bounded dependency batch. Count recovery memory and useful factors too.

Sweep small `k` values under an explicit workspace cap; include tables for
row tails and any provenance masks, not just the pivot stripe. Handle missing
pivots, short last stripes, rank deficiency and empty matrices exactly.
Consider sparse elimination followed by a dense trailing core only when
observed fill/density and total conversion cost justify it. Strassen-style
recursion is a later dense-kernel option, conditional on sufficiently large
dense blocks and a demonstrated base-case bottleneck.

### Stronger filtering without losing dependencies

Keep P3.3 singleton-only/weight-two reductions as the control. New P3.8 arms
can compare bounded higher-weight merges, minimum-weight spanning-tree merge
choices and fill-aware independent pivot batches. [Bouvier][bouvier] explains
why higher-way merges trade fewer dimensions for more nonzeros;
[the published SGE method][sge] and [current CADO merge][cado-merge] provide
batch/threshold designs. Recompute affected costs, cap merge degree/fill and
record reversible maps; ordinary numerical pivot tolerance is irrelevant.

Distinguish two graph ideas: disconnected components of the full constraint/
relation bipartite graph can be solved independently and lifted exactly;
filtering's weight-two "cliques" are components of a different graph.
Deleting a surplus clique intentionally discards candidate dependencies.
Such pruning is a separate conditional arm that must retain adequate usable
excess and measure complete factor recovery; it is not an exact reduction
preserving every original dependency. A smaller matrix alone is insufficient.

Compare `dimension * nonzeros` as a preliminary iterative-work proxy with
actual solve, lifting and factor time. Include oversieving/repeated filtering
and allocation costs; do not tune an average density in isolation. Valenta's
[historical factoring study][valenta] illustrates those whole-pipeline
trade-offs. Its 350-minute CADO versus 140-minute Msieve result used the same
53-million-relation set on their EC2 configuration, but different preprocessing
and solver stacks; it is not an algorithm-only or current-version comparison.

## Sparse products, block Lanczos and block Wiedemann

### Shared packed kernels

For a block `V` of `b` binary vectors, store `V[j]` as `b` independent bit
lanes. Compute `(M V)[i]` by XORing `V[j]` for each nonzero in row `i`.
This packs GF(2) vectors; it does not turn them into GF(2**b) scalars.
Implement the transpose consistently. Verify products with an independent
entry-wise GF(2) oracle and the bilinear identity
`uᵀ(M v) = (Mᵀ u)ᵀv`, including duplicate cancellation and tails.

Compare list/array indices, CSR/column storage and tiled coordinates, with a
separate dense heavy-row block where useful. [Yamaquasi][yama-gf2] provides a
particularly QS-specific design; [CADO buckets][cado-buckets] and
[Wouters][blanczos] demonstrate density-dependent layouts. Delta indices or
compressed gaps are conditional ideas: reduced storage can lose to Python
decoding and temporary allocation. Build and charge both orientations only
if their savings justify their aggregate memory and conversion cost.

Sweep block widths around representative 16/32/64/128-bit choices and relevant
PyPy integer representation boundaries, rather than inheriting a C word width.
Count XORs, index visits, auxiliary dense products, initial layout creation,
resident bytes, scratch storage and complete solve time. The GPU source's
[per-block autotuner][gpu-tuner] is a useful example of heterogeneous choices,
but Factor should freeze tuning before held-out evaluation and charge tuning
work. No kernel or device claim here has been reproduced.

### Block Lanczos challenger

Implement the actual finite-field block recurrence from [Montgomery][montgomery]
with rank-revealing selection/partial inverses for singular Gram blocks,
finite seeded restarts and original-matrix correction. Keep a small exact
reference before optimizing the packed recurrence. Singular/self-orthogonal
blocks are expected cases in characteristic two, not permission to divide by
zero, silently lose columns or port a floating-point Lanczos routine.

After the recurrence is correct, [Kruppa][kruppa] identifies two further
experiments: the homogeneous formulation can avoid a solution accumulator,
and a low-rank factorization of an auxiliary update can reduce block work.
Neither optimization changes the required original-kernel check. Compare
auxiliary operations and retained vectors as well as sparse iterations; wider
blocks reduce iteration counts while increasing dense auxiliary costs.

### Block Wiedemann challenger

Define a square operator `B` and its exact map to `M`. With seeded projection
blocks `X` and `Y`, generate the matrix sequence `S_i = Xᵀ B**i Y`; compute a
matrix-polynomial generator, reconstruct candidate vectors, then lift and
verify them. An entry-by-entry scalar Berlekamp-Massey computation is not the
block generator. Validate the generator against extra unused sequence terms
as well as the final operator check. [Coppersmith][coppersmith] supplies the
block algorithm; [CADO BWC][cado-bwc] supplies a concrete stage/checkpoint
decomposition: prep/secure, Krylov, collection, lingen, mksol and gather.

The approximate sequence length `N/m + N/n` for operator dimension `N` and
projection widths `m,n` is a planning model with overshoot/recovery, not a
guaranteed stop for every GF(2) matrix. Measure projection failures, missing
generators, zero outputs, dependent vectors and exhausted allowances. Bound
sequence extension and restarts, preserving the same total time/CPU/work
budget; exhausted random projections do not establish full rank or primality.

Include generator and reconstruction cost, not just Krylov products. For
square width `b`, a sequence of order `N/b` terms has order `N*b` bits of
matrix coefficients before checkpoints/object overhead. Larger blocks can
therefore increase memory/I/O and generator work. After a validated base-case
generator, compare [Thomé-style divide-and-conquer][thome] only if that stage
dominates. A faster generator does not imply a faster complete factorization.

Independent right-projection sequence blocks and reconstruction partitions
can reduce communication across separated machines; the recurrence within a
given sequence is still sequential. [Kruppa][kruppa] supports that topology
distinction, while the [EC2 comparison][valenta] shows it need not win on a
particular cluster. First establish a serial control on PyPy; any process
follow-up must compare identical assignments and total resource budgets.

Checkpoints must identify the matrix/atomic-relation hashes, field,
orientation, permutations/lifting history, operator, projection seeds/blocks,
widths, sequence extent, generator/reconstruction stage and consumed budget.
Retain the Krylov states reconstruction actually needs; [CADO's retention
documentation][cado-bwc] describes combinations that otherwise break mksol.
Detect stale, truncated and corrupted state before continuing. Cap disk,
RAM and retained checkpoints; mathematical spot checks supplement hashes and
the mandatory final original-matrix verification.

## Acceptance, measurements and deferrals

P3.8 records separate adopt/defer/reject decisions for its tracks. Research
completion is M25; all new implementation checkboxes remain open. Establish
the P3.3 exact control, a matrix corpus and at least one justified bounded
challenger before claiming an implemented optimization milestone. Defer
larger solvers with explicit matrix-cost/memory evidence when scaling does
not yet justify them. No mandatory implementation of every candidate or
automatic default change follows from this report.

Use hand-constructed and independently solved tiny matrices: zero/empty/full
rank, rectangular, repeated coordinates, equal-parity relations, singleton
cascades, disconnected components, singular Gram blocks, repeated invariant
factors, unlucky projections and the normal-equation counterexample. Check
known rank/nullity for exact dense solvers; require sound independent batches
and bounded failure outcomes for iterative solvers. Inject faulty lifting,
bad permutations, bit flips, truncated sequences and stale checkpoints.

Freeze genuine QS/SIQS matrices after P3.3/P3.4 and retain hashes, seeds,
dimensions, nonzeros, heavy-row weights, filtered excess and reduction history.
Synthetic sparse matrices test scaling and failure behavior; they cannot by
themselves demonstrate a SIQS crossover or justify a production default.
Measure filtering, representation conversion, products, generator/recurrence,
reconstruction, original verification, lifting and GCD extraction separately,
then measure complete factorizations with matched collectors/inputs/budgets.

Use the same supported PyPy release and language version for old/new arms,
default JIT, several seconds of validated warmup and stable repeated samples
(at least the existing three-second/nine-sample runner baseline). Separate
cold startup/JIT from warmed execution. Record median/spread/confidence
intervals, successful completion, exhausted budgets, peak/aggregate RSS,
CPU time and disk/I/O. Tune on training matrices, freeze parameters/selection,
then use held-out matrices and complete-factorization inputs. Apply the
unchanged [end-to-end promotion policy][promotion]. Report matrix capacity
gains separately when the dense control cannot fit; a failed baseline has
no successful timing ratio.

Defer extension-field preconditioning until small-field failures and its
extra arithmetic/lifting contract justify it. Do not substitute ordinary
NumPy/SciPy floating-point solve, eigensolvers, SVD or real-valued Lanczos for
GF(2) dependencies. Do not form sparse matrix powers or apply dense Strassen
indiscriminately. Native M4RI/LinBox/GPU acceleration would be a separate
backend decision; this milestone keeps Factor's algorithms in Python on
PyPy and uses upstream projects as design evidence.

## Verification and remaining uncertainty

[M25 content/provenance verification][verification] checks insertion-only
backlog changes, preservation of existing report/log prefixes and milestone
text, local links/anchors, fetched-source hashes, before/current runtime/support
hashes and historical evidence. Concurrent QS implementation/support changes
are recorded separately; M25 edits documentation only. Temporary source
downloads and extraction files are digest-checked and removed after review;
the manifest retains retrieval and review metadata. No runtime tests or
performance measurements apply to M25.

Unknown until measured: Factor's matrix bottleneck, the best layout/table/
block sizes, failure rates, useful merge stopping rule, dense/sparse crossover,
solver capacity, complete-factor benefit and valuable parallel topology.
Next: finish P3.3/P3.4's existing pipeline and use P3.8's distinct experiments
when its prerequisites and observed matrix costs make them worthwhile.

[p38]: TODOS.md#p38--evaluate-gf2-matrix-and-filtering-optimizations
[promotion]: TODOS.md#how-to-use-the-gates
[m20]: phase_three_m20_p31_summary.json
[m22]: phase_three_m22_p32_summary.json
[m24]: m24_p32_carryover_verification.json
[manifest]: m25_gf2_matrix_sources.json
[verification]: m25_gf2_matrix_verification.json
[wiedemann]: https://www.csd.uwo.ca/~mmorenom/CS424/Ressources/WIEDEMANN-IEEE-1986.pdf
[coppersmith]: https://www.ams.org/journals/mcom/1994-62-205/S0025-5718-1994-1192970-7/S0025-5718-1994-1192970-7.pdf
[montgomery]: https://link.springer.com/content/pdf/10.1007%2F3-540-49264-X_9.pdf
[kaltofen]: https://users.cs.duke.edu/~elk27/bibliography/91/KaSa91.pdf
[thome]: https://members.loria.fr/EThome/files/jsc.pdf
[m4rm]: https://arxiv.org/pdf/0811.1714
[ple]: https://arxiv.org/pdf/1111.6549
[bertolazzi]: https://arxiv.org/pdf/1209.5198
[bouvier]: https://inria.hal.science/hal-00734654v3/document
[sge]: https://journals.flvc.org/mathcryptology/article/download/126033/127787
[valenta]: https://ifca.ai/fc16/preproceedings/19_Valenta.pdf
[kruppa]: https://event.cwi.nl/wcnt2011/slides/kruppa.pdf
[blog-tuning]: https://martinralbrecht.wordpress.com/2012/12/21/m4ri-20121224/
[blog-hybrid]: https://martinralbrecht.wordpress.com/2012/06/24/linear-algebra-for-grobner-bases-over-gf2-m4ri/
[blog-lu]: https://martinralbrecht.wordpress.com/2013/08/17/enrico-bertolazzis-linear-algebra-code-over-gf2-available/
[cado]: https://github.com/cado-nfs/cado-nfs/tree/692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b
[cado-bwc]: https://github.com/cado-nfs/cado-nfs/blob/692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b/linalg/bwc/README
[cado-buckets]: https://github.com/cado-nfs/cado-nfs/blob/692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b/linalg/bwc/matmul-bucket.cpp
[cado-merge]: https://github.com/cado-nfs/cado-nfs/blob/692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b/filter/merge.cpp
[msieve]: https://github.com/radii/msieve/tree/c8727d91305bdbe0972d160ef0ce61dd02ce9193
[msieve-lanczos]: https://github.com/radii/msieve/blob/c8727d91305bdbe0972d160ef0ce61dd02ce9193/common/lanczos/lanczos.c
[yama]: https://github.com/remyoudompheng/yamaquasi/tree/3f95f43682ed15d8c1ed206a9a702dd655d7c8ad
[yama-gf2]: https://github.com/remyoudompheng/yamaquasi/blob/3f95f43682ed15d8c1ed206a9a702dd655d7c8ad/src/matrix/gf2.rs
[flint]: https://github.com/flintlib/flint/tree/00cc19b350b5302876b4fa8c877df5a676c2061c
[flint-lanczos]: https://github.com/flintlib/flint/blob/00cc19b350b5302876b4fa8c877df5a676c2061c/src/qsieve/block_lanczos.c
[m4ri]: https://github.com/malb/m4ri/tree/35f06e132363d12592cdd94c98473d36e736cef9
[m4ri-echelon]: https://github.com/malb/m4ri/blob/35f06e132363d12592cdd94c98473d36e736cef9/m4ri/echelonform.c
[m4ri-solve]: https://github.com/malb/m4ri/blob/35f06e132363d12592cdd94c98473d36e736cef9/m4ri/solve.c
[blanczos]: https://github.com/SebWouters/blanczos/tree/232a57090e07e95a436ae03f4e6536e30d6f727d
[gpu]: https://github.com/drjanosch42/block-wiedemann/tree/509fb79bdeb231bce4d0e22c36b814d8d3bbd390
[gpu-home]: https://math.uni-paderborn.de/ag/arbeitsgruppe-algebra-und-zahlentheorie/forschung/software
[gpu-tuner]: https://github.com/drjanosch42/block-wiedemann/blob/509fb79bdeb231bce4d0e22c36b814d8d3bbd390/cuda_spmm/src/autotuner.cpp
