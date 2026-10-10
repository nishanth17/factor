# C1 DLP research reconciliation

Read alongside the earlier [QS/GNFS review](../qs_gnfs_research.md), the
[frozen feasibility protocol](c1_protocol.md), and the C1 result in the
[benchmark guide](../../docs/studies.md). This is selective source review, not a complete
external audit. No external implementation was copied or executed. The local
`results/c1/research/` manifests preserve retrieved bytes, hashes and failures.

## Primary mathematics and economic limits

[Lenstra–Manasse, *Factoring with Two Large Primes*, §§1–3](https://www.researchgate.net/publication/266349860_Factoring_with_Two_Large_Primes)
was read through the author-uploaded text; the DOI/publisher fetch failed.
A partial is an edge from 1 to p; two-prime residuals are edges p–q. Even
incidence cancels large-prime parity in every component. A repeated endpoint
p=p contributes a square immediately. Their fundamental-cycle count is
edges minus vertices plus components; cycles still require factor-base
elimination. Section 3 explicitly charges false reports and cofactor splitting,
warns about small-input losses, and suggests asymmetric prime caps. Its
historical speedups/crossover are specific to the implementation and parameters.

[Boender–te Riele, §§5–7](https://ir.cwi.nl/pub/1367/1367D.pdf)
compares SLP/DLP and limits empirical prediction to similar input sizes and
fixed parameters. Its native SGI/Cray experiments do not predict a PyPy
crossover. Occupancy and collection duration must be observed together.

[Gower–Wagstaff, *Square Form Factorization*, §§3–5](https://homes.cerias.purdue.edu/~ssw/squfof.pdf)
provides an independently described alternative residual splitter, including
square/proper-form handling, bounded queues and multiplier considerations.
Its average-cost analysis makes assumptions; C1 does not implement SQUFOF or
transfer its native timing claims. The existing finite Brent routine is enough
for the declared feasibility screen, which is not a best-splitter search.

## Pinned implementations and licenses

Pins deliberately reuse the prior review. Source/license inspection happened
before choosing adaptations; C1 adapts no external code.

| Source | Inspected behavior and implication | Inspected terms |
| --- | --- | --- |
| [YAFU `8110dfb`, `factor/qs/tdiv.c`](https://github.com/bbuhrow/yafu/blob/8110dfbd8c6f9486d93b1a02de6eb7b180e55a80/factor/qs/tdiv.c), [`filter.c`](https://github.com/bbuhrow/yafu/blob/8110dfbd8c6f9486d93b1a02de6eb7b180e55a80/factor/qs/filter.c) | DLP intake uses separate product/endpoint bounds, base-2 PRP rejection and micro-ECM. Filtering counts components and constructs cycles, then sorts by length. PRP rejection can lose yield; it is not endpoint certification. Resieve and large-prime trial-division files were also inspected. | Public-domain file notices, with upstream acknowledgements; dependencies separate. |
| [Historic msieve `c8727d9`, `mpqs/sieve.c`](https://github.com/radii/msieve/blob/c8727d91305bdbe0972d160ef0ce61dd02ce9193/mpqs/sieve.c) | SLP intake precedes product screening; base-2 rejection precedes SQUFOF. Both endpoints are bounded and squares get explicit handling. Its factor-base-coverage assumptions cannot be inferred from arbitrary v2 `FactorBase` objects. Readme and common filter entry point also inspected. | Public-domain file notice. This is the previously pinned historic mirror, not a current official release. |
| [Yamaquasi `3f95f43`, `src/relations.rs`](https://github.com/remyoudompheng/yamaquasi/blob/3f95f43682ed15d8c1ed206a9a702dd655d7c8ad/src/relations.rs) | Greedy double-partial merging restricts cycles to the SLP-connected component, explicitly as a heuristic. Its comments on empirical rarity do not prove completeness. C1's offline incidence oracle includes disconnected triangles, self-loops and parallel edges. SIQS source and README also inspected. | BSD-3-Clause. |
| [FLINT `a4c9750`, `collect_relations.c`](https://github.com/flintlib/flint/blob/a4c9750d0d3d67bb01cf6d18c187591b313451c3/src/qsieve/collect_relations.c), [official documentation at the same pin](https://github.com/flintlib/flint/blob/a4c9750d0d3d67bb01cf6d18c187591b313451c3/doc/source/qsieve.rst) | Intake bounds small cofactors and checks multiplier GCD. The documentation describes equal-cofactor partial merging and singleton removal. No explicit intake primality call was found in that path; this does not extend v2's prime-residual promise. | LGPL-3.0-or-later QS file notice. |
| [JavaMath `088d01f`, `TDiv_QS_2LP`](https://github.com/TilmanNeumann/java-math-library/blob/088d01fa97e7d0412c6cbbb8f463fa7e5f78ab97/src/main/java/de/tilman_neumann/jml/factor/siqs/tdiv/TDiv_QS_2LP.java) | Size-based Hart/tiny-ECM/rho/QS splitting, probable-prime checks, and `Smooth1LargeSquare` distinguish repeated factors. Floating size gates and probable-prime semantics do not transfer to exact proven-residual admission. | GPL-3.0-or-later. |
| [SymPy `2f22a5f`, `ntheory/qs.py`](https://github.com/sympy/sympy/blob/2f22a5f81e2f4124380be3739a092e9ff20128de/sympy/ntheory/qs.py) | Python SLP matching provides a useful comparison of exact recovery and classification costs. Optional native arithmetic prevents treating its runtime as pure PyPy evidence. | BSD-3-Clause core, with separate component notices. |

Factor has no explicit top-level license in this checkout. Future copying must
check compatibility and notices; the mathematical ideas above do not authorize
copying GPL/LGPL implementations into the project. C1 reuses existing v2 routines and independently written graph, provenance
and diagnostic code; no upstream implementation was adapted.

## Release cross-check and current research

The official [msieve 1.53 source release](https://sourceforge.net/projects/msieve/files/msieve/Msieve%20v1.53/msieve153_src.tar.gz/download)
was downloaded and matched the publisher's SHA-256
`c5fcbaaff266a43aa8bca55239d5b087d3e3f138d1a95d75b776c04ce4d93bb4`.
Its public-domain `mpqs/sieve.c` enables DLP at 282 bits, approximately 85
digits, and uses the SLP limit raised to 1.8 instead of squaring it. The
comment explicitly attributes this to rare, poorly surviving large-large
pairs. Its residual path uses base-2 rejection then SQUFOF; `relation.c`
initializes roots for every component, but its cycle enumerator skips paths
longer than 100 edges per side. That finite heuristic is not a completeness
guarantee. C1's offline oracle instead checks the full E−V+C dimension and
includes a test for a 2,001-edge disconnected cycle. `Readme.qs`
describes initially slow then accelerating combined-relation accumulation. These inspected release sources
cross-check the historic mirror without implying that 85 digits is a PyPy
cutoff. Release source files were read, not compiled or executed.

The October 2026 [CUDA-MPQS preprint, v1](https://arxiv.org/html/2610.07126v1)
reports a GPU SIQS pipeline using SLP, explicitly leaving multiple-large-prime
variants unimplemented. Its useful-dependency failures and truncated-tuning
limitations reinforce the need to measure extraction and adequate collection
duration. Those are implementation-specific observations; its reported
small-input SLP failure is not a universal theorem and does not describe v2's
verified completing SLP controls. GPU timing is not evidence of PyPy speed.
The paper lists LGPL-3.0-only code with a CUDA exception; the attempted direct
license-file fetch failed, so no code adaptation or independent license audit
of that project is claimed. The primary-source abstract and implementation
sections were inspected; the reported record was not independently reproduced.

## Blogs and technical discussions

[Programming Praxis, SQUFOF continued-fraction exercise](https://programmingpraxis.com/2014/07/08/squfof-continued-fraction-version/)
led to the Gower–Wagstaff primary paper. Its erratum discussion reinforces the
need for independent exact tests; its code was not adapted.
[FLINT developer discussion, GSoC 2015](https://groups.google.com/g/flint-devel/c/qDqu17-J_O8)
was readable on this pass. The maintainers describe parameter, storage,
cofactoring and duplicate-removal costs; recommendations in that old native
context are hypotheses, not instructions or v2 performance evidence. Mersenne
Forum searches did not yield a readable relevant thread; no forum guarantee
is claimed. The review therefore grounds arithmetic guarantees in the primary
papers and pinned source, not excerpts or empirical comments.

## Implications for the user-authorized follow-up

The first screen used product cap 10,000 B², although the pinned Yamaquasi
`double_large_factor` uses about 100 B² near 200 bits. Its comments explain
that endpoint limits and product limits serve different purposes: expanding
the product increases trial-division and cofactoring pressure. The follow-up
therefore tests nested 64/128 B² products with endpoints capped at 100B.
This is an independently chosen finite experiment, not a copied parameter
selector or an assumption of equal performance.

Yamaquasi enables DLP by default only above 256 bits (roughly 77 decimal
digits); explicit preferences can enable it earlier. It also halves factor
base cardinality in its DLP preset. This supports the user's point that
small-input losses cannot settle larger-input usefulness. Native cutoffs do
not transfer to PyPy. Our 40/50/60-digit complete-pipeline feasibility cases
cover larger workloads than the original completing 30/40-digit controls;
they cannot establish a universal 80–100-digit crossover. No conclusion about
that uncalibrated range follows from a failure below it. Parameter selection
and complete-factor confirmation remain distinct from mathematical validity.

Lenstra–Manasse's delayed cycle accumulation is why the new protocol provides
longer finite collection windows rather than requiring early endpoint
collisions to earn an extension. The first screen's dense oracle refused
40-digit records before measuring their complete cycle space. A bounded
offline spanning forest removes that diagnostic artifact and is checked
against independently implemented generic GF(2) elimination. It does not
weaken production reservations or silently change accepted R3 provenance.

## Contracts preserved during feasibility

During the diagnostic studies, the SLP collector, admitted residual type,
R3 stable mixed order, complete immutable cache identities, original-row
lifting, exact square corrections, charged checkpoint replay and finite
stores remained byte-for-byte unchanged. Later opt-in implementation and its
separate acceptance are recorded in [the implementation report](c1_implementation_results.md).
Diagnostic records are not `AtomicRelation` objects and never enter a runtime
checkpoint. Full trial division on randomly chosen positions cross-checks the
root-guided diagnostic recovery, including rejected candidates. Offline
kernels are checked against the original incidence, then exact exponent sums
and both GCD signs. Diagnostic completeness alone does not establish a production graph's
eviction, retained ownership, cancellation or resume acceptance. The later
implementation report records those separate passed gates and explicitly
bounds production cycle length.
