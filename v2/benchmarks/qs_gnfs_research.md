# QS and GNFS research reconciliation — 9 October 2026

This review inspected 11 immutable repository pins, 76 selected source/license
files and primary papers against A10 commit `6794af6`. No external code was
adapted or executed, and no new performance measurements were made. The
arithmetic examples below were checked independently on PyPy 7.3.23 implementing
Python 3.11.15. Inspection was selective, not a complete correctness audit.

A10 establishes verified wider primality ranges and a bounded candidate
decision. It does not establish that current v2 is optimal. The following
experiments remain conditional work under the existing [roadmap](../ROADMAP.md)
owners; this research closes no implementation or performance gate.

## Candidate decisions and owners

| Owner | Candidate and trigger | Required comparison and correctness contract |
| --- | --- | --- |
| C1 / E3, coordinated relation ownership | Reduce duplicate residual classification if its complete cost matters. FLINT's intake motivates checked factor-base coverage or opaque cofactor matching. | First prove all possible primes below B were excluded before inferring primality for `1 < r < B**2`. Alternatively, introduce a separate exactly verified opaque-cofactor type: equal composite cofactors also pair into a square. Preserve independent relation verification, terminal certainty, finite occupancy and versioned checkpoints; current v2 promises prime residuals. |
| C1: R4 + P5.4 | Double-large-prime (DLP) collection if accepted B1 calibration still has insufficient useful yield. | Freeze calibrated single-large-prime (SLP) control; charge splitting, separate prime/product caps, stores, filtering and extraction. Use a full graph-cycle oracle, including disconnected components, repeated primes, self-loops, duplicates and eviction. Relation counts alone cannot establish a win. |
| A8 → B6; B7 / C4 separately | Four Russians table elimination on genuine costly post-filter matrices. | Python-int XOR already runs natively. Train bounded widths such as 4/6/8; include `2**k` rows, provenance width, pivot search, conversion, recovery and simultaneous storage. Test rectangular/rank-deficient matrices and verify every dependency against the original operator. Storage/provenance and stronger filtering retain their own owners. |
| E6: P8.6 | Bounded score-table reuse and bulk threshold-mask extraction if a changed workload exposes allocation/scanning cost. | v2 already has bytearray slice translation, Gray root reuse, root-guided division, buckets/resieving and prime-power plans. Preserve score saturation/clipped thresholds, exact refinement, position charges, cancellation and resume cursors. Existing cache losses stand until the declared trigger changes. |
| A7 / R5; acceptance in E1 | Refresh existing bounded SSS/SSSf comparisons after calibration. | Reuse repaired adapters and forced-factor/checkpoint contracts. Compare complete factors and certainty under matched budgets; account for loss, conversions, startup and recovery. The original upstream comparator is not current calibrated SIQS. |
| F5: P6.2 | Bradford–Monagan–Percival `A=A0*q` reuse if polynomial/root setup dominates. | Carry q's full exponent and partial-relation/provenance role when it lies above the factor base; bound duplicates and setup. The method changes relation semantics. Its score-subtraction shortcut is not a proof with v2's saturated/conservative scores. |
| F6 / F7: P6.2 | General multi-LP incidence after viable DLP economics; batch smooth parts only after a changed size/backend makes candidate division dominant. | Three large primes need general GF(2) incidence, not a two-endpoint graph. Product/remainder trees need complete valuation recovery and bounded nodes/storage/latency. Existing small-workload batching losses remain the control. |
| A9 onward; F1 after scaled integration | GNFS remains later coverage work after H1's coverage review. | Establish field/norm/ideal/character and algebraic-root contracts before scaling. Measure a crossover only after both calibrated SIQS and a correct scaled GNFS engine exist. Native support thresholds and asymptotic estimates do not establish a PyPy crossover. |

## External implementations and licenses

Each link pins the inspected commit. Native implementations supply design ideas
and contextual references; their runtime thresholds are not Python defaults.
License observations concern inspected files, not every bundled dependency.

| Project / immutable pin | Relevant implementation evidence and transfer limits | License observations |
| --- | --- | --- |
| [YAFU, `8110dfbd8c6f9486d93b1a02de6eb7b180e55a80`](https://github.com/bbuhrow/yafu/tree/8110dfbd8c6f9486d93b1a02de6eb7b180e55a80) | C/GMP/SIMD; `factor/qs/tdiv_resieve.c` and `factor/batch_factor.c` inform residual recovery, root reuse and batching. Limb/SIMD gains need separate Python evidence. | Public-domain notices in inspected files; dependencies separate. |
| [Historic msieve mirror, `c8727d91305bdbe0972d160ef0ce61dd02ce9193`](https://github.com/radii/msieve/tree/c8727d91305bdbe0972d160ef0ce61dd02ce9193) | 2011 C reference for DLP relations/filtering and Block Lanczos; not a current official release. | Public-domain notices in inspected files. |
| [CADO-NFS, `692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b`](https://github.com/cado-nfs/cado-nfs/tree/692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b) | C/C++, GMP/OpenMP/MPI; `filter/merge.cpp`, `polyselect/ropt_main.cpp`, lattice sieving and Block Wiedemann. README targets >85 digits and gives no support <60; these are project support choices. Some developer notes warn they are obsolete. | LGPL-2.1-or-later; embedded/optional components separate. |
| [Yamaquasi 0.3, `3f95f43682ed15d8c1ed206a9a702dd655d7c8ad`](https://github.com/remyoudompheng/yamaquasi/tree/3f95f43682ed15d8c1ed206a9a702dd655d7c8ad) | Rust fixed-width arithmetic and packed GF(2). `src/relations.rs` restricts cycles to the SLP-connected component; that is a yield heuristic, not general completeness. Python bindings retain native execution. | BSD-3-Clause. |
| [FLINT, `a4c9750d0d3d67bb01cf6d18c187591b313451c3`](https://github.com/flintlib/flint/tree/a4c9750d0d3d67bb01cf6d18c187591b313451c3) | C/fmpz/GMP; `src/qsieve/collect_relations.c` intake uses small-cofactor/size/GCD checks without an explicit primality call there. Supports investigating a proved shortcut or opaque matching, not silently changing v2's prime-residual API. | LGPL-3.0-or-later in inspected QS files. |
| [SymPy, `2f22a5f81e2f4124380be3739a092e9ff20128de`](https://github.com/sympy/sympy/tree/2f22a5f81e2f4124380be3739a092e9ff20128de) | Python QS with optional GMP-backed ground types, heuristic scaled logs and SLP `isprime` calls. `qs_factor` documents that returned keys need not all be prime; reconstruct and independently classify reference results. | BSD-3-Clause core plus component notices. |
| [JavaMath, `088d01fa97e7d0412c6cbbb8f463fa7e5f78ab97`](https://github.com/TilmanNeumann/java-math-library/tree/088d01fa97e7d0412c6cbbb8f463fa7e5f78ab97) | BigInteger/JVM/optional Unsafe; SIQS variants explore small-prime skipping, equal-log groups and wide scans. `TDiv_QS_2LP` uses probable-prime testing, not v2's proven-residual guarantee. | GPL-3.0-or-later. |
| [numthy 0.2.1, `1496df97faec010c909e54e14d6bf9cc4399e5fe`](https://github.com/ini/numthy/tree/1496df97faec010c909e54e14d6bf9cc4399e5fe) | Pure Python `numthy.py`: score translation is already in v2; threshold-mask extraction is a narrower candidate. Python-int incidence is useful, but unbounded partial retention and enumerating all possible LP primes do not transfer directly. | MIT root license; header also requests author acknowledgement. |
| [SSS, `8dbaf6d39ab88a40380965d25ec2c363d7f27358`](https://github.com/sbaresearch/smoothsubsumsearch/tree/8dbaf6d39ab88a40380965d25ec2c363d7f27358) | Python/SymPy with vendored historical sSIQS comparator; upstream has unbounded loops/exit paths. v2 already supplies bounded SSS/SSSf, so refresh comparison rather than duplicate the adapter. | No explicit root license found at this pin; clarify permission before copying. Vendored terms separate. |
| [M4RI, `ee69c3935da14a06757d37247a51c5428203ec3a`](https://github.com/malb/m4ri/tree/ee69c3935da14a06757d37247a51c5428203ec3a) | C dense GF(2), echelon/PLE/Four Russians. Small pivot tables are a plausible Python-int transfer, with rank/provenance/cache costs measured. Use the paper rather than the README's suspect typeset complexity expression. | GPL-2.0-or-later. |
| [CUDA-MPQS README v1.0.8, `b1a9c4500da15a2bea56a59ecad7ef61cf054756`](https://github.com/drjanosch42/cuda-mpqs/tree/b1a9c4500da15a2bea56a59ecad7ef61cf054756) | CUDA/C++; bucket layouts, saturation and matrix/provenance costs are leads. Its [paper](https://arxiv.org/html/2610.07126v1) reports H100 RSA-100 at 29.2 s versus 106–122 s on 96-core EPYC YAFU. These unreproduced cross-hardware author results establish no PyPy or Apple-GPU win. | LGPL-3.0 with NVIDIA linking permission; separately versioned algebra/dependencies require their own audit. |

No external implementation was adapted. Future copying must satisfy the relevant
notices and distribution terms; no explicit top-level Factor license was found
in the inspected checkout, so compatibility must not be assumed. Independent
mathematical implementations should cite the algorithm's authors.

## Mathematical checks and primary sources

- **SLP has no universal 85-digit failure boundary.** For `N=77`,
  `10**2-N=23` and `13**2-N=4*23`. Pairing gives `X=53`, `Y=46` modulo N,
  and the two GCDs give 7 and 11. This refutes the universal interpretation of
  the CUDA README statement, not its measured configurations or its paper's
  configuration-dependent genus observations. Measure useful divisor yield.
- **Opaque composite matching can preserve the square identity.** For
  `N=469`, `22**2-N=15` and `23**2-N=4*15`. The composite unit cofactor 15
  yields `X=37`, `Y=30`, and GCDs 7 and 67. This verifies an identity, not a
  performance win or compatibility with the current residual type.
- **The residual bound needs coverage and strict inequality.** If all possible
  prime factors below B were excluded, a composite survivor is at least B².
  Equality is unsafe (`B=7`, `r=49`). A sparse base containing only 2 with bound
  100 leaves `31**2-77=4*13*17`, residual 221 < 100². v2's public FactorBase
  validates supplied entries/roots, not completeness; cover A/multiplier primes
  as well before using this shortcut.
- **General DLP cycles need every component.** The edges `(p,q), (q,r), (r,p)`
  cancel all LP parity without reaching the SLP vertex 1. Compare heuristic
  collectors against a complete cycle oracle.
- **Smooth parts require full prime powers.** For input 72 and prime product
  30, `gcd(72,30)=6`, not 72. [Bernstein, Algorithm 2.1/Theorem 2.2](https://cr.yp.to/factorization/smoothparts-20040510.pdf)
  uses a sufficient modular power; valuation recovery still has a cost.
- **Gram-kernel membership is insufficient over GF(2).** For `M=[1;1]` and
  `d=1`, `M.T*M*d=0` while `M*d` is nonzero. Verify against the original
  operator, including after filtering/lifting and finite sparse-solver retries.

[Alford–Pomerance, equations 3.1–4.1](https://math.dartmouth.edu/~carlp/implementing.pdf)
and [Carrier–Wagstaff, §3.2](https://homes.cerias.purdue.edu/~ssw/qs4.pdf)
support normalized SIQS identities and CRT/Gray reuse; v2 already implements
the core techniques. Preserve A in the square-difference identity and handle
exceptional primes exactly. [Lenstra–Manasse (1994)](https://doi.org/10.1090/S0025-5718-1994-1250773-9)
provides the classic DLP method; its historical speed factors are not PyPy
predictions. [Albrecht–Bard–Pernet](https://arxiv.org/html/1111.6549v1) explains
table-based elimination, pivot/rank conditions and cache tradeoffs.
[Bouillaguet–Zimmermann structured elimination](https://journals.flvc.org/mathcryptology/article/download/126033/127787)
is already referenced by upstream CADO's merge implementation; OpenMP scaling
does not establish Python-thread scaling.

[SSS §4](https://arxiv.org/html/2301.10529v2#S4) measures full factorization for
30–70 digits, but 75–100 digits compare one-hour relation counts on two inputs
per size. Its 5.1–6.8× collection ratios and completion extrapolation are not
measured full-factor speedups. [Bradford–Monagan–Percival §2.3–2.4](https://www.cecm.sfu.ca/~mmonagan/papers/NT4.pdf)
describes the external-q polynomial variant and its setup tradeoffs.

[CADO's pipeline](https://cado-nfs.gitlabpages.inria.fr/),
[Bai–Brent–Thomé root optimization](https://arxiv.org/abs/1212.1958),
[the refined NFS analysis](https://arxiv.org/html/2007.02730v2) and
[alternative NFS sieving](https://eprint.iacr.org/2023/801) motivate later norm,
root-quality, special-q and two-sided cofactoring experiments. Smoothness
quality and NFS complexity estimates have heuristic assumptions; lower-order
terms preclude inferring a finite digit crossover. Equal rational primes can
name different roots/sides of GNFS ideals. QS edges are not that contract.

Blogs and implementation discussions were discovery leads. The original AMS
DLP PDF and Mersenne Forum pages returned 403; author-uploaded text and quoted
YAFU discussion excerpts were inspected without claiming full publisher/forum
access. FLINT discussion search excerpts did not provide a readable full page.
No downloaded implementation was locally benchmarked.

## Experiment contract for the next owner

Freeze the accepted integrated B1/B2/A10 control, runtime/backend/source pins,
versioned training/untouched-confirmation inputs, expected factors/certainty,
seeds, finite budgets/storage caps and selection rules before timing. Profile
separately and stratify by actual bottleneck. Preserve unsuccessful runs,
capacity refusals and unresolved cofactors; distinguish factor-one from complete
factorization. Recalibrate changed collectors on training inputs only.

Use the shared machine-wide performance lock, PyPy implementing Python 3.11,
at least three seconds validated warmup and nine interleaved samples; extend
unstable measurements or retain an inconclusive result. Keep cold startup and
instrumented profiles separate. The existing promotion policy requires a
prespecified complete-time/completion improvement with uncertainty and no
correctness, certainty or resource regression. A coverage/capacity gate must be
declared explicitly; primitive throughput alone cannot promote a default.

Charge setup, certification/splitting, retained stores, filtering, provenance,
recovery and final classification. Validate cancellation/exhaustion, corrupt
checkpoints, cumulative resume/RNG/work and exact reconstruction. Pin feasible
external builds and disclose backend, hardware and threads; add proof work to
match certainty requirements. Native timings are context, not proof of a faster
PyPy implementation. B14 retains prime certificates; F2/F3 retain complementary
ECM work after B2. E1/H1 combined acceptance remains open.
