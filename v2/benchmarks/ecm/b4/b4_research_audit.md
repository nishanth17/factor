# B4 research coverage audit — 9 October 2026

This follow-up checks the breadth of the completed
[B4 research record](b4_research.md). It adds source coverage and clarifies
claims after the frozen comparison; it changes no candidates, controls,
timing results or production code.

The supported claim is a source-informed comparison of five portable
Montgomery x-coordinate candidates on the declared PyPy workloads.
Neither the literature review nor the experiment proves an optimal kernel
or establishes the fastest ECM implementation across representations and
architectures. EFD's smallest listed costs are not an arithmetic lower-bound
proof, and the fastest measured arm is the fastest of the frozen candidates.

## Original source-level coverage

The original inspection read arithmetic functions and license files, with
immutable revisions, source URLs and SHA256 manifests retained locally.

| Source | Functions/mechanisms reviewed | Consequence for the bounded candidates |
| --- | --- | --- |
| [GMP-ECM, 8ea5e214](https://github.com/sethtroisi/gmp-ecm/blob/8ea5e214fdf2f0ddf9415141b8dc039ed6f5874e/ecm.c#L148) | add3, duplicate, ecm_mul; modular squaring, canonical temporaries and alias-safe output | Explicit squares and selected reductions need Python measurements; native residue primitives are not Python operators. |
| [CADO-NFS, 692ecb7e](https://github.com/cado-nfs/cado-nfs/blob/692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b/sieve/ecm/ec_arith_Montgomery.hpp#L80) | Checked affine conversion, montgomery_dbl_inl, montgomery_dadd_inl, exceptional-state comments and backend-dependent inlining | Unit checks and shared expressions transfer; C++ inlining results and zero sentinels do not establish PyPy speed or valid projective equality. |
| [AVX-ECM, 118e015c](https://github.com/bbuhrow/avx-ecm/blob/118e015caba54747c6ac9f0e293bfd9f8d58b3a0/ecm.c#L407) | vec_add, vec_duplicate; shared sums/differences, explicit squares, affine-difference saving and vector limb arithmetic | Fusion and fixed-difference normalization transfer algebraically; SIMD lanes and native carry/reduction scheduling require a different execution path. |
| [Yamaquasi, 3f95f436](https://github.com/remyoudompheng/yamaquasi/blob/3f95f43682ed15d8c1ed206a9a702dd655d7c8ad/src/ecm128.rs#L349) | Extended/projective Edwards dbladd and fixed-width arithmetic | This is a different representation and scalar execution package; its multiplication counts do not describe a replacement for the current x-only ladder step. |

The original record covers primary formula sources, the EFD, RFC 7748,
mixed-representation ECM, official gmpy2 documentation, and technical
blog/discussion leads checked against primary formulas. It did not name all
the additional papers and implementations below.

## Additional implementation inspection

YAFU is now pinned at
[8110dfbd8c6f9486d93b1a02de6eb7b180e55a80](https://github.com/bbuhrow/yafu/tree/8110dfbd8c6f9486d93b1a02de6eb7b180e55a80)
(commit dated 5 October 2026). Focused source review covers:

- [microecm.c, uecm_uadd and uecm_udup](https://github.com/bbuhrow/yafu/blob/8110dfbd8c6f9486d93b1a02de6eb7b180e55a80/factor/gmp-ecm/microecm.c#L616):
  the standard differential-addition and doubling polynomials, explicit
  modular squares, precomputed sums/differences for doubling, scalar native
  mulredc primitives and an eight-lane vector path. The header offers a choice
  of BSD-2-Clause/FreeBSD terms or MPL-2.0.
- [tinyecm.c, add and duplicate](https://github.com/bbuhrow/yafu/blob/8110dfbd8c6f9486d93b1a02de6eb7b180e55a80/factor/gmp-ecm/tinyecm.c#L195):
  two-word modular arithmetic, explicit squares, alias handling, and the
  explicit warning that the input difference cannot simply be assumed to
  have Z=1. The file header carries BSD-2-Clause/FreeBSD terms.

These functions confirm the same core x-only arithmetic and the
normalization precondition already represented in B4. Native word-size
bounds, assembly, SIMD and mulredc implementations are not equivalent to
rearranging PyPy integer operators. Compact PRAC programs belong to C6/B3.
No code was adapted. Downloaded copies in other lasieve directories are
provenance captures, not a claim of equally detailed review of every variant.

GMP-ECM's inspected mirror commit is dated 9 March 2026. Official
[Sage package documentation](https://doc.sagemath.org/html/en/reference/spkg/ecm.html)
identifies upstream Inria and release 7.0.6. The audit retrieved Sage's
[release checksum metadata](https://raw.githubusercontent.com/sagemath/sage/develop/build/pkgs/ecm/checksums.ini),
but available release download attempts failed retrieval or checksum
validation. No unverified archive was accepted, and no equivalence between
the inspected mirror functions and the release archive is asserted.
The immutable inspected commit remains the reproducible source reference;
the review does not certify the newest upstream release or every later patch.

## Formula costs and newer literature

The [EFD XZ table](https://hyperelliptic.org/EFD/g1p/auto-montgom-xz.html)
lists the shared projective ladder step at 6M+4S+1*a24 and its
unit-normalized-difference form at 5M+4S+1*a24. Sharing removes two
add/subtract operations compared with separate point operations.
These are the formulas exercised by B4's fusion and normalization arms.
For random Suyama curves, multiplying by a24 is a full-width operation;
it cannot be silently priced as multiplication by a small cryptographic
constant. EFD costs also omit Python allocation, JIT and remainder costs.
They justify candidate selection, not an optimality claim.

[Farashahi, Fadavi and Sabbaghian, Faster Complete Addition Laws for
Montgomery Curves, TCHES 2024(4), 737–762](https://ches.iacr.org/2024/papers-issue-4/4_111.pdf)
was missing from the original bibliography. Sections 4–5 and Theorems 2–3
use four extended coordinates (U:V:S:W) and multiplication-free maps to
extended Edwards coordinates. Their improved complete-addition costs are
for that representation, with stated finite-field conditions; some scaled
forms require square roots and quadratic-character assumptions.
They are not a lower-cost x-only differential-addition formula. Applying
them to ECM would require a separately proved representation/setup/recovery
package over composite moduli; field completeness alone is insufficient.

[Bernstein, Cottaar and Lange, Searching for differential addition chains,
Research in Number Theory 11, 45 (2025)](https://link.springer.com/article/10.1007/s40993-024-00604-8)
was also missing from the original record. Sections 1.1 and 3–5 distinguish
chain length from weighted M/S/constant/addition costs and introduce improved
pruned and meet-in-the-middle searches. Section 1.1 explicitly leaves
lower-level performance analysis to future work. This belongs to C6's chain
comparison; it does not prove a faster replacement for an xDBLADD step.

The primary [ECM using Edwards curves paper](https://eprint.iacr.org/2008/016)
and [EECM-MPFQ software overview](https://eecm.cr.yp.to/mpfq.html) identify
extended Edwards coordinates, signed windows, batched primes, small
parameters/base points and torsion as a combined optimization package.
[ECM at Work](https://eprint.iacr.org/2012/089) further studies chain
combination and GPU memory/performance. These sources were reviewed here as
literature/architecture context; EECM-MPFQ source was not inspected.
Their historical native comparisons cannot be transplanted into a PyPy
speed claim. Representation and curve-selection changes remain separate F2
work, and scalar program changes remain C6/B3.

## Remaining limits and decision

B4 does not exhaust reduction placements or candidate combinations, implement
setup-inversion reuse, search all differential chains, or evaluate every
native/SIMD/GPU ECM package. Prime-specific cryptographic reducers and
small-constant assumptions do not hold for arbitrary factoring moduli and
random Suyama parameters. No proof of global arithmetic optimality or
globally fastest PyPy implementation is claimed.

The follow-up sources supply useful coverage and explain why published
native or representation-changing improvements are not direct evidence for
a faster drop-in kernel on the frozen controls. They establish no new PyPy
winner. The bounded comparison and its positive native/GMP measurements
stand; any additional candidate requires a separate prespecified study.

New downloads, file hashes, repository tree metadata and rejected-release
provenance stay local under
results/b4/research-audit-20261009/. Required controls/corpora and accepted
timing sources remain unchanged. No benchmarks, heavy checks, production
edits or merges are part of this coverage audit.
