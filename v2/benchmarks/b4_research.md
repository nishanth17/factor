# B4 / P4.2: bounded Montgomery kernel study

Research inspected 9 October 2026. The control is committed mainline
`bcf5f3d1e57304694b48ba6e7ef8b4ea2ffd0db0`; experimental arithmetic,
independent controls and the comparison protocol are frozen at `08ccac2`.
Production arithmetic, defaults, schedules, checkpoint formats and `v1/`
are unchanged. This study compares five candidates and stops after fresh
confirmation. It does not reopen the stopped P4.3 backend study, C5 reducers,
C6 chain search or F2 curve-family work.

## Research and provenance

The arithmetic source is Montgomery's 1987 paper, pp. 260–261,
[Speeding the Pollard and elliptic curve methods of factorization](https://wstein.org/edu/124/misc/montgomery.pdf).
[Costello–Smith, arXiv:1703.01863v1, Algorithms 1–2](https://arxiv.org/pdf/1703.01863v1)
explains differential addition and doubling on the x-line, including the
excluded differences at infinity and the rational two-torsion point.
[Bernstein–Lange, 2017/293](https://eprint.iacr.org/2017/293) gives a
completeness proof for its cryptographic ladder; those field/encoding
conditions do not license accepting degenerate ECM points. The [EFD XZ formulas](https://hyperelliptic.org/EFD/g1p/auto-montgom-xz.html)
provide fused and normalized forms with explicit assumptions and verification
scripts. Their multiplication counts are algebraic costs, with parameter
multiplication accounted for separately; a random Suyama a24 is a full-width
integer. They are not measured PyPy costs.

The following immutable source revisions were retrieved, hashed locally,
and read alongside their licenses. No upstream code was copied or translated;
the Python candidates are independent rewrites of our readable control.

| Implementation / pinned source | Inspected mechanism | License inspection / transfer decision |
| --- | --- | --- |
| [GMP-ECM `8ea5e214fdf2f0ddf9415141b8dc039ed6f5874e`](https://github.com/sethtroisi/gmp-ecm/blob/8ea5e214fdf2f0ddf9415141b8dc039ed6f5874e/ecm.c#L148) | `add3`, `duplicate`, `ecm_mul`: explicit modular squares, canonical modular temporaries, adjacent-multiple ladder; b=(A+2)/4 | `ecm.c` and `COPYING.LIB`: LGPL-3.0-or-later. Mathematical formulas inform independent Python experiments; native residue storage and reducers are outside B4. Its zero-coordinate conventions are not imported. |
| [CADO-NFS `692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b`](https://github.com/cado-nfs/cado-nfs/blob/692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b/sieve/ecm/ec_arith_Montgomery.hpp) | `montgomery_dbl_inl`, `montgomery_dadd_inl`, `montgomery_point_to_affine`; checked inversion; comments distinguish fixed-size inlining from mpz costs | Root `COPYING`: LGPL-2.1; other bundled components have their own licenses. Portable common-subexpression sharing and unit checks merit measurement. C++ inlining percentages do not predict PyPy JIT behavior. Its optional `ELLM_SAFE_ADD` and (0:0) sentinel are not our projective validity contract. |
| [AVX-ECM `118e015caba54747c6ac9f0e293bfd9f8d58b3a0`](https://github.com/bbuhrow/avx-ecm/blob/118e015caba54747c6ac9f0e293bfd9f8d58b3a0/ecm.c#L407) | `vec_add` consumes shared sums/differences; `vec_duplicate` uses explicit vector squares; fixed-difference Z=1 comment identifies a saved multiply. `vecarith52.c` uses native vector limb arithmetic | Root `LICENSE`: BSD-2-Clause, but `ecm.c` also carries upstream ECM LGPL notices and `vecarith52.c` has Mayo Apache-2.0 notices. Do not infer a uniform license from the root alone. AVX-512, SIMD lanes, limb packing and carry scheduling do not transfer to Python operators on Apple ARM64. |
| [Yamaquasi `3f95f43682ed15d8c1ed206a9a702dd655d7c8ad`](https://github.com/remyoudompheng/yamaquasi/blob/3f95f43682ed15d8c1ed206a9a702dd655d7c8ad/src/ecm128.rs) | Fused `dbladd` and extended/projective Edwards operations, plus native fixed-width arithmetic in the small-integer path | Root `LICENSE`: BSD-3-Clause. Useful evidence that operation fusion is representation-specific; these Edwards formulas are not a Montgomery xDBLADD replacement. Curve changes and native reducers are deferred. |

[Bouvier–Imbert, Faster cofactorization with ECM using mixed representations](https://eprint.iacr.org/2018/669.pdf)
(PKC 2020) combines prime batches, double-base/Lucas chains and mixed curve
representations for NFS cofactorization. Its authors explicitly tie square
cost assumptions to their native 64/96/128-bit arithmetic. The paper supports
separate ownership of chain and representation experiments, not substituting
those algorithms into this kernel comparison. The
[authors' software page](https://www.lirmm.fr/eco/softwares.php) links their
chain-generation implementation; C6 owns that avenue.

[Monocypher author Loup Vaillant's technical tutorial](https://loup-vaillant.fr/tutorials/fast-scalarmult)
emphasizes setup/table costs and limits to precomputation benefits. It is
context for amortization rather than a proof of our x-only formulas.
The [EFD/RFC convention discussion](https://crypto.stackexchange.com/questions/67942/difference-on-montgomery-curve-equation-between-efd-and-rfc7748)
was checked against the primary [RFC 7748, section 5](https://www.rfc-editor.org/rfc/rfc7748.html#section-5):
the RFC uses (A−2)/4 and **AA**, whereas our implementation uses (A+2)/4 and
**BB**. Discussion threads are leads; primary formulas and the proof below
settle correctness. [Official gmpy2 integer documentation](https://gmpy2.readthedocs.io/en/stable/mpz.html)
was inspected together with the installed PyPy extension. The square arm
uses integer `**2` on both supported types; no per-operation backend callback,
mutable xmpz, new dependency or floating arithmetic is introduced.

A later [research coverage audit](b4_research_audit.md) adds YAFU, newer
literature and explicit coverage/optimality limits without changing this
frozen comparison.

## Formula and normalization proof

Write S=X+Z, D=X−Z, AA=S², BB=D², E=AA−BB=4XZ, and
c=(A+2)/4 modulo odd n. The unchanged doubling is

```
X(2P) = AA*BB
Z(2P) = E*(BB+c*E)                       (mod n)
```

Because AA=BB+E, `BB+c*E = AA+(c−1)*E`; the RFC coefficient is c−1.
Replacing BB by AA while retaining c is incorrect. All five candidates keep
c and BB; tests also exhibit the incorrect-convention counterexample.

For differential addition, let U=(XP−ZP)(XQ+ZQ),
V=(XP+ZP)(XQ−ZQ), with R=P−Q up to sign. The output is
`(ZR*(U+V)², XR*(U−V)²) mod n`. A fused step shares S and D between this
addition and doubling, saving two sums/differences, without changing the
polynomials. Whole-ladder fusion inlines that step. The pre/post bit swaps
select which adjacent multiple is doubled: a zero bit maps (k,k+1) to
(2k,2k+1); a one bit maps it to (2k+1,2k+2). The fixed difference remains P.
Explicit squares replace multiplication by the same integer with `**2`.
Selected reductions replace AA, BB, U and V by their residues before their
next products. Polynomial congruence proves exactly the same canonical
output pair, including exceptional inputs, for all unnormalized arms.

Fixed-difference normalization is legal only when ZR is a unit. Compute
`t=ZR⁻¹ mod n`, `R'=(t*XR,1)`. Differential addition with R' multiplies
both old output coordinates by the same unit t, so the X-coordinate no
longer needs multiplication by ZR. If working points have unit scales a,b,
the next addition has scale a²b²t and doubling has scale a⁴. Induction proves
that every ladder output differs from the readable control by a unit.
Thus its projective action, denominator GCD and coordinate ideal agree,
also over composite rings and prime powers. This argument does not justify
normalization by a nonunit or accepting (0:0) as a point.

`normalize_difference` retains the exact GCD in
`NormalizationNonunitError.divisor`, with a separately validated proper
`factor` or a saturated retry. The private stage adapter terminates the
already-reserved arithmetic action with that proper factor or curve retry;
it commits no partial point. Scalars 0/1 and infinity preserve the original
shortcuts. The public production engine is untouched.

Suyama setup inverts d=16*u³*v. When d is a unit, one can derive
`v⁻¹=16*u³*d⁻¹`, then `(v³)⁻¹=(v⁻¹)³`. Reusing d⁻¹ directly as the Z
inverse is wrong. This study deliberately pays a separate checked inverse
at each scalar entry, including chunk and giant initialization; later chunk
points have new Z values. Reusing setup inversion or caching normalized
chunks is not an uncharged saving or a sixth tuned candidate.

## Width, finite-resource and checkpoint contracts

Assume all incoming coordinates and c are in [0,n), as produced by setup
and every kernel exit, with b=bit_length(n). Before reduction,
`|S|<2n`, `|D|<n`, `AA<4n²`, `BB<n²`, `E<4n²`.
Late doubling has `X<4n⁴`, `|Z|<20n⁵`; late addition has
`|U|,|V|<2n²` and coordinate numerators below 16n⁵. The conservative
peak is **5b+5 bits**, with signed differences handled by exact integers.
Fusion and explicit squares preserve these bounds.

The selected-reduction loop has AA,BB,U,V in [0,n), `|E|<n`,
`|U+V|<2n`, and numerators below 4n³: **3b+2 bits** in its steady loop.
It adds four remainders per bit. Its one initial doubling still uses the
late 5b+5 bound. Normalization reduces the addition X numerator below
16n⁴ but retains the late bound for other operations. No encoded residues,
Barrett/Montgomery reduction context or native backend is introduced.

All candidates use constant-count arithmetic temporaries and the original
finite scalar bit string. The largest diagnostic modulus is 1024 bits;
full portfolios cap inputs at 331 bits, workspace at 16 MiB, work at
5,000,000 and wall/CPU at 30 seconds per input/seed. Chunk boundaries,
reservations, cancellation polling, saturation replay and result validation
come from the immutable engine. New normalization temporaries are bounded
by these widths and fit its conservative workspace reservation. No table
or scalar-dependent search is added. The harness additionally caps each
worker at 600 seconds and preserves uncensored/stable requirements.

Every experimental exit is canonical X:Z with the existing a24 convention.
Normalization retains no private representation in checkpoints. Tests cover
JSON roundtrips, cumulative work, cancellation, same-arm resume, and resuming
candidate checkpoints with the unchanged ladder engine. The standard
checkpoint schema and backend identity remain intact.

Independent affine group laws use full (x,y) coordinates; CRT controls check
both component fields, and Hensel lifts provide prime-power controls.
Primitive pairs require `gcd(X,Z,n)=1`; affine comparison also checks that Z
is a unit when an affine output is expected. Low-order and shared-factor
pairs have explicit degeneracy assertions. The benchmark labels matching
nonprimitive ideals as exceptional outcomes, never point-equality successes.
