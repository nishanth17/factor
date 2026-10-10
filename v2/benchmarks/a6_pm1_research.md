# A6 / P5.3: p−1 and exact bound continuation

Research and contract frozen 9 October 2026. Control: committed mainline
`bcf5f3d1e57304694b48ba6e7ef8b4ea2ffd0db0`. This tranche concerns ordinary
modular powering on PyPy Python 3.11. Williams p+1, Lucas optimization,
ECM pairing, ECM continuation and portfolio allocation have separate owners.

## Mathematical basis and independent derivation

[Pollard (1974), §4](https://www.cambridge.org/core/journals/mathematical-proceedings-of-the-cambridge-philosophical-society/article/abs/theorems-on-factorization-and-primality-testing/6762E84DBD34AEF13E6B1D1A8334A989)
is the original p−1 paper. The publisher exposes its bibliographic record;
access to its complete text was unavailable here. The operative stage-one
and stage-two descriptions were checked against [Kruppa's primary thesis,
§§2.2.2, 4.4, 4.7](https://docnum.univ-lorraine.fr/public/SCD_T_2010_0054_KRUPPA.pdf).
They use an LCM exponent, element orders and a single additional large prime.
This is a sufficient condition for divisibility by a prime factor, not a
promise that a proper divisor will survive a mixed-factor GCD.

The following ratio argument is an independent derivation used by the tests.
For inclusive integer B ≥ 1, define M(B)=lcm(1,…,B). Unique factorization gives
M(B)=∏ p^e(B,p), where e(B,p) is the largest integer e with p^e ≤ B (zero
when p>B). For 1 ≤ L ≤ U, M(L) divides M(U), so

```
R(L,U) = M(U)/M(L) = ∏[p ≤ U] p**(e(U,p)-e(L,p)).
```

Every prime ≤ U is considered, including primes already ≤ L. For example,
R(7,10)=2·3, with no new primes; R(15,17)=2·17. Computation uses integer
multiplication and exact division, never logarithms or floating roots.
For any residue a, including a nonunit, `(a**M(L))**R(L,U) = a**M(U)` modulo n.
The order/smoothness interpretation additionally requires gcd(a,n)=1.
A smaller order than p−1 can give incidental success below the sufficient
bound. Tests therefore construct and check exact base orders, rather than
asserting failure solely from the factorization of p−1.

For stage two, with A=a**M(B1), every prime B1<q≤B2 contributes A**q−1.
For successive q,r, A**r=A**q·A**(r−q); this is the exact gap-cache invariant.
The first relation starts at exponent zero and must also be included.
[Montgomery (1987)](https://wstein.org/edu/124/misc/montgomery.pdf) develops
faster continuations. [Montgomery and Kruppa's ANTS VIII presentation
(2008)](https://antsmath.org/ANTSVIII/files/kruppa.pdf) describes a polynomial
stage two. Those larger algorithms are reference material and deferred here.

## Pinned implementation and license inspection

| Source pin and files inspected | Useful distinction and decision |
| --- | --- |
| [GMP-ECM `8ea5e214`](https://github.com/sethtroisi/gmp-ecm/blob/8ea5e214fdf2f0ddf9415141b8dc039ed6f5874e/pm1.c): `pm1.c`, `pm1fs2.c`, `README`, `COPYING.LIB` | Stage one multiplies prime powers into bounded exponents/cascades; `r > B1done` includes increased powers of old primes. Its GMP sliding-window, limb thresholds and NTT stage two are native optimizations, not PyPy speed evidence. Library source is LGPL-3.0-or-later. No code copied. |
| [CADO-NFS `692ecb7e`](https://github.com/cado-nfs/cado-nfs/blob/692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b/sieve/ecm/pm1.cpp): `pm1.cpp`, `pm1.hpp`, `stage2.c`, `COPYING` | Bound-owned plans contain the stage-one exponent and separate powers of two for backtracking. The implementation converts x to x+1/x for its p+1-style second stage; its REDC layers and fixed-word cofactor workloads differ from this Python path. The repository's COPYING is LGPL-2.1; inspect individual files before any future adaptation. No code copied. |
| [FLINT `17950040` / v3.3.1](https://github.com/flintlib/flint/blob/17950040404e6ed797a4becd8a866fb3f62b5c5e/src/fmpz_factor/factor_pp1.c): `factor_pp1.c` in `fmpz_factor` and `ulong_extras`, `COPYING.LESSER` | These are p+1 implementations with Lucas recurrences, recovery and specialized multiprecision/word arithmetic. They support keeping shared integer schedules separate from group actions; their recurrence cannot replace ordinary p−1 powering. File headers specify LGPL-3.0-or-later. No code copied. |
| [PARI official 2.19.0 release](https://pari.math.u-bordeaux.fr/download.html): `src/basemath/ifactor1.c`, `COPYING` | Archive SHA-256 `f317b9722eb5d9094a60303774f066f3a83e3ec1f170be8546c44d7583f30b6d` verified before inspection. Its factoring driver provides ECM and other methods; no standalone ordinary p−1 continuation API was found in the inspected driver. This is not evidence that PARI has no use of p−1 elsewhere, especially primality proofs. COPYING is GPL-2.0; the source header permits GPL-2.0-or-later. No code copied. |
| [SymPy `fe935ceb` / 1.14.0](https://github.com/sympy/sympy/blob/fe935ceb303891d1f8bea4c03b19fd9ec9464b02/sympy/ntheory/factor_.py): `pollard_pm1`, `LICENSE` | Clear LCM/power-smoothness and base-dependent saturation examples, with a first-stage reference rather than this two-stage bounded campaign. BSD-3-Clause. No code copied. |

Source text and license captures are local, ignored research evidence. Repository
license ambiguity is avoided by independently implementing the integer
identities and using the project's existing sieve, budget and result tools.

The [GMP-ECM README's saved-residue example](https://github.com/sethtroisi/gmp-ecm/blob/8ea5e214fdf2f0ddf9415141b8dc039ed6f5874e/README)
confirms that increased B1 is an explicit resume operation with bound metadata.
Its native timings and automatic B2 selection are not adopted. The
[GIMPS mathematical overview](https://www.mersenne.org/various/math.php)
has specialized known-order/Mersenne structure, which does not transfer as a
general-composite allocation rule.

Two secondary leads were also inspected: [Programming Praxis's continuation
exercise and discussion](https://programmingpraxis.com/2010/03/19/extending-pollards-p-1-factorization-algorithm/)
and [this probability discussion](https://crypto.stackexchange.com/questions/75033/estimating-the-probability-of-sucess-of-pollards-p-1).
Their useful leads are checked against the thesis and implementations above;
neither supplies mathematical guarantees, licenses for adapted code, or
performance acceptance evidence for this tranche.

## Reusable continuation contract

`RATIO_VERSION="inclusive-lcm-ratio-v1"` identifies the integer schedule.
`prime_power_ratio(p,L,U)` assumes a caller-verified prime p; it returns the
exact ratio, including 1. `prime_power_ratios(L,U,**sieve_options)` verifies
primes through the bounded sieve and yields only nonidentity pairs. Old bound
1 denotes M(1)=1. This streaming helper itself follows the existing sieve's
work-accounting convention; a consumer must reserve its arithmetic separately.

Reuse requires the *complete* M(L) action for the same modulus and base/group
assignment. A last-prime cursor from an unfinished stage is not such a
certificate. The p−1 campaign retains one chunk start until its GCD succeeds;
saturation replays prime units within the finite recovery allowance. A fully
saturated residue is 1 modulo n and remains 1 under every positive ratio, so
this base stops. An incompatible checkpoint is rejected, never converted into
a purported completed-bound residue. A caller can explicitly start a fresh
campaign/base, but repeated bases do not create independent random group
orders as ECM curves do.

With increased B1, only the verified stage-one residue and modulus/base
identity carry forward. Prime generation/ratio compilation is repeated over
all primes ≤ new B1. Stage-two primes in the *new* interval (new B1,new B2]
must be processed using the new residue, including overlap with the previous
stage-two interval. Previous gap powers, products, terms and running powers
are invalidated. With unchanged B1 and increased B2, checked stage-two batches
can be retained as completed coverage; the running power, previous prime and
gap cache append only the new interval. No unchecked batch crosses a rung.

The ratio also specifies the integer scalar needed by later ECM extensions:
[M(U)]P=[R(L,U)]([M(L)]P). This does not certify an ECM checkpoint, migrate a
campaign, implement saturation recovery for points, or authorize new bounds
or curve counts. Those actions remain with the ECM continuation owner.

## Coverage audit and next bakeoffs — 9 October 2026

The accepted tranche is a bounded comparison, not an exhaustive claim about
state-of-the-art p−1 implementations or literature. Its correctness evidence
stands independently of its candidate selection. The following audit adds
important implementation families that were absent from that selection.
FLINT's inspected routines are p+1, and the inspected PARI driver is not a
standalone p−1 benchmark: neither is counted as an independent p−1 speed
baseline. All candidates below remain unmeasured in this project.

### Additional pinned source inspection

| Source and inspected files | Findings, limits and license |
| --- | --- |
| [Yamaquasi `3f95f436`](https://github.com/remyoudompheng/yamaquasi/blob/3f95f43682ed15d8c1ed206a9a702dd655d7c8ad/src/pollard_pm1.rs), `src/pollard_pm1.rs`, `LICENSE` | Both small-integer reusable factor bases and generic p−1 are present. Generic stage one uses bounded machine-word / 1024-bit exponent blocks. Classical stage two grows even gap powers by recurrence; above its configured threshold it builds a root polynomial and evaluates it by chirp-z/NTT convolution. The native crossover is not a PyPy crossover. Prime-power comparisons and approximate/rounded stage-two parameter selection differ from our inclusive exact-bound contract; normalize coverage before comparison. BSD-3-Clause. |
| [YAFU `8110dfbd`](https://github.com/bbuhrow/yafu/blob/8110dfbd8c6f9486d93b1a02de6eb7b180e55a80/factor/gmp-ecm/micropm1.c), `factor/gmp-ecm/micropm1.c`, `factor/avx-ecm/avxppm1.c` | Micro p−1 has precomputed stage-one window plans and D=60 square-exponent BSGS pairing, with fixed small B1 choices and B2=25·B1. That is an additional ordinary p−1 algorithm family, not merely a wrapper around GMP-ECM. The micro file permits FreeBSD two-clause or MPL-2.0 licensing; the AVX file contains both a two-clause notice and LGPL-3.0-or-later material. Review each adapted region's provenance. Fixed-word REDC and AVX execution do not transfer directly to Python integers. |
| [Prime95 mirror `027cb137`](https://github.com/primesearch/Prime95/blob/027cb13799d46bbcc5dc5fd208cb71ebe4efb348/ecm.cpp), `ecm.cpp`, `pair.cpp`, `license.txt` | The p−1 sections select pairing or polynomial stage two under memory/cost estimates. `calc_exp2` constructs bit-sized exponent batches with balanced multiplication and includes every prime-power threshold above the completed B1. Stage one limits uninterruptible work. The FFT arithmetic and its cost model target very different sizes/hardware. The source has Mersenne Research copyright/all-rights-reserved notices and a custom GIMPS EULA; public availability does not supply a general adaptation license. Use independently derived mathematics, not copied planner code. |
| [GPUOwl v6 `98ff9c78`](https://github.com/preda/gpuowl/blob/98ff9c78068543674866c333c7e46c5c92212686/Pm1Plan.cpp), `Pm1Plan.cpp`, `LICENSE` | This older execution-plan source contains D=30030 pairing and relocation of small multiples of eligible primes. GPL-3.0. Current repository HEAD [`4d0e7590`](https://github.com/preda/gpuowl/blob/4d0e75902d12e57dcea36f23b56bcfe364ad9df2/pm1/pm1.cpp) was separately checked: its `pm1/pm1.cpp` is a Mersenne probability/bound calculator, not a p−1 execution kernel. Do not count the two as independent execution baselines or transfer Mersenne allocation probabilities to generic composites. |
| [PrMers `af6f9965`](https://github.com/cherubrock-seb/PrMers/blob/af6f99659082105c1eb077ccf56a7b3674693fc1/src/modes/RunPM1.cpp), `src/modes/RunPM1.cpp`, `LICENSE`, `LICENSES/README.md` | Current Mersenne-specific code exposes classical BSGS, scalar-trace stage two, a Pair95 option and Prime95 handoff/checkpoint handling. This is a useful recent source lead, not evidence of superiority or an independently validated algorithm here. Authored source is MIT; its optional GPUOwl-derived Aevum engine is GPL-3.0 and files retain their own notices. Whole-program licensing cannot be inferred from the root MIT notice. |

Pins were resolved on the audit date; source/license captures and SHA-256
manifests remain in ignored `results/a6/research/coverage-audit/`. No upstream
code was copied. Mlucas 20.1.1 and CUDAPm1 are additional special-form leads,
not inspected execution baselines in this audit. Proprietary implementations,
all historical forks and every recent paper have not been exhaustively audited.

### Literature coverage and a recent claim

[Montgomery–Silverman (1990), *An FFT Extension to the p−1 Factoring
Algorithm*](https://cr.yp.to/bib/1990/montgomery.pdf) is the primary predecessor
of the later polynomial continuation. Its PDF was captured; the available web
text exposed the abstract and classical gap recurrence, but a complete local
text extraction was unavailable in this audit. [Montgomery–Kruppa (2008),
*Improved Stage 2 to P±1 Factoring Algorithms*](https://inria.hal.science/inria-00188192/file/pm1fft-final.pdf)
is the primary space-efficient reciprocal-polynomial continuation paper.
Its full HAL download was access-blocked; the original tranche inspected the
authors' ANTS presentation and Kruppa's thesis. These access limits must not be
reported as a full reading of those two papers.

The author-hosted [Brent–Kruppa–Zimmermann chapter, §8.3](https://members.loria.fr/PZimmermann/papers/Chap8.pdf)
was inspected in this audit. It supplies the p−1-specific scaling/product
construction and geometric-progression evaluation by Bluestein/chirp-z,
distinct from ECM's generic product/remainder tree. It explicitly notes
possible missing primes near interval edges for some set choices. A bounded
implementation must certify those edges and exceptional wheel primes rather
than infer exact coverage from a convolution size. Polynomial multiplication,
coefficient representation and evaluation costs remain part of stage two;
there is no justified standalone O(log B2) total-stage claim here.

[Xia–Wang–Gu (13 August 2026), *Dynamic Scaling Pollard's P-1
Algorithm*](https://www.mdpi.com/2410-387X/10/4/57) is a recent primary lead.
The publisher's indexed abstract and algorithm discussion were available;
direct HTML/XML/full-PDF retrieval failed. It studies dynamic exponent/prime
scaling, reuse of prime products and balanced multiplication. Its experiment
constructs smooth p−1 inputs, and its named comparisons do not establish a
win over our retained bounded chunk/gap implementation or GMP-ECM. Its
resistance assumptions for the other factor also require review for mixed
saturation. Treat it as a separate schedule hypothesis requiring complete
paper/pseudocode inspection and an exact exponent/coverage oracle before any
matched experiment; do not transfer its speedup claims or silently replace
M(B)=lcm(1,...,B).

### Ranked hypotheses for a separate bounded follow-up

1. **Ordinary p−1 wheel/± BSGS stage two.** This has the largest algorithmic
   opportunity because stage two accounts for substantial inclusive CPU in
   the A6 diagnostic profile. Compare a small bounded wheel against retained
   gap execution, including all table setup, singleton/tail work, GCD and
   finite replay costs. YAFU's square-exponent method and CADO/PrMers' trace
   representation are distinct possible arms; start with one. For unit A,
   writing T=A^(kD), H=A^r gives
   `T + T^-1 - H - H^-1 = (T-H)*(T*H-1)/(T*H)`.
   This independently explains the ± coverage modulo each prime factor.
   Check units before inversion, certify every eligible prime, retain
   exceptional wheel primes, and account explicitly for any incidental
   coverage from singletons or square exponents. Mixed factors can saturate
   a product, so direct-prime replay remains necessary. This is a separate
   p−1 executor; no Williams p+1, Lucas optimizer or ECM pairing changes.

2. **Exponent chunks capped by bit length.** Compare a few frozen bit caps
   with the accepted fixed-prime-count chunk 64 and retained chunk 16.
   Prime counts hide changing exponent lengths across bounds. Include
   sequential versus balanced product construction only where profiling
   shows setup matters. Preserve chunk-start recovery, latency limits and
   cumulative charges. Built-in `pow` remains the control; a Python window
   interpreter is not assumed to beat it because native code uses windows.

3. **Bounded even-gap power recurrence.** Compare the current per-distinct-gap
   `pow(A,gap,n)` cache with growing `A^2, A^4, ...` by multiplication and
   with a sparse compiled gap plan. Keep the first prime exponent separate,
   cap retained powers, and charge unused setup. The existing cache already
   removes repeated gap exponentiation, so the remaining opportunity may be
   small. Measure complete stages rather than just table construction.

4. **Compiled integer schedules reused across different inputs.** Measure
   packed prime-power/chunk/gap/ratio plans with construction, decoding,
   memory and resume rebuilding charged. The A6 generic prime-cache study
   does not settle this richer representation. The existing program store's
   `program_segment`/`power_values` interface is a coordination reference;
   avoid changing the ECM-owned store during a p−1 experiment. Amortize over
   genuinely different inputs at matched bounds. Repeated bases for one n
   are correlated smooth-order trials and cannot justify ECM-like success
   claims.

5. **Polynomial/chirp-z stage-two crossover.** Use GMP-ECM and Yamaquasi as
   architectural references after the classical stage-two alternatives.
   Establish a feasible exact polynomial kernel and finite coefficient/
   workspace/recovery contract before building a new transform backend.
   Compare complete stages over predeclared larger B2 ranges, including the
   project's production B2 when feasible. Native NTT/FFT thresholds and
   GPU/Mersenne speedups supply no Python-integer promotion evidence.

The cheapest first screen is bit-capped chunks plus even-gap recurrence;
the next substantive new algorithm is p−1 wheel/± stage two. A polynomial
implementation and dynamic-scaling schedule remain conditional. Resume
verification is another possible later study: the current replay cost is
measured, but accepting an unverified cached residue is not an optimization.

Freeze a fresh protocol and untouched confirmation inputs before any new
timing. Keep Python integers, exact bound/coverage identities and existing
allowances as controls. Include nonsplitting, independently certified general
inputs and marginal portfolio completion per CPU-second alongside constructed
boundary/recovery cases. Charge fresh construction and distinguish legitimate
batch amortization from warm-only reuse. Require PyPy Python 3.11, at least
three seconds of validated warmup, nine samples with stability extensions,
paired uncertainty estimates, and an exclusive window coordinated with B4/C6.
The revised roadmap permits gains below 10% when sustained beyond noise and
confirmed independently; practical regressions, storage, reconstruction and
completion evidence still govern promotion. This audit performs no new timing
and changes no API, defaults, allocation or roadmap completion status.

### Follow-up wheel contract (implementation pending acceptance)

The frozen follow-up independently implements the trace identity with D=30
or 210. For the nearest center c=kD and eligible primes c-r and c+r,

```
A**c + A**(-c) - A**r - A**(-r)
    = A**(-c) * (A**(c-r)-1) * (A**(c+r)-1)  (mod n).
```

Setup verifies gcd(A,n)=1 before inversion. Thus A**(-c) is a unit modulo n,
so the trace has the same GCD as the two ordinary relations, also over
composite rings and prime powers. A singleton evaluates its one direct
relation; an absent partner outside the inclusive interval is not introduced.
For prime q=c±r, gcd(r,D)>1 implies q divides D, so tables retain coprime
offsets plus the prime divisors of D. This compact representation requires
explicit small-prime exceptions, retained in the implementation and oracle.

Pending centers survive prime-segment boundaries. Records contain their
original eligible primes and are consumed in finite GCD batches. If different
factors satisfy different primes within one trace, the trace can saturate;
recovery therefore replays each original q rather than only the trace. The
recovery allowance counts these attempts. Increasing B1 clears all retained
wheel/even-gap state, then repeats the complete new stage-two interval.
Equal-B1 B2 append preserves checked coverage; it need not pair across rungs.

The executor has its own checkpoint identity and deterministic work ledger;
it does not copy CADO, YAFU, Prime95 or PrMers code. The existing group-neutral
LCM ratio remains the only shared schedule interface. Future ECM point
continuation and allocation remain separate. The frozen follow-up protocol
and source control live in versioned `inputs/`; final acceptance still requires
correctness, fresh confirmation and complete-stage/portfolio measurements.
