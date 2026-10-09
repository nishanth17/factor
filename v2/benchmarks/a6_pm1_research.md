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
