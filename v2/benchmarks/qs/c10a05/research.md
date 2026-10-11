# C10a-05: collection boundaries and downstream pacing

Starting point: `master`, `16fe8a6180b3b86553cb2d2e56e7c27c9098e2ca`;
managed branch `codex/c10a-05-pacing`. The original checkout's uncommitted
roadmap, benchmark index and two gap reports were read as local research;
they are not part of this branch. v1 and historical experiments are unchanged.

## Executable diagnosis

`SIQSJob.run` invokes `QSJob.run(batch_limit=1)`. The latter passes
`next_position >= hi` to `_solve(final)`. In `_solve`, `final` bypasses both
`filter_row_growth` and the post-filter `row_excess` requirement. SIQS then
marks its polynomial inactive and activates another one, preserving the store
but clearing the engine's `final_solve_done`. Thus a changed polynomial window
is effectively an unconditional downstream retry even during shared collection.
An unchanged already-solved store is skipped. This matters when each polynomial
adds fewer rows than the existing filter cadence; it is not a new cadence-32
experiment. Both selected 40/50-digit configurations already set both fields
to 32. Historical R3's cadence-32 win at 30 digits remains a separate result.

The nine retained C10 fresh `plan` receipts in the original checkout were
summed independently. At 50 digits: 1169.093374 s complete wall, 751.799106 s
collection, 185.076547 s preparation, 84.596251 s filtering, 84.809477 s
elimination, 6203 filters and 6194 solves. The three downstream stages sum to
354.482275 s (30.321%). At 40 digits they sum to 4.188817 / 61.289816 s
(6.834%). The old follow-up ledger's 85.020/83.518 values are not these
selected-plan sums; its historical text remains intact. Receipts and hashes
are recorded separately in the new source manifest. These are opportunities
for attribution, not predicted savings.

New bounded, instrumented control calls on the first previously revealed C10
input per band find 52/640 forced polynomial-boundary solves at 40/50 digits;
all have fewer than 32 newly retained rows. Their first nonempty solver output
also produces a proper factor at the last observed attempt. These observations
bound what was seen; they do not prove that no untested earlier prefix had a
useful kernel. The detailed event traces retain post-filter columns, rows,
dependencies, growth, positions and factor witnesses. They are not timing
acceptance samples, nor additional independent numbers.

## Primary source reconciliation

The [manifest](../../inputs/controls/c10a05_sources.json) records immutable
revisions, archive/file hashes, licenses and inspected anchors. Local captures
were checked against the prior C9/C10 pins; executable caller reads were
repeated for this task. A targeted web search additionally revisited official
msieve, the papers below and implementer notes. No upstream executable was
timed, no source code was copied, and this is not exhaustive SOTA coverage.

| Source | What the inspected caller does | Transfer boundary |
| --- | --- | --- |
| [YAFU, `8110dfb`](https://github.com/bbuhrow/yafu/blob/8110dfbd8c6f9486d93b1a02de6eb7b180e55a80/factor/qs/SIQS.c) `SIQS.c:172, 640, 1023–1181, 2853` | Initializes base-plus-extra target, collects across polynomial sets, filters, solves and extracts. Distinct branches for solver errors −1/−2 and no factors increase the target by 5%; retry counter is finite. | Target feedback is C9a. Do not equate insufficient matrix yield with nontrivial-kernel extraction failure, or import 5%/64 as PyPy optima. Public-domain notices in inspected files. |
| [Official msieve 1.53](https://sourceforge.net/projects/msieve/files/msieve/Msieve%20v1.53/) `mpqs/sieve.c:427–487, 580–680`; `mpqs/mpqs.c:97–125` | Uses base count plus 1.5 times its native excess; collection progress counts full relations plus graph cycle rank, not every raw partial. Collection loops across polynomial blocks; after the target and flushed savefile, filtering runs, followed by one solve/extraction in `mpqs`. | This caller has no YAFU-style target-retry loop. `Readme.qs` explains restart/duplicate records can falsely satisfy the collection estimate. Public domain; dependencies separate. |
| [FLINT, `a4c9750`](https://github.com/flintlib/flint/blob/a4c9750d0d3d67bb01cf6d18c187591b313451c3/src/qsieve/factor.c#L284-L418) `qsieve_factor_with_tune`; `large_prime_variant.c:457–567` | Raw full-plus-cycle target precedes relation processing. Processing can report too few useful relations and return to collection (adjusting its edge estimate). Successful processing precedes matrix reduction, repeated Lanczos and square roots; failure or A exhaustion can enlarge the base. | Distinguish collection estimates, processed yield and failed factor extraction. Base growth and the unbounded native Lanczos retry are not transferred. LGPL-3.0-or-later. |
| [Yamaquasi, `3f95f43`](https://github.com/remyoudompheng/yamaquasi/blob/3f95f43682ed15d8c1ed206a9a702dd655d7c8ad/src/siqs.rs#L119-L234) `siqs`, `sieve_a`, `SieveSIQS::new`; `relations.rs:125–153` | Polynomial end checks the target before scanning relation coverage. `gap` counts distinct odd-exponent primes and asks for 48 extra cycles; target starts at 0.8 of base count and is reset using measured gap. Final processing follows collection, with finite A population and an insufficient-yield panic. | `gap` is a coverage estimate, not recursive singleton-filter surplus or measured rank. Feedback and truncation assumptions do not transfer. BSD-3-Clause. The later C9-checked `181e282` change is not substituted for this immutable pin. |
| [JavaMath, `088d01f`](https://github.com/TilmanNeumann/java-math-library/blob/088d01fa97e7d0412c6cbbb8f463fa7e5f78ab97/src/main/java/de/tilman_neumann/jml/factor/base/congruence/CongruenceCollector01.java#L99-L208) `initialize`, `collectAndProcessAQPair`, `addSmooth:351`; `SIQS.java:265–301` | SIQS repeatedly generates/sieves polynomials; collector tests exact squares immediately. Otherwise it waits for base count plus fixed extra congruences; no factor advances its target by the extra amount (default 10). | Motivates preserving early square opportunities. Yield-triggered target changes remain C9a. We retain original rows, including trivial squares, unlike its discard path. GPL-3.0-or-later. |
| [Alpertron, `93c5c81`](https://github.com/alpertron/calculators/blob/93c5c8189cb2f149c8dae898eb11997c5f2e7980/siqs.c#L2342-L2345) `FactoringSIQS`, `InsertNewRelation:2644–2765`, `sieveThread:2816–2900`; `siqsLA.c:LinearAlgebraPhase` | Matrix target is max(base×33/32, base+50). Polynomial-set completion enters linear algebra only after that relation target; `while (!LinearAlgebraPhase(...))` retries downstream. Admission excludes empty parity rows and duplicate matrix rows. | Neither that loss policy nor unbounded retry is acceptable here. Native fixed arrays and 32-vector layout do not choose a PyPy target. GPL-3.0-or-later. |

## Literature and implementer notes

- [Ekkelkamp, *On the Amount of Sieving in Factorization Methods*](https://math.leidenuniv.nl/scripties/EkkelkampThesis.pdf),
  §4.3 and §§4.6.3–4.6.4: stopping on post-singleton surplus differs from
  raw relation volume. Additional collection changes matrix choice and weight;
  the smallest relation count need not give the shortest total time. These
  particular experiments are NFS, not SIQS constants or PyPy speed predictions.
- [Cavallar, *Strategies in Filtering in the Number Field Sieve*](https://ir.cwi.nl/pub/2078/2078D.pdf),
  §§2.1–2.2: pruning and merging consume surplus to trade matrix size against
  weight. A singleton excluded from a completed matrix need not remain a
  singleton after more collection. Our original store therefore stays intact.
- [Boender–te Riele, *Factoring Integers with Large-Prime Variations of the Quadratic Sieve*](https://ir.cwi.nl/pub/1367/1367D.pdf),
  §§5–7, and [Leyland et al., *MPQS with Three Large Primes*](https://infoscience.epfl.ch/bitstreams/cf37274e-b9f3-4986-8615-c751dab46f2e/download),
  §§2.3, 3, as mapped in C9: partial/cycle yield and the eventual matrix matter;
  no universal cycle growth law or upstream large-input constant is assumed.
- Official msieve's `Readme.qs` distributed/restart discussion was read with
  its executable target caller. [Percival's implementer note](https://www.daemonology.net/blog/2005-11-30-quadratic-sieve-constant.html)
  was revisited as a constant-cost motivation; it supplies no applicable
  fixed retry policy. Paper copyrights remain with their authors/publishers;
  only links, hashes and paraphrased findings are committed.

**Mathematical argument:** for m retained parity rows on c live columns,
nullity is m−rank, at least m−c. Positive surplus is sufficient for kernels,
not necessary; a zero row is already a kernel, and a kernel may give only
trivial GCDs. Provenance lifting and proper-divisor tests remain necessary.

**Engineering inference:** preserve the current cadence and surplus at
shared-store polynomial boundaries, while keeping immediate zero-parity and
genuine terminal opportunities. This changes scheduling alone. An unconditional
base-count target could miss early useful factors; adaptive yield/rate estimates
would cross into C9a. No extra fixed constants are justified by this diagnosis.
The selected small set is therefore the current `window` policy and one
`collection` challenger. Measurements must price its additional collection.

## C9a interface follow-up

The new pending terminal boundary is independent of any future target
controller. A later C9a controller would need separate observations for
raw atoms, retained full/cycle rows, post-filter live rows/columns, solver
kernel count, tested trivial dependencies and proper factors; it must not
interpret all failed downstream attempts as insufficient yield. Persist its
target and observation epoch, charge estimation/rebuild, and retain the
terminal drain defined here. No feedback, residual-domain or graph change is
implemented or accepted by this tranche.
