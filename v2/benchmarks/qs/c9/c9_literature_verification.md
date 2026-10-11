# C9.2 literature-transfer verification

Verified 11 October 2026, after the accepted C9.2 study. **C9.2 adapts a
selected, demonstrated subset of the literature. It does not adapt every
relevant mechanism or establish an optimal large-prime policy.** The
[accepted result](c9_growth_results.md) remains valid for its stated scope.
This audit records remaining transfers explicitly; it changes no runtime,
collector configuration, corpus, frozen protocol or timing verdict.

## What was rechecked

Re-read the residual, graph and downstream callers in YAFU, msieve,
Yamaquasi, JavaMath and FLINT against v2's [intake](../../../qs/large_prime_collector.py),
[forest](../../../qs/large_primes.py), [pipeline](../../../qs/pipeline.py)
and [SIQS caller](../../../qs/siqs.py). All **43 previously pinned inspected-file hashes** still match:
39 in `c9_source_pins.json`, four in `c9_growth_source_pins.json`.
GitHub default-branch HEADs still match the earlier pins except Yamaquasi.
Its new `181e282cbfcace30bd0d99c2b4bae181dee7ba69` adds three commits dated
11 October: Python conversion, warning cleanup and build/dependency updates.
The inspected source delta changes no LP admission, splitter, retention,
cycle or relation-target algorithm. This is a dated check, not a promise that
upstream HEAD will remain unchanged. The historical pins remain intact.

The versioned [verification manifest](../../inputs/controls/c9_verification_sources.json)
records these revisions, additional inspected files and primary-paper hashes.
Source archives, API receipts, PDFs and extracted text remain in ignored
`benchmarks/results/c9/literature-verification/` or the original local research
archive. Mathematical ideas were independently implemented; no upstream code
was copied. The original [comparison](c9_research.md) records licenses,
backends, certainty, units and size-dependent branches. Newly inspected
YAFU `batch_factor.c` has a public-domain notice; `microecm.c` is dual-licensed
under its FreeBSD/BSD-style terms or MPL-2.0. Yamaquasi `params.rs` uses its
BSD-3-Clause license. Their file hashes and terms are in the new manifest.

The primary cross-check includes [Lenstra–Manasse §§1–3](https://www.researchgate.net/publication/266349860_Factoring_with_Two_Large_Primes),
[Boender–te Riele §§5–7](https://ir.cwi.nl/pub/1367/1367D.pdf),
[Cavallar §§2.1–2.2](https://ir.cwi.nl/pub/2078/2078D.pdf) and
[Gower–Wagstaff §§3–5](https://homes.cerias.purdue.edu/~ssw/squfof.pdf). It also adds substantive
inspection of:

- [Ekkelkamp, *On the Amount of Sieving in Factorization Methods*](https://math.leidenuniv.nl/scripties/EkkelkampThesis.pdf),
  §§3.3.2–3.3.3: unequal bounds on the larger and smaller residual primes,
  with fixed product allowance; highest relation density and lowest time
  need not coincide. §§4.3 and 4.6.3–4.6.4 study post-singleton surplus,
  calibrated progress and extra collection versus matrix cost. The latter
  experiments are NFS-specific, not QS presets. **The asymmetric-domain
  mechanism was missing from C9's explicit transfer map.**
- [Bernstein, *How to find smooth parts of integers*, 10 May 2004 draft](https://cr.yp.to/factorization/smoothparts-20040510.pdf),
  Algorithm 2.1/Theorems 2.2–2.3 and §3: product/remainder trees amortize
  smooth-part extraction across inputs. Repeated squaring recovers prime
  powers; a single GCD with a square-free prime product is insufficient.
  This is a batch arithmetic algorithm, not a free replacement for a rho
  call or a measured PyPy benefit. The previous report mentioned batching
  without this algorithm-to-contract account.
- [Leyland et al., *MPQS with Three Large Primes*, ANTS 2002](https://infoscience.epfl.ch/bitstreams/cf37274e-b9f3-4986-8615-c751dab46f2e/download),
  §§2.3 and 3: cycle density affects downstream work; extra collection can
  allow selection among less dense cycles. Their growth data do not fit
  one universal power law. The 135-digit/TLP result cannot select v2's
  durations, cycle cap or DLP thresholds. Higher-arity implementation
  remains F6.

Author notes were also retraced: msieve's `Readme.qs` warns that duplicate
savefile records can falsely satisfy the collection target; Yamaquasi's
comments distinguish score targets from intake; [YAFU's batch source][yafu-batch] cites
Gerbicz's remainder-tree discussion and its memory tradeoff. The operative
caller, rather than a blog or comment alone, determines the policy below.

## Mechanism coverage and remaining gaps

B is v2's numeric factor-base bound; L is an endpoint ceiling. A product
limit, score allowance and endpoint limit are separate domains. “Existing”
means the capability exists; it does not mean all its settings are calibrated.

| Source mechanism | v2 and C9.2 disposition | Remaining issue / transfer cost |
| --- | --- | --- |
| Grow original history as well as graph tables — YAFU `tdiv.c:1117`, `filter.c:1177`; msieve `relation.c`; FLINT `large_prime_variant.c` | **Implemented and adopted:** allocator-grown partial/atom/row containers, no implicit FIFO/count ceiling in growing mode, explicit caller quotas and checked schema-five recovery. | Native packed arrays and PyPy reservations differ. Growing storage consumes more memory; disk spill and a tighter representation remain unimplemented. |
| Common SLP/DLP endpoint — YAFU and msieve intake | **Evaluated individually in training, adopted in the selected bundle:** SLP ceiling 3,000,000 instead of 900,000,000 at the same DLP product/score domain. | Deliberately loses rare high SLP matches. The fresh bundle result does not prove a universally optimal ceiling or separately quantify the causal share of alignment. |
| Nominal score target with wider recovered intake — [Yamaquasi `fbase::cofactor`][yama-fbase] | **Existing interface, untested policy:** `candidate_bound` can differ from `large_product_bound`; current selected values resolve to the same 115,200,000,000 allowance. | Yamaquasi can retain recovered products up to L² when endpoints fit. This transfer was documented but not evaluated. It can add useful atoms, splitting, certification and storage cost. |
| Unequal ordered endpoint bounds — Ekkelkamp §§3.3.2–3.3.3 | **Previously unmapped; missing exact policy:** v2 has a common endpoint cap and hyperbolic product cap, not an independent cap on the smaller endpoint. | A sorted p≤q policy with p≤L₂, q≤L₁ is not equivalent to pq≤P, q≤L. This is a C9 residual-policy candidate requiring calibration and an execution/checkpoint identity if implemented. |
| Size-dependent residual effort / fallback — YAFU micro-ECM, msieve SQUFOF, JavaMath Hart/ECM/rho/nested SIQS, Yamaquasi rho/ECM | **Deferred algorithm transfer:** v2 retains exact-square recovery, two bounded Brent walks and proven endpoint checks. Growing mode removes the unrelated cumulative split-count quota. | Near-zero current-domain failures support retaining success coverage, not declaring the splitter fastest. Failed walks, cost and rejected larger products require representative measurement before a new schedule or splitter. |
| k−1 independent equal-endpoint combinations | **Already supported and evaluated:** C1's graph-SLP route reuses a representative; ordinary SLP uses disjoint pairs. Growing graph-SLP was screened, not selected. | Two training inputs with process/JIT confounding do not reject reusable SLP generally. No new fresh graph-SLP win was established. |
| Independent partial basis / duplicate handling — JavaMath Collector01; msieve duplicate purge | **Already partly supported:** DLP publishes fundamental cycles and retains independent connectivity edges plus exact original provenance; repeated original IDs are rejected. | `duplicates=0` measures repeated original identity, not every algebraically equivalent row. No separate census of equivalent norm/square-root content was done. Equal parity alone does not make a row safe to discard. |
| Breadth-first cycle construction, short representatives and sparse-first cycle selection — Lenstra–Manasse §2; [msieve `relation.c`][msieve-relation]; [Yamaquasi `combine_single`/`combine_double`][yama-relations] | **Correctness supported, representation optimization deferred:** v2 handles every component, loops and parallel edges; union-by-size rerooting does not minimize cycle length or factor-base row weight. | Zero 256-atom-cap refusals does not show shorter cycles have no value. They could reduce provenance, verification and matrix fill-in. Replacement/rebuilding must preserve the original-atom kernel, ownership and checked restore. |
| Singleton pruning / merging — Cavallar; native filters | **Already supported for the current factor-base matrix:** recursive singletons, weight-two merging, lifted dependencies and exact reconstruction. | LP leaves must not be permanently deleted merely because a prefix has no cycle. Higher-order merging, clique selection and mixing LP columns into the final matrix remain matrix/provenance-owner work. |
| Downstream-feedback targets and amortized solve cadence — [YAFU retry caller][yafu-siqs], [JavaMath Collector01][java-collector], [FLINT caller][flint-factor] | **Partial support; important untested transfer:** v2 continues after trivial dependencies, but selected growth/excess 32 is overridden by changed-store polynomial-window finalization. | It lacks the source-style target/pacing policy. YAFU raises its target 5% on failed downstream attempts with a finite retry loop; JavaMath adds 10 after no factor. Those constants are not calibrated PyPy values. Scheduling is coordinated with C10a; yield-dependent target feedback remains C9. |
| Batch smooth-part extraction / cofactoring | **Missing algorithm, lower priority at the current measured domain:** no product/remainder-tree batch intake. | [YAFU `tdiv.c:341`][yafu-tdiv] disables its DLP-only batch branch; higher-arity branches differ. A PyPy transfer must account for tree setup, simultaneous storage, repeated powers, interruption and certification. It is not established SOTA DLP behavior simply because batching code exists. |
| Higher arity, probable-prime coverage shortcuts, opaque equal-cofactor matches | **Explicitly not transferred:** arity >2 remains F6; opaque cofactors require a new compatible relation contract; current endpoints stay proven. | FLINT's opaque-equality cancellation is exact in its contract, but is not v2's prime endpoint representation. Yamaquasi's disconnected-cycle omissions and all-even-row discard, or a one-base PRP loss, must not silently narrow supported v2 cases. |

Concrete recovered-product example: with the selected B=30,000 and
L=3,000,000, L² is 9,000,000,000,000, **78.125 times** the selected product
allowance. An exact residual 1,000,003×1,000,033 = 1,000,036,000,099 has two
endpoints below L, but current intake rejects it on product size. If candidate
recovery happened to report it under the nominal score target, a wider intake
could retain it. This is an arithmetic opportunity, not evidence that such
norms occur often, split cheaply or produce an earlier proper factor.

Likewise, a product cap allows balanced endpoint pairs that an explicit
smaller-endpoint cap can exclude. An asymmetric experiment therefore cannot
be described as merely changing `large_product_bound`. Both admitted and
rejected residual populations, including late prefixes, must be sampled.

[YAFU's micro-ECM caller and dispatch][yafu-microecm] pass `is_arbitrary=0`; the generic preliminary
small-factor curves are therefore not this QS path. Its scalar dispatch tiers
residual bit length and uses B2=25×B1. Optional p−1 calls are disabled by the
local macro. Native Montgomery/SIMD arithmetic and probable-prime yield
assumptions do not establish that ECM is cheaper than proven PyPy rho here.

## Revised investment order

1. **Coordinate solve scheduling and yield-dependent retry targets.** The
   accepted growth calls spent 157.662 seconds in preparation, 74.168 in
   filtering and 70.100 in elimination across nine calls, with 5,411 filters
   and 5,402 solves. Their approximately 301.93 seconds of downstream work is
   not wholly avoidable. Nevertheless, “sieving dominates” does not justify
   dismissing this transfer. Compare a frozen scheduling/feedback policy
   with unchanged residual collection; preserve immediate supported square
   dependencies, true final/storage-cap extraction, cancellation and charged
   continuation. Complete cost can worsen if useful early extraction is delayed.
2. **Attribute recovered-product and ordered-endpoint policies.** These are
   distinct hypotheses. First reuse retained atoms and representative rejected
   norms; then use a separately frozen small policy study only if it predicts
   useful dependency value. Keep candidate scoring fixed where possible;
   any changed recovery is a coordinated C9/C10 interface. Reject a transfer
   on complete cost or resource regression, not raw admission count.
3. **Measure cycle/provenance weight and algebraic duplication.** Replay the
   same original atoms through independent small-graph kernel and arithmetic
   oracles before trying shorter representations. Measure total path length,
   verified row support, post-filter kernels and extraction; no claimed speedup
   follows from shorter cycles alone. Matrix/provenance redesign retains its gates.
4. **Revisit split effort, batching or disk retention when the population or
   resource regime changes.** The instrumented 50-second training profiles
   measured only 0.31/0.53 seconds in splitting for fixed/aligned DLP, nested
   within 1.47/1.33 seconds in residual processing. Fresh fixed DLP subsequently
   had one failed split in 249,354 attempts, growth zero in 106,686. Those are
   observations at their admitted domains; they do not cover rejected larger
   products or prove that upstream effort schedules cannot improve cost.

This order is an inference from source mechanisms and v2 evidence, not a
newly timed selection. No additional experiment allowance is silently created.
C10's separate frozen training retains the shared machine; this verification
uses source/document reads and metadata checks, not heavy QA or timing.

## Acceptance and verification boundary

The accepted growing-DLP bundle still demonstrates 38.10% lower complete-call
wall cost than calibrated SLP and 53.08% lower than fixed DLP on three fresh
balanced 50-digit inputs, nine matched groups and 9/9 valid completions per
arm. Its higher memory and charged-work costs, generation bias and limited
independent-input count remain explicit in the original report. This audit
neither supplies new timing evidence nor changes the 30/40/default decisions.

All **209 accepted frozen source digests** match after this audit. Runtime,
benchmark runners, accepted timing controls and confirmation corpus are unchanged. The prior
596-test, lint and committed-only test/import/loaders receipt remains the
applicable QA evidence; the documentation/manifest change does not claim a
new test run. JSON parsing, local-link and hash checks pass, as does `git diff --check`
for this addendum; no heavy QA was run alongside C10's timings.

Only C9.2's named accepted tranche is closed. The remaining C9 transfers above
are open hypotheses or explicitly gated work. The research omission on
asymmetric domains and the incomplete attribution of cadence/row weight are
now visible; “all literature adapted” and “most optimal strategy” would be
unsupported descriptions of this task.

[yama-fbase]: https://github.com/remyoudompheng/yamaquasi/blob/3f95f43682ed15d8c1ed206a9a702dd655d7c8ad/src/fbase.rs#L480-L567
[yama-relations]: https://github.com/remyoudompheng/yamaquasi/blob/3f95f43682ed15d8c1ed206a9a702dd655d7c8ad/src/relations.rs#L440-L520
[msieve-relation]: https://github.com/MersenneForum/msieve/blob/03e0ab9e13c4862e5bb2c51a4a39718aa2500e19/mpqs/relation.c#L1089-L1190
[yafu-siqs]: https://github.com/bbuhrow/yafu/blob/8110dfbd8c6f9486d93b1a02de6eb7b180e55a80/factor/qs/SIQS.c#L1037-L1181
[yafu-tdiv]: https://github.com/bbuhrow/yafu/blob/8110dfbd8c6f9486d93b1a02de6eb7b180e55a80/factor/qs/tdiv.c#L326-L341
[java-collector]: https://github.com/TilmanNeumann/java-math-library/blob/088d01fa97e7d0412c6cbbb8f463fa7e5f78ab97/src/main/java/de/tilman_neumann/jml/factor/base/congruence/CongruenceCollector01.java#L186-L208
[flint-factor]: https://github.com/flintlib/flint/blob/a4c9750d0d3d67bb01cf6d18c187591b313451c3/src/qsieve/factor.c#L284-L325
[yafu-batch]: https://github.com/bbuhrow/yafu/blob/8110dfbd8c6f9486d93b1a02de6eb7b180e55a80/factor/batch_factor.c#L40-L77
[yafu-microecm]: https://github.com/bbuhrow/yafu/blob/8110dfbd8c6f9486d93b1a02de6eb7b180e55a80/factor/gmp-ecm/microecm.c#L3594-L3705
