# C3 continuing policy search: mechanisms beyond quick8

The first six-bundle study is an immutable negative result. Its 10 training
inputs selected quick8, which lost on 15 fresh inputs. That does not complete
the search for an economical ECM/SIQS policy. The first study underrepresented
intermediate smaller factors, did not fit stage-two depth jointly, and tested
little of the cost/probability scheduling described below. All revealed inputs
are now training/regression evidence. A later selection needs new confirmation.

## Cost and survival, rather than imported curve counts

For a surviving cofactor, let C_E be the cost of the next ECM block, q its
conditional chance of resolving the expensive search, C_Q the alternative
relation cost, and C_R the expected recursive cost after an ECM success.
The block is promising when C_E + q*C_R < q*C_Q. A complete model also includes
partial successes, finite-budget censoring and continuation on failure.
For illustration, a 0.02-second block that avoids a 4-second relation call
needs about 0.5% useful success before recursive cost; against a 0.2-second
relation call it needs about 10%. These are explanatory values, not fitted
v2 probabilities or a production rule. Costs must be CPU/wall observations;
ECM and SIQS accounting units cannot supply this comparison.

The success probability after failures is conditional: repeated failures
change the mixture of possible smaller factors. Known factor sizes can
stratify offline calibration; a dispatcher sees the cofactor's size, completed
attempts, prior expenditure and remaining resources. It never sees the hidden
factor, certificate or fixture ID.

## Additional primary and implementation evidence

Pinned downloaded bytes and failure receipts are in
[c3_round2_research_sources.json](../../inputs/controls/c3_round2_research_sources.json).
No upstream implementation is copied or executed.

- [Kleinjung, SHARCS 2006](https://www.hyperelliptic.org/tanja/SHARCS/talks06/kleinjung.pdf),
  slides 14–23, compares method sequences using expected time and success,
  including a p−1-then-MPQS example and a yield/time frontier. This is NFS
  cofactorization of two sifted residues, whose distribution and objective
  differ from complete factoring. The transferable idea is to measure method
  cost and conditional yield before choosing a sequence; the native timings,
  residue priors and smoothness-rejection rules do not transfer.
- [Kruppa's 2010 thesis](https://docnum.univ-lorraine.fr/public/SCD_T_2010_0054_KRUPPA.pdf),
  chapter 4 and conclusion on printed page 108, explicitly identifies the
  effect of failed attempts on subsequent factor-size expectations. This is
  motivation for recording surviving attempts, not a claim that the thesis
  supplies a ready-made optimal ECM/SIQS dispatcher.
- [CADO-NFS 692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b](https://github.com/cado-nfs/cado-nfs/tree/692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b),
  LGPL-2.1-or-later, `sieve/strategies/README`, `generate_strategies.cpp`,
  `finding_good_strategy.cpp` and `sieve/ecm/facul.cpp`: the strategy tooling
  benchmarks success by factor bits and cost by modulus size, weights later
  method costs by survival, and chooses a yield/time frontier using a workload
  distribution. `facul_all` selects a method sequence from observable residue
  sizes; `facul_aux` continues children at the next method. This is a strategy
  generator and executor, not evidence that every deployed default is freshly
  optimized. Native modular arithmetic, multiple curve families and the
  sifted-residue/smoothness objective differ from v2. Transfer measured cost,
  survival and completed-method credit; retain v2's exact validation and full
  factoring objective. Repeated p−1 attempts are correlated, so independent
  ECM-curve assumptions cannot silently be applied to all engines.
- [Alpertron 93c5c8189cb2f149c8dae898eb11997c5f2e7980](https://github.com/alpertron/calculators/blob/93c5c8189cb2f149c8dae898eb11997c5f2e7980/ecm.c),
  GPL-3.0-or-later header, `ecm.c:ecmCurve` and `factor.c`: a separate
  input-size table limits the number of ECM curves before `CHANGE_TO_SIQS`;
  the expected-curves/bound table alone does not establish dispatch. Its
  `NumberLength * 9` size proxy, C/WASM arithmetic and continuation differ.
  Transfer observable-size-dependent effort as a hypothesis, not the table
  entries or a claim of current v2 speed.

The already pinned YAFU caller adds current-cofactor target depth and completed
bound/count credit. It skips trivial ECM targets and reduces work for a usable
SNFS form; the latter has no matching v2 engine. Its tuned time estimates
select QS versus NFS. The inspected modern caller uses t-level ECM stopping;
older timed-ECM ratios must not be misrepresented as its current policy.
FLINT's current-cofactor goal and PARI's ordinary/deeper disjoint campaigns
remain useful residual-aware and explicit-mode mechanisms. Yamaquasi's eight
small curves were one transfer hypothesis, not the limit of this research.

## Mechanisms to implement and test

| Family | Expected benefit and failure mode | Bounded next experiment |
| --- | --- | --- |
| Cost-limited pretesting | Tie optional CPU/wall expenditure to an estimated relation cost, with a finite curve cap and protected fallback. A cheap relation engine warrants less investment; a costly fallback can justify more. A poor workload prior or noisy estimate can stop useful work. | Economic25/economic50 use historical 0.2/4.5-second relation estimates and cumulative 25%/50% ceilings, capped at 64 existing-bound curves. They are hypotheses, not upstream defaults. |
| Joint stage-two depth/count | Cheaper stage two may permit useful coverage at lower total cost; too short a continuation loses semismooth orders. | Compare 32 curves at B1/B2 2000/50000 with the 2000/147396 control. Keep accepted arithmetic, chaining and memory contracts. |
| Survival and finite escalation | Failed inexpensive attempts can justify a stronger independent tier when the remaining relation cost is high; a balanced workload may make escalation wasteful. | Compare 64 existing-bound curves and 16 existing-bound curves followed by four independent 11000/250000 curves. Record actual first-hit positions and per-curve costs in attribution only. |
| Residual-aware effort and prior credit | Recompute the appropriate finite effort for a smaller cofactor and skip already completed methods. Resetting budgets or replaying prior attempts creates false gains. | Use the broadened recursive/regression cohort, then fit a finite size/attempt decision table if the measured conditional frontier supports one. This requires a separately frozen implementation/identity before acceptance. |
| Reservation efficiency and funded admission | Policy overhead can erase small arithmetic gains. Unfunded reservations currently suppress useful cheap work under a 2M grant. | Pair unchanged 32-curve coverage with protected32, profile separately, and test default-grant behavior. Change the wrapper/admission only when attribution establishes a mechanism; freeze any changed source before training. |

Curve-family changes, scalar-chain searches and same-curve bound extension
would change arithmetic or reopen completed work. This search consumes the
accepted engines. A new result from those mechanisms is not assumed.
