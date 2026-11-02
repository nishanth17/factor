# C10a-01 revised combined-arm nomination

Decision: 2 November 2026. The four extensions are now implemented and
independently validated; see the [implementation handoff](c10a01_coupled_implementation.md).
The original nomination rationale follows. The earlier v2 candidate and its
evidence remain unchanged. The revised bundle has not been timed and cannot
execute under the old source freeze.

## One combined challenger

C10a-01 tests a coupled polynomial-population and multiplier policy. Research
findings that make that policy coherent belong in the same experimental arm.
They do not need independent speedup evidence before they can be tested
jointly. The earlier argument that an unmeasured component must stay outside
the arm applied an adoption requirement to experiment admission and was wrong.
A changed bundle requires a new source/protocol identity; an existing freeze
is historical evidence, not a reason to exclude a justified improvement.

| Component | Membership and concrete proposed change |
| --- | --- |
| Exact parity/discriminant representation | Retain 2P, the exact original relation identity, the D=1 mod 4 eligibility rule and explicit standard fallback for even D or D=3 mod 4. |
| Wider representation-aware multiplier selection | Include all 61 squarefree k from 1 through 97, including 20 even choices, and the first 300 primes through 1987. Keep exact derived characteristic-two and singular-prime weights, the half-log size penalty, deterministic ties and proper-GCD recovery. These finite limits are a source-motivated experiment specification, not a tuned optimum. |
| Reuse of multiplier residue calculations | Include bounded reuse within selection so the wider search does not repeat N's modular/Legendre calculations for every k. Reuse into base setup only if identities and ownership are proved. A separate performance arm is not necessary for this supporting implementation. |
| Normalized A target and actual-product construction | Include target-aware product selection, not only the already implemented halved odd-a target. Extend the existing one-flyer construction with a bounded two-compensator population, choosing the closest actual product from disjoint low/high prime sets around the core pool. This addresses actual target misses without repeating the retired 16-family reorder or quota/face grids. |
| Guaranteed powers of two | Include a specialized exact score initialization for normalized polynomials: stored F=2P has at least two factors of 2 for D=1 mod 8 and exactly one for D=5 mod 8. Preserve hit metadata and full relation exponents; account for higher powers exactly once. This is part of implementing the representation efficiently. |

The old literal standard control stays unchanged. The revised challenger is
a bundle comparison: any improvement would initially be attributed to that
bundle. Small component checks and ablations can explain it; component-level
attribution is not an admission prerequisite for a useful complete-factor
comparison. The proposed third standard-k=31 arm is optional explanatory work,
not the main answer to which mechanisms belong in the challenger.

The wider selector must be evaluated as its declared policy at all four sizes.
Do not assume it chooses the old k values, carry over the forced k=31 label,
or choose a new k after complete-factor outcomes. If it selects standard
fallback, report that result honestly and distinguish any separately declared
eligible-form diagnostic. Both residue-1 and residue-5 correctness support
remain; this nomination does not silently impose a mod-8-only constraint.

## A-product scope and obligations

The existing flyer already chooses the nearest available single prime for
each core product. Its disjoint core pool can exclude precisely the prime
range needed for a close completion. Sorting the same assignments again is
not the missing change. The proposed population uses r-2 core factors and
one compensator below and one above the core pool. All factors must be
nonsingular, distinct, odd factor-base primes. Disjoint sets make the core
recoverable from A's factorization, preserving uniqueness without an
unbounded seen set. A monotone two-pointer scan can minimize absolute product
error over the finite low/high Cartesian product without storing all pairs.

Before implementation admission, independently verify that minimum against
Cartesian enumeration, prove deterministic ties and finite exhaustion, and
bind the new population to resume/cache identities. Charge scans and rejected
or impossible constructions, including interrupted reconstruction. Establish
reachability for the literal 60/62/65/70 configurations and measure setup cost
as part of complete cost. Do not claim close A products guarantee useful
relations. Generic A diversification and alternative mean-norm objectives
remain separate; this change addresses the target population of this arm.

## Kept outside this bundle for specific reasons

Batch smoothness is a credible separate collector experiment: it changes
candidate-processing strategy and simultaneous radical/tree/pending-candidate
ownership, and the current SmoothBatch base-prime ceiling does not cover the
larger declared bases. Its benefit is not specific to normalization. That is
a concrete integration boundary, not a requirement to prove a speedup before
testing it. Relaxed admission, 3LP, new solvers, cross-polynomial batching and
hardware-specific lanes also remain separate changes.

Direct P storage enumerates the same trials and residuals as 2P. Keep 2P and
optimize its known characteristic-two work first, preserving the proven
extraction contract. Do not add native empirical weights or an alternate
continuous A target merely because another implementation uses them.

## Evidence and execution status

The membership decision draws on the retained expanded review and actual
Factor callers. Detailed provenance remains in local audit research. Existing
calibration is evidence for the v2 prototype only; its 3,124 charged seconds,
failures, source hashes and archives are unchanged. No extra calibration or
complete-factor pilot was run for this decision.

The implementation belongs to this same C10a-01 challenger, with schema 9,
a new source identity, independent oracles and resource accounting. After
bounded new calibration, freeze one current control
and this combined challenger at 60/62/65/70. Preserve the funded two-arm
complete-factor envelope as the design constraint; adding components does
not itself multiply arms or require the optional 42,600-second C3 proposal.
Re-derive work/setup reservations and reconfirm matched memory tiers before
execution. The old 35,200-second plan is not evidence that this revised bundle
has passed those checks, and no new numerical allocation is implied here.


The [combined v3 implementation](../studies/c10a01_coupled_implementation.md)
adds wider coupled ranking/residue reuse, paired A completion and guaranteed
powers-of-two initialization. `configuration_v3.py` constructs literal matched
control/challenger configurations. This bundle is correctness-validated but
unmeasured. Earlier calibration and source freezes remain historical; the old
pilot loader rejects changed sources. No new performance call is admitted by
the configuration helper or compatibility-manifest update.
