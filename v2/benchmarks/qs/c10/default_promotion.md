# C10 bounded tranche: scoped automatic plan selection

The user requested automatic use of the accepted C10 plan in the size bands
where it helped. This is a routing decision after the frozen first-tranche
experiment, not a new tuning search or a claim of optimal collection settings.
The [frozen controls](../../inputs/controls/c10_training.json), [fresh result](tranche1_results.md)
and [candidate reconciliation](candidate_reconciliation.md) retain their
original identities and verdicts.

The selected plan adds exactly `power_plan_bytes=2**21` to the 40- and
50-digit growing-DLP bundles, leaving the factor base, interval, A family,
residual admission, storage and stopping fields fixed. In nine matched fresh
complete calls per band, it reduced wall cost by 21.83% at 40 digits and
17.79% at 50 digits, with 9/9 validated completions in each arm and band.
The measured inputs were three independent certified balanced products per
band, conditioned on a large factor of `p-1`; arbitrary inputs of the same
digit count were not measured. Candidate, relation and final dependency
counts were identical to control in these comparisons. Owned workspace rose
by about 2 MiB; process RSS was a shared high-water mark.

At 40 digits C1 accepted its **fixed-store DLP** bundle on a different
balanced cohort (31.73% less complete cost than calibrated SLP). C9's later
growth-versus-fixed-DLP timing was unstable and showed no new growth gain;
both collectors did identical arithmetic work in that band. Its descriptive
growth-versus-SLP estimate was 31.72% lower complete cost. C10 then confirmed
the plan on top of the growing bundle. Automatic 40-digit routing is therefore
a user-directed composition of those bounded results, not a claim that growing
storage itself improved 40-digit timing. C9 accepted the 50-digit growing-DLP
bundle on its fresh balanced cohort
(38.10% less complete cost than calibrated SLP). These are separate matched
comparisons; their percentages must not be multiplied or added into a claimed
end-to-end default speedup. The 30-digit C10 plan did not pass selection and
remains off. No result covers 41–49 digits, other arithmetic backends, changed
score policies, lower owned memory, SSS workers or a combined portfolio policy.

The automatic scope is the Python-integer 40- and 50-digit bands with the
confirmed 1 GiB owned SIQS allowance. Direct `SIQSJob` creation, omitted
portfolio SIQS and no-override CLI SIQS/automatic calls use the same concrete
selection. Explicit `SIQSConfig`, QS flags, `siqs=None`, custom budget or
memory allowances and saved checkpoint configurations retain precedence.
The finite 10^13-work, 1,200-second wall/CPU and 1 GiB SIQS grant mirrors
the confirmation envelope; the outer portfolio reserves additional
coexistence headroom. Lower grants can return unresolved cofactors without
changing the correctness or resume contract. `SIQSConfig()` remains the
legacy reference configuration with `power_plan_bytes=0`.

The [C10a gates](c10a_followups.md) hold the distinct untested mechanisms,
including normalized/mod-8 polynomials, broader A-product and diversity
selection, per-job threshold feedback, coupled size policies and solve
cadence. The failed fixed 16-family ordering is a verdict on that arm alone.
A1a, C8, C9a, the matrix/solver owners and E1 retain their own boundaries.

The [integration validation](validation.md#subsequent-c9c10-default-integration-verification)
records 721 passing committed-only PyPy/gmpy2 tests, lint, 144 live benchmark
imports and complete default-path factoring of one certified fresh input per
band. These functional calls do not enlarge the frozen timing population.
