# C1 source-informed larger-input feasibility follow-up

This is a new, bounded protocol authorized by the user's challenge to the
initial screen. Preserve `c1_protocol.md`, its controls and all captures.
The initial screen did not establish DLP infeasibility. Its 40-digit matrix
refusal was a diagnostic representation limit; 512 blocks at 60 digits were
insufficient to assess delayed cycle formation. No completion is claimed.

## Research-driven changes

The pinned Yamaquasi `src/siqs.rs` separates endpoint and product caps and uses
about 100 B² for products near 60 digits, rather than the first screen's
10,000 B². Lenstra–Manasse discusses nonlinear, delayed cycle accumulation,
asymmetric limits and the expense of false candidates/cofactoring. A finite
longer collection with tighter products is justified. These are hypotheses
for this PyPy implementation, not imported native speedups or digit cutoffs.

First analyze every retained initial-screen record using an offline spanning
forest and complete fundamental-cycle basis. Include components disconnected
from vertex 1, loops and parallel edges. Independently compare with the
existing generic GF(2) incidence oracle on small multigraphs. This is a
benchmark-only analyzer, not a production graph collector. It removes the
unnecessary dense large-prime matrix while retaining worst-case reservations
for actual cycle/provenance and factor-base matrices. Do not weaken R3's
production matrix reservation. Six offline cells, at most 30 wall/CPU seconds
and 512 MiB per analysis, 180 seconds total.

## Fixed workloads and candidate domain

Use two independent already-certified training inputs in each band 40, 50
and 60 digits. At 40/60 digits use B1's balanced training input plus index 0
from `phase_three_p34_large_v3_corpus.json`; at 50 use indices 0 and 1 from
that corpus. All remain training; no previously inspected input becomes
fresh confirmation. Use seeds 7 and 29 respectively, not two purportedly
independent observations of the same polynomial prefix.

SLP uses B1's selected 40-digit configuration, a declared uncalibrated
50-digit interpolation (B=30,000, half-width=65,536, six A factors, flyer),
and B1's 60-digit diagnostic configuration. The latter two are feasibility
controls, not promoted complete-factor presets. Use 8,192 rows/partials and
32,768 atoms at 50/60; retain all original 40-digit limits. Production memory
is 256 MiB; production work is 10^13 units per run.

Audit the unchanged SLP collector with endpoint cap 100B and product cap
128 B². A nested 64 B² cohort uses exactly the same positions and records;
its splitting costs must be measured separately, not inferred from edge
counts. SLP retains its B² prime bound. Do not assume cofactor <B² is prime:
classify every admitted endpoint PROVEN under the existing certainty contract.
Use the existing two-attempt, 2,048-evaluation Brent splitter, batch 32,
recovery 64, exact-square shortcut, exact proper-divisor and product checks.
No splitter or polynomial-parameter search in this follow-up.

## Finite allowances and representative sampling

Six uninstrumented SLP runs: at most 30/120/300 wall and CPU seconds per
40/50/60-digit input. Six diagnostic runs: at most 180/480/900 wall and CPU
seconds, 131,072 blocks and 10^13 production work units each. Diagnostic
classification/splitting/recovery has a separate 10^13 work allowance within
the same total wall/CPU deadline; at most 131,072 splitting calls per run.
At most 65,536 exact records, 256 MiB conservative record reservation, and
64 MiB serialized captures per run. Entire follow-up retains at most 512 MiB
on disk and consumes at most 4,500 active wall/CPU seconds, including offline
analysis. Stop immediately on any exhausted allowance and report censoring.
Do not silently widen these limits or add workloads/bounds.

Every position has independent sampling probability 1/256, implemented by
geometric skipping. This gives equal inclusion probability to short polynomial
tails and full blocks. Full factor-base trial division independently checks
sampled residuals, including block/refinement rejects. Retain at most 128
reservoir samples per rejection/shape stratum and count their populations;
retain at most 32,768 sampled positions per run, then stop sampling and mark
its coverage prefix explicitly. Candidate census continues within other caps.
No threshold-loss rate is extrapolated beyond that sampled prefix.

Analyze cumulative prefixes at 512, 2,048, 8,192, 32,768 and 131,072 blocks,
and at terminal collection. Analysis costs remain inside the diagnostic
run's total allowance. Report both nested policies, ideal all-record incidence
and actual SLP matching/eviction. Record LP cycle nullity, degree/occupancy,
component sizes, factor-base singleton losses, post-filter core, independent
lifted dependencies, exact square congruences, proper factors and unresolved
cofactors. Verify every reported dependency, including filtering discoveries.
Every terminal result reconstructs the original n. These are instrumented
feasibility diagnostics, not accepted performance timings.

## Investment decision and stopping

A production implementation is justified if at least one candidate policy
on both independent inputs of any declared band produces independently
verified proper divisors from the diagnostic relations, with at least one
input achieving this while the same-position SLP-only incidence has no
proper divisor. Require at least 16 additional independent LP-cancelled
constraints, no arithmetic/certainty failure, and no record/analysis capacity
censoring before that prefix. Require splitting plus certification plus
complete offline processing CPU to fit within 50% of the uninstrumented
SLP control's CPU allowance or observed completion CPU, whichever is larger.
This is an investment gate, deliberately distinct from a factoring speed win.
A diagnostic collection deadline after a qualifying complete prefix does not
invalidate that prefix; retain the terminal censoring as well.

If no band qualifies, finish with the specific measured bottleneck and a
scoped decision about these policies/workloads. A capacity/time-censored
larger case is inconclusive, never proof that DLP is useless at larger sizes.
No repeated expansion after this follow-up. If a band qualifies, freeze the
implementation/training/fresh-confirmation protocol before production work;
opt-in support can pass independently of default promotion. It must preserve
all R3 ownership, lifting, certainty, bounded storage, cancellation and charged
resume contracts. Complete-factor comparisons then include all setup,
splitting, certification, graph, filtering, provenance, matrix and extraction,
with >=3 seconds validated PyPy 3.11 warmup and >=9 samples, extending unstable
measurements under the revised roadmap promotion policy.

Hold the machine-wide `/private/tmp/factor-performance.lock` for all runs and
heavy checks. B3 currently has the next measurement window. Source review and
implementation preparation may continue without measurement overlap. Do not
change B3/A7 code or shared schemas without coordinating. Generated evidence
stays ignored; version protocols, required controls/corpora, runners and tests.
