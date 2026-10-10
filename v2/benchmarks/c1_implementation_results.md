# C1 bounded DLP implementation and complete-factor evaluation

Status: production correctness checks passed; frozen training is running.
Fresh confirmation and promotion are pending. This report will retain failed
and unresolved outcomes along with complete factorizations.

## Why implementation proceeded

The [research review](c1_research.md) verifies all-component cycle mathematics
against primary literature and distinguishes pinned implementation policy from
mathematical guarantees. YAFU/msieve/Yamaquasi provide residual-bound and
smaller-base hypotheses, with native cutoffs explicitly excluded as PyPy
predictions. FLINT and JavaMath provide further certainty/ownership comparisons;
no external implementation code was adapted. Licenses and retrieved source
identities are recorded in that review.

The initial512-block screen could not decide DLP economics. Its final defer
was withdrawn. The longer bounded40/50/60-digit study retained rejected
residuals and sampled positions as well as survivors. Both40- and both50-digit
inputs produced independently verified DLP complete factors at prefixes where
SLP incidence had no dependency. Entire-census splitting/certification cost
was under1s per40-digit input and under3.3s per50-digit input, instrumented.
Those figures exclude collection and are not accepted performance evidence.
At60digits, both runs remained unresolved: many degree-one vertices survived,
and every LP-cancelled row was singleton-filtered. Raising a solver budget
would not repair those matrices.

The original literal follow-up gate charged both competing policy analyses,
SLP analysis and exhaustive extraction after a factor. It failed and remains
preserved. The explicitly post-observation attribution repair used retained
records only, stayed inside the original4,500-second envelope, and charged
all processing for one candidate through its first complete factor. That
supports an investment go at40/50, not a retrospective preregistered pass or
a factoring-speed claim. Full data and limitations are in
[the feasibility report](c1_followup_results.md).

## Implementation and acceptance

Production source `ebd049e` adds explicit `DoubleLargeSieveConfig`, separately
bounded proven endpoint pairs, bounded rho/square splitting, and a finite
all-component forest. The scalar residual contract remains unit-or-proven-prime.
Every combined row retains original atoms, exact factor-base exponents and
square corrections. Repeated primes, self-loops and parallel/disconnected
cycles are supported. A256-atom cycle limit is an explicit yield loss.

Only unowned forest edges can be evicted. Planning is read-only, followed by
verification and a final resource check before publication. Owned paths survive
rebuilding. SIQS checkpoint4 reconstructs the forest and verifies every row,
closing edge, atom reference and polynomial ownership, preserving original
mixed row order. Prior work/wall/CPU, replay and started split attempts stay
charged. Ordinary SLP uses checkpoint3. The same explicit SIQS configuration
works through portfolio resume; no shared portfolio/stage schema changed.

Independent tests enumerate the entire boolean incidence kernel of small
multigraphs and compare its span with emitted cycles, including eviction.
Exact N=91 fixtures expose proper7/13 splits in a disconnected self-loop and
a disconnected parallel-edge cycle whose known square shares a large prime.
Further checks cover all four division backends, conservative and deliberately
narrowed candidate domains, corrupted corrections/stores, cap refusals,
cancellation, native/GMP arithmetic, and charged nested portfolio resume.

The committed-only `b5066ad` archive passed525 tests, full `make -C v2 lint`,
107 benchmark imports and required C1 loaders. The final orphan-polynomial
bound in `ebd049e` additionally passed the15 DLP tests. Full final-source
committed-only validation remains scheduled after confirmation.

## Frozen evaluation

The [implementation protocol](c1_implementation_protocol.md) fixes five arms:
committed calibrated SLP, graph-only SLP reuse, full DLP, narrowed DLP candidate
allowance, and DLP with half the numeric base bound. The [control notes](c1_implementation_controls.md)
explain the separate owned SLP source package, finite graph reservations and
charging. Original40-digit input IDs collide across the two source corpora;
the `run-40-0-*` and `run-40-1-*` filenames and independently certified factors
identify the two distinct cases.

Training selects by completion, then total capped cost. Fresh settings must
freeze before generating new30/40/50-digit inputs. Confirmation includes at
least3s validated PyPy3.11 warmup and9 paired samples, extending to18/27 on the
prespecified instability checks, within10,800 wall/CPU seconds. All setup,
collection, splitting, certification, graph/eviction, provenance, filtering,
matrix, extraction, classification and final output validation are timed.
Cold processes and profiles are separate. RSS is the process high-water mark,
including prior arms in the same process; owned reservations are per-job.

## Cross-owner comparison boundary

A7's integrated reconciliation (`018236f`, read from mainline without changing
this branch) prepares calibrated30-digit SSS comparison arms for E1. It
explicitly identifies the absence of a calibrated40-digit SSS arm and assigns
that training/confirmation to E1 or an explicit challenger-class deferral.
C1's frozen five-arm study therefore does not invent larger SSS settings or
claim a DLP-versus-SSS crossover. P5.4's broader SSS comparison remains an E1
integration prerequisite. B3 released the performance window at07:12:46UTC;
C1's committed-only QA and subsequent study each hold the machine-wide lock.
No B3/A7 algorithm or shared portfolio/stage-job file was edited.
