# C9.1 endpoint admission and finite retention: frozen study

Status: training protocol, 10 October 2026. Parent: integrated `b3b3cfb`.
This closes only a named bounded C9 tranche if its gates pass. C10 collector
calibration, E1 combined acceptance, F6 higher arity and matrix/provenance
redesign remain separate. Nothing authorizes merging or pushing.

## Question and source rationale

Compare large-prime **endpoint alignment** with increased retention, starting
with unchanged C1 collectors. In the 50-digit control, SLP admits primes through
900,000,000 while DLP endpoints stop at 3,000,000. A prime above the latter
cannot link a DLP edge. Keeping it competes with repeatable smaller endpoints
in the finite FIFO forest. YAFU/msieve and Yamaquasi's common endpoint ceilings
motivate aligning the two, while their growable partial stores motivate a
capacity control. Neither proves a win. Rejecting larger SLP primes deliberately
loses some supported matches; retain the existing unrestricted option.

C1's proven endpoint contract and all-component forest already exist. Its
fresh comparison found no failed split and its profile was sieve dominated.
Do not implement a new splitter or weaken certification without new evidence.
The historical 40-digit 31.73% result, 50-digit 0/18 per arm at 120 seconds,
21,504/28,672 median evictions and separate 141.063-second half-base resume
remain their original observations. They do not isolate a 50-digit cause.

## Finite envelope and advance rules (before training)

- Total active study allowance: **14,400 wall seconds and 14,400 CPU seconds**
  across training, diagnostics, warmup, cold calls, confirmation and validation.
  Waiting for the shared lock and ordinary source research are excluded. Record
  every launched run, including failed/censored runs and harness repairs.
- Initial training: the two already inspected C1 50-digit training inputs,
  seeds 7 and 29 respectively; five arms once each, maximum 240 seconds each.
  Observe cumulative **60, 120, 240 seconds**, or the earlier proper factor or
  finite resource refusal. No termination just because an early matrix is empty.
- Each call: work `10**13`, wall/CPU 240 seconds, owned memory **256 MiB**,
  checkpoint configuration 1 MiB (no checkpoint in attribution timings),
  at most 32,768 atoms and 8,192 rows. The 16,384-partial arm changes only that
  capacity. A capacity refusal is an outcome; do not raise memory to rescue it.
  RSS is reported separately, with a campaign abort at 1 GiB RSS. Scratch and
  simultaneous graph rebuilding use existing reservations. Capture disk cap
  256 MiB; checkpoint diagnostics may use a separately declared 8 MiB checkpoint
  and 512 MiB memory, never substituted into comparisons.
- Arms: C1 SLP; SLP with residual limit `100*base_bound`; current C1 DLP;
  DLP with SLP limit equal to its DLP endpoint; DLP with 16,384 unowned edges.
  DLP candidate/product bounds remain unchanged in endpoint/capacity arms.
  The SLP low-endpoint arm necessarily narrows its score domain and is separately
  labelled as SLP calibration, not residual-only attribution.
- Select the SLP control and DLP challenger by number of complete outcomes,
  then sum of complete-call/capped cost; stable listed order breaks ties.
  Graph counts cannot select a winner over complete outcomes.
- After individual evaluation, at most **two** additional configurations across
  the same two inputs (four 240-second calls) are allowed under a new explicit
  source/configuration freeze: endpoint alignment plus capacity; or a selected
  policy with the already historical half-base collector. The latter requires
  recording the C10 interface boundary; it is a bundle, not isolated C9 credit.
  Advance only if an individual arm completes or has a surviving filtered
  matrix at the long prefix. If every arm has zero survivors at 240 seconds,
  stop training and retain an inconclusive feasibility verdict. No 60–100-digit
  campaign is authorized by this protocol.
- At most 600 active seconds for additional attribution/profile/census work,
  360 for charged checkpoint/resume, 120 for corpus generation, 600 for heavy
  QA. The aggregate cap dominates every suballowance. No unbounded search.

## Measurement and confirmation freeze

Freeze runtime hashes, runner, actual configurations, selection and corpus
generation rules **before generating or revealing** new inputs. Use three
new independently generated balanced 50-digit inputs and three balanced
40-digit regression inputs, each with seeds 7/29/47: **nine matched samples
per class**, three independent inputs, not nine independent inputs. A new
30-digit balanced input supplies nine cheap regression samples. Generate
Pocklington-certified primes with a fixed new seed; disclose the large p−1
factor conditioning. Compare all new n against every versioned corpus.
Expected factors/certificates enter output validation only.

Confirmation compares selected SLP, current accepted C1 DLP (half-base at
30/40, full-base at 50), and the selected challenger. All arms receive identical
total allowances: 240 seconds at 50, 30 at 40, 5 at 30; 256 MiB memory,
`10**13` work, identical row/atom/checkpoint caps. Changes at smaller bands
must be explicitly frozen; retaining C1 unchanged there is a valid narrow scope.
At least three seconds of validated PyPy 3.11 warmup precedes measurements
for each arm/band. Counterbalance arm order deterministically and rotate seeds.
Instrumented prefixes/profiles and one cold process per arm/band are separate.

Use paired capped-total-cost reduction and a 10,000-resample hierarchical
bootstrap over independent inputs and within-input paired samples. Report
individual-input and seed effects, completion and all censored costs. Extend
from nine to eighteen matched class samples once if within-input/seed paired
benefits reverse sign, or repeated fixed-seed first/last medians drift by >10%,
or relative IQR exceeds 10% for a repeated cell; a single observation is not
a stability estimate. Fresh input heterogeneity alone is not temporal drift.
If the aggregate cap cannot accommodate the complete extension, keep all
observations and mark timing inconclusive. No tuning on confirmation.

Any sustained complete-run gain may qualify; there is no 10% benefit floor.
Require a positive 95% interval, no correctness failures, unchanged finite
contracts, acceptable complexity and at most five percentage points completion
regression in each declared class. Narrow scope if necessary. Equal capped
failures earn no promotion. More duration/memory by itself earns no speedup.

## Attribution, correctness and ownership

Report candidate/refinement and residual counters, splitting/certification,
admitted atoms, endpoints/degree histogram, retained graph components and
cycle/core development, eviction, singleton survival, last solver dependencies,
verified extraction/proper factors, CPU/wall/work, stage times, owned memory,
RSS and checkpoint costs. A retained forest alone has cycle rank zero; compute
graph diagnostics over retained original atoms including closing edges and
distinguish cumulative emitted cycles from retained incidence rank. Core and
matrix counts are diagnostics, never success criteria. Preserve raw captures.

Use C1's retained complete census and independently sampled rejected positions
with explicit coverage cutoff. Extend sampling only if the policy diagnosis
needs an unobserved domain. Do not infer rejected residual factorization from
surviving atoms. Source pinning, licenses and rationale precede experiments.

No polynomial, multiplier, factor-base generator, score or recovery edits in
this freeze. C9 owns new benchmark files and residual-policy configurations;
potential shared runtime files are `sieve_collector.py`, `large_primes.py`,
`large_prime_collector.py`, `large_prime_checkpoint.py`, `siqs.py`. Any runtime
change requires a new identity/compatibility review before accepted timing.
README/roadmap/changelog additions are isolated commits for later integration.
C3 owns portfolio/handoff; C10 owns collector strategy. Inspect active chats
and worktrees; record any new collector-owner coordination before joint arms.

All timings and heavy checks hold `/private/tmp/factor-performance.lock`,
publish `/private/tmp/factor-performance-owner.json`, and refuse competing
benchmark/test processes. All outcomes independently reconstruct; every divisor
is proper and terminal factors agree with independently verified certificates.
Preserve independent small-graph kernel, arithmetic, corruption, ownership,
publication, refusal, cancellation and repeated-restore tests. Run
`make -C v2 test`, `make -C v2 lint`, then committed-files-only suites and all
benchmark imports/loaders before publishing reviewable commits.
