# C10a-01 v4 experiment design

Status: design only, 2 November 2026. No new numerical calls are authorized by
running the design checker. This is a separately identified proposal; the
original calibration is closed. Keep all previous freezes and receipts.

Implementation: `4b9a616`; documentation handoff: `17e58c6`. Current committed
mainline control: `0076e256871722d52d849ab344152e0a7535c2ce`, retained as an
owned source snapshot with independent module/type identities. The control
has the current combined queue/sign defaults and historical standard
arithmetic. The challenger is `configuration_v4.py`; no automatic factory
selection is allowed. Mainline advanced only in documentation since the
previous comparison. Do not rebase, merge, push or promote defaults.

Use Astra xhigh throughout this task, PyPy implementing Python 3.11, native
integers, one active performance worker, and `/private/tmp/factor-performance.lock`.
No other sessions or agents are needed. Ordinary correctness/harness QA is
separate from the finite experimental allocation and remains serialized.

## Questions and order

1. Does the full combined arm perform feasible, useful collection on current
   representative 60/62/65/70 inputs under matched resources?
2. Which root construction pays for itself over the actual 64-polynomial
   family prefix? Does the new A population improve retained relation yield
   enough to justify its total cost? Does norm ranking change anything besides
   selection setup on these exposed inputs?
3. Can one frozen nominee complete balanced instances under the separately
   funded 60/62/65/70 grants? These are feasibility pilots, not adoption tests.

Run clean primary prefixes first, beginning with the balanced 60 pair. A valid
slow, zero-row or censored control does not veto a valid candidate or larger
size. An incorrect output, impossible configuration or harness failure stops
the affected stage for audit. Independent sizes may continue once the failure
is shown local; never automatically retry the failed cell. Diagnose before
spending diagnostics if the prefix failed before reaching the mechanism.

Do not require a prefix speedup, a solver event, or 65-digit completion before
an informative complete-factor comparison at another size. No broad quota,
face, queue, bound or parameter grids are reopened.

## Frozen inputs and literal arms

Reuse the eight exposed, certified, distinct-factor semiprimes. Do not generate
new exploratory moduli or regenerate failed cells. The committed corpus hash
and exact fixtures are in `experiment_v4.json`. The algorithm receives N, seed
7 and explicit configuration only; factors/certificates are validation data.
These generator-conditioned inputs are not uniform RSA samples and cannot
become untouched confirmation. Balance classes remain separate in reports.

| Size | Balanced N mod 8 | Uneven N mod 8 | Base count / A factors | Owned / RSS tier |
| --- | ---: | ---: | --- | --- |
| 60 | 1 | 3 | 5000 / 8 | 1 / 3 GiB |
| 62 | 5 | 7 | 6000 / 9 | 2 / 5 GiB |
| 65 | 3 | 1 | 8500 / 9 | 2 / 5 GiB |
| 70 | 7 | 5 | 10500 / 10 | 4 / 8 GiB |

C is the isolated committed mainline control. F uses parity2, parity2_norm,
windowed_pair, batch root tables and parity_two marking. Both use M=65536,
block=4096, 64 polynomials per A, pool=64, at most 100000 families, the same
LP domains, queues/plans, solver/pacing and declared capacities. The entire
literal mappings, including parent memory subtraction, are frozen in JSON.

Other labels are diagnostics, not additional complete-factor arms. Root
variants change only standard/cached/batch root updating at fixed k. A variants
change only windowed_pair/paired_flyer at fixed k. Ranking variants change only
parity2_norm/parity2_wide with automatic k. Forced standard/parity forms use
the challenger runtime and are not mislabeled as mainline control.

## New short calibration proposal: 3400 wall AND CPU seconds

| Phase | Calls | Fixed work / scope | Per-call cap + margin | Reserved seconds |
| --- | ---: | --- | --- | ---: |
| Primary | 16 | Two inputs x C/F at each size; 512 polynomials | 30/45/60/90 + 5 | 980 |
| Root construction | 6 | Three policies at balanced 62 and 70; 9 samples each | 45 + 5 | 300 |
| A population | 4 | Two policies at balanced 62 and 70; 512 polynomials | 45 + 5 | 200 |
| Ranking setup | 4 | Two policies at balanced 62 and 70; 9 setup samples | 45 + 5 | 200 |
| Forced multipliers | 4 | Declared forms/residues below | 45 + 5 | 200 |
| Observer | 4 | Two on/off pairs, 32 polynomials | 15 + 5 | 80 |
| Warmup/drift | At most 40 contexts | At least 3 validated wall and CPU seconds | 25 per context | 1000 |
| Analysis/validation | Bounded receipts and selection | No new factoring | Aggregate | 180 |
| Orchestration | Startup, locks, watchdog and accounting | No new cells | Aggregate | 260 |
| Total | 38 measured calls maximum | No new matrix campaign or generation | | 3400 |

All phases are maxima, not mandatory spend. Reuse warmed processes for
compatible cells; the warmup pool permits fresh processes after native RSS
censoring without repeating a measured call. Unspent reserves cannot fund
extra arms, larger caps, additional inputs or retries. Root/ranking sample
loops must finish all nine samples within their existing call caps; incomplete
samples remain recorded and cannot be called a successful timing comparison.

The 512-polynomial prefix visits eight A families, 16896 blocks and exactly
67109376 positions, including every one-position interval tail. It quadruples
the original two-family prefix to expose more A variation. Retained 128-
polynomial observations took about 4.8–7.3 seconds, implying a rough 19–29-
second baseline for this work; new setup and policy costs are unmeasured.
The overall 12–20-minute estimate is provisional, not a scaling law or deadline.
If a cap binds, retain it as incomplete work; do not compare its raw relation
count or elapsed time to a completed prefix as equal-work throughput.

Primary order is C/F on balanced 60, F/C on balanced 62, C/F on balanced 65,
F/C on balanced 70; reverse each order for the corresponding uneven input.
For diagnostic pairs, use opposite arm orders at 62 and 70. Root and ranking
sample orders are counterbalanced within a cell group: rotate the three root
orders across nine repetitions; alternate the two ranking orders and reverse
the starting order at the second modulus. Reserve the whole group before its
first sample; all individual sample charges roll into the named call receipts.
No order is chosen after seeing a timing. For interleaved root/setup cells,
the cell wall cap applies to accumulated active elapsed intervals, including
its validation and receipts; another policy's interval is not charged twice.
The group also has a continuous wall deadline equal to the sum of its cell
reservations. CPU and shared setup are recorded once in the group ledger;
per-policy receipts identify that shared charge separately. Pre/post warmup
is separately reserved, and a pending group retains its full charge if lost.

### Root diagnostic

At balanced 62 fix k=13 (parity eligible); at balanced 70 fix k=1 (standard
fallback). These values come from pre-timing structural QA. Build the actual
base and first eight A assignments once per compatible group, accounting for
shared setup separately. Each of nine samples constructs each family anew
and consumes all 64 Gray steps. Include CRT/inverse/table construction and
all transitions; do not time only the cached subtraction loop. Base and
assignment setup are not falsely credited to any one root method.

Record inclusive call cost and the sum of family construction/traversal
intervals. Stream final-state checks after each interval, retaining at most
one family and its current/previous roots. Record identical A/B and final
root fingerprints, step counts and direct-integer identities. The independent
preflight verifies all transitions against CRT/brute-root oracles; timed
fingerprints alone do not constitute that proof. No sieving occurs in this
microdiagnostic; it cannot establish an end-to-end speedup.

### A and ranking diagnostics

The A comparison uses the same fixed k=13/1 values and all other F policies;
it compares two deliberately different populations over 67109376 positions.
The two selected inputs exposed large old paired-core misses before timing.
Report this exploratory conditioning. Include assignment/setup cost, actual
A errors, sampled norms, candidates, refined rejections, admitted/retained
rows, partials, graph/provenance growth and RSS. A reduction in norm size is
not a useful-yield or complete-factor result.

The ranking comparison measures nine fresh `SIQSJob._setup()` operations per
policy, including the selected base and A setup. Both policies currently choose
the same k on all eight exposed inputs. Record exact choices, score receipts,
base/assignment identities, work and both clocks. Do not spend four extra long
prefixes re-sieving identical populations merely to rediscover setup overhead.
If a choice changes, do not silently convert the setup experiment into a long
run; primary F remains measured and the ambiguity is retained for nomination.

### Forced coverage and observer checks

Use four full sieve prefixes: balanced 60, k=5, standard and parity forms,
768 polynomials each; uneven 62, k=2 (D mod 8=6), 640 polynomials; balanced 70,
k=31 (D mod 8=1), 512 polynomials. These target approximately 25–30 seconds
of retained reference work within 45-second caps. They exercise polynomial,
root, scoring/refinement and relation paths, not just multiplier selection.
The fixed-k form pair keeps the other policies equal. Its actual-A targets
and physical trial populations differ deliberately; exact coordinate
reparameterizations are separately verified by existing arithmetic oracles.
Missing D classes and exceptional k cases remain covered by correctness QA.

Observer pairs use F on balanced 60 and 70: parity and standard fallback.
Each on/off call covers 32 polynomials / 4194336 positions / 1056 blocks.
Keep the existing 8192-record ceiling; never enable it for the 16896-block
primary prefix. Compare semantic state AND charged work on identical work.
A record failure is `recorder_exhaustion`, not a native memory outcome.
Two pairs are an overhead/passivity screen, not a precision timing claim.
All primary and pilot calls remain uninstrumented regardless of this result.

## Metrics, passivity and clocks

Start complete-cost timing before multiplier selection and stop after exact
result reconstruction/classification and the minimal bounded receipt. Include
base/A setup, collection, failed splitting, preparation, filtering, solving
and extraction. Record cold startup/import, warmup and analysis separately.
Population scans, full semantic digests and offline norm summaries are charged
to analysis, outside the uninstrumented timer. They may not retain the job
through a subsequent call or change its cursors/charges.

Record wall AND process CPU, total and stage charged work, actual positions,
A families and polynomials, scored/forced k, D mod 8, active representation,
retained atoms/rows/partials, candidate/refinement counters and explicit
termination. Record process RSS high water with scope including warmup and
prior compatible calls. Peak reservations and summed peak counters are distinct.
Report norm sampling as sampling, not the full distribution.

For corresponding root or same-k ranking trials, compare semantic digests
without budget totals or policy-specific family IDs; compare accounting
separately. Observer on/off must also preserve work. The old state digest
includes work and cannot be used unchanged for root-policy equivalence.
Offline A/CRT descriptions cache each of at most 64 family records once;
never rescan a large paired completion for every polynomial in the family.

Every return reconstructs N, including unresolved cofactors. Proper divisors
are checked independently against exposed factors only after the algorithm
returns. Preserve proven/probable/unclassified labels. Early complete factoring
is a valid outcome, but no longer an equal-position prefix observation.

## Warmup, drift, resources and failures

A warm context has 25 reserved seconds: prewarm/sentinels at most 17 including
2 seconds of margin, post-sentinels at most 8 including 2 seconds of margin.
Prewarm must establish at least three validated wall AND CPU seconds, followed
by three fixed eight-polynomial sentinels; post uses the same three sentinels.
Root/setup-only diagnostics use their own semantic work for warmup/sentinels.
At most 32 fresh semantic invocations fit a context; stop if either time or
work guard refuses. Reuse compatible warm contexts, not stale warmup receipts
from another process. Missing post-drift after native resource refusal is
censoring, not successful readiness. Sentinel ratios/spread are descriptive;
flag a >20% median shift or >20% relative range as drift. Preserve the result;
do not automatically repeat a drifting call. These are readiness tolerances,
not a minimum acceptable performance gain.

Keep the smallest justified matched tiers: 1/3 GiB at 60, 2/5 GiB at 62 and
65, 4/8 GiB at 70. Parent reserves 32 MiB owned and 256 MiB RSS inside each
tier; literal native job grants already subtract the bounded harness/assignment
reserve. No grant may silently alter score, queue, plan, backend or pacing.
The root cache adds about 4–8.4 MiB; the retained conditional central ownership
estimates including parent are about 0.769/0.924/1.346/1.728 GiB. They assume
5F atoms, 4F graph edges, 1.2F rows and bounded exponent/provenance populations;
they are not RAM requirements. Row-stress cases still refuse at 60/65 under
their smaller grants. Existing dense-fill/matrix bounds remain, and new native
RSS is unmeasured. No new matrix fixture is needed merely for this root cache.

Verify actual RAM/headroom before every process admission: total tier RSS plus
4 GiB host reserve. Last check was 24 GiB physical RAM; never assume it remains
available. Cooperative checks may overshoot during an atomic operation.
After native RSS/memory censoring, terminate that worker before the next cell;
continue the counterpart in a fresh warmed process within pre-reserved slots.
No emergency full checkpoint serialization: its existing ownership limitation
remains. Small charged resume is correctness-tested, not a license to resume
a censored experiment with larger limits.

The new allocation permits 64 MiB of output INCLUDING its evidence archive:
40 MiB raw, 24 MiB archive. It also must fit the shared 512 MiB study ceiling;
about 383 MB was owned at the implementation handoff. Recount before launch.
Receipts <=512 KiB; worker stdout <=64 KiB, streamed with a bound before parsing;
parent summaries stay bounded. Do not capture unbounded stdout and check its
size afterward. Limit polynomial metadata to 4096 and family metadata to 64.
Recorder failures stop the affected stage for repair under a new identity.

Reserve each call before launch. Retain pending/interrupted reservations at
full declared charge when no trustworthy receipt exists. A lost or failed
measurement is not retried. Stop an uninformative diagnostic group rather
than use its remaining allowance for a larger grid. Unexpected exceptions,
invalid arithmetic or policy/source drift require audit. Native time, CPU,
work, family and memory exhaustion remain legitimate outcomes. A validated
candidate can continue even if its control is censored.

## Work allowances and source freeze

`plan_v4.py` derives every integer limit in `experiment_v4.json` and checks
all phase sums. Use the retained upper sustained rate of 8 billion units/CPU
second with factor-four operational headroom, plus 1.1 trillion separately
bounded setup units. Add the full declared schedule's root-cache, A-pair and
window charges and the existing analytic preparation/filter/solve/extraction
bounds. The last bounds are intentionally very loose and allow dense fill;
these huge work guards are not runtime forecasts or digit-only memory models.
They prevent an early collection rate being mistaken for a mature solver rate.
Actual production charges enter Budget once; allowances do not add synthetic
work to reported consumption. No grant is increased within a measured call.

Root/setup diagnostics have no matrix phase: their work limits are 4*rate*time
plus one/nine setup bounds and 100 billion units covering all declared family
traversals. Each warm context reserves 32 times the largest applicable
20-second full-call guard, covering its finite fresh invocations. Its time
ceiling remains 25 seconds. Analysis has 1 trillion units per stage; orchestration
has a time allowance and no algorithmic work. Freeze CPU/rate clocks and final
units exactly; JSON integers are not floating-point estimates.

Record new sustained rates for setup, collection, root diagnostics and any
naturally reached postprocessing. Use completed intervals of at least 0.25 CPU
seconds or aggregated compatible samples; an instantaneous atomic spike is not
a sustained rate. Before pilot admission use max(8 billion, observed sustained
maximum), retain fourfold headroom and analytic charges. If the existing
numerical pilot guards are inadequate, create a separately identified
prospective release before any pilot; do not overwrite this design or failed
outcomes. A measured work-limit termination is retained even if time remains.

## Nominee selection before complete pilots

The initial full arm is tested as implemented. Nomination may remove components;
it is not default adoption. Apply these rules once, using all valid outcomes,
and commit the chosen configuration/source hashes before any complete pilot.
There is one policy bundle across all four sizes, not a hidden per-input table.

- Root policy: require complete nine-sample groups and equivalent physical
  trials. If one policy has the lowest inclusive family cost at both anchor
  inputs and its advantage exceeds within-input uncertainty, nominate it.
  Use paired 10000-resample 99% intervals for median CPU differences across
  the nine counterbalanced repetitions (analysis seed 20261102104), only as
  repeatability diagnostics. Require the upper endpoint below zero against
  each alternative at both anchors; flag contrary wall-cost/drift evidence.
  If there is no common resolved winner, retain batch as the original hypothesis and
  record the uncertainty; do not launch more cells. Two moduli do not provide
  population-level speed evidence.
- A policy: retain windowed_pair unless the fixed-k comparisons consistently
  show more CPU with no retained-row/partial improvement at BOTH anchors,
  outside the contextual noise band, or it has a correctness/capacity defect.
  Define e as the larger relative pre/post sentinel median change or relative
  MAD across the six sentinel values, taking the maximum across both policy
  contexts at that anchor; require CPU_window/CPU_paired > 1+2e
  at each anchor, with both retained rows and partials no greater. This is a
  conservative exploratory switch rule, not a formal single-prefix speed test.
  In that case select paired_flyer. Mixed or censored evidence is unresolved,
  not a reason to exclude a feasible whole-job test. Rows are a proxy only.
- Ranking: if norm and wide choose the same k on ALL eight primary inputs,
  read the wide winner from each full-arm density-score receipt and use
  parity2_wide for the feasibility nominee: norm preview changes no
  subsequent trials here and adds discarded setup work. Preserve and report
  the measured setup cost. If any k differs, retain parity2_norm as the
  predeclared hypothesis; do not tune a new score on observed completions.

The selection menu is bounded to the implemented root choices, the two A
policies and the two wide/norm policies. Evaluating their configuration validity
is not executing their Cartesian product. Correctness or infeasibility can
reject a component; performance uncertainty cannot be disguised as rejection.
Validate the chosen combination and literal matched grants before release.
No claim that component gains add, that the chosen bundle is globally optimal,
or that prefix selection establishes a complete-factor speedup.

## Later complete-factor stages: keep all four sizes

| Stage | Balanced inputs / calls | Per-call wall AND CPU cap | Measured reservations | Warmup / analysis / orchestration | Total |
| --- | --- | --- | ---: | --- | ---: |
| A | 60 and 62; 4 | 1800 / 2400 | 8420 | 600 / 120 / 60 | 9200 |
| B | 65; 2 | 3600 | 7210 | 300 / 60 / 30 | 7600 |
| C | 70; 2 | 7200 | 14410 | 300 / 60 / 30 | 14800 |

Each stage has its own resource/arithmetic/readiness decision. B need not
complete before C. All pilot calls are clean, fresh-process, setup-to-validated-
result runs with the same matched tiers. Preserve a control failure and still
run a valid nominee. No larger-tier retry or second pair. Rewarm after process
creation; all parent/worker CPU and wall reservations must fit before launch.

Balanced 70 currently chooses k=1 and standard fallback. The new A/root policies
still differ from mainline, so this is informative about the whole scored arm.
The forced k=31 diagnostic tests parity at 70 separately; it does not silently
replace the scored pilot or add a third complete-factor call.

The 60 grant is informed by retained current-policy 426–796-second observations.
The old 65 compact 1800-second timeout preceded collector improvements and is
not a scaling curve. No calibrated expected complete runtime exists at 62/65/70.
The 70 grant funds an attempt, predicting neither completion nor a speedup.

The revised whole envelope is 3600 original + 3400 new + 31600 pilots = 38600
wall AND CPU seconds (10 h 43 min 20 s). With original actual reserved charges
3124, remaining prospective grants total 35000 and the corresponding charged
ceiling is 38124 seconds. No part of this proposal has run. Confirmation and
new storage/solver/native-reference campaigns are excluded and unfunded.

## Analysis, stopping and later confirmation

Report every input separately: finite-budget completion, actual CPU/wall cost,
termination, work, memory and retained population. For unfinished calls, capped
service cost assigns the full cap even to early resource refusal. Report actual
cost alongside it. If both arms finish, a descriptive complete-cost ratio is
valid for that input. If only one finishes, report completion difference and
capped cost; never claim successful-factor speed against an unfinished control.
One balanced input per size is feasibility evidence only. No pooled ratio may
hide a band or balance-class regression.

Early failure, recorder exhaustion, drift and source changes have distinct
labels. Missing post-drift means uncertainty, not automatic deletion or success.
Do not select only successful pairs, repeat favorable seeds as new moduli, or
stop on the first attractive ratio. Freeze all analysis rules before execution.

A promising nominee can justify a separately budgeted fresh comparison after
selection is committed. At least nine independent moduli per claimed band is
a floor, not automatic statistical sufficiency. Choose the sample count,
residue/balance strata, caps, primary capped-CPU-cost endpoint, paired input-level
uncertainty and multiplicity handling before revealing inputs. Resample whole
modulus pairs, not repeated seeds. No performance-percentage threshold is
universal; evidence must resolve benefit and unacceptable completion harm.

If a claim requires a one-sided 95% bound below 5 percentage points of completion
harm, nine zero-harm pairs cannot establish it. Even zero harmful discordances
need at least 59 independent pairs for the simple exact binomial bound
(0.95^59 < 0.05). That cost is not automatically funded. Do not turn a small
pilot into an adoption claim to evade this uncertainty. A broader fastest-Python
claim additionally needs the separate current-competitor protocol and matched
runtime/backends/resources; this workstream alone cannot establish it.

## Executable design check and execution prerequisites

From the task worktree, run:

    /opt/homebrew/bin/pypy3 -m v2.benchmarks.qs.c10a01.plan_v4 check

This verifies inputs, configurations, source/runtime hashes and exact sums;
it has no run action and creates no factoring job. The source snapshot loader
can be imported/tested from committed files without Git history or old worktrees.
The old pilot runner remains bound to its old source and must keep rejecting
this runtime. Do not repoint it or edit a historical JSON freeze.

Before execution, implement bounded adapters for the explicit cell kinds and
ledger/worker contract above, then pin them in a new execution release. Required
adapter tests cover position/tail targets, early factors, semantic/work digests,
shared setup accounting, nine-sample interleaving, recorder exhaustion, native
RSS reset, child-output bounds, warmup clocks, interrupted reservations and
nominee/source substitution. Validate these with stubs/small correctness cases,
not a hidden exploratory timing run. After calibration, emit the chosen nominee
and pilot releases with exact hashes and updated rate admission. These are
explicit remaining harness tasks, not missing arithmetic implementations.
