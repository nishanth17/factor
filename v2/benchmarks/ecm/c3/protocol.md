# C3 finite allocation/handoff protocol

Frozen before pilot timing and before confirmation generation. The JSON
control pins the mainline source snapshot, prepared implementation, runner,
relation configurations, historical training corpora and this protocol.
Changes after freezing require a separately named revision and renewed freeze;
never replace an earlier protocol or use confirmation to select a new winner.

## Source and objective

Control is integrated mainline b3b3cfbea6105f08db0a8484ec09c8260d310e28.
Run native integers on PyPy implementing Python 3.11. Use accepted ECM,
p−1 and QS machinery. No kernel, chain search, collector, splitter or graph
change belongs to this experiment. A6 p−1 and B3 chain routing remain intact.

Rank by complete factorizations, then summed capped complete-call wall cost,
then the declared arm order. Select one policy per observable size band
(<=35 digits uses the 30-digit bundle, otherwise the 40-digit bundle).
Factor-size/structure labels are available only for generation, validation
and reporting. Input fixtures and certificates never enter the dispatcher.
Do not claim optimality or general automatic dispatch from these two bands.

## Fixed candidate set

Arm order: control, no_ecm, quick8, short4, tiered, wide1. Bounds/counts are
pinned in the JSON and runner. Control retains 32 curves at 2000/147396.
No-ECM retains preprocessing/rho/p−1. Quick8 tests Yamaquasi's *mechanism*
with eight 200/7700 curves. Short4 tests four existing-bound curves.
Tiered uses eight existing-bound curves then two independent 11000/1873422
curves. Wide1 tests one 50000/12746592 curve. The latter two use explicit
campaign mode. A changed bound always starts a new curve; it never claims
same-curve continuation. These are six bundles, not a Cartesian grid.

Pretests have absolute cumulative ceilings of 2M work and 0.5 wall/CPU
seconds, including earlier preprocessing and children. Classification and
exact power checks remain mandatory under the total allowance even with a
zero pretest ceiling; optional-stage admission then observes all their cost.
Explicit policies prepare schedules lazily and charge snapshot verification
bytes on resume, in addition to elapsed decoding/rebuilding time. Fallback reservations
are 500M work/1 wall/CPU second at 30 digits and 15B work/10 wall/CPU seconds
at 40 digits. They are conservative admission floors, not completion promises.
Campaign arms have the same reserves but no extra pretest ceiling. All arms
receive identical total allowances: 10^13 work, 5/30 wall and CPU seconds,
288 MiB owned workspace, one process/core. The relation engine is the frozen
C1 30-digit calibrated SLP or the accepted 40-digit half-base DLP bundle.
This larger service envelope is a separately labelled policy regime; the
control receives it too. Default 2M-work/2-second diagnostics remain censored
resource probes, not accepted speed evidence or a raised production default.
Separate default-grant and explicit enlarged-total resume diagnostics have a
600-second total cap; report their single-run costs only as functional data.

## Training, pilots and stop criteria

Previously inspected C1 balanced 30/40-digit fixtures and B1 uneven-10,
power, smooth p−1 and close controls are training only. Use the first input
per kind and size, seeds 7/29. A balanced-only one-sample feasibility pilot
has at most 900 active wall/CPU seconds. It may reject an infeasible study,
but cannot change this candidate set or numeric settings. Training has
three matched samples after >=3 validated warmup seconds per arm/input and
at most 2100 active wall/CPU seconds. The stopped revision-4 training
window consumed 258.062 wall seconds; charge 300 seconds including an active
call/CPU allowance. Combined training remains bounded by the original 2400
seconds. Its completed partial rows are diagnostic only. Training estimates select settings;
they are not accepted performance comparisons.

Confirmation has at most 7200 active wall/CPU seconds, inclusive of warmup.
Check the envelope before each call; a last call can overshoot by at most its
finite 30-second cap plus cooperative atomic/output overhead. Save stopped
runs and declare incomplete acceptance if the study allowance is reached.
An engine failure, invalid result or competing heavy process fails the study.
No concurrent heavy QA, timing, compression or profile is permitted under
/private/tmp/factor-performance.lock. Publish C3 owner PID/cwd and inspect
the process inventory before the window; release it in finally.

## Selection before fresh generation

Commit the selected policy and training capture digest before generating
confirmation with seed 2026101103. Generator has 200000 random-draw and
1000 decimal-band attempt caps. Each 30/40-digit band gets two balanced
inputs and one each with approximately 8/12/14-digit smaller factors at 30 digits and
8/12/16-digit smaller factors at 40 digits.
Add recursion (three factors), a square, a prime, close small factors and a
smooth p−1 control. Independent recursive Pocklington certificates and exact
trial-division leaves validate expected factors. Large-q p=kq+1 generation
biases p−1; disclose this and the small number of independent inputs.
Rejected generation attempts are not factoring observations.

Compare control only with the selected policy, seeds 47/71/101, at least
three seconds of validated warmup and nine matched samples. Rotate/reverse
arm order deterministically. Per input/arm/seed relative IQR >15% extends both
arms to 18 samples after five seconds warmup, then at most 27 after eight.
Stable inconclusive outcomes do not justify more samples. If selected=control,
report retain without running duplicate indistinguishable arms.

## Metrics, uncertainty and decision

Capture complete/proper-factor yield, capped cost and complete-call time,
CPU/wall/work, finished curve counts/bounds, partial-curve handoff frontier,
fallback entry and execution, recursive events, recovery counters when
published, checkpoint bytes, configured owned reserves and process RSS
high-water. A high-water RSS includes earlier arms and JIT; it is not a
per-call heap measure or an OS-enforced cap. Setup, failed searches, output
packing and independent output validation are inside the measured call;
configuration decoding and certificate verification are outside, identically
for both arms. User-seed matching does not force identical SIQS assignment
seeds after different counts of ECM draws; report this source of variation.

Use paired medians per input/seed, retaining censored calls at their cap for
the cost objective. Bootstrap input subjects with all seeds nested, 4000
resamples and seed 193003; report a conditional paired 95% interval and
per-class absolute completion/time. Repeated timings do not add independent
inputs. Fewer than two independent inputs in a class support regression
controls only. Zero correctness failures, a positive improvement interval,
stable sampling and <=5 percentage-point completion loss in every declared
class are required for scoped promotion. Retain current defaults when
inconclusive; do not retune on confirmation. Any supported opt-in scope must
state its workload and service allowance. Broader G1/E1 remains open.

Recovery/resume and legacy compatibility are acceptance tests independent of
this speed study. Cold process startup and instrumented profiles, if run, are
separate diagnostics. Upper-band 50–99-digit general dispatch, native upstream
execution, same-curve extension and new stage-two geometry remain deferred;
reuse existing negative/censored evidence instead of asserting a universal
30/120-second feasibility ceiling.
