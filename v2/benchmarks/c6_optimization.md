# C6 executor optimization follow-up

The user reopened C6 on 9 October 2026 after the first bounded comparison.
The original captures, source hashes and negative verdict remain historical
controls. The new task is to optimize the executor and recovery policy before
repeating the chain comparison. No production kernel, default, stage job or
checkpoint is changed. The no-universal-percentage-floor policy applies.

## Why the first executor lost

Its full B1=2,000 Lucas stage has 4,301 point operations and, on the ordinary
unit path, 9,511 GCD calls excluding curve setup. The pinned ladder has 5,755
point operations and one final GCD. These are source-derived counts, not a
profile assigning time to individual operations. The earlier study establishes
that the conservative implementation loses; it does not isolate arithmetic
quality from that safety-policy cost.

The pinned GMP-ECM [top-level README, section 8](https://github.com/sethtroisi/gmp-ecm/blob/8ea5e214fdf2f0ddf9415141b8dc039ed6f5874e/README)
describes precomputed near-optimal codes as an optional `-param 0` stage-one
path replacing PRAC for primes >=11. It is not a promise that every GMP-ECM
configuration always uses that path. Its native dispatch and arithmetic
selection also differ from this Python/GMP-wrapper implementation. The
common-executor experiment separates chain quality from these costs.

## Exact factor-coverage certificate

Write a point as `(X,Z)`. All expressions below are modulo the odd modulus n.
For the pinned doubling formula,

`Z_double = 4 X Z ((X-Z)^2 + 4 a24 X Z)`.

Thus every prime divisor of either input coordinate that divides n also
divides output Z. In differential addition with known difference `(Xd,Zd)`,

`X_sum = Zd * u^2` and `Z_sum = Xd * v^2`.

Each difference coordinate's factors therefore propagate into the indicated
output. These polynomial divisibility identities do not assume the unknown
prime factors, curve order, squarefreeness or successful earlier arithmetic.

The planner constructs a DAG with one vertex per versioned coordinate,
including overwritten registers. Doubling adds edges from both input
coordinates to output Z. Addition adds edges from difference Z to output X
and difference X to output Z. Coordinates with no successor become guards,
except for the final output pair, which the caller promises to check. No
assumption is made that arbitrary addition operands propagate their factors.

An independent verifier propagates ancestor bitsets forward through the
record and requires the guards plus the final output coordinates to cover
**every** coordinate vertex. The scalar/differential record is independently
verified as well. Immutable certified records are checked once at load;
execution cannot silently replace their code or guards.

A raw executor returns its output and the product of its guard coordinates.
For a finite block of consecutive prime-power records, the stage multiplies
guard products modulo n and includes the block's final X and Z in one GCD.
Each next record certifies the previous record's output: by reverse induction,
a unit aggregate certifies all coordinate factors throughout the block.
All intermediate coordinates are then units, so the verified differential
identities apply without exceptional differences. This also proves projective
validity over composite and prime-power moduli under the existing nonsingular
curve precondition.

If the aggregate is nonunit, the block is replayed from its saved input using
the unchanged strict executor. The product may be saturated because different
coordinates expose different factors; replay checks those coordinates
individually and retains the first proper factor. A globally zero coordinate
also triggers this path. The existing finite ladder retry and prime-unit
replay remain available. A bad block is replayed once, without recursively
calling the optimized executor. Thus factor opportunities checked by the old
executor are covered, including factors in coordinates no longer live.

## Execution candidates and bounds

The shared certified representation supports a decoded tuple interpreter,
straight-line calls to the pinned kernels, and straight-line code containing
the exact pinned kernel expressions in the same operation order. Generated
source contains only fixed templates and independently validated numeric
register indices. Source text is never accepted from a catalog. PRAC, upstream
Lucas and binary prime-power controls use the same executors and proof rules.

B1 remains <=2,000; each catalog family has <=512 immutable records,
<=512 operations per record, <=16 point registers and <=1 MiB catalog input.
Generated source is capped at 128 KiB per function and 8 MiB per owned program.
A block contains at most 64 records. Guard products inside a record have at
most 2*(512+1)*bit_length(n) bits; the cross-record accumulator is reduced
modulo n each time. There is no global point or unbounded code cache.
Original block points, scalar products and arithmetic temporaries are charged
separately from the register limit. Cold construction, compilation, certificate
verification, code size and reuse are measured independently.

The comparison screens the three chain families, three executors,
and batches 1/16/64 on frozen training inputs. It uses fresh held-out
certified inputs for confirmation of qualifying candidates, retaining the
pinned whole-lcm ladder and original strict executor. Profiles remain separate
from performance evidence. New CF search stays conditional on a credible
opportunity in the optimized comparison; original negative results do not
veto this explicitly authorized follow-up.

## Frozen first screen

Commit `17aca0b` freezes the optimized common-executor sources and independent
held-out corpus. The GMP-enabled PyPy suite passes 421 tests and full lint.
The completed screen contains 58 groups: three chain families, three execution
modes, three batch sizes, plus both old strict controls, on each backend.
Every group has >=3 seconds of validated warmup per arm and >=9 paired samples.
One GMP inline-Lucas group extends to 18 samples; all final groups meet the
frozen stability rule. The raw capture hash is
`3f23585d99a37fc4cff2c1feb83d90a8828293436adc31cefd07d1f22bc734b4`.

The selection rule favors the simpler executor and smaller batch within 1%
of the best stable paired ratio. Selected stage-reuse training ratios are:

| Family | int choice / ratio (95% interval) | GMP choice / ratio (95% interval) |
| --- | --- | --- |
| Binary prime powers | tuple/16: 1.244 [1.216, 1.254] | tuple/16: 1.114 [1.072, 1.164] |
| PRAC | tuple/16: 0.879 [0.871, 0.888] | tuple/16: 0.888 [0.844, 0.909] |
| GMP-ECM Lucas | tuple/64: 0.867 [0.864, 0.873] | tuple/16: 0.855 [0.838, 0.880] |

Ratios use the pinned whole-lcm ladder as denominator; smaller is faster.
These are training-stage results, not production or complete-campaign claims.
The selected binary prime-power family remains a losing diagnostic control.
Straight-line expansion does not win this screen. Profiles must be kept
separate before attributing its loss to JIT behavior or call overhead.
The largest observed worker RSS is 96,223,232 bytes; that is process peak
memory, not the size of a retained program or an extra-memory estimate.

The input catalog is 202,461 bytes and holds 333 records per family, covering
all prime powers <=2,000. The complete catalog contains 6,025/4,507/4,483 point
operations and 6,048/1,912/2,046 guard coordinates for binary/PRAC/Lucas.
Those are catalog totals, not one stage's executed schedule; unused powers
and recovery-unit records must not be counted as ordinary-stage arithmetic.

## Conditional continued-fraction extension

The stable optimized PRAC and Lucas gains open the new bounded search gate.
The versioned `c6_cf_gate.json` records the evidence and finite selection rule
before any new search is executed. This explicitly supersedes the original
conservative executor's empirical search stop, without editing its captures.

The paper links its own [dacbench-20240609 release](
https://cr.yp.to/2024/dacbench-20240609.tar.gz). Its SHA-256 is
`9319a21b30425d68363c0a2f1a9f375a4745e9fd5274a9c4d6aaf942468ce2bf`.
The upstream README offers several permissive alternatives; this experiment
uses CC0-1.0 and retains the authors' attribution and license text in
`inputs/upstream/c6_dacbench/`. No license is inferred from the paper alone.

The adapter uses the published incremental-length, Fibonacci-pruned search
from Section 3.5. The immutable source is copied to isolated scratch; its
single floating-point floor expression is replaced by exact integer division.
External guards cap target primes at 2,000, depth at 18 bits, search nodes at
10 million, wall/CPU time at 60 seconds, process RSS at 512 MiB and output at
16 MiB. Upstream's threaded benchmark driver is never invoked. At this bound,
meet-in-the-middle tables do not have a demonstrated need and remain unrun.

The decoder independently reconstructs the differential instructions from
upstream integer chains. A different verifier enumerates every coprime
terminal pair `a<b`, `a+b=p` and runs the unique reverse Euclidean path to
`(1,2)`, checking the claimed minimum without trusting forward-search pruning.
This establishes minimum length only inside the defined CF family for the
bounded prime records. Repeating a prime chain for a prime power is verified
composition; it is not an optimality claim for that power or for the stage lcm.

Both the common tuple executor and Algorithm 1's specialized three-point
executor use identical CF arithmetic and the same guard certificate. The
specialized metadata is checked against the certified compact record before
execution. The known-difference identities, all intermediate coordinate
factors, overwritten points, final scalar and finite strict recovery remain
covered. Three persistent working points exclude the saved block/recovery
point, scalar guard accumulator and arithmetic temporaries.

## CF search and common-executor screen results

The first bounded upstream search used 239,243 nodes; independent reverse
verification used 1,923,503 nodes. The initial fresh search process took
0.119 seconds and independent verification 0.056 seconds; repeated generation
costs are reported separately below. The generated output is pinned at SHA-256
`ad3b413d76675158d4ccd7d7604a22a2c21d6aa1e4f5cb185c7fc818ebd1c2d0`.
The 303 prime records compose into 333 independently certified prime-power
records. No search is performed during a factoring attempt.

All twelve CF-screen groups settled at nine samples. The frozen choices are
`tuple/16` for int, ratio 0.867 [0.865, 0.874], and `tuple/64` for GMP,
0.875 [0.833, 0.887], against the original ladder. The distinct three-point
executor did not win the prespecified within-1% simplicity rule. These are
training-stage results. The raw capture hash is
`21d63d361d3764c7b11fdc8b2d1de98afdbb637a397a9c390d9b4f7677f883ee`.

CF uses 3,982 additions and 333 doublings per stage, 4,315 operations total.
GMP-ECM Lucas uses 3,437 additions and 864 doublings, 4,301 operations total.
Under the illustrative 6A+5D weighting, CF costs 25,557 versus 24,942, or 2.47%
more. This weighting is not a calibrated PyPy cost model. Exact minimality
within CF therefore supplies no arithmetic or runtime dominance over Lucas.
No larger search, meet-in-the-middle table or full-lcm search is justified by
this bounded chain-quality comparison.

## Combined B4 comparison, frozen before confirmation

The user's newly integrated B4 control is pinned separately at `a521573`
(native integration `0b568ec`). Its complete `ecm.py` is an immutable versioned
input loaded under a separate module name; C6 production sources remain at
`bcf5f3d`. The B4 change fuses ladder addition/doubling and reduces four
selected intermediates. Its historical 6.79% complete-factoring saving cannot
be added to C6 stage-one percentages.

The bounded new screen crosses PRAC, GMP-ECM Lucas and CF records with late
kernels, early-reduction kernels, and an adjacent independent D/A fusion pass,
using batches 16/64. All variants retain the same certified records and guard
masks. Fusion requires that the second instruction does not read the first
instruction's destination, the destinations differ, and the double input is
one addition input. Both outputs are computed before writes and both masks
are charged. Thus exact canonical coordinate residues and the coverage proof
are unchanged, including over nonfields. Strict replay retains its original
kernels. Independent tests cover every retained record, affine action,
nonunits, prime squares, split saturation and the published false-infinity
input. The new stage-two control uses B4's scalar actions; early point
reductions apply only to experimental stage one.

The paired control is B4 alone. The original ladder is also a confirmation
arm. Each family's stable screen winner is frozen before the unused heldout
corpus is timed; ties within 1% prefer late, reduced, then fused, followed by
the smaller batch. `c6_b4_protocol.json` pins sources and the unchanged
bounds, inputs, seeds, budgets and sampling policy. GMP stays in the separate
original-kernel comparison because B4 retained its readable GMP control.

Separate JIT diagnostics compare three seconds and an additional twenty
seconds of validated warmup. Tuple execution remains faster in both phases;
extra warmup does not rescue generated inline code. The initial optional
snapshot API aborted on this PyPy build and is excluded. The revised hooks
capture compile roots and abort reasons, not every inlined function, and
instrumented times are not acceptance evidence. GMP inline tracing reports
`ABORT_TOO_LONG`; this is diagnostic evidence, not proof of the entire cost
breakdown. See [PyPy's JIT-hook documentation](https://doc.pypy.org/jit-hooks.html).

The combined screen finishes all twenty groups; two instability extensions
use eighteen samples and all final captures pass the frozen spread rule.
The raw capture hash is
`5eaad3a29bf287358769e3dc1cc21503d442c50ab619dcac3c9e69504ece2c41`.
Commit `8ba1d60` freezes the following choices before heldout timing:

| Family | Chosen kernel / batch | Training ratio to B4 (95% interval) |
| --- | --- | --- |
| PRAC | reduced / 16 | 0.900 [0.890, 0.902] |
| GMP-ECM Lucas | reduced / 64 | 0.888 [0.879, 0.904] |
| CF | fused / 64 | 0.876 [0.847, 0.900] |

The full PyPy/GMP suite passes 431 tests and full lint. The initial lint run
found a formatter/pycodestyle slice disagreement in the report utility; an
explicit midpoint variable resolves it. No timed source changed. Neither
the screen nor the correctness pass closes B3's production acceptance gates.

## Fresh-input confirmation and complete-campaign decision

All 48 heldout groups finish and meet the frozen stability rule. Forty-seven
settle at nine paired samples; GMP CF reused-stage timing extends through
18 to 27. Both rejected shorter captures remain local. Each final arm has
at least three seconds of validated warmup. Choices and timed sources remain
unchanged after training; confidence intervals are the frozen 4,000-resample
paired-median intervals. Intervals describe repeat timing on this bounded
cohort, not an estimate over arbitrary 40–80-digit composites.

The following ratios use the named control as denominator; smaller is faster.
Fresh construction is inside each attempt. Reuse construction is separately
charged below. Absolute times from different workers are not interchangeable.

### Complete stage one

| Control / backend / candidate | Reused ratio [95% interval] | Fresh ratio [95% interval] |
| --- | --- | --- |
| Original / int / binary powers | 1.2161 [1.1340, 1.2413] | 1.8332 [1.8002, 1.9174] |
| Original / int / PRAC tuple/16 | 0.8742 [0.8640, 0.8801] | 1.4206 [1.4112, 1.4352] |
| Original / int / Lucas tuple/64 | 0.8700 [0.8609, 0.9201] | 1.4501 [1.4044, 1.4845] |
| Original / int / CF tuple/16 | 0.8655 [0.8621, 0.8669] | 2.4163 [2.2918, 2.4479] |
| Original / GMP / binary powers | 1.1149 [1.0745, 1.1736] | 1.2833 [1.2487, 1.3402] |
| Original / GMP / PRAC tuple/16 | 0.8723 [0.8500, 0.9155] | 1.0422 [1.0307, 1.0670] |
| Original / GMP / Lucas tuple/16 | 0.8427 [0.8004, 0.8880] | 1.0486 [0.9664, 1.0582] |
| Original / GMP / CF tuple/64 | 0.8472 [0.8369, 0.8762] | 1.3312 [1.2696, 1.3620] |
| B4 / int / original ladder | 1.1519 [1.1478, 1.1541] | 1.1400 [1.1028, 1.1566] |
| B4 / int / PRAC reduced/16 | 0.8894 [0.8708, 0.9000] | 2.1153 [1.9455, 2.1424] |
| B4 / int / Lucas reduced/64 | 0.8790 [0.8647, 0.8816] | 2.0662 [1.9696, 2.1428] |
| B4 / int / CF fused/64 | 0.8703 [0.8695, 0.8773] | 3.1751 [3.0045, 3.2110] |

### Complete two-stage ECM attempts

| Control / backend / candidate | Reused ratio [95% interval] | Fresh ratio [95% interval] |
| --- | --- | --- |
| Original / int / binary powers | 1.1024 [1.0970, 1.1134] | 1.1677 [1.1599, 1.1721] |
| Original / int / PRAC tuple/16 | 0.9265 [0.9211, 0.9278] | 0.9855 [0.9537, 0.9970] |
| Original / int / Lucas tuple/64 | 0.9380 [0.8920, 0.9686] | 0.9912 [0.9743, 0.9965] |
| Original / int / CF tuple/16 | 0.9341 [0.9122, 0.9445] | 1.0799 [1.0695, 1.0961] |
| Original / GMP / binary powers | 1.0413 [0.9959, 1.0458] | 1.0610 [1.0485, 1.0977] |
| Original / GMP / PRAC tuple/16 | 0.9272 [0.9041, 0.9514] | 0.9432 [0.9284, 0.9613] |
| Original / GMP / Lucas tuple/16 | 0.9158 [0.8954, 0.9324] | 0.9403 [0.9267, 0.9572] |
| Original / GMP / CF tuple/64 | 0.9356 [0.8889, 0.9579] | 0.9665 [0.9441, 0.9845] |
| B4 / int / original ladder | 1.0676 [1.0670, 1.0751] | 1.0669 [1.0579, 1.0708] |
| B4 / int / PRAC reduced/16 | 0.9327 [0.9274, 0.9402] | 1.0568 [1.0502, 1.0695] |
| B4 / int / Lucas reduced/64 | 0.9472 [0.9315, 1.0123] | 1.0588 [1.0423, 1.0661] |
| B4 / int / CF fused/64 | 0.9453 [0.9334, 0.9880] | 1.1667 [1.1571, 1.2386] |

Every campaign arm retains the same 10/20 proper splits and ten unresolved
balanced inputs, with no timeout and 103 curves. Each size/shape class has
identical completion in both seeds; no class loses a completion. The controls
make 98 stage-two calls; optimized programs make 97 and perform six finite
block replays. One earlier factor discovery accounts for that stage-two-call
difference and is part of the measured behavior. Every returned factor and
cofactor reconstructs its certified input. Every unsuccessful attempt remains
in the time totals and raw evidence; nine timing repetitions are not 180
independent factor trials. Stage-only cohorts retain all ten unsplit inputs.

The B4-relative PRAC reused-campaign ratio is 0.9327: 6.73% less time with
95% interval 5.98–7.26%. Its CPU ratio is 0.9328; chronological-half medians
are 0.9386/0.9324, seed medians 0.9343/0.9314 and balanced/small10 medians
0.9463/0.8743. Thus the sub-10% gain is not supported solely by one seed,
one input class or one chronological half. Its candidate/control cohort
medians are 0.8503/0.9096 seconds. This is additional benefit against B4,
not the sum of B4's historical portfolio percentage and a C6 stage percentage.

B4-relative Lucas has a favorable campaign median but its interval crosses
one, so that complete-campaign confirmation gate stays unpassed. CF's reused
interval is positive but its first-use and fresh-attempt costs are worse.
Neither family is declared globally faster from these separate paired arms.
The selected CF `fused/64` executor has zero eligible fused D/A pairs in
this schedule: its result reflects reduced kernels and that dispatch layout,
not a realized CF arithmetic-fusion saving. No selected C6 arithmetic-fusion
pass wins over the reduced-only PRAC/Lucas choices.
All combined fresh-attempt candidates lose to B4. Original-kernel int PRAC
and Lucas have small positive aggregate fresh-campaign intervals, but their
small10 class medians regress (1.0203/1.0568), and Lucas has one seed median
above one. Those narrow aggregate results do not justify a broad fresh-run
recommendation, especially against the improved native control.

GMP remains separate: reused Lucas takes 8.42% less campaign time
[6.76%, 10.46%], and fresh-per-attempt Lucas takes 5.97% less
[4.28%, 7.33%]. Both class and seed medians improve. Its reused cohort
medians are 3.4956/3.8663 seconds. This supplies a GMP candidate for B3;
it does not make GMP faster than native integers or change backend defaults.

### First-use charge and decision boundary

The report also adds all recorded program preparation, including the control's
preparation, to the candidate's first reused cohort. This conservative
accounting combines warmed execution with an observed first-use charge; it is
not a cold-process timing or a claim that preparation cost is constant.
For native B4-relative PRAC, 25.068 ms yields ratio 0.9602
[0.9547, 0.9679], retaining a 3.98% first-cohort saving. For GMP Lucas,
18.604 ms yields 0.9207 [0.9002, 0.9371]. Native combined CF's 75.865 ms yields
1.0283 [1.0150, 1.0703]; even its first complete reused cohort does not
amortize that observed preparation. Construction and cold repeats follow below.

Accept the bounded verified reusable executor as the C6 candidate. Recommend
native reduced-kernel PRAC with batch 16 as B3's first integrated-B4 candidate,
and tuple Lucas with batch 16 as the separate GMP candidate. Retain the
ladder for fresh native attempts and all current production defaults. Keep
Lucas/CF records as verified alternatives; defer combined Lucas promotion
and any CF search expansion. Stop binary-power dispatch, generated calls/
inline code, larger/MitM/full-lcm search and three-point promotion at this
bound. B3 must still pass real work-ledger, checkpoint/resume and complete
recursive-portfolio gates; the research program is not production integration.

## Construction, storage, layouts and cold processes

Nine fresh bounded CF generator processes reproduce the same output hash and
239,243 search nodes. Median generator time is 78.149 ms (75.126–79.057 ms);
the parent process's independent minimality verification is 13.685 ms median
(13.421–28.439 ms), always 1,923,503 nodes. The initial separately recorded
first search/verification remains 118.566/56.333 ms. These are offline costs,
not free work inside a factoring attempt. The original GMP-ECM reproduction
remains the separately measured nine-run 526.836 ms median plus a one-time
1.188 s C build; no second large upstream search is run or implied here.

Warmed catalog precomputation from pinned prime chains takes 23.912 ms;
loading and independently verifying all three compact families takes 2.127 ms.
CF catalog decoding/composition/verification takes 6.026 ms and does not rerun
the optimality search. All twenty construction/load/layout groups meet the
spread rule after the prescribed 9/18/27-sample extensions. Ordinary factoring
loads committed records; regenerating catalogs also charges the offline costs.

| Candidate | Warm build (ms) | Estimated full stages to repay warm build / observed first preparation |
| --- | ---: | ---: |
| Original int PRAC tuple/16 | 1.806 | 4 / 32 |
| Original int Lucas tuple/64 | 1.795 | 3 / 28 |
| Original int CF tuple/16 | 6.957 | 11 / 106 |
| GMP PRAC tuple/16 | 1.980 | 1 / 9 |
| GMP Lucas tuple/16 | 1.969 | 1 / 7 |
| GMP CF tuple/64 | 6.783 | 3 / 39 |
| B4 int PRAC reduced/16 | 3.614 | 8 / 56 |
| B4 int Lucas reduced/64 | 3.700 | 8 / 51 |
| B4 int CF fused/64 | 8.102 | 16 / 137 |

These estimates divide each separately measured construction cost by the
median paired saving per full balanced-input stage, then round up. They are
cohort-specific accounting estimates, not measured policy thresholds or cold
factoring times. Construction interleaved with short attempts, partial stages,
JIT effects and stage-two work can change the result; the directly measured
fresh-stage/campaign losses take precedence. The first-preparation column
conservatively includes both arms' preparation. Reusing a catalog across
program construction differs from reusing an already verified program.

The compact JSON catalog is 202,461 bytes and the CF prime catalog 16,292
bytes. Full-stage PRAC/Lucas/CF bytecode occupies 17,268/17,204/17,260 bytes,
with 4,620/4,604/4,618 guard-mask bytes. Their persistent point-register
requirements are 5/4/3; their guarded-coordinate counts are 1,816/1,928/1,280.
All selected variants generate zero source bytes. Reduced programs retain
the extra bounded decoded tuple; CF's plan has at most 17 entries per record
and zero fused pairs. These serialized quantities exclude Python object
headers, strict-recovery objects, saved points and arithmetic temporaries.
Maximum confirmation worker RSS is 91,275,264 bytes; the shared cost process
peaks at 117,637,120 bytes and B4 costs at 89,653,248 bytes. RSS is total process
memory including the JIT, not retained-program size or incremental allocation.

On all 45 eligible upstream primes with identical CF arithmetic over F_1009,
compact/rolling16/three-point medians are 13.417/14.041/11.708 microseconds for
int and 1.598/1.526/1.548 ms for GMP. The GMP compact diagnostic extends to
27 samples. This small-field layout diagnostic includes independent affine
validation and recovery; it does not establish a 40–80-digit full-stage gain.
The full-stage CF training selection still favors the common tuple executor,
so neither the rolling buffer nor the three-point specialization is promoted.

There are 153 fresh-process captures: nine per arm/control, including process
startup, certificate/oracle preparation, ten one-curve attempts, construction
inside each attempt and result validation. Native original ladder/PRAC/Lucas/CF
medians are 0.587/0.650/0.655/0.790 seconds (CF control 0.590 s); GMP values are
0.756/0.792/0.806/0.961 seconds (CF's separate control is 0.747 s).
B4 ladder/PRAC/Lucas/CF medians are 0.578/0.689/0.715/0.846 seconds. Cold runs
are not mixed with warmed samples and do not support a cold-start promotion.

Twelve separate profiles preserve all selected arms. Native reduced PRAC's
profile records 36,880 additions, 6,290 doublings and 270 `math.gcd` calls
across ten stages, including setup/validation. Kernel arithmetic, compact
dispatch and guard products remain visible costs after removing redundant
GCDs. GMP profiles expose wrapper/context-variable activity, with substantial
profiler perturbation; their instrumented times are not a calibrated cost
model or acceptance evidence. The separate longer-warmup JIT diagnostic
also fails to rescue generated inline code. No new optimization or selection
is made after examining heldout profiles.

The 48-group confirmation summary SHA-256 is
`373848b49ae7603f5bbc4c2996b25e17fa8da0ba72210d413ecebc0b5e43177e`;
the diagnostic summary is
`64a5d413d76c223b6f63871c88cec6dab8ce9854b15c6337384dfd8e7c2cfcb8`.
Raw captures, rejected shorter attempts, cold runs, profiles, failed optional
JIT diagnostics, factor/cofactor rows and a per-file SHA-256 manifest remain
local in `v2/benchmarks/results/c6-fast/`. Required generators, protocols,
records, upstream notices and certified corpora are versioned; no raw capture
is required by a test or benchmark loader.

## Reproduction and accounting

Use a GMP-enabled PyPy implementing Python 3.11, the committed inputs and
new output paths. Coordinate the entire sequence with other experiment
owners; each parent runner takes the shared lock and checks for competing
benchmark/test interpreters. Training selections are already frozen; never
regenerate or replace them after viewing confirmation results.

```sh
v2/.venv/bin/python -B -u -m v2.benchmarks.c6_fast_study --confirmation --scope stage_reuse --output v2/benchmarks/results/c6-fast/fast-stage-reuse-new.json
v2/.venv/bin/python -B -u -m v2.benchmarks.c6_cf_study --confirmation --scope stage_reuse --output v2/benchmarks/results/c6-fast/cf-stage-reuse-new.json
v2/.venv/bin/python -B -u -m v2.benchmarks.c6_b4_study --confirmation --scope stage_reuse --output v2/benchmarks/results/c6-fast/b4-stage-reuse-new.json
```

Repeat each runner for `stage_fresh`, `campaign_reuse`, and `campaign_fresh`.
A stage cohort has ten one-curve attempts; a campaign cohort has twenty
complete bounded two-stage ECM attempts, including failures. These are not
recursive portfolio factorizations. The ten certified heldout composites
are distinct from training, with seeds 56839/64758. Input sizes are
40/50/60/70/80 digits; target factors are 10 digits or 20/25/30/35/40 digits.
Each attempt has B1=2,000, B2=147,396, eight curves and twenty seconds wall/CPU.
Balanced cases intentionally test feasible failures, not balanced-80 success.

Reuse means one verified program is owned across the cohort; its construction
is separately reported. Fresh means construct/load/verify once per attempt,
inside its timer, then reuse across that attempt's curves. Timed work includes
conversion, setup, dispatch, intermediate guards, recovery and result
validation. Stage outputs are also compared with independently precomputed
affine targets inside the timed attempt. Certificate checking and computing
the affine targets happen before warmed timing and are included in fresh
process totals. Import/startup and instrumented profiles remain separate.

`c6_fast_costs` repeats bounded CF generation nine times and measures catalog
precomputation, load/verification, construction and identical-arithmetic
compact/ring/three-point layouts. `c6_b4_costs` measures the combined choices.
Use `--confirmation --scope stage_fresh --cold` with each study runner for
nine fresh processes per selected arm and control. `--profiles --scope
stage_reuse` produces separate profile files. `c6_fast_report` validates every
factor/cofactor again and reports repeat-timing intervals, CPU, seeds,
chronological halves, each input-size/shape class, completion and unresolved
values. Repeated timing samples do not create new independent factor trials.

## Additional research boundary

Neill Clift's author-maintained [Lucas-chain notes](https://additionchains.com/Lucas.html)
report extensive length-table enumeration and a corrected integer-overflow
pruning bug. This is an additional primary technical-blog lead, not an
independently certified catalog for this experiment. No multi-gigabyte table,
unspecified-license search implementation or claimed general optimum is
imported. Its claims do not replace the explicit CF-family proof above.
The blog and the GMP-ECM maintainer discussion reinforce the need to verify
pruning assumptions and integer overflow separately from successful sample
chains. The literature/software inspection is bounded, not a claim to have
audited every differential-chain implementation.

## Reusable contract for B3

The candidate is a program of independently certified scalar records, not a
replacement for a production stage job. B3 can build on the following explicit
boundary without importing the research runner's signal timers:

- Own immutable verified records per bounded schedule and backend; never cache
  curve points globally. Identity must cover the scalar schedule, record and
  guard digests, interpreter/kernel version and batch policy. A catalog is
  data requiring verification, not trusted executable code.
- Reserve an entire block's fast arithmetic, guard multiplications/GCD and
  worst-case strict replay before starting it. Retain one original block
  point. The research counters describe executed recovery; they are not the
  current portfolio's scalar-bit work currency and cannot replace its ledger.
- Commit a stage cursor and point only after the aggregate certifies the block
  or strict replay produces a valid continuation. A pending unchecked point
  is not a checkpoint. On nonunit aggregates, replay the exact saved block
  once; preserve a proper factor even when the aggregate itself is n.
- Preserve cancellation/time checks at a documented finite block boundary,
  cumulative allowances, curve/RNG identity and unresolved cofactors. Reject
  incompatible program/checkpoint identities; use an intentional new state
  version rather than reinterpreting an existing ladder checkpoint.
- This API computes a fresh M(B1) stage. Applying it to an already completed
  stage multiplies by M(B1) again; increased-B1 continuation requires the
  exact schedule ratio and its separately verified multiplicities. Coordinate
  that contract with A6; this study implements no ECM bound continuation.
- Keep construction and verification charged on cache misses and resumed
  reconstruction. Bound retained programs and metadata explicitly. Reduced
  execution retains at most two bounded decoded-operation tuples per action;
  fusion adds at most one plan row per original operation. Shared immutable
  records are not duplicated curve states. Serialized code/mask byte counts
  and process RSS are reported separately from Python object memory.
- Recheck complete recursive portfolio behavior against integrated B4 with
  real work reservations and checkpoint/resume enabled. This study's complete
  two-stage ECM attempts establish only the documented candidate scope.

No source or result here changes production defaults, portfolios or checkpoint
formats. The study was frozen on a separate C6 branch before its experimental
code and records were integrated into mainline.

## Committed-only acceptance

A `git archive` of code/input commit `a7642c2` passes 431 tests under the
GMP-enabled PyPy 7.3.23 / Python 3.11.15 environment (26.548 s), full lint
(164 Python files) and all 84 benchmark imports. Required original, optimized,
CF and B4 source hashes, 303 upstream Lucas records, every 333-record catalog
and both certified corpora pass from committed files only. The archive uses
only an external environment symlink; no ignored evidence supplies an import
or test input. Subsequent changes in the acceptance commit are documentation.

Commands are `make -C v2 test PYTHON="$PWD/v2/.venv/bin/python"` and
`make -C v2 lint`, followed by the committed-only import/hash/corpus checks.
All timing and heavy validation is serial and completed before releasing the
shared timing window. C6 closes only its bounded research/experimental-candidate
gate. `v1/`, production arithmetic, stage jobs and original immutable controls
remain at the original C6 base; the integrated B4 source is a separate pinned
benchmark input. Those statements describe the accepted C6 checkout before
integration. Mainline now contains the experimental PRAC, Lucas and CF code
and immutable records, while its production ECM path remains the B4 ladder.
The historical runners deliberately reject changed mainline source hashes;
reproduce their frozen measurements from commit `65cf845`, or freeze a new
protocol against current mainline for a new comparison. Do not relabel the
historical ratios as fresh mainline timing evidence.
