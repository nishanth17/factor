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

The first comparison will screen the three chain families, three executors,
and batches 1/16/64 on frozen training inputs. It will use fresh held-out
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
formats. The C6 branch remains separate and unmerged.
