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
