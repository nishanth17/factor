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
