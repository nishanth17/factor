# C10a-01 opt-in arithmetic contract

This is an arithmetic implementation, not a default or performance promotion.
`polynomial_policy="parity2"` adds an odd-coordinate SIQS family; the existing
`"standard"` policy remains the literal control. `multiplier_policy="parity2"`
requires that family policy and is selected only when `multiplier=0`.
A fixed positive multiplier bypasses ranking, including supported even values.

## Representation and coordinates

The existing polynomial is F(x) = A x² + 2B x + C, where C = (B² − D)/A,
D = kN, U = Ax + B, and U² − D = AF. Division by A already existed.
The new policy applies when D ≡ 1 (mod 4). Let a be a squarefree product of
odd, nonsingular factor-base primes. Choose odd B with B² ≡ D (mod 4a),
centered in [−a, a). Store A = 2a and F = 2P, where

    P(x) = a*x*x + B*x + (B*B-D)/(4*a)
    U(x) = 2*a*x+B
    U(x)**2-D = 4*a*P(x) = A*F(x).

All divisions are exact. The family CRT terms modulo 2a are even; adding a
to the initial signed sum supplies its odd lift. Changing a Gray sign changes
B by an even amount. Reduction modulo 2a preserves oddness and all odd-prime
roots. Fixing one sign chooses one of each pair B, −B; the existing adjacent
coordinate permutation changes traversal only. No square root or inversion
uses floating point.

For the standard root b0 modulo a, write B = b0 + aδ. The exact mapping is
y = 2x + δ. Standard F at y equals 4P(x), while stored parity F equals 2P(x);
U agrees. Thus corresponding physical trials have equal residuals after
removing factor-base powers. The selected parity class is only half of the
standard coordinates; the full streams are not identical.

Both policies visit the integer interval [−M, M], including its final point:
2M + 1 positions, divided into half-open blocks. A target uses the actual
half-width, not a block length. Standard A targets floor(sqrt(2D)/M); parity
odd a targets floor(sqrt(2D)/(2M)). Stored 2a has the former target magnitude.
Changing a changes the population deliberately. With b small, balancing the
endpoint and central polynomial magnitudes yields these targets; it does not
prove an optimal smoothness rate. Physical U strides are A and 2a.

## Domain, roots and exceptional cases

The core supports odd N ≥ 3 of at most 4096 bits and bounded positive k.
Primes and powers remain diagnostic inputs; production preprocessing and
certainty labels retain their prior contract. D ≡ 3 (mod 4) and even D use
the standard family. Even k is never passed to the odd-discriminant formula.
The new selector enumerates the squarefree integers at most 31 (20 choices,
including 2, 6, 10, 14, 22, 26 and 30); the old 13 odd choices are unchanged.
Forced nonsquarefree k keeps the existing explicit-multiplier semantics but
is outside the new score's squarefree heuristic derivation.

For odd p not dividing A, roots are (±sqrt(D) − B)/A modulo p. For p dividing
D this is one repeated root; when p divides squarefree k but not N it does
not lift to p². For p dividing a, B is a unit and the root is −C/(2B).
A primes dividing D are excluded before inversion. Characteristic two is
enumerated, including all-position roots, and prime powers are lifted with
explicit singular branches. No inverse is assumed at a nonunit.

For parity F, v2(F) = 1 at D ≡ 5 (mod 8); at D ≡ 1 (mod 8), v2(F) ≥ 2.
The two factors of 2 in AF are retained, one from A and at least one from F.
They are not discarded or represented by an inverse modulo N. All exponents
include the complete A contribution; the relation sign is the sign of F.
Negative values use absolute-value division and a separate sign column.
Existing known square coefficients and paired-large-prime corrections retain
their exact meaning. Dependency extraction reconstructs original atoms,
checks even totals and X² ≡ Y² (mod N), then validates both GCD signs.
Every actually exposed proper GCD is recovered; no claim is made that the
collector computes a GCD for every rejected nonsmooth position.

Zero polynomial values use the existing U/N GCD recovery. Multiplier and
base GCDs precede family construction. A nonunit with no proper split is
refused, not inverted. Recentring changes the coordinate offset in the exact
mapping; boundary/tail positions and negative windows are included in tests.
Invalid roots, unavailable A sets and exhausted finite families refuse or
terminate through the existing checked-cofactor paths.

## Coupled multiplier heuristic

Use the existing integer approximation floor(256 log2(p)) and half-log k
penalty. Odd nonsingular prime contribution is 2 log(p)/(p−1); a prime
solely dividing squarefree k contributes log(p)/p because its root cannot
lift. For standard F the expected 2-adic valuation is 2, 1, 1/2, 1/2 for
D mod 8 = 1, 5, 3, 7. Squarefree even k gives 1/2. Conditioning U to be odd
and subtracting the one factor 2 in stored A gives parity F weights 3 and 1
for residues 1 and 5. The new score uses those weights plus standard fallback
weights. It does not copy a native candidate table or a native tuning result.
Finite windows, A-prime selection and downstream graph utility make this a
ranking hypothesis, not a theorem of complete-factor performance.

## Derived state and compatibility

Polynomial identity already binds N, k, A, B and square coefficient. Family
identity additionally binds the parity convention; base identity binds ordered
primes and roots. Root, plan, presieve and sparse-queue caches are collector
local, reset on polynomial switches, and keyed by polynomial/window/policy
as appropriate. Backend and full configuration are checked on job restoration.
No coefficient-free cache is transferred between representations.

New jobs serialize as SIQS schema 8, canonical parity families as schema 3,
and swapped parity families as `c10a01-coordinate-family-v2`. Standard jobs
keep their previous schema and omit the two new default-valued fields.
Legacy schema injection and policy substitution during resume are rejected,
including monotone resource extensions. Interrupted private construction
leaves the unpublished family index/pending step intact; rebuilding is charged.
Existing storage and work reservations remain binding.

Independent tests use Cartesian CRT, direct integer norm evaluation, the
existing discriminant oracle, exhaustive small modular roots and high powers,
direct trial division, original-atom extraction and interrupted checkpoints.
They include both odd residue branches, even/odd standard fallback, queues,
presieving, plans, tails, cancellation and exact factor reconstruction.
