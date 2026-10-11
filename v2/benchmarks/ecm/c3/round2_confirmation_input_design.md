# C3 confirmation input preparation

This code prepares a stronger generation mechanism while the current frozen
calibration runs. No new confirmation integers have been generated or revealed.
Do not use them until the fitted selection and a separate confirmation protocol
have been committed. This helper does not change the 57 training subjects or
their source freeze.

The existing p=kq+1 construction is nonuniform and can overrepresent convenient
p-1 structure at smaller factor sizes. UniformPrimeSource instead draws odd
integers uniformly within the requested exact decimal band, filters composites,
and constructs a recursive full n-1 primality proof for each survivor. It covers
up to 20-digit factors, including both primes of the balanced 30/40-digit
strata. Larger partners may still use the existing Pocklington generator; that
remaining bias must be disclosed in the eventual corpus and acceptance report.
Decimal rejection on the final product also conditions the joint distribution.
Do not claim a fully uniform semiprime population.

Proof generation reuses the existing bounded p-1/rho portfolio without ECM.
It has a shared 150-billion-work / 180-wall-second / 180-CPU-second allowance,
one million random draws, 10,000 stored proof nodes and a finite 256-base
witness search. These are generation allowances, not benchmark service grants.
An unfinished factorization/proof ends generation. It never discards that
candidate and silently samples an easier one; a failed frozen generation must
be reported, not retried with an undeclared seed or increased allowance.

Probable-prime filtering and terminal factorer labels do not prove primality.
The independent checker uses only integer products, exact trial division,
modular powers and gcds. It verifies all recursively required prime factors,
full reconstruction of n-1, a common witness of order n-1, or the existing
partial Pocklington inequality and witness. Unknown proof types are rejected.
No expected factor or certificate enters portfolio dispatch.

The mathematical basis is the authors' Handbook of Applied Cryptography,
1996 edition, section 4.3.2, Facts 4.38 and 4.40, printed pages 143-144;
section 4.7 notes explains certificate verification and witness refinements:
https://cacr.uwaterloo.ca/hac/about/chap4.pdf
This implementation was written from those integer conditions; no upstream
code or book text is copied into the repository.

Before any generation, the confirmation freeze must still declare its seed,
exact strata/counts, disjointness checks, common resource regime, matching,
instability extensions, total expense cap and adopt/retain/defer criteria.
This preparation is not that freeze and is not evidence of a policy benefit.

Verification before any fresh generation: all 621 PyPy tests and full lint
pass. A fixed 20-digit candidate, 18446744073709551557, independently verifies
through eleven proof nodes and 228,125 generation work units within a
predeclared 15-second feasibility cap. This is a functional proof check, not
an accepted performance measurement or a new confirmation subject. The raw
proof and QA receipt remain in results/c3/round2/proof-qa/.

The functional probe can be repeated with PyPy 3.11:

```python
from v2.benchmarks.ecm.c3.round2_prime_proofs import (
    UniformPrimeSource, verify_proofs,
)
from v2.execution.budget import Budget
source = UniformPrimeSource(
    17, budget=Budget(work_limit=150_000_000_000, seconds=15, cpu_seconds=15)
)
source.prove(18446744073709551557)
verify_proofs(source.certificates)
```
