# Shared exact arithmetic

Arithmetic, preprocessing, prime generation and primality classification used
by the factoring engines. Run on PyPy implementing Python 3.11.

| Module | Responsibility |
| --- | --- |
| [arithmetic.py](arithmetic.py) | Native/GMP integer boundary, exact helpers and canonical outputs |
| [preprocessing.py](preprocessing.py) | Perfect powers, integer roots and reusable classifications |
| [prime_sieve.py](prime_sieve.py) | Prime enumeration with half-open endpoints |
| [utils.py](utils.py) | GCD/inversion, divisors, RNG helpers and primality certainty |

## Arithmetic backends

`python-int` is the dependency-free default. Explicit `gmpy2-mpz` selection
keeps long-lived arithmetic values as `mpz`; indices, seeds and accounting
remain Python integers. High-level results and checkpoint payloads use
canonical Python integers. No callback is added per multiplication.

Install the optional dependency into the chosen PyPy environment:

```sh
v2/.venv/bin/python -m pip install -r v2/requirements-gmp.txt
v2/.venv/bin/python -m v2.factor 626100403 --bounded --backend gmpy2-mpz
```

A missing dependency is an explicit error. Native calls do not import GMP;
there is no silent backend or CPython fallback and no automatic size crossover.
A portfolio and its relation fallback must select the same backend.
Checkpoint identities bind the selected backend and, for GMP, its build;
incompatible restores are rejected. Faster isolated operations do not imply
faster complete factoring. See [results](../benchmarks/RESULTS.md).

For low-level calls, use `arithmetic.get_backend(name).integer(value)` once
at the boundary. `pow` supports exact and modular powers; `divexact` checks
divisibility. `utils.gcd`, `extended_gcd`, `modular_inverse`, `isqrt` and
`preprocessing.integer_root` preserve exact inputs. Failed inversion raises
`NonInvertibleError`, whose `divisor` retains the GCD, including saturation.
Factoring arithmetic does not use floating division.

## Primality and certainty

`classify_prime` returns composite/probable/proven classification;
`is_prime` and `is_prime_fast` are Boolean convenience functions.
`deterministic_bases(n)` returns the supported witnesses or `None`; callers
handle small-input and divisibility cases first.

| Deterministic range | Witness contract |
| --- | --- |
| Below `2**64` | Existing smaller-range tests and seven-base word-domain test |
| Below `318665857834031151167461` | First 12 prime bases, through 37 |
| Below `3317044064679887385961981` | First 13 prime bases, through 41 |

These are strict upper bounds; the last two endpoints are composite strong
pseudoprimes to the preceding witness sets. Their guarantees use the
[recorded exhaustive computational result](https://arxiv.org/abs/1509.00864).
Above the supported range, surviving tests remain probable. Certificate
proofs remain separately owned work in the [roadmap](../ROADMAP.md).

Explicit `use_probabilistic=True` draws the requested `tolerance` witnesses
for nontrivial survivors, with early exit on a composite witness. Fixed ranges
otherwise consume their deterministic witness count and no RNG draws.
`PortfolioConfig.primality_rounds` controls random tests outside its fixed
policy, not a request to override deterministic tests.

Fresh bounded jobs pin `mr13-strict-v1`. Older saved jobs retain their policy,
RNG/work identity and labels; unknown policy/witness identities are rejected.
No completed probable factor is silently upgraded by resume. Result
completion and proof are distinct; see [API contracts](../docs/API.md).

## Endpoint and verification contracts

Prime iterators use `[lo, hi)` and emit no upper endpoint. p−1/ECM bounds are
inclusive, so callers translate endpoints explicitly. Perfect-power tests
verify integer roots exactly. Every accepted divisor is strictly between one
and n and divides n; every factorization reconstructs unresolved pieces too.
