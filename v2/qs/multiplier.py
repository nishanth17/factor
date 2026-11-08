"""Finite integer Knuth-Schroeppel-style scoring; no native digit table."""

from dataclasses import dataclass
from fractions import Fraction

from ..common import prime_sieve, utils
from ..common.arithmetic import gcd, pow
from ..execution.budget import Budget
from ..execution.resources import removed_arguments
from .factor_base import checked_target

DEFAULT_MULTIPLIERS = (1, 3, 5, 7, 11, 13, 15, 17, 19, 21, 23, 29, 31)
PARITY_MULTIPLIERS = (
    1,
    2,
    3,
    5,
    6,
    7,
    10,
    11,
    13,
    14,
    15,
    17,
    19,
    21,
    22,
    23,
    26,
    29,
    30,
    31,
)


@dataclass(frozen=True)
class MultiplierChoice:
    """Selected multiplier, proper GCD split or finite scored controls."""

    multiplier: int
    divisor: int | None
    scores: tuple


def _log_weight(value):
    """Floor 256*log2(value) using bounded integer powering."""
    return (value**256).bit_length() - 1


def select_multiplier(
    n, *, candidates=DEFAULT_MULTIPLIERS, budget=None, policy="standard"
):
    """Score supported residues/modulo-8 with a half-log size penalty.

    This is an integer-scaled finite heuristic. It is independently implemented
    and calibrated locally; it does not reproduce a native implementation's
    floating score or import its parameter tables. h=1 remains a control.
    """
    checked_target(n, 1)
    if policy not in ("standard", "parity2"):
        raise ValueError("unknown multiplier scoring policy")
    if not isinstance(candidates, tuple) or not 1 <= len(candidates) <= 64:
        raise ValueError("multiplier candidates must be a finite tuple")
    if len(set(candidates)) != len(candidates):
        raise ValueError("multiplier candidates must be distinct")
    for value in candidates:
        utils.require_integer(value, "multiplier", 1)
        if value > 255 or (policy == "standard" and value % 2 == 0):
            raise ValueError("unsupported scored multiplier")
        if any(value % (p * p) == 0 for p in (2, 3, 5, 7, 11, 13)):
            raise ValueError("scored multipliers must be squarefree")

    budget = budget if budget is not None else Budget()
    budget.consume(len(candidates) * (n.bit_length() + 8192))
    primes = prime_sieve.small_sieve(97)
    weights = {prime: _log_weight(prime) for prime in primes}
    scores = []

    for multiplier in candidates:
        divisor = gcd(n, multiplier)
        if utils.valid_divisor(divisor, n):
            return MultiplierChoice(multiplier, divisor, tuple(scores))
        if divisor != 1:
            continue
        target = n * multiplier
        residue8 = target % 8
        two_weight = {1: 2, 5: 1, 3: Fraction(1, 2), 7: Fraction(1, 2)}
        if policy == "parity2":
            # Stored 2P has the same target magnitude as standard F. At
            # D=1 mod 8 its expected 2-valuation is 3; at D=5 it is 1.
            # Squarefree even k uses standard F, with valuation 1 on half
            # the positions and zero elsewhere (no roots modulo 4).
            two_weight.update({1: 3, 2: Fraction(1, 2), 6: Fraction(1, 2)})
        score = Fraction(-_log_weight(multiplier), 2)
        score += 256 * two_weight[residue8]
        for prime in primes:
            if prime == 2:
                continue
            residue = target % prime
            if residue == 0:
                score += Fraction(weights[prime], prime)
            elif pow(residue, (prime - 1) // 2, prime) == 1:
                score += Fraction(2 * weights[prime], prime - 1)

        scores.append((multiplier, score.numerator // score.denominator))

    if not scores:
        raise ValueError("no nonsingular scored multiplier")
    chosen = max(scores, key=lambda item: (item[1], -item[0]))[0]
    return MultiplierChoice(chosen, None, tuple(scores))


# A separate policy preserves the original candidate list, arithmetic and
# charges. The sieve endpoint is exclusive, so 1988 includes prime 1987.
WIDE_MULTIPLIERS = (
    1,
    2,
    3,
    5,
    6,
    7,
    10,
    11,
    13,
    14,
    15,
    17,
    19,
    21,
    22,
    23,
    26,
    29,
    30,
    31,
    33,
    34,
    35,
    37,
    38,
    39,
    41,
    42,
    43,
    46,
    47,
    51,
    53,
    55,
    57,
    58,
    59,
    61,
    62,
    65,
    66,
    67,
    69,
    70,
    71,
    73,
    74,
    77,
    78,
    79,
    82,
    83,
    85,
    86,
    87,
    89,
    91,
    93,
    94,
    95,
    97,
)
WIDE_PRIME_BOUND = 1988
WIDE_MULTIPLIER_WORKSPACE = 8 * 2**20


def _legendre(value, prime):
    residue = pow(value % prime, (prime - 1) // 2, prime)
    return 0 if residue == 0 else 1 if residue == 1 else -1


@removed_arguments(("memory_bytes",))
def select_coupled_multiplier(n, *, budget=None):
    """Rank the fixed wider parity population with local exact-symbol reuse.

    All fractions share an exact denominator. The returned integer scores
    equal flooring the direct rational score, including negative results.
    No cache survives this call or crosses input/backend identities.
    """
    checked_target(n, 1)
    budget = budget if budget is not None else Budget()
    budget.consume(WIDE_PRIME_BOUND + 32 * len(WIDE_MULTIPLIERS))
    primes = prime_sieve.small_sieve(WIDE_PRIME_BOUND)[1:]
    factors = prime_sieve.small_sieve(98)
    support = tuple(
        tuple(i for i, prime in enumerate(factors) if k % prime == 0)
        for k in WIDE_MULTIPLIERS
    )

    # 2 and every p,p-1 divide scale. Its bit length is at most 6601,
    # independently of N. Bound allocations before computing the denominator.
    scale = 2
    for prime in primes:
        budget.consume(n.bit_length() + 4 * 6601)
        divisor = gcd(n, prime)
        if utils.valid_divisor(divisor, n):
            return MultiplierChoice(1, divisor, ())
        for denominator in (prime, prime - 1):
            scale *= denominator // gcd(scale, denominator)

    cache = []
    for prime in primes:
        budget.consume((len(factors) + 1) * prime.bit_length() ** 2 + 4096)
        weight = _log_weight(prime)
        cache.append(
            (
                _legendre(n, prime),
                tuple(_legendre(q, prime) for q in factors),
                weight * (scale // prime),
                2 * weight * (scale // (prime - 1)),
            )
        )

    scores = []
    words = (scale.bit_length() + 63) // 64
    for multiplier, indices in zip(WIDE_MULTIPLIERS, support):
        budget.consume(n.bit_length() + 8192)
        divisor = gcd(n, multiplier)
        if utils.valid_divisor(divisor, n):
            return MultiplierChoice(multiplier, divisor, tuple(scores))
        if divisor != 1:
            continue
        residue = n * multiplier % 8
        two_numerator = 6 if residue == 1 else 2 if residue == 5 else 1
        score = -_log_weight(multiplier) * (scale // 2)
        score += 128 * two_numerator * scale
        for symbol, factor_symbols, singular, nonsingular in cache:
            budget.consume(words + len(indices) + 1)
            value = symbol
            for index in indices:
                value *= factor_symbols[index]
            if value == 0:
                score += singular
            elif value == 1:
                score += nonsingular
        scores.append((multiplier, int(score // scale)))

    if not scores:
        raise ValueError("no nonsingular scored multiplier")
    chosen = max(scores, key=lambda item: (item[1], -item[0]))[0]
    return MultiplierChoice(chosen, None, tuple(scores))
