"""Finite integer Knuth-Schroeppel-style scoring; no native digit table."""

from dataclasses import dataclass
from fractions import Fraction

from ..common import prime_sieve, utils
from ..common.arithmetic import gcd, pow
from ..execution.budget import Budget
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
