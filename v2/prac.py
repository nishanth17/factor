"""Bounded, independently verified differential chains for experimental ECM.

PRAC's nine reductions follow Montgomery's Table 4 (also in GMP-ECM's
ecm.c). Records describe arithmetic, not the reduction algorithm: a separate
integer interpreter proves every differential relation and the final scalar.
See ROADMAP.md P4.1 for references, exceptional states and promotion limits.
"""

from dataclasses import dataclass
from functools import lru_cache

from . import utils
from .arithmetic import gcd

MAX_SCALAR_BITS = 32
MAX_CHAIN_STEPS = 512
CACHE_SIZE = 512
RATIO_DENOMINATOR = 10**17
# Exact rational approximations to GMP-ECM's continued-fraction choices.
RATIO_NUMERATORS = (
    61803398874989485,
    72360679774997897,
    58017872829546410,
    63283980608870629,
    61242994950949500,
    62018198080741576,
    61721461653440386,
    61834711965622806,
    61791440652881789,
    61807966846989581,
)


@dataclass(frozen=True)
class Chain:
    """Immutable SSA record; node 0 is P, each instruction appends a node.

    (i, i, -1) doubles node i. (i, j, h) applies differential addition with
    known third node h; h can represent either the sum or the difference.
    Only the zero scalar has output -1, denoting infinity.
    """

    scalar: int
    instructions: tuple
    output: int
    method: str
    split: int = 0

    def cost(self, add_cost=6, double_cost=5):
        """Abstract 4M+2S addition / 3M+2S doubling costs, not timings."""
        return sum(
            double_cost if difference == -1 else add_cost
            for _, _, difference in self.instructions
        )


class NonunitPointError(ArithmeticError):
    """A proper factor is available, or a curve must be retried.

    A factor result is never returned as a projective point. Callers can use
    ``factor`` after their usual divisor validation; None means retry.
    """

    def __init__(self, factor=None):
        self.factor = factor
        super().__init__(
            f"projective arithmetic exposed factor {factor}"
            if factor is not None
            else "degenerate projective arithmetic; retry the curve"
        )


class _ExceptionalDifferenceError(Exception):
    """A chain needs one retry from its original point using the ladder."""


def verify_chain(chain):
    """Verify bounded records without trusting PRAC reductions or metadata.

    X coordinates identify signed multiples. Given positive integers a,b,c,
    f(aP,bP,cP) represents a+b if c=|a-b|, or |a-b| if c=a+b. Checking this
    at every instruction is stronger than checking just the final scalar.
    """
    if not isinstance(chain, Chain):
        raise ValueError("expected a Chain record")
    utils.require_integer(chain.scalar, "scalar", 0)
    if chain.scalar.bit_length() > MAX_SCALAR_BITS:
        raise ValueError("scalar exceeds the chain limit")
    if (
        type(chain.instructions) is not tuple
        or len(chain.instructions) > MAX_CHAIN_STEPS
        or type(chain.output) is not int
        or type(chain.split) is not int
        or chain.split < 0
        or type(chain.method) is not str
        or chain.method not in ("binary", "prac", "lucas")
    ):
        raise ValueError("invalid chain structure")
    if chain.split:
        odd = (
            chain.scalar >> ((chain.scalar & -chain.scalar).bit_length() - 1)
            if chain.scalar
            else 0
        )
        if not odd // 2 < chain.split < odd or gcd(odd, chain.split) != 1:
            raise ValueError("invalid PRAC split")
    if chain.scalar == 0:
        if chain.instructions or chain.output != -1 or chain.split:
            raise ValueError("invalid zero chain")
        return True

    values = [1]
    for instruction in chain.instructions:
        if (
            type(instruction) is not tuple
            or len(instruction) != 3
            or any(type(index) is not int for index in instruction)
        ):
            raise ValueError("invalid chain instruction")
        left, right, difference = instruction
        if not 0 <= left < len(values) or not 0 <= right < len(values):
            raise ValueError("chain references an unavailable node")
        a, b = values[left], values[right]
        if difference == -1:
            if left != right:
                raise ValueError("doubling must use the same node")
            value = 2 * a
        elif 0 <= difference < len(values):
            c = values[difference]
            if c == abs(a - b):
                value = a + b
            elif c == a + b:
                value = abs(a - b)
            else:
                raise ValueError("invalid differential relation")
        else:
            raise ValueError("chain references an unavailable difference")
        if not 0 < value <= 2 * chain.scalar + 2:
            raise ValueError("chain exceeds its integer workspace bound")
        values.append(value)

    if not 0 <= chain.output < len(values):
        raise ValueError("invalid chain output")
    if values[chain.output] != chain.scalar:
        raise ValueError("chain computes the wrong scalar")
    return True


class _Builder:
    def __init__(self):
        self.instructions = []
        self.multiples = [1]

    def add(self, left, right, difference):
        if len(self.instructions) >= MAX_CHAIN_STEPS:
            raise ValueError("chain instruction limit exceeded")
        a, b = self.multiples[left], self.multiples[right]
        if difference == -1:
            value = 2 * a
        else:
            known = self.multiples[difference]
            total, delta = a + b, abs(a - b)
            if known == total:
                value = delta
            elif known == delta:
                value = total
            else:
                raise ValueError("invalid differential instruction")
        self.instructions.append((left, right, difference))
        self.multiples.append(value)
        return len(self.instructions)

    def double(self, point):
        return self.add(point, point, -1)


def _binary_chain(scalar):
    builder = _Builder()
    if scalar <= 1:
        return Chain(scalar, (), scalar - 1, "binary")
    q, r = 0, builder.double(0)
    for bit in bin(scalar)[3:]:
        total = builder.add(q, r, 0)
        if bit == "1":
            q, r = total, builder.double(r)
        else:
            q, r = builder.double(q), total
    return Chain(scalar, tuple(builder.instructions), q, "binary")


def _prac_chain(scalar, split):
    twos = (scalar & -scalar).bit_length() - 1
    odd = scalar >> twos
    builder = _Builder()
    if odd == 1:
        output = 0
    else:
        d, e = odd - split, 2 * split - odd
        a, b, c = builder.double(0), 0, 0
        # d*a + e*b = odd and c=|a-b| (in integer multiples of P).
        # Each reduction preserves these invariants and gcd(d,e), and
        # strictly decreases d+e. A finite cap also bounds implementation
        # mistakes independently of the mathematical termination argument.
        for _ in range(MAX_CHAIN_STEPS):
            if d == e:
                break
            previous_sum = d + e
            if d < e:
                d, e, a, b = e, d, b, a
            if 4 * d <= 5 * e and (d + e) % 3 == 0:
                d, e = (2 * d - e) // 3, (2 * e - d) // 3
                total = builder.add(a, b, c)
                a, b = builder.add(total, a, b), builder.add(b, total, a)
            elif 4 * d <= 5 * e and (d - e) % 6 == 0:
                d = (d - e) // 2
                b, a = builder.add(a, b, c), builder.double(a)
            elif d <= 4 * e:
                d -= e
                b, c = builder.add(b, a, c), b
            elif (d + e) % 2 == 0:
                d = (d - e) // 2
                b, a = builder.add(b, a, c), builder.double(a)
            elif d % 2 == 0:
                d //= 2
                c, a = builder.add(c, a, b), builder.double(a)
            elif d % 3 == 0:
                d = d // 3 - e
                twice = builder.double(a)
                total = builder.add(a, b, c)
                a, b, c = (
                    builder.add(twice, a, a),
                    builder.add(twice, total, c),
                    b,
                )
            elif (d + e) % 3 == 0:
                d = (d - 2 * e) // 3
                total = builder.add(a, b, c)
                b = builder.add(total, a, b)
                a = builder.add(a, builder.double(a), a)
            elif (d - e) % 3 == 0:
                d = (d - e) // 3
                total = builder.add(a, b, c)
                c, b = builder.add(c, a, b), total
                a = builder.add(a, builder.double(a), a)
            else:
                e //= 2
                c, b = builder.add(c, b, a), builder.double(b)
            if min(d, e) < 1 or d + e >= previous_sum:
                raise ValueError("invalid PRAC reduction")
            a_value, b_value = builder.multiples[a], builder.multiples[b]
            if d * a_value + e * b_value != odd or builder.multiples[c] != abs(
                a_value - b_value
            ):
                raise ValueError("PRAC reduction broke its integer invariant")
        if d != 1 or e != 1:
            raise ValueError("PRAC must terminate at d=e=1")
        output = builder.add(a, b, c)
    for _ in range(twos):
        output = builder.double(output)
    chain = Chain(scalar, tuple(builder.instructions), output, "prac", split)
    verify_chain(chain)
    return chain


def get_chain(scalar, *, add_cost=6, double_cost=5):
    """Get a verified record, or None above 32 bits for ladder fallback.

    At most 30 rational candidate splits, 512 steps per candidate, and 512
    cached records. No search for full-lcm scalars; no points in the cache.
    Costs are positive integer weights so selection is deterministic.
    """
    utils.require_integer(scalar, "scalar", 0)
    utils.require_integer(add_cost, "add_cost", 1)
    utils.require_integer(double_cost, "double_cost", 1)
    if max(add_cost, double_cost).bit_length() > 32:
        raise ValueError("cost weights exceed the cache key limit")
    if scalar.bit_length() > MAX_SCALAR_BITS:
        return None
    return _cached_chain(scalar, add_cost, double_cost)


@lru_cache(maxsize=CACHE_SIZE)
def _cached_chain(scalar, add_cost, double_cost):
    best = _binary_chain(scalar)
    verify_chain(best)
    if scalar < 2:
        return best
    odd = scalar >> ((scalar & -scalar).bit_length() - 1)
    if odd == 1:
        return _prac_chain(scalar, 0)
    candidates = set()
    for numerator in RATIO_NUMERATORS:
        center = (
            odd * numerator + RATIO_DENOMINATOR // 2
        ) // RATIO_DENOMINATOR
        for split in (center - 1, center, center + 1):
            if odd // 2 < split < odd and gcd(odd, split) == 1:
                candidates.add(split)
    for split in sorted(candidates):
        candidate = _prac_chain(scalar, split)
        if candidate.cost(add_cost, double_cost) < best.cost(
            add_cost, double_cost
        ):
            best = candidate
    return best


def clear_cache():
    """Release all immutable chain records (also useful for cold timings)."""
    _cached_chain.cache_clear()


def cache_info():
    """Report the finite record cache's hit/miss/size statistics."""
    return _cached_chain.cache_info()


def _check_point(point, n):
    x, z = point
    divisor = gcd(z, n)
    if divisor == 1:
        return point
    if divisor < n:
        raise NonunitPointError(divisor)
    divisor = gcd(x, n)
    if divisor == 1:
        return point
    if divisor < n:
        raise NonunitPointError(divisor)
    raise _ExceptionalDifferenceError


def _execute(chain, point, n, a24, add, double):
    points = [point]
    try:
        for left, right, difference in chain.instructions:
            if difference == -1:
                result = double(*points[left], n, a24)
            else:
                known = points[difference]
                # Infinity or the order-two point does not disambiguate P+Q
                # from P-Q. Even a nonzero output can be a false infinity.
                if known[0] == 0 or known[1] == 0:
                    raise _ExceptionalDifferenceError
                result = add(*points[left], *points[right], *known, n)
            points.append(_check_point(result, n))
    except _ExceptionalDifferenceError:
        # Preserve factor information before discarding the chain workspace.
        # Every retained Z has already been checked; at most 513 extra GCDs
        # inspect the X coordinates, including the original point.
        for x, _ in points:
            divisor = gcd(x, n)
            if 1 < divisor < n:
                raise NonunitPointError(divisor) from None
        raise
    return (1, 0) if chain.output == -1 else points[chain.output]


def _checked_ladder(scalar, point, n, a24, add, double):
    if scalar == 0 or point[1] == 0:
        return 1, 0
    if point[0] == 0:
        return point if scalar % 2 else (1, 0)
    if scalar == 1:
        return point
    q = point
    r = _check_point(double(*point, n, a24), n)
    for bit in bin(scalar)[3:]:
        total = _check_point(add(*q, *r, *point, n), n)
        if bit == "1":
            q, r = total, _check_point(double(*r, n, a24), n)
        else:
            q, r = _check_point(double(*q, n, a24), n), total
    return q


def multiply(scalar, point, n, a24, add, double, *, chain=None):
    """Execute a verified chain, with at most one checked ladder retry.

    Supplying an external record always re-verifies it. Cached compiler
    records were verified at construction. The curve must be nonsingular,
    with a24=(A+2)/4; callers constructing curves own that precondition.
    """
    utils.require_integer(scalar, "scalar", 0)
    utils.require_integer(n, minimum=2)
    point = point[0] % n, point[1] % n
    if point == (0, 0):
        raise ValueError("(0, 0) is not a projective point")
    _check_point(point, n)
    if chain is not None:
        verify_chain(chain)
        if chain.scalar != scalar:
            raise ValueError("chain scalar does not match request")
    else:
        chain = get_chain(scalar)
    try:
        if chain is not None and point[0] and point[1]:
            try:
                return _execute(chain, point, n, a24, add, double)
            except _ExceptionalDifferenceError:
                pass
        return _checked_ladder(scalar, point, n, a24, add, double)
    except _ExceptionalDifferenceError:
        raise NonunitPointError() from None
