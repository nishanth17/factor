"""Bounded precomputed chains for C6; production routing belongs to B3.

SSA records are compiled into four-byte register instructions. A separate
integer machine verifies those instructions, including overwritten slots.
Arithmetic uses the pinned, readable Montgomery kernels. No search occurs
in the executor. See c6_research.md for upstream provenance and scope.
"""

import json
from dataclasses import dataclass
from pathlib import Path

from .. import prac, prime_sieve, utils

MAX_STEPS = 512
MAX_POINTS = 16
MAX_RECORDS = 512
MAX_FILE_BYTES = 1024 * 1024
MAX_BOUND = 2000
DATA = Path(__file__).parent / "inputs/controls/c6_lucas_records.json"
DOUBLE = 255


@dataclass(frozen=True)
class Record:
    """Immutable register bytecode; slot zero starts with P."""

    scalar: int
    code: bytes
    output: int
    slots: int


def verify(record):
    """Prove scalar action without trusting SSA or generator metadata."""
    if type(record) is not Record:
        raise ValueError("expected a compact Record")
    if (
        type(record.scalar) is not int
        or not 1 <= record.scalar < 2**32
        or type(record.code) is not bytes
        or len(record.code) % 4
        or len(record.code) > 4 * MAX_STEPS
        or type(record.slots) is not int
        or not 1 <= record.slots <= MAX_POINTS
        or type(record.output) is not int
        or not 0 <= record.output < record.slots
    ):
        raise ValueError("compact record exceeds structural limits")

    values = [None] * record.slots
    values[0] = 1
    for offset in range(0, len(record.code), 4):
        stop = offset + 4
        dest, left, right, difference = record.code[offset:stop]
        if max(dest, left, right) >= record.slots:
            raise ValueError("register out of bounds")
        a, b = values[left], values[right]
        if a is None or b is None:
            raise ValueError("uninitialized register")
        if difference == DOUBLE:
            if left != right:
                raise ValueError("invalid doubling")
            value = 2 * a
        else:
            if difference >= record.slots or values[difference] is None:
                raise ValueError("uninitialized difference")
            known = values[difference]
            if known == abs(a - b):
                value = a + b
            elif known == a + b:
                value = abs(a - b)
            else:
                raise ValueError("invalid differential identity")
        if not 0 < value <= 2 * record.scalar + 2:
            raise ValueError("integer workspace exceeded")
        values[dest] = value
    if values[record.output] != record.scalar:
        raise ValueError("wrong scalar action")
    return True


def compact(chain, *, rolling=False):
    """Allocate last-use registers, or the upstream 16-point rolling ring."""
    prac.verify_chain(chain)
    if chain.scalar < 1:
        raise ValueError("zero uses the finite ladder fallback")
    uses = [0] * (len(chain.instructions) + 1)
    uses[chain.output] += 1
    for instruction in chain.instructions:
        for node in instruction:
            if node >= 0:
                uses[node] += 1
    slots, free, code = {0: 0}, [], bytearray()
    next_slot = 1
    for index, (left, right, difference) in enumerate(chain.instructions, 1):
        operands = (left, right, difference)
        if rolling and any(
            index - node > 16 for node in operands if node >= 0
        ):
            raise ValueError("record does not fit the rolling ring")
        encoded = [slots[node] if node >= 0 else DOUBLE for node in operands]
        for node in operands:
            if node >= 0:
                uses[node] -= 1
                if uses[node] == 0:
                    free.append(slots.pop(node))
        if rolling:
            dest = index % 16
            next_slot = min(16, index + 1)
        elif free:
            dest = free.pop()
        else:
            dest, next_slot = next_slot, next_slot + 1
        if next_slot > MAX_POINTS:
            raise ValueError("live-point limit exceeded")
        slots[index] = dest
        code.extend((dest, *encoded))
    result = Record(chain.scalar, bytes(code), slots[chain.output], next_slot)
    verify(result)
    return result


def decode_rows(rows):
    """Translate upstream offsets; integer proofs ignore its value column."""
    chains = {}
    if len(rows) > MAX_RECORDS:
        raise ValueError("too many upstream records")
    for row in rows:
        scalar, elements = row["prime"], row["elements"]
        if scalar in chains or len(elements) > 64:
            raise ValueError("duplicate or oversized upstream record")
        instructions = []
        values = [1]
        for index, (value, left, right, difference) in enumerate(elements, 1):
            parent = index - 1
            if any(
                type(v) is not int for v in (value, left, right, difference)
            ) or not (
                0 <= left <= parent
                and 0 <= right <= parent
                and 0 <= difference <= parent
            ):
                raise ValueError("invalid upstream offsets")
            operands = (parent - left, parent - right)
            known = -1 if difference == 0 else parent - difference
            instructions.append((*operands, known))
            # Check the upstream value independently, not only the endpoint.
            a, b = (values[node] for node in operands)
            if known == -1:
                actual = 2 * a
            else:
                if values[known] != abs(a - b):
                    raise ValueError("upstream differential mismatch")
                actual = a + b
            if actual != value:
                raise ValueError("upstream integer value mismatch")
            values.append(actual)
        chain = prac.Chain(scalar, tuple(instructions), len(elements), "lucas")
        prac.verify_chain(chain)
        verify(compact(chain))
        chains[scalar] = chain
    return chains


def load_lucas():
    """Read a finite, pinned upstream-decoded catalog; no persistent cache."""
    with DATA.open("rb") as source:
        raw = source.read(MAX_FILE_BYTES + 1)
    if len(raw) > MAX_FILE_BYTES:
        raise ValueError("catalog exceeds byte cap")
    data = json.loads(raw)
    chains = decode_rows(data["records"])
    if set(chains) != set(prime_sieve.prime_sieve(MAX_BOUND + 1)):
        raise ValueError("catalog does not cover the declared bound")
    return chains


def compose(chain, power):
    """Repeat a prime action to obtain its power, without full-lcm search."""
    remaining, count = power, 0
    while remaining > 1:
        remaining, remainder = divmod(remaining, chain.scalar)
        if remainder or count >= 31:
            raise ValueError("not a bounded power of this prime")
        count += 1
    instructions, output = [], 0
    for _ in range(count):
        mapping = [output]
        for left, right, difference in chain.instructions:
            instructions.append(
                (
                    mapping[left],
                    mapping[right],
                    mapping[difference] if difference >= 0 else -1,
                )
            )
            mapping.append(len(instructions))
        output = mapping[chain.output]
    result = prac.Chain(power, tuple(instructions), output, "lucas")
    prac.verify_chain(result)
    return result


class Executor:
    """Own a verified immutable record and its fixed arithmetic backend.

    Every X and Z is checked before a point can be overwritten. This extra
    X check preserves factors a rolling buffer would otherwise discard.
    Infinity/order-two differences receive at most one checked ladder retry.
    """

    def __init__(self, record, backend):
        verify(record)
        self.record = record
        self.backend = backend

    def __call__(self, point, n, a24):
        record, backend = self.record, self.backend
        check = backend.prac._check_point
        point = point[0] % n, point[1] % n
        if point == (0, 0):
            raise ValueError("invalid projective point")
        check(point, n)
        add, double, gcd = (
            backend.ecm.point_add,
            backend.ecm.point_double,
            backend.gcd,
        )
        divisor = gcd(point[0], n)
        if 1 < divisor < n:
            raise prac.NonunitPointError(divisor)
        try:
            if point[0] and point[1]:
                points = [None] * record.slots
                points[0] = point
                code = record.code
                for offset in range(0, len(code), 4):
                    stop = offset + 4
                    dest, left, right, difference = code[offset:stop]
                    if difference == DOUBLE:
                        result = double(*points[left], n, a24)
                    else:
                        known = points[difference]
                        if known[0] == 0 or known[1] == 0:
                            raise prac._ExceptionalDifferenceError
                        result = add(*points[left], *points[right], *known, n)
                    result = check(result, n)
                    divisor = gcd(result[0], n)
                    if 1 < divisor < n:
                        raise prac.NonunitPointError(divisor)
                    points[dest] = result
                return points[record.output]
        except prac._ExceptionalDifferenceError:
            pass
        try:
            return checked_ladder(record.scalar, point, n, a24, backend)
        except prac._ExceptionalDifferenceError:
            raise prac.NonunitPointError() from None


def build_program(bound, method, backend):
    """Construct and verify a bound-owned program, with <=512 records."""
    utils.require_integer(bound, "bound", 2)
    if bound > MAX_BOUND:
        raise ValueError("C6 supports B1 <= 2000")
    if method not in ("checked", "compact", "lucas", "rolling"):
        raise ValueError("unknown C6 method")
    lucas = load_lucas() if method in ("lucas", "rolling") else None
    program, retained = [], {}
    prac.clear_cache()
    for prime in prime_sieve.prime_sieve(bound + 1):
        power = utils.prime_power(prime, bound)
        actions = []
        for scalar in (power, prime):
            if scalar not in retained:
                if len(retained) >= MAX_RECORDS:
                    raise ValueError("program record cap exceeded")
                chain = (
                    compose(lucas[prime], scalar)
                    if lucas
                    else prac.get_chain(scalar)
                )
                retained[scalar] = (
                    chain
                    if method == "checked"
                    else Executor(
                        compact(chain, rolling=method == "rolling"), backend
                    )
                )
            actions.append(retained[scalar])
        program.append((prime, power, *actions))
    prac.clear_cache()
    return tuple(program)


def apply(action, scalar, point, n, a24, backend):
    if isinstance(action, Executor):
        return action(point, n, a24)
    return backend.prac.multiply(
        scalar,
        point,
        n,
        a24,
        backend.ecm.point_add,
        backend.ecm.point_double,
        chain=action,
    )


def continued_fraction_bits(chain):
    """Recognize the paper's restricted family without any chain search."""
    prac.verify_chain(chain)
    values = [1]
    for left, right, difference in chain.instructions:
        if difference == -1:
            # Algorithm 1 has only its initial doubling. The same integer
            # sequence can also use extra doublings; exclude those records
            # so this diagnostic changes dispatch, not arithmetic cost.
            if len(values) != 1:
                return None
            values.append(2 * values[left])
        elif values[difference] == abs(values[left] - values[right]):
            values.append(values[left] + values[right])
        else:
            return None
    if values[:3] != [1, 2, 3] or chain.output != len(values) - 1:
        return None
    a, b, c, bits = 1, 2, 3, []
    for value in values[3:]:
        if value == b + c:
            bits.append(0)
            a, b, c = b, c, value
        elif value == a + c:
            bits.append(1)
            a, b, c = a, c, value
        else:
            return None
    return bytes(bits)


class ThreePointExecutor(Executor):
    """Algorithm 1's three live points, for recognized CF records only."""

    def __init__(self, chain, backend):
        bits = continued_fraction_bits(chain)
        if bits is None:
            raise ValueError("chain is outside the continued-fraction family")
        super().__init__(compact(chain), backend)
        self.bits = bits

    def __call__(self, point, n, a24):
        backend = self.backend
        add, double = backend.ecm.point_add, backend.ecm.point_double

        def check(point):
            point = backend.prac._check_point(point, n)
            divisor = backend.gcd(point[0], n)
            if 1 < divisor < n:
                raise prac.NonunitPointError(divisor)
            return point

        def addition(left, right, known):
            if not known[0] or not known[1]:
                raise prac._ExceptionalDifferenceError
            return check(add(*left, *right, *known, n))

        point = point[0] % n, point[1] % n
        if point == (0, 0):
            raise ValueError("invalid projective point")
        check(point)
        try:
            if point[0] and point[1]:
                a = point
                b = check(double(*point, n, a24))
                c = addition(a, b, a)
                for bit in self.bits:
                    if bit:
                        a, b, c = a, c, addition(a, c, b)
                    else:
                        a, b, c = b, c, addition(b, c, a)
                return c
        except prac._ExceptionalDifferenceError:
            pass
        try:
            return checked_ladder(self.record.scalar, point, n, a24, backend)
        except prac._ExceptionalDifferenceError:
            raise prac.NonunitPointError() from None


def checked_ladder(scalar, point, n, a24, backend):
    """One finite retry, retaining X nonunits before any point is discarded."""

    def check(result):
        result = backend.prac._check_point(result, n)
        divisor = backend.gcd(result[0], n)
        if 1 < divisor < n:
            raise prac.NonunitPointError(divisor)
        return result

    if scalar == 0 or point[1] == 0:
        return 1, 0
    if point[0] == 0:
        return point if scalar % 2 else (1, 0)
    if scalar == 1:
        return point
    double, add = backend.ecm.point_double, backend.ecm.point_add
    q, r = point, check(double(*point, n, a24))
    for bit in bin(scalar)[3:]:
        total = check(add(*q, *r, *point, n))
        if bit == "1":
            q, r = total, check(double(*r, n, a24))
        else:
            q, r = check(double(*q, n, a24)), total
    return q
