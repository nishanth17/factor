"""C6 scalar/coordinate certificates and strict finite recovery, reused by B3.

These interpreters retain C6's arithmetic and independent verification. No
chain generator, allocator or search is part of production preparation.
"""

from dataclasses import dataclass

from . import prac

MAX_STEPS = 512
MAX_POINTS = 16
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


def instructions(record):
    code = record.code
    result = []
    for start in range(0, len(code), 4):
        stop = start + 4
        result.append(tuple(code[start:stop]))
    return tuple(result)


def verify_frontier(record, masks):
    """Independently prove coverage using forward ancestor bitsets.

    Guards together with the two output coordinates must cover every
    intermediate coordinate, including discarded points and the input.
    Checking the next record therefore also certifies this record's output.
    """
    verify(record)
    count = len(record.code) // 4
    if type(masks) is not bytes or len(masks) != count + 1:
        raise ValueError("wrong coverage certificate size")
    if any(mask > 3 for mask in masks):
        raise ValueError("invalid coordinate mask")
    states = [None] * record.slots
    states[0] = (1, 2)
    covered = (1 if masks[0] & 1 else 0) | (2 if masks[0] & 2 else 0)
    for index, (dest, left, right, difference) in enumerate(
        instructions(record), 1
    ):
        x, z = 1 << (2 * index), 1 << (2 * index + 1)
        if difference == DOUBLE:
            previous_x, previous_z = states[left]
            z |= previous_x | previous_z
        else:
            known_x, known_z = states[difference]
            x |= known_z
            z |= known_x
        if masks[index] & 1:
            covered |= x
        if masks[index] & 2:
            covered |= z
        states[dest] = x, z
    output_x, output_z = states[record.output]
    covered |= output_x | output_z
    if covered != (1 << (2 * (count + 1))) - 1:
        raise ValueError("certificate loses an intermediate factor")
    return True


def _factor(value, point, mask):
    if mask & 1:
        value *= point[0]
    if mask & 2:
        value *= point[1]
    return value
