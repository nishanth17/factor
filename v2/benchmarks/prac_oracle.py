"""Independent full-coordinate controls for the P4.1 chain experiment."""

from math import gcd


def affine_add(point, other, modulus, curve_a, curve_b=1):
    """Group law on B*y^2=x^3+A*x^2+x; None denotes infinity."""
    if point is None:
        return other
    if other is None:
        return point
    x, y = point
    u, v = other
    if x == u and (y + v) % modulus == 0:
        return None
    if point == other:
        numerator = 3 * x * x + 2 * curve_a * x + 1
        denominator = 2 * curve_b * y
    else:
        numerator, denominator = v - y, u - x
    slope = numerator * pow(denominator, -1, modulus) % modulus
    result_x = (curve_b * slope * slope - curve_a - x - u) % modulus
    return result_x, (slope * (x - result_x) - y) % modulus


def affine_multiply(scalar, point, modulus, curve_a, curve_b=1):
    """Right-to-left full-coordinate multiplication, independent of PRAC."""
    result = None
    while scalar:
        if scalar & 1:
            result = affine_add(result, point, modulus, curve_a, curve_b)
        scalar >>= 1
        if scalar:
            point = affine_add(point, point, modulus, curve_a, curve_b)
    return result


def matches(point, expected, modulus):
    """Reject invalid projective pairs before testing an affine equality."""
    x, z = point
    if gcd(gcd(x, z), modulus) != 1:
        return False
    if expected is None:
        return z % modulus == 0
    return gcd(z, modulus) == 1 and (x - expected[0] * z) % modulus == 0


def historical_points():
    """The original audit's first four points on each of four small curves."""
    for modulus, curve_a in ((1009, 6), (1013, 11), (10007, 6), (10009, 17)):
        roots = {}
        for y in range(1, modulus):
            roots.setdefault(y * y % modulus, y)
        count = 0
        for x in range(2, modulus):
            square = (x**3 + curve_a * x * x + x) % modulus
            if square in roots:
                yield modulus, curve_a, (x, roots[square])
                count += 1
                if count == 4:
                    break


def twist_point(point, a24, modulus):
    """Lift a unit Suyama X:Z point to (x,1) on its independent twist."""
    x = point[0] * pow(point[1], -1, modulus) % modulus
    curve_a = (4 * a24 - 2) % modulus
    curve_b = (x**3 + curve_a * x * x + x) % modulus
    if gcd(curve_b, modulus) != 1:
        raise ValueError("the point has no nonsingular y=1 twist")
    return (x, 1), curve_a, curve_b
