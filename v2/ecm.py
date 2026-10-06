"""Two-stage Montgomery ECM with modular Suyama setup and exact limits.

The readable ladder remains the production baseline. The experimental
multiply_prac entry point uses bounded, verified records and checked recovery.
"""

from dataclasses import dataclass
from math import gcd, isqrt, prod
from typing import Optional, Tuple

from . import constants, prac, prime_sieve, utils

Point = Tuple[int, int]


@dataclass(frozen=True)
class CurveSetup:
    """Exactly one of a usable curve, a proper factor, or an explicit retry."""

    point: Optional[Point] = None
    a24: Optional[int] = None
    factor: Optional[int] = None
    retry: bool = False


@dataclass
class EcmStats:
    """Curve attempts and transitions, including unusable setup attempts."""

    curves: int = 0
    setup_retries: int = 0
    stage_one_calls: int = 0
    stage_one_saturations: int = 0
    stage_two_calls: int = 0
    stage_two_saturations: int = 0


def compute_bounds(n):
    """Use modest defaults until Phase 2 calibrates factor-size tiers."""
    utils.require_integer(n, minimum=1)
    return constants.ECM_B1, constants.ECM_B2


def setup_curve(n, sigma):
    """Construct Suyama X:Z coordinates without integer/float division.

    Modular division is legal only for a unit. A nonunit setup denominator
    or singularity discriminant may itself reveal an unknown factor.
    """
    utils.require_integer(n, minimum=3)
    utils.require_integer(sigma, "sigma", 6)
    if n % 2 == 0:
        return CurveSetup(factor=2)

    u = (sigma * sigma - 5) % n
    v = 4 * sigma % n
    u_cubed = pow(u, 3, n)
    v_cubed = pow(v, 3, n)
    denominator = 16 * u_cubed * v % n
    divisor = gcd(denominator, n)

    if utils.valid_divisor(divisor, n):
        return CurveSetup(factor=divisor)
    if divisor == n:
        return CurveSetup(retry=True)

    inverse = utils.modular_inverse(denominator, n)
    a24 = pow(v - u, 3, n) * (3 * u + v) * inverse % n
    curve_a = (4 * a24 - 2) % n
    divisor = gcd(curve_a * curve_a - 4, n)

    if utils.valid_divisor(divisor, n):
        return CurveSetup(factor=divisor)
    if divisor == n:
        return CurveSetup(retry=True)

    return CurveSetup(point=(u_cubed, v_cubed), a24=a24)


def point_add(px, pz, qx, qz, rx, rz, n):
    """Differential addition: R must equal P-Q up to sign.

    X-only coordinates do not support arbitrary point addition. An output
    (0, 0) is degenerate; factoring callers must GCD-check and retry/recover.
    """
    u = (px - pz) * (qx + qz)
    v = (px + pz) * (qx - qz)
    total, difference = u + v, u - v
    x = rz * total * total % n
    z = rx * difference * difference % n
    return x, z


def point_double(px, pz, n, a24):
    """Double X:Z on y²=x³+A*x²+x with a24=(A+2)/4 modulo n."""
    total, difference = px + pz, px - pz
    sum_squared = total * total
    difference_squared = difference * difference
    delta = sum_squared - difference_squared
    x = sum_squared * difference_squared % n
    z = delta * (difference_squared + a24 * delta) % n
    return x, z


def scalar_multiply(scalar, px, pz, n, a24):
    """Montgomery ladder for nonnegative scalars; 0 returns infinity (1:0)."""
    utils.require_integer(scalar, "scalar", 0)
    utils.require_integer(n, minimum=2)
    px, pz = px % n, pz % n
    if px == 0 and pz == 0:
        raise ValueError("(0, 0) is not a projective point")
    if scalar == 0 or pz == 0:
        return 1, 0
    if scalar == 1:
        return px, pz

    # Q and R remain adjacent multiples, so their difference is always P.
    qx, qz = px, pz
    rx, rz = point_double(px, pz, n, a24)

    for bit in bin(scalar)[3:]:
        if bit == "1":
            qx, qz = point_add(rx, rz, qx, qz, px, pz, n)
            rx, rz = point_double(rx, rz, n, a24)
        else:
            rx, rz = point_add(qx, qz, rx, rz, px, pz, n)
            qx, qz = point_double(qx, qz, n, a24)

    return qx, qz


def multiply_prac(scalar, px, pz, n, a24):
    """Experimental verified PRAC; return X:Z or raise prac.NonunitPointError.

    The exception carries a proper factor, or factor=None for curve retry.
    Large scalars and exceptional differences use a checked ladder. This
    entry point is opt-in; factoring jobs retain scalar_multiply until B3.
    """
    return prac.multiply(scalar, (px, pz), n, a24, point_add, point_double)


def stage_one_scalar(b1):
    """Exact lcm(1, ..., B1); moderate default bounds keep it manageable."""
    utils.require_integer(b1, "b1", 2)
    return prod(
        utils.prime_power(prime, b1)
        for prime in prime_sieve.prime_sieve(b1 + 1)
    )


def _stage_two_terms(point, n, a24, b1, primes):
    """Yield the original even-baby-step cross differences for every q.

    Comparing [r]Q and [2d]Q covers q=r+2d because X(P)=X(-P).
    Future prime pairing can exploit the other sign. Here every eligible
    prime is traversed once, and initialization never uses a negative scalar.
    """
    px, pz = point
    center = b1 if b1 % 2 else b1 - 1
    distance = min(isqrt(primes[-1]), (center - 1) // 2)

    if distance < 2:
        # Tiny B1 has no legal positive preceding giant-step scalar.
        for prime in primes:
            yield scalar_multiply(prime, px, pz, n, a24)[1]
        return

    baby_steps = [None, point_double(px, pz, n, a24)]
    baby_steps.append(point_double(*baby_steps[1], n, a24))
    for index in range(3, distance + 1):
        baby_steps.append(
            point_add(
                *baby_steps[index - 1],
                *baby_steps[1],
                *baby_steps[index - 2],
                n,
            )
        )

    step = 2 * distance
    giant = scalar_multiply(center, px, pz, n, a24)
    previous = scalar_multiply(center - step, px, pz, n, a24)

    for prime in primes:
        while prime > center + step:
            previous, giant = (
                giant,
                point_add(*giant, *baby_steps[distance], *previous, n),
            )
            center += step

        index = (prime - center) // 2
        bx, bz = baby_steps[index]
        # Cross multiplication avoids affine conversion or another inversion.
        yield (giant[0] * bz - bx * giant[1]) % n


def stage_two(point, n, a24, b1, primes, batch_size=constants.GCD_BATCH_SIZE):
    """Return (factor_or_None, saturation_seen), replaying local batches."""
    if not primes:
        return None, False
    terms = []

    for term in _stage_two_terms(point, n, a24, b1, primes):
        terms.append(term)
        if len(terms) == batch_size:
            divisor, saturated = utils.batch_factor(terms, n)
            if divisor is not None or saturated:
                return divisor, saturated
            terms.clear()

    # The trailing partial batch still belongs to the inclusive B2 schedule.
    return utils.batch_factor(terms, n) if terms else (None, False)


def factorize_ecm(
    n,
    verbose=False,
    *,
    b1=None,
    b2=None,
    seed=None,
    rng=None,
    max_curves=constants.MAX_CURVES_ECM,
    batch_size=constants.GCD_BATCH_SIZE,
    stats=None,
    _known_composite=False,
):
    """Return a proper divisor or None after exactly at most max_curves.

    An unusable curve counts as an attempt. Fully saturated stage 1 retries
    immediately; stage 2 can recover a mixed-factor batch locally. The internal
    _known_composite hint avoids repeating the dispatcher's classification.
    """
    utils.require_integer(n, minimum=1)
    default_b1, default_b2 = compute_bounds(n)
    b1 = default_b1 if b1 is None else b1
    b2 = default_b2 if b2 is None else b2
    utils.require_integer(b1, "b1", 2)
    utils.require_integer(b2, "b2", b1)
    utils.require_integer(max_curves, "max_curves", 0)
    utils.require_integer(batch_size, "batch_size", 1)
    generator = utils.resolve_rng(seed, rng)
    if n == 1 or not max_curves:
        return None
    if n % 2 == 0:
        return 2 if n > 2 else None
    if not _known_composite and utils.is_prime(n, rng=generator):
        return None

    scalar = stage_one_scalar(b1)
    work = stats if stats is not None else EcmStats()
    stage_two_primes = None

    for _ in range(max_curves):
        work.curves += 1
        sigma = generator.randint(6, constants.MAX_RANDOM_ECM)
        setup = setup_curve(n, sigma)
        if setup.factor is not None:
            return setup.factor
        if setup.retry:
            work.setup_retries += 1
            continue
        if verbose:
            print(f"ECM curve {work.curves}, sigma={sigma}, bounds={b1}/{b2}")

        work.stage_one_calls += 1
        point = scalar_multiply(scalar, *setup.point, n, setup.a24)
        divisor = gcd(point[1], n)
        if utils.valid_divisor(divisor, n):
            return divisor
        if divisor == n:
            work.stage_one_saturations += 1
            continue

        # Stage-1 success or saturation never pays for the B2 prime list.
        if stage_two_primes is None:
            stage_two_primes = prime_sieve.segmented_sieve(b1 + 1, b2 + 1)
        work.stage_two_calls += 1
        divisor, saturated = stage_two(
            point, n, setup.a24, b1, stage_two_primes, batch_size
        )
        if saturated:
            work.stage_two_saturations += 1
        if utils.valid_divisor(divisor, n):
            return divisor

    return None
