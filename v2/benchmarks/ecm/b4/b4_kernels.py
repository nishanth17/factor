"""Bounded B4 candidates; production retains the readable Montgomery ladder.

All doubling uses a24=(A+2)/4 with the squared difference. These independent
implementations derive from the existing formulas; no upstream code is copied.
See b4_research.md for proofs, exceptional states and intermediate bounds.
"""

from ....common import arithmetic, utils
from . import b4_common as common

ARMS = (
    "baseline",
    "squares",
    "fused_step",
    "whole_ladder",
    "reductions",
    "normalized",
)


class NormalizationNonunitError(ValueError):
    """Retain the exact failed-inversion GCD, including full saturation."""

    def __init__(self, divisor, n):
        super().__init__("fixed-difference denominator is not a unit")
        self.divisor = int(divisor)
        self.factor = int(divisor) if utils.valid_divisor(divisor, n) else None


def normalize_difference(px, pz, n):
    """Return (X/Z, 1) only for unit Z; never discard inversion factors."""
    try:
        inverse = arithmetic.invert(pz, n)
    except arithmetic.NonInvertibleError as error:
        raise NormalizationNonunitError(error.divisor, n) from error
    return px * inverse % n, pz * 0 + 1


def square_add(px, pz, qx, qz, rx, rz, n):
    u = (px - pz) * (qx + qz)
    v = (px + pz) * (qx - qz)
    return rz * (u + v) ** 2 % n, rx * (u - v) ** 2 % n


def square_double(px, pz, n, a24):
    aa, bb = (px + pz) ** 2, (px - pz) ** 2
    delta = aa - bb
    return aa * bb % n, delta * (bb + a24 * delta) % n


def dbladd(px, pz, qx, qz, dx, dz, n, a24):
    """Return (2P, P+Q), sharing P's sum/difference without new assumptions."""
    total, difference = px + pz, px - pz
    aa, bb = total * total, difference * difference
    delta = aa - bb
    u, v = difference * (qx + qz), total * (qx - qz)
    added_total, added_difference = u + v, u - v
    return (
        aa * bb % n,
        delta * (bb + a24 * delta) % n,
        dz * added_total * added_total % n,
        dx * added_difference * added_difference % n,
    )


def _prepare(scalar, px, pz, n):
    utils.require_integer(scalar, "scalar", 0)
    utils.require_integer(n, minimum=2)
    px, pz = px % n, pz % n
    if px == 0 and pz == 0:
        raise ValueError("(0, 0) is not a projective point")
    return px, pz


def multiply_step(scalar, px, pz, n, a24):
    px, pz = _prepare(scalar, px, pz, n)
    if scalar == 0 or pz == 0:
        return px * 0 + 1, pz * 0
    if scalar == 1:
        return px, pz

    qx, qz = px, pz
    rx, rz = square_double(px, pz, n, a24)
    for bit in bin(scalar)[3:]:
        if bit == "1":
            rx, rz, qx, qz = dbladd(rx, rz, qx, qz, px, pz, n, a24)
        else:
            qx, qz, rx, rz = dbladd(qx, qz, rx, rz, px, pz, n, a24)
    return qx, qz


def multiply_whole(scalar, px, pz, n, a24):
    """Inline dbladd in the unchanged adjacent-multiple binary ladder."""
    px, pz = _prepare(scalar, px, pz, n)
    if scalar == 0 or pz == 0:
        return px * 0 + 1, pz * 0
    if scalar == 1:
        return px, pz

    qx, qz = px, pz
    rx, rz = square_double(px, pz, n, a24)
    for bit in bin(scalar)[3:]:
        if bit == "1":
            qx, rx, qz, rz = rx, qx, rz, qz

        total, difference = qx + qz, qx - qz
        aa, bb = total * total, difference * difference
        delta = aa - bb
        u = difference * (rx + rz)
        v = total * (rx - rz)
        added_total, added_difference = u + v, u - v
        qx, qz, rx, rz = (
            aa * bb % n,
            delta * (bb + a24 * delta) % n,
            pz * added_total * added_total % n,
            px * added_difference * added_difference % n,
        )
        if bit == "1":
            qx, rx, qz, rz = rx, qx, rz, qz
    return qx, qz


def multiply_reduced(scalar, px, pz, n, a24):
    """Bound wide temporaries by reducing AA, BB and the two cross products."""
    px, pz = _prepare(scalar, px, pz, n)
    if scalar == 0 or pz == 0:
        return px * 0 + 1, pz * 0
    if scalar == 1:
        return px, pz

    qx, qz = px, pz
    rx, rz = square_double(px, pz, n, a24)
    for bit in bin(scalar)[3:]:
        if bit == "1":
            qx, rx, qz, rz = rx, qx, rz, qz

        total, difference = qx + qz, qx - qz
        aa, bb = total * total % n, difference * difference % n
        delta = aa - bb
        u = difference * (rx + rz) % n
        v = total * (rx - rz) % n
        added_total, added_difference = u + v, u - v
        qx, qz, rx, rz = (
            aa * bb % n,
            delta * (bb + a24 * delta) % n,
            pz * added_total * added_total % n,
            px * added_difference * added_difference % n,
        )
        if bit == "1":
            qx, rx, qz, rz = rx, qx, rz, qz
    return qx, qz


def multiply_normalized(scalar, px, pz, n, a24):
    """Normalize only the fixed difference, once per scalar application."""
    px, pz = _prepare(scalar, px, pz, n)
    if scalar == 0 or pz == 0:
        return px * 0 + 1, pz * 0
    if scalar == 1:
        return px, pz

    dx, _ = normalize_difference(px, pz, n)
    qx, qz = px, pz
    rx, rz = square_double(px, pz, n, a24)
    for bit in bin(scalar)[3:]:
        if bit == "1":
            qx, rx, qz, rz = rx, qx, rz, qz

        total, difference = qx + qz, qx - qz
        aa, bb = total * total, difference * difference
        delta = aa - bb
        u = difference * (rx + rz)
        v = total * (rx - rz)
        added_total, added_difference = u + v, u - v
        qx, qz, rx, rz = (
            aa * bb % n,
            delta * (bb + a24 * delta) % n,
            added_total * added_total % n,
            dx * added_difference * added_difference % n,
        )
        if bit == "1":
            qx, rx, qz, rz = rx, qx, rz, qz
    return qx, qz


def engine(arm):
    """Bind a candidate in a private engine, preserving work ledgers."""
    if arm not in ARMS:
        raise ValueError("unknown B4 arm")
    result = common.load_control("_b4_" + arm)
    package = __import__(result.__package__, fromlist=["ecm", "stage_jobs"])
    ecm = package.ecm
    if arm == "squares":
        ecm.point_add, ecm.point_double = square_add, square_double
    elif arm != "baseline":
        ecm.scalar_multiply = {
            "fused_step": multiply_step,
            "whole_ladder": multiply_whole,
            "reductions": multiply_reduced,
            "normalized": multiply_normalized,
        }[arm]

    if arm == "normalized" and not getattr(result, "_b4_bound", False):
        jobs = package.stage_jobs
        original = jobs.advance_job

        def advance(job, *args, **kwargs):
            try:
                return original(job, *args, **kwargs)
            except NormalizationNonunitError as error:
                # Arithmetic reservations already precede each scalar action.
                # Failed normalization publishes a proper split or curve retry,
                # without committing a partial point or new checkpoint form.
                jobs._finish(job, error.factor)

        jobs.advance_job = advance
        result.advance_job = advance
        result._b4_bound = True
    return result


def ecm_module(engine):
    if engine.__package__ == "v2":
        from ....ecm import core

        return core

    return __import__(engine.__package__ + ".ecm", fromlist=["ecm"])
