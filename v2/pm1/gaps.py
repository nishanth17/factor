"""Finite validation of residue-specific p-1 recurrence tables."""

from ..common import utils
from ..common.arithmetic import pow


def verify_powers(job, budget, mode):
    """Reserve reconstruction of even and exceptional cached powers."""
    powers = job.get("even_powers", [])
    if not isinstance(powers, list) or len(powers) > 64:
        raise ValueError("invalid p-1 recurrence table")
    if powers and (mode != "recurrence" or job["phase"] != "stage_two"):
        # Saturated batches retain the same residue/table during finite replay.
        if mode != "recurrence" or job["phase"] != "term_replay":
            raise ValueError("incompatible p-1 recurrence phase")
    cache = job.get("gap_powers", {})
    if not isinstance(cache, dict) or len(cache) > 64:
        raise ValueError("invalid p-1 gap cache")
    gaps = []
    for key in cache:
        if (
            type(key) is not str
            or not key.isascii()
            or not key.isdecimal()
            or len(key) > len(str(job["b2"]))
        ):
            raise ValueError("invalid p-1 gap key")
        gap = int(key)
        if str(gap) != key or not 1 <= gap <= job["b2"]:
            raise ValueError("invalid p-1 gap key")
        gaps.append(gap)
    if not powers and not gaps:
        return
    budget.consume(
        (len(powers) + 1 if powers else 0)
        + sum(gap.bit_length() + 1 for gap in gaps)
    )
    n, value = job["n"], job["value"]
    square, expected = value * value % n, 1
    for power in powers:
        utils.require_integer(power, "p-1 even power", 0)
        expected = expected * square % n
        if power != expected:
            raise ValueError("corrupt p-1 recurrence table")
    for gap in gaps:
        power = cache[str(gap)]
        utils.require_integer(power, "p-1 gap power", 0)
        if power != pow(value, gap, n):
            raise ValueError("corrupt p-1 gap cache")
