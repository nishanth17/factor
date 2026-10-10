"""Finite validation of residue-specific p-1 recurrence tables."""

from . import utils


def verify_powers(job, budget, mode):
    """Reserve reconstruction before using any deserialized even powers."""
    powers = job.get("even_powers", [])
    if not isinstance(powers, list) or len(powers) > 64:
        raise ValueError("invalid p-1 recurrence table")
    if powers and (mode != "recurrence" or job["phase"] != "stage_two"):
        # Saturated batches retain the same residue/table during finite replay.
        if mode != "recurrence" or job["phase"] != "term_replay":
            raise ValueError("incompatible p-1 recurrence phase")
    if not powers:
        return
    budget.consume(len(powers) + 1)
    n, value = job["n"], job["value"]
    square, expected = value * value % n, 1
    for power in powers:
        utils.require_integer(power, "p-1 even power", 0)
        expected = expected * square % n
        if power != expected:
            raise ValueError("corrupt p-1 recurrence table")
