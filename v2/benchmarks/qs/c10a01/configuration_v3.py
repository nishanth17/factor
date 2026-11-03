"""Literal combined-arm configuration; no timed protocol is admitted here."""

from dataclasses import replace

from .configuration import configuration as previous_configuration


def configuration(digits, arm, *, multiplier=0):
    """Keep matched grants/capacities; opt into the four coupled changes."""
    if arm not in ("control", "coupled"):
        raise ValueError("unknown v3 arm")
    config = previous_configuration(
        digits,
        "control" if arm == "control" else "parity2",
        multiplier=multiplier,
    )
    if arm == "control":
        return config
    return replace(
        config,
        multiplier_policy="parity2_wide",
        assignment_policy="paired_flyer",
        collector=replace(config.collector, marking="parity_two"),
    )
