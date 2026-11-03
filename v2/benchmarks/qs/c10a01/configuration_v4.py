"""Literal v4 bundle and explicit component ablations; no timing admission."""

from dataclasses import replace

from .configuration_v3 import configuration as previous_configuration


def configuration(
    digits,
    arm,
    *,
    multiplier=0,
    root_update_policy="batch",
    assignment_policy="windowed_pair",
    multiplier_policy="parity2_norm",
):
    """Preserve matched capacities; select all four new components."""
    if arm not in ("control", "coupled"):
        raise ValueError("unknown v4 arm")
    config = previous_configuration(digits, arm, multiplier=multiplier)
    if arm == "control":
        return config
    return replace(
        config,
        root_update_policy=root_update_policy,
        assignment_policy=assignment_policy,
        multiplier_policy=multiplier_policy,
    )
