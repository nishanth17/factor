"""Size-scoped numerical plans adopted by C3 and automatic ECM selection."""

from dataclasses import replace

from ..common import utils
from .ecm_presets import ECM_PRESETS

COMPACT64 = ((2000, 50000, 64),)
ESCALATING140 = ECM_PRESETS["escalating140"]


def resolve_defaults(config, n):
    """Choose measured size bands without expanding any resource allowance.

    C3 selected 30/40 digits; automatic ECM adds 60/80/100. Other sizes and
    explicitly selected paired executors retain the previous numerical plan.
    A tight existing memory cap also retains that plan rather than reducing
    the relation engine's storage or enlarging the parent's allowance.
    """
    utils.require_integer(n)
    magnitude = abs(n)
    tiers = config.ecm_tiers
    if config.ecm_pair_distance is None and config.ecm_pair_wheel is None:
        if 10**29 <= magnitude < 10**30:
            tiers = COMPACT64
        elif any(
            10 ** (digits - 1) <= magnitude < 10**digits
            for digits in (40, 60, 80, 100)
        ):
            tiers = ESCALATING140
    try:
        return replace(config, ecm_tiers=tiers)
    except MemoryError:
        # Admission remains conservative even when ECM and SIQS coexist.
        # Freeze the old plan too, so children cannot trigger a new choice.
        return replace(config, ecm_tiers=config.ecm_tiers)
