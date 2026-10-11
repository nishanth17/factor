"""Size-scoped ECM defaults promoted from C3's revealed-input screen."""

from dataclasses import replace

from ..common import utils
from .ecm_presets import ECM_PRESETS

COMPACT64 = ((2000, 50000, 64),)
ESCALATING140 = ECM_PRESETS["yamaquasi_ecm64"]


def resolve_defaults(config, n):
    """Choose measured size bands without expanding any resource allowance.

    The human requested promotion before fresh confirmation. Other sizes and
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
        elif 10**39 <= magnitude < 10**40:
            tiers = ESCALATING140
    try:
        return replace(config, ecm_tiers=tiers)
    except MemoryError:
        # Admission remains conservative even when ECM and SIQS coexist.
        # Freeze the old plan too, so children cannot trigger a new choice.
        return replace(config, ecm_tiers=config.ecm_tiers)
