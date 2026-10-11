"""Explicit source-derived ECM schedules with protected SIQS admission."""

from dataclasses import replace
from types import MappingProxyType

from ..portfolio import PortfolioConfig

# Numerical schedules are facts transcribed from the pinned upstream review;
# no upstream arithmetic or scheduling code is copied. Native curve families
# and stage-two implementations differ, so these are explicit transfer presets,
# not claims of equal success probabilities or accepted PyPy speedups.
ECM_PRESETS = MappingProxyType(
    {
        "yamaquasi_auto160": ((200, 7700, 8),),
        "yamaquasi_ecm64": (
            (200, 7700, 10),
            (2000, 81000, 30),
            (10000, 554000, 100),
        ),
        "alpertron20": ((2000, 200000, 25), (11000, 1100000, 90)),
        "alpertron25": (
            (2000, 200000, 25),
            (11000, 1100000, 90),
            (50000, 5000000, 300),
        ),
        "gmp_ecm20": ((11000, 1900000, 74),),
        "gmp_ecm25": ((50000, 13000000, 214),),
    }
)


def with_ecm_preset(config, name):
    """Replace only ECM tiers, retaining explicit SIQS and admission bounds.

    The caller supplies the existing cumulative allocation with positive
    work/wall/CPU fallback floors. The normal PortfolioConfig constructor
    checks coexistence of the larger point/schedule workspace with unchanged
    SIQS storage; an unfunded preset raises rather than shrinking SIQS.
    Pretest ceilings, total run budgets and checkpoint behavior stay binding.
    """
    if not isinstance(config, PortfolioConfig):
        raise TypeError("ECM presets require a PortfolioConfig")
    if name not in ECM_PRESETS:
        raise ValueError("unknown source-derived ECM preset: " + str(name))
    policy = config.allocation
    if config.siqs is None or policy is None:
        raise ValueError("ECM presets require SIQS and protected allocation")
    if (
        min(
            policy.fallback_work,
            policy.fallback_seconds,
            policy.fallback_cpu_seconds,
        )
        <= 0
    ):
        raise ValueError(
            "ECM presets require positive SIQS work/wall/CPU floors"
        )

    return replace(config, ecm_tiers=ECM_PRESETS[name])
