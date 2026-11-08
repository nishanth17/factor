"""Prospective combined challenger; the historical v2 freeze is unchanged.

The primary comparison has two arms. Grouped runtime switches are correctness
and diagnostic controls, not a required all-combinations campaign.
"""

from math import isqrt

from v2.benchmarks.infrastructure.resource_control import (
    adapt_resource_configuration,
)
from v2.benchmarks.infrastructure.resource_control import (
    historical_asdict as asdict,
)
from v2.benchmarks.infrastructure.resource_control import (
    historical_replace as replace,
)

from ....qs import CombinedYieldSieveConfig
from . import plan as historical
from .configuration import SIZES
from .configuration import configuration as baseline

IDENTITY = "c9a-combined-yield-v3"
CURRENT_MAINLINE = "99a2943aa9ba681091184007d6dd8aa7171969a2"


def envelope(bound, incumbent):
    """Preserve the incumbent domain and expose four independent boundaries.

    The proposed ordinary DLP geometry is B <= p <= 4B, p <= q <= 16B,
    p*q <= 32B². It gives both ordered endpoints and product a distinct role.
    Retain old limits where they are larger. The constants are an exploratory
    logarithmic-domain hypothesis, not empirically optimized parameters.
    Reporting remains fixed to isolate admission of already reported values.
    """
    if type(bound) is not int or not 3 <= bound <= 1000000:
        raise ValueError("admission needs a resolved factor-base bound")
    old_product = incumbent.large_product_bound
    old_large = incumbent.large_prime_bound
    large = max(old_large, 16 * bound)
    small = min(large, max(isqrt(old_product), 4 * bound))
    product = max(old_product, 32 * bound**2)
    return dict(
        candidate_bound=incumbent.candidate_bound
        or max(incumbent.residual_bound, old_product),
        residual_bound=max(incumbent.residual_bound, large),
        large_prime_bound=large,
        small_prime_bound=small,
        large_product_bound=min(large * small, product),
    )


def configuration(digits, arm="C1", *, admission_base_bound=None):
    if arm not in ("C0", "C1"):
        raise ValueError("combined primary plan has exactly two arms")
    config = baseline(digits, "C0")
    # Out-of-band routing checks preserve the established routing policy.
    if arm == "C0" or digits not in SIZES:
        return config
    values = asdict(config.collector)
    values.update(envelope(admission_base_bound, config.collector))
    return replace(
        config,
        collector=adapt_resource_configuration(
            CombinedYieldSieveConfig, **values
        ),
    )


def convert(cell):
    result = dict(cell)
    result["configuration"] = asdict(
        configuration(
            cell["digits"],
            cell["arm"],
            admission_base_bound=cell["resolved"]["bound"],
        )
    )
    result["protocol"] = IDENTITY
    return result


def build(corpus, resolved):
    plan = historical.build(corpus, resolved)
    cells = [convert(c) for c in plan["cells"] if c["arm"] in ("C0", "C1")]
    # First stop after four clean representative prefixes. They are retained
    # as the first four A observations if the larger early screen is useful.
    first = [
        c
        for c in cells
        if c["stage"] == "A"
        and c["fixture"] in ("exposed_60_0", "exposed_65_0")
    ]
    cells = first + [c for c in cells if c not in first]
    plan.update(
        identity=IDENTITY,
        cells=cells,
        initial_cells=4,
        early_cells=32,
        current_mainline=CURRENT_MAINLINE,
        ablation_reuse=(
            "none required; targeted diagnosis needs a new decision"
        ),
        feasibility_question=(
            "Does the clean combined arm make useful progress "
            "within matched resources?"
        ),
        parameter_status=(
            "domain-derived starting hypothesis; not a measured optimum"
        ),
    )
    return plan


def confirmation_cells(corpus, resolved, band):
    return [
        convert(c)
        for c in historical.confirmation_cells(corpus, resolved, band)
    ]


def validate(plan):
    counts = {
        stage: sum(c["stage"] == stage for c in plan["cells"])
        for stage in historical.STAGES
    }
    if counts != dict(A=32, B1=4, B2=20, M=2, R=24, F=2, C=0):
        raise ValueError("combined plan cell count changed")
    if any(c["arm"] not in ("C0", "C1") for c in plan["cells"]):
        raise ValueError("unexpected primary arm")
    # Reuse the resource and routing checks without changing the historical
    # identity or fabricating extra measured observations.
    for stage, values in historical.STAGES.items():
        cells = [c for c in plan["cells"] if c["stage"] == stage]
        if sum(c["warm_seconds"] for c in cells) > values["warm"]:
            raise ValueError("combined warmup reservation exceeded")
        if (
            sum(c["seconds"] + c["margin"] + c["warm_seconds"] for c in cells)
            > values["seconds"]
        ):
            raise ValueError("combined lifecycle reservation exceeded")
    for cell in plan["cells"]:
        if cell["digits"] == 80 and cell["kind"] != "prefix":
            raise ValueError("complete 80-digit cell forbidden")
        if (
            cell["configuration"]["memory_bytes"] + cell["helper_bytes"]
            != cell["owned_bytes"]
        ):
            raise ValueError("combined owned memory is unmatched")
        if cell["owned_bytes"] > 4 * 2**30:
            raise ValueError("native owned ceiling exceeded")
    return counts
