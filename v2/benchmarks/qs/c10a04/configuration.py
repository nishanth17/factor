"""Exact, charged translation from eligible base counts to prime endpoints."""

from math import gcd

from v2.benchmarks.infrastructure.resource_control import (
    adapt_resource_configuration,
)
from v2.benchmarks.infrastructure.resource_control import (
    historical_asdict as asdict,
)
from v2.benchmarks.infrastructure.resource_control import (
    historical_replace as replace,
)

from ....common import prime_sieve, utils
from ....qs.defaults import default_siqs_config
from ....qs.factor_base import checked_target
from ....qs.multiplier import select_multiplier

ARMS = ("control", "compact", "wide")
BANDS = (55, 60, 65)
COUNTS = {55: 2500, 60: 5000, 65: 8500}
FACTORS = {55: 7, 60: 8, 65: 9}
WIDTHS = {55: 131072, 60: 131072, 65: 262144}
ENDPOINT_CEILING = 250000
HELPER_BYTES = 8192 + 64 * ENDPOINT_CEILING


class ConfigurationRefusedError(ValueError):
    """A finite count policy could not reach its requested configuration."""


def template(digits, arm):
    """Use the accepted 50/51-digit collector in every configuration arm."""
    if digits not in BANDS or arm not in ARMS:
        raise ValueError("unfrozen band or arm")
    control = adapt_resource_configuration(default_siqs_config, 10**49)
    # Preserve this historical collector when automatic marking advances.
    control = replace(
        control,
        multiplier=1,
        assignment_policy="flyer",
        polynomial_policy="standard",
        multiplier_policy="standard",
        root_update_policy="standard",
        collector=replace(control.collector, marking="sparse"),
    )
    if arm == "control":
        return control

    width = 65536 if arm == "compact" else WIDTHS[digits]
    return replace(
        control,
        multiplier=0,
        half_width=width,
        max_half_width=width,
        factor_count=FACTORS[digits],
        pool_size=64,
        polynomials_per_family=64,
    )


def endpoint_for_count(
    n, multiplier, count, budget, *, ceiling=ENDPOINT_CEILING
):
    """Find the minimal exclusive endpoint, or return a proper divisor.

    Include 2 and singular multiplier primes exactly as build_factor_base does.
    The finite sieve allocation is charged before construction. Euler's test is
    sufficient for primes; SIQS builds the roots later.
    """
    target = checked_target(n, multiplier)
    utils.require_integer(count, "eligible count", 1)
    utils.require_integer(ceiling, "endpoint ceiling", 3)
    if ceiling > ENDPOINT_CEILING:
        raise ValueError("endpoint exceeds the frozen helper allowance")
    budget.consume(n.bit_length())
    divisor = gcd(n, multiplier)
    if utils.valid_divisor(divisor, n):
        return None, divisor
    if divisor != 1:
        raise ValueError("multiplier is divisible by n")

    budget.consume(ceiling)
    found = 0
    for prime in prime_sieve.prime_sieve(ceiling):
        budget.consume(target.bit_length())
        divisor = gcd(prime, n)
        if utils.valid_divisor(divisor, n):
            return None, divisor
        residue = target % prime
        if (
            prime == 2
            or residue == 0
            or pow(residue, (prime - 1) // 2, prime) == 1
        ):
            found += 1
            if found == count:
                return prime + 1, None

    raise ConfigurationRefusedError(
        "eligible count is unreachable under endpoint allowance"
    )


def resolve(n, digits, arm, budget):
    """Resolve the joint policy inside the inclusive factoring allowance."""
    config = template(digits, arm)
    if arm == "control":
        return config, None
    choice = select_multiplier(n, budget=budget)
    if choice.divisor is not None:
        return config, choice.divisor
    bound, divisor = endpoint_for_count(
        n, choice.multiplier, COUNTS[digits], budget
    )
    if divisor is not None:
        return config, divisor
    return replace(
        config, multiplier=choice.multiplier, base_bound=bound
    ), None


def registry():
    """Serialize policy templates; actual endpoints are input dependent."""
    return {str(d): {a: asdict(template(d, a)) for a in ARMS} for d in BANDS}
