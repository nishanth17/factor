"""Explicit SIQS capacity diagnostics and monotone resume allowances."""

from dataclasses import asdict
from math import prod

from .linear_algebra import matrix_workspace
from .polynomial import a_target


def preparation_cache_reserve(config):
    """Return the maximum verified-row cache reservation used by QSJob."""
    return min(2 * 2**20, 2048 * config.collector.max_relations)


def validate_extension(previous, proposed):
    """Permit only larger resource/search quotas under the same mathematics.

    Reference assignment lists depend on their length and cannot be extended.
    Streaming identities exclude the total quota and preserve checked prefixes.
    Checkpoint and cache growth must not reduce remaining live workspace.
    """
    if type(proposed) is not type(previous):
        raise TypeError("extension requires the same configuration type")
    old, new = asdict(previous), asdict(proposed)
    quotas = {"max_stalled", "max_trivial", "memory_bytes", "checkpoint_bytes"}
    if previous.streaming:
        quotas.add("family_count")
    if previous.external_coefficients:
        quotas.add("coefficient_trials")
    for name in quotas:
        if new[name] < old[name]:
            raise ValueError(f"extension cannot reduce {name}")
        new[name] = old[name]
    for name in ("max_atoms", "max_partials", "max_relations"):
        if new["collector"][name] < old["collector"][name]:
            raise ValueError(f"extension cannot reduce {name}")
        new["collector"][name] = old["collector"][name]
    if new != old:
        raise ValueError(
            "extension changes mathematical or execution configuration"
        )
    if (
        proposed.memory_bytes
        - proposed.metadata_reserve
        - preparation_cache_reserve(proposed)
        < previous.memory_bytes
        - previous.metadata_reserve
        - preparation_cache_reserve(previous)
    ):
        raise ValueError(
            "extension reduces live workspace after checkpoint/cache reserves"
        )


def capacity_report(base, config):
    """Report exact reachability bounds and owned-storage reservations.

    Matrix storage alone is a necessary lower bound on the full job allowance;
    it excludes coexisting collector/provenance storage and is not process RSS.
    An attainable interval does not prove the selected A pool reaches a target.
    """
    eligible = tuple(p for p in base.primes if p != 2 and base.n_prime % p)
    count = config.factor_count
    target = a_target(base.n_prime, config.half_width)
    enough = len(eligible) >= count
    low = prod(eligible[:count]) if enough else None
    high = prod(eligible[-count:]) if enough else None
    matrix_bytes = matrix_workspace(
        config.collector.max_relations, len(base.entries) + 1
    )
    available = (
        config.memory_bytes
        - config.metadata_reserve
        - base.workspace_bytes
        - preparation_cache_reserve(config)
    )
    return dict(
        base_bound=base.bound,
        base_cardinality=len(base.entries),
        a_target=target,
        a_minimum=low,
        a_maximum=high,
        factor_base_product_envelope_applies=not config.external_coefficients,
        target_in_product_envelope=enough and low <= target <= high,
        family_allowance=config.family_count,
        polynomials_per_family=config.gray_limit,
        polynomial_allowance=config.polynomial_limit,
        residual_bound=config.collector.residual_bound,
        atom_allowance=config.collector.max_atoms,
        partial_allowance=config.collector.max_partials,
        relation_allowance=config.collector.max_relations,
        matrix_reserve_bytes=matrix_bytes,
        metadata_reserve_bytes=config.metadata_reserve,
        preparation_cache_reserve_bytes=preparation_cache_reserve(config),
        base_reserve_bytes=base.workspace_bytes,
        matrix_alone_fits=matrix_bytes <= available,
        owned_cap_bytes=config.memory_bytes,
    )
