"""Exact early rejection and shared recovery for the combined challenger."""

from math import isqrt

from ..common import utils
from .cofactor_policy import split_cofactor
from .relations import AtomicRelation


def divide_yield_residual(
    collector, position, value, remaining, exponents, stats, *, hint=None
):
    if getattr(collector.config, "extension_policy", None):
        from .extended_collector import recover

        return recover(
            collector, position, value, remaining, exponents, stats, hint=hint
        )
    config, checked, budget = (
        collector.config,
        collector._checked,
        collector.budget,
    )
    slp, product, small = (
        config.residual_bound,
        config.large_product_bound,
        config.small_prime_bound,
    )
    if remaining > max(slp, product):
        stats["product_rejections"] += 1
        return None, None
    endpoint_product = config.large_prime_bound * (
        small or config.large_prime_bound
    )
    if remaining > max(slp, endpoint_product):
        stats["endpoint_product_rejections"] = (
            stats.get("endpoint_product_rejections", 0) + 1
        )
        return None, None
    if (
        config.recovery_policy == "checked-v2"
        and remaining > slp
        and (
            remaining < collector._residual_floor**2
            or (small is not None and small < collector._residual_floor)
        )
    ):
        budget.consume(1)
        stats["domain_hole_rejections"] = (
            stats.get("domain_hole_rejections", 0) + 1
        )
        return None, None

    pair, scalar = (), remaining
    if remaining != 1:
        known_prime = (
            config.recovery_policy == "checked-v2"
            and checked.primes.get(remaining) is True
        )
        if hint is not None and not utils.valid_divisor(hint, remaining):
            raise ValueError("invalid shared cofactor divisor")
        if hint is None and not known_prime and config.batch_policy != "off":
            hint = collector._endpoints.divisor(remaining, budget)
            stats["known_endpoint_splits"] = stats.get(
                "known_endpoint_splits", 0
            ) + (hint is not None)
        elif hint is not None:
            stats["shared_endpoint_splits"] = (
                stats.get("shared_endpoint_splits", 0) + 1
            )
        square = False
        if hint is None and not known_prime and remaining % 16 in (0, 1, 4, 9):
            budget.consume(remaining.bit_length() ** 2)
            root = isqrt(remaining)
            square = root * root == remaining
            if square:
                hint = root
        # A proper divisor or an exact square already proves compositeness;
        # do not repeat a whole-cofactor primality test before certification.
        if hint is None and checked.prime(remaining, budget):
            if remaining > slp:
                stats["prime_rejections"] += 1
                return None, None
        else:
            stats["composite_residuals"] += 1
            if remaining > min(product, endpoint_product):
                stats["product_rejections"] += 1
                return None, None
            if (
                config.split_limit is not None
                and collector._split_calls >= config.split_limit
            ):
                stats["split_quota_rejections"] += 1
                return None, None
            collector._split_calls += 1
            stats["split_attempts"] += 1
            pair, reason, evaluations = split_cofactor(
                remaining,
                config.large_prime_bound,
                budget,
                checked,
                config.split_policy,
                hint=hint,
                square_checked=True,
            )
            stats["rho_evaluations"] += evaluations
            if not pair:
                stats[reason] += 1
                return None, None
            if small is not None and pair[0] > small:
                stats["small_endpoint_rejections"] = (
                    stats.get("small_endpoint_rejections", 0) + 1
                )
                return None, None
            stats["double_residuals"] += 1
            stats["square_residuals"] += pair[0] == pair[1]
            scalar = 1
    atom = AtomicRelation(
        collector.polynomial,
        position,
        -1 if value < 0 else 1,
        exponents
        if getattr(config, "indexed_recovery", False)
        else tuple(exponents),
        scalar,
        large_primes=pair,
    )
    checked.atom(atom, budget)
    for prime in atom.residual_primes + (atom.square_correction,):
        divisor = checked.unit(prime, budget)
        if utils.valid_divisor(divisor, collector.factor_base.n):
            return atom, divisor
        if divisor != 1:
            raise ValueError("nonunit relation has no proper split")
    if config.batch_policy != "off":
        collector._endpoints.add(atom.residual_primes, budget)
    return atom, None


def prepare_admission(collector, atom, stats):
    checked, budget = collector._checked, collector.budget
    checked.atom(atom, budget)
    if getattr(collector.config, "extension_policy", None) == "extended-v2":
        from .opaque_staging import stage_unmatched

        if stage_unmatched(collector, atom, stats):
            return False
    key = None
    if collector.config.duplicate_policy == "congruence-v2":
        key, divisor = checked.congruence(atom, budget)
        if divisor is not None:
            raise ValueError(
                "proper divisor must be returned before admission"
            )
        if key is not None and checked.duplicate(
            key, collector._atoms, budget
        ):
            stats["congruence_duplicates"] = (
                stats.get("congruence_duplicates", 0) + 1
            )
            return False
    collector._incoming_key = key
    return True


def publish_auxiliary(collector, atom, row, slot, evicted):
    """Part of the collector's poisoned-on-allocation-failure transaction."""
    for old in evicted:
        collector._checked.forget(old)
        if collector._arena is not None:
            collector._arena.forget(old.relation_id)
    if collector._arena is not None:
        collector._arena.commit(atom, slot)
    collector._checked.remember(
        atom, row, collector._incoming_key, collector._atoms
    )
