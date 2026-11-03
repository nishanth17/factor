"""Transactional generalized residual admission for the expanded policy."""

from .storage import RelationPublicationError, require_published


def recover(
    collector, position, value, residual, factors, stats, *, hint=None
):
    from math import isqrt

    from ..common import utils
    from .cofactor_policy import split_residual
    from .relations import AtomicRelation
    from .residual_relations import ResidualAtom

    config, checked, budget = (
        collector.config,
        collector._checked,
        collector.budget,
    )
    if not config.indexed_recovery:
        factors = tuple(factors)
    if residual > max(config.residual_bound, config.large_product_bound):
        stats["product_rejections"] += 1
        return None, None
    known_prime = checked.primes.get(residual) is True
    if residual > 1 and hint is None and not known_prime:
        if config.batch_policy != "off":
            hint = collector._endpoints.divisor(residual, budget)
        if hint is None and residual % 16 in (0, 1, 4, 9):
            budget.consume(residual.bit_length() ** 2)
            root = isqrt(residual)
            if root * root == residual:
                hint = root
    if hint is not None and not utils.valid_divisor(hint, residual):
        raise ValueError("invalid expanded cofactor hint")
    sign = -1 if value < 0 else 1
    if residual == 1:
        atom = AtomicRelation(collector.polynomial, position, sign, factors)
    elif hint is None and checked.prime(residual, budget):
        if residual > config.residual_bound:
            stats["prime_rejections"] += 1
            return None, None
        atom = ResidualAtom(
            collector.polynomial, position, sign, factors, primes=(residual,)
        )
    else:
        stats["composite_residuals"] += 1
        if residual > config.large_product_bound:
            stats["product_rejections"] += 1
            return None, None
        if residual < collector._residual_floor**2:
            stats["domain_hole_rejections"] = (
                stats.get("domain_hole_rejections", 0) + 1
            )
            return None, None
        primes, reason, evaluations = (), "split_failed", 0
        small = config.small_prime_bound or config.large_prime_bound
        split_eligible = residual <= small * config.large_prime_bound ** (
            config.residual_arity - 1
        )
        if config.residual_arity == 3 and (
            small * config.large_prime_bound
            < residual
            < collector._residual_floor**3
        ):
            split_eligible = False
        if not split_eligible:
            stats["split_domain_skips"] = (
                stats.get("split_domain_skips", 0) + 1
            )
            reason = "endpoint_bound"
        elif collector._split_calls < config.split_limit:
            collector._split_calls += 1
            stats["split_attempts"] += 1
            primes, reason, evaluations = split_residual(
                residual,
                config.large_prime_bound,
                budget,
                checked,
                config.split_policy,
                arity=config.residual_arity,
                hint=hint,
                dispatch=config.split_dispatch,
                memo=collector._split_memo,
            )
        else:
            stats["split_quota_rejections"] += 1
        stats["rho_evaluations"] += evaluations
        if (
            primes
            and config.small_prime_bound is not None
            and (primes[0] > config.small_prime_bound)
        ):
            primes, reason = (), "small_endpoint_rejections"
        if not primes:
            stats[reason] = stats.get(reason, 0) + 1
            if not config.opaque_pairs:
                return None, None
        atom = ResidualAtom(
            collector.polynomial,
            position,
            sign,
            factors,
            primes=primes,
            opaque=None if primes else residual,
        )
        counter = (
            "opaque_residuals"
            if not primes
            else "triple_residuals"
            if len(primes) == 3
            else "double_residuals"
        )
        stats[counter] = stats.get(counter, 0) + 1
    checked.atom(atom, budget)
    for token in atom.residual_primes + (atom.square_correction,):
        divisor = checked.unit(token, budget)
        if utils.valid_divisor(divisor, collector.factor_base.n):
            return atom, divisor
        if divisor != 1:
            raise ValueError("nonunit residual without proper divisor")
    if getattr(atom, "opaque", None) is None and config.batch_policy != "off":
        collector._endpoints.add(atom.residual_primes, budget)
    return atom, None


def atom_reserve(atom, config, base, plan=None):
    """Retained ownership; construction scratch is reserved separately."""
    if config.extension_policy == "extended-v2":
        from .opaque_staging import context_bits
        from .yield_memory import graph_reserve

        # One worst-case 4 KiB slab per record plus payload/map overhead;
        # incidence ownership is separate from temporary lifting workspace.
        return (
            8192
            + (0 if config.spill_bytes else 512 * len(atom.exponents))
            + 32 * context_bits(atom, base)
            + (
                0
                if plan is None
                else graph_reserve(plan, config.cycle_reduction_choices)
            )
        )
    return (
        32768
        + (0 if config.spill_bytes else 512 * len(atom.exponents))
        + 32 * base.n_prime.bit_length()
    )


def admit(collector, atom, stats):
    if collector.config.extension_policy == "extended-v2":
        from .opaque_staging import route

        handled, refusal = route(collector, atom, stats)
        if handled:
            return refusal
    return admit_retained(collector, atom, stats)


def admit_retained(collector, atom, stats):
    require_published(collector)
    config, budget, arena = (
        collector.config,
        collector.budget,
        collector._arena,
    )
    if atom.relation_id in collector._atoms:
        stats["duplicates"] += 1
        return None
    if getattr(atom, "opaque", None) is None and atom.residual_primes:
        from .opaque_conversion import convert

        if convert(
            collector, atom.residual_product, atom.residual_primes, budget
        ):
            stats["opaque_conversions"] = (
                stats.get("opaque_conversions", 0) + 1
            )
    if len(collector._atoms) >= config.atom_limit:
        return "atom_limit"
    if (
        not atom.residual_primes
        and len(collector._rows) >= config.relation_limit
    ):
        return "relation_limit"
    # Covers incidence vectors/lifting up to512 handles, simultaneous packing,
    # optional recomputation descriptors and a maximum-sized flat lazy node.
    reserve = atom_reserve(atom, config, collector.factor_base)
    construction = 512 * len(atom.exponents) if config.spill_bytes else 0
    if (
        collector._workspace
        + collector._scratch_bytes
        + reserve
        + construction
        > config.memory_limit
    ):
        return "memory_limit"
    try:
        plan = collector._graph.plan_atom(
            atom, arena.next_slot, config, budget
        )
    except MemoryError:
        return "memory_limit"
    reserve = atom_reserve(atom, config, collector.factor_base, plan)
    if (
        collector._workspace
        + collector._scratch_bytes
        + reserve
        + construction
        > config.memory_limit
    ):
        return "memory_limit"
    row = None
    if plan.cycle is not None:
        if len(collector._rows) >= config.relation_limit:
            return "relation_limit"
        leaves = tuple(
            atom if i == arena.next_slot else arena.atoms[i]
            for i in plan.cycle
        )
        extra = 8192 + 512 * sum(len(a.exponents) for a in leaves)
        if (
            collector._workspace + reserve + extra + construction
            > config.memory_limit
        ):
            return "memory_limit"
        result = collector._checked.combine(
            leaves, budget, config.memory_limit
        )
        if result.divisor is not None:
            raise ArithmeticError("checked unit residual changed")
        row = result.relation
        reserve += extra
    packed, row, slot = arena.plan(atom, row, budget)
    budget.consume(len(atom.exponents) + len(plan.provenance) + 1)
    collector._publication_failed = True
    try:
        arena.commit(packed, slot)
        collector._atoms[packed.relation_id] = packed
        collector._atom_bytes[packed.relation_id] = reserve
        collector._graph.commit_atom(
            plan, packed, config.cycle_reduction_choices
        )
        collector._checked.remember(
            packed, row, collector._incoming_key, collector._atoms
        )
        collector._workspace += reserve
        stats["admitted_atoms"] += 1
        if row is not None:
            collector._combined.append(row)
            collector._rows.append(row)
            stats["matches"] += 1
            stats["graph_cycles"] += 1
            stats["cycle_atoms"] += len(plan.cycle)
            collector._reporting.record(0, 0, 1)
    except (MemoryError, OSError) as error:
        raise RelationPublicationError(
            "expanded relation publication failed"
        ) from error
    else:
        collector._publication_failed = False
    return None
