"""Opt-in DLP intake and atomic publication into the relation store."""

from dataclasses import replace

from ..common import utils
from ..common.arithmetic import gcd
from .large_primes import CycleTooLongError, split_two_primes
from .relations import (
    AtomicRelation,
    _combination_workspace,
    combine_relations,
    combined_storage_reserve,
    verify_atomic,
)
from .storage import RelationPublicationError, require_published

COUNTERS = (
    "double_residuals",
    "square_residuals",
    "split_attempts",
    "rho_evaluations",
    "split_failed",
    "endpoint_bound",
    "not_two_primes",
    "product_rejections",
    "prime_rejections",
    "split_quota_rejections",
    "graph_cycles",
    "dlp_cycles",
    "cycle_atoms",
    "long_cycle_rejections",
    "graph_rebuilds",
)


def divide_large_residual(
    collector, position, value, remaining, exponents, stats
):
    """Certify both endpoint identities before any graph/store mutation."""
    config, budget, base = (
        collector.config,
        collector.budget,
        collector.factor_base,
    )
    pair, scalar = (), remaining
    small_bound = getattr(config, "small_prime_bound", None)
    endpoint_product = (
        None if small_bound is None else small_bound * config.large_prime_bound
    )
    if remaining > max(config.residual_bound, config.large_product_bound):
        stats["product_rejections"] += 1
        return None, None
    if remaining != 1:
        if endpoint_product is not None:
            budget.consume(1)
            # p<=S and q<=L implies pq<=SL. Preserve the independent SLP
            # allowance: a prime above SL may still be a valid scalar atom.
            if remaining > max(config.residual_bound, endpoint_product):
                stats["endpoint_product_rejections"] = (
                    stats.get("endpoint_product_rejections", 0) + 1
                )
                return None, None
        budget.consume(remaining.bit_length() ** 2)
        prime = utils.classify_prime(remaining) is utils.Primality.PROVEN
        if prime:
            if remaining > config.residual_bound:
                stats["prime_rejections"] += 1
                return None, None
        else:
            stats["composite_residuals"] += 1
            if remaining > config.large_product_bound:
                stats["product_rejections"] += 1
                return None, None
            if endpoint_product is not None and remaining > endpoint_product:
                stats["endpoint_product_rejections"] = (
                    stats.get("endpoint_product_rejections", 0) + 1
                )
                return None, None
            if (
                config.split_limit is not None
                and collector._split_calls >= config.split_limit
            ):
                stats["split_quota_rejections"] += 1
                return None, None
            # A started attempt remains charged even if a later cancellation
            # prevents publication. Resume cannot obtain free cofactoring.
            collector._split_calls += 1
            stats["split_attempts"] += 1
            pair, reason, evaluations = split_two_primes(
                remaining,
                config.large_prime_bound,
                budget,
                known_composite=True,
            )
            stats["rho_evaluations"] += evaluations
            if not pair:
                stats[reason] += 1
                return None, None
            # Reporting and recovered-product admission are separate. An
            # ordered policy pays for the exact split before rejecting its
            # smaller endpoint; it cannot publish an unchecked square loop.
            if small_bound is not None and pair[0] > small_bound:
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
        tuple(exponents),
        scalar,
        large_primes=pair,
    )
    verify_atomic(
        atom,
        base,
        residual_bound=config.residual_bound,
        large_prime_bound=config.large_prime_bound,
        large_product_bound=config.large_product_bound,
        small_prime_bound=getattr(config, "small_prime_bound", None),
        budget=budget,
    )
    for prime in atom.residual_primes:
        budget.consume(base.n.bit_length())
        divisor = gcd(prime, base.n)
        if utils.valid_divisor(divisor, base.n):
            return atom, divisor
        if divisor != 1:
            raise ValueError("nonunit large prime has no proper split")
    return atom, None


def admit_large_relation(collector, atom, stats):
    """Plan first, verify exact provenance, then commit a bounded change."""
    require_published(collector)
    config, graph, budget = (
        collector.config,
        collector._graph,
        collector.budget,
    )
    if atom.relation_id in collector._atoms:
        stats["duplicates"] += 1
        return None
    reserve = 4096 + 256 * len(atom.exponents)
    reserve += 16 * (
        abs(atom.position).bit_length()
        + collector.factor_base.n_prime.bit_length()
    )
    # The existing graph reservation covers two simultaneous map sets. Allow
    # one new forest edge and incoming atom before planning a staged rebuild;
    # old atoms remain live until commit, so do not subtract evictions here.
    endpoints = atom.large_primes or (1, atom.residual)
    shortcuts = collector._short_cycles
    try:
        direct = None
        if (
            shortcuts is not None
            and shortcuts.strategy == "short-v2"
            and shortcuts.edge_visits
        ):
            budget.consume(1)
            direct = shortcuts.neighbors.get(endpoints[0], {}).get(
                endpoints[1]
            )
        if direct or (
            shortcuts is not None
            and shortcuts.strategy == "short-v2"
            and endpoints[0] == endpoints[1]
        ):
            from .large_primes import GraphPlan

            path = (
                ()
                if endpoints[0] == endpoints[1]
                else (shortcuts._representative(direct, budget),)
            )
            # This edge already closes a cycle. Avoid walking the original
            # forest merely to rediscover a direct prior edge or square loop.
            plan = GraphPlan(cycle=path)
            stats["direct_cycle_shortcuts"] = (
                stats.get("direct_cycle_shortcuts", 0) + 1
            )
        else:
            plan = graph.plan(
                *endpoints, atom.relation_id, config.max_partials, budget
            )
    except CycleTooLongError:
        if shortcuts is None:
            stats["long_cycle_rejections"] += 1
            return None
        from .large_primes import GraphPlan

        plan = GraphPlan(cycle=None)
    if shortcuts is not None and (
        plan.cycle is not None
        or (plan.link is None and plan.rebuilt is None and not plan.dropped)
    ):
        path, visits = shortcuts.shorter_path(*endpoints, plan.cycle, budget)
        stats["shortcut_visits"] = stats.get("shortcut_visits", 0) + visits
        if path is None:
            stats["long_cycle_rejections"] += 1
            return None
        saved = len(plan.cycle) - len(path) if plan.cycle is not None else 0
        stats["shortcut_atoms_saved"] = (
            stats.get("shortcut_atoms_saved", 0) + saved
        )
        stats["shortcut_cycles"] = stats.get("shortcut_cycles", 0) + (
            path != plan.cycle
        )
        plan = replace(plan, cycle=path)
    if plan.dropped:
        stats["dropped_partials"] += 1
        return None
    accepted = plan.cycle is not None
    old_edges = len(graph.edges)
    new_edges = old_edges + int(not accepted) - len(plan.evicted)
    graph_delta = collector._graph_reserve(
        new_edges
    ) - collector._graph_reserve(old_edges)
    if shortcuts is not None:
        graph_delta += 2048 * (1 - len(plan.evicted))
    combination, atoms, scratch = None, (), 0
    if accepted:
        atoms = tuple(collector._atoms[i] for i in plan.cycle) + (atom,)
        reserve += combined_storage_reserve(
            sum(len(a.exponents) for a in atoms),
            collector.factor_base.n.bit_length(),
        )
        scratch = (
            _combination_workspace(atoms, collector.factor_base)
            - collector.factor_base.workspace_bytes
        )
    if accepted:
        collector._scratch_peak_bytes = max(
            collector._scratch_peak_bytes, collector._scratch_bytes + scratch
        )
        combination = (
            collector._checked.combine(atoms, budget)
            if collector._checked is not None
            and config.recovery_policy == "checked-v2"
            else combine_relations(
                atoms,
                collector.factor_base,
                budget=budget,
            )
        )
        if combination.divisor is not None:
            raise ArithmeticError("previously checked unit residual changed")
    compact_slot = None
    if collector._arena is not None:
        atom, row, compact_slot = collector._arena.plan(
            atom, None if combination is None else combination.relation, budget
        )
        if combination is not None:
            combination = replace(combination, relation=row)
    elif getattr(config, "record_policy", "objects") == "packed-v1":
        from .packed_relations import pack_relation

        atom = pack_relation(atom, budget)
        if combination is not None:
            combination = replace(
                combination,
                relation=pack_relation(combination.relation, budget),
            )
    budget.consume(len(atom.exponents) + len(atoms) + len(plan.evicted) + 1)
    if shortcuts is not None:
        budget.consume(4 * (1 + len(plan.evicted)))
    released = sum(collector._atom_bytes[i] for i in plan.evicted)
    # Containers allocate during insertion. If the allocator fails midway,
    # old external checkpoints remain valid, but this live store must never
    # masquerade as a resumable resource refusal or publish a partial row.
    collector._publication_failed = True
    collector.budget.publication()
    try:
        evicted_atoms = tuple(collector._atoms[i] for i in plan.evicted)
        graph.commit(plan)
        if shortcuts is not None:
            for identity in plan.evicted:
                shortcuts.remove(identity)
            shortcuts.add(
                atom.relation_id, *endpoints, weight=(len(atom.exponents),)
            )
        for identity in plan.evicted:
            del collector._atoms[identity]
            del collector._atom_bytes[identity]
        collector._atoms[atom.relation_id] = atom
        if collector._checked is not None:
            from .yield_collector import publish_auxiliary

            publish_auxiliary(
                collector,
                atom,
                None if combination is None else combination.relation,
                compact_slot,
                evicted_atoms,
            )
        collector._atom_bytes[atom.relation_id] = reserve
        collector._workspace += reserve - released + graph_delta
        stats["admitted_atoms"] += 1
        stats["evictions"] += len(plan.evicted)
        stats["graph_rebuilds"] += plan.rebuilt is not None
        if accepted:
            collector._combined.append(combination.relation)
            collector._rows.append(combination.relation)
            stats["matches"] += 1
            stats["graph_cycles"] += 1
            stats["dlp_cycles"] += any(a.large_primes for a in atoms)
            stats["cycle_atoms"] += len(atoms)
    except MemoryError as error:
        raise RelationPublicationError(
            "allocation failed during relation publication"
        ) from error
    else:
        collector._publication_failed = False
    return None
