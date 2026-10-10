"""Opt-in DLP intake and atomic publication into the relation store."""

from ..common import utils
from ..common.arithmetic import gcd
from .large_primes import CycleTooLongError, graph_reserve, split_two_primes
from .relations import (
    AtomicRelation,
    _combination_workspace,
    combine_relations,
    combined_storage_reserve,
    verify_atomic,
)

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
    if remaining > max(config.residual_bound, config.large_product_bound):
        stats["product_rejections"] += 1
        return None, None
    if remaining != 1:
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
            if collector._split_calls >= config.split_call_limit:
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
    if (
        collector._workspace + collector._scratch_bytes + reserve + 4096
        > config.memory_bytes
    ):
        return "memory_limit"
    endpoints = atom.large_primes or (1, atom.residual)
    try:
        plan = graph.plan(
            *endpoints, atom.relation_id, config.max_partials, budget
        )
    except MemoryError:
        return (
            "atom_limit"
            if len(collector._atoms) >= config.max_atoms
            else "memory_limit"
        )
    except CycleTooLongError:
        stats["long_cycle_rejections"] += 1
        return None
    if plan.dropped:
        stats["dropped_partials"] += 1
        return None
    accepted = plan.cycle is not None
    if accepted and len(collector._rows) >= config.max_relations:
        return "relation_limit"
    if len(collector._atoms) + 1 - len(plan.evicted) > config.max_atoms:
        return "atom_limit"
    old_edges = len(graph.edges)
    new_edges = old_edges + int(not accepted) - len(plan.evicted)
    graph_delta = graph_reserve(new_edges) - graph_reserve(old_edges)
    combination, atoms, scratch = None, (), 0
    if accepted:
        atoms = tuple(collector._atoms[i] for i in plan.cycle) + (atom,)
        reserve += combined_storage_reserve(
            sum(len(a.exponents) for a in atoms),
            collector.factor_base.n.bit_length(),
        )
        scratch = (
            _combination_workspace(
                atoms, collector.factor_base, config.memory_bytes
            )
            - collector.factor_base.workspace_bytes
        )
    if (
        collector._workspace
        + collector._scratch_bytes
        + reserve
        + scratch
        + max(0, graph_delta)
        > config.memory_bytes
    ):
        return "memory_limit"
    if accepted:
        collector._scratch_peak_bytes = max(
            collector._scratch_peak_bytes, collector._scratch_bytes + scratch
        )
        combination = combine_relations(
            atoms,
            collector.factor_base,
            budget=budget,
            memory_bytes=config.memory_bytes,
        )
        if combination.divisor is not None:
            raise ArithmeticError("previously checked unit residual changed")
    budget.consume(len(atom.exponents) + len(atoms) + len(plan.evicted) + 1)
    released = sum(collector._atom_bytes[i] for i in plan.evicted)
    graph.commit(plan)
    for identity in plan.evicted:
        del collector._atoms[identity]
        del collector._atom_bytes[identity]
    collector._atoms[atom.relation_id] = atom
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
    return None
