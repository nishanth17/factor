"""Finite transactional conversion of opaque factors into endpoints."""

from .residual_kernel import ResidualKernel


def convert(collector, residual, primes, budget, *, replay=False):
    config = collector.config
    if residual in collector._graph.resolutions:
        return False
    if (
        not replay
        and collector._conversion_attempts >= config.opaque_conversions
    ):
        return False
    budget.consume(1)
    if frozenset(((True, residual),)) not in collector._graph.direct:
        return False
    collector._conversion_attempts += 1
    product = 1
    for prime in primes:
        if (
            type(prime) is not int
            or not 2
            <= prime
            <= max(config.residual_bound, config.large_prime_bound)
            or not collector._checked.prime(prime, budget)
        ):
            raise ValueError("opaque conversion requires proven endpoints")
        product *= prime
    if product != residual or not 1 <= len(primes) <= config.residual_arity:
        raise ValueError("opaque conversion does not reconstruct")
    bounded = config.extension_policy == "extended-v2"
    scratch = 2**21 if bounded else 2**20 + 32768 * len(collector._atoms)
    if collector._workspace + scratch > config.memory_limit:
        if replay:
            raise MemoryError("opaque conversion replay exceeds its cap")
        return False
    kernel = ResidualKernel()
    kernel.resolutions = dict(collector._graph.resolutions)
    kernel.resolutions[residual] = tuple(primes)
    existing = {
        tuple(sorted(collector._arena.slots[i] for i in row.atom_ids))
        for row in collector._rows
    }
    added, reserve = [], 0
    for handle, atom in collector._arena.atoms.items():
        plan = kernel.plan_atom(atom, handle, config, budget)
        if bounded:
            from .yield_memory import graph_reserve

            reserve += graph_reserve(plan, config.cycle_reduction_choices)
            if collector._workspace + scratch + reserve > config.memory_limit:
                if replay:
                    raise MemoryError(
                        "opaque rebuilt incidence exceeds its cap"
                    )
                return False
        if plan.cycle is not None and plan.cycle not in existing:
            leaves = tuple(collector._arena.atoms[i] for i in plan.cycle)
            reserve += 8192 + 512 * sum(len(a.exponents) for a in leaves)
            if (
                collector._workspace + scratch + reserve > config.memory_limit
                or (len(collector._rows) + len(added) >= config.relation_limit)
            ):
                if replay:
                    raise MemoryError(
                        "opaque conversion rows exceed their cap"
                    )
                return False
            row = collector._checked.combine(
                leaves, budget, config.memory_limit
            ).relation
            if row is None:
                raise ValueError("opaque conversion encountered a nonunit")
            added.append(row)
            existing.add(plan.cycle)
        kernel.commit_atom(plan, atom, config.cycle_reduction_choices)
    budget.consume(len(added) + len(collector._atoms) + 1)
    # Allocate replacement lists before publishing the coupled graph/row state.
    rows = collector._rows + added
    combined = collector._combined + added
    events = collector._conversion_events + [
        [len(collector._atoms), residual, list(primes)]
    ]
    observed_rows = None if replay else collector._reporting.rows + len(added)
    collector._graph, collector._rows, collector._combined = (
        kernel,
        rows,
        combined,
    )
    collector._conversion_events = events
    if not replay:
        # Conversion commits independently of the incoming atom. Preserve its
        # yield even if that atom subsequently refuses or needs a new budget.
        collector._reporting.rows = observed_rows
    collector._workspace += reserve
    collector._conversion_bytes = (
        getattr(collector, "_conversion_bytes", 0) + reserve
    )
    return True
