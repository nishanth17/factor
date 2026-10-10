"""Checked DLP stores; the default SLP encoding remains unchanged."""

from ..common import arithmetic, utils
from .large_primes import CycleTooLongError, LargePrimeForest, graph_reserve
from .polynomial import Polynomial, checked_position
from .relations import (
    MAX_COMBINED_ATOMS,
    AtomicRelation,
    CombinedRelation,
    _combination_workspace,
    combined_storage_reserve,
    verify_atomic,
    verify_combined,
)


def pack_store(collector):
    atoms = list(collector._atoms.values())
    indices = {a.relation_id: i for i, a in enumerate(atoms)}
    polynomials, lookup, encoded = [], {}, []
    for atom in atoms:
        polynomial = atom.polynomial
        key = polynomial.a, polynomial.b, polynomial.square_coefficient
        if key not in lookup:
            lookup[key] = len(polynomials)
            polynomials.append(key)
        encoded.append(
            [
                lookup[key],
                atom.position,
                atom.sign,
                atom.exponents,
                atom.residual,
                atom.large_primes,
            ]
        )
    rows = collector._full + collector._combined
    order = {id(row): i for i, row in enumerate(rows)}
    return dict(
        polynomials=polynomials,
        atoms=encoded,
        full=[indices[a.relation_id] for a in collector._full],
        combined=[
            [
                [indices[i] for i in row.atom_ids],
                row.u,
                row.sign,
                row.exponents,
                row.square_correction,
            ]
            for row in collector._combined
        ],
        row_order=[order[id(row)] for row in collector._rows],
        forest=[indices[i] for i in collector._graph.edges],
        split_calls=collector._split_calls,
    )


def _index(value, size):
    utils.require_integer(value, "atomic index", 0)
    if value >= size:
        raise ValueError("checkpoint refers outside retained provenance")
    return value


def restore_store(payload, collector, budget):
    """Rebuild only after bounded shape checks; publish a fully checked store.

    Every combined row has exactly one distinct closing edge, and its other
    IDs must equal the unique path in the retained forest. Owned paths cannot
    have been evicted. This also reconstructs FIFO order without trusting an
    encoded union-find cache, and verifies disconnected cycles and loops.
    """
    if collector._atoms or collector._rows or collector._graph.edges:
        raise ValueError("DLP restoration requires an empty collector")
    if not isinstance(payload, dict) or set(payload) != {
        "polynomials",
        "atoms",
        "full",
        "combined",
        "row_order",
        "forest",
        "split_calls",
    }:
        raise ValueError("invalid DLP checkpoint store shape")
    config, base = collector.config, collector.factor_base
    for name, cap in (
        ("polynomials", config.max_atoms),
        ("atoms", config.max_atoms),
        ("full", config.max_relations),
        ("combined", config.max_relations),
        ("row_order", config.max_relations),
        ("forest", config.max_atoms),
    ):
        if not isinstance(payload[name], list) or len(payload[name]) > cap:
            raise ValueError("DLP checkpoint exceeds its store cap")
    calls = payload["split_calls"]
    utils.require_integer(calls, "split_calls", 0)
    if calls > config.split_call_limit:
        raise ValueError("DLP checkpoint exceeds its splitting quota")
    count = len(payload["atoms"])
    row_count = len(payload["full"]) + len(payload["combined"])
    order = payload["row_order"]
    if row_count > config.max_relations or len(order) != row_count:
        raise ValueError("invalid DLP mixed row count")
    for index in order:
        _index(index, row_count)
    if len(set(order)) != row_count:
        raise ValueError("DLP row order repeats an original row")
    reserves, used_polynomials = [], set()
    for record in payload["atoms"]:
        budget.consume(1)
        if not isinstance(record, list) or len(record) != 6:
            raise ValueError("invalid DLP atomic record")
        poly_index, position, _, exponents, _, pair = record
        _index(poly_index, len(payload["polynomials"]))
        used_polynomials.add(poly_index)
        checked_position(position)
        if (
            not isinstance(exponents, list)
            or len(exponents) > len(base.entries)
            or not isinstance(pair, list)
            or len(pair) not in (0, 2)
        ):
            raise ValueError("unbounded DLP atomic payload")
        reserves.append(
            4096
            + 256 * len(exponents)
            + 16 * (abs(position).bit_length() + base.n_prime.bit_length())
        )
    if len(used_polynomials) != len(payload["polynomials"]):
        raise ValueError("DLP checkpoint retains an unreferenced polynomial")
    combined_bytes = []
    for record in payload["combined"]:
        budget.consume(1)
        if not isinstance(record, list) or len(record) != 5:
            raise ValueError("invalid DLP combined record")
        indices, u, _, exponents, correction = record
        if (
            not isinstance(indices, list)
            or not 1 <= len(indices) <= MAX_COMBINED_ATOMS
        ):
            raise ValueError("DLP cycle provenance exceeds its cap")
        for index in indices:
            _index(index, count)
        if len(set(indices)) != len(indices):
            raise ValueError("DLP cycle reuses an atomic position")
        if not isinstance(exponents, list) or len(exponents) > len(
            base.entries
        ):
            raise ValueError("unbounded DLP combined exponents")
        for value in (u, correction):
            utils.require_integer(value, "modular payload", 0)
            if value >= base.n:
                raise ValueError("noncanonical DLP modular payload")
        combined_bytes.append(
            combined_storage_reserve(
                sum(len(payload["atoms"][i][3]) for i in indices),
                base.n.bit_length(),
            )
        )
    forest_indices = payload["forest"]
    previous = -1
    for index in forest_indices:
        _index(index, count)
        if index <= previous:
            raise ValueError("DLP forest is not in unique admission order")
        previous = index
    reserve = collector._workspace + sum(reserves) + sum(combined_bytes)
    reserve += graph_reserve(len(forest_indices)) - graph_reserve(0)
    if reserve > config.memory_bytes:
        raise MemoryError(
            "restored DLP provenance and graph exceed memory cap"
        )
    polynomials = []
    for record in payload["polynomials"]:
        budget.consume(1)
        if not isinstance(record, list) or len(record) != 3:
            raise ValueError("invalid DLP polynomial record")
        polynomials.append(Polynomial(base.n, base.multiplier, *record))
    if len({p.identity for p in polynomials}) != len(polynomials):
        raise ValueError("duplicate DLP polynomial")
    atoms, store, sizes = [], {}, {}
    for i, (poly_index, position, sign, exponents, scalar, pair) in enumerate(
        payload["atoms"]
    ):
        _index(poly_index, len(polynomials))
        atom = AtomicRelation(
            polynomials[poly_index],
            position,
            sign,
            tuple(tuple(item) for item in exponents),
            scalar,
            large_primes=tuple(pair),
        )
        verify_atomic(
            atom,
            base,
            residual_bound=config.residual_bound,
            large_prime_bound=config.large_prime_bound,
            large_product_bound=config.large_product_bound,
            budget=budget,
        )
        if atom.relation_id in store:
            raise ValueError("duplicate DLP atomic identity")
        for prime in atom.residual_primes:
            budget.consume(base.n.bit_length())
            if arithmetic.gcd(prime, base.n) != 1:
                raise ValueError("nonunit residual cannot enter a DLP store")
        atoms.append(atom)
        store[atom.relation_id], sizes[atom.relation_id] = atom, reserves[i]
    graph, forest = LargePrimeForest(), set(forest_indices)
    for index in forest_indices:
        atom = atoms[index]
        if atom.residual_product == 1:
            raise ValueError("full atom cannot be a forest edge")
        endpoints = atom.large_primes or (1, atom.residual)
        graph.commit_link(
            graph.plan_link(*endpoints, atom.relation_id, budget)
        )
    full, combined, used, closing_edges = [], [], set(forest), set()
    for index in payload["full"]:
        _index(index, count)
        if index in used or atoms[index].residual_product != 1:
            raise ValueError("invalid full DLP checkpoint atom")
        used.add(index)
        full.append(atoms[index])
    scratch_peak = 0
    for record, row_bytes in zip(payload["combined"], combined_bytes):
        indices, u, sign, exponents, correction = record
        closing = set(indices) - forest
        if len(closing) != 1:
            raise ValueError("DLP cycle must have one closing edge")
        closing = closing.pop()
        if (
            closing in closing_edges
            or closing != max(indices)
            or closing in used
        ):
            raise ValueError("DLP closing edge is reused or out of order")
        atom = atoms[closing]
        if atom.residual_product == 1:
            raise ValueError("full atom is marked as a DLP cycle")
        endpoints = atom.large_primes or (1, atom.residual)
        try:
            path = graph.path(*endpoints, budget)
        except CycleTooLongError as error:
            raise ValueError("checkpoint DLP cycle exceeds its cap") from error
        expected = {atoms[i].relation_id for i in indices if i != closing}
        if path is None or set(path) != expected:
            raise ValueError("DLP cycle differs from its retained forest path")
        selected = tuple(atoms[i] for i in indices)
        scratch = (
            _combination_workspace(selected, base, config.memory_bytes)
            - base.workspace_bytes
        )
        scratch_peak = max(scratch_peak, scratch)
        if reserve + scratch > config.memory_bytes:
            raise MemoryError("restored DLP verification scratch exceeds cap")
        item = CombinedRelation(
            tuple(a.relation_id for a in selected),
            arithmetic.backend_for(base.n).integer(u),
            sign,
            tuple(tuple(pair) for pair in exponents),
            arithmetic.backend_for(base.n).integer(correction),
        )
        verify_combined(
            item, base, store, budget=budget, memory_bytes=config.memory_bytes
        )
        for identity in path:
            graph.unowned.pop(identity, None)
        used.update(indices)
        closing_edges.add(closing)
        sizes[atom.relation_id] += row_bytes
        combined.append(item)
    if len(used) != count or len(graph.unowned) > config.max_partials:
        raise ValueError("DLP checkpoint has unowned overflow or orphan atoms")
    rows = full + combined
    budget.consume(count + row_count + 1)
    collector._atoms, collector._atom_bytes = store, sizes
    collector._full, collector._combined = full, combined
    collector._rows = [rows[i] for i in order]
    collector._graph, collector._split_calls = graph, calls
    collector._workspace = reserve
    collector._scratch_peak_bytes = scratch_peak
