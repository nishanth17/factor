"""Schema13 structural checkpoints with independently checked restore."""

import base64
from types import SimpleNamespace

from .cold_factors import ColdFactors, RecomputedFactors, SpillStore
from .indexed_relations import (
    FactorRecord,
    IndexedAtomicRelation,
    IndexedResidualAtom,
    RecordArena,
    encode_factors,
    import_factors,
)
from .polynomial import Polynomial
from .relations import AtomicRelation
from .residual_relations import ResidualAtom
from .storage import require_published


def pack_factors(factors, base):
    if isinstance(factors, RecomputedFactors):
        return {
            "recompute": [
                list(factors.small),
                len(factors),
                pack_factors(factors.retained, base),
            ]
        }
    if isinstance(factors, ColdFactors):
        return {"cold": [list(factors.reference), len(factors)]}
    data = (
        factors.encoded()
        if isinstance(factors, FactorRecord)
        else encode_factors(factors, base)
    )
    return {
        "index-v2": base64.b64encode(data).decode("ascii"),
        "count": len(factors),
    }


def pack_store(collector):
    require_published(collector)
    atoms = list(collector._atoms.values())
    indices = {atom.relation_id: i for i, atom in enumerate(atoms)}
    encoded = []
    for atom in atoms:
        poly = atom.polynomial
        encoded.append(
            [
                [poly.a, poly.b, poly.square_coefficient],
                atom.position,
                atom.sign,
                pack_factors(atom.exponents, collector.factor_base),
                list(atom.residual_primes)
                if getattr(atom, "opaque", None) is None
                else [],
                getattr(atom, "opaque", None),
            ]
        )
    spill = collector._arena.spill
    return {
        "extended_store": 1,
        "atoms": encoded,
        "rows": [
            [indices[i] for i in row.atom_ids] for row in collector._rows
        ],
        "row_order": list(range(len(collector._rows))),
        "kernel": collector._graph.state(),
        "conversion_events": collector._conversion_events,
        "conversion_attempts": collector._conversion_attempts,
        "reporting": collector._reporting.state(),
        "split_calls": collector._split_calls,
        "endpoints": list(collector._endpoints.primes),
        "memo": [
            [list(key), list(result[0]), result[1]]
            for key, result in collector._split_memo.values.items()
            if result[0]
        ],
        "canonical": [indices[i] for i in collector._checked.canonical_ids],
        "spill": None if spill is None else spill.manifest(),
    }


def restore_factors(
    payload, base, poly, position, config, spill, budget, *, nested=False
):
    if not isinstance(payload, dict):
        raise ValueError("invalid expanded factor record")
    if set(payload) == {"recompute"}:
        if nested or not config.recompute_prime_limit:
            raise ValueError("unexpected recomputation descriptor")
        small, count, retained = payload["recompute"]
        if (
            not isinstance(small, list)
            or small != sorted(set(small))
            or any(
                type(p) is not int
                or p not in base._columns
                or p > config.recompute_prime_limit
                for p in small
            )
        ):
            raise ValueError("invalid omitted factor set")
        record = restore_factors(
            retained, base, poly, position, config, spill, budget, nested=True
        )
        if type(count) is not int or count != len(small) + len(record):
            raise ValueError("recomputed factor count differs")
        if any(p <= config.recompute_prime_limit for p, _ in record):
            raise ValueError("recomputed factors overlap retained factors")
        return RecomputedFactors(poly, position, small, record, count)
    if set(payload) == {"cold"}:
        if spill is None:
            raise ValueError("cold record without an owned spill store")
        reference, count = payload["cold"]
        if type(count) is not int or not 0 <= count <= len(base.entries):
            raise ValueError("cold factor count exceeds base")
        record = ColdFactors(spill, reference, base.primes, count)
        # Canonical decoding and complete arithmetic validation follow below.
        tuple(record)
        return record
    pairs = import_factors(payload, base, budget)
    return tuple(tuple(pair) for pair in pairs)


def restore_store(payload, collector, budget):
    config, base = collector.config, collector.factor_base
    keys = {
        "extended_store",
        "atoms",
        "rows",
        "row_order",
        "kernel",
        "reporting",
        "split_calls",
        "endpoints",
        "memo",
        "canonical",
        "spill",
        "conversion_events",
        "conversion_attempts",
    }
    if (
        not isinstance(payload, dict)
        or set(payload) != keys
        or payload["extended_store"] != 1
    ):
        raise ValueError("invalid expanded store identity")
    if collector._atoms or collector._rows:
        raise ValueError("expanded restore requires an empty collector")
    atoms, rows = payload["atoms"], payload["rows"]
    if (
        not isinstance(atoms, list)
        or not isinstance(rows, list)
        or (
            len(atoms) > config.atom_limit or len(rows) > config.relation_limit
        )
    ):
        raise ValueError("expanded checkpoint exceeds its retained caps")
    if payload["row_order"] != list(range(len(rows))):
        raise ValueError("invalid expanded row order")

    # Encoded input, decoded pairs and published records coexist. Reserve the
    # largest allowed factor payload per atom before opening or decoding data.
    def factor_count(record, nested=False):
        if not isinstance(record, dict):
            raise ValueError("invalid expanded factors")
        if set(record) == {"recompute"} and not nested:
            small, count, retained = record["recompute"]
            if not isinstance(small, list) or count != len(
                small
            ) + factor_count(retained, True):
                raise ValueError("invalid expanded recomputation count")
        elif set(record) == {"cold"}:
            _, count = record["cold"]
        elif set(record) == {"index-v2", "count"}:
            count = record["count"]
            if type(count) is not int or not 0 <= count <= len(base.entries):
                raise ValueError("expanded factor count exceeds base")
            if not isinstance(record["index-v2"], str) or len(
                record["index-v2"]
            ) > 4 * ((10 * count + 2) // 3):
                raise ValueError("expanded encoding exceeds declared count")
        else:
            raise ValueError("invalid expanded factor shape")
        if type(count) is not int or not 0 <= count <= len(base.entries):
            raise ValueError("expanded factor count exceeds base")
        return count

    scratch = 65536
    for record in atoms:
        if not isinstance(record, list) or len(record) != 6:
            raise ValueError("invalid expanded atom shape")
        scratch += 32768 + 512 * factor_count(record[3])
    if collector._workspace + scratch > config.memory_limit:
        raise MemoryError("expanded restore coexistence exceeds memory cap")
    spill = None
    if config.spill_bytes:
        spill = SpillStore(
            config.spill_directory,
            config.spill_bytes,
            config.spill_cache_bytes,
            budget,
        )
        spill.restore_manifest(payload["spill"])
    elif payload["spill"] is not None:
        raise ValueError("unexpected production spill")
    arena = RecordArena(base, config=config, budget=budget)
    arena.spill = spill
    store, sizes, rebuilt_rows = {}, {}, []
    kernel = type(collector._graph)()
    checked = type(collector._checked)(base, config)
    reserve = collector._workspace
    events = payload["conversion_events"]
    attempts = payload["conversion_attempts"]
    if (
        not isinstance(events, list)
        or len(events) > config.opaque_conversions
        or (
            type(attempts) is not int
            or not len(events) <= attempts <= config.opaque_conversions
        )
    ):
        raise ValueError("invalid opaque conversion allowance")
    previous_prefix = 0
    for event in events:
        if (
            not isinstance(event, list)
            or len(event) != 3
            or (
                type(event[0]) is not int
                or not previous_prefix <= event[0] <= len(atoms)
                or not isinstance(event[2], list)
            )
        ):
            raise ValueError("invalid opaque conversion event")
        previous_prefix = event[0]
    event_index = 0

    def convert_prefix():
        nonlocal event_index, kernel, reserve, rebuilt_rows
        from .opaque_conversion import convert

        while event_index < len(events) and events[event_index][0] == len(
            store
        ):
            _, residual, primes = events[event_index]
            staged = SimpleNamespace(
                config=config,
                _graph=kernel,
                _checked=checked,
                _atoms=store,
                _arena=arena,
                _rows=rebuilt_rows,
                _combined=rebuilt_rows,
                _workspace=reserve,
                _conversion_attempts=event_index,
                _conversion_events=events[:event_index],
            )
            old_count = len(rebuilt_rows)
            if not convert(
                staged, residual, tuple(primes), budget, replay=True
            ):
                raise ValueError("opaque event has no retained source")
            kernel, reserve, rebuilt_rows = (
                staged._graph,
                staged._workspace,
                staged._rows,
            )
            for index in range(old_count, len(rebuilt_rows)):
                expected = sorted(
                    arena.slots[i] for i in rebuilt_rows[index].atom_ids
                )
                if index >= len(rows) or rows[index] != expected:
                    raise ValueError(
                        "opaque converted row differs from replay"
                    )
            event_index += 1

    for record in atoms:
        convert_prefix()
        if not isinstance(record, list) or len(record) != 6:
            raise ValueError("invalid expanded atom shape")
        coefficients, position, sign, factors, primes, opaque = record
        poly = Polynomial(base.n, base.multiplier, *coefficients)
        decoded = restore_factors(
            factors, base, poly, position, config, spill, budget
        )
        source = (
            ResidualAtom(
                poly, position, sign, (), primes=tuple(primes), opaque=opaque
            )
            if primes or opaque is not None
            else AtomicRelation(poly, position, sign, ())
        )
        cls = (
            IndexedResidualAtom
            if isinstance(source, ResidualAtom)
            else IndexedAtomicRelation
        )
        atom = cls(source, decoded)
        checked.atom(atom, budget)
        for token in atom.residual_primes + (atom.square_correction,):
            if checked.unit(token, budget) != 1:
                raise ValueError("nonunit expanded checkpoint residual")
        if atom.relation_id in store:
            raise ValueError("duplicate expanded atom")
        handle = arena.next_slot
        plan = kernel.plan_atom(atom, handle, config, budget)
        from .extended_collector import atom_reserve

        cost = atom_reserve(atom, config, base)
        if plan.cycle is not None:
            leaves = tuple(
                atom if i == handle else arena.atoms[i] for i in plan.cycle
            )
            row = checked.combine(leaves, budget, config.memory_limit).relation
            if row is None:
                raise ValueError("checkpoint combination exposed a nonunit")
            expected = list(plan.cycle)
            if (
                len(rebuilt_rows) >= len(rows)
                or rows[len(rebuilt_rows)] != expected
            ):
                raise ValueError(
                    "expanded cycle differs from exact incidence replay"
                )
            rebuilt_rows.append(row)
            cost += 8192 + 512 * sum(len(a.exponents) for a in leaves)
        reserve += cost
        if reserve + scratch > config.memory_limit:
            raise MemoryError("expanded restored ownership exceeds memory cap")
        arena.commit(atom, handle)
        kernel.commit_atom(plan, atom, config.cycle_reduction_choices)
        store[atom.relation_id], sizes[atom.relation_id] = atom, cost
    convert_prefix()
    if len(rows) != len(rebuilt_rows) or kernel.state() != payload["kernel"]:
        raise ValueError("expanded kernel state differs from checked replay")
    for row in rebuilt_rows:
        checked._record(row, atoms=tuple(store[i] for i in row.atom_ids))
    calls = payload["split_calls"]
    if type(calls) is not int or not 0 <= calls <= config.split_limit:
        raise ValueError("invalid expanded split charges")
    endpoints = type(collector._endpoints)()
    saved = payload["endpoints"]
    if (
        not isinstance(saved, list)
        or len(saved) > endpoints.maximum
        or len(set(saved)) != len(saved)
    ):
        raise ValueError("invalid expanded endpoint cache")
    for prime in saved:
        if type(prime) is not int or not 2 <= prime <= max(
            config.residual_bound, config.large_prime_bound
        ):
            raise ValueError("invalid expanded endpoint bound")
        if (
            not checked.prime(prime, budget)
            or checked.unit(prime, budget) != 1
        ):
            raise ValueError("invalid expanded endpoint proof")
    endpoints.add(saved, budget)
    memo = type(collector._split_memo)(config.split_cache_entries)
    if (
        not isinstance(payload["memo"], list)
        or len(payload["memo"]) > memo.capacity
    ):
        raise ValueError("invalid expanded split memo")
    for key, primes, reason in payload["memo"]:
        if len(key) != 5 or key[1:] != [
            config.large_prime_bound,
            config.split_policy,
            config.residual_arity,
            config.split_dispatch,
        ]:
            raise ValueError("split memo uses a different policy")
        product = 1
        if not 1 <= len(primes) <= config.residual_arity or primes != sorted(
            primes
        ):
            raise ValueError("invalid split memo arity")
        for prime in primes:
            if (
                type(prime) is not int
                or not 2 <= prime <= config.large_prime_bound
                or not checked.prime(prime, budget)
            ):
                raise ValueError("invalid split memo endpoint")
            product *= prime
        if (
            type(key[0]) is not int
            or product != key[0]
            or reason not in ("double", "square")
        ):
            raise ValueError("split memo reconstruction differs")
        memo.remember(tuple(key), (tuple(primes), reason))
    canonical = payload["canonical"]
    if (
        not isinstance(canonical, list)
        or len(canonical) > len(atoms)
        or len(set(canonical)) != len(canonical)
    ):
        raise ValueError("invalid expanded canonical index")
    for index in canonical:
        if type(index) is not int or not 0 <= index < len(atoms):
            raise ValueError("canonical handle outside retained atoms")
        atom = arena.atoms[index]
        key, divisor = checked.congruence(atom, budget)
        if divisor is not None or key in checked.canonical:
            raise ValueError("invalid canonical representative")
        checked.remember(atom, key=key)
    if len(checked.canonical) != len(canonical):
        raise MemoryError("saved canonical cache exceeds its allowance")
    original_reporting = collector._reporting
    reporting = type(original_reporting)(
        original_reporting.lower,
        original_reporting.upper,
        original_reporting.exposure,
        original_reporting.epochs,
    )
    reporting.restore(payload["reporting"])
    published_rows = list(rebuilt_rows)
    budget.consume(len(atoms) + len(rows) + 1)
    collector._atoms, collector._atom_bytes = store, sizes
    collector._arena, collector._graph, collector._checked = (
        arena,
        kernel,
        checked,
    )
    collector._endpoints, collector._split_memo = endpoints, memo
    collector._combined = rebuilt_rows
    collector._rows = published_rows
    collector._full = []
    collector._split_calls = calls
    collector._conversion_events = events
    collector._conversion_attempts = attempts
    collector._reporting = reporting
    collector._candidate_bound = collector._reporting.bound
    collector._conversion_bytes = (
        reserve - collector._workspace - sum(sizes.values())
    )
    collector._workspace = reserve
    cache = getattr(collector, "_preparation_cache", None)
    if cache is not None:
        cache.checked = checked
