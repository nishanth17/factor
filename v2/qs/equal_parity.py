"""Immediate checked two-row dependencies without deleting equal parity."""

from .extraction import (
    PreparedRelations,
    _prepare_relations,
    extract_dependency,
)
from .relations import AtomicRelation, parity_bits


def _weight(row):
    return (
        1 if isinstance(row, AtomicRelation) else len(row.atom_ids),
        len(row.exponents),
    )


class EqualParityIndex:
    """Bounded representatives; failed/unchecked pairs remain ordinary rows."""

    def __init__(self, memory_bytes):
        self.memory_bytes = memory_bytes
        self.used = self.checked_count = 0
        self.representatives = {}

    def clear(self):
        self.representatives.clear()
        self.used = self.checked_count = 0

    def _remember(self, bits, index, rows):
        if not bits:
            return
        old = self.representatives.get(bits)
        if old is not None:
            if _weight(rows[index]) < _weight(rows[old]):
                self.representatives[bits] = index
            return
        reserve = 512 + 8 * ((bits.bit_length() + 7) // 8)
        if self.used + reserve <= self.memory_bytes:
            self.representatives[bits] = index
            self.used += reserve

    def rebuild(self, rows, count, base, budget):
        """Rebuild derived keys without repeating checked dependencies."""
        if type(count) is not int or not 0 <= count <= len(rows):
            raise ValueError("equal-parity cursor exceeds retained rows")
        for index in range(count):
            budget.consume(len(rows[index].exponents) + 1)
            self._remember(parity_bits(rows[index], base), index, rows)
        self.checked_count = count

    def probe(self, rows, collector, budget, stats):
        """Complete one checked pair before advancing its resumable cursor."""
        while self.checked_count < len(rows):
            index = self.checked_count
            row = rows[index]
            base = collector.factor_base
            budget.consume(
                (len(row.exponents) + 1) * (1 + len(base.entries) // 64)
            )
            bits = (
                collector._checked.row(
                    row,
                    collector._atoms,
                    budget,
                    collector.config.memory_limit,
                )
                if collector._checked is not None
                else parity_bits(row, base)
            )
            old = self.representatives.get(bits) if bits else None
            divisor = None
            if old is not None:
                selected = (rows[old], row)
                try:
                    if collector._checked is not None:
                        scratch = 65536 + 512 * sum(
                            len(r.exponents) for r in selected
                        )
                        if (
                            collector._workspace + scratch
                            > collector.config.memory_limit
                        ):
                            raise MemoryError(
                                "checked pair extraction scratch"
                            )
                        parities = tuple(
                            collector._checked.row(
                                r,
                                collector._atoms,
                                budget,
                                collector.config.memory_limit,
                            )
                            for r in selected
                        )
                        prepared = PreparedRelations(
                            selected,
                            parities,
                            (old, index),
                            (),
                            base,
                            (),
                            scratch,
                        )
                        result = extract_dependency(
                            prepared,
                            3,
                            budget=budget,
                            product_policy=collector.config.extraction_policy,
                        )
                        divisor = result.divisor
                    else:
                        divisor = self._checked_pair(
                            selected, base, collector, budget
                        )
                except MemoryError:
                    stats["pair_scratch_refusals"] = (
                        stats.get("pair_scratch_refusals", 0) + 1
                    )
                else:
                    stats["pair_dependencies"] = (
                        stats.get("pair_dependencies", 0) + 1
                    )
                    stats["pair_trivial"] = stats.get("pair_trivial", 0) + (
                        divisor is None
                    )
            self._remember(bits, index, rows)
            self.checked_count += 1
            if divisor is not None:
                return divisor
        return None

    @staticmethod
    def _checked_pair(selected, base, collector, budget):
        # Two cycles own at most 512 original IDs. Reserve the
        # temporary set/map before constructing either container;
        # preparation accounts for its other coexisting scratch.
        scratch = 65536 + 512 * sum(_weight(r)[0] for r in selected)
        if collector._workspace + scratch > collector.config.memory_limit:
            raise MemoryError("equal-parity provenance scratch")
        identities = set()
        for item in selected:
            budget.consume(_weight(item)[0])
            identities.update(
                (item.relation_id,)
                if isinstance(item, AtomicRelation)
                else item.atom_ids
            )
        store = {key: collector._atoms[key] for key in identities}
        prepared = _prepare_relations(
            selected,
            base,
            store,
            budget=budget,
            memory_bytes=collector.config.memory_limit,
            retained_workspace_bytes=collector._workspace + scratch,
            max_rows=2,
            max_atoms=512,
        )
        if len(prepared.rows) == 2:
            result = extract_dependency(prepared, 3, budget=budget)
            divisor = result.divisor
        return divisor if len(prepared.rows) == 2 else None
