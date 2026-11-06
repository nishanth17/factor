"""Bounded sparse incidence, tree elimination and panel-table kernels."""

from array import array
from bisect import bisect_left


def incidence_add(mask, index):
    """Toggle a row in a compact sorted incidence segment."""
    result = array("I", mask)
    position = bisect_left(result, index)
    if position < len(result) and result[position] == index:
        result.pop(position)
    else:
        result.insert(position, index)
    return result


def spanning_merges(indices, active, budget):
    """Minimum spanning tree of exact pair-XOR weights on at most eight rows.

    The k-1 edge vectors form a basis of even selections of the eliminated
    column. Their original-row masks provide exact back-substitution.
    """
    if not 3 <= len(indices) <= 8:
        raise ValueError("structured merge degree exceeds its cap")
    seen, remaining, edges = {indices[0]}, set(indices[1:]), []
    while remaining:
        choices = []
        for a in sorted(seen):
            for b in sorted(remaining):
                left, right = active[a], active[b]
                budget.consume(
                    left[0].bit_length() + right[0].bit_length() + 1
                )
                merged = left[0] ^ right[0]
                choices.append(
                    (merged.bit_count(), a, b, merged, left[1] ^ right[1])
                )
        best = min(choices)
        _, _, endpoint, row, mask = best
        edges.append((row, mask))
        remaining.remove(endpoint)
        seen.add(endpoint)
    return tuple(edges)


def panel_reduce(solver, row, mask):
    """Eliminate a bounded triangular pivot panel with one table lookup."""
    width = solver.matrix.stats.get("panel_bits", 0)
    if not width:
        return None
    leading = (
        1 << (row.bit_length() - 1)
        if solver.pivot == "highest"
        else row & -row
    )
    # A pivot key is one column bit. The incrementally maintained presence
    # mask selects ordered keys directly, without sorting or scanning pivots.
    solver.budget.consume(width * (solver.matrix.stats["working_columns"] + 1))
    bits = solver._panel_pivots
    if solver.pivot == "highest":
        bits &= (leading << 1) - 1
    else:
        bits &= -leading
    selected = []
    while bits and len(selected) < width:
        bit = (
            1 << (bits.bit_length() - 1)
            if solver.pivot == "highest"
            else bits & -bits
        )
        selected.append(bit)
        bits ^= bit
    keys = tuple(selected)
    if len(keys) < 2:
        return None
    # New pivots outside this panel cannot invalidate its immutable table.
    identity = keys
    cache = getattr(solver, "_panel_cache", None)
    if cache is None or cache[0] != identity:
        cost = (1 << len(keys)) * (
            solver.matrix.stats["working_columns"]
            + len(solver.matrix.original_rows)
            + len(keys)
            + 1
        )
        solver.budget.consume(cost)
        vectors = [(0, 0)]
        for key in keys:
            pivot_row, pivot_mask = solver.pivots[key]
            vectors.extend(
                tuple((r ^ pivot_row, m ^ pivot_mask) for r, m in vectors)
            )
        table = {}
        for value, provenance in vectors:
            projection = sum(
                (bool(value & bit) << i) for i, bit in enumerate(keys)
            )
            if projection in table:
                raise ArithmeticError("singular elimination panel")
            table[projection] = value, provenance
        cache = identity, table
        solver._panel_cache = cache
    solver.budget.consume(row.bit_length() + mask.bit_length() + len(keys))
    projection = sum((bool(row & bit) << i) for i, bit in enumerate(keys))
    value, provenance = cache[1][projection]
    return row ^ value, mask ^ provenance
