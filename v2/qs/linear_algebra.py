"""Bounded row-oriented GF(2) filtering with original-row provenance."""

from array import array
from dataclasses import dataclass
from heapq import heapify, heappop, heappush

from ..common import arithmetic, utils
from ..execution.resources import removed_arguments, resolve_budget

MAX_MATRIX_ROWS = 65536
MAX_MATRIX_COLUMNS = 100_001


@dataclass(frozen=True)
class FilteredMatrix:
    """Prime/sign columns in row bitsets; masks lift to original rows.

    zero_dependencies are kernels discovered during filtering, never dropped.
    Exact duplicate payloads must be removed before calling filter_matrix;
    equal parity alone is not a reason to discard a distinct relation.
    """

    original_rows: tuple
    rows: tuple
    masks: tuple
    zero_dependencies: tuple
    stats: dict
    workspace_bytes: int


def matrix_workspace(row_count, column_count):
    """Reserve worst-case fill-in, lift masks, pivots and temporary copies."""
    return (
        32768
        + column_count * 512
        + row_count
        * (512 + 8 * ((row_count + 7) // 8 + (column_count + 7) // 8))
    )


@removed_arguments(("max_atoms", "max_rows", "memory_bytes"))
def filter_matrix(
    rows,
    *,
    weight_two=False,
    budget=None,
    merge_policy="reference",
    light_rows=False,
    surplus=32,
    merge_weight=2,
    merge_fill=0,
    sparse_incidence=False,
    panel_bits=0,
    survivor_mask=None,
    _incremental=None,
    _component_contractions=True,
):
    """Remove singleton constraints and optionally eliminate weight two.

    A column in exactly one row forces that row out of every kernel. A
    column in exactly two rows forces equal selection of those rows, so
    replace them by their XOR and XOR their original-row masks. All kernels
    of the surviving matrix lift exactly to kernels of the original matrix.
    Refusal raises before publishing a result; caller input is never mutated.
    """
    if merge_policy not in ("reference", "fill-v2"):
        raise ValueError("unknown weight-two merge policy")
    for name, value, maximum in (
        ("merge_weight", merge_weight, 8),
        ("merge_fill", merge_fill, 4096),
        ("panel_bits", panel_bits, 8),
    ):
        utils.require_integer(value, name, 0)
        if value > maximum:
            raise ValueError(name + " exceeds its finite cap")
    if merge_weight < 2 or type(sparse_incidence) is not bool:
        raise ValueError("invalid structured filtering policy")
    if type(light_rows) is not bool:
        raise TypeError("light_rows must be Boolean")
    utils.require_integer(surplus, "matrix surplus", 0)
    if surplus > 4096:
        raise ValueError("matrix surplus exceeds its bound")
    if not isinstance(rows, tuple):
        raise TypeError("rows must be an immutable tuple")
    if survivor_mask is not None and (
        type(survivor_mask) is not int
        or survivor_mask < 0
        or survivor_mask.bit_length() > len(rows)
    ):
        raise ValueError("invalid incremental singleton mask")
    for row in rows:
        utils.require_integer(row, "row", 0)
    backend = arithmetic.backend_for(
        next((row for row in rows if arithmetic.is_mpz(row)), 0)
    )
    zero_bit, one_bit = backend.integer(0), backend.integer(1)
    union = zero_bit
    for row in rows:
        union |= row
    original_columns = union.bit_length()
    active_columns = union.bit_count()
    columns = original_columns
    compact = columns > 2 * union.bit_count()
    if compact:
        columns = union.bit_count()
    reserve = matrix_workspace(len(rows), columns)
    if panel_bits:
        # Pivot presence and coexisting masks used for ordered panel lookup.
        reserve += 512 + 32 * ((columns + 7) // 8)
        reserve += (1 << panel_bits) * (
            512 + 8 * ((len(rows) + columns + 7) // 8)
        )
    if compact:
        # The original wide bitsets coexist with remapped rows and the
        # index map. Reserve them explicitly, without dense wide fill-in.
        reserve += 256 * columns + 16 * ((original_columns + 7) // 8)
        reserve += sum(128 + 8 * ((row.bit_length() + 7) // 8) for row in rows)

    budget = resolve_budget(budget)
    if _incremental is not None:
        survivor_mask = _incremental.observe(rows, budget)
    working_rows = rows
    if compact:
        nonzeros = sum(row.bit_count() for row in rows)
        source_words = (original_columns + 63) // 64
        target_words = (columns + 63) // 64
        budget.consume(
            (source_words + 1) * (columns + nonzeros)
            + target_words * nonzeros
            + len(rows)
        )
        indices, bits = {}, union

        while bits:
            bit = bits & -bits
            indices[bit.bit_length() - 1] = len(indices)
            bits ^= bit
            if len(indices) % 64 == 0:
                budget.consume(0)

        working_rows = []

        for row in rows:
            budget.consume(0)
            remapped, bits = zero_bit, row
            while bits:
                bit = bits & -bits
                remapped |= one_bit << indices[bit.bit_length() - 1]
                bits ^= bit
            working_rows.append(remapped)

    active = {
        index: (row, one_bit << index)
        for index, row in enumerate(working_rows)
        if survivor_mask is None or survivor_mask & (1 << index)
    }
    selected_count = len(rows)
    if light_rows and len(active) > active_columns + surplus:
        budget.consume(
            len(active) * (1 + len(active).bit_length() + (columns + 63) // 64)
        )
        kept_count = max(
            sum(not row for row, _ in active.values()),
            active_columns + surplus,
        )
        selected_count = len(rows) - len(active) + kept_count
        selected = sorted(active, key=lambda i: (active[i][0].bit_count(), i))[
            :kept_count
        ]
        active = {i: active[i] for i in selected}
    zero, singletons, merges, rounds = [], len(rows) - len(active), 0, 0
    incremental_core = False
    if (
        _incremental is not None
        and _component_contractions
        and weight_two
        and selected_count == len(rows)
    ):
        core = _incremental.weight_two_core(budget, merge_policy)
        if core is not None:
            incremental_core = True
            pairs, zeros, singletons, merges = core
            active = {}
            for index, (row, mask) in enumerate(pairs):
                budget.consume(row.bit_count() + mask.bit_count() + 1)
                if compact:
                    mapped, bits = zero_bit, row
                    while bits:
                        bit = bits & -bits
                        mapped |= one_bit << indices[bit.bit_length() - 1]
                        bits ^= bit
                    row = mapped
                active[index] = row, mask
            zero = list(zeros)
    if compact:
        del indices
    input_nonzeros = sum(row.bit_count() for row in rows)
    peak_nonzeros = input_nonzeros
    # Column incidence uses row-index bitsets. The existing worst-case
    # matrix reservation covers these masks and the bounded queues.
    # Charge initial incidence by actual nonzeros and bounded 64-bit words.
    # Each nonzero removes a column bit and updates a row-index bitset.
    # Dense worst-case storage remains reserved; sparse work need not be dense.
    word_cost = 1 + (columns + 63) // 64 + (len(rows) + 63) // 64
    budget.consume(len(rows) + input_nonzeros * word_cost)
    incidence = {}

    if sparse_incidence:
        # Lighter-row selection orders by payload weight. Sparse segments
        # must instead start in row-number order for later bisect toggles.
        budget.consume(len(active) * (len(active).bit_length() + 1))
        initial_rows = sorted(active.items())
    else:
        initial_rows = active.items()
    for index, (row, _) in initial_rows:
        budget.consume(0)
        bits = row
        while bits:
            bit = bits & -bits
            column = bit.bit_length() - 1
            if sparse_incidence:
                incidence.setdefault(column, array("I")).append(index)
            else:
                incidence[column] = incidence.get(column, zero_bit) | (
                    one_bit << index
                )
            bits ^= bit

    def count(mask):
        return len(mask) if sparse_incidence else mask.bit_count()

    def members(mask):
        if sparse_incidence:
            return tuple(mask)
        result = []
        while mask:
            bit = mask & -mask
            result.append(bit.bit_length() - 1)
            mask ^= bit
        return tuple(result)

    empty = array("I") if sparse_incidence else zero_bit
    single_columns = {
        column for column, mask in incidence.items() if count(mask) == 1
    }
    pair_columns = [
        column
        for column, mask in incidence.items()
        if weight_two and count(mask) == 2
    ]
    heapify(pair_columns)
    queued_pairs = set(pair_columns)
    pending_zeros = [index for index, (row, _) in active.items() if not row]
    nonzeros = sum(row.bit_count() for row, _ in active.values())

    def toggle(index, bits):
        """Update affected columns; queues contain at most one copy each."""
        row_bit = one_bit << index

        while bits:
            bit = bits & -bits
            column = bit.bit_length() - 1
            budget.consume(word_cost)
            if sparse_incidence:
                from .structured_matrix import incidence_add

                previous = incidence.get(column, empty)
                budget.consume(len(previous) + 1)
                mask = incidence_add(previous, index)
            else:
                mask = incidence.get(column, zero_bit) ^ row_bit
            if mask:
                incidence[column] = mask
            else:
                incidence.pop(column, None)
            degree = count(mask)
            if degree == 1:
                single_columns.add(column)
            else:
                single_columns.discard(column)
            if weight_two and degree == 2 and column not in queued_pairs:
                heappush(pair_columns, column)
                queued_pairs.add(column)
            bits ^= bit

    structured_merges = 0
    while active:
        budget.consume(0)
        rounds += 1
        for index in pending_zeros:
            zero.append(active.pop(index)[1])
        pending_zeros = []
        forced = {members(incidence[column])[0] for column in single_columns}
        if forced:
            for index in forced:
                row, _ = active.pop(index)
                toggle(index, row)
                nonzeros -= row.bit_count()
            singletons += len(forced)
            continue

        pair = None

        proposals = []
        while pair_columns:
            column = heappop(pair_columns)
            queued_pairs.remove(column)
            mask = incidence.get(column, empty)
            if count(mask) == 2:
                candidate = members(mask)
                if merge_policy == "reference":
                    pair = candidate
                    break
                budget.consume(3 * word_cost)
                left, right = (active[i][0] for i in candidate)
                weight = (left ^ right).bit_count()
                proposals.append(
                    (
                        weight - left.bit_count() - right.bit_count(),
                        weight,
                        column,
                        candidate,
                    )
                )
                if len(proposals) == 8:
                    break

        if proposals:
            chosen = min(proposals)
            pair = chosen[3]
            for proposal in proposals:
                if proposal is not chosen:
                    heappush(pair_columns, proposal[2])
                    queued_pairs.add(proposal[2])

        if pair is None:
            from .structured_matrix import spanning_merges

            proposals = []
            if merge_weight > 2:
                for column in sorted(incidence):
                    budget.consume(1)
                    degree = count(incidence[column])
                    if not 3 <= degree <= merge_weight:
                        continue
                    indices = members(incidence[column])
                    edges = spanning_merges(indices, active, budget)
                    old_weight = sum(active[i][0].bit_count() for i in indices)
                    delta = sum(r.bit_count() for r, _ in edges) - old_weight
                    if delta <= merge_fill:
                        proposals.append((delta, column, indices, edges))
                    if len(proposals) >= 8:
                        break
                if proposals:
                    delta, _, indices, edges = min(proposals)
                    for index in indices:
                        row, _ = active.pop(index)
                        toggle(index, row)
                    for index, (row, mask) in zip(indices, edges):
                        active[index] = row, mask
                        toggle(index, row)
                        if not row:
                            pending_zeros.append(index)
                    nonzeros += delta
                    peak_nonzeros = max(peak_nonzeros, nonzeros)
                    structured_merges += 1
                    continue
            break
        first, second = pair
        left, left_mask = active[first]
        right, right_mask = active.pop(second)
        budget.consume(word_cost)
        merged = left ^ right
        # Removing the right row and XORing it into the left row changes
        # precisely the right row's columns, twice.
        toggle(second, right)
        toggle(first, right)
        active[first] = merged, left_mask ^ right_mask
        nonzeros += merged.bit_count() - left.bit_count() - right.bit_count()
        merges += 1
        peak_nonzeros = max(peak_nonzeros, nonzeros)
        if not merged:
            pending_zeros.append(first)

    remaining = tuple(active.values())
    output_rows = tuple(row for row, _ in remaining)
    union = 0
    for row in output_rows:
        union |= row
    stats = {
        "input_rows": len(rows),
        "input_columns": original_columns,
        "working_columns": columns,
        "input_nonzeros": input_nonzeros,
        "output_rows": len(output_rows),
        "output_columns": union.bit_count(),
        "output_nonzeros": sum(row.bit_count() for row in output_rows),
        "peak_nonzeros": peak_nonzeros,
        "singletons_removed": singletons,
        "weight_two_merges": merges,
        "zero_dependencies": len(zero),
        "rounds": rounds,
    }
    if merge_policy != "reference" or light_rows:
        stats.update(
            input_active_columns=active_columns, selected_rows=selected_count
        )
    if merge_weight > 2 or sparse_incidence or panel_bits:
        stats.update(
            structured_merges=structured_merges,
            panel_bits=panel_bits,
            sparse_incidence=sparse_incidence,
        )
    if _incremental is not None:
        stats.update(
            incremental_survivors=survivor_mask is not None,
            incremental_core=incremental_core,
        )
    return FilteredMatrix(
        rows,
        output_rows,
        tuple(mask for _, mask in remaining),
        tuple(zero),
        stats,
        reserve,
    )


def verify_dependency(mask, original_rows):
    """Independently XOR selected original rows; reject empty/invalid masks."""
    utils.require_integer(mask, "dependency mask", 1)
    if mask.bit_length() > len(original_rows):
        raise ValueError("dependency mask refers outside original rows")
    parity = 0
    bits = mask
    while bits:
        bit = bits & -bits
        parity ^= original_rows[bit.bit_length() - 1]
        bits ^= bit
    if parity:
        raise ValueError("dependency is not in the original-row kernel")
    return True


class DependencySolver:
    """Resumable bitset elimination; refusal retains the pending row.

    lowest/highest select deterministic pivots. Every XOR is charged before
    state mutation. The matrix reservation already includes dense fill-in.
    Dependencies refer to original rows, including filtering transformations.
    """

    def __init__(self, matrix, *, pivot="highest", budget=None):
        """Retain a filtered matrix and verify preexisting zero kernels."""
        if pivot not in ("highest", "lowest"):
            raise ValueError("unknown pivot strategy")
        self.matrix, self.pivot = matrix, pivot
        self.budget = resolve_budget(budget)
        self.dependencies = list(matrix.zero_dependencies)
        for mask in self.dependencies:
            verify_dependency(mask, matrix.original_rows)
        self.pivots = {}
        self._panel_pivots = 0
        self.next_row = 0
        self.pending = None
        self.xors = 0
        self.peak_nonzeros = 0
        self.pivot_nonzeros = 0

    def step(self):
        """Commit one elimination action; keep the pending row on refusal."""
        if self.next_row >= len(self.matrix.rows):
            return
        if self.pending is None:
            self.pending = (
                self.matrix.rows[self.next_row],
                self.matrix.masks[self.next_row],
            )

        row, mask = self.pending
        self.budget.consume(row.bit_length() + mask.bit_length() + 1)
        if row == 0:
            verify_dependency(mask, self.matrix.original_rows)
            self.dependencies.append(mask)
            self.pending = None
            self.next_row += 1
            return

        bit = (
            1 << (row.bit_length() - 1)
            if self.pivot == "highest"
            else row & -row
        )
        previous = self.pivots.get(bit)
        if previous is None:
            if self.matrix.stats.get("panel_bits", 0):
                self.budget.consume(self.matrix.stats["working_columns"] + 1)
                self._panel_pivots |= bit
            self.pivots[bit] = row, mask
            self.pivot_nonzeros += row.bit_count()
            self.peak_nonzeros = max(self.peak_nonzeros, self.pivot_nonzeros)
            self.pending = None
            self.next_row += 1
        else:
            # Carry the XOR through provenance so a reduced zero row still
            # selects a kernel of the original, unfiltered matrix.
            from .structured_matrix import panel_reduce

            reduced = panel_reduce(self, row, mask)
            self.pending = (
                reduced
                if reduced is not None
                else (row ^ previous[0], mask ^ previous[1])
            )
            self.xors += 1

    def run(self):
        """Finish elimination or raise with resumable state."""
        while self.next_row < len(self.matrix.rows):
            self.step()
        return tuple(self.dependencies)
