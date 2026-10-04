"""Bounded row-oriented GF(2) filtering with original-row provenance."""

from dataclasses import dataclass
from heapq import heapify, heappop, heappush

from .. import utils
from ..budget import Budget
from .factor_base import DEFAULT_MEMORY_BYTES

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


def filter_matrix(
    rows,
    *,
    weight_two=False,
    budget=None,
    memory_bytes=DEFAULT_MEMORY_BYTES,
):
    """Remove singleton constraints and optionally eliminate weight two.

    A column in exactly one row forces that row out of every kernel. A
    column in exactly two rows forces equal selection of those rows, so
    replace them by their XOR and XOR their original-row masks. All kernels
    of the surviving matrix lift exactly to kernels of the original matrix.
    Refusal raises before publishing a result; caller input is never mutated.
    """
    if not isinstance(rows, tuple):
        raise TypeError("rows must be an immutable tuple")
    if len(rows) > MAX_MATRIX_ROWS:
        raise ValueError("matrix row cap exceeded")
    for row in rows:
        utils.require_integer(row, "row", 0)
        if row.bit_length() > MAX_MATRIX_COLUMNS:
            raise ValueError("matrix column cap exceeded")
    utils.require_integer(memory_bytes, "memory_bytes", 0)
    union = 0
    for row in rows:
        union |= row
    original_columns = union.bit_length()
    columns = original_columns
    compact = columns > 2 * union.bit_count()
    if compact:
        columns = union.bit_count()
    reserve = matrix_workspace(len(rows), columns)
    if compact:
        # The original wide bitsets coexist with remapped rows and the
        # index map. Reserve them explicitly, without dense wide fill-in.
        reserve += 256 * columns + 16 * ((original_columns + 7) // 8)
        reserve += sum(128 + 8 * ((row.bit_length() + 7) // 8) for row in rows)
    if reserve > memory_bytes:
        raise MemoryError("matrix fill-in/provenance exceeds memory_bytes")
    budget = budget if budget is not None else Budget()
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
            remapped, bits = 0, row
            while bits:
                bit = bits & -bits
                remapped |= 1 << indices[bit.bit_length() - 1]
                bits ^= bit
            working_rows.append(remapped)
        del indices
    active = {
        index: (row, 1 << index) for index, row in enumerate(working_rows)
    }
    zero, singletons, merges, rounds = [], 0, 0, 0
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
    for index, (row, _) in active.items():
        budget.consume(0)
        bits = row
        while bits:
            bit = bits & -bits
            column = bit.bit_length() - 1
            incidence[column] = incidence.get(column, 0) | (1 << index)
            bits ^= bit
    single_columns = {
        column for column, mask in incidence.items() if mask.bit_count() == 1
    }
    pair_columns = [
        column
        for column, mask in incidence.items()
        if weight_two and mask.bit_count() == 2
    ]
    heapify(pair_columns)
    queued_pairs = set(pair_columns)
    pending_zeros = [index for index, (row, _) in active.items() if not row]
    nonzeros = input_nonzeros

    def toggle(index, bits):
        """Update affected columns; queues contain at most one copy each."""
        row_bit = 1 << index
        while bits:
            bit = bits & -bits
            column = bit.bit_length() - 1
            budget.consume(word_cost)
            mask = incidence.get(column, 0) ^ row_bit
            if mask:
                incidence[column] = mask
            else:
                incidence.pop(column, None)
            count = mask.bit_count()
            if count == 1:
                single_columns.add(column)
            else:
                single_columns.discard(column)
            if weight_two and count == 2 and column not in queued_pairs:
                heappush(pair_columns, column)
                queued_pairs.add(column)
            bits ^= bit

    while active:
        budget.consume(0)
        rounds += 1
        for index in pending_zeros:
            zero.append(active.pop(index)[1])
        pending_zeros = []
        forced = {
            incidence[column].bit_length() - 1 for column in single_columns
        }
        if forced:
            for index in forced:
                row, _ = active.pop(index)
                toggle(index, row)
                nonzeros -= row.bit_count()
            singletons += len(forced)
            continue
        pair = None
        while pair_columns:
            column = heappop(pair_columns)
            queued_pairs.remove(column)
            mask = incidence.get(column, 0)
            if mask.bit_count() == 2:
                first_bit = mask & -mask
                pair = (
                    first_bit.bit_length() - 1,
                    (mask ^ first_bit).bit_length() - 1,
                )
                break
        if pair is None:
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
        self.budget = budget if budget is not None else Budget()
        self.dependencies = list(matrix.zero_dependencies)
        for mask in self.dependencies:
            verify_dependency(mask, matrix.original_rows)
        self.pivots = {}
        self.next_row = 0
        self.pending = None
        self.xors = 0
        self.peak_nonzeros = 0

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
            self.pivots[bit] = row, mask
            self.peak_nonzeros = max(
                self.peak_nonzeros,
                sum(value.bit_count() for value, _ in self.pivots.values()),
            )
            self.pending = None
            self.next_row += 1
        else:
            self.pending = row ^ previous[0], mask ^ previous[1]
            self.xors += 1

    def run(self):
        """Finish elimination or raise with resumable state."""
        while self.next_row < len(self.matrix.rows):
            self.step()
        return tuple(self.dependencies)
