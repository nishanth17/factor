"""Bounded row-oriented GF(2) filtering with original-row provenance."""

from dataclasses import dataclass

from .. import utils
from ..budget import Budget
from .factor_base import DEFAULT_MEMORY_BYTES

MAX_MATRIX_ROWS = 4096
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
    columns = max((row.bit_length() for row in rows), default=0)
    reserve = matrix_workspace(len(rows), columns)
    if reserve > memory_bytes:
        raise MemoryError("matrix fill-in/provenance exceeds memory_bytes")
    budget = budget if budget is not None else Budget()
    active = {index: (row, 1 << index) for index, row in enumerate(rows)}
    zero, singletons, merges, rounds = [], 0, 0, 0
    input_nonzeros = sum(row.bit_count() for row in rows)
    peak_nonzeros = input_nonzeros
    while active:
        budget.consume(len(active) * (columns + len(rows) + 1))
        rounds += 1
        incidence = {}
        for index, (row, mask) in tuple(active.items()):
            if row == 0:
                zero.append(mask)
                del active[index]
                continue
            bits = row
            while bits:
                bit = bits & -bits
                column = bit.bit_length() - 1
                previous = incidence.get(column)
                if previous is None:
                    incidence[column] = [1, index, None]
                else:
                    previous[0] += 1
                    if previous[0] == 2:
                        previous[2] = index
                bits ^= bit
        forced = {
            first for count, first, _ in incidence.values() if count == 1
        }
        if forced:
            for index in forced:
                del active[index]
            singletons += len(forced)
            continue
        pair = (
            next(
                (
                    (first, second)
                    for _, (count, first, second) in sorted(incidence.items())
                    if count == 2
                ),
                None,
            )
            if weight_two
            else None
        )
        if pair is None:
            break
        first, second = pair
        left, left_mask = active[first]
        right, right_mask = active.pop(second)
        active[first] = left ^ right, left_mask ^ right_mask
        merges += 1
        peak_nonzeros = max(
            peak_nonzeros, sum(row.bit_count() for row, _ in active.values())
        )
    remaining = tuple(active.values())
    output_rows = tuple(row for row, _ in remaining)
    union = 0
    for row in output_rows:
        union |= row
    stats = {
        "input_rows": len(rows),
        "input_columns": columns,
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
    for index, row in enumerate(original_rows):
        if mask & (1 << index):
            parity ^= row
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

    def run(self):
        """Finish elimination or raise with resumable state."""
        while self.next_row < len(self.matrix.rows):
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
                continue
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
                    sum(
                        value.bit_count() for value, _ in self.pivots.values()
                    ),
                )
                self.pending = None
                self.next_row += 1
            else:
                self.pending = row ^ previous[0], mask ^ previous[1]
                self.xors += 1
        return tuple(self.dependencies)
