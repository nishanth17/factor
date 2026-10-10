"""Bounded filtering challengers, including conversion and exact lifting."""

from dataclasses import replace


def live_compaction(matrix, budget, memory_bytes):
    """Remap surviving columns and retain their original-label inverse map."""
    union = 0
    for row in matrix.rows:
        union |= row
    labels, bits = [], union
    while bits:
        bit = bits & -bits
        labels.append(bit.bit_length() - 1)
        bits ^= bit
    reserve = matrix.workspace_bytes + 32768 + 256 * len(labels)
    reserve += sum(
        128 + 8 * ((row.bit_length() + 7) // 8) for row in matrix.rows
    )
    if reserve > memory_bytes:
        raise MemoryError("live compaction coexistence exceeds memory cap")
    budget.consume(
        (union.bit_length() // 64 + 1)
        * (len(labels) + sum(r.bit_count() for r in matrix.rows))
    )
    lookup = {column: index for index, column in enumerate(labels)}
    remapped = []

    for row in matrix.rows:
        budget.consume(0)
        bits, value = row, 0
        while bits:
            bit = bits & -bits
            value |= 1 << lookup[bit.bit_length() - 1]
            bits ^= bit
        remapped.append(value)

    original_union = 0
    for row in matrix.original_rows:
        original_union |= row
    initial_labels, bits = [], original_union
    if original_union.bit_length() > 2 * original_union.bit_count():
        while bits:
            bit = bits & -bits
            initial_labels.append(bit.bit_length() - 1)
            bits ^= bit

    inverse = tuple(initial_labels[i] if initial_labels else i for i in labels)
    stats = dict(matrix.stats, inverse_columns=inverse)
    return replace(
        matrix, rows=tuple(remapped), stats=stats, workspace_bytes=reserve
    )


def history_filter(
    rows, module, budget, memory_bytes, batch_size=1, use_history=True
):
    """Evaluate disjoint pivot batches and immutable trees with deferred lift.

    This bounded prototype deliberately reserves the existing dense workspace
    plus tree storage. It provides no claim of lower solver capacity. A batch
    uses disjoint incident rows from one incidence snapshot; every merge still
    updates touched columns. Lifting is iterative and independently checked.
    """
    if (
        len(rows) > 4096
        or max((r.bit_length() for r in rows), default=0) > 4096
    ):
        raise ValueError("history prototype is limited to 4096 rows/columns")

    columns = max((r.bit_length() for r in rows), default=0)
    reserve = module.matrix_workspace(len(rows), columns)
    reserve += 65536
    if use_history:
        reserve += 1024 * (2 * len(rows))
    if reserve > memory_bytes:
        raise MemoryError("merge history coexistence exceeds memory cap")
    active = {
        index: (row, index if use_history else 1 << index)
        for index, row in enumerate(rows)
    }
    history = [(index,) for index in range(len(rows))] if use_history else []
    zeros, incidence = [], {}
    words = 1 + (columns + 63) // 64 + (len(rows) + 63) // 64
    budget.consume(len(rows) + sum(r.bit_count() for r in rows) * words)
    for index, row in enumerate(rows):
        bits = row
        while bits:
            bit = bits & -bits
            incidence[bit] = incidence.get(bit, 0) | (1 << index)
            bits ^= bit

    def toggle(index, bits):
        while bits:
            bit = bits & -bits
            budget.consume(words)
            value = incidence.get(bit, 0) ^ (1 << index)
            if value:
                incidence[bit] = value
            else:
                incidence.pop(bit, None)
            bits ^= bit

    merges = rounds = 0

    while active:
        budget.consume(0)
        rounds += 1
        empty = [i for i, (row, _) in active.items() if row == 0]
        for index in empty:
            zeros.append(active.pop(index)[1])
        forced = {
            mask.bit_length() - 1
            for mask in incidence.values()
            if mask.bit_count() == 1
        }
        if forced:
            for index in forced:
                toggle(index, active.pop(index)[0])
            continue
        selected, used = [], set()
        budget.consume(len(incidence) * words)
        for _, mask in sorted(incidence.items()):
            if mask.bit_count() != 2:
                continue
            bit = mask & -mask
            first, second = bit.bit_length() - 1, (mask ^ bit).bit_length() - 1
            if first in used or second in used:
                continue
            selected.append((first, second))
            used.update((first, second))
            if len(selected) == batch_size:
                break

        if not selected:
            break
        for first, second in selected:
            left, left_node = active[first]
            right, right_node = active.pop(second)
            budget.consume(words)
            toggle(second, right)
            toggle(first, right)
            if use_history:
                history.append((left_node, right_node))
                if len(history) > 2 * len(rows):
                    raise AssertionError("history node bound")
                node = len(history) - 1
            else:
                node = left_node ^ right_node

            active[first] = left ^ right, node
            merges += 1

    def lift(node):
        if not use_history:
            return node
        stack, mask = [node], 0

        while stack:
            budget.consume(words)
            value = history[stack.pop()]
            if len(value) == 1:
                mask ^= 1 << value[0]
            else:
                stack.extend(value)

        return mask

    masks = tuple(lift(node) for _, node in active.values())
    kernels = tuple(lift(node) for node in zeros)
    for mask in kernels:
        module.verify_dependency(mask, rows)
    out = tuple(row for row, _ in active.values())
    union = 0
    for row in out:
        union |= row
    return module.FilteredMatrix(
        rows,
        out,
        masks,
        kernels,
        dict(
            input_rows=len(rows),
            output_rows=len(out),
            output_columns=union.bit_count(),
            rounds=rounds,
            weight_two_merges=merges,
            history_nodes=len(history),
        ),
        reserve,
    )
