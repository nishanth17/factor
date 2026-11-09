"""Reversible component-local singleton filtering over append-only rows."""


class IncrementalStructure:
    """Retain originals and recompute only components touched by new rows.

    A new relation can revive previously peeled rows. Connected components
    therefore own original rows, never just the last surviving core. This is
    derived state: restore rebuilds it with charged work from checked parity.
    """

    def __init__(self, memory_bytes):
        self.limit = memory_bytes
        self.rows, self.parent, self.members, self.columns = [], {}, {}, {}
        self.dirty = set()
        self.alive = self.nonzeros = self.row_bytes = 0
        self.disabled = False
        self.cores, self.core_sizes = {}, {}
        self.core_policy = None
        self.owned = 65536

    def _root(self, row):
        while self.parent[row] != row:
            row = self.parent[row]
        return row

    @staticmethod
    def _columns(bits):
        while bits:
            bit = bits & -bits
            yield bit.bit_length() - 1
            bits ^= bit

    def observe(self, rows, budget):
        try:
            return self._observe(rows, budget)
        except MemoryError:
            # Derived metadata can be discarded; original rows remain owned.
            self.rows.clear()
            self.parent.clear()
            self.members.clear()
            self.columns.clear()
            self.dirty.clear()
            self.cores.clear()
            self.core_sizes.clear()
            self.disabled = True
            return None

    def _observe(self, rows, budget):
        if self.disabled:
            return None
        budget.consume(len(rows) + 1)
        if len(rows) < len(self.rows) or any(
            a != b for a, b in zip(rows, self.rows)
        ):
            raise ValueError("incremental parity prefix changed")
        needed = (
            65536
            + 512 * len(rows)
            + self.row_bytes
            + sum(
                4 * ((r.bit_length() + 7) // 8) for r in rows[len(self.rows) :]
            )
            + 192
            * (
                self.nonzeros
                + sum(r.bit_count() for r in rows[len(self.rows) :])
            )
        )
        if needed > self.limit:
            self.disabled = True
            self.rows.clear()
            self.parent.clear()
            self.members.clear()
            self.columns.clear()
            self.dirty.clear()
            self.cores.clear()
            self.core_sizes.clear()
            return None
        self.owned = needed
        if needed + sum(self.core_sizes.values()) > self.limit:
            self.cores.clear()
            self.core_sizes.clear()
        for index in range(len(self.rows), len(rows)):
            row = rows[index]
            columns = tuple(self._columns(row))
            roots = {
                self._root(self.columns[c])
                for c in columns
                if c in self.columns
            }
            touched = 1 + sum(len(self.members[root]) for root in roots)
            budget.consume(
                touched + len(columns) * (len(rows).bit_length() + 1)
            )
            root = (
                max(roots, key=lambda r: (len(self.members[r]), -r))
                if roots
                else index
            )
            group = {index}
            for old in roots:
                group.update(self.members[old])
            # Prepare the merged membership before publishing the new prefix.
            self.parent[index] = root
            if root == index:
                self.parent[root] = root
            for old in roots:
                self.parent[old] = root
                self.members.pop(old)
                self.cores.pop(old, None)
                self.core_sizes.pop(old, None)
                self.dirty.discard(old)
            self.members[root] = group
            for column in columns:
                self.columns[column] = root
            self.rows.append(row)
            self.nonzeros += len(columns)
            self.row_bytes += 4 * ((row.bit_length() + 7) // 8)
            self.dirty.add(root)
        for root in sorted(tuple(self.dirty)):
            members = self.members[root]
            incidence = {}
            for index in sorted(members):
                budget.consume(self.rows[index].bit_count() + 1)
                for column in self._columns(self.rows[index]):
                    incidence.setdefault(column, set()).add(index)
            queue = [
                column
                for column, values in incidence.items()
                if len(values) == 1
            ]
            active = set(members)
            while queue:
                column = queue.pop()
                values = incidence[column]
                if len(values) != 1:
                    continue
                index = next(iter(values))
                budget.consume(self.rows[index].bit_count() + 1)
                active.remove(index)
                for affected in self._columns(self.rows[index]):
                    incidence[affected].remove(index)
                    if len(incidence[affected]) == 1:
                        queue.append(affected)
            removed_mask = sum(1 << i for i in members)
            kept_mask = sum(1 << i for i in active)
            budget.consume(len(members) * (1 + len(rows) // 64))
            self.alive = (self.alive & ~removed_mask) | kept_mask
            self.dirty.remove(root)
        return self.alive

    def weight_two_core(self, budget, merge_policy):
        """Reuse exact component contractions; new rows revive originals.

        Only unchanged components keep their contractions. Every changed
        component is filtered again from original rows, with local masks
        lifted back to the append-only global row numbering.
        """
        from .linear_algebra import filter_matrix, matrix_workspace

        if self.disabled:
            return None
        if self.core_policy != merge_policy:
            self.cores.clear()
            self.core_sizes.clear()
            self.core_policy = merge_policy
        try:
            plan, total_output, peak_scratch = [], 0, 0
            for root in sorted(self.members):
                members = tuple(sorted(self.members[root]))
                rows = tuple(self.rows[i] for i in members)
                budget.consume(len(rows) + sum(r.bit_count() for r in rows))
                union = 0
                for row in rows:
                    union |= row
                # Conservatively reserve the wide matrix even when the
                # batch filter will compact columns, plus lifted outputs.
                output_bytes = len(rows) * (
                    512 + 8 * ((len(self.rows) + union.bit_length() + 7) // 8)
                )
                scratch = matrix_workspace(len(rows), union.bit_length())
                total_output += output_bytes
                peak_scratch = max(peak_scratch, scratch)
                plan.append((root, members, rows, union, output_bytes))
            self.limit - self.owned - total_output
            for root, members, rows, union, output_bytes in plan:
                self.core_sizes[root] = output_bytes
                if root in self.cores:
                    continue
                matrix = filter_matrix(
                    rows,
                    weight_two=True,
                    merge_policy=merge_policy,
                    budget=budget,
                )
                columns = tuple(self._columns(union))
                compact = union.bit_length() > 2 * len(columns)

                def lift(mask):
                    budget.consume(len(members) * (1 + len(self.rows) // 64))
                    return sum(1 << members[i] for i in self._columns(mask))

                def expand(row):
                    if not compact:
                        return row
                    budget.consume(len(columns) + row.bit_count())
                    return sum(1 << columns[i] for i in self._columns(row))

                pairs = tuple(
                    (expand(row), lift(mask))
                    for row, mask in zip(matrix.rows, matrix.masks)
                )
                zeros = tuple(lift(mask) for mask in matrix.zero_dependencies)
                entry = (
                    pairs,
                    zeros,
                    matrix.stats["singletons_removed"],
                    matrix.stats["weight_two_merges"],
                )
                # Publish a complete component only; a work refusal leaves
                # it absent and retries with already-consumed charges intact.
                self.cores[root], self.core_sizes[root] = entry, output_bytes
            entries = tuple(self.cores[root] for root in sorted(self.members))
            budget.consume(len(self.rows) + 1)
            return (
                tuple(pair for entry in entries for pair in entry[0]),
                tuple(mask for entry in entries for mask in entry[1]),
                sum(entry[2] for entry in entries),
                sum(entry[3] for entry in entries),
            )
        except MemoryError:
            self.cores.clear()
            self.core_sizes.clear()
            return None
