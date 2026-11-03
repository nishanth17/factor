"""Integer-handle incidence elimination for two/three/opaque residuals."""

from collections import Counter
from dataclasses import dataclass
from heapq import nsmallest

from .relations import MAX_COMBINED_ATOMS


@dataclass(frozen=True)
class IncidencePlan:
    handle: int
    tokens: frozenset
    vector: frozenset
    provenance: frozenset
    cycle: tuple | None
    replacement: object = None


def tokens(atom, resolutions=None):
    opaque = getattr(atom, "opaque", None)
    if (
        opaque is not None
        and resolutions is not None
        and opaque in resolutions
    ):
        return frozenset(
            (False, p)
            for p, count in Counter(resolutions[opaque]).items()
            if count % 2
        )
    values = getattr(atom, "incidence_tokens", None)
    if values is None:
        values = tuple((False, p) for p in atom.residual_primes)
    return frozenset(p for p, count in Counter(values).items() if count % 2)


class ResidualKernel:
    """All-component sparse GF(2) incidence basis with bounded lifting.

    Every pivot carries original integer handles. New kernels contain their
    newest handle, so retain independence even after earlier-kernel reduction.
    No component is discarded for lacking an SLP anchor. Refusal is explicit.
    """

    def __init__(self, *, edge_limit=None):
        self.edge_limit = edge_limit
        self.edges, self.unowned, self.pivots = {}, {}, {}
        self.cycles, self.weights = [], {}
        self.direct, self.cycle_index = {}, {}
        self.resolutions = {}

    def plan_atom(self, atom, handle, config, budget):
        if type(handle) is not int or not 0 <= handle < 2**64:
            raise ValueError("invalid integer incidence handle")
        if handle in self.edges:
            raise ValueError("duplicate incidence handle")
        original = tokens(atom, self.resolutions)
        if original and len(self.pivots) >= config.partial_limit:
            # A closure does not consume another independent partial slot.
            if (
                original not in self.direct
                and max(original) not in self.pivots
            ):
                raise MemoryError("incidence partial cap")
        vector, provenance = original, frozenset((handle,))
        direct = self.direct.get(original) if original else None
        if direct is not None:
            budget.consume(1)
            vector, provenance = frozenset(), frozenset((direct, handle))
        while vector:
            pivot = max(vector)
            old = self.pivots.get(pivot)
            if old is None:
                break
            budget.consume(
                len(vector) + len(old[0]) + len(provenance) + len(old[1])
            )
            vector, provenance = vector ^ old[0], provenance ^ old[1]
            if len(provenance) > 2 * MAX_COMBINED_ATOMS:
                raise MemoryError("incidence lifting scratch cap")
        cycle = None
        if not vector:

            def weight(ids):
                budget.consume(len(ids) + 1)
                return sum(
                    len(atom.exponents) if i == handle else self.weights[i]
                    for i in ids
                )

            # Bounded postconstruction reduction, retaining original rows.
            # XOR with earlier kernels cannot remove the newest handle.
            # One-atom kernels are minimal. Two-atom direct matches already
            # use the cheapest equivalent representative; searching longer
            # history cannot improve their selected length/weight objective.
            passes = (
                config.cycle_reduction_passes if len(provenance) > 2 else 0
            )
            for _ in range(passes):
                before = (len(provenance), weight(provenance))
                candidates = set(self.cycles)
                visits = 0
                for item in sorted(provenance):
                    for previous in self.cycle_index.get(item, ()):
                        budget.consume(1)
                        candidates.add(previous)
                        visits += 1
                        if visits >= config.cycle_search_edges:
                            break
                    if visits >= config.cycle_search_edges:
                        break
                choices = nsmallest(
                    config.cycle_reduction_choices,
                    candidates,
                    key=lambda ids: (
                        len(ids),
                        weight(ids),
                        tuple(sorted(ids)),
                    ),
                )
                for old in choices:
                    budget.consume(len(provenance) + len(old))
                    candidate = provenance ^ old
                    if (len(candidate), weight(candidate)) < (
                        len(provenance),
                        weight(provenance),
                    ):
                        provenance = candidate
                if (len(provenance), weight(provenance)) == before:
                    break
            if len(provenance) > MAX_COMBINED_ATOMS:
                raise MemoryError("incidence cycle provenance cap")
            cycle = tuple(sorted(provenance))
        elif len(provenance) > MAX_COMBINED_ATOMS:
            raise MemoryError("incidence pivot provenance cap")
        if vector and len(self.pivots) >= config.partial_limit:
            raise MemoryError("incidence partial cap")
        replacement = None
        if original:
            old = self.pivots.get(max(original))
            if old is not None and old[0] == original and len(old[1]) == 1:
                prior = next(iter(old[1]))
                if len(atom.exponents) < self.weights[prior]:
                    replacement = (
                        max(original),
                        (original, frozenset((handle,))),
                    )
        return IncidencePlan(
            handle, original, vector, provenance, cycle, replacement
        )

    def commit_atom(self, plan, atom, choice_limit):
        self.edges[plan.handle] = plan.tokens
        self.weights[plan.handle] = len(atom.exponents)
        old = self.direct.get(plan.tokens)
        if plan.tokens and (
            old is None or self.weights[plan.handle] < self.weights[old]
        ):
            self.direct[plan.tokens] = plan.handle
        if plan.vector:
            self.pivots[max(plan.vector)] = plan.vector, plan.provenance
            self.unowned[plan.handle] = None
        else:
            self.cycles.append(plan.provenance)
            for handle in plan.provenance:
                entries = self.cycle_index.setdefault(handle, [])
                entries.append(plan.provenance)
                if len(entries) > choice_limit:
                    del entries[0]
            if len(self.cycles) > choice_limit:
                del self.cycles[0]
        if plan.replacement is not None:
            key, value = plan.replacement
            self.pivots[key] = value

    def state(self):
        return {
            "pivots": [
                [list(key), [list(v) for v in sorted(vector)], sorted(ids)]
                for key, (vector, ids) in sorted(self.pivots.items())
            ],
            "cycles": [sorted(ids) for ids in self.cycles],
            "unowned": list(self.unowned),
        }
