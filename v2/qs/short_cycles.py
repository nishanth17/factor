"""Bounded exact shortcuts through all retained large-prime edges."""

from collections import deque


class ShortCycleIndex:
    """An auxiliary index; the original forest still owns connectivity.

    Closing edges remain available as shortcuts. Each emitted cycle has a
    distinct newest edge, so using earlier closing edges preserves the
    triangular independent basis. Plans never mutate this index.
    """

    def __init__(self, edge_visits):
        self.edge_visits = edge_visits
        self.edges, self.neighbors = {}, {}

    @property
    def scratch_reserve(self):
        # At most one discovered vertex per inspected adjacency plus the
        # initial vertex; queue, parent/depth maps and returned path coexist.
        return 65536 + 1024 * self.edge_visits

    @staticmethod
    def edge_reserve(count):
        # Edges and two endpoint maps, including parallel-edge ownership.
        return 32768 + 2048 * count

    def add(self, identity, left, right):
        if identity in self.edges:
            raise ValueError("duplicate shortcut edge identity")
        self.edges[identity] = (left, right)
        for a, b in (
            ((left, right),)
            if left == right
            else ((left, right), (right, left))
        ):
            self.neighbors.setdefault(a, {}).setdefault(b, {})[identity] = None

    def remove(self, identity):
        left, right = self.edges.pop(identity)
        for a, b in (
            ((left, right),)
            if left == right
            else ((left, right), (right, left))
        ):
            if a not in self.neighbors or b not in self.neighbors[a]:
                continue
            edges = self.neighbors[a][b]
            del edges[identity]
            if not edges:
                del self.neighbors[a][b]
            if not self.neighbors[a]:
                del self.neighbors[a]

    def shorter_path(self, left, right, fallback, budget):
        """Return a strictly shorter path or the original bounded fallback.

        A visit cap bounds unsuccessful exploration; exhausting it is not a
        graph/resource failure. Budget interruption propagates without a
        publication. A None fallback permits recovery of an overlong tree
        path only when this search finds at most 255 original edges.
        """
        limit = 256 if fallback is None else len(fallback)
        if left == right:
            return (), 0
        if limit <= 1 or self.edge_visits == 0:
            return fallback, 0
        parent, depth, queue = {left: None}, {left: 0}, deque([left])
        visits = 0
        while queue:
            vertex = queue.popleft()
            if depth[vertex] + 1 >= limit:
                continue
            for neighbor, edges in self.neighbors.get(vertex, {}).items():
                if visits >= self.edge_visits:
                    return fallback, visits
                budget.consume(1)
                visits += 1
                if neighbor in parent:
                    continue
                parent[neighbor] = (vertex, next(iter(edges)))
                depth[neighbor] = depth[vertex] + 1
                if neighbor == right:
                    path = []
                    while neighbor != left:
                        budget.consume(1)
                        neighbor, identity = parent[neighbor]
                        path.append(identity)
                    return tuple(reversed(path)), visits
                queue.append(neighbor)
        return fallback, visits


def check_prior_path(indices, atoms, closing, budget):
    """Check a simple original-edge path without trusting forest caches."""
    atom = atoms[closing]
    left, right = atom.large_primes or (1, atom.residual)
    adjacency = {}
    for index in indices:
        budget.consume(1)
        if index == closing:
            continue
        if index >= closing:
            raise ValueError("cycle shortcut is not an earlier original atom")
        other = atoms[index]
        a, b = other.large_primes or (1, other.residual)
        if a == b:
            raise ValueError("shortcut path contains a separate loop")
        adjacency.setdefault(a, []).append(b)
        adjacency.setdefault(b, []).append(a)
    if left == right:
        if adjacency:
            raise ValueError("square loop has extraneous provenance")
        return
    if any(
        len(v) != (1 if k in (left, right) else 2)
        for k, v in adjacency.items()
    ) or (left not in adjacency or right not in adjacency):
        raise ValueError("shortcut provenance is not a simple endpoint path")
    seen, pending = set(), [left]
    while pending:
        vertex = pending.pop()
        if vertex in seen:
            continue
        seen.add(vertex)
        pending.extend(adjacency[vertex])
    if len(seen) != len(adjacency):
        raise ValueError("shortcut provenance has a disconnected cycle")
