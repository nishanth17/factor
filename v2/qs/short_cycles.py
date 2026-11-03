"""Bounded exact shortcuts through all retained large-prime edges."""

from collections import deque
from itertools import islice


class ShortCycleIndex:
    """An auxiliary index; the original forest still owns connectivity.

    Closing edges remain available as shortcuts. Each emitted cycle has a
    distinct newest edge, so using earlier closing edges preserves the
    triangular independent basis. Plans never mutate this index.
    """

    def __init__(self, edge_visits, *, strategy="short-v1"):
        if strategy not in ("short-v1", "short-v2"):
            raise ValueError("unknown bounded cycle search")
        self.edge_visits = edge_visits
        self.strategy = strategy
        self.edges, self.neighbors = {}, {}
        self.weights = {}

    @property
    def scratch_reserve(self):
        # At most one discovered vertex per inspected adjacency plus the
        # initial vertex; queue, parent/depth maps and returned path coexist.
        return 65536 + 1024 * self.edge_visits

    @staticmethod
    def edge_reserve(count):
        # Edges and two endpoint maps, including parallel-edge ownership.
        return 32768 + 2048 * count

    def add(self, identity, left, right, *, weight=()):
        if identity in self.edges:
            raise ValueError("duplicate shortcut edge identity")
        self.edges[identity] = (left, right)
        if self.strategy == "short-v2":
            self.weights[identity] = weight
        for a, b in (
            ((left, right),)
            if left == right
            else ((left, right), (right, left))
        ):
            self.neighbors.setdefault(a, {}).setdefault(b, {})[identity] = None

    def remove(self, identity):
        left, right = self.edges.pop(identity)
        self.weights.pop(identity, None)
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
        if self.strategy == "short-v2":
            return self._bidirectional(left, right, fallback, budget)
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

    def _representative(self, edges, budget):
        # A small deterministic prefix bounds tie-breaking even for a very
        # large parallel class. This proxy is not the final XOR row weight.
        candidates = tuple(islice(edges, 8))
        budget.consume(len(candidates))
        return min(candidates, key=lambda edge: self.weights[edge])

    def _bidirectional(self, left, right, fallback, budget):
        """Return a shorter simple path with a shared two-frontier cap."""
        limit = 256 if fallback is None else len(fallback)
        if left == right:
            return (), 0
        if limit <= 1 or not self.edge_visits:
            return fallback, 0

        budget.consume(1)
        direct = self.neighbors.get(left, {}).get(right)
        if direct:
            return (self._representative(direct, budget),), 1
        # With no direct edge, no path can improve a two-edge fallback.
        if limit == 2:
            return fallback, 1

        parents = [{left: None}, {right: None}]
        depths = [{left: 0}, {right: 0}]
        frontiers = [[left], [right]]
        visits = 1
        while all(frontiers) and visits < self.edge_visits:
            budget.consume(sum(map(len, frontiers)))
            side = min(
                (0, 1),
                key=lambda i: sum(
                    len(self.neighbors.get(v, ())) for v in frontiers[i]
                ),
            )
            other = 1 - side
            following = []
            for vertex in frontiers[side]:
                if depths[side][vertex] + 1 >= limit:
                    continue
                for neighbor, edges in self.neighbors.get(vertex, {}).items():
                    if visits >= self.edge_visits:
                        return fallback, visits
                    budget.consume(1)
                    visits += 1
                    if neighbor in parents[side]:
                        continue
                    parents[side][neighbor] = (
                        vertex,
                        self._representative(edges, budget),
                    )
                    depths[side][neighbor] = depths[side][vertex] + 1
                    following.append(neighbor)
                    if neighbor not in parents[other] or (
                        depths[side][neighbor] + depths[other][neighbor]
                        >= limit
                    ):
                        continue
                    paths = []
                    for tree in parents:
                        cursor, path = neighbor, []
                        while tree[cursor] is not None:
                            budget.consume(1)
                            cursor, edge = tree[cursor]
                            path.append(edge)
                        paths.append(path)
                    # The first frontier intersection gives disjoint parent
                    # paths except at the meeting vertex. Preserve original
                    # edges; the caller appends the distinct newest edge.
                    return tuple(reversed(paths[0])) + tuple(paths[1]), visits
            frontiers[side] = following
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
