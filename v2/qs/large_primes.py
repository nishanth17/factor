"""Finite residual splitting and transactional all-component cycle forests."""

from dataclasses import dataclass
from itertools import islice
from math import isqrt

from .. import utils
from ..pollard_rho import RhoStats, factorize_rho
from .relations import MAX_COMBINED_ATOMS, MAX_RESIDUAL

MAX_GRAPH_EDGES = 65536


def split_two_primes(residual, bound, budget, *, known_composite=False):
    """Return two proven bounded primes, or an explicit finite rejection.

    Worst-case rho work is reserved before either walk. A walk has at most
    2,048 evaluations including recovery; cancellation is polled before and
    after this bounded operation. No failed split implies primality.
    """
    utils.require_integer(residual, "residual", 4)
    utils.require_integer(bound, "large prime bound", 2)
    if residual > MAX_RESIDUAL**2 or bound > MAX_RESIDUAL:
        raise ValueError("double-large-prime input exceeds its finite bound")
    bits = residual.bit_length()
    budget.consume(bits**2)
    if not known_composite and utils.classify_prime(residual) is (
        utils.Primality.PROVEN
    ):
        return (), "prime", 0
    root = isqrt(residual)
    if root * root == residual:
        divisor, evaluations = root, 0
    else:
        budget.consume(4096 * bits**2)
        stats = RhoStats()
        divisor = factorize_rho(
            residual,
            seed=7,
            max_attempts=2,
            max_evaluations=2048,
            batch_size=32,
            recovery_limit=64,
            stats=stats,
            _known_composite=True,
        )
        evaluations = stats.evaluations
        budget.consume(0)
    if divisor is None:
        return (), "split_failed", evaluations
    if not utils.valid_divisor(divisor, residual):
        raise ArithmeticError("residual splitter returned an improper divisor")
    pair = tuple(sorted((int(divisor), int(residual // divisor))))
    if pair[-1] > bound:
        return (), "endpoint_bound", evaluations
    budget.consume(sum(prime.bit_length() ** 2 for prime in pair))
    if any(
        utils.classify_prime(p) is not utils.Primality.PROVEN for p in pair
    ):
        return (), "not_two_primes", evaluations
    if pair[0] * pair[1] != residual:
        raise ArithmeticError("residual split failed reconstruction")
    return pair, "square" if pair[0] == pair[1] else "double", evaluations


def graph_reserve(edges):
    """Cover two forests plus traversal/planning scratch, not shared atoms.

    A forest has at most 2E vertices, each a <=40-bit label. Maps retain
    component/parent/depth and adjacency; keys/atomic IDs are shared. 4KiB per
    edge reserves simultaneous old/new maps during a staged FIFO rebuild.
    """
    utils.require_integer(edges, "graph edges", 0)
    if edges > MAX_GRAPH_EDGES:
        raise MemoryError("large-prime graph edge cap")
    return 32768 + 4096 * edges


class CycleTooLongError(Exception):
    """The unique forest path cannot fit the retained provenance contract."""


@dataclass(frozen=True)
class LinkPlan:
    identity: str
    left: int
    right: int
    root: int
    old_root: int
    size: int
    updates: tuple


@dataclass(frozen=True)
class GraphPlan:
    cycle: tuple | None = None
    link: LinkPlan | None = None
    rebuilt: object = None
    evicted: tuple = ()
    dropped: bool = False


class LargePrimeForest:
    """Spanning forest with pinned cycle paths and FIFO unowned edges.

    Closing edges belong to immutable combined rows, not connectivity. Every
    retained closing edge therefore has a unique independent fundamental
    cycle, including in components without vertex1. Planning never mutates
    this forest; commit has no arithmetic, allocations of large scratch, or
    cancellation points after its caller's final resource check.
    """

    def __init__(self):
        self.edges, self.adjacency = {}, {}
        self.parent, self.depth, self.component, self.sizes = {}, {}, {}, {}
        self.unowned = {}

    @staticmethod
    def _vertices(left, right):
        for vertex in (left, right):
            utils.require_integer(vertex, "large-prime vertex", 1)
            if vertex > MAX_RESIDUAL:
                raise ValueError("large-prime vertex exceeds its bound")

    def _identity(self, identity):
        if not isinstance(identity, str) or not 1 <= len(identity) <= 64:
            raise ValueError("unbounded atomic graph identity")
        if identity in self.edges:
            raise ValueError("duplicate forest edge identity")

    def path(self, left, right, budget):
        """Return a bounded path, None for disconnected vertices, or refuse."""
        self._vertices(left, right)
        budget.consume(1)
        if left == right:
            return ()
        if left not in self.component or right not in self.component:
            return None
        if self.component[left] != self.component[right]:
            return None
        path = []
        while left != right:
            budget.consume(1)
            if len(path) >= MAX_COMBINED_ATOMS - 1:
                raise CycleTooLongError()
            if self.depth[left] < self.depth[right]:
                left, right = right, left
            left, identity = self.parent[left]
            path.append(identity)
        return tuple(path)

    def plan_link(self, left, right, identity, budget):
        """Precompute a union by rerooting only the smaller component."""
        self._vertices(left, right)
        self._identity(identity)
        if len(self.edges) >= MAX_GRAPH_EDGES:
            raise MemoryError("large-prime graph edge cap")
        a, b = self.component.get(left, left), self.component.get(right, right)
        if a == b:
            raise ValueError("link would make a forest cycle")
        size_a, size_b = self.sizes.get(a, 1), self.sizes.get(b, 1)
        if size_a < size_b:
            left, right, a, b = right, left, b, a
        stack = [(right, left, identity, self.depth.get(left, 0) + 1)]
        updates = []
        while stack:
            vertex, parent, edge, depth = stack.pop()
            budget.consume(1 + len(self.adjacency.get(vertex, ())))
            updates.append((vertex, parent, edge, depth))
            for child, child_edge in self.adjacency.get(vertex, {}).items():
                if child != parent:
                    stack.append((child, vertex, child_edge, depth + 1))
        return LinkPlan(
            identity, left, right, a, b, size_a + size_b, tuple(updates)
        )

    def commit_link(self, plan):
        if plan.left not in self.adjacency:
            self.adjacency[plan.left] = {}
            self.parent[plan.left] = (plan.left, None)
            self.depth[plan.left] = 0
            self.component[plan.left] = plan.root
        if plan.right not in self.adjacency:
            self.adjacency[plan.right] = {}
        for vertex, parent, edge, depth in plan.updates:
            self.parent[vertex] = parent, edge
            self.depth[vertex] = depth
            self.component[vertex] = plan.root
        self.sizes.pop(plan.old_root, None)
        self.sizes[plan.root] = plan.size
        self.adjacency[plan.left][plan.right] = plan.identity
        self.adjacency[plan.right][plan.left] = plan.identity
        self.edges[plan.identity] = plan.left, plan.right
        self.unowned[plan.identity] = None

    def plan(self, left, right, identity, max_partials, budget):
        self._identity(identity)
        utils.require_integer(max_partials, "max_partials", 0)
        if max_partials > MAX_GRAPH_EDGES:
            raise ValueError("partial cap exceeds graph limit")
        path = self.path(left, right, budget)
        if path is not None:
            return GraphPlan(cycle=path)
        if not max_partials:
            return GraphPlan(dropped=True)
        if len(self.unowned) < max_partials:
            return GraphPlan(
                link=self.plan_link(left, right, identity, budget)
            )
        evicted = tuple(islice(self.unowned, max(1, max_partials // 8)))
        removed = set(evicted)
        rebuilt = LargePrimeForest()
        for edge, (a, b) in self.edges.items():
            budget.consume(1)
            if edge in removed:
                continue
            rebuilt.commit_link(rebuilt.plan_link(a, b, edge, budget))
            if edge not in self.unowned:
                del rebuilt.unowned[edge]
        rebuilt.commit_link(rebuilt.plan_link(left, right, identity, budget))
        return GraphPlan(rebuilt=rebuilt, evicted=evicted)

    def commit(self, plan):
        if plan.cycle is not None:
            for identity in plan.cycle:
                self.unowned.pop(identity, None)
        elif plan.rebuilt is not None:
            self.__dict__ = plan.rebuilt.__dict__
        elif plan.link is not None:
            self.commit_link(plan.link)
