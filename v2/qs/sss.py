"""Experimental SSS/SSSf collision collection with common QS extraction."""

import random
from dataclasses import dataclass, field, replace
from functools import partial
from math import prod

from ..common import arithmetic, utils
from ..common.arithmetic import pow
from ..execution.budget import Budget, BudgetExhaustedError
from .factor_base import build_factor_base, checked_target
from .pipeline import QSJob, QSResult
from .polynomial import qs_polynomial
from .sieve_collector import SieveCollector, SieveConfig, SieveResult
from .smooth_batch import SmoothBatch


@dataclass(frozen=True)
class SSSConfig:
    """Finite serial search; an SSSf cutoff intentionally loses yield.

    small_bound splits the factor base by prime magnitude. Zero uses its
    first fifth. A positive filter_bound rejects candidates whose non-small-
    base part is >= that value (except fully smooth values). Zero disables
    that rejection while retaining SSSf's two-stage smoothness processing.
    These choices are adaptations, not upstream parameter defaults.
    """

    backend: str = field(default="python-int", kw_only=True)
    mode: str = "sss"
    base_bound: int = 1000
    small_bound: int = 0
    selection_size: int = 6
    search_rounds: int = 256
    collision_min: int = 3
    max_candidates: int = 4096
    max_tree_bits: int = 262144
    max_tree_nodes: int = 8192
    filter_divisor: int = 10
    filter_bound: int = 0
    weight_two: bool = True
    memory_bytes: int = 64 * 1024 * 1024
    checkpoint_bytes: int = 256 * 1024
    collector: SieveConfig = field(
        default_factory=lambda: SieveConfig(
            division="roots",
            residual_bound=10000,
            max_atoms=4096,
            max_relations=2048,
            max_partials=512,
        )
    )

    def __post_init__(self):
        """Refuse oversized schedules before setup or random assignments."""
        if self.backend not in ("python-int", "gmpy2-mpz"):
            raise ValueError("unknown arithmetic backend")

        limits = {
            "base_bound": (3, 100000),
            "small_bound": (0, 100000),
            "selection_size": (1, 16),
            "search_rounds": (1, 65536),
            "collision_min": (1, 32),
            "max_candidates": (1, 4096),
            "max_tree_bits": (1, 1048576),
            "max_tree_nodes": (1, 8192),
            "filter_divisor": (1, 1024),
            "memory_bytes": (0, 1024 * 1024 * 1024),
            "checkpoint_bytes": (4096, 1024 * 1024),
        }
        for name, (low, high) in limits.items():
            value = utils.require_integer(getattr(self, name), name, low)
            if value > high:
                raise ValueError(f"{name} exceeds the SSS limit")
        utils.require_integer(self.filter_bound, "filter_bound", 0)
        if self.filter_bound.bit_length() > 4096:
            raise ValueError("filter_bound exceeds the input bit limit")
        if self.mode not in ("sss", "sssf"):
            raise ValueError("unknown SSS mode")
        if type(self.weight_two) is not bool:
            raise TypeError("weight_two must be Boolean")
        if not isinstance(self.collector, SieveConfig):
            raise TypeError("collector must be a SieveConfig")

    @property
    def metadata_reserve(self):
        """Reserve coexistence with bounded checkpoint encoding/decoding."""
        return 65536 + 8 * self.checkpoint_bytes


def collision_candidates(
    polynomial,
    roots,
    small_count,
    selected,
    coefficients,
    *,
    budget,
    minimum=3,
    max_candidates=4096,
):
    """Generate one deterministic assignment using distinct-prime collisions.

    Cumulative root switches and M/q reuse follow the pinned SSS search.
    Count each prime once per signed shift, even at a singular root. Every
    candidate includes a checked known divisor m of the original F(x).
    Exceeding the candidate cap refuses the whole unpublished assignment.
    """
    primes = tuple(roots[index].prime for index in selected)
    modulus = prod(
        primes, start=arithmetic.backend_for(polynomial.n).integer(1)
    )
    budget.consume(sum(value.bit_length() for value in coefficients) + 1)
    position = sum(coefficients[i] * roots[i].roots[0] for i in selected)
    position %= modulus
    remaining_roots = roots[small_count:]
    inverses = []

    # Each unpublished metadata chunk has at most 64 primes. Reserve the
    # identical logical work once, then keep invariant arithmetic in locals.
    for start in range(0, len(remaining_roots), 64):
        end = start + 64
        chunk = remaining_roots[start:end]
        budget.consume(
            sum(
                modulus.bit_length() + r.prime.bit_length() ** 2 for r in chunk
            )
        )
        inverses.extend(pow(modulus, -1, r.prime) for r in chunk)

    candidates, seen = [], set()

    for changed in selected:
        entry = roots[changed]
        if len(entry.roots) != 2:
            continue
        budget.consume(modulus.bit_length())
        position = (
            position
            + coefficients[changed] * (entry.roots[1] - entry.roots[0])
        ) % modulus
        affine = []

        for start in range(0, len(remaining_roots), 64):
            end = start + 64
            chunk = remaining_roots[start:end]
            budget.consume(
                sum(
                    position.bit_length() + 2 * r.prime.bit_length()
                    for r in chunk
                )
            )
            for entry, inverse in zip(chunk, inverses[start:end]):
                shifts = tuple(
                    (root - position) * inverse % entry.prime
                    for root in entry.roots
                )
                if len(shifts) == 2 and shifts[0] == shifts[1]:
                    shifts = shifts[:1]
                elif len(shifts) > 2:
                    shifts = tuple(dict.fromkeys(shifts))
                affine.append(shifts)

        for dropped in (None,) + selected:
            if dropped == changed:
                continue
            prime = 1 if dropped is None else roots[dropped].prime
            step = modulus // prime
            counts = {}

            for start in range(0, len(remaining_roots), 64):
                end = start + 64
                chunk = remaining_roots[start:end]
                budget.consume(sum(4 * r.prime.bit_length() for r in chunk))
                for entry, shifts in zip(chunk, affine[start:end]):
                    # Certified roots are distinct. Multiplication by the
                    # dropped small prime is a unit at every remaining prime;
                    # positive and negative representatives cannot overlap.
                    for shift in shifts:
                        residue = shift * prime % entry.prime
                        counts[residue] = counts.get(residue, 0) + 1
                        negative = residue - entry.prime
                        counts[negative] = counts.get(negative, 0) + 1

            shifts = sorted(counts)

            for start in range(0, len(shifts), 64):
                end = start + 64
                chunk = shifts[start:end]
                budget.consume(len(chunk))
                for shift in chunk:
                    if counts[shift] < minimum:
                        continue
                    argument = position + shift * step
                    if argument in seen:
                        continue
                    if len(candidates) >= max_candidates:
                        raise MemoryError("SSS candidate limit reached")
                    budget.consume(2 * (abs(argument).bit_length() + 1))
                    value = polynomial.value(argument)
                    if value % step:
                        raise ArithmeticError(
                            "SSS CRT divisor invariant failed"
                        )
                    candidates.append(
                        (int(argument), arithmetic.divexact(abs(value), step))
                    )
                    seen.add(argument)

    return tuple(candidates)


class SSSCollector(SieveCollector):
    """Reuse checked admission/provenance; search indices replace sieve x.

    Refusal retains the first uncommitted candidate inside an assignment.
    Only completed assignments advance next_position. The generated batch,
    collision counters, trees and retained rows all have reserved storage.
    """

    def __init__(
        self,
        polynomial,
        factor_base,
        *,
        config=None,
        budget=None,
        search_config=None,
        seed=7,
    ):
        self.search = (
            search_config if search_config is not None else SSSConfig()
        )
        utils.require_integer(seed, "seed", 0)
        if seed.bit_length() > 4096:
            raise ValueError("seed exceeds the input bit limit")
        if polynomial.a != 1 or polynomial.multiplier != 1:
            raise ValueError("SSS currently requires A=1 and multiplier=1")
        config = config if config is not None else self.search.collector
        config = replace(
            config,
            division="roots",
            threshold_extra=0,
            score_policy="adaptive",
            memory_bytes=self.search.memory_bytes,
        )
        count = len(factor_base.entries)
        bits, nodes = self.search.max_tree_bits, self.search.max_tree_nodes
        self.search_workspace = 65536 + 2048 * count
        self.search_workspace += 1024 * self.search.max_candidates
        self.search_workspace += 384 * nodes + 8 * bits * (
            nodes.bit_length() + 2
        )
        if self.search_workspace > config.memory_bytes:
            raise MemoryError(
                "SSS search/tree reservation exceeds memory_bytes"
            )
        # Reserve search alongside inherited setup before either allocates.
        super().__init__(
            polynomial,
            factor_base,
            budget=budget,
            config=replace(
                config,
                memory_bytes=config.memory_bytes - self.search_workspace,
            ),
        )
        self.config = config
        self._workspace += self.search_workspace
        self.seed = seed
        self.small_count = (
            sum(entry.prime < self.search.small_bound for entry in self._roots)
            if self.search.small_bound
            else max(1, count // 5)
        )
        if not 1 <= self.small_count < count:
            raise ValueError("SSS needs nonempty small and remaining bases")
        self.budget.consume(sum(r.prime.bit_length() for r in self._roots))
        small_primes = tuple(r.prime for r in self._roots[: self.small_count])
        small_product = arithmetic.backend_for(self.factor_base.n).integer(
            prod(small_primes)
        )
        coefficients = []

        for prime in small_primes:
            self.budget.consume(
                small_product.bit_length() + prime.bit_length() ** 2
            )
            quotient = arithmetic.divexact(small_product, prime)
            # The CRT coefficient is one at this prime and zero at the rest.
            coefficients.append(quotient * pow(quotient, -1, prime))

        self.coefficients = tuple(coefficients)
        cut = max(1, count // self.search.filter_divisor)
        self.batch = self._batch(factor_base.primes)
        self.first_batch = (
            self._batch(factor_base.primes[:cut])
            if (self.search.mode == "sssf")
            else None
        )
        self._assignment, self._cursor, self._round = None, 0, None
        self._last_stop = 0

    def _batch(self, primes):
        return SmoothBatch(
            primes,
            backend=arithmetic.backend_for(self.factor_base.n).name,
            budget=self.budget,
            max_bits=self.search.max_tree_bits,
            max_nodes=self.search.max_tree_nodes,
            memory_bytes=self.search_workspace,
        )

    def _candidate_passes(self, value, offset):
        """SSS has no sieve-score threshold; tree results select candidates."""
        return True

    def assignment(self, index):
        """Seed each index independently; retries repeat the identical work."""
        generator = random.Random(f"factor-sss:{self.seed}:{index}")
        self.budget.consume(self.search.selection_size)
        selected = tuple(
            sorted(
                set(
                    generator.choices(
                        range(self.small_count),
                        k=self.search.selection_size,
                    )
                )
            )
        )
        return collision_candidates(
            self.polynomial,
            self._roots,
            self.small_count,
            selected,
            self.coefficients,
            budget=self.budget,
            minimum=self.search.collision_min,
            max_candidates=self.search.max_candidates,
        )

    def _prepare(self, index, stats):
        candidates = self.assignment(index)
        values = tuple(value for _, value in candidates if value)
        self.batch.budget = self.budget
        residuals = None
        if self.first_batch is not None:
            self.first_batch.budget = self.budget
            first = self.first_batch.residuals(values)
            eligible = tuple(
                i
                for i, value in enumerate(first)
                if (
                    value == 1
                    or self.search.filter_bound == 0
                    or value < self.search.filter_bound
                )
            )
            second = self.batch.residuals(tuple(first[i] for i in eligible))
            residuals = dict(zip(eligible, second))
            stats["filter_rejections"] += len(values) - len(eligible)

        if residuals is None:
            residuals = dict(enumerate(self.batch.residuals(values)))
        admissible, value_index = [], 0

        for position, value in candidates:
            if value == 0:
                admissible.append(position)
                continue
            if residuals.get(value_index, self.config.residual_bound + 1) <= (
                self.config.residual_bound
            ):
                admissible.append(position)
            value_index += 1

        stats["generated_candidates"] += len(candidates)
        stats["tree_rejections"] += len(candidates) - len(admissible)
        self._assignment, self._cursor, self._round = (
            tuple(admissible),
            0,
            index,
        )

    def collect(self, lo, hi):
        """Collect assignments; retain an unfinished candidate prefix."""
        utils.require_integer(lo, "lo", 0)
        utils.require_integer(hi, "hi", lo)
        if hi > self.search.search_rounds or lo != self._last_stop:
            raise ValueError(
                "SSS collection must resume at its next assignment"
            )
        stats = dict.fromkeys(
            (
                "scanned",
                "assignments",
                "generated_candidates",
                "filter_rejections",
                "tree_rejections",
                "candidates",
                "zeros",
                "division_primes",
                "division_steps",
                "composite_residuals",
                "admitted_atoms",
                "matches",
                "duplicates",
                "evictions",
                "dropped_partials",
            ),
            0,
        )
        index, divisor, reason = lo, None, "complete"

        try:
            while index < hi:
                if self._assignment is None:
                    self._prepare(index, stats)
                while self._cursor < len(self._assignment):
                    position = self._assignment[self._cursor]
                    atom, divisor = self._divide(position, 0, stats)
                    if atom is not None and divisor is None:
                        refusal = self._admit(atom, stats)
                        if refusal:
                            reason = refusal
                            break

                    self._cursor += 1
                    stats["scanned"] += 1
                    if divisor is not None:
                        reason = "factor_found"
                        break

                if reason != "complete":
                    break
                self._assignment, self._round = None, None
                index += 1
                stats["assignments"] += 1
        except BudgetExhaustedError:
            reason = self.budget.reason
        except MemoryError:
            reason = "memory_limit"

        self._last_stop = index
        return SieveResult(
            tuple(self._atoms.values()),
            tuple(self._full),
            tuple(self._combined),
            tuple(self._pending.values()),
            divisor,
            index,
            reason,
            stats,
            self._workspace,
        )


class SSSJob:
    """Experimental serial factoring challenger; no dispatcher promotion.

    run(batch_limit=...) pauses between assignments. Extend the same Budget
    in place to resume refused work without resetting cumulative expenditure.
    Serialized resume rebuilds and verifies the retained store, assignment
    and solver prefix under the cumulative allowance.
    """

    def __init__(self, n, *, seed=7, config=None, budget=None):
        checked_target(n, 1)
        utils.require_integer(seed, "seed", 0)
        if seed.bit_length() > 4096:
            raise ValueError("seed exceeds the input bit limit")
        self.n, self.seed = n, seed
        self.config = config if config is not None else SSSConfig()
        if not isinstance(self.config, SSSConfig):
            raise TypeError("config must be an SSSConfig")
        self.n = arithmetic.get_backend(self.config.backend).integer(n)
        self.budget = (
            budget if budget is not None else Budget(work_limit=200_000_000)
        )
        self._original_budget = self.budget
        self._setup_memory_refused = False
        self.pipeline, self.divisor, self.base = None, None, None

    def _setup(self):
        """Reserve checkpoint storage alongside collection before setup."""
        memory = self.config.memory_bytes - self.config.metadata_reserve
        if memory <= 0:
            raise MemoryError("SSS checkpoint reserve exceeds memory_bytes")
        base = build_factor_base(
            self.n,
            bound=self.config.base_bound,
            budget=self.budget,
            memory_bytes=memory,
        )
        if base.divisor is not None:
            self.divisor = base.divisor
            return
        search = replace(self.config, memory_bytes=memory)
        constructor = partial(
            SSSCollector, search_config=search, seed=self.seed
        )
        self.pipeline = QSJob(
            qs_polynomial(base.factor_base),
            base.factor_base,
            0,
            self.config.search_rounds,
            budget=self.budget,
            config=replace(self.config.collector, memory_bytes=memory),
            weight_two=self.config.weight_two,
            batch_width=1,
            collector_class=constructor,
        )
        self.base = base.factor_base

    def checkpoint(self):
        """Return a compact, byte-capped full-job checkpoint."""
        from .sss_checkpoint import pack_job

        return pack_job(self)

    @classmethod
    def from_checkpoint(cls, checkpoint, *, budget, config=None):
        """Restore checked state while retaining and charging resources."""
        from .sss_checkpoint import restore_job

        try:
            return restore_job(checkpoint, budget=budget, config=config)
        except (
            KeyError,
            TypeError,
            IndexError,
            AttributeError,
            RecursionError,
        ) as error:
            raise ValueError("malformed SSS checkpoint") from error

    def run(self, *, batch_limit=None):
        """Charge setup, collection and common QS factor extraction."""
        if self.budget is not self._original_budget:
            raise ValueError(
                "resume by extending the original Budget in place"
            )
        if batch_limit is not None:
            utils.require_integer(batch_limit, "batch_limit", 1)
        try:
            if self._setup_memory_refused:
                return QSResult(
                    None,
                    self.n,
                    "memory_limit",
                    0,
                    {"work_used": self.budget.used},
                )

            if self.pipeline is None and self.divisor is None:
                self._setup()
            if self.divisor is not None:
                return QSResult(
                    self.divisor,
                    self.n // self.divisor,
                    "factor_found",
                    0,
                    {"work_used": self.budget.used},
                )

            self.pipeline.budget = self.budget
            result = self.pipeline.run(batch_limit=batch_limit)
            result.stats["workspace_bytes"] = (
                result.stats["workspace_bytes"] + self.config.metadata_reserve
            )
            if result.reason == "window_exhausted":
                result = replace(result, reason="search_exhausted")
            return result
        except BudgetExhaustedError:
            reason = self.budget.reason
        except MemoryError:
            reason = "memory_limit"
            self._setup_memory_refused = self.pipeline is None

        position = self.pipeline.next_position if self.pipeline else 0
        stats = dict(self.pipeline.stats) if self.pipeline else {}
        stats["work_used"] = self.budget.used
        if "workspace_bytes" in stats:
            stats["workspace_bytes"] += self.config.metadata_reserve
        return QSResult(None, self.n, reason, position, stats)
