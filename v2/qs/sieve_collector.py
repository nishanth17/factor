"""Bounded conservative score sieving and single-large-prime matching."""

from array import array
from dataclasses import dataclass
from math import gcd

from .. import utils
from ..budget import Budget, BudgetExhaustedError
from .factor_base import DEFAULT_MEMORY_BYTES
from .polynomial import checked_position, polynomial_roots
from .relations import (
    AtomicRelation,
    checked_residual_bound,
    combine_relations,
    verify_atomic,
)

MAX_BLOCK_WIDTH = 4096
MAX_WINDOW_WIDTH = 1_000_000
MAX_STORED_ATOMS = 4096
MAX_METADATA_CHUNK = 1024


@dataclass(frozen=True)
class SieveConfig:
    """Finite collector controls; threshold_extra > 0 is intentionally lossy.

    Scores upper-bound factor-base log contributions, including all prime
    powers. Zero extra guarantees candidate coverage, not admission: exact
    division and the verifier decide admission. bytearray scores saturate at
    255; thresholds saturate too, allowing extra false candidates, never
    missed admissible values. Cutoffs omit marking only, never exact recovery.
    """

    block_width: int = 256
    metadata_chunk: int = 64
    score_backend: str = "list"
    marking: str = "sparse"
    division: str = "roots"
    small_prime_cutoff: int = 0
    threshold_extra: int = 0
    score_policy: str = "adaptive"
    residual_bound: int = 500
    max_partials: int = 256
    max_relations: int = 1024
    max_atoms: int = 2048
    memory_bytes: int = DEFAULT_MEMORY_BYTES

    def __post_init__(self):
        """Validate bounded choices before allocating buffers or roots."""
        limits = {
            "block_width": (1, MAX_BLOCK_WIDTH),
            "metadata_chunk": (1, MAX_METADATA_CHUNK),
            "small_prime_cutoff": (0, 100_000),
            "threshold_extra": (0, 65536),
            "max_partials": (0, MAX_STORED_ATOMS),
            "max_relations": (0, MAX_STORED_ATOMS),
            "max_atoms": (0, MAX_STORED_ATOMS),
        }
        for name, (minimum, maximum) in limits.items():
            value = getattr(self, name)
            utils.require_integer(value, name, minimum)
            if value > maximum:
                raise ValueError(f"{name} exceeds the collector limit")
        utils.require_integer(self.memory_bytes, "memory_bytes", 0)
        checked_residual_bound(self.residual_bound)
        if self.score_backend not in ("list", "bytearray", "array"):
            raise ValueError("unknown score backend")
        if self.marking not in ("dense", "sparse", "bucket"):
            raise ValueError("unknown marking policy")
        if self.division not in ("full", "roots", "bucket", "resieve"):
            raise ValueError("unknown division policy")
        if self.score_policy not in ("adaptive", "conservative", "candidate"):
            raise ValueError("unknown score policy")


@dataclass(frozen=True)
class SieveResult:
    """Cumulative checked store and per-call diagnostics for [lo,hi).

    atoms includes full/combined provenance and unmatched partials. FIFO
    eviction touches unmatched atoms only. next_position is the first
    uncommitted position; resuming on the same collector retains matches.
    Returned snapshots own references; release old snapshots before resuming
    if their memory is to be excluded from the collector's reservation.
    """

    atoms: tuple
    full_relations: tuple
    combined_relations: tuple
    partial_ids: tuple
    divisor: int | None
    next_position: int
    reason: str
    stats: dict
    workspace_bytes: int


class SieveCollector:
    """Cache one polynomial's roots/logs and reuse bounded working blocks.

    Constructor setup can raise ValueError, MemoryError or budget exhaustion.
    Collection refusals return a checked prefix. Single-thread use only;
    state is in-memory, not the serialized SIQS checkpoint of P3.4.
    """

    def __init__(self, polynomial, factor_base, *, config=None, budget=None):
        """Reserve metadata/buffers and validate A's factor-base support."""
        self.config = config if config is not None else SieveConfig()
        self.polynomial, self.factor_base = polynomial, factor_base
        self.budget = budget if budget is not None else Budget()
        if (polynomial.n, polynomial.multiplier) != (
            factor_base.n,
            factor_base.multiplier,
        ):
            raise ValueError("polynomial and factor-base identity mismatch")
        count = len(factor_base.entries)
        # Include roots, Python object overhead, score integers, hit bitsets,
        # candidate temporaries, and two bounded translation slices. Account
        # for combination scratch before retaining any atom.
        self._workspace = factor_base.workspace_bytes + 32768 + count * 512
        self._workspace += self.config.block_width * (128 + (count + 7) // 8)
        if self.config.division == "resieve":
            # Candidate-only exact recovery retains bounded values and
            # sparse exponent lists. Reserve the dense worst case as well.
            self._workspace += self.config.block_width * (512 + 128 * count)
        self._workspace += 128 * (polynomial.n_prime.bit_length() + 16384)
        if self._workspace > self.config.memory_bytes:
            raise MemoryError("collector setup exceeds memory_bytes")
        remaining, a_exponents, roots, logs = polynomial.a, [], [], []
        for entry in factor_base.entries:
            self.budget.consume(polynomial.a.bit_length() + 1)
            exponent = 0
            while remaining % entry.prime == 0:
                remaining //= entry.prime
                exponent += 1
            a_exponents.append(exponent)
            roots.append(
                polynomial_roots(
                    polynomial, factor_base, entry, budget=self.budget
                )
            )
            # ceil(log2 p) uses exact bit lengths, with 2 exactly one bit.
            logs.append((entry.prime - 1).bit_length())
        if remaining != 1:
            raise ValueError("A must factor completely over the factor base")
        self._roots, self._logs = tuple(roots), tuple(logs)
        self._a_exponents = tuple(a_exponents)
        width = self.config.block_width
        if self.config.score_backend == "bytearray":
            self._scores = bytearray(width)
        elif self.config.score_backend == "array":
            self._scores = array("I", [0]) * width
        else:
            self._scores = [0] * width
        self._hits = [0] * width
        self._resieved = {}
        self._scratch_bytes = 0
        self._scratch_peak_bytes = 0
        self._omitted_allowance = 0
        self._atoms, self._pending = {}, {}
        self._full, self._combined = [], []
        self._atom_bytes = {}

    def _bounds(self, lo, hi):
        """Bound |F| with endpoints and the two integer vertex neighbors.

        A positive quadratic's extrema on integer positions occur here.
        A sign crossing gets lower bound zero, conservatively selecting all
        positions. This evaluates only O(1) big integers per working block.
        """
        vertex = -self.polynomial.b // self.polynomial.a
        positions = {lo, hi - 1}
        for position in (vertex, vertex + 1):
            if lo <= position < hi:
                positions.add(position)
        values = [self.polynomial.value(position) for position in positions]
        minimum, maximum = min(values), max(values)
        lower = 0 if minimum <= 0 <= maximum else min(map(abs, values))
        return lower, max(map(abs, values))

    def _weights(self, maximum):
        """Upper-bound all p-power log contributions using exact powers.

        If p divides F, its valuation is at most floor(log_p(max |F|)).
        Multiplying that bound by ceil(log2 p) cannot underestimate any
        smooth contribution. No Hensel depth cutoff can lose high powers.
        """
        weights = []
        for entry, log in zip(self.factor_base.entries, self._logs):
            self.budget.consume(maximum.bit_length() + 1)
            exponent, power = 0, entry.prime
            while power <= maximum:
                exponent += 1
                power *= entry.prime
            weights.append(exponent * log)
        return weights

    def _add(self, indices, weight, stats):
        """Update scores; byte translations allocate at most block width."""
        if self.config.score_backend == "bytearray":
            if isinstance(indices, range) and indices:
                table = bytes(min(255, value + weight) for value in range(256))
                start, stop, step = indices.start, indices.stop, indices.step
                segment = self._scores[start:stop:step]
                self._scores[start:stop:step] = segment.translate(table)
                stats["slice_bytes"] += 2 * len(segment) + 256
                stats["slice_allocations"] += 3
            else:
                for index in indices:
                    self._scores[index] = min(
                        255, self._scores[index] + weight
                    )
        else:
            for index in indices:
                self._scores[index] += weight

    def _sieve(self, lo, hi, stats):
        """Mark complete root hits in bounded prime chunks, without F(x)."""
        width = hi - lo
        self._omitted_allowance = 0
        self.budget.consume(width)
        for index in range(width):
            self._scores[index] = 0
            self._hits[index] = 0
        lower, maximum = self._bounds(lo, hi)
        lower_threshold = max(0, lower.bit_length() - 1)
        lower_threshold -= (self.config.residual_bound - 1).bit_length()
        if (
            self.config.score_policy == "adaptive"
            and lower_threshold <= 0
            and self.config.threshold_extra == 0
            and self.config.division != "bucket"
        ):
            # Every score is nonnegative. A zero lower threshold selects
            # every position, so weight calculation/marking cannot affect
            # coverage. Root division or candidate-only resieving recovers
            # exact exponents later; buckets still need their hit metadata.
            stats["blocks"] += 1
            stats["threshold_min"] = 0
            stats["skipped_score_blocks"] = (
                stats.get("skipped_score_blocks", 0) + 1
            )
            return 0
        weights = self._weights(maximum)
        omitted = 0
        buckets = (
            self.config.marking == "bucket" or self.config.division == "bucket"
        )
        chunk = self.config.metadata_chunk
        for start in range(0, len(self._roots), chunk):
            for index in range(start, min(start + chunk, len(self._roots))):
                roots, weight = self._roots[index], weights[index]
                prime = roots.prime
                skipped = prime < self.config.small_prime_cutoff
                if skipped:
                    omitted += weight
                residues = (0,) if roots.all_positions else roots.roots
                step = 1 if roots.all_positions else prime
                for root in residues:
                    offset = (root - lo) % step
                    hits = range(offset, width, step)
                    self.budget.consume(len(hits) + 1)
                    stats["root_hits"] += len(hits)
                    if buckets:
                        for hit in hits:
                            self._hits[hit] |= 1 << index
                    if skipped or self.config.marking == "bucket":
                        continue
                    if self.config.marking == "sparse" and step >= width:
                        hits = (offset,) if offset < width else ()
                    self._add(hits, weight, stats)
        if self.config.marking == "bucket":
            for offset in range(width):
                bits = self._hits[offset]
                while bits:
                    bit = bits & -bits
                    index = bit.bit_length() - 1
                    if (
                        self._roots[index].prime
                        >= self.config.small_prime_cutoff
                    ):
                        self._add((offset,), weights[index], stats)
                    bits ^= bit
        # log2 |F| - log2 residual <= log2 factor-base part. Lower-bound
        # the first term and upper-bound the second; omitted primes get a
        # proved universal allowance. Positive extra deliberately loses it.
        threshold = max(0, lower.bit_length() - 1)
        self._omitted_allowance = omitted
        threshold -= (self.config.residual_bound - 1).bit_length() + omitted
        threshold = max(0, threshold + self.config.threshold_extra)
        if self.config.score_backend == "bytearray":
            threshold = min(255, threshold)
        stats["blocks"] += 1
        stats["threshold_min"] = min(stats["threshold_min"], threshold)
        stats["threshold_max"] = max(stats["threshold_max"], threshold)
        return threshold

    def _candidate_passes(self, value, offset):
        """Refine a safe threshold using an already evaluated candidate norm.

        floor(log2 |F(x)|) replaces the block-wide lower bound. Scores still
        upper-bound all prime-power contributions; omitted primes retain
        their full allowance. Thus zero-extra refinement cannot reject an
        admissible norm. Clip both sides for byte scores as in block scoring.
        """
        if self.config.score_policy != "candidate":
            return True
        threshold = abs(value).bit_length() - 1
        threshold -= (self.config.residual_bound - 1).bit_length()
        threshold -= self._omitted_allowance
        threshold = max(0, threshold + self.config.threshold_extra)
        if self.config.score_backend == "bytearray":
            threshold = min(255, threshold)
        return self._scores[offset] >= threshold

    def _resieve(self, lo, hi, threshold, stats):
        """Recover exact candidate exponents in a separate root-hit pass.

        Only selected candidate values are formed. Iterate prime/root hits
        across this block, divide every valuation to completion, and retain
        A's exponents separately. Refusal changes scratch only; no position
        or checked store is committed before the full recovery pass returns.
        """
        self._resieved.clear()
        _, maximum = self._bounds(lo, hi)
        scratch = (hi - lo) * (512 + 4 * maximum.bit_length())
        if self._workspace + scratch > self.config.memory_bytes:
            return False
        self._scratch_bytes = scratch
        self._scratch_peak_bytes = max(self._scratch_peak_bytes, scratch)
        values, remaining, recovered = {}, {}, {}
        for offset in range(hi - lo):
            if self._scores[offset] >= threshold:
                self.budget.consume(self.polynomial.n_prime.bit_length() + 1)
                value = self.polynomial.value(lo + offset)
                if not self._candidate_passes(value, offset):
                    continue
                values[offset] = value
                remaining[offset] = abs(value)
                recovered[offset] = []
        for roots, a_exponent in zip(self._roots, self._a_exponents):
            self.budget.consume(len(values) + 1)
            prime = roots.prime
            residues = (0,) if roots.all_positions else roots.roots
            step = 1 if roots.all_positions else prime
            hit_exponents = {}
            for root in residues:
                for offset in range((root - lo) % step, hi - lo, step):
                    if offset not in values or remaining[offset] == 0:
                        continue
                    self.budget.consume(remaining[offset].bit_length() + 1)
                    exponent = 0
                    while remaining[offset] % prime == 0:
                        remaining[offset] //= prime
                        exponent += 1
                    hit_exponents[offset] = exponent
                    stats["division_primes"] += 1
                    stats["division_steps"] += exponent
            for offset in values:
                exponent = a_exponent + hit_exponents.get(offset, 0)
                if exponent:
                    recovered[offset].append((prime, exponent))
        self._resieved = {
            offset: (values[offset], remaining[offset], recovered[offset])
            for offset in values
        }
        return True

    def _divide(self, position, offset, stats):
        """Recover every exponent exactly; roots/buckets certify hit coverage.

        Root exclusion is exact, independent of scores and omitted primes.
        All valuations at every included prime are divided to completion.
        Full factor-base division remains an explicit reference control.
        """
        if self.config.division == "resieve":
            if offset not in self._resieved:
                self.budget.consume()
                stats["candidates"] += 1
                stats["refined_rejections"] = (
                    stats.get("refined_rejections", 0) + 1
                )
                return None, None
            value, remaining, exponents = self._resieved[offset]
        else:
            value = self.polynomial.value(position)
            remaining, exponents = abs(value), []
        self.budget.consume(
            (len(self._roots) + 1) * (abs(value).bit_length() + 1) + 1
        )
        stats["candidates"] += 1
        if value == 0:
            stats["zeros"] += 1
            divisor = gcd(
                self.polynomial.u_value(position), self.factor_base.n
            )
            return None, divisor if utils.valid_divisor(
                divisor, self.factor_base.n
            ) else None
        if not self._candidate_passes(value, offset):
            stats["refined_rejections"] = (
                stats.get("refined_rejections", 0) + 1
            )
            return None, None
        roots_to_divide = (
            ()
            if self.config.division == "resieve"
            else (enumerate(self._roots))
        )
        for index, roots in roots_to_divide:
            prime = roots.prime
            if self.config.division == "full":
                hit = True
            elif self.config.division == "bucket":
                hit = bool(self._hits[offset] & (1 << index))
            else:
                hit = roots.all_positions or position % prime in roots.roots
            exponent = self._a_exponents[index]
            if hit:
                stats["division_primes"] += 1
                while remaining % prime == 0:
                    remaining //= prime
                    exponent += 1
                    stats["division_steps"] += 1
            if exponent:
                exponents.append((prime, exponent))
        if remaining > self.config.residual_bound:
            return None, None
        self.budget.consume(remaining.bit_length() ** 2)
        if remaining != 1 and utils.classify_prime(remaining) != (
            utils.Primality.PROVEN
        ):
            stats["composite_residuals"] += 1
            return None, None
        atom = AtomicRelation(
            self.polynomial,
            position,
            -1 if value < 0 else 1,
            tuple(exponents),
            remaining,
        )
        verify_atomic(
            atom,
            self.factor_base,
            residual_bound=self.config.residual_bound,
            budget=self.budget,
        )
        if remaining != 1:
            self.budget.consume(self.factor_base.n.bit_length())
            divisor = gcd(remaining, self.factor_base.n)
            if utils.valid_divisor(divisor, self.factor_base.n):
                return atom, divisor
            if divisor != 1:
                raise ValueError("nonunit residual has no proper split")
        return atom, None

    def _admit(self, atom, stats):
        """Commit only after verification/combination and all cap checks.

        FIFO eviction cannot reach pinned atoms. Refusal leaves store state
        untouched, so the same position can be retried with a new Budget.
        """
        if atom.relation_id in self._atoms:
            stats["duplicates"] += 1
            return None
        partner = (
            self._pending.get(atom.residual) if atom.residual != 1 else None
        )
        accepted = atom.residual == 1 or partner is not None
        if accepted and len(self._full) + len(self._combined) >= (
            self.config.max_relations
        ):
            return "relation_limit"
        evicted = None
        if not accepted:
            if self.config.max_partials == 0:
                stats["dropped_partials"] += 1
                return None
            if len(self._pending) >= self.config.max_partials:
                evicted = next(iter(self._pending.values()))
        released = self._atom_bytes[evicted] if evicted else 0
        if (
            len(self._atoms) + 1 - int(evicted is not None)
            > self.config.max_atoms
        ):
            return "atom_limit"
        reserve = 4096 + 256 * len(atom.exponents)
        reserve += 16 * (
            abs(atom.position).bit_length()
            + self.polynomial.n_prime.bit_length()
        )
        if partner is not None:
            reserve += 2048 + 256 * len(self.factor_base.entries)
        if (
            self._workspace + self._scratch_bytes + reserve - released
            > self.config.memory_bytes
        ):
            return "memory_limit"
        combination = None
        if partner is not None:
            combination = combine_relations(
                (self._atoms[partner], atom),
                self.factor_base,
                budget=self.budget,
                memory_bytes=self.config.memory_bytes,
            )
            if combination.divisor is not None:
                raise ArithmeticError(
                    "previously checked unit residual changed"
                )
        self.budget.consume(len(atom.exponents) + 1)
        if evicted:
            del self._pending[self._atoms[evicted].residual]
            del self._atoms[evicted]
            del self._atom_bytes[evicted]
            stats["evictions"] += 1
        self._atoms[atom.relation_id] = atom
        self._atom_bytes[atom.relation_id] = reserve
        self._workspace += reserve - released
        stats["admitted_atoms"] += 1
        if atom.residual == 1:
            self._full.append(atom)
        elif partner is not None:
            del self._pending[atom.residual]
            self._combined.append(combination.relation)
            stats["matches"] += 1
        else:
            self._pending[atom.residual] = atom.relation_id
        return None

    def collect(self, lo, hi):
        """Collect a bounded window using reusable half-open working blocks.

        A refused block or position returns its first uncommitted x. Scoring
        scratch may change on refusal; accepted store state never does.
        Eviction and thresholds may lose yield only as documented in config.
        """
        checked_position(lo)
        checked_position(hi)
        if not 0 <= hi - lo <= MAX_WINDOW_WIDTH:
            raise ValueError("window must be ordered and at most 1000000 wide")
        stats = dict.fromkeys(
            (
                "blocks",
                "scanned",
                "candidates",
                "root_hits",
                "division_primes",
                "division_steps",
                "zeros",
                "composite_residuals",
                "admitted_atoms",
                "matches",
                "duplicates",
                "evictions",
                "dropped_partials",
                "slice_bytes",
                "slice_allocations",
                "threshold_max",
            ),
            0,
        )
        stats["threshold_min"] = 2**63
        position, reason, divisor = lo, "complete", None
        try:
            while position < hi:
                block_lo = position
                block_hi = min(hi, block_lo + self.config.block_width)
                threshold = self._sieve(block_lo, block_hi, stats)
                if self.config.division == "resieve":
                    if not self._resieve(block_lo, block_hi, threshold, stats):
                        reason = "memory_limit"
                        break
                while position < block_hi:
                    offset = position - block_lo
                    if self._scores[offset] >= threshold:
                        # Candidate division reserves the position unit and
                        # recovery work in one atomic check. The former two
                        # adjacent polls did no intervening published work.
                        atom, divisor = self._divide(position, offset, stats)
                        if divisor is not None:
                            position += 1
                            stats["scanned"] += 1
                            reason = "factor_found"
                            break
                        if atom is not None:
                            refusal = self._admit(atom, stats)
                            if refusal:
                                reason = refusal
                                break
                    else:
                        self.budget.consume()
                    position += 1
                    stats["scanned"] += 1
                if reason != "complete":
                    break
        except BudgetExhaustedError:
            reason = self.budget.reason
        if stats["blocks"] == 0:
            stats["threshold_min"] = 0
        self._resieved.clear()
        self._scratch_bytes = 0
        return SieveResult(
            tuple(self._atoms.values()),
            tuple(self._full),
            tuple(self._combined),
            tuple(self._pending.values()),
            divisor,
            position,
            reason,
            stats,
            self._workspace + self._scratch_peak_bytes,
        )
