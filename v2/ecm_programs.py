"""Bound-owned immutable ECM schedules with finite, run-local retention."""

from dataclasses import dataclass
from math import gcd
from struct import iter_unpack, pack

from . import utils

PROGRAM_VERSION = "ecm-packed-blocks-v1"
PAIRED_VERSION = "ecm-packed-pairs-v1"
WORD_LIMIT = 2**64


@dataclass(frozen=True)
class ProgramBlock:
    """One half-open prime block, optionally with inclusive-B1 powers.

    Bytes keep published records immutable and independent of a modulus.
    Decoded values belong to the consumer's already-reserved cursor workspace.
    """

    lo: int
    hi: int
    bound: int | None
    primes: bytes
    powers: bytes

    def prime_values(self):
        """Decode exact primes without exposing mutable packed storage."""
        return [value for (value,) in iter_unpack("<Q", self.primes)]


class ECMPrograms:
    """Retain completed blocks across curves; regenerate when the cap fills.

    There is no global cache or eviction. The fixed reserve covers temporary
    packing, one unretained block and small object headers. Retained records
    include a conservative per-key allowance. Generation, power compilation
    and every decode are charged before publication. A new run or resume owns
    a new store: rebuilding missing blocks consumes the cumulative allowance.
    """

    def __init__(self, context, *, memory_bytes):
        utils.require_integer(memory_bytes, "program memory", 4096)
        self.context = context
        self.segment_size = context.segment_size
        self.base_primes = context.base_primes
        self.memory_bytes = memory_bytes
        self.used_bytes = 4096 + 256 * self.segment_size
        if self.used_bytes > memory_bytes:
            raise MemoryError("ECM program scratch exceeds configured cap")
        self.blocks = {}
        self.last_key = None
        self.last_block = None
        self.hits = 0
        self.misses = 0
        self.unretained = 0
        self.coverage_blocks = {}
        self.coverage_hits = 0
        self.coverage_misses = 0

    def coverage(self, cursor, *, b1, b2, distance, budget):
        """Reuse integer certificates; decoded records belong to one curve.

        The caller reserves coverage construction and decoded workspace in
        addition to this store's cap. Checkpoints retain the current decoded
        block, so resuming never needs to regenerate a consumed prefix.
        """
        if (
            cursor["next"] - cursor["left"] > 2 * self.segment_size
            or len(cursor["values"]) > self.segment_size
        ):
            raise ValueError("paired coverage exceeds one prime segment")
        key = (cursor["left"], cursor["next"], b1, b2, distance)
        coverage = self.coverage_blocks.get(key)
        if coverage is None:
            values = cursor["values"]
            budget.consume(len(values))
            block = ProgramBlock(
                key[0],
                key[1],
                None,
                b"".join(pack("<Q", prime) for prime in values),
                b"",
            )
            coverage = pair_coverage(
                block,
                b1=b1,
                b2=b2,
                distance=distance,
                memory_bytes=4096 + 512 * len(values),
                budget=budget,
            )
            # Reserve decoding before publishing anything to the store.
            budget.consume(len(coverage.data) // 32)
            reserve = 512 + len(coverage.data)
            if self.used_bytes + reserve <= self.memory_bytes:
                self.coverage_blocks[key] = coverage
                self.used_bytes += reserve
            else:
                self.unretained += 1
            self.coverage_misses += 1
        else:
            budget.consume(len(coverage.data) // 32)
            self.coverage_hits += 1
        return [list(record) for record in coverage.records()]

    def program_segment(self, lo, hi, budget, *, bound=None):
        """Return a complete block, never publishing an exhausted prefix."""
        utils.require_integer(lo, "lo", 2)
        utils.require_integer(hi, "hi", lo)
        if hi > WORD_LIMIT or hi - lo > 2 * self.segment_size:
            raise ValueError("ECM program exceeds packed block bounds")
        if bound is not None:
            utils.require_integer(bound, "bound", 2)
            if bound >= WORD_LIMIT or hi > bound + 1:
                raise ValueError("power program exceeds inclusive B1")
        key = (lo, hi, bound)
        block = self.blocks.get(key)

        if block is None:
            # Refused reservations leave both the cursor and store unchanged.
            budget.consume(self.segment_size + len(self.base_primes))
            segment = getattr(self.context, "prime_segment", None)
            values = (
                segment(lo, hi)
                if segment is not None
                else list(self.context.primes(lo, hi))
            )
            budget.consume(len(values) * (2 if bound is not None else 1))
            prime_data = b"".join(pack("<Q", value) for value in values)
            power_data = (
                b"".join(
                    pack("<Q", utils.prime_power(prime, bound))
                    for prime in values
                )
                if bound is not None
                else b""
            )
            block = ProgramBlock(lo, hi, bound, prime_data, power_data)
            reserve = 512 + len(prime_data) + len(power_data)
            if self.used_bytes + reserve <= self.memory_bytes:
                self.blocks[key] = block
                self.used_bytes += reserve
            else:
                self.unretained += 1
            self.misses += 1
        else:
            budget.consume(len(block.primes) // 8)
            values = block.prime_values()
            self.hits += 1

        self.last_key, self.last_block = key, block
        return values

    def power_values(self, lo, hi, bound, start, count):
        """Read a power slice; resumed buffers may need fallback."""
        key = (lo, hi, bound)
        block = (
            self.last_block if key == self.last_key else self.blocks.get(key)
        )
        if block is None:
            # A checkpoint retains its prime buffer, not this run-local store.
            # Recompute that buffer's powers instead of resieving old work.
            return None
        first_byte = 8 * start
        final_byte = 8 * (start + count)
        data = block.powers[first_byte:final_byte]
        if len(data) != 8 * count:
            raise ValueError("power slice disagrees with prime cursor")
        return [value for (value,) in iter_unpack("<Q", data)]


@dataclass(frozen=True)
class CoverageBlock:
    """Bounded +/- coverage records: (center, distance, minus, plus).

    A zero distance denotes direct scalar evaluation of the recorded prime.
    Nonzero records certify one or two eligible primes, independent of points.
    These certificates are also consumed by the opt-in paired executor.
    """

    lo: int
    hi: int
    b1: int
    b2: int
    distance: int
    data: bytes

    def records(self):
        """Iterate immutable integer coverage certificates."""
        return iter_unpack("<QQQQ", self.data)


def pair_coverage(block, *, b1, b2, distance, memory_bytes, budget):
    """Compile bounded coverage for a positive existing giant recurrence.

    Centers are odd and spaced by 2*D, so odd-prime distances are even.
    The initial previous scalar must remain positive. Wheel-divisor and
    center-prime exceptions use direct scalars; no pruning is assumed.
    A block boundary can separate a pair, but never omit either prime.
    """
    utils.require_integer(b1, "B1", 2)
    utils.require_integer(b2, "B2", b1)
    utils.require_integer(distance, "distance", 0)
    utils.require_integer(memory_bytes, "coverage memory", 4096)
    origin = b1 if b1 % 2 else b1 - 1
    if b2 >= WORD_LIMIT or block.lo < b1 + 1 or block.hi > b2 + 1:
        raise ValueError("coverage block lies outside inclusive bounds")
    if distance and (distance < 2 or distance % 2 or 2 * distance >= origin):
        raise ValueError("D needs even distances and positive initialization")
    count = len(block.primes) // 8
    if 4096 + 512 * count > memory_bytes:
        raise MemoryError("coverage construction exceeds configured cap")
    budget.consume(count)
    records = {}
    previous = block.lo - 1

    for (prime,) in iter_unpack("<Q", block.primes):
        if not previous < prime < block.hi or prime < block.lo:
            raise ValueError("coverage primes must be ordered in the block")
        previous = prime
        if not distance or gcd(prime, 2 * distance) != 1:
            records[(prime, 0)] = [prime, 0]
            continue
        center = origin + ((prime - origin + distance) // (2 * distance)) * (
            2 * distance
        )
        offset = abs(prime - center)
        if center >= WORD_LIMIT:
            raise ValueError("coverage center exceeds packed word bounds")
        if offset == 0:
            records[(prime, 0)] = [prime, 0]
            continue

        pair = records.setdefault((center, offset), [0, 0])
        side = int(prime > center)
        if pair[side]:
            raise ValueError("duplicate coverage prime")
        pair[side] = prime

    data = b"".join(
        pack("<QQQQ", center, offset, *pair)
        for (center, offset), pair in sorted(records.items())
    )
    return CoverageBlock(block.lo, block.hi, b1, b2, distance, data)
