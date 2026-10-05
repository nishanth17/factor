"""Half-open prime streams with bounded packed storage and private scratch."""

from array import array
from collections import OrderedDict
from math import isqrt

from . import utils


class SieveContext:
    """Reuse packed base primes up to max_hi, with one active iterator.

    memory_bytes caps owned array/bytearray payloads and temporary marking
    slices, with a 4096-byte reserve for small object headers. It is not a
    process RSS limit. Consumer-retained output is outside this workspace cap.
    Endpoints are Python integers; output above 2**32 is never narrowed.
    """

    def __init__(
        self,
        max_hi,
        *,
        memory_bytes=1_048_576,
        segment_size=4096,
        rolling=False,
        budget=None,
    ):
        utils.require_integer(max_hi, "max_hi", 3)
        utils.require_integer(memory_bytes, "memory_bytes", 8192)
        utils.require_integer(segment_size, "segment_size", 1)
        self.max_hi = max_hi
        self.memory_bytes = memory_bytes
        self.segment_size = segment_size
        self.rolling = rolling
        self.active = False
        root = isqrt(max_hi - 1)
        size = (root + 1) // 2
        # Bound construction conservatively before allocating. Worst-case
        # base storage includes every odd candidate, plus rolling offsets.
        required = 4096 + 18 * size + 3 * segment_size
        if root >= 2**64 or required > memory_bytes:
            raise MemoryError("base-prime workspace exceeds configured cap")
        if budget is not None:
            budget.consume(size)
        flags = bytearray(b"\x01") * size
        if flags:
            flags[0] = 0
        for prime in range(3, isqrt(root) + 1, 2):
            if flags[prime // 2]:
                start = prime * prime // 2
                count = (size - 1 - start) // prime + 1
                flags[start::prime] = b"\x00" * count
        self.base_primes = array(
            "Q", (2 * i + 1 for i in range(1, size) if flags[i])
        )
        self._flags = bytearray(segment_size)
        self.payload_cap = required

    def prime_segment(self, lo, hi):
        """Materialize one bounded segment without per-prime resumptions.

        This shares the context's private marking buffer with primes(); the
        same single-consumer rule and half-open endpoints apply. The caller
        reserves generation work before invoking this cursor-oriented path.
        """
        utils.require_integer(lo, "lo")
        utils.require_integer(hi, "hi")
        if hi > self.max_hi or hi - lo > 2 * self.segment_size:
            raise ValueError("segment exceeds context bounds")
        if self.active:
            raise RuntimeError("context already has an active iterator")
        self.active = True
        try:
            if hi <= max(lo, 2):
                return []
            left = max(lo, 3) | 1
            size = max(0, (hi - left + 1) // 2)
            self._flags[:size] = b"\x01" * size
            for prime in self.base_primes:
                if prime * prime >= hi:
                    break
                first = max(
                    prime * prime, ((left + prime - 1) // prime) * prime
                )
                if first % 2 == 0:
                    first += prime
                index = (first - left) // 2
                if index < size:
                    count = (size - 1 - index) // prime + 1
                    self._flags[index:size:prime] = b"\x00" * count
            values = [left + 2 * i for i in range(size) if self._flags[i]]
            if lo <= 2 < hi:
                values.insert(0, 2)
            return values
        finally:
            self.active = False

    def primes(self, lo, hi, *, budget=None):
        """Yield lo <= p < hi; release ownership on exhaustion or close.

        Rolling strikes restart from exact multiples on every call. A caller
        must close an abandoned iterator before reusing this context.
        """
        utils.require_integer(lo, "lo")
        utils.require_integer(hi, "hi")
        if hi > self.max_hi:
            raise ValueError("endpoint exceeds context maximum")
        if self.active:
            raise RuntimeError("context already has an active iterator")
        self.active = True
        try:
            if hi <= max(lo, 2):
                return
            if lo <= 2 < hi:
                if budget is not None:
                    budget.consume()
                yield 2
            start = max(lo, 3) | 1
            strikes = array("Q")
            if self.rolling:
                for prime in self.base_primes:
                    first = max(
                        prime * prime, ((start + prime - 1) // prime) * prime
                    )
                    if first % 2 == 0:
                        first += prime
                    # Absolute strikes may exceed 64 bits on high intervals.
                    # Store relative offsets instead, bounded by 2*p.
                    strikes.append((first - start) // 2)
            for left in range(start, hi, 2 * self.segment_size):
                right = min(hi, left + 2 * self.segment_size)
                size = (right - left + 1) // 2
                if budget is not None:
                    budget.consume(size + len(self.base_primes))
                self._flags[:size] = b"\x01" * size
                for position, prime in enumerate(self.base_primes):
                    if prime * prime >= right and not self.rolling:
                        break
                    if self.rolling:
                        index = strikes[position]
                    else:
                        first = max(
                            prime * prime,
                            ((left + prime - 1) // prime) * prime,
                        )
                        if first % 2 == 0:
                            first += prime
                        index = (first - left) // 2
                    if index < size:
                        count = (size - 1 - index) // prime + 1
                        stop = index + count * prime
                        self._flags[index:size:prime] = b"\x00" * count
                        index = stop
                    if self.rolling:
                        strikes[position] = index - size
                for index in range(size):
                    if self._flags[index]:
                        yield left + 2 * index
        finally:
            self.active = False


def iter_primes(lo, hi, **options):
    """Yield a bounded fresh stream; **options configure SieveContext."""
    utils.require_integer(lo, "lo")
    utils.require_integer(hi, "hi")
    budget = options.pop("budget", None)
    if hi <= max(lo, 2):
        return
    context = SieveContext(max(hi, 3), budget=budget, **options)
    yield from context.primes(lo, hi, budget=budget)


def prime_powers(bound, **options):
    """Yield (prime, largest exact prime power <= bound), inclusively."""
    utils.require_integer(bound, "bound", 2)
    for prime in iter_primes(2, bound + 1, **options):
        yield prime, utils.prime_power(prime, bound)


class ScheduleCache:
    """Optional bounded cache of integer primes, gaps, or exact powers.

    Keys include half-open endpoints and representation. Modulus-dependent
    residues/points never enter this cache. Only one active consumer is allowed
    so an evicted array cannot remain invisibly live through another iterator.
    Large schedules stream without caching. Arrays use unsigned 64-bit values;
    wider outputs remain Python integers and are deliberately uncached.
    """

    def __init__(self, context, *, cache_bytes=65_536, max_entries=8):
        utils.require_integer(cache_bytes, "cache_bytes", 4096)
        utils.require_integer(max_entries, "max_entries", 1)
        self.context = context
        self.cache_bytes = cache_bytes
        self.max_entries = max_entries
        self.entries = OrderedDict()
        self.used_bytes = 0
        self.active = False
        self.hits = 0
        self.misses = 0

    @property
    def segment_size(self):
        """Expose the underlying stream's finite segment width."""
        return self.context.segment_size

    @property
    def base_primes(self):
        """Expose packed base primes for work accounting."""
        return self.context.base_primes

    def primes(self, lo, hi):
        """Provide the prime-cursor interface without caching curve state."""
        yield from self.values(lo, hi)

    def values(self, lo, hi, *, kind="primes", bound=None):
        """Yield a schedule; gaps are relative to lo, then preceding primes."""
        utils.require_integer(lo, "lo")
        utils.require_integer(hi, "hi")
        if kind not in ("primes", "gaps", "powers"):
            raise ValueError("unknown integer schedule representation")
        if kind != "powers" and bound is not None:
            raise ValueError("bound applies only to prime-power schedules")
        if kind == "powers":
            utils.require_integer(bound, "bound", max(lo, 2))
            if hi > bound + 1:
                raise ValueError("power schedule includes primes above bound")
        if self.active:
            raise RuntimeError("cache already has an active consumer")
        self.active = True
        key = (lo, hi, "half-open", kind, bound)
        try:
            if key in self.entries:
                self.hits += 1
                self.entries.move_to_end(key)
                yield from self.entries[key]
                return
            self.misses += 1
            # Worst-case packed size, with capacity-growth/header allowance.
            reserve = 1024 + 16 * max(0, (hi - max(lo, 2) + 1) // 2 + 1)
            retain = (
                reserve <= self.cache_bytes
                and max(hi, bound or 0) <= (2**64)
                and (kind != "powers" or bound < 2**64)
            )
            if retain:
                while self.entries and (
                    self.used_bytes + reserve > self.cache_bytes
                    or len(self.entries) >= self.max_entries
                ):
                    _, removed = self.entries.popitem(last=False)
                    self.used_bytes -= 1024 + 16 * len(removed)
            values = array("Q") if retain else None
            previous = lo
            stream = self.context.primes(lo, hi)
            try:
                for prime in stream:
                    if kind == "gaps":
                        value, previous = prime - previous, prime
                    elif kind == "powers":
                        value = utils.prime_power(prime, bound)
                    else:
                        value = prime
                    if values is not None:
                        values.append(value)
                    yield value
            finally:
                stream.close()
            if values is not None:
                self.entries[key] = values
                self.used_bytes += 1024 + 16 * len(values)
        finally:
            self.active = False
