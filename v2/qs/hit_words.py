"""Bounded native-width recovery bits with lazy candidate materialization."""

from ..common import utils

WORD_BITS = 60
WORD_MASK = (1 << WORD_BITS) - 1
MAX_HIT_BYTES = 32 * 2**20
MAX_EXTENDED_HIT_BYTES = 128 * 2**20


class WordHits:
    """Keep marking in unboxed integer lists; assemble masks only on reads.

    Epochs give untouched words their block's initial value without clearing
    every factor-base word at every position. This is derived scratch: resume
    reconstructs it, and no checkpoint or canonical relation owns its bytes.
    """

    @staticmethod
    def reservation(count, width):
        groups = (count + WORD_BITS - 1) // WORD_BITS
        # Include both word arrays, descriptors, candidate-cache stamps,
        # capacity growth and temporary slices during periodic initialization.
        return 65536 + 512 * count + 32 * width * (groups + 1)

    @staticmethod
    def construction_work(count, width):
        return (
            2 * ((count + WORD_BITS - 1) // WORD_BITS) * width + count + width
        )

    def __init__(self, count, cache, budget, *, byte_limit=MAX_HIT_BYTES):
        utils.require_integer(count, "hit count", 1)
        utils.require_integer(byte_limit, "hit word bytes", 1)
        if byte_limit > MAX_EXTENDED_HIT_BYTES:
            raise ValueError("hit word allowance exceeds its finite ceiling")
        width = len(cache)
        if not width or self.reservation(count, width) > byte_limit:
            raise ValueError("hit word scratch exceeds its bounded allowance")
        work = self.construction_work(count, width)
        budget.consume(work)
        self.work = work
        self.width = width
        self.words = [[0] * width for _ in range((count + 59) // 60)]
        self.stamps = self._new_stamps(width)
        self.descriptors = [
            (self.words[i // 60], self.stamps[i // 60], 1 << (i % 60))
            for i in range(count)
        ]
        self.epoch = self.cache_epoch = 0
        self.default = 0
        self.cache = cache
        self.cache_stamps = [0] * width
        self.sealed = False
        # Private buffer-poisoning oracles use negative integers. Production
        # recovery masks are nonnegative; keep that rare view operation exact.
        self.negative = {}

    def _charge(self, budget, amount):
        budget.consume(amount)
        self.work += amount

    def _invalidate(self, budget=None):
        if self.cache_epoch == WORD_MASK:
            if budget is not None:
                self._charge(budget, self.width + 1)
            self.cache_stamps[:] = [0] * self.width
            self.cache_epoch = 0
        self.cache_epoch += 1

    def begin(self, initial, budget):
        """Reserve each reset before changing its corresponding scratch."""
        if not 0 <= initial <= WORD_MASK:
            raise ValueError("block initialization must fit the low word")
        if self.epoch == WORD_MASK:
            # Retain the old epoch until all clears succeed. A refusal can
            # safely retry; no stale epoch becomes current after wraparound.
            for stamps in self.stamps:
                self._charge(budget, self.width + 1)
                stamps[:] = [0] * self.width
            self.epoch = 0
        self._invalidate(budget)
        self._charge(budget, 1)
        self.epoch += 1
        self.default = initial
        self.sealed = False
        self.negative.clear()

    def seal(self):
        """Permit read caching after dense and sparse marking both finish."""
        self.sealed = True

    def _unseal(self):
        if self.sealed:
            self._invalidate()
            self.sealed = False

    def add(self, index, offsets, scores=None, weight=0):
        """Fuse scores and first-level metadata on one certified range."""
        self._unseal()
        words, stamps, bit = self.descriptors[index]
        epoch = self.epoch
        initial = self.default if index < WORD_BITS else 0
        if self.negative:
            for offset in offsets:
                self[offset] = self[offset] | (1 << index)
                if scores is not None:
                    scores[offset] += weight
        elif scores is None:
            for offset in offsets:
                if stamps[offset] != epoch:
                    words[offset], stamps[offset] = initial | bit, epoch
                else:
                    words[offset] |= bit
        else:
            for offset in offsets:
                if stamps[offset] != epoch:
                    words[offset], stamps[offset] = initial | bit, epoch
                else:
                    words[offset] |= bit
                scores[offset] += weight

    def add_one(self, index, offset):
        """Apply a queued first-level event without constructing a wide bit."""
        self._unseal()
        if self.negative:
            self[offset] = self[offset] | (1 << index)
            return
        words, stamps, bit = self.descriptors[index]
        if stamps[offset] != self.epoch:
            words[offset] = (self.default if index < WORD_BITS else 0) | bit
            stamps[offset] = self.epoch
        else:
            words[offset] |= bit

    def assign_prefix(self, values, budget):
        """Copy a periodic initializer; retain a general multiword fallback."""
        width = len(values)
        if width > self.width:
            raise ValueError("initializer exceeds the configured block")
        low = all(0 <= value <= WORD_MASK for value in values)
        self._charge(budget, width + 1 if low else width * len(self.words) + 1)
        self.begin(0, budget)
        if low:
            self.words[0][:width] = values
            self.stamps[0][:width] = [self.epoch] * width
        else:
            for offset, value in enumerate(values):
                self[offset] = value

    def mask(self, offset, budget=None):
        if offset < 0:
            offset += self.width
        if not 0 <= offset < self.width:
            raise IndexError("hit offset outside the configured block")
        if offset in self.negative:
            return self.negative[offset]
        if self.sealed and self.cache_stamps[offset] == self.cache_epoch:
            return self.cache[offset]
        if budget is not None:
            self._charge(budget, len(self.words) + 1)
        value = 0
        for group, words in enumerate(self.words):
            word = (
                words[offset]
                if self.stamps[group][offset] == self.epoch
                else (self.default if group == 0 else 0)
            )
            value |= word << (WORD_BITS * group)
        if self.sealed:
            self.cache[offset] = value
            self.cache_stamps[offset] = self.cache_epoch
        return value

    def __len__(self):
        return self.width

    def __iter__(self):
        return (self.mask(i) for i in range(self.width))

    def __eq__(self, other):
        return list(self) == other

    def __getitem__(self, offset):
        if isinstance(offset, slice):
            return [self.mask(i) for i in range(*offset.indices(self.width))]
        return self.mask(offset)

    def __setitem__(self, offset, value):
        if isinstance(offset, slice):
            indices = range(*offset.indices(self.width))
            if len(indices) != len(value):
                raise ValueError("hit buffer has a fixed width")
            for i, item in zip(indices, value):
                self[i] = item
            return
        if offset < 0:
            offset += self.width
        if not 0 <= offset < self.width:
            raise IndexError("hit offset outside the configured block")
        self._invalidate()
        self.sealed = False
        if value < 0:
            self.negative[offset] = value
            return
        if value.bit_length() > WORD_BITS * len(self.words):
            raise ValueError("hit mask exceeds the allocated words")
        self.negative.pop(offset, None)
        for group, words in enumerate(self.words):
            words[offset] = (value >> (WORD_BITS * group)) & WORD_MASK
            self.stamps[group][offset] = self.epoch

    def _new_stamps(self, width):
        """Keep the default native integer representation unchanged."""
        return [[0] * width for _ in self.words]
