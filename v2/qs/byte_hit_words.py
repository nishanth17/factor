"""Optional byte-sized epochs for the existing exact recovery-word buffers."""

from .hit_words import WORD_MASK, WordHits


class ByteStampWordHits(WordHits):
    """Shrink bulk stamps; retain the original conservative owned allowance.

    Native words, descriptors and candidate-read caches keep their existing
    representation. Only per-word stamps use bytes. A funded clear before
    epoch reuse prevents stale bits after 255 blocks. The inherited hot
    marking/read methods and their logical hit charges remain unchanged.
    """

    def _new_stamps(self, width):
        self.epoch_resets = 0
        return [bytearray(width) for _ in self.words]

    def begin(self, initial, budget):
        if not 0 <= initial <= WORD_MASK:
            raise ValueError("block initialization must fit the low word")
        if self.epoch >= 255:
            # Charge each bounded clear before touching its buffer. Keep
            # the old epoch until every clear succeeds, as in WordHits.
            for stamps in self.stamps:
                self._charge(budget, self.width + 1)
                stamps[:] = bytes(self.width)
            self.epoch = 0
            self.epoch_resets += 1
        super().begin(initial, budget)

    def assign_prefix(self, values, budget):
        """Retain the low-word fast path and exact general-mask fallback."""
        width = len(values)
        if width > self.width:
            raise ValueError("initializer exceeds the configured block")
        low = all(0 <= value <= WORD_MASK for value in values)
        self._charge(budget, width + 1 if low else width * len(self.words) + 1)
        self.begin(0, budget)
        if low:
            self.words[0][:width] = values
            self.stamps[0][:width] = bytes([self.epoch]) * width
        else:
            for offset, value in enumerate(values):
                self[offset] = value
