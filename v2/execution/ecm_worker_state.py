"""Worker-private polling and immutable decoding caches; no new arithmetic."""

from collections import OrderedDict

from .budget import Budget
from .worker_contract import integer


class ActionBudget(Budget):
    """Exact work checks with clocks/cancellation every bounded N charges.

    Stage changes and publication force real checks. Deferred checks
    are cooperative, bounded by N charged atomic actions, not by milliseconds.
    """

    def __init__(self, *, poll_every=1, **kwargs):
        super().__init__(**kwargs)
        self.poll_every = integer(poll_every, "poll interval", 1, 64)
        self.until_poll = 0
        self.phase = None

    def poll(self):
        super().consume(0)
        self.until_poll = self.poll_every

    def enter_phase(self, phase):
        if phase != self.phase:
            self.poll()
            self.phase = phase

    def consume(self, amount=1):
        integer(amount, "atomic work")
        if self.until_poll <= 0 or amount > self.work_limit - self.used:
            # A refused work reservation still checks cancel/time precedence.
            super().consume(amount)
            self.until_poll = self.poll_every
        else:
            self.used += amount
        self.until_poll -= 1


class DecodedPrograms:
    """Bounded immutable prime/power tuples shared only inside one worker.

    Live cursors own their references separately. Never store a point/modulus.
    Refusal leaves the ordinary packed program path available unchanged.
    """

    def __init__(self, programs, memory_bytes):
        self.programs = programs
        self.memory_bytes = integer(memory_bytes, "decoded memory")
        self.blocks = OrderedDict()
        self.used_bytes = 0
        self.hits = self.misses = self.refusals = 0

    def __getattr__(self, name):
        return getattr(self.programs, name)

    def program_segment(self, lo, hi, budget, *, bound=None):
        key = (lo, hi, bound)
        if key in self.blocks:
            budget.consume(1)
            self.blocks.move_to_end(key)
            self.hits += 1
            return self.blocks[key][0]

        values = self.programs.program_segment(lo, hi, budget, bound=bound)
        self.misses += 1
        # Includes Python integer/references/tuple/key ownership, not RSS.
        reserve = 512 + 64 * len(values) * (2 if bound is not None else 1)
        if reserve > self.memory_bytes:
            self.refusals += 1
            return values
        budget.consume(len(values) * (2 if bound is not None else 1) + 1)
        powers = (
            None
            if bound is None
            else self.programs.power_values(lo, hi, bound, 0, len(values))
        )
        if bound is not None and powers is None:
            self.refusals += 1
            return values
        while self.blocks and (
            self.used_bytes + reserve > self.memory_bytes
            or len(self.blocks) >= 64
        ):
            _, (_, _, owned) = self.blocks.popitem(last=False)
            self.used_bytes -= owned
        entry = (
            tuple(values),
            None if powers is None else tuple(powers),
            reserve,
        )
        self.blocks[key] = entry
        self.used_bytes += reserve
        return entry[0]

    def power_values(self, lo, hi, bound, start, count):
        cached = self.blocks.get((lo, hi, bound))
        if cached is None or cached[1] is None:
            return self.programs.power_values(lo, hi, bound, start, count)
        stop = start + count
        values = cached[1][start:stop]
        if len(values) != count:
            raise ValueError("power slice disagrees with prime cursor")
        return values
