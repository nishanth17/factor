"""Worker-private polling and immutable decoding caches; no new arithmetic."""

from ..common import utils
from ..execution.budget import Budget
from ..execution.worker_contract import integer


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
        if self.until_poll <= 0:
            # A refused work reservation still checks cancel/time precedence.
            super().consume(amount)
            self.until_poll = self.poll_every
        else:
            self.used += amount
        self.until_poll -= 1


class ProspectiveActionBudget(ActionBudget):
    """Same checked contract with one prospective exact work sum."""

    def consume(self, amount=1):
        integer(amount, "atomic work")
        prospective_used = self.used + amount
        if self.until_poll <= 0:
            Budget.consume(self, amount)
            self.until_poll = self.poll_every
        else:
            self.used = prospective_used
        self.until_poll -= 1


class BoundContext:
    """Bind compatible schedule references once; never retain curve state."""

    def __init__(self, source):
        from . import core

        self.source = source
        self.segment_size = source.segment_size
        self.base_primes = source.base_primes
        self.chains = getattr(source, "chains", None)
        self.program_segment = getattr(source, "program_segment", None)
        self.power_values = getattr(source, "power_values", None)
        self.prime_segment = getattr(source, "prime_segment", None)
        self.point_add = getattr(source, "point_add", core.point_add)
        self.point_double = getattr(source, "point_double", core.point_double)
        self.stage_two_schedule = getattr(source, "stage_two_schedule", None)

    def __getattr__(self, name):
        return getattr(self.source, name)


class DecodedPrograms:
    """Bounded immutable prime/power tuples shared only inside one worker.

    Live cursors own their references separately. Never store a point/modulus.
    Byte admission and a 4096-block ceiling bound even very sparse schedules.
    Refusal leaves the ordinary packed program path available unchanged.
    """

    def __init__(self, programs, memory_bytes):
        self.programs = programs
        self.memory_bytes = integer(memory_bytes, "decoded memory")
        self.blocks = {}
        self.used_bytes = 0
        self.hits = self.misses = self.refusals = 0

    def __getattr__(self, name):
        return getattr(self.programs, name)

    def program_segment(self, lo, hi, budget, *, bound=None):
        key = (lo, hi, bound)
        if key in self.blocks:
            budget.consume(1)
            self.hits += 1
            return self.blocks[key][0]

        read = getattr(self.programs, "program_segment", None)
        if read is None:
            # Streamed bounds still share immutable primes across curves.
            budget.consume(self.segment_size + len(self.base_primes))
            values = self.programs.prime_segment(lo, hi)
        else:
            values = read(lo, hi, budget, bound=bound)
        self.misses += 1
        # Includes Python integer/references/tuple/key ownership, not RSS.
        reserve = 512 + 64 * len(values) * (2 if bound is not None else 1)
        if (
            self.used_bytes + reserve > self.memory_bytes
            or len(self.blocks) >= 4096
        ):
            # Curves scan schedules from the beginning. Evicting on a scan
            # larger than the cap can miss on every block of every curve.
            # Keep the funded prefix; stream the uncached tail exactly.
            self.refusals += 1
            return values
        budget.consume(len(values) * (2 if bound is not None else 1) + 1)
        read_powers = getattr(self.programs, "power_values", None)
        if bound is None:
            powers = None
        elif read_powers is None:
            powers = [utils.prime_power(prime, bound) for prime in values]
        else:
            powers = read_powers(lo, hi, bound, 0, len(values))
        if bound is not None and powers is None:
            self.refusals += 1
            return values
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
            read = getattr(self.programs, "power_values", None)
            return None if read is None else read(lo, hi, bound, start, count)
        stop = start + count
        values = cached[1][start:stop]
        if len(values) != count:
            raise ValueError("power slice disagrees with prime cursor")
        return values
