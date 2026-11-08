"""Exact work reservations with finite cooperative polling intervals."""

from ..common import utils


class PollingBudget:
    """Poll before the first action, every interval actions, and on consume(0).

    Every work reservation is checked and charged before its action. Only
    deadline/cancellation checks are amortized. At most interval - 1 bounded
    atomic actions can follow an external stop before the next poll. Callers
    must force consume(0) before publishing a completed batch. This adapter
    is private to one active execution; it shares the original Budget ledger.
    """

    def __init__(self, budget, *, interval=64, poll=None):
        utils.require_integer(interval, "poll interval", 1)
        if interval > 64:
            raise ValueError("poll interval exceeds 64 atomic actions")
        self.budget, self.interval, self.poll = budget, interval, poll
        self._remaining = 0

    def __getattr__(self, name):
        return getattr(self.budget, name)

    def __setattr__(self, name, value):
        # Public ledger updates, including cancellation installed during a
        # callback, must reach the shared Budget rather than shadow its fields.
        if name in ("budget", "interval", "poll", "_remaining"):
            object.__setattr__(self, name, value)
        else:
            setattr(self.budget, name, value)

    def consume(self, amount=1):
        """Refuse an overdraw immediately, even between external polls."""
        utils.require_integer(amount, "amount", 0)
        if (
            amount == 0
            or self._remaining == 0
        ):
            if self.poll is not None:
                self.poll(self.budget)
            self.budget.consume(amount)
            self._remaining = self.interval - 1
        else:
            # Skip only the external poll; the reservation still costs work.
            self.budget.used += amount
            self._remaining -= 1

    def consume_span(self, amount):
        """Reserve up to one interval of certified unit-cost actions.

        Unlike consume(amount), this counts each unit as a separate action.
        Poll early when the whole span would cross the next polling boundary.
        A failed reservation leaves the shared ledger and caller cursor intact.
        """
        utils.require_integer(amount, "span amount", 1)
        if amount > self.interval:
            raise ValueError("span exceeds the cooperative polling interval")
        if (
            amount > self._remaining
        ):
            if self.poll is not None:
                self.poll(self.budget)
            self.budget.consume(amount)
            self._remaining = self.interval - amount
        else:
            self.budget.used += amount
            self._remaining -= amount

    def consume_actions(self, count, work_each):
        """Reserve a group while counting every original logical action.

        Derived queue construction groups already counted events. Each event
        keeps its charge and contributes to the polling gap. A refused group
        changes neither the ledger nor any caller-owned event buffer.
        """
        utils.require_integer(count, "action count", 1)
        utils.require_integer(work_each, "work per action", 1)
        if count > self.interval:
            raise ValueError("group exceeds the cooperative polling interval")
        amount = count * work_each
        if (
            count > self._remaining
        ):
            if self.poll is not None:
                self.poll(self.budget)
            self.budget.consume(amount)
            self._remaining = self.interval - count
        else:
            self.budget.used += amount
            self._remaining -= count
