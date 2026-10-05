"""Exact work reservations with finite cooperative polling intervals."""

from . import utils


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

    def consume(self, amount=1):
        """Refuse an overdraw immediately, even between external polls."""
        utils.require_integer(amount, "amount", 0)
        if (
            amount == 0
            or self._remaining == 0
            or amount > self.budget.work_limit - self.budget.used
        ):
            if self.poll is not None:
                self.poll(self.budget)
            self.budget.consume(amount)
            self._remaining = self.interval - 1
        else:
            # Skip only the external poll; the reservation still costs work.
            self.budget.used += amount
            self._remaining -= 1
