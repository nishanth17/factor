"""Private packet lease counter; not a shared or checkpointed Budget view."""

from ..common import utils
from ..execution.work_budget import PollingBudget


class ParentCPUSnapshot:
    """Read worker counters under one lock at the unchanged parent cadence."""

    def __init__(self, pool):
        from .worker_packet import _process_cpu_reader

        self.pool = pool
        # begin() has fixed this run's baseline; finish() drains workers
        # before another begin can replace it. Capture only this run's view.
        self.worker_cpu = _process_cpu_reader(pool.state)

    def __call__(self, budget):
        pool = self.pool
        current = self.worker_cpu()
        budget.prior_cpu += max(0, current - pool.accounted_cpu)
        pool.accounted_cpu = current
        pool.state["parent_cpu"].value = budget.cpu_used - current
        budget.consume(0)


class PacketLeaseBudget(PollingBudget):
    """Keep prepaid work in plain fields and poll at most every 64 actions.

    Only this private packet may access the backing ledger. Synchronization
    before a poll gives callbacks the same paid prefix as PollingBudget;
    reconciliation after a packet copies that prefix once to its report.
    The raw Budget retains deadline/cancellation precedence and exact refusal.
    """

    def __init__(self, budget, *, poll=None):
        self.budget = budget
        self.work_limit = budget.work_limit
        self.used = budget.used
        self.interval = 64
        self.poll = poll
        self._actions_left = 0

    # This private view has no shared-ledger forwarding on its hot path.
    __setattr__ = object.__setattr__

    @property
    def reason(self):
        return self.budget.reason

    @reason.setter
    def reason(self, value):
        self.budget.reason = value

    @property
    def _remaining(self):
        return self._actions_left

    def _charge_poll(self, amount, actions):
        self.budget.used = self.used
        if self.poll is not None:
            self.poll(self.budget)
        self.budget.consume(amount)
        self.work_limit = self.budget.work_limit
        self.used = self.budget.used
        self._actions_left = self.interval - actions

    def consume(self, amount=1):
        utils.require_integer(amount, "amount", 0)
        if (
            amount == 0
            or not self._actions_left
            or amount > self.work_limit - self.used
        ):
            self._charge_poll(amount, 1)
        else:
            self.used += amount
            self._actions_left -= 1

    def consume_span(self, amount):
        utils.require_integer(amount, "span amount", 1)
        if amount > self.interval:
            raise ValueError("span exceeds cooperative polling interval")
        if amount > self._actions_left or amount > self.work_limit - self.used:
            self._charge_poll(amount, amount)
        else:
            self.used += amount
            self._actions_left -= amount

    def consume_actions(self, count, work_each):
        utils.require_integer(count, "action count", 1)
        utils.require_integer(work_each, "work per action", 1)
        if count > self.interval:
            raise ValueError("group exceeds cooperative polling interval")
        amount = count * work_each
        if count > self._actions_left or amount > self.work_limit - self.used:
            self._charge_poll(amount, count)
        else:
            self.used += amount
            self._actions_left -= count

    def reconcile(self):
        """Retain every paid action even when publication is cancelled."""
        self.budget.used = self.used
        return self.used
