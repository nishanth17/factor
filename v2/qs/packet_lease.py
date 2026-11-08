"""Private packet lease counter; not a shared or checkpointed Budget view."""

from ..common import utils
from ..execution.work_budget import PollingBudget


class ParentCPUSnapshot:
    """Read worker counters under one lock at the unchanged parent cadence."""

    def __init__(self, pool):
        self.pool = pool
        # begin() has fixed this run's baseline; finish() drains workers
        # before another begin can replace it. Capture only this run's view.
        baseline = tuple(pool.state["cpu_base"])
        counters = pool.state["cpu"]
        values, lock = counters.get_obj(), counters.get_lock()
        samples = [0.0] * len(baseline)

        def used():
            # This reader has one owner. Retain a private consistent copy and
            # release the shared lock before subtraction, clipping and sum;
            # workers can publish while the parent reduces its snapshot.
            with lock:
                for index in range(len(samples)):
                    samples[index] = values[index]
            return sum(
                max(0, value - prior)
                for value, prior in zip(samples, baseline)
            )

        self.worker_cpu = used

    def _account_cpu(self, budget):
        """Publish one snapshot before PollingBudget checks its reservation."""
        pool = self.pool
        current = self.worker_cpu()
        budget.prior_cpu += max(0, current - pool.accounted_cpu)
        pool.accounted_cpu = current
        pool.state["parent_cpu"].value = budget.cpu_used - current

    def __call__(self, budget):
        self._account_cpu(budget)
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
        self.marking_pairs = 0

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

    def _precharge_marking_pair(self, first, second):
        """Fund two certified hit ranges without crossing a poll boundary.

        The internal caller supplies positive native range lengths plus one.
        On a boundary return false without any mutation: the original scalar
        path must retain its first affordable action and exact refusal prefix.
        """
        amount = first + second
        if self._actions_left < 2 or amount > self.work_limit - self.used:
            return False
        self.used += amount
        self._actions_left -= 2
        self.marking_pairs += 1
        return True
