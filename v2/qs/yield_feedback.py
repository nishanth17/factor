"""Finite useful-outcome feedback and reporting adaptation with hysteresis."""

from ..common import utils
from .yield_support import UsefulTarget


class OutcomeTarget(UsefulTarget):
    def __init__(self):
        super().__init__()
        self.attempts = self.trivial = self.streak = 0

    def extracted(self, count, trials, skips, divisor, growth, limit):
        self.attempts += 1
        self.trivial += skips + sum(trial.divisor is None for trial in trials)
        self.streak = min(8, self.streak + 1) if divisor is None else 0
        if divisor is None:
            # Failed square roots require new arithmetic information even
            # when nullity is sufficient. Bound escalation by one quarter.
            extra = min(
                max(growth, (count + 3) // 4), growth * (1 + self.streak)
            )
            desired = max(self.next_target, count + extra)
            self.next_target = (
                desired if limit is None else min(limit, desired)
            )

    def state(self):
        return dict(
            super().state(),
            attempts=self.attempts,
            trivial=self.trivial,
            streak=self.streak,
        )

    def restore(self, state, rows, collector, budget):
        if not isinstance(state, dict) or set(state) != set(self.state()):
            raise ValueError("invalid outcome-target state")
        original = UsefulTarget()
        original.restore(
            {key: state[key] for key in original.state()},
            rows,
            collector,
            budget,
        )
        for name in ("attempts", "trivial", "streak"):
            utils.require_integer(state[name], name, 0)
            if state[name] > (8 if name == "streak" else 2**64 - 1):
                raise ValueError("outcome feedback exceeds its finite cap")
        if state["streak"] > state["attempts"]:
            raise ValueError("outcome streak exceeds attempts")
        self.__dict__.update(original.__dict__)
        self.attempts, self.trivial, self.streak = (
            state["attempts"],
            state["trivial"],
            state["streak"],
        )


class ReportingController:
    """Explore a bounded reporting interval only at completed block boundaries.

    This is an exploratory policy, not an optimum inferred from raw relations.
    Two consecutive low-useful-yield windows are required before widening.
    It never changes endpoint admission or touches held-out confirmation data.
    """

    def __init__(self, lower, upper, exposure, epochs):
        if not 1 <= lower <= upper:
            raise ValueError("invalid reporting interval")
        self.lower, self.upper = lower, upper
        self.exposure, self.epochs = exposure, epochs
        self.bound = lower
        self.positions = self.candidates = self.rows = self.bad = (
            self.epoch
        ) = 0

    def observe(self, positions, candidates, rows):
        self.record(positions, candidates, rows)
        return self.decide()

    def record(self, positions, candidates, rows):
        """Persist committed progress across refusal and resume."""
        self.positions += positions
        self.candidates += candidates
        self.rows += rows

    def decide(self):
        if self.positions < self.exposure or self.epoch >= self.epochs:
            return self.bound
        # A row is a closed exact relation, not an unmatched residual edge.
        low = self.rows * 64 < max(1, self.candidates)
        self.bad = self.bad + 1 if low else 0
        if self.bad >= 2:
            self.bound = min(self.upper, max(self.bound + 1, 2 * self.bound))
            self.bad = 0
        self.epoch += 1
        self.positions = self.candidates = self.rows = 0
        return self.bound

    def state(self):
        return {
            name: getattr(self, name)
            for name in (
                "bound",
                "positions",
                "candidates",
                "rows",
                "bad",
                "epoch",
            )
        }

    def restore(self, state):
        if not isinstance(state, dict) or set(state) != set(self.state()):
            raise ValueError("invalid reporting progress")
        for key, value in state.items():
            utils.require_integer(value, key, 0)
            if value >= 2**64:
                raise ValueError("reporting progress exceeds its cap")
        if not self.lower <= state["bound"] <= self.upper or (
            state["bad"] > 1 or state["epoch"] > self.epochs
        ):
            raise ValueError("reporting state differs from the policy")
        self.__dict__.update(state)
