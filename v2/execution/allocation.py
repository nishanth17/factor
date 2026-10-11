"""Explicit finite pretesting and protected relation-engine allowances."""

import math
from dataclasses import dataclass

from ..common import utils

ALLOCATION_VERSION = "cumulative-pretest/reserved-fallback-v1"


@dataclass(frozen=True)
class ECMAllocation:
    """Bound optional searches against the original cumulative ledger.

    None retains the legacy portfolio policy. A pretest ceiling includes all
    prior work, setup, classification and earlier recursive/fallback calls;
    it is never a fresh allowance per child or resume. Campaigns instead run
    the caller's finite tiers, subject to the same optional fallback reserves.
    These are service limits, not estimated factor sizes or success guarantees.
    """

    mode: str
    pretest_work: int | None = None
    pretest_seconds: float | None = None
    pretest_cpu_seconds: float | None = None
    fallback_work: int = 0
    fallback_seconds: float = 0.0
    fallback_cpu_seconds: float = 0.0

    def __post_init__(self):
        if self.mode not in ("pretest", "campaign"):
            raise ValueError("ECM allocation mode must be pretest or campaign")
        utils.require_integer(self.fallback_work, "fallback_work", 0)
        if self.mode == "pretest":
            utils.require_integer(self.pretest_work, "pretest_work", 0)
        elif any(
            value is not None
            for value in (
                self.pretest_work,
                self.pretest_seconds,
                self.pretest_cpu_seconds,
            )
        ):
            raise ValueError(
                "campaigns use finite tiers, not pretest ceilings"
            )

        for name in (
            "pretest_seconds",
            "pretest_cpu_seconds",
            "fallback_seconds",
            "fallback_cpu_seconds",
        ):
            value = getattr(self, name)
            if value is None and name.startswith("pretest"):
                continue
            if (
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not math.isfinite(value)
                or value < 0
            ):
                raise ValueError(f"{name} must be finite and nonnegative")


class HandoffRequiredError(Exception):
    """An optional action would consume the selected engine's allowance."""


class PretestBudget:
    """Refuse before mutation while delegating every charge to one Budget."""

    def __init__(self, budget, policy, *, fallback):
        self.budget = budget
        self.policy = policy
        self.fallback = fallback
        self.needs_wall = policy.pretest_seconds is not None
        self.needs_cpu = policy.pretest_cpu_seconds is not None

    def __getattr__(self, name):
        return getattr(self.budget, name)

    def consume(self, amount=1):
        self.budget._consume_reserved(amount, self)

    def check(self, amount, wall, cpu):
        """Refuse optional work against the base ledger's sampled clocks."""
        policy = self.policy
        if (
            policy.pretest_work is not None
            and amount > policy.pretest_work - self.budget.used
        ):
            raise HandoffRequiredError("pretest_work")
        if (
            policy.pretest_seconds is not None
            and wall >= policy.pretest_seconds
        ):
            raise HandoffRequiredError("pretest_wall")
        if (
            policy.pretest_cpu_seconds is not None
            and cpu >= policy.pretest_cpu_seconds
        ):
            raise HandoffRequiredError("pretest_cpu")

        if self.fallback:
            if (
                amount
                > self.budget.work_limit
                - self.budget.used
                - policy.fallback_work
            ):
                raise HandoffRequiredError("reserved_work")
            if (
                self.budget.seconds is not None
                and self.budget.seconds - wall <= policy.fallback_seconds
            ):
                raise HandoffRequiredError("reserved_wall")
            if (
                self.budget.cpu_seconds is not None
                and self.budget.cpu_seconds - cpu
                <= policy.fallback_cpu_seconds
            ):
                raise HandoffRequiredError("reserved_cpu")


def fallback_refusal(budget, policy):
    """Check a first admission; resumed collection retains its admission."""
    budget.consume(0)
    if budget.work_limit - budget.used < policy.fallback_work:
        return "insufficient_fallback_work"
    if (
        budget.seconds is not None
        and budget.seconds - budget.wall_used < policy.fallback_seconds
    ):
        return "insufficient_fallback_wall"
    if (
        budget.cpu_seconds is not None
        and budget.cpu_seconds - budget.cpu_used < policy.fallback_cpu_seconds
    ):
        return "insufficient_fallback_cpu"
    return None
