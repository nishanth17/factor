"""Shared cooperative work, wall-time, CPU-time, and cancellation limits."""

import math
import time
from dataclasses import dataclass, field
from typing import Callable, Optional

from . import utils


class BudgetExhaustedError(Exception):
    """An atomic operation was refused; its input state remains resumable."""


@dataclass
class Budget:
    """Charge work before mutation across every stage and recursive child.

    Work units are algorithmic allowances, not seconds or bigint operations.
    Deadlines are cooperative: a running native integer operation cannot be
    interrupted. Input sizes and chunk sizes therefore also have finite caps.
    A resumed run includes previously consumed work, wall time, and CPU time.
    """

    work_limit: int = 2_000_000
    seconds: Optional[float] = 30.0
    cpu_seconds: Optional[float] = 30.0
    cancelled: Optional[Callable[[], bool]] = None
    used: int = 0
    prior_wall: float = 0.0
    prior_cpu: float = 0.0
    reason: Optional[str] = None
    _wall_start: float = field(default_factory=time.monotonic, repr=False)
    _cpu_start: float = field(default_factory=time.process_time, repr=False)

    def __post_init__(self):
        """Reject invalid limits before any work or clock comparison."""
        utils.require_integer(self.work_limit, "work_limit", 0)
        utils.require_integer(self.used, "used", 0)
        for name in ("seconds", "cpu_seconds", "prior_wall", "prior_cpu"):
            value = getattr(self, name)
            if value is not None and (
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not math.isfinite(value)
                or value < 0
            ):
                raise ValueError(f"{name} must be finite and nonnegative")
        if self.used > self.work_limit:
            raise ValueError("work_limit is below work already consumed")

    @property
    def wall_used(self):
        """Total active-run wall seconds; time spent paused is excluded."""
        return self.prior_wall + time.monotonic() - self._wall_start

    @property
    def cpu_used(self):
        """Total process CPU seconds across all resumptions."""
        return self.prior_cpu + time.process_time() - self._cpu_start

    def consume(self, amount=1):
        """Reserve an atomic action or raise without charging it."""
        utils.require_integer(amount, "amount", 0)
        reason = None
        if self.cancelled is not None and self.cancelled():
            reason = "cancelled"
        elif self.seconds is not None and self.wall_used >= self.seconds:
            reason = "wall_limit"
        elif (
            self.cpu_seconds is not None and self.cpu_used >= self.cpu_seconds
        ):
            reason = "cpu_limit"
        elif amount > self.work_limit - self.used:
            reason = "work_limit"
        if reason is not None:
            self.reason = reason
            raise BudgetExhaustedError(reason)
        self.used += amount
