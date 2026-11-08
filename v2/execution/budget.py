"""One observed resource owner shared by stages and recursive children."""

import os
import sys
import time

from ..common import utils
from .resources import (
    DEFAULT_RSS,
    RESOURCE_VERSION,
    RemovedResourceLimitError,
    default_provider,
    validate_seconds,
)


class BudgetExhaustedError(Exception):
    """A cooperative boundary stopped before publishing its next action."""

    def __init__(self, reason="cancelled", payload=None):
        self.payload = payload
        super().__init__(reason)


class Budget:
    """Cumulative active wall/CPU and sampled current owned-tree RSS.

    Omitted wall/CPU are unlimited. Omitted RSS is floor(0.8 * physical RAM);
    explicit None disables it. Work is an exact diagnostic counter, not a
    refusal allowance. Checkpoint resumes inherit saved charges; unsaved
    consumption survives only when the caller retains this live owner.
    Paused/idle intervals are excluded by close()/start(). No durable ledger
    is created. The RSS guard can overshoot between cooperative samples.
    """

    def __init__(
        self,
        *,
        seconds=None,
        cpu_seconds=None,
        rss_bytes=DEFAULT_RSS,
        cancelled=None,
        used=0,
        prior_wall=0.0,
        prior_cpu=0.0,
        wall_clock=None,
        cpu_clock=None,
        provider=None,
        scope="factorization",
        **removed,
    ):
        if removed:
            name = next(iter(removed))
            if name == "work_limit":
                raise RemovedResourceLimitError(
                    name, "seconds/cpu_seconds; used remains diagnostic"
                )
            raise TypeError("unexpected Budget argument: " + name)
        for name, value in (
            ("seconds", seconds),
            ("cpu_seconds", cpu_seconds),
            ("prior_wall", prior_wall),
            ("prior_cpu", prior_cpu),
        ):
            validate_seconds(value, name)
        utils.require_integer(used, "used", 0)
        if rss_bytes is not DEFAULT_RSS and rss_bytes is not None:
            utils.require_integer(rss_bytes, "rss_bytes", 0)
        self.seconds, self.cpu_seconds = seconds, cpu_seconds
        self.cancelled, self.scope = cancelled, scope
        self.used, self.prior_wall, self.prior_cpu = (
            used,
            prior_wall,
            prior_cpu,
        )
        self.reason = self.stop = None
        self._wall_clock = wall_clock or time.monotonic
        self._cpu_clock = cpu_clock or time.process_time
        self.provider = provider
        self.roots = {os.getpid()}
        self._children = {}
        self._child_cpu = 0.0
        self._last_sample = None
        self.observed_rss = 0
        self._active = False
        if self.provider is None and (
            rss_bytes is not None or sys.platform == "darwin"
        ):
            self.provider = default_provider()
        self.rss_bytes = (
            4 * self.provider.physical_memory() // 5
            if rss_bytes is DEFAULT_RSS
            else rss_bytes
        )
        self.start()

    def start(self):
        """Begin an active interval, including retained-context startup."""
        if self._active:
            return self
        self._wall_start = self._wall_clock()
        self._cpu_start = self._cpu_clock()
        self._active = True
        self._last_sample = None
        self._sample(baseline=True)
        return self

    def close(self):
        """Settle an active interval after owned children have drained."""
        if self._active:
            self._sample(force=True)
            self.prior_wall = self.wall_used
            self.prior_cpu = self.cpu_used
            self._child_cpu = 0.0
            self._active = False
        return self

    def __enter__(self):
        return self.start()

    def __exit__(self, *exception):
        self.close()

    @property
    def wall_used(self):
        return self.prior_wall + (
            self._wall_clock() - self._wall_start if self._active else 0
        )

    @property
    def cpu_used(self):
        return (
            self.prior_cpu
            + self._child_cpu
            + (self._cpu_clock() - self._cpu_start if self._active else 0)
        )

    def register(self, pid):
        """Include an owned helper/worker even after it is reparented."""
        self.roots.add(pid)

    def observe_child(self, identity, cpu_total, *, baseline=False):
        """Reconcile OS samples, acknowledgments and final wait4 totals once.

        identity is (pid, birth), never a reusable PID or a packet identity.
        Reports use process-lifetime CPU. Out-of-order smaller observations
        cannot erase charges; a duplicate acknowledgment cannot add charges.
        """
        validate_seconds(cpu_total, "child CPU")
        previous = self._children.get(identity)
        if baseline or previous is None:
            previous = cpu_total if baseline else 0.0
        current = max(previous, cpu_total)
        self._child_cpu += current - previous
        self._children[identity] = current

    def _sample(self, *, force=False, baseline=False):
        if self.provider is None:
            return
        now = self._wall_clock()
        if (
            not force
            and self._last_sample is not None
            and now - self._last_sample < 0.05
        ):
            return
        samples = self.provider.sample(self.roots)
        self.observed_rss = sum(sample.rss for sample in samples)
        for sample in samples:
            if sample.pid != os.getpid():
                self.observe_child(
                    (sample.pid, sample.birth), sample.cpu, baseline=baseline
                )
        self._last_sample = now

    def _check(self, *, publication=False):
        self._sample(force=publication)
        if self.cancelled is not None and self.cancelled():
            self.reason, self.stop = "cancelled", None
            raise BudgetExhaustedError(self.reason)
        for name, observed, limit in (
            ("wall", self.wall_used, self.seconds),
            ("cpu", self.cpu_used, self.cpu_seconds),
            ("rss", self.observed_rss, self.rss_bytes),
        ):
            if limit is not None and observed >= limit:
                self.reason = name + "_limit"
                self.stop = dict(
                    reason=self.reason,
                    limit=name,
                    observed=observed,
                    limit_value=limit,
                    scope=self.scope,
                )
                raise BudgetExhaustedError(self.reason, self.stop)

    def consume(self, amount=1):
        utils.require_integer(amount, "amount", 0)
        self._check()
        self.used += amount

    def publication(self):
        """Force a current tree sample before atomic state publication."""
        self._check(publication=True)

    def _consume_reserved(self, amount, reservation):
        utils.require_integer(amount, "amount", 0)
        self._check()
        reservation.check(amount, self.wall_used, self.cpu_used)
        self.used += amount

    def restore_charges(self, *, work, wall, cpu):
        """Restore saved history without erasing a retained owner's charges."""
        utils.require_integer(work, "saved work", 0)
        validate_seconds(wall, "saved wall")
        validate_seconds(cpu, "saved CPU")
        self.used = max(self.used, work)
        self.prior_wall = max(self.prior_wall, wall)
        self.prior_cpu = max(self.prior_cpu, cpu)

    def charges(self):
        return dict(
            resource_version=RESOURCE_VERSION,
            work=self.used,
            wall=self.wall_used,
            cpu=self.cpu_used,
        )

    def acknowledge(self, pid, cpu_total):
        """Reconcile a process-lifetime acknowledgment with OS samples."""
        self._sample(force=True)
        identities = [key for key in self._children if key[0] == pid]
        if not identities:
            raise RuntimeError("CPU report from an unregistered process")
        self.observe_child(identities[-1], cpu_total)

    def track_process(self, process):
        """Retain final wait4 usage through this owned Popen's sole reaper.

        PyPy 3.11 Popen.wait delegates reaping to poll. This instance adapter
        uses the same flags, status conversion and returncode cache as poll,
        retaining the final CPU that waitpid discards. Rebinding a retained
        process changes only its current accounting owner; no global process
        or multiprocessing runtime is modified.
        """
        self.register(process.pid)
        self._sample(force=True)
        identities = [key for key in self._children if key[0] == process.pid]
        identity = identities[-1] if identities else (process.pid, id(process))
        self._children.setdefault(identity, 0.0)
        popen = process._popen
        existing = getattr(popen, "_factor_accounting", None)
        if existing is not None:
            existing[0], existing[1] = self, identity
            return
        binding = [self, identity]
        popen._factor_accounting = binding

        def poll(flag=os.WNOHANG):
            if popen.returncode is None:
                try:
                    pid, status, usage = os.wait4(popen.pid, flag)
                except OSError as error:
                    # Match Popen's pre-creation/ECHILD handling. Other OS
                    # Other OS failures keep their original error identity.
                    if error.errno == 10:
                        return None
                    raise
                if pid == popen.pid:
                    owner, key = binding
                    owner.observe_child(key, usage.ru_utime + usage.ru_stime)
                    popen.returncode = os.waitstatus_to_exitcode(status)
            return popen.returncode

        popen.poll = poll
