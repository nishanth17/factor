"""One observed resource owner shared by stages and recursive children."""

import itertools
import math
import os
import sys
import threading
import time
from functools import wraps

from ..common import utils
from .resources import (
    DEFAULT_RSS,
    RESOURCE_VERSION,
    RemovedResourceLimitError,
    ResourceProviderError,
    default_provider,
    validate_seconds,
)

_PROCESS_REGISTRATIONS = itertools.count()


def read_clock(clock, name):
    value = clock()
    try:
        valid = type(value) in (int, float) and math.isfinite(value)
    except OverflowError:
        valid = False
    if not valid:
        raise ResourceProviderError("invalid " + name + " clock observation")
    return value


def synchronized(function):
    """Serialize parent polling with an owned executor's reaping thread."""

    @wraps(function)
    def call(self, *args, **kwargs):
        with self._accounting_lock:
            return function(self, *args, **kwargs)

    return call


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
            if name.startswith("prior_") and value is None:
                raise ValueError(name + " must be finite and nonnegative")
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
        self._history_work = used
        self._history_wall = prior_wall
        self._history_cpu = prior_cpu
        self.reason = self.stop = None
        self.resource_events = []
        self._accounting_lock = threading.RLock()
        self._wall_clock = wall_clock or time.monotonic
        self._cpu_clock = cpu_clock or time.process_time
        self.provider = provider
        self.roots = {os.getpid(): None}
        self._process_bindings = {}
        self._process_identities = {}
        self._live_identities = {}
        self._children = {}
        self._child_cpu = 0.0
        self._last_sample = None
        self.observed_rss = 0
        self._active = False
        if self.provider is None and (
            rss_bytes is not None or sys.platform == "darwin"
        ):
            self.provider = default_provider()
        if rss_bytes is DEFAULT_RSS:
            physical = self.provider.physical_memory()
            if type(physical) is not int or physical <= 0:
                raise ResourceProviderError("invalid physical memory size")
        self.rss_bytes = (
            4 * physical // 5 if rss_bytes is DEFAULT_RSS else rss_bytes
        )
        self.start()

    @synchronized
    def start(self):
        """Begin an active interval, including retained-context startup."""
        if self._active:
            return self
        self._wall_start = read_clock(self._wall_clock, "wall")
        self._cpu_start = read_clock(self._cpu_clock, "CPU")
        self._active = True
        self._last_sample = None
        self._sample(baseline=True)
        return self

    @synchronized
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
        elapsed = (
            read_clock(self._wall_clock, "wall") - self._wall_start
            if self._active
            else 0
        )
        if elapsed < 0:
            raise ResourceProviderError("wall clock moved backwards")
        return self.prior_wall + elapsed

    @property
    def cpu_used(self):
        elapsed = (
            read_clock(self._cpu_clock, "CPU") - self._cpu_start
            if self._active
            else 0
        )
        if elapsed < 0:
            raise ResourceProviderError("CPU clock moved backwards")
        return self.prior_cpu + self._child_cpu + elapsed

    def register(self, pid):
        """Include an owned helper/worker even after it is reparented."""
        self.roots[pid] = None

    @synchronized
    def observe_child(self, identity, cpu_total, *, baseline=False):
        """Reconcile OS samples, acknowledgments and final wait4 totals once.

        identity is (pid, birth), never a reusable PID or a packet identity.
        Reports use process-lifetime CPU. Out-of-order smaller observations
        cannot erase charges; a duplicate acknowledgment cannot add charges.
        """
        validate_seconds(cpu_total, "child CPU")
        if cpu_total is None:
            raise ResourceProviderError("invalid observed child CPU")
        previous = self._children.get(identity)
        if baseline or previous is None:
            previous = cpu_total if baseline else 0.0
        current = max(previous, cpu_total)
        self._child_cpu += current - previous
        self._children[identity] = current

    @synchronized
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
        self._live_identities = {
            sample.pid: (sample.pid, sample.birth) for sample in samples
        }
        for sample in samples:
            if type(sample.rss) is not int or sample.rss < 0:
                raise ResourceProviderError("invalid observed resident size")
            if sample.pid in self.roots:
                self.roots[sample.pid] = sample.birth
            binding = self._process_bindings.get(sample.pid)
            identity = (sample.pid, sample.birth)
            if binding is not None and binding[1] != identity:
                # A child may finish exec before libproc first observes it.
                # Move its report counter to the now-known birth identity;
                # otherwise report -> sample -> exit could bill it twice.
                previous = self._children.pop(binding[1], 0.0)
                self._children[identity] = max(
                    previous, self._children.get(identity, 0.0)
                )
                binding[1] = identity
                self._process_identities[sample.pid] = identity
        self.observed_rss = sum(sample.rss for sample in samples)
        for sample in samples:
            if sample.pid != os.getpid():
                self.observe_child(
                    (sample.pid, sample.birth), sample.cpu, baseline=baseline
                )
        self._last_sample = now

    def _check(self, *, publication=False):
        if self.cancelled is not None and self.cancelled():
            self.reason, self.stop = "cancelled", None
            raise BudgetExhaustedError(self.reason)
        wall = self.wall_used
        if self.seconds is not None and wall >= self.seconds:
            self.reason = "wall_limit"
            self.stop = dict(
                reason=self.reason,
                limit="wall",
                observed=wall,
                limit_value=self.seconds,
                scope=self.scope,
            )
            raise BudgetExhaustedError(self.reason, self.stop)
        self._sample(force=publication)
        for name, observed, limit in (
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
        self.reason = self.stop = None

    def consume(self, amount=1):
        utils.require_integer(amount, "amount", 0)
        self._check()
        self.used += amount

    @synchronized
    def record_incurred_work(self, amount):
        """Record completed foreign work before any cancellable publication poll."""
        utils.require_integer(amount, "amount", 0)
        self.used += amount

    def publication(self):
        """Force a current tree sample before atomic state publication."""
        self._check(publication=True)

    def _consume_reserved(self, amount, reservation):
        utils.require_integer(amount, "amount", 0)
        self._check()
        reservation.check(amount, self.wall_used, self.cpu_used)
        self.used += amount

    @synchronized
    def restore_charges(self, *, work, wall, cpu):
        """Restore saved history without erasing a retained owner's charges."""
        utils.require_integer(work, "saved work", 0)
        validate_seconds(wall, "saved wall")
        validate_seconds(cpu, "saved CPU")
        if wall is None or cpu is None:
            raise ValueError("saved charges must be finite and nonnegative")
        self.used += max(0, work - self._history_work)
        self.prior_wall += max(0, wall - self._history_wall)
        self.prior_cpu += max(0, cpu - self._history_cpu)
        self.record_charges(work=work, wall=wall, cpu=cpu)

    @synchronized
    def record_charges(self, *, work, wall, cpu):
        """Remember an emitted/restored prefix without a durable ledger."""
        self._history_work = max(self._history_work, work)
        self._history_wall = max(self._history_wall, wall)
        self._history_cpu = max(self._history_cpu, cpu)

    @synchronized
    def charges(self):
        values = dict(work=self.used, wall=self.wall_used, cpu=self.cpu_used)
        self.record_charges(**values)
        return dict(resource_version=RESOURCE_VERSION, **values)

    def acknowledge(self, pid, cpu_total):
        """Reconcile a process-lifetime acknowledgment with OS samples."""
        self._sample(force=True)
        identity = self._process_identities.get(pid)
        if identity is None:
            identity = self._live_identities.get(pid)
        if identity is None:
            raise RuntimeError("CPU report from an unregistered process")
        self.observe_child(identity, cpu_total)

    @synchronized
    def track_process(self, process):
        """Retain final wait4 usage through this owned Popen's sole reaper.

        PyPy 3.11 Popen.wait delegates reaping to poll. This instance adapter
        uses the same flags, status conversion and returncode cache as poll,
        retaining the final CPU that waitpid discards. Rebinding a retained
        process changes only its current accounting owner; no global process
        or multiprocessing runtime is modified.
        """
        popen = process._popen
        existing = getattr(popen, "_factor_accounting", None)
        self.register(process.pid)
        self._sample(force=True)
        identity = (
            existing[1]
            if existing is not None
            else self._live_identities.get(
                process.pid,
                (process.pid, ("unobserved", next(_PROCESS_REGISTRATIONS))),
            )
        )
        if existing is not None:
            previous_owner = existing[0]
            if previous_owner is not self:
                previous_owner.roots.pop(process.pid, None)
                previous_owner._process_bindings.pop(process.pid, None)
            if popen.returncode is not None:
                # An idle exit predates this new active interval.
                self._children.setdefault(
                    identity, previous_owner._children.get(identity, 0.0)
                )
            existing[0], existing[1] = self, identity
            binding = existing
        else:
            self._children.setdefault(identity, 0.0)
            binding = [self, identity]
            popen._factor_accounting = binding
        self._process_identities[process.pid] = identity
        self._process_bindings[process.pid] = binding
        if popen.returncode is not None:
            self.roots.pop(process.pid, None)
        if existing is not None:
            return

        reap_lock = threading.RLock()

        def poll(flag=os.WNOHANG):
            with reap_lock:
                if popen.returncode is not None:
                    return popen.returncode
                owner, key = binding
                peek = getattr(owner.provider, "peek_exit", None)
                try:
                    final = None if peek is None else peek(popen.pid, flag)
                    if peek is not None and final is None:
                        return None
                    pid, status, usage = os.wait4(popen.pid, flag)
                except OSError as error:
                    if error.errno == 10:
                        return popen.returncode
                    raise
                if pid == popen.pid:
                    with owner._accounting_lock:
                        if final is not None:
                            identity = (final.pid, final.birth)
                            if identity != key:
                                previous = owner._children.pop(key, 0.0)
                                owner._children[identity] = max(
                                    previous,
                                    owner._children.get(identity, 0.0),
                                )
                                binding[1] = key = identity
                            total = final.cpu
                        else:
                            total = usage.ru_utime + usage.ru_stime
                        owner.observe_child(
                            key, total, baseline=not owner._active
                        )
                        owner.roots.pop(popen.pid, None)
                        owner._process_bindings.pop(popen.pid, None)
                        popen.returncode = os.waitstatus_to_exitcode(status)
                return popen.returncode

        popen.poll = poll
