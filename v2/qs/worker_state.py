"""One identity-bound, disposable current-QS context per owned worker."""

import time
from concurrent.futures import wait
from dataclasses import asdict

from ..execution.budget import Budget, BudgetExhaustedError
from .factor_base import FactorBase
from .worker_packet import digest


def context_key(payload):
    record = payload["base"]
    return digest(
        [
            "qs-worker-state-v1",
            payload["job"],
            payload["execution"],
            asdict(payload["config"]),
            payload["settings"].record(),
            [
                record[0],
                record[1],
                record[2],
                [[e.prime, list(e.square_roots)] for e in record[3]],
            ],
        ]
    )


class WorkerState:
    """Keep bounded derived scratch, independent of checkpoints and stores."""

    def __init__(self, payload, budget):
        if payload["key"] != context_key(payload):
            raise ValueError("worker context identity mismatch")
        self.key, self.job = payload["key"], payload["job"]
        self.config, self.settings = payload["config"], payload["settings"]
        self.record = payload["base"]
        budget.consume(payload["work"])
        self.base = FactorBase(*self.record)
        # Includes immutable entries, indices, config and in-flight metadata.
        self.bytes = (
            self.base.workspace_bytes + 65536 + 512 * len(self.base.entries)
        )
        if self.bytes > self.settings.state_cache_bytes:
            raise MemoryError(
                "immutable worker context exceeds cache allowance"
            )
        budget.consume(0)
        self.family = self.collector = self.recovery = None
        self.primes = None
        self.family_hits = self.collector_hits = 0

    def reset_recovery(self, collector):
        """Retain exact reusable proofs, not packet-dependent split hints."""
        reset_packet_hints(collector)
        collector._split_calls = 0
        if self.recovery is not None:
            collector._divide_double = self.recovery


def reset_packet_hints(collector):
    """Invalidate trajectory-dependent hints without changing split leases."""
    checked = collector._checked
    if checked is not None:
        checked.records.clear()
        checked.canonical.clear()
        checked.canonical_ids.clear()
        checked.canonical_used = 0
        checked.used = 512 * (len(checked.primes) + len(checked.units))
    endpoints = getattr(collector, "_endpoints", None)
    if endpoints is not None:
        endpoints.primes.clear()
        endpoints.seen.clear()
        endpoints.product = 1
    memo = getattr(collector, "_split_memo", None)
    if memo is not None:
        memo.values.clear()
    if getattr(collector.config, "extension_policy", None):
        from .yield_feedback import ReportingController

        collector._reporting = ReportingController(
            collector._candidate_bound,
            max(
                collector._candidate_bound,
                collector.config.large_product_bound,
            ),
            collector.config.threshold_exposure,
            collector.config.threshold_epochs,
        )


def install_context(payload, state=None):
    """Install on every slot before acknowledging the shared barrier."""
    from . import parallel

    process = state is None
    slot = None
    if process:
        state, slot = parallel._PROCESS_STATE
    # Retire old derived objects before constructing a different context.
    state.pop("_qs_context", None)
    budget = Budget(
        work_limit=payload["work"],
        seconds=payload["seconds"],
        cpu_seconds=payload["cpu"],
        cancelled=state["stop"].is_set,
    )
    try:
        context = WorkerState(payload, budget)
        state["_qs_context"] = context
        if process:
            state["barrier"].wait(timeout=payload["seconds"])
        return dict(key=context.key, bytes=context.bytes)
    finally:
        if process:
            state["cpu"][slot] = time.process_time()
            state["rss"][slot] = parallel._rss()


def ensure_context(job, pool):
    """Fund installation before dispatch and require every acknowledgement."""
    base = job.source.base
    context_bytes = base.workspace_bytes + 65536 + 512 * len(base.entries)
    if context_bytes > job.workers.state_cache_bytes:
        raise MemoryError("immutable worker context exceeds cache allowance")
    work = 1 + len(base.entries) * base.bound.bit_length() ** 2
    cpu = job.workers.assignment_cpu
    if job.budget.cpu_seconds is not None:
        cpu = min(
            cpu,
            (
                job.budget.cpu_seconds
                - job.budget.cpu_used
                - job.workers.parent_cpu_reserve
            )
            / pool.workers,
        )
    if cpu <= 0:
        job.budget.reason = "cpu_limit"
        raise BudgetExhaustedError("cpu_limit")
    payload = dict(
        job=job.identity,
        execution=job.execution_identity,
        config=job.config,
        settings=job.workers,
        base=(base.n, base.multiplier, base.bound, base.entries),
        work=work,
        cpu=cpu,
        seconds=job.workers.assignment_seconds,
    )
    payload["key"] = context_key(payload)
    if getattr(pool, "_qs_context_key", None) == payload["key"]:
        return payload["key"]
    job.budget.consume(work * pool.workers)
    job.stats["state_install_work"] += work * pool.workers
    started = time.monotonic()
    try:
        if pool.executor is None:
            reports = [install_context(payload, pool.state)]
        else:
            futures = {
                pool.executor.submit(install_context, payload)
                for _ in range(pool.workers)
            }
            reports = []
            while futures:
                pool.poll(job.budget)
                if time.monotonic() - started > payload["seconds"]:
                    raise TimeoutError(
                        "worker context acknowledgement timed out"
                    )
                done, futures = wait(futures, timeout=0.01)
                reports.extend(f.result() for f in done)
        if len(reports) != pool.workers or any(
            r["key"] != payload["key"]
            or not 0 < r["bytes"] <= job.workers.state_cache_bytes
            for r in reports
        ):
            raise ValueError("worker context acknowledgement mismatch")
    except BaseException:
        # Unknown initialization CPU is conservatively charged; paid setup
        # work is never refunded. Only this pool's processes are retired.
        job.budget.prior_cpu += cpu * pool.workers
        pool.abort()
        raise
    pool._qs_context_key = payload["key"]
    pool._qs_state_settings = job.workers.record()
    job.stats["state_install_seconds"] += time.monotonic() - started
    return payload["key"]


def compact_task(task, key):
    """Send a leased delta after every worker acknowledged its context."""
    return {
        k: v
        for k, v in dict(task, context_key=key).items()
        if k not in ("base", "config", "settings")
    }


def expand_task(task, state=None):
    if "context_key" not in task:
        return task
    if state is None:
        from .parallel import _PROCESS_STATE

        state = _PROCESS_STATE[0]
    context = state.get("_qs_context")
    if (
        context is None
        or context.key != task["context_key"]
        or context.job != task["job"]
    ):
        raise ValueError("uninstalled or stale worker context delta")
    return dict(
        task,
        config=context.config,
        settings=context.settings,
        base=context.record,
        _runtime=context,
    )
