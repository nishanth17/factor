"""Opt-in bounded independent ECM curves with explicit spawn ownership."""

import multiprocessing
import os
import socket
import sys
import threading
import time
from dataclasses import asdict, dataclass, replace

from .worker_contract import (
    VERSION,
    Decoder,
    EngagementRule,
    LeaseRefusedError,
    Ledger,
    Manifest,
    Search,
    WorkerSettings,
    canonical,
    digest,
    finite,
    frame,
    integer,
)


@dataclass(frozen=True)
class WorkerRun:
    divisor: int | None
    reason: str
    checkpoint: dict
    stats: dict


def arithmetic_options(config, settings):
    """Keep full tier eligibility; remove only unrelated engine ownership."""
    selected = replace(
        config,
        siqs=None,
        sss=None,
        allocation=None,
        memory_bytes=settings.worker_memory,
    )
    return asdict(selected)


def _base_budget(budget):
    return getattr(budget, "budget", budget)


def _allowances(budget):
    """Child CPU/work cannot borrow protected relation-engine reservations."""
    base = _base_budget(budget)
    work, cpu = base.work_limit - base.used, base.cpu_seconds
    if cpu is None or base.seconds is None:
        raise ValueError("workers require finite parent wall and CPU limits")
    cpu -= base.cpu_used
    wall = base.seconds - base.wall_used
    policy = getattr(budget, "policy", None)
    if policy is not None:
        if policy.pretest_work is not None:
            work = min(work, policy.pretest_work - base.used)
        if policy.pretest_cpu_seconds is not None:
            cpu = min(cpu, policy.pretest_cpu_seconds - base.cpu_used)
        if policy.pretest_seconds is not None:
            wall = min(wall, policy.pretest_seconds - base.wall_used)
        if budget.fallback:
            work = min(
                work, base.work_limit - base.used - policy.fallback_work
            )
            cpu = min(
                cpu,
                base.cpu_seconds - base.cpu_used - policy.fallback_cpu_seconds,
            )
            wall = min(
                wall, base.seconds - base.wall_used - policy.fallback_seconds
            )
    return max(0, work), max(0, cpu), max(0, wall)


class ECMExecutor:
    """Explicit context; lazy workers, one active caller, explicit close.

    A rule without calibration uses equivalent-manifest serial execution.
    force_workers is for diagnostics, never an implicit engagement policy.
    Lost reports retain full grants; no automatic in-call retry/replacement.
    """

    def __init__(
        self, *, seed=0, settings=None, rule=None, force_workers=False
    ):
        integer(seed, "method seed", 0, 2**256 - 1)
        if type(force_workers) is not bool:
            raise ValueError("force_workers must be Boolean")
        self.seed = seed
        self.settings = settings or WorkerSettings()
        self.rule = rule or EngagementRule()
        self.force_workers = force_workers
        self.slots = []
        self.closed = False
        self.lock = threading.Lock()
        self.serial_engine = None
        self.last_stats = {}
        self.stats_history = []
        self.lifetime_stats = {}
        self.identity = digest(
            [
                VERSION,
                seed,
                asdict(self.settings),
                asdict(self.rule),
                force_workers,
            ]
        )

    def _record(self, stats):
        self.last_stats = stats
        self.stats_history = (self.stats_history + [stats])[-64:]
        for key in (
            "startup_seconds",
            "cold_pool_start_seconds",
            "charged_work",
            "charged_child_cpu",
            "reported_child_cpu",
            "parent_cpu",
            "drain_seconds",
            "elapsed",
            "lost_reports",
        ):
            self.lifetime_stats[key] = self.lifetime_stats.get(
                key, 0
            ) + stats.get(key, 0)
        self.lifetime_stats["worker_hwm_sum_peak"] = max(
            self.lifetime_stats.get("worker_hwm_sum_peak", 0),
            stats.get("worker_peak_rss_sum", 0),
        )

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def ready(self, manifest):
        return len(self.slots) == self.settings.workers and all(
            s["schedule"] == manifest.schedule and s["process"].is_alive()
            for s in self.slots
        )

    def _send(self, slot, message):
        slot["channel"].sendall(frame(message, self.settings.frame_bytes))

    def _read(self, slot):
        from .ecm_worker import read_message

        return read_message(slot["channel"], slot["decoder"])

    def _destroy(self):
        """Only owned children/channels; every join has a finite timeout."""
        errors = []
        for slot in self.slots:
            try:
                self._send(slot, dict(kind="close"))
            except (OSError, ValueError):
                pass
        for slot in self.slots:
            slot["survivor"] = False
            process = slot["process"]
            process.join(timeout=self.settings.drain_seconds)
            if process.is_alive():
                process.terminate()
                process.join(timeout=self.settings.kill_seconds)
            if process.is_alive():
                process.kill()
                process.join(timeout=self.settings.kill_seconds)
            slot["channel"].close()
            if process.is_alive():
                errors.append(process.pid)
                slot["survivor"] = True
            else:
                process.close()
        self.slots = [s for s in self.slots if s.get("survivor", False)]
        if errors:
            raise RuntimeError("owned workers did not exit: " + repr(errors))

    def close(self):
        if not self.lock.acquire(blocking=False):
            raise ValueError("cannot close an active ECM call")
        try:
            started = time.monotonic()
            self._destroy()
            self.serial_engine = None
            self.closed = True
            self.last_stats["close_seconds"] = time.monotonic() - started
            self.lifetime_stats["close_seconds"] = (
                self.lifetime_stats.get("close_seconds", 0)
                + self.last_stats["close_seconds"]
            )
        finally:
            self.lock.release()

    def _spawn(self, ledger, deadline, poll):
        from .ecm_worker import worker_main

        if threading.current_thread() is not threading.main_thread():
            raise ValueError("spawn preparation requires the main thread")
        context = multiprocessing.get_context("spawn")
        if context.get_start_method() != "spawn":
            raise RuntimeError("explicit spawn context unavailable")
        for index in range(self.settings.workers):
            poll()
            startup = ledger.reserve(
                "startup",
                self.settings.startup_work,
                self.settings.startup_seconds,
            )
            shutdown = ledger.reserve(
                "shutdown", 0, self.settings.shutdown_cpu
            )
            # Pay unreported teardown now, including retained-pool calls.
            ledger.settle(shutdown, "shutdown", 0, self.settings.shutdown_cpu)
            parent, child = socket.socketpair()
            parent.settimeout(0.01)
            process = context.Process(
                target=worker_main,
                args=(
                    child,
                    self.settings.frame_bytes,
                    self.settings.startup_work,
                ),
                daemon=True,
            )
            names = (
                "OMP_NUM_THREADS",
                "OPENBLAS_NUM_THREADS",
                "MKL_NUM_THREADS",
                "VECLIB_MAXIMUM_THREADS",
                "NUMEXPR_NUM_THREADS",
            )
            previous = {name: os.environ.get(name) for name in names}
            try:
                for name in names:
                    os.environ[name] = "1"
                process.start()
            except BaseException:
                parent.close()
                child.close()
                ledger.lose(startup)
                raise
            finally:
                for name, value in previous.items():
                    if value is None:
                        os.environ.pop(name, None)
                    else:
                        os.environ[name] = value
            child.close()
            slot = dict(
                process=process,
                channel=parent,
                decoder=Decoder(self.settings.frame_bytes),
                schedule=None,
                active=None,
                rss=0,
                reported_cpu=0.0,
                index=index,
            )
            self.slots.append(slot)
            end = min(
                deadline, time.monotonic() + self.settings.startup_seconds
            )
            try:
                while True:
                    poll()
                    if time.monotonic() >= end:
                        raise LeaseRefusedError("startup_timeout")
                    message = self._read(slot)
                    if message is not None:
                        if message.get("kind") != "ready":
                            raise ValueError("missing worker ready report")
                        integer(message["rss"], "startup RSS")
                        ledger.settle(
                            startup, "startup", message["work"], message["cpu"]
                        )
                        slot["rss"] = message["rss"]
                        slot["reported_cpu"] = message["cpu"]
                        break
                    if not process.is_alive():
                        raise EOFError("worker exited during startup")
                    time.sleep(0.001)
            except BaseException:
                if startup in ledger.active:
                    ledger.lose(startup)
                raise

    def start(self, *, budget):
        """Measure cold pool readiness separately; no curve or plan is run."""
        if sys.implementation.name != "pypy" or sys.version_info[:2] != (
            3,
            11,
        ):
            raise RuntimeError("PyPy implementing Python 3.11 required")
        if self.closed or not self.lock.acquire(blocking=False):
            raise ValueError("ECM context closed or active")
        base = _base_budget(budget)
        started, cpu_start = time.monotonic(), time.process_time()
        ledger = None
        reported_before = sum(s["reported_cpu"] for s in self.slots)
        try:
            work, cpu, wall = _allowances(budget)
            ledger = Ledger(work, cpu)
            deadline = started + wall

            def poll():
                budget.consume(0)
                if time.monotonic() >= deadline:
                    raise LeaseRefusedError("wall_limit")
                ledger.charge(parent_cpu=time.process_time() - cpu_start)

            poll()
            if not self.slots:
                self._spawn(ledger, deadline, poll)
            poll()
        except BaseException:
            self._destroy()
            raise
        finally:
            if ledger is not None:
                for token in list(ledger.active):
                    ledger.lose(token)
                base.used += ledger.work
                base.prior_cpu += ledger.cpu
            self.lock.release()
        self.last_stats = dict(
            cold_pool_start_seconds=time.monotonic() - started,
            parent_cpu=time.process_time() - cpu_start,
            reported_child_cpu=sum(s["reported_cpu"] for s in self.slots)
            - reported_before,
            charged_child_cpu=ledger.cpu,
            owned_bytes=self.settings.owned_bytes,
        )
        self._record(self.last_stats)
        return dict(self.last_stats)

    def run(
        self,
        n,
        *,
        config,
        budget,
        lineage="root",
        checkpoint=None,
        fixed_work=False,
    ):
        """Return a checked split and a drained, charged checkpoint."""
        if sys.implementation.name != "pypy" or sys.version_info[:2] != (
            3,
            11,
        ):
            raise RuntimeError(
                "ECM executor requires PyPy implementing Python 3.11"
            )
        if self.closed or not self.lock.acquire(blocking=False):
            raise ValueError("ECM context closed or already active")
        try:
            return self._run(
                n, config, budget, lineage, checkpoint, fixed_work
            )
        finally:
            self.lock.release()

    def _run(self, n, config, budget, lineage, checkpoint, fixed_work):
        from .budget import BudgetExhaustedError
        from .ecm_worker import CurveEngine

        call_started = time.monotonic()
        base = _base_budget(budget)
        if base.cpu_seconds is None or base.seconds is None:
            raise ValueError("finite wall/CPU grants required")
        manifest = Manifest(
            n, self.seed, lineage, arithmetic_options(config, self.settings)
        )
        search = Search(manifest, self.settings, fixed_work=fixed_work)
        if checkpoint is not None:
            payload = dict(checkpoint)
            overhead = payload.pop("serialization_overhead")
            overhead_hash = payload.pop("overhead_sha256")
            if digest(overhead) != overhead_hash or set(overhead) != {
                "wall",
                "cpu",
            }:
                raise ValueError("corrupt serialization overhead")
            for key, value in overhead.items():
                finite(value, "serialization " + key)
            checksum = payload.pop("sha256")
            if digest(payload) != checksum:
                raise ValueError("corrupt cumulative checkpoint")
            if checkpoint["executor"] != self.identity:
                raise ValueError("executor/rule changed on resume")
            saved = dict(checkpoint["budget"])
            saved["cpu"] += overhead["cpu"]
            saved["wall"] += overhead["wall"]
            for key in ("work", "cpu", "wall"):
                finite(saved[key], "saved " + key)
            integer(saved["work"], "saved work")
            if (
                saved["work"] > base.work_limit
                or saved["cpu"] > base.cpu_seconds
                or saved["wall"] > base.seconds
            ):
                raise ValueError("resume grant below prior charges")
            # A new Budget must retain all prior active-run charges.
            if base.used == 0 and base.prior_cpu == 0 and base.prior_wall == 0:
                base.used = saved["work"]
                base.prior_cpu, base.prior_wall = saved["cpu"], saved["wall"]
            elif (
                base.used < saved["work"]
                or base.cpu_used < saved["cpu"]
                or base.wall_used < saved["wall"]
            ):
                raise ValueError("resume would reset cumulative charges")
            search = Search.restore(
                checkpoint["search"],
                manifest,
                self.settings,
                fixed_work=fixed_work,
            )
        if search.runs >= self.settings.max_runs:
            raise LeaseRefusedError("run_limit")
        search.runs += 1
        work, cpu, wall = _allowances(budget)
        output_reserve = min(work, self.settings.checkpoint_bytes)
        ledger = Ledger(work - output_reserve, cpu)
        start_cpu, started = time.process_time(), time.monotonic()
        deadline = started + wall
        reported_before = sum(s["reported_cpu"] for s in self.slots)
        ready = self.ready(manifest)
        decision_started = time.monotonic()
        remaining = tuple(
            sum(
                position not in search.completed and job[0] == tier
                for position, job in enumerate(manifest.jobs)
            )
            for tier in range(len(manifest.tiers))
        )
        engaged = self.force_workers or self.rule.engages(
            manifest,
            ready=ready,
            remaining=remaining,
        )
        funded = ready or (
            cpu
            > self.settings.workers
            * (self.settings.startup_seconds + self.settings.shutdown_cpu)
            and work - output_reserve
            > self.settings.workers
            * (self.settings.startup_work + 2 * self.settings.frame_bytes)
        )
        if not self.force_workers and not funded:
            engaged = False
        stats = dict(
            decision_seconds=time.monotonic() - decision_started,
            planning_seconds=decision_started - call_started,
            worker_dispatches=0,
            funded=funded,
            version=VERSION,
            manifest=manifest.identity,
            schedule=manifest.schedule,
            engaged=engaged,
            pool_ready=ready,
            owned_bytes=self.settings.owned_bytes,
            startup_seconds=0.0,
            drain_seconds=0.0,
            worker_peak_rss_sum=0,
            completed_before=len(search.completed),
        )
        reason = "curves_exhausted"
        stopped = False

        def poll():
            budget.consume(0)
            if time.monotonic() >= deadline:
                raise LeaseRefusedError("wall_limit")
            ledger.charge(parent_cpu=time.process_time() - start_cpu)

        def receive(slot, message, publish):
            if message.get("kind") != "result" or slot["active"] is None:
                raise ValueError("unsolicited worker result")
            position, token = slot["active"]
            identity, grant_work, grant_cpu = ledger.active[token]
            used = integer(message["work"], "worker work")
            consumed = finite(message["cpu"], "worker CPU")
            integer(message["rss"], "worker RSS")
            if used > grant_work or consumed > grant_cpu:
                # Preserve demonstrated overruns instead of silently clipping.
                ledger.work += max(0, used - grant_work)
                ledger.cpu += max(0, consumed - grant_cpu)
                ledger.receipts[token - 1].update(
                    overrun_work=used, overrun_cpu=consumed
                )
                raise ValueError("worker exceeded reservation")
            if (
                message.get("reason") == "complete"
                or message.get("sigma") is not None
            ) and message.get("sigma") != manifest.assignment(position)[
                "sigma"
            ]:
                raise ValueError("curve sigma differs from manifest")
            search.receive(position, token, message, publish=publish)
            ledger.settle(token, identity, used, consumed)
            ledger.receipts[token - 1].update(
                reason=message["reason"], failure=message.get("failure")
            )
            slot["active"] = None
            slot["schedule"] = message.get("prepared_schedule")
            slot["rss"] = max(slot["rss"], message["rss"])
            slot["reported_cpu"] += consumed
            if message["reason"] != "complete":
                raise LeaseRefusedError(message["reason"])

        try:
            poll()
            ledger.charge(
                len(canonical(manifest.options)) + 64 * len(manifest.jobs)
            )
            if not search.retired and engaged and not self.slots:
                before = time.monotonic()
                self._spawn(ledger, deadline, poll)
                stats["startup_seconds"] = time.monotonic() - before
            while not search.retired and len(search.completed) < len(
                manifest.jobs
            ):
                poll()
                if engaged:
                    for slot in self.slots:
                        if slot["active"] is not None:
                            continue
                        position = search.next(slot["index"])
                        if position is None:
                            continue
                        transport_work = 2 * self.settings.frame_bytes
                        if (
                            ledger.work_limit - ledger.work <= transport_work
                            or ledger.cpu_limit
                            - ledger.cpu
                            - ledger.parent_cpu
                            <= 0
                        ):
                            if search.active:
                                break
                            raise LeaseRefusedError("lease_refused")
                        ledger.charge(transport_work)
                        grant_work = min(
                            self.settings.assignment_work,
                            ledger.work_limit - ledger.work,
                        )
                        grant_cpu = min(
                            self.settings.assignment_cpu,
                            ledger.cpu_limit - ledger.cpu - ledger.parent_cpu,
                        )
                        if grant_work <= 0 or grant_cpu <= 0:
                            raise LeaseRefusedError("lease_refused")
                        token = ledger.reserve(
                            manifest.assignment(position)["id"],
                            grant_work,
                            grant_cpu,
                        )
                        search.dispatch(position, token)
                        slot["active"] = (position, token)
                        stats["worker_dispatches"] += 1
                        task = self._task(
                            manifest,
                            search,
                            position,
                            token,
                            deadline,
                            grant_work,
                            grant_cpu,
                        )
                        self._send(slot, task)
                    for slot in self.slots:
                        message = self._read(slot)
                        if message is not None:
                            poll()
                            receive(slot, message, publish=True)
                        if not slot["process"].is_alive():
                            raise EOFError(
                                "worker exited without acknowledgement"
                            )
                    time.sleep(0.001)
                else:
                    position = search.next(0)
                    # Equivalent-manifest serial ignores static slot ownership.
                    if position is None:
                        pending = [
                            i
                            for i in range(len(search.attempts))
                            if i not in search.completed
                        ]
                        position = pending[0]
                    grant_work = min(
                        self.settings.assignment_work,
                        ledger.work_limit - ledger.work,
                    )
                    grant_cpu = min(
                        self.settings.assignment_cpu,
                        ledger.cpu_limit - ledger.parent_cpu,
                    )
                    if grant_work <= 0 or grant_cpu <= 0:
                        raise LeaseRefusedError("lease_refused")
                    token = ledger.reserve(
                        manifest.assignment(position)["id"], grant_work, 0
                    )
                    search.dispatch(position, token)
                    if self.serial_engine is None:
                        self.serial_engine = CurveEngine()
                    result = self.serial_engine.execute(
                        self._task(
                            manifest,
                            search,
                            position,
                            token,
                            deadline,
                            grant_work,
                            grant_cpu,
                        ),
                        cancelled=lambda: bool(
                            base.cancelled and base.cancelled()
                        ),
                    )
                    poll()
                    search.receive(position, token, result)
                    ledger.settle(token, result["id"], result["work"], 0)
                    if result["reason"] != "complete":
                        raise LeaseRefusedError(result["reason"])
        except (BudgetExhaustedError, LeaseRefusedError) as error:
            reason, stopped = str(error), True
        except (
            OSError,
            EOFError,
            ValueError,
            MemoryError,
            KeyError,
            TypeError,
        ) as error:
            reason, stopped = "worker_failure", True
            stats["failure"] = type(error).__name__ + ": " + str(error)[:256]
        except BaseException:
            stopped = True
            raise
        finally:
            before = time.monotonic()
            for slot in self.slots:
                if slot["active"] is not None:
                    try:
                        self._send(slot, dict(kind="cancel"))
                    except OSError:
                        stopped = True
            drain_end = time.monotonic() + self.settings.drain_seconds
            while (
                any(s["active"] is not None for s in self.slots)
                and time.monotonic() < drain_end
            ):
                for slot in self.slots:
                    if slot["active"] is None:
                        continue
                    try:
                        message = self._read(slot)
                        if message is not None:
                            receive(slot, message, publish=False)
                    except (
                        OSError,
                        EOFError,
                        ValueError,
                        KeyError,
                        TypeError,
                        LeaseRefusedError,
                    ):
                        stopped = True
                time.sleep(0.001)
            stats["worker_peak_rss_sum"] = sum(s["rss"] for s in self.slots)
            for slot in self.slots:
                if slot["active"] is not None:
                    stopped = True
            stats["reported_child_cpu"] = (
                sum(s["reported_cpu"] for s in self.slots) - reported_before
            )
            if stopped:
                try:
                    self._destroy()
                except RuntimeError as error:
                    stats["cleanup_failure"] = str(error)
            for position, (token, _) in list(search.active.items()):
                if token in ledger.active:
                    ledger.lose(token)
                del search.active[position]
            # Any pre-dispatch/startup reservation without a trustworthy report
            # remains spent. No missing report can manufacture a refund.
            for token in list(ledger.active):
                ledger.lose(token)
            base.used += ledger.work
            base.prior_cpu += ledger.cpu
            stats.update(
                drain_seconds=time.monotonic() - before,
                charged_work=ledger.work,
                charged_child_cpu=ledger.cpu,
                parent_cpu=time.process_time() - start_cpu,
                lost_reports=ledger.lost,
                elapsed=time.monotonic() - started,
            )
        stats["engagement_selected"] = stats["engaged"]
        stats["engaged"] = stats["worker_dispatches"] > 0
        stats["total_call_seconds"] = time.monotonic() - call_started
        stats["completed_after"] = len(search.completed)
        sigmas = {}
        for position, message in search.completed.items():
            sigmas.setdefault(message.get("sigma"), []).append(position)
        stats["observed_sigma_collisions"] = [
            positions
            for sigma, positions in sigmas.items()
            if sigma is not None and len(positions) > 1
        ]
        snapshot = dict(
            executor=self.identity,
            leases=([] if checkpoint is None else checkpoint["leases"])
            + ledger.receipts,
            search=search.snapshot(),
            budget=dict(
                work=base.used, cpu=base.cpu_used, wall=base.wall_used
            ),
        )
        if len(canonical(snapshot)) > self.settings.checkpoint_bytes:
            raise MemoryError("complete checkpoint exceeds encoding allowance")
        # Serialization is finite work even when search exhausted its grant.
        # Reserve output before work leases by withholding checkpoint_bytes.
        # Fixed header reserve covers the two integrity/overhead envelopes.
        output_work = len(canonical(snapshot)) + 512
        if output_work > output_reserve:
            raise LeaseRefusedError("checkpoint_work_limit")
        base.used += output_work
        snapshot["budget"]["work"] = base.used
        snapshot["sha256"] = digest(snapshot)
        self._record(stats)
        overhead = dict(
            cpu=base.cpu_used - snapshot["budget"]["cpu"],
            wall=base.wall_used - snapshot["budget"]["wall"],
        )
        # Match the incumbent checkpoint convention: encoding is paid on a
        # fresh-Budget resume too, rather than resetting active-run clocks.
        snapshot["serialization_overhead"] = overhead
        snapshot["overhead_sha256"] = digest(overhead)
        if stats.get("cleanup_failure"):
            raise RuntimeError(stats["cleanup_failure"])
        return WorkerRun(
            search.factor,
            "factor_found" if search.factor is not None else reason,
            snapshot,
            stats,
        )

    def _task(self, manifest, search, position, token, deadline, work, cpu):
        return dict(
            kind="curve",
            assignment=manifest.assignment(position),
            schedule=manifest.schedule,
            options=manifest.options,
            token=token,
            attempt=search.attempts[position],
            work=work,
            cpu=cpu,
            deadline=deadline,
        )
