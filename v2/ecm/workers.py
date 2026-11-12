"""Opt-in bounded independent ECM curves with explicit spawn ownership."""

import os
import select
import socket
import sys
import threading
import time
from collections import OrderedDict
from dataclasses import asdict, dataclass, replace

from ..execution.resources import (
    RESOURCE_VERSION,
    UnsupportedEncodingError,
    accounted_entry,
    require_resource_version,
)
from ..execution.worker_contract import (
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
    resource_stop: dict | None = None


def arithmetic_options(config, settings):
    """Keep full tier eligibility; remove only unrelated engine ownership."""
    selected = replace(
        config,
        siqs=None,
        sss=None,
        allocation=None,
    )
    options = asdict(selected)
    if selected.ecm_kernel == "python":
        options.pop("ecm_kernel")
    return options


def _base_budget(budget):
    return getattr(budget, "budget", budget)


def _deadline(budget):
    """Root/stage ordering deadline; no worker receives a CPU/work lease."""
    base = _base_budget(budget)
    remaining = (
        None if base.seconds is None else max(0, base.seconds - base.wall_used)
    )
    return None if remaining is None else time.monotonic() + remaining


class ECMExecutor:
    """Explicit context; lazy workers, one active caller, explicit close.

    A rule without calibration uses equivalent-manifest serial execution.
    force_workers is for diagnostics, never an implicit engagement policy.
    Lost reports retain full grants; no automatic in-call retry/replacement.
    """

    def __init__(
        self,
        *,
        seed=0,
        settings=None,
        rule=None,
        force_workers=False,
        metadata_cache=None,
    ):
        integer(seed, "method seed", 0, 2**256 - 1)
        if type(force_workers) is not bool:
            raise ValueError("force_workers must be Boolean")
        self.seed = seed
        self.settings = settings or WorkerSettings()
        if metadata_cache is not None and (
            self.settings.workers != 1 or force_workers
        ):
            raise ValueError(
                "metadata cache requires the serial whole-native owner"
            )
        self.metadata_cache = metadata_cache
        self.rule = rule or EngagementRule()
        self.force_workers = force_workers
        self.slots = []
        self.packet_sequence = 0
        self.helper_funded = False
        self.startup_in_progress = False
        self.cost_identity = digest(
            [VERSION, "cost-policy", asdict(self.settings)]
        )
        self.closed = False
        self.lock = threading.Lock()
        self.serial_engine = None
        self.manifests = OrderedDict()
        self.manifest_bytes = 0
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
        if metadata_cache is not None:
            self.identity = digest(
                [self.identity, "native-metadata-cache", str(metadata_cache)]
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
            "helper_cpu_reserved",
            "preparations",
            "cache_hits",
            "cache_evictions",
            "decode_hits",
            "decode_misses",
            "decode_refusals",
            "stage_two_hits",
            "stage_two_misses",
            "stage_two_refusals",
            "stage_two_build_bytes",
            "manifest_cache_hits",
            "manifest_cache_misses",
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

    def _check_configuration(self):
        actual = digest(
            [
                VERSION,
                self.seed,
                asdict(self.settings),
                asdict(self.rule),
                self.force_workers,
            ]
        )
        if self.metadata_cache is not None:
            actual = digest(
                [actual, "native-metadata-cache", str(self.metadata_cache)]
            )
        if actual != self.identity:
            raise ValueError(
                "executor settings/rule mutated after construction"
            )

    def pool_started(self):
        return len(self.slots) == self.settings.workers and all(
            slot["process"].is_alive() for slot in self.slots
        )

    def ready(self, manifest):
        return self.pool_started() and all(
            manifest.metadata in slot["metadata_cache"] for slot in self.slots
        )

    def _readable(self, timeout):
        active = [slot for slot in self.slots if slot["active"]]
        buffered = [slot for slot in active if slot["decoder"].complete]
        if buffered or not active:
            return buffered
        channels, _, _ = select.select(
            [slot["channel"] for slot in active], [], [], max(0, timeout)
        )
        return [slot for slot in active if slot["channel"] in channels]

    def _send(self, slot, message):
        slot["channel"].sendall(frame(message))

    def _read(self, slot):
        from .worker import read_message

        return read_message(slot["channel"], slot["decoder"])

    def _destroy(self):
        """Broadcast each stop phase; one shared deadline for all its joins."""
        for slot in self.slots:
            try:
                self._send(slot, dict(kind="close"))
            except (OSError, ValueError):
                pass
        for action, allowance in (
            (None, self.settings.drain_seconds),
            ("terminate", self.settings.kill_seconds),
            ("kill", self.settings.kill_seconds),
        ):
            live = [
                s["process"] for s in self.slots if s["process"].is_alive()
            ]
            if action is not None:
                for process in live:
                    try:
                        getattr(process, action)()
                    except ProcessLookupError:
                        pass
            end = time.monotonic() + allowance
            for process in live:
                process.join(timeout=max(0, end - time.monotonic()))
        survivors = []
        for slot in self.slots:
            slot["channel"].close()
            if slot["process"].is_alive():
                survivors.append(slot)
            else:
                slot["process"].close()
        self.slots = survivors
        if not survivors:
            self.startup_in_progress = False
        if survivors:
            raise RuntimeError(
                "owned workers did not exit: "
                + repr([slot["process"].pid for slot in survivors])
            )

    @accounted_entry
    def close(self, *, budget=None):
        if not self.lock.acquire(blocking=False):
            raise ValueError("cannot close an active ECM call")
        try:
            self._budget = budget
            for slot in self.slots:
                budget.track_process(slot["process"])
            started = time.monotonic()
            self._destroy()
            if self.serial_engine is not None:
                self.serial_engine.close()
            self.serial_engine = None
            self.manifests.clear()
            self.manifest_bytes = 0
            self.closed = True
            self.last_stats["close_seconds"] = time.monotonic() - started
            self.lifetime_stats["close_seconds"] = (
                self.lifetime_stats.get("close_seconds", 0)
                + self.last_stats["close_seconds"]
            )
        finally:
            self.lock.release()

    def _spawn(self, ledger, deadline, poll):
        from .worker import worker_main

        if threading.current_thread() is not threading.main_thread():
            raise ValueError("spawn preparation requires the main thread")
        import multiprocessing

        context = multiprocessing.get_context("spawn")
        if context.get_start_method() != "spawn":
            raise RuntimeError("explicit spawn context unavailable")
        self.helper_funded = True
        self.startup_in_progress = True
        startup_tokens, pending = [], []
        try:
            # Start all funded children before collecting ready reports.
            # Waiting after each start serializes interpreter/import costs.
            for index in range(self.settings.workers):
                poll()
                startup = ledger.reserve("startup")
                startup_tokens.append(startup)
                parent, child = socket.socketpair()
                names = (
                    "OMP_NUM_THREADS",
                    "OPENBLAS_NUM_THREADS",
                    "MKL_NUM_THREADS",
                    "VECLIB_MAXIMUM_THREADS",
                    "NUMEXPR_NUM_THREADS",
                )
                previous = {name: os.environ.get(name) for name in names}
                end = time.monotonic() + self.settings.startup_seconds
                if deadline is not None:
                    end = min(end, deadline)
                try:
                    parent.settimeout(0.01)
                    process = context.Process(
                        target=worker_main,
                        args=(
                            child,
                            asdict(self.settings),
                        ),
                        daemon=True,
                    )
                    for name in names:
                        os.environ[name] = "1"
                    process.start()
                except BaseException:
                    parent.close()
                    child.close()
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
                    decoder=Decoder(),
                    schedule=None,
                    active=[],
                    metadata_cache=[],
                    packet=None,
                    counters={},
                    rss=0,
                    reported_cpu=0.0,
                    index=index,
                )
                self.slots.append(slot)
                # A failed CPU preflight still leaves an owned worker to drain.
                self._budget.track_process(process)
                pending.append((slot, startup, end))

            while pending:
                progress = False
                for entry in tuple(pending):
                    slot, startup, end = entry
                    poll()
                    # Each child keeps its launch deadline; scanning another
                    # child or receiving a fragment cannot extend that grant.
                    if time.monotonic() >= end:
                        raise LeaseRefusedError("startup_timeout")
                    message = self._read(slot)
                    if message is not None:
                        if message.get("kind") != "ready":
                            raise ValueError("missing worker ready report")
                        integer(message["rss"], "startup RSS")
                        used = integer(message["work"], "startup work")
                        consumed = finite(message["cpu"], "startup CPU")
                        ledger.settle(startup, "startup", used, consumed)
                        slot["rss"] = message["rss"]
                        slot["reported_cpu"] = message["cpu"]
                        self._budget.acknowledge(
                            slot["process"].pid, message["cpu"]
                        )
                        pending.remove(entry)
                        progress = True
                    elif not slot["process"].is_alive():
                        raise EOFError("worker exited during startup")
                if pending and not progress:
                    time.sleep(0.001)
            self.startup_in_progress = False
        except BaseException:
            # Every launched child without a checked report keeps its full
            # lease, even if a different child's startup caused the failure.
            for startup in startup_tokens:
                if startup in ledger.active:
                    ledger.lose(startup)
            raise

    @accounted_entry
    def start(self, *, budget):
        """Measure cold pool readiness separately; no curve or plan is run."""
        if sys.implementation.name != "pypy" or sys.version_info[:2] != (
            3,
            11,
        ):
            raise RuntimeError("PyPy implementing Python 3.11 required")
        self._check_configuration()
        if self.closed or not self.lock.acquire(blocking=False):
            raise ValueError("ECM context closed or active")
        base = _base_budget(budget)
        self._budget = base
        for slot in self.slots:
            base.track_process(slot["process"])
        started, cpu_start = time.monotonic(), time.process_time()
        ledger = None
        reported_before = sum(s["reported_cpu"] for s in self.slots)
        helper_before = self.helper_funded
        try:
            ledger = Ledger()
            deadline = _deadline(budget)

            def poll():
                budget.consume(0)
                if deadline is not None and time.monotonic() >= deadline:
                    budget.consume(0)
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
            self.lock.release()
        self.last_stats = dict(
            cold_pool_start_seconds=time.monotonic() - started,
            parent_cpu=time.process_time() - cpu_start,
            reported_child_cpu=sum(s["reported_cpu"] for s in self.slots)
            - reported_before,
            charged_child_cpu=ledger.cpu,
            owned_bytes=self.settings.owned_bytes,
            helper_cpu_reserved=0 if helper_before else 0.0,
            helper_cpu_observed=None,
        )
        self._record(self.last_stats)
        return dict(self.last_stats)

    @accounted_entry
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
        self._check_configuration()
        if self.closed or not self.lock.acquire(blocking=False):
            raise ValueError("ECM context closed or already active")
        try:
            return self._run(
                n, config, budget, lineage, checkpoint, fixed_work
            )
        finally:
            self.lock.release()

    def _manifest(self, n, config, lineage):
        """Reuse identities only, with independent per-call search state."""
        options = arithmetic_options(config, self.settings)
        cap = self.settings.manifest_cache_bytes
        if not cap:
            return Manifest(n, self.seed, lineage, options), False
        key = (n, lineage, canonical(options))
        if key in self.manifests:
            self.manifests.move_to_end(key)
            return self.manifests[key][0], True

        # Reserve every possible assignment before retaining a lazy manifest.
        # This includes its duplicate canonical key/options and integer IDs.
        reserve = (
            65536
            + 2048 * sum(tier[2] for tier in options["ecm_tiers"])
            + 4 * len(key[2])
            + 4 * len(lineage)
            + n.bit_length()
        )
        retain = reserve <= cap
        if retain:
            while self.manifests and (
                self.manifest_bytes + reserve > cap
                or len(self.manifests) >= 16
            ):
                _, (_, owned) = self.manifests.popitem(last=False)
                self.manifest_bytes -= owned
        result = Manifest(n, self.seed, lineage, options)
        if retain:
            self.manifests[key] = (result, reserve)
            self.manifest_bytes += reserve
        return result, False

    def _run(self, n, config, budget, lineage, checkpoint, fixed_work):
        from ..execution.budget import BudgetExhaustedError
        from .worker import CurveEngine

        call_started = time.monotonic()
        base = _base_budget(budget)
        self._budget = base
        for slot in self.slots:
            base.track_process(slot["process"])
        if n.bit_length() > config.max_input_bits:
            raise ValueError("cofactor exceeds configured input bound")
        manifest, manifest_hit = self._manifest(n, config, lineage)
        search = Search(manifest, self.settings, fixed_work=fixed_work)
        if checkpoint is not None:
            require_resource_version(checkpoint)
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
            base.restore_charges(
                work=saved["work"], cpu=saved["cpu"], wall=saved["wall"]
            )
            search = Search.restore(
                checkpoint["search"],
                manifest,
                self.settings,
                fixed_work=fixed_work,
            )
            search.check_leases(checkpoint["leases"], saved)
        search.runs += 1
        ledger = Ledger()
        start_cpu, started = time.process_time(), time.monotonic()
        deadline = _deadline(budget)
        reported_before = sum(s["reported_cpu"] for s in self.slots)
        helper_before = self.helper_funded
        ready = self.ready(manifest)
        decision_started = time.monotonic()
        remaining = tuple(
            sum(
                position not in search.completed
                and position not in search.boundaries
                and job[0] == tier
                for position, job in enumerate(manifest.jobs)
            )
            for tier in range(len(manifest.tiers))
        )
        continuations = tuple(
            sum(
                manifest.jobs[position][0] == tier
                for position in search.boundaries
            )
            for tier in range(len(manifest.tiers))
        )
        decision = self.rule.decision(
            manifest,
            ready=ready,
            pool_started=self.pool_started(),
            remaining=remaining,
            remaining_stage_two=continuations,
            execution=self.cost_identity,
        )
        engaged = self.force_workers or decision["engage"]
        funded = True
        stats = dict(
            decision_seconds=time.monotonic() - decision_started,
            planning_seconds=decision_started - call_started,
            manifest_cache_hits=int(manifest_hit),
            manifest_cache_misses=int(
                bool(self.settings.manifest_cache_bytes) and not manifest_hit
            ),
            manifest_cache_bytes=self.manifest_bytes,
            worker_dispatches=0,
            packets_dispatched=0,
            funded=funded,
            version=VERSION,
            manifest=manifest.identity,
            schedule=manifest.schedule,
            metadata=manifest.metadata,
            cost_identity=self.cost_identity,
            engaged=engaged,
            pool_ready=ready,
            pool_started=self.pool_started(),
            metadata_ready_workers=sum(
                manifest.metadata in slot["metadata_cache"]
                and slot["process"].is_alive()
                for slot in self.slots
            ),
            engagement=decision,
            owned_bytes=self.settings.owned_bytes,
            startup_seconds=0.0,
            drain_seconds=0.0,
            worker_peak_rss_sum=0,
            completed_before=len(search.completed),
            residue_queue_peak=len(search.boundaries),
        )
        counter_names = (
            "preparations",
            "cache_hits",
            "cache_evictions",
            "decode_hits",
            "decode_misses",
            "decode_refusals",
            "stage_two_hits",
            "stage_two_misses",
            "stage_two_refusals",
            "stage_two_build_bytes",
        )
        serial_before = {
            key: getattr(self.serial_engine, key, 0) for key in counter_names
        }
        reason = "curves_exhausted"
        damaged = False

        def poll():
            budget.consume(0)
            if deadline is not None and time.monotonic() >= deadline:
                budget.consume(0)
            ledger.charge(parent_cpu=time.process_time() - start_cpu)

        def receive(slot, message, publish):
            if message.get("kind") != "result" or not slot["active"]:
                raise ValueError("unsolicited worker result")
            if message.get("packet") != slot["packet"]:
                raise ValueError("stale packet result")
            position, token = slot["active"][0]
            if message.get("stage") != search.stage(position):
                raise ValueError("missing or stale wire stage")
            identity = ledger.active[token]
            used = integer(message["work"], "worker work")
            consumed = finite(message["cpu"], "worker CPU")
            integer(message["rss"], "worker RSS")
            keys = message.get("metadata_cache")
            if (
                not isinstance(keys, list)
                or len(keys) > self.settings.metadata_slots
                or len(set(keys)) != len(keys)
            ):
                raise ValueError("invalid metadata cache report")
            if any(not isinstance(key, str) or len(key) != 64 for key in keys):
                raise ValueError("invalid metadata cache identity")
            if message["reason"] in ("complete", "boundary") and (
                message.get("sigma") != manifest.assignment(position)["sigma"]
                or message.get("prepared_schedule") != manifest.metadata
                or manifest.metadata not in keys
            ):
                raise ValueError(
                    "worker mathematical/metadata identity changed"
                )
            counter_updates = {}
            for key in (
                "preparations",
                "cache_hits",
                "cache_evictions",
                "decode_hits",
                "decode_misses",
                "decode_refusals",
                "stage_two_hits",
                "stage_two_misses",
                "stage_two_refusals",
                "stage_two_build_bytes",
            ):
                value = integer(message.get(key, 0), key)
                previous = slot["counters"].get(key, 0)
                if value < previous:
                    raise ValueError("worker counter decreased")
                counter_updates[key] = (value, value - previous)
            search.receive(position, token, message, publish=publish)
            ledger.settle(token, identity, used, consumed)
            ledger.receipts[token - 1].update(
                reason=message["reason"],
                stage=message["stage"],
                failure=message.get("failure"),
            )
            for key, (value, delta) in counter_updates.items():
                stats[key] = stats.get(key, 0) + delta
                slot["counters"][key] = value
            if message["reason"] == "boundary":
                ledger.receipts[token - 1]["boundary_sha256"] = message[
                    "residue"
                ]["sha256"]
            slot["active"].pop(0)
            slot["schedule"] = message.get("prepared_schedule")
            slot["metadata_cache"] = keys
            slot["rss"] = max(slot["rss"], message["rss"])
            slot["reported_cpu"] += consumed
            base.acknowledge(slot["process"].pid, slot["reported_cpu"])
            stats["residue_queue_peak"] = max(
                stats["residue_queue_peak"], len(search.boundaries)
            )
            if publish and message["reason"] not in ("complete", "boundary"):
                raise LeaseRefusedError(message["reason"])

        try:
            poll()
            ledger.charge(
                len(canonical(manifest.options)) + 64 * len(manifest.jobs)
            )
            if engaged and self.slots and not self.pool_started():
                self._destroy()
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
                        if slot["active"]:
                            continue
                        packet, issued = self._packet(
                            manifest,
                            search,
                            ledger,
                            deadline,
                            slot,
                            peers=sum(
                                not other["active"] for other in self.slots
                            ),
                        )
                        if not issued:
                            continue
                        slot["active"], slot["packet"] = (
                            issued,
                            packet["packet"],
                        )
                        stats["packets_dispatched"] += 1
                        stats["worker_dispatches"] += len(issued)
                        self._send(slot, packet)
                    if not search.active:
                        raise LeaseRefusedError("attempts_exhausted")
                    wait = self.settings.parent_poll_seconds
                    if deadline is not None:
                        wait = min(wait, max(0, deadline - time.monotonic()))
                    for slot in self._readable(wait):
                        message = self._read(slot)
                        if message is not None:
                            poll()
                            receive(slot, message, publish=True)
                            if search.retired:
                                break
                    for slot in self.slots:
                        if not slot["process"].is_alive():
                            raise EOFError(
                                "worker exited without acknowledgement"
                            )
                else:
                    if self.serial_engine is None:
                        self.serial_engine = CurveEngine(
                            self.settings, metadata_cache=self.metadata_cache
                        )
                    slot = dict(
                        index=None,
                        metadata_cache=list(self.serial_engine.contexts),
                    )
                    packet, issued = self._packet(
                        manifest, search, ledger, deadline, slot, serial=True
                    )
                    if not issued:
                        raise LeaseRefusedError("attempts_exhausted")
                    cancelled_after_factor = False
                    for index, (task, (position, token)) in enumerate(
                        zip(packet["jobs"], issued)
                    ):
                        task = dict(
                            task,
                            schedule=manifest.metadata,
                            options=packet.get("options")
                            if index == 0
                            else None,
                            packet=packet["packet"],
                        )
                        if config.ecm_kernel == "c-gmp-full":
                            from .native_resources import native_bulk_allowed

                            task["native_bulk"] = native_bulk_allowed(budget)
                        result = self.serial_engine.execute(
                            task,
                            cancelled=lambda: (
                                cancelled_after_factor
                                or bool(base.cancelled and base.cancelled())
                            ),
                        )
                        poll()
                        search.receive(position, token, result)
                        ledger.settle(token, result["id"], result["work"], 0)
                        ledger.receipts[token - 1].update(
                            reason=result["reason"], stage=result["stage"]
                        )
                        if result["reason"] == "boundary":
                            ledger.receipts[token - 1]["boundary_sha256"] = (
                                result["residue"]["sha256"]
                            )
                        if result["factor"] is not None and not fixed_work:
                            cancelled_after_factor = True
                        if (
                            result["reason"] not in ("complete", "boundary")
                            and not cancelled_after_factor
                        ):
                            raise LeaseRefusedError(result["reason"])
        except (BudgetExhaustedError, LeaseRefusedError) as error:
            reason = str(error)
            damaged = reason in ("error", "allocation_failure")
        except (MemoryError, UnsupportedEncodingError):
            damaged = True
            raise
        except (
            OSError,
            EOFError,
            ValueError,
            KeyError,
            TypeError,
        ) as error:
            reason, damaged = "worker_failure", True
            stats["failure"] = type(error).__name__ + ": " + str(error)[:256]
        except BaseException:
            damaged = True
            raise
        finally:
            before = time.monotonic()
            for slot in self.slots:
                if slot["active"]:
                    try:
                        self._send(
                            slot, dict(kind="cancel", packet=slot["packet"])
                        )
                    except OSError:
                        damaged = True
            drain_end = time.monotonic() + self.settings.drain_seconds
            while (
                any(s["active"] for s in self.slots)
                and time.monotonic() < drain_end
            ):
                for slot in self._readable(drain_end - time.monotonic()):
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
                        damaged = True
                        break
                if damaged:
                    break
            stats["worker_peak_rss_sum"] = sum(s["rss"] for s in self.slots)
            stats["reported_child_cpu"] = (
                sum(s["reported_cpu"] for s in self.slots) - reported_before
            )
            damaged = (
                damaged
                or self.startup_in_progress
                or any(s["active"] for s in self.slots)
            )
            if damaged:
                self.serial_engine = None
                try:
                    self._destroy()
                except RuntimeError as error:
                    stats["cleanup_failure"] = str(error)
            if reason == "worker_failure":
                for position in search.active:
                    search.record_failure(position)
            search.active.clear()
            for token in list(ledger.active):
                ledger.lose(token)
            base.used += ledger.work
            if not engaged and self.serial_engine is not None:
                for key in counter_names:
                    stats[key] = (
                        getattr(self.serial_engine, key) - serial_before[key]
                    )
            stats["helper_cpu_reserved"] = (
                0.0 if self.helper_funded and not helper_before else 0
            )
            stats["helper_cpu_observed"] = None
            stats.update(
                drain_seconds=time.monotonic() - before,
                charged_work=ledger.work,
                charged_child_cpu=ledger.cpu,
                parent_cpu=time.process_time() - start_cpu,
                lost_reports=ledger.lost,
                elapsed=time.monotonic() - started,
            )
        try:
            base.publication()
        except BudgetExhaustedError as error:
            reason = str(error)
        # A proper factor may coexist with an interrupted fixed-work search.
        # Keep the public factor outcome and expose the actual traversal stop.
        stats["termination_reason"] = reason
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
            resource_version=RESOURCE_VERSION,
            executor=self.identity,
            leases=([] if checkpoint is None else checkpoint["leases"])
            + ledger.receipts,
            search=search.snapshot(),
            budget=dict(
                work=base.used, cpu=base.cpu_used, wall=base.wall_used
            ),
        )
        output_work = len(canonical(snapshot)) + 512
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
            resource_stop=base.stop,
        )

    def _task(self, manifest, search, position, token, deadline):
        task = dict(
            assignment=manifest.assignment(position),
            token=token,
            attempt=search.attempt(position) + 1,
            stage=search.stage(position),
            deadline=deadline,
        )
        if position in search.boundaries:
            task["residue"] = search.boundaries[position]
        return task

    def _packet(
        self,
        manifest,
        search,
        ledger,
        deadline,
        slot,
        *,
        serial=False,
        peers=1,
    ):
        """Own every queued job before sending; shrink to actual IPC/funds."""
        integer(self.packet_sequence + 1, "packet sequence", 1)
        packet = dict(
            kind="packet",
            packet=self.packet_sequence + 1,
            schedule=manifest.metadata,
            fixed_work=search.fixed_work,
            jobs=[],
        )
        if manifest.metadata not in slot["metadata_cache"]:
            packet["options"] = manifest.options
        issued = []
        integer(peers, "idle peers", 1, self.settings.workers)
        # Share grants and ready work among idle slots before filling packets.
        # Batching must not serialize an otherwise parallel short manifest.
        count = search.dispatchable_count()
        limit = min(
            self.settings.packet_curves, max(1, (count + peers - 1) // peers)
        )
        if (
            not serial
            and self.settings.packet_policy == "tail_one"
            and count <= 2 * self.settings.workers
        ):
            limit = 1
        for _ in range(limit):
            position = search.next(slot["index"])
            if position is None:
                break
            token = ledger.reserve(manifest.assignment(position)["id"])
            task = self._task(manifest, search, position, token, deadline)
            search.dispatch(position, token)
            ledger.receipts[token - 1].update(
                stage=task["stage"], attempt=task["attempt"]
            )
            if task["stage"] == "stage_two":
                ledger.receipts[token - 1]["boundary_sha256"] = task[
                    "residue"
                ]["sha256"]
            packet["jobs"].append(task)
            issued.append((position, token))
        if issued:
            self.packet_sequence += 1
        return packet, issued


_LAZY_IMPORTS = frozenset(("multiprocessing",))


def __getattr__(name):
    if name in _LAZY_IMPORTS:
        module = __import__(name)
        globals()[name] = module
        return module
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return sorted(set(globals()) | _LAZY_IMPORTS)


__all__ = sorted(
    {name for name in globals() if not name.startswith("_")} | _LAZY_IMPORTS
)
