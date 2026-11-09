"""Opt-in current SIQS collection on the retained A7 executor."""

import hashlib
import json
import math
import time
from collections import defaultdict
from concurrent.futures import wait
from concurrent.futures.process import BrokenProcessPool
from dataclasses import asdict, dataclass, replace

from ..common import utils
from ..execution.budget import BudgetExhaustedError
from ..execution.resources import (
    AllocationFailureError,
    accounted_job,
    removed_arguments,
)
from .parallel import CollectionPool, _ParentBudget
from .pipeline import QSJob, QSResult
from .polynomial import qs_polynomial
from .siqs import SIQSConfig, SIQSJob
from .storage import require_published
from .worker_packet import (
    VERSION,
    checked_int,
    decode_atom,
    digest,
    execute_packet,
    polynomial_at,
)


@removed_arguments(
    (
        "assignment_work",
        "assignment_seconds",
        "assignment_cpu",
        "parent_work_reserve",
        "parent_cpu_reserve",
        "startup_work",
        "shutdown_work",
        "shutdown_cpu",
        "helper_cpu",
        "worker_memory_bytes",
        "runtime_memory_bytes",
        "memory_bytes",
        "coexistence_bytes",
        "checkpoint_bytes",
        "context_limit",
        "receipt_limit",
        "max_attempts",
        "max_runs",
    )
)
@dataclass(frozen=True)
class QSWorkerConfig:
    """Retained scheduling/cache capacities under one observed root owner."""

    mode: str = "process"
    workers: int = 1
    polynomials_per_task: int = 4
    packet_bytes: int = 4 * 2**20
    max_batch_atoms: int = 8192
    split_calls_per_task: int = 4096
    drain_seconds: float = 2.0
    observe: bool = False
    lease_budget: bool = False
    reuse_state: bool = False
    state_cache_bytes: int = 64 * 2**20
    packed_transport: bool = False
    verification: str = "parent"
    verification_batch_atoms: int = 64
    lookahead: int = 0
    whole_family: bool = False
    arrival_order: bool = False
    private_stop: bool = False
    cpu_snapshot: bool = False
    stop_backend: str = "event"
    coefficient_cache: bool = False
    family_affinity: bool = False
    grouped_lease_charging: bool = False

    def __post_init__(self):
        if type(self.observe) is not bool:
            raise ValueError("observe must be Boolean")
        if type(self.lease_budget) is not bool:
            raise ValueError("lease_budget must be Boolean")
        if type(self.reuse_state) is not bool:
            raise ValueError("reuse_state must be Boolean")
        if type(self.packed_transport) is not bool:
            raise ValueError("packed_transport must be Boolean")
        if self.verification not in ("parent", "worker", "batch"):
            raise ValueError("unknown relation verification policy")
        if type(self.whole_family) is not bool:
            raise ValueError("whole_family must be Boolean")
        if type(self.arrival_order) is not bool:
            raise ValueError("arrival_order must be Boolean")
        if type(self.private_stop) is not bool:
            raise ValueError("private_stop must be Boolean")
        if self.private_stop and (
            self.mode != "process" or not self.lease_budget
        ):
            raise ValueError(
                "private stop flags require process packet leases"
            )
        if type(self.cpu_snapshot) is not bool:
            raise ValueError("cpu_snapshot must be Boolean")
        if type(self.coefficient_cache) is not bool:
            raise ValueError("coefficient_cache must be Boolean")
        if type(self.family_affinity) is not bool or (
            self.family_affinity
            and (self.mode != "process" or not self.reuse_state)
        ):
            raise ValueError("family affinity requires reused process state")
        if type(self.grouped_lease_charging) is not bool or (
            self.grouped_lease_charging and not self.lease_budget
        ):
            raise ValueError("grouped charging requires packet leases")
        if self.stop_backend not in ("event", "semaphore") or (
            self.stop_backend == "semaphore" and not self.private_stop
        ):
            raise ValueError("invalid private stop backend")
        if self.cpu_snapshot and (
            self.mode != "process" or not self.lease_budget
        ):
            raise ValueError("CPU snapshots require process packet leases")
        if self.mode not in ("serial", "process"):
            raise ValueError(
                "current workers support serial or spawn processes"
            )
        if type(self.workers) is not int or not 1 <= self.workers <= 20:
            raise ValueError("current preparation supports 1 to 20 workers")
        checked_int(self.lookahead, "lookahead", 0, 4 * self.workers)
        if self.mode == "serial" and self.workers != 1:
            raise ValueError("serial execution requires one worker")
        for name, low, high in (
            ("polynomials_per_task", 1, 64),
            ("packet_bytes", 4096, 16 * 2**20),
            ("max_batch_atoms", 1, 65536),
            ("split_calls_per_task", 1, 1048576),
            ("state_cache_bytes", 4096, 256 * 2**20),
            ("verification_batch_atoms", 1, 256),
        ):
            checked_int(getattr(self, name), name, low, high)
        if (
            type(self.drain_seconds) not in (int, float)
            or not math.isfinite(self.drain_seconds)
            or not 0 < self.drain_seconds <= 3600
        ):
            raise ValueError("drain_seconds must be finite and positive")

    def record(self):
        """Serialize retained capacities and enabled execution options."""
        values = asdict(self)
        if not self.lease_budget:
            values.pop("lease_budget")
        if not self.reuse_state:
            values.pop("reuse_state")
            values.pop("state_cache_bytes")
        if not self.packed_transport:
            values.pop("packed_transport")
        if self.verification == "parent":
            values.pop("verification")
        if self.verification != "batch":
            values.pop("verification_batch_atoms")
        if not self.lookahead:
            values.pop("lookahead")
        if not self.whole_family:
            values.pop("whole_family")
        if not self.arrival_order:
            values.pop("arrival_order")
        if not self.private_stop:
            values.pop("private_stop")
        if not self.cpu_snapshot:
            values.pop("cpu_snapshot")
        if self.stop_backend == "event":
            values.pop("stop_backend")
        if not self.coefficient_cache:
            values.pop("coefficient_cache")
        if not self.family_affinity:
            values.pop("family_affinity")
        if not self.grouped_lease_charging:
            values.pop("grouped_lease_charging")
        return values

    @property
    def credit_limit(self):
        """Bound both executing leases and paid queued/pending responses."""
        return self.workers + self.lookahead


class CurrentParallelSIQSJob:
    """Ordered bounded packets and one authoritative current relation store.

    Mathematical identity excludes executor count. Execution identity includes
    count/mode and is immutable across resume. Interrupted prefixes retain all
    reported atoms, a first-unpublished cursor, and their consumed charges.
    """

    def __init__(self, n, *, seed=7, config=None, workers=None, budget=None):
        self.workers = workers if workers is not None else QSWorkerConfig()
        if not isinstance(self.workers, QSWorkerConfig):
            raise TypeError("workers must be a QSWorkerConfig")
        self.source = SIQSJob(n, seed=seed, config=config, budget=budget)
        self.config, self.budget = self.source.config, self.source.budget
        config = self.config
        if not isinstance(config, SIQSConfig) or (
            config.mode != "siqs"
            or config.growth_steps
            or not config.shared_relations
            or not config.diverse
            or config.backend != "python-int"
        ):
            raise ValueError(
                "current workers require fixed-width native shared SIQS"
            )
        if getattr(config.collector, "recovery_workers", 0):
            raise ValueError("worker recovery must be resident and non-nested")
        self.n, self.seed = self.source.n, seed
        self.identity = digest(
            [
                VERSION,
                "qs_collect",
                self.n,
                seed,
                asdict(config),
                (
                    config.gray_limit
                    if self.workers.whole_family
                    else self.workers.polynomials_per_task
                ),
            ]
        )
        self.execution_identity = digest(
            [self.identity, self.workers.record()]
        )
        self._original = (self.config, self.workers, self.budget)
        self.engine = None
        self.next_assignment = self.cursor = self.runs = 0
        self.pending, self.progress, self.attempts, self.contexts = (
            {},
            {},
            {},
            {},
        )
        self.receipts = []
        self.ahead = {}
        self.split_used = 0
        self.split_unknown = False
        self.divisor = self.deferred_divisor = None
        self.fixed_work = None
        self._running = False
        self.stats = dict(
            scanned=0,
            worker_work=0,
            cancelled_work=0,
            lost_work=0,
            lost_packets=0,
            admitted_atoms=0,
            useful_rows=0,
            received_atoms=0,
            received_packets=0,
            assignments=0,
            parent_cpu_seconds=0.0,
            worker_cpu_seconds=0.0,
            worker_seconds=0.0,
            transport_seconds=0.0,
            credit_wait_seconds=0.0,
            lifecycle_work=0,
            drain_work_upper=0,
            validation_seconds=0.0,
            admission_seconds=0.0,
            solve_seconds=0.0,
            decode_seconds=0.0,
            wait_seconds=0.0,
            ordered_stall_seconds=0.0,
            drain_seconds=0.0,
            pending_bytes_high=0,
            pending_atoms_high=0,
            pending_packets_high=0,
            backlog_byte_seconds=0.0,
            backlog_atom_seconds=0.0,
            observed_seconds=0.0,
        )
        self._backlog_time = None
        if self.workers.lease_budget:
            self.stats["lease_reconciliations"] = 0
        if self.workers.grouped_lease_charging:
            self.stats["grouped_marking_pairs"] = 0
        if self.workers.reuse_state:
            self.stats.update(
                state_install_work=0,
                state_install_seconds=0.0,
                state_context_hits=0,
                state_family_reuses=0,
                state_collector_reuses=0,
            )

    @property
    def chunks_per_family(self):
        width = self.packet_width
        return (self.config.gray_limit + width - 1) // width

    @property
    def packet_width(self):
        return (
            self.config.gray_limit
            if self.workers.whole_family
            else self.workers.polynomials_per_task
        )

    @property
    def assignment_count(self):
        return (
            0
            if self.source.assignments is None
            else len(self.source.assignments) * self.chunks_per_family
        )

    def _setup(self):
        self.source._setup()
        self.divisor = self.source.divisor
        if self.divisor is not None or self.engine is not None:
            return
        options = replace(self.config.collector)
        self.engine = QSJob(
            qs_polynomial(self.source.base),
            self.source.base,
            0,
            0,
            budget=self.budget,
            config=options,
            weight_two=self.config.weight_two,
            row_excess=self.config.row_excess,
            batch_width=self.config.batch_width,
            filter_row_growth=self.config.filter_row_growth,
            tested_dependencies=self.config.tested_dependencies,
            collection_continues=self.config.paced_collection,
        )
        if self.workers.reuse_state:
            from .worker_state import context_key

            base = self.source.base
            self._state_key = context_key(
                dict(
                    config=self.config,
                    settings=self.workers,
                    base=(base.n, base.multiplier, base.bound, base.entries),
                    job=self.identity,
                    execution=self.execution_identity,
                )
            )

    def _spec(self, index):
        family, chunk = divmod(index, self.chunks_per_family)
        start = chunk * self.packet_width
        end = min(self.config.gray_limit, start + self.packet_width)
        primes = self.source.assignments[family]
        return primes, start, end

    def _task(self, index, lease, cpu, seconds):
        primes, start, end = self._spec(index)
        attempt = self.attempts.get(index, 0) + 1
        cursor = self.progress.get(index, [start, -self.config.half_width])
        limit = getattr(self.config.collector, "split_call_limit", None)
        reserved = sum(
            r.get("split_grant", 0) or 0
            for r in self.receipts
            if r["status"] == "issued"
        )
        grant = (
            None
            if limit is None
            else min(
                self.workers.split_calls_per_task,
                max(0, limit - self.split_used - reserved),
            )
        )
        token = digest(
            [self.execution_identity, index, attempt, cursor, grant]
        )
        base = self.source.base
        return dict(
            config=self.config,
            settings=self.workers,
            base=(base.n, base.multiplier, base.bound, base.entries),
            job=self.identity,
            id=index,
            attempt=attempt,
            token=token,
            primes=primes,
            start=start,
            end=end,
            cursor=list(cursor),
            work=lease,
            cpu=cpu,
            seconds=seconds,
            fixed_work=self.fixed_work,
            split_grant=grant,
            split_soft=grant is not None and grant > 0,
        )

    def _task_for_restore(self, index, report, receipt):
        primes, start, end = self._spec(index)
        return dict(
            id=index,
            job=self.identity,
            token=report["token"],
            attempt=report["attempt"],
            cursor=receipt["start"],
            end=end,
            start=start,
            primes=primes,
            work=receipt["work"],
            split_grant=receipt["split_grant"],
        )

    def _backlog(self):
        now = time.monotonic()
        atoms = sum(
            len(r["rows"]) - r["admitted"] for r in self.pending.values()
        )
        size = sum(r["owned_bytes"] for r in self.pending.values())
        if self.workers.arrival_order:
            size += sum(r["owned_bytes"] for r in self.ahead.values())
        if self._backlog_time is not None:
            elapsed = now - self._backlog_time
            self.stats["backlog_byte_seconds"] += elapsed * self._backlog_size
            self.stats["backlog_atom_seconds"] += elapsed * self._backlog_atoms
            self.stats["observed_seconds"] += elapsed
        self._backlog_time, self._backlog_size, self._backlog_atoms = (
            now,
            size,
            atoms,
        )
        for name, value in (
            ("bytes", size),
            ("atoms", atoms),
            ("packets", len(self.pending) + len(self.ahead)),
        ):
            key = "pending_" + name + "_high"
            self.stats[key] = max(self.stats[key], value)

    def _validate_report(self, report, task):
        if self.workers.packed_transport and (
            type(report) is not dict or report.get("encoding") != "packed-v1"
        ):
            raise ValueError("packed transport encoding mismatch")
        if self.workers.reuse_state:
            if (
                type(report) is not dict
                or report.get("context_key") != self._state_key
            ):
                raise ValueError("worker context report identity mismatch")
            state_stats = report.get("state_stats")
            if type(state_stats) is not dict or set(state_stats) != {
                "context_hits",
                "family_reuses",
                "collector_reuses",
            }:
                raise ValueError("invalid worker state counters")
            for value in state_stats.values():
                checked_int(value, "state counter", 0, 1)
            if max(state_stats.values()) > state_stats["context_hits"]:
                raise ValueError("uninstalled worker state reports reuse")
        if self.workers.lease_budget and (
            type(report) is not dict
            or report.get("lease_policy") != "packet-parent-v1"
            or report.get("publication_stop")
            not in (None, "cancelled", "wall_limit", "cpu_limit", "work_limit")
        ):
            raise ValueError("packet lease guard policy mismatch")
        if type(report) is not dict or any(
            report.get(k) != v
            for k, v in (
                ("version", VERSION),
                ("kind", "qs_collect"),
                ("job", self.identity),
                ("id", task["id"]),
                ("token", task["token"]),
                ("attempt", task["attempt"]),
                ("start", task["cursor"]),
            )
        ):
            raise ValueError("stale or malformed QS assignment report")
        work = checked_int(report["work"], "reported work", 0)
        if self.workers.grouped_lease_charging:
            checked_int(
                report.get("grouped_marking_pairs"), "marking groups", 0, work
            )
        checked_int(
            report["split_calls"],
            "reported splitting",
            0,
            task["split_grant"]
            if task["split_grant"] is not None
            else 2**63 - 1,
        )
        for name in ("cpu_seconds", "seconds", "finished_at"):
            value = report[name]
            if (
                type(value) not in (float, int)
                or not math.isfinite(value)
                or value < 0
            ):
                raise ValueError("invalid worker timing")
        raw = report["packet"]
        checked_int(
            report["workspace_bytes"],
            "worker ownership",
            0,
        )
        if type(raw) is not bytes or len(raw) > self.workers.packet_bytes:
            raise ValueError("oversized QS packet")
        if hashlib.sha256(raw).hexdigest() != report["sha256"]:
            raise ValueError("QS packet digest mismatch")
        cursor = report["cursor"]
        lo, hi = -self.config.half_width, self.config.half_width + 1
        if not isinstance(cursor, list) or len(cursor) != 2:
            raise ValueError("invalid publication cursor")
        checked_int(cursor[0], "Gray cursor", task["cursor"][0], task["end"])
        if (
            type(cursor[1]) is not int
            or not lo <= cursor[1] <= hi
            or (cursor[0] == task["end"] and cursor[1] != lo)
            or tuple(cursor) < tuple(task["cursor"])
        ):
            raise ValueError("publication cursor outside assignment")
        span = 2 * self.config.half_width + 1
        coverage = (
            (cursor[0] - task["cursor"][0]) * span
            + cursor[1]
            - task["cursor"][1]
        )
        if checked_int(report["scanned"], "scanned", 0) != coverage:
            raise ValueError("reported coverage disagrees with cursor")
        reason = report["reason"]
        if reason not in (
            "complete",
            "packet_limit",
            "work_limit",
            "wall_limit",
            "cpu_limit",
            "cancelled",
            "memory_limit",
            "factor_found",
            "split_lease_limit",
        ):
            raise ValueError("invalid worker outcome")
        if (reason == "complete") != (cursor[0] == task["end"]):
            raise ValueError("complete report has incomplete coverage")
        if report["divisor"] is not None and not utils.valid_divisor(
            report["divisor"], self.n
        ):
            raise ValueError("invalid worker divisor")
        started = time.perf_counter()
        if self.workers.packed_transport:
            from .worker_wire import decode, norm_bits

            bits = norm_bits(self.config, self.n)
            rows = decode(
                raw,
                max_rows=self.workers.max_batch_atoms,
                max_exponents=min(
                    len(self.source.base.entries), bits // 2 + 1
                ),
                max_bits=bits,
                legacy_limit=self.workers.packet_bytes,
            )
        else:
            rows = json.loads(raw)
        if (
            not isinstance(rows, list)
            or len(rows) > self.workers.max_batch_atoms
        ):
            raise ValueError("packet atom cap exceeded")
        previous = None
        for row in rows:
            if not isinstance(row, list) or len(row) != 8:
                raise ValueError("invalid packet row")
            gray = checked_int(
                row[0], "atom Gray index", task["cursor"][0], task["end"] - 1
            )
            position = row[1]
            if type(position) is not int or not lo <= position < hi:
                raise ValueError("atom outside window")
            key = (gray, position)
            if (
                key < tuple(task["cursor"])
                or key >= tuple(cursor)
                or (previous is not None and key <= previous)
            ):
                raise ValueError("atom outside published prefix/order")
            if not isinstance(row[3], list) or len(row[3]) > len(
                self.source.base.entries
            ):
                raise ValueError("atom exponent dimension exceeds its base")
            previous = key
        self.stats["decode_seconds"] += time.perf_counter() - started
        # Bound decoded objects before retaining them. Byte JSON itself is
        # bounded; full atom expansion was reserved prospectively in settings.
        decoded_bound = sum(
            8192 + 512 * len(row[3])
            for row in rows
            if isinstance(row[3], list)
        )
        return dict(
            report,
            rows=rows,
            admitted=0,
            owned_bytes=4 * len(raw) + decoded_bound,
        ), work

    def _receive(self, report, task, receipt):
        if receipt["status"] != "issued" or not (
            self.next_assignment
            <= task["id"]
            < self.next_assignment + self.workers.credit_limit
        ):
            raise ValueError("repeated or foreign lease acknowledgement")
        self._backlog()
        checked, work = self._validate_report(report, task)
        if task["id"] in self.pending:
            raise ValueError("repeated assignment publication")
        self.budget.used += work
        self.split_used += report["split_calls"]
        self.engine.collector._split_calls = self.split_used
        receipt.update(
            status="reported",
            charged_work=work,
            reason=report["reason"],
            cpu=report["cpu_seconds"],
            cursor=report["cursor"],
            split_grant=task["split_grant"],
            split_calls=report["split_calls"],
        )
        self.pending[task["id"]] = checked
        if self.workers.arrival_order:
            checked["arrival_seq"] = self.stats["received_packets"]
            receipt["arrival_seq"] = checked["arrival_seq"]
        if self.workers.lease_budget:
            self.stats["lease_reconciliations"] += 1
        if self.workers.grouped_lease_charging:
            self.stats["grouped_marking_pairs"] += report[
                "grouped_marking_pairs"
            ]
        checked["received_at"] = time.monotonic()
        self.stats["received_packets"] += 1
        self.stats["received_atoms"] += len(checked["rows"])
        self.stats["worker_work"] += work
        self.stats["worker_cpu_seconds"] += report["cpu_seconds"]
        self.stats["worker_seconds"] += report["seconds"]
        if self.workers.reuse_state:
            for key, value in report["state_stats"].items():
                self.stats["state_" + key] += value
        self.stats["transport_seconds"] += max(
            0, time.monotonic() - report["finished_at"]
        )
        self.stats["scanned"] += report["scanned"]
        if report["reason"] != "complete":
            self.stats["cancelled_work"] += work
        self._backlog()

    def _ready_assignment(self):
        if not self.workers.arrival_order:
            return (
                self.next_assignment
                if self.next_assignment in self.pending
                else None
            )
        eligible = [
            i
            for i, r in self.pending.items()
            if i == self.next_assignment or r["reason"] == "complete"
        ]
        if not eligible:
            return None
        return min(
            eligible,
            key=lambda i: (
                not bool(self.pending[i]["admitted"]),
                self.pending[i]["arrival_seq"],
            ),
        )

    def _merge(self, index=None):
        index = self.next_assignment if index is None else index
        report = self.pending[index]
        if index != self.next_assignment and (
            not self.workers.arrival_order or report["reason"] != "complete"
        ):
            raise ValueError("unfunded or incomplete out-of-order admission")
        primes, _, _ = self._spec(index)
        polynomials = {}
        coefficient_cache = None
        if self.workers.coefficient_cache:
            from .coefficient_cache import CoefficientCache

            coefficient_cache = CoefficientCache(
                self.source.base, primes, self.config
            )
        collector = self.engine.collector
        stats = defaultdict(int)
        prepared_atoms = {}
        while report["admitted"] < len(report["rows"]):
            row = report["rows"][report["admitted"]]
            gray = row[0]
            if gray not in polynomials:
                polynomials[gray] = polynomial_at(
                    self.source.base,
                    primes,
                    gray,
                    self.config,
                    self.engine.budget,
                    coefficient_cache=coefficient_cache,
                )
            polynomial = polynomials[gray]
            started = time.perf_counter() if self.workers.observe else 0
            checked = collector._checked
            if self.workers.verification == "batch" and checked is not None:
                if not prepared_atoms:
                    stop = min(
                        len(report["rows"]),
                        report["admitted"]
                        + self.workers.verification_batch_atoms,
                    )
                    for row_index in range(report["admitted"], stop):
                        member = report["rows"][row_index]
                        if member[0] not in polynomials:
                            polynomials[member[0]] = polynomial_at(
                                self.source.base,
                                primes,
                                member[0],
                                self.config,
                                self.engine.budget,
                                coefficient_cache=coefficient_cache,
                            )
                        candidate = decode_atom(
                            member,
                            polynomials[member[0]],
                            self.config.collector,
                            self.source.base,
                            self.engine.budget,
                            verify=False,
                        )
                        checked.atom(candidate, self.engine.budget)
                        prepared_atoms[row_index] = candidate
                atom = prepared_atoms.pop(report["admitted"])
            else:
                atom = decode_atom(
                    row,
                    polynomial,
                    self.config.collector,
                    self.source.base,
                    self.engine.budget,
                    verify=not (
                        self.workers.verification == "worker"
                        and checked is not None
                    ),
                )
            if self.workers.observe:
                self.stats["validation_seconds"] += (
                    time.perf_counter() - started
                )
            started = time.perf_counter() if self.workers.observe else 0
            before = len(collector._rows)
            refusal = collector._admit(atom, stats)
            if self.workers.observe:
                self.stats["admission_seconds"] += (
                    time.perf_counter() - started
                )
            if refusal:
                raise MemoryError(refusal)
            self.contexts[polynomial.identity] = [index, gray]
            report["admitted"] += 1
            self.stats["admitted_atoms"] += 1
            self.stats["useful_rows"] += len(collector._rows) - before
            if self.workers.observe:
                self._backlog()
        self._backlog()
        if report["divisor"] is not None:
            self.deferred_divisor = report["divisor"]
            if not self.fixed_work:
                self.divisor = report["divisor"]
        self.progress[index] = report["cursor"]
        complete = report["reason"] == "complete"
        reason = report["reason"]
        self.stats["credit_wait_seconds"] += (
            time.monotonic() - report["received_at"]
        )
        del self.pending[index]
        if complete:
            del self.progress[index]
            self.attempts.pop(index, None)
            if self.workers.arrival_order:
                self.ahead[index] = report
                while self.next_assignment in self.ahead:
                    del self.ahead[self.next_assignment]
                    self.next_assignment += 1
            else:
                self.next_assignment += 1
            self.stats["assignments"] += 1
        self._backlog()
        if not self.fixed_work and self.divisor is None:
            self._solve(False)
        return (
            None
            if complete or reason in ("packet_limit", "split_lease_limit")
            else reason
        )

    def _solve(self, final):
        started = time.perf_counter()
        try:
            self.divisor = self.engine._solve(final)
        finally:
            self.stats["solve_seconds"] += time.perf_counter() - started

    def _lost(self, task, receipt):
        if receipt["status"] != "issued":
            return
        receipt.update(status="lost", charged_work=0, cpu=None)
        self.stats["lost_packets"] += 1
        # Missing private work is unknown; root CPU retains OS observations.
        grant = task["split_grant"]
        receipt["split_grant"] = grant
        receipt["split_calls"] = grant
        if grant is None:
            self.split_unknown = True
        else:
            self.split_used += grant
        if self.engine is not None:
            self.engine.collector._split_calls = self.split_used

    @accounted_job
    def run(self, *, pool=None, max_assignments=None, fixed_work=False):
        """Drain before return; extend the same Budget to continue a pause."""
        if (
            self.config,
            self.workers,
            self.budget,
        ) != self._original or self.budget is not self._original[2]:
            raise ValueError("job configuration/budget identity changed")
        if type(fixed_work) is not bool or (
            self.fixed_work is not None and self.fixed_work != fixed_work
        ):
            raise ValueError("fixed-work mode changed")
        if self.workers.arrival_order and fixed_work:
            raise ValueError(
                "arrival-order admission requires a trajectory comparison"
            )
        if self._running:
            raise ValueError("job already running")
        if self.engine is not None:
            require_published(self.engine.collector)
        if max_assignments is not None:
            checked_int(max_assignments, "max assignments", 1)
        if self.divisor is not None:
            return self._result("factor_found")
        self.fixed_work = fixed_work
        own_pool = pool is None
        pool = (
            pool
            if pool is not None
            else CollectionPool(
                self.workers.mode,
                self.workers.workers,
                family_affinity=self.workers.family_affinity,
            )
        )
        if (
            pool.mode != self.workers.mode
            or pool.workers != self.workers.workers
            or pool.family_affinity != self.workers.family_affinity
        ):
            raise ValueError("pool differs from pinned execution identity")
        from .worker_stop import configure_stop

        configure_stop(
            pool, self.workers.private_stop, self.workers.stop_backend
        )
        retained_settings = getattr(pool, "_qs_state_settings", None)
        if (
            retained_settings is not None
            and retained_settings != self.workers.record()
        ):
            raise ValueError(
                "retained worker settings changed; use a new pool"
            )
        self._running = True
        self.runs += 1
        inflight, begun, broken = {}, False, False
        reason, merged = "families_exhausted", 0
        parent_started = time.process_time()
        try:
            pool.begin(self.budget)
            begun = True
            self._setup()
            state_key = None
            if self.workers.reuse_state and self.engine is not None:
                from .worker_state import ensure_context

                state_key = ensure_context(self, pool)
            self._backlog_time = None
            self._backlog()
            if self.engine is not None:
                if self.workers.cpu_snapshot:
                    from ..execution.work_budget import PollingBudget
                    from .packet_lease import ParentCPUSnapshot

                    parent_budget = PollingBudget(
                        self.budget, poll=ParentCPUSnapshot(pool)._account_cpu
                    )
                else:
                    parent_budget = _ParentBudget(self.budget, pool)
                self.engine.budget = self.engine.collector.budget = (
                    parent_budget
                )
            while (
                not self.divisor
                and self.next_assignment < self.assignment_count
            ):
                pool.poll(self.budget)
                if (
                    self.engine.solver is not None
                    or self.engine.extractor is not None
                ):
                    self._solve(False)
                    if self.divisor:
                        break
                ready = self._ready_assignment()
                if ready is not None:
                    previous = self.stats["assignments"]
                    stopped = self._merge(ready)
                    merged += self.stats["assignments"] - previous
                    if stopped is not None:
                        reason = stopped
                        break
                    if (
                        max_assignments is not None
                        and merged >= max_assignments
                    ):
                        reason = "paused"
                        break
                    continue
                occupied = (
                    set(self.pending)
                    | set(self.ahead)
                    | {task["id"] for task, _ in inflight.values()}
                )
                upper = min(
                    self.assignment_count,
                    self.next_assignment + self.workers.credit_limit,
                )
                if max_assignments is not None:
                    upper = min(
                        upper, self.next_assignment + max_assignments - merged
                    )
                for index in range(self.next_assignment, upper):
                    if index in occupied:
                        continue
                    self.budget.consume(0)
                    seconds = (
                        None
                        if self.budget.seconds is None
                        else max(
                            0, self.budget.seconds - self.budget.wall_used
                        )
                    )
                    task = self._task(index, None, None, seconds)
                    split_limit = getattr(
                        self.config.collector, "split_call_limit", None
                    )
                    if (
                        task["split_grant"] == 0
                        and split_limit is not None
                        and self.split_used < split_limit
                    ):
                        # Other leases own the remaining splitter calls.
                        # Wait for settlement rather than discard DLP intake.
                        break
                    self.attempts[index] = task["attempt"]
                    receipt = dict(
                        id=index,
                        attempt=task["attempt"],
                        token=task["token"],
                        start=task["cursor"],
                        work=None,
                        cpu_grant=None,
                        status="issued",
                        charged_work=0,
                        controller_work=0,
                        split_grant=task["split_grant"],
                    )
                    self.receipts.append(receipt)
                    wire_task = task
                    if state_key is not None:
                        from .worker_state import compact_task

                        wire_task = compact_task(task, state_key)
                    if pool.executor is None:
                        try:
                            self._receive(
                                execute_packet(wire_task, pool.state),
                                task,
                                receipt,
                            )
                        except BaseException:
                            self._lost(task, receipt)
                            raise
                    else:
                        try:
                            if self.workers.family_affinity:
                                slot = (
                                    index // self.chunks_per_family
                                ) % self.workers.workers
                                future = pool.executor.submit_to(
                                    slot, execute_packet, wire_task
                                )
                            else:
                                future = pool.executor.submit(
                                    execute_packet, wire_task
                                )
                            # BrokenProcessPool's manager may clear its table
                            # before abort runs. Keep our own bounded handles
                            # from successful issuance until pool retirement.
                            pool._a7a_processes = tuple(
                                {
                                    process.pid: process
                                    for process in (
                                        *getattr(pool, "_a7a_processes", ()),
                                        *(
                                            pool.executor._processes or {}
                                        ).values(),
                                    )
                                }.values()
                            )
                        except BaseException:
                            self._lost(task, receipt)
                            broken = True
                            raise
                        inflight[future] = task, receipt
                if self._ready_assignment() is not None:
                    continue
                if not inflight:
                    break
                started = time.perf_counter()
                done, _ = wait(inflight, timeout=0.01)
                elapsed = time.perf_counter() - started
                self.stats["wait_seconds"] += elapsed
                if self.pending:
                    self.stats["ordered_stall_seconds"] += elapsed
                for future in done:
                    task, receipt = inflight.pop(future)
                    try:
                        self._receive(future.result(), task, receipt)
                    except AllocationFailureError:
                        self._lost(task, receipt)
                        broken = True
                        raise
                    except Exception:
                        self._lost(task, receipt)
                        broken = True
                        reason = "worker_failure"
                        break
                if broken:
                    break
            if (
                not self.divisor
                and self.next_assignment == self.assignment_count
            ):
                self._solve(True)
                if self.divisor is None:
                    self.divisor = self.deferred_divisor
        except BrokenProcessPool:
            # Submission can fail before a future reaches the tracking table.
            # A dead worker may own the stop lock: retire without broadcasting.
            broken = True
            reason = "worker_failure"
        except BudgetExhaustedError:
            reason = self.budget.reason
        except MemoryError as error:
            # A failed parent publication can leave the packet cursor behind
            # an inserted atom. Never serialize or resume that uncertain store.
            if self.engine is not None:
                self.engine.collector._publication_failed = True
            raise AllocationFailureError(
                "current SIQS worker collection", "qs_worker"
            ) from error
        finally:
            try:
                if begun:
                    started = time.monotonic()
                    # A detected dead worker may own the stop lock. Preserve
                    # missing-report uncertainty and use hard retirement.
                    if not broken:
                        pool.state["stop"].set()
                    deadline = started + self.workers.drain_seconds
                    while (
                        inflight and not broken and time.monotonic() < deadline
                    ):
                        done, _ = wait(
                            inflight,
                            timeout=min(
                                0.01, max(0, deadline - time.monotonic())
                            ),
                        )
                        for future in done:
                            task, receipt = inflight.pop(future)
                            try:
                                self._receive(future.result(), task, receipt)
                                self.stats["drain_work_upper"] += receipt[
                                    "charged_work"
                                ]
                            except Exception:
                                self._lost(task, receipt)
                                broken = True
                    if inflight or broken:
                        for task, receipt in inflight.values():
                            self._lost(task, receipt)
                        pool.abort()
                        reason = "worker_failure"
                    else:
                        try:
                            pool.finish(self.budget)
                        except BudgetExhaustedError:
                            reason = self.budget.reason
                        except Exception:
                            pool.abort()
                            reason = "worker_failure"
                    self.stats["drain_seconds"] += time.monotonic() - started
            finally:
                if begun:
                    pool.lock.release()
                if self.engine is not None:
                    self.engine.budget = self.engine.collector.budget = (
                        self.budget
                    )
                    for continuation in (
                        self.engine.solver,
                        self.engine.extractor,
                    ):
                        if continuation is not None:
                            continuation.budget = self.budget
                if own_pool:
                    pool.close(budget=self.budget)
                self._backlog()
                self._running = False
                self.stats["parent_cpu_seconds"] += (
                    time.process_time() - parent_started
                )
        return self._result("factor_found" if self.divisor else reason)

    def _result(self, reason):
        if self.divisor is not None and not utils.valid_divisor(
            self.divisor, self.n
        ):
            raise ArithmeticError("invalid current-worker split")
        stats = dict(
            self.stats,
            work_used=self.budget.used,
            cpu_used=self.budget.cpu_used,
            job=self.identity,
            execution=self.execution_identity,
            pending_assignments=sorted(self.pending),
            schedule_complete=self.next_assignment == self.assignment_count,
            observed_rss_bytes=self.budget.observed_rss,
            split_calls_charged=self.split_used,
            split_calls_unknown=self.split_unknown,
        )
        if self.engine is not None:
            stats["engine"] = dict(self.engine.stats)
        return QSResult(
            self.divisor,
            self.n // self.divisor if self.divisor else self.n,
            reason,
            self.next_assignment,
            stats,
            resource_stop=self.budget.stop,
        )

    def checkpoint(self):
        from .current_checkpoint import pack

        return pack(self)

    @classmethod
    def from_checkpoint(cls, checkpoint, *, budget, workers=None):
        from .current_checkpoint import restore

        return restore(cls, checkpoint, budget=budget, workers=workers)
