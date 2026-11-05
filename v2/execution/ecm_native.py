"""Opt-in single-process ECM owner with C-only OpenMP stage-one chunks."""

import hashlib
import sys
import threading
import time
from dataclasses import asdict
from math import prod
from pathlib import Path

from ..common import utils
from ..ecm.native_packet import poll_budget, source_identity
from ..ecm.native_team import (
    SCHEDULES,
    TEAM_OVERHEAD,
    NativeTeam,
    TeamRefusedError,
)
from .budget import BudgetExhaustedError
from .ecm_worker import CurveEngine
from .ecm_workers import WorkerRun, _base_budget, arithmetic_options
from .worker_contract import (
    Manifest,
    WorkerSettings,
    canonical,
    digest,
    finite,
    integer,
)

VERSION = "b4a-openmp-owner-v1"


class _BudgetView:
    """Share the parent ledger; process CPU already contains native threads."""

    def __init__(self, budget):
        self.budget = budget

    def consume(self, amount=1):
        self.budget.consume(amount)

    def poll(self):
        poll_budget(self.budget)


class NativeECMExecutor:
    """Experimental owner retaining immutable metadata within its life.

    Stable B12 assignments are evaluated in tier order, at most capacity curves
    live at once. Setup, strict recovery and classical stage two run in the
    Python owner. Only funded stage-one packets enter OpenMP. Pauses drain C,
    retain completed curves and charges, and restart unfinished curve IDs on
    resume; partial native/curve state is deliberately never serialized.
    """

    def __init__(
        self, *, seed=0, settings=None, schedule="dynamic", capacity=None
    ):
        self.seed = integer(seed, "method seed", 0, 2**256 - 1)
        self.settings = settings or WorkerSettings(workers=2)
        if self.settings.workers not in (1, 2, 4) or schedule not in SCHEDULES:
            raise ValueError("unsupported OpenMP team configuration")
        if (
            self.settings.stage_policy != "whole"
            or self.settings.assignment_policy != "dynamic"
        ):
            raise ValueError(
                "native owner requires whole-curve dynamic admission"
            )
        self.capacity = integer(
            2 * self.settings.workers if capacity is None else capacity,
            "team capacity",
            1,
            32,
        )
        self.schedule = schedule
        # Reserve simultaneous Python jobs, metadata, transport, OpenMP stacks
        # and parent ownership before any native allocation or worker launch.
        self.owned_bytes = (
            self.settings.parent_memory
            + TEAM_OVERHEAD
            + self.capacity * (self.settings.worker_memory + 2 * 2**20)
            + self.settings.metadata_slots
            * (
                self.settings.worker_memory
                + self.settings.decoded_cache_bytes
                + self.settings.packed_cache_bytes
                + self.settings.stage_two_cache_bytes
                + 4096 * self.settings.bind_context
            )
            + 3 * self.settings.checkpoint_bytes
        )
        if self.owned_bytes > self.settings.memory_bytes:
            raise MemoryError("simultaneous native-team ownership refused")
        self.engine = CurveEngine(self.settings)
        self.team = None
        self.lock = threading.Lock()
        self.closed = False
        self.lifetime_stats = {}
        self.identity = self._identity()

    def _identity(self):
        source = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        team_source = hashlib.sha256(
            Path(__file__)
            .parents[1]
            .joinpath("ecm/native_team.py")
            .read_bytes()
        ).hexdigest()
        return digest(
            [
                VERSION,
                source_identity(),
                source,
                team_source,
                self.seed,
                asdict(self.settings),
                self.schedule,
                self.capacity,
            ]
        )

    def _team(self):
        if self.team is None:
            self.team = NativeTeam(
                workers=self.settings.workers,
                schedule=self.schedule,
                capacity=self.capacity,
                memory_bytes=TEAM_OVERHEAD + self.capacity * 2 * 2**20,
            )
        return self.team

    def start(self, *, budget):
        """Explicit startup diagnostic, including first OpenMP region entry."""
        if self.closed or not self.lock.acquire(blocking=False):
            raise ValueError("native executor closed or active")
        try:
            if self.identity != self._identity():
                raise ValueError("native executor settings/source changed")
            budget.consume(2)
            self._team().execute(
                [("ladder", 1, (2, 1), 1013, 2)],
                budget,
                lambda index, result: False,
            )
        finally:
            self.lock.release()

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
        if sys.implementation.name != "pypy" or sys.version_info[:2] != (
            3,
            11,
        ):
            raise RuntimeError("native executor requires PyPy Python 3.11")
        if self.closed or not self.lock.acquire(blocking=False):
            raise ValueError("native executor closed or active")
        try:
            if self.identity != self._identity():
                raise ValueError("native executor settings/source changed")
            return self._run(
                n, config, budget, lineage, checkpoint, fixed_work
            )
        finally:
            self.lock.release()

    def _restore(self, checkpoint, manifest, base, fixed_work):
        if (
            not isinstance(checkpoint, dict)
            or len(canonical(checkpoint)) > self.settings.checkpoint_bytes
        ):
            raise ValueError("oversized native checkpoint")
        payload = checkpoint["payload"]
        overhead = checkpoint["serialization_overhead"]
        if checkpoint["sha256"] != digest(payload) or checkpoint[
            "overhead_sha256"
        ] != digest(overhead):
            raise ValueError("corrupt native checkpoint")
        if set(payload) != {
            "version",
            "executor",
            "manifest",
            "fixed_work",
            "completed",
            "attempts",
            "runs",
            "budget",
        }:
            raise ValueError("invalid native checkpoint fields")
        if (
            payload["version"],
            payload["executor"],
            payload["manifest"],
            payload["fixed_work"],
        ) != (VERSION, self.identity, manifest.identity, fixed_work):
            raise ValueError(
                "native source/team/schedule/search changed on resume"
            )
        if (
            type(payload["fixed_work"]) is not bool
            or type(fixed_work) is not bool
            or set(overhead) != {"wall", "cpu"}
        ):
            raise ValueError("invalid native checkpoint policy/overhead")
        completed, attempts = payload["completed"], payload["attempts"]
        if (
            not isinstance(completed, dict)
            or not isinstance(attempts, list)
            or len(attempts) != len(manifest.jobs)
        ):
            raise ValueError("invalid native checkpoint traversal")
        for count in attempts:
            integer(count, "curve attempts", 0, self.settings.max_attempts)
        for key, result in completed.items():
            position = integer(
                int(key), "completed curve", 0, len(attempts) - 1
            )
            assignment = manifest.assignment(position)
            if (
                key != str(position)
                or set(result) != {"id", "sigma", "factor", "reason"}
                or not attempts[position]
            ):
                raise ValueError("invalid native completed record")
            if (
                result["id"] != assignment["id"]
                or result["sigma"] != assignment["sigma"]
                or result["reason"] != "complete"
            ):
                raise ValueError("native completed assignment changed")
            if result["factor"] is not None and not utils.valid_divisor(
                result["factor"], manifest.n
            ):
                raise ValueError("invalid native checkpoint divisor")
        integer(payload["runs"], "native runs", 1, self.settings.max_runs)
        saved = dict(payload["budget"])
        if set(saved) != {"work", "wall", "cpu"}:
            raise ValueError("invalid native resource record")
        integer(saved["work"], "saved work")
        for name in ("wall", "cpu"):
            saved[name] = finite(saved[name], "saved " + name) + finite(
                overhead[name], "serialization " + name
            )
        if (
            saved["work"] > base.work_limit
            or saved["wall"] > base.seconds
            or saved["cpu"] > base.cpu_seconds
        ):
            raise ValueError("resume grant below native prior charges")
        if base.used == base.prior_wall == base.prior_cpu == 0:
            base.used = saved["work"]
            base.prior_wall, base.prior_cpu = saved["wall"], saved["cpu"]
        elif (
            base.used < saved["work"]
            or base.wall_used < saved["wall"]
            or base.cpu_used < saved["cpu"]
        ):
            raise ValueError("native resume would reset cumulative charges")
        return dict(completed), list(attempts), payload["runs"]

    def _run(self, n, config, budget, lineage, checkpoint, fixed_work):
        from . import stage_jobs as jobs

        integer(
            n, "native cofactor", 3, 2 ** min(1024, config.max_input_bits) - 1
        )
        if (
            not n % 2
            or config.backend != "python-int"
            or config.ecm_kernel != "c-gmp"
        ):
            raise ValueError("native team requires odd Python-int/c-gmp input")
        if (
            config.ecm_pair_distance is not None
            or config.ecm_pair_wheel is not None
        ):
            raise ValueError(
                "native owner currently supports classical continuation"
            )
        if type(fixed_work) is not bool:
            raise ValueError("fixed_work must be Boolean")
        base = _base_budget(budget)
        if base.seconds is None or base.cpu_seconds is None:
            raise ValueError("finite native wall/CPU grants required")
        if (
            self.settings.assignment_cpu < base.cpu_seconds
            or self.settings.assignment_work < base.work_limit
        ):
            raise ValueError(
                "assignment caps must cover the native global grant"
            )
        view = _BudgetView(budget)
        options = arithmetic_options(config, self.settings)
        manifest = Manifest(n, self.seed, lineage, options)
        output_bound = 8192 + 2048 * len(manifest.jobs)
        if output_bound > self.settings.checkpoint_bytes:
            raise MemoryError(
                "native checkpoint ownership refused before dispatch"
            )
        completed, attempts, runs = {}, [0] * len(manifest.jobs), 0
        if checkpoint is not None:
            completed, attempts, runs = self._restore(
                checkpoint, manifest, base, fixed_work
            )
        if runs >= self.settings.max_runs:
            raise ValueError("native resume run limit")
        runs += 1
        # Output and manifest preparation remain charged even on cancellation.
        # There are no refundable child leases: all CPU is this process's CPU.
        active, curve_work = {}, {}
        reason, factor = "curves_exhausted", None
        start_wall, start_cpu = time.monotonic(), time.process_time()
        slices_before = 0 if self.team is None else self.team.slices
        lost_chunks, chunks = 0, 0

        def finish(position, job):
            nonlocal factor
            if not job["done"]:
                return False
            assignment = manifest.assignment(position)
            if job.get("sigma") != assignment["sigma"]:
                raise ValueError("native curve identity changed")
            divisor = job["factor"]
            if divisor is not None and not utils.valid_divisor(divisor, n):
                raise ValueError("native owner produced an invalid divisor")
            completed[str(position)] = dict(
                id=assignment["id"],
                sigma=assignment["sigma"],
                factor=divisor,
                reason="complete",
            )
            active.pop(position)
            curve_work.pop(position)
            if divisor is not None and factor is None:
                factor = divisor
            return factor is not None and not fixed_work

        try:
            view.consume(output_bound + len(canonical(options)))
            if (
                self.team is not None
                and self.engine.identity != manifest.metadata
            ):
                # A changed immutable program family starts a fresh useful
                # native context/store. Parked libomp workers keep their DSO.
                self.team.close()
                self.team = None
                slices_before = 0
            self.engine._prepare(
                dict(schedule=manifest.metadata, options=options), view
            )
            context = self.engine.programs or self.engine.context
            selected = self.engine.config
            advance = jobs.advance_job
            if self.settings.stage_two_cache_bytes:
                from .ecm_stage_two import advance_prepared_job

                advance = advance_prepared_job
            # A resumed successful search needs no further arithmetic.
            factors = [
                r["factor"]
                for r in completed.values()
                if r["factor"] is not None
            ]
            if factors:
                factor = factors[0]
            while len(completed) < len(manifest.jobs) and not (
                factor is not None and not fixed_work
            ):
                view.poll()
                remaining = [
                    i
                    for i in range(len(manifest.jobs))
                    if str(i) not in completed
                ]
                tier = manifest.jobs[remaining[0]][0]
                for position in remaining:
                    if len(active) >= self.capacity:
                        break
                    if (
                        manifest.jobs[position][0] != tier
                        or position in active
                    ):
                        continue
                    if attempts[position] >= self.settings.max_attempts:
                        raise BudgetExhaustedError("attempt_limit")
                    view.consume(1)
                    assignment = manifest.assignment(position)
                    active[position] = jobs.new_job(
                        "ecm",
                        n,
                        assignment["seed"],
                        assignment["b1"],
                        assignment["b2"],
                    )
                    attempts[position] += 1
                    curve_work[position] = 0
                requests, bindings = [], []
                for position, job in list(active.items()):
                    before = base.used
                    if job["phase"] == "stage_one" and jobs.prepare_stage_one(
                        job, view, context, selected
                    ):
                        plan = jobs.stage_one_plan(job, context, view)
                        start = list(job["value"])
                        if plan is None:
                            payload = prod(power for _, power in job["powers"])
                            view.consume(payload.bit_length() + 1)
                            kind = "ladder"
                        else:
                            payload = plan.prepare(job["powers"], view)
                            kind = "records"
                        requests.append((kind, payload, start, n, job["a24"]))
                        bindings.append((position, job, plan, start, payload))
                    elif job["phase"] != "stage_one":
                        advance(job, view, context, selected)
                    curve_work[position] += base.used - before
                    if curve_work[position] > self.settings.assignment_work:
                        raise BudgetExhaustedError("assignment_work_limit")
                    if job["done"] and finish(position, job):
                        break
                if factor is not None and not fixed_work:
                    lost_chunks += len(requests)
                    break
                if requests:
                    published = 0

                    def publish(index, result):
                        nonlocal published, chunks
                        position, job, plan, start, payload = bindings[index]
                        if plan is None:
                            jobs.commit_ladder_chunk(job, start, list(result))
                        else:
                            value, divisor, replayed = plan.recover(
                                payload, start, *result, n, job["a24"]
                            )
                            jobs.commit_chain_chunk(
                                job, plan, start, value, divisor, replayed
                            )
                        published += 1
                        chunks += 1
                        return finish(position, job)

                    try:
                        self._team().execute(requests, view, publish)
                    finally:
                        lost_chunks += len(requests) - published
            if factor is not None and not fixed_work:
                reason = "factor"
            view.poll()
        except BudgetExhaustedError as error:
            reason, factor = str(error), None
        except TeamRefusedError:
            reason, factor = "native_team_refused", None
        # Active curves are deliberately restarted at their stable assignment
        # boundary on resume. All work for discarded chunks stays consumed.
        stats = dict(
            manifest=manifest.identity,
            native_threads=self.settings.workers,
            schedule=self.schedule,
            capacity=self.capacity,
            chunks=chunks,
            lost_chunks=lost_chunks,
            restarted_on_resume=len(active),
            completed=len(completed),
            owned_bytes=self.owned_bytes,
            native_slices=(0 if self.team is None else self.team.slices)
            - slices_before,
            elapsed=time.monotonic() - start_wall,
            process_cpu=time.process_time() - start_cpu,
            charged_child_cpu=0,
            reported_child_cpu=0,
        )
        for key in (
            "chunks",
            "lost_chunks",
            "native_slices",
            "elapsed",
            "process_cpu",
        ):
            self.lifetime_stats[key] = (
                self.lifetime_stats.get(key, 0) + stats[key]
            )
        serialization_wall, serialization_cpu = (
            time.monotonic(),
            time.process_time(),
        )
        payload = dict(
            version=VERSION,
            executor=self.identity,
            manifest=manifest.identity,
            fixed_work=fixed_work,
            completed=completed,
            attempts=attempts,
            runs=runs,
            budget=dict(
                work=base.used, wall=base.wall_used, cpu=base.cpu_used
            ),
        )
        checkpoint = dict(payload=payload, sha256=digest(payload))
        overhead = dict(
            wall=time.monotonic() - serialization_wall,
            cpu=time.process_time() - serialization_cpu,
        )
        checkpoint.update(
            serialization_overhead=overhead, overhead_sha256=digest(overhead)
        )
        if len(canonical(checkpoint)) > output_bound:
            raise MemoryError("native checkpoint exceeded its funded bound")
        return WorkerRun(factor, reason, checkpoint, stats)

    def close(self):
        if not self.lock.acquire(blocking=False):
            raise ValueError("cannot close an active native executor")
        try:
            if self.team is not None:
                self.team.close()
            self.engine.close()
            self.closed = True
        finally:
            self.lock.release()
