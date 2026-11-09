"""Stable assignment/publication owner for the whole-ECM native coroutine."""

import time

from ..common import utils
from ..ecm.native_curve_team import NativeCurveTeam
from ..ecm.native_resources import poll_native_budget
from .budget import BudgetExhaustedError
from .ecm_native import VERSION, _BudgetView
from .ecm_workers import WorkerRun, _base_budget, arithmetic_options
from .resources import RESOURCE_VERSION
from .stage_jobs import new_job
from .worker_contract import Manifest, canonical, digest, integer


def run_full(owner, n, config, budget, lineage, checkpoint, fixed_work):
    """Retain the manifest, history and stable ready publication."""
    integer(n, "native cofactor", 3, 2 ** min(1024, config.max_input_bits) - 1)
    if (
        not n % 2
        or config.backend != "python-int"
        or config.ecm_kernel != "c-gmp-full"
    ):
        raise ValueError(
            "whole native ECM requires odd Python-int/c-gmp-full input"
        )
    if (
        config.ecm_pair_distance is not None
        or config.ecm_pair_wheel is not None
    ):
        raise ValueError(
            "whole native ECM supports classical continuation only"
        )
    if type(fixed_work) is not bool:
        raise ValueError("fixed_work must be Boolean")
    base, view = _base_budget(budget), _BudgetView(budget)
    options = arithmetic_options(config, owner.settings)
    manifest = Manifest(n, owner.seed, lineage, options)
    completed, attempts, runs = {}, [0] * len(manifest.jobs), 0
    if checkpoint is not None:
        completed, attempts, runs = owner._restore(
            checkpoint, manifest, base, fixed_work
        )
    runs += 1
    active = {}
    reason, factor = "curves_exhausted", None
    elapsed, cpu = time.monotonic(), time.process_time()
    counters = dict(setup_ops=0, stage_one_ops=0, stage_two_ops=0)
    phase_ns = [0, 0, 0]
    slices_before = 0
    lost = 0

    def capture(position):
        packet = owner.team.bindings.get(position)
        if packet is None:
            return
        info = packet._read()
        for key in counters:
            counters[key] += int(getattr(info, key))
        for index, value in enumerate(info.phase_ns):
            phase_ns[index] += int(value)

    try:
        view.consume(len(canonical(options)))
        if owner.team is not None and (
            not isinstance(owner.team, NativeCurveTeam)
            or owner.engine.identity != manifest.metadata
        ):
            owner.team.close()
            owner.team = None
        owner.engine._prepare(
            dict(schedule=manifest.metadata, options=options), view
        )
        context = owner.engine.programs or owner.engine.context
        selected = owner.engine.config
        if owner.team is None:
            owner.team = NativeCurveTeam(
                workers=owner.settings.workers,
                schedule=owner.schedule,
                capacity=owner.capacity,
                chunk_size=owner.chunk_size,
            )
        slices_before = owner.team.slices
        factors = [
            record["factor"]
            for record in completed.values()
            if record["factor"] is not None
        ]
        if factors:
            factor = factors[0]
        while len(completed) < len(manifest.jobs) and not (
            factor is not None and not fixed_work
        ):
            poll_native_budget(budget, owner.team)
            remaining = [
                position
                for position in range(len(manifest.jobs))
                if str(position) not in completed
            ]
            tier = manifest.jobs[remaining[0]][0]
            for position in remaining:
                if len(active) >= owner.capacity:
                    break
                if position in active or manifest.jobs[position][0] != tier:
                    continue
                if attempts[position] >= owner.settings.max_attempts:
                    raise BudgetExhaustedError("attempt_limit")
                view.consume(1)
                assignment = manifest.assignment(position)
                active[position] = new_job(
                    "ecm",
                    n,
                    assignment["seed"],
                    assignment["b1"],
                    assignment["b2"],
                )
                attempts[position] += 1
            ready = owner.team.advance_jobs(
                list(active.items()), view, context, selected
            )
            for position, job in ready:
                assignment = manifest.assignment(position)
                if job.get("sigma") != assignment["sigma"]:
                    raise ValueError("whole native curve identity changed")
                divisor = job["factor"]
                if divisor is not None and not utils.valid_divisor(divisor, n):
                    raise ValueError(
                        "whole native curve produced invalid divisor"
                    )
            # All ready exports have been checked. Queue order, rather than
            # completion time or worker rank, resolves simultaneous factors.
            for position, job in ready:
                assignment = manifest.assignment(position)
                capture(position)
                completed[str(position)] = dict(
                    id=assignment["id"],
                    sigma=assignment["sigma"],
                    factor=job["factor"],
                    reason="complete",
                )
                active.pop(position)
                if factor is None and job["factor"] is not None:
                    factor = job["factor"]
                if factor is not None and not fixed_work:
                    break
        if factor is not None and not fixed_work:
            reason = "factor"
        poll_native_budget(budget, owner.team)
    except BudgetExhaustedError as error:
        reason, factor = str(error), None
    except Exception as error:
        from ..ecm.native_team import TeamRefusedError

        if isinstance(error, TeamRefusedError):
            reason, factor = "native_team_refused", None
        else:
            raise
    finally:
        if isinstance(owner.team, NativeCurveTeam):
            for position in active:
                capture(position)
            lost = len(active)
            owner.team.clear_jobs()

    stats = dict(
        manifest=manifest.identity,
        native_threads=owner.settings.workers,
        schedule=owner.schedule,
        chunk_size=owner.chunk_size,
        capacity=owner.capacity,
        chunks=len(completed),
        lost_chunks=lost,
        restarted_on_resume=lost,
        completed=len(completed),
        owned_bytes=owner.owned_bytes,
        native_slices=(0 if owner.team is None else owner.team.slices)
        - slices_before,
        elapsed=time.monotonic() - elapsed,
        process_cpu=time.process_time() - cpu,
        charged_child_cpu=0,
        reported_child_cpu=0,
        native_coverage=counters,
        native_phase_seconds=[value / 1e9 for value in phase_ns],
    )
    for key in (
        "chunks",
        "lost_chunks",
        "native_slices",
        "elapsed",
        "process_cpu",
    ):
        owner.lifetime_stats[key] = (
            owner.lifetime_stats.get(key, 0) + stats[key]
        )
    wall_start, cpu_start = time.monotonic(), time.process_time()
    payload = dict(
        version=VERSION,
        executor=owner.identity,
        manifest=manifest.identity,
        fixed_work=fixed_work,
        completed=completed,
        attempts=attempts,
        runs=runs,
        budget=dict(work=base.used, wall=base.wall_used, cpu=base.cpu_used),
    )
    base.record_charges(**payload["budget"])
    checkpoint = dict(
        resource_version=RESOURCE_VERSION,
        payload=payload,
        sha256=digest(payload),
    )
    overhead = dict(
        wall=time.monotonic() - wall_start, cpu=time.process_time() - cpu_start
    )
    checkpoint.update(
        serialization_overhead=overhead, overhead_sha256=digest(overhead)
    )
    return WorkerRun(
        factor, reason, checkpoint, stats, resource_stop=base.stop
    )
