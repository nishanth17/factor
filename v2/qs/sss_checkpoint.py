"""Compact SSS provenance and charged assignment/solver reconstruction."""

import hashlib
import json
from dataclasses import asdict
from math import isfinite

from .. import utils
from .checkpoint import _restore_store, _solver_digest, _store
from .extraction import (
    DependencyExtractor,
    extract_dependency,
    prepare_relations,
)
from .families import _checked_resources, _checksum, _identity
from .linear_algebra import MAX_MATRIX_ROWS, DependencySolver, filter_matrix
from .sieve_collector import SieveConfig

VERSION = 2
MAX_BLOB_BYTES = 1024 * 1024


def pack_job(job):
    """Encode a checked store once; regenerate candidate batches on resume."""
    if job.config.metadata_reserve > job.config.memory_bytes:
        raise MemoryError("SSS checkpoint reserve exceeds memory_bytes")
    engine = job.pipeline
    divisor = job.divisor or (engine.divisor if engine else None)
    if divisor is not None:
        engine = None
    progress = store = None
    if engine is not None:
        collector, solver = engine.collector, engine.solver
        store = _store(collector)
        progress = dict(
            next_position=engine.next_position,
            final_solve_done=engine.final_solve_done,
            storage_reason=engine.storage_reason,
            storage_solve_done=engine.storage_solve_done,
            last_count=engine.last_count,
            last_solved_count=engine.last_solved_count,
            stats=engine.stats,
            workspace_peak=engine._workspace_peak,
            assignment=None
            if collector._assignment is None
            else dict(
                index=collector._round,
                cursor=collector._cursor,
                identity=_checksum(collector._assignment),
            ),
            solver=None
            if solver is None
            else dict(
                actions=solver.next_row + solver.xors,
                digest=_solver_digest(solver),
                pending=solver.pending is not None,
                extract_index=None
                if engine.extractor is None
                else engine.extractor.next_dependency,
            ),
        )
    payload = dict(
        version=VERSION,
        n=job.n,
        seed=job.seed,
        config=asdict(job.config),
        divisor=divisor,
        setup_memory_refused=job._setup_memory_refused,
        base_identity=_identity(job.base) if engine else None,
        engine=progress,
        store=store,
    )
    parts, size = [], 0
    encoder = json.JSONEncoder(sort_keys=True, separators=(",", ":"))
    for part in encoder.iterencode(payload):
        size += len(part.encode())
        if size + 1024 > job.config.checkpoint_bytes:
            raise MemoryError("SSS checkpoint exceeds checkpoint_bytes")
        parts.append(part)
    blob = "".join(parts)
    resources = _checked_resources(
        dict(
            work_used=job.budget.used,
            wall_used=job.budget.wall_used,
            cpu_used=job.budget.cpu_used,
        )
    )
    return dict(
        version=VERSION,
        blob=blob,
        sha256=hashlib.sha256(blob.encode()).hexdigest(),
        resources=resources,
        resources_sha256=_checksum(resources),
    )


def _restore_solver(engine, prefix, resources):
    """Rebuild the verified matrix and replay only its retained prefix."""
    budget, collector = engine.budget, engine.collector
    relations = collector.matrix_relations
    prepared = prepare_relations(
        relations,
        collector.factor_base,
        collector._atoms,
        budget=budget,
        memory_bytes=engine.config.memory_bytes,
        retained_workspace_bytes=collector._workspace,
    )
    live = collector._workspace + prepared.workspace_bytes
    live -= prepared.shared_workspace_bytes
    matrix = filter_matrix(
        prepared.rows,
        weight_two=engine.weight_two,
        budget=budget,
        memory_bytes=max(0, engine.config.memory_bytes - live),
    )
    solver = DependencySolver(matrix, budget=budget)
    actions = utils.require_integer(prefix["actions"], "solver actions", 0)
    if actions > min(
        resources["work_used"], MAX_MATRIX_ROWS**2 + MAX_MATRIX_ROWS
    ):
        raise ValueError("SSS solver prefix exceeds its bound")
    for _ in range(actions):
        if solver.next_row == len(matrix.rows):
            raise ValueError("SSS solver prefix exceeds its matrix")
        solver.step()
    if type(prefix["pending"]) is not bool:
        raise ValueError("invalid SSS solver pending flag")
    if (
        prefix["pending"]
        and solver.pending is None
        and solver.next_row < len(matrix.rows)
    ):
        solver.pending = (
            matrix.rows[solver.next_row],
            matrix.masks[solver.next_row],
        )
    if _solver_digest(solver) != prefix["digest"]:
        raise ValueError("SSS elimination prefix mismatch")
    engine.prepared, engine.solver = prepared, solver
    index = prefix["extract_index"]
    if index is not None:
        utils.require_integer(index, "extraction index", 0)
        if solver.next_row != len(matrix.rows) or index > len(
            solver.dependencies
        ):
            raise ValueError("SSS extraction prefix exceeds its kernels")
        extractor = DependencyExtractor(
            prepared, tuple(solver.dependencies), budget=budget
        )
        for mask in solver.dependencies[:index]:
            trial = extract_dependency(prepared, mask, budget=budget)
            if trial.divisor is not None:
                raise ValueError("SSS checkpoint skipped a proper factor")
            extractor.trials.append(trial)
        extractor.next_dependency = index
        engine.extractor = extractor
    engine._workspace_peak = max(
        engine._workspace_peak, live + matrix.workspace_bytes
    )


def restore_job(checkpoint, *, budget, config=None):
    """Validate all arithmetic and retain prior resources before rebuilding."""
    from .sss import SSSConfig, SSSJob

    if not isinstance(checkpoint, dict) or set(checkpoint) != {
        "version",
        "blob",
        "sha256",
        "resources",
        "resources_sha256",
    }:
        raise ValueError("invalid SSS checkpoint envelope")
    if type(checkpoint["version"]) is not int or (
        checkpoint["version"] not in (1, VERSION)
    ):
        raise ValueError("unsupported SSS checkpoint version")
    blob = checkpoint["blob"]
    if (
        not isinstance(blob, str)
        or len(blob) > MAX_BLOB_BYTES
        or (len(blob.encode()) > MAX_BLOB_BYTES)
    ):
        raise ValueError("SSS checkpoint exceeds its hard byte cap")
    if hashlib.sha256(blob.encode()).hexdigest() != checkpoint["sha256"]:
        raise ValueError("SSS checkpoint integrity mismatch")
    resources = _checked_resources(checkpoint["resources"])
    if _checksum(resources) != checkpoint["resources_sha256"]:
        raise ValueError("SSS checkpoint resource integrity mismatch")
    payload = json.loads(blob)
    if (
        type(payload["version"]) is not int
        or payload["version"] != checkpoint["version"]
    ):
        raise ValueError("invalid SSS checkpoint payload version")
    if payload["version"] >= 2 and payload["store"] is not None:
        if not isinstance(payload["store"], dict) or (
            "row_order" not in payload["store"]
        ):
            raise ValueError("mixed-order checkpoint lacks row_order")
    values = dict(payload["config"])
    values["collector"] = SieveConfig(**values["collector"])
    saved_config = SSSConfig(**values)
    if config is not None and config != saved_config:
        raise ValueError("SSS checkpoint configuration mismatch")
    config = saved_config
    if len(blob.encode()) + 1024 > config.checkpoint_bytes:
        raise ValueError("SSS checkpoint exceeds configured byte cap")
    if budget.work_limit < resources["work_used"]:
        raise ValueError("resume allowance is below consumed work")
    if budget.used == 0 and budget.prior_wall == 0 and budget.prior_cpu == 0:
        budget.used = resources["work_used"]
        budget.prior_wall = resources["wall_used"]
        budget.prior_cpu = resources["cpu_used"]
    elif (
        budget.used < resources["work_used"]
        or budget.wall_used < resources["wall_used"]
        or budget.cpu_used < resources["cpu_used"]
    ):
        raise ValueError("resume must retain consumed resources")
    budget.consume(0)
    job = SSSJob(
        payload["n"], seed=payload["seed"], config=config, budget=budget
    )
    job.divisor = payload["divisor"]
    if job.divisor is not None and not utils.valid_divisor(job.divisor, job.n):
        raise ValueError("invalid SSS checkpoint divisor")
    if type(payload["setup_memory_refused"]) is not bool:
        raise ValueError("invalid SSS setup refusal flag")
    job._setup_memory_refused = payload["setup_memory_refused"]
    saved = payload["engine"]
    if saved is None:
        if (
            payload["store"] is not None
            or payload["base_identity"] is not None
        ):
            raise ValueError("SSS checkpoint store has no engine")
        return job
    if job.divisor is not None or job._setup_memory_refused:
        raise ValueError("finished SSS checkpoint retains an active engine")
    job._setup()
    if job.base is None or _identity(job.base) != payload["base_identity"]:
        raise ValueError("SSS checkpoint base identity mismatch")
    engine, collector = job.pipeline, job.pipeline.collector
    _restore_store(payload["store"], collector, budget)
    for name in ("next_position", "last_count", "last_solved_count"):
        value = utils.require_integer(
            saved[name], name, -1 if name.startswith("last") else 0
        )
        limit = (
            config.search_rounds
            if name == "next_position"
            else (len(collector._full) + len(collector._combined))
        )
        if value > limit:
            raise ValueError("SSS checkpoint progress exceeds its bound")
        setattr(engine, name, value)
    collector._last_stop = engine.next_position
    for name in ("final_solve_done", "storage_solve_done"):
        if type(saved[name]) is not bool:
            raise ValueError("invalid SSS completion flag")
        setattr(engine, name, saved[name])
    if engine.final_solve_done and engine.next_position != engine.hi:
        raise ValueError("SSS final solve precedes assignment exhaustion")
    engine.storage_reason = saved["storage_reason"]
    if engine.storage_reason not in (
        None,
        "relation_limit",
        "atom_limit",
        "memory_limit",
    ) or (engine.storage_solve_done and engine.storage_reason is None):
        raise ValueError("invalid SSS storage stop")
    stats = saved["stats"]
    required = {
        "filter_calls",
        "solve_calls",
        "trivial_dependencies",
        "collected_positions",
        "stage_seconds",
        "collector_counts",
    }
    if not isinstance(stats, dict) or not required <= stats.keys():
        raise ValueError("invalid SSS engine statistics")
    for name in required - {"stage_seconds", "collector_counts"}:
        utils.require_integer(stats[name], name, 0)
    for name in ("stage_seconds", "collector_counts"):
        values = stats[name]
        if not isinstance(values, dict) or len(values) > 32:
            raise ValueError("invalid SSS statistic map")
        for key, value in values.items():
            if not isinstance(key, str) or len(key) > 80:
                raise ValueError("invalid SSS statistic name")
            if name == "collector_counts":
                utils.require_integer(value, key, 0)
            elif (
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not isfinite(value)
                or value < 0
            ):
                raise ValueError("invalid SSS stage duration")
    engine.stats = stats
    peak = utils.require_integer(saved["workspace_peak"], "workspace peak", 0)
    if peak > engine.config.memory_bytes:
        raise ValueError("SSS workspace peak exceeds its bound")
    engine._workspace_peak = max(peak, collector._workspace)
    assignment = saved["assignment"]
    if assignment is not None:
        if assignment["index"] != engine.next_position:
            raise ValueError("SSS assignment disagrees with collection cursor")
        stats = dict.fromkeys(
            ("filter_rejections", "generated_candidates", "tree_rejections"),
            0,
        )
        collector._prepare(engine.next_position, stats)
        if _checksum(collector._assignment) != assignment["identity"]:
            raise ValueError("SSS assignment identity mismatch")
        cursor = utils.require_integer(assignment["cursor"], "cursor", 0)
        if cursor > len(collector._assignment):
            raise ValueError("SSS candidate cursor exceeds its assignment")
        collector._cursor = cursor
    if saved["solver"] is not None:
        _restore_solver(engine, saved["solver"], resources)
    return job
