"""Quiescent current-worker snapshots with independent arithmetic replay."""

import base64
import hashlib
import json
import time

from ..common import utils
from .checkpoint import _restore_store, _store
from .current_parallel import QSWorkerConfig
from .families import _checked_resources, _checksum, _identity
from .worker_configuration import decode_config, encode_config
from .worker_packet import checked_int, decode_atom, digest, polynomial_at

SCHEMA = "a7a-current-checkpoint-v2"


def _paid_prefix(job, index, receipts, pending):
    """Prove contiguous paid coverage before a merged partial cursor."""
    _, start, end = job._spec(index)
    lo, hi = -job.config.half_width, job.config.half_width + 1
    cursor = [start, lo]
    for receipt in receipts:
        if receipt["id"] != index:
            continue
        if receipt["start"] != cursor:
            raise ValueError("paid prefix receipt chain is discontinuous")
        if receipt["status"] == "lost":
            continue
        if (index, receipt["attempt"]) in pending:
            break
        following = receipt.get("cursor")
        if (
            not isinstance(following, list)
            or len(following) != 2
            or any(type(v) is not int for v in following)
            or not start <= following[0] <= end
            or not lo <= following[1] <= hi
            or following[0] == end
            and following[1] != lo
            or tuple(following) < tuple(cursor)
        ):
            raise ValueError("invalid acknowledged prefix cursor")
        cursor = following
    if cursor != job.progress.get(index):
        raise ValueError("merged prefix lacks contiguous paid receipts")
    return tuple(cursor)


def _pack_report(report, settings):
    record = {
        k: v
        for k, v in report.items()
        if k not in ("rows", "owned_bytes", "received_at")
    }
    if settings.packed_transport:
        record["packet"] = base64.b64encode(record["packet"]).decode("ascii")
        record["packet_encoding"] = "base64"
    else:
        record["packet"] = record["packet"].decode("ascii")
    return record


def _unpack_report(report, settings):
    report = dict(report)
    if settings.packed_transport:
        if report.pop("packet_encoding", None) != "base64":
            raise ValueError("packed checkpoint encoding mismatch")
        raw = base64.b64decode(report["packet"], validate=True)
    else:
        raw = report["packet"].encode("ascii")
    return dict(report, packet=raw)


def pack(job):
    if job._running or any(r["status"] == "issued" for r in job.receipts):
        raise ValueError("checkpoint requires drained and settled workers")
    pending = [
        [i, _pack_report(r, job.workers)]
        for i, r in sorted(job.pending.items())
    ]
    payload = dict(
        schema=SCHEMA,
        n=job.n,
        seed=job.seed,
        config=encode_config(job.config),
        workers=job.workers.record(),
        identity=job.identity,
        execution=job.execution_identity,
        base=None if job.source.base is None else _identity(job.source.base),
        next_assignment=job.next_assignment,
        progress=sorted(job.progress.items()),
        attempts=sorted(job.attempts.items()),
        contexts=job.contexts,
        split_used=job.split_used,
        split_unknown=job.split_unknown,
        pending=pending,
        receipts=job.receipts,
        stats=job.stats,
        runs=job.runs,
        divisor=job.divisor,
        deferred_divisor=job.deferred_divisor,
        fixed_work=job.fixed_work,
        store=None if job.engine is None else _store(job.engine.collector),
    )
    if job.workers.arrival_order:
        payload["ahead"] = [
            [i, _pack_report(r, job.workers)]
            for i, r in sorted(job.ahead.items())
        ]
    # Eight coexisting representations are reserved in QSWorkerConfig before
    # execution. Stop encoding prospectively; never build an unbounded blob.
    chunks, size = [], 0
    for chunk in json.JSONEncoder(
        sort_keys=True, separators=(",", ":"), allow_nan=False
    ).iterencode(payload):
        size += len(chunk.encode("ascii"))
        if size > job.workers.checkpoint_bytes:
            raise MemoryError("current-worker checkpoint exceeds cap")
        chunks.append(chunk)
    blob = "".join(chunks)
    resources = dict(
        work_used=job.budget.used,
        wall_used=job.budget.wall_used,
        cpu_used=job.budget.cpu_used,
    )
    return dict(
        blob=blob,
        sha256=hashlib.sha256(blob.encode("ascii")).hexdigest(),
        resources=resources,
        resources_sha256=_checksum(resources),
    )


def _restore(cls, checkpoint, budget, workers):
    blob = checkpoint["blob"]
    limit = workers.checkpoint_bytes if workers is not None else 256 * 2**20
    if type(blob) is not str or len(blob) > limit or not blob.isascii():
        raise ValueError("current-worker checkpoint size/encoding")
    if (
        hashlib.sha256(blob.encode("ascii")).hexdigest()
        != checkpoint["sha256"]
    ):
        raise ValueError("current-worker checkpoint digest mismatch")
    resources = _checked_resources(checkpoint["resources"])
    if _checksum(resources) != checkpoint["resources_sha256"]:
        raise ValueError("checkpoint resource digest mismatch")
    saved = json.loads(blob)
    if saved["schema"] != SCHEMA:
        raise ValueError("unsupported current-worker checkpoint")
    settings = QSWorkerConfig(**saved["workers"])
    if workers is not None and workers != settings:
        raise ValueError(
            "resume execution settings changed, including workers"
        )
    if len(blob) > settings.checkpoint_bytes:
        raise ValueError("checkpoint exceeds its declared cap")
    budget.used = max(budget.used, resources["work_used"])
    budget.prior_wall = max(budget.prior_wall, resources["wall_used"])
    budget.prior_cpu = max(budget.prior_cpu, resources["cpu_used"])
    if budget.used > budget.work_limit:
        raise ValueError("restored charges exceed work allowance")
    job = cls(
        saved["n"],
        seed=saved["seed"],
        config=decode_config(saved["config"]),
        workers=settings,
        budget=budget,
    )
    if (
        saved["identity"] != job.identity
        or saved["execution"] != job.execution_identity
    ):
        raise ValueError("checkpoint job identity mismatch")
    job._setup()
    if saved["base"] != (
        None if job.source.base is None else _identity(job.source.base)
    ):
        raise ValueError("checkpoint factor base differs from reconstruction")
    next_index = checked_int(
        saved["next_assignment"], "next assignment", 0, job.assignment_count
    )
    job.next_assignment = next_index
    job.runs = checked_int(saved["runs"], "runs", 0, settings.max_runs)
    if (
        saved["fixed_work"] is not None
        and type(saved["fixed_work"]) is not bool
    ):
        raise ValueError("invalid fixed-work checkpoint")
    job.fixed_work = saved["fixed_work"]
    if settings.arrival_order and job.fixed_work:
        raise ValueError("arrival-order checkpoint cannot be fixed-work")
    raw_ahead = saved.get("ahead", []) if settings.arrival_order else []
    if (
        not isinstance(raw_ahead, list)
        or len(raw_ahead) >= settings.credit_limit
    ):
        raise ValueError("arrival frontier witness exceeds credit")
    ahead_ids = set()
    for index, _ in raw_ahead:
        checked_int(
            index,
            "ahead assignment",
            next_index + 1,
            min(job.assignment_count, next_index + settings.credit_limit) - 1,
        )
        if index in ahead_ids:
            raise ValueError("duplicate arrival frontier witness")
        ahead_ids.add(index)
    for name in ("divisor", "deferred_divisor"):
        value = saved[name]
        if value is not None and not utils.valid_divisor(value, job.n):
            raise ValueError("invalid checkpoint divisor")
        setattr(job, name, value)
    receipts = saved["receipts"]
    if (
        not isinstance(receipts, list)
        or len(receipts) > settings.receipt_limit
    ):
        raise ValueError("checkpoint receipt limit")
    charged, attempts, acknowledged = 0, {}, {}
    for receipt in receipts:
        index = checked_int(
            receipt["id"],
            "receipt assignment",
            0,
            max(0, job.assignment_count - 1),
        )
        attempt = checked_int(
            receipt["attempt"], "receipt attempt", 1, settings.max_attempts
        )
        if attempt != attempts.get(index, 0) + 1:
            raise ValueError("receipt attempt sequence changed")
        if receipt["token"] != digest(
            [
                job.execution_identity,
                index,
                attempt,
                receipt["start"],
                receipt["split_grant"],
            ]
        ):
            raise ValueError("receipt token differs from assignment identity")
        attempts[index] = attempt
        grant = checked_int(
            receipt["work"], "lease work", 1, settings.assignment_work
        )
        spent = checked_int(receipt["charged_work"], "charged work", 0, grant)
        if receipt["status"] not in ("reported", "lost") or (
            receipt["status"] == "lost" and spent != grant
        ):
            raise ValueError("unsettled or refunded lost lease")
        charged += spent + checked_int(
            receipt.get("controller_work", 0),
            "controller tail",
            0,
            settings.packet_bytes + 1024,
        )
        if receipt["status"] == "reported":
            acknowledged[(index, attempt)] = receipt
    if charged > resources["work_used"]:
        raise ValueError("checkpoint loses retained work charges")
    if settings.arrival_order:
        sequences = [
            checked_int(
                r["arrival_seq"],
                "arrival sequence",
                0,
                settings.receipt_limit - 1,
            )
            for r in acknowledged.values()
        ]
        if sorted(sequences) != list(range(len(sequences))) or saved["stats"][
            "received_packets"
        ] != len(sequences):
            raise ValueError("arrival history is incomplete or repeated")
    split_total = sum(r["split_calls"] or 0 for r in receipts)
    if (
        split_total != saved["split_used"]
        or type(saved["split_unknown"]) is not bool
    ):
        raise ValueError("checkpoint splitting charges changed")
    limit = getattr(job.config.collector, "split_call_limit", None)
    if limit is not None and split_total > limit:
        raise ValueError("checkpoint exceeds cumulative split quota")
    for row in receipts:
        grant = row["split_grant"]
        if (grant is None) != (limit is None):
            raise ValueError("split lease differs from collector quota")
        if grant is not None:
            checked_int(
                grant, "split grant", 0, job.workers.split_calls_per_task
            )
            checked_int(row["split_calls"], "split charge", 0, grant)
            if row["status"] == "lost" and row["split_calls"] != grant:
                raise ValueError("lost splitter grant was refunded")
    unknown = any(
        r["status"] == "lost" and r["split_grant"] is None for r in receipts
    )
    if unknown != saved["split_unknown"]:
        raise ValueError("missing splitter report uncertainty was reset")
    job.split_used, job.split_unknown = split_total, unknown
    if job.engine is not None:
        job.engine.collector._split_calls = split_total
    job.receipts = receipts
    for name in ("progress", "attempts"):
        records = saved[name]
        if (
            not isinstance(records, list)
            or len(records) > settings.credit_limit
        ):
            raise ValueError("unbounded checkpoint traversal")
        values = dict(records)
        if len(values) != len(records):
            raise ValueError("duplicate traversal key")
        for index in values:
            checked_int(
                index,
                "pending assignment",
                next_index,
                min(job.assignment_count, next_index + settings.credit_limit)
                - 1,
            )
        setattr(job, name, values)
    if job.attempts != {
        i: count
        for i, count in attempts.items()
        if i >= next_index and i not in ahead_ids
    }:
        raise ValueError("attempt cursors disagree with paid history")
    for index, cursor in job.progress.items():
        _, start, end = job._spec(index)
        if (
            not isinstance(cursor, list)
            or len(cursor) != 2
            or type(cursor[1]) is not int
        ):
            raise ValueError("invalid resumed traversal cursor")
        checked_int(cursor[0], "resumed Gray", start, end)
        if (
            not -job.config.half_width
            <= cursor[1]
            <= job.config.half_width + 1
        ):
            raise ValueError("resumed position outside window")
        pending_keys = {(i, r["attempt"]) for i, r in saved["pending"]}
        prior = [
            r
            for r in receipts
            if r["id"] == index
            and r["status"] == "reported"
            and (index, r["attempt"]) not in pending_keys
        ]
        if not prior or prior[-1]["cursor"] != cursor:
            raise ValueError("cursor has no paid publication")
    pending = saved["pending"]
    if not isinstance(pending, list) or len(pending) > settings.credit_limit:
        raise ValueError("checkpoint pending bound")
    for index, report in pending:
        checked_int(
            index,
            "pending ID",
            next_index,
            min(job.assignment_count, next_index + settings.credit_limit) - 1,
        )
        if index in job.pending:
            raise ValueError("duplicate pending packet")
        receipt = acknowledged.get((index, report["attempt"]))
        if (
            receipt is None
            or receipt["token"] != report["token"]
            or receipt["charged_work"] != report["work"]
        ):
            raise ValueError("pending packet lacks paid receipt")
        task = job._task_for_restore(index, report, receipt)
        report = _unpack_report(report, settings)
        checked, _ = job._validate_report(report, task)
        checked["admitted"] = checked_int(
            report["admitted"], "admission cursor", 0, len(checked["rows"])
        )
        if (
            index != next_index
            and checked["admitted"]
            and not (settings.arrival_order and report["reason"] == "complete")
        ):
            raise ValueError("unordered checkpoint admission")
        if (
            settings.arrival_order
            and checked_int(
                report.get("arrival_seq"),
                "arrival sequence",
                0,
                settings.receipt_limit - 1,
            )
            != receipt["arrival_seq"]
        ):
            raise ValueError("pending packet arrival sequence changed")
        job.pending[index] = checked
        checked["received_at"] = time.monotonic()
    if len(job.pending) + len(raw_ahead) > settings.credit_limit:
        raise ValueError("arrival witnesses exceed aggregate credit")
    for index, report in raw_ahead:
        receipt = acknowledged.get((index, report["attempt"]))
        if (
            index in job.pending
            or index in job.progress
            or receipt is None
            or receipt["token"] != report["token"]
            or receipt["charged_work"] != report["work"]
            or checked_int(
                report.get("arrival_seq"),
                "arrival sequence",
                0,
                settings.receipt_limit - 1,
            )
            != receipt["arrival_seq"]
        ):
            raise ValueError("arrival witness lacks its paid acknowledgement")
        task = job._task_for_restore(index, report, receipt)
        report = _unpack_report(report, settings)
        checked, _ = job._validate_report(report, task)
        if report["reason"] != "complete" or checked_int(
            report["admitted"], "ahead admission", 0, len(checked["rows"])
        ) != len(checked["rows"]):
            raise ValueError("unfinished arrival witness")
        checked["admitted"] = len(checked["rows"])
        checked["received_at"] = time.monotonic()
        job.ahead[index] = checked
    if job.engine is not None and saved["store"] is not None:
        _restore_store(saved["store"], job.engine.collector, budget)
        if job.engine.collector._split_calls != job.split_used:
            raise ValueError("store resets cumulative splitting charges")
    contexts = saved["contexts"]
    if (
        not isinstance(contexts, dict)
        or len(contexts) > settings.context_limit
    ):
        raise ValueError("checkpoint context bound")
    verified = {}
    for identity, pair in contexts.items():
        if not isinstance(pair, list) or len(pair) != 2:
            raise ValueError("invalid polynomial context")
        index, gray = pair
        checked_int(
            index,
            "context assignment",
            0,
            min(
                job.assignment_count - 1,
                next_index
                + (settings.credit_limit - 1 if settings.arrival_order else 0),
            ),
        )
        primes, start, end = job._spec(index)
        checked_int(gray, "context Gray", start, end - 1)
        poly = polynomial_at(job.source.base, primes, gray, job.config, budget)
        if poly.identity != identity:
            raise ValueError("polynomial context binding failed")
        verified[identity] = poly
    job.contexts = contexts
    prefix_ids = set()
    for index, report in list(job.pending.items()) + list(job.ahead.items()):
        primes, _, _ = job._spec(index)
        for offset, row in enumerate(report["rows"]):
            poly = polynomial_at(
                job.source.base, primes, row[0], job.config, budget
            )
            atom = decode_atom(
                row, poly, job.config.collector, job.source.base, budget
            )
            if (
                index == next_index or settings.arrival_order
            ) and offset < report["admitted"]:
                prefix_ids.add(atom.relation_id)
    if job.engine is not None:
        collector = job.engine.collector
        atoms = list(collector._atoms.values())
        staging = getattr(collector, "_opaque_staging", None)
        if staging is not None:
            atoms.extend(staging.atoms.values())
        extension = any(
            name in settings.record()
            for name in (
                "lease_budget",
                "reuse_state",
                "packed_transport",
                "verification",
                "lookahead",
                "whole_family",
                "arrival_order",
            )
        )
        pending_receipts = {(i, r["attempt"]) for i, r in saved["pending"]}
        paid = (
            {
                i: _paid_prefix(job, i, receipts, pending_receipts)
                for i in job.progress
            }
            if extension
            else {}
        )
        for atom in atoms:
            poly = verified.get(atom.polynomial.identity)
            if poly != atom.polynomial:
                raise ValueError("store lacks checked assignment provenance")
            index, _ = contexts[poly.identity]
            if index >= next_index and atom.relation_id not in prefix_ids:
                _, start, _ = job._spec(index)
                if (
                    index not in paid
                    or not (start, -job.config.half_width)
                    <= (contexts[poly.identity][1], atom.position)
                    < paid[index]
                ):
                    raise ValueError(
                        "store includes unpublished current assignment"
                    )
    if set(saved["stats"]) != set(job.stats):
        raise ValueError("checkpoint statistics shape")
    for name, value in saved["stats"].items():
        if type(value) not in (int, float) or not math_is_nonnegative(value):
            raise ValueError("invalid checkpoint statistics")
    if saved["stats"]["assignments"] != next_index + len(job.ahead):
        raise ValueError("checkpoint progress disagrees with statistics")
    job.stats = saved["stats"]
    return job


def math_is_nonnegative(value):
    return 0 <= value < float("inf")


def restore(cls, checkpoint, *, budget, workers=None):
    try:
        return _restore(cls, checkpoint, budget, workers)
    except (
        KeyError,
        TypeError,
        IndexError,
        AttributeError,
        RecursionError,
    ) as error:
        raise ValueError("malformed current-worker checkpoint") from error
