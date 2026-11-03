"""Bounded atom/boundary journals with distinct recorder failures."""

import hashlib
import importlib
import json
from contextlib import contextmanager

STREAM_BYTES = 256 * 2**20
RECORD_BYTES = 2**20
FOOTER_BYTES = 4096


class RecorderExhaustedError(RuntimeError):
    """The diagnostic recorder stopped; no production resource inference."""


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode()


class Journal:
    def __init__(
        self, path, header, *, byte_limit=STREAM_BYTES, work_limit=2**32
    ):
        if not FOOTER_BYTES + 1024 <= byte_limit <= STREAM_BYTES:
            raise ValueError("invalid journal byte allowance")
        self.byte_limit, self.work_limit = byte_limit, work_limit
        self.stream = path.open("xb", buffering=65536)
        self.bytes = self.work = self.sequence = self.admissions = 0
        self.previous = "0" * 64
        self.exhausted = False
        self.emit("header", header)

    def emit(self, kind, data, *, terminal=False):
        record = dict(
            sequence=self.sequence,
            previous=self.previous,
            kind=kind,
            data=data,
        )
        encoded = canonical(record)
        digest = hashlib.sha256(encoded).hexdigest()
        wire = canonical(dict(record, sha256=digest)) + b"\n"
        limit = self.byte_limit if terminal else self.byte_limit - FOOTER_BYTES
        if (
            len(wire) > RECORD_BYTES
            or self.bytes + len(wire) > limit
            or self.work + 8 * len(wire)
            > self.work_limit - (0 if terminal else 8 * FOOTER_BYTES)
        ):
            self.exhausted = True
            raise RecorderExhaustedError("recorder byte/event-work allowance")
        self.stream.write(wire)
        self.bytes += len(wire)
        self.work += 8 * len(wire)
        self.sequence += 1
        self.previous = digest
        return record["sequence"]

    def close(self, result, *, complete):
        try:
            self.emit(
                "footer",
                dict(result=result, complete=complete and not self.exhausted),
                terminal=True,
            )
            self.stream.flush()
        finally:
            self.stream.close()


def records(path, *, require_complete=True, budget=None):
    """Check order, hashes, record bounds and terminal framing."""
    sequence, previous, size, footer = 0, "0" * 64, 0, False
    pending = []
    with path.open("rb") as stream:
        while True:
            line = stream.readline(RECORD_BYTES + 1)
            if not line:
                break
            size += len(line)
            if footer or len(line) > RECORD_BYTES or size > STREAM_BYTES:
                raise ValueError("journal length or trailing data")
            if budget is not None:
                budget.consume(8 * len(line))
            record = json.loads(line)
            if set(record) != {
                "sequence",
                "previous",
                "kind",
                "data",
                "sha256",
            }:
                raise ValueError("invalid journal record shape")
            digest = record.pop("sha256")
            if (
                record["sequence"] != sequence
                or record["previous"] != previous
                or hashlib.sha256(canonical(record)).hexdigest() != digest
            ):
                raise ValueError("journal order or hash mismatch")
            kind = record["kind"]
            if kind not in {
                "header",
                "base",
                "window",
                "atom_begin",
                "atom_end",
                "graph",
                "solve_begin",
                "solve_end",
                "footer",
            } or not isinstance(record["data"], dict):
                raise ValueError("unknown journal event or payload")
            if sequence == 0 and kind != "header":
                raise ValueError("journal lacks header")
            if sequence and kind == "header":
                raise ValueError("duplicate journal header")
            if kind in ("atom_begin", "solve_begin"):
                if pending:
                    raise ValueError("nested collector operations")
                pending.append((kind, sequence))
            elif kind in ("atom_end", "solve_end"):
                expected = (
                    kind.replace("_end", "_begin"),
                    record["data"]["begin"],
                )
                if not pending or pending.pop() != expected:
                    raise ValueError("journal publication/boundary order")
            elif kind == "graph":
                if not pending or pending[-1][0] != "atom_begin":
                    raise ValueError("graph event outside atomic admission")
            elif kind in ("window", "base") and pending:
                raise ValueError("metadata changed during publication")
            elif kind == "footer":
                footer = True
                if require_complete and (
                    not record["data"]["complete"] or pending
                ):
                    raise ValueError(
                        "journal is an incomplete observation window"
                    )
            previous, sequence = digest, sequence + 1
            yield record
    if not footer:
        raise ValueError("journal terminal event missing")


def atom_record(atom):
    polynomial = atom.polynomial
    return dict(
        id=atom.relation_id,
        polynomial=[
            int(polynomial.a),
            int(polynomial.b),
            int(polynomial.square_coefficient),
        ],
        position=int(atom.position),
        sign=atom.sign,
        exponents=atom.exponents,
        residual=int(atom.residual),
        large_primes=atom.large_primes,
    )


@contextmanager
def observe(qs, writer):
    """Capture every eligible pre-graph atom and every relevant solve boundary.

    Empty unchanged blocks are coalesced: they cannot change replayed rows.
    Window/final/cap and pending-solver boundaries are always retained. This
    is a single-process diagnostic adapter, never a clean measurement path.
    """
    admit = qs.SieveCollector._admit
    solve = qs.QSJob._solve
    activate = qs.SIQSJob._activate
    graph_class = importlib.import_module(
        qs.__name__ + ".large_primes"
    ).LargePrimeForest
    commit = graph_class.commit
    boundary_serial = -1
    base_written = False

    def admission(collector, atom, stats):
        nonlocal base_written
        if not base_written:
            base = collector.factor_base
            writer.emit(
                "base",
                dict(
                    n=int(base.n),
                    multiplier=base.multiplier,
                    bound=base.bound,
                    primes=list(base.primes),
                ),
            )
            base_written = True
        begin = writer.emit("atom_begin", atom_record(atom))
        writer.admissions += 1
        before = dict(stats)
        try:
            result = admit(collector, atom, stats)
        except BaseException as error:
            writer.emit(
                "atom_end", dict(begin=begin, error=type(error).__name__)
            )
            raise
        writer.emit(
            "atom_end",
            dict(
                begin=begin,
                refusal=result,
                retained=atom.relation_id in collector._atoms,
                rows=len(collector._rows),
                counters={
                    k: v - before.get(k, 0)
                    for k, v in stats.items()
                    if v != before.get(k, 0)
                },
            ),
        )
        return result

    def preparation(job, final):
        nonlocal boundary_serial
        needed = (
            writer.admissions != boundary_serial
            or final
            or job.next_position >= job.hi
            or job.storage_reason is not None
            or job.solver is not None
            or job.extractor is not None
        )
        if not needed:
            return solve(job, final)
        begin = writer.emit(
            "solve_begin",
            dict(
                final=final,
                next_position=job.next_position,
                hi=job.hi,
                storage_reason=job.storage_reason,
                rows=len(job.collector._rows),
            ),
        )
        boundary_serial = writer.admissions
        try:
            result = solve(job, final)
        except BaseException as error:
            writer.emit(
                "solve_end", dict(begin=begin, error=type(error).__name__)
            )
            raise
        writer.emit(
            "solve_end",
            dict(
                begin=begin,
                divisor=result,
                last_count=job.last_count,
                filters=job.stats["filter_calls"],
                solves=job.stats["solve_calls"],
            ),
        )
        return result

    def activation(job):
        count = job.stats.get("polynomials", 0)
        result = activate(job)
        if job.stats.get("polynomials", 0) != count:
            writer.emit(
                "window",
                dict(
                    polynomial=atom_polynomial(job.engine.collector),
                    lo=job.engine.lo,
                    hi=job.engine.hi,
                    memory_bytes=job.engine.config.memory_limit,
                ),
            )
        return result

    def publication(graph, plan):
        # This event precedes mutation. atom_end distinguishes committed
        # ownership from a failed publication; incomplete streams are refused.
        writer.emit(
            "graph",
            dict(
                evicted=list(plan.evicted),
                pinned=list(plan.cycle or ()),
                closing=plan.cycle is not None,
                rebuilt=plan.rebuilt is not None,
            ),
        )
        return commit(graph, plan)

    qs.SieveCollector._admit = admission
    qs.QSJob._solve = preparation
    qs.SIQSJob._activate = activation
    graph_class.commit = publication
    try:
        yield
    finally:
        qs.SieveCollector._admit = admit
        qs.QSJob._solve = solve
        qs.SIQSJob._activate = activate
        graph_class.commit = commit


def atom_polynomial(collector):
    polynomial = collector.polynomial
    return [
        int(polynomial.a),
        int(polynomial.b),
        int(polynomial.square_coefficient),
    ]
