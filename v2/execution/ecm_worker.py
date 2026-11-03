"""Private arithmetic adapter and spawn entry point for the B12 ECM owner."""

import os
import resource
import select
import socket
import sys
import time

from .worker_contract import Decoder, frame


def rss_bytes():
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return value if sys.platform == "darwin" else value * 1024


def read_message(channel, decoder, *, block=False):
    """Read at most one bounded fragment; never wait on an incomplete frame."""
    buffered = decoder.feed(b"")
    if buffered is not None:
        return buffered
    if not block and not select.select([channel], [], [], 0)[0]:
        return None
    try:
        raw = channel.recv(decoder.limit + 4 - len(decoder.buffer))
    except socket.timeout:
        return None
    if not raw:
        raise EOFError("worker channel closed mid-message")
    return decoder.feed(raw)


class CurveEngine:
    """Retain one schedule's immutable metadata, never a curve's points."""

    def __init__(self):
        self.identity = None
        self.context = self.programs = self.config = None
        self.preparations = 0

    def execute(self, message, *, cancelled, previous_cpu=None):
        from ..portfolio import PortfolioConfig, _execution_context
        from .budget import Budget, BudgetExhaustedError
        from .stage_jobs import advance_job, new_job

        assignment = message["assignment"]
        deadline = message["deadline"]
        now = time.monotonic()
        budget = Budget(
            work_limit=message["work"],
            seconds=max(0, deadline - now),
            cpu_seconds=message["cpu"],
            cancelled=cancelled,
        )
        if previous_cpu is not None:
            budget._cpu_start = previous_cpu
        job = new_job(
            "ecm",
            assignment["n"],
            assignment["seed"],
            assignment["b1"],
            assignment["b2"],
        )
        reason, failure = "complete", None
        try:
            budget.consume(0)
            if self.identity != message["schedule"]:
                # Drop all old ownership before constructing a replacement.
                self.identity = None
                self.context = self.programs = self.config = None
                config = PortfolioConfig(**message["options"])
                context, programs = _execution_context({}, config, budget)
                self.config, self.context, self.programs = (
                    config,
                    context,
                    programs,
                )
                self.identity = message["schedule"]
                self.preparations += 1
            while not job["done"]:
                budget.consume(0)
                advance_job(
                    job, budget, self.programs or self.context, self.config
                )
            # Check stop precedence before returning a valid split.
            budget.consume(0)
        except BudgetExhaustedError as error:
            reason = str(error)
        except (MemoryError, ValueError, ArithmeticError) as error:
            reason = "error"
            failure = type(error).__name__ + ": " + str(error)[:192]

        return dict(
            id=assignment["id"],
            token=message["token"],
            attempt=message["attempt"],
            work=budget.used,
            factor=int(job["factor"])
            if reason == "complete" and job["factor"] is not None
            else None,
            reason=reason,
            failure=failure,
            phase=job["phase"],
            sigma=job.get("sigma"),
            chain_chunks=job.get("chain_chunks", 0),
            preparations=self.preparations,
            prepared_schedule=self.identity,
        )


def worker_main(channel, frame_bytes, startup_work):
    """One private socket, one caller, explicit spawn; no nested processes."""
    for name in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[name] = "1"
    channel.settimeout(0.01)
    decoder = Decoder(frame_bytes)
    try:
        # Include arithmetic imports in the startup report.
        from .. import portfolio  # noqa: F401

        if sys.implementation.name != "pypy" or sys.version_info[:2] != (
            3,
            11,
        ):
            raise RuntimeError(
                "ECM workers require PyPy implementing Python 3.11"
            )
        engine = CurveEngine()
        previous_cpu = time.process_time()
        channel.sendall(
            frame(
                dict(
                    kind="ready",
                    cpu=previous_cpu,
                    rss=rss_bytes(),
                    work=startup_work,
                ),
                frame_bytes,
            )
        )
        while True:
            message = read_message(channel, decoder, block=True)
            if message is None:
                continue
            if message["kind"] == "close":
                break
            if message["kind"] == "cancel":
                continue
            if message["kind"] != "curve":
                raise ValueError("invalid worker command")
            stopped = False

            def cancelled():
                nonlocal stopped
                command = read_message(channel, decoder)
                if command is not None:
                    if command.get("kind") not in ("cancel", "close"):
                        raise ValueError("assignment received while busy")
                    stopped = True
                return stopped

            result = engine.execute(
                message, cancelled=cancelled, previous_cpu=previous_cpu
            )
            current_cpu = time.process_time()
            result.update(
                kind="result", cpu=current_cpu - previous_cpu, rss=rss_bytes()
            )
            previous_cpu = current_cpu
            channel.sendall(frame(result, frame_bytes))
    finally:
        channel.close()
