"""Private arithmetic adapter and spawn entry point for the B12 ECM owner."""

import os
import resource
import select
import socket
import sys
import time
from collections import OrderedDict

from .worker_contract import (
    RESIDUE_VERSION,
    Decoder,
    WorkerSettings,
    checked_residue,
    digest,
    frame,
    metadata_identity,
)


def rss_bytes():
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return value if sys.platform == "darwin" else value * 1024


def read_message(channel, decoder, *, block=False):
    """Read bounded fragments; idle workers block without periodic wakeups."""
    buffered = decoder.feed(b"")
    if buffered is not None:
        return buffered
    if not select.select([channel], [], [], None if block else 0)[0]:
        return None
    try:
        raw = channel.recv(decoder.limit + 4 - len(decoder.buffer))
    except socket.timeout:
        return None
    if not raw:
        raise EOFError("worker channel closed mid-message")
    return decoder.feed(raw)


class CurveEngine:
    """Bounded immutable metadata LRU; every curve owns fresh state."""

    def __init__(self, settings=None):
        self.settings = settings or WorkerSettings()
        self.identity = None
        self.context = self.programs = self.config = None
        self.contexts = OrderedDict()
        self.preparations = self.cache_hits = self.cache_evictions = 0
        self.decode_hits = self.decode_misses = self.decode_refusals = 0

    def _prepare(self, message, budget):
        if self.settings.worker_imports == "ecm":
            from .portfolio_config import PortfolioConfig, _execution_context
        else:
            from ..portfolio import PortfolioConfig, _execution_context
        from .ecm_worker_state import DecodedPrograms

        identity = message["schedule"]
        options = message.get("options")
        if options is not None and metadata_identity(options) != identity:
            raise ValueError("metadata configuration identity changed")
        budget.poll()
        budget.consume(1)
        if identity in self.contexts:
            self.contexts.move_to_end(identity)
            self.cache_hits += 1
            config, context, programs = self.contexts[identity]
        else:
            if options is None:
                raise ValueError("missing uncached metadata configuration")
            # Release current aliases and LRU ownership before allocating.
            self.identity = None
            self.context = self.programs = self.config = None
            if len(self.contexts) >= self.settings.metadata_slots:
                self.contexts.popitem(last=False)
                self.cache_evictions += 1
            config = PortfolioConfig(**options)
            context, programs = _execution_context({}, config, budget)
            if self.settings.packed_cache_bytes:
                if programs is None:
                    from ..ecm.programs import ECMPrograms

                    programs = ECMPrograms(
                        context, memory_bytes=self.settings.packed_cache_bytes
                    )
                else:
                    # Add funded retention without replacing accepted chains.
                    programs.memory_bytes += self.settings.packed_cache_bytes
            if self.settings.decoded_cache_bytes:
                programs = DecodedPrograms(
                    programs or context, self.settings.decoded_cache_bytes
                )
            if self.settings.stage_two_kernel == "reduced":
                from ..ecm.chains import point_add, point_double

                active = programs or context
                active.point_add, active.point_double = point_add, point_double
            self.contexts[identity] = (config, context, programs)
            self.preparations += 1
        self.identity = identity
        self.config, self.context, self.programs = config, context, programs

    def execute(self, message, *, cancelled, previous_cpu=None):
        from .budget import BudgetExhaustedError
        from .ecm_worker_state import ActionBudget
        from .stage_jobs import advance_job, new_job, prime_cursor

        assignment = message["assignment"]
        stage = message.get("stage", "curve")
        if stage not in ("curve", "stage_one", "stage_two"):
            raise ValueError("invalid worker stage")
        budget = ActionBudget(
            poll_every=self.settings.poll_every,
            work_limit=message["work"],
            seconds=max(0, message["deadline"] - time.monotonic()),
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
        reason, failure, residue = "complete", None, None
        boundary_hash = None
        tracked, before_decode = None, (0, 0, 0)
        try:
            budget.poll()
            self._prepare(message, budget)
            if hasattr(self.programs, "programs"):
                tracked = self.programs
                before_decode = (
                    tracked.hits,
                    tracked.misses,
                    tracked.refusals,
                )
            if stage == "stage_two":
                budget.consume(self.settings.residue_bytes)
                boundary = checked_residue(
                    message.get("residue"),
                    assignment,
                    message["schedule"],
                    self.settings.residue_bytes,
                )
                boundary_hash = boundary["sha256"]
                job.update(
                    phase="stage_two_setup",
                    value=boundary["value"],
                    a24=boundary["a24"],
                    sigma=boundary["sigma"],
                    chain_chunks=boundary["chain_chunks"],
                    cursor=prime_cursor(job["b1"] + 1, job["b2"] + 1),
                )
            while not job["done"]:
                budget.enter_phase(job["phase"])
                if stage == "stage_one" and job["phase"] == "stage_two_setup":
                    budget.consume(self.settings.residue_bytes)
                    payload = dict(
                        version=RESIDUE_VERSION,
                        id=assignment["id"],
                        metadata=message["schedule"],
                        n=job["n"],
                        b1=job["b1"],
                        b2=job["b2"],
                        sigma=job["sigma"],
                        value=list(job["value"]),
                        a24=int(job["a24"]),
                        chain_chunks=job.get("chain_chunks", 0),
                    )
                    residue = checked_residue(
                        dict(payload, sha256=digest(payload)),
                        assignment,
                        message["schedule"],
                        self.settings.residue_bytes,
                    )
                    reason = "boundary"
                    break
                advance_job(
                    job, budget, self.programs or self.context, self.config
                )
            # Real checks always precede a factor or boundary publication.
            budget.poll()
        except BudgetExhaustedError as error:
            reason, residue = str(error), None
        except (MemoryError, ValueError, ArithmeticError) as error:
            reason, residue = "error", None
            failure = type(error).__name__ + ": " + str(error)[:192]
        if stage == "stage_two" and boundary_hash is None:
            # Cancellation before arithmetic still names its queued residue.
            raw_boundary = message.get("residue")
            boundary_hash = (
                raw_boundary.get("sha256")
                if isinstance(raw_boundary, dict)
                else None
            )
        if tracked is not None:
            self.decode_hits += tracked.hits - before_decode[0]
            self.decode_misses += tracked.misses - before_decode[1]
            self.decode_refusals += tracked.refusals - before_decode[2]
        result = dict(
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
            stage=stage,
            sigma=job.get("sigma"),
            chain_chunks=job.get("chain_chunks", 0),
            preparations=self.preparations,
            prepared_schedule=self.identity,
            metadata_cache=list(self.contexts),
            cache_hits=self.cache_hits,
            cache_evictions=self.cache_evictions,
            decode_hits=self.decode_hits,
            decode_misses=self.decode_misses,
            decode_refusals=self.decode_refusals,
        )
        if boundary_hash is not None:
            result["boundary_sha256"] = boundary_hash
        if residue is not None:
            result["residue"] = residue
        if "packet" in message:
            result["packet"] = message["packet"]
        return result


def worker_main(channel, frame_bytes, startup_work, settings_data=None):
    """One bounded packet, streamed per-job results, no nested processes."""
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
        if sys.implementation.name != "pypy" or sys.version_info[:2] != (
            3,
            11,
        ):
            raise RuntimeError(
                "ECM workers require PyPy implementing Python 3.11"
            )
        settings = WorkerSettings(**(settings_data or {}))
        if settings.worker_imports == "ecm":
            from . import portfolio_config  # noqa: F401
        else:
            from .. import portfolio  # noqa: F401
        engine = CurveEngine(settings)
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
            if message["kind"] != "packet":
                raise ValueError("invalid worker command")
            jobs = message["jobs"]
            if (
                not isinstance(jobs, list)
                or not 1 <= len(jobs) <= settings.packet_curves
            ):
                raise ValueError("invalid bounded packet")
            stopped, closing = False, False

            def cancelled():
                nonlocal stopped, closing
                command = read_message(channel, decoder)
                if command is not None:
                    if command.get("kind") == "close":
                        stopped = closing = True
                    elif command.get("kind") == "cancel":
                        if command.get("packet") == message["packet"]:
                            stopped = True
                    else:
                        raise ValueError("assignment received while busy")
                return stopped

            for index, task in enumerate(jobs):
                task = dict(
                    task,
                    packet=message["packet"],
                    schedule=message["schedule"],
                    options=message.get("options") if index == 0 else None,
                )
                result = engine.execute(
                    task, cancelled=cancelled, previous_cpu=previous_cpu
                )
                current_cpu = time.process_time()
                result.update(
                    kind="result",
                    cpu=current_cpu - previous_cpu,
                    rss=rss_bytes(),
                )
                previous_cpu = current_cpu
                channel.sendall(frame(result, frame_bytes))
                if result["factor"] is not None and not message["fixed_work"]:
                    stopped = True
                if result["reason"] not in ("complete", "boundary"):
                    stopped = True
            if closing:
                break
    finally:
        channel.close()
