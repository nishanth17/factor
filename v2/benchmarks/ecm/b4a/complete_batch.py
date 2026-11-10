"""Complete fixed-work ECM calls with scalar, team or pool ownership."""

import ctypes
import errno
import hashlib
import json
import os
import time

from .ownership import GroupObserver, Usage, cpu_seconds


class BatchSession:
    """Retain only explicitly keyed owners; never reset warm state silently."""

    def __init__(self):
        from ....execution.resources import MacResources

        self.owners = {}
        self.resources = MacResources()
        self.observer = GroupObserver()

    def cpu(self):
        """Count subtree CPU, including reaped children and native threads."""
        records = []
        for sample in self.resources.sample({os.getpid(): None}):
            usage = Usage()
            if self.observer.lib.proc_pid_rusage(
                sample.pid, 2, ctypes.byref(usage)
            ):
                if ctypes.get_errno() == errno.ESRCH:
                    continue
                raise OSError("batch CPU observation failed")
            records.append(usage)

        return cpu_seconds(records, self.observer.numer, self.observer.denom)

    @staticmethod
    def make_owner(request):
        from ....execution.ecm_native import NativeECMExecutor
        from ....execution.ecm_workers import ECMExecutor
        from ....execution.worker_contract import WorkerSettings

        role = request["role"]
        workers = request.get("workers", 1)
        settings = WorkerSettings(
            workers=workers,
            allow_wide_workers=workers not in (1, 2, 4),
        )
        if role == "team":
            return NativeECMExecutor(
                seed=request.get("seed", 73),
                settings=settings,
                schedule=request.get("schedule", "dynamic"),
                capacity=request.get("capacity"),
                chunk_size=request.get("queue_chunk", 1),
            )
        if role not in ("python", "scalar", "process"):
            raise ValueError("unknown complete-batch role")
        if role != "process" and workers != 1:
            raise ValueError("serial batch requires one worker")
        return ECMExecutor(
            seed=request.get("seed", 73),
            settings=settings,
            force_workers=role == "process",
        )

    @staticmethod
    def state(owner):
        team = getattr(owner, "team", None)
        if team is not None:
            packets = list(team.packets)
            store = team.blocks
            lifecycle = (
                list(team.lifecycle)
                if team.owner is None
                else [
                    int(team.lib.b4_pthread_counter(team.owner, index))
                    for index in range(4)
                ]
            )
        else:
            engine = getattr(owner, "serial_engine", None)
            source = (
                None if engine is None else engine.programs or engine.context
            )
            packet = getattr(source, "native_packet", None)
            packets = [] if packet is None else [packet]
            store = getattr(packet, "blocks", None)
            lifecycle = None
        return dict(
            contexts=len(packets),
            live_contexts=sum(p.context is not None for p in packets),
            returns=sum(getattr(p, "calls", 0) for p in packets),
            slices=getattr(team, "slices", None),
            blocks=len(getattr(store, "blocks", {})),
            retained_bytes=getattr(store, "used", 0),
            thread_lifecycle=lifecycle,
        )

    def call(self, request):
        from ....common import utils
        from ....execution.budget import Budget
        from ....execution.portfolio_config import PortfolioConfig

        role = request["role"]
        key = (
            role,
            request.get("workers", 1),
            request.get("schedule", "dynamic"),
            request.get("queue_chunk", 1),
            request.get("capacity"),
            request["b1"],
            request["b2"],
            request.get("seed", 73),
        )
        if key not in self.owners:
            self.owners[key] = self.make_owner(request)
        owner = self.owners[key]
        before_state = self.state(owner)
        before_cpu = self.cpu()
        started = time.monotonic()
        budget = Budget(
            seconds=request.get("seconds", 30),
            cpu_seconds=request.get("cpu_seconds", 120),
            rss_bytes=8 * 2**30,
        )
        try:
            config = PortfolioConfig(
                trial_bound=2,
                rho_attempts=0,
                pm1_attempts=0,
                siqs=None,
                ecm_tiers=((request["b1"], request["b2"], request["curves"]),),
                ecm_kernel="python" if role == "python" else "c-gmp-full",
                ecm_chain_mode="auto",
                max_input_bits=1024,
            )
            result = owner.run(
                request["n"], config=config, budget=budget, fixed_work=True
            )
            payload = (
                result.checkpoint["payload"]
                if role == "team"
                else result.checkpoint["search"]["payload"]
            )
            outcomes = [
                [key, value["id"], value["sigma"], value["factor"]]
                for key, value in sorted(payload["completed"].items())
            ]
            for _, _, _, factor in outcomes:
                if factor is not None and not utils.valid_divisor(
                    factor, request["n"]
                ):
                    raise ValueError("batch returned an improper divisor")
                if factor is not None:
                    assert factor * (request["n"] // factor) == request["n"]
            exhausted = result.reason == "curves_exhausted"
            if exhausted and len(outcomes) != request["curves"]:
                raise ValueError("complete batch omitted a numerical curve")
            if role == "process" and exhausted and not (
                result.stats.get("engaged")
                and result.stats.get("worker_dispatches", 0)
            ):
                raise ValueError("process control did not dispatch workers")
            state = self.state(owner)
            row = dict(
                request,
                reason=result.reason,
                completed=len(outcomes),
                exhausted=exhausted,
                outcomes=outcomes,
                signature=hashlib.sha256(
                    json.dumps(outcomes).encode()
                ).hexdigest(),
                work=budget.used,
                observed_rss=budget.observed_rss,
                state=state,
                native_returns=state["returns"] - before_state["returns"],
                stats=result.stats,
            )
        finally:
            budget.close()
        row.update(
            wall=time.monotonic() - started, cpu=self.cpu() - before_cpu
        )
        return row

    def close(self):
        for owner in self.owners.values():
            owner.close()
        self.owners.clear()


def serve(stream, output, *, once=False):
    """Process JSON requests; cold callers time the entire child process."""
    session = BatchSession()
    try:
        for line in stream:
            request = json.loads(line)
            if request.get("op") == "warm":
                useful, calls = 0.0, 0
                while useful < 3:
                    row = session.call(request)
                    if not row["exhausted"]:
                        raise ValueError("warmup did not finish fixed work")
                    useful += row["wall"]
                    calls += 1
                row = dict(request, useful=useful, calls=calls)
            else:
                row = session.call(request)
            if once:
                session.close()
            output.write(json.dumps(row) + "\n")
            output.flush()
            if once:
                break
    finally:
        session.close()
