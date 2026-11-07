"""Bounded outer-stage counters, installed only by the experiment child."""

import time
from collections import defaultdict
from functools import wraps


class TierA:
    """Time outer boundaries without wrapping arithmetic or atom loops."""

    def __init__(self):
        self.patches = []
        self.reset()

    def reset(self):
        self.costs = defaultdict(lambda: [0, 0.0, 0.0])
        self.jobs = {}
        self.transport_bytes = 0
        self.started_cpu = time.process_time()

    def patch(self, cls, name, category):
        original = getattr(cls, name)

        @wraps(original)
        def timed(owner, *args, **kwargs):
            key = category
            if name == "_timed":
                key += "." + args[0]
            wall, cpu = time.perf_counter(), time.process_time()
            try:
                return original(owner, *args, **kwargs)
            finally:
                row = self.costs[key]
                row[0] += 1
                row[1] += time.perf_counter() - wall
                row[2] += time.process_time() - cpu

        self.patches.append((cls, name, original))
        setattr(cls, name, timed)

    def install(self):
        from ....qs.current_parallel import CurrentParallelSIQSJob
        from ....qs.pipeline import QSJob
        from ....qs.siqs import SIQSJob

        for cls, name, category in (
            (SIQSJob, "_timed", "serial"),
            (QSJob, "_timed", "pipeline"),
            (CurrentParallelSIQSJob, "_setup", "parent.setup"),
            (
                CurrentParallelSIQSJob,
                "_validate_report",
                "parent.decode_check",
            ),
            (CurrentParallelSIQSJob, "_merge", "parent.admission_solve"),
            (CurrentParallelSIQSJob, "_solve", "parent.solve"),
        ):
            self.patch(cls, name, category)

        original_receive = CurrentParallelSIQSJob._receive
        original_run = CurrentParallelSIQSJob.run

        def receive(owner, report, *args, **kwargs):
            self.transport_bytes += len(report["packet"])
            return original_receive(owner, report, *args, **kwargs)

        def run(owner, *args, **kwargs):
            try:
                return original_run(owner, *args, **kwargs)
            finally:
                if len(self.jobs) >= 64 and owner.identity not in self.jobs:
                    raise RuntimeError("tier-A job observation cap exhausted")
                # Job statistics are cumulative across portfolio waves. Replace
                # each snapshot rather than counting the same work twice.
                self.jobs[owner.identity] = dict(
                    stats=dict(owner.stats), workers=owner.workers.workers
                )

        for name, wrapper in (("_receive", receive), ("run", run)):
            self.patches.append(
                (
                    CurrentParallelSIQSJob,
                    name,
                    getattr(CurrentParallelSIQSJob, name),
                )
            )
            setattr(CurrentParallelSIQSJob, name, wrapper)

    def report(self, seconds):
        return dict(
            stages={
                key: dict(calls=r[0], wall=r[1], parent_cpu=r[2])
                for key, r in self.costs.items()
            },
            parent_cpu=time.process_time() - self.started_cpu,
            parent_busy_fraction=(time.process_time() - self.started_cpu)
            / max(seconds, 1e-9),
            relation_payload_bytes=self.transport_bytes,
            jobs=list(self.jobs.values()),
            limitations=(
                "Inclusive stage times overlap. Payload bytes exclude "
                "executor envelopes and setup messages. Idle capacity uses "
                "packet wall and includes startup/teardown/wait; "
                "it is not an OS scheduler measurement. Parent busy uses "
                "CPU/wall, not an Amdahl serial fraction."
            ),
        )

    def close(self):
        for cls, name, original in reversed(self.patches):
            setattr(cls, name, original)
        self.patches.clear()
