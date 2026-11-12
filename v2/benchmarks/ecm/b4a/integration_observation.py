"""Runner-side phase, canonical-work and owned-thread handoff observations.

Hooks preserve the production calls, numerical arguments and return values.
Only the benchmark runner installs them; clean ECM comparisons do not import
this module. Per-atom prefix hashes are bounded evidence, never routing data.
"""

import hashlib
import json
import time
from contextlib import ExitStack
from dataclasses import asdict
from unittest.mock import patch

from .thread_observation import NativeThreadProbe


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


class IntegrationObserver:
    """Observe one portfolio lifetime with unchanged production calls."""

    def __init__(self, owner, prefix_path, *, prefix_bytes=8 * 2**20):
        self.owner = owner
        self.prefix_path = prefix_path
        self.prefix_limit = prefix_bytes
        self.stream = None
        self.stack = ExitStack()
        self.phases = []
        self.snapshots = []
        self.jobs = {}
        self.assignments = {}
        self.current_assignment = None
        self.probe = None
        self.qs_entry = None
        self.qs_early = None
        self.observer_wall = 0.0
        self.qs_calls = 0
        self.merge_calls = 0

    @staticmethod
    def cpu_tree():
        import os

        from ....execution.resources import MacResources

        rows = MacResources().sample({os.getpid(): None})
        return dict(
            cpu=sum(row.cpu + row.reaped_cpu for row in rows),
            rss=sum(row.rss for row in rows),
            processes=[
                dict(
                    pid=row.pid,
                    birth=row.birth,
                    cpu=row.cpu,
                    reaped_cpu=row.reaped_cpu,
                    rss=row.rss,
                )
                for row in rows
            ],
        )

    def owner_state(self):
        team = getattr(self.owner, "team", None)
        if team is None:
            return dict(
                native=False,
                closed=getattr(self.owner, "closed", None),
                ecm_worker_pids=[
                    slot["process"].pid
                    for slot in getattr(self.owner, "slots", ())
                ],
            )
        return dict(
            native=True,
            closed=team.closed,
            contexts=len(team.packets),
            live_contexts=sum(p.context is not None for p in team.packets),
            shared_native_blocks=len(team.blocks.blocks),
            metadata_reserve_estimate=team.blocks.used,
            activity=team.worker_activity(),
            lifecycle=team.lifecycle,
            ownership="explicit caller retains ECM owner until final close",
        )

    def snapshot(self, label):
        started = time.monotonic()
        state = self.owner_state()
        if state["native"] and not state["closed"]:
            if self.probe is None:
                self.probe = NativeThreadProbe(self.owner.team)
            threads = self.probe.snapshot(label)
        else:
            threads = None
        self.snapshots.append(
            dict(
                label=label,
                at=time.monotonic(),
                tree=self.cpu_tree(),
                native=state,
                native_threads=threads,
            )
        )
        self.observer_wall += time.monotonic() - started

    def _timed(self, phase, call, *args, **kwargs):
        before = self.cpu_tree()
        started = time.monotonic()
        try:
            return call(*args, **kwargs)
        finally:
            after = self.cpu_tree()
            self.phases.append(
                dict(
                    phase=phase,
                    wall=time.monotonic() - started,
                    cpu=after["cpu"] - before["cpu"],
                    start_rss=before["rss"],
                    end_rss=after["rss"],
                )
            )

    def _prefix(self, job, index, rows):
        started = time.monotonic()
        key = f"{job.identity}/{index}"
        record = self.assignments.get(key)
        if record is None:
            record = dict(
                job=job.identity,
                index=index,
                offset=self.stream.tell(),
                rows=0,
                digest=hashlib.sha256(key.encode()).hexdigest(),
            )
            self.assignments[key] = record
        if self.current_assignment not in (None, key):
            if record["rows"]:
                raise ValueError("QS canonical assignment order regressed")
        self.current_assignment = key
        if self.stream.tell() + len(rows) * 32 > self.prefix_limit:
            raise RuntimeError("bounded canonical-QS recorder exhausted")
        digest = bytes.fromhex(record["digest"])
        for row in rows:
            digest = hashlib.sha256(digest + canonical(row).encode()).digest()
            self.stream.write(digest)
        record["rows"] += len(rows)
        record["digest"] = digest.hex()
        self.observer_wall += time.monotonic() - started

    def __enter__(self):
        from .... import portfolio
        from ....ecm.workers import ECMExecutor
        from ....qs.current_parallel import CurrentParallelSIQSJob
        from ....qs.pipeline import QSJob

        self.stream = self.stack.enter_context(self.prefix_path.open("xb"))
        handoff = portfolio._handoff
        setup = CurrentParallelSIQSJob._setup
        merge = CurrentParallelSIQSJob._merge
        run = CurrentParallelSIQSJob.run
        b12_run = ECMExecutor.run
        native_run = getattr(type(self.owner), "run", None)
        qs_timed = QSJob._timed

        def observed_ecm(original):
            def execute(owner, *args, **kwargs):
                if owner is self.owner:
                    return self._timed(
                        "ecm_search_and_prepare",
                        original,
                        owner,
                        *args,
                        **kwargs,
                    )
                return original(owner, *args, **kwargs)

            return execute

        def observed_handoff(state, config, budget, reason):
            previous = state["current"]["stage"]
            result = self._timed(
                "handoff", handoff, state, config, budget, reason
            )
            current = state["current"]
            if previous == "ecm" and current and current["stage"] == "siqs":
                self.qs_entry = time.monotonic()
                self.snapshot("qs-entry")
            return result

        def observed_setup(job):
            first = job.engine is None
            result = (
                self._timed("qs_prepare", setup, job) if first else setup(job)
            )
            if first:
                self.jobs[job.identity] = dict(
                    mathematical=job.identity,
                    execution=job.execution_identity,
                    seed=job.seed,
                    config=asdict(job.config),
                    workers=job.workers.record(),
                )
            return result

        def observed_merge(job, index=None):
            index = job.next_assignment if index is None else index
            report = job.pending[index]
            before = report["admitted"]
            try:
                return merge(job, index)
            finally:
                self._prefix(
                    job, index, report["rows"][before : report["admitted"]]
                )
                self.merge_calls += 1
                if self.merge_calls == 1:
                    self.snapshot("qs-early-published-work")
                    self.qs_early = time.monotonic()
                elif (
                    self.merge_calls > 1
                    and self.qs_early is not None
                    and time.monotonic() - self.qs_early >= 3
                    and not any(
                        s["label"] == "qs-steady" for s in self.snapshots
                    )
                ):
                    self.snapshot("qs-steady")

        def observed_run(job, *args, **kwargs):
            result = self._timed("qs_run_inclusive", run, job, *args, **kwargs)
            self.qs_calls += 1
            return result

        def observed_qs_timed(job, stage, call):
            if stage in (
                "preparation",
                "filtering",
                "elimination",
                "extraction",
            ):
                return self._timed("qs_" + stage, qs_timed, job, stage, call)
            return qs_timed(job, stage, call)

        self.stack.enter_context(
            patch.object(portfolio, "_handoff", observed_handoff)
        )
        self.stack.enter_context(
            patch.object(ECMExecutor, "run", observed_ecm(b12_run))
        )
        if native_run is not None and type(self.owner) is not ECMExecutor:
            self.stack.enter_context(
                patch.object(type(self.owner), "run", observed_ecm(native_run))
            )
        self.stack.enter_context(
            patch.object(CurrentParallelSIQSJob, "_setup", observed_setup)
        )
        self.stack.enter_context(
            patch.object(CurrentParallelSIQSJob, "_merge", observed_merge)
        )
        self.stack.enter_context(
            patch.object(CurrentParallelSIQSJob, "run", observed_run)
        )
        self.stack.enter_context(
            patch.object(QSJob, "_timed", observed_qs_timed)
        )
        return self

    def finish(self):
        if self.qs_entry is not None:
            self.snapshot("qs-finished-before-ecm-close")
        self.stream.flush()
        return dict(
            phases=self.phases,
            snapshots=self.snapshots,
            jobs=self.jobs,
            assignments=list(self.assignments.values()),
            prefix_bytes=self.stream.tell(),
            observer_wall=self.observer_wall,
            merge_calls=self.merge_calls,
            qs_calls=self.qs_calls,
            scope="Inclusive phases overlap; never sum QS run with its parts. "
            "Thread CPU is explanatory, not added to H3 CPU. The recording "
            "timer covers prefix/snapshot work; matched tests measure all "
            "observer overhead including phase queries.",
        )

    def __exit__(self, *args):
        return self.stack.__exit__(*args)
