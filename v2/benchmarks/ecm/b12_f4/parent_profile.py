"""Bounded parent/metadata diagnostics with explicit observer overhead."""

import argparse
import json
import time
from pathlib import Path


class RecorderLimitError(RuntimeError):
    """Observer exhaustion, distinct from factorer resource refusal."""


class CoarseCounter:
    """Inclusive/exclusive clocks with a finite number of observed calls."""

    def __init__(self, maximum=200000, *, cpu_clock=None, wall_clock=None):
        self.maximum = maximum
        self.cpu_clock = cpu_clock or time.process_time
        self.wall_clock = wall_clock or time.monotonic
        self.enabled = False
        self.events = 0
        self.stack = []
        self.rows = {}

    def wrap(self, name, function):
        def observed(*args, **kwargs):
            if not self.enabled:
                return function(*args, **kwargs)
            if self.events >= self.maximum:
                raise RecorderLimitError("recorder_call_limit")
            self.events += 1
            cpu, wall = self.cpu_clock(), self.wall_clock()
            children = [0.0, 0.0]
            self.stack.append(children)
            try:
                return function(*args, **kwargs)
            finally:
                elapsed_cpu = self.cpu_clock() - cpu
                elapsed_wall = self.wall_clock() - wall
                self.stack.pop()
                if self.stack:
                    self.stack[-1][0] += elapsed_cpu
                    self.stack[-1][1] += elapsed_wall
                row = self.rows.setdefault(name, [0, 0.0, 0.0, 0.0, 0.0])
                row[0] += 1
                row[1] += elapsed_cpu
                row[2] += max(0.0, elapsed_cpu - children[0])
                row[3] += elapsed_wall
                row[4] += max(0.0, elapsed_wall - children[1])

        return observed


def child(path):
    """Patch this diagnostic interpreter; leave worker execution unchanged."""
    from ....common import utils
    from ....execution import ecm_workers, stage_jobs, worker_contract
    from . import profile_runner

    request = json.loads(path.read_text())
    marker = Path(request["capture"]) / "measuring"
    counter = CoarseCounter()
    patches = []

    def install(owner, name, wrapped):
        patches.append((owner, name, getattr(owner, name)))
        setattr(owner, name, wrapped)

    for owner, names in (
        (
            worker_contract.Search,
            ("next", "dispatchable_count", "receive", "snapshot"),
        ),
        (worker_contract.Manifest, ("assignment",)),
        (
            ecm_workers.ECMExecutor,
            ("_manifest", "_packet", "_send", "_readable"),
        ),
        (stage_jobs, ("prod",)),
        (utils, ("prime_power",)),
    ):
        for name in names:
            label = getattr(owner, "__name__", "unknown") + "." + name
            install(owner, name, counter.wrap(label, getattr(owner, name)))
    for name in ("canonical", "digest", "frame"):
        wrapped = counter.wrap(name, getattr(worker_contract, name))
        install(worker_contract, name, wrapped)
        install(ecm_workers, name, wrapped)
    original = ecm_workers.ECMExecutor.run

    def measured_run(*args, **kwargs):
        counter.enabled = marker.exists()
        try:
            return original(*args, **kwargs)
        finally:
            counter.enabled = False

    install(ecm_workers.ECMExecutor, "run", measured_run)
    try:
        profile_runner.child(path)
    finally:
        for owner, name, previous in reversed(patches):
            setattr(owner, name, previous)
    output = dict(
        instrumented=True,
        events=counter.events,
        columns=[
            "calls",
            "inclusive_cpu",
            "exclusive_cpu",
            "inclusive_wall",
            "exclusive_wall",
        ],
        rows=counter.rows,
        limitation="Wrapping perturbs JIT/clock costs; compare the clean arm.",
    )
    with (Path(request["capture"]) / "parent-profile.json").open(
        "x"
    ) as stream:
        json.dump(output, stream)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--request", type=Path, required=True)
    child(parser.parse_args().request)


if __name__ == "__main__":
    main()
