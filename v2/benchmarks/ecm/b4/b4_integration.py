"""Frozen production-source bridge for the selected B4 kernels."""

import argparse
import fcntl
import hashlib
import json
import os
import platform
import resource
import statistics
import subprocess
import sys
import time
from pathlib import Path

from ....common import arithmetic
from ...support.paths import (
    matches_source_pin,
    source_path,
)
from ..p52 import p52_realistic
from . import b4_bakeoff as bakeoff
from . import b4_common as common
from . import b4_kernels as kernels
from . import b4_study as study
from .build_b4_inputs import write_new

ROOT = bakeoff.ROOT
CONTROL = common.INPUTS / "baselines/b4_integrated_ecm.json"
FREEZE = common.INPUTS / "controls/b4_integration_freeze.json"


def source_data():
    """Reconstruct the immutable production package from its source delta."""
    original = json.loads(source_path(common.CONTROL).read_text())
    delta = json.loads(source_path(CONTROL).read_text())
    if (
        delta["parent_sha256"] != common.digest(common.CONTROL)
        or hashlib.sha256(delta["source"].encode()).hexdigest()
        != delta["sha256"]
    ):
        raise ValueError("integration control changed")
    sources = dict(original["source"])
    sources["v2/ecm.py"] = delta["source"]
    return sources


def verify_current():
    """Prove the measured package is exactly the current production source."""
    for path, source in source_data().items():
        if (source_path(ROOT / path)).read_text() != source:
            raise ValueError("production differs from frozen bridge: " + path)


def freeze():
    """Create a new one-shot integration freeze; preserve earlier controls."""
    bakeoff.verify_inputs()
    if CONTROL.exists() or FREEZE.exists():
        raise ValueError("integration freeze already exists")
    original = json.loads(source_path(common.CONTROL).read_text())
    changes = [
        path
        for path, source in original["source"].items()
        if (source_path(ROOT / path)).read_text() != source
    ]
    if changes != ["v2/ecm.py"]:
        raise ValueError("unexpected production changes: " + repr(changes))
    source = (source_path(ROOT / "v2/ecm.py")).read_text()
    write_new(
        CONTROL,
        dict(
            parent_sha256=common.digest(common.CONTROL),
            source=source,
            sha256=hashlib.sha256(source.encode()).hexdigest(),
        ),
    )
    paths = (Path(__file__), CONTROL, bakeoff.FREEZE)
    write_new(
        FREEZE,
        dict(
            scope="Production bridge on the already inspected confirmation "
            "cohort; no new independent-input claim or candidate tuning.",
            choices={"python-int": "reductions", "gmpy2-mpz": "whole_ladder"},
            source_hashes={
                str(path.relative_to(ROOT)): common.digest(path)
                for path in paths
            },
            protocol_sha256=common.digest(bakeoff.PROTOCOL),
            criteria="Require exact outcomes/work and the frozen revised "
            "timing gate: stable pairs, positive bridge interval/CPU "
            "effect and both chronological halves. No expansion on failure.",
        ),
    )
    verify()
    verify_current()


def verify():
    protocol, corpus, _ = bakeoff.verify_inputs()
    pins = json.loads(source_path(FREEZE).read_text())
    for path, expected in pins["source_hashes"].items():
        if not matches_source_pin(ROOT / path, expected):
            raise ValueError("integration source/input changed: " + path)
    if common.digest(bakeoff.PROTOCOL) != pins["protocol_sha256"]:
        raise ValueError("integration protocol changed")
    source_data()
    return protocol, corpus, pins


class IntegratedFinder(common.FrozenFinder):
    def exec_module(self, module):
        filename = self.path(module.__name__)
        module.__file__ = str(CONTROL) + ":" + filename
        exec(
            compile(self.sources[filename], module.__file__, "exec"),
            module.__dict__,
        )


def engine(arm):
    if arm == "baseline":
        return kernels.engine("baseline")
    name = "_b4_integrated"
    if name not in sys.modules:
        sys.meta_path.insert(0, IntegratedFinder(name, source_data()))
    return __import__(name + ".portfolio", fromlist=["portfolio"])


def worker(args):
    protocol, corpus, _ = verify()
    selected = engine(args.arm)
    fixtures = [f for f in corpus["fixtures"] if f["split"] == "confirmation"]
    expected = None

    def measure():
        nonlocal expected
        started, cpu = time.perf_counter(), time.process_time()
        rows = study.measure_full(selected, protocol, args.backend, fixtures)
        elapsed, used = (
            time.perf_counter() - started,
            time.process_time() - cpu,
        )
        current = study.signature(rows)
        if expected is not None and current != expected:
            raise AssertionError("integration outcome/work changed")
        expected = current
        return rows, elapsed, used

    start, warmups = time.perf_counter(), 0
    while time.perf_counter() - start < args.warmup:
        measure()
        warmups += 1
    warm_elapsed = time.perf_counter() - start
    samples, elapsed, cpu = [], 0.0, 0.0
    while not samples or cpu < protocol["sample_cpu_seconds"]:
        rows, wall_used, cpu_used = measure()
        samples.append(rows)
        elapsed += wall_used
        cpu += cpu_used
    averaged = []
    for index, row in enumerate(samples[0]):
        row = dict(row)
        row["seconds"] = statistics.mean(
            sample[index]["seconds"] for sample in samples
        )
        averaged.append(row)
    return dict(
        arm=args.arm,
        backend=args.backend,
        block=args.block,
        warmup_seconds=warm_elapsed,
        warmup_cohorts=warmups,
        measured_cohorts=len(samples),
        seconds=elapsed / len(samples),
        cpu_seconds=cpu / len(samples),
        rows=averaged,
        rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    )


def launch(arm, backend, block, warmup, deadline):
    p52_realistic.check_quiet({os.getpid()})
    command = [
        sys.executable,
        "-B",
        "-m",
        "v2.benchmarks.ecm.b4.b4_integration",
        "--worker",
        "--arm",
        arm,
        "--backend",
        backend,
        "--block",
        str(block),
        "--warmup",
        str(warmup),
    ]
    process = subprocess.Popen(
        command, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    started = time.perf_counter()
    try:
        while True:
            try:
                stdout, stderr = process.communicate(timeout=1)
                break
            except subprocess.TimeoutExpired:
                p52_realistic.check_quiet({os.getpid(), process.pid})
                if time.perf_counter() > min(deadline, started + 600):
                    raise TimeoutError("integration child/window expired")
        if process.returncode:
            raise RuntimeError(stderr)
        p52_realistic.check_quiet({os.getpid()})
        return json.loads(stdout)
    finally:
        if process.poll() is None:
            process.terminate()
        process.communicate()


def run(args):
    protocol, _, pins = verify()
    if args.current:
        verify_current()
    report = dict(
        freeze=pins,
        runtime=sys.version,
        executable=sys.executable,
        platform=platform.platform(),
        backend_identities={
            name: arithmetic.get_backend(name).identity
            for name in protocol["backends"]
        },
        verified_current_source=args.current,
        source_revision=subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        captures=[],
        summary=[],
        rounds=[],
    )
    p52_realistic.check_quiet({os.getpid()})
    with open(protocol["machine_lock"], "a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        owner = Path("/private/tmp/factor-performance-owner.json")
        owner.write_text(
            json.dumps(
                dict(
                    owner="B4 production integration",
                    pid=os.getpid(),
                    worktree=str(ROOT),
                    phase="bridge",
                )
            )
        )
        deadline = time.perf_counter() + protocol["window_seconds"]
        active = set(protocol["backends"])
        completed = 0
        try:
            with args.output.open("x") as output:
                bakeoff.save(output, report)
                for warmup, target in protocol["sampling"]:
                    for block in range(completed, target):
                        backends = list(protocol["backends"])
                        if block % 2:
                            backends.reverse()
                        for backend in backends:
                            if backend not in active:
                                continue
                            order = bakeoff.block_order(
                                ("baseline", "production"),
                                backend,
                                block,
                                protocol["order_seed"],
                            )
                            for arm in order:
                                if time.perf_counter() > deadline:
                                    raise TimeoutError("bridge window expired")
                                report["captures"].append(
                                    launch(
                                        arm, backend, block, warmup, deadline
                                    )
                                )
                                bakeoff.save(output, report)
                        print("integration block", block + 1, flush=True)
                    completed = target
                    report["summary"] = [
                        bakeoff.compare(
                            report["captures"], backend, "production", protocol
                        )
                        for backend in protocol["backends"]
                    ]
                    if not all(
                        row["matched_outcomes_and_work"]
                        for row in report["summary"]
                    ):
                        raise AssertionError("production bridge changed work")
                    report["rounds"].append(
                        dict(blocks=target, summary=report["summary"])
                    )
                    bakeoff.save(output, report)
                    active = {
                        row["backend"]
                        for row in report["summary"]
                        if not row["stable"]
                    }
                    if not active:
                        break
                report["accepted"] = all(
                    row["timing_evidence_passes"] for row in report["summary"]
                )
                verify()
                if args.current:
                    verify_current()
                bakeoff.save(output, report)
        finally:
            owner.unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--freeze", action="store_true")
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--current", action="store_true")
    parser.add_argument("--arm", choices=("baseline", "production"))
    parser.add_argument("--backend", choices=("python-int", "gmpy2-mpz"))
    parser.add_argument("--block", type=int, default=0)
    parser.add_argument("--warmup", type=float, default=3)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    common.require_runtime()
    if args.freeze:
        freeze()
    elif args.worker:
        print(json.dumps(worker(args)))
    else:
        if args.output is None:
            parser.error("output is required")
        run(args)


if __name__ == "__main__":
    main()
