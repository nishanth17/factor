"""Finite, source-pinned diagnosis of the retained 65-digit slowdown."""

import argparse
import copy
import gc
import hashlib
import json
import os
import resource
import subprocess
import sys
import time
from pathlib import Path

from ..c10a04.accounting import digest, output_lock
from . import execute_v4, nominee_v4, pilots_v4, plan_v4

ROOT = plan_v4.ROOT
PLAN = Path(__file__).with_name("regression65_v1.json")
OUTPUT = ROOT / "v2/benchmarks/results/c10a01/regression65-v1"
MODULE = "v2.benchmarks.qs.c10a01.regression65"
WARM_WORK = 290_000_000_000
Ledger = pilots_v4.BoundedLedger


def build_plan():
    """Freeze caps before measurement, using the existing analytic charges."""
    frozen = json.loads(nominee_v4.FREEZE.read_text())
    configurations = {
        row["arm"]: row["configuration"]
        for row in frozen["calls"]
        if row["digits"] == 65
    }
    rows = []
    # Reverse the observer order between arms. These are identical-work
    # overhead checks, not independent-input evidence for adoption.
    for arm, observations in (
        ("control", (False, True)),
        ("nominee", (True, False)),
    ):
        for observed in observations:
            rows.append(
                dict(
                    id=f"screen-{arm}-{'on' if observed else 'off'}",
                    stage="screen",
                    arm=arm,
                    observed=observed,
                    polynomials=1024,
                    chunk_polynomials=128 if observed else 1024,
                    seconds=90,
                    configuration=configurations[arm],
                )
            )
    for arm in ("nominee", "control"):
        rows.append(
            dict(
                id=f"long-{arm}",
                stage="long",
                arm=arm,
                observed=True,
                polynomials=8192,
                chunk_polynomials=1024,
                seconds=600,
                configuration=configurations[arm],
            )
        )
    for row in rows:
        row.update(digits=65, fixture="c10a01_65_balanced", tier=2)
        row["work"] = plan_v4.call_work(65, row["seconds"])
    files = subprocess.check_output(
        [
            "git",
            "ls-files",
            "v2/qs",
            "v2/execution",
            "v2/common",
            "v2/linear_algebra",
            "v2/benchmarks/qs/c10a01",
            "v2/benchmarks/qs/c10a04/accounting.py",
        ],
        cwd=ROOT,
        text=True,
    ).splitlines()
    files.extend(
        [
            str(Path(__file__).relative_to(ROOT)),
            str(plan_v4.CORPUS.relative_to(ROOT)),
            str(nominee_v4.FREEZE.relative_to(ROOT)),
            str(pilots_v4.control_v4.BASELINE.relative_to(ROOT)),
        ]
    )
    sources = {
        name: digest(ROOT / name)
        for name in sorted(set(files))
        if name.endswith(".py")
        or name.endswith(".json")
        and name in files[-3:]
    }
    phases = {
        stage: dict(
            wall=sum(r["seconds"] + 5 for r in rows if r["stage"] == stage),
            cpu=sum(r["seconds"] + 5 for r in rows if r["stage"] == stage),
            work=sum(r["work"] for r in rows if r["stage"] == stage),
        )
        for stage in ("screen", "long")
    }
    phases["warmup"] = dict(wall=420, cpu=420, work=12 * WARM_WORK)
    phases["analysis"] = dict(wall=120, cpu=120, work=1_000_000_000_000)
    phases["orchestration"] = dict(wall=60, cpu=60, work=0)
    return dict(
        identity="c10a01-regression65-v1",
        source_commit=subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        sources=sources,
        runtime=digest(Path(sys.executable).resolve()),
        rows=rows,
        limits=dict(
            **phases,
            total={
                key: sum(p[key] for p in phases.values())
                for key in ("wall", "cpu", "work")
            },
            disk_bytes=16 * 2**20,
        ),
        work_derivation="plan_v4.call_work: 4 * 8e9 * seconds + setup, "
        "root/assignment charges and conservative postprocessing; "
        "warmups retain the original finite 290e9 work per call.",
        decision="First validate matching work and observer overhead. "
        "Long calls require separate launch after review. No retries.",
    )


def check():
    plan = json.loads(PLAN.read_text())
    if sys.implementation.name != "pypy" or sys.version_info[:2] != (3, 11):
        raise RuntimeError("PyPy implementing Python 3.11 required")
    if digest(Path(sys.executable).resolve()) != plan["runtime"]:
        raise RuntimeError("runtime drift")
    for name, expected in plan["sources"].items():
        if digest(ROOT / name) != expected:
            raise RuntimeError("source drift: " + name)
    for key in ("wall", "cpu", "work"):
        assert (
            sum(
                p[key]
                for name, p in plan["limits"].items()
                if name not in ("total", "disk_bytes")
            )
            == plan["limits"]["total"][key]
        )
    return plan


def snapshot():
    """Read GC totals; this build's JIT snapshot API aborts outright."""
    import pypyjit

    stats = gc.get_stats(memory_pressure=False)
    return dict(
        gc_time_native=stats.total_gc_time,
        gc_used=stats.total_gc_memory,
        gc_allocated=stats.total_allocated_memory,
        gc_peak=stats.peak_memory,
        jit_bytes=pypyjit.get_stats_asmmemmgr(),
        rss_high_water_bytes=resource.getrusage(
            resource.RUSAGE_SELF
        ).ru_maxrss,
    )


def state_digest(job):
    """Bind semantic atoms independently of clocks and diagnostic snapshots."""
    state = hashlib.sha256()
    state.update(
        str(
            (job.multiplier, job.family_index, job.stats["polynomials"])
        ).encode()
    )
    if job.engine is not None:
        for key in sorted(job.engine.collector._atoms):
            state.update(str(key).encode())
    return state.hexdigest()


def prefix(row, fixture, modules):
    job, budget, _ = execute_v4.instantiate(
        row, fixture, modules, seconds=row["seconds"], work=row["work"]
    )
    started, cpu = time.monotonic(), time.process_time()
    chunks, previous_work, previous_positions = [], 0, 0
    total_blocks = row["polynomials"] * 33
    for offset in range(0, total_blocks, row["chunk_polynomials"] * 33):
        before = snapshot() if row["observed"] else None
        begin_wall, begin_cpu = time.monotonic(), time.process_time()
        result = job.run(
            max_blocks=min(
                row["chunk_polynomials"] * 33, total_blocks - offset
            )
        )
        positions = (
            0
            if job.engine is None
            else job.engine.stats["collected_positions"]
        )
        chunks.append(
            dict(
                wall=time.monotonic() - begin_wall,
                cpu=time.process_time() - begin_cpu,
                work=budget.used - previous_work,
                positions=positions - previous_positions,
                families=job.family_index,
                polynomials=job.stats["polynomials"],
                before=before,
                after=snapshot() if row["observed"] else None,
                stages=copy.deepcopy(result.stats["stage_seconds"]),
                engine_stages=copy.deepcopy(
                    result.stats.get("engine", {}).get("stage_seconds", {})
                ),
            )
        )
        previous_work, previous_positions = budget.used, positions
        if result.reason != "paused":
            break
    validated = execute_v4.validate_result(
        fixture, result, budget, modules[2], modules[1].BudgetExhaustedError
    )
    if (
        result.reason == "paused"
        and previous_positions != row["polynomials"] * 131073
    ):
        raise RuntimeError("fixed-position target missed")
    return dict(
        **validated,
        reason=result.reason,
        work=budget.used,
        wall=time.monotonic() - started,
        cpu=time.process_time() - cpu,
        positions=previous_positions,
        digest=state_digest(job),
        chunks=chunks,
        stats=copy.deepcopy(result.stats),
        memory=snapshot(),
    )


def worker(row, plan):
    execute_v4.OUTPUT = OUTPUT
    execute_v4.plan_v4.RSS = {
        key: value - pilots_v4.PARENT_RSS for key, value in plan_v4.RSS.items()
    }
    modules = execute_v4.runtime_modules(row)
    fixture = execute_v4.identity(row, json.loads(plan_v4.CORPUS.read_text()))
    ledger = Ledger(OUTPUT, digest(PLAN), plan["limits"])
    identity = dict(cell=row["id"], protocol=digest(PLAN))
    for phase in ("before", "measured", "after"):
        ledger.call(
            row["stage"] if phase == "measured" else "warmup",
            row["id"] + "-" + phase,
            row["seconds"] + 5 if phase == "measured" else 35,
            row["work"] if phase == "measured" else WARM_WORK,
            identity,
            lambda phase=phase: prefix(row, fixture, modules)
            if phase == "measured"
            else pilots_v4.sentinel(row, fixture, modules),
        )


def run(stage, plan):
    OUTPUT.mkdir(parents=True, exist_ok=True)
    with output_lock(OUTPUT):
        ledger = Ledger(OUTPUT, digest(PLAN), plan["limits"])
        for row in plan["rows"]:
            if row["stage"] != stage:
                continue
            if any(
                c["identity"].get("cell") == row["id"]
                for c in ledger.data["calls"]
            ):
                raise RuntimeError("cell already started; preserve its result")
            pilots_v4.host_admission(2)
            env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
            if row["observed"]:
                env["PYPYLOG"] = "jit-summary:" + str(
                    OUTPUT / (row["id"] + "-jit.log")
                )
            child = subprocess.Popen(
                [sys.executable, "-m", MODULE, "worker", "--cell", row["id"]],
                cwd=ROOT,
                env=env,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            from .execute_v4_stage import read_worker

            try:
                answer = read_worker(
                    child, time.monotonic() + row["seconds"] + 80
                )
            finally:
                if child.poll() is None:
                    child.kill()
                    child.wait()
            print(json.dumps(answer), flush=True)
            ledger = Ledger(OUTPUT, digest(PLAN), plan["limits"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "action", choices=("freeze", "check", "screen", "long", "worker")
    )
    parser.add_argument("--cell")
    args = parser.parse_args()
    if args.action == "freeze":
        if PLAN.exists():
            raise RuntimeError("freeze already exists")
        PLAN.write_text(
            json.dumps(build_plan(), indent=2, sort_keys=True) + "\n"
        )
        return
    plan = check()
    if args.action == "check":
        print(json.dumps(dict(ready=True, limits=plan["limits"])))
    elif args.action == "worker":
        row = next(r for r in plan["rows"] if r["id"] == args.cell)
        worker(row, plan)
        print(json.dumps(dict(cell=row["id"], complete=True)))
    else:
        run(args.action, plan)


if __name__ == "__main__":
    main()
