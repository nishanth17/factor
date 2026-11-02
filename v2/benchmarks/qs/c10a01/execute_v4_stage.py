"""Run the clean v4 primary cells with bounded receipts and durable charges."""

import argparse
import json
import os
import selectors
import subprocess
import sys
import time
from pathlib import Path

from ..c10a04.accounting import Ledger, digest
from . import execute_v4, plan_v4

RELEASE = Path(__file__).with_name("execution_v4_sources_r2.json")


def check_release(*, require_prior=False):
    plan = json.loads(plan_v4.PLAN.read_text())
    plan_v4.verify_plan(plan)
    sources = json.loads(plan_v4.SOURCES.read_text())
    for name, expected in sources["files"].items():
        if plan_v4.digest(plan_v4.ROOT / name) != expected:
            raise RuntimeError("parent source drift: " + name)
    manifest = json.loads(RELEASE.read_text())
    prior = plan_v4.ROOT / "v2/benchmarks/results/c10a01/calibration-v4"
    prior_available = True
    for name, expected in manifest["prior_receipts"].items():
        if not (prior / name).exists():
            prior_available = False
            if require_prior:
                raise RuntimeError(
                    "prior failed-release evidence missing: " + name
                )
            continue
        if digest(prior / name) != expected:
            raise RuntimeError("prior failed-release evidence drift: " + name)
    for name, expected in manifest["files"].items():
        if plan_v4.digest(plan_v4.ROOT / name) != expected:
            raise RuntimeError("execution source drift: " + name)
    if digest(plan_v4.PLAN) != manifest["plan_sha256"]:
        raise RuntimeError("plan drift")
    if sys.version != manifest["runtime_version"]:
        raise RuntimeError("runtime drift")
    if (
        digest(Path(sys.executable).resolve())
        != manifest["runtime_binary_sha256"]
    ):
        raise RuntimeError("runtime binary drift")
    return plan, manifest, prior_available


def host_admission(tier):
    ram = int(
        subprocess.check_output(["sysctl", "-n", "hw.memsize"], text=True)
    )
    pressure = subprocess.check_output(["memory_pressure", "-Q"], text=True)
    free_percent = int(
        pressure.split("System-wide memory free percentage: ")[1].split("%")[0]
    )
    available = ram * free_percent // 100
    required = plan_v4.RSS[tier] + 4 * 2**30
    if available < required:
        raise RuntimeError("host headroom refuses tier")
    return dict(
        ram=ram,
        free_percent=free_percent,
        available=available,
        required=required,
    )


def read_worker(process, deadline):
    selector = selectors.DefaultSelector()
    selector.register(process.stdout, selectors.EVENT_READ, "stdout")
    selector.register(process.stderr, selectors.EVENT_READ, "stderr")
    buffers = {"stdout": bytearray(), "stderr": bytearray()}
    while selector.get_map():
        if time.monotonic() >= deadline:
            process.kill()
            process.wait()
            raise TimeoutError("worker hard wall watchdog")
        for key, _ in selector.select(timeout=0.1):
            chunk = os.read(key.fileobj.fileno(), 4096)
            if not chunk:
                selector.unregister(key.fileobj)
                continue
            name = key.data
            buffers[name].extend(chunk)
            if len(buffers[name]) > 64 * 1024:
                process.kill()
                process.wait()
                raise RuntimeError("recorder_exhaustion: worker pipe")
    code = process.wait(timeout=1)
    if code:
        raise RuntimeError(
            "worker_exit_"
            + str(code)
            + ": "
            + buffers["stderr"].decode(errors="replace")[-1024:]
        )
    rows = buffers["stdout"].splitlines()
    if len(rows) != 1:
        raise RuntimeError("worker protocol requires one bounded result")
    return json.loads(rows[0])


def run_worker(row):
    command = [
        sys.executable,
        "-u",
        "-m",
        "v2.benchmarks.qs.c10a01.execute_v4",
        row["id"],
    ]
    started = time.monotonic()
    process = subprocess.Popen(
        command,
        cwd=plan_v4.ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        env=dict(os.environ, PYTHONDONTWRITEBYTECODE="1"),
    )
    try:
        result = read_worker(process, started + row["seconds"] + 35)
    finally:
        if process.poll() is None:
            process.kill()
            process.wait()
    before, after, measured = (
        result["before"],
        result["after"],
        result["result"],
    )
    if before["wall"] > 17 or before["cpu"] > 17:
        raise RuntimeError("prewarm overrun")
    if min(before["wall"], before["cpu"]) < 3:
        raise RuntimeError("warmup too short")
    if after["wall"] > 8 or after["cpu"] > 8:
        raise RuntimeError("postwarm overrun")
    if measured["seconds"] > row["seconds"] + 5:
        raise RuntimeError("measured wall overrun")
    if measured["cpu_seconds"] > row["seconds"] + 5:
        raise RuntimeError("measured CPU overrun")
    if measured["work"] > row["work"]:
        raise RuntimeError("measured work overrun")
    if (
        measured.get("reason") == "paused"
        and measured["positions"] != row["positions"]
    ):
        raise RuntimeError("fixed-position target missed")
    if result["total_work"] != (
        before["work"] + measured["work"] + after["work"]
    ):
        raise RuntimeError("worker work total mismatch")
    if result["total_work"] > row["work"] + warm_slot_work():
        raise RuntimeError("total work bound")
    # Ledger.call settles result['work']; preserve the separate components.
    result["work"] = result["total_work"]
    return result


def warm_slot_work():
    plan = json.loads(plan_v4.PLAN.read_text())
    return plan["calibration"]["limits"]["warm_slot_work"]


def limits(plan):
    phases = {}
    calls = plan["calibration"]["calls"]
    slot = plan["calibration"]["limits"]["warm_slot_work"]
    for phase in {row["phase"] for row in calls}:
        group = [row for row in calls if row["phase"] == phase]
        phases[phase] = dict(
            wall=sum(row["reservation_seconds"] + 25 for row in group),
            cpu=sum(row["reservation_seconds"] + 25 for row in group),
            work=sum(row["work"] + slot for row in group),
        )
    return dict(
        **phases,
        total=dict(
            wall=3400,
            cpu=3400,
            work=plan["calibration"]["limits"]["total"]["work"],
        ),
        disk_bytes=64 * 2**20,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("check", "primary"))
    args = parser.parse_args()
    plan, release, prior_available = check_release(
        require_prior=args.action == "primary"
    )
    if args.action == "check":
        print(
            json.dumps(
                dict(
                    execution_source_verified=True,
                    primary_calls=16,
                    performance_calls=0,
                )
            )
        )
        return
    if sys.implementation.name != "pypy":
        raise RuntimeError("PyPy required")
    output = execute_v4.OUTPUT
    output.mkdir(parents=True, exist_ok=True)
    ledger = Ledger(output, digest(RELEASE), limits(plan))
    for row in plan["calibration"]["calls"]:
        if row["phase"] != "primary":
            continue
        tag = row["id"]
        if any(c["receipt"] == tag + ".json" for c in ledger.data["calls"]):
            continue
        host = host_admission(row["tier"])
        identity = dict(
            cell=tag,
            plan=digest(plan_v4.PLAN),
            release=digest(RELEASE),
            fixture=row["fixture"],
            source=release["source_commit"],
            host=host,
        )
        result = ledger.call(
            "primary",
            tag,
            row["reservation_seconds"] + 25,
            row["work"] + warm_slot_work(),
            identity,
            lambda: run_worker(row),
        )
        print(
            json.dumps(
                dict(
                    cell=tag,
                    reason=result["result"]["reason"],
                    positions=result["result"]["positions"],
                    seconds=result["result"]["seconds"],
                )
            ),
            flush=True,
        )


if __name__ == "__main__":
    main()
