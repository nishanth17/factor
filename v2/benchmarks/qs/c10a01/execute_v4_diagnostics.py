"""Run the frozen diagnostic cells under separate receipts and phase guards."""

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

from ..c10a04.accounting import Ledger, digest
from . import execute_v4_stage, plan_v4

OUTPUT = plan_v4.ROOT / "v2/benchmarks/results/c10a01/diagnostics-v4"
PRIMARY = plan_v4.ROOT / "v2/benchmarks/results/c10a01/calibration-v4-r2"
RELEASE = Path(__file__).with_name("diagnostic_v4_sources.json")
GROUPS = {"roots", "ranking"}


def check_release(*, require_prior=False):
    plan = json.loads(plan_v4.PLAN.read_text())
    plan_v4.verify_plan(plan)
    release = json.loads(RELEASE.read_text())
    if digest(plan_v4.PLAN) != release["plan_sha256"]:
        raise RuntimeError("diagnostic plan drift")
    if (
        sys.version != release["runtime_version"]
        or digest(Path(sys.executable).resolve())
        != release["runtime_binary_sha256"]
    ):
        raise RuntimeError("diagnostic runtime drift")
    for name, expected in release["files"].items():
        if digest(plan_v4.ROOT / name) != expected:
            raise RuntimeError("diagnostic source drift: " + name)
    prior_available = True
    for name, expected in release["prior"].items():
        path = (
            PRIMARY / name
            if name == "ledger.json"
            else (plan_v4.ROOT / "v2/audit/c10a01" / name)
        )
        if not path.exists():
            prior_available = False
            if require_prior:
                raise RuntimeError("primary evidence missing: " + name)
            continue
        if digest(path) != expected:
            raise RuntimeError("primary evidence drift: " + name)
    if require_prior:
        analysis = json.loads(
            (
                plan_v4.ROOT / "v2/audit/c10a01/primary-r2-analysis.json"
            ).read_text()
        )
        for entry in analysis["rows"]:
            receipt = PRIMARY / (entry["cell"] + ".json")
            if digest(receipt) != entry["sha256"]:
                raise RuntimeError("primary cell receipt drift")
    return plan, release, prior_available


def schedule(plan):
    rows = plan["calibration"]["calls"]
    calls = []
    for digits in (62, 70):
        roots = [
            r for r in rows if r["phase"] == "roots" and r["digits"] == digits
        ]
        ranking = [
            r
            for r in rows
            if r["phase"] == "ranking" and r["digits"] == digits
        ]
        population = [
            r
            for r in rows
            if r["phase"] == "population" and r["digits"] == digits
        ]
        if digits == 70:
            population.reverse()
        calls.extend(
            (
                dict(id=f"roots-{digits}", rows=roots),
                *[dict(id=r["id"], rows=[r]) for r in population],
                dict(id=f"ranking-{digits}", rows=ranking),
            )
        )
    for phase in ("forced", "observer"):
        for row in rows:
            if row["phase"] == phase:
                calls.append(dict(id=row["id"], rows=[row]))
    assert sum(len(c["rows"]) for c in calls) == 22
    return calls


def limits(plan):
    rows = [r for r in plan["calibration"]["calls"] if r["phase"] != "primary"]
    slot = plan["calibration"]["limits"]["warm_slot_work"]
    phase_limits = {}
    for phase in {r["phase"] for r in rows}:
        group = [r for r in rows if r["phase"] == phase]
        phase_limits[phase] = dict(
            wall=sum(r["reservation_seconds"] + 25 for r in group),
            cpu=sum(r["reservation_seconds"] + 25 for r in group),
            work=sum(r["work"] + slot for r in group),
        )
    parent = plan["calibration"]["limits"]["total"]
    primary = execute_v4_stage.limits(plan)["primary"]
    total = {f: parent[f] - primary[f] for f in ("wall", "cpu", "work")}
    return dict(**phase_limits, total=total, disk_bytes=64 * 2**20)


def run_child(call):
    rows = call["rows"]
    group = len(rows) > 1
    if group:
        command = [
            sys.executable,
            "-u",
            "-m",
            "v2.benchmarks.qs.c10a01.execute_v4_group",
            rows[0]["phase"],
            str(rows[0]["digits"]),
        ]
    else:
        command = [
            sys.executable,
            "-u",
            "-m",
            "v2.benchmarks.qs.c10a01.execute_v4",
            rows[0]["id"],
        ]
    maximum = sum(r["reservation_seconds"] + 25 for r in rows) + 5
    start = time.monotonic()
    process = subprocess.Popen(
        command,
        cwd=plan_v4.ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        env=dict(
            os.environ,
            C10A01_V4_OUTPUT=str(OUTPUT),
            PYTHONDONTWRITEBYTECODE="1",
        ),
    )
    try:
        result = execute_v4_stage.read_worker(process, start + maximum)
    finally:
        if process.poll() is None:
            process.kill()
            process.wait()
    if group:
        if result["reason"] != "diagnostic_complete":
            raise RuntimeError("incomplete diagnostic group")
        if [m["id"] for m in result["members"]] != [r["id"] for r in rows]:
            raise RuntimeError("diagnostic member mismatch")
        for member, row in zip(result["members"], rows):
            if len(member["samples"]) != 9:
                raise RuntimeError("diagnostic samples incomplete")
            if min(member["before"]["wall"], member["before"]["cpu"]) < 3:
                raise RuntimeError("diagnostic warmup too short")
            if (
                member["measured_wall"] > row["seconds"]
                or member["measured_cpu"] > row["seconds"]
            ):
                raise RuntimeError("diagnostic measured cap")
    else:
        result = execute_v4_stage.run_worker_result(result, rows[0])
    if result["work"] > sum(
        r["work"] + execute_v4_stage.warm_slot_work() for r in rows
    ):
        raise RuntimeError("diagnostic group work cap")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("check", "run"))
    args = parser.parse_args()
    plan, release, available = check_release(
        require_prior=args.action == "run"
    )
    calls = schedule(plan)
    if args.action == "check":
        print(
            json.dumps(
                dict(
                    source_verified=True,
                    prior_evidence_available=available,
                    calls=sum(len(c["rows"]) for c in calls),
                    groups=len(calls),
                    performance_calls=0,
                )
            )
        )
        return
    OUTPUT.mkdir(parents=True, exist_ok=True)
    ledger = Ledger(OUTPUT, digest(RELEASE), limits(plan))
    for call in calls:
        name = call["id"]
        if any(r["receipt"] == name + ".json" for r in ledger.data["calls"]):
            continue
        rows = call["rows"]
        for tier in {r["tier"] for r in rows}:
            execute_v4_stage.host_admission(tier)
        seconds = sum(r["reservation_seconds"] + 25 for r in rows)
        work = sum(r["work"] + execute_v4_stage.warm_slot_work() for r in rows)
        identity = dict(
            group=name,
            cells=[r["id"] for r in rows],
            release=digest(RELEASE),
            primary=digest(PRIMARY / "ledger.json"),
        )
        result = ledger.call(
            rows[0]["phase"],
            name,
            seconds,
            work,
            identity,
            lambda: run_child(call),
        )
        if rows[0]["phase"] == "observer":
            other_flag = int(not rows[0]["observer"])
            other = OUTPUT / f"observer-{rows[0]['digits']}-{other_flag}.json"
            if other.exists():
                previous = json.loads(other.read_text())["result"]["result"]
                current = result["result"]
                if (
                    previous["semantic_digest"] != current["semantic_digest"]
                    or previous["work"] != current["work"]
                    or previous["positions"] != current["positions"]
                ):
                    raise RuntimeError("observer changed semantic work")
        print(
            json.dumps(
                dict(
                    group=name,
                    cells=len(rows),
                    reason=result.get(
                        "reason", result.get("result", {}).get("reason")
                    ),
                    work=result["work"],
                )
            ),
            flush=True,
        )


if __name__ == "__main__":
    main()
