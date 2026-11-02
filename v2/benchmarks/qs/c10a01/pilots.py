"""Frozen complete-factor feasibility pilots; run only in a later session."""

import argparse
import fcntl
import json
import statistics
import subprocess
import sys
import time
from pathlib import Path

from ..c10a04.accounting import Ledger, digest, output_lock
from . import runner
from .configuration import TIERS, configuration
from .work_model import charges, work_limit

PLAN = Path(__file__).with_name("pilots.json")
SOURCES = Path(__file__).with_name("pilot_sources.json")


class BoundedLedger(Ledger):
    """Keep each receipt below 1 MiB with wrapper/identity headroom."""

    def call(self, phase, name, seconds, work, identity, function):
        def bounded():
            result = function()
            if len(json.dumps(result).encode()) > 900 * 1024:
                raise RuntimeError("recorder_exhaustion: receipt bound")
            return result

        return super().call(phase, name, seconds, work, identity, bounded)


def verify_plan(plan):
    """Independently reconcile every call reservation and stage sum."""
    for stage in plan["stages"].values():
        calls = stage["calls"]
        measured = sum(c["seconds"] + 5 for c in calls)
        limits = stage["limits"]
        assert measured == limits["measured"]["wall"]
        assert sum(c["work"] for c in calls) == limits["measured"]["work"]
        for call in calls:
            assert call["work"] == work_limit(call["digits"], call["seconds"])
            assert call["tier"] == TIERS[call["digits"]]
            from dataclasses import asdict

            assert call["configuration"] == asdict(
                configuration(
                    call["digits"], call["arm"], multiplier=call["multiplier"]
                )
            )
            assert call["postprocessing"] == charges(call["digits"])
        for field in ("wall", "cpu", "work"):
            assert (
                sum(
                    limits[k][field]
                    for k in (
                        "measured",
                        "warmup",
                        "analysis",
                        "orchestration",
                    )
                )
                == (limits["total"][field])
            )
        assert limits["warmup"]["wall"] >= 70 * len(calls)
        assert limits["warmup"]["work"] >= 580000000000 * len(calls)
    assert (
        3600
        + sum(s["limits"]["total"]["wall"] for s in plan["stages"].values())
        == 35200
    )
    return True


def readiness(stage_name):
    runner.runtime()
    plan = json.loads(PLAN.read_text())
    verify_plan(plan)
    current_runtime = dict(
        version=sys.version,
        binary_sha256=digest(Path(sys.executable).resolve()),
    )
    if current_runtime != plan["runtime"]:
        raise RuntimeError("runtime drift requires a new pilot identity")
    sources = json.loads(SOURCES.read_text())
    for name, expected in sources["files"].items():
        if digest(runner.ROOT / name) != expected:
            raise RuntimeError("source drift: " + name)
    corpus = json.loads(runner.CORPUS.read_text())
    runner.validate(corpus)
    if digest(runner.CORPUS) != plan["corpus_sha256"]:
        raise RuntimeError("input drift")
    subprocess.run(
        ["git", "diff", "--exit-code", "HEAD", "--", "v2"],
        cwd=runner.ROOT,
        check=True,
        stdout=subprocess.DEVNULL,
    )
    stage = plan["stages"][stage_name]
    hosts = {
        d: runner.host_admission(d)
        for d in {c["digits"] for c in stage["calls"]}
    }
    return plan, corpus, hosts


def sentinel(fixture, call):
    """Warm both clocks, then retain nine short samples for drift evidence."""
    start, cpu = time.monotonic(), time.process_time()
    warming = runner.warm(
        fixture["n"],
        call["digits"],
        call["arm"],
        multiplier=call["multiplier"],
    )
    work, samples = warming["work"], []
    for _ in range(9):
        remaining = 30 - max(
            time.monotonic() - start, time.process_time() - cpu
        )
        if remaining <= 0:
            raise RuntimeError("sentinel reservation exhausted")
        row = runner.run_one(
            fixture["n"],
            call["digits"],
            call["arm"],
            remaining,
            290000000000 - work,
            264,
            multiplier=call["multiplier"],
        )
        runner.verify_outcome(row, fixture)
        if row["reason"] != "paused" or row["positions"] != 8 * 131073:
            raise RuntimeError("sentinel did not finish fixed work")
        work += row["work"]
        samples.append(
            dict(
                wall=row["seconds"],
                cpu=row["cpu_seconds"],
                work=row["work"],
                digest=row["state_digest"],
            )
        )
    assert len({r["digest"] for r in samples}) == 1
    return dict(
        work=work,
        warmup=warming,
        samples=samples,
        wall=time.monotonic() - start,
        cpu=time.process_time() - cpu,
    )


def complete_call(fixture, call):
    start, cpu = time.monotonic(), time.process_time()
    try:
        row = runner.run_one(
            fixture["n"],
            call["digits"],
            call["arm"],
            call["seconds"],
            call["work"],
            None,
            multiplier=call["multiplier"],
        )
    except MemoryError as error:
        # No emergency checkpoint: serialization has a separate ownership gap.
        job = runner.LAST_JOB
        row = dict(
            n=fixture["n"],
            complete=False,
            divisor=None,
            cofactor=fixture["n"],
            factors=[],
            certainty=[],
            unresolved=[fixture["n"]],
            reason="native_memory_refusal",
            detail=str(error),
            work=0 if job is None else job.budget.used,
            seconds=time.monotonic() - start,
            cpu_seconds=time.process_time() - cpu,
            rss_high_water_bytes=runner.rss(),
        )
    runner.verify_outcome(row, fixture)
    row["capped_service_wall"] = (
        row["seconds"] if row["complete"] else call["seconds"]
    )
    row["capped_service_cpu"] = (
        row["cpu_seconds"] if row["complete"] else call["seconds"]
    )
    row["cold_import"] = dict(
        wall=runner.IMPORT_COST[0], cpu=runner.IMPORT_COST[1]
    )
    return row


def analyze(rows):
    """Report every pair, censoring and pre/post drift; no adoption test."""
    pairs = []
    for digits in sorted({r["call"]["digits"] for r in rows}):
        selected = [r for r in rows if r["call"]["digits"] == digits]
        by_arm = {r["call"]["arm"]: r for r in selected}
        control, candidate = by_arm["control"], by_arm["parity2"]
        pair = dict(
            digits=digits,
            calls=selected,
            both_complete=all(r["result"]["complete"] for r in selected),
        )
        if pair["both_complete"]:
            pair["descriptive_wall_ratio"] = (
                candidate["result"]["seconds"] / control["result"]["seconds"]
            )
            pair["descriptive_cpu_ratio"] = (
                candidate["result"]["cpu_seconds"]
                / control["result"]["cpu_seconds"]
            )
        pairs.append(pair)
    return dict(
        work=1000000 * len(rows),
        pairs=pairs,
        independent_inputs=len(pairs),
        adoption=False,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("check", "run"))
    parser.add_argument("--stage", choices=("A", "B", "C"), required=True)
    args = parser.parse_args()
    plan, corpus, hosts = readiness(args.stage)
    stage = plan["stages"][args.stage]
    if args.action == "check":
        print(
            json.dumps(
                dict(
                    stage=args.stage,
                    ready=True,
                    hosts=hosts,
                    limits=stage["limits"],
                )
            )
        )
        return
    output = runner.ROOT / (
        "v2/benchmarks/results/c10a01/pilot-v1-" + args.stage
    )
    identity = dict(
        plan=digest(PLAN),
        sources=digest(SOURCES),
        corpus=digest(runner.CORPUS),
        stage=args.stage,
    )
    with open("/private/tmp/factor-performance.lock", "a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        with output_lock(output):
            ledger = BoundedLedger(output, digest(PLAN), stage["limits"])
            if ledger.data["calls"]:
                raise RuntimeError(
                    "stage already started; preserve outcomes, no retry"
                )
            rows = []
            for index, call in enumerate(stage["calls"]):
                fixture = next(
                    r for r in corpus["fixtures"] if r["id"] == call["fixture"]
                )
                runner.host_admission(call["digits"])
                name = f"{index}-{call['digits']}-{call['arm']}"
                cell = dict(identity, **call)
                before = ledger.call(
                    "warmup",
                    name + "-before",
                    35,
                    290000000000,
                    cell,
                    lambda: sentinel(fixture, call),
                )
                assert (
                    before["warmup"]["last"]["multiplier"]
                    == (call["expected_multiplier"])
                )
                assert before["warmup"]["last"]["stats"]["capacity"][
                    "target_in_product_envelope"
                ]
                result = ledger.call(
                    "measured",
                    name,
                    call["seconds"] + 5,
                    call["work"],
                    cell,
                    lambda: complete_call(fixture, call),
                )
                after = ledger.call(
                    "warmup",
                    name + "-after",
                    35,
                    290000000000,
                    cell,
                    lambda: sentinel(fixture, call),
                )
                drift = {
                    key: statistics.median(r[key] for r in after["samples"])
                    / statistics.median(r[key] for r in before["samples"])
                    for key in ("wall", "cpu")
                }
                rows.append(
                    dict(
                        call=call,
                        result=result,
                        drift=drift,
                        drift_pass=all(
                            abs(v - 1) <= 0.15 for v in drift.values()
                        ),
                    )
                )
                print(
                    json.dumps(
                        dict(
                            call=name,
                            reason=result["reason"],
                            complete=result["complete"],
                            drift=drift,
                        )
                    ),
                    flush=True,
                )
            ledger.call(
                "analysis",
                "analysis",
                stage["limits"]["analysis"]["wall"],
                stage["limits"]["analysis"]["work"],
                identity,
                lambda: analyze(rows),
            )
            ledger.call(
                "orchestration",
                "orchestration",
                stage["limits"]["orchestration"]["wall"],
                0,
                identity,
                lambda: dict(work=0, scope="conservative overhead charge"),
            )


if __name__ == "__main__":
    main()
