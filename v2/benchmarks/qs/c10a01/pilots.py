"""Frozen complete-factor feasibility pilots; run only in a later session."""

import argparse
import fcntl
import json
import os
import statistics
import subprocess
import sys
import time
from dataclasses import replace
from pathlib import Path

from ..c10a04.accounting import Ledger, digest, output_lock
from . import runner
from .configuration import TIERS, configuration
from .work_model import charges, work_limit

PLAN = Path(__file__).with_name("pilots_v2.json")
SOURCES = Path(__file__).with_name("pilot_sources_v2.json")
PARENT_OWNED = 32 * 2**20
PARENT_RSS = 256 * 2**20
CELL_RECEIPT_BYTES = 64 * 1024


def pilot_configuration(digits, arm, *, multiplier=0):
    """Suballocate parent ownership without changing any algorithm policy."""
    config = configuration(digits, arm, multiplier=multiplier)
    return replace(config)


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
    assert plan["identity"] == "c10a01-pilot-v2"
    assert plan["parent_owned_reserve"] == PARENT_OWNED
    assert plan["parent_rss_reserve"] == PARENT_RSS
    assert plan["source_manifest"] == SOURCES.name
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

            current = asdict(
                pilot_configuration(
                    call["digits"], call["arm"], multiplier=call["multiplier"]
                )
            )
            # The frozen v2 config predates this disabled opt-in. Enabled
            # policies still fail the equality; source pins remain strict.
            assert current.pop("root_update_policy") == "standard"
            assert call["configuration"] == current
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


def cell_result(call, result, before, after, name):
    """Bound parent data; keep full native counters in their receipt."""
    drift = None
    if after is not None and not after.get("missing"):
        drift = {
            key: statistics.median(r[key] for r in after["samples"])
            / statistics.median(r[key] for r in before["samples"])
            for key in ("wall", "cpu")
        }
    keys = (
        "n",
        "complete",
        "divisor",
        "cofactor",
        "factors",
        "certainty",
        "unresolved",
        "reason",
        "detail",
        "work",
        "seconds",
        "cpu_seconds",
        "rss_high_water_bytes",
        "capped_service_wall",
        "capped_service_cpu",
        "cold_import",
        "state_digest",
    )
    row = dict(
        call=call,
        result={key: result[key] for key in keys if key in result},
        native_receipt=name + ".json",
        process_id=os.getpid(),
        drift=drift,
        drift_pass=drift is not None
        and all(abs(v - 1) <= 0.15 for v in drift.values()),
        drift_missing_reason=(
            None if drift is not None else "native memory/RSS censoring"
        ),
        worker_rss_cap=runner.RSS[TIERS[call["digits"]]],
        parent_rss_reservation=PARENT_RSS,
        parent_owned_reservation=PARENT_OWNED,
    )
    if len(json.dumps(row).encode()) > CELL_RECEIPT_BYTES:
        raise RuntimeError("recorder_exhaustion: parent cell receipt bound")
    return row


def post_sentinel(fixture, call):
    """Preserve missing drift when the worker cannot recover after a call."""
    try:
        return sentinel(fixture, call)
    except (MemoryError, RuntimeError) as error:
        if (
            not isinstance(error, MemoryError)
            and runner.rss() <= (runner.RSS[TIERS[call["digits"]]])
        ):
            raise
        # Charge the whole reservation when the sentinel did not return its
        # cumulative work. This is censoring, never a successful drift check.
        return dict(
            work=290000000000,
            missing=True,
            reason="native_memory_or_rss",
            detail=str(error),
        )


def execute_cell(ledger, fixture, call, name, identity):
    """Run prewarm, one measured call and optional drift in a fresh worker."""
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
        before["warmup"]["last"]["multiplier"] == (call["expected_multiplier"])
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
    after = None
    if result["reason"] not in ("rss_limit", "native_memory_refusal"):
        after = ledger.call(
            "warmup",
            name + "-after",
            35,
            290000000000,
            cell,
            lambda: post_sentinel(fixture, call),
        )
    return cell_result(call, result, before, after, name)


def output_path(stage):
    return runner.ROOT / ("v2/benchmarks/results/c10a01/pilot-v2-" + stage)


def plan_identity(stage):
    return dict(
        plan=digest(PLAN),
        sources=digest(SOURCES),
        corpus=digest(runner.CORPUS),
        stage=stage,
    )


def worker(stage_name, index, lock_fd):
    """Child entry point; the parent owns both locks across every cell."""
    inherited = os.fstat(lock_fd)
    expected = os.stat("/private/tmp/factor-performance.lock")
    if (inherited.st_dev, inherited.st_ino) != (
        expected.st_dev,
        expected.st_ino,
    ):
        raise RuntimeError("worker requires the parent's performance lock")
    # The inherited open-file description shares the parent's lock. No second
    # output lock is taken: all child/parent ledger access is sequential.
    fcntl.flock(lock_fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
    plan, corpus, _ = readiness(stage_name)
    stage = plan["stages"][stage_name]
    call = stage["calls"][index]
    fixture = next(r for r in corpus["fixtures"] if r["id"] == call["fixture"])
    runner.configuration = pilot_configuration
    runner.RSS = {tier: cap - PARENT_RSS for tier, cap in runner.RSS.items()}
    ledger = BoundedLedger(
        output_path(stage_name), digest(PLAN), stage["limits"]
    )
    name = f"{index}-{call['digits']}-{call['arm']}"
    row = execute_cell(ledger, fixture, call, name, plan_identity(stage_name))
    print(json.dumps(row), flush=True)


def launch_worker(stage_name, index, call, lock_fd):
    """Reset lifetime RSS with exec; never retry a failed worker."""
    if runner.rss() > PARENT_RSS:
        raise RuntimeError("pilot parent exceeds its RSS suballocation")
    process = subprocess.run(
        [
            sys.executable,
            "-m",
            "v2.benchmarks.qs.c10a01.pilots",
            "worker",
            "--stage",
            stage_name,
            "--cell-index",
            str(index),
            "--lock-fd",
            str(lock_fd),
        ],
        cwd=runner.ROOT,
        pass_fds=(lock_fd,),
        capture_output=True,
        text=True,
        timeout=call["seconds"] + 75,
        check=True,
    )
    if len(process.stdout.encode()) > CELL_RECEIPT_BYTES:
        raise RuntimeError("recorder_exhaustion: worker output bound")
    row = json.loads(process.stdout)
    if row["call"] != call or row["process_id"] == os.getpid():
        raise RuntimeError("worker cell identity mismatch")
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("check", "run", "worker"))
    parser.add_argument("--stage", choices=("A", "B", "C"), required=True)
    parser.add_argument("--cell-index", type=int, help=argparse.SUPPRESS)
    parser.add_argument("--lock-fd", type=int, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.action == "worker":
        if args.cell_index is None or args.lock_fd is None:
            parser.error("worker requires an inherited parent lock and index")
        worker(args.stage, args.cell_index, args.lock_fd)
        return
    plan, corpus, hosts = readiness(args.stage)
    stage = plan["stages"][args.stage]
    if args.action == "check":
        if runner.rss() > PARENT_RSS:
            raise RuntimeError("pilot parent exceeds its RSS suballocation")
        print(
            json.dumps(
                dict(
                    stage=args.stage,
                    ready=True,
                    hosts=hosts,
                    limits=stage["limits"],
                    parent_rss=runner.rss(),
                    parent_rss_reservation=PARENT_RSS,
                    parent_owned_reservation=PARENT_OWNED,
                )
            )
        )
        return
    output = output_path(args.stage)
    identity = plan_identity(args.stage)
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
                runner.host_admission(call["digits"])
                row = launch_worker(args.stage, index, call, lock.fileno())
                fixture = next(
                    r for r in corpus["fixtures"] if r["id"] == call["fixture"]
                )
                runner.verify_outcome(row["result"], fixture)
                rows.append(row)
                print(
                    json.dumps(
                        dict(
                            call=index,
                            reason=row["result"]["reason"],
                            complete=row["result"]["complete"],
                            drift=row["drift"],
                        )
                    ),
                    flush=True,
                )
            # Workers mutated the durable ledger. Reload before parent writes.
            ledger = BoundedLedger(output, digest(PLAN), stage["limits"])
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
