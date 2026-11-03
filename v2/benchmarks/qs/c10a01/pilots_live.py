"""Four-hour 70-digit restart with acknowledged live grant extensions."""

import argparse
import hashlib
import json
import os
import selectors
import statistics
import subprocess
import sys
import time
from pathlib import Path

from ..c10a04.accounting import digest, output_lock, receipt
from . import execute_v4, live_limits, nominee_v4, plan_v4
from . import pilots_v4 as original

ROOT = plan_v4.ROOT
FREEZE = Path(__file__).with_name("pilot_live.json")
SOURCES = Path(__file__).with_name("pilot_live_sources.json")
OUTPUT = ROOT / "v2/benchmarks/results/c10a01/pilot-live-70"
INTERRUPTED = ROOT / "v2/benchmarks/results/c10a01/pilot-v4-C/ledger.json"
INTERRUPTED_SHA = (
    "9f9404abd77a09279b0c630c9a7c411a99f3026518ec52dc40f957d6cf91265f"
)
MODULE = "v2.benchmarks.qs.c10a01.pilots_live"
DISK_BYTES = 16 * 2**20


def build_freeze():
    nominee, plan = original.source_check()
    calls = []
    for source in plan["pilots"]["C"]["calls"]:
        row = dict(source)
        selected = next(c for c in nominee["calls"] if c["id"] == row["id"])
        row.update(
            seconds=14400,
            reservation_seconds=14405,
            work=plan_v4.call_work(70, 14400),
            configuration=selected["configuration"],
            configuration_status="frozen_live_nominee_or_control",
        )
        calls.append(row)
    limits = plan_v4.totals(calls, 12, 60, 30)
    assert limits["total"]["wall"] == 29200
    return dict(
        identity="c10a01-live-four-hour-restart-1",
        calls=calls,
        limits=limits,
        previous_ledger_sha256=INTERRUPTED_SHA,
        previous_wall_cpu=7240,
        previous_work=81438260391796876608,
        original_nominee_sha256=digest(nominee_v4.FREEZE),
        extension_semantics=1,
        source_policies="unchanged",
        disk_bytes=DISK_BYTES,
        live_state_bytes=live_limits.STATE_BYTES,
        maximum_extensions=live_limits.MAX_REVISIONS,
    )


def source_check(require_receipts=False):
    freeze = json.loads(FREEZE.read_text())
    if freeze != build_freeze():
        raise RuntimeError("live restart input/configuration/grant drift")
    sources = json.loads(SOURCES.read_text())
    for name, expected in sources["files"].items():
        if digest(ROOT / name) != expected:
            raise RuntimeError("live runner source drift: " + name)
    if sources["freeze_sha256"] != digest(FREEZE):
        raise RuntimeError("live release drift")
    if require_receipts:
        nominee_v4.verify(require_receipts=True)
        if digest(INTERRUPTED) != INTERRUPTED_SHA:
            raise RuntimeError("interrupted original reservation changed")
    return freeze


def disk_check():
    used = sum(p.stat().st_size for p in OUTPUT.rglob("*") if p.is_file())
    if used + 2**20 > DISK_BYTES:
        raise RuntimeError("recorder_exhaustion: live stage storage")


def publish(name, result):
    if len(json.dumps(result).encode()) > original.RECEIPT_LIMIT - 4096:
        raise RuntimeError("recorder_exhaustion: live receipt")
    disk_check()
    path = OUTPUT / (name + ".json")
    receipt(path, result)
    return dict(name=path.name, sha256=digest(path))


def phase_receipt(identity, name, phase, function):
    started, cpu = time.monotonic(), time.process_time()
    result = function()
    wall, used_cpu = time.monotonic() - started, time.process_time() - cpu
    reference = publish(
        name + "-" + phase,
        dict(
            result=result,
            actual_wall=wall,
            actual_cpu=used_cpu,
            charge_wall=max(35, wall),
            charge_cpu=max(35, used_cpu),
            work=result["work"],
        ),
    )
    reference.update(
        wall=max(35, wall), cpu=max(35, used_cpu), work=result["work"]
    )
    with live_limits.transaction(OUTPUT, identity) as state:
        state["calls"][name]["receipts"][phase] = reference
    if wall > 35 or used_cpu > 35 or result["work"] > 290_000_000_000:
        raise RuntimeError("sentinel reservation exceeded; receipt retained")
    return result


def measured(row, fixture, modules, identity, *, blocks=None):
    """Attach the live grant to one job; never recreate it for an extension."""
    started, cpu = time.monotonic(), time.process_time()
    job, budget, cfg = execute_v4.instantiate(
        row, fixture, modules, seconds=row["seconds"], work=row["work"]
    )
    live_limits.LiveBudget(OUTPUT, identity, row["id"], budget)
    try:
        result = job.run(max_blocks=blocks)
    except MemoryError:
        # No emergency serialization: preserve native refusal and exact work.
        from types import SimpleNamespace

        result = SimpleNamespace(
            divisor=None,
            cofactor=fixture["n"],
            reason="native_memory_refusal",
            stats=job.stats,
        )
    validated = execute_v4.validate_result(
        fixture, result, budget, modules[2], modules[1].BudgetExhaustedError
    )
    positions = (
        0 if job.engine is None else job.engine.stats["collected_positions"]
    )
    state_digest = hashlib.sha256()
    state_digest.update(
        str((job.multiplier, job.family_index, positions)).encode()
    )
    if job.engine is not None:
        for key in sorted(job.engine.collector._atoms):
            state_digest.update(str(key).encode())
    reason = result.reason
    if reason == "cancelled" and execute_v4.rss() > plan_v4.RSS[row["tier"]]:
        reason = "rss_limit"
    return dict(
        **validated,
        reason=reason,
        work=budget.used,
        seconds=time.monotonic() - started,
        cpu_seconds=time.process_time() - cpu,
        positions=positions,
        polynomials=job.stats["polynomials"],
        families=job.family_index,
        multiplier=job.multiplier,
        d_mod8=None if job.base is None else job.base.n_prime % 8,
        stats=result.stats,
        rss_high_water_bytes=execute_v4.rss(),
        semantic_digest=state_digest.hexdigest(),
        final_seconds=budget.seconds,
        final_work_limit=budget.work_limit,
        complete=bool(validated["factors"]) and not validated["unresolved"],
    )


def worker(index):
    freeze = source_check(require_receipts=True)
    identity = digest(FREEZE)
    state = live_limits.read(OUTPUT)
    if state["identity"] != identity:
        raise RuntimeError("unexpected live stage")
    name = freeze["calls"][index]["id"]
    row = dict(freeze["calls"][index])
    for key in ("seconds", "work", "reservation_seconds"):
        row[key] = state["calls"][name][key]
    original.host_admission(row["tier"])
    execute_v4.OUTPUT = OUTPUT
    plan_v4.RSS = {k: v - original.PARENT_RSS for k, v in plan_v4.RSS.items()}
    fixture = execute_v4.identity(row, json.loads(plan_v4.CORPUS.read_text()))
    modules = execute_v4.runtime_modules(row)
    live_limits.mark_phase(OUTPUT, identity, name, "warming")
    before = phase_receipt(
        identity,
        name,
        "before",
        lambda: original.sentinel(row, fixture, modules),
    )
    current = live_limits.read(OUTPUT)["calls"][name]
    for key in ("seconds", "work", "reservation_seconds"):
        row[key] = current[key]
    result = measured(row, fixture, modules, identity)
    reference = publish(name, result)
    live_limits.close_measurement(
        OUTPUT, identity, name, reference["name"], reference["sha256"], result
    )
    granted = live_limits.read(OUTPUT)["calls"][name]["seconds"]
    result["capped_service_wall"] = (
        result["seconds"] if result["complete"] else granted
    )
    result["capped_service_cpu"] = (
        result["cpu_seconds"] if result["complete"] else granted
    )
    after = None
    if result["reason"] not in ("rss_limit", "native_memory_refusal"):
        after = phase_receipt(
            identity,
            name,
            "after",
            lambda: original.sentinel(row, fixture, modules),
        )
    drift = (
        None
        if after is None
        else {
            key: statistics.median(s[key] for s in after["samples"])
            / statistics.median(s[key] for s in before["samples"])
            for key in ("wall", "cpu")
        }
    )
    summary = dict(
        name=name,
        result={k: v for k, v in result.items() if k != "stats"},
        drift=drift,
        drift_pass=drift is not None
        and all(abs(value - 1) <= 0.15 for value in drift.values()),
    )
    reference = publish(name + "-summary", summary)
    with live_limits.transaction(OUTPUT, identity) as state:
        state["calls"][name]["status"] = "complete"
        state["calls"][name]["receipts"]["summary"] = reference
    print(json.dumps(summary), flush=True)


def launch(index, name, identity):
    if execute_v4.rss() > original.PARENT_RSS:
        raise RuntimeError("parent RSS exceeds its suballocation")
    process = subprocess.Popen(
        [sys.executable, "-m", MODULE, "worker", "--cell-index", str(index)],
        cwd=ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    selector = selectors.DefaultSelector()
    selector.register(process.stdout, selectors.EVENT_READ, "stdout")
    selector.register(process.stderr, selectors.EVENT_READ, "stderr")
    captured = {"stdout": bytearray(), "stderr": bytearray()}
    started, next_read = time.monotonic(), 0.0
    deadline = (
        started + live_limits.read(OUTPUT)["calls"][name]["seconds"] + 80
    )
    try:
        while selector.get_map():
            now = time.monotonic()
            if now >= next_read or now >= deadline - 1:
                state = live_limits.read(OUTPUT)
                if state["identity"] != identity:
                    raise RuntimeError("supervisor live identity drift")
                deadline = live_limits.watchdog_deadline(state, name, started)
                next_read = now + 1
            if now >= deadline:
                raise TimeoutError("live worker outer grant exhausted")
            for key, _ in selector.select(timeout=0.1):
                chunk = os.read(key.fileobj.fileno(), 4096)
                if not chunk:
                    selector.unregister(key.fileobj)
                    continue
                captured[key.data].extend(chunk)
                if len(captured[key.data]) > original.OUTPUT_LIMIT:
                    raise RuntimeError(
                        "recorder_exhaustion: live child output"
                    )
        if process.wait(timeout=1):
            raise RuntimeError(
                captured["stderr"].decode(errors="replace")[-2048:]
            )
        summary = json.loads(captured["stdout"])
        if summary["name"] != name:
            raise RuntimeError("live child identity mismatch")
        return summary
    except BaseException as error:
        if process.poll() is None:
            process.kill()
        process.wait()
        live_limits.interrupt(OUTPUT, identity, name, str(error))
        raise
    finally:
        selector.close()


def run():
    freeze = source_check(require_receipts=True)
    identity = digest(FREEZE)
    for call in freeze["calls"]:
        original.host_admission(call["tier"])
    with output_lock(OUTPUT):
        live_limits.initialize(OUTPUT, identity, freeze["calls"])
        rows = []
        for index, call in enumerate(freeze["calls"]):
            original.host_admission(call["tier"])
            live_limits.check_capacity(live_limits.read(OUTPUT))
            summary = launch(index, call["id"], identity)
            rows.append(summary)
            print(json.dumps(summary), flush=True)
        live_limits.check_capacity(live_limits.read(OUTPUT))
        reference = publish("analysis", dict(rows=rows, adoption=False))
        with live_limits.transaction(OUTPUT, identity) as state:
            state["analysis"] = reference
            state["status"] = "complete"
            overhead = 90 + live_limits.EXTENSION_OVERHEAD * state["revision"]
            state["overhead_charged"] = dict(
                wall=overhead, cpu=overhead, work=0
            )


def acknowledgement(names, revision, *, wait_seconds=5):
    """Report application separately from a durably reserved request."""
    deadline = time.monotonic() + wait_seconds
    while True:
        state = live_limits.read(OUTPUT)
        replies = {}
        waiting = False
        for name in names:
            call = state["calls"][name]
            if call["acknowledged"] >= revision:
                status = "applied"
            elif call["status"] in ("pending", "warming"):
                status = "reserved_for_start"
            elif call["status"] == "measured":
                status = "awaiting_worker"
                waiting = True
            else:
                status = "closed_without_application"
            replies[name] = dict(
                status=status,
                acknowledged=call["acknowledged"],
                seconds=call["seconds"],
                pid=call["pid"],
                used_at_ack=call["used_at_ack"],
            )
        if not waiting or time.monotonic() >= deadline:
            return replies
        time.sleep(0.05)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "action", choices=("check", "run", "worker", "status", "extend")
    )
    parser.add_argument("--cell-index", type=int)
    parser.add_argument("--seconds", type=int)
    parser.add_argument("--expect-revision", type=int)
    parser.add_argument("--reason")
    parser.add_argument(
        "--arm", choices=("nominee", "control", "both"), default="both"
    )
    args = parser.parse_args()
    if args.action == "check":
        freeze = source_check()
        print(
            json.dumps(
                dict(
                    ready=True,
                    initial_limits=freeze["limits"]["total"],
                    hosts={3: original.host_admission(3)},
                    performance_calls=0,
                )
            )
        )
    elif args.action == "worker":
        if args.cell_index not in (0, 1):
            parser.error("worker requires cell index 0 or 1")
        worker(args.cell_index)
    elif args.action == "run":
        run()
    elif args.action == "status":
        state = live_limits.read(OUTPUT)
        fields = (
            "status",
            "seconds",
            "work",
            "pid",
            "started",
            "revision",
            "acknowledged",
            "used_at_ack",
            "reason",
        )
        print(
            json.dumps(
                dict(
                    status=state["status"],
                    revision=state["revision"],
                    limits=state["limits"]["total"],
                    calls={
                        name: {k: call[k] for k in fields if k in call}
                        for name, call in state["calls"].items()
                    },
                ),
                indent=2,
            )
        )
    else:
        source_check(require_receipts=True)
        state = live_limits.read(OUTPUT)
        names = [
            name
            for name, call in state["calls"].items()
            if args.arm == "both" or call["arm"] == args.arm
        ]
        result = live_limits.extend(
            OUTPUT,
            digest(FREEZE),
            names,
            args.seconds,
            args.expect_revision,
            args.reason,
        )
        result["acknowledgements"] = acknowledgement(names, result["revision"])
        result["stage_limits"] = live_limits.read(OUTPUT)["limits"]["total"]
        print(json.dumps(result))


if __name__ == "__main__":
    main()
