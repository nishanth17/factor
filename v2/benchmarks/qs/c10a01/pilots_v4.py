"""Source-pinned complete-factor feasibility pilots for the v4 nominee."""

import argparse
import json
import os
import selectors
import statistics
import subprocess
import sys
import time
from pathlib import Path

from ..c10a04.accounting import Ledger, digest, output_lock
from . import control_v4, execute_v4, nominee_v4, plan_v4

ROOT = plan_v4.ROOT
SOURCES = Path(__file__).with_name("pilot_v4_sources.json")
PARENT_RSS = 256 * 2**20
SENTINEL_BLOCKS = 8 * 33
OUTPUT_LIMIT = 64 * 1024
RECEIPT_LIMIT = 512 * 1024


class BoundedLedger(Ledger):
    """Refuse an oversized receipt as a recorder failure."""

    def call(self, phase, name, seconds, work, identity, function):
        def bounded():
            result = function()
            if len(json.dumps(result).encode()) > RECEIPT_LIMIT - 4096:
                raise RuntimeError("recorder_exhaustion: cell receipt")
            return result

        return super().call(phase, name, seconds, work, identity, bounded)


def source_check():
    """Check committed-only runtime sources and the exact selected policies."""
    if sys.implementation.name != "pypy" or sys.version_info[:2] != (3, 11):
        raise RuntimeError("PyPy implementing Python 3.11 required")
    nominee = nominee_v4.verify()
    sources = json.loads(SOURCES.read_text())
    if sources["nominee_sha256"] != digest(nominee_v4.FREEZE):
        raise RuntimeError("nominee release drift")
    if sources["runtime_binary_sha256"] != digest(
        Path(sys.executable).resolve()
    ):
        raise RuntimeError("PyPy binary drift")
    if sources["runtime_version"] != sys.version:
        raise RuntimeError("PyPy version drift")
    for name, expected in sources["files"].items():
        if digest(ROOT / name) != expected:
            raise RuntimeError("pilot source drift: " + name)
    inherited = json.loads(plan_v4.SOURCES.read_text())
    for name, expected in inherited["files"].items():
        if digest(ROOT / name) != expected:
            raise RuntimeError("inherited production source drift: " + name)
    plan = json.loads(plan_v4.PLAN.read_text())
    assert plan_v4.verify_plan(plan)
    assert control_v4.checked_baseline()["commit"] == nominee["control_commit"]
    for stage in ("A", "B", "C"):
        frozen = [
            c
            for c in nominee["calls"]
            if c["id"] in {r["id"] for r in plan["pilots"][stage]["calls"]}
        ]
        assert len(frozen) == len(plan["pilots"][stage]["calls"])
        assert nominee["stages"][stage] == plan["pilots"][stage]["limits"]
    return nominee, plan


def host_admission(tier):
    # Reuse the already exercised calibration host admission calculation.
    from .execute_v4_stage import host_admission as admit

    return admit(tier)


def stage_check(stage, *, require_receipts=False):
    nominee, plan = source_check()
    nominee_v4.verify(require_receipts=require_receipts)
    hosts = {
        row["tier"]: host_admission(row["tier"])
        for row in nominee["calls"]
        if row["id"] in {c["id"] for c in plan["pilots"][stage]["calls"]}
    }
    return nominee, plan, hosts


def semantic_run(row, fixture, modules, *, seconds, work, blocks):
    result = execute_v4.sieve(
        row, fixture, modules, seconds=seconds, work=work, blocks=blocks
    )
    if result["positions"] < 0 or result["work"] > work:
        raise RuntimeError("invalid semantic work")
    if (
        result["rss_high_water_bytes"] > plan_v4.RSS[row["tier"]]
        and result["reason"] == "cancelled"
    ):
        result["reason"] = "rss_limit"
    return result


def sentinel(row, fixture, modules):
    """Warm both clocks and measure nine identical semantic prefixes."""
    started, cpu = time.monotonic(), time.process_time()
    used, warmed, samples = 0, 0, []
    while min(time.monotonic() - started, time.process_time() - cpu) < 3:
        if max(time.monotonic() - started, time.process_time() - cpu) > 18:
            raise RuntimeError("warmup reservation exhausted")
        result = semantic_run(
            row,
            fixture,
            modules,
            seconds=20,
            work=290_000_000_000 - used,
            blocks=SENTINEL_BLOCKS,
        )
        if result["reason"] not in ("paused", "factor_found"):
            raise RuntimeError("warmup semantic refusal")
        used += result["work"]
        warmed += 1
    for _ in range(9):
        result = semantic_run(
            row,
            fixture,
            modules,
            seconds=20,
            work=290_000_000_000 - used,
            blocks=SENTINEL_BLOCKS,
        )
        if result["reason"] != "paused" or result["positions"] != 8 * 131073:
            raise RuntimeError("sentinel fixed work did not finish")
        used += result["work"]
        samples.append(
            {
                "wall": result["seconds"],
                "cpu": result["cpu_seconds"],
                "work": result["work"],
                "digest": result["semantic_digest"],
            }
        )
    if len({s["digest"] for s in samples}) != 1:
        raise RuntimeError("sentinel state drift")
    return {
        "work": used,
        "wall": time.monotonic() - started,
        "cpu": time.process_time() - cpu,
        "warmed_calls": warmed,
        "samples": samples,
    }


def summarize(result, call):
    complete = bool(result["factors"]) and not result["unresolved"]
    return {
        key: result[key]
        for key in (
            "reason",
            "work",
            "seconds",
            "cpu_seconds",
            "positions",
            "polynomials",
            "families",
            "multiplier",
            "d_mod8",
            "rss_high_water_bytes",
            "factors",
            "certainty",
            "unresolved",
            "divisor",
            "cofactor",
            "semantic_digest",
        )
    } | {
        "complete": complete,
        "capped_wall": result["seconds"] if complete else call["seconds"],
        "capped_cpu": result["cpu_seconds"] if complete else call["seconds"],
    }


def run_cell(stage, index, nominee, plan):
    call = plan["pilots"][stage]["calls"][index]
    frozen = next(c for c in nominee["calls"] if c["id"] == call["id"])
    corpus = json.loads(plan_v4.CORPUS.read_text())
    fixture = execute_v4.identity(call, corpus)
    if (
        frozen["configuration"] != call["configuration"]
        and call["arm"] == "control"
    ):
        raise RuntimeError("control configuration drift")
    row = dict(call, configuration=frozen["configuration"])
    execute_v4.OUTPUT = ROOT / (
        "v2/benchmarks/results/c10a01/pilot-v4-" + stage
    )
    execute_v4.plan_v4.RSS = {
        tier: cap - PARENT_RSS for tier, cap in plan_v4.RSS.items()
    }
    modules = execute_v4.runtime_modules(row)
    out = execute_v4.OUTPUT
    ledger = BoundedLedger(
        out, digest(nominee_v4.FREEZE), plan["pilots"][stage]["limits"]
    )
    identity = {
        "stage": stage,
        "cell": call["id"],
        "nominee": digest(nominee_v4.FREEZE),
    }
    before = ledger.call(
        "warmup",
        call["id"] + "-before",
        35,
        290_000_000_000,
        identity,
        lambda: sentinel(row, fixture, modules),
    )
    result = ledger.call(
        "measured",
        call["id"],
        call["seconds"] + 5,
        call["work"],
        identity,
        lambda: semantic_run(
            row,
            fixture,
            modules,
            seconds=call["seconds"],
            work=call["work"],
            blocks=None,
        ),
    )
    after = None
    if result["reason"] != "rss_limit":
        after = ledger.call(
            "warmup",
            call["id"] + "-after",
            35,
            290_000_000_000,
            identity,
            lambda: sentinel(row, fixture, modules),
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
    output = {
        "call": call["id"],
        "result": summarize(result, call),
        "drift": drift,
        "drift_pass": drift is not None
        and all(abs(value - 1) <= 0.15 for value in drift.values()),
        "before_work": before["work"],
        "after_work": None if after is None else after["work"],
    }
    if len(json.dumps(output).encode()) > OUTPUT_LIMIT:
        raise RuntimeError("recorder_exhaustion: worker summary")
    return output


def stage_run(stage):
    nominee, plan, _ = stage_check(stage, require_receipts=True)
    out = ROOT / ("v2/benchmarks/results/c10a01/pilot-v4-" + stage)
    with output_lock(out):
        ledger = BoundedLedger(
            out, digest(nominee_v4.FREEZE), plan["pilots"][stage]["limits"]
        )
        if ledger.data["calls"]:
            raise RuntimeError("stage already started; preserve all outcomes")
        rows = []
        for index, call in enumerate(plan["pilots"][stage]["calls"]):
            host_admission(call["tier"])
            command = [
                sys.executable,
                "-m",
                __name__,
                "worker",
                "--stage",
                stage,
                "--cell-index",
                str(index),
            ]
            child = subprocess.Popen(
                command,
                cwd=ROOT,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            selector = selectors.DefaultSelector()
            selector.register(child.stdout, selectors.EVENT_READ, "stdout")
            selector.register(child.stderr, selectors.EVENT_READ, "stderr")
            captured = {"stdout": bytearray(), "stderr": bytearray()}
            deadline = time.monotonic() + call["seconds"] + 80
            try:
                while selector.get_map():
                    if time.monotonic() >= deadline:
                        raise TimeoutError("pilot worker outer reservation")
                    for key, _ in selector.select(timeout=0.1):
                        chunk = os.read(key.fileobj.fileno(), 4096)
                        if not chunk:
                            selector.unregister(key.fileobj)
                            continue
                        captured[key.data].extend(chunk)
                        if len(captured[key.data]) > OUTPUT_LIMIT:
                            raise RuntimeError(
                                "recorder_exhaustion: child output"
                            )
                if child.wait(timeout=1):
                    raise RuntimeError(
                        f"worker failed at {call['id']}: "
                        + captured["stderr"].decode(errors="replace")[-2000:]
                    )
            except BaseException:
                child.kill()
                child.wait()
                raise
            finally:
                selector.close()
            row = json.loads(captured["stdout"])
            rows.append(row)
            print(
                json.dumps({"cell": call["id"], "result": row["result"]}),
                flush=True,
            )
        # Workers settle measured and warmup reservations in the ledger.
        ledger = BoundedLedger(
            out, digest(nominee_v4.FREEZE), plan["pilots"][stage]["limits"]
        )
        ledger.call(
            "analysis",
            "analysis",
            plan["pilots"][stage]["limits"]["phases"]["analysis"]["wall"],
            1_000_000_000_000,
            {"stage": stage},
            lambda: {"work": 0, "rows": rows, "adoption": False},
        )
        ledger.call(
            "orchestration",
            "orchestration",
            plan["pilots"][stage]["limits"]["phases"]["orchestration"]["wall"],
            0,
            {"stage": stage},
            lambda: {"work": 0},
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("check", "run", "worker"))
    parser.add_argument("--stage", choices=("A", "B", "C"), required=True)
    parser.add_argument("--cell-index", type=int)
    args = parser.parse_args()
    if args.action == "worker":
        nominee, plan, _ = stage_check(args.stage, require_receipts=True)
        if args.cell_index is None:
            parser.error("worker needs a cell index")
        print(json.dumps(run_cell(args.stage, args.cell_index, nominee, plan)))
    elif args.action == "run":
        stage_run(args.stage)
    else:
        nominee, plan, hosts = stage_check(args.stage)
        print(
            json.dumps(
                {
                    "stage": args.stage,
                    "ready": True,
                    "hosts": hosts,
                    "limits": plan["pilots"][args.stage]["limits"],
                    "nominee": nominee["identity"],
                    "execution": "not started",
                }
            )
        )


if __name__ == "__main__":
    main()
