"""Prospective A1b gap repair: only the six unvisited 55-digit cells."""

import argparse
import builtins
import fcntl
import json
import os
import statistics
import sys
import time
from math import exp, log
from unittest.mock import patch

from . import a1b
from .c10a04.accounting import Ledger, digest, output_lock, receipt

ROOT = a1b.ROOT
CONTROL = ROOT / "v2/benchmarks/inputs/controls/a1b_continuation.json"
FREEZE = ROOT / "v2/benchmarks/inputs/controls/a1b_continuation_source.json"
OUTPUT = ROOT / "v2/benchmarks/results/a1b_continuation_v1"


def geometric_mean(values):
    values = tuple(values)
    return exp(sum(log(value) for value in values) / len(values))


def can_warm(completed, remaining_ends, maximum):
    """Reserve only endpoints that this process has not yet attempted."""
    return completed + remaining_ends < maximum


def deep_bytes(value, seen=None):
    """Bound the small JSON-shaped helper graph without double counting."""
    if seen is None:
        seen = set()
    if id(value) in seen:
        return 0
    seen.add(id(value))
    # PyPy deliberately cannot report object sizes. A 1024-byte node
    # allowance plus character/bit/container storage bounds this small,
    # finite JSON-shaped graph conservatively rather than returning zero.
    size = 1024
    if isinstance(value, str):
        size += 4 * len(value)
    elif isinstance(value, int):
        size += value.bit_length() // 8 + 1
    if isinstance(value, dict):
        size += sum(
            deep_bytes(k, seen) + deep_bytes(v, seen) for k, v in value.items()
        )
    elif isinstance(value, (list, tuple)):
        size += 16 * len(value)
        size += sum(deep_bytes(v, seen) for v in value)
    return size


def bounded_sorted(entries):
    """Refuse validation scratch before sorting more than 8192 atoms."""
    if len(entries) > 8192:
        raise MemoryError("A1b validation-sort helper refusal")
    return builtins.sorted(entries)


def helper_guard(values, row):
    """Reserve globals and two serialization buffers inside the 8 MiB grant."""
    encoded = len(json.dumps(row).encode())
    used = 2 * 2**20 + deep_bytes(values) + 2 * encoded
    if encoded > 2**20 or used > 8 * 2**20:
        raise MemoryError("A1b helper observation/serialization refusal")
    return used


def load(check_evidence=True):
    frozen = json.loads(FREEZE.read_text())
    for name, expected in frozen["files"].items():
        if digest(ROOT / name) != expected:
            raise ValueError("continuation source changed: " + name)
    parent = a1b.load()
    protocol = json.loads(CONTROL.read_text())
    for name, expected in (
        protocol["parent_evidence"].items() if check_evidence else ()
    ):
        if digest(a1b.OUTPUT / name) != expected:
            raise ValueError("parent evidence changed: " + name)
    if protocol["schedule"] != [
        c for c in parent["schedule"] if c["digits"] == 55
    ]:
        raise ValueError("continuation may visit only original 55 cells")
    for cell in protocol["schedule"]:
        if (a1b.OUTPUT / f"cell-{cell['cell']:02d}.json").exists():
            raise ValueError("parent 55 cell is already visited")
    return protocol, parent


def call_prefix(
    protocol, ledger, fixture, config, blocks, seconds, work, retained=()
):
    # Only benchmark validation gets a scratch bound; native SIQS is untouched.
    with patch.object(a1b, "sorted", bounded_sorted, create=True):
        row = a1b.run_one(fixture, 7, config, blocks, seconds, work)
    row["helper_owned_bytes"] = helper_guard(
        [protocol, ledger.data, row, retained], row
    )
    return row


def measure(protocol, parent, ledger):
    """One fresh process, one band, durable starts and all attempted ends."""
    process = f"{os.getpid()}-{time.time_ns()}"
    rules = parent["warm"]
    fixtures = [f for f in parent["fixtures"] if f["digits"] == 55]
    configs = {
        arm: a1b._decode_siqs_config(values)
        for arm, values in parent["configurations"]["55"].items()
    }
    starts = {}
    warm_count = sum(r["phase"] == "warm" for r in ledger.data["calls"])
    receipt(
        OUTPUT / f"process-{process}.json",
        dict(
            process=process,
            interpreter=sys.executable,
            version=sys.version,
            freeze=digest(FREEZE),
            scope="55 only; fresh process readiness",
        ),
    )
    try:
        for arm, config in configs.items():
            rows = []
            while not a1b.ready(rows, rules):
                if not can_warm(warm_count, 2, 10):
                    raise RuntimeError("finite continuation readiness stop")
                identity = dict(
                    process=process, arm=arm, kind="warm", index=warm_count
                )
                row = ledger.call(
                    "warm",
                    f"warm-{process}-{warm_count}",
                    5,
                    10**10,
                    identity,
                    lambda: call_prefix(
                        protocol,
                        ledger,
                        fixtures[0],
                        config,
                        1056,
                        4,
                        10**10,
                        retained=[parent, rows],
                    ),
                )
                rows.append(row)
                warm_count += 1
                if row["reason"] != "paused":
                    raise RuntimeError(
                        "native warm outcome recorded; "
                        "physical prefix incomplete"
                    )
            starts[arm] = {
                key: statistics.median(r[key] for r in rows[-3:])
                for key in ("seconds", "cpu_seconds")
            }
            receipt(
                OUTPUT / f"start-{process}-{arm}.json",
                dict(
                    process=process,
                    arm=arm,
                    calls=len(rows),
                    start=starts[arm],
                ),
            )
        for cell in protocol["schedule"]:
            name = f"cell-{cell['cell']:02d}"
            if any(
                r["receipt"] == name + ".json" for r in ledger.data["calls"]
            ):
                continue
            fixture = next(f for f in fixtures if f["id"] == cell["id"])
            identity = dict(
                sample=cell["cell"],
                fixture=fixture["id"],
                n=fixture["n"],
                digits=55,
                seed=7,
                arm=cell["arm"],
                cap=20,
                process=process,
            )
            row = ledger.call(
                "measured",
                name,
                25,
                2 * 10**10,
                identity,
                lambda: call_prefix(
                    protocol,
                    ledger,
                    fixture,
                    configs[cell["arm"]],
                    2146,
                    20,
                    2 * 10**10,
                    retained=[parent, rows],
                ),
            )
            print(
                cell["cell"],
                cell["arm"],
                row["reason"],
                round(row["seconds"], 6),
                flush=True,
            )
            if row["reason"] != "paused":
                raise RuntimeError(
                    "native measured outcome recorded; "
                    "physical prefix incomplete"
                )
    finally:
        remaining_ends = 2
        for arm in configs:
            name = f"end-{process}-{arm}"
            receipt(
                OUTPUT / (name + "-attempt.json"),
                dict(process=process, arm=arm, start_available=arm in starts),
            )
            if arm not in starts:
                receipt(
                    OUTPUT / (name + "-missing.json"),
                    dict(reason="no validated start in this process"),
                )
                remaining_ends -= 1
                continue
            try:
                if warm_count >= 10:
                    raise RuntimeError("finite endpoint allowance")
                row = ledger.call(
                    "warm",
                    name,
                    5,
                    10**10,
                    dict(process=process, arm=arm, kind="end"),
                    lambda: call_prefix(
                        protocol,
                        ledger,
                        fixtures[0],
                        configs[arm],
                        1056,
                        4,
                        10**10,
                        retained=[parent, rows],
                    ),
                )
                warm_count += 1
                ratios = {
                    key: row[key] / starts[arm][key] for key in starts[arm]
                }
                receipt(
                    OUTPUT / (name + "-stability.json"),
                    dict(
                        process=process,
                        arm=arm,
                        ratios=ratios,
                        passed=row["reason"] == "paused"
                        and all(0.8 <= v <= 1.2 for v in ratios.values()),
                    ),
                )
            except BaseException as error:
                receipt(
                    OUTPUT / (name + "-failure.json"),
                    dict(error=type(error).__name__ + ": " + str(error)),
                )
                raise
            finally:
                remaining_ends -= 1


def summarize_band(parent, digits, source):
    """Use the original input medians and advancement rule without retuning."""
    inputs = []
    for fixture in (f for f in parent["fixtures"] if f["digits"] == digits):
        values = {key: [] for key in ("wall", "cpu")}
        for seed in (7, 29, 47) if digits == 60 else (7,):
            rows = {}
            for cell in parent["schedule"]:
                if cell["id"] != fixture["id"] or cell["seed"] != seed:
                    continue
                path = source / f"cell-{cell['cell']:02d}.json"
                if not path.exists():
                    return dict(
                        qualifies=False, reason="missing measured cell"
                    )
                rows[cell["arm"]] = json.loads(path.read_text())["result"]
            if any(r["reason"] != "paused" for r in rows.values()):
                return dict(qualifies=False, reason="incomplete native prefix")
            control, fixed = rows["control"], rows["fixed"]
            if control["physical"] != fixed["physical"]:
                raise ArithmeticError("physical prefix mismatch")
            for key, clock in (("wall", "seconds"), ("cpu", "cpu_seconds")):
                values[key].append(fixed[clock] / control[clock])
        inputs.append(
            dict(
                id=fixture["id"],
                pair_ratios=values,
                wall=statistics.median(values["wall"]),
                cpu=statistics.median(values["cpu"]),
            )
        )
    ends = [
        json.loads(p.read_text())
        for p in source.glob(
            "end-*-stability.json"
            if digits == 55
            else f"end-*-{digits}-*-stability.json"
        )
    ]
    starts = list(
        source.glob(
            "start-*.json" if digits == 55 else f"start-*-{digits}-*.json"
        )
    )
    stability = len(ends) == len(starts) == 2 and all(
        e["passed"] for e in ends
    )
    wall = geometric_mean(r["wall"] for r in inputs)
    cpu = geometric_mean(r["cpu"] for r in inputs)
    wins = sum(r["wall"] < 1 for r in inputs)
    return dict(
        inputs=inputs,
        independent_inputs=3,
        geometric_wall=wall,
        geometric_cpu=cpu,
        wall_wins=wins,
        stability=stability,
        qualifies=stability and wins >= 2 and wall < 1 and cpu <= 1,
    )


def analyze(protocol, parent, ledger):
    result = dict(identity=protocol["identity"], work=0, bands={})
    for digits in (50, 55, 60):
        source = OUTPUT if digits == 55 else a1b.OUTPUT
        result["bands"][str(digits)] = summarize_band(parent, digits, source)
    result["decision"] = (
        "pilot design needed"
        if any(b["qualifies"] for b in result["bands"].values())
        else "retain"
    )
    result["parent_charges_preserved"] = 1240
    result["continuation_charges_before_analysis"] = {
        key: ledger.used(None, key) for key in ("wall", "cpu", "work")
    }
    result["source_sha256"] = digest(FREEZE)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("verify", "run", "analyze"))
    args = parser.parse_args()
    protocol, parent = load(check_evidence=args.command != "verify")
    if args.command == "verify":
        print(protocol["identity"], digest(FREEZE))
        return
    with (
        output_lock(OUTPUT),
        open("/private/tmp/factor-performance.lock", "a+") as lock,
    ):
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        ledger = Ledger(OUTPUT, digest(FREEZE), protocol["limits"])
        if args.command == "run":
            ledger.call(
                "legacy_analysis",
                "legacy-analysis-reconciliation",
                10,
                0,
                dict(kind="preserve and debit prior omission"),
                lambda: dict(
                    work=0,
                    known_legacy_wall=0.024204791989,
                    legacy_cpu=None,
                    other_legacy_analysis_actual=None,
                    note="Conservative repair charge; missing "
                    "actual clocks remain missing. Parent "
                    "failed60s charge is unchanged.",
                ),
            )
            ledger.call(
                "orchestration",
                "orchestration-reservation",
                5,
                0,
                dict(kind="continuation overhead"),
                lambda: dict(work=0),
            )
            measure(protocol, parent, ledger)
        else:
            result = ledger.call(
                "analysis",
                "analysis",
                15,
                0,
                dict(kind="unchanged all-band screen rule"),
                lambda: analyze(protocol, parent, ledger),
            )
            print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
