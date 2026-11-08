"""Frozen existing fixed/list screen against integrated SIQS collection."""

import argparse
import fcntl
import gc
import hashlib
import json
import os
import platform
import resource
import statistics
import sys
import time
from math import exp, log
from pathlib import Path

from v2.benchmarks.infrastructure.resource_control import (
    RecordingBudget as Budget,
)

from ...common.utils import classify_prime
from ...portfolio import _decode_siqs_config
from ...qs import SIQSJob, verify_atomic
from .c10a04.accounting import Ledger, digest, output_lock, receipt

ROOT = Path(__file__).resolve().parents[3]
CONTROL = ROOT / "v2/benchmarks/inputs/controls/a1b_screen.json"
FREEZE = ROOT / "v2/benchmarks/inputs/controls/a1b_source.json"
OUTPUT = ROOT / "v2/benchmarks/results/a1b_screen_v1"


def load():
    """Refuse executable, source, interpreter or protocol substitutions."""
    frozen = json.loads(FREEZE.read_text())
    for name, expected in frozen["files"].items():
        if digest(ROOT / name) != expected:
            raise ValueError("A1b frozen source changed: " + name)
    if (
        platform.python_implementation() != "PyPy"
        or sys.version_info[:2] != (3, 11)
        or digest(Path(sys.executable).resolve()) != frozen["interpreter"]
    ):
        raise ValueError("A1b requires the frozen PyPy 3.11 interpreter")
    return json.loads(CONTROL.read_text())


def limits(protocol):
    """Keep every phase separate; unused funds never finance retries."""
    result = {
        key: dict(wall=value, cpu=value, work=0)
        for key, value in protocol["accounts"].items()
    }
    for key, value in protocol["work_limits"].items():
        result[key]["work"] = value
    result["total"]["work"] = sum(protocol["work_limits"].values())
    result["disk_bytes"] = protocol["memory"]["output"]
    return result


def run_one(fixture, seed, config, blocks, seconds, work):
    """Time native setup-to-prefix work and exact output validation."""
    gc.collect()
    wall, cpu = time.perf_counter(), time.process_time()
    budget = Budget(work_limit=work, seconds=seconds, cpu_seconds=seconds)
    job = SIQSJob(fixture["n"], seed=seed, config=config, budget=budget)
    result = job.run(max_blocks=blocks)
    native_wall = time.perf_counter() - wall
    native_cpu = time.process_time() - cpu
    n = fixture["n"]
    if result.divisor is not None:
        if not 1 < result.divisor < n or result.divisor * result.cofactor != n:
            raise ArithmeticError("invalid A1b divisor/reconstruction")
    elif result.cofactor != n:
        raise ArithmeticError("invalid unresolved reconstruction")
    certainty = [classify_prime(result.cofactor).value]
    if result.divisor:
        certainty.append(classify_prime(result.divisor).value)
    collector = job.engine.collector if job.engine else None
    state = hashlib.sha256()
    queue_peak = 0
    if collector:
        # Stream validation and hashes: do not retain a second relation store.
        for key, atom in sorted(collector._atoms.items()):
            verify_atomic(
                atom,
                job.base,
                residual_bound=collector.config.residual_bound,
                large_prime_bound=collector.config.large_prime_bound,
                large_product_bound=collector.config.large_product_bound,
                budget=budget,
            )
            state.update(
                repr(
                    (key, atom.sign, atom.exponents, atom.residual_product)
                ).encode()
            )
        queue_peak = getattr(collector, "_sparse_queue_peak", 0)
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform != "darwin":
        rss *= 1024
    if rss > 3 * 2**30:
        raise MemoryError("A1b RSS ceiling")
    if result.stats.get("workspace_bytes", 0) > config.memory_limit:
        raise MemoryError("A1b native owned ceiling")
    physical = dict(
        blocks=result.stats["blocks"],
        polynomials=result.stats.get("polynomials"),
        family=job.family_index,
        polynomial=collector.polynomial.identity if collector else None,
        position=job.engine.next_position if job.engine else None,
        scanned=result.stats.get("engine", {}).get("collected_positions"),
    )
    return dict(
        seconds=time.perf_counter() - wall,
        cpu_seconds=time.process_time() - cpu,
        native_wall=native_wall,
        native_cpu=native_cpu,
        validation_included=True,
        work=budget.used,
        reason=result.reason,
        divisor=result.divisor,
        cofactor=result.cofactor,
        certainty=certainty,
        complete=bool(result.divisor)
        and all(label != "composite" for label in certainty),
        physical=physical,
        stats=result.stats,
        state_sha256=state.hexdigest(),
        queue_peak=queue_peak,
        rss_bytes=rss,
        process=os.getpid(),
    )


def ready(rows, rules):
    """Interpret readiness prospectively, including ordinary native stops."""
    if len(rows) < rules["min_calls"]:
        return False
    if not all(row["reason"] == "paused" for row in rows):
        return False
    return all(
        sum(row[key] for row in rows) >= rules["min_wall_cpu"]
        and (
            max(row[key] for row in rows[-3:])
            - min(row[key] for row in rows[-3:])
        )
        / statistics.median(row[key] for row in rows[-3:])
        <= rules["last_three_range_over_median"]
        for key in ("seconds", "cpu_seconds")
    )


def measure(protocol, ledger):
    """Persist readiness starts and attempted ends; visit each cell once."""
    warm_rules = protocol["warm"]
    warm_calls = sum(r["phase"] == "warm" for r in ledger.data["calls"])
    process = str(os.getpid()) + "-" + str(time.time_ns())
    for digits in (60, 50, 55):
        band = [f for f in protocol["fixtures"] if f["digits"] == digits]
        configs = {
            arm: _decode_siqs_config(values)
            for arm, values in protocol["configurations"][str(digits)].items()
        }
        starts = {}
        try:
            for arm, config in configs.items():
                rows = []
                # Reserve the six defined end observations before readiness.
                while not ready(rows, warm_rules):
                    if warm_calls >= warm_rules["max_calls"] - 6:
                        raise RuntimeError("A1b finite readiness allowance")
                    identity = dict(
                        process=process,
                        digits=digits,
                        arm=arm,
                        kind="warm",
                        index=warm_calls,
                    )
                    row = ledger.call(
                        "warm",
                        f"warm-{process}-{warm_calls}",
                        10,
                        10**10,
                        identity,
                        lambda: run_one(
                            band[0], 7, config, warm_rules["blocks"], 5, 10**10
                        ),
                    )
                    rows.append(row)
                    warm_calls += 1
                    if row["reason"] != "paused":
                        raise RuntimeError("A1b warm prefix not completed")
                starts[arm] = {
                    key: statistics.median(row[key] for row in rows[-3:])
                    for key in ("seconds", "cpu_seconds")
                }
                receipt(
                    OUTPUT / f"start-{process}-{digits}-{arm}.json",
                    dict(
                        process=process,
                        digits=digits,
                        arm=arm,
                        start=starts[arm],
                        validated_calls=len(rows),
                    ),
                )
            for cell in protocol["schedule"]:
                if cell["digits"] != digits:
                    continue
                fixture = next(f for f in band if f["id"] == cell["id"])
                identity = dict(
                    sample=cell["cell"],
                    fixture=fixture["id"],
                    n=fixture["n"],
                    digits=digits,
                    arm=cell["arm"],
                    seed=cell["seed"],
                    cap=cell["seconds"],
                    process=process,
                )
                name = f"cell-{cell['cell']:02d}"
                # A prior cell remains visited, including lost/censored work.
                if any(
                    r["receipt"] == name + ".json"
                    for r in ledger.data["calls"]
                ):
                    continue
                row = ledger.call(
                    "measured",
                    name,
                    cell["seconds"] + 5,
                    cell["work"],
                    identity,
                    lambda: run_one(
                        fixture,
                        cell["seed"],
                        configs[cell["arm"]],
                        cell["blocks"],
                        cell["seconds"],
                        cell["work"],
                    ),
                )
                print(
                    digits,
                    cell["cell"],
                    cell["arm"],
                    row["reason"],
                    round(row["seconds"], 4),
                    flush=True,
                )
                if row["reason"] != "paused":
                    raise RuntimeError("A1b prefix incomplete; audit required")
        finally:
            # An endpoint failure is durable and cannot be replaced later.
            for arm, start in starts.items():
                name = f"end-{process}-{digits}-{arm}"
                attempted = dict(process=process, digits=digits, arm=arm)
                receipt(OUTPUT / (name + "-attempt.json"), attempted)
                try:
                    if warm_calls >= warm_rules["max_calls"]:
                        raise RuntimeError("A1b end allowance exhausted")
                    row = ledger.call(
                        "warm",
                        name,
                        10,
                        10**10,
                        attempted,
                        lambda: run_one(
                            band[0],
                            7,
                            configs[arm],
                            warm_rules["blocks"],
                            5,
                            10**10,
                        ),
                    )
                    warm_calls += 1
                    ratios = {key: row[key] / start[key] for key in start}
                    receipt(
                        OUTPUT / (name + "-stability.json"),
                        dict(
                            ratios=ratios,
                            passed=row["reason"] == "paused"
                            and all(
                                0.8 <= ratio <= 1.2
                                for ratio in ratios.values()
                            ),
                        ),
                    )
                except BaseException as error:
                    receipt(
                        OUTPUT / (name + "-failure.json"),
                        dict(error=type(error).__name__ + ": " + str(error)),
                    )
                    raise


def geometric_mean(values):
    return exp(sum(map(log, values)) / len(values))


def analyze(protocol, ledger):
    """Use every input and the frozen band rule; censoring cannot qualify."""
    bands = {}
    for digits in (50, 55, 60):
        pairs = []
        complete = True
        for fixture in (
            f for f in protocol["fixtures"] if f["digits"] == digits
        ):
            ratios = {key: [] for key in ("wall", "cpu")}
            for seed in (7, 29, 47) if digits == 60 else (7,):
                rows = {}
                for cell in protocol["schedule"]:
                    if cell["id"] == fixture["id"] and cell["seed"] == seed:
                        path = OUTPUT / f"cell-{cell['cell']:02d}.json"
                        if path.exists():
                            rows[cell["arm"]] = json.loads(path.read_text())[
                                "result"
                            ]
                if set(rows) != {"control", "fixed"} or any(
                    row["reason"] != "paused" for row in rows.values()
                ):
                    complete = False
                    continue
                if rows["control"]["physical"] != rows["fixed"]["physical"]:
                    raise ArithmeticError("A1b physical prefix changed")
                for key, clock in (
                    ("wall", "seconds"),
                    ("cpu", "cpu_seconds"),
                ):
                    ratios[key].append(
                        rows["fixed"][clock] / rows["control"][clock]
                    )
            pairs.append(
                dict(
                    id=fixture["id"],
                    ratios=ratios,
                    wall=statistics.median(ratios["wall"])
                    if ratios["wall"]
                    else None,
                    cpu=statistics.median(ratios["cpu"])
                    if ratios["cpu"]
                    else None,
                )
            )
        sentinels = [
            json.loads(p.read_text())
            for p in OUTPUT.glob(f"end-*-{digits}-*-stability.json")
        ]
        starts = list(OUTPUT.glob(f"start-*-{digits}-*.json"))
        stability = len(starts) == len(sentinels) == 2 and all(
            s["passed"] for s in sentinels
        )
        wall = geometric_mean(p["wall"] for p in pairs) if complete else None
        cpu = geometric_mean(p["cpu"] for p in pairs) if complete else None
        wins = sum(p["wall"] < 1 for p in pairs) if complete else 0
        bands[str(digits)] = dict(
            inputs=pairs,
            independent_inputs=3,
            complete_prefix_pairs=complete,
            geometric_wall=wall,
            geometric_cpu=cpu,
            wall_wins=wins,
            stability=stability,
            qualifies=complete
            and stability
            and wins >= 2
            and wall < 1
            and cpu <= 1,
        )
    return dict(
        bands=bands,
        work=0,
        decision="pilot design needed"
        if any(b["qualifies"] for b in bands.values())
        else "retain",
        meaning="Exploratory prefixes, no complete-factor speed claim",
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("verify", "run", "analyze"))
    args = parser.parse_args()
    protocol = load()
    if args.command == "verify":
        print(protocol["identity"], digest(FREEZE))
        return
    with (
        output_lock(OUTPUT),
        open("/private/tmp/factor-performance.lock", "a+") as lock,
    ):
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        ledger = Ledger(OUTPUT, digest(FREEZE), limits(protocol))
        if args.command == "run":
            ledger.call(
                "diagnostic",
                "prior-diagnostic",
                90,
                5 * 10**11,
                dict(kind="imported diagnostic"),
                lambda: dict(
                    work=5098787305,
                    actual_wall=13.051687584,
                    actual_cpu=13.035224,
                    note="Preserved audit receipts, including v1",
                ),
            )
            ledger.call(
                "orchestration",
                "orchestration-reservation",
                90,
                0,
                dict(kind="finite overhead reservation"),
                lambda: dict(work=0),
            )
            measure(protocol, ledger)
        else:
            result = ledger.call(
                "analysis",
                "analysis",
                60,
                0,
                dict(kind="frozen band analysis"),
                lambda: analyze(protocol, ledger),
            )
            print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
