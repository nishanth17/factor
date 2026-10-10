"""Long-batch extension of inconclusive PyPy portfolio measurements."""

import argparse
import hashlib
import json
import math
import random
import statistics
import time
from pathlib import Path

from ...support.paths import (
    BENCHMARK_ROOT,
    source_path,
)
from .a6_pm1 import assert_quiet, performance_window, require_runtime
from .a6_pm1_followup import relative_iqr
from .a6_production import confidence, identity
from .a6_production_portfolio import run

PROTOCOL = BENCHMARK_ROOT / "inputs/controls/a6_long_protocol.json"


def measure(functions, confirm):
    samples = {arm: [] for arm in functions}
    warmups, loops = {}, {}
    generator = random.Random(2026100929)
    for arm, function in functions.items():
        start, cpu, count = time.monotonic(), time.process_time(), 0
        while time.monotonic() - start < 10:
            function()
            count += 1
        elapsed = time.process_time() - cpu
        warmups[arm] = dict(
            seconds=time.monotonic() - start, iterations=count, cpu=elapsed
        )
        loops[arm] = max(1, math.ceil(count / elapsed))
    target = 9
    for index in range(63):
        order = list(functions)
        generator.shuffle(order)
        for arm in order:
            assert_quiet()
            start = time.process_time()
            for _ in range(loops[arm]):
                functions[arm]()
            samples[arm].append((time.process_time() - start) / loops[arm])
            assert_quiet()
        if index + 1 < target:
            continue
        unstable = any(relative_iqr(v) > 0.05 for v in samples.values())
        lower = [
            confidence(samples["control"], values)["ci95"][0]
            for arm, values in samples.items()
            if arm not in ("control", "no_pm1")
        ]
        unresolved = any(v <= 0 for v in lower) if confirm else max(lower) <= 0
        if target < 63 and (unstable or unresolved):
            target = 27 if target < 27 else 63
        else:
            break
    return {
        arm: dict(
            cpu_samples=values,
            median_cpu=statistics.median(values),
            relative_iqr=relative_iqr(values),
            loops=loops[arm],
            warmup=warmups[arm],
        )
        for arm, values in samples.items()
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("phase", choices=("screen", "confirmation"))
    parser.add_argument("--arm", default="recurrence64")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    require_runtime()
    before = identity()
    protocol = json.loads(source_path(PROTOCOL).read_text())
    arms = (
        protocol["screen_arms"]
        if args.phase == "screen"
        else ["no_pm1", "control", args.arm]
    )
    functions = {
        arm: lambda arm=arm: run(
            arm, "python-int", 0 if args.phase == "screen" else 12
        )
        for arm in arms
    }
    reference = functions["control"]()
    for arm, function in functions.items():
        if arm != "no_pm1" and [dict(v, work=0) for v in function()] != [
            dict(v, work=0) for v in reference
        ]:
            raise AssertionError("portfolio outcomes changed")
    with performance_window():
        records = measure(functions, args.phase == "confirmation")
    for arm in records:
        records[arm]["outcome"] = functions[arm]()
    if identity() != before:
        raise ValueError("source changed during long-batch capture")
    comparisons = {
        arm: confidence(
            records["control"]["cpu_samples"], records[arm]["cpu_samples"]
        )
        for arm in arms
        if arm not in ("control", "no_pm1")
    }
    result = dict(
        phase=args.phase,
        source=before,
        arms=records,
        comparisons=comparisons,
        runner_sha256=hashlib.sha256(
            source_path(Path(__file__)).read_bytes()
        ).hexdigest(),
        protocol_sha256=hashlib.sha256(
            source_path(PROTOCOL).read_bytes()
        ).hexdigest(),
    )
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(comparisons, indent=2))


if __name__ == "__main__":
    main()
