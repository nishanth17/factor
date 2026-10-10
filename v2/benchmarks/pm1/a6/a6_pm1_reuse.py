"""One frozen A6 follow-up separating schedule reuse from LRU churn."""

import argparse
import hashlib
import json
import platform
import sys
from pathlib import Path

from ....execution.schedules import ScheduleCache, SieveContext
from ...support.paths import (
    BENCHMARK_ROOT,
    source_path,
)
from .a6_pm1 import (
    _measure_arms,
    performance_window,
    require_runtime,
    stage,
    stage_fixtures,
    verify_inputs,
)

PROTOCOL = BENCHMARK_ROOT / "inputs/corpora/a6_pm1_reuse_protocol.json"


def run(fixture, arm):
    context = SieveContext(20001, segment_size=256)
    if arm != "stream":
        context = ScheduleCache(
            context,
            cache_bytes=262144,
            max_entries=8 if arm == "entries8" else 64,
        )
    outcomes = [
        stage(fixture["n"], 2000, 20000, base=base, context=context)
        for base in (2, 3, 4)
    ]
    hits = getattr(context, "hits", 0)
    if arm == "entries64" and hits != 80:
        raise AssertionError("complete schedule did not reuse all 40 blocks")
    return {
        "stages": outcomes,
        "hits": hits,
        "retained_bytes": getattr(context, "used_bytes", 0),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    require_runtime()
    records = []
    with performance_window():
        verify_inputs()
        for fixture in stage_fixtures("confirmation"):
            functions = {
                arm: lambda arm=arm: run(fixture, arm)
                for arm in ("stream", "entries8", "entries64")
            }
            expected = run(fixture, "stream")["stages"]
            for arm in functions:
                if run(fixture, arm)["stages"] != expected:
                    raise AssertionError("reused outcomes or work differ")
            measurements = _measure_arms(functions)
            for arm, function in functions.items():
                records.append(
                    {
                        "id": fixture["id"],
                        "arm": arm,
                        "outcome": function(),
                        **measurements[arm],
                    }
                )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(
            {
                "runtime": sys.version,
                "platform": platform.platform(),
                "protocol_sha256": hashlib.sha256(
                    source_path(PROTOCOL).read_bytes()
                ).hexdigest(),
                "runner_sha256": hashlib.sha256(
                    source_path(Path(__file__)).read_bytes()
                ).hexdigest(),
                "parent_runner_sha256": hashlib.sha256(
                    source_path(
                        source_path(BENCHMARK_ROOT / "a6_pm1.py")
                    ).read_bytes()
                ).hexdigest(),
                "records": records,
            },
            indent=2,
        )
        + "\n"
    )
    print(args.output)


if __name__ == "__main__":
    main()
