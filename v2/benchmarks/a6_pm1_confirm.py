"""Frozen repeat extension for A6's below-ten-percent chunk candidate."""

import argparse
import hashlib
import json
import platform
import sys
from pathlib import Path

from .a6_pm1 import (
    _measure_arms,
    performance_window,
    require_runtime,
    stage,
    stage_fixtures,
    verify_inputs,
)

PROTOCOL = Path(__file__).parent / "inputs/corpora/a6_pm1_confirm_protocol.json"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--samples", type=int, choices=(27, 63), default=27)
    args = parser.parse_args()
    require_runtime()
    records = []

    with performance_window():
        verify_inputs()
        for fixture in stage_fixtures("confirmation"):
            for b1, b2 in ((2000, 20000), (11000, 100000)):
                functions = {
                    str(chunk): lambda chunk=chunk: stage(
                        fixture["n"], b1, b2, chunk=chunk
                    )
                    for chunk in (16, 64)
                }
                measurements = _measure_arms(functions, samples=args.samples)
                for arm, function in functions.items():
                    outcome = function()
                    if outcome["factor"] is not None:
                        raise AssertionError("confirmation must be nonsplitting")
                    records.append(
                        {
                            "id": fixture["id"],
                            "bounds": [b1, b2],
                            "arm": arm,
                            "outcome": outcome,
                            **measurements[arm],
                        }
                    )

    source_paths = [
        Path(__file__),
        Path(__file__).with_name("a6_pm1.py"),
        PROTOCOL,
        Path(__file__).parents[1] / "stage_jobs.py",
        Path(__file__).parents[1] / "schedules.py",
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(
            {
                "runtime": sys.version,
                "platform": platform.platform(),
                "source_hashes": {
                    path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                    for path in source_paths
                },
                "records": records,
            },
            indent=2,
        )
        + "\n"
    )
    print(args.output)


if __name__ == "__main__":
    main()
