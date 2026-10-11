"""At most two frozen C9 joint configurations after individual attribution."""

import argparse
import hashlib
import json
import random
import time
from dataclasses import asdict, replace
from pathlib import Path

from ..b1.b1_calibration import save
from . import c9_study as study

CONTROL = study.INPUTS / "controls/c9_combinations.json"


def hashes():
    result = study.hashes()
    path = Path(__file__)
    result[str(path.relative_to(study.ROOT))] = hashlib.sha256(
        path.read_bytes()
    ).hexdigest()
    return result


def freeze(training):
    initial = json.loads(study.CONTROL.read_text())
    if json.loads((training / "manifest.json").read_text()) != initial:
        raise ValueError("training manifest differs from initial freeze")
    rows = [
        json.loads(p.read_text()) for p in sorted(training.glob("run-*.json"))
    ]
    if len(rows) != 10 or not any(r["complete"] for r in rows):
        raise ValueError("C9 individual-arm advancement gate did not pass")
    base = study.decode(initial["configurations"]["slp"])
    half = study.arm_config(base, "dlp_half")
    aligned = replace(
        half,
        collector=replace(
            half.collector, residual_bound=half.collector.large_prime_bound
        ),
    )
    save(
        CONTROL,
        dict(
            schema=1,
            source=hashes(),
            fixtures=initial["fixtures"],
            seeds=initial["seeds"],
            configurations={
                "half": asdict(half),
                "half_aligned": asdict(aligned),
            },
            individual_results=[
                {
                    k: row[k]
                    for k in (
                        "id",
                        "arm",
                        "complete",
                        "seconds",
                        "capped_seconds",
                        "reason",
                    )
                }
                for row in rows
            ],
            gate="Individual controls completed; capacity alone refused at "
            "256 MiB. Test historical half-base DLP and its aligned endpoint. "
            "This is a joint bundle; C9 does not claim isolated credit.",
            ownership="No active C10 chat was found. "
            "No shared collector source edits. "
            "Consume the historical C1 half-base configuration. "
            "C3 was informed. "
            "C10 retains broader calibration.",
            seconds=240,
            prefixes=[60, 120, 240],
            memory_bytes=256 * 2**20,
        ),
    )


def run(output):
    control = json.loads(CONTROL.read_text())
    if control["source"] != hashes():
        raise ValueError("joint configurations changed after freeze")
    if len(control["configurations"]) > 2 or len(control["fixtures"]) != 2:
        raise ValueError("C9 joint-grid cap")
    save(output / "manifest.json", control)
    start, cpu = time.monotonic(), time.process_time()
    for index, fixture in enumerate(control["fixtures"]):
        arms = list(control["configurations"])
        random.Random(9910 + index).shuffle(arms)
        for arm in arms:
            row = study.run_one(
                fixture,
                control["seeds"][index],
                study.decode(control["configurations"][arm]),
                arm,
                240,
                prefixes=True,
            )
            save(output / f"run-{index}-{arm}.json", row)
            print(
                json.dumps(
                    {
                        k: row[k]
                        for k in (
                            "id",
                            "arm",
                            "complete",
                            "reason",
                            "seconds",
                        )
                    }
                ),
                flush=True,
            )
    save(
        output / "finished.json",
        dict(
            seconds=time.monotonic() - start,
            cpu_seconds=time.process_time() - cpu,
        ),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("freeze", "run"))
    parser.add_argument("--output", type=Path)
    parser.add_argument("--training", type=Path)
    args = parser.parse_args()
    study.runtime()
    with study.machine_window():
        if args.mode == "freeze":
            if args.training is None:
                parser.error("freeze requires training directory")
            freeze(args.training)
        else:
            if args.output is None:
                parser.error("run requires output directory")
            run(args.output)


if __name__ == "__main__":
    main()
