"""Direct portfolio default selection and independent held-out confirmation."""

import argparse
import hashlib
import json
from pathlib import Path

from .a6_pm1 import performance_window, require_runtime, verify_inputs
from .a6_pm1_followup import measure
from .a6_production import confidence, config, identity

INPUT = (
    Path(__file__).parent
    / "inputs/corpora/phase_two_m15_independent_corpus.json"
)
PROTOCOL = Path(__file__).parent / "inputs/controls/a6_portfolio_protocol.json"


def run(arm, backend, start):
    data = json.loads(INPUT.read_text())
    fixtures = [f for f in data["fixtures"] if len(str(f["n"])) == 20]
    stop = start + 12
    fixtures = fixtures[start:stop]
    module, configuration = config(
        "control" if arm == "no_pm1" else arm,
        backend,
        pm1_attempts=int(arm != "no_pm1"),
    )
    outcomes = []
    for fixture in fixtures:
        from v2.budget import Budget

        outcome = module.factorize_bounded(
            fixture["n"],
            seed=7,
            config=configuration,
            budget=Budget(work_limit=500000, seconds=5, cpu_seconds=5),
        )
        if outcome.result.reconstruct() != fixture["n"]:
            raise AssertionError("portfolio reconstruction failed")
        outcomes.append(
            {
                "complete": outcome.result.complete,
                "factors": [
                    (int(f.value), f.exponent, f.certainty.value)
                    for f in outcome.result.factors
                ],
                "remaining": list(map(int, outcome.result.remaining)),
                "work": outcome.work_used,
            }
        )
    return outcomes


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("phase", choices=("screen", "confirmation"))
    parser.add_argument("--backend", default="python-int")
    parser.add_argument("--arm", default="recurrence64")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    require_runtime()
    verify_inputs()
    protocol = json.loads(PROTOCOL.read_text())
    arms = (
        protocol["screen_arms"]
        if args.phase == "screen"
        else ["no_pm1", "control", args.arm]
    )
    start = 0 if args.phase == "screen" else 12
    source = {
        "main": identity(),
        "runner": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "protocol": hashlib.sha256(PROTOCOL.read_bytes()).hexdigest(),
    }
    functions = {
        arm: lambda arm=arm: run(arm, args.backend, start) for arm in arms
    }
    reference = functions["control"]()
    for arm, function in functions.items():
        result = function()
        if arm != "no_pm1" and [dict(v, work=0) for v in result] != [
            dict(v, work=0) for v in reference
        ]:
            raise AssertionError("portfolio outcomes changed")
    with performance_window():
        measurements = measure(functions, 63)
    records = {
        arm: {**measurements[arm], "outcome": function()}
        for arm, function in functions.items()
    }
    if identity() != source["main"]:
        raise ValueError("source changed during capture")
    result = {
        "source": source,
        "phase": args.phase,
        "start": start,
        "backend": args.backend,
        "arms": records,
        "comparisons": {
            arm: confidence(
                records["control"]["cpu_samples"], records[arm]["cpu_samples"]
            )
            for arm in arms
            if arm not in ("control", "no_pm1")
        },
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result["comparisons"], indent=2))


if __name__ == "__main__":
    main()
