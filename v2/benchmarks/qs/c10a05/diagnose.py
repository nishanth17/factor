"""Instrument caller boundaries separately from accepted timing samples."""

import json
import time
from unittest.mock import patch

from ....qs.pipeline import QSJob
from ..c9 import c9_growth_study as c9
from ..c10.study import machine_window


def trace(fixture, config, seconds):
    events = []
    original = QSJob._solve
    started = time.monotonic()

    def solve(job, final):
        before = job.stats["filter_calls"], job.stats["solve_calls"]
        count = len(job.collector._rows)
        prior = job.last_count
        divisor = original(job, final)
        after = job.stats["filter_calls"], job.stats["solve_calls"]
        if after != before:
            events.append(
                dict(
                    seconds=time.monotonic() - started,
                    final=final,
                    count=count,
                    growth=count - prior,
                    filter=after[0] - before[0],
                    solve=after[1] - before[1],
                    position=job.next_position,
                    divisor=divisor,
                    matrix=dict(job.stats.get("last_filter", {})),
                    solver=dict(job.stats.get("last_solver", {})),
                )
            )
        return divisor

    with patch.object(QSJob, "_solve", solve):
        result = c9.run_one(fixture, 7, config, "diagnostic", seconds)
    return dict(result=result, events=events, instrumented=True)


def main():
    c9.runtime()
    inputs = c9.INPUTS
    corpus = json.loads((inputs / "corpora/c10_fresh_corpus.json").read_text())
    configs = json.loads((inputs / "controls/c10_training.json").read_text())
    output = c9.ROOT / "v2/benchmarks/results/c10a05/feasibility"
    output.mkdir(parents=True, exist_ok=True)
    with machine_window():
        for digits in (40, 50):
            path = output / f"baseline-{digits}.json"
            if path.exists():
                continue
            fixture = next(
                f for f in corpus["fixtures"] if f["digits"] == digits
            )
            config = c9.decode(configs["configurations"][str(digits)]["plan"])
            value = trace(fixture, config, 240)
            path.write_text(json.dumps(value, indent=2) + "\n")
            print(
                digits,
                value["result"]["seconds"],
                len(value["events"]),
                flush=True,
            )


if __name__ == "__main__":
    main()
