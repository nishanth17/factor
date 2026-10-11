"""Censored default-grant and charged-resume diagnostics; no speed claims."""

import argparse
import importlib
import json
import time
from pathlib import Path

from .... import portfolio
from .c3_study import (
    SELECTED,
    band_of,
    baseline,
    config_for,
    machine_window,
    run_one,
    runtime,
    save,
    training,
    validate,
    verify_protocol,
)


def diagnose(output):
    runtime()
    verify_protocol()
    selected = json.loads(SELECTED.read_text())["bands"]
    control = baseline()
    rows, resumes = [], []
    started = time.monotonic()
    with machine_window():
        for fixture in training():
            arms = list(
                dict.fromkeys(
                    ["control", selected[str(band_of(fixture["n"]))]]
                )
            )
            for arm in arms:
                if time.monotonic() - started > 600:
                    raise RuntimeError("diagnostic allowance exhausted")
                rows.append(
                    run_one(fixture, 7, arm, control, regime="default-probe")
                )
                if fixture["kind"] != "balanced":
                    continue
                module = control if arm == "control" else portfolio
                budget_module = importlib.import_module(
                    module.__package__ + ".execution.budget"
                )
                config = config_for(fixture["n"], arm, module)
                first = module.factorize_bounded(
                    fixture["n"],
                    seed=7,
                    config=config,
                    budget=budget_module.Budget(
                        work_limit=2_000_000, seconds=30, cpu_seconds=30
                    ),
                )
                validate(first.result, fixture)
                restored = module.factorize_bounded(
                    fixture["n"],
                    config=config,
                    checkpoint=first.checkpoint,
                    budget=budget_module.Budget(
                        work_limit=10**13, seconds=30, cpu_seconds=30
                    ),
                )
                validate(restored.result, fixture)
                if (
                    restored.work_used < first.work_used
                    or restored.wall_seconds < first.wall_seconds
                    or restored.cpu_seconds < first.cpu_seconds
                ):
                    raise AssertionError("resume lost prior expenditure")
                resumes.append(
                    dict(
                        id=fixture["id"],
                        arm=arm,
                        first_reason=first.reason,
                        first_work=first.work_used,
                        first_wall=first.wall_seconds,
                        first_cpu=first.cpu_seconds,
                        first_stage=(
                            first.checkpoint["payload"]["state"]["current"]
                            or {}
                        ).get("stage"),
                        cumulative_work=restored.work_used,
                        cumulative_wall=restored.wall_seconds,
                        cumulative_cpu=restored.cpu_seconds,
                        complete=restored.result.complete,
                        reason=restored.reason,
                        events=restored.events,
                    )
                )
    save(
        output,
        dict(
            rows=rows,
            resumes=resumes,
            seconds=time.monotonic() - started,
            limitation="Single-run diagnostics; no performance claim.",
        ),
    )
    print("default probes", len(rows), "resumes", len(resumes), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    diagnose(parser.parse_args().output)


if __name__ == "__main__":
    main()
