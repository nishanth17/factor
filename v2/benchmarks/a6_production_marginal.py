"""Supplemental frozen marginal portfolio completion and CPU comparison."""

import hashlib
import json
import sys
from dataclasses import replace
from pathlib import Path
from unittest.mock import patch

from . import a6_production as study
from .a6_pm1 import performance_window, require_runtime, verify_inputs
from .a6_pm1_followup import measure


def call(arm):
    """Match wrapper overhead and change only predeclared p-1 presence."""
    original = study.config

    def configuration(name, backend, **changes):
        module, config = original(name, backend, **changes)
        return module, replace(config, pm1_attempts=int(arm != "no_pm1"))

    with patch.object(study, "config", configuration):
        return study.portfolio_call(
            "control" if arm == "no_pm1" else arm, "python-int"
        )


def main():
    require_runtime()
    verify_inputs()
    identity = study.identity()
    source = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    functions = {
        name: lambda name=name: call(name)
        for name in ("no_pm1", "control", "recurrence64")
    }
    with performance_window():
        measured = measure(functions, 27)
        records = {
            name: {**measured[name], "outcome": function()}
            for name, function in functions.items()
        }
    if study.identity() != identity:
        raise ValueError("source changed during marginal capture")
    control = records["control"]["cpu_samples"]
    result = {
        "source_identity": identity,
        "supplement_source_sha256": source,
        "criteria": "Fixed 27 samples, extend all arms to63 for IQR>15%; "
        "no new arm selection; same 12 inputs, budgets and seed "
        "as primary portfolio confirmation.",
        "arms": records,
        "comparison": study.confidence(
            control, records["recurrence64"]["cpu_samples"]
        ),
    }
    Path(sys.argv[1]).write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result["comparison"], indent=2))


if __name__ == "__main__":
    main()
