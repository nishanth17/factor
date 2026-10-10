"""Separate construction/storage costs for frozen combined C6/B4 choices."""

import argparse
import hashlib
import json
import resource
from pathlib import Path

from . import c6_b4_study, c6_fast_study, c6_study
from .build_c6_inputs import require_runtime
from .c6_fast_costs import measured


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    require_runtime()
    c6_b4_study.install()
    settings = c6_fast_study.protocol()
    selected = json.loads(c6_b4_study.SELECTION.read_text())
    report = dict(protocol=settings, selection=selected, costs={}, storage={})
    with c6_study.performance_window("combined C6/B4 construction costs"):
        backend = c6_fast_study.p41_campaign.PYTHON_BACKEND
        for arm in selected["candidates"]["int"]:
            if arm == "old-ladder":
                continue

            def construct():
                program = c6_b4_study.construct(arm, backend)
                if len(program.entries) != 303:
                    raise AssertionError("wrong stage schedule")
                return program

            report["costs"][arm] = measured(construct)
            program = construct()
            actions = [row[2] for row in program.entries]
            report["storage"][arm] = dict(
                registers=max(action.record.slots for action in actions),
                bytecode=sum(len(action.record.code) for action in actions),
                masks=sum(len(action.masks) for action in actions),
                fused_pairs=sum(
                    row[0] == 2
                    for action in actions
                    for row in getattr(action, "plan", ())
                ),
                max_plan_entries=max(
                    len(getattr(action, "plan", ())) for action in actions
                ),
            )
        report["maxrss_bytes_macos"] = resource.getrusage(
            resource.RUSAGE_SELF
        ).ru_maxrss
        report["source_sha256"] = hashlib.sha256(
            Path(__file__).read_bytes()
        ).hexdigest()
        with args.output.open("x") as output:
            json.dump(report, output, indent=2)
            output.write("\n")
        c6_fast_study.protocol()


if __name__ == "__main__":
    main()
