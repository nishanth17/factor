"""Host-audited v4 rechecks with the existing live-budget supervisor."""

import argparse
import copy
import json
import os
import sys
from pathlib import Path

from ..c10a04.accounting import digest
from . import host_state, nominee_v4, pilots_live, plan_v4, regression65

MODULE = "v2.benchmarks.qs.c10a01.regression65_recheck"
ORIGINAL_SENTINEL = pilots_live.original.sentinel
ORIGINAL_MEASURED = pilots_live.measured
STAGE = "screen"
FREEZE = None


def build(stage, seconds):
    """A new prospective freeze; never enlarge an existing measured cell."""
    previous = regression65.check()
    configurations = {
        row["arm"]: row["configuration"] for row in previous["rows"]
    }
    arms = (
        ("control", "nominee", "nominee", "control")
        if stage == "screen"
        else ("nominee", "control")
    )
    calls = []
    for index, arm in enumerate(arms):
        calls.append(
            dict(
                id=f"{stage}65-{arm}-{index}",
                arm=arm,
                digits=65,
                tier=2,
                fixture="c10a01_65_balanced",
                phase="measured",
                seconds=seconds,
                reservation_seconds=seconds + 5,
                work=plan_v4.call_work(65, seconds),
                blocks=128 * 33 if stage == "screen" else 0,
                configuration=copy.deepcopy(configurations[arm]),
            )
        )
    files = dict(previous["sources"])
    for path in (Path(__file__), Path(host_state.__file__)):
        files[str(path.relative_to(plan_v4.ROOT))] = digest(path)
    return dict(
        identity=f"c10a01-regression65-host-{stage}-v2",
        sources=files,
        runtime=previous["runtime"],
        calls=calls,
        original_nominee_sha256=digest(nominee_v4.FREEZE),
        expected_power=host_state.power_state(),
        limits=plan_v4.totals(calls, 12, 60, 30),
        decision="No changed arithmetic. Stable matched feasibility only; "
        "keep the original regression and censored diagnostics.",
        native_canary_work_per_sentinel=9 * host_state.CANARY_ITERATIONS,
    )


def source_check(require_receipts=False):
    value = json.loads(FREEZE.read_text())
    if digest(Path(sys.executable).resolve()) != value["runtime"]:
        raise RuntimeError("PyPy binary drift")
    if sys.implementation.name != "pypy" or sys.version_info[:2] != (3, 11):
        raise RuntimeError("PyPy implementing Python 3.11 required")
    for name, expected in value["sources"].items():
        if digest(plan_v4.ROOT / name) != expected:
            raise RuntimeError("source drift: " + name)
    if digest(nominee_v4.FREEZE) != value["original_nominee_sha256"]:
        raise RuntimeError("literal policies changed")
    if require_receipts:
        nominee_v4.verify(require_receipts=True)
    assert value["limits"] == plan_v4.totals(value["calls"], 12, 60, 30)
    return value


def audited_sentinel(row, fixture, modules):
    value = ORIGINAL_SENTINEL(row, fixture, modules)
    value["host"] = dict(
        power=host_state.power_state(), runtime=host_state.runtime_state()
    )
    value["native_canary"] = host_state.native_canary()
    value["work"] += value["native_canary"]["work"]
    return value


def measured(row, fixture, modules, identity, *, blocks=None):
    expected = json.loads(FREEZE.read_text())["expected_power"]
    if host_state.power_state() != expected:
        raise RuntimeError("host power changed before measurement")
    return ORIGINAL_MEASURED(
        row,
        fixture,
        modules,
        identity,
        blocks=row["blocks"] or blocks,
    )


def configure(stage):
    global STAGE, FREEZE

    STAGE = stage
    FREEZE = Path(__file__).with_name(f"regression65_{stage}_v2.json")
    os.environ["C10A01_REGRESSION_STAGE"] = stage
    pilots_live.FREEZE = FREEZE
    pilots_live.OUTPUT = plan_v4.ROOT / (
        f"v2/benchmarks/results/c10a01/regression65-{stage}-v2"
    )
    pilots_live.MODULE = MODULE
    pilots_live.source_check = source_check
    pilots_live.original.sentinel = audited_sentinel
    pilots_live.measured = measured


def main():
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument(
        "--stage",
        choices=("screen", "pair"),
        default=os.environ.get("C10A01_REGRESSION_STAGE", "screen"),
    )
    parser.add_argument("--freeze-seconds", type=int)
    args, rest = parser.parse_known_args()
    configure(args.stage)
    if args.freeze_seconds is not None:
        if FREEZE.exists() or args.freeze_seconds <= 0:
            raise RuntimeError("existing freeze or invalid finite cap")
        FREEZE.write_text(
            json.dumps(
                build(args.stage, args.freeze_seconds),
                indent=2,
                sort_keys=True,
            )
            + "\n"
        )
        return
    if rest and rest[0] == "check":
        value = source_check()
        print(
            json.dumps(
                dict(
                    ready=True,
                    expected_power=value["expected_power"],
                    limits=value["limits"]["total"],
                    host=pilots_live.original.host_admission(2),
                )
            )
        )
    elif rest and rest[0] == "worker":
        worker_parser = argparse.ArgumentParser()
        worker_parser.add_argument("--cell-index", type=int, required=True)
        worker_args = worker_parser.parse_args(rest[1:])
        value = source_check(require_receipts=True)
        if not 0 <= worker_args.cell_index < len(value["calls"]):
            raise ValueError("invalid worker index")
        pilots_live.worker(worker_args.cell_index)
    else:
        sys.argv = [sys.argv[0]] + rest
        pilots_live.main()


if __name__ == "__main__":
    main()
