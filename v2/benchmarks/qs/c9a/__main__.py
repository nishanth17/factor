"""Prepare C9a or execute one selected stage with its recorded decision."""

import argparse
import fcntl
import json
import os
import time
from pathlib import Path

from ...ecm.p52.p52_realistic import check_quiet
from ..c10a04.accounting import output_lock, receipt
from .accounting import StudyLedger
from .analysis import intake_summary, observer_summary, summarize
from .confirmation import freeze_confirmation, fresh, select
from .gates import admit, outcomes
from .inputs import validate
from .launcher import launch
from .plan import observer_cells
from .preparation import INPUTS, ROOT, freeze, load_freeze, prepare
from .replay import replay
from .worker import run, runtime


def parser():
    result = argparse.ArgumentParser(description=__doc__)
    commands = result.add_subparsers(dest="command", required=True)
    worker = commands.add_parser("worker", help=argparse.SUPPRESS)
    worker.add_argument("request", type=Path)
    preparation = commands.add_parser("prepare")
    preparation.add_argument("--exclusions", type=Path, required=True)
    frozen = commands.add_parser("freeze")
    frozen.add_argument("--output", type=Path, required=True)
    frozen.add_argument("--combined", action="store_true")
    frozen.add_argument("--extended", action="store_true")
    frozen.add_argument("--integration", type=Path)
    frozen.add_argument(
        "--profile",
        default="three-prime",
        choices=(
            "fixed-report",
            "wide-report",
            "wide-slp",
            "three-prime",
            "adaptive-report",
            "balanced-endpoints",
            "bounded-report",
            "cost-report",
        ),
    )
    selection = commands.add_parser("select-confirmation")
    selection.add_argument("--freeze", type=Path, required=True)
    selection.add_argument("--output", type=Path, required=True)
    selection.add_argument(
        "--band", type=int, choices=(60, 62, 65, 70), required=True
    )
    selection.add_argument("--decision", type=Path, required=True)
    fresh_command = commands.add_parser("prepare-confirmation")
    fresh_command.add_argument("--selection", type=Path, required=True)
    fresh_command.add_argument("--exclusions", type=Path, required=True)
    fresh_freeze = commands.add_parser("freeze-confirmation")
    fresh_freeze.add_argument("--selection", type=Path, required=True)
    fresh_freeze.add_argument("--output", type=Path, required=True)
    for name in ("run", "observer", "replay", "analyze"):
        command = commands.add_parser(name)
        command.add_argument("--freeze", type=Path, required=True)
        command.add_argument("--output", type=Path, required=True)
        if name == "run":
            command.add_argument("--continue-early", action="store_true")
            command.add_argument(
                "--stage",
                choices=("A", "B1", "B2", "M", "R", "F", "C"),
                required=True,
            )
            command.add_argument("--decision", type=Path, required=True)
        elif name == "replay":
            command.add_argument("--stream", type=Path, required=True)
            command.add_argument(
                "--arm",
                choices=("C0", "C1", "noP", "noG", "noD", "noS"),
                required=True,
            )
            command.add_argument(
                "--slot", type=int, choices=range(8), required=True
            )
        elif name == "analyze":
            command.add_argument(
                "--stage",
                choices=("A", "B1", "B2", "R", "F", "C"),
                required=True,
            )
    return result


def execution_cells(plan, stage, continue_early, decision, ledger):
    """Keep initial names stable so continuation reuses exact receipts."""
    selected = [cell for cell in plan["cells"] if cell["stage"] == stage]
    initial = plan.get("initial_cells")
    initial_only = initial and stage == "A" and not continue_early
    if initial_only:
        selected = selected[:initial]
    elif initial and stage == "A":
        if len(outcomes(ledger, "A")) < initial or not decision.get(
            "feasibility_review"
        ):
            raise RuntimeError(
                "review the four clean initial cells before enlarging A"
            )
    return selected, "-initial" if initial_only else ""


def execute(args):
    started, cpu = time.monotonic(), time.process_time()
    frozen = load_freeze(args.freeze)
    corpus_path = (
        ROOT / frozen["corpus_path"]
        if "corpus_path" in frozen
        else INPUTS / "corpus.json"
    )
    corpus = json.loads(corpus_path.read_text())
    validate(corpus, band=frozen.get("confirmation_band"))
    fixtures = {row["id"]: row for row in corpus["fixtures"]}
    with output_lock(args.output):
        ledger = StudyLedger(args.output, frozen["sha256"], frozen["limits"])
        stage = args.stage if args.command in ("run", "analyze") else "A"
        grant = max(10, time.monotonic() - started, time.process_time() - cpu)
        # Source/interpreter hashing and input validation precede launch but
        # still consume resources. Charge their entire invocation allowance,
        # including repeated controller starts after an interruption.
        read_bytes = sum(
            (ROOT / name).stat().st_size for name in frozen["source"]["files"]
        )
        read_bytes += Path(frozen["source"]["interpreter"]).stat().st_size
        ledger.call_controller(
            stage,
            grant,
            10**10,
            lambda: dict(work=8 * read_bytes, checked_source=frozen["sha256"]),
        )
        if args.command == "run":
            if (args.output / "confirmation-selection.json").exists():
                raise RuntimeError(
                    "exposed execution closed by confirmation selection"
                )
            decision = json.loads(args.decision.read_text())
            admit(args.stage, decision, frozen, ledger)
            selected, suffix = execution_cells(
                frozen["plan"],
                args.stage,
                args.continue_early,
                decision,
                ledger,
            )
            receipt(
                args.output / (args.stage + suffix + "-decision.json"),
                decision,
            )
            if not selected:
                raise ValueError(
                    "no frozen cells; confirmation needs its fresh freeze"
                )
            for index, cell in enumerate(selected):
                name = (
                    f"{args.stage}-{index:03}-{cell['fixture']}-{cell['arm']}"
                )
                ledger.call_cell(
                    args.stage,
                    name,
                    cell,
                    lambda: launch(
                        cell, fixtures[cell["fixture"]], args.output, name
                    ),
                )
        elif args.command == "observer":
            if len(outcomes(ledger, "A")) < frozen["plan"].get(
                "early_cells", 40
            ):
                raise RuntimeError(
                    "clean representative prefixes precede instrumentation"
                )
            results = []
            for index, cell in enumerate(observer_cells(frozen["plan"])):
                name = f"A-observer-{index:02}"
                results.append(
                    ledger.call_cell(
                        "A",
                        name,
                        cell,
                        lambda: launch(
                            cell, fixtures[cell["fixture"]], args.output, name
                        ),
                    )
                )
            ledger.call_operation(
                "A",
                "A-observer-analysis",
                10,
                10**7,
                {"kind": "observer_analysis"},
                lambda: observer_summary(results),
            )
        elif args.command == "replay":
            if len(outcomes(ledger, "B1")) != 4:
                raise RuntimeError("mature replay allocation follows B1")
            streams = {
                row.get("identity", {}).get("stream")
                for row in ledger.data["calls"]
                if row.get("identity", {}).get("kind") == "mature_replay"
            }
            streams.discard(None)
            if str(args.stream.resolve()) not in streams and len(streams) >= 2:
                raise RuntimeError(
                    "mature replay pool permits at most two streams"
                )
            ledger.call_operation(
                "A",
                f"A-replay-{args.slot:02}",
                190,
                10**13,
                dict(
                    kind="mature_replay",
                    stream=str(args.stream.resolve()),
                    arm=args.arm,
                ),
                lambda: replay(
                    args.stream,
                    args.arm,
                    combined=frozen["plan"].get("early_cells") == 32,
                    extended=frozen["plan"].get("extension_policy")
                    in ("extended-v1", "extended-v2"),
                ),
            )
        else:
            results = outcomes(ledger, args.stage)
            if args.stage == "B2":
                results = outcomes(ledger, "B1") + results
            ledger.call_operation(
                args.stage,
                args.stage + "-analysis",
                120,
                10**9,
                {"kind": "analysis"},
                lambda: intake_summary(results)
                if args.stage == "A"
                else summarize(results),
            )
        return dict(
            calls=len(ledger.data["calls"]),
            charged_wall=ledger.used(),
            charged_cpu=ledger.used(field="cpu"),
            charged_work=ledger.used(field="work"),
        )


def main():
    args = parser().parse_args()
    runtime()
    if args.command == "worker":
        request = json.loads(args.request.read_bytes())
        capture = request["capture_path"]
        result = run(
            request["cell"],
            capture_path=None if capture is None else Path(capture),
        )
    elif args.command == "freeze":
        result = freeze(
            args.output,
            combined=args.combined,
            extended=args.extended,
            integration=args.integration,
            profile=args.profile,
        )
        result = dict(freeze=result["sha256"], output=str(args.output))
    elif args.command == "freeze-confirmation":
        result = freeze_confirmation(args.selection, args.output)
    else:
        with Path("/private/tmp/factor-performance.lock").open("a+") as lease:
            fcntl.flock(lease, fcntl.LOCK_EX | fcntl.LOCK_NB)
            check_quiet({os.getpid()})
            if args.command == "prepare":
                result = prepare(args.exclusions)
            elif args.command == "prepare-confirmation":
                result = fresh(args.selection, args.exclusions)
            elif args.command == "select-confirmation":
                result = select(
                    args.freeze,
                    args.output,
                    args.band,
                    json.loads(args.decision.read_text()),
                )
            else:
                result = execute(args)
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
