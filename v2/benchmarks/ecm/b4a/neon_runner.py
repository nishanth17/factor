"""Own the one optional NEON screen inside the common Stage A ledger."""

import fcntl
import json
import time
from pathlib import Path

from .neon_screen import COLUMNS, check_screen_cell, validate_report
from .runner import (
    execution_environment,
    reserve,
    run_owned,
    sha,
    write_new,
)


def read_raw(path, report, cell):
    """Reconstruct measured observations from bounded retained JSONL rows."""
    if (
        path.stat().st_size > cell["raw_bytes"] + 4096
        or sha(path) != report["raw_sha256"]
    ):
        raise ValueError("missing/changed NEON raw batches")
    measured, rows = [], 0
    with path.open() as stream:
        header = json.loads(next(stream))
        if tuple(header["columns"]) != COLUMNS:
            raise ValueError("NEON raw schema changed")
        for line in stream:
            values = json.loads(line)
            if len(values) != len(COLUMNS):
                raise ValueError("invalid NEON raw row")
            row = dict(zip(COLUMNS, values))
            row["arm"] = header["arms"][row["arm"]]
            row["operation"] = header["operation_names"][row["operation"]]
            rows += 1
            if rows > cell["max_rows"]:
                raise ValueError("NEON raw recorder exhausted")
            if row["phase"] == "measure":
                measured.append(row)
    if rows != report["raw_rows"] or measured != report["measurements"]:
        raise ValueError("NEON measurements do not match retained raw rows")


def execute_screen(freeze, cell, directory):
    """No retry/refund: failure retains the same Stage A reservation."""
    check_screen_cell(cell)
    admission = freeze.get("neon_admission", {})
    correct_stage = admission.get("stage") == "NEON-correctness-admission" or (
        admission.get("stage") == "Q" and admission.get("neon_passed") is True
    )
    if (
        not correct_stage
        or admission.get("passed") is not True
        or admission.get("commit") != freeze["commit"]
        or admission.get("timing_executed") is not False
    ):
        raise ValueError("exact NEON build/correctness admission required")
    directory = Path(directory).resolve()
    directory.mkdir(parents=True, exist_ok=True)
    with open("/private/tmp/factor-performance.lock", "r+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        reserve(directory / "ledger.json", freeze, cell)
        output = directory / cell["id"]
        raw = output / "batches.jsonl"
        request = dict(
            cell=cell,
            stage="A",
            admission_passed=True,
            entry=time.monotonic(),
            raw_path=str(raw),
        )
        request_path = directory / (cell["id"] + ".request.json")
        write_new(request_path, request)
        record = run_owned(
            [
                freeze["pypy"],
                "-B",
                "-m",
                "v2.benchmarks.ecm.b4a.neon_screen",
                str(request_path),
            ],
            output,
            cell,
            env=execution_environment(),
            extra_outputs=(raw,),
        )
        started, cpu_started = time.monotonic(), time.process_time()
        try:
            report = json.loads((output / "stdout").read_text())
            validate_report(report, cell)
            read_raw(raw, report, cell)
            status = (
                "screen_complete"
                if report["complete"]
                else "screen_unadmitted"
            )
            outcome = dict(status=status, child=report)
        except (
            ValueError,
            KeyError,
            TypeError,
            OSError,
            IndexError,
            StopIteration,
        ) as error:
            outcome = dict(
                status="recorder_or_validation_failure", failure=repr(error)
            )
        outcome.update(
            validation_seconds=time.monotonic() - started,
            validation_cpu=time.process_time() - cpu_started,
        )
        outcome["observed_cpu"] = (
            record["observed_cpu"] + outcome["validation_cpu"]
        )
        outcome["complete_lifecycle_seconds"] = (
            record["wall"] + outcome["validation_seconds"]
        )
        write_new(output / "validated.json", outcome)
        if (
            record["exit"] != 0
            or record["failure"]
            or not record["group_empty"]
            or record["guardian"].get("reason")
            not in (None, "controller_closed")
            or outcome["observed_cpu"] > cell["cpu"]
            or outcome["complete_lifecycle_seconds"] > cell["wall"]
        ):
            raise RuntimeError(
                "NEON screen stopped; full reservation retained"
            )
        if outcome["status"] != "screen_complete":
            raise RuntimeError(
                "NEON screen unadmitted; diagnose this prototype"
            )
        # A valid negative/inconclusive complete-path result returns normally;
        # it does not prevent the independent scalar/OpenMP Stage A cells.
        return outcome
