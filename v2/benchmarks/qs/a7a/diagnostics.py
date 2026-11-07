"""Bounded CPU/heap capture, kept separate from clean performance evidence."""

import argparse
import hashlib
import json
import os
import signal
import subprocess
import sys
import time
from pathlib import Path

from .protocol import INPUTS, ROOT, corpus, schedule
from .runner import digest, source, stage_lock, validate_output, write

IDENTITY = "a7a-diagnostics-20261105-v1"


def plan():
    complete = schedule("H6")["cells"][0]
    cells = []
    for index, mode in enumerate(("off", "cpu", "cpu", "off", "heap")):
        spec = dict(
            complete,
            id=f"D-{index:03}",
            stage="D",
            arm="P4",
            endpoint="prefix",
            polynomials=768,
            wall=45 if mode == "heap" else 30,
            cpu=450 if mode == "heap" else 300,
            storage=3 * 2**30 if mode == "heap" else 256 * 2**20,
        )
        cells.append(dict(cell=spec, mode=mode, cpu_window=[3, 13], heap_at=6))
    cells.append(
        dict(
            cell=dict(
                complete,
                id="D-005",
                stage="D",
                wall=330,
                cpu=3300,
                storage=3 * 2**30,
            ),
            mode="mature",
            cpu_window=[150, 165],
            heap_at=180,
        )
    )
    return dict(
        identity=IDENTITY,
        wall=600,
        cpu=6000,
        storage=8 * 2**30,
        reserve_wall=105,
        reserve_cpu=1050,
        expected_wall=330,
        question=(
            "Identify worker CPU cost and mature parent retention/control; "
            "calibrate CPU recording on identical 768-polynomial work first."
        ),
        decision=(
            "Stop on recorder, validation or resource failure. Retain useful "
            "captures; select a named follow-up, without adopting a speedup."
        ),
        cells=cells,
    )


def freeze(path):
    payload = dict(
        source=source(),
        plan=plan(),
        inputs={
            p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in INPUTS.glob("*.json")
        },
    )
    write(path, dict(payload=payload, sha256=digest(payload)))


def _events(folder, parent):
    owned = {}
    for path in folder.glob("*/events.json"):
        rows = json.loads(path.read_text())
        ready = rows[0]
        if ready["pid"] != int(path.parent.name) or not (
            (ready["role"] == "parent" and ready["pid"] == parent)
            or (ready["role"] == "worker" and ready["ppid"] == parent)
        ):
            raise RuntimeError("diagnostic process ownership mismatch")
        if any(r["kind"] == "recorder_failure" for r in rows):
            raise RuntimeError("diagnostic recorder failed; preserve attempt")
        owned[ready["pid"]] = rows
    return owned


def _signal_live(pid, signum):
    try:
        os.kill(pid, signum)
        return True
    except ProcessLookupError:
        return False


def _capture(process, folder, spec):
    started, active, stopped, heaps = time.monotonic(), set(), set(), set()
    mode = spec["mode"]
    while process.poll() is None:
        elapsed = time.monotonic() - started
        if elapsed >= spec["cell"]["wall"]:
            raise TimeoutError("diagnostic cell wall reservation exhausted")
        owned = _events(folder, process.pid)
        if mode in ("cpu", "mature"):
            start, stop = spec["cpu_window"]
            for pid in owned:
                if start <= elapsed < stop and pid not in active:
                    if _signal_live(pid, signal.SIGUSR1):
                        active.add(pid)
                if elapsed >= stop and pid in active and pid not in stopped:
                    if _signal_live(pid, signal.SIGUSR1):
                        stopped.add(pid)
        if mode == "heap" and elapsed >= spec["heap_at"] and not heaps:
            workers = [p for p, r in owned.items() if r[0]["role"] == "worker"]
            if workers:
                pid = min(workers)
                if _signal_live(pid, signal.SIGUSR2):
                    heaps.add(pid)
        if mode == "mature" and elapsed >= spec["heap_at"] and not heaps:
            if _signal_live(process.pid, signal.SIGUSR2):
                heaps.add(process.pid)
        time.sleep(0.1)
    if process.returncode:
        raise RuntimeError("diagnostic cell exited unsuccessfully")
    owned = _events(folder, process.pid)
    if mode in ("cpu", "mature"):
        # Ordinary worker shutdown can bypass Python atexit. Explicit window
        # stop acknowledgements, rather than nonempty files, establish flush.
        if len(active) != (5 if mode == "cpu" else 1):
            raise RuntimeError("not every expected process was profiled")
        for pid in active:
            kinds = [r["kind"] for r in owned[pid]]
            if (
                kinds.count("cpu_started") != 1
                or kinds.count("cpu_stopped") != 1
            ):
                raise RuntimeError("CPU window lacks an acknowledged flush")
    if mode in ("heap", "mature"):
        if len(heaps) != 1 or not any(
            r["kind"] == "heap_finished" for p in heaps for r in owned[p]
        ):
            raise RuntimeError("live heap capture missing")
    return dict(active=sorted(active), heaps=sorted(heaps), processes=owned)


def run(frozen_path, output):
    if sys.implementation.name != "pypy" or sys.version_info[:2] != (3, 11):
        raise RuntimeError("PyPy implementing Python 3.11 required")
    frozen = json.loads(frozen_path.read_text())
    payload = frozen["payload"]
    if (
        frozen["sha256"] != digest(payload)
        or payload["source"] != source()
        or payload["plan"] != plan()
    ):
        raise ValueError("diagnostic freeze/source/plan mismatch")
    if any(
        hashlib.sha256((INPUTS / name).read_bytes()).hexdigest() != expected
        for name, expected in payload["inputs"].items()
    ):
        raise ValueError("diagnostic frozen input changed")
    output.mkdir(parents=True, exist_ok=False)
    ledger = dict(freeze=frozen, attempts=[])
    fixtures = {r["id"]: r for r in corpus()["fixtures"]}
    signatures = set()
    started, cpu = time.monotonic(), time.process_time()
    with stage_lock("D"):
        for spec in plan()["cells"]:
            if (
                time.monotonic() - started
                > plan()["wall"] - spec["cell"]["wall"]
            ):
                raise RuntimeError("insufficient cumulative diagnostic wall")
            cell = spec["cell"]
            folder = output / cell["id"]
            folder.mkdir()
            capture = folder / "captures"
            capture.mkdir()
            request = dict(cell=cell, n=fixtures[cell["fixture"]]["n"])
            write(folder / "request.json", request)
            receipt = dict(
                id=cell["id"],
                status="issued",
                wall=cell["wall"],
                cpu=cell["cpu"],
            )
            ledger["attempts"].append(receipt)
            (output / "ledger.json").write_text(json.dumps(ledger, indent=2))
            before = time.monotonic()
            module = "cell" if spec["mode"] == "off" else "diagnostic_hook"
            with (
                (folder / "stdout.txt").open("x") as stdout,
                (folder / "stderr.txt").open("x") as stderr,
            ):
                process = subprocess.Popen(
                    [
                        sys.executable,
                        "-B",
                        "-m",
                        f"v2.benchmarks.qs.a7a.{module}",
                        "--request",
                        str((folder / "request.json").resolve()),
                        "--result",
                        str((folder / "result.json").resolve()),
                    ],
                    cwd=ROOT,
                    env=dict(
                        os.environ,
                        PYTHONPATH=str(ROOT),
                        PYTHONDONTWRITEBYTECODE="1",
                        FACTOR_A7A_DIAGNOSTIC_DIR=str(capture.resolve()),
                    ),
                    stdout=stdout,
                    stderr=stderr,
                    start_new_session=True,
                )
                try:
                    receipt["capture"] = _capture(process, capture, spec)
                    result = json.loads((folder / "result.json").read_text())
                    validate_output(result, request)
                    if result["process_tree_cpu"] > cell["cpu"]:
                        raise RuntimeError("cell aggregate CPU cap exceeded")
                    if cell["endpoint"] == "prefix":
                        row = result["samples"][0]
                        signatures.add(row["signature"])
                        if row["stats"].get("scanned") != 768 * 131073:
                            raise RuntimeError(
                                "fixed-work traversal incomplete"
                            )
                        if len(signatures) != 1:
                            raise RuntimeError(
                                "diagnostic changed exact atoms"
                            )
                    else:
                        if not result["samples"][0]["complete"]:
                            raise RuntimeError(
                                "diagnostic complete call unresolved"
                            )
                        actual = sorted(
                            p
                            for p, exponent, _ in result["samples"][0][
                                "factors"
                            ]
                            for _ in range(exponent)
                        )
                        if actual != sorted(
                            fixtures[cell["fixture"]]["factors"]
                        ):
                            raise RuntimeError(
                                "complete fixture oracle mismatch"
                            )
                    size = sum(
                        p.stat().st_size
                        for p in folder.rglob("*")
                        if p.is_file()
                    )
                    if size > cell["storage"]:
                        raise RuntimeError("diagnostic cell storage exhausted")
                    receipt.update(
                        status="reported",
                        wall=time.monotonic() - before,
                        cpu=result["process_tree_cpu"],
                        storage=size,
                    )
                except Exception as error:
                    receipt.update(status="lost", failure=repr(error))
                    if process.poll() is None:
                        os.killpg(process.pid, signal.SIGTERM)
                        try:
                            process.wait(timeout=1)
                        except subprocess.TimeoutExpired:
                            os.killpg(process.pid, signal.SIGKILL)
                            process.wait()
            receipt["actual_wall"] = time.monotonic() - before
            write(folder / "receipt.json", receipt)
            (output / "ledger.json").write_text(json.dumps(ledger, indent=2))
            if receipt["status"] != "reported":
                break
    actual_wall = time.monotonic() - started
    ledger.update(
        actual_wall=actual_wall,
        charged_wall=actual_wall
        + sum(
            max(0, r["wall"] - r.get("actual_wall", r["wall"]))
            for r in ledger["attempts"]
            if r["status"] != "reported"
        ),
        charged_cpu=sum(r["cpu"] for r in ledger["attempts"])
        + time.process_time()
        - cpu,
        all_reported=len(ledger["attempts"]) == len(plan()["cells"])
        and all(r["status"] == "reported" for r in ledger["attempts"]),
        storage=sum(
            p.stat().st_size for p in output.rglob("*") if p.is_file()
        ),
    )
    ledger["limits_exceeded"] = any(
        ledger[field] > plan()[limit]
        for field, limit in (
            ("charged_wall", "wall"),
            ("charged_cpu", "cpu"),
            ("storage", "storage"),
        )
    )
    write(output / "terminal.json", ledger)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("freeze", "run"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--freeze", type=Path)
    args = parser.parse_args()
    if args.command == "freeze":
        freeze(args.output)
    else:
        if args.freeze is None:
            parser.error("run requires a committed-source freeze")
        run(args.freeze, args.output)


if __name__ == "__main__":
    main()
