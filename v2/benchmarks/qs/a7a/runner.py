"""Freeze and execute finite A7a stages; execution requires a later phase."""

import argparse
import fcntl
import hashlib
import json
import os
import resource
import signal
import subprocess
import sys
import time
from contextlib import ExitStack, contextmanager
from pathlib import Path
from tempfile import TemporaryDirectory

from .protocol import (
    CPU_READ_CONTROL,
    IDENTITY,
    INPUTS,
    JOINT_CONTROL,
    LIMITS,
    PACKET_CONTROL,
    ROOT,
    STATISTICS_CONTROL,
    confirmation_tranche,
    corpus,
    schedule,
)


def digest(value):
    return hashlib.sha256(
        json.dumps(value, sort_keys=True).encode()
    ).hexdigest()


def write(path, value):
    with Path(path).open("x") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")


def source():
    commit = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
    ).strip()
    dirty = subprocess.check_output(
        ["git", "status", "--porcelain", "--untracked-files=normal"],
        cwd=ROOT,
        text=True,
    )
    if dirty:
        raise ValueError(
            "commit implementation and required inputs before freezing"
        )
    return commit


def freeze(path):
    payload = dict(
        identity=IDENTITY,
        implementation=source(),
        control="62c474dee93c7597fc578491a1819051fe09434f",
        b12="62c474dee93c7597fc578491a1819051fe09434f",
        packet_control=PACKET_CONTROL,
        statistics_control=STATISTICS_CONTROL,
        cpu_read_control=CPU_READ_CONTROL,
        joint_control=JOINT_CONTROL,
        inputs={
            p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in INPUTS.glob("*.json")
        },
        stages={name: schedule(name) for name in LIMITS},
        note=(
            "H0 is the first current-mainline screen after Q. Historical "
            "A-E schedules require separate admission; "
            "no automatic advancement."
        ),
    )
    write(path, dict(payload=payload, sha256=digest(payload)))


def checked_freeze(path):
    record = json.loads(path.read_text())
    payload = record["payload"]
    if (
        record["sha256"] != digest(payload)
        or payload["implementation"] != source()
    ):
        raise ValueError("freeze/source identity mismatch")
    for name, expected in payload["inputs"].items():
        if (
            hashlib.sha256((INPUTS / name).read_bytes()).hexdigest()
            != expected
        ):
            raise ValueError("frozen corpus/policy changed")
    return record


def validate_output(result, request):
    from .cell import validate_result

    if result["cell"] != request["cell"] or result["n"] != request["n"]:
        raise ValueError("result differs from issued cell")
    for row in result["samples"]:
        validate_result(row, request["n"])
    if result["work"] > request["cell"]["work"]:
        raise ValueError("cell work cap exceeded")


@contextmanager
def stage_lock(stage):
    """Own the shared window and an optional immutable packet control."""
    with ExitStack() as stack:
        lock = stack.enter_context(
            open("/private/tmp/factor-performance.lock", "r+")
        )
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        controls = {}
        commits = {
            "H3": {"packet-control": PACKET_CONTROL},
            "H7": {"statistics-control": STATISTICS_CONTROL},
            "H8": {"cpu-read-control": CPU_READ_CONTROL},
            "H9": {
                "cpu-read-control": CPU_READ_CONTROL,
                "joint-control": JOINT_CONTROL,
            },
            "H12": {"budget-control": CPU_READ_CONTROL},
        }.get(stage, {})
        for implementation, control_commit in commits.items():
            control = Path(
                stack.enter_context(
                    TemporaryDirectory(prefix="factor-a7a-packet-control-")
                )
            )
            archive = control / "source.tar"
            subprocess.run(
                [
                    "git",
                    "archive",
                    "--format=tar",
                    "--output",
                    str(archive),
                    control_commit,
                ],
                cwd=ROOT,
                check=True,
                timeout=30,
            )
            if archive.stat().st_size > 64 * 2**20:
                raise ValueError("packet control export exceeds 64 MiB")
            subprocess.run(
                ["tar", "-xf", str(archive), "-C", str(control)],
                check=True,
                timeout=30,
            )
            archive.unlink()
            controls[implementation] = control
        yield controls


def thermal_status():
    """Retain platform telemetry, including unavailable readings."""
    started = time.monotonic()
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    try:
        result = subprocess.run(
            ["pmset", "-g", "therm"],
            capture_output=True,
            text=True,
            timeout=3,
        )
        after = resource.getrusage(resource.RUSAGE_CHILDREN)
        return dict(
            monotonic=started,
            load_average=os.getloadavg(),
            returncode=result.returncode,
            stdout=result.stdout,
            stderr=result.stderr,
            cpu=after.ru_utime
            + after.ru_stime
            - before.ru_utime
            - before.ru_stime,
        )
    except (OSError, subprocess.TimeoutExpired) as error:
        return dict(monotonic=started, unavailable=repr(error))


def _host_command_status(command):
    """Bound telemetry duration and retained output without hiding failures."""
    started = time.monotonic()
    try:
        result = subprocess.run(
            command, capture_output=True, text=True, timeout=1
        )
        maximum = 128 * 1024
        return dict(
            monotonic=started,
            returncode=result.returncode,
            stdout=result.stdout[:maximum],
            stderr=result.stderr[:maximum],
            truncated=len(result.stdout) > maximum
            or len(result.stderr) > maximum,
        )
    except (OSError, subprocess.TimeoutExpired) as error:
        return dict(monotonic=started, unavailable=repr(error))


def host_pressure_status():
    """Read whole-host pressure outside the measured cell interval."""
    return dict(
        vm=_host_command_status(["vm_stat"]),
        cpu=_host_command_status(["ps", "-axo", "pid=,ppid=,pcpu=,comm="]),
    )


def run_stage(args):
    frozen = checked_freeze(args.freeze)
    stage_started, controller_cpu = time.monotonic(), time.process_time()
    plan = schedule(args.stage, args.nominee)
    if args.stage == "E":
        plan = confirmation_tranche(plan, args.tranche)
    admission = json.loads(args.admission.read_text())
    if (
        admission.get("freeze") != frozen["sha256"]
        or admission.get("stage") != args.stage
        or admission.get("prior_stage_passed") is not True
    ):
        raise ValueError("stage needs a source-bound prior-stage admission")
    if args.stage not in ("A", "H0") and not admission.get(
        "recalibrated_numeric_plan"
    ):
        raise ValueError("settle downstream numerical reservations first")
    order = admission.get("order", [row["id"] for row in plan["cells"]])
    indexed = {row["id"]: row for row in plan["cells"]}
    if len(order) != len(indexed) or set(order) != set(indexed):
        raise ValueError("admitted order changed cell membership")
    plan["cells"] = [indexed[key] for key in order]
    if admission.get("plan_sha256") != digest(plan):
        raise ValueError(
            "admission does not bind the exact numeric ordered plan"
        )
    if (
        args.stage == "E"
        and args.tranche == "E1"
        and not admission.get("e0_ledger_sha256")
    ):
        raise ValueError("E1 requires an explicit review of E0")
    fixtures = {r["id"]: r for r in corpus()["fixtures"]}
    if args.stage == "E":
        if args.confirmation is None:
            raise ValueError("fresh corpus must be frozen before confirmation")
        fresh = json.loads(args.confirmation.read_text())
        training = {r["n"] for r in fixtures.values()}
        if any(r["n"] in training for r in fresh["fixtures"]):
            raise ValueError("confirmation overlaps exposed inputs")
        fixtures.update({r["id"]: r for r in fresh["fixtures"]})
    args.output.mkdir(parents=True, exist_ok=False)
    ledger = dict(
        freeze=frozen["sha256"],
        plan=plan,
        admission=admission,
        attempts=[],
        wall=0,
        cpu=0,
    )
    signatures = {}
    with stage_lock(args.stage) as control_roots:
        for spec in plan["cells"]:
            fixture = fixtures[spec["fixture"]]
            request = dict(cell=spec, n=fixture["n"])
            folder = args.output / spec["id"]
            folder.mkdir()
            write(folder / "request.json", request)
            # Issuance is durable before starting. Missing reports retain the
            # entire reservation; a changed protocol needs a new identity.
            receipt = dict(
                id=spec["id"],
                status="issued",
                wall=spec["wall"],
                cpu=spec["cpu"],
                work=spec["work"],
            )
            ledger["attempts"].append(receipt)
            (args.output / "ledger.json").write_text(
                json.dumps(ledger, indent=2) + "\n"
            )
            started = time.monotonic()
            module = (
                "parent_guard_probe"
                if spec.get("cpu_guard") == "parent"
                else "cell"
            )
            profile = spec.get("cpu_profile")
            capture = None
            if profile is not None:
                if spec.get("cpu_guard") != "peer":
                    raise ValueError("profile must use the incumbent guard")
                module = "diagnostic_hook"
                capture = folder / "captures"
                capture.mkdir()
            command = [
                sys.executable,
                "-B",
                "-m",
                f"v2.benchmarks.qs.a7a.{module}",
                "--request",
                str((folder / "request.json").resolve()),
                "--result",
                str((folder / "result.json").resolve()),
            ]
            implementation = spec.get("implementation", "candidate")
            if implementation not in (
                "candidate",
                "packet-control",
                "statistics-control",
                "cpu-read-control",
                "joint-control",
                "budget-control",
            ):
                raise ValueError("unknown cell implementation")
            execution_root = (
                ROOT
                if implementation == "candidate"
                else control_roots.get(implementation)
            )
            if execution_root is None:
                raise ValueError(
                    "packet control is not admitted for this stage"
                )
            receipt["thermal_start"] = thermal_status()
            measured_started = time.monotonic()
            with (
                (folder / "stdout.txt").open("x") as stdout,
                (folder / "stderr.txt").open("x") as stderr,
            ):
                process = subprocess.Popen(
                    command,
                    cwd=execution_root,
                    env=dict(
                        os.environ,
                        PYTHONPATH=str(execution_root),
                        PYTHONDONTWRITEBYTECODE="1",
                        **(
                            {
                                "FACTOR_A7A_DIAGNOSTIC_DIR": str(
                                    capture.resolve()
                                )
                            }
                            if capture is not None
                            else {}
                        ),
                    ),
                    stdout=stdout,
                    stderr=stderr,
                    start_new_session=True,
                )
                try:
                    if capture is None:
                        process.wait(timeout=spec["wall"])
                    else:
                        from .diagnostics import _capture

                        receipt["capture"] = _capture(
                            process,
                            capture,
                            dict(cell=spec, mode="cpu", cpu_window=profile),
                        )
                    if process.returncode:
                        raise RuntimeError("cell exited unsuccessfully")
                    if (
                        sum(
                            p.stat().st_size
                            for p in folder.rglob("*")
                            if p.is_file()
                        )
                        > spec["storage"]
                    ):
                        raise RuntimeError("recorder storage exhausted")
                    result = json.loads((folder / "result.json").read_text())
                    validate_output(result, request)
                    if args.stage in ("H9", "H10"):
                        for row in result["samples"]:
                            expected = signatures.setdefault(
                                spec["fixture"], row["signature"]
                            )
                            if (
                                row["signature"] != expected
                                or row["stats"].get("scanned")
                                != spec["polynomials"] * 131073
                            ):
                                raise ArithmeticError(
                                    "changed fixed-work traversal or atoms"
                                )
                    if result.get("observer_failure") or result.get(
                        "observer_exhausted"
                    ):
                        raise RuntimeError(
                            "observer harness or sampled RSS stop"
                        )
                    elapsed = time.monotonic() - measured_started
                    if (
                        elapsed > spec["wall"]
                        or result["process_tree_cpu"] > spec["cpu"]
                    ):
                        raise RuntimeError("cell resource allowance exceeded")
                    receipt.update(
                        status="reported",
                        measured_wall=elapsed,
                        wall=elapsed,
                        cpu=result["process_tree_cpu"],
                        work=result["work"],
                    )
                except Exception as error:
                    receipt.update(status="lost", failure=str(error))
                    if process.poll() is None:
                        os.killpg(process.pid, signal.SIGTERM)
                        try:
                            process.wait(timeout=1)
                        except subprocess.TimeoutExpired:
                            os.killpg(process.pid, signal.SIGKILL)
                            process.wait()
            receipt["thermal_end"] = thermal_status()
            receipt["actual_wall"] = time.monotonic() - started
            if receipt["status"] == "reported":
                receipt["wall"] = receipt["actual_wall"]
                receipt["cpu"] += sum(
                    receipt[name].get("cpu", 0)
                    for name in ("thermal_start", "thermal_end")
                )
            write(folder / "receipt.json", receipt)
            ledger["wall"] += receipt["wall"]
            ledger["cpu"] += receipt["cpu"]
            (args.output / "ledger.json").write_text(
                json.dumps(ledger, indent=2) + "\n"
            )
            if receipt["status"] != "reported":
                break
    ledger["all_scheduled_reported"] = len(ledger["attempts"]) == len(
        plan["cells"]
    ) and all(r["status"] == "reported" for r in ledger["attempts"])
    ledger["controller_cpu"] = time.process_time() - controller_cpu
    ledger["controller_wall"] = max(
        0,
        time.monotonic()
        - stage_started
        - sum(r["actual_wall"] for r in ledger["attempts"]),
    )
    ledger["charged_wall_total"] = ledger["wall"] + ledger["controller_wall"]
    ledger["charged_cpu_total"] = ledger["cpu"] + ledger["controller_cpu"]
    ledger["limits_exceeded"] = (
        ledger["charged_wall_total"] > plan["wall"]
        or ledger["charged_cpu_total"] > plan["cpu"]
    )
    write(args.output / "terminal.json", ledger)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("plan", "freeze", "run", "qa"))
    parser.add_argument("--stage", choices=tuple(LIMITS), default="H0")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--freeze", type=Path)
    parser.add_argument("--admission", type=Path)
    parser.add_argument("--nominee", choices=("P1", "P2", "P4"), default="P2")
    parser.add_argument("--confirmation", type=Path)
    parser.add_argument("--tranche", choices=("E0", "E1"), default="E0")
    parser.add_argument("--execute-phase3", action="store_true")
    args = parser.parse_args()
    if args.command == "plan":
        write(
            args.output,
            {stage: schedule(stage, args.nominee) for stage in LIMITS},
        )
    elif args.command == "freeze":
        freeze(args.output)
    else:
        if not args.execute_phase3:
            parser.error("performance and full QA execution belong to Phase 3")
        if sys.implementation.name != "pypy" or sys.version_info[:2] != (
            3,
            11,
        ):
            parser.error("PyPy implementing Python 3.11 is required")
        if args.command == "qa":
            # Historical quiet checks recognize benchmark module commands.
            # This process only supervises checks under the shared lock; give
            # it an explicit controller invocation rather than exempting real
            # competing benchmark processes in the historical checker.
            os.execv(
                sys.executable,
                [
                    sys.executable,
                    "-B",
                    "-c",
                    (
                        "from pathlib import Path; import sys; "
                        "from v2.benchmarks.qs.a7a.qa import run; "
                        "run(Path(sys.argv[1]))"
                    ),
                    str(args.output),
                ],
            )
        else:
            if (
                args.stage == "Q"
                or args.freeze is None
                or args.admission is None
            ):
                parser.error("run needs a freeze, admission, and A-E stage")
            run_stage(args)


if __name__ == "__main__":
    main()
