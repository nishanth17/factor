"""B4a frozen-cell runner. Inspection never admits performance execution."""

import argparse
import fcntl
import hashlib
import json
import math
import os
import random
import re
import resource
import select
import signal
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
ARMS = {
    "S0": ("python", 0),
    "M0": ("python", 1),
    "N1": ("c-gmp", 1),
    "Native-S0": ("c-gmp", 0),
    "B12-P2": ("python", 2),
    "B12-P4": ("python", 4),
    "Native-P2": ("c-gmp", 2),
    "Native-P4": ("c-gmp", 4),
    "Native-T1": ("c-gmp", 1),
    "Native-T2": ("c-gmp", 2),
    "Native-T4": ("c-gmp", 4),
}


def execution_environment():
    """Control one parallelism layer without inherited runtime overrides."""
    environment = {
        key: value
        for key, value in os.environ.items()
        if not key.startswith(("OMP_", "KMP_", "GOMP_", "LIBOMP_", "DYLD_"))
    }
    environment.update(
        PYTHONDONTWRITEBYTECODE="1",
        PYTHONPATH=str(ROOT),
        OMP_NUM_THREADS="1",
        OMP_WAIT_POLICY="PASSIVE",
        KMP_BLOCKTIME="0",
        OMP_DYNAMIC="FALSE",
        OMP_MAX_ACTIVE_LEVELS="1",
        OMP_THREAD_LIMIT="4",
        KMP_STACKSIZE="4M",
        LIBOMP_NUM_HIDDEN_HELPER_THREADS="0",
        OPENBLAS_NUM_THREADS="1",
        MKL_NUM_THREADS="1",
        VECLIB_MAXIMUM_THREADS="1",
        NUMEXPR_NUM_THREADS="1",
    )
    return environment


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def canonical(value):
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False
    )


def write_new(path, value):
    with Path(path).open("x") as stream:
        stream.write(
            json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n"
        )


def save(path, value):
    temporary = Path(path).with_suffix(".pending")
    temporary.write_text(canonical(value) + "\n")
    temporary.replace(path)


def check_cell(cell):
    if not re.fullmatch(r"[a-zA-Z0-9_-]{1,128}", cell["id"]):
        raise ValueError("invalid bounded cell identity")
    if cell["arm"] not in ARMS:
        raise ValueError("unimplemented arm")
    if cell.get("schedule", "dynamic") not in ("dynamic", "static", "guided"):
        raise ValueError("unsupported native schedule")
    if (
        type(cell.get("capacity", 8)) is not int
        or not 1 <= cell.get("capacity", 8) <= 32
    ):
        raise ValueError("unsupported native queue capacity")
    if not cell["arm"].startswith("Native-T") and (
        cell.get("schedule", "dynamic") != "dynamic" or "capacity" in cell
    ):
        raise ValueError("native schedule/queue options require a thread arm")
    if cell["mode"] not in ("cold", "warm") or cell["endpoint"] not in (
        "complete",
        "fixed",
        "startup",
    ):
        raise ValueError("invalid endpoint/lifetime")
    for name in (
        "wall",
        "cpu",
        "shutdown",
        "rss",
        "output",
        "work",
        "processes",
    ):
        value = cell[name]
        if (
            type(value) not in (int, float)
            or not math.isfinite(value)
            or value <= 0
        ):
            raise ValueError("invalid finite allowance: " + name)
    if (
        cell["shutdown"] > 3
        or cell["processes"] > 6
        or cell["output"] > 16 * 2**20
    ):
        raise ValueError("cell outside runner ownership bounds")
    if cell["wall"] <= cell["shutdown"] + 1 or cell["rss"] > 4 * 2**30:
        raise ValueError("invalid lifecycle/memory allowance")
    for name in ("work", "output", "processes"):
        if type(cell[name]) is not int:
            raise ValueError("integer ownership allowance required")
    if (
        cell.get("assignment_work") != cell["work"]
        or cell.get("assignment_cpu") != cell["cpu"]
    ):
        raise ValueError("common global-grant assignment caps required")
    if cell["mode"] == "warm":
        for name in (
            "warmup_wall",
            "warmup_cpu",
            "measurement_wall",
            "measurement_cpu",
            "measurement_work",
        ):
            if (
                type(cell.get(name)) not in (int, float)
                or not math.isfinite(cell[name])
                or cell[name] <= 0
            ):
                raise ValueError(
                    "finite matched warmup/measurement grants required"
                )
        if (
            cell["warmup_wall"] + cell["measurement_wall"] + cell["shutdown"]
            > cell["wall"]
            or cell["warmup_cpu"] + cell["measurement_cpu"] + 1 > cell["cpu"]
        ):
            raise ValueError("warmup/measurement partition exceeds cell")


def validate_freeze(freeze, root):
    if freeze.get("status") != "execution-admitted" or not freeze.get(
        "qa", {}
    ).get("passed"):
        raise ValueError("prospective plan is not an execution grant")
    if freeze.get("integration", {}).get("reviewed") is not True:
        raise ValueError("settled B12/mainline integration review required")
    commit = subprocess.check_output(
        ["git", "-C", str(root), "rev-parse", "HEAD"], text=True
    ).strip()
    if commit != freeze["commit"]:
        raise ValueError("source differs from freeze")
    if subprocess.check_output(
        ["git", "-C", str(root), "status", "--porcelain"], text=True
    ):
        raise ValueError(
            "execution checkout must contain only committed sources"
        )
    for name, pin in freeze["sha256"].items():
        path = (root / name).resolve()
        if not path.is_relative_to(root.resolve()) or sha(path) != pin:
            raise ValueError("changed/outside frozen source or input: " + name)
    identities = set()
    for cell in freeze["cells"]:
        check_cell(cell)
        if cell["id"] in identities:
            raise ValueError("duplicate frozen cell")
        identities.add(cell["id"])
    if (
        sum(c["wall"] for c in freeze["cells"]) > freeze["wall"]
        or sum(c["cpu"] for c in freeze["cells"]) > freeze["cpu"]
    ):
        raise ValueError("cells exceed cumulative reservation")


def reserve(path, freeze, cell):
    """Book the full attempt before launch; interruptions never refund it."""
    identity = hashlib.sha256(canonical(freeze).encode()).hexdigest()
    ledger = (
        json.loads(path.read_text())
        if path.exists()
        else dict(identity=identity, wall=0, cpu=0, cells=[])
    )
    if ledger["identity"] != identity or cell["id"] in ledger["cells"]:
        raise ValueError(
            "changed freeze or previously issued cell; no selective retry"
        )
    if (
        ledger["wall"] + cell["wall"] > freeze["wall"]
        or ledger["cpu"] + cell["cpu"] > freeze["cpu"]
    ):
        raise ValueError("cumulative reservation exhausted")
    ledger["wall"] += cell["wall"]
    ledger["cpu"] += cell["cpu"]
    ledger["cells"].append(cell["id"])
    save(path, ledger)


def run_owned(command, directory, cell, *, env=None, extra_outputs=()):
    """Own each child with finite resources and process-group cleanup."""
    from .ownership import GroupObserver

    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=False)
    stdout, stderr = directory / "stdout", directory / "stderr"
    guard_path, request_path = (
        directory / "guardian.json",
        directory / "ownership.json",
    )
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    own, started = time.process_time(), time.monotonic()
    observer = None
    read_pipe, write_pipe = os.pipe()
    gate_read, gate_write = os.pipe()
    process = guardian = None
    failure = None

    def signal_group(sig):
        if process is not None:
            try:
                os.killpg(process.pid, sig)
            except ProcessLookupError:
                pass

    def limits():
        resource.setrlimit(
            resource.RLIMIT_CPU,
            (math.ceil(cell["cpu"]), math.ceil(cell["cpu"]) + 1),
        )
        resource.setrlimit(
            resource.RLIMIT_FSIZE, (cell["output"], cell["output"])
        )

    try:
        observer = GroupObserver()
        with stdout.open("xb") as out, stderr.open("xb") as err:
            # The shell cannot exec arithmetic until its guardian is ready.
            # Parent death closes this gate; a wedged pre-guardian launcher
            # also expires after three seconds. Arguments never enter shell
            # source: exec receives them through the quoted positional array.
            gate = (
                f"IFS= read -r -t 3 factor_start <&{gate_read} || exit 125; "
                '[ "$factor_start" = go ] || exit 126; exec "$@"'
            )
            process = subprocess.Popen(
                ["/bin/sh", "-c", gate, "b4a-launch", *command],
                cwd=ROOT,
                env=env,
                stdout=out,
                stderr=err,
                start_new_session=True,
                preexec_fn=limits,
                pass_fds=(gate_read,),
            )
            os.close(gate_read)
            gate_read = None
            ready = directory / "guardian.ready"
            request = dict(
                pgid=process.pid,
                deadline=started + cell["wall"] - cell["shutdown"],
                cpu=max(0.01, cell["cpu"] - 1),
                rss=cell["rss"],
                processes=cell["processes"],
                threads=cell.get("threads", 32),
                outputs=[str(stdout), str(stderr), *map(str, extra_outputs)],
                output=cell["output"],
                shutdown=cell["shutdown"],
                ready=str(ready),
            )
            write_new(request_path, request)
            guardian = subprocess.Popen(
                [
                    sys.executable,
                    "-B",
                    "-m",
                    "v2.benchmarks.ecm.b4a.runner",
                    "guard",
                    str(request_path),
                    str(read_pipe),
                    str(guard_path),
                ],
                cwd=ROOT,
                env=env,
                pass_fds=(read_pipe,),
                start_new_session=True,
                stdout=subprocess.DEVNULL,
                stderr=err,
            )
            os.close(read_pipe)
            read_pipe = None
            readiness_deadline = time.monotonic() + 1.5
            while not ready.exists():
                if (
                    guardian.poll() is not None
                    or time.monotonic() >= readiness_deadline
                ):
                    raise RuntimeError("guardian did not become ready")
                select.select([], [], [], 0.005)
            os.write(gate_write, b"go\n")
            os.close(gate_write)
            gate_write = None
            try:
                process.wait(timeout=cell["wall"] - cell["shutdown"])
            except subprocess.TimeoutExpired:
                signal_group(signal.SIGTERM)
                try:
                    process.wait(timeout=cell["shutdown"] / 2)
                except subprocess.TimeoutExpired:
                    signal_group(signal.SIGKILL)
                    process.wait(timeout=cell["shutdown"] / 2)
    except BaseException as error:
        failure = repr(error)
    finally:
        # Guardian also cleans up on controller failure/death.
        os.close(write_pipe)
        for descriptor in (gate_read, gate_write):
            if descriptor is not None:
                os.close(descriptor)
        if read_pipe is not None:
            os.close(read_pipe)
        if process is not None and process.poll() is None:
            signal_group(signal.SIGKILL)
            try:
                process.wait(timeout=1)
            except subprocess.TimeoutExpired:
                failure = failure or "owned root did not exit after SIGKILL"
        if guardian is not None:
            try:
                guardian.wait(timeout=cell["shutdown"] + 1)
            except subprocess.TimeoutExpired:
                if process is not None:
                    try:
                        signal_group(signal.SIGKILL)
                    except ProcessLookupError:
                        pass
                guardian.kill()
                try:
                    guardian.wait(timeout=1)
                except subprocess.TimeoutExpired:
                    failure = (
                        failure or "owned guardian did not exit after kill"
                    )

    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    cpu = (
        after.ru_utime
        + after.ru_stime
        - before.ru_utime
        - before.ru_stime
        + time.process_time()
        - own
    )
    group_empty = process is None
    if process is not None:
        try:
            group_empty = not observer.sample(process.pid)["pids"]
        except RuntimeError as error:
            failure = repr(error)
    report = (
        json.loads(guard_path.read_text())
        if guard_path.exists()
        else dict(reason="recorder_failure: missing guardian receipt")
    )
    record = dict(
        command=command,
        exit=None if process is None else process.returncode,
        wall=time.monotonic() - started,
        observed_cpu=cpu,
        cpu_reservation=cell["cpu"],
        group_empty=group_empty,
        guardian=report,
        failure=failure,
        stdout_sha256=sha(stdout) if stdout.exists() else None,
        stderr_sha256=sha(stderr) if stderr.exists() else None,
    )
    write_new(directory / "receipt.json", record)
    return record


def validate_child(data, cell, fixture):
    """Check reconstruction and completion from raw mathematical output."""
    n = fixture["n"]
    if type(n) is not int or not 2 <= n < 2**4096:
        raise ValueError("input outside frozen arithmetic domain")
    if (
        data["n"] != n
        or data["arm"] != cell["arm"]
        or data["mode"] != cell["mode"]
    ):
        raise ValueError("child input/arm/lifetime changed")
    if (
        type(data["complete"]) is not bool
        or type(data["total_work"]) is not int
        or not 0 <= data["total_work"] <= cell["work"]
    ):
        raise ValueError("invalid completion/work report")
    if cell["mode"] == "warm" and data["warmup_seconds"] < 3:
        raise ValueError("warm endpoint lacks validated three-second warmup")
    if data["endpoint"] != cell["endpoint"]:
        raise ValueError("child endpoint changed")
    if cell["endpoint"] == "complete":
        from ....common import utils

        if len(data["factors"]) + len(data["remaining"]) > n.bit_length():
            raise ValueError("unbounded output factor count")
        value = 1
        generator = random.Random(2026110423)
        for prime, exponent, certainty in data["factors"]:
            if (
                type(prime) is not int
                or prime < 2
                or type(exponent) is not int
                or not 1 <= exponent <= 4096
                or certainty not in ("proven_prime", "probable_prime")
            ):
                raise ValueError("invalid factor/certainty record")
            if (prime.bit_length() - 1) * exponent + 1 > n.bit_length():
                raise ValueError("factor power exceeds input")
            power = prime**exponent
            if value > n // power:
                raise ValueError("factors exceed input reconstruction")
            classification = utils.classify_prime(prime, rng=generator)
            if classification is utils.Primality.COMPOSITE or (
                certainty == "proven_prime"
                and classification is not utils.Primality.PROVEN
            ):
                raise ValueError(
                    "factor or certainty failed independent recheck"
                )
            value *= power
        for cofactor in data["remaining"]:
            if type(cofactor) is not int or cofactor < 2:
                raise ValueError("invalid unresolved cofactor")
            if value > n // cofactor:
                raise ValueError("unresolved cofactors exceed input")
            value *= cofactor
        if (
            value != n
            or data["reconstruct"] != n
            or data["complete"] != (not data["remaining"])
        ):
            raise ValueError("raw reconstruction/completion mismatch")
        if data["complete"] and "factors" in fixture:
            if sorted([row[:2] for row in data["factors"]]) != sorted(
                fixture["factors"]
            ):
                raise ValueError(
                    "completed factorization differs from independent fixture"
                )
    elif cell["endpoint"] == "fixed":
        identities = set()
        for position, row in data["completed"].items():
            expected = cell.get("expected_assignments", {}).get(position)
            if (
                expected is None
                or row["id"] != expected["id"]
                or row["sigma"] != expected["sigma"]
            ):
                raise ValueError(
                    "curve ID/sigma differs from pre-frozen assignment"
                )
            if row["id"] in identities:
                raise ValueError("duplicate accepted curve publication")
            identities.add(row["id"])
            divisor = row["factor"]
            if divisor is not None and not (
                type(divisor) is int and 1 < divisor < n and n % divisor == 0
            ):
                raise ValueError("invalid curve divisor")
            if row["reason"] != "complete":
                raise ValueError("incomplete curve labeled accepted")
        if data["complete"] and len(identities) != cell["curves"]:
            raise ValueError("fixed-work manifest coverage incomplete")


def cli_output(path, cell, fixture):
    """Decode raw output from the public CLI's canonical snapshot."""
    snapshot = json.loads(Path(path).read_text())
    payload = snapshot["payload"]
    if (
        hashlib.sha256(canonical(payload).encode()).hexdigest()
        != snapshot["sha256"]
    ):
        raise ValueError("literal CLI checkpoint checksum mismatch")
    expected_kernel = ARMS[cell["arm"]][0]
    if payload["config"].get("ecm_kernel", "python") != expected_kernel:
        raise ValueError("literal CLI executor differs from frozen arm")
    state = payload["state"]
    if state["original"] != fixture["n"]:
        raise ValueError("literal CLI checkpoint belongs to another input")
    remaining = list(state["remaining"])
    pending = list(state["pending"])
    if state["current"] is not None:
        pending.append([state["current"]["n"], state["current"]["mult"]])
    for value, multiplicity in pending:
        if (
            type(multiplicity) is not int
            or not 1 <= multiplicity <= 4096
            or len(remaining) + multiplicity > 4096
        ):
            raise ValueError("unbounded literal CLI cofactor multiplicity")
        remaining.extend([value] * multiplicity)
    factors = [
        [int(value), exponent, certainty]
        for value, (exponent, certainty) in state["factors"].items()
    ]
    return dict(
        n=fixture["n"],
        arm=cell["arm"],
        mode="cold",
        endpoint="complete",
        complete=not remaining,
        total_work=payload["work_used"],
        factors=factors,
        remaining=remaining,
        reconstruct=fixture["n"],
        call_seconds=payload["wall_used"],
        warmup_seconds=0,
        source_endpoint="literal-v2.factor-CLI",
    )


def execute(freeze_path, cell_id, directory):
    freeze = json.loads(Path(freeze_path).read_text())
    validate_freeze(freeze, ROOT)
    cell = next(c for c in freeze["cells"] if c["id"] == cell_id)
    directory = Path(directory).resolve()
    directory.mkdir(parents=True, exist_ok=True)
    with open("/private/tmp/factor-performance.lock", "r+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        reserve(directory / "ledger.json", freeze, cell)
        profile_path = directory / cell_id / "observer.profile"
        request = dict(
            cell,
            input=freeze["inputs"][cell["input"]],
            config=freeze["configurations"][cell["config"]],
            settings=freeze["settings"],
            entry=time.monotonic(),
            profile_path=str(profile_path),
        )
        request_path = directory / (cell_id + ".request.json")
        write_new(request_path, request)
        command = [
            freeze["pypy"],
            "-B",
            "-m",
            "v2.benchmarks.ecm.b4a.child",
            str(request_path),
        ]
        cli_path = directory / cell_id / "result.checkpoint.json"
        if cell.get("literal_cli"):
            if (
                cell["mode"] != "cold"
                or cell["arm"] not in ("S0", "Native-S0")
                or request["config"]
            ):
                raise ValueError("literal CLI arm/configuration changed")
            command = [
                freeze["pypy"],
                "-B",
                "-m",
                "v2.factor",
                str(request["input"]["n"]),
                "--bounded",
                "--seed",
                str(cell["seed"]),
                "--work-limit",
                str(cell["work"]),
                "--seconds",
                str(cell["wall"] - cell["shutdown"] - 2),
                "--cpu-seconds",
                str(cell["cpu"] - 2),
                "--checkpoint",
                str(cli_path),
            ]
            if cell["arm"] == "Native-S0":
                command.extend(["--ecm-kernel", "c-gmp"])
        environment = execution_environment()
        record = run_owned(
            command,
            directory / cell_id,
            cell,
            env=environment,
            extra_outputs=(profile_path, cli_path),
        )
        validation_started, validation_cpu = (
            time.monotonic(),
            time.process_time(),
        )
        try:
            data = (
                cli_output(cli_path, cell, request["input"])
                if cell.get("literal_cli")
                else json.loads((directory / cell_id / "stdout").read_text())
            )
            validate_child(data, cell, request["input"])
            outcome = dict(
                status="complete" if data["complete"] else "incomplete",
                child=data,
            )
        except (
            ValueError,
            KeyError,
            TypeError,
            OverflowError,
            OSError,
        ) as error:
            outcome = dict(
                status="production_resource_stop"
                if record["guardian"].get("reason")
                in ("wall_limit", "cpu_limit", "rss_limit")
                else "recorder_or_output_failure",
                failure=repr(error),
            )
        outcome["validation_seconds"] = time.monotonic() - validation_started
        outcome["validation_cpu"] = time.process_time() - validation_cpu
        outcome["complete_lifecycle_seconds"] = (
            record["wall"] + outcome["validation_seconds"]
        )
        outcome["observed_cpu"] = (
            record["observed_cpu"] + outcome["validation_cpu"]
        )
        write_new(directory / cell_id / "validated.json", outcome)
        expected_exit = (
            1
            if cell.get("literal_cli")
            and outcome.get("status") == "incomplete"
            else 0
        )
        if (
            not record["group_empty"]
            or record["exit"] != expected_exit
            or record["failure"]
            or outcome["complete_lifecycle_seconds"] > cell["wall"]
            or outcome["observed_cpu"] > cell["cpu"]
            or record["guardian"].get("reason")
            not in (None, "controller_closed")
        ):
            raise RuntimeError(
                "cell failed/censored; retain charge and diagnose"
            )
        if outcome["status"] == "recorder_or_output_failure":
            raise RuntimeError(
                "invalid/missing output; diagnose the retained receipt"
            )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("inspect", "execute", "guard"))
    parser.add_argument("path", type=Path)
    parser.add_argument("extra", nargs="*")
    args = parser.parse_args()
    if args.action == "inspect":
        plan = json.loads(args.path.read_text())
        print(
            canonical(
                dict(
                    identity=plan["identity"],
                    status=plan["status"],
                    cells=len(plan.get("cells", [])),
                    executable=plan["status"] == "execution-admitted",
                )
            )
        )
    elif args.action == "guard":
        from .ownership import guard

        guard(
            json.loads(args.path.read_text()),
            int(args.extra[0]),
            args.extra[1],
        )
    else:
        if sys.implementation.name != "pypy" or sys.version_info[:2] != (
            3,
            11,
        ):
            raise RuntimeError("PyPy Python 3.11 required")
        if len(args.extra) != 2:
            parser.error("execute requires CELL_ID RECEIPT_DIRECTORY")
        execute(args.path, *args.extra)


if __name__ == "__main__":
    main()
