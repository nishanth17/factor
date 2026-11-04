"""Single-worker lifecycle bounds, host admission and parent validation."""

import json
import math
import os
import platform
import resource
import selectors
import signal
import subprocess
import sys
import time
from pathlib import Path

from ...ecm.p52.p52_realistic import check_quiet
from .journal import canonical

OUTPUT_BYTES = 128 * 2**10


def host_admission(required):
    """Observe host capacity before reserving a large native worker."""
    if platform.system() == "Darwin":
        raw = subprocess.check_output(["vm_stat"], text=True)
        page = int(raw.split("page size of ", 1)[1].split()[0])
        values = {
            line.split(":")[0]: int(line.split(":")[1].strip().strip("."))
            for line in raw.splitlines()[1:]
            if ":" in line
        }
        available = page * sum(
            values.get(key, 0)
            for key in ("Pages free", "Pages inactive", "Pages purgeable")
        )
    elif Path("/proc/meminfo").exists():
        values = dict(
            line.split(":", 1)
            for line in Path("/proc/meminfo").read_text().splitlines()
        )
        available = int(values["MemAvailable"].split()[0]) * 1024
    else:
        raise RuntimeError("host admission has no verified memory observer")
    if available < required + 2 * 2**30:
        raise RuntimeError(
            "host available memory cannot cover RSS plus 2 GiB headroom"
        )
    return dict(
        available_bytes=available,
        rss_allowance=required,
        host_headroom=2 * 2**30,
    )


def terminate(process):
    if process.poll() is None:
        os.killpg(process.pid, signal.SIGKILL)
    process.wait()


def collect_output(process, out, err, deadline, rss_limit):
    """Drain pipes into bounded transcripts and enforce lifecycle limits."""
    reason, next_host = None, time.monotonic() + 1
    sizes = {process.stdout: 0, process.stderr: 0}
    with selectors.DefaultSelector() as selector:
        selector.register(process.stdout, selectors.EVENT_READ, out)
        selector.register(process.stderr, selectors.EVENT_READ, err)
        while selector.get_map():
            now = time.monotonic()
            if now >= deadline:
                reason = "lifecycle_timeout"
            if now >= next_host and process.poll() is None:
                sampled = subprocess.run(
                    ["ps", "-o", "rss=", "-p", str(process.pid)],
                    capture_output=True,
                    text=True,
                    check=False,
                )
                if (
                    sampled.stdout.strip()
                    and int(sampled.stdout) * 1024 > rss_limit
                ):
                    reason = "process_rss_limit"
                try:
                    check_quiet({os.getpid(), process.pid})
                except RuntimeError:
                    reason = "host_overlap"
                next_host = now + 5
            if reason is not None:
                terminate(process)
                break
            for key, _ in selector.select(timeout=0.2):
                chunk = os.read(key.fileobj.fileno(), 65536)
                if not chunk:
                    selector.unregister(key.fileobj)
                    continue
                room = OUTPUT_BYTES - sizes[key.fileobj]
                key.data.write(chunk[:room])
                sizes[key.fileobj] += min(room, len(chunk))
                if len(chunk) > room:
                    reason = "output_limit"
                    terminate(process)
                    break
    # EOF can precede process termination. Preserve the same lifecycle cap.
    try:
        process.wait(timeout=max(0, deadline - time.monotonic()))
    except subprocess.TimeoutExpired:
        reason = "lifecycle_timeout"
        terminate(process)
    process.stdout.close()
    process.stderr.close()
    return reason


def validate_result(result, row, cell):
    """Keep independent certificates in the parent, outside algorithms."""
    if (result["n"], result["seed"], result["arm"]) != (
        row["n"],
        cell["seed"],
        cell["arm"],
    ):
        raise ValueError("worker output belongs to another cell")
    divisor = result["divisor"]
    expected = sorted(row["factors"])
    if divisor is not None:
        if sorted(result["factors"]) != expected or result["remaining"]:
            raise ArithmeticError(
                "worker factorization disagrees with independent proofs"
            )
        if not (1 < divisor < row["n"] and row["n"] % divisor == 0):
            raise ArithmeticError("worker divisor is improper")
    elif result["factors"] or result["remaining"] != [row["n"]]:
        raise ArithmeticError("unresolved input was not reconstructed")
    result["complete"] = divisor is not None
    result["certainty"] = ["proven_prime"] * len(result["factors"])
    result["uneven"] = row["uneven"]
    result["cap_seconds"] = cell["seconds"]
    result["capped_seconds"] = (
        min(cell["seconds"], result["seconds"])
        if divisor is not None
        else cell["seconds"]
    )
    # Resource compliance is independent of the historical, uncalibrated
    # canary flag. New analysis retains that flag as a descriptive warning.
    result["resource_valid"] = bool(
        result["rss_bytes"] <= cell["rss_bytes"]
        and result["seconds"] <= cell["seconds"]
        and result["cpu_seconds"] <= cell["seconds"]
        and result.get("reason") != "recorder_limit"
    )
    result["timing_valid"] = bool(
        result["timing_valid"] and result["resource_valid"]
    )
    return result


def launch(cell, row, directory, name, *, execution_root=None):
    """No automatic retry; a missing output becomes a charged censored cell."""
    host = host_admission(cell["rss_bytes"])
    directory = Path(directory)
    request = directory / (name + ".request.json")
    output = directory / (name + ".worker.json")
    errors = directory / (name + ".stderr.txt")
    capture = (
        directory / (name + ".atoms.jsonl")
        if cell["kind"] == "capture"
        else None
    )
    payload = dict(
        cell=cell, capture_path=None if capture is None else str(capture)
    )
    with request.open("xb") as stream:
        stream.write(canonical(payload))
    allowance = cell["seconds"] + cell["warm_seconds"] + cell["margin"]

    def child_limits():
        resource.setrlimit(
            resource.RLIMIT_CPU, (math.ceil(allowance), math.ceil(allowance))
        )
        resource.setrlimit(resource.RLIMIT_FSIZE, (257 * 2**20, 257 * 2**20))

    started = time.monotonic()
    with output.open("xb") as out, errors.open("xb") as err:
        process = subprocess.Popen(
            [
                sys.executable,
                "-m",
                "v2.benchmarks.qs.c9a",
                "worker",
                str(request),
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            start_new_session=True,
            preexec_fn=child_limits,
            cwd=execution_root,
            env={
                **os.environ,
                "PYTHONDONTWRITEBYTECODE": "1",
                **(
                    {"PYTHONPATH": str(Path(execution_root).resolve())}
                    if execution_root is not None
                    else {}
                ),
            },
        )
        reason = collect_output(
            process, out, err, started + allowance, cell["rss_bytes"]
        )
    if reason is None and process.returncode != 0:
        # A Python exception is a harness/correctness failure and stops the
        # stage. Hard resource signals are retained as censored outcomes.
        if process.returncode in (
            -signal.SIGKILL,
            -signal.SIGXCPU,
            -signal.SIGXFSZ,
        ):
            reason = "worker_resource_signal"
        else:
            raise RuntimeError(
                "worker failed; inspect retained stderr, no retry"
            )
    if reason is not None:
        measured, warm = cell["work_limit"], cell["warm_work_limit"]
        observer = cell.get("observer_work_limit", 0)
        return dict(
            fixture=cell["fixture"],
            n=cell["n"],
            seed=cell["seed"],
            arm=cell["arm"],
            digits=cell["digits"],
            uneven=row["uneven"],
            complete=False,
            divisor=None,
            factors=[],
            certainty=[],
            remaining=[cell["n"]],
            reason=reason,
            seconds=cell["seconds"],
            cpu_seconds=cell["seconds"],
            capped_seconds=cell["seconds"],
            cap_seconds=cell["seconds"],
            work=measured + warm + observer,
            measured_work=measured,
            warm_work=warm,
            observer_work=observer,
            timing_valid=False,
            resource_censored=True,
            stats={},
            host=host,
            observation="capture" if capture else "clean",
        )
    if output.stat().st_size > OUTPUT_BYTES:
        raise RuntimeError("worker output exceeded serialization allowance")
    result = validate_result(json.loads(output.read_bytes()), row, cell)
    result.update(host=host, lifecycle_wall=time.monotonic() - started)
    return result
