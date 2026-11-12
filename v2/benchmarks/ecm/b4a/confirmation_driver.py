"""Run a frozen finite ECM panel with cold and retained-owner costs."""

import argparse
import hashlib
import json
import os
import select
import subprocess
import sys
import time
from pathlib import Path

RUNTIME_FIELDS = {
    "role",
    "workers",
    "schedule",
    "chunk_size",
    "capacity",
    "seed",
    "n",
    "b1",
    "b2",
    "curves",
    "seconds",
    "cpu_seconds",
    "worker_settings",
    "op",
}
METADATA_FIELDS = {
    "owner",
    "tag",
    "lane",
    "arm",
    "case",
    "digits",
    "stratum",
    "repeat",
    "b1",
    "b2",
    "curves",
}


def thermal():
    started = time.monotonic()
    try:
        result = subprocess.run(
            ["/usr/bin/pmset", "-g", "therm"],
            capture_output=True,
            text=True,
            timeout=5,
        )
        return dict(
            exit=result.returncode,
            stdout=result.stdout,
            stderr=result.stderr,
            wall=time.monotonic() - started,
        )
    except subprocess.TimeoutExpired as error:
        return dict(error=repr(error), wall=time.monotonic() - started)


def reap(process):
    """Get one interpreter's self-plus-reaped-descendant CPU exactly once."""
    pid, status, usage = os.wait4(process.pid, 0)
    if pid != process.pid:
        raise RuntimeError("wait4 returned another owner")
    process.returncode = os.waitstatus_to_exitcode(status)
    return usage.ru_utime + usage.ru_stime


def response(process, timeout):
    ready, _, _ = select.select([process.stdout], [], [], timeout)
    if not ready:
        raise TimeoutError("finite confirmation response expired")
    line = process.stdout.readline(2**20 + 1)
    if not line or len(line) > 2**20:
        raise RuntimeError("missing or oversized confirmation response")
    return json.loads(line)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("freeze", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    spec = json.loads(args.freeze.read_text())
    output = args.output.resolve()
    tree = Path(spec["tree"]).resolve()
    if not output.is_dir():
        raise ValueError("declared output directory is missing")
    head = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=tree, text=True
    ).strip()
    if head != spec["source"]:
        raise ValueError("confirmation source changed")
    for name, digest in spec["files"].items():
        if hashlib.sha256(Path(name).read_bytes()).hexdigest() != digest:
            raise ValueError("frozen confirmation input changed: " + name)
    env = dict(os.environ, PYTHONPATH=str(tree), PYTHONDONTWRITEBYTECODE="1")
    command = [
        sys.executable,
        "-B",
        "-m",
        "v2.benchmarks.ecm.b4a.confirmation_runtime",
    ]
    servers, markers, signatures, manifests = {}, {}, {}, {}
    rows = (output / (spec["identity"] + "-rows.jsonl")).open("x")

    def record(row, index):
        rows.write(json.dumps(row) + "\n")
        rows.flush()
        print(
            json.dumps(
                dict(
                    index=index + 1,
                    total=len(spec["requests"]),
                    op=row["op"],
                    arm=row["arm"],
                    case=row["case"],
                    reason=row.get("reason"),
                )
            ),
            flush=True,
        )

    try:
        for index, request in enumerate(spec["requests"]):
            key = request["owner"]
            if request["op"] == "close":
                process, error_file = servers.pop(key)
                started = time.monotonic()
                process.stdin.close()
                cpu = reap(process)
                if process.returncode:
                    raise RuntimeError("retained owner failed during close")
                marker = markers.pop(key)
                if cpu < marker["cpu"] - 0.002:
                    raise ValueError("owner CPU went backwards")
                row = dict(
                    request,
                    close_wall=time.monotonic() - started,
                    owner_cpu=max(0, cpu - marker["cpu"]),
                    owner_calls=marker["calls"],
                    owner_scope_wall=time.monotonic() - marker["at"],
                    reaped_lifetime_cpu=cpu,
                )
                process.stdout.close()
                error_file.close()
                record(row, index)
                continue
            # Fixture/certificate/stratum labels are parent metadata. Only
            # mathematical/runtime parameters cross into production calls.
            body = {k: v for k, v in request.items() if k in RUNTIME_FIELDS}
            before = thermal() if request["op"] == "call" else None
            if request["lane"] == "cold":
                error_path = output / f"{spec['identity']}-cold-{index}.stderr"
                with error_path.open("x") as error_file:
                    own, started = time.process_time(), time.monotonic()
                    process = subprocess.Popen(
                        command + ["cold"],
                        cwd=tree,
                        env=env,
                        stdin=subprocess.PIPE,
                        stdout=subprocess.PIPE,
                        stderr=error_file,
                        text=True,
                    )
                    process.stdin.write(json.dumps(body) + "\n")
                    process.stdin.close()
                    row = response(process, request["seconds"] + 30)
                    cpu = reap(process)
                    wall = time.monotonic() - started
                    process.stdout.close()
                    if process.returncode:
                        raise RuntimeError(
                            "cold confirmation interpreter failed"
                        )
                    row.update(
                        wall=wall,
                        total_cpu=cpu + time.process_time() - own,
                        reaped_lifetime_cpu=cpu,
                    )
            else:
                if key not in servers:
                    error_file = (
                        output / f"{spec['identity']}-{key}.stderr"
                    ).open("x")
                    process = subprocess.Popen(
                        command,
                        cwd=tree,
                        env=env,
                        stdin=subprocess.PIPE,
                        stdout=subprocess.PIPE,
                        stderr=error_file,
                        text=True,
                    )
                    servers[key] = process, error_file
                process = servers[key][0]
                started = time.monotonic()
                process.stdin.write(json.dumps(body) + "\n")
                process.stdin.flush()
                row = response(process, request.get("seconds", 20) + 60)
                row["wall"] = time.monotonic() - started
                if request["op"] == "mark":
                    markers[key] = dict(
                        cpu=row["marker"]["cpu"], calls=0, at=time.monotonic()
                    )
                elif request["op"] == "call" and key in markers:
                    markers[key]["calls"] += 1
            row.update(
                {k: v for k, v in request.items() if k in METADATA_FIELDS}
            )
            if request["op"] == "call":
                row["thermal"] = [before, thermal()]
                if row.get("measured"):
                    row["censored"] = bool(
                        row["resource_stop"]
                        or row["reason"]
                        in (
                            "wall_limit",
                            "cpu_limit",
                            "rss_limit",
                            "work_limit",
                            "cancelled",
                            "native_team_refused",
                            "worker_failure",
                        )
                    )
                    if row["publications"]:
                        first = row["publications"][0]["at"] - started
                        if not 0 <= first <= row["wall"]:
                            raise ValueError(
                                "publication clock is outside owner call"
                            )
                        row["publication_latency"] = first
                        row["post_publication_wall"] = row["wall"] - first
                    if not row["censored"]:
                        identity = (
                            request["n"],
                            request["b1"],
                            request["b2"],
                            request["curves"],
                            request["seed"],
                        )
                        signature = (
                            row["factors"],
                            row["remaining"],
                            row["complete"],
                        )
                        if signature != signatures.setdefault(
                            identity, signature
                        ):
                            raise ArithmeticError(
                                "matched ECM outcomes differ"
                            )
                        if row["manifests"] != manifests.setdefault(
                            identity, row["manifests"]
                        ):
                            raise ArithmeticError(
                                "matched ECM curve streams differ"
                            )
                        if (
                            not row["complete"]
                            and row["completed_curves"] != request["curves"]
                        ):
                            raise ArithmeticError(
                                "failed search did not exhaust its quota"
                            )
            record(row, index)
    finally:
        cleanup = []
        for key, (process, error_file) in servers.items():
            process.stdin.close()
            cpu = reap(process)
            cleanup.append(
                dict(
                    owner=key, exit=process.returncode, reaped_lifetime_cpu=cpu
                )
            )
            process.stdout.close()
            error_file.close()
        rows.close()
        if cleanup:
            (output / (spec["identity"] + "-unplanned-close.json")).write_text(
                json.dumps(cleanup, indent=2) + "\n"
            )


if __name__ == "__main__":
    main()
