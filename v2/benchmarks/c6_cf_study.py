"""Run the frozen CF extension through the unchanged paired C6 harness.

Only protocol, candidate construction, worker module and selection change.
All timed attempts, validations, kernels, budgets and sampling stay shared.
"""

import hashlib
import json
import os
import subprocess
import sys
import time
from pathlib import Path

from . import c6_cf
from . import c6_fast_study as study
from .build_c6_cf_inputs import PROTOCOL

SELECTION = Path(__file__).parent / "inputs/controls/c6_cf_selection.json"
ORIGINAL_CONSTRUCT = study.construct


def construct(arm, backend):
    if arm.startswith("fast/cf/"):
        _, _, mode, batch = arm.split("/")
        return c6_cf.build_program(2000, mode, backend, int(batch))
    return ORIGINAL_CONSTRUCT(arm, backend)


def select(path):
    data, settings = json.loads(path.read_text()), study.protocol()
    result = dict(
        screen_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        protocol_sha256=hashlib.sha256(PROTOCOL.read_bytes()).hexdigest(),
        candidates={},
    )
    for backend in settings["backends"]:
        rows = [
            row
            for row in data["captures"]
            if row["backend"] == backend
            and row["arm"].startswith("fast/cf/")
            and row["summary"]["stable"]
        ]
        if not rows:
            raise RuntimeError("no stable CF candidate")
        best = min(row["summary"]["ratio"] for row in rows)
        winner = min(
            (row for row in rows if row["summary"]["ratio"] <= best * 1.01),
            key=lambda row: (
                ("tuple", "three").index(row["arm"].split("/")[2]),
                int(row["arm"].split("/")[3]),
            ),
        )
        result["candidates"][backend] = [winner["arm"]]
    with SELECTION.open("x") as output:
        json.dump(result, output, indent=2)
        output.write("\n")
    print(json.dumps(result, indent=2))


def launch(args, arm, backend, warmup, count, *, single=False, profile=None):
    command = [
        sys.executable,
        "-B",
        "-m",
        __spec__.name,
        "--worker",
        "--arm",
        arm,
        "--backend",
        backend,
        "--scope",
        args.scope,
        "--warmup",
        str(warmup),
        "--samples",
        str(count),
    ]
    if args.confirmation:
        command.append("--confirmation")
    if single:
        command.append("--single")
    if profile is not None:
        command += ["--profile", str(profile)]
    started = time.perf_counter()
    process = subprocess.Popen(
        command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    try:
        while True:
            try:
                stdout, stderr = process.communicate(timeout=1)
                break
            except subprocess.TimeoutExpired:
                study.p52_realistic.check_quiet({os.getpid(), process.pid})
                if time.perf_counter() - started > 1200:
                    raise RuntimeError("CF worker exceeded finite limit")
        study.p52_realistic.check_quiet({os.getpid(), process.pid})
        if process.returncode:
            raise RuntimeError(stderr)
    except BaseException:
        process.terminate()
        process.communicate()
        raise
    return json.loads(stdout), time.perf_counter() - started


def install():
    study.PROTOCOL, study.SELECTION = PROTOCOL, SELECTION
    study.construct, study.select, study.launch = construct, select, launch


if __name__ == "__main__":
    install()
    study.main()
