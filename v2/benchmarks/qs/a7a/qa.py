"""Phase 3 Stage Q only: full suites, source-only checkout and imports."""

import fcntl
import importlib
import json
import os
import pkgutil
import resource
import signal
import subprocess
import sys
import tempfile
import time
from pathlib import Path

from .protocol import ROOT


def benchmark_imports():
    import v2.benchmarks

    names = [
        m.name
        for m in pkgutil.walk_packages(
            v2.benchmarks.__path__, "v2.benchmarks."
        )
    ]
    for name in names:
        importlib.import_module(name)
    return len(names)


def run(output):
    output.mkdir(parents=True, exist_ok=False)
    receipts = []
    # Includes source checkout and import setup; subprocess limits leave the
    # remaining 480 seconds / 1,920 CPU for controller/failure/cleanup work.
    jobs = [
        (
            "tests",
            480,
            ["make", "-C", "v2", "test", "PYTHON=" + sys.executable],
        ),
        ("lint", 120, ["make", "-C", "v2", "lint"]),
        (
            "fixtures",
            120,
            [
                sys.executable,
                "-B",
                "-m",
                "unittest",
                "v2.tests.test_current_qs_workers",
                "v2.tests.test_a7a_protocol",
                "v2.tests.test_c1_implementation.ImplementationHarnessTests."
                "test_owned_control_and_current_dlp_reconstruct_same_input",
                "-v",
            ],
        ),
    ]
    jobs.sort(key=lambda job: {"fixtures": 0, "lint": 1, "tests": 2}[job[0]])
    with open("/private/tmp/factor-performance.lock", "r+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        with tempfile.TemporaryDirectory(
            prefix="factor-a7a-committed-"
        ) as folder:
            checkout = Path(folder)
            commit = subprocess.check_output(
                ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
            ).strip()
            archive = subprocess.Popen(
                ["git", "archive", commit], cwd=ROOT, stdout=subprocess.PIPE
            )
            subprocess.run(
                ["tar", "-x", "-C", str(checkout)],
                stdin=archive.stdout,
                check=True,
            )
            archive.stdout.close()
            if archive.wait():
                raise RuntimeError("committed source export failed")
            # Development executables may be shared; project imports resolve
            # exclusively from the committed export, with no source symlinks.
            os.symlink(ROOT / "v2/.venv", checkout / "v2/.venv")
            jobs += [
                (
                    "committed-tests",
                    480,
                    ["make", "-C", "v2", "test", "PYTHON=" + sys.executable],
                ),
                (
                    "committed-imports",
                    120,
                    [
                        sys.executable,
                        "-B",
                        "-c",
                        (
                            "from v2.benchmarks.qs.a7a.qa "
                            "import benchmark_imports; "
                            "print(benchmark_imports())"
                        ),
                    ],
                ),
            ]
            for name, wall, command in jobs:
                cwd = checkout if name.startswith("committed") else ROOT
                started = time.monotonic()
                before = resource.getrusage(resource.RUSAGE_CHILDREN)
                with (output / (name + ".txt")).open("x") as stream:
                    environment = dict(
                        os.environ,
                        PYTHONPATH=str(cwd),
                        PYTHONDONTWRITEBYTECODE="1",
                    )
                    process = subprocess.Popen(
                        command,
                        env=environment,
                        cwd=cwd,
                        stdout=stream,
                        stderr=subprocess.STDOUT,
                        start_new_session=True,
                    )
                    try:
                        code = process.wait(timeout=wall)
                    except subprocess.TimeoutExpired:
                        os.killpg(process.pid, signal.SIGTERM)
                        try:
                            process.wait(timeout=1)
                        except subprocess.TimeoutExpired:
                            os.killpg(process.pid, signal.SIGKILL)
                            process.wait()
                        code = process.returncode
                after = resource.getrusage(resource.RUSAGE_CHILDREN)
                cpu = (
                    after.ru_utime
                    + after.ru_stime
                    - before.ru_utime
                    - before.ru_stime
                )
                receipts.append(
                    dict(
                        name=name,
                        command=command,
                        commit=commit,
                        wall=time.monotonic() - started,
                        cpu=cpu,
                        returncode=code,
                    )
                )
                (output / "receipts.json").write_text(
                    json.dumps(receipts, indent=2) + "\n"
                )
                if (
                    code
                    or sum(r["wall"] for r in receipts) > 1800
                    or sum(r["cpu"] for r in receipts) > 7200
                ):
                    raise RuntimeError(
                        "Stage Q failed or exhausted its allowance"
                    )
    (output / "passed.json").write_text(
        json.dumps(
            dict(passed=True, implementation=commit, receipts=receipts),
            indent=2,
        )
        + "\n"
    )
