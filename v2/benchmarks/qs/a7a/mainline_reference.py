"""Literal source-export reference; this file is only a benchmark adapter."""

import argparse
import json
import math
import os
import resource
import signal
import subprocess
import sys
import tempfile
import time
from dataclasses import asdict
from pathlib import Path


def child(request):
    """Import production only from the exported mainline beside this file."""
    from v2.execution.budget import Budget
    from v2.portfolio import PortfolioConfig, factorize_bounded
    from v2.qs.defaults import default_siqs_config

    n, spec = request["n"], request["cell"]
    if sys.implementation.name != "pypy" or sys.version_info[:2] != (3, 11):
        raise RuntimeError("PyPy 3.11 required")
    resource.setrlimit(resource.RLIMIT_CPU, (math.ceil(spec["cpu"]),) * 2)
    config = default_siqs_config(n, memory_bytes=spec["qs_memory_bytes"])
    options = PortfolioConfig(
        siqs=config,
        rho_attempts=0,
        pm1_attempts=0,
        ecm_tiers=(),
        memory_bytes=spec["memory"],
    )
    budget = Budget(
        work_limit=spec["work"],
        seconds=spec["wall"] - 10,
        cpu_seconds=spec["cpu"],
    )
    samples = []

    def call():
        started = time.monotonic()
        run = factorize_bounded(
            n, seed=spec["seed"], config=options, budget=budget
        )
        factors = [
            (f.value, f.exponent, f.certainty.value)
            for f in run.result.factors
        ]
        remaining = list(run.result.remaining)
        if math.prod(p**e for p, e, _ in factors) * math.prod(remaining) != n:
            raise ArithmeticError("reference reconstruction failed")
        return dict(
            complete=run.result.complete,
            factors=factors,
            remaining=remaining,
            reason=run.reason,
            library_seconds=time.monotonic() - started,
        )

    warm_seconds = 0
    if spec["lifetime"] == "warm":
        started = time.monotonic()
        while time.monotonic() - started < 3:
            call()
        warm_seconds = time.monotonic() - started
    for _ in range(spec["samples"]):
        samples.append(call())
    return dict(
        n=n,
        samples=samples,
        warm_seconds=warm_seconds,
        configuration=asdict(config),
        work=budget.used,
        cpu=budget.cpu_used,
        rss_hwm=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    )


def run(plan_path, output):
    """Export the exact pin; no production files come from the task branch."""
    import fcntl

    from .performance_lock import open_lock
    from .protocol import ROOT
    from .runner import write

    plan = json.loads(Path(plan_path).read_text())
    if len(plan["cells"]) != 1:
        raise ValueError("reference batches contain one owned actor")
    spec = plan["cells"][0]
    if spec["arm"] != "MAINLINE" or spec["endpoint"] != "complete":
        raise ValueError("reference endpoint must be complete serial QS")
    if spec.get("features") or spec.get("collector_overrides"):
        raise ValueError("reference cannot carry experimental options")
    if spec["lifetime"] == "cold" and spec["samples"] != 1:
        raise ValueError("cold reference needs a fresh interpreter")
    commit = subprocess.check_output(
        ["git", "rev-parse", plan["reference_source"] + "^{commit}"],
        cwd=ROOT,
        text=True,
    ).strip()
    if commit != plan["reference_source"]:
        raise ValueError("reference needs a full immutable source pin")
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    with open_lock() as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        started = time.monotonic()
        before = resource.getrusage(resource.RUSAGE_CHILDREN)
        status = "issued"
        process = None
        with tempfile.TemporaryDirectory(prefix="a7a-mainline-") as directory:
            root = Path(directory)
            archive = root / "source.tar"
            subprocess.run(
                ["git", "archive", "-o", str(archive), commit],
                cwd=ROOT,
                check=True,
                timeout=20,
            )
            subprocess.run(
                ["tar", "-xf", str(archive), "-C", directory],
                check=True,
                timeout=20,
            )
            adapter = root / "reference_adapter.py"
            adapter.write_bytes(Path(__file__).read_bytes())
            from .combined import fixtures

            n = fixtures(plan["corpus_kind"])[spec["fixture"]]["n"]
            request = root / "request.json"
            write(request, dict(n=n, cell=spec))
            try:
                with (
                    (output / "stdout.json").open("x") as stdout,
                    (output / "stderr.txt").open("x") as stderr,
                ):
                    process = subprocess.Popen(
                        [sys.executable, "-B", str(adapter), str(request)],
                        cwd=root,
                        env=dict(os.environ, PYTHONPATH=str(root)),
                        stdout=stdout,
                        stderr=stderr,
                        start_new_session=True,
                    )
                    code = process.wait(
                        timeout=min(spec["wall"], plan["wall"] - 45)
                    )
                    if code:
                        raise RuntimeError(f"reference child exited {code}")
                body = json.loads((output / "stdout.json").read_text())
                write(output / "result.json", body)
                status = "reported"
            finally:
                if process is not None and process.poll() is None:
                    os.killpg(process.pid, signal.SIGTERM)
                    try:
                        process.wait(timeout=2)
                    except subprocess.TimeoutExpired:
                        os.killpg(process.pid, signal.SIGKILL)
                        process.wait()
                after = resource.getrusage(resource.RUSAGE_CHILDREN)
                write(
                    output / "ledger.json",
                    dict(
                        wall=time.monotonic() - started,
                        cpu=after.ru_utime
                        + after.ru_stime
                        - before.ru_utime
                        - before.ru_stime,
                        reference_source=commit,
                        attempts=[dict(status=status)],
                    ),
                )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("request", type=Path)
    args = parser.parse_args()
    print(json.dumps(child(json.loads(args.request.read_text()))))
