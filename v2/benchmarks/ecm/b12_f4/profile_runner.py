"""Finite CPU/GC/heap diagnostics; execute only under an inherited lock."""

import argparse
import fcntl
import hashlib
import json
import os
import resource
import sys
import time
from pathlib import Path

from .diagnostics import (
    Capture,
    RSSSampler,
    captured_worker,
    install_engine_hooks,
)


def signature(checkpoint, n, curves):
    """Independently check every result and retain exact replay identities."""
    raw = checkpoint["search"]["payload"]["completed"]
    completed = {}
    for position, result in raw.items():
        if type(position) is int:
            key = str(position)
        elif isinstance(position, str) and position.isdecimal():
            key = str(int(position))
            if key != position:
                raise ValueError("noncanonical curve position")
        else:
            raise ValueError("invalid curve position")
        if key in completed:
            raise ValueError("duplicate canonical curve position")
        completed[key] = result
    if set(completed) != {str(i) for i in range(curves)}:
        raise ValueError("fixed-work curve set is incomplete")
    rows, ids = [], set()
    for position in range(curves):
        result = completed[str(position)]
        divisor = result["factor"]
        if divisor is not None and not (
            type(divisor) is int and 1 < divisor < n and n % divisor == 0
        ):
            raise ValueError("invalid returned divisor")
        if result["reason"] != "complete" or result["id"] in ids:
            raise ValueError("unfinished or duplicate assignment")
        ids.add(result["id"])
        rows.append(
            [
                result["id"],
                result["sigma"],
                divisor,
                result["phase"],
                result["chain_chunks"],
            ]
        )
    return rows


def component_arm(arm, workers):
    """Validate forced component labels without changing production routing."""
    from .runner import forced_arm

    if arm == "S0":
        raise ValueError("incumbent stream requires the incumbent runner")
    return forced_arm(arm, workers)


def child(path):
    """Run a frozen cell; observer output is never performance evidence."""
    if sys.implementation.name != "pypy" or sys.version_info[:2] != (3, 11):
        raise RuntimeError("PyPy implementing Python 3.11 required")
    inherited = int(os.environ["FACTOR_DIAGNOSTIC_LOCK_FD"])
    expected = os.stat("/private/tmp/factor-performance.lock")
    actual = os.fstat(inherited)
    if (actual.st_dev, actual.st_ino) != (expected.st_dev, expected.st_ino):
        raise RuntimeError("inherited lock identity changed")
    fcntl.flock(inherited, fcntl.LOCK_EX | fcntl.LOCK_NB)
    request = json.loads(path.read_text())
    if request["observer"] not in ("none", "cpu", "counters", "heap", "rss"):
        raise ValueError("unknown observer")
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    output_cap = (32 if request["observer"] == "heap" else 4) * 2**20
    resource.setrlimit(resource.RLIMIT_FSIZE, (output_cap, output_cap))
    resource.setrlimit(resource.RLIMIT_CPU, (100, 101))
    directory = Path(request["capture"])
    directory.mkdir(exist_ok=False)
    capture = None
    if request["observer"] in ("cpu", "counters", "heap"):
        capture = Capture(request, "parent")
    from ....execution.budget import Budget
    from ....execution.ecm_workers import ECMExecutor
    from ....execution.worker_contract import WorkerSettings
    from ....portfolio import PortfolioConfig

    if capture is not None:
        from ....execution import ecm_worker

        os.environ["FACTOR_ECM_DIAGNOSTIC"] = str(path.resolve())
        ecm_worker.worker_main = captured_worker
        install_engine_hooks(capture)
    config = PortfolioConfig(**request["config"])
    executor = ECMExecutor(
        seed=request["seed"],
        settings=WorkerSettings(**request["settings"]),
        force_workers=component_arm(
            request["arm"], request["settings"]["workers"]
        ),
    )
    sampler = None
    if request["observer"] in ("counters", "heap", "rss"):
        sampler = RSSSampler(executor)
        sampler.thread.start()
    charged_children, work_used, warmup_calls = 0.0, 0, 0
    curves = sum(tier[2] for tier in config.ecm_tiers)
    reference, calls, stats = None, [], []

    def invoke():
        nonlocal charged_children, work_used, reference
        wall = request["wall"] - (time.monotonic() - request["entry"]) - 4
        cpu = request["cpu"] - time.process_time() - charged_children
        if wall <= 0 or cpu <= 0:
            raise RuntimeError("diagnostic cell allowance exhausted")
        budget = Budget(
            work_limit=10**12 - work_used,
            seconds=wall,
            cpu_seconds=cpu,
            cancelled=None if sampler is None else lambda: sampler.over_limit,
        )
        started = time.monotonic()
        run = executor.run(
            request["n"], config=config, budget=budget, fixed_work=True
        )
        duration = time.monotonic() - started
        work_used += budget.used
        charged_children += budget.prior_cpu
        if run.reason not in ("curves_exhausted", "factor_found"):
            raise RuntimeError("fixed work incomplete: " + run.reason)
        try:
            rows = signature(run.checkpoint, request["n"], curves)
        except ValueError:
            with (directory / "failure-checkpoint.json").open("x") as stream:
                json.dump(run.checkpoint, stream)
            raise
        if reference is None:
            reference = rows
        elif rows != reference:
            raise ValueError("warmup/repetition changed search results")
        return duration, dict(run.stats)

    started = time.monotonic()
    try:
        if request["mode"] == "warm":
            for _ in range(128):
                invoke()
                warmup_calls += 1
                if time.monotonic() - started >= 3:
                    break
            else:
                raise RuntimeError(
                    "validated warmup failed to reach 3 seconds"
                )
        warmup = time.monotonic() - started if warmup_calls else 0.0
        (directory / "measuring").write_text("measured fixed work\n")
        if capture is not None:
            if request["observer"] == "cpu":
                capture.start_profile()
            elif (
                request["observer"] == "heap"
                and request["settings"]["workers"] > 1
            ):
                capture.dump("parent_retained_pool")
        measured = time.monotonic()
        for _ in range(request["repetitions"]):
            duration, last = invoke()
            calls.append(duration)
            stats.append(last)
        measured = time.monotonic() - measured
        if capture is not None:
            capture.stop_profile()
    finally:
        executor.close()
        if sampler is not None:
            sampler.close()
            sampler.executor = None
        lifetime = dict(executor.lifetime_stats)
        executor = None
        if capture is not None:
            if request["observer"] == "heap":
                capture.dump("parent_after_close")
            capture.finish()
    usage = resource.getrusage(resource.RUSAGE_SELF)
    children = resource.getrusage(resource.RUSAGE_CHILDREN)
    result = dict(
        complete=True,
        n=request["n"],
        arm=request["arm"],
        mode=request["mode"],
        observer=request["observer"],
        rows=reference,
        calls=calls,
        stats=stats,
        measured_wall=measured,
        warmup_seconds=warmup,
        warmup_calls=warmup_calls,
        lifetime=lifetime,
        work_used=work_used,
        aggregate_cpu=usage.ru_utime
        + usage.ru_stime
        + children.ru_utime
        + children.ru_stime,
        aggregate_cpu_charged=usage.ru_utime
        + usage.ru_stime
        + charged_children,
        helper_coverage="CPU reserved; helper CPU/RSS unobserved",
        parent_hwm=usage.ru_maxrss,
        rss=None if sampler is None else sampler.samples,
    )
    with (directory / "result.json").open("x") as stream:
        json.dump(result, stream, allow_nan=False)
    digest = hashlib.sha256(json.dumps(reference).encode()).hexdigest()
    print(
        json.dumps(
            dict(
                complete=True,
                signature=digest,
                measured_wall=measured,
                calls=calls,
                aggregate_cpu=result["aggregate_cpu"],
                aggregate_cpu_charged=result["aggregate_cpu_charged"],
            )
        )
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--request", type=Path, required=True)
    args = parser.parse_args()
    child(args.request)


if __name__ == "__main__":
    main()
