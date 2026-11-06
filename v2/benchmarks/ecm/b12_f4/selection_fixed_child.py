"""Source-isolated fixed curves with stabilized retained metadata."""

import argparse
import fcntl
import hashlib
import json
import math
import os
import resource
import sys
import time
from pathlib import Path
from types import SimpleNamespace


def retained_growth(stats):
    """Separate retained inserts from repeated bounded-cache refusals."""
    decoded = stats.get("decode_misses", 0) - stats.get("decode_refusals", 0)
    if decoded < 0:
        raise ValueError("invalid decoded-cache counters")
    return (
        decoded
        + stats.get("stage_two_build_bytes", 0)
        + stats.get("preparations", 0)
    )


def endpoint(request):
    """Keep this runner engine-only; reject configured fallback paths."""
    value = request.get("endpoint", "fixed")
    if value not in ("fixed", "first"):
        raise ValueError("ECM-only endpoint must be fixed or first")
    if value == "first" and (
        request["config"].get("siqs", "missing") is not None
        or request["config"].get("sss", "missing") is not None
    ):
        raise ValueError("ECM-only runs require explicit disabled fallbacks")
    return value


def checked_divisor(n, divisor):
    """Check the public result independently of the ECM implementation."""
    if divisor is not None and (
        type(divisor) is not int or not 1 < divisor < n or n % divisor
    ):
        raise ValueError("invalid ECM divisor")
    return divisor


class MatchedSerialECM:
    """Run a supplied manifest through the selected source's serial ECM.

    Only curve seeds are adapted. Full-schedule program/chain eligibility,
    setup, both arithmetic stages and recovery remain incumbent code. A warm
    interpreter still pays per-call preparation, as the incumbent API does.
    """

    def __init__(self, request):
        self.assignments = request["assignments"]
        self.manifest = request["manifest"]
        self.lifetime_stats = {}

    def run(self, n, *, config, budget, fixed_work=False, **kwargs):
        from v2.execution.budget import BudgetExhaustedError
        from v2.execution.stage_jobs import advance_job, new_job
        from v2.portfolio import _execution_context

        if fixed_work:
            raise ValueError("matched serial adapter requires first-factor")
        expected = [
            (b1, b2)
            for b1, b2, count in config.ecm_tiers
            for _ in range(count)
        ]
        actual = [(a["b1"], a["b2"]) for a in self.assignments]
        if actual != expected or any(a["n"] != n for a in self.assignments):
            raise ValueError("serial manifest differs from numerical schedule")
        divisor, completed, reason = None, 0, "curves_exhausted"
        try:
            context, programs = _execution_context({}, config, budget)
            for assignment in self.assignments:
                job = new_job(
                    "ecm",
                    n,
                    assignment["seed"],
                    assignment["b1"],
                    assignment["b2"],
                )
                while not job["done"]:
                    advance_job(job, budget, programs or context, config)
                if job["sigma"] != assignment["sigma"]:
                    raise ValueError("incumbent generated a different curve")
                completed += 1
                divisor = checked_divisor(n, job["factor"])
                if divisor is not None:
                    reason = "factor_found"
                    break
        except BudgetExhaustedError:
            reason = budget.reason
        return SimpleNamespace(
            divisor=divisor,
            reason=reason,
            stats=dict(
                manifest=self.manifest,
                engaged=False,
                completed_curves=completed,
                startup_seconds=0.0,
                drain_seconds=0.0,
                reported_child_cpu=0.0,
            ),
        )

    def close(self):
        pass


class SelectedECMOnly:
    """Apply the adopted engagement policy with a matched serial decline."""

    def __init__(self, request):
        from v2.execution.ecm_auto import selected_context

        self.owner = selected_context(
            seed=request["seed"],
            prefer_warm=request["mode"] == "warm",
            parent_memory=request["parent_memory"],
            memory_bytes=request["ecm_memory_bytes"],
        )
        self.serial = MatchedSerialECM(request)

    @property
    def lifetime_stats(self):
        return self.owner.lifetime_stats

    def run(self, n, **kwargs):
        result = self.owner.run_portfolio(n, **kwargs)
        if result is not None:
            return result
        result = self.serial.run(n, **kwargs)
        result.stats.update(self.owner.last_stats)
        result.stats["manifest"] = self.serial.manifest
        return result

    def close(self):
        self.owner.close()


def child(request_path, source_root):
    if sys.implementation.name != "pypy" or sys.version_info[:2] != (3, 11):
        raise RuntimeError("PyPy implementing Python 3.11 required")
    descriptor = int(os.environ["FACTOR_DIAGNOSTIC_LOCK_FD"])
    expected = os.stat("/private/tmp/factor-performance.lock")
    actual = os.fstat(descriptor)
    if (actual.st_dev, actual.st_ino) != (expected.st_dev, expected.st_ino):
        raise RuntimeError("performance lock identity changed")
    fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
    request = json.loads(request_path.read_text())
    measured_endpoint = endpoint(request)
    if request["mode"] not in ("cold", "warm"):
        raise ValueError("unknown pool condition")
    if request["mode"] == "cold" and request["repetitions"] != 1:
        raise ValueError("cold samples require separate interpreters")
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    resource.setrlimit(resource.RLIMIT_FSIZE, (4 * 2**20, 4 * 2**20))
    resource.setrlimit(
        resource.RLIMIT_CPU,
        (math.ceil(request["cpu"]), math.ceil(request["cpu"]) + 1),
    )
    sys.path.insert(0, source_root)
    from v2.execution.budget import Budget
    from v2.portfolio import PortfolioConfig

    config = PortfolioConfig(**request["config"])
    if request["arm"] == "MASTER-MATCHED":
        owner = MatchedSerialECM(request)
    elif request["arm"] == "SELECTED":
        owner = SelectedECMOnly(request)
    else:
        from v2.benchmarks.ecm.b12_f4.profile_runner import component_arm
        from v2.execution.ecm_workers import ECMExecutor
        from v2.execution.worker_contract import WorkerSettings

        owner = ECMExecutor(
            seed=request["seed"],
            settings=WorkerSettings(**request["settings"]),
            force_workers=component_arm(
                request["arm"], request["settings"]["workers"]
            ),
        )
    directory = Path(request["capture"])
    directory.mkdir(exist_ok=False)
    work_used, charged_children = 0, 0.0
    reference, calls, statistics = None, [], []
    outcomes, warmup_outcomes = [], []
    warmup_calls, warmup_seconds, stable_calls = 0, 0.0, 0
    measured_started = None
    measured_wall = 0.0
    stop_reason = None

    def invoke():
        nonlocal work_used, charged_children, reference
        seconds = request["wall"] - (time.monotonic() - request["entry"]) - 4
        cpu = request["cpu"] - time.process_time() - charged_children
        if seconds <= 0 or cpu <= 0:
            raise TimeoutError("whole-cell allowance exhausted")
        budget = Budget(
            work_limit=request.get("work", 10**12) - work_used,
            seconds=seconds,
            cpu_seconds=min(cpu, request.get("call_cpu", cpu)),
        )
        if request.get("call_seconds") is not None:
            budget.seconds = min(budget.seconds, request["call_seconds"])
        started, cpu_started = time.monotonic(), time.process_time()
        try:
            result = owner.run(
                request["n"],
                config=config,
                budget=budget,
                fixed_work=measured_endpoint == "fixed",
                lineage=request.get("lineage", "root"),
            )
        finally:
            work_used += budget.used
            charged_children += budget.prior_cpu
        duration = time.monotonic() - started
        if measured_endpoint == "first":
            divisor = checked_divisor(request["n"], result.divisor)
            if result.reason == "worker_failure":
                raise RuntimeError("worker failure: " + str(result.stats))
            manifest = result.stats["manifest"]
            if reference is None:
                reference = manifest
            elif reference != manifest:
                raise ValueError("retained context changed curve manifest")
            outcome = dict(
                divisor=divisor,
                found=divisor is not None,
                reason=result.reason,
                seconds=duration,
                parent_cpu=time.process_time() - cpu_started,
                reported_child_cpu=result.stats.get("reported_child_cpu", 0),
                charged_child_cpu=budget.prior_cpu,
                work=budget.used,
                manifest=manifest,
            )
            return duration, dict(result.stats), outcome
        reason = result.stats.get("termination_reason", result.reason)
        if reason != "curves_exhausted":
            with (directory / "stopped-checkpoint.json").open("x") as stream:
                json.dump(result.checkpoint, stream)
            if reason == "worker_failure":
                raise RuntimeError("worker failure: " + str(result.stats))
            raise TimeoutError(reason)
        from v2.benchmarks.ecm.b12_f4.profile_runner import signature

        curves = sum(tier[2] for tier in config.ecm_tiers)
        rows = signature(result.checkpoint, request["n"], curves)
        if reference is None:
            reference = rows
        elif reference != rows:
            raise ValueError("retained context changed curve results")
        return duration, dict(result.stats), None

    try:
        if request["mode"] == "warm":
            started = time.monotonic()
            for _ in range(4096):
                _, stats, outcome = invoke()
                if outcome is not None:
                    warmup_outcomes.append(outcome)
                warmup_calls += 1
                stable_calls = (
                    stable_calls + 1 if not retained_growth(stats) else 0
                )
                warmup_seconds = time.monotonic() - started
                if warmup_seconds >= 3 and stable_calls >= 2:
                    break
            else:
                raise TimeoutError("finite metadata warmup did not stabilize")
        measured_started = time.monotonic()
        for _ in range(request["repetitions"]):
            duration, stats, outcome = invoke()
            calls.append(duration)
            statistics.append(stats)
            if outcome is not None:
                outcomes.append(outcome)
        measured_wall = time.monotonic() - measured_started
    except TimeoutError as error:
        if measured_started is not None:
            measured_wall = time.monotonic() - measured_started
        stop_reason = str(error)
    finally:
        close_started = time.monotonic()
        owner.close()
        close_seconds = time.monotonic() - close_started
    usage = resource.getrusage(resource.RUSAGE_SELF)
    children = resource.getrusage(resource.RUSAGE_CHILDREN)
    result = dict(
        contract="b12-f4-fixed-selection-v2",
        complete=stop_reason is None,
        stop_reason=stop_reason,
        n=request["n"],
        arm=request["arm"],
        mode=request["mode"],
        observer="none",
        endpoint=measured_endpoint,
        qs_disabled=measured_endpoint == "first",
        lineage=request.get("lineage", "root"),
        rows=reference,
        outcomes=outcomes,
        warmup_outcomes=warmup_outcomes,
        close_seconds=close_seconds,
        calls=calls,
        stats=statistics,
        measured_wall=measured_wall,
        warmup_seconds=warmup_seconds,
        warmup_calls=warmup_calls,
        warmup_stable_calls=stable_calls,
        lifetime=owner.lifetime_stats,
        work_used=work_used,
        aggregate_cpu=usage.ru_utime
        + usage.ru_stime
        + children.ru_utime
        + children.ru_stime,
        aggregate_cpu_charged=usage.ru_utime
        + usage.ru_stime
        + charged_children,
        helper_coverage="unpublished helper/teardown CPU retained",
        parent_hwm=usage.ru_maxrss,
        rss=None,
    )
    with (directory / "result.json").open("x") as stream:
        json.dump(result, stream, allow_nan=False)
    print(
        json.dumps(
            dict(
                complete=result["complete"],
                stop_reason=stop_reason,
                calls=calls,
                aggregate_cpu=result["aggregate_cpu"],
                aggregate_cpu_charged=result["aggregate_cpu_charged"],
                signature=hashlib.sha256(
                    json.dumps(reference).encode()
                ).hexdigest(),
            )
        )
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--request", type=Path, required=True)
    parser.add_argument("--source-root", required=True)
    arguments = parser.parse_args()
    child(arguments.request, arguments.source_root)


if __name__ == "__main__":
    main()
