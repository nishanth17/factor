"""Quiet, source-isolated complete-call measurements for F4 selection.

Execute as a file so a committed mainline S0 checkout never imports candidate
factorer modules. The request contains input/configuration, never certificates
or known factors. Parent-side validation remains independent of these checks.
"""

import argparse
import fcntl
import json
import math
import os
import resource
import sys
import time
from pathlib import Path


def validate_reconstruction(row):
    """Reject malformed raw output before retaining a timing observation."""
    if (
        type(row["n"]) is not int
        or row["n"] < 2
        or type(row["complete"]) is not bool
        or type(row["factorization_complete"]) is not bool
    ):
        raise ValueError("invalid input or completion flags")
    value = 1
    for prime, exponent, certainty in row["factors"]:
        if (
            type(prime) is not int
            or prime < 2
            or type(exponent) is not int
            or not 1 <= exponent <= 4096
            or certainty not in ("proven_prime", "probable_prime")
        ):
            raise ValueError("invalid factor record")
        value *= prime**exponent
    for cofactor in row["remaining"]:
        if type(cofactor) is not int or cofactor < 2:
            raise ValueError("invalid unresolved cofactor")
        value *= cofactor
    if value != row["n"] or row["reconstruct"] != row["n"]:
        raise ValueError("raw output does not reconstruct the input")
    if row["factorization_complete"] != (not row["remaining"]):
        raise ValueError("completion flag contradicts unresolved cofactors")
    if row["endpoint"] == "first" and row["complete"]:
        if not any(
            1 < divisor < row["n"] and row["n"] % divisor == 0
            for divisor in [factor[0] for factor in row["factors"]]
            + row["remaining"]
        ):
            raise ValueError("first-factor result has no split")


def run_accounted(function, budget, totals, **kwargs):
    """Retain work and child charges even when result construction refuses."""
    try:
        return function(budget=budget, **kwargs)
    finally:
        totals["work"] += budget.used
        totals["children"] += budget.prior_cpu


def summarize_stages(events, dropped, seconds):
    """Use existing bounded events; missing attribution remains unknown."""
    ecm = [
        event
        for event in events
        if event.get("stage") in ("ecm", "ecm_workers")
    ]
    success = any(
        event.get("outcome") in ("factor", "factor_found") for event in ecm
    )
    engaged = any(event.get("stats", {}).get("engaged") for event in ecm)
    return dict(
        ecm_reached=bool(ecm)
        or any(
            key.startswith("ecm") and value > 0
            for key, value in seconds.items()
        ),
        ecm_success="success" if success else "unknown" if dropped else "none",
        engaged=True if engaged else None if dropped else False,
        fallback_reached=any(
            event.get("stage") in ("siqs", "sss") for event in events
        )
        or any(
            key.startswith(("siqs", "sss")) and value > 0
            for key, value in seconds.items()
        ),
        dropped_events=dropped,
        coarse_stage_seconds=seconds,
    )


def child(request_path, source_root):
    """Measure cold CLI or repeated calls with one retained executor."""
    if sys.implementation.name != "pypy" or sys.version_info[:2] != (3, 11):
        raise RuntimeError("PyPy implementing Python 3.11 required")
    descriptor = int(os.environ["FACTOR_DIAGNOSTIC_LOCK_FD"])
    expected = os.stat("/private/tmp/factor-performance.lock")
    actual = os.fstat(descriptor)
    if (actual.st_dev, actual.st_ino) != (expected.st_dev, expected.st_ino):
        raise RuntimeError("inherited performance lock identity changed")
    fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
    request = json.loads(Path(request_path).read_text())
    if request["endpoint"] not in ("first", "complete"):
        raise ValueError("portfolio endpoint required")
    if request["mode"] not in ("cold", "warm"):
        raise ValueError("unknown pool condition")
    if request["mode"] == "cold" and request["repetitions"] != 1:
        raise ValueError("each cold CLI sample requires a fresh interpreter")
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    resource.setrlimit(resource.RLIMIT_FSIZE, (4 * 2**20, 4 * 2**20))
    resource.setrlimit(
        resource.RLIMIT_CPU,
        (math.ceil(request["cpu"]), math.ceil(request["cpu"]) + 1),
    )
    sys.path.insert(0, source_root)
    from v2.execution.budget import Budget
    from v2.portfolio import PortfolioConfig, factorize_bounded

    config = PortfolioConfig(**request["config"])
    executor = None
    if request["arm"] == "POLICY":
        from v2.execution.ecm_selection import SelectedECMContext

        executor = SelectedECMContext(
            request["selection_profile"],
            request["selection_profile_sha256"],
            seed=request["seed"],
            prefer_warm=request.get("prefer_warm", False),
        )
    elif request["arm"] != "S0":
        from v2.execution.ecm_workers import ECMExecutor
        from v2.execution.worker_contract import EngagementRule, WorkerSettings

        executor = ECMExecutor(
            seed=request["seed"],
            settings=WorkerSettings(**request["settings"]),
            rule=EngagementRule(**(request.get("rule") or {})),
            force_workers=request["arm"] != "M0" and request["arm"] != "N",
        )
    totals = dict(children=0.0, work=0)
    warmup_calls = 0
    warmup_seconds = 0.0
    observations = []
    stop_reason = None

    def invoke():
        remaining_wall = (
            request["wall"] - (time.monotonic() - request["entry"]) - 3
        )
        remaining_cpu = (
            request["cpu"] - time.process_time() - totals["children"]
        )
        if remaining_wall <= 0 or remaining_cpu <= 0:
            raise TimeoutError("whole-cell allowance exhausted")
        call_seconds = request["call_seconds"]
        if request["mode"] == "cold":
            # A cold call pays imports and startup from the same entry clock.
            call_seconds -= time.monotonic() - request["entry"]
        if call_seconds <= 0:
            raise TimeoutError("cold CLI startup exhausted its allowance")
        budget = Budget(
            work_limit=request["work"] - totals["work"],
            seconds=min(remaining_wall, call_seconds),
            cpu_seconds=remaining_cpu,
        )
        before_executor_work = (
            0
            if executor is None
            else executor.lifetime_stats.get("charged_work", 0)
        )
        started = time.monotonic()
        kwargs = {} if executor is None else {"ecm_executor": executor}
        try:
            run = run_accounted(
                factorize_bounded,
                budget,
                totals,
                n=request["n"],
                seed=request["seed"],
                config=config,
                stop_after_split=request["endpoint"] == "first",
                **kwargs,
            )
        except MemoryError as error:
            if not request.get("record_memory_refusals", False):
                raise
            # No native result was returned. Preserve the original cofactor
            # and all charges; this is a censored refusal, never a fast result.
            row = dict(
                n=request["n"],
                endpoint=request["endpoint"],
                reconstruct=request["n"],
                factorization_complete=False,
                complete=False,
                factors=[],
                remaining=[request["n"]],
                reason="memory_refusal",
                error=str(error)[:512],
                seconds=time.monotonic() - started,
                work=budget.used,
                charged_child_cpu=budget.prior_cpu,
                executor_stats={}
                if executor is None
                else dict(executor.last_stats),
                stage_summary=dict(
                    ecm_reached=None,
                    ecm_success="unknown",
                    engaged=None,
                    fallback_reached=None,
                    dropped_events=0,
                    coarse_stage_seconds={},
                ),
            )
            validate_reconstruction(row)
            return row
        duration = time.monotonic() - started
        factors = [
            [factor.value, factor.exponent, factor.certainty.value]
            for factor in run.result.factors
        ]
        row = dict(
            n=request["n"],
            endpoint=request["endpoint"],
            reconstruct=run.result.reconstruct(),
            factorization_complete=run.result.complete,
            complete=any(
                1 < divisor < request["n"] and request["n"] % divisor == 0
                for divisor in [factor[0] for factor in factors]
                + list(run.result.remaining)
            )
            if request["endpoint"] == "first"
            else run.result.complete,
            factors=factors,
            remaining=list(run.result.remaining),
            reason=run.reason,
            seconds=duration,
            work=budget.used,
            charged_child_cpu=budget.prior_cpu,
            executor_stats={}
            if executor is None
            or executor.lifetime_stats.get("charged_work", 0)
            == before_executor_work
            else dict(executor.last_stats),
            stage_summary=summarize_stages(
                run.events,
                run.dropped_events,
                run.checkpoint["payload"]["state"]["stage_seconds"],
            ),
        )
        validate_reconstruction(row)
        return row

    try:
        if request["mode"] == "warm":
            warm_started = time.monotonic()
            for _ in range(500_000):
                invoke()
                warmup_calls += 1
                warmup_seconds = time.monotonic() - warm_started
                if warmup_seconds >= 3:
                    break
            else:
                raise TimeoutError("finite validated warmup cap exhausted")
        for _ in range(request["repetitions"]):
            observations.append(invoke())
    except TimeoutError as error:
        stop_reason = str(error)
    finally:
        if executor is not None:
            executor.close()
    usage = resource.getrusage(resource.RUSAGE_SELF)
    children = resource.getrusage(resource.RUSAGE_CHILDREN)
    result = dict(
        contract="b12-f4-portfolio-selection-v1",
        portfolio_source=factorize_bounded.__code__.co_filename,
        source_root=source_root,
        n=request["n"],
        observations=observations,
        complete_observation_set=stop_reason is None,
        stop_reason=stop_reason,
        warmup_calls=warmup_calls,
        warmup_seconds=warmup_seconds,
        total_work=totals["work"],
        aggregate_cpu=usage.ru_utime
        + usage.ru_stime
        + children.ru_utime
        + children.ru_stime,
        aggregate_cpu_charged=usage.ru_utime
        + usage.ru_stime
        + totals["children"],
        helper_coverage="unpublished helper/teardown CPU retained in leases",
        parent_hwm=usage.ru_maxrss,
        child_hwm=children.ru_maxrss,
        lifetime={} if executor is None else executor.lifetime_stats,
    )
    print(json.dumps(result, allow_nan=False))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--request", type=Path, required=True)
    parser.add_argument("--source-root", required=True)
    arguments = parser.parse_args()
    child(arguments.request, arguments.source_root)


if __name__ == "__main__":
    main()
