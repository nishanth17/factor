"""One cold CLI driver or one explicitly retained warm context lifetime."""

import argparse
import cProfile
import json
import resource
import sys
import time
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("request", type=Path)
    args = parser.parse_args()
    if sys.implementation.name != "pypy" or sys.version_info[:2] != (3, 11):
        raise RuntimeError("PyPy implementing Python 3.11 required")
    request = json.loads(args.request.read_text())

    from .runner import ARMS
    from ....execution.budget import Budget
    from ....execution.ecm_workers import ECMExecutor
    from ....execution.portfolio_config import PortfolioConfig
    from ....execution.worker_contract import WorkerSettings
    from ....portfolio import factorize_bounded

    kernel, workers = ARMS[request["arm"]]
    options = dict(request["config"], ecm_kernel=kernel)
    # Leaving ecm_tiers/siqs omitted preserves automatic root routing intent.
    config = PortfolioConfig(**options)
    n = request["input"]["n"]
    executor = None
    settings = dict(
        request["settings"],
        workers=workers or 1,
        assignment_work=request["assignment_work"],
        assignment_cpu=request["assignment_cpu"],
    )
    if request["arm"].startswith("Native-T"):
        from ....execution.ecm_native import NativeECMExecutor

        executor = NativeECMExecutor(
            seed=request["seed"],
            settings=WorkerSettings(**settings),
            schedule=request.get("schedule", "dynamic"),
            capacity=request.get("capacity"),
        )
    elif workers:
        executor = ECMExecutor(
            seed=request["seed"],
            settings=WorkerSettings(**settings),
            force_workers=workers > 1,
        )
    used_work, warmup_calls, warmup_seconds, charged_children = 0, 0, 0, 0
    first_useful = time.monotonic() - request["entry"]

    def invoke(*, warming=False):
        nonlocal used_work, charged_children
        usage = resource.getrusage(resource.RUSAGE_SELF)
        children = resource.getrusage(resource.RUSAGE_CHILDREN)
        observed = (
            usage.ru_utime
            + usage.ru_stime
            + max(charged_children, children.ru_utime + children.ru_stime)
        )
        remaining_wall = (
            request["wall"]
            - (time.monotonic() - request["entry"])
            - request["shutdown"]
        )
        remaining_cpu = request["cpu"] - observed - 1
        work = max(0, request["work"] - used_work)
        if request["mode"] == "warm":
            if warming:
                remaining_wall = min(
                    remaining_wall - request["measurement_wall"],
                    request["warmup_wall"]
                    - (time.monotonic() - started_warmup),
                )
                remaining_cpu = min(
                    remaining_cpu - request["measurement_cpu"],
                    request["warmup_cpu"] - (observed - warmup_cpu_start),
                )
                work -= request["measurement_work"]
            else:
                if (
                    remaining_wall < request["measurement_wall"]
                    or remaining_cpu < request["measurement_cpu"]
                    or work < request["measurement_work"]
                ):
                    raise RuntimeError(
                        "matched warm measurement reservation refused"
                    )
                remaining_wall = request["measurement_wall"]
                remaining_cpu = request["measurement_cpu"]
                work = request["measurement_work"]
        budget = Budget(
            work_limit=max(0, work),
            seconds=max(0, remaining_wall),
            cpu_seconds=max(0, remaining_cpu),
        )
        started = time.monotonic()
        if request["endpoint"] == "startup":
            if executor is None or (
                workers < 2 and not request["arm"].startswith("Native-T")
            ):
                raise ValueError(
                    "startup endpoint requires an explicit process arm"
                )
            executor.start(budget=budget)
            result = dict(endpoint="startup", n=n, complete=True)
        elif request["endpoint"] == "fixed":
            if executor is None:
                raise ValueError(
                    "fixed work requires equivalent-manifest execution arms"
                )
            run = executor.run(
                n, config=config, budget=budget, fixed_work=True
            )
            search = (
                run.checkpoint["payload"]
                if request["arm"].startswith("Native-T")
                else run.checkpoint["search"]["payload"]
            )
            result = dict(
                endpoint="fixed",
                n=n,
                complete=run.reason == "curves_exhausted",
                reason=run.reason,
                completed=search["completed"],
                manifest=search["manifest"],
                factor=run.divisor,
            )
        else:
            run = factorize_bounded(
                n,
                seed=request["seed"],
                config=config,
                budget=budget,
                ecm_executor=executor,
            )
            result = dict(
                endpoint="complete",
                n=n,
                complete=run.result.complete,
                reconstruct=run.result.reconstruct(),
                reason=run.reason,
                factors=[
                    [f.value, f.exponent, f.certainty.value]
                    for f in run.result.factors
                ],
                remaining=list(run.result.remaining),
                stage_seconds=run.checkpoint["payload"]["state"][
                    "stage_seconds"
                ],
            )
        used_work += budget.used
        charged_children += budget.prior_cpu
        result.update(
            call_seconds=time.monotonic() - started,
            work=budget.used,
            call_charged_cpu=budget.cpu_used,
        )
        return result

    profile = (
        cProfile.Profile() if request.get("observer") == "cprofile" else None
    )
    try:
        if request["mode"] == "warm":
            started_warmup = time.monotonic()
            usage = resource.getrusage(resource.RUSAGE_SELF)
            warmup_cpu_start = usage.ru_utime + usage.ru_stime
            # Both work and total cell clocks are cumulative across warmup.
            for warmup_calls in range(1, 4097):
                checked = invoke(warming=True)
                if not checked["complete"]:
                    raise RuntimeError("warmup is incomplete/censored")
                warmup_seconds = time.monotonic() - started_warmup
                if warmup_seconds >= 3:
                    break
            else:
                raise RuntimeError(
                    "finite warmup did not reach three validated seconds"
                )
        if profile is not None:
            profile.enable()
        result = invoke()
        if profile is not None:
            profile.disable()
            profile.dump_stats(request["profile_path"])
    finally:
        if executor is not None:
            executor.close()

    usage = resource.getrusage(resource.RUSAGE_SELF)
    children = resource.getrusage(resource.RUSAGE_CHILDREN)
    result.update(
        arm=request["arm"],
        mode=request["mode"],
        warmup_seconds=warmup_seconds,
        warmup_calls=warmup_calls,
        total_work=used_work,
        first_dispatch_ready_seconds=first_useful,
        self_cpu=usage.ru_utime + usage.ru_stime,
        reaped_child_cpu=children.ru_utime + children.ru_stime,
        parent_hwm=usage.ru_maxrss,
        child_hwm=children.ru_maxrss,
        hwm_units="bytes",
        stats={} if executor is None else executor.lifetime_stats,
        observed_instrumented=profile is not None,
    )
    print(json.dumps(result, allow_nan=False))


if __name__ == "__main__":
    main()
