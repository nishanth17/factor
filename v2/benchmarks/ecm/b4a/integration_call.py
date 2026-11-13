"""One cold, matched ECM-to-QS portfolio call for integration guards."""

import json
import sys
import time
from pathlib import Path

from .confirmation_runtime import close_spawn_helper
from .integration_observation import IntegrationObserver


def run_request(request):
    from ....execution.budget import Budget
    from ....execution.portfolio_config import PortfolioConfig
    from ....execution.worker_contract import WorkerSettings
    from ....portfolio import factorize_bounded
    from ....qs.current_parallel import QSWorkerConfig
    from ....qs.worker_configuration import decode_config

    role = request["role"]
    if role == "automatic":
        from ....ecm.auto import AutomaticECMOwner

        owner = AutomaticECMOwner(seed=request["seed"])
    elif role == "incumbent-auto":
        from ....ecm.auto import selected_context

        owner = selected_context(seed=request["seed"])
    elif role == "b12":
        from ....ecm.workers import ECMExecutor

        owner = ECMExecutor(
            seed=request["seed"],
            settings=WorkerSettings(workers=4),
            force_workers=True,
        )
    elif role == "pthread":
        from ....execution.ecm_native import NativeECMExecutor

        owner = NativeECMExecutor(
            seed=request["seed"],
            settings=WorkerSettings(
                workers=request.get("workers", 8), allow_wide_workers=True
            ),
            schedule="dynamic",
            chunk_size=1,
            capacity=32,
        )
    else:
        raise ValueError("unknown integration ECM role")
    values = dict(
        trial_bound=2,
        rho_attempts=0,
        pm1_attempts=0,
        ecm_tiers=tuple(tuple(t) for t in request["ecm_tiers"]),
        ecm_chain_mode="auto",
        siqs=decode_config(request["qs_config"]),
        max_input_bits=1024,
    )
    if role == "pthread":
        values["ecm_kernel"] = "c-gmp-full"
    if role in ("automatic", "incumbent-auto"):
        # Exercise ordinary pretests and numerical recipes. Automatic intent
        # is omitted here exactly as in a normal library/CLI configuration.
        for name in (
            "trial_bound",
            "rho_attempts",
            "pm1_attempts",
            "max_input_bits",
        ):
            values.pop(name)
    config = PortfolioConfig(**values)
    workers = QSWorkerConfig(**request["qs_workers"])
    budget = Budget(
        seconds=request["seconds"],
        cpu_seconds=request["cpu_seconds"],
        rss_bytes=8 * 2**30,
    )
    observer = None
    started = time.monotonic()
    try:
        if request.get("observe", True):
            observer = IntegrationObserver(
                owner,
                Path(request["prefix_path"]),
                prefix_bytes=request.get("prefix_bytes", 8 * 2**20),
            )
            with observer:
                run = factorize_bounded(
                    request["n"],
                    seed=request["seed"],
                    config=config,
                    budget=budget,
                    ecm_executor=owner,
                    qs_workers=workers,
                )
                evidence = observer.finish()
                observer._timed("ecm_final_close", owner.close)
                observer.snapshot("after-ecm-close")
        else:
            run = factorize_bounded(
                request["n"],
                seed=request["seed"],
                config=config,
                budget=budget,
                ecm_executor=owner,
                qs_workers=workers,
            )
            owner.close()
            evidence = None
        if run.result.reconstruct() != request["n"]:
            raise ArithmeticError(
                "integration factorization did not reconstruct"
            )
        row = dict(
            role=role,
            n=request["n"],
            reason=run.reason,
            complete=run.result.complete,
            factors=[
                [int(p.value), p.exponent, p.certainty.value]
                for p in run.result.factors
            ],
            remaining=[int(n) for n in run.result.remaining],
            wall=time.monotonic() - started,
            cpu_through_owner_close=budget.cpu_used,
            rss=budget.observed_rss,
            resource_stop=run.resource_stop,
            work=run.work_used,
            events=run.events,
            ecm_execution=getattr(run, "ecm_execution", None),
            portfolio_stage_seconds=run.checkpoint["payload"]
            .get("state", {})
            .get("stage_seconds", {}),
            observation=evidence,
            owner_closed=owner.closed,
        )
        if request.get("require_native"):
            from ....ecm.execution_identity import require_native

            try:
                expected = request.get("expected_workers")
                require_native(
                    run,
                    workers=expected,
                    variant=request.get(
                        "expected_variant",
                        ("scalar" if expected == 1 else "pthread"),
                    ),
                )
            except ValueError as error:
                row["native_assertion_failure"] = str(error)
        return row
    finally:
        owner.close()
        budget.close()


def main():
    try:
        request = json.loads(sys.stdin.readline())
        row = run_request(request)
        print(json.dumps(row), flush=True)
        return 1 if row.get("native_assertion_failure") else 0
    finally:
        close_spawn_helper()


if __name__ == "__main__":
    raise SystemExit(main())
