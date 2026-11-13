"""One cold default-path call, preserving unexpected fallback observations.

Run as a file with one JSON request on stdin. source_root may point to a
committed control checkout; imports are delayed until after that root is set.
The ordinary automatic arm never constructs or preloads a native owner.
"""

import json
import sys
import time


def run_request(request):
    started, cpu_started = time.monotonic(), time.process_time()
    from v2.execution.budget import Budget
    from v2.portfolio import PortfolioConfig, factorize_bounded

    values = dict(request.get("config", {}))
    if "ecm_tiers" in values:
        values["ecm_tiers"] = tuple(map(tuple, values["ecm_tiers"]))
    config = PortfolioConfig(**values)
    budget = Budget(
        seconds=request["seconds"],
        cpu_seconds=request["cpu_seconds"],
        rss_bytes=request.get("rss_bytes", 8 * 2**30),
    )
    owner = None
    role = request.get("role", "automatic")
    row = dict(
        role=role,
        source=request.get("source"),
        n=request["n"],
        seed=request["seed"],
        requested_kernel=values.get("ecm_kernel", "omitted"),
        numerical_tiers=[
            list(t) for t in config.for_input(request["n"]).ecm_tiers
        ],
        require_native=bool(request.get("require_native", False)),
        valid=False,
    )
    try:
        executor = "auto"
        if role == "serial":
            executor = None
        elif role == "explicit-native":
            from v2.execution.worker_contract import WorkerSettings

            workers = request["workers"]
            settings = WorkerSettings(
                workers=workers,
                allow_wide_workers=workers > 4,
                worker_imports="ecm",
                budget_mode="prospective",
            )
            if workers == 1:
                from v2.ecm.workers import ECMExecutor

                owner = ECMExecutor(seed=request["seed"], settings=settings)
            else:
                from v2.execution.ecm_native import NativeECMExecutor

                owner = NativeECMExecutor(
                    seed=request["seed"],
                    settings=settings,
                    capacity=32,
                    schedule="dynamic",
                    chunk_size=1,
                )
            executor = owner
        elif role not in ("automatic", "incumbent"):
            raise ValueError("unknown automatic comparison role")
        run = factorize_bounded(
            request["n"],
            seed=request["seed"],
            config=config,
            budget=budget,
            ecm_executor=executor,
            qs_workers=None,
            stop_after_split=request.get("stop_after_split", True),
        )
        row.update(
            reason=run.reason,
            complete=run.result.complete,
            reconstructed=run.result.reconstruct() == request["n"],
            factors=[
                [int(p.value), p.exponent, p.certainty.value]
                for p in run.result.factors
            ],
            remaining=[int(n) for n in run.result.remaining],
            work=run.work_used,
            resource_stop=run.resource_stop,
            execution=getattr(run, "ecm_execution", None),
            events=run.events,
            checkpoint=run.checkpoint
            if request.get("save_checkpoint")
            else None,
        )
        if row["execution"] is None:
            # A committed pre-activation control has no new receipt schema.
            # Its unchanged events establish Python serial/process execution;
            # absent worker counts remain unknown rather than inferred.
            jobs = []
            for event in run.events:
                if event["stage"] not in ("ecm", "ecm_workers"):
                    continue
                stats = event.get("stats", {})
                jobs.append(
                    dict(
                        backend="python",
                        implementation="python",
                        topology="scalar"
                        if event["stage"] == "ecm"
                        else "processes",
                        workers=1
                        if event["stage"] == "ecm"
                        else stats.get("configured_workers"),
                        native_entries=0,
                        event=event,
                    )
                )
            row["execution"] = dict(
                schema="historical-python-events-v1",
                requested="omitted",
                backend="python" if jobs else "not-engaged",
                native_entries=0,
                jobs=jobs,
                reason="historical_committed_default",
                build=None,
            )
        if not row["reconstructed"]:
            raise ArithmeticError("automatic comparison did not reconstruct")
        # Save every paid observation before applying a non-forcing assertion.
        if request.get("observation_path"):
            _save(request["observation_path"], row)
        if request.get("require_native"):
            from v2.ecm.execution_identity import require_native

            require_native(run, workers=request.get("expected_workers"))
        row["valid"] = True
    except Exception as error:
        row.update(failure=type(error).__name__, detail=str(error))
    finally:
        if owner is not None:
            owner.close()
        row.update(
            wall=time.monotonic() - started,
            parent_cpu=time.process_time() - cpu_started,
            aggregate_cpu=budget.cpu_used,
            observed_rss=budget.observed_rss,
        )
        budget.close()
        if request.get("observation_path"):
            _save(request["observation_path"], row)
    return row


def _save(path, row):
    from pathlib import Path

    Path(path).write_text(json.dumps(row, sort_keys=True) + "\n")


def main():
    request = json.loads(sys.stdin.readline())
    if request.get("source_root"):
        sys.path.insert(0, request["source_root"])
    row = run_request(request)
    print(json.dumps(row), flush=True)
    return 0 if row["valid"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
