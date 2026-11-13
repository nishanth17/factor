"""Retained automatic API protocol with validated warmup and paid final close.

This benchmark adapter also runs against a pinned pre-activation source root.
It never warms an observation labelled cold; use automatic_call.py for cold.
"""

import importlib.util
import json
import os
import subprocess
import sys
import time


def marker():
    from v2.execution.resources import MacResources

    rows = MacResources().sample({os.getpid(): None})
    return dict(
        cpu=sum(r.cpu + r.reaped_cpu for r in rows),
        rss=sum(r.rss for r in rows),
        pids=[r.pid for r in rows],
    )


def thermal(request):
    if not request.get("collect_thermal"):
        return None
    try:
        result = subprocess.run(
            ["pmset", "-g", "therm"], capture_output=True, text=True, timeout=2
        )
        return dict(
            exit=result.returncode,
            stdout=result.stdout[:4096],
            stderr=result.stderr[:4096],
        )
    except (OSError, subprocess.TimeoutExpired) as error:
        return dict(missing=type(error).__name__, detail=str(error))


def state(owner):
    selected = getattr(owner, "native_context", None) or owner
    team = getattr(selected, "team", None)
    engine = getattr(selected, "serial_engine", None)
    context = None if engine is None else engine.programs or engine.context
    packet = getattr(context, "native_packet", None)
    native = bool(team is not None or getattr(packet, "full_ecm", False))
    return dict(
        closed=owner.closed,
        native_owner_retained=native and not owner.closed,
        native_contexts=(
            len(team.packets)
            if team is not None
            else int(packet is not None and packet.context is not None)
        ),
        native_lifecycle=None if team is None else team.lifecycle,
    )


def main():
    first = json.loads(sys.stdin.readline())
    sys.path.insert(0, first["source_root"])
    spec = importlib.util.spec_from_file_location(
        "b4a_automatic_observation_adapter", first["adapter_path"]
    )
    adapter = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(adapter)
    role = first["role"]
    if role == "automatic":
        from v2.ecm.auto import AutomaticECMOwner

        owner = AutomaticECMOwner(seed=first["seed"])
    elif role == "incumbent":
        from v2.ecm.auto import selected_context

        owner = selected_context(seed=first["seed"])
    elif role == "explicit-native":
        from v2.execution.worker_contract import WorkerSettings

        workers = first["workers"]
        settings = WorkerSettings(
            workers=workers,
            allow_wide_workers=True,
            worker_imports="ecm",
            budget_mode="prospective",
        )
        if workers == 1:
            from v2.ecm.workers import ECMExecutor

            owner = ECMExecutor(seed=first["seed"], settings=settings)
        else:
            from v2.execution.ecm_native import NativeECMExecutor

            owner = NativeECMExecutor(
                seed=first["seed"], settings=settings, capacity=32
            )
    else:
        raise ValueError("unknown retained comparison role")
    key = (role, first["seed"], first.get("workers"))
    try:
        request = first
        while True:
            if (
                request["role"],
                request["seed"],
                request.get("workers"),
            ) != key:
                raise ValueError("retained comparison owner changed")
            operation = request["op"]
            if operation == "warm":
                useful, calls = 0.0, []
                before_thermal = thermal(request)
                before = marker()
                started = time.monotonic()
                while useful < 3 and len(calls) < 5000:
                    latest = adapter.run_request(request, retained=owner)
                    calls.append(
                        dict(
                            wall=latest["wall"],
                            valid=latest["valid"],
                            reason=latest.get("reason"),
                            execution=latest.get("execution"),
                        )
                    )
                    useful += latest["wall"]
                    if not latest["valid"] or not latest.get(
                        "endpoint_reached"
                    ):
                        break
                    if time.monotonic() - started > 45:
                        break
                row = dict(
                    op="warm",
                    useful=useful,
                    calls=len(calls),
                    warmup_complete=useful >= 3
                    and latest["valid"]
                    and latest.get("endpoint_reached", False),
                    before=before,
                    after=marker(),
                    evidence=calls,
                    failure=(
                        None
                        if latest["valid"] and latest.get("endpoint_reached")
                        else latest
                    ),
                    state=state(owner),
                    thermal_before=before_thermal,
                    thermal_after=thermal(request),
                )
            elif operation == "call":
                before_thermal = thermal(request)
                before = marker()
                row = adapter.run_request(request, retained=owner)
                row.update(
                    op="call",
                    before=before,
                    after=marker(),
                    state=state(owner),
                    thermal_before=before_thermal,
                    thermal_after=thermal(request),
                )
            elif operation == "close":
                before = marker()
                started = time.monotonic()
                owner.close()
                close_wall = time.monotonic() - started
                row = dict(
                    op="close",
                    wall=close_wall,
                    before=before,
                    after=marker(),
                    state=state(owner),
                )
                print(json.dumps(row), flush=True)
                break
            else:
                raise ValueError("unknown retained operation")
            print(json.dumps(row), flush=True)
            line = sys.stdin.readline()
            if not line:
                break
            request = json.loads(line)
    finally:
        owner.close()
        module = sys.modules.get("multiprocessing.resource_tracker")
        if module is not None:
            tracker = module._resource_tracker
            if tracker._pid is not None:
                tracker._stop()


if __name__ == "__main__":
    main()
