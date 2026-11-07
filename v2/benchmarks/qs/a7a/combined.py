"""Source-pinned, finite combined-candidate cells and paired diagnostics."""

import argparse
import fcntl
import hashlib
import json
import os
import resource
import signal
import stat
import statistics
import subprocess
import sys
import tempfile
import time
import traceback
from pathlib import Path

from .diagnostics import _events
from .performance_lock import open_lock
from .protocol import INPUTS, schedule
from .runner import (
    host_pressure_status,
    source,
    thermal_status,
    validate_output,
    write,
)

FEATURES = dict(
    lease_budget=True,
    reuse_state=True,
    packed_transport=True,
    verification="worker",
    lookahead=4,
    whole_family=True,
)


def fixtures(kind="combined"):
    if kind not in ("combined", "all", "regimes", "confirmation"):
        raise ValueError("unknown source-pinned input panel")
    data = json.loads((INPUTS / "combined_corpus.json").read_text())
    result = {r["id"]: r for r in data["fixtures"]}
    names = [] if kind == "combined" else ["corpus.json"]
    if kind in ("regimes", "confirmation"):
        names.append("regime_corpus.json")
    if kind == "confirmation":
        names.append("confirmation_g4.json")
    for name in names:
        other = json.loads((INPUTS / name).read_text())
        for row in other["fixtures"]:
            if row["id"] in result:
                if row["n"] != result[row["id"]]["n"]:
                    raise ValueError("conflicting fixture identity")
            else:
                result[row["id"]] = row
    return result


def corpus_identity(kind="combined"):
    """Bind broad searches to both panels without changing old freezes."""
    if kind == "combined":
        return hashlib.sha256(
            (INPUTS / "combined_corpus.json").read_bytes()
        ).hexdigest()
    if kind not in ("all", "regimes", "confirmation"):
        raise ValueError("unknown source-pinned input panel")
    digest = hashlib.sha256()
    names = ["combined_corpus.json", "corpus.json"]
    if kind in ("regimes", "confirmation"):
        names.append("regime_corpus.json")
    if kind == "confirmation":
        names.append("confirmation_g4.json")
    for name in names:
        digest.update(name.encode() + b"\0")
        digest.update((INPUTS / name).read_bytes())
    return digest.hexdigest()


def cell(arm, fixture="exposed_60_0", **changes):
    row = dict(schedule("J0")["cells"][0])
    row.update(
        arm=arm,
        fixture=fixture,
        wall=300,
        cpu=3000,
        stage="J1a",
        tier_a=True,
        features=dict(FEATURES) if arm.startswith("C") else {},
    )
    row.update(changes)
    return row


def storage_bytes(folder):
    """Count captures while atomic metadata replacements may be in flight."""
    size = 0
    for path in folder.rglob("*"):
        try:
            details = path.stat()
        except FileNotFoundError:
            if path.suffix == ".tmp":
                # events.tmp is renamed over events.json atomically. Losing
                # this transient name is expected, not lost capture evidence.
                continue
            raise
        if stat.S_ISREG(details.st_mode):
            size += details.st_size
    return size


def _charged_cpu(receipt, reservation, *, launched):
    """Retain a failed launch's CPU grant when child usage may be missing."""
    observed = receipt["cpu"]
    if launched and receipt["status"] == "failed":
        # Killing a parent before it reaps workers can omit orphan CPU from
        # RUSAGE_CHILDREN. Keep the unknown portion funded, not measured.
        return max(observed, reservation)
    return observed


def _capture(process, capture, elapsed, window, state, *, live_heaps=True):
    owned = _events(capture, process.pid)
    if state["targets"] is None and process.pid in owned:
        workers = [
            p for p, events in owned.items() if events[0]["role"] == "worker"
        ]
        if not state["needs_worker"] or workers:
            state["targets"] = [process.pid] + (
                [min(workers)] if workers else []
            )
    targets = state["targets"] or []
    for target in targets:
        events = owned[target]
        kinds = [e["kind"] for e in events]
        if "recorder_failure" in kinds:
            raise RuntimeError("diagnostic recorder failure")
        if window[0] <= elapsed < window[1] and target not in state["started"]:
            os.kill(target, signal.SIGUSR1)
            state["started"].add(target)
        if (
            elapsed >= window[1]
            and target in state["started"]
            and target not in state["stopped"]
        ):
            os.kill(target, signal.SIGUSR1)
            state["stopped"].add(target)
    # Acknowledged CPU flush precedes sequential live heaps. A deferred signal
    # writes no dump; at most three bounded attempts observe a useful point.
    if (
        live_heaps
        and elapsed >= window[2]
        and all(
            any(e["kind"] == "cpu_stopped" for e in owned[p]) for p in targets
        )
    ):
        for target in targets:
            kinds = [e["kind"] for e in owned[target]]
            if "heap_finished" in kinds:
                continue
            attempts, last = state["heaps"].get(target, (0, 0))
            if attempts == 0 or (
                kinds.count("heap_deferred") == attempts
                and elapsed - last >= 0.5
            ):
                if attempts >= 3:
                    raise RuntimeError("no useful live-heap capture point")
                os.kill(target, signal.SIGUSR2)
                state["heaps"][target] = (attempts + 1, elapsed)
            break


def _capture_full(process, capture, elapsed, window, state):
    """Trace owned processes throughout; request two labeled live heaps."""
    owned = _events(capture, process.pid)
    if process.pid not in owned:
        return
    targets = [process.pid]
    workers = [
        pid for pid, rows in owned.items() if rows[0]["role"] == "worker"
    ]
    if workers:
        targets.append(min(workers))
    if state["targets"] != targets:
        state["targets"] = targets
        temporary = capture / "heap-targets.tmp"
        temporary.write_text(json.dumps(targets))
        temporary.replace(capture / "heap-targets.json")
    if elapsed < window[2]:
        return
    for target in targets:
        kinds = [event["kind"] for event in owned[target]]
        if "heap_finished" in kinds or target in state["heaps"]:
            continue
        try:
            # The private session/group binds these PIDs to this invocation.
            # Never send a heap request to a process outside the owned group.
            if os.getpgid(target) != process.pid:
                continue
            os.kill(target, signal.SIGUSR2)
        except ProcessLookupError:
            continue
        state["heaps"][target] = (1, elapsed)
        break


def run(plan_path, output):
    plan_path, output = plan_path.resolve(), output.resolve()
    plan_raw = plan_path.read_bytes()
    plan = json.loads(plan_raw)
    commit = source()
    if plan["source"] != commit:
        raise ValueError("plan is not bound to current committed source")
    panel = plan.get("corpus_kind", "combined")
    if plan["corpus_sha256"] != corpus_identity(panel):
        raise ValueError("panel identity mismatch")
    if plan.get("analyze_profiles", False):
        # Admission fails before spending a factoring window if the optional
        # analysis reader is missing. It is never sent into worker imports.
        import vmprof  # noqa: F401
    output.mkdir(parents=True, exist_ok=False)
    ledger = dict(
        plan=plan,
        plan_sha256=hashlib.sha256(plan_raw).hexdigest(),
        attempts=[],
    )
    started, cpu_start = time.monotonic(), time.process_time()
    child_start = resource.getrusage(resource.RUSAGE_CHILDREN)
    prior_wall = sum(r["wall"] for r in plan.get("prior_charges", []))
    prior_cpu = sum(r["cpu"] for r in plan.get("prior_charges", []))
    profiles = Path(plan["profile_folder"])
    profiles.mkdir(parents=True, exist_ok=False)
    sealed_capture_bytes = 0
    unobserved_cpu_reserve = 0.0

    def save():
        children = resource.getrusage(resource.RUSAGE_CHILDREN)
        observed_cpu = (
            time.process_time()
            - cpu_start
            + children.ru_utime
            + children.ru_stime
            - child_start.ru_utime
            - child_start.ru_stime
        )
        ledger.update(
            wall=time.monotonic() - started,
            observed_cpu=observed_cpu,
            cpu=observed_cpu + unobserved_cpu_reserve,
            unobserved_cpu_reserve=unobserved_cpu_reserve,
        )
        ledger["cumulative_wall"] = prior_wall + ledger["wall"]
        ledger["cumulative_cpu"] = prior_cpu + ledger["cpu"]
        temporary = output / "ledger.tmp"
        temporary.write_text(json.dumps(ledger, indent=2) + "\n")
        temporary.replace(output / "ledger.json")

    with (
        open_lock() as lock,
        tempfile.TemporaryDirectory(prefix="a7j-") as temporary,
    ):
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        archive = Path(temporary) / "source.tar"
        subprocess.run(
            [
                "git",
                "archive",
                "--format=tar",
                "--output",
                str(archive),
                commit,
            ],
            check=True,
            timeout=20,
        )
        ledger["archive_sha256"] = hashlib.sha256(
            archive.read_bytes()
        ).hexdigest()
        subprocess.run(
            ["tar", "-xf", str(archive), "-C", temporary],
            check=True,
            timeout=20,
        )
        archive.unlink()
        save()
        for ordinal, planned in enumerate(plan["cells"]):
            save()
            spec = dict(planned)
            mode = spec.pop("mode", "clean")
            if mode not in ("clean", "profile", "packet", "cpu"):
                raise ValueError("unknown observation mode")
            packet_counters = spec.get("packet_counters", False)
            if type(packet_counters) is not bool or (
                packet_counters and not spec["arm"].startswith("C")
            ):
                raise ValueError("packet counters need process workers")
            spec["id"] = f"{plan['identity']}-{ordinal:03}"
            reserved_wall, reserved_cpu = spec["wall"], spec["cpu"]
            if spec.get("profile_required", False):
                following = plan["cells"][ordinal + 1 : ordinal + 2]
                if (
                    not following
                    or following[0].get("mode") != "profile"
                    or following[0].get("diagnostic_for") != spec.get("pair")
                ):
                    raise ValueError(
                        "clean run lacks its required diagnostic pair"
                    )
                reserved_wall += following[0]["wall"] + 5
                reserved_cpu += following[0]["cpu"]
            if (
                ledger["cumulative_wall"] + reserved_wall + 5 > plan["wall"]
                or ledger["cumulative_cpu"] + reserved_cpu > plan["cpu"]
            ):
                ledger["stopped"] = (
                    "Next whole reservation exceeds cumulative allowance"
                )
                break
            if (
                sealed_capture_bytes
                + plan.get("prior_storage", 0)
                + (
                    3 * 2**30
                    if mode == "profile" or spec.get("profile_required", False)
                    else 160 * 2**20
                    if mode == "cpu"
                    else 0
                )
                > plan["storage"]
            ):
                ledger["stopped"] = (
                    "Next recorder reservation exceeds capture allowance"
                )
                break
            fixture = fixtures(panel)[spec["fixture"]]
            request = dict(n=fixture["n"], cell=spec)
            folder = output / spec["id"]
            folder.mkdir()
            capture = profiles / spec["id"]
            capture.mkdir()
            request_path, result_path = (
                folder / "request.json",
                folder / "result.json",
            )
            write(request_path, request)
            receipt = dict(
                id=spec["id"],
                arm=spec["arm"],
                fixture=spec["fixture"],
                mode=mode,
                status="issued",
                pair=spec.get("pair"),
                diagnostic_for=spec.get("diagnostic_for"),
                thermal_start=thermal_status(),
                host_pressure_start=host_pressure_status(),
            )
            ledger["attempts"].append(receipt)
            save()
            begin = time.monotonic()
            before = resource.getrusage(resource.RUSAGE_CHILDREN)
            process = None
            state = dict(
                targets=None,
                needs_worker=spec["arm"] != "S0",
                started=set(),
                stopped=set(),
                heaps={},
            )
            reference = [
                r["body"]
                for r in ledger["attempts"]
                if r.get("status") == "reported"
                and r["mode"] == "clean"
                and r["arm"] == spec["arm"]
                and r["fixture"] == spec["fixture"]
                and (
                    "diagnostic_for" not in spec
                    or r.get("pair") == spec["diagnostic_for"]
                )
            ]
            estimate = (
                statistics.median(reference)
                if reference
                else spec.get("profile_reference_seconds", 60)
            )
            if spec.get("full_profile", False):
                # CPU covers the complete execution. Place the single live
                # heap within the repeated workload, including later warmed
                # cache/GC state, rather than always inside its first sample.
                estimate *= spec["samples"]
            window = [
                max(1, 0.20 * estimate),
                max(3, 0.40 * estimate),
                max(4, 0.50 * estimate),
            ]
            if spec.get("full_profile", False):
                matched = next(
                    (
                        r
                        for r in ledger["attempts"]
                        if r.get("pair") == spec.get("diagnostic_for")
                        and r.get("mode") == "clean"
                        and r.get("status") == "reported"
                    ),
                    {},
                )
                window[2] = max(
                    0.1,
                    0.5 * (estimate + matched.get("warm_seconds", 0) + 0.2),
                )
            if mode == "cpu" and "cpu_window" in spec:
                window = spec["cpu_window"]
                if (
                    not isinstance(window, list)
                    or len(window) != 3
                    or not 0 < window[0] < window[1] < spec["wall"] - 5
                    or window[2] != window[1]
                ):
                    raise ValueError("invalid bounded CPU-only window")
            try:
                env = dict(
                    os.environ,
                    PYTHONPATH=temporary,
                    PYTHONDONTWRITEBYTECODE="1",
                )
                env.pop("FACTOR_A7A_PACKET_DIR", None)
                if "nursery_mb" in spec:
                    nursery = spec["nursery_mb"]
                    if type(nursery) is not int or nursery not in (
                        4,
                        16,
                        64,
                        128,
                    ):
                        raise ValueError("unadmitted nursery setting")
                    env["PYPY_GC_NURSERY"] = f"{nursery}MB"
                module = "v2.benchmarks.qs.a7a.cell"
                if mode in ("profile", "cpu"):
                    if packet_counters:
                        packet_folder = capture / "packet-counters"
                        packet_folder.mkdir()
                        env["FACTOR_A7A_PACKET_DIR"] = str(packet_folder)
                    env.update(
                        FACTOR_A7A_DIAGNOSTIC_DIR=str(capture),
                        FACTOR_A7A_LIVE_HEAPS="1"
                        if mode == "profile"
                        else "0",
                    )
                    if spec.get("full_profile", False):
                        env.update(
                            FACTOR_A7A_PROFILE_ALL="1",
                            FACTOR_A7A_LIVE_HEAPS="0",
                        )
                    module = "v2.benchmarks.qs.a7a.diagnostic_hook"
                    if mode == "profile":
                        spec["storage"] = max(spec["storage"], 2 * 2**30)
                    request_path.write_text(json.dumps(request))
                elif mode == "packet":
                    env["FACTOR_A7A_PACKET_DIR"] = str(capture)
                    module = "v2.benchmarks.qs.a7a.packet_observer"
                with (
                    (folder / "stdout.txt").open("x") as stdout,
                    (folder / "stderr.txt").open("x") as stderr,
                ):
                    process = subprocess.Popen(
                        [
                            sys.executable,
                            "-B",
                            "-m",
                            module,
                            "--request",
                            str(request_path),
                            "--result",
                            str(result_path),
                        ],
                        cwd=temporary,
                        env=env,
                        stdout=stdout,
                        stderr=stderr,
                        start_new_session=True,
                    )
                    while process.poll() is None:
                        elapsed = time.monotonic() - begin
                        if (
                            elapsed > spec["wall"]
                            or prior_wall + time.monotonic() - started
                            > plan["wall"]
                        ):
                            raise TimeoutError(
                                "finite cell/stage wall exhausted"
                            )
                        if mode in ("profile", "cpu"):
                            if spec.get("full_profile", False):
                                _capture_full(
                                    process, capture, elapsed, window, state
                                )
                            else:
                                _capture(
                                    process,
                                    capture,
                                    elapsed,
                                    window,
                                    state,
                                    live_heaps=mode == "profile",
                                )
                        if mode in ("profile", "packet", "cpu"):
                            if (
                                sealed_capture_bytes
                                + storage_bytes(capture)
                                + plan.get("prior_storage", 0)
                                > plan["storage"]
                            ):
                                raise RuntimeError(
                                    "recorder storage exhausted"
                                )
                        time.sleep(0.1)
                if process.returncode:
                    raise RuntimeError(f"cell exited {process.returncode}")
                result = json.loads(result_path.read_text())
                validate_output(result, request)
                from .observer import validate_observation

                validate_observation(result)
                if packet_counters and mode in ("profile", "cpu"):
                    from .packet_observer import (
                        summarize_capture,
                        validate_capture,
                    )

                    validate_capture(packet_folder, int(spec["arm"][1:]))
                    receipt["packet_counters"] = summarize_capture(
                        packet_folder,
                        [s["measurement_interval"] for s in result["samples"]],
                    )
                if mode == "packet":
                    from .packet_observer import validate_capture

                    validate_capture(capture, int(spec["arm"][1:]))
                    if any(
                        sample["reason"] == "worker_failure"
                        for sample in result["samples"]
                    ):
                        raise RuntimeError(
                            "diagnostic worker failure; inspect recorder "
                            "before attributing a production failure"
                        )
                for sample in result["samples"]:
                    if (
                        sample["complete"]
                        and "factors" in fixture
                        and sorted(
                            p
                            for p, e, _ in sample["factors"]
                            for _ in range(e)
                        )
                        != sorted(fixture["factors"])
                    ):
                        raise ArithmeticError(
                            "result disagrees with certified factor oracle"
                        )
                if mode in ("profile", "cpu"):
                    owned = _events(capture, process.pid)
                    if spec.get("full_profile", False):
                        workers = [
                            pid
                            for pid, events in owned.items()
                            if events[0]["role"] == "worker"
                        ]
                        state["targets"] = [process.pid] + (
                            [min(workers)] if workers else []
                        )
                        if any(
                            not all(
                                any(e["kind"] == kind for e in events)
                                for kind in ("cpu_started", "cpu_stopped")
                            )
                            for events in owned.values()
                        ):
                            raise RuntimeError(
                                "owned process lacks its full CPU trace"
                            )
                    required = ["cpu_started", "cpu_stopped"]
                    if mode == "profile":
                        required.append("heap_finished")
                    if not state["targets"] or any(
                        not all(
                            any(e["kind"] == kind for e in owned[p])
                            for kind in required
                        )
                        for p in state["targets"]
                    ):
                        raise RuntimeError(
                            "profile lacks required flushed captures"
                        )
                    if mode == "profile" and plan.get(
                        "analyze_profiles", False
                    ):
                        from .capture_review import review

                        analysis = review(
                            capture,
                            result,
                            source=commit,
                            input_id=spec["fixture"],
                            expected_pids=(
                                owned
                                if spec.get("full_profile", False)
                                else state["targets"]
                            ),
                        )
                        receipt["analysis"] = dict(
                            file=str(capture / "analysis.json"),
                            wall=analysis["wall"],
                            cpu=analysis["cpu"],
                        )
                receipt.update(
                    status="reported",
                    body=result["samples"][0]["library_seconds"],
                    cpu_body=result["process_tree_cpu"],
                    complete=result["samples"][0]["complete"],
                    reason=result["samples"][0]["reason"],
                    result=str(result_path),
                    warm_seconds=result["warm_seconds"],
                    targets=state["targets"],
                    profile_window=window
                    if mode in ("profile", "cpu")
                    else None,
                )
            except BaseException as error:
                if process is not None and process.poll() is None:
                    os.killpg(process.pid, signal.SIGTERM)
                    try:
                        process.wait(timeout=2)
                    except subprocess.TimeoutExpired:
                        os.killpg(process.pid, signal.SIGKILL)
                        process.wait(timeout=2)
                receipt.update(
                    status="failed",
                    error=repr(error),
                    error_detail=str(error),
                    error_file=getattr(error, "filename", None),
                    error_traceback=traceback.format_exc(),
                )
                ledger["stopped"] = (
                    "Harness/validation failure: reassess before new issuance"
                )
            finally:
                after = resource.getrusage(resource.RUSAGE_CHILDREN)
                receipt.update(
                    actual_wall=time.monotonic() - begin,
                    cpu=after.ru_utime
                    + after.ru_stime
                    - before.ru_utime
                    - before.ru_stime,
                    thermal_end=thermal_status(),
                    host_pressure_end=host_pressure_status(),
                )
                receipt["cpu_charged"] = _charged_cpu(
                    receipt, spec["cpu"], launched=process is not None
                )
                receipt["unobserved_cpu_reserve"] = (
                    receipt["cpu_charged"] - receipt["cpu"]
                )
                unobserved_cpu_reserve += receipt["unobserved_cpu_reserve"]
                # Finished capture folders are immutable. Revisit only the
                # current folder while recording, avoiding an ever-growing
                # directory scan in every100ms observation tick.
                sealed_capture_bytes += storage_bytes(capture)
                save()
                print(json.dumps(receipt), flush=True)
            if receipt["status"] == "failed":
                break
        index = []
        for capture in profiles.iterdir():
            receipt = next(
                r for r in ledger["attempts"] if r["id"] == capture.name
            )
            request = json.loads(
                (output / capture.name / "request.json").read_text()
            )
            for path in capture.rglob("*"):
                if path.is_file():
                    index.append(
                        dict(
                            path=str(path),
                            bytes=path.stat().st_size,
                            sha256=hashlib.sha256(
                                path.read_bytes()
                            ).hexdigest(),
                            source=commit,
                            arm=receipt["arm"],
                            flags=request["cell"].get("features", {}),
                            collector_overrides=request["cell"].get(
                                "collector_overrides", {}
                            ),
                            siqs_overrides=request["cell"].get(
                                "siqs_overrides", {}
                            ),
                            nursery_mb=request["cell"].get("nursery_mb"),
                            pair=receipt.get("pair"),
                            diagnostic_for=receipt.get("diagnostic_for"),
                            input=receipt["fixture"],
                            repetition=capture.name,
                            capture_point=receipt.get("profile_window"),
                        )
                    )
        write(profiles / "index.json", dict(files=index))
        ledger["capture_bytes"] = storage_bytes(profiles)
        save()
    print(
        json.dumps(
            {k: v for k, v in ledger.items() if k not in ("plan", "attempts")}
        ),
        flush=True,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.plan, args.output)


if __name__ == "__main__":
    main()
