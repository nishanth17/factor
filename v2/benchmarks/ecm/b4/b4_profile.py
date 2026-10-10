"""Separately instrumented baseline costs, never performance evidence."""

import argparse
import cProfile
import fcntl
import json
import os
import pstats
import time
from pathlib import Path

from ..p52.p52_realistic import check_quiet
from . import b4_common as common


def profile(output):
    common.require_runtime()
    protocol, corpus = common.inputs()
    engine = common.load_control()
    package = __import__(engine.__package__, fromlist=["stage_jobs"])
    jobs = package.stage_jobs
    original = jobs.advance_job
    phases = {}
    calls = {}
    components = {}
    active_phase = "outside"
    stack = []

    def timed(function, label):
        def wrapped(*args, **kwargs):
            started = time.perf_counter()
            stack.append(0.0)
            try:
                return function(*args, **kwargs)
            finally:
                elapsed = time.perf_counter() - started
                nested = stack.pop()
                if stack:
                    stack[-1] += elapsed
                key = active_phase + "/" + label
                components[key] = components.get(key, 0) + elapsed - nested

        return wrapped

    for name in (
        "point_add",
        "point_double",
        "scalar_multiply",
        "setup_curve",
    ):
        setattr(package.ecm, name, timed(getattr(package.ecm, name), name))
    for name in ("peek_prime", "gcd", "_batch_check"):
        setattr(jobs, name, timed(getattr(jobs, name), name))

    def advance(job, *args, **kwargs):
        nonlocal active_phase
        phase = job["phase"]
        active_phase = phase
        started = time.perf_counter()
        try:
            return original(job, *args, **kwargs)
        finally:
            phases[phase] = (
                phases.get(phase, 0) + time.perf_counter() - started
            )
            calls[phase] = calls.get(phase, 0) + 1

    # The private package owns every replacement, isolating production.
    jobs.advance_job = advance
    engine.advance_job = advance
    rows = []
    for backend in protocol["backends"]:
        for fixture in corpus["fixtures"]:
            if fixture["split"] != "stages":
                continue
            phases.clear()
            components.clear()
            calls.clear()
            config = common.config(engine, protocol, fixture["case"], backend)
            profiler = cProfile.Profile()
            started = time.perf_counter()
            profiler.enable()
            run = engine.factorize_bounded(
                fixture["n"],
                seed=7,
                config=config,
                budget=engine.Budget(
                    work_limit=protocol["work_limit"],
                    seconds=30,
                    cpu_seconds=30,
                ),
            )
            profiler.disable()
            elapsed = time.perf_counter() - started
            result = common.validate_run(run, fixture)
            stats = pstats.Stats(profiler)
            functions = []
            for (filename, line, name), values in stats.stats.items():
                primitive, total, own, cumulative, _ = values
                if "b4_mainline.json:" in filename:
                    functions.append(
                        dict(
                            file=filename.split("b4_mainline.json:")[1],
                            name=name,
                            calls=total,
                            own_seconds=own,
                            cumulative_seconds=cumulative,
                        )
                    )
            rows.append(
                dict(
                    backend=backend,
                    fixture=fixture["id"],
                    bits=fixture["n"].bit_length(),
                    seconds=elapsed,
                    phases=dict(phases),
                    phase_calls=dict(calls),
                    phase_components=dict(components),
                    functions=functions,
                    result=result,
                )
            )
            print(
                backend,
                fixture["id"],
                round(elapsed, 3),
                dict(phases),
                flush=True,
            )
    with output.open("x") as stream:
        json.dump(
            dict(
                kind="instrumented; overhead included; diagnostic only",
                protocol_sha256=common.digest(common.PROTOCOL),
                rows=rows,
            ),
            stream,
            indent=2,
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    check_quiet({os.getpid()})
    with open("/private/tmp/factor-performance.lock", "a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        owner = Path("/private/tmp/factor-performance-owner.json")
        owner.write_text(
            json.dumps(dict(owner="B4 baseline profile", pid=os.getpid()))
        )
        try:
            profile(args.output)
        finally:
            owner.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
