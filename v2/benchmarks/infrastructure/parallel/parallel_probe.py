"""Early fixed-work probe: serial/thread/process ECM candidate execution.

This measures throughput, including cold spawn and warm worker reuse. It
intentionally completes every assignment; first-factor cancellation and an
aggregate RSS cap must be tested before promoting production parallelism.
"""

import argparse
import json
import multiprocessing
import resource
import statistics
import sys
import sysconfig
import time
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from dataclasses import asdict
from pathlib import Path

from ....common import utils
from ....ecm import core as ecm
from ...suites.phase_one import environment

# Known balanced input; workers receive only n, bounds, and their seed.
INPUT = 1000000000039 * 1000000000061


def _curve_job(assignment):
    seed, b1, b2 = assignment
    work = ecm.EcmStats()
    started = time.perf_counter()
    cpu_started = time.thread_time()
    divisor = ecm.factorize_ecm(
        INPUT, b1=b1, b2=b2, seed=seed, max_curves=1, stats=work
    )
    return {
        "seed": seed,
        "factor": divisor,
        "wall_seconds": time.perf_counter() - started,
        "job_cpu_seconds": time.thread_time() - cpu_started,
        "stats": asdict(work),
        "worker_peak_rss": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    }


def _cpu_usage():
    main = resource.getrusage(resource.RUSAGE_SELF)
    children = resource.getrusage(resource.RUSAGE_CHILDREN)
    return (
        main.ru_utime + main.ru_stime + children.ru_utime + children.ru_stime
    )


def _check(results, assignments):
    if [row["seed"] for row in results] != [job[0] for job in assignments]:
        raise AssertionError("assignment was duplicated, lost, or reordered")
    for row in results:
        if row["factor"] is not None and not utils.valid_divisor(
            row["factor"], INPUT
        ):
            raise AssertionError("invalid worker factor")
        if row["stats"]["curves"] != 1:
            raise AssertionError("worker did not execute exactly one curve")


def _run_configuration(mode, workers, warm, assignments, repetitions):
    samples = []
    rounds = []
    cpu_started = _cpu_usage()
    lifecycle_started = time.perf_counter()
    warmup_seconds = 0
    context = multiprocessing.get_context("spawn")

    def make_executor():
        if mode == "threads":
            return ThreadPoolExecutor(max_workers=workers)
        return ProcessPoolExecutor(max_workers=workers, mp_context=context)

    executor = None
    if warm and mode != "serial":
        executor = make_executor()
        start = time.perf_counter()
        warmup = list(executor.map(_curve_job, assignments[:workers]))
        _check(warmup, assignments[:workers])
        warmup_seconds = time.perf_counter() - start

    try:
        for _ in range(repetitions):
            start = time.perf_counter()
            if mode == "serial":
                results = [_curve_job(job) for job in assignments]
            elif warm:
                results = list(executor.map(_curve_job, assignments))
            else:
                # Cold samples include pool construction and full shutdown.
                with make_executor() as cold_executor:
                    results = list(cold_executor.map(_curve_job, assignments))

            samples.append(time.perf_counter() - start)
            _check(results, assignments)
            rounds.append(results)
    finally:
        if executor is not None:
            executor.shutdown()

    return {
        "mode": mode,
        "workers": workers,
        "pool": "warm" if warm else "cold",
        "samples_seconds": samples,
        "median_seconds": statistics.median(samples),
        "warmup_seconds": warmup_seconds,
        "lifecycle_wall_seconds": time.perf_counter() - lifecycle_started,
        "lifecycle_cpu_seconds": _cpu_usage() - cpu_started,
        "results": rounds,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--jobs", type=int, default=24)
    parser.add_argument("--repetitions", type=int, default=3)
    args = parser.parse_args()
    if args.jobs < 4 or args.repetitions < 1:
        parser.error("at least four jobs and one repetition are required")

    configurations = []

    for b1, b2 in ((200, 5000), (2000, 147396)):
        assignments = [(seed, b1, b2) for seed in range(args.jobs)]
        baseline = None

        for mode, workers, warm in (
            ("serial", 1, False),
            ("threads", 4, False),
            ("processes", 2, False),
            ("processes", 4, False),
            ("processes", 2, True),
            ("processes", 4, True),
        ):
            row = _run_configuration(
                mode, workers, warm, assignments, args.repetitions
            )
            row["b1"], row["b2"] = b1, b2
            signature = [result["factor"] for result in row["results"][0]]
            if baseline is None:
                baseline = signature
                serial_time = row["median_seconds"]
            if signature != baseline:
                raise AssertionError(
                    "execution modes changed candidate results"
                )
            row["serial_over_candidate_time_ratio"] = (
                serial_time / row["median_seconds"]
            )
            configurations.append(row)
            print(
                b1,
                mode,
                workers,
                row["pool"],
                round(row["median_seconds"] * 1000, 3),
                "ms",
            )

    gil_enabled = getattr(sys, "_is_gil_enabled", lambda: None)()
    args.output.write_text(
        json.dumps(
            {
                "environment": environment(),
                "gil_enabled": gil_enabled,
                "free_threaded_build": sysconfig.get_config_var(
                    "Py_GIL_DISABLED"
                ),
                "start_method": "spawn",
                "input": INPUT,
                "jobs_per_sample": args.jobs,
                "repetitions": args.repetitions,
                "worker_rss_units": "bytes on macOS; platform units elsewhere",
                "limitations": [
                    "Fixed-work throughput; no early-factor cancellation test",
                    "Worker RSS peaks are not aggregate concurrent peak RSS",
                    "Warm pools do not yet share prime/curve schedules",
                    "No held-out corpus or default worker count promotion",
                ],
                "configurations": configurations,
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
