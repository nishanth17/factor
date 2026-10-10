"""Matched fixed-work and first-factor probes with cooperative cancellation."""

import argparse
import json
import multiprocessing
import os
import queue
import resource
import statistics
import subprocess
import sys
import sysconfig
import threading
import time
from concurrent.futures import (
    ProcessPoolExecutor,
    ThreadPoolExecutor,
    as_completed,
)
from pathlib import Path

from ....common import utils
from ....execution.budget import Budget, BudgetExhaustedError
from ....execution.schedules import SieveContext
from ....execution.stage_jobs import advance_job, new_job
from ....portfolio import PortfolioConfig
from ...suites.phase_one import environment

_STOP = None
_FOUND_AT = None
_CPU_GATE = None


def _initialize(stop, found_at, reports=None, cpu_gate=None):
    """Give each worker the shared cancellation flag at spawn time."""
    global _STOP, _FOUND_AT, _CPU_GATE
    _STOP = stop
    _FOUND_AT = found_at
    _CPU_GATE = cpu_gate
    if reports is not None:
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        reports.put(
            (os.getpid(), rss if sys.platform == "darwin" else rss * 1024)
        )


def _cpu_snapshot():
    """Sample process CPU and require every pool worker to participate."""
    snapshot = (os.getpid(), time.process_time())
    # Blocking ensures a fast worker cannot consume multiple sample tasks,
    # leaving another worker unaccounted for. Snapshot work is not a candidate.
    _CPU_GATE.wait(timeout=30)
    return snapshot


def _worker_cpu(executor, workers):
    """Read a unique CPU counter for each live worker or fail explicitly."""
    futures = [executor.submit(_cpu_snapshot) for _ in range(workers)]
    readings = [future.result() for future in futures]
    if len({pid for pid, _ in readings}) != workers:
        raise RuntimeError("CPU snapshot did not cover every worker")
    return dict(readings)


def _job(assignment):
    """Execute one unique candidate; never give it known input factors."""
    config = PortfolioConfig(
        trial_bound=100,
        max_input_bits=256,
        rho_evaluations=5000,
        segment_size=1024,
        pm1_b1=assignment["b1"],
        pm1_b2=assignment["b2"],
        ecm_tiers=((assignment["b1"], assignment["b2"], 1),),
    )
    candidate = new_job(
        assignment["kind"],
        assignment["n"],
        assignment["seed"],
        assignment["b1"],
        assignment["b2"],
    )
    started = time.perf_counter()
    cpu_started = time.thread_time()
    budget = Budget(
        work_limit=250000,
        seconds=None,
        cpu_seconds=None,
        cancelled=_STOP.is_set if assignment["early"] else None,
    )
    reason = "exhausted"

    try:
        context = SieveContext(config.max_hi, segment_size=config.segment_size)
        while not candidate["done"]:
            advance_job(candidate, budget, context, config)
    except BudgetExhaustedError as error:
        reason = str(error)

    divisor = candidate["factor"]
    if divisor is not None:
        if not utils.valid_divisor(divisor, assignment["n"]):
            raise AssertionError("worker returned an invalid factor")
        reason = "factor"
        if assignment["early"]:
            with _FOUND_AT.get_lock():
                if not _FOUND_AT.value:
                    _FOUND_AT.value = time.monotonic()
            _STOP.set()

    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return {
        "id": assignment["id"],
        "seed": assignment["seed"],
        "pid": os.getpid(),
        "factor": divisor,
        "reason": reason,
        "work": budget.used,
        "seconds": time.perf_counter() - started,
        "cpu_seconds": time.thread_time() - cpu_started,
        "peak_rss_bytes": rss if sys.platform == "darwin" else rss * 1024,
        "finished_at": time.monotonic(),
    }


class RssObserver:
    """Sample aggregate live parent/worker RSS through the OS ps interface.

    The measured peak is a lower bound between samples. Individual process
    high-water marks provide a separate conservative upper bound where every
    worker reports; neither is described as an OS-enforced memory limit.
    """

    def __init__(self, executor, cap_bytes, stop):
        self.executor = executor
        self.cap_bytes = cap_bytes
        self.stop = stop
        self.peak = 0
        self.samples = []
        self.finished = threading.Event()
        self.thread = threading.Thread(target=self._observe, daemon=True)

    def _observe(self):
        """Read numeric process data only; sampler overhead is disclosed."""
        while not self.finished.is_set():
            pids = [os.getpid()]
            if isinstance(self.executor, ProcessPoolExecutor):
                processes = getattr(self.executor, "_processes", None)
                if processes:
                    try:
                        snapshot = list(processes.values())
                    except RuntimeError:
                        snapshot = []
                    pids.extend(process.pid for process in snapshot)

            result = subprocess.run(
                ["ps", "-o", "rss=", "-p", ",".join(map(str, pids))],
                capture_output=True,
                text=True,
                check=False,
            )
            values = [
                int(value) * 1024
                for value in result.stdout.split()
                if value.isdigit()
            ]
            total = sum(values)
            self.peak = max(self.peak, total)
            self.samples.append(total)
            if total > self.cap_bytes:
                self.stop.set()
            self.finished.wait(0.02)

    def start(self):
        """Start observation before submitting candidates."""
        self.thread.start()

    def close(self):
        """Collect the last observer result before destroying the executor."""
        self.finished.set()
        self.thread.join()


def _cpu_usage():
    """Account for terminated workers as well as the parent process."""
    main = resource.getrusage(resource.RUSAGE_SELF)
    children = resource.getrusage(resource.RUSAGE_CHILDREN)
    return (
        main.ru_utime + main.ru_stime + children.ru_utime + children.ru_stime
    )


def _configuration(mode, workers, warm, assignments, args):
    """Use the identical fixed candidate list in every execution mode."""
    context = multiprocessing.get_context("spawn")
    configuration_cpu_started = _cpu_usage()
    stop = context.Event() if mode == "processes" else threading.Event()
    found_at = context.Value("d", 0.0)
    warmup = {"seconds": 0, "calls": 0}
    reports = context.Queue()
    initializer_rss = {}
    cpu_gate = context.Barrier(workers) if mode == "processes" else None

    def create_executor():
        """Spawn isolated processes; the thread arm records the actual GIL."""
        if mode == "serial":
            _initialize(stop, found_at)
            return None
        factory = (
            ProcessPoolExecutor if mode == "processes" else ThreadPoolExecutor
        )
        options = {
            "max_workers": workers,
            "initializer": _initialize,
            "initargs": (stop, found_at, reports, cpu_gate),
        }
        if mode == "processes":
            options["mp_context"] = context
        return factory(**options)

    def execute(executor, jobs):
        """Keep cancelled assignments and validate successful factors."""
        if executor is None:
            return [_job(job) for job in jobs]
        return [
            future.result()
            for future in as_completed(
                [executor.submit(_job, assignment) for assignment in jobs]
            )
        ]

    executor = create_executor() if warm or mode == "serial" else None
    if warm:
        started = time.perf_counter()
        jobs = [{**job, "early": False} for job in assignments]
        while time.perf_counter() - started < args.warmup_seconds:
            execute(executor, jobs)
            warmup["calls"] += len(jobs)
        warmup["seconds"] = time.perf_counter() - started

    samples = []

    try:
        for _ in range(args.repetitions):
            stop.clear()
            found_at.value = 0.0
            started = time.perf_counter()
            monotonic_started = time.monotonic()
            cpu_started = _cpu_usage()
            current = executor
            if current is None and mode != "serial":
                initializer_rss.clear()
                current = create_executor()
            cpu_before = (
                _worker_cpu(current, workers)
                if warm and mode == "processes"
                else {}
            )
            observer = RssObserver(current, args.rss_mib * 2**20, stop)
            observer.start()
            try:
                results = execute(current, assignments)
            finally:
                observer.close()
                if not warm and current is not None:
                    current.shutdown()

            wall = time.perf_counter() - started
            cpu = _cpu_usage() - cpu_started
            if warm and mode == "processes":
                cpu_after = _worker_cpu(current, workers)
                if cpu_after.keys() != cpu_before.keys():
                    raise RuntimeError("pool membership changed during sample")
                cpu = (
                    _cpu_usage()
                    - cpu_started
                    + sum(
                        cpu_after[pid] - value
                        for pid, value in cpu_before.items()
                    )
                )

            results.sort(key=lambda result: result["id"])
            if [row["id"] for row in results] != [
                job["id"] for job in assignments
            ]:
                raise AssertionError(
                    "lost or duplicated candidate assignments"
                )

            while True:
                try:
                    pid, rss = reports.get_nowait()
                except queue.Empty:
                    break
                initializer_rss[pid] = rss

            active_pids = {result["pid"] for result in results}
            # Include late-starting workers even if another worker found the
            # factor before they received a candidate assignment.
            rss_by_pid = dict(initializer_rss)
            active_pids.update(rss_by_pid)
            for result in results:
                pid = result["pid"]
                rss_by_pid[pid] = max(
                    rss_by_pid.get(pid, 0), result["peak_rss_bytes"]
                )

            parent_rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            if sys.platform != "darwin":
                parent_rss *= 1024
            if os.getpid() not in rss_by_pid:
                rss_by_pid[os.getpid()] = parent_rss
            rss_upper = sum(rss_by_pid.values())
            reporting_complete = (
                mode != "processes" or len(active_pids) >= workers
            )
            factor_results = [result for result in results if result["factor"]]
            samples.append(
                {
                    "wall_seconds": wall,
                    "cpu_seconds": cpu,
                    "worker_cpu_start": cpu_before,
                    "worker_cpu_end": cpu_after if cpu_before else {},
                    "work": sum(result["work"] for result in results),
                    "sampled_aggregate_peak_rss_bytes": observer.peak,
                    "sum_process_high_water_rss_bytes": rss_upper,
                    "rss_observations": observer.samples,
                    "cpu_cap_pass": cpu <= args.cpu_seconds,
                    "rss_cap_pass": reporting_complete
                    and (rss_upper <= args.rss_mib * 2**20),
                    "reported_worker_count": len(active_pids),
                    "rss_reporting_complete": reporting_complete,
                    "success": bool(factor_results),
                    "first_factor_seconds": found_at.value - monotonic_started
                    if found_at.value
                    else None,
                    "cancellation_tail_seconds": max(
                        (
                            max(0, result["finished_at"] - found_at.value)
                            for result in results
                            if result["reason"] == "cancelled"
                        ),
                        default=0,
                    ),
                    "results": results,
                }
            )
    finally:
        if executor is not None:
            executor.shutdown()
        reports.close()
        reports.join_thread()

    return {
        "mode": mode,
        "core_budget": workers if mode == "processes" else 1,
        "workers": workers,
        "pool": "warm" if warm else "cold",
        "warmup": warmup,
        "configuration_cpu_seconds": _cpu_usage() - configuration_cpu_started,
        "samples": samples,
        "median_wall_seconds": statistics.median(
            sample["wall_seconds"] for sample in samples
        ),
        "median_cpu_seconds": statistics.median(
            sample["cpu_seconds"] for sample in samples
        ),
    }


def main():
    """Capture feasibility evidence without changing worker defaults."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--jobs", type=int, default=8)
    parser.add_argument("--repetitions", type=int, default=9)
    parser.add_argument("--warmup-seconds", type=float, default=3)
    parser.add_argument("--cpu-seconds", type=float, default=20)
    parser.add_argument("--rss-mib", type=int, default=1024)
    args = parser.parse_args()
    if args.output.exists() or args.jobs < 4 or args.repetitions < 1:
        parser.error("use fresh output, >=4 jobs, and positive repetitions")
    rows = []
    measured_environment = environment()
    measured_environment["core_budget"] = 4

    for kind, n, b1, b2 in (
        ("ecm", 1000000000039 * 1000000000061, 200, 5000),
        ("ecm", 1000000000039 * 1000000000061, 2000, 147396),
        ("ecm", 1009 * 1000000000039, 2000, 147396),
        ("rho", 25013 * 25031, 2, 2),
    ):
        for early in (False, True):
            assignments = [
                {
                    "id": index,
                    "seed": 20261003 + index,
                    "kind": kind,
                    "n": n,
                    "b1": b1,
                    "b2": b2,
                    "early": early,
                }
                for index in range(args.jobs)
            ]
            baseline = None

            for mode, workers, warm in (
                ("serial", 1, True),
                ("threads", 4, True),
                ("processes", 2, False),
                ("processes", 4, False),
                ("processes", 2, True),
                ("processes", 4, True),
            ):
                row = _configuration(mode, workers, warm, assignments, args)

                if not early:
                    signatures = [
                        [result["factor"] for result in sample["results"]]
                        for sample in row["samples"]
                    ]
                    if baseline is None:
                        baseline = signatures[0]
                    if any(signature != baseline for signature in signatures):
                        raise AssertionError(
                            "execution mode changed fixed candidate results"
                        )

                row.update(kind=kind, n=n, b1=b1, b2=b2, early=early)
                rows.append(row)
                print(kind, b1, early, mode, workers, warm, "done", flush=True)

    if environment()["source_sha256"] != measured_environment["source_sha256"]:
        raise RuntimeError("source changed during measurements")
    args.output.write_text(
        json.dumps(
            {
                "milestone": "phase_two_parallel_feasibility",
                "environment": measured_environment,
                "backend": "python-int",
                "gil_enabled": getattr(sys, "_is_gil_enabled", lambda: None)(),
                "gil_build_flag": sysconfig.get_config_var("Py_GIL_DISABLED"),
                "start_method": "spawn",
                "jobs": args.jobs,
                "candidate_work_allowance": 250000,
                "total_work_cap": args.jobs * 250000,
                "cpu_cap_seconds": args.cpu_seconds,
                "rss_cap_bytes": args.rss_mib * 2**20,
                "repetitions": args.repetitions,
                "configurations": rows,
                "decision": "retain serial pending held-out evidence",
                "limitations": [
                    "RSS sampling gaps can miss peaks",
                    "RSS upper bound requires every worker to report",
                    "Observer overhead needs an uninstrumented timing control",
                    "Cancellation latency uses a shared monotonic timestamp",
                    "CPU cap is measured; candidate work caps are enforced",
                    "Warm CPU snapshots include dispatch/IPC between readings",
                    "CPU snapshot instrumentation adds accounting overhead",
                    "Small workloads do not establish broad benefit",
                ],
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
