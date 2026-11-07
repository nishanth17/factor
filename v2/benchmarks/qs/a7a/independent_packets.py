"""G2H fixed disjoint-packet throughput, without a shared admission graph."""

import argparse
import hashlib
import json
import multiprocessing
import os
import resource
import sys
import time
from dataclasses import replace
from pathlib import Path

from .cell import configuration, settings
from .combined import fixtures

POLYNOMIALS = 64
PACKET = 4
SEEDS = (104729, 130363, 155921)


def partition(tasks, workers):
    """Keep the same distinct ordered packets at every concurrency level."""
    if workers not in (1, 4, 8) or len(tasks) < workers:
        raise ValueError("unsupported fixed-packet allocation")
    return [tasks[index::workers] for index in range(workers)]


def stream_packets(tasks, payload):
    """Retain private state and exact verification, returning only hashes."""
    from ....execution.budget import Budget
    from ....qs.worker_packet import collect_packet
    from ....qs.worker_state import WorkerState
    from ....qs.worker_wire import PackedRows

    started, cpu = time.monotonic(), time.process_time()
    budget = Budget(work_limit=96 * 10**12, seconds=15, cpu_seconds=150)
    runtime = WorkerState(payload, budget)
    writer = PackedRows(payload["settings"].packet_bytes)
    reports = []
    for task in tasks:
        budget.consume(0)
        report = collect_packet(
            dict(task, _runtime=runtime, _transport=writer)
        )
        if report["reason"] != "complete" or report["publication_stop"]:
            raise RuntimeError("independent fixed packet did not complete")
        budget.consume(report["work"])
        reports.append(
            dict(
                index=task["id"],
                sha256=report["sha256"],
                scanned=report["scanned"],
                work=report["work"],
                family_reuses=report["state_stats"]["family_reuses"],
            )
        )
    return dict(
        pid=os.getpid(),
        wall=time.monotonic() - started,
        cpu=time.process_time() - cpu,
        rss_hwm=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        reports=reports,
    )


def _child(connection, tasks, payload, barrier):
    resource.setrlimit(resource.RLIMIT_CPU, (150, 150))
    try:
        barrier.wait(timeout=15)
        connection.send(dict(result=stream_packets(tasks, payload)))
    except BaseException as error:
        connection.send(dict(error=repr(error)))
    finally:
        connection.close()


def summarize(streams, expected_packets):
    reports = sorted(
        (r for stream in streams for r in stream["reports"]),
        key=lambda row: row["index"],
    )
    if [r["index"] for r in reports] != list(range(expected_packets)):
        raise ValueError("missing or duplicate disjoint packet")
    signature = hashlib.sha256(
        json.dumps([(r["index"], r["sha256"]) for r in reports]).encode()
    ).hexdigest()
    return dict(
        signature=signature,
        scanned=sum(r["scanned"] for r in reports),
        work=sum(r["work"] for r in reports),
        family_reuses=sum(r["family_reuses"] for r in reports),
        worker_cpu=sum(s["cpu"] for s in streams),
        longest_stream=max(s["wall"] for s in streams),
    )


def execute(fixture, workers):
    """Run local (zero) or one private stream per spawned owned process."""
    from ....execution.budget import Budget
    from ....qs.current_parallel import CurrentParallelSIQSJob
    from ....qs.worker_state import context_key

    if workers not in (0, 1, 4, 8):
        raise ValueError("unsupported concurrency")
    digits, index = map(
        int, fixture.removeprefix("regime_calibration_").split("_")
    )
    if digits not in (50, 65) or index not in range(3):
        raise ValueError("G2H uses only its calibration inputs")
    n = fixtures("regimes")[fixture]["n"]
    started = time.monotonic()
    config = replace(
        configuration(n, "Y0"),
        memory_bytes=4 * 2**30,
        checkpoint_bytes=64 * 2**20,
    )
    options = settings(
        f"C{workers or 1}",
        22 * 2**30,
        False,
        polynomials_per_task=PACKET,
        max_batch_atoms=2048,
        features=dict(
            lease_budget=True,
            reuse_state=True,
            packed_transport=True,
            verification="worker",
            lookahead=0,
            private_stop=True,
            stop_backend="semaphore",
            coefficient_cache=True,
        ),
    )
    job = CurrentParallelSIQSJob(
        n,
        seed=SEEDS[index],
        config=config,
        workers=options,
        budget=Budget(work_limit=96 * 10**12, seconds=15, cpu_seconds=150),
    )
    job.fixed_work = True
    job._setup()
    if job.divisor is not None:
        raise ValueError("input factors during preparation; no throughput row")
    payload = dict(
        config=config,
        settings=options,
        base=(
            job.source.base.n,
            job.source.base.multiplier,
            job.source.base.bound,
            job.source.base.entries,
        ),
        job=job.identity,
        execution=job.execution_identity,
        work=1
        + len(job.source.base.entries)
        * job.source.base.bound.bit_length() ** 2,
    )
    payload["key"] = context_key(payload)
    tasks = [
        job._task(i, 2 * 10**12, 150, 15) for i in range(POLYNOMIALS // PACKET)
    ]
    prepared = time.monotonic()
    if not workers:
        streams = [stream_packets(tasks, payload)]
    else:
        context = multiprocessing.get_context("spawn")
        barrier = context.Barrier(workers)
        processes, connections, streams = [], [], []
        try:
            for subset in partition(tasks, workers):
                parent, child = context.Pipe(duplex=False)
                process = context.Process(
                    target=_child,
                    args=(child, subset, payload, barrier),
                )
                processes.append(process)
                connections.append(parent)
                process.start()
                child.close()
            deadline = started + 18
            for connection in connections:
                if not connection.poll(max(0, deadline - time.monotonic())):
                    raise TimeoutError("independent actor deadline")
                row = connection.recv()
                if "error" in row:
                    raise RuntimeError(row["error"])
                streams.append(row["result"])
            for process in processes:
                process.join(timeout=max(0, deadline - time.monotonic()))
                if process.exitcode != 0:
                    raise RuntimeError("independent worker exit failure")
        finally:
            for process in processes:
                if process.pid is not None and process.is_alive():
                    process.terminate()
            for process in processes:
                if process.pid is not None:
                    process.join(timeout=0.2)
                    if process.is_alive():
                        process.kill()
                        process.join()
            for connection in connections:
                connection.close()
    result = summarize(streams, len(tasks))
    expected = POLYNOMIALS * (2 * config.half_width + 1)
    if result["scanned"] != expected:
        raise ValueError("fixed polynomial traversal incomplete")
    return dict(
        fixture=fixture,
        workers=workers,
        seed=SEEDS[index],
        polynomials=POLYNOMIALS,
        packet=PACKET,
        setup_wall=prepared - started,
        actor_wall=time.monotonic() - started,
        collection_wall=time.monotonic() - prepared,
        summary=result,
        streams=streams,
        limitations=(
            "Private streams omit shared CPU publication, parent admission "
            "and graph work. Worker verification and packet encoding remain. "
            "Cold startup is included; state reuse varies with allocation. "
            "This is an empirical packet-throughput reference, not a proven "
            "upper bound or complete-factor speedup."
        ),
    )


def main():
    if sys.implementation.name != "pypy" or sys.version_info[:2] != (3, 11):
        raise RuntimeError("PyPy implementing Python 3.11 required")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fixture")
    parser.add_argument("workers", type=int)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise ValueError("never overwrite a throughput receipt")
    result = execute(args.fixture, args.workers)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")


if __name__ == "__main__":
    main()
