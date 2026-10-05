"""Instrument collection separately from validated transport micro-controls.

Profiles identify inclusive costs; their durations are never speed evidence.
Echo controls include dispatch/serialization/wait, not factoring or startup.
They do not isolate the GIL. RSS values are lifetime high-water observations.
"""

import argparse
import cProfile
import json
import pickle
import sys
import time
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

from ..budget import Budget
from ..qs import parallel
from ..qs.relations import verify_atomic
from .performance_audit import CORPUS, fingerprint, measure, verify_corpus


def echo(payload):
    """Return a bounded immutable payload for a measured round trip."""
    return payload


def validate(result, base):
    if result["reason"] != "complete":
        raise AssertionError("transport control needs a completed batch")
    for atom in result["atoms"]:
        verify_atomic(atom, base, residual_bound=10000)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("choose a new capture path")
    if sys.implementation.name != "pypy" or sys.version_info[:2] != (3, 11):
        raise RuntimeError("PyPy implementing Python 3.11 is required")
    corpus = json.loads(CORPUS.read_text())
    verify_corpus(corpus)
    fixture = next(
        f
        for f in corpus["fixtures"]
        if f["split"] == "training" and f["band"] == "medium"
    )
    root = Path(__file__).resolve().parents[2]
    data = dict(
        source_sha256=fingerprint(root),
        fixture=fixture["id"],
        diagnostics={},
        transport={},
    )
    control = parallel.ParallelConfig(
        base_bound=400,
        factor_count=3,
        half_width=512,
        family_count=4,
        pool_size=16,
        max_batch_atoms=1024,
    )

    for width, interval in ((0, 1), (0, 64), (256, 64)):
        config = replace(control, batch_width=width, poll_interval=interval)
        job = parallel.ParallelSIQSJob(
            fixture["n"], config=config, budget=Budget(work_limit=10**12)
        )
        job._setup()
        with parallel.CollectionPool() as pool:
            task = (
                (
                    job.base.n,
                    job.base.multiplier,
                    job.base.bound,
                    job.base.entries,
                ),
                job._assignment_primes(0),
                0,
                config,
                50_000_000,
                None,
                True,
            )

            def collect():
                result = parallel._collect(task, pool.state)
                validate(result, job.base)
                return result

            started = time.perf_counter()
            while time.perf_counter() - started < 3:
                collect()
            profile = cProfile.Profile()
            result = profile.runcall(collect)
            name = f"width{width}_poll{interval}"
            profile.dump_stats(str(args.output.with_suffix(f".{name}.prof")))
            data["diagnostics"][name] = dict(
                atoms=len(result["atoms"]),
                scanned=result["scanned"],
                work=result["work"],
                workspace=result["workspace_bytes"],
                pickle_bytes=len(pickle.dumps(result, protocol=5)),
            )
            if width == 256:
                payload = result
                base = job.base

    settings = SimpleNamespace(warmup_seconds=3, repetitions=9)

    for mode, workers in (
        ("serial", 1),
        ("thread", 1),
        ("thread", 2),
        ("thread", 4),
        ("process", 1),
        ("process", 2),
        ("process", 4),
    ):
        with parallel.CollectionPool(mode, workers) as pool:

            def roundtrip():
                if pool.executor is None:
                    results = [
                        pickle.loads(pickle.dumps(payload, protocol=5))
                        for _ in range(4)
                    ]
                else:
                    futures = [
                        pool.executor.submit(echo, payload) for _ in range(4)
                    ]
                    results = [future.result() for future in futures]

                for result in results:
                    # Equality checks the full payload; separate verification
                    # follows timings to avoid attributing it to transport.
                    if result != payload:
                        raise AssertionError("transport corrupted the batch")

                return []

            data["transport"][f"{mode}_{workers}"] = measure(
                roundtrip, settings
            )

    validate(payload, base)
    if data["source_sha256"] != fingerprint(root):
        raise AssertionError("source changed during diagnostic capture")
    args.output.write_text(json.dumps(data, indent=2) + "\n")


if __name__ == "__main__":
    main()
