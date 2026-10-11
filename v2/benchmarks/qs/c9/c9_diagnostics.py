"""Separate bounded profile and charged repeated-resume C9 diagnostics."""

import argparse
import cProfile
import gc
import json
import pstats
import time
from dataclasses import replace
from pathlib import Path

from ....common import utils
from ....execution.budget import Budget
from ....qs import SIQSJob
from ..b1.b1_calibration import save
from . import c9_study as study


def profile(config, fixture, output):
    profiler = cProfile.Profile()
    started, cpu = time.monotonic(), time.process_time()
    row = profiler.runcall(study.run_one, fixture, 7, config, "profile", 60)
    profiler.dump_stats(str(output.with_suffix(".prof")))
    stats = pstats.Stats(profiler)
    top = sorted(
        stats.stats.items(), key=lambda item: item[1][3], reverse=True
    )[:60]
    save(
        output,
        dict(
            outcome=row,
            seconds=time.monotonic() - started,
            cpu_seconds=time.process_time() - cpu,
            functions=[
                dict(
                    file=key[0],
                    line=key[1],
                    name=key[2],
                    primitive_calls=value[0],
                    calls=value[1],
                    own_seconds=value[2],
                    cumulative_seconds=value[3],
                )
                for key, value in top
            ],
            limitation="Instrumented profile; no accepted timing inference",
        ),
    )


def resume(config, fixture, output):
    """Release old stores before two restores; charge encoding and replay."""
    # Checkpoint verification reserves eight simultaneous encoded copies.
    # This is a correctness diagnostic, separately bounded from timing arms.
    config = replace(
        config,
        memory_bytes=512 * 2**20,
        checkpoint_bytes=8 * 2**20,
        collector=replace(config.collector, memory_bytes=512 * 2**20),
    )
    start, cpu = time.monotonic(), time.process_time()
    budget = Budget(work_limit=10**13, seconds=300, cpu_seconds=300)
    job = SIQSJob(fixture["n"], seed=7, config=config, budget=budget)
    receipts = []
    for tier in (10, 20):
        # Pause on public block boundaries without resetting a deadline.
        while budget.wall_used < tier:
            result = job.run(max_blocks=64)
            if result.reason != "paused":
                break
        if result.divisor is not None:
            break
        checkpoint = job.checkpoint()
        encoded = json.dumps(checkpoint)
        if len(encoded.encode()) > 8 * 2**20:
            raise MemoryError("C9 checkpoint byte cap")
        path = output.parent / f"checkpoint-{tier}.json"
        with path.open("x") as stream:
            stream.write(encoded)
        before = budget.used
        saved = dict(checkpoint["resources"])
        del encoded, checkpoint, result, job
        gc.collect()
        checkpoint = json.loads(path.read_text())
        job = SIQSJob.from_checkpoint(checkpoint, config=config, budget=budget)
        del checkpoint
        if budget.used <= before or budget.wall_used < saved["wall_used"]:
            raise AssertionError("resume omitted earlier work or wall time")
        receipts.append(
            dict(
                tier=tier,
                bytes=path.stat().st_size,
                replay_work=budget.used - before,
                cumulative_wall=budget.wall_used,
                cumulative_cpu=budget.cpu_used,
            )
        )
    result = job.run()
    if (result.divisor or 1) * result.cofactor != fixture["n"]:
        raise AssertionError("resume reconstruction failed")
    if result.divisor is not None:
        if not utils.valid_divisor(result.divisor, fixture["n"]):
            raise AssertionError("resume returned improper divisor")
        if sorted((result.divisor, result.cofactor)) != fixture["factors"]:
            raise AssertionError("resume disagrees with independent corpus")
    save(
        output,
        dict(
            seconds=time.monotonic() - start,
            cpu_seconds=time.process_time() - cpu,
            receipts=receipts,
            reason=result.reason,
            divisor=result.divisor,
            remaining=[] if result.divisor else [fixture["n"]],
            reconstructed=(result.divisor or 1) * result.cofactor,
            stats=result.stats,
            work=budget.used,
            config=study.asdict(config),
            limitation="512-MiB/8-MiB repeated-resume correctness only; "
            "not a matched speed estimate",
        ),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("profile", "resume"))
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--arm", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    study.runtime()
    control = json.loads(args.config.read_text())
    config = study.decode(control["configurations"][args.arm])
    fixture = control["fixtures"][0]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with study.machine_window():
        (profile if args.mode == "profile" else resume)(
            config, fixture, args.output
        )


if __name__ == "__main__":
    main()
