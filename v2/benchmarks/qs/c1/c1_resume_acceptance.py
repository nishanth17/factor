"""One finite DLP checkpoint/resume acceptance, outside timing arms."""

import argparse
import gc
import hashlib
import json
import platform
import time
from collections import Counter
from dataclasses import asdict, replace
from math import prod
from pathlib import Path

from ....common import utils
from ....execution.budget import Budget, BudgetExhaustedError
from ....qs import SIQSJob
from ...support.paths import (
    source_path,
)
from ..b1.b1_calibration import save
from ..phase_three.phase_three_reference import _rss_bytes
from .c1_feasibility import machine_window
from .c1_followup import HERE, load_followup
from .c1_implementation import arm_config, source_hashes

CONTROL = HERE / "inputs/controls/c1_resume_acceptance.json"


def run(output):
    control = json.loads(source_path(CONTROL).read_text())
    _, fixtures, configs = load_followup()
    fixture = next(
        f
        for f in fixtures
        if f["digits"] == control["digits"]
        and f["c1_index"] == control["training_index"]
    )
    config = arm_config(configs[control["digits"]], control["arm"])
    config = replace(
        config,
        memory_bytes=control["memory_bytes"],
        checkpoint_bytes=control["checkpoint_bytes"],
        collector=replace(
            config.collector, memory_bytes=control["memory_bytes"]
        ),
    )
    save(
        output / "manifest.json",
        dict(
            control=control,
            control_sha256=hashlib.sha256(
                source_path(CONTROL).read_bytes()
            ).hexdigest(),
            driver_sha256=hashlib.sha256(
                source_path(Path(__file__)).read_bytes()
            ).hexdigest(),
            source=source_hashes(),
            fixture=fixture,
            config=asdict(config),
        ),
    )
    started, cpu = time.monotonic(), time.process_time()
    prior = Budget(
        work_limit=control["work"],
        seconds=control["pause_seconds"],
        cpu_seconds=control["pause_seconds"],
    )
    job = SIQSJob(
        fixture["n"], seed=control["seed"], config=config, budget=prior
    )
    first = job.run()
    assert (first.divisor or 1) * first.cofactor == fixture["n"]
    first_stats = json.loads(json.dumps(first.stats))
    checkpoint = job.checkpoint()
    assert checkpoint["version"] == 4
    checkpoint_path = output / "checkpoint.json"
    save(checkpoint_path, checkpoint)
    assert checkpoint_path.stat().st_size <= control["retained_byte_limit"]
    saved_resources = dict(checkpoint["resources"])
    del checkpoint, first, job
    gc.collect()

    # Charge encoding, disk roundtrip and collection of the released old job.
    # No old live atom store is intentionally retained beside the restored one.
    checkpoint = json.loads(source_path(checkpoint_path).read_text())
    budget = Budget(
        work_limit=control["work"],
        seconds=control["cumulative_seconds"],
        cpu_seconds=control["cumulative_seconds"],
        used=prior.used,
        prior_wall=prior.wall_used,
        prior_cpu=prior.cpu_used,
    )
    restored = SIQSJob.from_checkpoint(
        checkpoint, config=config, budget=budget
    )
    assert budget.used > saved_resources["work_used"]
    assert budget.prior_wall >= saved_resources["wall_used"]
    assert budget.prior_cpu >= saved_resources["cpu_used"]
    restore_work = budget.used - saved_resources["work_used"]
    del checkpoint
    result = restored.run()
    assert (result.divisor or 1) * result.cofactor == fixture["n"]
    assert result.divisor is None or utils.valid_divisor(
        result.divisor, fixture["n"]
    )
    factors, certainty, remaining = [], [], [fixture["n"]]
    if result.divisor is not None:
        remaining = [result.divisor, result.cofactor]
        try:
            for child in tuple(remaining):
                budget.consume(child.bit_length() ** 2)
                label = utils.classify_prime(child)
                if label is utils.Primality.COMPOSITE:
                    break
                factors.append(child)
                certainty.append(label.value)
                remaining.pop(0)
        except BudgetExhaustedError:
            pass
    assert prod(factors) * prod(remaining) == fixture["n"]
    assert not Counter(factors) - Counter(fixture["factors"])
    assert result.stats["workspace_bytes"] <= control["memory_bytes"]
    assert budget.used <= control["work"]
    seconds = time.monotonic() - started
    assert seconds <= control["harness_seconds"]
    save(
        output / "result.json",
        dict(
            n=fixture["n"],
            first_stats=first_stats,
            stats=result.stats,
            complete=not remaining,
            factors=factors,
            certainty=certainty,
            divisor=result.divisor,
            remaining=remaining,
            reason=result.reason,
            seconds=seconds,
            cpu_seconds=time.process_time() - cpu,
            work=budget.used,
            restore_work=restore_work,
            saved_resources=saved_resources,
            checkpoint_bytes=checkpoint_path.stat().st_size,
            rss_bytes=_rss_bytes(),
            classification=(
                "correctness/resume only; larger declared memory "
                "and deadline than timed arms"
            ),
        ),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if (
        platform.python_implementation() != "PyPy"
        or platform.python_version_tuple()[:2] != ("3", "11")
    ):
        raise RuntimeError("C1 requires PyPy implementing Python 3.11")
    args.output.mkdir(parents=True, exist_ok=True)
    with machine_window():
        run(args.output)


if __name__ == "__main__":
    main()
