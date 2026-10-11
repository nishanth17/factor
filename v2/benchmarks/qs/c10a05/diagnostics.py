"""C10a-05 attribution and replay evidence, explicitly separate from timing."""

import cProfile
import gc
import json
import pstats
import time
from dataclasses import replace

from ....execution.budget import Budget
from ....qs import SIQSJob
from ..b1.b1_calibration import save
from . import study


def profile(fixture, config, arm, path):
    profiler = cProfile.Profile()
    row = profiler.runcall(study.c9.run_one, fixture, 7, config, arm, 30)
    profiler.dump_stats(str(path.with_suffix(".prof")))
    stats = pstats.Stats(profiler)
    top = sorted(
        stats.stats.items(), key=lambda pair: pair[1][3], reverse=True
    )[:80]
    return dict(
        outcome=row,
        functions=[
            dict(
                file=k[0],
                line=k[1],
                name=k[2],
                own_seconds=v[2],
                cumulative_seconds=v[3],
                calls=v[1],
            )
            for k, v in top
        ],
        meaning="Instrumented attribution; nested times are not additive "
        "or accepted speed evidence",
    )


def checkpoint(fixture, config):
    # This separate diagnostic shares a larger encoding allowance across arms.
    config = replace(config, checkpoint_bytes=16 * 2**20)
    budget = Budget(work_limit=10**13, seconds=60, cpu_seconds=60)
    job = SIQSJob(fixture["n"], seed=7, config=config, budget=budget)
    records = []
    for _ in range(2):
        result = job.run(max_blocks=2)
        if result.divisor is not None:
            if not (
                1 < result.divisor < fixture["n"]
                and result.divisor * result.cofactor == fixture["n"]
            ):
                raise AssertionError("checkpoint diagnostic divisor")
            break
        if result.reason != "paused":
            records.append(dict(reason=result.reason))
            break
        started, work = time.monotonic(), budget.used
        payload = job.checkpoint()
        encoded = json.dumps(payload)
        if len(encoded.encode()) > 16 * 2**20:
            raise MemoryError("C10a-05 diagnostic checkpoint cap")
        encode_seconds, encode_work = (
            time.monotonic() - started,
            budget.used - work,
        )
        identity = job.assignments.identity
        del job, result
        gc.collect()
        started, work = time.monotonic(), budget.used
        job = SIQSJob.from_checkpoint(json.loads(encoded), budget=budget)
        if job.assignments.identity != identity or budget.used <= work:
            raise AssertionError(
                "replay must preserve identity and charge reconstruction"
            )
        records.append(
            dict(
                bytes=len(encoded.encode()),
                encode_seconds=encode_seconds,
                encode_work=encode_work,
                restore_seconds=time.monotonic() - started,
                restore_work=budget.used - work,
                cumulative_work=budget.used,
                stats=job.stats,
            )
        )
        del payload, encoded
    return dict(
        records=records,
        seconds=budget.wall_used,
        cpu_seconds=budget.cpu_used,
        rss_bytes=study.c9._rss_bytes(),
        work=budget.used,
        meaning="Charged replay diagnostic; larger checkpoint allowance, "
        "no timing inference",
    )


def run(output):
    control = study.load_control()
    expense = study.expenses(output)
    fixtures, registry = control["fixtures"], control["configurations"]
    fixture = next(f for f in fixtures if f["digits"] == 50)
    for arm in study.ARMS:
        config = study.c9.decode(registry["50"][arm])
        for kind, cap in (("profile", 32), ("checkpoint", 62)):
            path = output / f"diagnostic-{kind}-{arm}.json"
            if path.exists():
                continue
            result = study.charged(
                expense,
                "diagnostic",
                f"{kind}-{arm}",
                cap,
                lambda: profile(fixture, config, arm, path)
                if kind == "profile"
                else checkpoint(fixture, config),
            )
            save(path, result)
