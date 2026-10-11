"""C10 attribution and replay evidence, explicitly separate from timing."""

import cProfile
import gc
import json
import pstats
import time
from dataclasses import replace

from ....execution.budget import Budget
from ....qs import SIQSJob, build_factor_base
from ....qs.assignment_stream import AssignmentStream
from ..b1.b1_calibration import save
from . import study
from .oracles import product_census


def census(fixtures, registry):
    records = []
    budget = Budget(work_limit=10**9, seconds=20, cpu_seconds=20)
    for fixture in fixtures:
        config = study.c9.decode(registry[str(fixture["digits"])]["control"])
        base = build_factor_base(
            fixture["n"],
            bound=config.base_bound,
            multiplier=config.multiplier,
            budget=budget,
        ).factor_base
        for policy in ("flyer", "quality"):
            stream = AssignmentStream(
                base,
                config.half_width,
                factor_count=config.factor_count,
                family_count=128,
                pool_size=config.pool_size,
                seed=7,
                policy=policy,
                budget=budget,
                memory_bytes=config.memory_limit,
            )
            families = [stream[i] for i in range(min(128, len(stream)))]
            records.append(
                dict(
                    id=fixture["id"],
                    policy=policy,
                    prefixes={
                        str(size): product_census(
                            families[:size], stream.target
                        )
                        for size in (4, 8, 16, 64, 128)
                        if size <= len(families)
                    },
                )
            )
    return dict(
        records=records,
        work=budget.used,
        meaning="Exact product/overlap diagnostics, "
        "not performance or useful-yield evidence",
    )


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
            raise MemoryError("C10 diagnostic checkpoint cap")
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
    path = output / "diagnostic-products.json"
    if not path.exists():
        result = study.charged(
            expense,
            "diagnostic",
            "products",
            22,
            lambda: census(fixtures, registry),
        )
        save(path, result)
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
