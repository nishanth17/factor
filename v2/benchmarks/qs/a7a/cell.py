"""One-shot experiment CLI; supervisor owns source pin, lock and total caps."""

import argparse
import hashlib
import json
import resource
import sys
import time
from dataclasses import asdict, replace
from math import prod
from pathlib import Path

from .protocol import INPUTS


def configuration(n, policy):
    from ....qs.defaults import default_siqs_config
    from ....qs.worker_configuration import decode_config

    digits = len(str(n))
    if policy == "G0" or digits < 50:
        return default_siqs_config(n, memory_bytes=2**30)
    templates = json.loads((INPUTS / "policies.json").read_text())["policies"]
    config = decode_config(templates["60" if digits <= 65 else "70"])
    if digits in (60, 70):
        return config
    counts = {50: 1800, 55: 2500, 62: 6000, 65: 8500}
    factors = {50: 6, 55: 7, 62: 9, 65: 9, 75: 11, 80: 12}
    if digits in counts:
        return replace(
            config, base_count=counts[digits], factor_count=factors[digits]
        )
    from ....qs.siqs import SIQSConfig

    values = asdict(config)
    values.pop("base_count")
    values.update(
        base_bound=400000 if digits == 75 else 600000,
        factor_count=factors[digits],
        collector=config.collector,
    )
    return SIQSConfig(**values)


def settings(
    arm,
    memory,
    observe,
    polynomials_per_task=4,
    max_batch_atoms=8192,
    features=None,
):
    from ....qs.current_parallel import QSWorkerConfig

    if arm == "S0":
        if features:
            raise ValueError("serial control cannot carry worker features")
        return None
    features = {} if features is None else dict(features)
    allowed = {
        "lease_budget",
        "reuse_state",
        "packed_transport",
        "verification",
        "verification_batch_atoms",
        "lookahead",
        "whole_family",
        "arrival_order",
        "state_cache_bytes",
        "private_stop",
        "cpu_snapshot",
        "stop_backend",
        "coefficient_cache",
        "worker_memory_bytes",
        "packet_bytes",
    }
    if features.keys() - allowed:
        raise ValueError("unknown experimental worker features")
    return QSWorkerConfig(
        mode="serial" if arm == "M0" else "process",
        workers=1 if arm == "M0" else int(arm[1:]),
        memory_bytes=memory,
        observe=observe,
        polynomials_per_task=polynomials_per_task,
        max_batch_atoms=max_batch_atoms,
        **features,
    )


def validate_result(result, n):
    factors = result.get("factors", [])
    remaining = result.get("remaining", [])
    if prod(int(p) ** int(e) for p, e, _ in factors) * prod(remaining) != n:
        raise ArithmeticError("complete result fails exact reconstruction")
    if result["complete"] != (not remaining):
        raise ArithmeticError("completion metadata disagrees with cofactors")
    for p, e, certainty in factors:
        if (
            p < 2
            or e < 1
            or certainty not in ("proven_prime", "probable_prime")
        ):
            raise ArithmeticError("invalid terminal factor/certainty")


def prefix(
    n, config, worker_settings, budget, tasks, pool, *, polynomials=None
):
    from ....qs.siqs import SIQSJob

    if polynomials is not None:
        if type(polynomials) is not int or not 1 <= polynomials <= 8192:
            raise ValueError(
                "finite positive prefix polynomial quota required"
            )
        if polynomials > config.family_count * config.gray_limit:
            raise ValueError("prefix quota exceeds the available traversal")
        if worker_settings is not None:
            width = (
                config.gray_limit
                if worker_settings.whole_family
                else worker_settings.polynomials_per_task
            )
            families, tail = divmod(polynomials, config.gray_limit)
            if tail % width:
                raise ValueError("prefix quota splits a worker packet")
            chunks = (config.gray_limit + width - 1) // width
            tasks = families * chunks + tail // width

    if worker_settings is None:
        job = SIQSJob(n, seed=budget.a7a_seed, config=config, budget=budget)
        blocks = (
            2 * config.half_width + config.batch_width
        ) // config.batch_width
        count = tasks * 4 if polynomials is None else polynomials
        result = job.run(max_blocks=count * blocks)
    else:
        from ....qs.current_parallel import CurrentParallelSIQSJob

        job = CurrentParallelSIQSJob(
            n,
            seed=budget.a7a_seed,
            config=config,
            workers=worker_settings,
            budget=budget,
        )
        result = job.run(pool=pool, max_assignments=tasks, fixed_work=True)
    atoms = (
        []
        if job.engine is None
        else list(job.engine.collector._atoms.values())
    )
    signature = hashlib.sha256()
    for atom in sorted(atoms, key=lambda a: a.relation_id):
        signature.update(
            str(
                (
                    atom.relation_id,
                    atom.sign,
                    tuple(atom.exponents),
                    atom.residual,
                    atom.large_primes,
                    getattr(atom, "opaque", None),
                )
            ).encode()
        )
    if result.divisor is not None and not (
        1 < result.divisor < n and n % result.divisor == 0
    ):
        raise ArithmeticError("prefix divisor fails exact validation")
    return dict(
        complete=False,
        factors=[],
        remaining=[n],
        reason=result.reason,
        divisor=result.divisor,
        stats=result.stats,
        signature=signature.hexdigest(),
        endpoint="fixed_work_collection"
        if worker_settings is not None
        else "serial_prefix",
    )


def complete(
    n,
    config,
    worker_settings,
    budget,
    pool,
    endpoint,
    *,
    memory_bytes=16 * 2**30,
):
    from ....portfolio import PortfolioConfig, factorize_bounded
    from ....qs.defaults import uses_siqs_defaults

    if endpoint == "portfolio":
        # Preserve the unconditioned production dispatch population. A bigger
        # global worker grant must not retune the child's collector policy.
        child = configuration(n, "G0") if uses_siqs_defaults(n) else None
        options = PortfolioConfig(siqs=child, memory_bytes=memory_bytes)
    else:
        options = PortfolioConfig(
            siqs=config,
            rho_attempts=0,
            pm1_attempts=0,
            ecm_tiers=(),
            memory_bytes=memory_bytes,
        )
    run = factorize_bounded(
        n,
        seed=budget.a7a_seed,
        config=options,
        budget=budget,
        qs_workers=worker_settings,
        qs_pool=pool,
        stop_after_split=endpoint == "first",
    )
    state = run.checkpoint["payload"]["state"]
    current = state["current"]
    result = dict(
        complete=run.result.complete,
        factors=[
            (f.value, f.exponent, f.certainty.value)
            for f in run.result.factors
        ],
        remaining=list(run.result.remaining),
        reason=run.reason,
        events=list(run.events),
        dropped_events=run.dropped_events,
        portfolio_stage_seconds=state["stage_seconds"],
        active_stage=None if current is None else current["stage"],
        endpoint=endpoint,
    )
    validate_result(result, n)
    return result


def complete_call(n, config, workers, budget, pool, spec, *, seconds_cap=None):
    """Use the real entry point and reconcile each warm/sample child ledger."""
    from ....execution.budget import Budget

    seconds = max(0, budget.seconds - budget.wall_used)
    if seconds_cap is not None:
        seconds = min(seconds, seconds_cap)
    child = Budget(
        work_limit=budget.work_limit - budget.used,
        seconds=seconds,
        cpu_seconds=max(0, budget.cpu_seconds - budget.cpu_used),
    )
    child.a7a_seed = spec["seed"]
    child.cancelled = budget.cancelled
    try:
        return complete(
            n,
            config,
            workers,
            child,
            pool,
            spec["endpoint"],
            memory_bytes=spec["memory"],
        )
    finally:
        # Parent CPU is already on the encompassing process clock. Publish
        # only this child's worker CPU once, including on exceptional exit.
        budget.used += child.used
        budget.prior_cpu += child.prior_cpu


def execute(request):
    from ....execution.budget import Budget

    if sys.implementation.name != "pypy" or sys.version_info[:2] != (3, 11):
        raise RuntimeError("PyPy implementing Python 3.11 is required")
    n, spec = request["n"], request["cell"]
    warm_target = spec.get("warmup_seconds", 3)
    if type(warm_target) is not int or warm_target not in (3, 10, 20):
        raise ValueError("unadmitted warmup duration")
    if spec["lifetime"] == "warm" and spec["wall"] <= warm_target + 10:
        raise ValueError("warmup must leave timed work and cleanup funded")
    config = configuration(n, spec["policy"])
    if "checkpoint_bytes" in spec:
        # A complete pilot needs a funded resume snapshot beyond tiny prefixes;
        # this grant changes no base, graph, queue or recovery selection.
        config = replace(config, checkpoint_bytes=spec["checkpoint_bytes"])
    if "qs_memory_bytes" in spec:
        config = replace(config, memory_bytes=spec["qs_memory_bytes"])
    numerical = spec.get("siqs_overrides", {})
    numerical_fields = {
        "base_count",
        "base_bound",
        "factor_count",
        "half_width",
        "max_half_width",
        "polynomials_per_family",
        "pool_size",
        "batch_width",
    }
    if type(numerical) is not dict or numerical.keys() - numerical_fields:
        raise ValueError("unknown experimental numerical overrides")
    if numerical:
        config = replace(config, **numerical)
    overrides = spec.get("collector_overrides", {})
    allowed = {
        "marking",
        "sparse_queue_bytes",
        "power_plan_bytes",
        "metadata_chunk",
        "batch_policy",
        "batch_size",
        "opaque_cache_entries",
        "candidate_bound",
        "large_prime_bound",
        "small_prime_bound",
        "large_product_bound",
        "residual_bound",
        "split_call_limit",
        "max_atoms",
        "max_partials",
        "max_relations",
    }
    if type(overrides) is not dict or overrides.keys() - allowed:
        raise ValueError("unknown experimental collector overrides")
    if overrides:
        config = replace(
            config, collector=replace(config.collector, **overrides)
        )
    startup_gc = None
    if "nursery_mb" in spec:
        import gc
        import os
        import re

        nursery = spec["nursery_mb"]
        if type(nursery) is not int or nursery not in (4, 16, 64, 128):
            raise ValueError("unadmitted nursery setting")
        if os.environ.get("PYPY_GC_NURSERY") != f"{nursery}MB":
            raise ValueError("nursery must be selected before process startup")
        startup_gc = str(gc.get_stats())
        observed = re.search(r"nursery:\s+([\d.]+)MB", startup_gc)
        if observed is None or float(observed.group(1)) != nursery:
            raise ValueError(
                "runtime nursery differs from the requested setting"
            )
    workers = settings(
        spec["arm"],
        spec["memory"],
        spec["observe"],
        spec.get("polynomials_per_task", 4),
        spec.get("max_batch_atoms", 8192),
        spec.get("features"),
    )
    # Imports, warmup, samples, validation and teardown share the supervisor's
    # hard reservation. Inner work/time stop leaves a ten-second cleanup tail.
    budget = Budget(
        work_limit=spec["work"],
        seconds=max(0, spec["wall"] - 10),
        cpu_seconds=max(0, spec["cpu"] - 10),
    )
    budget.a7a_seed = spec["seed"]
    if workers is not None:
        from ....qs.parallel import CollectionPool

    pool = (
        None
        if workers is None
        else CollectionPool(workers.mode, workers.workers)
    )
    observer = None
    if spec["observe"]:
        from .observer import OwnedRSS

        observer = OwnedRSS(pool, limit=spec["memory"])
        observer.__enter__()
        budget.cancelled = (
            lambda: observer.exceeded or observer.error is not None
        )
    from .tier_a import TierA

    counters = TierA() if spec.get("tier_a", False) else None
    if counters is not None:
        counters.install()
    samples, warm_seconds, warm_work = [], 0.0, 0
    try:
        if spec["lifetime"] == "warm":
            started, before = time.monotonic(), budget.used
            while time.monotonic() - started < warm_target:
                budget.consume(0)
                if spec["endpoint"] == "prefix":
                    prefix(
                        n,
                        config,
                        workers,
                        budget,
                        max(1, workers.workers if workers else 1),
                        pool,
                    )
                else:
                    # The portfolio adds coexistence ownership to its QS
                    # settings. Warm that same entry point, not a bare prefix
                    # with another pool identity. A bounded partial call is
                    # sufficient; warmup need not solve the input first.
                    complete_call(
                        n,
                        config,
                        workers,
                        budget,
                        pool,
                        spec,
                        seconds_cap=max(
                            0, warm_target - (time.monotonic() - started)
                        ),
                    )
            warm_seconds, warm_work = (
                time.monotonic() - started,
                budget.used - before,
            )
        for _ in range(spec["samples"]):
            budget.consume(0)
            if counters is not None:
                counters.reset()
            interval_start = time.monotonic()
            started = time.perf_counter()
            if spec["endpoint"] == "prefix":
                result = prefix(
                    n,
                    config,
                    workers,
                    budget,
                    spec["tasks"],
                    pool,
                    polynomials=spec.get("polynomials"),
                )
            else:
                result = complete_call(n, config, workers, budget, pool, spec)
            result["library_seconds"] = time.perf_counter() - started
            result["measurement_interval"] = [
                interval_start,
                time.monotonic(),
            ]
            if counters is not None:
                result["tier_a"] = counters.report(result["library_seconds"])
            samples.append(result)
    finally:
        if observer is not None:
            observer.__exit__()
        if pool is not None:
            pool.close()
        if counters is not None:
            counters.close()
    children = resource.getrusage(resource.RUSAGE_CHILDREN)
    own = resource.getrusage(resource.RUSAGE_SELF)
    return dict(
        cell=spec,
        n=n,
        samples=samples,
        warm_seconds=warm_seconds,
        warm_work=warm_work,
        warm_endpoint=spec["endpoint"],
        warmup_target=warm_target,
        work=budget.used,
        charged_cpu=budget.cpu_used,
        process_tree_cpu=own.ru_utime
        + own.ru_stime
        + children.ru_utime
        + children.ru_stime,
        rss_self_hwm=own.ru_maxrss,
        rss_children_hwm=children.ru_maxrss,
        memory_note="Lifetime HWM fields are not simultaneous aggregate RSS.",
        configuration=asdict(config),
        startup_gc=startup_gc,
        rss_samples=[] if observer is None else observer.rows,
        rss_backend=None if observer is None else observer.backend,
        observer_failure=None if observer is None else observer.error,
        observer_exhausted=False if observer is None else observer.exceeded,
        worker_configuration=None if workers is None else asdict(workers),
        transition="new QS pool"
        if workers and spec["endpoint"] == "portfolio"
        else "QS pool"
        if workers
        else "serial",
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--request", type=Path, required=True)
    parser.add_argument("--result", type=Path, required=True)
    args = parser.parse_args()
    if args.result.exists():
        parser.error("preserve the existing result")
    request = json.loads(args.request.read_text())
    resource.setrlimit(
        resource.RLIMIT_FSIZE, (request["cell"]["storage"],) * 2
    )
    result = execute(request)
    with args.result.open("x") as stream:
        json.dump(result, stream, sort_keys=True)
        stream.write("\n")


if __name__ == "__main__":
    main()
