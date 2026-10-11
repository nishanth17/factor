"""Finite attribution and feasibility before the next C3 selection freeze."""

import argparse
import contextlib
import json
import os
import resource
import sys
import time
from dataclasses import asdict
from importlib import import_module
from unittest.mock import patch

from .... import portfolio
from ....execution.allocation import ECMAllocation
from ...qs.c1.c1_implementation import decode_config
from ..p52.p52_realistic import check_quiet
from . import c3_study as first
from .round2_corpus import load_training, write_training

CONTROL = first.INPUTS / "controls/c3_round2_pilot_v2.json"
# These are mechanism bundles, not a Cartesian sweep. Economic ceilings are
# cumulative fractions of historical C1-sized relation-call estimates.
POLICIES = {
    "control": dict(tiers=((2000, 147396, 32),), mode=None),
    "protected32": dict(tiers=((2000, 147396, 32),), mode="campaign"),
    "compact32": dict(tiers=((2000, 50000, 32),), mode="campaign"),
    "expanded64": dict(tiers=((2000, 147396, 64),), mode="campaign"),
    "deep4": dict(
        tiers=((2000, 147396, 16), (11000, 250000, 4)), mode="campaign"
    ),
    "economic25": dict(
        tiers=((2000, 147396, 64),), mode="pretest", fraction=0.25
    ),
    "economic50": dict(
        tiers=((2000, 147396, 64),), mode="pretest", fraction=0.5
    ),
}


def configuration(n, arm, engine):
    """Use only input size, public policy and frozen relation settings."""
    band = first.band_of(n)
    selected = json.loads(
        (first.INPUTS / "controls/c1_selected.json").read_text()
    )
    values = selected["configurations"][str(band)][
        "slp" if band == 30 else "dlp_half"
    ]
    relation = decode_config(values, import_module(engine.__package__ + ".qs"))
    settings = POLICIES[arm]
    options = dict(
        ecm_tiers=settings["tiers"],
        siqs=relation,
        memory_bytes=288 * 2**20,
        trace_limit=512,
    )
    mode = settings["mode"]
    if mode:
        limit = settings.get("fraction", 0) * (0.2 if band == 30 else 4.5)
        options["allocation"] = ECMAllocation(
            mode,
            pretest_work=10**13 if mode == "pretest" else None,
            pretest_seconds=limit if mode == "pretest" else None,
            pretest_cpu_seconds=limit if mode == "pretest" else None,
            fallback_work=500_000_000 if band == 30 else 15_000_000_000,
            fallback_seconds=1 if band == 30 else 10,
            fallback_cpu_seconds=1 if band == 30 else 10,
        )
    return engine.PortfolioConfig(**options)


@contextlib.contextmanager
def attribution(engine, enabled):
    """Instrument curve cost separately; never call this acceptance timing."""
    curves, active = [], {}
    if not enabled:
        yield curves
        return
    original = engine.advance_job

    def advance(job, budget, context, config):
        began, cpu = time.perf_counter(), time.process_time()
        try:
            return original(job, budget, context, config)
        finally:
            key = id(job)
            costs = active.setdefault(key, [0.0, 0.0])
            costs[0] += time.perf_counter() - began
            costs[1] += time.process_time() - cpu
            if job["done"]:
                curves.append(
                    dict(
                        stage=job["kind"],
                        n=int(job["n"]),
                        seed=job["seed"],
                        b1=job["b1"],
                        b2=job["b2"],
                        outcome="factor" if job["factor"] else "exhausted",
                        wall=costs[0],
                        cpu=costs[1],
                        work=budget.used - job["start_work"],
                    )
                )
                active.pop(key)

    with patch.object(engine, "advance_job", advance):
        yield curves


def run_one(
    fixture, seed, arm, control, *, instrument=False, regime="service"
):
    engine = control if arm == "control" else portfolio
    config = configuration(fixture["n"], arm, engine)
    band = first.band_of(fixture["n"])
    seconds = (5 if band == 30 else 30) if regime == "service" else 2
    work = 10**13 if regime == "service" else 2_000_000
    ledger = import_module(engine.__package__ + ".execution.budget")
    began, cpu = time.perf_counter(), time.process_time()
    with attribution(engine, instrument) as curve_costs:
        run = engine.factorize_bounded(
            fixture["n"],
            seed=seed,
            config=config,
            budget=ledger.Budget(
                work_limit=work, seconds=seconds, cpu_seconds=seconds
            ),
        )
    first.validate(run.result, fixture)
    elapsed, cpu_elapsed = (
        time.perf_counter() - began,
        time.process_time() - cpu,
    )
    current = run.checkpoint["payload"]["state"]["current"] or {}
    return dict(
        id=fixture["id"],
        kind=fixture["kind"],
        band=band,
        seed=seed,
        arm=arm,
        instrumented=instrument,
        complete=run.result.complete,
        proper_factor=any(
            first.utils.valid_divisor(value, abs(fixture["n"]))
            for value in (
                *[p.value for p in run.result.factors],
                *run.result.remaining,
            )
        ),
        reason=run.reason,
        wall=elapsed,
        cpu=cpu_elapsed,
        work=run.work_used,
        cap_seconds=seconds,
        events=run.events,
        curve_costs=curve_costs,
        stage_seconds=run.checkpoint["payload"]["state"]["stage_seconds"],
        bounds=[list(t) for t in config.ecm_tiers],
        allocation=(
            asdict(config.allocation)
            if getattr(config, "allocation", None)
            else None
        ),
        active_curve=current.get("job"),
        fallback_started=any(e["stage"] == "siqs" for e in run.events)
        or "siqs_seed" in current,
        memory_cap=config.memory_bytes,
        workspace_reserve=config.workspace_reserve,
        fallback_owned_peak=max(
            (e.get("stats", {}).get("workspace_bytes", 0) for e in run.events),
            default=0,
        ),
        fallback_recoveries=sum(
            e.get("stats", {}).get("recoveries", 0) for e in run.events
        ),
        rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        checkpoint_bytes=len(json.dumps(run.checkpoint).encode()),
        curves=sum(
            e["stage"] == "ecm" and e["outcome"] != "handoff"
            for e in run.events
        ),
    )


def verify_control():
    protocol = json.loads(CONTROL.read_text())
    if protocol["policies"] != json.loads(json.dumps(POLICIES)):
        raise ValueError("round-two pilot policies changed")
    for name, expected in protocol["sha256"].items():
        if first.digest(first.ROOT / name) != expected:
            raise ValueError("round-two pilot source changed: " + name)
    return protocol


def measure(output, *, own_window=True):
    """At most 480 active seconds, leaving 120 seconds for focused QA."""
    first.runtime()
    protocol = verify_control()
    fixtures = load_training()
    generated = [f for f in fixtures if f["id"].startswith("r2_")]
    # The first independent input in each selected generation stratum is fixed
    # before any timing. Expected factors stay in the validator.
    cases = [
        f
        for f in generated
        if f["id"].endswith("_0")
        and f["kind"] in ("balanced", "uneven_12", "uneven_14", "uneven_16")
    ]
    rows = []
    started, cpu = time.monotonic(), time.process_time()
    stopped = None
    window = first.machine_window() if own_window else contextlib.nullcontext()
    with window:
        control = first.baseline()
        for fixture in cases:
            for seed_index, seed in enumerate(protocol["seeds"]):
                arms = list(POLICIES)
                offset = seed_index % len(arms)
                arms = arms[offset:] + arms[:offset]
                if seed_index % 2:
                    arms.reverse()
                for arm in arms:
                    # Reserve the full next call; a partial pilot must not
                    # silently overrun the lease while a finite call is active.
                    cap = 5 if first.band_of(fixture["n"]) == 30 else 30
                    if (
                        max(
                            time.monotonic() - started,
                            time.process_time() - cpu,
                        )
                        + cap
                        > protocol["active_seconds"]
                    ):
                        stopped = "pilot_allowance"
                        break
                    check_quiet({os.getpid()})
                    row = run_one(fixture, seed, arm, control, instrument=True)
                    first.save(
                        output.parent
                        / (output.stem + "-rows")
                        / f"{fixture['id']}-{seed}-{arm}.json",
                        row,
                    )
                    rows.append(row)
                    print(
                        fixture["id"],
                        seed,
                        arm,
                        row["reason"],
                        round(row["wall"], 4),
                        flush=True,
                    )
                if stopped:
                    break
            if stopped:
                break
    first.save(
        output,
        dict(
            mode="round2-instrumented-pilot",
            rows=rows,
            stopped=stopped,
            wall=time.monotonic() - started,
            cpu=time.process_time() - cpu,
            protocol_sha256=first.digest(CONTROL),
            runtime=sys.version,
            timing_status="Instrumented feasibility; no performance claim",
        ),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("generate-training", "pilot"))
    parser.add_argument("--output", type=first.Path)
    args = parser.parse_args()
    first.runtime()
    if args.mode == "generate-training":
        with first.machine_window():
            print(write_training())
    else:
        if args.output is None:
            parser.error("pilot requires --output")
        measure(args.output)


if __name__ == "__main__":
    main()
