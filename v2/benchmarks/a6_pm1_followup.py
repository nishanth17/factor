"""Frozen A6 execution follow-up with complete API and portfolio costs."""

import argparse
import builtins
import contextlib
import hashlib
import importlib.util
import json
import math
import platform
import random
import resource
import statistics
import subprocess
import sys
import time
from pathlib import Path
from unittest.mock import patch

from v2 import pm1_bounded, portfolio, stage_jobs, utils
from v2.benchmarks.a6_pm1 import (
    assert_quiet,
    performance_window,
    require_runtime,
    stage_fixtures,
    verify_inputs,
)
from v2.benchmarks.build_phase_two_corpus import verify_certificates
from v2.budget import Budget
from v2.pm1_tuning import PM1TuningConfig

ROOT = Path(__file__).parents[1]
INPUTS = Path(__file__).parent / "inputs/corpora"
PROTOCOL = INPUTS / "a6_pm1_followup_protocol.json"
SELECTION = INPUTS / "a6_pm1_followup_selection.json"
FRESH = INPUTS / "a6_followup_corpus.json"
BASELINE = (
    Path(__file__).parent / "inputs/baselines/a6_followup_pm1_32b3c65.py"
)


def frozen_control():
    name = "v2._a6_followup_control"
    if name not in sys.modules:
        spec = importlib.util.spec_from_file_location(name, BASELINE)
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
    return sys.modules[name]


def identity():
    sources = (
        Path(__file__),
        PROTOCOL,
        BASELINE,
        ROOT / "pm1_bounded.py",
        ROOT / "pm1_tuning.py",
        FRESH,
    )
    return {
        str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sources
    }


def verify_fresh():
    data = json.loads(FRESH.read_text())
    verify_certificates(data["certificates"])
    for fixture in data["fixtures"]:
        if math.prod(fixture["factors"]) != fixture["n"]:
            raise AssertionError("fresh reconstruction failed")
    return data["fixtures"]


def selected():
    manifest = json.loads(SELECTION.read_text())
    if manifest["source_identity"] != identity():
        raise ValueError("selected source/input identity changed")
    return manifest


def options_for(arm, selection=None):
    if arm == "legacy_hook64":
        return {"chunk_size": 64, "legacy_hook": True}
    if arm.startswith("control"):
        return {"chunk_size": int(arm.removeprefix("control"))}
    if arm.startswith("bits"):
        return {"chunk_size": 256, "chunk_bits": int(arm.removeprefix("bits"))}
    if arm == "recurrence64":
        return {"chunk_size": 64, "gap_mode": "recurrence"}
    selection = selected() if selection is None else selection
    arithmetic = options_for(selection["arithmetic"])
    if arm in ("arithmetic_control", "selected_arithmetic"):
        return arithmetic
    if arm in ("wheel30", "wheel210"):
        return {**arithmetic, "wheel": int(arm.removeprefix("wheel"))}
    if arm == "selected_pairing":
        return (
            {**arithmetic, "wheel": selection["wheel"]}
            if selection["wheel"]
            else arithmetic
        )
    raise ValueError("unknown frozen arm")


def configuration(options, bounds):
    options = dict(options)
    legacy_hook = options.pop("legacy_hook", False)
    common = dict(
        bounds=tuple(tuple(p) for p in bounds), gcd_batch=64, segment_size=256
    )
    if set(options) == {"chunk_size"}:
        module = pm1_bounded if legacy_hook else frozen_control()
        return module, module.PM1Config(**common, **options)
    return pm1_bounded, PM1TuningConfig(**common, **options)


def grant(work=2_000_000, seconds=30):
    return Budget(work_limit=work, seconds=seconds, cpu_seconds=seconds)


def run(n, bounds, options, *, checkpoint=None, max_actions=None):
    module, cfg = configuration(options, bounds)
    outcome = module.factorize_pm1_bounded(
        n,
        config=cfg,
        budget=grant(),
        checkpoint=checkpoint,
        max_actions=max_actions,
    )
    if outcome.result.reconstruct() != n:
        raise AssertionError("complete API reconstruction failed")
    if outcome.divisor is not None and not utils.valid_divisor(
        outcome.divisor, n
    ):
        raise AssertionError("invalid split")
    if outcome.reason not in (
        "factor_found",
        "exhausted",
        "saturated",
        "paused",
    ):
        raise AssertionError("unexpected finite refusal: " + outcome.reason)
    return outcome


def summary(outcome):
    return {
        "divisor": outcome.divisor,
        "reason": outcome.reason,
        "work": outcome.work_used,
        "verification_work": outcome.verification_work,
        "checkpoint_bytes": len(json.dumps(outcome.checkpoint).encode()),
        "steps": outcome.checkpoint["payload"]["state"]["steps"],
    }


def measure(functions, samples):
    generator = random.Random(2026100917)
    loops, warmups = {}, {}
    values = {name: [] for name in functions}
    for name, function in functions.items():
        start, count = time.monotonic(), 0
        while time.monotonic() - start < 3:
            function()
            count += 1
        warmups[name] = {
            "iterations": count,
            "seconds": time.monotonic() - start,
        }
        loops[name] = 1
        while True:
            start = time.process_time()
            for _ in range(loops[name]):
                function()
            if time.process_time() - start >= 0.08:
                break
            loops[name] *= 2
    target = samples
    for index in range(63):
        order = list(functions)
        generator.shuffle(order)
        for name in order:
            assert_quiet()
            start = time.process_time()
            for _ in range(loops[name]):
                functions[name]()
            values[name].append((time.process_time() - start) / loops[name])
            assert_quiet()
        if index + 1 < target:
            continue
        unstable = any(relative_iqr(v) > 0.15 for v in values.values())
        if unstable and target < 63:
            target = 27 if target < 27 else 63
        else:
            break
    return {
        name: {
            "cpu_samples": values[name],
            "median_cpu": statistics.median(values[name]),
            "relative_iqr": relative_iqr(values[name]),
            "loops": loops[name],
            "warmup": warmups[name],
        }
        for name in functions
    }


def relative_iqr(values):
    q = statistics.quantiles(values, n=4)
    return (q[2] - q[0]) / statistics.median(values)


def stages(phase, samples):
    protocol = json.loads(PROTOCOL.read_text())
    selection = None if phase == "arithmetic" else selected()
    if phase == "arithmetic":
        arms = protocol["chunk_screen"]["arms"]
    elif phase == "paired":
        arms = protocol["paired_screen"]["arms"]
    else:
        arms = protocol["confirmation"]["arms"]
    fixtures = (
        verify_fresh() if phase == "confirmation" else stage_fixtures("screen")
    )
    bounds = (
        protocol["confirmation_bounds"]
        if phase == "confirmation"
        else protocol["screen_bounds"]
    )
    records = []
    for fixture in fixtures:
        for b1, b2 in bounds:
            functions = {
                arm: lambda arm=arm: run(
                    fixture["n"], ((b1, b2),), options_for(arm, selection)
                )
                for arm in arms
            }
            reference = functions[arms[0]]()
            for arm, function in functions.items():
                result = function()
                if bool(result.divisor) != bool(reference.divisor):
                    raise AssertionError("matched stage split outcome changed")
            measurements = measure(functions, samples)
            records.append(
                {
                    "id": fixture["id"],
                    "digits": fixture["digits"],
                    "bounds": [b1, b2],
                    "arms": {
                        arm: {
                            "options": options_for(arm, selection),
                            "outcome": summary(function()),
                            **measurements[arm],
                        }
                        for arm, function in functions.items()
                    },
                }
            )
            print(phase, fixture["id"], b1, b2, flush=True)
    return records


def paired_interval(cells, arm, reference):
    generator = random.Random(2026100918)
    ratios = [
        cell["arms"][arm]["median_cpu"] / cell["arms"][reference]["median_cpu"]
        for cell in cells
    ]
    gain = 1 - math.exp(statistics.mean(math.log(r) for r in ratios))
    resamples = []
    for _ in range(10000):
        logs = []
        for cell in cells:
            candidate = cell["arms"][arm]["cpu_samples"]
            control = cell["arms"][reference]["cpu_samples"]
            indices = [
                generator.randrange(min(len(candidate), len(control)))
                for _ in range(min(len(candidate), len(control)))
            ]
            ratio = statistics.median(
                candidate[i] for i in indices
            ) / statistics.median(control[i] for i in indices)
            logs.append(math.log(ratio))
        resamples.append(1 - math.exp(statistics.mean(logs)))
    resamples.sort()
    return {
        "gain": gain,
        "ci95": [resamples[250], resamples[9749]],
        "all_positive": all(r < 1 for r in ratios),
    }


def select(data, phase):
    if data["source_identity"] != identity():
        raise ValueError("capture identity changed")
    cells = data["records"]
    reference = "control64" if phase == "arithmetic" else "arithmetic_control"
    arms = (
        ("bits256", "bits512", "bits1024", "recurrence64")
        if phase == "arithmetic"
        else ("wheel30", "wheel210")
    )
    decisions = {arm: paired_interval(cells, arm, reference) for arm in arms}
    eligible = [
        arm
        for arm in arms
        if decisions[arm]["all_positive"] and decisions[arm]["ci95"][0] > 0
    ]
    scope = "all_screen_bounds" if eligible else "retain"
    if phase == "paired" and not eligible:
        larger = [c for c in cells if c["bounds"] == [11000, 100000]]
        for arm in arms:
            decisions[arm]["larger_bounds"] = paired_interval(
                larger, arm, reference
            )
        eligible = [
            arm
            for arm in arms
            if decisions[arm]["larger_bounds"]["all_positive"]
            and decisions[arm]["larger_bounds"]["ci95"][0] > 0
        ]
        if eligible:
            scope = "larger_bounds"
    winner = (
        max(eligible, key=lambda arm: decisions[arm]["gain"])
        if eligible
        else reference
    )
    manifest = (
        {"source_identity": identity(), "arithmetic": winner, "wheel": 0}
        if phase == "arithmetic"
        else selected()
    )
    if phase == "paired":
        manifest["wheel_scope"] = scope
        manifest["wheel"] = (
            int(winner.removeprefix("wheel"))
            if winner.startswith("wheel")
            else 0
        )
    manifest[phase + "_capture_sha256"] = hashlib.sha256(
        json.dumps(data, sort_keys=True).encode()
    ).hexdigest()
    SELECTION.write_text(json.dumps(manifest, indent=2) + "\n")
    return {"winner": winner, "comparisons": decisions}


def bridge(options):
    original = stage_jobs.advance_job
    _, cfg = configuration(options, ((2000, 200000),))
    module = (
        frozen_control()
        if type(cfg) is frozen_control().PM1Config
        else pm1_bounded
    )

    def advance(job, budget, context, config):
        if job["kind"] != "pm1":
            original(job, budget, context, config)
            return
        if "a6_state" not in job:
            base = 2 + job["seed"] % max(1, job["n"] - 3)
            job["a6_state"] = module._initial(job["n"], base, cfg)
        module._advance(job["a6_state"], budget, context, cfg)
        if job["a6_state"]["done"]:
            job.update(done=True, factor=job["a6_state"]["factor"])

    return advance


def portfolio_run(arm):
    fixtures = json.loads(
        (INPUTS / "phase_two_m15_independent_corpus.json").read_text()
    )["fixtures"]
    fixtures = [f for f in fixtures if len(str(f["n"])) == 20][:12]
    config = portfolio.PortfolioConfig(
        trial_bound=100,
        chunk_size=16,
        gcd_batch=64,
        rho_attempts=1,
        rho_evaluations=512,
        pm1_attempts=0 if arm == "no_pm1" else 1,
        pm1_b1=2000,
        pm1_b2=200000,
        ecm_tiers=((50, 1000, 2),),
        segment_size=256,
    )
    completed = splits = 0
    for fixture in fixtures:
        outcome = portfolio.factorize_bounded(
            fixture["n"], seed=7, config=config, budget=grant(500000, 5)
        )
        if outcome.result.reconstruct() != fixture["n"]:
            raise AssertionError("portfolio reconstruction failed")
        completed += outcome.result.complete
        splits += len(outcome.result.remaining) != 1 or bool(
            outcome.result.factors
        )
    return {"complete": completed, "splits": splits, "inputs": len(fixtures)}


def portfolio_call(arm, chosen_options=None):
    if arm in ("no_pm1", "retained_portfolio"):
        return portfolio_run(arm)
    options = (
        {"chunk_size": 16}
        if arm == "bounded_control_bridge"
        else chosen_options
    )
    with patch.object(portfolio, "advance_job", bridge(options)):
        return portfolio_run(arm)


def portfolio_measure(samples):
    arms = json.loads(PROTOCOL.read_text())["portfolio"]["arms"]
    chosen_options = options_for("selected_pairing")
    functions = {
        arm: lambda arm=arm: portfolio_call(arm, chosen_options)
        for arm in arms
    }
    measurements = measure(functions, samples)
    return {
        arm: {"outcome": function(), **measurements[arm]}
        for arm, function in functions.items()
    }


def continuation_call(n, options, arm, pause_steps=None):
    bounds = ((500, 5000), (2000, 20000))
    if arm == "fresh_final":
        return summary(run(n, (bounds[-1],), options))
    module, cfg = configuration(options, bounds)
    if arm == "fresh_each_shared_allowance":
        prior_work = prior_wall = prior_cpu = 0
        for pair in bounds:
            _, single = configuration(options, (pair,))
            outcome = module.factorize_pm1_bounded(
                n,
                config=single,
                budget=Budget(
                    work_limit=2000000 - prior_work,
                    seconds=30 - prior_wall,
                    cpu_seconds=30 - prior_cpu,
                ),
            )
            prior_work += outcome.work_used
            prior_wall += outcome.wall_seconds
            prior_cpu += outcome.cpu_seconds
        result = summary(outcome)
        result["work"] = prior_work
        return result
    if arm == "in_memory":
        return summary(run(n, bounds, options))
    if arm != "verified_checkpoint":
        raise ValueError("unknown continuation arm")
    paused = run(n, bounds, options, max_actions=pause_steps)
    resumed = run(
        n,
        bounds,
        options,
        checkpoint=json.loads(json.dumps(paused.checkpoint)),
    )
    return summary(resumed)


def continuation_measure(samples):
    records = []
    arms = json.loads(PROTOCOL.read_text())["continuation"]["arms"]
    for fixture in verify_fresh():
        for name, options in (
            ("control64", {"chunk_size": 64}),
            ("selected", options_for("selected_pairing")),
        ):
            module, cfg = configuration(options, ((500, 5000), (2000, 20000)))
            state = module._initial(fixture["n"], 2, cfg)
            from v2.schedules import SieveContext

            context, ledger = SieveContext(20001, segment_size=256), grant()
            while state["rung"] == 0 and not state["done"]:
                module._advance(state, ledger, context, cfg)
            pause_steps = state["steps"]
            functions = {
                arm: lambda arm=arm: continuation_call(
                    fixture["n"], options, arm, pause_steps
                )
                for arm in arms
            }
            measurements = measure(functions, samples)
            records.append(
                {
                    "id": fixture["id"],
                    "configuration": name,
                    "options": options,
                    "arms": {
                        arm: {"outcome": function(), **measurements[arm]}
                        for arm, function in functions.items()
                    },
                }
            )
    return records


def diagnostic_profile():
    """Instrument complete calls separately from accepted timing samples."""
    from v2 import pm1_tuning

    fixture = next(f for f in stage_fixtures("screen") if f["digits"] == 50)
    records = []
    for arm in (
        "control16",
        "control64",
        "bits256",
        "bits512",
        "bits1024",
        "recurrence64",
        "wheel30",
        "wheel210",
    ):
        options = options_for(arm, {"arithmetic": "control64"})
        module, cfg = configuration(options, ((2000, 20000),))
        phases, counts = (
            {},
            {
                "pow": 0,
                "exponent_bits": 0,
                "inversions": 0,
                "gcd": 0,
                "gcd_argument_bits": 0,
            },
        )
        original = module._advance

        def advance(state, budget, context, config):
            phase = state["phase"]
            start = time.process_time()
            original(state, budget, context, config)
            record = phases.setdefault(phase, {"calls": 0, "cpu": 0})
            record["calls"] += 1
            record["cpu"] += time.process_time() - start

        def powering(base, exponent, modulus):
            counts["pow"] += 1
            counts["inversions"] += exponent < 0
            counts["exponent_bits"] += max(0, exponent).bit_length()
            return builtins.pow(base, exponent, modulus)

        def gcd(left, right):
            counts["gcd"] += 1
            counts["gcd_argument_bits"] += max(
                left.bit_length(), right.bit_length()
            )
            return math.gcd(left, right)

        with contextlib.ExitStack() as patches:
            patches.enter_context(patch.object(module, "_advance", advance))
            for target in {module, pm1_bounded, pm1_tuning}:
                patches.enter_context(
                    patch.object(target, "pow", powering, create=True)
                )
                patches.enter_context(
                    patch.object(target, "gcd", gcd, create=True)
                )
            start = time.process_time()
            outcome = run(fixture["n"], ((2000, 20000),), options)
            total = time.process_time() - start
        records.append(
            {
                "arm": arm,
                "options": options,
                "total_instrumented_cpu": total,
                "phases": phases,
                "operation_counts": counts,
                "outcome": summary(outcome),
                "owned_workspace_reserve": cfg.workspace_reserve,
                "configured_memory_cap": cfg.memory_bytes,
            }
        )
    return {
        "instrumented_not_performance_evidence": True,
        "records": records,
        "process_peak_rss": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        "rss_units": "bytes on macOS; KiB on Linux",
    }


def cold_calls():
    records = []
    for arm in ("control64", "selected_pairing"):
        values = []
        for _ in range(9):
            assert_quiet()
            start = time.monotonic()
            output = subprocess.check_output(
                [
                    sys.executable,
                    "-m",
                    "v2.benchmarks.a6_pm1_followup",
                    "cold-one",
                    "--arm",
                    arm,
                    "--output",
                    "/private/tmp/a6-cold-unused.json",
                ],
                text=True,
            )
            values.append(time.monotonic() - start)
            if json.loads(output)["reason"] not in (
                "factor_found",
                "exhausted",
                "saturated",
            ):
                raise AssertionError("cold validation failed")
        records.append(
            {
                "arm": arm,
                "wall_samples_including_startup": values,
                "median_wall": statistics.median(values),
            }
        )
    return records


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "phase",
        choices=(
            "arithmetic",
            "paired",
            "confirmation",
            "portfolio",
            "continuation",
            "select-arithmetic",
            "select-paired",
            "profile",
            "cold",
            "cold-one",
        ),
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--input", type=Path)
    parser.add_argument("--samples", type=int, choices=(9, 27, 63), default=9)
    parser.add_argument("--arm", choices=("control64", "selected_pairing"))
    args = parser.parse_args()
    require_runtime()
    if args.phase == "cold-one":
        selected()  # Match source/input verification in both cold arms.
        fixture = next(
            f for f in stage_fixtures("screen") if f["digits"] == 20
        )
        print(
            json.dumps(
                summary(
                    run(fixture["n"], ((2000, 200000),), options_for(args.arm))
                )
            )
        )
        return
    if args.phase.startswith("select-"):
        report = select(
            json.loads(args.input.read_text()),
            args.phase.removeprefix("select-"),
        )
        args.output.write_text(json.dumps(report, indent=2) + "\n")
        print(json.dumps(report, indent=2))
        return
    with performance_window():
        verify_inputs()
        verify_fresh()
        if args.phase in ("arithmetic", "paired", "confirmation"):
            records = stages(args.phase, args.samples)
        elif args.phase == "profile":
            records = diagnostic_profile()
        elif args.phase == "cold":
            records = cold_calls()
        elif args.phase == "portfolio":
            records = portfolio_measure(args.samples)
        else:
            records = continuation_measure(args.samples)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(
            {
                "phase": args.phase,
                "runtime": sys.version,
                "platform": platform.platform(),
                "source_identity": identity(),
                "records": records,
            },
            indent=2,
        )
        + "\n"
    )
    print(args.output)


if __name__ == "__main__":
    main()
