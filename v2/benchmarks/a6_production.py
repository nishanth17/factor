"""Frozen direct-production A6 default promotion comparison."""

import argparse
import hashlib
import importlib.util
import json
import math
import random
import statistics
import sys
from pathlib import Path

from v2 import portfolio, stage_jobs
from v2.benchmarks.a6_pm1 import performance_window, require_runtime
from v2.benchmarks.a6_pm1_followup import measure, verify_fresh
from v2.benchmarks.build_phase_two_corpus import verify_certificates
from v2.budget import Budget
from v2.schedules import SieveContext

ROOT = Path(__file__).parents[1]
CONTROLS = Path(__file__).parent / "inputs/controls/a6_production_319d5c6"
PROTOCOL = CONTROLS.parent / "a6_production_protocol.json"
MANIFEST = CONTROLS.parent / "a6_production_identity.json"
ARMS = ("control", "recurrence16", "chunk64", "recurrence64")


def load_control(name):
    qualified = "v2._a6_production_" + name
    if qualified not in sys.modules:
        from importlib.machinery import SourceFileLoader

        path = CONTROLS / (name + ".py.txt")
        loader = SourceFileLoader(qualified, str(path))
        spec = importlib.util.spec_from_loader(qualified, loader)
        module = importlib.util.module_from_spec(spec)
        sys.modules[qualified] = module
        loader.exec_module(module)
        if name == "portfolio":
            module.advance_job = load_control("stage_jobs").advance_job
    return sys.modules[qualified]


def identity():
    files = [
        *sorted(ROOT.glob("*.py")),
        Path(__file__),
        PROTOCOL,
        *sorted(CONTROLS.glob("*.txt")),
        ROOT / "benchmarks/a6_pm1.py",
        ROOT / "benchmarks/a6_pm1_followup.py",
        ROOT / "benchmarks/build_phase_two_corpus.py",
        ROOT / "benchmarks/inputs/corpora/a6_followup_corpus.json",
        ROOT / "benchmarks/inputs/corpora/p43_size_corpus.json",
        ROOT
        / "benchmarks/inputs/corpora/phase_two_m15_independent_corpus.json",
    ]
    return {
        str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in files
    }


def options(arm):
    return dict(
        pm1_gap_mode="recurrence" if "recurrence" in arm else "cached",
        pm1_chunk_size=64 if "64" in arm else None,
    )


def config(arm, backend, **changes):
    module = load_control("portfolio") if arm == "control" else portfolio
    common = dict(
        trial_bound=100,
        rho_attempts=1,
        rho_evaluations=512,
        chunk_size=16,
        gcd_batch=64,
        segment_size=256,
        pm1_b1=2000,
        pm1_b2=200000,
        ecm_tiers=((50, 1000, 2),),
        backend=backend,
    )
    if arm != "control":
        common.update(options(arm))
    common.update(changes)
    return module, module.PortfolioConfig(**common)


def grant():
    return Budget(work_limit=2_000_000, seconds=30, cpu_seconds=30)


def stage(arm, fixture, bounds, backend):
    b1, b2 = bounds
    _, cfg = config(arm, backend, pm1_b1=b1, pm1_b2=b2, ecm_tiers=())
    from v2.arithmetic import get_backend

    n = get_backend(backend).integer(fixture["n"])
    jobs = load_control("stage_jobs") if arm == "control" else stage_jobs
    context = SieveContext(b2 + 1, segment_size=256)
    ledger = grant()
    ledger.consume((math.isqrt(b2) + 1) // 2)
    job = jobs.new_job("pm1", n, 0, b1, b2)
    while not job["done"]:
        jobs.advance_job(job, ledger, context, cfg)
    if job["factor"] is not None:
        raise AssertionError("complete-stage fixture split unexpectedly")
    return {
        "work": ledger.used,
        "factor": None,
        "workspace": cfg.workspace_reserve,
        "table_entries": len(job.get("even_powers", [])),
        "job_bytes": len(json.dumps(job, default=int)),
    }


def fixtures(kind):
    if kind == "confirmation":
        return verify_fresh()
    path = ROOT / "benchmarks/inputs/corpora/p43_size_corpus.json"
    data = json.loads(path.read_text())
    verify_certificates(data["certificates"])
    return [
        next(
            f
            for f in data["fixtures"]
            if f["digits"] == digits and f["split"] == "screen"
        )
        for digits in (20, 50, 100)
    ]


def portfolio_call(arm, backend):
    path = (
        ROOT
        / "benchmarks/inputs/corpora/phase_two_m15_independent_corpus.json"
    )
    data = json.loads(path.read_text())
    selected = [f for f in data["fixtures"] if len(str(f["n"])) == 20][:12]
    module, cfg = config(arm, backend)
    outcomes = []
    for fixture in selected:
        outcome = module.factorize_bounded(
            fixture["n"],
            seed=7,
            config=cfg,
            budget=Budget(work_limit=500000, seconds=5, cpu_seconds=5),
        )
        if outcome.result.reconstruct() != fixture["n"]:
            raise AssertionError("portfolio reconstruction failed")
        outcomes.append(
            {
                "complete": outcome.result.complete,
                "factors": [
                    (int(f.value), f.multiplicity, f.certainty.value)
                    for f in outcome.result.factors
                ],
                "remaining": list(map(int, outcome.result.remaining)),
                "work": outcome.work_used,
            }
        )
    return outcomes


def confidence(control, candidate):
    generator = random.Random(2026100927)
    gains = []
    for _ in range(10000):
        a = statistics.median(generator.choices(control, k=len(control)))
        b = statistics.median(generator.choices(candidate, k=len(candidate)))
        gains.append(1 - b / a)
    gains.sort()
    return {
        "gain": 1 - statistics.median(candidate) / statistics.median(control),
        "ci95": [gains[250], gains[9749]],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "phase", choices=("freeze", "screen", "confirmation", "portfolio")
    )
    parser.add_argument("--output", type=Path)
    parser.add_argument("--backend", default="python-int")
    parser.add_argument("--arm", choices=ARMS, default="recurrence64")
    parser.add_argument("--samples", type=int, choices=(9, 27, 63), default=27)
    args = parser.parse_args()
    require_runtime()
    if args.phase == "freeze":
        MANIFEST.write_text(json.dumps(identity(), indent=2) + "\n")
        return
    frozen = json.loads(MANIFEST.read_text())
    if identity() != frozen:
        raise ValueError("production source/input identity changed")
    records = []
    with performance_window():
        if args.phase == "portfolio":
            arms = ("control", args.arm)
            functions = {
                arm: lambda arm=arm: portfolio_call(arm, args.backend)
                for arm in arms
            }
            reference = functions["control"]()
            for function in functions.values():
                outcome = function()
                if [dict(v, work=0) for v in outcome] != [
                    dict(v, work=0) for v in reference
                ]:
                    raise AssertionError("portfolio outcomes changed")
            measurements = measure(functions, args.samples)
            records.append(
                {
                    "arms": {
                        arm: {**measurements[arm], "outcome": function()}
                        for arm, function in functions.items()
                    }
                }
            )
        else:
            arms = ARMS if args.phase == "screen" else ("control", args.arm)
            for fixture in fixtures(args.phase):
                for bounds in ((2000, 20000), (11000, 100000), (2000, 200000)):
                    functions = {
                        arm: lambda arm=arm: stage(
                            arm, fixture, bounds, args.backend
                        )
                        for arm in arms
                    }
                    measured = measure(functions, args.samples)
                    records.append(
                        {
                            "id": fixture["id"],
                            "digits": fixture["digits"],
                            "bounds": bounds,
                            "arms": {
                                arm: {**measured[arm], "outcome": function()}
                                for arm, function in functions.items()
                            },
                        }
                    )
    if identity() != frozen:
        raise ValueError("source changed during measurement")
    for record in records:
        control = record["arms"]["control"]["cpu_samples"]
        record["comparisons"] = {
            arm: confidence(control, data["cpu_samples"])
            for arm, data in record["arms"].items()
            if arm != "control"
        }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(
            {
                "phase": args.phase,
                "backend": args.backend,
                "identity": frozen,
                "records": records,
            },
            indent=2,
        )
        + "\n"
    )
    print(
        json.dumps(
            [{k: v for k, v in r.items() if k != "arms"} for r in records],
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
