"""Matched P5.2 program controls, finite campaign probes and cold startup."""

import argparse
import hashlib
import json
import math
import platform
import random
import resource
import statistics
import subprocess
import sys
import time
import types
from collections import Counter
from pathlib import Path
from struct import pack

from .... import portfolio
from ....ecm.programs import ECMPrograms
from ....execution.stage_jobs import advance_job, new_job
from ...suites.build_phase_two_corpus import verify_certificates
from ...support.paths import (
    BENCHMARK_ROOT,
    REPOSITORY_ROOT,
    source_path,
)

INPUTS = BENCHMARK_ROOT / "inputs"
BASELINE = INPUTS / "baselines/p52_a3_baseline.json"
DEPENDENCIES = INPUTS / "baselines/p52_a3_dependencies.json"
CORPUS = INPUTS / "corpora/p52_a3_corpus.json"
ARMS = ("baseline", "default", "programs", "regenerated")
SCHEDULE_CASES = {
    "default_schedule": (2000, 147396, 3),
    "middle_schedule": (11000, 1900000, 3),
    "large_schedule": (50000, 5000000, 3),
}
CASES = ("small", "medium", "campaign", "large_campaign", *SCHEDULE_CASES)


def load_control():
    """Load a hash-checked private snapshot for ECM-only control arms.

    Inactive QS fallback types support the original dataclass annotations;
    control experiments always disable those fallbacks. Arithmetic, sieve,
    budget and result helpers are frozen, independent of later v2 changes.
    """
    data = json.loads(source_path(BASELINE).read_text())
    dependencies = json.loads(source_path(DEPENDENCIES).read_text())
    if (
        hashlib.sha256(source_path(BASELINE).read_bytes()).hexdigest()
        != dependencies["original_baseline_sha256"]
    ):
        raise ValueError("P5.2 original baseline changed")
    for name, expected in data["dependency_sha256"].items():
        if dependencies["source_sha256"][name] != expected:
            raise ValueError("P5.2 dependency does not match original pin")
    sources = {**dependencies["source"], **data["source"]}
    hashes = {**dependencies["source_sha256"], **data["source_sha256"]}
    for name, source in sources.items():
        if hashlib.sha256(source.encode()).hexdigest() != hashes[name]:
            raise ValueError("corrupt P5.2 frozen source: " + name)

    package_name = "_p52_control"
    package = types.ModuleType(package_name)
    package.__path__ = []
    sys.modules[package_name] = package
    from .... import qs
    from ....qs import sss

    # No fallback is executed by these ECM-only arms. Sharing its type names
    # avoids copying unrelated QS engines into an arithmetic control snapshot.
    sys.modules[package_name + ".qs"] = qs
    sys.modules[package_name + ".qs.sss"] = sss
    package.qs = qs
    modules = {}
    for name in (
        "constants",
        "utils",
        "prime_sieve",
        "budget",
        "preprocessing",
        "ecm",
        "pollard_rho",
        "factor",
        "schedules",
        "stage_jobs",
        "portfolio",
    ):
        path = "v2/" + name + ".py"
        qualified = package_name + "." + name
        module = types.ModuleType(qualified)
        module.__package__ = package_name
        snapshot = BASELINE if path in data["source"] else DEPENDENCIES
        module.__file__ = str(snapshot) + ":" + path
        sys.modules[qualified] = module
        setattr(package, name, module)
        exec(compile(sources[path], module.__file__, "exec"), module.__dict__)
        modules[name] = module
    return modules["portfolio"], modules["stage_jobs"]


def source_hashes(*runners):
    """Pin active helpers and validation runners with each capture."""
    root = REPOSITORY_ROOT
    names = (
        "constants",
        "utils",
        "prime_sieve",
        "budget",
        "preprocessing",
        "ecm",
        "pollard_rho",
        "factor",
        "schedules",
        "stage_jobs",
        "portfolio",
        "ecm_programs",
    )
    paths = ["v2/" + name + ".py" for name in names]
    paths.extend(
        "v2/benchmarks/" + name + ".py"
        for name in ("p52_a3", "build_phase_two_corpus", *runners)
    )
    return {
        name: hashlib.sha256(
            (source_path(root / name)).read_bytes()
        ).hexdigest()
        for name in paths
    }


def load_corpus():
    """Check independent prime proofs and exact input reconstruction."""
    corpus = json.loads(source_path(CORPUS).read_text())
    verify_certificates(corpus["certificates"])
    for fixture in corpus["fixtures"]:
        n = 1
        for prime, exponent in fixture["factors"]:
            n *= prime**exponent
        if n != fixture["n"]:
            raise ValueError("corpus does not reconstruct")
    return corpus


def options(case, arm):
    """Keep seeds, bounds and total caps matched across every arm."""
    b1, b2, curves = {
        "small": (50, 2000, 8),
        "medium": (200, 20000, 16),
        "campaign": (2000, 147396, 4),
        "large_campaign": (11000, 1900000, 3),
        **SCHEDULE_CASES,
    }[case]
    values = dict(
        trial_bound=5,
        rho_attempts=0,
        pm1_attempts=0,
        ecm_tiers=((b1, b2, curves),),
        max_input_bits=329,
        memory_bytes=16 * 2**20,
    )
    if arm == "programs":
        values["ecm_program_bytes"] = 8 * 2**20
    elif arm == "regenerated":
        values["ecm_program_bytes"] = 4096 + 256 * 1024
    return values


def measure(case, arm, corpus, engine, stages):
    """Time validated whole cohorts or explicitly exhausted curve campaigns."""
    if case in SCHEDULE_CASES:
        return measure_schedule(case, arm, corpus, engine)
    started = time.perf_counter()
    cpu_started = time.process_time()
    config = engine.PortfolioConfig(**options(case, arm))
    rows = []
    for fixture in corpus["fixtures"]:
        band = "campaign" if case == "large_campaign" else case
        if fixture["band"] != band:
            continue
        if case == "large_campaign" and fixture["id"] != "campaign_0":
            continue
        seeds = (
            corpus["seeds"][:2]
            if case == "large_campaign"
            else corpus["seeds"]
        )
        for seed in seeds:
            budget = engine.Budget(
                work_limit=50_000_000, seconds=30, cpu_seconds=30
            )
            n = fixture["n"]
            if case in ("campaign", "large_campaign"):
                context = engine.SieveContext(
                    config.max_hi,
                    segment_size=config.segment_size,
                    budget=budget,
                )
                programs = (
                    ECMPrograms(context, memory_bytes=config.ecm_program_bytes)
                    if arm in ("programs", "regenerated")
                    else None
                )
                b1, b2, curves = config.ecm_tiers[0]
                outcomes = []
                for index in range(curves):
                    job = stages.new_job("ecm", n, seed + index, b1, b2)
                    while not job["done"]:
                        if arm in ("programs", "regenerated"):
                            advance_job(job, budget, programs, config)
                        else:
                            stages.advance_job(job, budget, context, config)
                    divisor = job["factor"]
                    if divisor is not None:
                        if not 1 < divisor < n or n % divisor:
                            raise AssertionError("invalid campaign divisor")
                        if divisor * (n // divisor) != n:
                            raise AssertionError("campaign reconstruction")
                    outcomes.append(divisor)
                rows.append(
                    dict(
                        fixture=fixture["id"],
                        seed=seed,
                        factors=outcomes,
                        complete=False,
                        work=budget.used,
                        retained_bytes=programs.used_bytes if programs else 0,
                        hits=programs.hits if programs else 0,
                        misses=programs.misses if programs else 0,
                        unretained=programs.unretained if programs else 0,
                    )
                )
            else:
                run = engine.factorize_bounded(
                    n, seed=seed, config=config, budget=budget
                )
                actual = Counter(
                    {f.value: f.exponent for f in run.result.factors}
                )
                expected = Counter(dict(fixture["factors"]))
                if run.result.reconstruct() != n:
                    raise AssertionError("factorization reconstruction")
                if actual - expected or (
                    run.result.complete and actual != expected
                ):
                    raise AssertionError("unexpected terminal factors")
                rows.append(
                    dict(
                        fixture=fixture["id"],
                        seed=seed,
                        complete=run.result.complete,
                        factors=sorted(actual.items()),
                        remaining=run.result.remaining,
                        reason=run.reason,
                        work=run.work_used,
                    )
                )
    return dict(
        seconds=time.perf_counter() - started,
        cpu_seconds=time.process_time() - cpu_started,
        rows=rows,
    )


def reference_schedule(b1, b2):
    """Full integer-index Eratosthenes oracle, outside retained timings."""
    if b2 > 5_000_000:
        raise ValueError("schedule oracle has a finite five-million bound")
    flags = bytearray(b"\x01") * (b2 + 1)
    flags[:2] = b"\x00\x00"
    for prime in range(2, b2 + 1):
        if flags[prime] and prime * prime <= b2:
            for composite in range(prime * prime, b2 + 1, prime):
                flags[composite] = 0
    digests = [hashlib.sha256(), hashlib.sha256()]
    counts = [0, 0]
    for prime in range(2, b2 + 1):
        if not flags[prime]:
            continue
        index = int(prime > b1)
        counts[index] += 1
        if index:
            digests[index].update(pack("<Q", prime))
        else:
            power = prime
            while power <= b1 // prime:
                power *= prime
            digests[index].update(pack("<QQ", prime, power))
    return dict(
        counts=counts, digests=[digest.hexdigest() for digest in digests]
    )


def measure_schedule(case, arm, corpus, engine):
    """Consume whole integer programs; this is not factorization evidence."""
    started = time.perf_counter()
    cpu_started = time.process_time()
    config = engine.PortfolioConfig(**options(case, arm))
    budget = engine.Budget(work_limit=50_000_000, seconds=30, cpu_seconds=30)
    context = engine.SieveContext(
        config.max_hi, segment_size=config.segment_size, budget=budget
    )
    programs = (
        ECMPrograms(context, memory_bytes=config.ecm_program_bytes)
        if arm in ("programs", "regenerated")
        else None
    )
    b1, b2, curves = config.ecm_tiers[0]
    rows = []
    for curve in range(curves):
        digests = [hashlib.sha256(), hashlib.sha256()]
        counts = [0, 0]
        for index, (lo, hi, bound) in enumerate(
            ((2, b1 + 1, b1), (b1 + 1, b2 + 1, None))
        ):
            for left in range(lo, hi, 2 * context.segment_size):
                right = min(hi, left + 2 * context.segment_size)
                if programs is None:
                    budget.consume(
                        context.segment_size + len(context.base_primes)
                    )
                    primes = context.prime_segment(left, right)
                    powers = (
                        [
                            engine.utils.prime_power(prime, b1)
                            for prime in primes
                        ]
                        if bound
                        else None
                    )
                else:
                    primes = programs.program_segment(
                        left, right, budget, bound=bound
                    )
                    powers = (
                        programs.power_values(left, right, b1, 0, len(primes))
                        if bound
                        else None
                    )
                budget.consume(len(primes))
                counts[index] += len(primes)
                for position, prime in enumerate(primes):
                    data = (
                        pack("<QQ", prime, powers[position])
                        if bound
                        else pack("<Q", prime)
                    )
                    digests[index].update(data)
        actual = dict(
            counts=counts, digests=[digest.hexdigest() for digest in digests]
        )
        if actual != corpus["schedule_oracle"]:
            raise AssertionError(
                "schedule differs from independent full sieve"
            )
        rows.append(
            dict(
                fixture=case,
                seed=0,
                curve=curve,
                complete=False,
                work=budget.used,
                **actual,
                retained_bytes=programs.used_bytes if programs else 0,
                hits=programs.hits if programs else 0,
                misses=programs.misses if programs else 0,
                unretained=programs.unretained if programs else 0,
            )
        )
    return dict(
        seconds=time.perf_counter() - started,
        cpu_seconds=time.process_time() - cpu_started,
        rows=rows,
    )


def campaign_probes():
    """Measure deterministic work/storage feasibility, without speed ratios."""
    corpus = load_corpus()
    fixture = next(
        item for item in corpus["fixtures"] if item["band"] == "campaign"
    )
    n = fixture["n"]
    control, stages = load_control()
    rows = []
    for case in SCHEDULE_CASES:
        for arm in ("baseline", "programs"):
            engine = control if arm == "baseline" else portfolio
            config = engine.PortfolioConfig(**options(case, arm))
            budget = engine.Budget(
                work_limit=50_000_000, seconds=120, cpu_seconds=120
            )
            context = engine.SieveContext(
                config.max_hi, segment_size=config.segment_size, budget=budget
            )
            programs = (
                ECMPrograms(context, memory_bytes=config.ecm_program_bytes)
                if arm == "programs"
                else None
            )
            b1, b2, curves = config.ecm_tiers[0]
            curve_work = []
            outcomes = []
            for curve in range(curves):
                start_work = budget.used
                initializer = stages.new_job if arm == "baseline" else new_job
                job = initializer("ecm", n, 7 + curve, b1, b2)
                while not job["done"]:
                    if programs is None:
                        stages.advance_job(job, budget, context, config)
                    else:
                        advance_job(job, budget, programs, config)
                divisor = job["factor"]
                if divisor is not None and (
                    not 1 < divisor < n or n % divisor
                ):
                    raise AssertionError("invalid probe divisor")
                outcomes.append(divisor)
                curve_work.append(budget.used - start_work)
            rows.append(
                dict(
                    case=case,
                    arm=arm,
                    b1=b1,
                    b2=b2,
                    curves=curves,
                    input_bits=n.bit_length(),
                    max_input_bits=329,
                    curve_work=curve_work,
                    total_work=budget.used,
                    outcomes=outcomes,
                    workspace_reserve=config.workspace_reserve,
                    program_used_bytes=programs.used_bytes if programs else 0,
                    hits=programs.hits if programs else 0,
                    misses=programs.misses if programs else 0,
                )
            )
    for index in range(0, len(rows), 2):
        if rows[index]["outcomes"] != rows[index + 1]["outcomes"]:
            raise AssertionError("campaign probes changed curve outcomes")
    return rows


def worker(args):
    corpus = load_corpus()
    if args.case in SCHEDULE_CASES:
        corpus["schedule_oracle"] = reference_schedule(
            *SCHEDULE_CASES[args.case][:2]
        )
    if args.arm == "baseline":
        engine, stages = load_control()
    else:
        from ....execution import stage_jobs

        engine, stages = portfolio, stage_jobs
    started = time.perf_counter()
    count = 0
    while time.perf_counter() - started < args.warmup_seconds:
        measure(args.case, args.arm, corpus, engine, stages)
        count += 1
    warmup = dict(seconds=time.perf_counter() - started, cohorts=count)
    samples = [
        measure(args.case, args.arm, corpus, engine, stages)
        for _ in range(args.repetitions)
    ]
    # macOS reports ru_maxrss in bytes, unlike Linux's KiB.
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform != "darwin":
        rss *= 1024
    return dict(
        arm=args.arm,
        case=args.case,
        warmup=warmup,
        samples=samples,
        peak_rss_bytes=rss,
    )


def interval(before, after):
    """Unpaired bootstrap median reduction; independent worker processes."""
    generator = random.Random(52003)
    reductions = []
    for _ in range(3000):
        old = statistics.median(generator.choices(before, k=len(before)))
        new = statistics.median(generator.choices(after, k=len(after)))
        reductions.append(100 * (1 - new / old))
    reductions.sort()
    return [reductions[75], reductions[2924]]


def launch(case, arm, warmup, repetitions):
    command = [
        sys.executable,
        "-B",
        "-m",
        "v2.benchmarks.ecm.p52.p52_a3",
        "--worker",
        "--case",
        case,
        "--arm",
        arm,
        "--warmup-seconds",
        str(warmup),
        "--repetitions",
        str(repetitions),
    ]
    started = time.perf_counter()
    completed = subprocess.run(
        command, check=True, text=True, capture_output=True
    )
    return json.loads(completed.stdout), time.perf_counter() - started


def run(args):
    captures = []
    summaries = []
    for case in args.cases:
        times = {}
        signatures = []
        for arm in ARMS:
            attempts = []
            for attempt in range(3):
                warmup = (
                    args.warmup_seconds
                    if attempt == 0
                    else (5, 8)[attempt - 1]
                )
                repetitions = (
                    args.repetitions if attempt == 0 else (31, 63)[attempt - 1]
                )
                result, _ = launch(case, arm, warmup, repetitions)
                values = [sample["seconds"] for sample in result["samples"]]
                quartiles = statistics.quantiles(values, n=4)
                spread = (quartiles[2] - quartiles[0]) / statistics.median(
                    values
                )
                attempts.append(result)
                if spread <= 0.15:
                    break
            result = {
                **result,
                "prior_attempts": attempts[:-1],
                "relative_iqr": spread,
                "stable": spread <= 0.15,
            }
            captures.append(result)
            times[arm] = values
            first_rows = result["samples"][0]["rows"]
            if any(
                sample["rows"] != first_rows for sample in result["samples"]
            ):
                raise AssertionError(
                    "repeated samples changed outcomes or work"
                )
            signatures.append(
                [
                    {
                        key: value
                        for key, value in row.items()
                        if key
                        not in (
                            "work",
                            "retained_bytes",
                            "hits",
                            "misses",
                            "unretained",
                        )
                    }
                    for row in result["samples"][0]["rows"]
                ]
            )
            print(case, arm, round(statistics.median(values), 6), flush=True)
        if any(signature != signatures[0] for signature in signatures[1:]):
            raise AssertionError("matched arms changed outcomes or seeds")
        for arm in ARMS[1:]:
            summaries.append(
                dict(
                    case=case,
                    arm=arm,
                    baseline_seconds=statistics.median(times["baseline"]),
                    seconds=statistics.median(times[arm]),
                    reduction_percent=100
                    * (
                        1
                        - statistics.median(times[arm])
                        / statistics.median(times["baseline"])
                    ),
                    reduction_interval=interval(times["baseline"], times[arm]),
                )
            )
    cold = []
    if args.cold:
        for arm in ARMS:
            for _ in range(9):
                result, seconds = launch("small", arm, 0, 1)
                cold.append(
                    dict(arm=arm, seconds=seconds, sample=result["samples"][0])
                )
    output = dict(
        schema=1,
        runtime=sys.version,
        platform=platform.platform(),
        corpus_sha256=hashlib.sha256(
            source_path(CORPUS).read_bytes()
        ).hexdigest(),
        baseline_sha256=hashlib.sha256(
            source_path(BASELINE).read_bytes()
        ).hexdigest(),
        dependency_sha256=hashlib.sha256(
            source_path(DEPENDENCIES).read_bytes()
        ).hexdigest(),
        source_sha256=source_hashes(),
        probes=campaign_probes() if args.probes else [],
        summaries=summaries,
        captures=captures,
        cold=cold,
        limitations=[
            "Campaign timing measures finite curves, not factoring success",
            "RSS includes interpreter/JIT, corpus validation and warmup",
            "Separate processes; uncertainty uses an unpaired bootstrap",
            "Owned program workspace is distinct from process RSS",
        ],
    )
    with Path(args.output).open("x") as stream:
        json.dump(output, stream, indent=2)
        stream.write("\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--arm", choices=ARMS, default="programs")
    parser.add_argument("--case", choices=CASES, default="small")
    parser.add_argument(
        "--cases", choices=CASES, nargs="+", default=list(CASES)
    )
    parser.add_argument("--warmup-seconds", type=float, default=3)
    parser.add_argument("--repetitions", type=int, default=9)
    parser.add_argument("--cold", action="store_true")
    parser.add_argument("--probes", action="store_true")
    parser.add_argument("--output")
    args = parser.parse_args()
    if platform.python_implementation() != "PyPy" or sys.version_info[:2] != (
        3,
        11,
    ):
        parser.error("P5.2 measurements require PyPy implementing Python 3.11")
    if (
        not math.isfinite(args.warmup_seconds)
        or args.warmup_seconds < 0
        or args.repetitions < 1
    ):
        parser.error("warmup and repetitions must be finite and nonnegative")
    if not args.worker and (
        args.warmup_seconds < 3 or args.repetitions < 9 or not args.output
    ):
        parser.error("require >=3 seconds warmup, >=9 samples and an output")
    if args.worker:
        print(json.dumps(worker(args)))
    else:
        run(args)


if __name__ == "__main__":
    main()
