"""Full ECM-only portfolios on certified 30--80 digit workload strata."""

import argparse
import hashlib
import json
import math
import os
import platform
import resource
import shlex
import statistics
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path

from v2 import portfolio
from v2.benchmarks import p52_a3
from v2.benchmarks.build_phase_two_corpus import verify_certificates

CORPUS = p52_a3.INPUTS / "corpora/p52_realistic_corpus.json"
ARMS = ("baseline", "default", "programs")
# Pretests share total allowances; campaigns share declared curve counts.
POLICIES = {
    "pretest": ((2000, 147396, 32), 2_000_000),
    "middle": ((11000, 1900000, 8), 50_000_000),
    "large": ((50000, 5000000, 4), 50_000_000),
    "deep": ((11000, 1900000, 256), 400_000_000),
}
DEFAULT_POLICIES = ("pretest", "middle", "large")


def load_corpus():
    """Verify proofs, factor references, reconstruction and declared strata."""
    corpus = json.loads(CORPUS.read_text())
    verify_certificates(corpus["certificates"])
    for fixture in corpus["fixtures"]:
        product = 1
        for prime, exponent in fixture["factors"]:
            if str(prime) not in corpus["certificates"]:
                raise ValueError("factor lacks an independent certificate")
            product *= prime**exponent
        if (
            product != fixture["n"]
            or len(str(product)) != fixture["digits"]
            or len(str(min(p for p, _ in fixture["factors"])))
            != fixture["small_digits"]
        ):
            raise ValueError("corpus reconstruction or stratum mismatch")
    return corpus


def options(policy, arm):
    """Match all controls except the explicitly requested program cap."""
    return dict(
        trial_bound=5,
        rho_attempts=0,
        pm1_attempts=0,
        ecm_tiers=(POLICIES[policy][0],),
        max_input_bits=329,
        memory_bytes=16 * 2**20,
        **({"ecm_program_bytes": 8 * 2**20} if arm == "programs" else {}),
    )


def signature(rows, *, work=True):
    """Keep deterministic results and certainty apart from timing evidence."""
    ignored = {"seconds", "cpu_seconds"}
    if not work:
        ignored |= {"work", "events"}
    return [{k: v for k, v in row.items() if k not in ignored} for row in rows]


def measure(policy, arm, corpus, engine):
    """Include preprocessing, classification, ECM and recursive validation."""
    config = engine.PortfolioConfig(**options(policy, arm))
    rows = []
    cohort_started = time.perf_counter()
    for fixture in corpus["fixtures"]:
        for seed in corpus["seeds"]:
            budget = engine.Budget(
                work_limit=POLICIES[policy][1], seconds=120, cpu_seconds=120
            )
            started = time.perf_counter()
            cpu_started = time.process_time()
            run = engine.factorize_bounded(
                fixture["n"], seed=seed, config=config, budget=budget
            )
            elapsed = time.perf_counter() - started
            cpu = time.process_time() - cpu_started
            expected = Counter(dict(fixture["factors"]))
            actual = Counter({f.value: f.exponent for f in run.result.factors})
            if run.result.reconstruct() != fixture["n"]:
                raise AssertionError("whole portfolio does not reconstruct")
            if actual - expected or (
                run.result.complete and actual != expected
            ):
                raise AssertionError("terminal factors disagree with proofs")
            if run.reason in ("wall_limit", "cpu_limit"):
                raise AssertionError(
                    "time censoring prevents matched evidence"
                )
            if run.dropped_events:
                raise AssertionError(
                    "curve trace exceeded its finite capacity"
                )
            events = [
                [e["seed"], e["outcome"]]
                for e in run.events
                if e["stage"] == "ecm"
            ]
            current = run.checkpoint["payload"]["state"]["current"]
            active_job = current.get("job") if current else None
            active_ecm = (
                active_job
                if active_job and active_job["kind"] == "ecm"
                else None
            )
            rows.append(
                dict(
                    fixture=fixture["id"],
                    digits=fixture["digits"],
                    shape=fixture["shape"],
                    small_digits=fixture["small_digits"],
                    seed=seed,
                    complete=run.result.complete,
                    reason=run.reason,
                    factors=[
                        [f.value, f.exponent, f.certainty.value]
                        for f in run.result.factors
                    ],
                    remaining=run.result.remaining,
                    events=events,
                    active_seed=active_ecm["seed"] if active_ecm else None,
                    paused_phase=active_ecm["phase"] if active_ecm else None,
                    work=run.work_used,
                    seconds=elapsed,
                    cpu_seconds=cpu,
                )
            )
    return dict(seconds=time.perf_counter() - cohort_started, rows=rows)


def worker(args):
    corpus = load_corpus()
    if args.fixtures:
        available = {fixture["id"] for fixture in corpus["fixtures"]}
        if set(args.fixtures) - available:
            raise ValueError("unknown workload fixture")
        corpus["fixtures"] = [
            fixture
            for fixture in corpus["fixtures"]
            if fixture["id"] in args.fixtures
        ]
    if args.seeds is not None:
        corpus["seeds"] = args.seeds
    engine = p52_a3.load_control()[0] if args.arm == "baseline" else portfolio
    warmup_started = time.perf_counter()
    warmed = []
    while time.perf_counter() - warmup_started < args.warmup_seconds:
        warmed.append(
            signature(measure(args.policy, args.arm, corpus, engine)["rows"])
        )
    warmup_elapsed = time.perf_counter() - warmup_started
    samples = [
        measure(args.policy, args.arm, corpus, engine)
        for _ in range(args.repetitions)
    ]
    expected = signature(samples[0]["rows"])
    if any(rows != expected for rows in warmed) or any(
        signature(sample["rows"]) != expected for sample in samples
    ):
        raise AssertionError(
            "repeat changed results, work or curve assignments"
        )
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return dict(
        policy=args.policy,
        arm=args.arm,
        seeds=corpus["seeds"],
        config=options(args.policy, args.arm),
        work_limit=POLICIES[args.policy][1],
        warmup_seconds=warmup_elapsed,
        warmup_cohorts=len(warmed),
        samples=samples,
        peak_rss_bytes=rss if sys.platform == "darwin" else 1024 * rss,
    )


def spread(values):
    quartiles = statistics.quantiles(values, n=4)
    return (quartiles[2] - quartiles[0]) / statistics.median(values)


def competing_processes(listing, allowed):
    """Select other benchmark/test interpreters from a scoped process list."""
    competing = []
    for line in listing.splitlines():
        fields = line.strip().split(None, 2)
        if len(fields) != 3:
            continue
        process_id, _, command = fields
        try:
            process_id = int(process_id)
            arguments = shlex.split(command)
        except ValueError:
            continue
        if process_id in allowed or not arguments:
            continue
        interpreter = Path(arguments[0]).name
        if not interpreter.startswith(("python", "pypy")):
            continue
        if "-m" not in arguments:
            continue
        position = arguments.index("-m") + 1
        if position >= len(arguments):
            continue
        module = arguments[position]
        if module.startswith("v2.benchmarks.") or module == "unittest":
            competing.append(dict(pid=process_id, command=command))
    return competing


def check_quiet(allowed):
    """Fail closed if the process inventory is unavailable or overlapping."""
    listing = subprocess.check_output(
        ["ps", "-Ao", "pid=,ppid=,command="], text=True
    )
    competing = competing_processes(listing, allowed)
    if competing:
        raise RuntimeError("benchmark/test overlap: " + json.dumps(competing))


def launch(
    policy, arm, warmup, repetitions, fixtures=None, quiet=False, seeds=None
):
    command = [
        sys.executable,
        "-B",
        "-m",
        "v2.benchmarks.p52_realistic",
        "--worker",
        "--policy",
        policy,
        "--arm",
        arm,
        "--warmup-seconds",
        str(warmup),
        "--repetitions",
        str(repetitions),
    ]
    if fixtures:
        command.extend(("--fixtures", *fixtures))
    if seeds is not None:
        command.extend(("--seeds", *(str(seed) for seed in seeds)))
    if quiet:
        check_quiet({os.getpid()})
    process = subprocess.Popen(
        command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    try:
        while True:
            try:
                stdout, stderr = process.communicate(timeout=1)
                break
            except subprocess.TimeoutExpired:
                if quiet:
                    check_quiet({os.getpid(), process.pid})
        if quiet:
            check_quiet({os.getpid(), process.pid})
        if process.returncode:
            raise RuntimeError("ECM worker failed: " + stderr)
    except BaseException:
        # Only this runner's worker is terminated; partial captures stay local.
        process.terminate()
        process.communicate()
        raise
    return json.loads(stdout)


def summaries(captures):
    """Report finite-cohort costs and completion without imputed success."""
    output = []
    for policy in POLICIES:
        arms = {c["arm"]: c for c in captures if c["policy"] == policy}
        if "baseline" not in arms:
            continue
        for arm, capture in arms.items():
            baseline = arms["baseline"]
            before = [s["seconds"] for s in baseline["samples"]]
            after = [s["seconds"] for s in capture["samples"]]
            rows = capture["samples"][0]["rows"]
            grouped = []
            for digits in (30, 40, 60, 80):
                for shape in ("small10", "target", "balanced"):
                    selected = [
                        r
                        for r in rows
                        if r["digits"] == digits and r["shape"] == shape
                    ]
                    if not selected:
                        continue
                    indices = [i for i, r in enumerate(rows) if r in selected]
                    values = [
                        sum(sample["rows"][i]["seconds"] for i in indices)
                        for sample in capture["samples"]
                    ]
                    control = [
                        sum(sample["rows"][i]["seconds"] for i in indices)
                        for sample in baseline["samples"]
                    ]
                    grouped.append(
                        dict(
                            digits=digits,
                            shape=shape,
                            small_digits=selected[0]["small_digits"],
                            complete=sum(r["complete"] for r in selected),
                            attempts=len(selected),
                            seconds=statistics.median(values),
                            baseline_seconds=statistics.median(control),
                            relative_iqr=spread(values),
                            stable=spread(values) <= 0.15
                            and spread(control) <= 0.15,
                            work=[r["work"] for r in selected],
                            curves=[len(r["events"]) for r in selected],
                            started_curves=[
                                len(r["events"])
                                + int(r["active_seed"] is not None)
                                for r in selected
                            ],
                            reasons=[r["reason"] for r in selected],
                        )
                    )
            output.append(
                dict(
                    policy=policy,
                    arm=arm,
                    complete=sum(r["complete"] for r in rows),
                    attempts=len(rows),
                    seconds=statistics.median(after),
                    reduction_percent=100
                    * (
                        1
                        - statistics.median(after) / statistics.median(before)
                    ),
                    reduction_interval=p52_a3.interval(before, after),
                    stable=capture["stable"] and baseline["stable"],
                    work=sum(r["work"] for r in rows),
                    curves=sum(len(r["events"]) for r in rows),
                    started_curves=sum(
                        len(r["events"]) + int(r["active_seed"] is not None)
                        for r in rows
                    ),
                    grouped=grouped,
                )
            )
    return output


def run(args):
    captures = []
    for policy in args.policies:
        for arm in ARMS:
            attempts = []
            for warmup, repetitions in (
                (args.warmup_seconds, args.repetitions),
                (5, 31),
                (8, 63),
            ):
                capture = launch(
                    policy,
                    arm,
                    warmup,
                    repetitions,
                    args.fixtures,
                    args.check_quiet,
                    args.seeds,
                )
                values = [s["seconds"] for s in capture["samples"]]
                attempts.append(capture)
                if spread(values) <= 0.15:
                    break
            capture.update(
                prior_attempts=attempts[:-1],
                relative_iqr=spread(values),
                stable=spread(values) <= 0.15,
            )
            captures.append(capture)
            print(
                policy,
                arm,
                round(statistics.median(values), 4),
                "complete",
                sum(r["complete"] for r in capture["samples"][0]["rows"]),
                flush=True,
            )
            partial = Path(str(args.output) + f".{policy}-{arm}.json")
            with partial.open("x") as stream:
                json.dump(capture, stream, indent=2)
                stream.write("\n")
        current = {c["arm"]: c for c in captures if c["policy"] == policy}
        baseline_rows = current["baseline"]["samples"][0]["rows"]
        if signature(baseline_rows) != signature(
            current["default"]["samples"][0]["rows"]
        ):
            raise AssertionError(
                "default path changed actions, certainty or work"
            )
        # Equal grants can fund different coverage under a new work ledger.
        # Fixed-curve comparisons must still preserve every assignment/outcome.
        if policy != "pretest":
            program_rows = current["programs"]["samples"][0]["rows"]
            if signature(baseline_rows, work=False) != signature(
                program_rows, work=False
            ):
                raise AssertionError("fixed-curve campaigns changed outcomes")
            if [r["events"] for r in baseline_rows] != [
                r["events"] for r in program_rows
            ]:
                raise AssertionError(
                    "fixed-curve campaigns changed assignments"
                )
    output = dict(
        schema=1,
        runtime=sys.version,
        platform=platform.platform(),
        corpus_sha256=hashlib.sha256(CORPUS.read_bytes()).hexdigest(),
        baseline_sha256=hashlib.sha256(
            p52_a3.BASELINE.read_bytes()
        ).hexdigest(),
        dependency_sha256=hashlib.sha256(
            p52_a3.DEPENDENCIES.read_bytes()
        ).hexdigest(),
        source_sha256=p52_a3.source_hashes("p52_realistic"),
        selected_fixtures=args.fixtures,
        selected_seeds=args.seeds,
        quiet_monitor=args.check_quiet,
        captures=captures,
        summaries=summaries(captures),
        limitations=[
            "One certified input per stratum, two seeds; fixed exploration.",
            "Pretest uses equal work grants; programs can fund more curves.",
            "Campaigns compare equal declared curves, not factoring success.",
            "Unresolved elapsed times are search costs, not time to factor.",
            "RSS includes runtime/JIT, proof validation and warmup.",
        ],
    )
    with Path(args.output).open("x") as stream:
        json.dump(output, stream, indent=2)
        stream.write("\n")


def coverage_sweep(args):
    """Compare seeded completion, with no single-sample speed estimate."""
    captures = []
    for policy in args.policies:
        for arm in ("baseline", "programs"):
            capture = launch(
                policy, arm, 0, 1, args.fixtures, args.check_quiet, args.seeds
            )
            captures.append(capture)
            with Path(str(args.output) + f".{policy}-{arm}.json").open(
                "x"
            ) as stream:
                json.dump(capture, stream, indent=2)
                stream.write("\n")
            rows = capture["samples"][0]["rows"]
            print(
                policy,
                arm,
                sum(r["complete"] for r in rows),
                len(rows),
                flush=True,
            )
        before, after = (c["samples"][0]["rows"] for c in captures[-2:])
        if signature(before, work=False) != signature(after, work=False):
            raise AssertionError("coverage sweep changed outcomes")
        if [r["events"] for r in before] != [r["events"] for r in after]:
            raise AssertionError("coverage sweep changed assignments")
    output = dict(
        schema=1,
        kind="seeded-completion-only",
        runtime=sys.version,
        quiet_monitor=args.check_quiet,
        selected_fixtures=args.fixtures,
        selected_seeds=args.seeds,
        corpus_sha256=hashlib.sha256(CORPUS.read_bytes()).hexdigest(),
        baseline_sha256=hashlib.sha256(
            p52_a3.BASELINE.read_bytes()
        ).hexdigest(),
        dependency_sha256=hashlib.sha256(
            p52_a3.DEPENDENCIES.read_bytes()
        ).hexdigest(),
        source_sha256=p52_a3.source_hashes("p52_realistic"),
        captures=captures,
        limitations=[
            "No warmup or repeated samples; no speed claims.",
            "Additional seeds share fixed certified inputs.",
        ],
    )
    with Path(args.output).open("x") as stream:
        json.dump(output, stream, indent=2)
        stream.write("\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--policy", choices=POLICIES, default="middle")
    parser.add_argument(
        "--policies",
        choices=POLICIES,
        nargs="+",
        default=list(DEFAULT_POLICIES),
    )
    parser.add_argument("--arm", choices=ARMS, default="programs")
    parser.add_argument("--fixtures", nargs="+")
    parser.add_argument("--seeds", nargs="+", type=int)
    parser.add_argument("--coverage-only", action="store_true")
    parser.add_argument("--check-quiet", action="store_true")
    parser.add_argument("--warmup-seconds", type=float, default=3)
    parser.add_argument("--repetitions", type=int, default=9)
    parser.add_argument("--output")
    args = parser.parse_args()
    if platform.python_implementation() != "PyPy" or sys.version_info[:2] != (
        3,
        11,
    ):
        parser.error("require PyPy implementing Python 3.11")
    if (
        not math.isfinite(args.warmup_seconds)
        or args.warmup_seconds < 0
        or args.repetitions < 1
    ):
        parser.error("invalid warmup or repetitions")
    if (
        not args.worker
        and not args.coverage_only
        and (
            args.warmup_seconds < 3 or args.repetitions < 9 or not args.output
        )
    ):
        parser.error("require >=3 seconds warmup, >=9 samples and an output")
    if args.coverage_only and (
        args.worker or not args.output or "pretest" in args.policies
    ):
        parser.error("coverage sweep requires fixed-curve policies and output")
    if args.worker:
        print(json.dumps(worker(args)))
    elif args.coverage_only:
        coverage_sweep(args)
    else:
        run(args)


if __name__ == "__main__":
    main()
