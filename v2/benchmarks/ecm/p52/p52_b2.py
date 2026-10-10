"""Frozen paired ECM comparisons; run only in an approved exclusive window."""

import argparse
import fcntl
import hashlib
import json
import os
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

from .... import portfolio
from ...suites.build_phase_two_corpus import verify_certificates
from ...support.paths import (
    REPOSITORY_ROOT,
    source_path,
)
from . import p52_a3, p52_realistic
from .build_p52_b2_inputs import INPUTS, MODULES

PROTOCOL = INPUTS / "controls/p52_b2_protocol.json"
CONTROL = INPUTS / "baselines/p52_b2_mainline.json"
CORPUS = INPUTS / "corpora/p52_b2_corpus.json"


def digest(path):
    return hashlib.sha256(source_path(path).read_bytes()).hexdigest()


def inputs():
    protocol = json.loads(source_path(PROTOCOL).read_text())
    if (
        digest(CONTROL) != protocol["controls_sha256"]
        or digest(CORPUS) != protocol["corpus_sha256"]
    ):
        raise ValueError("frozen B2 inputs changed")
    corpus = json.loads(source_path(CORPUS).read_text())
    verify_certificates(corpus["certificates"])
    for fixture in corpus["fixtures"]:
        product = 1
        for prime, exponent in fixture["factors"]:
            if str(prime) not in corpus["certificates"]:
                raise ValueError("missing independent factor proof")
            product *= prime**exponent
        if product != fixture["n"]:
            raise ValueError("corpus reconstruction failed")
    return protocol, corpus


def load_control():
    """Private A2/A3 engine; inactive QS types are annotations only."""
    if (
        digest(CONTROL)
        != json.loads(source_path(PROTOCOL).read_text())["controls_sha256"]
    ):
        raise ValueError("frozen B2 control changed")
    data = json.loads(source_path(CONTROL).read_text())
    package_name = "_p52_b2_control"
    package = types.ModuleType(package_name)
    package.__path__ = []
    sys.modules[package_name] = package
    from .... import qs
    from ....qs import sss

    sys.modules[package_name + ".qs"] = qs
    sys.modules[package_name + ".qs.sss"] = sss
    package.qs = qs
    for name in MODULES:
        source = data["source"][name]
        if hashlib.sha256(source.encode()).hexdigest() != data["sha256"][name]:
            raise ValueError("corrupt committed B2 control")
        qualified = package_name + "." + name
        module = types.ModuleType(qualified)
        module.__package__ = package_name
        module.__file__ = str(CONTROL) + ":" + name
        sys.modules[qualified] = module
        setattr(package, name, module)
        exec(compile(source, module.__file__, "exec"), module.__dict__)
    return package.portfolio


def source_hashes():
    root = REPOSITORY_ROOT
    names = [f"v2/{name}.py" for name in (*MODULES, "ecm_paired")]
    names += [
        "v2/benchmarks/" + name + ".py"
        for name in (
            "p52_b2",
            "build_p52_b2_inputs",
            "p52_a3",
            "p52_realistic",
            "build_phase_two_corpus",
        )
    ]
    return {name: digest(root / name) for name in names}


def options(protocol, case, arm):
    record = protocol["cases"][case]
    result = dict(
        trial_bound=5,
        rho_attempts=0,
        pm1_attempts=0,
        ecm_tiers=(tuple(record["tier"]),),
        backend="python-int",
        memory_bytes=protocol["memory_bytes"],
        max_input_bits=329,
        segment_size=protocol["segment_size"],
    )
    if arm != "streamed":
        result["ecm_program_bytes"] = protocol["program_bytes"]
    if arm.startswith("paired_") or arm == "regenerated":
        index = int(arm[-1]) if arm.startswith("paired_") else 1
        result["ecm_pair_distance"] = record["distances"][index]
    if arm == "regenerated":
        result["ecm_program_bytes"] = 4096 + 256 * protocol["segment_size"]
    return result


def measure(protocol, fixtures, engine, case, arm):
    """Whole factoring, including setup, table/recovery and output costs."""
    start, cpu_start = time.perf_counter(), time.process_time()
    config = engine.PortfolioConfig(**options(protocol, case, arm))
    rows = []
    for fixture in fixtures:
        for seed in protocol["seeds"]:
            run = engine.factorize_bounded(
                fixture["n"],
                seed=seed,
                config=config,
                budget=engine.Budget(
                    work_limit=protocol["work_limit"],
                    seconds=protocol["seconds"],
                    cpu_seconds=protocol["cpu_seconds"],
                ),
            )
            expected = Counter(dict(fixture["factors"]))
            actual = Counter({f.value: f.exponent for f in run.result.factors})
            if run.result.reconstruct() != fixture["n"] or actual - expected:
                raise AssertionError(
                    "incorrect complete/partial factorization"
                )
            if run.result.complete and actual != expected:
                raise AssertionError(
                    "terminal factors differ from certificates"
                )
            # Reconstruction includes unresolved cofactors and every split.
            for event in run.events:
                if event["n"] < 2 or fixture["n"] % event["n"]:
                    raise AssertionError("invalid recursive cofactor")
            if run.dropped_events:
                raise AssertionError(
                    "campaign exceeds reserved trace capacity"
                )
            rows.append(
                dict(
                    fixture=fixture["id"],
                    seed=seed,
                    complete=run.result.complete,
                    reason=run.reason,
                    factors=[
                        [f.value, f.exponent, f.certainty.value]
                        for f in run.result.factors
                    ],
                    remaining=list(run.result.remaining),
                    work=run.work_used,
                    curves=[
                        [e["seed"], e["outcome"]]
                        for e in run.events
                        if e["stage"] == "ecm"
                    ],
                )
            )
    return dict(
        seconds=time.perf_counter() - start,
        cpu_seconds=time.process_time() - cpu_start,
        rows=rows,
    )


def worker(args):
    protocol, corpus = inputs()
    engine = (
        load_control() if args.arm in ("streamed", "programs") else portfolio
    )
    fixtures = [
        f
        for f in corpus["fixtures"]
        if f["case"] == args.case and f["split"] == args.phase
    ]
    if not fixtures:
        raise ValueError("empty workload")
    start, count, signature = time.perf_counter(), 0, None
    deterministic = True
    while time.perf_counter() - start < args.warmup:
        sample = measure(protocol, fixtures, engine, args.case, args.arm)
        signature = sample["rows"] if signature is None else signature
        deterministic &= sample["rows"] == signature
        count += 1
    warmup = time.perf_counter() - start
    samples = [
        measure(protocol, fixtures, engine, args.case, args.arm)
        for _ in range(args.samples)
    ]
    signature = samples[0]["rows"] if signature is None else signature
    deterministic &= all(s["rows"] == signature for s in samples)
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return dict(
        case=args.case,
        arm=args.arm,
        phase=args.phase,
        config=options(protocol, args.case, args.arm),
        warmup_seconds=warmup,
        warmup_cohorts=count,
        deterministic=deterministic,
        samples=samples,
        peak_rss_bytes=rss if sys.platform == "darwin" else 1024 * rss,
    )


def launch(args, case, arm, warmup, samples):
    p52_realistic.check_quiet({os.getpid()})
    command = [
        sys.executable,
        "-B",
        "-m",
        "v2.benchmarks.ecm.p52.p52_b2",
        "--worker",
        "--phase",
        args.phase,
        "--case",
        case,
        "--arm",
        arm,
        "--warmup",
        str(warmup),
        "--samples",
        str(samples),
    ]
    started = time.perf_counter()
    process = subprocess.Popen(
        command, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    try:
        while True:
            try:
                stdout, stderr = process.communicate(timeout=1)
                break
            except subprocess.TimeoutExpired:
                p52_realistic.check_quiet({os.getpid(), process.pid})
        if process.returncode:
            raise RuntimeError(stderr)
        return json.loads(stdout), time.perf_counter() - started
    except BaseException:
        if process.poll() is None:
            process.terminate()
        process.communicate()
        raise


def stable(capture, protocol):
    samples = capture["samples"]
    return (
        capture["deterministic"]
        and p52_realistic.spread([s["seconds"] for s in samples])
        <= protocol["relative_iqr_limit"]
        and all(
            r["reason"] not in ("wall_limit", "cpu_limit")
            for s in samples
            for r in s["rows"]
        )
    )


def completion(capture):
    rows = capture["samples"][0]["rows"]
    return sum(r["complete"] for r in rows) / len(rows)


def completion_interval(control, candidate):
    """Cluster seed outcomes by fixture; reruns are not new inputs."""
    outcomes = {}
    for sign, capture in ((-1, control), (1, candidate)):
        for row in capture["samples"][0]["rows"]:
            key = row["fixture"]
            outcomes.setdefault(key, [0, 0])
            outcomes[key][0] += sign * int(row["complete"])
            outcomes[key][1] += int(sign == 1)
    differences = [100 * value / count for value, count in outcomes.values()]
    generator = random.Random(52009)
    samples = sorted(
        statistics.mean(generator.choices(differences, k=len(differences)))
        for _ in range(3000)
    )
    return [samples[75], samples[2924]]


def select(captures, protocol):
    """Freeze one D per tested bound tier before opening held-out timings."""
    selected = {}
    for case in protocol["cases"]:
        if case == "structured":
            continue
        controls = [
            c
            for c in captures
            if c["case"] == case and c["arm"] in ("streamed", "programs")
        ]
        candidates = [
            c
            for c in captures
            if c["case"] == case
            and c["arm"].startswith("paired_")
            and stable(c, protocol)
        ]
        if len(controls) != 2 or not all(
            stable(c, protocol) for c in controls
        ):
            selected[case] = None
            continue
        candidates = [
            c
            for c in candidates
            if completion(c)
            >= max(completion(control) for control in controls)
        ]
        selected[case] = (
            min(
                candidates,
                key=lambda c: (
                    statistics.median(s["seconds"] for s in c["samples"]),
                    c["arm"],
                ),
            )["arm"]
            if candidates
            else None
        )
    selected["structured"] = selected.get("small")
    return selected


def summaries(captures, protocol):
    rows = []
    for candidate in captures:
        case, arm = candidate["case"], candidate["arm"]
        if not arm.startswith("paired_"):
            continue
        for baseline in ("streamed", "programs"):
            control = next(
                c
                for c in captures
                if c["case"] == case and c["arm"] == baseline
            )
            before = [s["seconds"] for s in control["samples"]]
            after = [s["seconds"] for s in candidate["samples"]]
            rows.append(
                dict(
                    case=case,
                    arm=arm,
                    baseline=baseline,
                    median_reduction_percent=100
                    * (
                        1
                        - statistics.median(after) / statistics.median(before)
                    ),
                    conditional_timing_interval=p52_a3.interval(before, after),
                    completion_change_points=100
                    * (completion(candidate) - completion(control)),
                    fixture_cluster_completion_interval=completion_interval(
                        control, candidate
                    ),
                    stable=stable(candidate, protocol)
                    and stable(control, protocol),
                    scope="finite campaign"
                    if case == "campaign"
                    else "whole factoring",
                )
            )
    return rows


def diagnostics(protocol, corpus):
    """Instrument finite curves separately, without time ratios."""
    rows = []
    control = load_control()
    for case in protocol["cases"]:
        fixture = next(
            f
            for f in corpus["fixtures"]
            if f["case"] == case and f["split"] == "training"
        )
        for arm in protocol["training_arms"]:
            engine = control if arm in ("streamed", "programs") else portfolio
            config = engine.PortfolioConfig(**options(protocol, case, arm))
            budget = engine.Budget(
                work_limit=protocol["work_limit"],
                seconds=protocol["seconds"],
                cpu_seconds=protocol["cpu_seconds"],
            )
            context = engine.SieveContext(
                config.max_hi, segment_size=config.segment_size, budget=budget
            )
            programs = (
                engine.ECMPrograms(
                    context, memory_bytes=config.ecm_program_bytes
                )
                if config.ecm_program_bytes
                else None
            )
            b1, b2, curves = config.ecm_tiers[0]
            for curve in range(min(3, curves)):
                start = budget.used
                seed = protocol["seeds"][0] + curve
                job = engine.new_job("ecm", fixture["n"], seed, b1, b2)
                counts = Counter()
                peak_table = 0
                while not job["done"]:
                    phase = job["phase"]
                    terms, center = len(job["terms"]), job.get("center")
                    engine.advance_job(
                        job, budget, programs or context, config
                    )
                    counts[phase + "_actions"] += 1
                    added = max(0, len(job["terms"]) - terms)
                    counts["term_products"] += added
                    if (
                        phase in ("stage_two", "pair_terms")
                        and job.get("center") != center
                    ):
                        counts["giant_advances"] += 1
                    peak_table = max(peak_table, len(job.get("baby", [])))
                    if job.get("term_records") and added:
                        counts["certified_primes"] += sum(
                            bool(record[2]) + bool(record[3])
                            for record in job["term_records"][-added:]
                        )
                divisor = job["factor"]
                if divisor is not None and not (
                    1 < divisor < fixture["n"] and fixture["n"] % divisor == 0
                ):
                    raise AssertionError("invalid diagnostic divisor")
                rows.append(
                    dict(
                        case=case,
                        arm=arm,
                        fixture=fixture["id"],
                        seed=seed,
                        curve=curve,
                        divisor=divisor,
                        counts=dict(counts),
                        curve_work=budget.used - start,
                        cumulative_work=budget.used,
                        peak_table_slots=peak_table,
                        workspace_reserve=config.workspace_reserve,
                        program_used_bytes=programs.used_bytes
                        if programs
                        else 0,
                        hits=programs.hits if programs else 0,
                        coverage_hits=getattr(programs, "coverage_hits", 0),
                        coverage_misses=getattr(
                            programs, "coverage_misses", 0
                        ),
                    )
                )
    return rows


def run(args):
    protocol, corpus = inputs()
    pins = source_hashes()
    selection = None
    if args.phase == "held_out":
        if args.selection is None:
            raise ValueError("held-out runs require frozen training selection")
        training = json.loads(source_path(args.selection).read_text())
        if (
            training["phase"] != "training"
            or training["source_sha256"] != pins
            or training["protocol_sha256"] != digest(PROTOCOL)
        ):
            raise ValueError("training selection is incompatible")
        selection = training["selected"]
    report = dict(
        phase=args.phase,
        source_sha256=pins,
        protocol_sha256=digest(PROTOCOL),
        protocol=protocol,
        runtime=sys.version,
        machine=platform.platform(),
        captures=[],
        attempts=[],
        cold=[],
        status="running",
    )
    with open(protocol["machine_lock"], "a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        owner = Path("/private/tmp/factor-performance-owner.json")
        owner.write_text(json.dumps(dict(owner="B2 ECM", pid=os.getpid())))
        try:
            with args.output.open("x") as output:

                def save():
                    output.seek(0)
                    json.dump(report, output, indent=2)
                    output.truncate()
                    output.flush()

                save()
                for case in protocol["cases"]:
                    arms = list(protocol["training_arms"])
                    if selection is not None:
                        arms = ["streamed", "programs"]
                        if selection.get(case):
                            arms.append(selection[case])
                    random.Random(2026100952).shuffle(arms)
                    for arm in arms:
                        for warmup, samples in protocol["sampling"]:
                            capture, _ = launch(
                                args, case, arm, warmup, samples
                            )
                            report["attempts"].append(capture)
                            save()
                            if stable(capture, protocol):
                                break
                        report["captures"].append(capture)
                        save()
                        print(
                            case,
                            arm,
                            "stable"
                            if stable(capture, protocol)
                            else "inconclusive",
                            flush=True,
                        )
                        if args.cold:
                            for _ in range(9):
                                cold, elapsed = launch(args, case, arm, 0, 1)
                                report["cold"].append(
                                    dict(seconds=elapsed, capture=cold)
                                )
                                save()
                if args.diagnostics:
                    p52_realistic.check_quiet({os.getpid()})
                    report["instrumented_diagnostics"] = diagnostics(
                        protocol, corpus
                    )
                    save()
                if pins != source_hashes():
                    raise RuntimeError("sources changed during comparisons")
                report.update(
                    selected=select(report["captures"], protocol)
                    if selection is None
                    else selection,
                    comparisons=summaries(report["captures"], protocol),
                    status="complete",
                    decision="Retain defaults pending held-out gate review.",
                )
                save()
        finally:
            owner.unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--phase", choices=("training", "held_out"), default="training"
    )
    parser.add_argument("--output", type=Path)
    parser.add_argument("--selection", type=Path)
    parser.add_argument("--cold", action="store_true")
    parser.add_argument("--diagnostics", action="store_true")
    parser.add_argument(
        "--worker", action="store_true", help=argparse.SUPPRESS
    )
    parser.add_argument("--case")
    parser.add_argument("--arm")
    parser.add_argument("--warmup", type=float, default=3)
    parser.add_argument("--samples", type=int, default=9)
    args = parser.parse_args()
    if platform.python_implementation() != "PyPy" or sys.version_info[:2] != (
        3,
        11,
    ):
        parser.error("requires PyPy implementing Python 3.11")
    if args.worker:
        print(json.dumps(worker(args)))
    elif args.output is None:
        parser.error("an unused output path is required")
    else:
        run(args)


if __name__ == "__main__":
    main()
