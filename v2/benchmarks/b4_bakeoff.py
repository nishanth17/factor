"""Fresh B4 bakeoff with counterbalanced independent PyPy process blocks."""

import argparse
import fcntl
import json
import os
import platform
import random
import resource
import statistics
import subprocess
import sys
import time
from pathlib import Path

from . import b4_common as common
from . import b4_kernels as kernels
from . import b4_study as study
from . import p52_realistic
from .build_phase_two_corpus import verify_certificates

ROOT = Path(__file__).parents[2]
PROTOCOL = common.INPUTS / "controls/b4_bakeoff_protocol.json"
CORPUS = common.INPUTS / "corpora/b4_bakeoff_corpus.json"
FREEZE = common.INPUTS / "controls/b4_bakeoff_freeze.json"


def verify_inputs():
    """Preserve the first study and verify the new independent freeze."""
    study.verify_freeze()
    original, old_corpus = common.inputs()
    pins = json.loads(FREEZE.read_text())
    for name, expected in pins["source_hashes"].items():
        if common.digest(ROOT / name) != expected:
            raise ValueError("bakeoff source/input changed: " + name)
    protocol = json.loads(PROTOCOL.read_text())
    corpus = json.loads(CORPUS.read_text())
    verify_certificates(corpus["certificates"])
    old_numbers = {f["n"] for f in old_corpus["fixtures"]}
    numbers, identifiers = set(), set()
    for fixture in corpus["fixtures"]:
        product = 1
        for prime, exponent in fixture["factors"]:
            if str(prime) not in corpus["certificates"]:
                raise ValueError("missing independent certificate")
            product *= prime**exponent
        if (
            product != fixture["n"]
            or product in old_numbers | numbers
            or fixture["id"] in identifiers
            or product.bit_length() > protocol["max_input_bits"]
        ):
            raise ValueError("invalid or overlapping bakeoff fixture")
        numbers.add(product)
        identifiers.add(fixture["id"])
    if (
        protocol["revision"] != original["revision"]
        or pins["arms"] != list(kernels.ARMS)
        or {f["split"] for f in corpus["fixtures"]}
        != {"screen", "confirmation"}
    ):
        raise ValueError("bakeoff control identity changed")
    return protocol, corpus, pins


def block_order(arms, backend, block, seed):
    """Rotate a fixed shuffled order; reverse after each complete cycle."""
    order = list(arms)
    generator = random.Random(seed + sum(map(ord, backend)))
    generator.shuffle(order)
    offset = block % len(order)
    order = order[offset:] + order[:offset]
    return order if block // len(order) % 2 == 0 else list(reversed(order))


def paired_gain(before, after):
    """Bootstrap matched process blocks, conditional on the fixed cohort."""
    if len(before) != len(after) or len(before) < 9:
        raise ValueError("need at least nine complete paired blocks")
    ratios = [new / old for old, new in zip(before, after)]
    generator = random.Random(2026100957)
    boot = sorted(
        100 * (1 - statistics.median(generator.choices(ratios, k=len(ratios))))
        for _ in range(5000)
    )
    return 100 * (1 - statistics.median(ratios)), [boot[125], boot[4874]]


def capture_map(captures, backend, arm):
    return {
        capture["block"]: capture
        for capture in captures
        if capture["backend"] == backend and capture["arm"] == arm
    }


def compare(captures, backend, arm, protocol):
    """Assess a fixed candidate without changing the declared gate."""
    before = capture_map(captures, backend, "baseline")
    after = capture_map(captures, backend, arm)
    blocks = sorted(after)
    if not set(blocks) <= set(before) or len(blocks) < 9:
        raise ValueError("unpaired bakeoff captures")
    old_times = [before[b]["seconds"] for b in blocks]
    new_times = [after[b]["seconds"] for b in blocks]
    gain, interval = paired_gain(old_times, new_times)
    cpu_gain, cpu_interval = paired_gain(
        [before[b]["cpu_seconds"] for b in blocks],
        [after[b]["cpu_seconds"] for b in blocks],
    )
    ratios = [new / old for old, new in zip(old_times, new_times)]
    spreads = [
        p52_realistic.spread(values)
        for values in (old_times, new_times, ratios)
    ]
    stable = max(spreads) <= protocol["relative_iqr_limit"]
    split = len(blocks) // 2
    half_gains = [
        100 * (1 - statistics.median(values))
        for values in (ratios[:split], ratios[split:])
    ]
    old_rows, new_rows = before[blocks[0]]["rows"], after[blocks[0]]["rows"]
    # Validate deterministic results across fresh processes as well as within
    # each process. Comparing only the first block could hide a varying exit.
    for captures_by_block, expected in ((before, old_rows), (after, new_rows)):
        for block in blocks:
            if study.signature(captures_by_block[block]["rows"]) != (
                study.signature(expected)
            ):
                raise AssertionError("cross-process bakeoff outcome changed")
    classes = []
    for case in sorted({row["case"] for row in old_rows}):
        old_class = [
            sum(r["seconds"] for r in before[b]["rows"] if r["case"] == case)
            for b in blocks
        ]
        new_class = [
            sum(r["seconds"] for r in after[b]["rows"] if r["case"] == case)
            for b in blocks
        ]
        class_gain, class_interval = paired_gain(old_class, new_class)
        old_done = [r["complete"] for r in old_rows if r["case"] == case]
        new_done = [r["complete"] for r in new_rows if r["case"] == case]
        classes.append(
            dict(
                case=case,
                reduction_percent=class_gain,
                conditional_interval=class_interval,
                completion_change_points=100
                * (statistics.mean(new_done) - statistics.mean(old_done)),
            )
        )
    completion_ok = all(
        row["completion_change_points"] >= -5 for row in classes
    )
    # Nine independent process blocks support only fixed-cohort timing
    # uncertainty. Held-out inputs, not additional repeats of screening
    # inputs, protect the eventual selection from a lucky training winner.
    timing_pass = (
        stable
        and completion_ok
        and gain > 0
        and interval[0] > 0
        and cpu_gain > 0
        and min(half_gains) > 0
    )
    return dict(
        backend=backend,
        arm=arm,
        samples=len(blocks),
        control_seconds=statistics.median(old_times),
        candidate_seconds=statistics.median(new_times),
        reduction_percent=gain,
        conditional_paired_interval=interval,
        cpu_reduction_percent=cpu_gain,
        cpu_conditional_paired_interval=cpu_interval,
        first_second_half_reductions=half_gains,
        relative_iqr=dict(
            baseline=spreads[0], candidate=spreads[1], paired_ratio=spreads[2]
        ),
        stable=stable,
        completion_ok=completion_ok,
        timing_evidence_passes=timing_pass,
        complete=[
            sum(r["complete"] for r in rows) for rows in (old_rows, new_rows)
        ],
        trials=len(old_rows),
        matched_outcomes_and_work=(
            study.signature(old_rows) == study.signature(new_rows)
        ),
        classes=classes,
    )


def select(rows):
    """Select one positive stable timing candidate per backend."""
    selected = {}
    for backend in {row["backend"] for row in rows}:
        eligible = [
            row
            for row in rows
            if row["backend"] == backend
            and row["stable"]
            and row["completion_ok"]
            and row["reduction_percent"] > 0
        ]
        selected[backend] = (
            max(eligible, key=lambda row: row["reduction_percent"])["arm"]
            if eligible
            else None
        )
    return selected


def worker(args):
    protocol, corpus, _ = verify_inputs()
    engine = kernels.engine(args.arm)
    fixtures = [f for f in corpus["fixtures"] if f["split"] == args.phase]
    expected = None

    def measure():
        nonlocal expected
        started, cpu = time.perf_counter(), time.process_time()
        rows = study.measure_full(engine, protocol, args.backend, fixtures)
        elapsed, used = (
            time.perf_counter() - started,
            time.process_time() - cpu,
        )
        current = study.signature(rows)
        if expected is not None and current != expected:
            raise AssertionError("deterministic bakeoff outcome changed")
        expected = current
        return rows, elapsed, used

    start, warmups = time.perf_counter(), 0
    while time.perf_counter() - start < args.warmup:
        measure()
        warmups += 1
    warm_elapsed = time.perf_counter() - start
    samples, elapsed, cpu = [], 0.0, 0.0
    while not samples or cpu < protocol["sample_cpu_seconds"]:
        rows, wall_used, cpu_used = measure()
        samples.append(rows)
        elapsed += wall_used
        cpu += cpu_used
    averaged = []
    for index, row in enumerate(samples[0]):
        row = dict(row)
        row["seconds"] = statistics.mean(
            sample[index]["seconds"] for sample in samples
        )
        averaged.append(row)
    return dict(
        phase=args.phase,
        arm=args.arm,
        backend=args.backend,
        block=args.block,
        warmup_seconds=warm_elapsed,
        warmup_cohorts=warmups,
        measured_cohorts=len(samples),
        seconds=elapsed / len(samples),
        cpu_seconds=cpu / len(samples),
        rows=averaged,
        rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    )


def launch(args, arm, backend, block, warmup, deadline):
    p52_realistic.check_quiet({os.getpid()})
    command = [
        sys.executable,
        "-B",
        "-m",
        "v2.benchmarks.b4_bakeoff",
        "--worker",
        "--phase",
        args.phase,
        "--arm",
        arm,
        "--backend",
        backend,
        "--block",
        str(block),
        "--warmup",
        str(warmup),
    ]
    process = subprocess.Popen(
        command, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    started = time.perf_counter()
    try:
        while True:
            try:
                stdout, stderr = process.communicate(timeout=1)
                break
            except subprocess.TimeoutExpired:
                p52_realistic.check_quiet({os.getpid(), process.pid})
                now = time.perf_counter()
                if now - started > 600 or now > deadline:
                    raise TimeoutError("bounded bakeoff worker/window expired")
        if process.returncode:
            raise RuntimeError(stderr)
        p52_realistic.check_quiet({os.getpid()})
        return json.loads(stdout)
    finally:
        if process.poll() is None:
            process.terminate()
        process.communicate()


def save(output, report):
    output.seek(0)
    json.dump(report, output, indent=2)
    output.truncate()
    output.flush()


def run(args):
    protocol, _, pins = verify_inputs()
    selected, selection = None, None
    if args.phase == "confirmation":
        selection = json.loads(args.selection.read_text())
        if (
            selection["freeze"] != pins
            or selection["phase"] != "screen"
            or selection["selected"] != select(selection["summary"])
        ):
            raise ValueError("selection does not match the frozen screen")
        selected = selection["selected"]
    arms = {
        backend: list(kernels.ARMS)
        if selected is None
        else (["baseline", selected[backend]] if selected[backend] else [])
        for backend in protocol["backends"]
    }
    report = dict(
        phase=args.phase,
        freeze=pins,
        protocol_sha256=common.digest(PROTOCOL),
        runtime=sys.version,
        executable=sys.executable,
        platform=platform.platform(),
        captures=[],
        summary=[],
        rounds=[],
        selection_sha256=common.digest(args.selection) if selection else None,
        source_revision=subprocess.check_output(
            ["git", "rev-parse", "HEAD"], text=True
        ).strip(),
    )
    p52_realistic.check_quiet({os.getpid()})
    with open(protocol["machine_lock"], "a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        owner = Path("/private/tmp/factor-performance-owner.json")
        owner.write_text(
            json.dumps(
                dict(
                    owner="B4 fresh bakeoff",
                    pid=os.getpid(),
                    phase=args.phase,
                    worktree=str(ROOT),
                )
            )
        )
        deadline = time.perf_counter() + protocol["window_seconds"]
        active = dict(arms)
        try:
            with args.output.open("x") as output:
                save(output, report)
                completed = 0
                for warmup, target in protocol["sampling"]:
                    for block in range(completed, target):
                        backends = list(protocol["backends"])
                        if block % 2:
                            backends.reverse()
                        for backend in backends:
                            if not active[backend]:
                                continue
                            order = block_order(
                                active[backend],
                                backend,
                                block,
                                protocol["order_seed"],
                            )
                            for position, arm in enumerate(order):
                                if time.perf_counter() > deadline:
                                    raise TimeoutError(
                                        "bakeoff window expired"
                                    )
                                capture = launch(
                                    args, arm, backend, block, warmup, deadline
                                )
                                capture["position"] = position
                                report["captures"].append(capture)
                                save(output, report)
                        print(args.phase, "block", block + 1, flush=True)
                    completed = target
                    report["summary"] = [
                        compare(report["captures"], backend, arm, protocol)
                        for backend in protocol["backends"]
                        for arm in arms[backend]
                        if arm != "baseline"
                    ]
                    report["rounds"].append(
                        dict(blocks=target, summary=report["summary"])
                    )
                    save(output, report)
                    for backend in protocol["backends"]:
                        unstable = [
                            row["arm"]
                            for row in report["summary"]
                            if row["backend"] == backend and not row["stable"]
                        ]
                        active[backend] = (
                            ["baseline", *unstable] if unstable else []
                        )
                    # Extend instability, never a stable interval merely
                    # because it includes zero. Otherwise repeated looks
                    # could turn a noise hunt into a reported positive result.
                    if not any(active.values()):
                        break
                if selected is None:
                    report["selected"] = select(report["summary"])
                else:
                    report["decisions"] = [
                        dict(
                            backend=row["backend"],
                            arm=row["arm"],
                            confirmed=row["timing_evidence_passes"]
                            and next(
                                old["reduction_percent"] > 0
                                for old in selection["summary"]
                                if old["backend"] == row["backend"]
                                and old["arm"] == row["arm"]
                            ),
                            default_promoted=False,
                        )
                        for row in report["summary"]
                    ]
                verify_inputs()
                save(output, report)
        finally:
            owner.unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", action="store_true")
    parser.add_argument(
        "--phase", choices=("screen", "confirmation"), default="screen"
    )
    parser.add_argument("--arm", choices=kernels.ARMS, default="baseline")
    parser.add_argument(
        "--backend",
        choices=("python-int", "gmpy2-mpz"),
        default="python-int",
    )
    parser.add_argument("--block", type=int, default=0)
    parser.add_argument("--warmup", type=float, default=3)
    parser.add_argument("--selection", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    common.require_runtime()
    if args.worker:
        print(json.dumps(worker(args)))
    else:
        if args.output is None or (
            args.phase == "confirmation" and args.selection is None
        ):
            parser.error("output and confirmation selection are required")
        run(args)


if __name__ == "__main__":
    main()
