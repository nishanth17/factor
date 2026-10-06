"""Finite larger-bound ECM probes without performance promotion."""

import argparse
import hashlib
import json
import os
import platform
import resource
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path

from v2 import portfolio
from v2.benchmarks import p52_a3, p52_realistic

PROBES = {
    "factor30": ("60d_balanced", (250000, 130000000, 4), 1_000_000_000),
    "factor40": ("80d_balanced", (3000000, 5700000000, 1), 50_000_000_000),
}
ARMS = ("baseline", "programs")


def worker(args):
    """Validate a cold capped search, retaining explicit partial progress."""
    corpus = p52_realistic.load_corpus()
    fixture_id, tier, work_limit = PROBES[args.probe]
    fixture = next(f for f in corpus["fixtures"] if f["id"] == fixture_id)
    engine = p52_a3.load_control()[0] if args.arm == "baseline" else portfolio
    options = dict(
        trial_bound=5,
        rho_attempts=0,
        pm1_attempts=0,
        ecm_tiers=(tier,),
        max_input_bits=329,
        memory_bytes=512 * 2**20,
    )
    if args.arm == "programs":
        options["ecm_program_bytes"] = 128 * 2**20
    config = engine.PortfolioConfig(**options)
    budget = engine.Budget(work_limit=work_limit, seconds=120, cpu_seconds=120)
    started = time.perf_counter()
    run = engine.factorize_bounded(
        fixture["n"], seed=7, config=config, budget=budget
    )
    elapsed = time.perf_counter() - started
    actual = Counter({f.value: f.exponent for f in run.result.factors})
    expected = Counter(dict(fixture["factors"]))
    if run.result.reconstruct() != fixture["n"] or actual - expected:
        raise AssertionError("wider-bound result does not reconstruct")
    if run.result.complete and actual != expected:
        raise AssertionError("wider-bound factors disagree with proofs")
    if run.dropped_events:
        raise AssertionError("wider-bound trace exceeded capacity")

    current = run.checkpoint["payload"]["state"]["current"]
    job = current.get("job") if current else None
    cursor = job.get("cursor") if job else None
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return dict(
        probe=args.probe,
        arm=args.arm,
        fixture=fixture_id,
        input_bits=fixture["n"].bit_length(),
        smaller_factor_digits=fixture["small_digits"],
        config=options,
        workspace_reserve=config.workspace_reserve,
        work_limit=work_limit,
        wall_limit=120,
        cpu_limit=120,
        complete=run.result.complete,
        reason=run.reason,
        elapsed_seconds=elapsed,
        consumed_wall=run.wall_seconds,
        consumed_cpu=run.cpu_seconds,
        work=run.work_used,
        completed_curves=len(run.events),
        phase=job["phase"] if job else None,
        next_endpoint=cursor["next"] if cursor else None,
        target_endpoint=cursor["hi"] if cursor else None,
        factors=[
            [f.value, f.exponent, f.certainty.value]
            for f in run.result.factors
        ],
        remaining=run.result.remaining,
        peak_rss_bytes=rss if sys.platform == "darwin" else 1024 * rss,
    )


def launch(probe, arm, quiet):
    """Monitor externally so the library's work ledger is untouched."""
    if quiet:
        p52_realistic.check_quiet({os.getpid()})
    process = subprocess.Popen(
        [
            sys.executable,
            "-B",
            "-m",
            "v2.benchmarks.p52_wider",
            "--worker",
            "--probe",
            probe,
            "--arm",
            arm,
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    try:
        while True:
            try:
                stdout, stderr = process.communicate(timeout=1)
                break
            except subprocess.TimeoutExpired:
                if quiet:
                    p52_realistic.check_quiet({os.getpid(), process.pid})
        if quiet:
            p52_realistic.check_quiet({os.getpid(), process.pid})
        if process.returncode:
            raise RuntimeError("wider-bound probe failed: " + stderr)
    except BaseException:
        process.terminate()
        process.communicate()
        raise
    return json.loads(stdout)


def run(args):
    rows = []
    for probe in PROBES:
        for arm in ARMS:
            row = launch(probe, arm, args.check_quiet)
            rows.append(row)
            with Path(str(args.output) + f".{probe}-{arm}.json").open(
                "x"
            ) as stream:
                json.dump(row, stream, indent=2)
                stream.write("\n")
            print(
                probe, arm, row["reason"], row["completed_curves"], flush=True
            )
    result = dict(
        schema=1,
        runtime=sys.version,
        platform=platform.platform(),
        quiet_monitor=args.check_quiet,
        corpus_sha256=hashlib.sha256(
            p52_realistic.CORPUS.read_bytes()
        ).hexdigest(),
        baseline_sha256=hashlib.sha256(
            p52_a3.BASELINE.read_bytes()
        ).hexdigest(),
        dependency_sha256=hashlib.sha256(
            p52_a3.DEPENDENCIES.read_bytes()
        ).hexdigest(),
        source_sha256=p52_a3.source_hashes("p52_realistic", "p52_wider"),
        rows=rows,
        limitations=[
            "One cold seeded probe per arm; not comparative speed evidence.",
            "Time/work exhaustion is retained, never imputed as completion.",
            "Native factor-size tiers are hypotheses for this PyPy engine.",
            "Owned workspace and process RSS are separate measures.",
        ],
    )
    with Path(args.output).open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--probe", choices=PROBES, default="factor30")
    parser.add_argument("--arm", choices=ARMS, default="programs")
    parser.add_argument("--check-quiet", action="store_true")
    parser.add_argument("--output")
    args = parser.parse_args()
    if platform.python_implementation() != "PyPy" or sys.version_info[:2] != (
        3,
        11,
    ):
        parser.error("require PyPy implementing Python 3.11")
    if not args.worker and not args.output:
        parser.error("require an exclusive output path")
    if args.worker:
        print(json.dumps(worker(args)))
    else:
        run(args)


if __name__ == "__main__":
    main()
