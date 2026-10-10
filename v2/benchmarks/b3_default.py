"""Frozen user-directed default bridge; historical C6/B3 evidence is intact."""

import argparse
import hashlib
import json
import resource
import signal
import subprocess
import sys
import time
from pathlib import Path

from .. import portfolio
from . import b3_coverage
from . import b3_production as study
from .a6_pm1 import performance_window, require_runtime
from .b4_common import FrozenFinder
from .build_phase_two_corpus import verify_certificates

ROOT = Path(__file__).resolve().parents[2]
INPUTS = Path(__file__).parent / "inputs"
PROTOCOL = INPUTS / "controls/b3_default_protocol.json"
BASELINE = INPUTS / "baselines/b3_pre_default.json"
CORPORA = {
    "training": INPUTS / "corpora/b3_corpus.json",
    "confirmation": INPUTS / "corpora/b3_coverage_confirmation.json",
}


def freeze():
    """Save current committed mainline and bridge controls before timing."""
    if PROTOCOL.exists() or BASELINE.exists():
        raise ValueError("default bridge is already frozen")
    ref = subprocess.check_output(
        ["git", "rev-parse", "master"], text=True
    ).strip()
    names = subprocess.check_output(
        ["git", "ls-tree", "-r", "--name-only", ref, "v2"], text=True
    ).splitlines()
    sources = {
        name: subprocess.check_output(
            ["git", "show", ref + ":" + name]
        ).decode()
        for name in names
        if name.endswith(".py")
        and not name.startswith(("v2/benchmarks/", "v2/tests/"))
    }
    BASELINE.write_text(
        json.dumps(
            dict(
                commit=ref,
                source=sources,
                sha256={
                    name: hashlib.sha256(value.encode()).hexdigest()
                    for name, value in sources.items()
                },
            ),
            indent=2,
        )
        + "\n"
    )
    pins = [ROOT / name for name in sources]
    pins += [
        ROOT / "v2/ecm_chains.py",
        ROOT / "v2/ecm_chain_records.py",
        ROOT / "v2/ecm_chain_options.py",
        Path(__file__),
        Path(study.__file__),
        Path(b3_coverage.__file__),
        BASELINE,
        INPUTS / "controls/c6_fast_records.json",
        INPUTS / "controls/b3_cf_records.json",
        *CORPORA.values(),
    ]
    PROTOCOL.write_text(
        json.dumps(
            dict(
                schema=1,
                control_commit=ref,
                source_commit=subprocess.check_output(
                    ["git", "rev-parse", "HEAD"], text=True
                ).strip(),
                sha256={
                    str(path.relative_to(ROOT)): study.digest(path)
                    for path in pins
                },
                training_seeds=[31001, 38920],
                confirmation_seeds=[62677, 70596],
                b1=2000,
                b2=147396,
                curves=32,
                work=32_000_000,
                seconds=20,
                cpu_seconds=20,
                memory_bytes=32 * 1024**2,
                program_bytes=512 * 1024,
                chain_bytes=8 * 1024**2,
                sampling=[[3, 9], [5, 18], [8, 27]],
                max_relative_iqr=0.15,
                bootstrap_seed=193001,
                arms=["ladder", "default", "lucas", "cf"],
                scopes=["portfolio", "resume_equal", "default_budget"],
                policy=(
                    "User-directed default; no new selection or population "
                    "confirmation. Previously used certified inputs supply "
                    "a regression bridge. Setup, proof, recovery, failed "
                    "curves and validation are timed. Native ints only; "
                    "optional Lucas/CF use batch16. Equal-search runs need "
                    "32 curves for each unresolved cofactor. The 2M "
                    "budget diagnostic is censored, with no equal-search "
                    "speed claim. User promotion implies no measured gain."
                ),
            ),
            indent=2,
        )
        + "\n"
    )


def inputs(split):
    settings = json.loads(PROTOCOL.read_text())
    for name, expected in settings["sha256"].items():
        if study.digest(ROOT / name) != expected:
            raise ValueError("default bridge pin changed: " + name)
    corpus = json.loads(CORPORA[split].read_text())
    verify_certificates(corpus["certificates"])
    fixtures = [
        f for f in corpus["fixtures"] if f.get("split", split) == split
    ]
    return settings, fixtures


def baseline():
    name = "_b3_pre_default"
    data = json.loads(BASELINE.read_text())
    for path, source in data["source"].items():
        if hashlib.sha256(source.encode()).hexdigest() != data["sha256"][path]:
            raise ValueError("corrupt pre-default baseline")
    if name not in sys.modules:
        sys.meta_path.insert(0, FrozenFinder(name, data["source"]))
    return __import__(name + ".portfolio", fromlist=["portfolio"])


def cohort(settings, fixtures, seeds, arm, scope):
    engine = baseline() if arm == "ladder" else portfolio
    options = dict(
        trial_bound=5,
        rho_attempts=0,
        pm1_attempts=0,
        ecm_tiers=((settings["b1"], settings["b2"], settings["curves"]),),
        chunk_size=16,
        segment_size=128,
        max_input_bits=512,
        memory_bytes=settings["memory_bytes"],
        ecm_program_bytes=settings["program_bytes"],
    )
    if arm in ("lucas", "cf"):
        options["ecm_chain_family"] = arm
    config = engine.PortfolioConfig(**options)
    rows = []
    for fixture in fixtures:
        if scope == "resume_equal" and not fixture["id"].startswith(
            "balanced"
        ):
            continue
        for seed in seeds:
            started = time.perf_counter()
            checkpoint = (
                b3_coverage.equal_prefix_pause(engine, fixture, config, seed)
                if scope == "resume_equal"
                else None
            )
            run = engine.factorize_bounded(
                fixture["n"],
                seed=seed,
                config=config,
                checkpoint=checkpoint,
                budget=engine.Budget(
                    work_limit=2_000_000
                    if scope == "default_budget"
                    else settings["work"],
                    seconds=settings["seconds"],
                    cpu_seconds=settings["cpu_seconds"],
                ),
            )
            row = study.validate(run, fixture)
            curves = {
                n: sum(
                    e.get("stage") == "ecm" and e.get("n") == n
                    for e in run.events
                )
                for n in run.result.remaining
            }
            if scope != "default_budget" and not run.result.complete:
                if run.reason != "exhausted" or any(
                    count != settings["curves"] for count in curves.values()
                ):
                    raise AssertionError("unequal unsuccessful curve coverage")
            rows.append(
                dict(
                    row,
                    seed=seed,
                    curves=curves,
                    seconds=time.perf_counter() - started,
                )
            )
    return rows


def worker(args):
    settings, fixtures = inputs(args.split)
    seeds = settings[args.split + "_seeds"]
    arms = (
        settings["arms"]
        if args.scope == "portfolio"
        else ["ladder", "default"]
    )
    signatures, samples, warmups = {}, [], []

    def measure(arm):
        started, cpu = time.perf_counter(), time.process_time()
        rows = cohort(settings, fixtures, seeds, arm, args.scope)
        observation = dict(
            wall=time.perf_counter() - started,
            cpu=time.process_time() - cpu,
            rows=rows,
        )
        signature = json.dumps(
            [{k: v for k, v in row.items() if k != "seconds"} for row in rows],
            sort_keys=True,
        )
        if arm in signatures and signatures[arm] != signature:
            raise AssertionError("default bridge outcome/work changed")
        signatures[arm] = signature
        return observation

    for seconds, count in settings["sampling"]:
        for arm in arms:
            started, rounds = time.perf_counter(), 0
            while time.perf_counter() - started < seconds:
                measure(arm)
                rounds += 1
            warmups.append(
                dict(
                    arm=arm,
                    seconds=time.perf_counter() - started,
                    cohorts=rounds,
                )
            )
        while len(samples) < count:
            order = list(range(len(arms)))
            shift = len(samples) % len(arms)
            order = order[shift:] + order[:shift]
            if len(samples) % 2:
                order.reverse()
            row = [None] * len(arms)
            for index in order:
                row[index] = measure(arms[index])
            samples.append(row)
        summaries = {
            arm: study.summarize(
                [[row[0], row[i]] for row in samples], settings
            )
            for i, arm in enumerate(arms[1:], 1)
        }
        if all(summary["stable"] for summary in summaries.values()):
            break
    return dict(
        scope=args.scope,
        split=args.split,
        backend="python-int",
        arms=arms,
        samples=samples,
        warmups=warmups,
        summaries=summaries,
        peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        protocol_sha256=study.digest(PROTOCOL),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--freeze", action="store_true")
    parser.add_argument("--split", choices=CORPORA, default="training")
    parser.add_argument(
        "--scope",
        choices=("portfolio", "resume_equal", "default_budget"),
        default="portfolio",
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    require_runtime()
    if args.freeze:
        freeze()
        return
    if args.output is None or args.output.exists():
        parser.error("choose a new output path")
    with performance_window():

        def expired(signum, frame):
            raise TimeoutError("default bridge worker exceeded 600 seconds")

        signal.signal(signal.SIGALRM, expired)
        signal.alarm(600)
        try:
            result = worker(args)
        finally:
            signal.alarm(0)
    args.output.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
