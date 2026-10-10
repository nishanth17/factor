"""Frozen C1 complete-factor training and fresh paired PyPy confirmation."""

import argparse
import hashlib
import importlib
import json
import os
import platform
import random
import statistics
import sys
import tempfile
import time
import types
from dataclasses import asdict, replace
from pathlib import Path

from .. import budget as current_budget
from .. import qs, utils
from .b1_calibration import save
from .build_phase_two_corpus import certified_prime, verify_certificates
from .c1_confirmation_support import cold_start, summarize
from .c1_feasibility import machine_window, validate_control
from .c1_followup import HERE, ROOT, load_followup
from .p52_realistic import check_quiet
from .performance_audit import checked_baseline, fingerprint, verify_corpus
from .phase_three_reference import _rss_bytes

BASELINE = HERE / "inputs/baselines/c1_slp_baseline.json"
PROTOCOL = HERE / "c1_implementation_protocol.md"
ARMS = ("slp", "graph_slp", "dlp", "dlp_narrow", "dlp_half")
CAPS = {30: 5, 40: 30, 50: 120}
SEEDS = (7, 29, 47)


def source_hashes():
    hashes = fingerprint(ROOT)
    for path in (
        Path(__file__),
        HERE / "c1_confirmation_support.py",
        HERE / "c1_implementation_controls.md",
        BASELINE,
        PROTOCOL,
    ):
        hashes[str(path.relative_to(ROOT))] = hashlib.sha256(
            path.read_bytes()
        ).hexdigest()
    return hashes


def load_slp(directory):
    """Use the committed pre-DLP implementation, with its own budget type."""
    for name, source in checked_baseline(BASELINE)["source"].items():
        path = Path(directory) / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(source)
    package = types.ModuleType("_c1_slp_control")
    package.__path__ = [str(Path(directory) / "v2")]
    package.__package__ = package.__name__
    sys.modules[package.__name__] = package
    return tuple(
        importlib.import_module(package.__name__ + "." + name)
        for name in ("qs", "budget", "utils")
    )


def arm_config(config, arm):
    if arm == "slp":
        return config
    bound = config.base_bound // 2 if arm == "dlp_half" else config.base_bound
    options = asdict(config.collector)
    options["residual_bound"] = bound**2
    options.update(
        large_prime_bound=2 if arm == "graph_slp" else 100 * bound,
        large_product_bound=4 if arm == "graph_slp" else 128 * bound**2,
        candidate_bound=bound**2 if arm == "dlp_narrow" else 0,
        split_call_limit=131072,
    )
    return replace(
        config,
        base_bound=bound,
        collector=qs.DoubleLargeSieveConfig(**options),
    )


def decode_config(value, module=qs):
    values = dict(value)
    name = (
        "DoubleLargeSieveConfig"
        if "large_prime_bound" in values["collector"]
        else "SieveConfig"
    )
    values["collector"] = getattr(module, name)(**values["collector"])
    return module.SIQSConfig(**values)


def run_one(fixture, seed, config, arm, seconds, baseline=None, *, cold=False):
    """Time all computation and final result validation, including failures."""
    allowed = {os.getpid(), os.getppid()} if cold else {os.getpid()}
    check_quiet(allowed)
    module, ledger, numbers = (
        baseline if arm == "slp" else (qs, current_budget, utils)
    )
    config = decode_config(asdict(config), module)
    last_poll = time.monotonic()

    def poll():
        nonlocal last_poll
        if time.monotonic() - last_poll >= 1:
            check_quiet(allowed)
            last_poll = time.monotonic()
        return False

    started, cpu = time.perf_counter(), time.process_time()
    budget = ledger.Budget(
        work_limit=10**13, seconds=seconds, cpu_seconds=seconds, cancelled=poll
    )
    job = module.SIQSJob(fixture["n"], seed=seed, config=config, budget=budget)
    result = job.run()
    if (result.divisor or 1) * result.cofactor != fixture["n"]:
        raise AssertionError("SIQS result does not reconstruct")
    factors, labels, remaining = [], [], [fixture["n"]]
    reason = result.reason
    if result.divisor is not None:
        remaining = [result.divisor, result.cofactor]
        try:
            for child in tuple(remaining):
                budget.consume(child.bit_length() ** 2)
                label = numbers.classify_prime(child)
                if label is numbers.Primality.COMPOSITE:
                    break
                factors.append(child)
                labels.append(label.value)
                remaining.pop(0)
            budget.consume(0)
        except ledger.BudgetExhaustedError:
            reason = "classification_" + budget.reason
    row = dict(
        id=fixture["id"],
        digits=fixture["digits"],
        seed=seed,
        arm=arm,
        complete=not remaining,
        divisor=result.divisor,
        factors=factors,
        remaining=remaining,
        certainty=labels,
        reason=reason,
        work=budget.used,
        cap_seconds=seconds,
        stats=result.stats,
    )
    validate_control(row, fixture)
    row["seconds"] = time.perf_counter() - started
    row["cpu_seconds"] = time.process_time() - cpu
    row["rss_bytes"] = _rss_bytes()
    row["capped_seconds"] = row["seconds"] if row["complete"] else seconds
    check_quiet(allowed)
    return row


def campaign_limit(start, cpu, seconds, reserve=0):
    if (
        max(time.monotonic() - start, time.process_time() - cpu) + reserve
        >= seconds
    ):
        raise RuntimeError("frozen aggregate campaign allowance reached")


def warmup(fixture, arm, config, baseline, start, cpu, cap):
    rows, elapsed = [], 0
    while elapsed < 3:
        campaign_limit(start, cpu, cap, CAPS[fixture["digits"]] + 1)
        row = run_one(
            fixture, 7, config, arm, CAPS[fixture["digits"]], baseline
        )
        rows.append(row)
        elapsed += row["seconds"]
    return rows


def train(output, baseline):
    _, fixtures, configs = load_followup()
    fixtures = [f for f in fixtures if f["digits"] in (40, 50)]
    metadata = dict(
        schema=1,
        phase="training",
        source=source_hashes(),
        runtime=platform.python_version(),
        configurations={
            str(d): {a: asdict(arm_config(configs[d], a)) for a in ARMS}
            for d in (40, 50)
        },
    )
    save(output / "manifest.json", metadata)
    start, cpu = time.monotonic(), time.process_time()
    for digits in (40, 50):
        cases = [f for f in fixtures if f["digits"] == digits]
        for arm in ARMS:
            config = arm_config(configs[digits], arm)
            save(
                output / f"warmup-{digits}-{arm}.json",
                warmup(cases[0], arm, config, baseline, start, cpu, 4500),
            )
        for i, fixture in enumerate(cases):
            order = list(ARMS)
            random.Random(731 + digits + i).shuffle(order)
            for arm in order:
                campaign_limit(start, cpu, 4500, CAPS[digits] + 1)
                row = run_one(
                    fixture,
                    SEEDS[i],
                    arm_config(configs[digits], arm),
                    arm,
                    CAPS[digits],
                    baseline,
                )
                save(output / f"run-{digits}-{i}-{arm}.json", row)
                print(
                    json.dumps(
                        {
                            k: row[k]
                            for k in (
                                "id",
                                "arm",
                                "complete",
                                "reason",
                                "seconds",
                            )
                        }
                    ),
                    flush=True,
                )
    save(
        output / "finished.json",
        dict(
            seconds=time.monotonic() - start,
            cpu_seconds=time.process_time() - cpu,
        ),
    )


class GenerationRandom(random.Random):
    """Bound certificate generation by random calls and elapsed time."""

    def __init__(self):
        super().__init__(202610100138)
        self.calls = 0
        self.started = time.monotonic()

    def randrange(self, *args, **kwargs):
        self.calls += 1
        if self.calls > 1000000 or time.monotonic() - self.started > 120:
            raise RuntimeError("fresh corpus generation allowance exhausted")
        return super().randrange(*args, **kwargs)


def fresh_corpus():
    generator, certificates, fixtures, seen = GenerationRandom(), {}, [], set()
    for digits in (30, 40, 50):
        bits = (digits * 3322 // 1000 + 1) // 2
        for index in range(2):
            for _ in range(10000):
                p = certified_prime(bits, generator, certificates)
                q = certified_prime(bits, generator, certificates)
                if p != q and p * q not in seen and len(str(p * q)) == digits:
                    break
            else:
                raise RuntimeError("fresh product generation cap")
            seen.add(p * q)
            fixtures.append(
                dict(
                    id=f"c1_fresh_{digits}_{index}",
                    n=p * q,
                    factors=sorted((p, q)),
                    kind="balanced",
                    digits=digits,
                    split="held_out",
                )
            )
    verify_certificates(certificates)
    corpus = dict(
        schema=1,
        fixtures=fixtures,
        certificates=certificates,
        seed=202610100138,
        sampling=(
            "Pocklington-certified balanced primes with a large p-1 "
            "factor; not an RSA sample."
        ),
        generation_calls=generator.calls,
    )
    verify_corpus(corpus)
    return corpus


def freeze(training, selected_path, corpus_path):
    manifest = json.loads((training / "manifest.json").read_text())
    if manifest["source"] != source_hashes():
        raise ValueError("source changed since training")
    if not (training / "finished.json").exists():
        raise ValueError("incomplete training cannot select a policy")
    _, _, configs = load_followup()
    choices, selected = {}, {}
    for digits in (40, 50):

        def rank(arm):
            rows = [
                json.loads(
                    (training / f"run-{digits}-{i}-{arm}.json").read_text()
                )
                for i in range(2)
            ]
            return (
                -sum(r["complete"] for r in rows),
                sum(r["capped_seconds"] for r in rows),
            )

        winner = min(ARMS[2:], key=rank)
        arms = ["slp", winner]
        if rank("graph_slp") < min(rank("slp"), rank(winner)):
            arms.append("graph_slp")
        choices[str(digits)] = arms
        selected[str(digits)] = {
            a: asdict(arm_config(configs[digits], a)) for a in arms
        }
    choices["30"] = choices["40"]
    selected["30"] = {
        a: asdict(arm_config(configs[30], a)) for a in choices["30"]
    }
    # Exclusive creation precedes fresh integer generation: no fresh tuning.
    save(
        selected_path,
        dict(
            schema=1,
            source=source_hashes(),
            choices=choices,
            configurations=selected,
            training_manifest=hashlib.sha256(
                (training / "manifest.json").read_bytes()
            ).hexdigest(),
        ),
    )
    save(corpus_path, fresh_corpus())


def unstable(rows):
    values = [r["seconds"] for r in rows]
    median = statistics.median(values)
    quartiles = statistics.quantiles(values, n=4, method="inclusive")
    third = len(values) // 3
    return (quartiles[2] - quartiles[0]) / median > 0.10 or abs(
        statistics.median(values[:third]) - statistics.median(values[-third:])
    ) / median > 0.10


def confirm(selected_path, corpus_path, output, baseline):
    selection = json.loads(selected_path.read_text())
    if selection["source"] != source_hashes():
        raise ValueError("frozen source changed before confirmation")
    corpus = json.loads(corpus_path.read_text())
    verify_corpus(corpus)
    save(
        output / "manifest.json",
        dict(
            source=selection["source"],
            selected_sha256=hashlib.sha256(
                selected_path.read_bytes()
            ).hexdigest(),
            corpus_sha256=hashlib.sha256(corpus_path.read_bytes()).hexdigest(),
            runtime=platform.python_version(),
        ),
    )
    start, cpu = time.monotonic(), time.process_time()
    for fixture in corpus["fixtures"]:
        digits, identity = fixture["digits"], fixture["id"]
        arms = selection["choices"][str(digits)]
        configs = {
            a: decode_config(selection["configurations"][str(digits)][a])
            for a in arms
        }
        samples = {a: [] for a in arms}
        if identity.endswith("_0"):
            for arm in arms:
                campaign_limit(start, cpu, 10800, CAPS[digits] + 31)
                cold_start(
                    fixture,
                    arm,
                    selected_path,
                    corpus_path,
                    output,
                    CAPS[digits],
                )
        for arm in arms:
            save(
                output / f"{identity}-warmup-{arm}.json",
                warmup(
                    fixture, arm, configs[arm], baseline, start, cpu, 10800
                ),
            )
        for target in (9, 18, 27):
            for index in range(len(samples[arms[0]]), target):
                order = list(arms)
                random.Random(
                    380410 + digits * 100 + index + int(identity[-1])
                ).shuffle(order)
                for arm in order:
                    campaign_limit(start, cpu, 10800, CAPS[digits] + 1)
                    row = run_one(
                        fixture,
                        SEEDS[index % 3],
                        configs[arm],
                        arm,
                        CAPS[digits],
                        baseline,
                    )
                    row["sample"] = index
                    samples[arm].append(row)
                    save(output / f"{identity}-{index:02d}-{arm}.json", row)
                    print(
                        json.dumps(
                            {
                                k: row[k]
                                for k in (
                                    "id",
                                    "arm",
                                    "sample",
                                    "complete",
                                    "reason",
                                    "seconds",
                                )
                            }
                        ),
                        flush=True,
                    )
            if not any(unstable(rows) for rows in samples.values()):
                break
        save(
            output / f"{identity}-stability.json",
            {
                a: dict(count=len(rows), unstable=unstable(rows))
                for a, rows in samples.items()
            },
        )
    save(
        output / "finished.json",
        dict(
            seconds=time.monotonic() - start,
            cpu_seconds=time.process_time() - cpu,
        ),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "mode", choices=("train", "freeze", "confirm", "cold-worker")
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--training", type=Path)
    parser.add_argument("--fixture")
    parser.add_argument("--arm", choices=ARMS)
    parser.add_argument(
        "--selected",
        type=Path,
        default=HERE / "inputs/controls/c1_selected.json",
    )
    parser.add_argument(
        "--corpus",
        type=Path,
        default=HERE / "inputs/corpora/c1_fresh_corpus.json",
    )
    args = parser.parse_args()
    if platform.python_implementation() != "PyPy" or sys.version_info[:2] != (
        3,
        11,
    ):
        raise RuntimeError("C1 timings require PyPy implementing Python 3.11")
    if args.mode == "cold-worker":
        selection = json.loads(args.selected.read_text())
        if selection["source"] != source_hashes():
            raise ValueError("cold worker source mismatch")
        corpus = json.loads(args.corpus.read_text())
        verify_corpus(corpus)
        fixture = next(
            f for f in corpus["fixtures"] if f["id"] == args.fixture
        )
        config = decode_config(
            selection["configurations"][str(fixture["digits"])][args.arm]
        )
        with tempfile.TemporaryDirectory(prefix="c1-cold-") as directory:
            row = run_one(
                fixture,
                7,
                config,
                args.arm,
                CAPS[fixture["digits"]],
                load_slp(directory),
                cold=True,
            )
            save(args.output, row)
        return
    with machine_window():
        if args.mode == "freeze":
            freeze(args.training, args.selected, args.corpus)
            return
        args.output.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix="c1-slp-") as directory:
            baseline = load_slp(directory)
            if args.mode == "train":
                train(args.output, baseline)
            else:
                confirm(args.selected, args.corpus, args.output, baseline)
                save(args.output / "summary.json", summarize(args.output))


if __name__ == "__main__":
    main()
