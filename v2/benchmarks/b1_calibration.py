"""Frozen B1 joint QS/MPQS/SIQS calibration on the integrated R2 control.

Run freeze, probe, train, select and confirm in that order. All evidence is
exclusive-create; confirmation generates its corpus only after selection.
"""

import argparse
import fcntl
import hashlib
import json
import os
import platform
import statistics
import subprocess
import sys
import time
from collections import Counter
from contextlib import contextmanager
from dataclasses import asdict, replace
from math import isqrt, prod
from pathlib import Path

from .. import utils
from ..budget import Budget, BudgetExhaustedError
from ..qs import SieveConfig, SIQSConfig, SIQSJob
from ..qs import siqs as siqs_module
from .build_p38_r1_corpus import BANDS, build
from .p38_r1_capacity import decode_config
from .p38_r2 import selected_config
from .p38_r3 import comparison, paired_measure
from .p52_realistic import check_quiet
from .performance_audit import fingerprint, verify_corpus
from .phase_three_reference import _rss_bytes

ROOT = Path(__file__).resolve().parents[2]
INPUTS = Path(__file__).parent / "inputs"
TRAINING = INPUTS / "corpora/p38_r1_training_corpus.json"
R2_TRAINING = INPUTS / "corpora/p38_r2_training_corpus.json"
WORK = 10**13
MEMORY = 256 * 2**20
SEEDS = (7, 29)
PROBE_SECONDS = {30: 5, 40: 30, 60: 30, 70: 15, 80: 15, 90: 15, 99: 15}
LOCK = Path("/private/tmp/factor-performance.lock")
OWNER = Path("/private/tmp/factor-performance-owner.json")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path, data):
    """Refuse to replace either a frozen input or a generated capture."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x") as stream:
        json.dump(data, stream, indent=2)
        stream.write("\n")


def hashes():
    result = fingerprint(ROOT)
    for name in (
        "b1_calibration.py",
        "build_p38_r1_corpus.py",
        "p38_r1_capacity.py",
        "p38_r2.py",
        "p38_r3.py",
        "p52_realistic.py",
        "performance_audit.py",
    ):
        path = Path(__file__).with_name(name)
        result[str(path.relative_to(ROOT))] = digest(path)
    return result


@contextmanager
def performance_window():
    """Share one machine-wide lease with B2; refuse overlap or missing ps."""
    with LOCK.open("a+") as lease:
        fcntl.flock(lease, fcntl.LOCK_EX | fcntl.LOCK_NB)
        check_quiet({os.getpid()})
        OWNER.write_text(
            json.dumps(
                {
                    "owner": "B1",
                    "pid": os.getpid(),
                    "cwd": str(ROOT),
                    "started_utc": time.strftime(
                        "%Y-%m-%dT%H:%M:%SZ", time.gmtime()
                    ),
                }
            )
        )
        try:
            yield
            check_quiet({os.getpid()})
        finally:
            OWNER.unlink(missing_ok=True)


def common_config(
    bound,
    width,
    count,
    *,
    policy="flyer",
    gray=64,
    residual=None,
    rows=2048,
    partials=2048,
    memory=MEMORY,
):
    return SIQSConfig(
        base_bound=bound,
        half_width=width,
        max_half_width=width,
        factor_count=count,
        assignment_policy=policy,
        family_count=100000,
        pool_size=max(32, 4 * count),
        polynomials_per_family=gray,
        max_stalled=100000,
        max_trivial=100000,
        row_excess=32,
        filter_row_growth=32,
        batch_width=4096,
        memory_bytes=memory,
        checkpoint_bytes=2**20,
        collector=SieveConfig(
            block_width=4096,
            score_policy="powers",
            division="bucket",
            residual_bound=residual or min(10**12, bound**2),
            max_relations=rows,
            max_partials=partials,
            max_atoms=min(65536, 2 * (rows + partials)),
            memory_bytes=memory,
        ),
    )


def training_configs():
    """Prespecified joint bundles, plus controls isolating reuse and stores."""
    control = selected_config(siqs_module, "30d", "current")
    base = common_config(3000, 32768, 4)
    configs = {
        "r2_control": control,
        "nearest_3k_wide": replace(base, assignment_policy="nearest"),
        "flyer_3k_wide": base,
        "flyer_3k_narrow": common_config(3000, 8192, 4),
        "flyer_1k": common_config(1000, 8192, 5),
        "flyer_10k_c3": common_config(10000, 8192, 3),
        "flyer_10k_c4": common_config(10000, 8192, 4),
        "flyer_10k_c5": common_config(10000, 8192, 5),
        "flyer_10k_wide": common_config(10000, 32768, 4),
        "gray_one": replace(base, polynomials_per_family=1),
        "gray_four": replace(base, polynomials_per_family=4),
        "residual_low": replace(
            base, collector=replace(base.collector, residual_bound=10**6)
        ),
        "residual_high": replace(
            base, collector=replace(base.collector, residual_bound=10**8)
        ),
        "store_small": common_config(
            3000, 32768, 4, rows=512, partials=512, memory=32 * 2**20
        ),
        "store_large": common_config(3000, 32768, 4, rows=8192, partials=8192),
        "mpqs_external": replace(
            base,
            mode="mpqs",
            assignment_policy="reference",
            external_coefficients=True,
            factor_count=1,
            polynomials_per_family=1,
        ),
        "qs_fixed": replace(
            control, mode="qs", factor_count=1, family_count=1
        ),
    }
    for name, bound, width in (
        ("3k_narrow", 3000, 8192),
        ("10k_narrow", 10000, 8192),
        ("10k_wide", 10000, 32768),
    ):
        config = common_config(bound, width, 1)
        configs["mpqs_" + name] = replace(
            config,
            mode="mpqs",
            assignment_policy="reference",
            external_coefficients=True,
            polynomials_per_family=1,
        )
        configs["qs_" + name] = replace(
            config,
            mode="qs",
            assignment_policy="reference",
            family_count=1,
            polynomials_per_family=0,
        )
    return configs


def upper_config(digits, mode, *, wide=False):
    bound = (30000 if wide else 10000) if digits == 40 else 100000
    width = 131072 if wide else 65536
    target = isqrt(2 * 10**digits) // width
    count = next(k for k in range(1, 33) if (bound // 3) ** k >= target)
    config = common_config(bound, width, count, rows=8192, partials=8192)
    if mode == "mpqs":
        return replace(
            config,
            mode=mode,
            assignment_policy="reference",
            external_coefficients=True,
            factor_count=1,
            polynomials_per_family=1,
        )
    if mode == "qs":
        return replace(
            config,
            mode=mode,
            assignment_policy="reference",
            family_count=1,
            factor_count=1,
            polynomials_per_family=0,
        )
    return config


def failure_class(reason, stats):
    """Keep the stopping resource separate from observed yield diagnostics."""
    if reason == "factor_found":
        return "complete"
    if "classification" in reason:
        return "classification"
    if reason in ("relation_limit", "atom_limit", "memory_limit"):
        engine = stats.get("engine", {})
        return "matrix" if engine.get("matrix_rows", 0) else "storage"
    if reason in ("wall_limit", "cpu_limit"):
        return "time"
    if reason == "work_limit":
        return "work"
    if reason in ("stalled_yield", "trivial_dependency_limit"):
        return "useful_yield"
    return "collection_schedule"


def validate(row, fixture):
    """Independently reconstruct resolved and unresolved parts and labels."""
    if prod(row["factors"]) * prod(row["remaining"]) != fixture["n"]:
        raise AssertionError("result lost a resolved or unresolved cofactor")
    expected = Counter(fixture["factors"])
    if Counter(row["factors"]) - expected:
        raise AssertionError("terminal factor disagrees with certified proof")
    if row["complete"] and Counter(row["factors"]) != expected:
        raise AssertionError("incorrect complete multiplicity")
    divisor = row["divisor"]
    if divisor is not None and not utils.valid_divisor(divisor, fixture["n"]):
        raise AssertionError("improper divisor")
    if len(row["certainty"]) != len(row["factors"]):
        raise AssertionError("factor labels missing")
    for factor, label in zip(row["factors"], row["certainty"]):
        wanted = (
            "proven_prime"
            if utils.deterministic_bases(factor) is not None
            else "probable_prime"
        )
        if label != wanted:
            raise AssertionError("runtime certainty was overstated or changed")
    if row["work"] > WORK or row["stats"].get("workspace_bytes", 0) > MEMORY:
        raise AssertionError("finite work/storage envelope exceeded")


def run_one(fixture, seed, config, seconds, *, job_type=SIQSJob):
    """Include setup/classification; R5 can supply the existing SSS job."""
    check_quiet({os.getpid()})
    last_poll = time.monotonic()

    def overlap_poll():
        nonlocal last_poll
        if time.monotonic() - last_poll >= 1:
            check_quiet({os.getpid()})
            last_poll = time.monotonic()
        return False

    budget = Budget(
        work_limit=WORK,
        seconds=seconds,
        cpu_seconds=seconds,
        cancelled=overlap_poll,
    )
    started, cpu = time.perf_counter(), time.process_time()
    job = job_type(fixture["n"], seed=seed, config=config, budget=budget)
    result = job.run()
    if (result.divisor or 1) * result.cofactor != fixture["n"]:
        raise AssertionError("split/cofactor reconstruction failed")
    factors, remaining, labels = [], [fixture["n"]], []
    reason = result.reason
    if result.divisor is not None:
        remaining = [result.divisor, result.cofactor]
        try:
            for child in tuple(remaining):
                budget.consume(child.bit_length() ** 2)
                label = utils.classify_prime(child)
                if label is utils.Primality.COMPOSITE:
                    break
                factors.append(child)
                labels.append(label.value)
                remaining.pop(0)
            budget.consume(0)
        except BudgetExhaustedError:
            reason = "classification_" + budget.reason

    elapsed = time.perf_counter() - started
    stats = result.stats
    row = dict(
        id=fixture["id"],
        kind=fixture["kind"],
        digits=fixture["digits"],
        seed=seed,
        complete=not remaining,
        divisor=result.divisor,
        factors=factors,
        remaining=remaining,
        certainty=labels,
        reason=reason,
        failure_class=failure_class(reason, stats),
        seconds=elapsed,
        cpu_seconds=time.process_time() - cpu,
        work=budget.used,
        cap_seconds=seconds,
        stats=stats,
        rss_bytes=_rss_bytes(),
        zero_useful_rows=not stats.get("relations", 0),
    )
    validate(row, fixture)
    check_quiet({os.getpid()})
    return row


def load_corpus(path):
    corpus = json.loads(path.read_text())
    verify_corpus(corpus)
    return corpus


def load_training():
    """Use only previously inspected inputs; add R2's two 30-digit cases."""
    corpus = load_corpus(TRAINING)
    extra = load_corpus(R2_TRAINING)
    for fixture in extra["fixtures"]:
        if fixture["band"] == "30d":
            corpus["fixtures"].append(
                dict(
                    fixture,
                    kind="balanced",
                    digits=30,
                    smaller_factor_digits=min(
                        len(str(p)) for p in fixture["factors"]
                    ),
                )
            )
    corpus["certificates"].update(extra["certificates"])
    verify_corpus(corpus)
    return corpus


def freeze(path):
    corpus = load_training()
    save(
        path,
        dict(
            schema=1,
            environment={
                "implementation": platform.python_implementation(),
                "python": platform.python_version(),
                "build": sys.version,
                "platform": platform.platform(),
                "machine": platform.machine(),
            },
            source_commit=subprocess.check_output(
                ["git", "rev-parse", "HEAD"], text=True
            ).strip(),
            source_sha256=hashes(),
            training_sha256=digest(TRAINING),
            r2_training_sha256=digest(R2_TRAINING),
            r2_integration_sha256=digest(
                INPUTS / "controls/p38_r2_integration_frozen.json"
            ),
            seeds=SEEDS,
            work_limit=WORK,
            owned_memory_limit=MEMORY,
            probe_seconds=PROBE_SECONDS,
            training_seconds=5,
            confirmation_seconds=5,
            held_out_generation_seed=202610090138,
            configurations={
                k: asdict(v) for k, v in training_configs().items()
            },
            upper_configurations={
                str(d): {
                    m: asdict(upper_config(d, m))
                    for m in ("qs", "mpqs", "siqs")
                }
                for d in BANDS
                if d != 30
            },
            additional_40d_wide=asdict(upper_config(40, "siqs", wide=True)),
            training_fixtures=[f["id"] for f in corpus["fixtures"]],
            selection={
                "order": "completion; median cost or useful rows; name",
                "affordable": "all six starts complete in <=2s each",
                "cohort": "three balanced 30-digit inputs and two seeds",
                "upper": "one start/mode/band; extra 40d wide SIQS",
                "extensions": "No post-sweep allowance expansion.",
                "training_screen": (
                    "Screen each bundle on six 5s starts; warm/sample "
                    "only if all complete <=2s/start; otherwise diagnostic"
                ),
                "promotion": (
                    "Stable fresh paired cohorts: >=10% complete median gain "
                    "or >=10pp completion; CI excludes zero; <=5pp regression"
                ),
                "uncertainty": (
                    "Repeat CI conditional on fixed inputs; tiny cohorts "
                    "cannot settle general default/dispatch promotion"
                ),
                "bounded_stop": (
                    "One joint sweep; select one configuration/mode at 30d; "
                    "defer upper bands and retain defaults"
                ),
            },
            confirmation={
                "new_inputs_per_class_band": 1,
                "seeds": SEEDS,
                "classes": [
                    "balanced",
                    "uneven_5",
                    "uneven_10",
                    "pm1_smooth",
                    "pp1_smooth",
                    "close",
                    "power",
                ],
                "bands": "30d only; larger bands remain feasibility probes",
                "measurement": (
                    "Nine samples; extend to 15 up to three blocks; >=3s "
                    "validated warmup/arm/block; held-out arms interleave"
                ),
                "cold_and_profiles": "excluded from warmed evidence",
                "screen": (
                    "Two-seed 5s screen/class; time arms <=2s/start; "
                    "other arms diagnostic; require the R2 control"
                ),
            },
            sampling=corpus["sampling"],
        ),
    )


def checked_protocol(path):
    protocol = json.loads(path.read_text())
    if protocol["source_sha256"] != hashes():
        raise ValueError("source changed after protocol freeze")
    if protocol["training_sha256"] != digest(TRAINING):
        raise ValueError("training corpus changed after freeze")
    if protocol["r2_training_sha256"] != digest(R2_TRAINING):
        raise ValueError("R2 training corpus changed after freeze")
    return protocol


def probe(protocol, output):
    corpus = load_corpus(TRAINING)
    capture = dict(protocol_sha256=digest(protocol), rows=[])
    with performance_window():
        for digits in BANDS:
            fixture = next(
                f
                for f in corpus["fixtures"]
                if f["kind"] == "balanced" and f["digits"] == digits
            )
            configs = (
                {
                    k: training_configs()[k]
                    for k in ("r2_control", "mpqs_external", "qs_fixed")
                }
                if digits == 30
                else {
                    m: upper_config(digits, m) for m in ("siqs", "mpqs", "qs")
                }
            )
            if digits == 40:
                configs["siqs_wide"] = upper_config(40, "siqs", wide=True)
            for name, config in configs.items():
                row = run_one(fixture, SEEDS[0], config, PROBE_SECONDS[digits])
                row.update(configuration=name, config=asdict(config))
                capture["rows"].append(row)
                print(
                    "probe",
                    digits,
                    name,
                    row["reason"],
                    row["stats"].get("relations", 0),
                    flush=True,
                )
                save(
                    output.with_name(output.stem + f"-{digits}-{name}.json"),
                    row,
                )
    save(output, capture)


def train(protocol, output):
    fixtures = [
        f
        for f in load_training()["fixtures"]
        if f["kind"] == "balanced" and f["digits"] == 30
    ]
    capture = dict(protocol_sha256=digest(protocol), arms={}, measurements={})
    with performance_window():
        for name, config in training_configs().items():

            def call():
                return [
                    run_one(f, seed, config, 5)
                    for f in fixtures
                    for seed in SEEDS
                ]

            screen = call()
            if not all(r["complete"] and r["seconds"] <= 2 for r in screen):
                capture["arms"][name] = screen
                print(
                    "infeasible",
                    name,
                    [r["reason"] for r in screen],
                    flush=True,
                )
                save(
                    output.with_name(output.stem + "-" + name + ".json"),
                    dict(rows=screen, measurement=None),
                )
                continue
            measurement = paired_measure({name: call})
            capture["measurements"][name] = measurement
            samples = measurement["attempts"][-1]["arms"][name]["samples"]
            rows = [row for sample in samples for row in sample["result"]]
            capture["arms"][name] = rows
            print(
                "train",
                name,
                sum(r["complete"] for r in rows),
                round(sum(r["seconds"] for r in rows), 3),
                flush=True,
            )
            save(
                output.with_name(output.stem + "-" + name + ".json"),
                dict(rows=rows, measurement=measurement),
            )
    save(output, capture)


def rank(arms, costs=None):
    """Rank completion first; censored costs never become speed ratios."""
    return sorted(
        arms,
        key=lambda name: (
            -sum(r["complete"] for r in arms[name]) / len(arms[name]),
            (
                costs[name]
                if costs
                else statistics.mean(r["seconds"] for r in arms[name])
            )
            if all(r["complete"] for r in arms[name])
            else float("inf"),
            -statistics.mean(
                r["stats"].get("relations", 0) for r in arms[name]
            ),
            name,
        ),
    )


def select(protocol_path, training_path, probe_path, output):
    protocol = checked_protocol(protocol_path)
    training = json.loads(training_path.read_text())
    probes = json.loads(probe_path.read_text())
    if any(
        c["protocol_sha256"] != digest(protocol_path)
        for c in (training, probes)
    ):
        raise ValueError("capture belongs to a different frozen protocol")
    arms = training["arms"]
    costs = {
        name: (
            training["measurements"][name]["attempts"][-1]["arms"][name][
                "median_seconds"
            ]
            if name in training["measurements"]
            else sum(r["seconds"] for r in rows)
        )
        for name, rows in arms.items()
    }
    chosen = {
        mode: rank(
            {
                k: rows
                for k, rows in arms.items()
                if protocol["configurations"][k]["mode"] == mode
                and k != "r2_control"
            },
            costs,
        )[0]
        for mode in ("siqs", "mpqs", "qs")
    }
    affordable = all(
        r["complete"] and r["seconds"] <= 2 for r in arms[chosen["siqs"]]
    )
    names = ["r2_control", *chosen.values()]
    save(
        output,
        dict(
            schema=1,
            protocol_sha256=digest(protocol_path),
            training_capture_sha256=digest(training_path),
            probe_capture_sha256=digest(probe_path),
            source_sha256=hashes(),
            selected=chosen,
            ranking=rank(arms, costs),
            configurations={k: protocol["configurations"][k] for k in names},
            confirmation_bands=[30] if affordable else [],
            decision="Freeze scoped challenger; retain production defaults.",
            upper_decision="Defer upper calibration; probes are diagnostic.",
        ),
    )


def confirm(protocol_path, selected_path, corpus_path, output):
    selected = json.loads(selected_path.read_text())
    if selected["protocol_sha256"] != digest(protocol_path):
        raise ValueError("selection/protocol identity mismatch")
    if selected["source_sha256"] != hashes():
        raise ValueError("source changed after selection freeze")
    # The generator seed and classes were prespecified. No held-out integers
    # exist until the selected configuration has an immutable checksum.
    corpus = build(
        202610090138,
        split="b1_held_out",
        count=1,
        frozen_sha256=digest(selected_path),
    )
    save(corpus_path, corpus)
    configs = {
        k: decode_config(v) for k, v in selected["configurations"].items()
    }
    capture = dict(
        selected_sha256=digest(selected_path),
        corpus_sha256=digest(corpus_path),
        classes={},
    )
    with performance_window():
        fixture = next(
            f
            for f in load_training()["fixtures"]
            if f["kind"] == "balanced" and f["digits"] == 30
        )
        checks = {}
        for name, config in configs.items():
            job = SIQSJob(
                fixture["n"],
                config=config,
                budget=Budget(work_limit=WORK, seconds=60, cpu_seconds=60),
            )
            first = job.run(max_blocks=1)
            checkpoint = job.checkpoint()
            resumed = SIQSJob.from_checkpoint(
                checkpoint,
                config=config,
                budget=Budget(work_limit=WORK, seconds=60, cpu_seconds=60),
            )
            if resumed.budget.used < first.stats["work_used"]:
                raise AssertionError("resume reset consumed work")
            if (
                resumed.budget.prior_wall
                < checkpoint["resources"]["wall_used"]
            ):
                raise AssertionError("resume reset consumed time")
            result = resumed.run(max_blocks=1)
            if (result.divisor or 1) * result.cofactor != fixture["n"]:
                raise AssertionError("resume lost unresolved cofactor")
            for bad_config, bad_budget in (
                (
                    replace(config, row_excess=config.row_excess + 1),
                    Budget(work_limit=WORK, seconds=60, cpu_seconds=60),
                ),
                (config, Budget(work_limit=0, seconds=60, cpu_seconds=60)),
            ):
                try:
                    SIQSJob.from_checkpoint(
                        checkpoint, config=bad_config, budget=bad_budget
                    )
                except ValueError:
                    pass
                else:
                    raise AssertionError("unchecked resume accepted")
            checks[name] = dict(
                reason=result.reason,
                prior_work=resumed.budget.used,
                checkpoint_bytes=len(checkpoint["blob"].encode()),
            )
        save(output.with_name(output.stem + "-resume.json"), checks)
        for kind in (
            "balanced",
            "uneven_5",
            "uneven_10",
            "pm1_smooth",
            "pp1_smooth",
            "close",
            "power",
        ):
            fixtures = [
                f
                for f in corpus["fixtures"]
                if f["kind"] == kind and f["digits"] == 30
            ]
            if not selected["confirmation_bands"]:
                break
            calls = {
                name: (
                    lambda config=config: [
                        run_one(f, seed, config, 5)
                        for f in fixtures
                        for seed in SEEDS
                    ]
                )
                for name, config in configs.items()
            }
            screen = {name: call() for name, call in calls.items()}
            calls = {
                name: call
                for name, call in calls.items()
                if all(
                    row["complete"] and row["seconds"] <= 2
                    for row in screen[name]
                )
            }
            if "r2_control" not in calls or len(calls) < 2:
                capture["classes"][kind] = dict(
                    screen=screen,
                    decision="Diagnostic feasibility only; no timing claim.",
                )
                save(
                    output.with_name(output.stem + "-" + kind + ".json"),
                    capture["classes"][kind],
                )
                continue
            measurement = paired_measure(calls)
            capture["classes"][kind] = dict(
                screen=screen,
                measurement=measurement,
                comparisons={
                    name: comparison(measurement, "r2_control", name)
                    for name in calls
                    if name != "r2_control"
                },
            )
            save(
                output.with_name(output.stem + "-" + kind + ".json"),
                capture["classes"][kind],
            )
    save(output, capture)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "phase", choices=("freeze", "probe", "train", "select", "confirm")
    )
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--training", type=Path)
    parser.add_argument("--probes", type=Path)
    parser.add_argument("--selected", type=Path)
    parser.add_argument("--corpus", type=Path)
    args = parser.parse_args()
    if sys.implementation.name != "pypy" or sys.version_info[:2] != (3, 11):
        parser.error("PyPy implementing Python 3.11 is required")
    if args.phase == "freeze":
        freeze(args.protocol)
        return
    checked_protocol(args.protocol)
    if args.output is None or args.output.exists():
        parser.error("a new --output path is required")
    if args.phase == "probe":
        probe(args.protocol, args.output)
    elif args.phase == "train":
        train(args.protocol, args.output)
    elif args.phase == "select":
        select(args.protocol, args.training, args.probes, args.output)
    else:
        confirm(args.protocol, args.selected, args.corpus, args.output)
    if hashes() != json.loads(args.protocol.read_text())["source_sha256"]:
        raise RuntimeError("source changed during capture")
    print(
        platform.python_implementation(),
        platform.python_version(),
        args.output,
    )


if __name__ == "__main__":
    main()
