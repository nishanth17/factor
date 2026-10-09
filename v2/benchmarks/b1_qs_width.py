"""Frozen wider-QS follow-up to B1's finite-window exhaustion controls."""

import argparse
import json
import sys
from dataclasses import asdict, replace
from pathlib import Path

from . import b1_calibration as b1
from .build_p38_r1_corpus import build
from .p38_r1_capacity import decode_config
from .p38_r3 import comparison, paired_measure


def configurations():
    result = {}
    for bound, width in (
        (3000, 499999),
        (10000, 131072),
        (10000, 499999),
        (30000, 499999),
    ):
        config = b1.common_config(bound, width, 1, rows=8192, partials=8192)
        result[f"qs_{bound}_{width}"] = replace(
            config,
            mode="qs",
            assignment_policy="reference",
            family_count=1,
            polynomials_per_family=0,
        )
    return result


def freeze(path, parent):
    b1.checked_protocol(parent)
    b1.save(
        path,
        dict(
            schema=1,
            parent_sha256=b1.digest(parent),
            source_sha256=b1.hashes(),
            driver_sha256=b1.digest(Path(__file__)),
            configurations={k: asdict(v) for k, v in configurations().items()},
            training_sha256=b1.digest(b1.TRAINING),
            r2_training_sha256=b1.digest(b1.R2_TRAINING),
            seeds=b1.SEEDS,
            cap_seconds=15,
            work_limit=b1.WORK,
            owned_memory=b1.MEMORY,
            cohort="Three known balanced 30-digit inputs, two seeds",
            timing_gate="All six starts complete <=5s; other arms diagnostic",
            selection="Completion; complete median or verified rows",
            measurement="Validated >=3s warmup; 9 samples; <=3 blocks",
            fresh_seed=202610093038,
            confirmation="One new balanced 30-digit input after selection",
            rationale="Fixed QS windows ended before resource allowances",
            decision="No global default/dispatch promotion",
        ),
    )


def run(protocol, selected_path, corpus_path, output, siqs_selected):
    data = json.loads(protocol.read_text())
    if (
        data["source_sha256"] != b1.hashes()
        or data["driver_sha256"] != b1.digest(Path(__file__))
        or data["training_sha256"] != b1.digest(b1.TRAINING)
        or data["r2_training_sha256"] != b1.digest(b1.R2_TRAINING)
    ):
        raise ValueError("wider-QS source/corpus freeze changed")
    fixtures = [
        f
        for f in b1.load_training()["fixtures"]
        if f["kind"] == "balanced" and f["digits"] == 30
    ]
    configs = {k: decode_config(v) for k, v in data["configurations"].items()}
    capture = dict(protocol_sha256=b1.digest(protocol), training={})
    with b1.performance_window():
        arms, costs = {}, {}
        for name, config in configs.items():

            def call():
                return [
                    b1.run_one(f, seed, config, 15)
                    for f in fixtures
                    for seed in b1.SEEDS
                ]

            screen = call()
            measurement = None
            arms[name] = screen
            costs[name] = sum(r["seconds"] for r in screen)
            if all(r["complete"] and r["seconds"] <= 5 for r in screen):
                measurement = paired_measure({name: call})
                final = measurement["attempts"][-1]["arms"][name]
                costs[name] = final["median_seconds"]
                arms[name] = [r for s in final["samples"] for r in s["result"]]
            capture["training"][name] = dict(
                screen=screen, measurement=measurement
            )
            b1.save(
                output.with_name(output.stem + "-train-" + name + ".json"),
                capture["training"][name],
            )
            print("QS width", name, [r["reason"] for r in screen], flush=True)
        chosen = b1.rank(arms, costs)[0]
        previous = json.loads(siqs_selected.read_text())
        if previous["source_sha256"] != b1.hashes():
            raise ValueError("SIQS selection source mismatch")
        siqs_name = previous["selected"]["siqs"]
        configs.update(
            {
                k: decode_config(previous["configurations"][k])
                for k in ("r2_control", siqs_name)
            }
        )
        names = ["r2_control", siqs_name, chosen]
        b1.save(
            selected_path,
            dict(
                protocol_sha256=b1.digest(protocol),
                selected=chosen,
                source_sha256=b1.hashes(),
                driver_sha256=b1.digest(Path(__file__)),
                siqs_selected_sha256=b1.digest(siqs_selected),
                configurations={k: asdict(configs[k]) for k in names},
            ),
        )
        corpus = build(
            data["fresh_seed"],
            split="b1_qs_held_out",
            count=1,
            frozen_sha256=b1.digest(selected_path),
        )
        b1.save(corpus_path, corpus)
        fixture = next(
            f
            for f in corpus["fixtures"]
            if f["kind"] == "balanced" and f["digits"] == 30
        )
        calls = {
            name: (
                lambda config=configs[name]: [
                    b1.run_one(fixture, seed, config, 15) for seed in b1.SEEDS
                ]
            )
            for name in names
        }
        capture["screen"] = {name: call() for name, call in calls.items()}
        calls = {
            name: call
            for name, call in calls.items()
            if all(
                r["complete"] and r["seconds"] <= 5
                for r in capture["screen"][name]
            )
        }
        if "r2_control" in calls and len(calls) > 1:
            measurement = paired_measure(calls)
            capture["confirmation"] = measurement
            capture["comparisons"] = {
                name: comparison(measurement, "r2_control", name)
                for name in calls
                if name != "r2_control"
            }
        capture.update(
            selected_sha256=b1.digest(selected_path),
            corpus_sha256=b1.digest(corpus_path),
        )
    if data["source_sha256"] != b1.hashes():
        raise RuntimeError("source changed during wider-QS capture")
    b1.save(output, capture)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("phase", choices=("freeze", "run"))
    for name in (
        "protocol",
        "parent",
        "selected",
        "corpus",
        "output",
        "siqs-selected",
    ):
        parser.add_argument(
            "--" + name, type=Path, required=name == "protocol"
        )
    args = parser.parse_args()
    if sys.implementation.name != "pypy" or sys.version_info[:2] != (3, 11):
        parser.error("PyPy implementing Python 3.11 is required")
    if args.phase == "freeze":
        freeze(args.protocol, args.parent)
    else:
        run(
            args.protocol,
            args.selected,
            args.corpus,
            args.output,
            args.siqs_selected,
        )


if __name__ == "__main__":
    main()
