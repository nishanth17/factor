"""Bounded 40-digit follow-up after B1's successful feasibility probes.

Freeze this protocol before training. Freeze mode selections before creating
fresh confirmation inputs. This does not reopen the original 30-digit freeze.
"""

import argparse
import json
import platform
import sys
from dataclasses import asdict, replace
from pathlib import Path

from ...support.paths import (
    source_path,
)
from ..p38.build_p38_r1_corpus import build
from ..p38.p38_r1_capacity import decode_config
from ..p38.p38_r3 import comparison, paired_measure
from . import b1_calibration as b1


def configurations():
    configs = {
        mode + suffix: b1.upper_config(40, mode, wide=wide)
        for mode in ("siqs", "mpqs", "qs")
        for suffix, wide in (("_10k", False), ("_30k", True))
    }
    configs["siqs_10k_c4"] = replace(configs["siqs_10k"], factor_count=4)
    configs["siqs_10k_gray1"] = replace(
        configs["siqs_10k"], polynomials_per_family=1
    )
    return configs


def freeze(protocol, parent):
    b1.checked_protocol(parent)
    b1.save(
        protocol,
        dict(
            schema=1,
            parent_protocol_sha256=b1.digest(parent),
            source_sha256=b1.hashes(),
            driver_sha256=b1.digest(Path(__file__)),
            training_sha256=b1.digest(b1.TRAINING),
            seeds=b1.SEEDS,
            work_limit=b1.WORK,
            memory_limit=b1.MEMORY,
            cap_seconds=30,
            source_control="94caf40 integrated R2; unchanged arithmetic",
            configurations={k: asdict(v) for k, v in configurations().items()},
            training="One inspected balanced 40-digit input, two seeds",
            timing_gate="All starts <=10s complete; otherwise diagnostic",
            measurement=(
                "Three seconds validated warmup, nine samples; extend "
                "unstable blocks to 15, up to three blocks"
            ),
            selection=(
                "Completion fraction, then complete cohort median; partials "
                "rank by useful rows; one configuration/mode"
            ),
            confirmation=(
                "One new balanced 40-digit input, both seeds; generate "
                "after selection; no structured/uneven timing at 40d"
            ),
            held_out_generation_seed=202610094038,
            default_decision=(
                "Retain defaults; no population/dispatch promotion from "
                "one held-out input"
            ),
            rationale=(
                "Original protocol capped affordable starts at 2s. The 40d "
                "probe completed at 4-6s, motivating this separate extension."
            ),
            environment=dict(python=sys.version, platform=platform.platform()),
        ),
    )


def checked(protocol):
    data = json.loads(source_path(protocol).read_text())
    if data["source_sha256"] != b1.hashes():
        raise ValueError("source changed after 40d freeze")
    if data["driver_sha256"] != b1.digest(Path(__file__)):
        raise ValueError("40d driver changed after freeze")
    if data["training_sha256"] != b1.digest(b1.TRAINING):
        raise ValueError("40d training corpus changed")
    return data


def run(protocol, selected_path, corpus_path, output):
    data = checked(protocol)
    fixture = next(
        f
        for f in b1.load_corpus(b1.TRAINING)["fixtures"]
        if f["kind"] == "balanced" and f["digits"] == 40
    )
    capture = dict(protocol_sha256=b1.digest(protocol), training={})
    configs = {k: decode_config(v) for k, v in data["configurations"].items()}
    with b1.performance_window():
        for name, config in configs.items():

            def call():
                return [
                    b1.run_one(fixture, seed, config, 30) for seed in b1.SEEDS
                ]

            screen = call()
            measurement = None
            if all(r["complete"] and r["seconds"] <= 10 for r in screen):
                measurement = paired_measure({name: call})
            capture["training"][name] = dict(
                screen=screen, measurement=measurement
            )
            b1.save(
                output.with_name(output.stem + "-train-" + name + ".json"),
                capture["training"][name],
            )
            print("40d train", name, [r["reason"] for r in screen], flush=True)

        costs, arms = {}, {}
        for name, result in capture["training"].items():
            measurement = result["measurement"]
            if measurement is None:
                arms[name] = result["screen"]
                costs[name] = sum(r["seconds"] for r in arms[name])
            else:
                final = measurement["attempts"][-1]["arms"][name]
                arms[name] = [r for s in final["samples"] for r in s["result"]]
                costs[name] = final["median_seconds"]
        selected = {
            mode: b1.rank(
                {k: v for k, v in arms.items() if configs[k].mode == mode},
                costs,
            )[0]
            for mode in ("siqs", "mpqs", "qs")
        }
        b1.save(
            selected_path,
            dict(
                schema=1,
                protocol_sha256=b1.digest(protocol),
                source_sha256=b1.hashes(),
                driver_sha256=b1.digest(Path(__file__)),
                selected=selected,
                configurations={
                    k: asdict(configs[k])
                    for k in set(selected.values()) | {"siqs_10k"}
                },
                default_decision=data["default_decision"],
            ),
        )
        corpus = build(
            data["held_out_generation_seed"],
            split="b1_40d_held_out",
            count=1,
            frozen_sha256=b1.digest(selected_path),
        )
        b1.save(corpus_path, corpus)
        fixture = next(
            f
            for f in corpus["fixtures"]
            if f["kind"] == "balanced" and f["digits"] == 40
        )
        names = list(dict.fromkeys(["siqs_10k", *selected.values()]))
        calls = {
            name: (
                lambda config=configs[name]: [
                    b1.run_one(fixture, seed, config, 30) for seed in b1.SEEDS
                ]
            )
            for name in names
        }
        screen = {name: call() for name, call in calls.items()}
        calls = {
            name: call
            for name, call in calls.items()
            if all(r["complete"] and r["seconds"] <= 10 for r in screen[name])
        }
        capture.update(
            selected_sha256=b1.digest(selected_path),
            corpus_sha256=b1.digest(corpus_path),
            screen=screen,
        )
        if "siqs_10k" in calls and len(calls) > 1:
            measurement = paired_measure(calls)
            capture["confirmation"] = measurement
            capture["comparisons"] = {
                name: comparison(measurement, "siqs_10k", name)
                for name in calls
                if name != "siqs_10k"
            }
    checked(protocol)
    b1.save(output, capture)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("phase", choices=("freeze", "run"))
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--parent", type=Path)
    parser.add_argument("--selected", type=Path)
    parser.add_argument("--corpus", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if sys.implementation.name != "pypy" or sys.version_info[:2] != (3, 11):
        parser.error("PyPy implementing Python 3.11 is required")
    if args.phase == "freeze":
        freeze(args.protocol, args.parent)
    else:
        run(args.protocol, args.selected, args.corpus, args.output)


if __name__ == "__main__":
    main()
