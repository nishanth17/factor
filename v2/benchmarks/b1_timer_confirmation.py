"""Extend B1 repeats using factoring timers, excluding process-audit costs."""

import argparse
import json
import sys
import time
from pathlib import Path

from . import b1_calibration as b1
from .p38_r1_capacity import decode_config
from .p38_r3 import comparison, stable


def run(selected_path, corpus_path, kind, output):
    selected = json.loads(selected_path.read_text())
    source = b1.hashes()
    driver_hash = b1.digest(Path(__file__))
    if selected["source_sha256"] != source:
        raise ValueError("selection/source freeze mismatch")
    corpus = b1.load_corpus(corpus_path)
    if corpus["frozen_sha256"] != b1.digest(selected_path):
        raise ValueError("corpus/selection freeze mismatch")
    fixtures = [
        f
        for f in corpus["fixtures"]
        if f["kind"] == kind and f["digits"] == 30
    ]
    if not fixtures:
        raise ValueError("empty confirmation class")
    names = [
        "r2_control",
        selected["selected"]["siqs"],
        selected["selected"]["mpqs"],
    ]
    configs = {
        name: decode_config(selected["configurations"][name]) for name in names
    }
    capture = dict(
        source_sha256=source,
        driver_sha256=driver_hash,
        kind=kind,
        selected_sha256=b1.digest(selected_path),
        corpus_sha256=b1.digest(corpus_path),
        seeds=b1.SEEDS,
        configurations={
            name: selected["configurations"][name] for name in names
        },
        measurement="Factoring warmup >=3s; 9 then 15 repeats; <=3 blocks",
        cap_seconds=5,
        work_limit=b1.WORK,
        memory_limit=b1.MEMORY,
        attempts=[],
    )
    b1.save(output.with_suffix(".frozen.json"), capture)

    def call(name):
        return [
            b1.run_one(f, seed, configs[name], 5)
            for f in fixtures
            for seed in b1.SEEDS
        ]

    with b1.performance_window():
        for block in range(3):
            warmups = {}
            for name in names:
                started = time.perf_counter()
                seconds, count = 0, 0
                while seconds < (3 if block == 0 else 5):
                    rows = call(name)
                    seconds += sum(r["seconds"] for r in rows)
                    count += 1
                warmups[name] = dict(
                    seconds=seconds,
                    validated_calls=count,
                    wall_seconds=time.perf_counter() - started,
                )
            samples = {name: [] for name in names}
            for index in range(9 if block == 0 else 15):
                for name in names if index % 2 == 0 else list(reversed(names)):
                    rows = call(name)
                    samples[name].append(
                        dict(
                            seconds=sum(r["seconds"] for r in rows),
                            result=rows,
                        )
                    )
            arms = {name: stable(values) for name, values in samples.items()}
            capture["attempts"].append(dict(warmups=warmups, arms=arms))
            b1.save(
                output.with_name(output.stem + f"-block{block}.json"),
                capture["attempts"][-1],
            )
            print(
                "factor-timer stability",
                {k: v["stable"] for k, v in arms.items()},
                flush=True,
            )
            if all(arm["stable"] for arm in arms.values()):
                break
    capture["stable"] = all(arm["stable"] for arm in arms.values())
    capture["comparisons"] = {
        name: comparison(capture, "r2_control", name)
        for name in names
        if name != "r2_control"
    }
    if source != b1.hashes() or driver_hash != b1.digest(Path(__file__)):
        raise RuntimeError("source changed during timer confirmation")
    b1.save(output, capture)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selected", type=Path, required=True)
    parser.add_argument("--corpus", type=Path, required=True)
    parser.add_argument("--kind", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if sys.implementation.name != "pypy" or sys.version_info[:2] != (3, 11):
        parser.error("PyPy implementing Python 3.11 is required")
    run(args.selected, args.corpus, args.kind, args.output)


if __name__ == "__main__":
    main()
