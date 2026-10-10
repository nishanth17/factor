"""Run pinned, bounded upstream CF search and independent minimum checking."""

import argparse
import hashlib
import json
import shutil
import sys
import time
from pathlib import Path

from . import c6_cf, c6_fast_study, c6_study
from .build_c6_inputs import require_runtime, run

UPSTREAM = Path(__file__).parent / "inputs/upstream/c6_dacbench"
PROTOCOL = Path(__file__).parent / "inputs/controls/c6_cf_protocol.json"

DRIVER = """import json
from math import isqrt
import dac, prune, tally
original_node, original_search = tally.node, prune.search

def bounded_node(kind):
    original_node(kind)
    if tally.nodes > 10000000:
        raise RuntimeError("CF search node cap")

def bounded_search(target, length, thresholds):
    if not 3 <= target <= 2000 or length > 18:
        raise RuntimeError("CF search depth/input cap")
    return original_search(target, length, thresholds)

tally.node, prune.search = bounded_node, bounded_search
primes = [p for p in range(2,2001)
          if all(p%d for d in range(2,1+isqrt(p)))]
records = {}
for p in primes:
    values = prune.chain(p)
    assert values[-1] == p and dac.is_dac_starting_from_1(values)
    records[str(p)] = values
print(json.dumps(dict(schema=1,bound=2000,primes=records,nodes=tally.nodes)))
"""


def generate(directory):
    if not directory.resolve().is_relative_to(c6_study.ROOT / "v2/audit"):
        raise ValueError("upstream Python scratch must be under v2/audit")
    began = time.perf_counter()
    directory.mkdir(parents=True, exist_ok=False)
    for name in ("prune", "ladder", "tally", "dac"):
        shutil.copyfile(
            UPSTREAM / (name + ".py.txt"), directory / (name + ".py")
        )
    source = directory / "prune.py"
    text = source.read_text()
    if text.count("int((n-1)/f1)") != 1:
        raise AssertionError("upstream exact-integer patch changed")
    source.write_text(text.replace("int((n-1)/f1)", "((n-1)//f1)"))
    (directory / "generate.py").write_text(DRIVER)
    seconds = run(
        [sys.executable, "-B", "generate.py"], directory, "output.json"
    )
    data = json.loads((directory / "output.json").read_text())
    started = time.perf_counter()
    nodes = 0
    for scalar, values in data["primes"].items():
        if time.perf_counter() - began > 60:
            raise RuntimeError("CF generation/verification wall cap")
        chain, _ = c6_cf.decode(values)
        if int(scalar) != chain.scalar:
            raise AssertionError("wrong upstream scalar")
        if chain.scalar > 2:
            nodes += c6_cf.verify_minimum(chain.scalar, len(values) - 1)
    return data, dict(
        generator_process_seconds=seconds,
        independent_verification_seconds=time.perf_counter() - started,
        search_nodes=data["nodes"],
        verification_nodes=nodes,
        output_sha256=hashlib.sha256(
            (directory / "output.json").read_bytes()
        ).hexdigest(),
    )


def freeze():
    settings = c6_fast_study.protocol().copy()
    paths = list(settings["sha256"])
    paths += [
        "v2/benchmarks/" + name
        for name in (
            "c6_cf.py",
            "c6_cf_study.py",
            "build_c6_cf_inputs.py",
            "inputs/controls/c6_cf_records.json",
            "inputs/controls/c6_cf_gate.json",
        )
    ]
    paths += [str(p.relative_to(c6_study.ROOT)) for p in UPSTREAM.iterdir()]
    settings["sha256"] = {
        name: hashlib.sha256((c6_study.ROOT / name).read_bytes()).hexdigest()
        for name in paths
    }
    settings["families"] = ["cf"]
    settings["modes"] = ["tuple", "three"]
    settings["batches"] = [16, 64]
    settings["cf_gate"] = json.loads(
        (PROTOCOL.parent / "c6_cf_gate.json").read_text()
    )
    with PROTOCOL.open("x") as output:
        json.dump(settings, output, indent=2)
        output.write("\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scratch", type=Path)
    parser.add_argument("--freeze", action="store_true")
    args = parser.parse_args()
    require_runtime()
    with c6_study.performance_window("bounded CF generation/verification"):
        if args.scratch:
            data, costs = generate(args.scratch.resolve())
            with c6_cf.DATA.open("x") as output:
                json.dump(data, output, separators=(",", ":"))
                output.write("\n")
            (args.scratch / "costs.json").write_text(
                json.dumps(costs, indent=2)
            )
            c6_cf.load_catalog()
            print(json.dumps(costs))
        if args.freeze:
            freeze()


if __name__ == "__main__":
    main()
