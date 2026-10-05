"""Run hash-checked, unchanged upstream SSS functions on supported PyPy."""

import argparse
import contextlib
import hashlib
import importlib.util
import io
import json
import platform
import random
import re
import resource
import statistics
import subprocess
import sys
import time
from math import isfinite
from pathlib import Path

from .build_phase_two_corpus import verify_certificates

HASHES = {
    "sss.py": (
        "0a63e67de313418c4e723668559eab9aeaea4eebadea45560c34aca8acf988ab"
    ),
    "sssf.py": (
        "3dca1ace4175dd1a1236f0c2983103884e71f247d427c85cb172c50b0400c1f2"
    ),
    "mstep.py": (
        "7c3888929bc3dca3c41b8ee89c1993e916ba75cebbb60c02c486fd4f8a9c0464"
    ),
}
PIN = "8dbaf6d39ab88a40380965d25ec2c363d7f27358"


def supervise(args):
    """Keep unchanged upstream code inside a wall/RSS-capped child."""
    command = [
        sys.executable,
        "-m",
        "v2.benchmarks.phase_three_sss_upstream",
        *sys.argv[1:],
        "--worker",
    ]
    started, peak = time.monotonic(), 0
    process = subprocess.Popen(
        command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    failure = None

    try:
        while process.poll() is None:
            if time.monotonic() - started >= args.wall_limit:
                failure = "upstream child exceeds the external wall limit"
                break
            reading = subprocess.run(
                ["ps", "-o", "rss=", "-p", str(process.pid)],
                capture_output=True,
                text=True,
                check=False,
            )
            if reading.returncode:
                if process.poll() is not None:
                    break
                failure = "cannot enforce the upstream RSS limit"
                break

            peak = max(peak, int(reading.stdout.strip()) * 1024)
            if peak > args.rss_limit_mib * 1024 * 1024:
                failure = "upstream child exceeds the external RSS limit"
                break
            time.sleep(0.1)
    finally:
        if process.poll() is None:
            process.kill()
        stdout, stderr = process.communicate()

    if failure or process.returncode:
        raise RuntimeError(
            failure or f"upstream child failed: {stderr[-2048:]}"
        )
    result = json.loads(args.output.read_text())
    result["external_monitor"] = dict(
        wall_seconds=time.monotonic() - started,
        wall_limit=args.wall_limit,
        rss_limit_bytes=args.rss_limit_mib * 1024 * 1024,
        observed_peak_rss_bytes=peak,
        exit_code=process.returncode,
    )
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(stdout, end="")


class OutputCheck(io.TextIOBase):
    """Consume progress and check splits without changing the algorithms."""

    def __init__(self, n):
        self.n, self.tail, self.divisor, self.relations = n, "", None, None

    def write(self, text):
        self.tail = (self.tail + text)[-2048:]
        matches = re.findall(
            r"Proper factors found: (\d+) \| (\d+)", self.tail
        )

        if matches:
            left, right = map(int, matches[-1])

            if not 1 < left < self.n or left * right != self.n:
                raise AssertionError("upstream output is not a proper split")
            self.divisor = left

        matches = re.findall(r"SSS finished: (\d+) relations found", self.tail)
        if matches:
            self.relations = int(matches[-1])
        return len(text)


def run_one(module, fixture, seed, filtered):
    random.seed(seed)
    output = OutputCheck(fixture["n"])
    started, cpu = time.perf_counter(), time.process_time()
    reason = "no_factor"
    with contextlib.redirect_stdout(output):
        try:
            if filtered:
                # Explicit small-input settings; upstream demo is (10, 5).
                module.SSS(fixture["n"], 2, 0)
            else:
                module.SSS(fixture["n"])
        except SystemExit as error:
            reason = "upstream_exit"
            output.write(str(error))

    factors = (
        sorted((output.divisor, fixture["n"] // output.divisor))
        if (output.divisor is not None)
        else []
    )
    if factors and factors != fixture["factors"]:
        raise AssertionError(
            "upstream factors differ from independent certificates"
        )
    return dict(
        id=fixture["id"],
        seed=seed,
        seconds=time.perf_counter() - started,
        cpu_seconds=time.process_time() - cpu,
        completed=bool(factors),
        factors=factors,
        remaining=[] if factors else [fixture["n"]],
        relations=output.relations,
        reason="factor_found" if factors else reason,
    )


def main():
    if platform.python_implementation() != "PyPy" or sys.version_info[:2] != (
        3,
        11,
    ):
        raise RuntimeError("upstream comparison requires PyPy Python 3.11")

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--mode", choices=("sss", "sssf"), default="sss")
    parser.add_argument("--repetitions", type=int, default=9)
    parser.add_argument("--warmup-seconds", type=float, default=3)
    parser.add_argument("--limit-inputs", type=int, default=4)
    parser.add_argument("--cpu-limit", type=int, default=180)
    parser.add_argument("--wall-limit", type=float, default=180)
    parser.add_argument("--rss-limit-mib", type=int, default=512)
    parser.add_argument(
        "--worker", action="store_true", help=argparse.SUPPRESS
    )
    args = parser.parse_args()
    if (
        args.repetitions < 9
        or args.warmup_seconds < 3
        or not isfinite(args.warmup_seconds)
        or args.cpu_limit <= 0
        or args.wall_limit <= 0
        or not isfinite(args.wall_limit)
        or args.rss_limit_mib <= 0
        or args.limit_inputs <= 0
    ):
        parser.error("require validated samples and positive finite limits")

    if not args.worker:
        supervise(args)
        return
    resource.setrlimit(
        resource.RLIMIT_CPU, (args.cpu_limit, args.cpu_limit + 1)
    )
    for name, expected in HASHES.items():
        if (
            hashlib.sha256((args.source_dir / name).read_bytes()).hexdigest()
            != expected
        ):
            raise ValueError("modified upstream source: " + name)

    import gmpy2
    import sympy

    corpus_path = (
        Path(__file__).parent / "inputs/corpora/phase_three_p34_corpus.json"
    )
    corpus = json.loads(corpus_path.read_text())
    verify_certificates(corpus["certificates"])
    fixtures = [
        f
        for f in corpus["fixtures"]
        if f["split"] == "held_out" and f["band"] == "small"
    ][: args.limit_inputs]
    sys.path.insert(0, str(args.source_dir.resolve()))
    spec = importlib.util.spec_from_file_location(
        "factor_upstream_" + args.mode,
        args.source_dir / (args.mode + ".py"),
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    def cohort():
        return [
            run_one(module, fixture, seed, args.mode == "sssf")
            for fixture in fixtures
            for seed in corpus["seeds"]
        ]

    attempts = []

    for attempt in range(3):
        started, calls = time.perf_counter(), 0
        while time.perf_counter() - started < max(
            args.warmup_seconds, 3 + 2 * attempt
        ):
            cohort()
            calls += 1

        warmup = time.perf_counter() - started
        samples = []

        for _ in range(max(args.repetitions, 9 if attempt == 0 else 15)):
            wall, cpu = time.perf_counter(), time.process_time()
            rows = cohort()
            samples.append(
                dict(
                    seconds=time.perf_counter() - wall,
                    cpu_seconds=time.process_time() - cpu,
                    rows=rows,
                )
            )

        times = [sample["seconds"] for sample in samples]
        median = statistics.median(times)
        q1, _, q3 = statistics.quantiles(times, n=4)
        drift = abs(
            statistics.median(times[:3]) / statistics.median(times[-3:]) - 1
        )
        iqr = (q3 - q1) / median
        stable = drift <= 0.15 and iqr <= 0.2
        attempts.append(
            dict(
                warmup_seconds=warmup,
                warmup_calls=calls,
                samples=samples,
                stable=stable,
                relative_iqr=iqr,
                drift=drift,
                median_seconds=median,
            )
        )
        if stable:
            break

    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak = int(peak if sys.platform == "darwin" else peak * 1024)
    result = dict(
        schema=1,
        upstream_pin=PIN,
        source_sha256=HASHES,
        algorithm_source_changed=False,
        mode=args.mode,
        arguments={"prop": 2, "digred": 0} if args.mode == "sssf" else {},
        seed_policy="upstream global random.seed before each input",
        environment=dict(
            runtime=sys.version,
            sympy=sympy.__version__,
            gmpy2=gmpy2.version(),
            gmp=gmpy2.mp_version(),
        ),
        corpus_sha256=hashlib.sha256(corpus_path.read_bytes()).hexdigest(),
        attempts=attempts,
        median_seconds=median,
        stable=stable,
        peak_rss_bytes=peak,
        completion_per_sample=[
            sum(r["completed"] for r in s["rows"]) for s in samples
        ],
        role=(
            "upstream reproduction; distinct backend/settings; "
            "no promotion ratio"
        ),
        limits=("external CPU/wall/RSS caps; no internal work/storage cap"),
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        args.mode, result["completion_per_sample"], median, stable, flush=True
    )


if __name__ == "__main__":
    main()
