"""Standalone matched legacy-SIQS control on an immutable runtime snapshot.

Run this file directly. It imports arithmetic only after selecting one runtime.
Use the same corpus, seeds and caps for the current and --baseline-json arms.
"""

import argparse
import hashlib
import importlib
import json
import statistics
import sys
import tempfile
import time
from math import gcd, isqrt, prod
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def verify_fixture(fixture, certificates):
    """Prove the selected corpus factors independently of the runtime."""
    verified = set()

    def prove(n):
        if n in verified:
            return
        proof = certificates[str(n)]
        if proof["kind"] == "trial":
            if not 2 <= n < 2**16 or any(
                n % d == 0 for d in range(2, isqrt(n) + 1)
            ):
                raise ValueError("invalid trial proof")
        elif proof["kind"] == "pocklington":
            q, witness = proof["q"], proof["witness"]
            if not 2 <= q < n or (n - 1) % q:
                raise ValueError("invalid certificate dependency")
            prove(q)
            if not (
                q * q > n
                and pow(witness, n - 1, n) == 1
                and gcd(pow(witness, (n - 1) // q, n) - 1, n) == 1
            ):
                raise ValueError("invalid Pocklington proof")
        else:
            raise ValueError("unknown proof kind")

        verified.add(n)

    if prod(fixture["factors"]) != fixture["n"]:
        raise ValueError("fixture does not reconstruct")
    for factor in fixture["factors"]:
        prove(factor)


def materialize(path, destination):
    """Check and materialize only bounded owned Python sources."""
    data = json.loads(path.read_text())
    if data["schema"] != 1 or len(data["source"]) > 100:
        raise ValueError("invalid R1 baseline")
    for name, source in data["source"].items():
        parts = name.split("/")
        if (
            len(parts) not in (2, 3)
            or parts[0] != "v2"
            or (len(parts) == 3 and parts[1] != "qs")
            or not parts[-1].endswith(".py")
            or any(part in ("", ".", "..") for part in parts)
        ):
            raise ValueError("invalid baseline source path")

        if (
            hashlib.sha256(source.encode()).hexdigest()
            != data["source_sha256"][name]
        ):
            raise ValueError("corrupt baseline source")

        output = destination / name
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(source)

    return data["source_sha256"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--corpus", type=Path, required=True)
    parser.add_argument("--baseline-json", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if sys.implementation.name != "pypy" or sys.version_info[:2] != (3, 11):
        parser.error("PyPy implementing Python 3.11 is required")
    if args.output.exists():
        parser.error("preserve existing results")
    if "v2" in sys.modules:
        parser.error("run the driver directly for import isolation")
    with tempfile.TemporaryDirectory() as directory:
        runtime = Path(directory) if args.baseline_json else ROOT
        if args.baseline_json:
            hashes = materialize(args.baseline_json, runtime)
        else:
            files = sorted((ROOT / "v2").glob("*.py"))
            files += sorted((ROOT / "v2/qs").glob("*.py"))
            hashes = {
                str(p.relative_to(ROOT)): hashlib.sha256(
                    p.read_bytes()
                ).hexdigest()
                for p in files
            }

        sys.path.insert(0, str(runtime))
        budget_module = importlib.import_module("v2.budget")
        qs = importlib.import_module("v2.qs")
        utils = importlib.import_module("v2.utils")
        config = qs.SIQSConfig(
            base_bound=3000,
            half_width=32768,
            max_half_width=32768,
            factor_count=4,
            family_count=64,
            pool_size=32,
            max_stalled=8192,
            max_trivial=4096,
            row_excess=32,
            batch_width=4096,
            memory_bytes=384 * 2**20,
            checkpoint_bytes=16 * 2**20,
            collector=qs.SieveConfig(
                block_width=4096,
                score_policy="powers",
                division="bucket",
                residual_bound=3000**2,
                max_atoms=8192,
                max_relations=2048,
                max_partials=2048,
            ),
        )
        corpus = json.loads(args.corpus.read_text())
        fixture = next(
            f
            for f in corpus["fixtures"]
            if f["kind"] == "balanced" and f["digits"] == 30
        )
        verify_fixture(fixture, corpus["certificates"])

        def cohort():
            rows = []

            for seed in corpus["seeds"]:
                allowance = budget_module.Budget(
                    work_limit=10**13, seconds=30, cpu_seconds=30
                )
                started = time.perf_counter()

                result = qs.SIQSJob(
                    fixture["n"], seed=seed, config=config, budget=allowance
                ).run()
                labels = []
                reason = result.reason
                if result.divisor:
                    if not utils.valid_divisor(result.divisor, fixture["n"]):
                        raise AssertionError("invalid split")
                    try:
                        for child in (result.divisor, result.cofactor):
                            allowance.consume(32 * child.bit_length())
                            label = utils.classify_prime(child)

                            if label is utils.Primality.COMPOSITE:
                                raise AssertionError(
                                    "certified child classified composite"
                                )
                            labels.append(label.value)
                    except budget_module.BudgetExhaustedError:
                        labels = []
                        reason = "classification_" + allowance.reason

                    if (
                        sorted((result.divisor, result.cofactor))
                        != fixture["factors"]
                    ):
                        raise AssertionError(
                            "split disagrees with certified factors"
                        )

                if (result.divisor or 1) * result.cofactor != fixture["n"]:
                    raise AssertionError(
                        "unresolved outcome lost its cofactor"
                    )
                rows.append(
                    dict(
                        seconds=time.perf_counter() - started,
                        seed=seed,
                        divisor=result.divisor,
                        cofactor=result.cofactor,
                        reason=reason,
                        complete=len(labels) == 2,
                        certainty=labels,
                        work_used=allowance.used,
                        stats=result.stats,
                    )
                )

            return rows

        attempts = []

        for attempt in range(3):
            start, warm_calls = time.perf_counter(), 0
            while time.perf_counter() - start < (3 if attempt == 0 else 5):
                cohort()
                warm_calls += 1
            warm_seconds = time.perf_counter() - start
            samples = [cohort() for _ in range(9 if attempt == 0 else 15)]
            times = [sum(r["seconds"] for r in sample) for sample in samples]
            median = statistics.median(times)
            mad = statistics.median(abs(t - median) for t in times)
            stable = (
                mad <= median * 0.1
                and max(times) - min(times) <= median * 0.35
            )
            attempts.append(
                dict(
                    samples=samples,
                    median_seconds=median,
                    mad_seconds=mad,
                    stable=stable,
                    warmup_seconds=warm_seconds,
                    warmup_cohorts=warm_calls,
                )
            )
            if stable:
                break

        for name, digest in hashes.items():
            if (
                hashlib.sha256((runtime / name).read_bytes()).hexdigest()
                != digest
            ):
                raise RuntimeError("runtime source changed during capture")

        result = dict(
            source_sha256=hashes,
            driver_sha256=hashlib.sha256(
                Path(__file__).read_bytes()
            ).hexdigest(),
            runtime=sys.version,
            attempts=attempts,
            corpus_sha256=hashlib.sha256(args.corpus.read_bytes()).hexdigest(),
            scope="One inspected 30-digit input at two fixed seeds. "
            "Setup, split and classification; no checkpoint serialization. "
            "No general scaling or promotion claim.",
        )
        args.output.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
