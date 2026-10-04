"""Final matched QS implementation and snapshot-ownership cost checks."""

import argparse
import hashlib
import json
import math
import resource
import subprocess
import sys
import time
from pathlib import Path

from .phase_one import environment
from .phase_three_pipeline import CHANGES, _complete, _config, _corpus, _record
from .phase_three_reference import _rss_bytes
from .qs_snapshot import ROOT, load_qs_arm

FREEZE = ROOT / "benchmarks/phase_three_m25_p33_refinement_pypy.frozen.json"
OLD_PIPELINE = ROOT / "audit/m25_pre_ownership_pipeline.json"
HELD_OUT_SEED = 335


def _arm(name):
    """Use known owned sources; only owner_before restores the old pipeline."""
    arm = load_qs_arm("_qs_final_" + name, () if name == "m22" else CHANGES)
    if name == "owner_before":
        snapshot = json.loads(OLD_PIPELINE.read_text())
        source = snapshot["source"]
        if hashlib.sha256(source.encode()).hexdigest() != snapshot["sha256"]:
            raise ValueError("corrupt owned pipeline snapshot")
        module = arm.qs.pipeline
        exec(compile(source, str(OLD_PIPELINE), "exec"), module.__dict__)
    return arm


def main():
    """Keep the prior freeze, then evaluate a fresh cohort on final sources."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--warmup-seconds", type=float, default=3)
    parser.add_argument("--repetitions", type=int, default=9)
    parser.add_argument("--cold-arm", choices=("m22", "after"))
    args = parser.parse_args()
    frozen = json.loads(FREEZE.read_text())
    controls = frozen["controls"]
    if args.cold_arm:
        arm = _arm(args.cold_arm)
        config = _config(arm, **frozen["config"])
        answers = _complete(
            arm, _corpus(HELD_OUT_SEED, 16), config, **controls
        )
        usage = resource.getrusage(resource.RUSAGE_SELF)
        print(
            json.dumps(
                {
                    "answers": answers,
                    "cpu_to_output_seconds": usage.ru_utime + usage.ru_stime,
                    "peak_rss_bytes": _rss_bytes(),
                }
            )
        )
        return
    if (
        not math.isfinite(args.warmup_seconds)
        or args.warmup_seconds < 3
        or args.repetitions < 9
    ):
        parser.error("require three seconds warmup and nine samples")
    if args.output.exists():
        parser.error("use a unique output capture")
    measured_environment = environment()
    corpus = _corpus(HELD_OUT_SEED, 16)
    expected = [fixture["factors"] for fixture in corpus]
    rows = []
    for name in ("m22", "owner_before", "after"):
        arm = _arm(name)
        config = _config(arm, **frozen["config"])
        rows.append(
            _record(
                name,
                lambda: _complete(arm, corpus, config, **controls),
                expected,
                args,
                details=_complete(
                    arm, corpus, config, details=True, **controls
                ),
            )
        )
    cold = []
    for name in ("m22", "after"):
        for _ in range(9):
            started = time.perf_counter()
            prior = resource.getrusage(resource.RUSAGE_CHILDREN)
            child = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "v2.benchmarks.phase_three_ownership",
                    "--output",
                    str(args.output),
                    "--cold-arm",
                    name,
                ],
                check=True,
                capture_output=True,
                text=True,
                timeout=30,
            )
            result = json.loads(child.stdout)
            assert result["answers"] == expected
            usage = resource.getrusage(resource.RUSAGE_CHILDREN)
            result.update(
                arm=name,
                lifecycle_seconds=time.perf_counter() - started,
                lifecycle_cpu_seconds=(
                    usage.ru_utime
                    + usage.ru_stime
                    - prior.ru_utime
                    - prior.ru_stime
                ),
            )
            cold.append(result)
    assert (
        environment()["source_sha256"]
        == (measured_environment["source_sha256"])
    )
    args.output.write_text(
        json.dumps(
            {
                "milestone": (
                    "M25 / P3.3 final ownership and implementation check"
                ),
                "environment": measured_environment,
                "command": sys.orig_argv,
                "locked_config": frozen,
                "locked_config_source": str(FREEZE),
                "locked_config_sha256": hashlib.sha256(
                    FREEZE.read_bytes()
                ).hexdigest(),
                "held_out_seed": HELD_OUT_SEED,
                "held_out": corpus,
                "results": rows,
                "cold_samples": cold,
                "scope": (
                    "same frozen bucket/filter settings; "
                    "M22 vs final collector"
                ),
                "owner_control": str(OLD_PIPELINE),
                "owner_control_sha256": hashlib.sha256(
                    OLD_PIPELINE.read_bytes()
                ).hexdigest(),
                "budgets": dict(
                    work=200_000_000,
                    wall_seconds=10,
                    cpu_seconds=10,
                    owned_bytes=32 * 1024 * 1024,
                ),
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
