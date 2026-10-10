"""Instrumented C6 JIT sensitivity; never use these times for promotion."""

import argparse
import hashlib
import json
import statistics
import time
from collections import Counter
from pathlib import Path

from . import c6_fast_study as study
from . import c6_study
from .build_c6_inputs import require_runtime


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--backend", choices=("int", "gmp"), required=True)
    parser.add_argument("--mode", choices=("tuple", "inline"), required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    require_runtime()
    import pypyjit

    settings = study.protocol()
    arm = f"fast/lucas/{args.mode}/16"
    backend = (
        study.p41_campaign.gmp_backend()
        if args.backend == "gmp"
        else study.p41_campaign.PYTHON_BACKEND
    )
    cases = [
        c for c in study.corpus()["fixtures"] if c["id"].startswith("balanced")
    ]
    compiled, aborted, generated_codes = Counter(), Counter(), set()

    def compiled_hook(info):
        if info.jitdriver_name != "pypyjit":
            return
        code = info.greenkey[2]
        key = Path(code.co_filename).name + ":" + code.co_name
        if len(compiled) < 512 or key in compiled:
            compiled[key] += 1
        if code.co_filename == "<verified-c6-chain>":
            generated_codes.add(code)

    def aborted_hook(driver, greenkey, reason, operations):
        if len(aborted) < 32 or reason in aborted:
            aborted[reason] += 1

    with c6_study.performance_window("instrumented C6 JIT sensitivity"):
        program = study.construct(arm, backend)
        controls = {
            (case["id"], seed): c6_study.affine_controls(
                case["n"], tuple(case["factors"]), seed, 2000
            )
            for seed in settings["seeds"]
            for case in cases
        }
        pypyjit.set_compile_hook(compiled_hook, operations=False)
        pypyjit.set_abort_hook(aborted_hook)
        started = time.perf_counter()

        def sample():
            if time.perf_counter() - started > 180:
                raise RuntimeError("JIT diagnostic wall cap")
            rows = [
                study.attempt(
                    case,
                    seed,
                    arm,
                    backend,
                    settings,
                    "stage_reuse",
                    program,
                    controls[case["id"], seed],
                )
                for seed in settings["seeds"]
                for case in cases
            ]
            return sum(row["seconds"] for row in rows)

        phases = []
        for warmup in (3, 20):
            elapsed, runs = 0.0, 0
            while elapsed < warmup:
                elapsed += sample()
                runs += 1
            samples = [sample() for _ in range(9)]
            snapshot = pypyjit.get_stats_snapshot()
            phases.append(
                dict(
                    additional_warmup_seconds=elapsed,
                    warmup_cohorts=runs,
                    instrumented_samples=samples,
                    instrumented_median=statistics.median(samples),
                    compiled_roots=dict(compiled),
                    aborted=dict(aborted),
                    generated_root_codes=len(generated_codes),
                    counters=snapshot.counters,
                    counter_times=snapshot.counter_times,
                    assembler_bytes=pypyjit.get_stats_asmmemmgr(),
                )
            )
        pypyjit.set_compile_hook(None)
        pypyjit.set_abort_hook(None)
        report = dict(
            arm=arm,
            backend=args.backend,
            phases=phases,
            defaults=pypyjit.defaults,
            source_sha256=hashlib.sha256(
                Path(__file__).read_bytes()
            ).hexdigest(),
            scope="Instrumented diagnostic, not paired performance evidence. "
            "Root hooks do not enumerate every inlined function.",
        )
        with args.output.open("x") as output:
            json.dump(report, output, indent=2)
            output.write("\n")


if __name__ == "__main__":
    main()
