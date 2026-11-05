"""Locked finite build/correctness checks; Stage Q is a future explicit run."""

import argparse
import fcntl
import json
import subprocess
import sys
from pathlib import Path

from .runner import (
    ROOT,
    execution_environment,
    run_owned,
    save,
    sha,
    write_new,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=("tiny", "Q"), required=True)
    parser.add_argument("--window-released", action="store_true")
    parser.add_argument("--receipts", type=Path, required=True)
    parser.add_argument("--identity", required=True)
    parser.add_argument(
        "--gmp-prefix", type=Path, default=Path("/opt/homebrew/opt/gmp")
    )
    parser.add_argument(
        "--openmp-prefix", type=Path, default=Path("/opt/homebrew/opt/libomp")
    )
    args = parser.parse_args()
    if not args.window_released:
        parser.error(
            "lock availability also requires no competing accepted timing"
        )
    if sys.implementation.name != "pypy" or sys.version_info[:2] != (3, 11):
        raise RuntimeError("PyPy implementing Python 3.11 required")
    args.receipts = args.receipts.resolve()
    args.receipts.mkdir(parents=True, exist_ok=True)
    ceiling = (600, 2400) if args.stage == "tiny" else (1800, 7200)
    reservation = (180, 600) if args.stage == "tiny" else ceiling
    commit = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
    ).strip()
    if args.stage == "Q" and subprocess.check_output(
        ["git", "status", "--porcelain"], cwd=ROOT, text=True
    ):
        raise RuntimeError(
            "Stage Q requires the exact combined committed source"
        )
    commands = [
        [
            sys.executable,
            "-B",
            "-m",
            "v2.ecm.native.build",
            "--gmp-prefix",
            str(args.gmp_prefix),
        ],
        [
            sys.executable,
            "-B",
            "-m",
            "unittest",
            "v2.tests.test_native_ecm",
            "v2.tests.test_native_team",
            "v2.tests.test_b4a_runner",
            "-v",
        ],
        ["make", "-C", "v2", "lint"],
    ]
    commands.insert(
        1,
        commands[0]
        + ["--openmp-prefix", str(args.openmp_prefix), "--vector-report"],
    )
    if args.stage == "Q":
        commands[2].insert(-1, "v2.tests.test_b4a_ownership")
        commands.append(
            ["make", "-C", "v2", "test", "PYTHON=" + sys.executable]
        )
        # A committed-source archive independently builds, tests and imports.
        # Dependencies are explicit; no ignored source or input is copied.
        commands.append(
            [
                sys.executable,
                "-B",
                "-m",
                "v2.benchmarks.ecm.b4a.clean_check",
                "--gmp-prefix",
                str(args.gmp_prefix),
                "--openmp-prefix",
                str(args.openmp_prefix),
            ]
        )
    with open("/private/tmp/factor-performance.lock", "r+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        ledger_path = args.receipts / (args.stage + "-ledger.json")
        ledger = (
            json.loads(ledger_path.read_text())
            if ledger_path.exists()
            else dict(wall=0, cpu=0, attempts=[])
        )
        if (
            args.identity in ledger["attempts"]
            or ledger["wall"] + reservation[0] > ceiling[0]
            or ledger["cpu"] + reservation[1] > ceiling[1]
        ):
            raise ValueError("QA identity/cumulative reservation refused")
        ledger["wall"] += reservation[0]
        ledger["cpu"] += reservation[1]
        ledger["attempts"].append(args.identity)
        save(ledger_path, ledger)
        records = []
        passed = False
        environment = execution_environment()
        if args.stage == "Q":
            environment["FACTOR_B4A_LIFECYCLE_QA"] = "1"
        # Per-command maxima partition the same reservation; unused slots do
        # not enlarge any command or authorize extra attempts.
        walls = (
            [20, 20, 90, 50]
            if args.stage == "tiny"
            else [40, 40, 120, 100, 600, 900]
        )
        cpus = (
            [70, 70, 300, 160]
            if args.stage == "tiny"
            else [160, 160, 480, 400, 2400, 3600]
        )
        try:
            for index, (command, wall, cpu) in enumerate(
                zip(commands, walls, cpus)
            ):
                record = run_owned(
                    command,
                    args.receipts / f"{args.identity}-{index}",
                    dict(
                        wall=wall,
                        cpu=cpu,
                        shutdown=2,
                        rss=4 * 2**30,
                        output=16 * 2**20,
                        processes=10 if args.stage == "Q" else 6,
                        threads=48 if args.stage == "Q" else 24,
                    ),
                    env=environment,
                )
                records.append(record)
                if (
                    record["exit"] != 0
                    or record["failure"]
                    or not record["group_empty"]
                    or record["wall"] > wall
                    or record["observed_cpu"] > cpu
                    or record["guardian"].get("reason")
                    not in (None, "controller_closed")
                ):
                    raise RuntimeError(
                        "QA command failed; later commands not launched"
                    )
            passed = True
            if args.stage == "Q" and (
                subprocess.check_output(
                    ["git", "status", "--porcelain"], cwd=ROOT, text=True
                )
                or subprocess.check_output(
                    ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
                ).strip()
                != commit
            ):
                passed = False
                raise RuntimeError("source changed during Stage Q")
        finally:
            receipt = dict(
                stage=args.stage,
                identity=args.identity,
                commit=commit,
                pypy=sys.executable,
                passed=passed,
                records=records,
                reservation=reservation,
                pins={
                    str(p.relative_to(ROOT)): sha(p)
                    for p in ROOT.glob("v2/ecm/native/packet.*")
                },
            )
            write_new(args.receipts / (args.identity + ".json"), receipt)


if __name__ == "__main__":
    main()
