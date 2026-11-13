"""Finite untimed build/correctness admission for the optional NEON screen."""

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
    parser.add_argument("--window-released", action="store_true")
    parser.add_argument("--identity", required=True)
    parser.add_argument("--receipts", type=Path, required=True)
    parser.add_argument("--gmp-prefix", default="/opt/homebrew/opt/gmp")
    args = parser.parse_args()
    if not args.window_released:
        parser.error("a released shared window is required for executable QA")
    if sys.implementation.name != "pypy" or sys.version_info[:2] != (3, 11):
        raise RuntimeError("PyPy Python 3.11 required")
    if subprocess.check_output(
        ["git", "status", "--porcelain"], cwd=ROOT, text=True
    ):
        raise RuntimeError(
            "commit the exact preparation before NEON admission"
        )
    commit = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
    ).strip()
    args.receipts = args.receipts.resolve()
    args.receipts.mkdir(parents=True, exist_ok=True)
    records, passed = [], False
    ledger_path = ROOT / "v2/audit/b4a/phase12-neon-v3-ledger.json"
    if (
        ledger_path.exists()
        and len(json.loads(ledger_path.read_text())["attempts"]) >= 2
    ):
        raise ValueError("two Phase1-2 lock decision points already used")
    with open("/private/tmp/factor-performance.lock", "r+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        ledger = (
            json.loads(ledger_path.read_text())
            if ledger_path.exists()
            else dict(wall=0, cpu=0, attempts=[])
        )
        if (
            args.identity in ledger["attempts"]
            or ledger["wall"] + 160 > 600
            or ledger["cpu"] + 320 > 2400
        ):
            raise ValueError("phase1-2 cumulative admission grant refused")
        ledger["wall"] += 160
        ledger["cpu"] += 320
        ledger["attempts"].append(args.identity)
        save(ledger_path, ledger)
        environment = execution_environment()
        environment["FACTOR_B4A_NEON_QA"] = "1"
        environment["RAYON_NUM_THREADS"] = "1"
        commands = [
            (
                [
                    sys.executable,
                    "-B",
                    "-m",
                    "v2.ecm.native.neon_build",
                    "--gmp-prefix",
                    args.gmp_prefix,
                ],
                40,
                80,
            ),
            (
                [
                    sys.executable,
                    "-B",
                    "-m",
                    "unittest",
                    "v2.tests.test_neon_ecm",
                    "v2.tests.test_b4a_neon_screen",
                    "-v",
                ],
                60,
                120,
            ),
            (["make", "-C", "v2", "lint"], 60, 120),
        ]
        try:
            for index, (command, wall, cpu) in enumerate(commands):
                record = run_owned(
                    command,
                    args.receipts / f"{args.identity}-{index}",
                    dict(
                        wall=wall,
                        cpu=cpu,
                        shutdown=2,
                        rss=4 * 2**30,
                        output=8 * 2**20,
                        processes=4,
                        threads=8,
                    ),
                    env=environment,
                )
                records.append(record)
                if (
                    record["exit"] != 0
                    or record["failure"]
                    or not record["group_empty"]
                    or record["observed_cpu"] > cpu
                    or record["wall"] > wall
                    or record["guardian"].get("reason")
                    not in (None, "controller_closed")
                ):
                    raise RuntimeError(
                        "untimed NEON admission failed; reservation retained"
                    )
            if (
                subprocess.check_output(
                    ["git", "status", "--porcelain"], cwd=ROOT, text=True
                )
                or subprocess.check_output(
                    ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
                ).strip()
                != commit
            ):
                raise RuntimeError("source changed during NEON admission")
            passed = True
        finally:
            manifest = ROOT / "v2/ecm/native/build/neon/manifest.json"
            receipt = dict(
                identity=args.identity,
                commit=commit,
                passed=passed,
                stage="NEON-correctness-admission",
                timing_executed=False,
                reservation=dict(wall=160, cpu=320),
                records=records,
                manifest_sha256=sha(manifest) if manifest.exists() else None,
                pypy=sys.executable,
            )
            write_new(args.receipts / (args.identity + ".json"), receipt)


if __name__ == "__main__":
    main()
