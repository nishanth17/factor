"""Restore the immutable measured R1 sources into a new scratch directory."""

import argparse
import hashlib
import json
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    if sys.implementation.name != "pypy" or sys.version_info[:2] != (3, 11):
        parser.error("PyPy implementing Python 3.11 is required")
    if args.destination.exists():
        parser.error("destination must not already exist")
    data = json.loads((HERE / "p38_r1_measured_sources.json").read_text())
    if data["schema"] != 1 or len(data["source"]) > 128:
        raise ValueError("invalid measured-source manifest")
    for name, content in data["source"].items():
        parts = name.split("/")
        if (
            len(parts) not in (2, 3)
            or parts[0] != "v2"
            or (
                len(parts) == 3
                and parts[1] not in ("qs", "benchmarks", "tests")
            )
            or any(part in ("", ".", "..") for part in parts)
            or not parts[-1].endswith(".py")
        ):
            raise ValueError("invalid measured-source path")
        if (
            hashlib.sha256(content.encode()).hexdigest()
            != data["source_sha256"][name]
        ):
            raise ValueError("corrupt measured source")
    for name, content in data["source"].items():
        output = args.destination / name
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(content)
    for name in (
        "p38_r1_training_corpus.json",
        "p38_r1_confirmation_corpus.json",
        "p38_r1_frozen.json",
        "p38_r1_baseline.json",
    ):
        shutil.copyfile(HERE / name, args.destination / "v2/benchmarks" / name)


if __name__ == "__main__":
    main()
