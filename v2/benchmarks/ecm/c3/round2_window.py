"""Own one finite C3 pilot/verification lease after C9's safe release."""

import argparse
import json
import subprocess
import time

from . import c3_study as first
from . import round2_pilot as pilot
from .round2_corpus import TRAINING, load_training, write_training


class WindowAllowanceError(RuntimeError):
    """A next bounded phase cannot fit the remaining active lease."""


def run(output):
    first.runtime()
    started, cpu = time.monotonic(), time.process_time()
    receipt = dict(
        mode="round2-window",
        limit_seconds=570,
        prior_charge_seconds=30,
        phases=[],
    )

    def remaining():
        return 570 - max(time.monotonic() - started, time.process_time() - cpu)

    with first.machine_window():
        try:
            if not TRAINING.exists():
                write_training()
            fixtures = load_training()
            if len(fixtures) != 57:
                raise ValueError("round-two training requires 57 subjects")
            receipt["phases"].append("generated and verified training")

            if remaining() < 120:
                raise WindowAllowanceError("focused QA does not fit")
            command = [
                str(first.ROOT / "v2/.venv/bin/python"),
                "-B",
                "-m",
                "unittest",
                "v2.tests.test_c3_allocation",
                "v2.tests.test_c3_round2",
            ]
            with (output.parent / "focused-qa.log").open("x") as stream:
                completed = subprocess.run(
                    command,
                    cwd=first.ROOT,
                    stdout=stream,
                    stderr=subprocess.STDOUT,
                    timeout=110,
                    check=False,
                )
            if completed.returncode:
                raise RuntimeError("C3 focused verification failed")
            receipt["phases"].append("focused tests passed")

            if not pilot.CONTROL.exists():
                previous = json.loads(first.CONTROL.read_text())
                names = set(previous["sha256"])
                names.update(
                    str(path.relative_to(first.ROOT))
                    for path in (first.ROOT / "v2/benchmarks/ecm/c3").glob(
                        "round2*.py"
                    )
                )
                names.update(
                    (
                        "v2/benchmarks/ecm/c3/round2_protocol_v2.md",
                        "v2/benchmarks/ecm/c3/round2_research.md",
                        str(TRAINING.relative_to(first.ROOT)),
                        (
                            "v2/benchmarks/inputs/controls/"
                            "c3_round2_research_sources.json"
                        ),
                        "v2/benchmarks/inputs/corpora/c3_confirmation.json",
                    )
                )
                first.save(
                    pilot.CONTROL,
                    dict(
                        source_commit=subprocess.check_output(
                            ["git", "rev-parse", "HEAD"],
                            cwd=first.ROOT,
                            text=True,
                        ).strip(),
                        seeds=[17, 43],
                        active_seconds=480,
                        prior_charge_seconds=30,
                        previous_protocol_sha256=first.digest(
                            first.INPUTS / "controls/c3_round2_pilot.json"
                        ),
                        policies=pilot.POLICIES,
                        sha256={
                            name: first.digest(first.ROOT / name)
                            for name in sorted(names)
                        },
                    ),
                )
                subprocess.run(
                    [
                        "git",
                        "add",
                        str(TRAINING.relative_to(first.ROOT)),
                        str(pilot.CONTROL.relative_to(first.ROOT)),
                    ],
                    cwd=first.ROOT,
                    check=True,
                )
                subprocess.run(
                    [
                        "git",
                        "commit",
                        "-m",
                        "Freeze C3 attribution corpus and policy pilot",
                    ],
                    cwd=first.ROOT,
                    check=True,
                )
            pilot.verify_control()
            receipt["phases"].append("committed pilot freeze")

            if remaining() < 480:
                raise WindowAllowanceError("full pilot allowance does not fit")
            pilot.measure(output.parent / "pilot.json", own_window=False)
            receipt["phases"].append("instrumented pilot finished")
        except Exception as error:
            receipt["error"] = str(error)
            raise
        finally:
            receipt["wall"] = time.monotonic() - started
            receipt["cpu"] = time.process_time() - cpu
            first.save(output, receipt)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=first.Path)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    run(args.output)


if __name__ == "__main__":
    main()
