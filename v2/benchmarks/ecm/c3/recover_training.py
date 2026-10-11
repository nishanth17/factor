"""Finish interrupted C3 training receipts within the original allowance."""

import argparse
import json
import os
import subprocess
import time
from pathlib import Path

from . import c3_study as study

RECOVERY = study.INPUTS / "controls/c3_training_recovery.json"


def recover(output):
    """Rewarm unfinished inputs and run only missing frozen assignments."""
    study.runtime()
    protocol = study.verify_protocol()
    recovery = json.loads(RECOVERY.read_text())
    if recovery["protocol_sha256"] != study.digest(study.CONTROL):
        raise ValueError("recovery protocol changed")
    source = Path(__file__)
    if recovery["runner_sha256"] != study.digest(source):
        raise ValueError("recovery runner changed")
    for path in (RECOVERY, source):
        committed = subprocess.check_output(
            ["git", "show", "HEAD:" + str(path.relative_to(study.ROOT))]
        )
        if committed != path.read_bytes():
            raise ValueError("commit the recovery freeze before execution")

    folder = output.parent / (output.stem + "-rows")
    rows = []
    observed = set()
    for name, expected in recovery["receipts_sha256"].items():
        path = folder / name
        if study.digest(path) != expected:
            raise ValueError("interrupted receipt changed")
        row = json.loads(path.read_text())
        key = (row["id"], row["seed"], row["arm"], row["repetition"])
        if key in observed:
            raise ValueError("duplicated interrupted receipt")
        observed.add(key)
        rows.append({k: v for k, v in row.items() if k != "events"})
    fixtures = study.training()
    seeds = protocol["training_seeds"]
    arms = list(study.POLICIES)
    expected = {
        (fixture["id"], seed, arm, repetition)
        for fixture in fixtures
        for seed in seeds
        for arm in arms
        for repetition in range(3)
    }
    if not observed <= expected or len(rows) != recovery["completed_calls"]:
        raise ValueError("unexpected interrupted assignments")
    if set(path.name for path in folder.glob("*.json")) != set(
        recovery["receipts_sha256"]
    ):
        raise ValueError("unexpected additional receipts")
    if output.exists():
        raise ValueError("training aggregate already exists")

    control = study.baseline()
    started, cpu_started = time.monotonic(), time.process_time()
    allowance = protocol["envelopes"]["train"] - recovery["charged_seconds"]

    def allowed():
        return (
            max(time.monotonic() - started, time.process_time() - cpu_started)
            <= allowance
        )

    def finish(stopped=None):
        capture = dict(
            mode="train",
            rows=rows,
            wall=recovery["charged_seconds"] + time.monotonic() - started,
            cpu=recovery["charged_seconds"]
            + time.process_time()
            - cpu_started,
            runtime=study.sys.version,
            protocol_sha256=study.digest(study.CONTROL),
            recovery_sha256=study.digest(RECOVERY),
            prior_elapsed_seconds=recovery["prior_elapsed_seconds"],
            prior_charge_seconds=recovery["charged_seconds"],
        )
        if stopped:
            capture["stopped"] = stopped
        study.save(output, capture)

    with study.machine_window():
        for fixture in fixtures:
            if all(
                key in observed for key in expected if key[0] == fixture["id"]
            ):
                continue
            # A process interruption loses JIT state. Rewarm every arm for the
            # unfinished input; this cost consumes the remaining study grant.
            for arm in arms:
                began = time.monotonic()
                while time.monotonic() - began < 3:
                    if not allowed():
                        finish("study_allowance")
                        return
                    study.run_one(fixture, seeds[0], arm, control)
            for repetition in range(3):
                for seed in seeds:
                    offset = repetition % len(arms)
                    order = arms[offset:] + arms[:offset]
                    if (repetition // len(arms)) % 2:
                        order = list(reversed(order))
                    for arm in order:
                        key = (fixture["id"], seed, arm, repetition)
                        if key in observed:
                            continue
                        if not allowed():
                            finish("study_allowance")
                            return
                        study.check_quiet({os.getpid()})
                        row = study.run_one(fixture, seed, arm, control)
                        row["repetition"] = repetition
                        name = (
                            f"{fixture['id']}-{seed}-{repetition}-{arm}.json"
                        )
                        study.save(folder / name, row)
                        rows.append(
                            {k: v for k, v in row.items() if k != "events"}
                        )
                        observed.add(key)
            print(
                "recovered",
                fixture["id"],
                "total calls",
                len(rows),
                flush=True,
            )
    if observed != expected:
        raise ValueError("recovery did not finish training assignments")
    finish()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    recover(parser.parse_args().output)
