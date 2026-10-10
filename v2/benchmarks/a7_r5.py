"""Pinned R5 comparison arms; preparation does not run E1 confirmation."""

import hashlib
import json
import sys
from dataclasses import asdict
from pathlib import Path

from ..qs import SieveConfig
from ..qs.parallel import ParallelConfig
from ..qs.sss import SSSJob
from . import b1_calibration as b1
from .p38_r1_capacity import decode_config
from .phase_three_sss import deserialize_config

ROOT = Path(__file__).resolve().parents[2]
PLAN = Path(__file__).parent / "inputs/controls/a7_r5_e1_arms.json"


def load_plan():
    """Reject changed runtime/adapter/input bytes before using frozen arms."""
    plan = json.loads(PLAN.read_text())
    if plan["schema"] != 1:
        raise ValueError("unknown R5 arm schema")
    for group in ("source_sha256", "input_sha256"):
        for name, expected in plan[group].items():
            path = ROOT / name
            if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
                raise ValueError("R5 arm pin changed: " + name)

    resources = plan["serial_resources"]
    if resources["work_limit"] != b1.WORK or (
        resources["memory_bytes"] != b1.MEMORY
    ):
        raise ValueError("R5 and B1 resource envelopes differ")
    return plan


def serial_configurations(band, *, plan=None):
    """Decode calibrated controls and separately labelled SSS challengers."""
    plan = load_plan() if plan is None else plan
    configs = {}
    for name, arm in plan["serial"][band].items():
        decoder = (
            deserialize_config if arm["engine"] == "sss" else decode_config
        )
        config = decoder(arm["config"])
        if config.memory_bytes != plan["serial_resources"]["memory_bytes"]:
            raise ValueError("R5 serial arm has a different memory envelope")
        configs[name] = config
    return configs


def worker_configuration(band, *, plan=None):
    """Retain the accepted fixed-family contract, independent of B1 flyer."""
    plan = load_plan() if plan is None else plan
    values = dict(plan["workers"]["configurations"][band])
    values["collector"] = SieveConfig(**values["collector"])
    return ParallelConfig(**values)


def run_serial(fixture, seed, band, arm):
    """Reuse B1's complete-call budget, validation and certainty accounting.

    The caller owns the exclusive performance window and sampling protocol.
    Known factors enter B1 validation only, never collector configuration.
    """
    if sys.implementation.name != "pypy" or sys.version_info[:2] != (3, 11):
        raise RuntimeError(
            "R5 comparisons require PyPy implementing Python 3.11"
        )
    plan = load_plan()
    if seed not in plan["seeds"] or fixture["digits"] != int(band[:-1]):
        raise ValueError("fixture or seed outside the prepared R5 scope")
    config = serial_configurations(band, plan=plan)[arm]
    job_type = (
        SSSJob if plan["serial"][band][arm]["engine"] == "sss" else b1.SIQSJob
    )
    return b1.run_one(
        fixture,
        seed,
        config,
        plan["serial_resources"]["seconds"][band],
        job_type=job_type,
    )


def main():
    plan = load_plan()
    for band in plan["serial"]:
        print(band, ", ".join(serial_configurations(band, plan=plan)))
    for band in plan["workers"]["configurations"]:
        # Construction checks the existing worker limits without collecting.
        asdict(worker_configuration(band, plan=plan))
    print(
        "Prepared only; fresh E1 inputs, timing and worker calibration "
        "remain open."
    )


if __name__ == "__main__":
    main()
