"""Pinned R5 comparison arms; preparation does not run E1 confirmation."""

import hashlib
import json
import sys
from dataclasses import asdict

from ....qs import SieveConfig
from ....qs.parallel import ParallelConfig
from ....qs.sss import SSSJob
from ...support.paths import (
    BENCHMARK_ROOT,
    REPOSITORY_ROOT,
    matches_source_pin,
    source_path,
)
from ..b1 import b1_calibration as b1
from ..p38.p38_r1_capacity import decode_config
from ..phase_three.phase_three_sss import deserialize_config

ROOT = REPOSITORY_ROOT
PLAN = BENCHMARK_ROOT / "inputs/controls/a7_r5_e1_arms.json"
CURRENT_SOURCES = BENCHMARK_ROOT / "inputs/controls/a7_r5_c3_v5_sources.json"
PREVIOUS_SOURCES = BENCHMARK_ROOT / "inputs/controls/a7_r5_c3_sources.json"


def load_plan():
    """Reject changed runtime/adapter/input bytes before using frozen arms."""
    plan_bytes = source_path(PLAN).read_bytes()
    plan = json.loads(plan_bytes)
    if plan["schema"] != 1:
        raise ValueError("unknown R5 arm schema")
    integrated = json.loads(source_path(CURRENT_SOURCES).read_text())
    if (
        integrated["schema"] != 1
        or integrated["base_control_sha256"]
        != hashlib.sha256(plan_bytes).hexdigest()
        or integrated["previous_manifest_sha256"]
        != hashlib.sha256(
            source_path(PREVIOUS_SOURCES).read_bytes()
        ).hexdigest()
        or not integrated["source_sha256"].keys()
        <= plan["source_sha256"].keys()
    ):
        raise ValueError("R5 integrated source pin changed")
    for group in ("source_sha256", "input_sha256"):
        for name, expected in plan[group].items():
            if group == "source_sha256":
                expected = integrated["source_sha256"].get(name, expected)
            path = ROOT / name
            if not matches_source_pin(path, expected):
                raise ValueError("R5 arm pin changed: " + name)

    for name, expected in integrated["additional_source_sha256"].items():
        if not matches_source_pin(ROOT / name, expected):
            raise ValueError("R5 additional runtime pin changed: " + name)

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
