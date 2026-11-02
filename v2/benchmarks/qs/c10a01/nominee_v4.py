"""Frozen calibration nominee and exact later-pilot configurations."""

import hashlib
import json
from pathlib import Path

from . import plan_v4

HERE = Path(__file__).parent
FREEZE = HERE / "nominee_v4.json"


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def build(*, receipt_hashes=None):
    """Bind the measured receipts and selected literal policies."""
    root = plan_v4.ROOT
    primary = root / "v2/benchmarks/results/c10a01/calibration-v4-r2"
    diagnostic = root / "v2/benchmarks/results/c10a01/diagnostics-v4"
    policies = dict(
        root_update_policy="batch",
        assignment_policy="windowed_pair",
        multiplier_policy="parity2_wide",
        polynomial_policy="parity2",
        marking_policy="parity_two",
    )
    calls = []
    plan = json.loads(plan_v4.PLAN.read_text())
    for stage in ("A", "B", "C"):
        for call in plan["pilots"][stage]["calls"]:
            row = {
                key: call[key]
                for key in (
                    "id",
                    "fixture",
                    "digits",
                    "arm",
                    "seconds",
                    "tier",
                    "work",
                )
            }
            options = (
                {
                    key: policies[key]
                    for key in (
                        "root_update_policy",
                        "assignment_policy",
                        "multiplier_policy",
                    )
                }
                if call["arm"] == "nominee"
                else {}
            )
            row["configuration"] = plan_v4.config_values(
                call["digits"], call["arm"], **options
            )
            calls.append(row)
    if receipt_hashes is None:
        receipt_hashes = (
            digest(primary / "ledger.json"),
            digest(diagnostic / "ledger.json"),
        )
    return dict(
        identity="c10a01-v4-calibrated-nominee-1",
        policies=policies,
        control_commit=plan["control_commit"],
        implementation_commit=plan["implementation_commit"],
        experiment_sha256=digest(plan_v4.PLAN),
        corpus_sha256=plan["corpus_sha256"],
        primary_ledger_sha256=receipt_hashes[0],
        diagnostic_ledger_sha256=receipt_hashes[1],
        calls=calls,
        stages={
            stage: plan["pilots"][stage]["limits"] for stage in ("A", "B", "C")
        },
        selection="Frozen v4 rules: batch, windowed_pair, wide ranking",
        pilot_status="unrun",
    )


def verify(*, require_receipts=False):
    frozen = json.loads(FREEZE.read_text())
    root = plan_v4.ROOT
    receipt_paths = (
        root / "v2/benchmarks/results/c10a01/calibration-v4-r2/ledger.json",
        root / "v2/benchmarks/results/c10a01/diagnostics-v4/ledger.json",
    )
    if require_receipts and not all(path.exists() for path in receipt_paths):
        raise FileNotFoundError("calibration receipts required for execution")
    for path, key in zip(
        receipt_paths,
        ("primary_ledger_sha256", "diagnostic_ledger_sha256"),
    ):
        if path.exists() and digest(path) != frozen[key]:
            raise ValueError("calibration receipt drift")
    expected = build(
        receipt_hashes=(
            frozen["primary_ledger_sha256"],
            frozen["diagnostic_ledger_sha256"],
        )
    )
    if frozen != expected:
        raise ValueError("nominee, pilot resources or receipts drifted")
    return frozen
