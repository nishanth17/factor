"""Recorded stage decisions, without requiring an unrelated control to win."""

import json
from pathlib import Path


def outcomes(ledger, stage):
    result = []
    for row in ledger.data["calls"]:
        if row["phase"] == stage and "sample" in row.get("identity", {}):
            if row["status"] != "complete":
                raise RuntimeError("unsettled or failed cell needs an audit")
            value = json.loads((ledger.output / row["receipt"]).read_text())[
                "result"
            ]
            # The shared ledger's preserved lost-cell receipt predates C9a's
            # fields. Normalize its analysis view, never rewrite the receipt.
            identity = row["identity"]
            value.setdefault("fixture", identity["fixture"])
            value.setdefault("uneven", identity.get("uneven", False))
            value.setdefault("timing_valid", False)
            value.setdefault(
                "observation",
                "capture" if identity.get("kind") == "capture" else "clean",
            )
            result.append(value)
    return result


def admit(stage, decision, freeze, ledger):
    """Record usefulness and cost without changing a freeze."""
    if (
        decision.get("freeze") != freeze["sha256"]
        or decision.get("stage") != stage
    ):
        raise ValueError("stage decision belongs to another freeze")
    for field in (
        "question",
        "evidence",
        "reason",
        "stop_rule",
        "expected_seconds",
    ):
        if not decision.get(field):
            raise ValueError("stage decision lacks " + field)
    if decision.get("advance") is not True:
        raise RuntimeError("recorded decision does not advance this stage")
    estimate = decision["expected_seconds"]
    if (
        type(estimate) is not int
        or not 0 < estimate <= freeze["limits"][stage]["wall"]
    ):
        raise ValueError("expected stage runtime exceeds its ceiling")
    if stage not in ("A", "C"):
        early = outcomes(ledger, "A")
        if (
            len(
                [
                    r
                    for r in early
                    if r.get("fixture", "").startswith("exposed")
                    and r.get("observation") == "clean"
                    and r.get("cap_seconds", 0) > 10
                ]
            )
            < 40
        ):
            raise RuntimeError("primary early cells are incomplete")
    if stage in ("B2", "M") and len(outcomes(ledger, "B1")) != 4:
        raise RuntimeError("clean B1 anchors precede this stage")
    if stage == "M":
        observer = ledger.output / "A-observer-analysis.json"
        if not observer.exists():
            raise RuntimeError("observer overhead check is missing")
        if not decision.get("missing_mature_question"):
            raise ValueError("capture needs a named remaining mature question")
        projection = decision.get("projected_capture_bytes")
        if type(projection) is not int or not 0 < projection <= 192 * 2**20:
            raise ValueError(
                "useful observation window lacks 25% recorder headroom"
            )
        if not decision.get("capacity_evidence"):
            raise ValueError("capture projection lacks atom/window accounting")
    if stage == "F":
        large = [r for r in outcomes(ledger, "B2") if r["digits"] == 70]
        if len(large) != 4 or not all(
            r["complete"]
            and r["timing_valid"]
            and max(r["seconds"], r["cpu_seconds"]) < 5000
            for r in large
        ):
            raise RuntimeError(
                "F requires all four clean 70 calls below 5000 s"
            )
        if not decision.get("own_75_feasibility"):
            raise ValueError(
                "F also requires its own 75-digit resource/yield evidence"
            )
    if stage == "C":
        if not decision.get("fresh_selection"):
            raise ValueError(
                "confirmation requires a prior frozen one-band selection"
            )
    # Referenced receipts must exist; their contents remain reviewer-visible.
    if not isinstance(decision["evidence"], list) or not all(
        Path(path).is_file() for path in decision["evidence"]
    ):
        raise ValueError("stage decision references missing evidence")
