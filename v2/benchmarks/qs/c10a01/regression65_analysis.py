"""Verify host-audited receipts without changing historical timing verdicts."""

import argparse
import hashlib
import json
import statistics
from pathlib import Path

from ..c10a04.accounting import digest, receipt
from . import host_state, live_limits, plan_v4


def timing_check(before, after, expected_power):
    """A passing SIQS sentinel alone is insufficient for a clean comparison."""
    checks = dict(
        host_stable=host_state.matches(before["host"], after["host"]),
        expected_power=all(
            value["host"]["power"] == expected_power
            for value in (before, after)
        ),
        native_digest=(
            before["native_canary"]["digest"]
            == after["native_canary"]["digest"]
        ),
        warmup=all(
            min(value["wall"], value["cpu"]) >= 3
            and len(value["samples"]) == 9
            and len(value["native_canary"]["samples"]) == 9
            for value in (before, after)
        ),
    )
    ratios = {}
    for kind, field in (("siqs", "samples"), ("native", "native_canary")):
        samples = [
            value[field] if kind == "siqs" else value[field]["samples"]
            for value in (before, after)
        ]
        ratios[kind] = {
            key: statistics.median(s[key] for s in samples[1])
            / statistics.median(s[key] for s in samples[0])
            for key in ("wall", "cpu")
        }
        checks[kind + "_drift"] = all(
            abs(ratio - 1) <= 0.15 for ratio in ratios[kind].values()
        )
    return dict(valid=all(checks.values()), checks=checks, ratios=ratios)


def referenced(output, reference):
    name = reference["name"]
    if Path(name).name != name:
        raise ValueError("unsafe receipt path")
    path = output / name
    if path.stat().st_size > 512 * 1024 or digest(path) != reference["sha256"]:
        raise ValueError("receipt size or hash mismatch")
    return json.loads(path.read_text())


def analyze(stage):
    freeze_path = Path(__file__).with_name(f"regression65_{stage}_v2.json")
    frozen = json.loads(freeze_path.read_text())
    output = plan_v4.ROOT / (
        f"v2/benchmarks/results/c10a01/regression65-{stage}-v2"
    )
    state = live_limits.read(output)
    if state["status"] == "running":
        raise RuntimeError("stage is still running; use the live status command")
    if state["identity"] != digest(freeze_path):
        raise ValueError("stage/freeze identity mismatch")
    corpus = json.loads(plan_v4.CORPUS.read_text())
    fixtures = {value["id"]: value for value in corpus["fixtures"]}
    rows, charges = [], dict(wall=0, cpu=0, work=0)
    for frozen_call in frozen["calls"]:
        name = frozen_call["id"]
        call = state["calls"][name]
        if call["configuration"] != frozen_call["configuration"]:
            raise ValueError("policy substitution")
        references = call["receipts"]
        values = {}
        for phase, reference in references.items():
            values[phase] = referenced(output, reference)
            for key in charges:
                charges[key] += reference.get(key, 0)
        for key in charges:
            charges[key] += call.get("reserved_interrupted", {}).get(key, 0)
        if "measured" not in values:
            rows.append(
                dict(call=name, status=call["status"], timing_valid=False)
            )
            continue
        result = values["measured"]
        fixture = fixtures[frozen_call["fixture"]]
        if (result["divisor"] or 1) * result["cofactor"] != fixture["n"]:
            raise ArithmeticError("factor reconstruction failed")
        complete = bool(result["factors"]) and not result["unresolved"]
        if complete and (
            sorted(result["factors"]) != sorted(fixture["factors"])
            or result["reason"] != "factor_found"
        ):
            raise ArithmeticError("complete factor validation failed")
        timing = None
        if "before" in values and "after" in values:
            timing = timing_check(
                values["before"]["result"],
                values["after"]["result"],
                frozen["expected_power"],
            )
        rows.append(
            dict(
                call=name,
                arm=call["arm"],
                complete=complete,
                timing_valid=bool(timing and timing["valid"]),
                timing=timing,
                result=result,
            )
        )
    for key in charges:
        charges[key] += state.get("overhead_charged", {}).get(key, 0)
        if charges[key] > state["limits"]["total"][key]:
            raise ArithmeticError("settled stage exceeds grant")
    return dict(
        identity=frozen["identity"],
        stage_status=state["status"],
        rows=rows,
        charges=charges,
        granted=state["limits"]["total"],
        state_sha256=digest(output / "live.json"),
        analyzer_sha256=hashlib.sha256(
            Path(__file__).read_bytes()
        ).hexdigest(),
        independent_moduli=1,
        adoption=False,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=("screen", "pair"), required=True)
    args = parser.parse_args()
    value = analyze(args.stage)
    path = plan_v4.ROOT / (
        f"v2/audit/c10a01/regression65-{args.stage}-v2-verified.json"
    )
    receipt(path, value)
    print(
        json.dumps(
            dict(
                report=str(path),
                charges=value["charges"],
                timing_valid=[row["timing_valid"] for row in value["rows"]],
            )
        )
    )


if __name__ == "__main__":
    main()
