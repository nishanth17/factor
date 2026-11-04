"""Inspect frozen fixed-work receipts without importing factoring code."""

import argparse
import json
import math
from pathlib import Path


def summarise(freeze, records):
    """Validate declared curves and compare complete cold/warm observations.

    Warm CPU counters describe the measured call. Whole-cell CPU also includes
    warmup and teardown; those values must not be divided by warm-call latency.
    """
    cells = {cell["id"]: cell for cell in freeze["cells"]}
    seen, signatures, rows = set(), {}, []
    failures = []
    for record in records:
        cell = record["cell"]
        identity = cell["id"]
        if (
            record["identity"] != freeze["identity"]
            or cells.get(identity) != cell
            or identity in seen
        ):
            raise ValueError("duplicate, foreign or changed cell receipt")
        seen.add(identity)
        expected_root = freeze[
            "control_root" if cell["arm"] == "S0" else "candidate_root"
        ]
        if record["source_root"] != expected_root:
            raise ValueError("receipt source differs from frozen arm")
        child = record.get("child")
        if (
            child is None
            or not child["complete"]
            or record["timed_out"]
            or record["returncode"]
            or record["failure"]
        ):
            failures.append(identity)
            continue
        n = freeze["inputs"][cell["input"]]["n"]
        if child["n"] != n or child["reconstruct"] != n:
            raise ValueError("receipt input or reconstruction changed")
        config = freeze.get("configurations", {}).get(
            cell["input"], freeze.get("config", {})
        )
        curves = sum(tier[2] for tier in config["ecm_tiers"])
        hits, chunks = None, None
        checkpoint = child.get("execution_checkpoint")
        if cell["arm"] != "S0":
            completed = checkpoint["search"]["payload"]["completed"]
            if set(completed) != {str(i) for i in range(curves)}:
                raise ValueError("declared curve set incomplete or duplicated")
            signature, ids, hits, chunks = [], set(), 0, 0
            for position in range(curves):
                result = completed[str(position)]
                divisor = result["factor"]
                if divisor is not None:
                    if type(divisor) is not int or not (
                        1 < divisor < n and n % divisor == 0
                    ):
                        raise ValueError("invalid proper divisor")
                    hits += 1
                if result["reason"] != "complete" or result["id"] in ids:
                    raise ValueError("unfinished or duplicated curve identity")
                ids.add(result["id"])
                chunks += result["chain_chunks"]
                signature.append(
                    (
                        result["id"],
                        result["sigma"],
                        divisor,
                        result["phase"],
                        result["chain_chunks"],
                    )
                )
            if signatures.setdefault(cell["input"], signature) != signature:
                raise ValueError("search differs between equivalent arms")
        stats = child["startup_and_shutdown_stats"]
        last = stats[-1] if stats else {}
        warm = cell["mode"] == "warm"
        if warm and child["warmup_seconds"] < 3:
            raise ValueError("warm context lacks validated warmup")
        duration = child["warm_call_seconds"] if warm else record["wall"]
        if not math.isfinite(duration) or duration <= 0:
            raise ValueError("invalid measured duration")
        rows.append(
            dict(
                cell=identity,
                input=cell["input"],
                arm=cell["arm"],
                mode=cell["mode"],
                seconds=duration,
                curves=curves,
                curves_with_factors=hits,
                chain_chunks=chunks,
                whole_cell_observed_cpu=child["aggregate_cpu"],
                call_parent_cpu=last.get("parent_cpu"),
                call_reported_child_cpu=last.get("reported_child_cpu"),
                pool_start_seconds=last.get("startup_seconds"),
                preparations=child["pool_lifetime_stats"].get("preparations"),
                worker_hwm_sum=child["pool_lifetime_stats"].get(
                    "worker_hwm_sum_peak"
                ),
                parent_hwm=child["parent_peak_rss"],
            )
        )
    comparisons = []
    for row in rows:
        if row["arm"] not in ("S1", "P2", "P4"):
            continue
        controls = [
            control
            for control in rows
            if control["input"] == row["input"]
            and control["mode"] == row["mode"]
            and control["arm"] == "M0"
        ]
        if controls:
            control = controls[0]
            comparisons.append(
                dict(
                    input=row["input"],
                    mode=row["mode"],
                    arm=row["arm"],
                    serial_seconds=control["seconds"],
                    candidate_seconds=row["seconds"],
                    exploratory_ratio=control["seconds"] / row["seconds"],
                )
            )
    return dict(
        identity=freeze["identity"],
        rows=rows,
        comparisons=comparisons,
        failures=failures,
        missing=sorted(cells.keys() - seen),
        inference="exploratory fixed work; no quantile or adoption claim",
        cpu_coverage="self/reaped children; helper incomplete; retain charges",
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("freeze", type=Path)
    parser.add_argument("receipts", type=Path)
    args = parser.parse_args()
    freeze = json.loads(args.freeze.read_text())
    records = [
        json.loads((args.receipts / (cell["id"] + ".json")).read_text())
        for cell in freeze["cells"]
        if (args.receipts / (cell["id"] + ".json")).exists()
    ]
    print(json.dumps(summarise(freeze, records), indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
