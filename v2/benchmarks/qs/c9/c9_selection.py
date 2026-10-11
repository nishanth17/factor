"""Select C9 policies from complete training costs before fresh generation."""

import argparse
import hashlib
import json
from pathlib import Path

from ..b1.b1_calibration import save
from . import c9_combinations as combinations
from . import c9_confirmation as confirmation
from . import c9_study as study


def ranked(rows, arms):
    """Every failure costs its full cap; ties retain the declared arm order."""
    result = []
    for order, arm in enumerate(arms):
        samples = [r for r in rows if r["arm"] == arm]
        if len(samples) != 2:
            raise ValueError("selection requires both training inputs per arm")
        result.append(
            dict(
                arm=arm,
                complete=sum(r["complete"] for r in samples),
                capped_seconds=sum(r["capped_seconds"] for r in samples),
                order=order,
            )
        )
    return sorted(
        result,
        key=lambda r: (-r["complete"], r["capped_seconds"], r["order"]),
    )


def freeze(training, joint, spent_wall, spent_cpu):
    if confirmation.CORPUS.exists():
        raise ValueError("fresh corpus already revealed; no reselection")
    if confirmation.SELECTED.exists():
        raise ValueError("selection already frozen")
    initial = json.loads(study.CONTROL.read_text())
    extra = json.loads(combinations.CONTROL.read_text())
    if initial["source"] != study.hashes():
        raise ValueError("initial study source changed")
    if extra["source"] != combinations.hashes():
        raise ValueError("combination study source changed")
    for directory, control in ((training, initial), (joint, extra)):
        if json.loads((directory / "manifest.json").read_text()) != control:
            raise ValueError("training manifest differs from control")
    paths = sorted(training.glob("run-*.json"))
    paths += sorted(joint.glob("run-*.json"))
    rows = [json.loads(path.read_text()) for path in paths]
    if len(rows) != 14:
        raise ValueError("C9 requires ten individual and four joint calls")
    slp_rank = ranked(rows, ("slp", "slp_aligned"))
    dlp_rank = ranked(
        rows, ("dlp", "aligned", "retention", "half", "half_aligned")
    )
    registry = initial["configurations"] | extra["configurations"]
    c1 = json.loads((study.INPUTS / "controls/c1_selected.json").read_text())[
        "configurations"
    ]
    configs = {
        str(digits): dict(
            slp=c1[str(digits)]["slp"],
            dlp=c1[str(digits)]["dlp_half"],
            challenger=c1[str(digits)]["dlp_half"],
        )
        for digits in (30, 40)
    }
    configs["50"] = dict(
        slp=registry[slp_rank[0]["arm"]],
        dlp=registry["dlp"],
        challenger=registry[dlp_rank[0]["arm"]],
    )
    if min(spent_wall, spent_cpu) <= 0 or max(spent_wall, spent_cpu) >= 14400:
        raise ValueError("invalid study expense ledger")
    save(
        confirmation.SELECTED,
        dict(
            schema=1,
            source=confirmation.hashes(),
            configurations=configs,
            choices={"slp": slp_rank, "dlp": dlp_rank},
            spent_before_confirmation=dict(wall=spent_wall, cpu=spent_cpu),
            evidence=[
                dict(
                    name=str(p),
                    sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
                    outcome={
                        k: row[k]
                        for k in (
                            "id",
                            "arm",
                            "complete",
                            "reason",
                            "seconds",
                            "cpu_seconds",
                            "capped_seconds",
                            "work",
                        )
                    },
                )
                for p, row in zip(paths, rows)
            ],
            scope="Selected 50-digit balanced candidate only. Current C1 "
            "30/40-digit comparison bundles remain unchanged. C1's 30-digit "
            "DLP comparison was not generally promoted; this is a regression "
            "control, not a new smaller-band recommendation.",
            selection_rule="Completion, then capped total cost, then stable "
            "predeclared arm order; no tuning on fresh inputs.",
        ),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--training", type=Path, required=True)
    parser.add_argument("--joint", type=Path, required=True)
    parser.add_argument("--spent-wall", type=float, required=True)
    parser.add_argument("--spent-cpu", type=float, required=True)
    args = parser.parse_args()
    study.runtime()
    with study.machine_window():
        freeze(args.training, args.joint, args.spent_wall, args.spent_cpu)


if __name__ == "__main__":
    main()
