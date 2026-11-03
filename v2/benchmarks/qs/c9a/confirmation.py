"""One selection before fresh sampling; preserve every exposed charge."""

import hashlib
import json
from pathlib import Path

from ..c10a04.accounting import atomic, digest, output_lock, receipt
from .accounting import StudyLedger
from .gates import outcomes
from .inputs import validate
from .plan import CONFIRMATION_SECONDS, confirmation_cells, limits
from .preparation import (
    INPUTS,
    ROOT,
    load_freeze,
    prepare,
    require_clean,
    source_manifest,
)


def select(freeze_path, output, band, decision):
    frozen = load_freeze(freeze_path)
    if band not in CONFIRMATION_SECONDS:
        raise ValueError("unfrozen confirmation band")
    output = Path(output)
    with output_lock(output):
        ledger = StudyLedger(output, frozen["sha256"], frozen["limits"])
        if (
            len(outcomes(ledger, "B1")) != 4
            or len(outcomes(ledger, "B2")) != 20
        ):
            raise RuntimeError(
                "broader exposed completion screen is incomplete"
            )
        if len(outcomes(ledger, "R")) != 24:
            raise RuntimeError("lower-band regression screen is incomplete")
        if not all(
            decision.get(key)
            for key in (
                "reason",
                "uncertainty",
                "regressions",
                "expected_seconds",
                "stop_rule",
            )
        ):
            raise ValueError(
                "selection lacks its exposed cost/completion assessment"
            )
        selection = dict(
            parent_freeze=frozen["sha256"],
            band=band,
            nominee="C1",
            decision=decision,
            ledger_sha256=digest(ledger.path),
            exposed_charges={
                k: ledger.used(field=k) for k in ("wall", "cpu", "work")
            },
            source=frozen["source"],
            parent_output=str(output.resolve()),
        )
        receipt(output / "confirmation-selection.json", selection)
        return selection


def fresh(selection_path, exclusions):
    selection = json.loads(Path(selection_path).read_text())
    parent = Path(selection["parent_output"])
    if digest(parent / "ledger.json") != selection["ledger_sha256"]:
        raise ValueError("exposed ledger changed after fresh selection")
    output = INPUTS / ("confirmation_" + str(selection["band"]))
    # Include the entire exposed corpus plus every certified factor. This
    # supplements, rather than replaces, the historical overlap inventory.
    excluded = json.loads(Path(exclusions).read_text())
    exposed = json.loads((INPUTS / "corpus.json").read_text())
    values = set(excluded["integers"])
    for row in exposed["fixtures"]:
        values.update([row["n"], *row["factors"]])
    combined = parent / "confirmation-exclusions.json"
    receipt(combined, dict(integers=sorted(values)))
    resources = dict(
        C=dict(wall=1200, cpu=1200, work=10**13, disk_bytes=128 * 2**20),
        total=dict(wall=1200, cpu=1200, work=10**13),
        disk_bytes=128 * 2**20,
    )
    directory = parent / "confirmation-preparation"
    with output_lock(directory):
        ledger = StudyLedger(
            directory, digest(Path(selection_path)), resources
        )
        return ledger.call_operation(
            "C",
            "fresh-inputs",
            1200,
            10**13,
            {"kind": "fresh_generation"},
            lambda: prepare(
                combined, band=selection["band"], output=output, seconds=990
            ),
        )


def freeze_confirmation(selection_path, output):
    """Bind committed fresh inputs and both unchanged arms."""
    selection = json.loads(Path(selection_path).read_text())
    require_clean()
    parent = Path(selection["parent_output"])
    if digest(parent / "ledger.json") != selection["ledger_sha256"]:
        raise ValueError("exposed charges changed after selection")
    current = source_manifest()
    old = selection["source"]
    if any(
        current["files"].get(name) != sha for name, sha in old["files"].items()
    ):
        raise ValueError(
            "confirmation cannot silently change either underlying arm"
        )
    if (current["interpreter_sha256"], current["version"]) != (
        old["interpreter_sha256"],
        old["version"],
    ):
        raise ValueError("confirmation interpreter changed")
    band = selection["band"]
    corpus_path = INPUTS / ("confirmation_" + str(band)) / "corpus.json"
    corpus = json.loads(corpus_path.read_text())
    validate(corpus, band=band)
    resolved = json.loads(
        corpus_path.with_name("structural.json").read_text()
    )["resolved"]
    plan = dict(
        identity="c9a-full-yield-v2-fresh-" + str(band),
        cells=confirmation_cells(corpus, resolved, band),
        confirmation_band=band,
        parent_selection=selection,
    )
    resources = limits(plan, confirmation_band=band)
    prep = json.loads(
        (parent / "confirmation-preparation/ledger.json").read_text()
    )
    if any(row["status"] != "complete" for row in prep["calls"]):
        raise RuntimeError("fresh preparation has unsettled or failed work")
    preparation_charges = {
        key: sum(row[key] for row in prep["calls"])
        for key in ("wall", "cpu", "work")
    }
    for key, charge in preparation_charges.items():
        resources["C"][key] -= charge
    resources["total"] = {
        key: resources["C"][key] for key in ("wall", "cpu", "work")
    }
    resources["disk_bytes"] = resources["C"]["disk_bytes"]
    data = dict(
        plan=plan,
        limits=resources,
        source=current,
        corpus_path=str(corpus_path.relative_to(ROOT)),
        corpus_sha256=digest(corpus_path),
        confirmation_band=band,
        preparation_charges=preparation_charges,
    )
    data["sha256"] = hashlib.sha256(
        json.dumps(data, sort_keys=True).encode()
    ).hexdigest()
    if Path(output).exists():
        raise ValueError("fresh freeze already exists")
    atomic(Path(output), data)
    return dict(freeze=data["sha256"], band=band)
