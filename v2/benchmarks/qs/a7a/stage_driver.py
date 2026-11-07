"""Small gated batch driver: immutable proposals, receipts and no successor."""

import hashlib
import json
import math
import resource
import time
from pathlib import Path

from .combined import corpus_identity, fixtures, run
from .confirmation import BANDS
from .runner import source, write

ROWS = {(d, "cold") for d in BANDS} | {(40, "warm"), (50, "warm")}
CELL_FIELDS = frozenset(
    (
        "arm checkpoint_bytes collector_overrides cpu cpu_window "
        "diagnostic_for endpoint features fixture full_profile id lifetime "
        "max_batch_atoms memory nursery_mb observe packet_counters pair "
        "policy polynomials polynomials_per_task profile_reference_seconds "
        "profile_required progress_path qs_memory_bytes rss_interval "
        "samples seed siqs_overrides storage tasks tier_a timeline_path "
        "wall warmup_seconds work digits mode stage variant purpose "
        "worker_observe "
    ).split()
)


def save(path, value):
    """Atomically replace mutable accounting, preserving exclusive freezes."""
    path = Path(path)
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


def digest(value):
    return hashlib.sha256(
        json.dumps(value, sort_keys=True).encode()
    ).hexdigest()


def validate(plan, approval):
    """Reject incomplete coverage and unfunded work before opening a window."""
    if plan["stage"] not in ("G2", "G2H", "G3", "G4"):
        raise ValueError("stage is not admitted")
    if approval != {
        "stage": plan["stage"],
        "plan_sha256": digest(plan),
        "approved": True,
    }:
        raise ValueError("explicit approval must bind this exact stage plan")
    if plan["source"] != source():
        raise ValueError("implementation source changed")
    if plan["corpus_sha256"] != corpus_identity(plan["corpus_kind"]):
        raise ValueError("corpus source changed")
    for grant in [plan] + plan["batches"]:
        for name in ("wall", "cpu"):
            value = grant[name]
            if (
                type(value) not in (int, float)
                or not math.isfinite(value)
                or value <= 0
            ):
                raise ValueError("grant must be positive and finite")
    rows = {(50, "cold"), (65, "cold")} if plan["stage"] == "G2H" else ROWS
    panel = fixtures(plan["corpus_kind"])
    present, ids, recipes = set(), set(), {}
    wall = cpu = 0
    for batch in plan["batches"]:
        if batch["identity"] in ids or not batch["cells"]:
            raise ValueError("duplicate or empty batch")
        ids.add(batch["identity"])
        wall += batch["wall"]
        cpu += batch["cpu"]
        if batch["source"] != plan["source"] or (
            batch["corpus_sha256"] != plan["corpus_sha256"]
            or batch["corpus_kind"] != plan["corpus_kind"]
        ):
            raise ValueError("batch source differs from stage")
        reserved_wall = reserved_cpu = 0
        for cell in batch["cells"]:
            if cell.keys() - CELL_FIELDS:
                raise ValueError("unsupported cell fields")
            if any(
                type(cell[k]) is not int or cell[k] <= 0
                for k in ("memory", "work")
            ):
                raise ValueError("integer work and memory grants required")
            row = (cell["digits"], cell["lifetime"])
            if len(str(panel[cell["fixture"]]["n"])) != cell["digits"]:
                raise ValueError("fixture differs from declared regime")
            diagnostic = cell.get("purpose") == "diagnostic"
            if row not in rows or cell["endpoint"] not in (
                ("complete", "prefix") if diagnostic else ("complete",)
            ):
                raise ValueError("missing complete endpoint or invalid regime")
            if not diagnostic:
                present.add(row)
            if cell["lifetime"] == "cold" and cell["samples"] != 1:
                raise ValueError("cold observations require fresh actors")
            if (
                min(cell["wall"], cell["cpu"], cell["memory"], cell["work"])
                <= 0
            ):
                raise ValueError("nonpositive allowance")
            if not all(math.isfinite(cell[k]) for k in ("wall", "cpu")):
                raise ValueError("cell grant must be finite")
            reserved_wall += cell["wall"] + 5
            reserved_cpu += cell["cpu"]
            # These fields identify execution, rather than an input/repetition.
            excluded = {
                "fixture",
                "seed",
                "digits",
                "samples",
                "pair",
                "diagnostic_for",
                "variant",
                "purpose",
                "endpoint",
                "progress_path",
                "timeline_path",
                "tier_a",
                "observe",
                "rss_interval",
                "full_profile",
                "profile_required",
                "profile_reference_seconds",
                "cpu_window",
                "packet_counters",
                "mode",
                "id",
                "stage",
            }
            recipe = {k: v for k, v in cell.items() if k not in excluded}
            recipes.setdefault(row, set()).add(digest(recipe))
        if reserved_wall > batch["wall"] or reserved_cpu > batch["cpu"]:
            raise ValueError("whole batch reservation exceeds allowance")
    if present != rows:
        raise ValueError("coverage barrier has missing rows")
    if wall > plan["wall"] or cpu > plan["cpu"]:
        raise ValueError("whole stage reservation exceeds allowance")
    for row, choices in recipes.items():
        cap = plan["recipe_caps"][f"{row[0]}:{row[1]}"]
        maximum = 28 if plan["stage"] == "G2H" else 24
        if not 1 <= cap <= maximum or len(choices) > cap:
            raise ValueError("configuration count exceeds stage cap")
    return recipes


def execute(plan, approval, output):
    """Resume untouched batches; interrupted attempts need review."""
    validate(plan, approval)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    freeze = output / "plan.json"
    if freeze.exists():
        if json.loads(freeze.read_text()) != plan:
            raise ValueError("never overwrite a historical freeze")
    else:
        write(freeze, plan)
    ledger_path = output / "stage.json"
    ledger = (
        json.loads(ledger_path.read_text()) if ledger_path.exists() else {}
    )
    for batch in plan["batches"]:
        identity = batch["identity"]
        if identity in ledger:
            if ledger[identity]["status"] == "reported":
                continue
            raise ValueError("unfinished attempt requires reconciliation")
        if (
            sum(r.get("wall", 0) for r in ledger.values()) + batch["wall"]
            > plan["wall"]
            or sum(r.get("cpu", 0) for r in ledger.values()) + batch["cpu"]
            > plan["cpu"]
        ):
            raise RuntimeError("remaining stage grant cannot fund next batch")
        started, cpu_start = time.monotonic(), time.process_time()
        before = resource.getrusage(resource.RUSAGE_CHILDREN)
        ledger[identity] = {"status": "issued"}
        save(ledger_path, ledger)
        folder = output / identity
        folder.mkdir(exist_ok=False)
        frozen = folder / "plan.json"
        write(frozen, batch)
        try:
            if batch.get("reference_source"):
                from .mainline_reference import run as reference_run

                reference_run(frozen, folder / "receipts")
            else:
                run(frozen, folder / "receipts")
        finally:
            after = resource.getrusage(resource.RUSAGE_CHILDREN)
            ledger[identity].update(
                wall=time.monotonic() - started,
                cpu=time.process_time()
                - cpu_start
                + after.ru_utime
                + after.ru_stime
                - before.ru_utime
                - before.ru_stime,
            )
            receipt = folder / "receipts/ledger.json"
            if receipt.exists():
                measured = json.loads(receipt.read_text())
                ledger[identity].update(
                    wall=max(measured["wall"], ledger[identity]["wall"]),
                    cpu=max(measured["cpu"], ledger[identity]["cpu"]),
                )
                attempts = measured["attempts"]
                if len(attempts) == len(batch["cells"]) and all(
                    r["status"] == "reported" for r in attempts
                ):
                    ledger[identity]["status"] = "reported"
            save(ledger_path, ledger)
        if ledger[identity]["status"] != "reported":
            raise RuntimeError("batch incomplete; stop at this decision")
    return ledger
