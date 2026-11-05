"""Prospective v3 amendment replacing twelve cells with one NEON screen."""

import argparse
import copy
import json
from pathlib import Path

from .neon_screen import CELL_ID, screen_cell
from .runner import ROOT, sha, validate_freeze, write_new


def amend(base, admission):
    """Prepare a new identity, preserving the source freeze and old ledger."""
    correct_stage = admission.get("stage") == "NEON-correctness-admission" or (
        admission.get("stage") == "Q" and admission.get("neon_passed") is True
    )
    if (
        not correct_stage
        or admission.get("passed") is not True
        or admission.get("commit") != base["commit"]
        or admission.get("timing_executed") is not False
    ):
        raise ValueError("exact successful untimed NEON admission required")
    result = copy.deepcopy(base)
    removed = [
        cell
        for cell in result["cells"]
        if cell["id"].startswith("schedule-")
        and cell.get("capacity") == int(cell["arm"][-1])
    ]
    if len(removed) != 12 or any(
        cell["id"] == CELL_ID for cell in result["cells"]
    ):
        raise ValueError("v3 expects the preserved 93-cell OpenMP v2 base")
    identities = {cell["id"] for cell in removed}
    result["cells"] = [
        cell for cell in result["cells"] if cell["id"] not in identities
    ] + [screen_cell()]
    result.update(
        identity="b4a-A-neon-v3/" + base["commit"],
        supersedes_identity=base["identity"],
        neon_admission=admission,
        removed_cells=[cell["id"] for cell in removed],
    )
    if (
        len(result["cells"]) != 82
        or sum(c["wall"] for c in result["cells"]) != 1695
        or sum(c["cpu"] for c in result["cells"]) != 6420
    ):
        raise ValueError("v3 cumulative allocation changed")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-freeze", type=Path, required=True)
    parser.add_argument("--admission", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    base = json.loads(args.base_freeze.read_text())
    validate_freeze(base, ROOT)
    admission = json.loads(args.admission.read_text())
    result = amend(base, admission)
    for name in (
        "manifest.json",
        "libfactor_ecm_neon.dylib",
        "disassembly.txt",
    ):
        path = "v2/ecm/native/build/neon/" + name
        result["sha256"][path] = sha(ROOT / path)
    result["neon_admission_sha256"] = sha(args.admission)
    result["base_freeze_sha256"] = sha(args.base_freeze)
    validate_freeze(result, ROOT)
    write_new(args.output, result)


if __name__ == "__main__":
    main()
