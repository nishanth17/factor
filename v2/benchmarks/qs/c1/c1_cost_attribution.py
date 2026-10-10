"""Post-observation retained-data cost audit; never accepted performance."""

import argparse
import hashlib
import json
import platform
import subprocess
import time
from pathlib import Path

from ....execution.budget import Budget, BudgetExhaustedError
from ...support.paths import (
    source_path,
)
from .c1_feasibility import HERE, machine_window
from .c1_followup import (
    AnalysisBudget,
    compact_report,
    load_followup,
    save_capture,
)

PROTOCOL = HERE / "c1_cost_attribution.md"


def audit_cell(report, records, control, *, parent=None):
    """Seek a bounded direct witness without collecting another position."""
    allowance_budget = AnalysisBudget(parent)
    budget = allowance_budget.local
    snapshots = report["snapshots"]
    witnessed = [
        p
        for p in snapshots
        if p["reports"]["128"]["complete"]
        and not p["reports"]["slp"]["complete"]
        and not p["stopped"]
    ]
    if witnessed:
        selected = min(witnessed, key=lambda p: p["blocks"])
        lower, upper, enclosing = None, selected["blocks"], selected
    else:
        lower, upper, enclosing = 0, None, None
        for prefix in snapshots:
            candidate = prefix["reports"]["128"]
            if "censored" in candidate or candidate["complete"]:
                upper, enclosing = prefix["blocks"], prefix
                break
            lower = prefix["blocks"]
    result = dict(n=report["n"], observations=[], witnesses=[])
    if upper is None:
        result["reason"] = "no_complete_or_censored_bracket"
        return result
    costs = enclosing["costs"]
    split_cpu = costs.get("split_128", 0) + costs.get("classify_128", 0)
    allowance = max(control["cap_seconds"], control["cpu_seconds"]) / 2
    try:
        for _ in range(1 if lower is None else 8):
            block = upper if lower is None else (lower + upper) // 2
            retained = [r for r in records if r["block"] <= block]
            candidate = compact_report(
                retained,
                report["n"],
                budget=allowance_budget,
                first_factor=True,
            )
            slp = compact_report(
                [r for r in retained if r["kind"] != "dlp"],
                report["n"],
                budget=allowance_budget,
                first_factor=True,
            )
            charged = candidate["cpu_seconds"] + split_cpu
            observation = dict(
                blocks=block,
                candidate=candidate,
                slp=slp,
                split_certification_cpu_upper_bound=split_cpu,
                candidate_processing_cpu=candidate["cpu_seconds"],
                candidate_total_cpu_upper_bound=charged,
                cost_allowance=allowance,
                enclosing_blocks=enclosing.get(
                    "enclosing_blocks", enclosing["blocks"]
                ),
            )
            result["observations"].append(observation)
            if "censored" in candidate or "censored" in slp:
                result["reason"] = "analysis_censored"
                break
            if (
                candidate["complete"]
                and candidate["proper_divisors"]
                and candidate["independent_lp_constraints"]
                - slp["independent_lp_constraints"]
                >= 16
                and charged <= allowance
                and (
                    not report["stopped"] or block < report["counts"]["blocks"]
                )
            ):
                result["witnesses"].append(
                    dict(
                        blocks=block,
                        slp_has_factor=slp["complete"],
                        charged_cpu=charged,
                    )
                )
            if lower is None or upper - lower <= 1:
                break
            if candidate["complete"]:
                upper = block
            else:
                lower = block
    except BudgetExhaustedError:
        result["reason"] = budget.reason
    result["search_cpu"] = budget.cpu_used
    result["search_wall"] = budget.wall_used
    result["search_work"] = budget.used
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if (
        platform.python_implementation(),
        platform.python_version_tuple()[:2],
    ) != ("PyPy", ("3", "11")):
        raise RuntimeError("C1 requires PyPy implementing Python 3.11")
    prior = json.loads(
        (source_path(args.source / "decision.json")).read_text()
    )
    # Conservatively reserve the earlier retained analysis's entire 180
    # seconds as well as this audit; never extend the 4,500-second study.
    if max(prior["seconds"], prior["cpu_seconds"]) + 360 > 4500:
        raise RuntimeError("no remaining study allowance for attribution")
    _, fixtures, _ = load_followup()
    with machine_window():
        started, cpu = time.monotonic(), time.process_time()
        study = Budget(work_limit=6 * 10**12, seconds=180, cpu_seconds=180)
        result = dict(
            protocol_sha256=hashlib.sha256(
                source_path(PROTOCOL).read_bytes()
            ).hexdigest(),
            original_go_policies=prior["go_policies"],
            commit=subprocess.check_output(
                ["git", "rev-parse", "HEAD"], text=True
            ).strip(),
            driver_sha256=hashlib.sha256(
                source_path(Path(__file__)).read_bytes()
            ).hexdigest(),
            analyzer_sha256=hashlib.sha256(
                (source_path(HERE / "c1_followup.py")).read_bytes()
            ).hexdigest(),
            runtime=platform.python_version(),
            cells=[],
        )
        for fixture in fixtures:
            if (
                max(time.monotonic() - started, time.process_time() - cpu)
                > 150
            ):
                result["stopped"] = "audit_study_limit"
                break
            name = f"{fixture['digits']}-{fixture['c1_index']}"
            paths = {
                kind: args.source / f"{name}-{kind}.json"
                for kind in ("probe", "records", "control")
            }
            if not all(path.exists() for path in paths.values()):
                continue
            values = {
                kind: json.loads(source_path(path).read_text())
                for kind, path in paths.items()
            }
            cell = audit_cell(
                values["probe"],
                values["records"],
                values["control"],
                parent=study,
            )
            for observation in cell["observations"]:
                for arm in ("candidate", "slp"):
                    outcome = observation[arm]
                    for divisor in outcome.get("proper_divisors", ()):
                        assert (
                            sorted((divisor, fixture["n"] // divisor))
                            == fixture["factors"]
                        )
            cell["name"] = name
            cell["sources"] = {
                kind: hashlib.sha256(
                    source_path(path).read_bytes()
                ).hexdigest()
                for kind, path in paths.items()
            }
            result["cells"].append(cell)
            print(name, cell["witnesses"], cell.get("reason"), flush=True)
        result["go_bands"] = []
        for digits in (40, 50, 60):
            witnesses = [
                cell["witnesses"]
                for cell in result["cells"]
                if cell["name"].startswith(str(digits) + "-")
            ]
            if (
                len(witnesses) == 2
                and all(witnesses)
                and any(
                    not row["slp_has_factor"]
                    for rows in witnesses
                    for row in rows
                )
            ):
                result["go_bands"].append(digits)
        result.update(
            seconds=time.monotonic() - started,
            cpu_seconds=time.process_time() - cpu,
        )
        data = json.dumps(result)
        if len(data.encode()) > 16 * 2**20:
            raise ValueError("attribution capture byte cap")
        if (
            sum(p.stat().st_size for p in args.source.glob("*.json"))
            + len(data.encode())
            > 512 * 2**20
        ):
            raise ValueError("combined study capture byte cap")
        save_capture(args.output, result)


if __name__ == "__main__":
    main()
