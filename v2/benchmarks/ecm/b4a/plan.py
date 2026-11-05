"""Finite prospective cells and calibrated admission rules; no timings."""

import argparse
import hashlib
import json
import math
import random
import subprocess
from pathlib import Path

from .runner import ARMS, ROOT, canonical, check_cell, sha, write_new

INPUT = ROOT / "v2/benchmarks/inputs/controls/b4a_preparation.json"


def uniform_panel(seed, count):
    """Unconditioned exact-digit integers, including primes and cheap exits."""
    generator = random.Random(seed)
    return [
        dict(
            id=f"u-{seed}-{digits}-{index}",
            n=generator.randrange(10 ** (digits - 1), 10**digits),
            digits=digits,
            population="uniform",
        )
        for digits in (60, 80, 100)
        for index in range(count)
    ]


def cell(
    identity,
    arm,
    fixture,
    config,
    *,
    wall,
    mode="cold",
    endpoint="fixed",
    curves=8,
    observer="none",
    schedule="dynamic",
    capacity=None,
):
    work = (
        32_000_000 * curves + 2_000_000
        if endpoint == "fixed"
        else 2_000_000_000
    )
    if mode == "warm":
        work *= 4097  # Finite 4096-call warmup plus one measured call.
    result = dict(
        id=identity,
        arm=arm,
        input=fixture,
        config=config,
        mode=mode,
        endpoint=endpoint,
        curves=curves,
        wall=wall,
        cpu=4 * wall,
        shutdown=2,
        rss=4 * 2**30,
        output=8 * 2**20,
        work=work,
        assignment_work=work,
        assignment_cpu=4 * wall,
        processes=6,
        seed=7919,
        observer=observer,
        literal_cli=(
            endpoint == "complete"
            and arm in ("S0", "Native-S0")
            and config == "automatic"
        ),
    )
    if arm.startswith("Native-T"):
        result["schedule"] = schedule
        if capacity is not None:
            result["capacity"] = capacity
    if mode == "warm":
        result.update(
            warmup_wall=8,
            warmup_cpu=32,
            measurement_wall=wall - 10,
            measurement_cpu=4 * (wall - 10),
            measurement_work=work // 4097,
        )
    check_cell(result)
    return result


def stage_a():
    cells = []
    # Current supported full-tier eligibility: 60-digit/2000 uses accepted
    # reused PRAC. 100-digit/2000 and 11000 retain the supported ladders.
    for width in (60, 100):
        for bound in (2000, 11000):
            for arm in ("M0", "N1"):
                cells.append(
                    cell(
                        f"scalar-{width}-{bound}-{arm}",
                        arm,
                        f"width{width}",
                        f"b{bound}",
                        wall=20,
                    )
                )
    for arm in ARMS:
        cells.append(
            cell(
                "cheap-" + arm,
                arm,
                "cheap30",
                "automatic",
                wall=15,
                endpoint="complete",
                curves=0,
            )
        )
    for digits in (60, 80, 100):
        for index in range(2):
            for arm in ("S0", "Native-S0"):
                cells.append(
                    cell(
                        f"uniform-{digits}-{index}-{arm}",
                        arm,
                        f"uniform{digits}-{index}",
                        "automatic",
                        wall=15,
                        endpoint="complete",
                        curves=0,
                    )
                )
    for fixture in ("width60", "width80"):
        for arm in ("S0", "Native-S0"):
            cells.append(
                cell(
                    f"useful-{fixture}-{arm}",
                    arm,
                    fixture,
                    "automatic",
                    wall=20,
                    endpoint="complete",
                    curves=0,
                )
            )
    for arm in (
        "B12-P2",
        "B12-P4",
        "Native-P2",
        "Native-P4",
        "Native-T1",
        "Native-T2",
        "Native-T4",
    ):
        cells.append(
            cell(
                "startup-" + arm,
                arm,
                "cheap30",
                "automatic",
                wall=15,
                endpoint="startup",
                curves=0,
            )
        )
    for arm in (
        "M0",
        "N1",
        "Native-P2",
        "Native-P4",
        "Native-T2",
        "Native-T4",
    ):
        cells.append(
            cell("warm-" + arm, arm, "width60", "b2000", wall=35, mode="warm")
        )
    for arm in ("B12-P2", "B12-P4", "Native-P2", "Native-P4"):
        cells.append(cell("process-" + arm, arm, "width60", "b11000", wall=25))
    for arm in ("M0", "N1"):
        for observer in ("none", "cprofile"):
            cells.append(
                cell(
                    f"observer-{arm}-{observer}",
                    arm,
                    "width60",
                    "b2000",
                    wall=20,
                    observer=observer,
                )
            )
    for digits in (39, 41, 79, 81):
        for arm in ("S0", "Native-S0"):
            cells.append(
                cell(
                    f"boundary-{digits}-{arm}",
                    arm,
                    f"cheap{digits}",
                    "automatic",
                    wall=10,
                    endpoint="complete",
                    curves=0,
                )
            )
    for index in range(2):
        cells.append(cell(f"aa-{index}", "M0", "width60", "b2000", wall=20))
    for arm in ("Native-T1", "Native-T2", "Native-T4"):
        cells.append(cell("thread-" + arm, arm, "width60", "b11000", wall=25))
    # Small selected schedule/queue pilot, not a whole-machine worker sweep.
    # Full-modulus PRAC and ladder widths use the same eight-curve manifest.
    for width in (60, 100):
        for workers in (2, 4):
            for capacity in (workers, 2 * workers):
                for schedule in ("dynamic", "static", "guided"):
                    identity = f"schedule-{width}-{workers}-{capacity}"
                    cells.append(
                        cell(
                            identity + "-" + schedule,
                            f"Native-T{workers}",
                            f"width{width}",
                            "b2000",
                            wall=20,
                            schedule=schedule,
                            capacity=capacity,
                        )
                    )
    if len(cells) != 93:
        raise AssertionError("Stage A allocation changed")
    return cells


def downstream_reservation(stage, observations, *, cells, curves):
    """Require informative clean A rows, then freeze literal downstream caps.

    This computes a proposal; it never grants execution. Removing regimes or
    arms must occur before fresh outcomes. Never reduce required matched
    counts, nine samples or independent-input requirements to fit a budget.
    """
    ceilings = {
        "B": (3600, 14400, 20),
        "C": (7200, 28800, 120),
        "D": (18000, 72000, 120),
    }
    wall_max, cpu_max, per_cell_max = ceilings[stage]
    clean = [
        r
        for r in observations
        if r.get("status") == "complete" and not r.get("instrumented")
    ]
    if len(clean) < 3 or any(r.get("failure") for r in observations):
        raise ValueError(
            "diagnose missing/failed/uninformative calibration before freezing"
        )
    largest = max(r["wall"] for r in clean)
    cap = max(12 if stage == "B" else 3, math.ceil(2 * largest + 3))
    if (
        cap > per_cell_max
        or cells * cap > wall_max
        or cells * 4 * cap > cpu_max
    ):
        raise ValueError(
            "proposed cells do not fit; narrow scope before execution or stop"
        )
    rates = [r["work"] / r["wall"] for r in clean if r["wall"] > 0]
    if not rates or max(rates) <= 0:
        raise ValueError("logical-work rate calibration required")
    # The analytical floor includes B12's retained 32M/curve grant and native
    # setup/copy work. The rate term prevents work truncating the time grant.
    work = max(
        32_000_000 * curves + 2_000_000, math.ceil(2 * max(rates) * cap)
    )
    if work >= 2**63:
        raise ValueError("work cap exceeds checked worker integer domain")
    return dict(
        wall=cap,
        cpu=4 * cap,
        work=work,
        assignment_work=work,
        assignment_cpu=4 * cap,
        cells=cells,
        reserved_wall=cells * cap,
        reserved_cpu=cells * 4 * cap,
        wall_ceiling=wall_max,
        cpu_ceiling=cpu_max,
        shutdown=2,
        rss=4 * 2**30,
        output=8 * 2**20,
        processes=6,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--qa", type=Path, required=True)
    parser.add_argument("--integration", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    qa = json.loads(args.qa.read_text())
    integration = json.loads(args.integration.read_text())
    commit = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
    ).strip()
    if (
        not qa.get("passed")
        or qa.get("stage") != "Q"
        or qa.get("commit") != commit
    ):
        raise ValueError(
            "Stage Q on this exact combined committed source is required"
        )
    if (
        integration.get("reviewed") is not True
        or integration.get("commit") != commit
    ):
        raise ValueError(
            "explicit settled B12/mainline integration receipt required"
        )
    plan = json.loads(INPUT.read_text())
    if plan["status"] != "prospective-not-executable":
        raise ValueError("historical preparation identity changed")
    cells = stage_a()
    from ....execution.ecm_workers import arithmetic_options
    from ....execution.portfolio_config import PortfolioConfig
    from ....execution.worker_contract import Manifest, WorkerSettings

    for selected in cells:
        actual = PortfolioConfig(
            **dict(
                plan["configurations"][selected["config"]],
                ecm_kernel=ARMS[selected["arm"]][0],
            )
        ).for_input(plan["inputs"][selected["input"]]["n"])
        if actual.memory_bytes > plan["settings"]["parent_memory"]:
            raise ValueError(
                "parent reserve does not include the automatic portfolio"
            )
        if selected["endpoint"] != "fixed":
            continue
        kernel, workers = ARMS[selected["arm"]]
        config = PortfolioConfig(
            **dict(
                plan["configurations"][selected["config"]], ecm_kernel=kernel
            )
        )
        settings = WorkerSettings(
            **dict(plan["settings"], workers=workers or 1)
        )
        manifest = Manifest(
            plan["inputs"][selected["input"]]["n"],
            selected["seed"],
            "root",
            arithmetic_options(config, settings),
        )
        selected["expected_assignments"] = {
            str(i): manifest.assignment(i) for i in range(selected["curves"])
        }
    paths = subprocess.check_output(
        ["git", "ls-files", "v2"], cwd=ROOT, text=True
    ).splitlines()
    paths += [
        "v2/ecm/native/build/manifest.json",
        "v2/ecm/native/build/libfactor_ecm.dylib",
        "v2/ecm/native/build/openmp/manifest.json",
        "v2/ecm/native/build/openmp/libfactor_ecm.dylib",
    ]
    freeze = dict(
        identity="b4a-A-openmp-v2/" + commit,
        status="execution-admitted",
        commit=commit,
        qa=qa,
        integration=integration,
        pypy=qa["pypy"],
        sha256={p: sha(ROOT / p) for p in paths if (ROOT / p).is_file()},
        cells=cells,
        inputs=plan["inputs"],
        configurations=plan["configurations"],
        settings=plan["settings"],
        wall=1800,
        cpu=7200,
        preparation_sha256=sha(INPUT),
        qa_sha256=sha(args.qa),
        integration_sha256=sha(args.integration),
    )
    freeze["prospective_identity"] = hashlib.sha256(
        canonical(plan).encode()
    ).hexdigest()
    write_new(args.output, freeze)


if __name__ == "__main__":
    main()
