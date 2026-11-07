"""Executable cell schedule; no imports or factoring on plan construction."""

import json
from pathlib import Path

IDENTITY = "a7a-current-workers-phase3-v18"
PACKET_CONTROL = "f38520e85398a60b115a8a34f9b7b9685ec998b6"
STATISTICS_CONTROL = "a90dfa6e4e8f46cfe7e19113d6a0fc8e6337130b"
CPU_READ_CONTROL = "7ff5135ebc16c0648df67a1bb2dfd9397f2dbb9c"
JOINT_CONTROL = "6baa1508f278d1b2bd69635d3f78ab3155e22d9b"
ROOT = Path(__file__).resolve().parents[4]
INPUTS = ROOT / "v2/benchmarks/inputs/a7a"
LIMITS = {
    "Q": (1800, 7200),
    "H0": (600, 6000),
    "H1": (600, 6000),
    "H2": (600, 6000),
    "H3": (600, 6000),
    "H4": (750, 7500),
    "H5": (900, 9000),
    "H6": (900, 9000),
    "H7": (1050, 10500),
    "H8": (1050, 10500),
    "H9": (660, 6600),
    "H10": (480, 4800),
    "H11": (1800, 18000),
    "H12": (1500, 15000),
    "H13": (2100, 21000),
    "H14": (1500, 15000),
    "J0": (900, 9000),
    "A": (3600, 14400),
    "B": (3600, 14400),
    "C": (18000, 72000),
    "D": (14400, 57600),
    "E": (108000, 432000),
}


def cell(
    stage,
    fixture,
    arm,
    *,
    wall,
    endpoint="prefix",
    lifetime="cold",
    samples=1,
    tasks=4,
    observe=False,
    policy="Y0",
    seed=104729,
):
    return dict(
        stage=stage,
        fixture=fixture,
        arm=arm,
        wall=wall,
        cpu=4 * wall,
        work=32 * 10**9 * 4 * wall,
        memory=16 * 2**30,
        storage=8 * 2**20,
        endpoint=endpoint,
        lifetime=lifetime,
        samples=samples,
        tasks=tasks,
        observe=observe,
        policy=policy,
        seed=seed,
    )


def schedule(stage, nominee="P2"):
    """Freeze limited questions without a Cartesian parameter sweep."""
    rows = []
    if stage == "H0":
        # One exposed modulus per size is a feasibility screen, not a claim.
        # Reverse the endpoint arm order and interleave the larger-work arms.
        for size, arms in ((30, ("S0", "P1")), (50, ("P1", "S0"))):
            fixture = "engage_30_0" if size == 30 else "exposed_50_0"
            for arm in arms:
                rows.append(
                    cell(
                        stage,
                        fixture,
                        arm,
                        wall=15,
                        endpoint="complete",
                        policy="G0",
                    )
                )
        for fixture, arms in (
            ("exposed_60_0", ("S0", "P4", "M0", "P2", "P1")),
            ("exposed_70_0", ("P1", "P2", "M0", "P4", "S0")),
        ):
            for arm in arms:
                rows.append(cell(stage, fixture, arm, wall=30))
        for row in rows:
            row["cpu"] = 10 * row["wall"]
            row["work"] = 32 * 10**9 * row["cpu"]
    elif stage == "H1":
        # H0 completed its small prefixes in about 1–2 seconds per arm.
        # Match 64-polynomial cold/warm work to isolate lifecycle cost, then
        # price the existing observer on exactly that same warm workload.
        for fixture, arms in (
            ("exposed_60_0", ("S0", "P4", "M0", "P2", "P1")),
            ("exposed_70_0", ("P1", "P2", "M0", "P4", "S0")),
        ):
            for arm in arms:
                rows.append(
                    cell(
                        stage,
                        fixture,
                        arm,
                        wall=30,
                        tasks=16,
                        lifetime="warm",
                    )
                )
                if arm == "P4":
                    observed = cell(
                        stage,
                        fixture,
                        arm,
                        wall=30,
                        tasks=16,
                        lifetime="warm",
                        observe=True,
                    )
                    if fixture == "exposed_60_0":
                        rows.append(observed)
                    else:
                        rows.insert(len(rows) - 1, observed)
            cold_arms = (
                ("P4", "P1", "S0")
                if fixture == "exposed_60_0"
                else ("S0", "P1", "P4")
            )
            for arm in cold_arms:
                rows.append(cell(stage, fixture, arm, wall=30, tasks=16))
        for row in rows:
            row["cpu"] = 10 * row["wall"]
            row["work"] = 32 * 10**9 * row["cpu"]
            row["storage"] = 2**20
    elif stage == "H2":
        for fixture, variants in (
            (
                "exposed_60_0",
                (
                    ("S0", 4),
                    ("M0", 64),
                    ("P1", 4),
                    ("M0", 16),
                    ("P1", 64),
                    ("M0", 4),
                    ("P1", 16),
                ),
            ),
            (
                "exposed_70_0",
                (
                    ("P1", 16),
                    ("M0", 4),
                    ("P1", 64),
                    ("M0", 16),
                    ("P1", 4),
                    ("M0", 64),
                    ("S0", 4),
                ),
            ),
        ):
            for arm, width in variants:
                row = cell(stage, fixture, arm, wall=30, lifetime="warm")
                row.update(
                    polynomials=128,
                    polynomials_per_task=width,
                    cpu=300,
                    work=32 * 10**9 * 300,
                    storage=2**20,
                )
                rows.append(row)
    elif stage == "H3":
        for fixture, variants in (
            (
                "exposed_60_0",
                (
                    ("S0", "candidate"),
                    ("P4", "packet-control"),
                    ("P4", "candidate"),
                    ("P1", "candidate"),
                    ("P1", "packet-control"),
                ),
            ),
            (
                "exposed_70_0",
                (
                    ("P1", "packet-control"),
                    ("P1", "candidate"),
                    ("P4", "candidate"),
                    ("P4", "packet-control"),
                    ("S0", "candidate"),
                ),
            ),
        ):
            for arm, implementation in variants:
                row = cell(stage, fixture, arm, wall=45, lifetime="warm")
                row.update(
                    polynomials=128,
                    polynomials_per_task=4,
                    implementation=implementation,
                    cpu=450,
                    work=32 * 10**9 * 450,
                    storage=2**20,
                )
                rows.append(row)
    elif stage == "H4":
        # Four full-family packets still occupy four workers. Equal work is
        # essential: equal task counts would give larger packets more work.
        variants = (
            ("S0", 4),
            ("P4", 64),
            ("P1", 4),
            ("P2", 16),
            ("P4", 4),
            ("P1", 64),
            ("P2", 4),
            ("P4", 16),
            ("P1", 16),
            ("P2", 64),
        )
        for fixture, ordered in (
            ("exposed_60_0", variants),
            ("exposed_70_0", tuple(reversed(variants))),
        ):
            for arm, width in ordered:
                row = cell(stage, fixture, arm, wall=30, lifetime="warm")
                row.update(
                    polynomials=256,
                    polynomials_per_task=width,
                    cpu=300,
                    work=32 * 10**9 * 300,
                    storage=2**20,
                )
                rows.append(row)
    elif stage == "J0":
        for arm in ("M0", "S0", "S0", "M0"):
            row = cell(
                stage, "exposed_60_0", arm, wall=210, endpoint="complete"
            )
            row.update(
                cpu=2100,
                work=32 * 10**9 * 3000,
                checkpoint_bytes=64 * 2**20,
                qs_memory_bytes=4 * 2**30,
                memory=20 * 2**30,
                polynomials_per_task=64,
                max_batch_atoms=2048,
            )
            rows.append(row)
    elif stage == "H14":
        for guard in ("peer", "parent", "parent", "peer"):
            row = cell(
                stage, "exposed_60_0", "P10", wall=300, endpoint="complete"
            )
            row.update(
                cpu_guard=guard,
                cpu=3000,
                work=32 * 10**9 * 3000,
                checkpoint_bytes=64 * 2**20,
                qs_memory_bytes=4 * 2**30,
                memory=20 * 2**30,
                polynomials_per_task=64,
                max_batch_atoms=2048,
            )
            rows.append(row)
    elif stage == "H13":
        # Bracket the wider-count screen to expose sustained-run drift.
        for arm in ("S0", "P4", "P6", "P10", "P8", "P4"):
            row = cell(
                stage, "exposed_60_0", arm, wall=300, endpoint="complete"
            )
            row.update(
                cpu=3000,
                work=32 * 10**9 * 3000,
                checkpoint_bytes=64 * 2**20,
                qs_memory_bytes=4 * 2**30,
                memory=20 * 2**30,
                polynomials_per_task=64,
                max_batch_atoms=2048,
            )
            rows.append(row)
    elif stage == "H12":
        for arm, implementation in (
            ("S0", "budget-control"),
            ("S0", "candidate"),
            ("P4", "candidate"),
            ("P4", "budget-control"),
        ):
            row = cell(
                stage, "exposed_60_0", arm, wall=300, endpoint="complete"
            )
            row.update(
                implementation=implementation,
                cpu=3000,
                work=32 * 10**9 * 3000,
                checkpoint_bytes=64 * 2**20,
                qs_memory_bytes=12 * 2**30,
                memory=20 * 2**30,
                polynomials_per_task=64,
            )
            rows.append(row)
    elif stage == "H11":
        for index, guard in enumerate(
            ("peer", "parent", "parent", "peer", "peer")
        ):
            row = cell(
                stage, "exposed_60_0", "P4", wall=300, endpoint="complete"
            )
            row.update(
                cpu_guard=guard,
                cpu=3000,
                work=32 * 10**9 * 3000,
                checkpoint_bytes=64 * 2**20,
                qs_memory_bytes=12 * 2**30,
                memory=20 * 2**30,
                polynomials_per_task=64,
            )
            if index == 4:
                row.update(cpu_profile=[60, 75], storage=384 * 2**20)
            rows.append(row)
    elif stage == "H10":
        for fixture in ("exposed_60_0", "exposed_70_0"):
            for guard in ("peer", "parent", "parent", "peer"):
                row = cell(stage, fixture, "P4", wall=45)
                row.update(
                    cpu_guard=guard,
                    polynomials=768,
                    polynomials_per_task=64,
                    cpu=450,
                    work=32 * 10**9 * 450,
                )
                rows.append(row)
    elif stage == "H9":
        variants = ("cpu-read-control", "candidate", "joint-control")
        for fixture, variants in (
            ("exposed_60_0", variants),
            ("exposed_70_0", tuple(reversed(variants))),
        ):
            for implementation in variants:
                row = cell(
                    stage, fixture, "P4", wall=90, lifetime="warm", samples=9
                )
                row.update(
                    implementation=implementation,
                    polynomials=256,
                    polynomials_per_task=64,
                    cpu=900,
                    work=32 * 10**9 * 900,
                )
                rows.append(row)
    elif stage in ("H5", "H6", "H7", "H8"):
        # Complete calls expose mature graph/solve costs hidden by fixed-work
        # prefixes. Failure to finish one arm does not forbid the other arm.
        arms = ("S0", "P4", "S0") if stage == "H7" else ("S0", "P4")
        if stage == "H8":
            arms = ("S0", "P4", "P4")
        for index, arm in enumerate(arms):
            row = cell(
                stage,
                "exposed_60_0",
                arm,
                wall=300,
                endpoint="complete",
            )
            row.update(
                polynomials_per_task=64,
                checkpoint_bytes=64 * 2**20,
                cpu=3000,
                work=32 * 10**9 * 3000,
            )
            if stage in ("H6", "H7", "H8"):
                # The 64-MiB snapshot reserves eight coexisting copies. Fund
                # the collector after that known reserve, without retuning it.
                row["qs_memory_bytes"] = 12 * 2**30
                row["memory"] = 20 * 2**30
            if stage == "H7":
                row["implementation"] = (
                    "statistics-control" if index == 0 else "candidate"
                )
            if stage == "H8":
                row["implementation"] = (
                    "cpu-read-control" if index == 1 else "candidate"
                )
            rows.append(row)
    elif stage == "A":
        for size in (30, 40, 50):
            for index in range(2):
                fixture = (
                    f"engage_{size}_{index}"
                    if size < 50
                    else f"exposed_50_{index}"
                )
                for arm in ("S0", "P1"):
                    rows.append(
                        cell(
                            stage,
                            fixture,
                            arm,
                            wall=15,
                            endpoint="complete",
                            policy="G0",
                        )
                    )
        for fixture in (
            "exposed_50_0",
            "exposed_60_0",
            "exposed_60_1",
            "exposed_70_0",
        ):
            for arm in ("S0", "M0", "P1", "P2", "P4"):
                rows.append(cell(stage, fixture, arm, wall=60))
        for fixture in ("exposed_60_0", "exposed_70_0"):
            for arm in ("S0", "M0", "P1", "P2", "P4"):
                rows.append(
                    cell(stage, fixture, arm, wall=60, lifetime="warm")
                )
        for arm in ("P1", "P2", "P4"):
            for observe in (False, True):
                for repeat in range(2):
                    rows.append(
                        cell(
                            stage,
                            "exposed_60_0",
                            arm,
                            wall=30,
                            tasks=4,
                            observe=observe,
                            seed=104729,
                        )
                    )
        # Three cold/warm pool/context probes, with matching current G0.
        for arm in ("P1", "P2", "P4"):
            for lifetime in ("cold", "warm"):
                rows.append(
                    cell(
                        stage,
                        "exposed_60_0",
                        arm,
                        wall=15,
                        tasks=1,
                        lifetime=lifetime,
                        policy="G0",
                    )
                )
    elif stage == "B":
        for fixture in ("exposed_60_0", "c10a01_60_balanced", "exposed_70_0"):
            for arm in ("P1", "P2", "P4"):
                for repeat in range(3):
                    rows.append(
                        cell(
                            stage,
                            fixture,
                            arm,
                            wall=60,
                            lifetime="warm",
                            seed=104729,
                        )
                    )
                rows.append(
                    cell(
                        stage,
                        fixture,
                        arm,
                        wall=60,
                        lifetime="warm",
                        observe=True,
                    )
                )
            for repeat in range(3):
                rows.append(
                    cell(
                        stage,
                        fixture,
                        "S0",
                        wall=60,
                        lifetime="warm",
                        seed=104729,
                    )
                )
        for arm in ("S0", nominee):
            rows.append(
                cell(stage, "exposed_60_0", arm, wall=360, endpoint="complete")
            )
        for size in (55, 62, 65):
            for arm in ("S0", nominee):
                rows.append(cell(stage, f"exposed_{size}_0", arm, wall=20))
    elif stage == "C":
        for size, wall in ((60, 600), (70, 3600)):
            for fixture in (f"exposed_{size}_0", f"c10a01_{size}_balanced"):
                for arm in ("S0", nominee):
                    rows.append(
                        cell(
                            stage, fixture, arm, wall=wall, endpoint="complete"
                        )
                    )
        for arm in ("S0", nominee):
            rows.append(
                cell(stage, "exposed_60_1", arm, wall=120, endpoint="first")
            )
    elif stage == "D":
        for size in (75, 80):
            for index in range(2):
                for arm in ("S0", nominee):
                    rows.append(
                        cell(stage, f"exposed_{size}_{index}", arm, wall=120)
                    )
        for size in (60, 80, 100):
            for index in range(12):
                for arm in ("S0", nominee):
                    rows.append(
                        cell(
                            stage,
                            f"random_{size}_{index:02}",
                            arm,
                            wall=60,
                            endpoint="portfolio",
                            policy="G0",
                        )
                    )
        for fixture in (
            "exposed_55_0",
            "exposed_65_0",
            "derived_60_12",
            "derived_70_21",
            "recursive_60",
            "repeated_60",
        ):
            for arm in ("S0", nominee):
                rows.append(
                    cell(stage, fixture, arm, wall=90, endpoint="complete")
                )
        # One frontier candidate pilot requires a separate prospectively pinned
        # admission and size choice. This reservation is never automatic.
    elif stage == "E":
        for size, wall in ((60, 600), (70, 3600)):
            for index in range(9):
                seeds = (104729, 130363) if index < 3 else (104729,)
                for seed in seeds:
                    for arm in ("S0", nominee):
                        rows.append(
                            cell(
                                stage,
                                f"fresh_{size}_{index:02}",
                                arm,
                                wall=wall,
                                endpoint="complete",
                                seed=seed,
                            )
                        )
    elif stage != "Q":
        raise ValueError("unknown stage")
    for index, row in enumerate(rows):
        row["id"] = f"{stage}-{index:03}"
    wall, cpu = LIMITS[stage]
    if (
        sum(r["wall"] for r in rows) > wall
        or sum(r["cpu"] for r in rows) > cpu
    ):
        raise ValueError("scheduled reservations exceed stage allowance")
    return dict(
        identity=IDENTITY,
        stage=stage,
        cells=rows,
        wall=wall,
        cpu=cpu,
        reserved_wall=sum(r["wall"] for r in rows),
        reserved_cpu=sum(r["cpu"] for r in rows),
        reserve_wall=wall - sum(r["wall"] for r in rows),
        reserve_cpu=cpu - sum(r["cpu"] for r in rows),
        machine_memory=(
            20
            if stage in ("H6", "H7", "H8", "H11", "H12", "H13", "H14", "J0")
            else 16
        )
        * 2**30,
        storage=32 * 2**20 if stage in ("H1", "H2", "H3", "H4") else 2**30,
    )


def corpus():
    return json.loads((INPUTS / "corpus.json").read_text())


def confirmation_tranche(plan, tranche):
    """Do not automatically spend a day on fresh confirmation."""
    if plan["stage"] != "E" or tranche not in ("E0", "E1"):
        raise ValueError("confirmation needs an E0/E1 tranche")
    early = {"fresh_60_00", "fresh_60_01", "fresh_60_02", "fresh_70_00"}
    cells = [
        row
        for row in plan["cells"]
        if (
            (row["fixture"] in early and row["seed"] == 104729)
            == (tranche == "E0")
        )
    ]
    wall = 12000 if tranche == "E0" else 96000
    reserved = sum(row["wall"] for row in cells)
    return dict(
        plan,
        cells=cells,
        tranche=tranche,
        wall=wall,
        cpu=4 * wall,
        reserved_wall=reserved,
        reserved_cpu=4 * reserved,
        reserve_wall=wall - reserved,
        reserve_cpu=4 * (wall - reserved),
    )
