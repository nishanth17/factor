"""Numerical finite cells; creating a plan never executes a campaign."""

from dataclasses import asdict

from .configuration import (
    COMPLETE_CAPS,
    HELPER_BYTES,
    IDENTITY,
    PREFIX_CAPS,
    SEEDS,
    WORK_CAPS,
    WORK_RATE_GUARD,
    configuration,
    owned_bytes,
    rss_bytes,
)

STAGES = {
    "A": dict(seconds=10200, warm=1800, disk=256 * 2**20),
    "B1": dict(seconds=9000, warm=300, disk=64 * 2**20),
    "M": dict(seconds=4000, warm=200, disk=640 * 2**20),
    "B2": dict(seconds=57000, warm=1500, disk=256 * 2**20),
    "R": dict(seconds=6000, warm=600, disk=128 * 2**20),
    "F": dict(seconds=29700, warm=600, disk=64 * 2**20),
    "C": dict(seconds=91440, warm=2400, disk=512 * 2**20),
}
CONFIRMATION_SECONDS = {60: 91440, 62: 120240, 65: 177840, 70: 350640}
PREFIX_BLOCKS = 512 * 33
WARM_WORK = 10**11
REPLAY_WORK = 10**13


def cell(row, resolved, stage, arm, *, kind="complete", seed=SEEDS[0]):
    digits = row["digits"]
    prefix = kind in ("prefix", "routing")
    if digits == 80 and not prefix:
        raise ValueError("no complete-factor 80-digit cell is authorized")
    seconds = (
        10
        if kind == "routing"
        else (PREFIX_CAPS[digits] if prefix else COMPLETE_CAPS[digits])
    )
    warm = 30 if stage == "A" else 20 if stage == "R" else 50
    config = configuration(digits, arm, admission_base_bound=resolved["bound"])
    return dict(
        fixture=row["id"],
        digits=digits,
        n=row["n"],
        uneven=row["uneven"],
        seed=seed,
        sample=SEEDS.index(seed),
        arm=arm,
        stage=stage,
        kind=kind,
        seconds=seconds,
        margin=5 if prefix else 15,
        warm_seconds=warm,
        end_reserve=10,
        warm_work_limit=WARM_WORK,
        work_limit=WORK_RATE_GUARD * seconds,
        observer_work_limit=2**32 if kind == "capture" else 0,
        observer_byte_limit=256 * 2**20 if kind == "capture" else 0,
        max_blocks=PREFIX_BLOCKS
        if kind == "prefix"
        else 33
        if kind == "routing"
        else None,
        owned_bytes=owned_bytes(digits),
        rss_bytes=rss_bytes(digits),
        helper_bytes=HELPER_BYTES,
        configuration=asdict(config),
        source="immutable" if arm == "C0" else "candidate",
        resolved=resolved,
    )


def build(corpus, resolved):
    """A's identical C0/C1 anchor prefixes are reused, never counted twice."""
    cells = []
    for row in corpus["fixtures"]:
        identity, digits = row["id"], row["digits"]
        info = resolved[identity]
        order = ("C0", "C1") if row["n"] % 4 == 1 else ("C1", "C0")
        if identity.startswith("exposed"):
            cells.extend(
                cell(row, info, "A", arm, kind="prefix") for arm in order
            )
            if digits in (62, 70):
                cells.extend(
                    cell(row, info, "A", arm, kind="prefix")
                    for arm in ("noE", "wide")
                )
            if digits <= 70:
                stage = (
                    "B1" if digits in (60, 62) and not row["uneven"] else "B2"
                )
                cells.extend(cell(row, info, stage, arm) for arm in order)
            if digits == 60 and not row["uneven"]:
                cells.extend(
                    cell(row, info, "M", arm, kind="capture") for arm in order
                )
            if digits == 75 and not row["uneven"]:
                cells.extend(cell(row, info, "F", arm) for arm in order)
        elif identity.startswith("regression"):
            cells.extend(cell(row, info, "R", arm) for arm in order)
        elif identity.startswith("routing"):
            cells.extend(
                cell(row, info, "R", arm, kind="routing") for arm in order
            )
    return dict(
        identity=IDENTITY,
        cells=cells,
        ablation_reuse="reuse eight primary cells; add eight noE/wide cells",
        mature_replay=dict(
            stage="A",
            calls=8,
            seconds=180,
            margin=10,
            work_limit=REPLAY_WORK,
            after="B1 and optional M",
        ),
        stages=STAGES,
        generated_bytes=2**31,
        reference_bytes=128 * 2**20,
        complete_80=False,
        confirmation=dict(
            bands=list(CONFIRMATION_SECONDS),
            default=60,
            moduli=12,
            balanced=9,
            uneven=3,
            seeds=list(SEEDS),
            calls=48,
            overhead=4320,
            ceilings=CONFIRMATION_SECONDS,
        ),
    )


def confirmation_cells(corpus, resolved, band, *, nominee="C1"):
    if band not in CONFIRMATION_SECONDS or nominee not in ("C1", "wide"):
        raise ValueError("unfrozen confirmation choice")
    result = []
    for i, row in enumerate(corpus["fixtures"]):
        if row["digits"] != band:
            raise ValueError("confirmation includes another band")
        for j, seed in enumerate(SEEDS):
            order = ("C0", nominee) if (i + j) % 2 == 0 else (nominee, "C0")
            result.extend(
                cell(row, resolved[row["id"]], "C", arm, seed=seed)
                for arm in order
            )
    return result


def limits(plan, *, confirmation_band=60):
    """Bound every alternative by a single cumulative, source-bound ledger."""
    result = {}
    for stage, values in STAGES.items():
        seconds = (
            CONFIRMATION_SECONDS[confirmation_band]
            if stage == "C"
            else values["seconds"]
        )
        measured = sum(
            c["work_limit"] + c["warm_work_limit"] + c["observer_work_limit"]
            for c in plan["cells"]
            if c["stage"] == stage
        )
        if stage == "C":
            measured = 48 * (WORK_CAPS[confirmation_band] + WARM_WORK)
        if stage == "A":
            measured += 8 * REPLAY_WORK
            measured += 18 * (10 * WORK_RATE_GUARD + WARM_WORK) + 9 * 2**32
        # Explicit generation, verification, analysis, observer and controller
        # work reserve. Call grants and failures still bind each operation.
        result[stage] = dict(
            wall=seconds,
            cpu=seconds,
            work=measured + 10**13,
            warm=values["warm"],
            disk_bytes=values["disk"],
        )
    result["total"] = {
        key: sum(values[key] for values in result.values())
        for key in ("wall", "cpu", "work")
    }
    result["disk_bytes"] = 2**31
    return result


def validate(plan):
    counts = {
        stage: sum(c["stage"] == stage for c in plan["cells"])
        for stage in STAGES
    }
    if counts != {"A": 40, "B1": 4, "B2": 20, "M": 2, "R": 24, "F": 2, "C": 0}:
        raise ValueError("plan cell count changed")
    for stage in STAGES:
        selected = [c for c in plan["cells"] if c["stage"] == stage]
        if sum(c["warm_seconds"] for c in selected) > STAGES[stage]["warm"]:
            raise ValueError("warmup/end reservations exceed stage")
        if (
            sum(
                c["seconds"] + c["margin"] + c["warm_seconds"]
                for c in selected
            )
            > STAGES[stage]["seconds"]
        ):
            raise ValueError("worker reservations exceed stage")
    for c in plan["cells"]:
        if c["digits"] == 80 and c["kind"] != "prefix":
            raise ValueError("full 80 cell forbidden")
        if (
            c["configuration"]["memory_bytes"] + HELPER_BYTES
            != c["owned_bytes"]
        ):
            raise ValueError("owned helper budget is not matched")
        if c["owned_bytes"] > 4 * 2**30:
            raise ValueError("native owned ceiling exceeded")
    return counts


def observer_cells(plan):
    """Nine matched one-block clean/captured checks within A's 640 s pool.

    Each reserves 10 s work + 5 s lifecycle + 20 s warm/end = 35 s;
    18 calls reserve 630 s, including 360 s of A's warm allocation.
    These are diagnostic overhead observations, never factoring speed claims.
    """
    anchor = next(
        c
        for c in plan["cells"]
        if c["stage"] == "A"
        and c["digits"] == 60
        and c["fixture"].endswith("_0")
        and c["arm"] == "C1"
    )
    result = []
    for sample in range(9):
        for capture in (False, True) if sample % 2 == 0 else (True, False):
            result.append(
                dict(
                    anchor,
                    sample=sample,
                    kind="capture" if capture else "observer",
                    max_blocks=1,
                    seconds=10,
                    margin=5,
                    warm_seconds=20,
                    warm_work_limit=WARM_WORK,
                    work_limit=10 * WORK_RATE_GUARD,
                    observer_work_limit=2**32 if capture else 0,
                    observer_byte_limit=8 * 2**20 if capture else 0,
                )
            )
    return result
