"""Combined C6/B4 paired study, using an isolated committed B4 control."""

import hashlib
import json
from pathlib import Path

from . import c6_b4, c6_cf_study
from . import c6_fast_study as study

PROTOCOL = Path(__file__).parent / "inputs/controls/c6_b4_protocol.json"
SELECTION = Path(__file__).parent / "inputs/controls/c6_b4_selection.json"
ORIGINAL_BACKEND = study.p41_campaign.PYTHON_BACKEND
ORIGINAL_ATTEMPT = study.attempt
ORIGINAL_CONSTRUCT = study.construct


def construct(arm, backend):
    if arm == "old-ladder":
        return ORIGINAL_CONSTRUCT("ladder", ORIGINAL_BACKEND)
    if arm.startswith("fast/"):
        _, family, mode, batch = arm.split("/")
        return c6_b4.build_program(family, mode, backend, int(batch))
    return ORIGINAL_CONSTRUCT(arm, backend)


def attempt(case, seed, arm, backend, settings, scope, program, oracle):
    if arm == "old-ladder":
        result = ORIGINAL_ATTEMPT(
            case,
            seed,
            "ladder",
            ORIGINAL_BACKEND,
            settings,
            scope,
            program,
            oracle,
        )
        result["arm"] = arm
        return result
    return ORIGINAL_ATTEMPT(
        case, seed, arm, backend, settings, scope, program, oracle
    )


def select(path):
    data, settings = json.loads(path.read_text()), study.protocol()
    choices = ["old-ladder"]
    for family in settings["families"]:
        rows = [
            row
            for row in data["captures"]
            if row["arm"].startswith("fast/" + family + "/")
            and row["summary"]["stable"]
        ]
        if not rows:
            raise RuntimeError("no stable combined candidate")
        best = min(row["summary"]["ratio"] for row in rows)
        winner = min(
            (row for row in rows if row["summary"]["ratio"] <= best * 1.01),
            key=lambda row: (
                settings["modes"].index(row["arm"].split("/")[2]),
                int(row["arm"].split("/")[3]),
            ),
        )
        choices.append(winner["arm"])
    result = dict(
        screen_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        protocol_sha256=hashlib.sha256(PROTOCOL.read_bytes()).hexdigest(),
        candidates={"int": choices},
    )
    with SELECTION.open("x") as output:
        json.dump(result, output, indent=2)
        output.write("\n")
    print(json.dumps(result, indent=2))


def launch(args, arm, backend, warmup, count, *, single=False, profile=None):
    # The shared launcher reads its module specification at call time. Give
    # its private imported module our identity; no production module changes.
    original = c6_cf_study.__spec__
    c6_cf_study.__spec__ = __spec__
    try:
        return c6_cf_study.launch(
            args, arm, backend, warmup, count, single=single, profile=profile
        )
    finally:
        c6_cf_study.__spec__ = original


def install():
    study.PROTOCOL, study.SELECTION = PROTOCOL, SELECTION
    study.construct, study.attempt = construct, attempt
    study.select, study.launch = select, launch
    study.p41_campaign.PYTHON_BACKEND = c6_b4.load_backend(ORIGINAL_BACKEND)


if __name__ == "__main__":
    install()
    study.main()
