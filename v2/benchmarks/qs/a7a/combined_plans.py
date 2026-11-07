"""Pure cell factories for the amended combined-candidate experiment."""

from .combined import FEATURES, cell

LIMITS = {
    "60": (7200, 72000, 8 * 2**30),
    "70": (9000, 90000, 8 * 2**30),
    "ablation": (2700, 27000, 128 * 2**20),
}


def stage1a_cells():
    """Three clean calls and one diagnostic per arm on three60 inputs."""
    inputs = ("exposed_60_0", "combined_60_1", "combined_60_2")
    arms = ("S0", "P4", "C4", "C6")
    rows = []
    for repetition in range(3):
        for index, fixture in enumerate(inputs):
            shift = (repetition + index) % len(arms)
            for arm in arms[shift:] + arms[:shift]:
                rows.append(cell(arm, fixture, rep=repetition + 1))

    for index, fixture in enumerate(inputs):
        for arm in arms[index:] + arms[:index]:
            rows.append(
                cell(arm, fixture, mode="profile", rep=1, wall=360, cpu=3600)
            )
    return rows


def stage1b_cells(candidate="C4"):
    """Two70 inputs, one clean call; candidate selection follows60 report."""
    if candidate not in ("C4", None):
        raise ValueError("70 confirmation accepts C4 or no candidate")
    arms = ("S0", "P4") + ((candidate,) if candidate else ())
    rows = []
    for index, fixture in enumerate(("exposed_70_0", "combined_70_1")):
        for arm in arms[index:] + arms[:index]:
            rows.append(
                cell(arm, fixture, stage="J1b", rep=1, wall=1300, cpu=13000)
            )
    # Conditional diagnostics require a separate source-bound plan that carries
    # these charges. Unused60 time never extends the70 allowance or arm set.
    return rows


def ablation_cells():
    """Two L/R/S leave-one-out calls with bracketing combined controls."""
    rows = [cell("C4", variant="all", stage="J2", wall=240, cpu=2400)]
    for repetition, components in ((1, "LRS"), (2, "SRL")):
        for component in components:
            features = dict(FEATURES)
            if component == "L":
                features["lease_budget"] = False
            elif component == "R":
                features["reuse_state"] = False
            else:
                features.update(lookahead=0, whole_family=False)
            rows.append(
                cell(
                    "C4",
                    features=features,
                    variant="without_" + component,
                    stage="J2",
                    rep=repetition,
                    wall=240,
                    cpu=2400,
                )
            )
    rows.append(cell("C4", variant="all", stage="J2", wall=240, cpu=2400))
    return rows
