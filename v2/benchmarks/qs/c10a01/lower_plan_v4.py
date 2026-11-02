"""Lower-band supplement to the frozen v4 design; no experiment runner."""

import argparse
import copy
import json
from pathlib import Path

from . import plan_v4

PLAN = Path(__file__).with_name("experiment_lower_v4.json")
SOURCE = Path(__file__).with_name("lower_sources_v4.json")
CORPORA = plan_v4.ROOT / "v2/benchmarks/inputs/corpora"
CAPS = {50: 120, 55: 300}


def configuration(digits, arm):
    """Preserve each calibrated marking bundle under literal matched caps."""
    if digits not in (50, 51, 52, 55, 56):
        raise ValueError("unplanned lower band")
    if arm not in ("control", "nominee"):
        raise ValueError("unknown lower-band arm")
    cfg = copy.deepcopy(plan_v4.config_values(60, "control"))
    if digits <= 51:
        cfg.pop("base_count")
        cfg.update(
            base_bound=30000, factor_count=6, pool_size=32, multiplier=1
        )
    elif digits <= 55:
        cfg.update(base_count=2500, factor_count=7)
    if digits <= 55:
        cfg["collector"]["marking"] = "density+presieve+scan_fused"
    if arm == "nominee":
        cfg.update(
            multiplier=0,
            polynomial_policy="parity2",
            multiplier_policy="parity2_norm",
            assignment_policy="windowed_pair",
            root_update_policy="batch",
        )
        marking = cfg["collector"]["marking"]
        cfg["collector"]["marking"] = (
            "parity_two" if marking == "sparse" else marking + "+parity_two"
        )
    return cfg


def build_plan():
    corpus = json.loads((CORPORA / "a1a_integration_fresh.json").read_text())
    fixtures = {row["id"]: row for row in corpus["fixtures"]}
    prior = json.loads((CORPORA / "c10a05_plus_training.json").read_text())
    fixtures.update({r["id"]: r for r in prior["fixtures"]})
    cells = []
    for digits in (50, 55):
        for ordinal, index in enumerate((0, 1, 4)):
            arms = ("control", "nominee")
            if (ordinal + (digits == 55)) % 2:
                arms = arms[::-1]
            for arm in arms:
                cells.append(
                    dict(
                        id=f"lower-{digits}-{index}-{arm}",
                        fixture=f"a1a_integration_{digits}_{index}",
                        digits=digits,
                        arm=arm,
                        kind="complete",
                        seconds=CAPS[digits],
                        polynomials=0,
                    )
                )
    boundaries = {
        51: "c10a05_plus_feas_51",
        52: "pending_lower_boundary_52",
        55: "a1a_integration_55_0",
        56: "pending_lower_boundary_56",
    }
    for digits, fixture in boundaries.items():
        arms = (
            ("control", "nominee")
            if digits in (51, 55)
            else ("nominee", "control")
        )
        for arm in arms:
            cells.append(
                dict(
                    id=f"boundary-{digits}-{arm}",
                    fixture=fixture,
                    digits=digits,
                    arm=arm,
                    kind="prefix",
                    seconds=30,
                    polynomials=128,
                )
            )
    for cell in cells:
        cell.update(
            reservation_seconds=cell["seconds"] + 5,
            work=plan_v4.call_work(60, cell["seconds"]),
            owned=2**30,
            rss=3 * 2**30,
            configuration=configuration(cell["digits"], cell["arm"]),
            configuration_status="frozen_control"
            if cell["arm"] == "control"
            else "template_requires_global_nominee_freeze",
        )
        if cell["polynomials"]:
            cell.update(positions=128 * 131073, blocks=128 * 33)
    selected = {c["fixture"] for c in cells if c["fixture"] in fixtures}
    inputs = []
    for name in sorted(selected):
        row = fixtures[name]
        p, q = row["factors"]
        assert p * q == row["n"] and p != q
        assert len(str(row["n"])) == row["digits"]
        inputs.append(
            dict(
                id=name,
                n=row["n"],
                digits=row["digits"],
                n_mod8=row["n"] % 8,
                factor_bits=[p.bit_length(), q.bit_length()],
            )
        )
    phases = {
        "complete": dict(
            seconds=2580, work=sum(c["work"] for c in cells[:12])
        ),
        "boundary": dict(seconds=280, work=sum(c["work"] for c in cells[12:])),
        "warmup": dict(seconds=500, work=20 * 32 * plan_v4.call_work(60, 20)),
        "generation_validation": dict(seconds=120, work=10**12),
        "analysis": dict(seconds=120, work=10**12),
        "orchestration": dict(seconds=100, work=0),
    }
    return dict(
        identity="c10a01-v4-lower-design-1",
        execution_ready=False,
        parent_plan_sha256=plan_v4.digest(plan_v4.PLAN),
        parent_sources_sha256=plan_v4.digest(plan_v4.SOURCES),
        seed=7,
        calls=cells,
        inputs=inputs,
        corpora={
            name: plan_v4.digest(CORPORA / name)
            for name in (
                "a1a_integration_fresh.json",
                "c10a05_plus_training.json",
            )
        },
        pending_generation=[
            dict(
                id=boundaries[d],
                digits=d,
                balance="balanced",
                n_mod8=residue,
                seed=202611020150 + d,
                attempt_limit=10000,
            )
            for d, residue in ((52, 5), (56, 7))
        ],
        generation_draw_limit=1000000,
        phases=phases,
        wall_seconds=3700,
        cpu_seconds=3700,
        work=sum(p["work"] for p in phases.values()),
        combined_ceiling_seconds=42300,
        prior_charged_plus_prospective=41824,
        new_output_bytes=8 * 2**20,
        remaining_prerequisites=[
            "Parent bounded worker/ledger adapters and execution release.",
            "Two generated, independently certified boundary inputs.",
            "Freeze global nominee; validate lower-band "
            "composition and resources.",
            "Refresh prospective work/rate admission "
            "after parent calibration.",
        ],
    )


def verify(plan):
    if plan != build_plan():
        raise ValueError("lower plan/input/configuration drift")
    assert len(plan["calls"]) == 20
    assert sum(p["seconds"] for p in plan["phases"].values()) == 3700
    assert 38600 + 3700 == plan["combined_ceiling_seconds"]
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("check",))
    parser.parse_args()
    plan_v4.main()
    verify(json.loads(PLAN.read_text()))
    sources = json.loads(SOURCE.read_text())
    for name, expected in sources.items():
        if plan_v4.digest(plan_v4.ROOT / name) != expected:
            raise RuntimeError("lower supplement source drift: " + name)
    print(
        json.dumps(
            dict(
                lower_plan_verified=True,
                execution_ready=False,
                performance_calls=0,
                lower_seconds=3700,
                combined_ceiling_seconds=42300,
            )
        )
    )


if __name__ == "__main__":
    main()
