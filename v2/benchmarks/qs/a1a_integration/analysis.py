"""One frozen input-cluster decision; repeats never create new inputs."""

import itertools
import random
import statistics


def stability(rows):
    """Check repeated cells without inspecting the treatment effect."""
    cells = {}
    for row in rows:
        cells.setdefault((row["id"], row["arm"]), []).append(row)
    reasons = []
    for key, values in cells.items():
        for field in ("seconds", "cpu_seconds"):
            times = [r[field] for r in values]
            if (
                len(times) >= 2
                and (max(times) - min(times)) / statistics.median(times) > 0.15
            ):
                reasons.append(f"repeated-cell drift: {key}/{field}")
        paths = [
            tuple(
                r["stats"].get(k)
                for k in ("blocks", "polynomials", "relations")
            )
            for r in values
        ]
        if all(r["complete"] for r in values) and len(set(paths)) != 1:
            reasons.append("repeated deterministic path changed: " + str(key))
    return reasons


def summarize(rows, *, stable_warmup):
    """Use all retained observations and a 99.5% input bootstrap interval."""
    groups = {}
    for row in rows:
        groups.setdefault(row["id"], {}).setdefault(row["arm"], []).append(row)
    if len(groups) != 9:
        raise ValueError("acceptance requires nine independent inputs")
    pairs = []
    for arms in groups.values():
        if set(arms) != {"control", "selected"}:
            raise ValueError("incomplete matched input")
        old, new = arms["control"], arms["selected"]
        if len(old) != len(new) or len(old) not in (2, 4):
            raise ValueError("incomplete opposite-order schedule")
        expected = set(range(len(old)))
        if {r["repeat"] for r in old} != expected or {
            r["repeat"] for r in new
        } != expected:
            raise ValueError("unmatched repeat identities")
        pairs.append(
            (
                statistics.mean(r["seconds"] for r in old),
                statistics.mean(r["seconds"] for r in new),
            )
        )
    control = sum(p[0] for p in pairs)
    selected = sum(p[1] for p in pairs)
    change = selected / control - 1
    generator = random.Random(2026103055)
    draws = []
    for _ in range(50000):
        sample = [generator.choice(pairs) for _ in pairs]
        draws.append(sum(p[1] for p in sample) / sum(p[0] for p in sample) - 1)
    draws.sort()
    interval = [draws[int((len(draws) - 1) * q)] for q in (0.0025, 0.9975)]
    differences = [b - a for a, b in pairs]
    observed = sum(differences)
    p_value = (
        sum(
            sum(s * d for s, d in zip(signs, differences)) <= observed
            for signs in itertools.product((-1, 1), repeat=9)
        )
        / 512
    )
    reasons = stability(rows)
    if not stable_warmup:
        reasons.append("warmup did not stabilize")
    if not all(r["complete"] for r in rows):
        reasons.append("incomplete/censored outcome; no complete-cost claim")
    if interval[1] >= 0 or p_value > 0.005:
        reasons.append("cost uncertainty gate failed")
    return dict(
        independent_inputs=9,
        pairs=len(rows) // 2,
        cost_change=change,
        interval_995=interval,
        exact_sign_flip_p=p_value,
        complete=sum(r["complete"] for r in rows),
        calls=len(rows),
        cost_gate_passed=not reasons,
        reasons=reasons,
        limitation=(
            "Fixed native SIQS bundles and conditioned certified population; "
            "no portfolio/E1 or additive gain claim."
        ),
    )
