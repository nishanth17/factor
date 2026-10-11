"""Frozen C10 selection, independent fresh inputs and paired summaries."""

import json
import random
import statistics
import time

from ...infrastructure.performance.performance_audit import verify_corpus
from ...suites.build_phase_two_corpus import (
    certified_prime,
    verify_certificates,
)
from ..b1.b1_calibration import save
from ..c1.c1_confirmation_support import interval
from . import study


def paired(rows, challenger):
    groups = {}
    for row in rows:
        groups.setdefault((row["id"], row["sample"]), {})[row["arm"]] = row
    return [
        (p["control"], p[challenger])
        for p in groups.values()
        if "control" in p and challenger in p
    ]


def instability(rows):
    """Detect paired direction changes and process/repetition drift."""
    reasons = []
    for arm in sorted({r["arm"] for r in rows} - {"control"}):
        pairs = paired(rows, arm)
        effects = [
            (a, b["capped_seconds"] / a["capped_seconds"] - 1)
            for a, b in pairs
        ]
        for identity in {r["id"] for r in rows}:
            values = [v for a, v in effects if a["id"] == identity]
            if values and min(values) < 0 < max(values):
                reasons.append(f"mixed paired direction: {identity}/{arm}")
        ordered = [
            v for a, v in sorted(effects, key=lambda item: item[0]["sample"])
        ]
        if (
            len(ordered) >= 9
            and abs(
                statistics.mean(ordered[:3]) - statistics.mean(ordered[-3:])
            )
            > 0.15
        ):
            reasons.append("early/late paired effect: " + arm)
    cells = {}
    for row in rows:
        cells.setdefault((row["id"], row["seed"], row["arm"]), []).append(row)
    for key, cell in cells.items():
        ordered = sorted(cell, key=lambda r: r["sample"])
        if len(ordered) > 1:
            values = [r["seconds"] for r in ordered]
            if abs(values[-1] - values[0]) / statistics.median(values) > 0.15:
                reasons.append("repeated-cell drift: " + str(key))
    return reasons


def rank(rows):
    results = []
    for arm in study.ARMS:
        group = [r for r in rows if r["arm"] == arm]
        results.append(
            dict(
                arm=arm,
                samples=len(group),
                complete=sum(r["complete"] for r in group),
                capped_seconds=sum(r["capped_seconds"] for r in group),
            )
        )
    return sorted(
        results,
        key=lambda r: (
            -r["complete"],
            r["capped_seconds"],
            study.ARMS.index(r["arm"]),
        ),
    )


def read_rows(output, phase, digits):
    return [
        json.loads(p.read_text())
        for p in sorted(output.glob(f"{phase}-{digits}-[0-9][0-9]-*.json"))
    ]


def select(output):
    if study.SELECTED.exists() or study.CORPUS.exists():
        raise ValueError("selection already exists or fresh inputs revealed")
    control = study.load_control()
    choices, evidence, ranking = {}, [], {}
    for digits in study.TRAIN_CAPS:
        rows = read_rows(output, "train", digits)
        groups = {}
        for row in rows:
            groups.setdefault(row["sample"], set()).add(row["arm"])
        if not groups or any(
            arms != set(study.ARMS) for arms in groups.values()
        ):
            raise ValueError("incomplete matched training groups")
        if len(groups) < 9:
            raise ValueError("selection needs at least nine training samples")
        ranked = rank(rows)
        ranking[str(digits)] = ranked
        best = ranked[0]["arm"]
        if best != "control":
            choices[str(digits)] = best
        for path in sorted(output.glob(f"train-{digits}-[0-9][0-9]-*.json")):
            evidence.append(
                dict(
                    path=str(path.relative_to(study.ROOT)),
                    sha256=study.digest(path),
                )
            )
    save(
        study.SELECTED,
        dict(
            schema=1,
            source=control["source"],
            training_sha256=study.digest(study.CONTROL),
            choices=choices,
            rankings=ranking,
            evidence=evidence,
            limitation="Training selection only; "
            "fresh complete outcomes required",
        ),
    )


class GenerationRandom(random.Random):
    def __init__(self):
        super().__init__(202610111010)
        self.calls = 0
        self.started = time.monotonic()

    def randrange(self, *args, **kwargs):
        self.calls += 1
        if self.calls > 1000000 or time.monotonic() - self.started > 115:
            raise RuntimeError("C10 corpus generation cap")
        return super().randrange(*args, **kwargs)


def generate(output):
    if study.CORPUS.exists():
        raise ValueError("fresh C10 inputs already exist; no regeneration")
    selection = study.load_control(study.SELECTED)
    if selection["training_sha256"] != study.digest(study.CONTROL):
        raise ValueError("training identity changed")
    expense = study.expenses(output)

    def create():
        seen = study.c9_confirmation.existing_inputs()
        generator, certificates, fixtures = GenerationRandom(), {}, []
        for digits in sorted(map(int, selection["choices"])):
            bits = (digits * 3322 // 1000 + 1) // 2
            for index in range(3):
                for _ in range(10000):
                    p = certified_prime(bits, generator, certificates)
                    q = certified_prime(bits, generator, certificates)
                    n = p * q
                    if p != q and n not in seen and len(str(n)) == digits:
                        break
                else:
                    raise RuntimeError("C10 product generation cap")
                seen.add(n)
                fixtures.append(
                    dict(
                        id=f"c10_fresh_{digits}_{index}",
                        n=n,
                        factors=sorted((p, q)),
                        digits=digits,
                        kind="balanced",
                        split="held_out",
                    )
                )
        verify_certificates(certificates)
        corpus = dict(
            schema=1,
            fixtures=fixtures,
            certificates=certificates,
            seed=202610111010,
            selection_sha256=study.digest(study.SELECTED),
            generation_seconds=time.monotonic() - generator.started,
            sampling="Certified balanced primes conditioned on a large "
            "p-1 factor; not RSA-distributed",
        )
        verify_corpus(corpus)
        save(study.CORPUS, corpus)

    expense.call("generation", 120, create)


def comparison(rows, arm):
    pairs = paired(rows, arm)
    by_input = {}
    for control, candidate in pairs:
        by_input.setdefault(control["id"], []).append((control, candidate))
    generator, changes = random.Random(101019), []
    identities = sorted(by_input)
    if not identities:
        return dict(samples=0)
    for _ in range(4000):
        sample = []
        for identity in generator.choices(identities, k=len(identities)):
            source = by_input[identity]
            sample.extend(generator.choices(source, k=len(source)))
        changes.append(
            sum(b["capped_seconds"] for a, b in sample)
            / sum(a["capped_seconds"] for a, b in sample)
            - 1
        )
    result = dict(
        samples=len(pairs),
        independent_inputs=len(identities),
        cost_change=sum(b["capped_seconds"] for a, b in pairs)
        / sum(a["capped_seconds"] for a, b in pairs)
        - 1,
        cost_change_interval=interval(changes),
        complete_control=sum(a["complete"] for a, b in pairs),
        complete_challenger=sum(b["complete"] for a, b in pairs),
        per_input={
            i: dict(
                samples=len(group),
                cost_change=sum(b["capped_seconds"] for a, b in group)
                / sum(a["capped_seconds"] for a, b in group)
                - 1,
                complete_control=sum(a["complete"] for a, b in group),
                complete_challenger=sum(b["complete"] for a, b in group),
            )
            for i, group in by_input.items()
        },
        unstable=instability([r for pair in pairs for r in pair]),
    )
    result["passes_scoped_gate"] = (
        len(pairs) >= 9
        and len(identities) >= 3
        and result["complete_challenger"]
        == result["complete_control"]
        == len(pairs)
        and result["cost_change_interval"][1] < 0
        and not result["unstable"]
    )
    return result


def summarize(output):
    report = {}
    for phase in ("train", "confirm"):
        report[phase] = {}
        for digits in study.TRAIN_CAPS:
            rows = read_rows(output, phase, digits)
            report[phase][str(digits)] = {
                arm: comparison(rows, arm)
                for arm in sorted({r["arm"] for r in rows} - {"control"})
            }
    report["expenses"] = study.expenses(output).data
    report["limitation"] = (
        "Training gates are descriptive only; adoption uses confirmation"
    )
    save(output / "summary.json", report)
