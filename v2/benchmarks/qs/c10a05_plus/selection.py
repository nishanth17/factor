"""Precommitted 51-digit selection, certified inputs and confirmation gate."""

import json
import random
import time

from ...infrastructure.performance.performance_audit import verify_corpus
from ...suites.build_phase_two_corpus import (
    certified_prime,
    verify_certificates,
)
from ..b1.b1_calibration import save
from ..c10a05.selection import comparison, instability
from . import study


class GenerationRandom(random.Random):
    def __init__(self, seed):
        super().__init__(seed)
        self.calls = 0
        self.started = time.monotonic()

    def randrange(self, *args, **kwargs):
        self.calls += 1
        if self.calls > 1000000 or time.monotonic() - self.started > 115:
            raise RuntimeError("larger-band certified generation cap")
        return super().randrange(*args, **kwargs)


def existing_inputs():
    seen = set()

    def inspect(value):
        if isinstance(value, dict):
            if type(value.get("n")) is int:
                seen.add(value["n"])
            for child in value.values():
                inspect(child)
        elif isinstance(value, list):
            for child in value:
                inspect(child)

    for path in study.INPUTS.rglob("*.json"):
        inspect(json.loads(path.read_text()))
    return seen


def products(generator, certificates, seen, prefix, count):
    fixtures = []
    bits = (51 * 3322 // 1000 + 1) // 2
    for index in range(count):
        for _ in range(10000):
            p = certified_prime(bits, generator, certificates)
            q = certified_prime(bits, generator, certificates)
            n = p * q
            if p != q and n not in seen and len(str(n)) == 51:
                break
        else:
            raise RuntimeError("51-digit product generation cap")
        seen.add(n)
        fixtures.append(
            dict(
                id=f"{prefix}_{index}",
                n=n,
                factors=sorted((p, q)),
                digits=51,
                kind="balanced",
                split="training" if prefix.endswith("train") else "held_out",
            )
        )
    return fixtures


def training():
    if study.TRAINING.exists():
        raise ValueError("training inputs already generated")
    source = json.loads(
        (study.INPUTS / "corpora/c10a05_plus_feasibility.json").read_text()
    )
    verify_corpus(source)
    first = next(f for f in source["fixtures"] if f["digits"] == 51)
    certificates = dict(source["certificates"])
    seen = existing_inputs()
    generator = GenerationRandom(202610111202)
    fixtures = [first] + products(
        generator, certificates, seen, "c10a05_plus_train", 2
    )
    verify_certificates(certificates)
    corpus = dict(
        schema=1,
        seed=202610111202,
        fixtures=fixtures,
        certificates=certificates,
        generation_seconds=time.monotonic() - generator.started,
        sampling="revealed, certified balanced products conditioned on a "
        "large p-1 factor; first input was inspected in feasibility",
    )
    verify_corpus(corpus)
    save(study.TRAINING, corpus)


def select(output):
    if study.SELECTED.exists() or study.CORPUS.exists():
        raise ValueError("selection exists or fresh inputs have been revealed")
    control = study.load()
    finished = json.loads((output / "train-finished.json").read_text())
    rows = study.read_rows(output, "train")
    by_group = {}
    for row in rows:
        by_group.setdefault(row["sample"], set()).add(row["arm"])
    count = len(by_group)
    if count not in (9, 18) or any(
        arms != set(study.ARMS) for arms in by_group.values()
    ):
        raise ValueError("training lacks complete matched groups")
    if finished["samples"] != count:
        raise ValueError("training finish marker disagrees")
    ranked = sorted(
        (
            dict(
                arm=arm,
                complete=sum(r["complete"] for r in rows if r["arm"] == arm),
                capped_seconds=sum(
                    r["capped_seconds"] for r in rows if r["arm"] == arm
                ),
            )
            for arm in study.ARMS
        ),
        key=lambda r: (
            -r["complete"],
            r["capped_seconds"],
            study.ARMS.index(r["arm"]),
        ),
    )
    chosen = "current" if finished["reasons"] else ranked[0]["arm"]
    evidence = [
        dict(path=str(path.relative_to(study.ROOT)), sha256=study.digest(path))
        for path in sorted(output.glob("train-51-[0-9][0-9]-*.json"))
    ]
    save(
        study.SELECTED,
        dict(
            schema=1,
            source=control["source"],
            training_sha256=study.digest(study.CONTROL),
            chosen=chosen,
            ranked=ranked,
            instability=finished["reasons"],
            evidence=evidence,
            limitation="training choice only; fresh complete outcomes needed",
        ),
    )


def generate(output):
    if study.CORPUS.exists():
        raise ValueError("fresh corpus already exists")
    selected = study.load(study.SELECTED)
    if selected["training_sha256"] != study.digest(study.CONTROL):
        raise ValueError("selection does not match control")
    if selected["chosen"] == "current":
        return
    generator = GenerationRandom(202610111203)
    certificates = {}
    fixtures = products(
        generator,
        certificates,
        existing_inputs(),
        "c10a05_plus_fresh",
        3,
    )
    verify_certificates(certificates)
    corpus = dict(
        schema=1,
        seed=202610111203,
        fixtures=fixtures,
        certificates=certificates,
        selection_sha256=study.digest(study.SELECTED),
        generation_seconds=time.monotonic() - generator.started,
        sampling="certified balanced products conditioned on a large "
        "p-1 factor, generated after committed selection",
    )
    verify_corpus(corpus)
    save(study.CORPUS, corpus)


def summarize(output):
    if (output / "summary.json").exists():
        raise ValueError("summary is immutable")
    report = {}
    for phase in ("train", "confirm"):
        rows = study.read_rows(output, phase)
        report[phase] = {}
        # The prior comparison assumes an arm named control.
        if rows:
            normalized = [
                dict(r, arm="control") if r["arm"] == "current" else r
                for r in rows
            ]
            report[phase] = {
                arm: comparison(normalized, arm)
                for arm in sorted({r["arm"] for r in rows} - {"current"})
            }
            for result in report[phase].values():
                # The reference can stop with no factor almost immediately.
                # Treat full fresh completion as a coverage result, never as
                # a timing speedup against that unsuccessful reference.
                result["coverage_gain_gate"] = (
                    result["sample_plan_complete"]
                    and result["complete_challenger"] == result["samples"]
                    and result["complete_control"] < result["samples"]
                    and not result["unstable"]
                )
        report[phase + "_unstable"] = (
            instability(
                [
                    dict(r, arm="control") if r["arm"] == "current" else r
                    for r in rows
                ]
            )
            if rows
            else []
        )
    report["expenses"] = study.Expenses(output).data
    report["limitation"] = (
        "Only fresh 51-digit results can support 51-digit selection; "
        "baseline early failure is not a timing speedup"
    )
    save(output / "summary.json", report)
