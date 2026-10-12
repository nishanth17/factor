"""Frozen corpora, diagnostic screen and fresh acceptance gates."""

import json
import random
import time

from ...infrastructure.performance.performance_audit import verify_corpus
from ...suites.build_phase_two_corpus import (
    certified_prime,
    verify_certificates,
)
from ..c1.c1_confirmation_support import interval
from ..c10a05.selection import comparison, instability, paired
from ..c10a05_plus.selection import GenerationRandom, existing_inputs
from . import study
from .accounting import Ledger, digest, receipt
from .configuration import ARMS, BANDS


def products(
    digits, count, generator, certificates, seen, used_primes, prefix
):
    fixtures = []
    bits = (digits * 3322 // 1000 + 1) // 2
    for index in range(count):
        for _ in range(10000):
            p = certified_prime(bits, generator, certificates)
            q = certified_prime(bits, generator, certificates)
            if len(certificates) > 4096:
                raise RuntimeError("certificate storage allowance exhausted")
            n = p * q
            if (
                p != q
                and n not in seen
                and not {p, q} & used_primes
                and len(str(n)) == digits
            ):
                break
        else:
            raise RuntimeError("product generation allowance exhausted")
        seen.add(n)
        used_primes.update((p, q))
        fixtures.append(
            dict(
                id=f"{prefix}_{digits}_{index}",
                n=n,
                factors=sorted((p, q)),
                digits=digits,
                kind="balanced",
                split=prefix,
            )
        )
    return fixtures


def training():
    """Reuse revealed 55/60 pilots in training and add a 65-digit band."""
    if study.TRAINING.exists():
        raise ValueError("training corpus already exists")
    prior = json.loads(
        (study.INPUTS / "corpora/c10a05_plus_feasibility.json").read_text()
    )
    verify_corpus(prior)
    certificates, fixtures = dict(prior["certificates"]), []
    seen, used_primes = existing_inputs(), set()
    generator = GenerationRandom(202610120401)
    for digits in BANDS:
        old = [
            dict(f, split="revealed_training")
            for f in prior["fixtures"]
            if f["digits"] == digits
        ]
        fixtures.extend(old)
        used_primes.update(p for f in old for p in f["factors"])
        fixtures.extend(
            products(
                digits,
                3 - len(old),
                generator,
                certificates,
                seen,
                used_primes,
                "c10a04_train",
            )
        )
    verify_certificates(certificates)
    corpus = dict(
        schema=1,
        fixtures=fixtures,
        certificates=certificates,
        seed=202610120401,
        sampling=(
            "certified balanced products conditioned on a large p-1 "
            "factor; revealed 55/60 pilots are training only"
        ),
        generation_seconds=time.monotonic() - generator.started,
    )
    verify_corpus(corpus)
    receipt(study.TRAINING, corpus)


def screen(output):
    """Diagnose one input/seed per band, separately from acceptance."""
    corpus = json.loads(study.TRAINING.read_text())
    verify_corpus(corpus)
    sources = study.hashes()
    identity = dict(
        source=sources, corpus=digest(study.TRAINING), kind="diagnostic"
    )
    limits = {
        "total": dict(wall=1800, cpu=1800, work=10**15),
        "screen": dict(wall=1800, cpu=1800, work=10**15),
        "disk_bytes": 64 * 2**20,
    }
    ledger = Ledger(output, identity, limits)
    for digits in BANDS:
        fixture = next(f for f in corpus["fixtures"] if f["digits"] == digits)
        for arm in ARMS:
            name = f"screen-{digits}-{arm}"
            row = ledger.call(
                "screen",
                name,
                50,
                study.WORK,
                dict(identity, digits=digits, arm=arm),
                lambda: study.run_one(fixture, 7, arm, 45),
            )
            print(
                json.dumps(
                    dict(
                        digits=digits,
                        arm=arm,
                        seconds=row["seconds"],
                        complete=row["complete"],
                        reason=row["reason"],
                        relations=row["stats"].get("relations", 0),
                        base=row["base"],
                    )
                ),
                flush=True,
            )
    # A second, separately charged matched grant probes whether the larger-band
    # first prefix merely precedes cycle formation. It never replaces that row.
    for digits in (60, 65):
        fixture = next(f for f in corpus["fixtures"] if f["digits"] == digits)
        for arm in ARMS:
            name = f"screen-long-{digits}-{arm}"
            row = ledger.call(
                "screen",
                name,
                125,
                study.WORK,
                dict(identity, digits=digits, arm=arm, tier="long"),
                lambda: study.run_one(fixture, 7, arm, 120),
            )
            print(
                json.dumps(
                    dict(
                        digits=digits,
                        arm=arm,
                        tier="long",
                        seconds=row["seconds"],
                        complete=row["complete"],
                        reason=row["reason"],
                        relations=row["stats"].get("relations", 0),
                    )
                ),
                flush=True,
            )


def grants(output, digits):
    """Evaluate a larger matched time grant outside acceptance."""
    study.load()
    checked_rows(output, "train", digits, ARMS)
    corpus = json.loads(study.TRAINING.read_text())
    verify_corpus(corpus)
    fixture = next(f for f in corpus["fixtures"] if f["digits"] == digits)
    ledger = study.expenses(output)
    study.warm(ledger, "extension", digits, fixture, ARMS)
    cap = study.EXTENSION_CAPS[digits]
    missing = [
        a for a in ARMS if not (output / f"grant-{digits}-{a}.json").exists()
    ]
    if not ledger.room(
        "extension", len(missing) * (cap + 5), len(missing) * study.WORK
    ):
        raise RuntimeError(
            "extension allowance cannot cover the matched grant"
        )
    order = list(ARMS)
    random.Random(1044026 + digits).shuffle(order)
    for arm in order:
        row = ledger.call(
            "extension",
            f"grant-{digits}-{arm}",
            cap + 5,
            study.WORK,
            dict(
                control=digest(study.CONTROL), digits=digits, arm=arm, cap=cap
            ),
            lambda arm=arm: study.run_one(fixture, 7, arm, cap),
        )
        print(
            json.dumps(
                dict(
                    arm=arm,
                    complete=row["complete"],
                    reason=row["reason"],
                    seconds=row["seconds"],
                )
            ),
            flush=True,
        )


def checked_rows(output, phase, digits, arms):
    """Reject missing cells, changed schedules and unfinished extensions."""
    finished = json.loads(
        (output / f"{phase}-{digits}-finished.json").read_text()
    )
    count = finished["samples"]
    if count not in (9, 18):
        raise ValueError("invalid matched observation count")
    rows = study.rows(output, phase, digits)
    if len(rows) != count * len(arms):
        raise ValueError("incomplete matched observations")
    corpus = json.loads(
        (study.TRAINING if phase == "train" else study.FRESH).read_text()
    )
    fixtures = [f for f in corpus["fixtures"] if f["digits"] == digits]
    for index in range(count):
        group = [r for r in rows if r["sample"] == index]
        fixture, seed = fixtures[index % 3], study.SEEDS[(index // 3) % 3]
        if {r["arm"] for r in group} != set(arms) or len(group) != len(arms):
            raise ValueError("unmatched group")
        for row in group:
            if row["id"] != fixture["id"] or row["seed"] != seed:
                raise ValueError("sample schedule changed")
            study.validate(row, fixture)
    if finished["reasons"] != instability(rows):
        raise ValueError("stability marker disagrees with outcomes")
    return rows, finished["reasons"]


def select(output):
    study.load()
    if study.SELECTED.exists() or study.FRESH.exists():
        raise ValueError(
            "selection exists or fresh inputs are already revealed"
        )
    choices, rankings, evidence = {}, {}, {}
    study.expenses(output)  # Reconcile receipts before trusting any sample.
    for digits in BANDS:
        rows, reasons = checked_rows(output, "train", digits, ARMS)
        ranked = sorted(
            [
                dict(
                    arm=a,
                    complete=sum(r["complete"] for r in rows if r["arm"] == a),
                    capped_seconds=sum(
                        r["capped_seconds"] for r in rows if r["arm"] == a
                    ),
                )
                for a in ARMS
            ],
            key=lambda r: (
                -r["complete"],
                r["capped_seconds"],
                ARMS.index(r["arm"]),
            ),
        )
        choices[str(digits)] = "control" if reasons else ranked[0]["arm"]
        rankings[str(digits)] = ranked
        for path in sorted(output.glob(f"train-{digits}-[0-9][0-9]-*.json")):
            evidence[path.name] = digest(path)
    receipt(
        study.SELECTED,
        dict(
            schema=1,
            control_sha256=digest(study.CONTROL),
            choices=choices,
            rankings=rankings,
            evidence=evidence,
            limitation=(
                "training selection only; commit before fresh generation"
            ),
        ),
    )


def generate(output):
    """Require committed selection and exclude prior inputs from fresh ones."""
    study.load()
    selected = json.loads(study.SELECTED.read_text())
    if selected["control_sha256"] != digest(study.CONTROL):
        raise ValueError("selection identity changed")
    study.require_committed(study.SELECTED)
    if study.FRESH.exists():
        corpus = json.loads(study.FRESH.read_text())
        if corpus["selection_sha256"] != digest(study.SELECTED):
            raise ValueError("fresh corpus belongs to another selection")
        return

    def create():
        generator = GenerationRandom(202610120402)
        seen, used_primes, certificates, fixtures = (
            existing_inputs(),
            set(),
            {},
            [],
        )
        for digits in BANDS:
            if selected["choices"][str(digits)] != "control":
                fixtures.extend(
                    products(
                        digits,
                        3,
                        generator,
                        certificates,
                        seen,
                        used_primes,
                        "c10a04_fresh",
                    )
                )
        verify_certificates(certificates)
        corpus = dict(
            schema=1,
            fixtures=fixtures,
            certificates=certificates,
            selection_sha256=digest(study.SELECTED),
            seed=202610120402,
            sampling=(
                "certified balanced products conditioned on a large p-1 "
                "factor; all versioned inputs excluded"
            ),
        )
        verify_corpus(corpus)
        receipt(study.FRESH, corpus)
        return dict(work=generator.calls, corpus_sha256=digest(study.FRESH))

    study.expenses(output).call(
        "confirm",
        "fresh-generation",
        120,
        study.WORK,
        dict(selection=digest(study.SELECTED)),
        create,
    )


def completion_interval(rows, arm):
    """Bootstrap input clusters, then completion observations within them."""
    by_input = {}
    for control, candidate in paired(rows, arm):
        by_input.setdefault(control["id"], []).append((control, candidate))
    identities = sorted(by_input)
    generator, differences = random.Random(104001), []
    for _ in range(4000):
        sample = []
        for identity in generator.choices(identities, k=len(identities)):
            group = by_input[identity]
            sample.extend(generator.choices(group, k=len(group)))
        differences.append(
            sum(int(b["complete"]) - int(a["complete"]) for a, b in sample)
            / len(sample)
        )
    return interval(differences)


def summarize(output):
    study.load()
    selected = json.loads(study.SELECTED.read_text())
    report = {}
    for digits in BANDS:
        arm = selected["choices"][str(digits)]
        if arm == "control":
            report[str(digits)] = dict(
                decision="retain control after training"
            )
            continue
        rows, reasons = checked_rows(
            output, "confirm", digits, ("control", arm)
        )
        result = comparison(rows, arm)
        # Completion gains require all fresh candidate calls to complete and
        # no individual input to lose a completion. Capped cost is not speed.
        gain_interval = completion_interval(rows, arm)
        coverage = (
            gain_interval[0] > 0
            and result["complete_challenger"] == result["samples"]
            and result["complete_control"] < result["samples"]
            and all(
                v["complete_challenger"] >= v["complete_control"]
                for v in result["per_input"].values()
            )
        )
        result.update(
            completion_gain_interval=gain_interval,
            coverage_gain_gate=coverage and not reasons,
            decision="accept scoped bundle"
            if result["passes_scoped_gate"] or (coverage and not reasons)
            else "retain control under this allowance",
        )
        report[str(digits)] = result
    receipt(
        output / "summary.json",
        dict(
            bands=report,
            accounting=study.expenses(output).data,
            limitation=(
                "matched configuration envelopes only; no E1/default "
                "promotion or higher-band extrapolation"
            ),
        ),
    )
