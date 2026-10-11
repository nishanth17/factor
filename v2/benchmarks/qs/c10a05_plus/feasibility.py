"""One revealed input per larger band; diagnostic, not acceptance timing."""

import hashlib
import json
import os
import random
import sys
import time
from dataclasses import asdict

from ....qs import SIQSConfig, default_siqs_config
from ...infrastructure.performance.performance_audit import verify_corpus
from ...suites.build_phase_two_corpus import (
    certified_prime,
    verify_certificates,
)
from ..c9 import c9_growth_study as c9
from ..c10a05 import study
from . import study as extension


def write_once(path, value):
    if path.exists():
        raise RuntimeError(f"refusing to overwrite {path}")
    path.write_text(json.dumps(value, indent=2) + "\n")


def generate(path):
    if path.exists():
        corpus = json.loads(path.read_text())
        verify_corpus(corpus)
        return corpus
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

    for source in study.INPUTS.rglob("*.json"):
        inspect(json.loads(source.read_text()))
    generator = random.Random(202610111201)
    certificates, fixtures = {}, []
    started = time.monotonic()
    for digits in (51, 55, 60):
        bits = (digits * 3322 // 1000 + 1) // 2
        for _ in range(10000):
            if time.monotonic() - started > 115:
                raise RuntimeError("feasibility input generation cap")
            p = certified_prime(bits, generator, certificates)
            q = certified_prime(bits, generator, certificates)
            n = p * q
            if p != q and n not in seen and len(str(n)) == digits:
                break
        else:
            raise RuntimeError("prime generation attempts exhausted")
        seen.add(n)
        fixtures.append(
            dict(
                id=f"c10a05_plus_feas_{digits}",
                n=n,
                factors=sorted((p, q)),
                digits=digits,
                kind="balanced",
                split="inspected_feasibility",
            )
        )
    verify_certificates(certificates)
    corpus = dict(
        schema=1,
        fixtures=fixtures,
        certificates=certificates,
        seed=202610111201,
        sampling="certified balanced primes conditioned on a large p-1 factor",
        generation_seconds=time.monotonic() - started,
    )
    verify_corpus(corpus)
    write_once(path, corpus)
    return corpus


def main():
    c9.runtime()
    source = study.INPUTS / "corpora/c10a05_plus_feasibility.json"
    output = study.ROOT / "v2/benchmarks/results/c10a05_plus/feasibility"
    output.mkdir(parents=True, exist_ok=True)
    corpus = generate(source)
    recorded = {
        "schema": 1,
        "source_commit": "3d4bf42",
        "corpus_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "arms": {
            "current_mainline": asdict(SIQSConfig()),
            "c9_c10_pacing": asdict(default_siqs_config(10**49)),
        },
        "caps": {"current_mainline": 60, "c9_c10_pacing": 240},
        "limit": (
            "one revealed input/seed per band; diagnostic, not acceptance"
        ),
        "max_total_wall": 1200,
    }
    manifest = output / "manifest.json"
    if manifest.exists():
        if json.loads(manifest.read_text()) != recorded:
            raise RuntimeError("feasibility manifest changed")
    else:
        write_once(manifest, recorded)
    started = time.monotonic()
    with study.machine_window():
        for index, fixture in enumerate(corpus["fixtures"]):
            arms = list(recorded["arms"])
            if index % 2:
                arms.reverse()
            for arm in arms:
                receipt = output / f"{fixture['id']}-{arm}.json"
                if receipt.exists():
                    continue
                cap = recorded["caps"][arm]
                if (
                    time.monotonic() - started + cap + 5
                    > recorded["max_total_wall"]
                ):
                    raise RuntimeError("feasibility aggregate allowance")
                config = (
                    SIQSConfig()
                    if arm == "current_mainline"
                    else default_siqs_config(10**49)
                )
                row = c9.run_one(fixture, 7, config, arm, cap)
                row["feasibility"] = True
                write_once(receipt, row)
                print(
                    json.dumps(
                        {
                            k: row[k]
                            for k in (
                                "id",
                                "arm",
                                "complete",
                                "reason",
                                "seconds",
                                "work",
                            )
                        }
                    ),
                    flush=True,
                )


def long55():
    """Probe one revealed 55-digit input under the omitted default grant."""
    c9.runtime()
    source = study.INPUTS / "corpora/c10a05_plus_feasibility.json"
    fixture = next(
        f
        for f in json.loads(source.read_text())["fixtures"]
        if f["digits"] == 55
    )
    output = study.ROOT / "v2/benchmarks/results/c10a05_plus/feasibility"
    config = default_siqs_config(10**49)
    identity = dict(
        source_commit="bdc2330",
        corpus_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        fixture=fixture,
        arm="c9_c10_pacing",
        config=asdict(config),
        seed=7,
        cap_seconds=1200,
        meaning="revealed feasibility; not a fresh acceptance sample",
    )
    manifest = output / "55-longer-identity.json"
    if manifest.exists():
        if json.loads(manifest.read_text()) != identity:
            raise ValueError("55-digit feasibility identity changed")
    else:
        write_once(manifest, identity)
    receipt = output / "55-longer-result.json"
    if receipt.exists():
        print(json.loads(receipt.read_text())["reason"])
        return
    reservation = output / "55-longer-reservation.json"
    if reservation.exists():
        raise RuntimeError("prior reserved call must be audited before retry")
    with study.machine_window():
        write_once(
            reservation,
            dict(
                cap_seconds=1202,
                started_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                pid=os.getpid(),
            ),
        )
        row = c9.run_one(fixture, 7, config, "c9_c10_pacing", 1200)
        row["feasibility"] = True
        write_once(receipt, row)
        print(
            json.dumps(
                {k: row[k] for k in ("complete", "reason", "seconds", "work")}
            ),
            flush=True,
        )


def long60():
    """Check whether the unchanged bundle completes in the default grant."""
    c9.runtime()
    source = study.INPUTS / "corpora/c10a05_plus_feasibility.json"
    fixture = next(
        f
        for f in json.loads(source.read_text())["fixtures"]
        if f["digits"] == 60
    )
    output = study.ROOT / "v2/benchmarks/results/c10a05_plus/feasibility"
    config = default_siqs_config(10**49)
    identity = dict(
        source_commit="bdc2330",
        corpus_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        fixture=fixture,
        arm="c9_c10_pacing",
        config=asdict(config),
        seed=7,
        cap_seconds=1200,
        meaning="revealed feasibility; not a fresh acceptance sample",
    )
    manifest = output / "60-longer-identity.json"
    if manifest.exists():
        if json.loads(manifest.read_text()) != identity:
            raise ValueError("60-digit feasibility identity changed")
    else:
        write_once(manifest, identity)
    receipt = output / "60-longer-result.json"
    if receipt.exists():
        print(json.loads(receipt.read_text())["reason"])
        return
    reservation = output / "60-longer-reservation.json"
    if reservation.exists():
        raise RuntimeError("prior reserved call must be audited before retry")
    with study.machine_window():
        write_once(
            reservation,
            dict(
                cap_seconds=1202,
                started_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                pid=os.getpid(),
            ),
        )
        row = c9.run_one(fixture, 7, config, "c9_c10_pacing", 1200)
        row["feasibility"] = True
        write_once(receipt, row)
        print(
            json.dumps(
                {k: row[k] for k in ("complete", "reason", "seconds", "work")}
            ),
            flush=True,
        )


def windows51():
    """Attribute the C10 plan and pacing on the revealed 51-digit input."""
    c9.runtime()
    source = study.INPUTS / "corpora/c10a05_plus_feasibility.json"
    fixture = next(
        f
        for f in json.loads(source.read_text())["fixtures"]
        if f["digits"] == 51
    )
    output = study.ROOT / "v2/benchmarks/results/c10a05_plus/feasibility"
    configs = extension.configurations()
    arms = ("c9_window", "c10_window")
    manifest = dict(
        corpus_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        fixture=fixture,
        seed=7,
        cap_seconds=360,
        configurations={arm: asdict(configs[arm]) for arm in arms},
        meaning="revealed feasibility, not accepted timing",
    )
    path = output / "51-window-identity.json"
    if path.exists():
        if json.loads(path.read_text()) != manifest:
            raise ValueError("51-digit window identity changed")
    else:
        write_once(path, manifest)
    with study.machine_window():
        for arm in arms:
            receipt = output / f"51-window-{arm}.json"
            if receipt.exists():
                continue
            reservation = output / f"51-window-{arm}.reservation.json"
            if reservation.exists():
                raise RuntimeError("prior reserved call must be audited")
            write_once(
                reservation,
                dict(
                    cap_seconds=362,
                    started_utc=time.strftime(
                        "%Y-%m-%dT%H:%M:%SZ", time.gmtime()
                    ),
                    pid=os.getpid(),
                ),
            )
            row = c9.run_one(fixture, 7, configs[arm], arm, 360)
            row["feasibility"] = True
            write_once(receipt, row)
            print(
                json.dumps(
                    {
                        k: row[k]
                        for k in (
                            "id",
                            "arm",
                            "complete",
                            "reason",
                            "seconds",
                            "work",
                        )
                    }
                ),
                flush=True,
            )


if __name__ == "__main__":
    actions = {
        "short": main,
        "long55": long55,
        "long60": long60,
        "windows51": windows51,
    }
    actions[sys.argv[1] if len(sys.argv) > 1 else "short"]()
