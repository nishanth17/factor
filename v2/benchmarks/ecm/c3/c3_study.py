"""Frozen finite allocation study; hidden factors enter validation only."""

import argparse
import contextlib
import fcntl
import hashlib
import json
import math
import os
import platform
import random
import resource
import statistics
import subprocess
import sys
import time
from dataclasses import asdict
from pathlib import Path

from .... import portfolio
from ....common import utils
from ....execution.allocation import ECMAllocation
from ...qs.c1.c1_implementation import decode_config
from ...suites.build_phase_two_corpus import (
    certified_prime,
    verify_certificates,
)
from ...support.paths import BENCHMARK_ROOT, REPOSITORY_ROOT
from ..b4.b4_common import FrozenFinder
from ..p52.p52_realistic import check_quiet

ROOT = REPOSITORY_ROOT
INPUTS = BENCHMARK_ROOT / "inputs"
CONTROL = INPUTS / "controls/c3_protocol.json"
SELECTED = INPUTS / "controls/c3_selected.json"
BASELINE = INPUTS / "baselines/c3_mainline.json"
FRESH = INPUTS / "corpora/c3_confirmation.json"
POLICIES = {
    "control": (((2000, 147396, 32),), None),
    "no_ecm": ((), ("pretest", 2_000_000)),
    "quick8": (((200, 7700, 8),), ("pretest", 2_000_000)),
    "short4": (((2000, 147396, 4),), ("pretest", 2_000_000)),
    "tiered": (((2000, 147396, 8), (11000, 1873422, 2)), ("campaign", None)),
    "wide1": (((50000, 12746592, 1),), ("campaign", None)),
}


def save(path, value):
    data = json.dumps(value, indent=2) + "\n"
    if len(data.encode()) > 16 * 2**20:
        raise ValueError("capture exceeds finite storage cap")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x") as stream:
        stream.write(data)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def runtime():
    if platform.python_implementation() != "PyPy" or sys.version_info[:2] != (
        3,
        11,
    ):
        raise RuntimeError("C3 requires PyPy implementing Python 3.11")


@contextlib.contextmanager
def machine_window():
    """Publish the C3 owner under the repository's machine-wide lock."""
    owner = Path("/private/tmp/factor-performance-owner.json")
    with open("/private/tmp/factor-performance.lock", "a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        owner.write_text(
            json.dumps(
                dict(
                    owner="C3",
                    pid=os.getpid(),
                    cwd=str(ROOT),
                    started=time.time(),
                )
            )
        )
        try:
            check_quiet({os.getpid()})
            yield
            check_quiet({os.getpid()})
        finally:
            owner.unlink(missing_ok=True)


def verify_protocol():
    protocol = json.loads(CONTROL.read_text())
    for name, expected in protocol["sha256"].items():
        if digest(ROOT / name) != expected:
            raise ValueError("frozen C3 source/input changed: " + name)
    return protocol


class SourceFinder(FrozenFinder):
    """Use real resource paths with immutable compiled production bytes."""

    def exec_module(self, module):
        filename = self.path(module.__name__)
        module.__file__ = str(ROOT / filename)
        exec(
            compile(self.sources[filename], module.__file__, "exec"),
            module.__dict__,
        )


def baseline():
    data = json.loads(BASELINE.read_text())
    for name, source in data["source"].items():
        if hashlib.sha256(source.encode()).hexdigest() != data["sha256"][name]:
            raise ValueError("corrupt C3 baseline")
    name = "_c3_mainline"
    if name not in sys.modules:
        sys.meta_path.insert(0, SourceFinder(name, data["source"]))
    return __import__(name + ".portfolio", fromlist=["portfolio"])


def band_of(n):
    """Only an observable size selects the already frozen relation bundle."""
    return 30 if len(str(abs(n))) <= 35 else 40


def config_for(n, arm, module):
    band = band_of(n)
    selected = json.loads((INPUTS / "controls/c1_selected.json").read_text())
    relation = decode_config(
        selected["configurations"][str(band)][
            "slp" if band == 30 else "dlp_half"
        ]
    )
    values = asdict(relation)
    # The historical package owns its unchanged collector types.
    from importlib import import_module

    qs_module = import_module(module.__package__ + ".qs")
    relation = decode_config(values, qs_module)
    tiers, policy = POLICIES[arm]
    options = dict(
        ecm_tiers=tiers,
        siqs=relation,
        memory_bytes=288 * 2**20,
        trace_limit=512,
    )
    if policy is not None:
        options["allocation"] = ECMAllocation(
            policy[0],
            pretest_work=policy[1],
            pretest_seconds=0.5 if policy[0] == "pretest" else None,
            pretest_cpu_seconds=0.5 if policy[0] == "pretest" else None,
            fallback_work=500_000_000 if band == 30 else 15_000_000_000,
            fallback_seconds=1 if band == 30 else 10,
            fallback_cpu_seconds=1 if band == 30 else 10,
        )
    return module.PortfolioConfig(**options)


def validate(result, fixture):
    if result.reconstruct() != fixture["n"]:
        raise AssertionError("factoring result does not reconstruct")
    expected = {int(p): int(e) for p, e in fixture["factors"]}
    for factor in result.factors:
        if (
            factor.value not in expected
            or factor.exponent > expected[factor.value]
        ):
            raise AssertionError("unexpected resolved factor")
        actual = utils.classify_prime(factor.value)
        if (
            actual is utils.Primality.COMPOSITE
            or factor.certainty.value != actual.value
        ):
            raise AssertionError("incorrect primality label")
    if (
        result.complete
        and {int(p.value): p.exponent for p in result.factors} != expected
    ):
        raise AssertionError("incomplete multiplicities")


def run_one(fixture, seed, arm, control, *, regime="service"):
    """Time setup, failed searches, recursion, packing and validation."""
    module = control if arm == "control" else portfolio
    config = config_for(fixture["n"], arm, module)
    band = band_of(fixture["n"])
    seconds = (5 if band == 30 else 30) if regime == "service" else 2
    work = 10**13 if regime == "service" else 2_000_000
    ledger_module = __import__(
        module.__package__ + ".execution.budget", fromlist=["Budget"]
    )
    started, cpu = time.perf_counter(), time.process_time()
    run = module.factorize_bounded(
        fixture["n"],
        seed=seed,
        config=config,
        budget=ledger_module.Budget(
            work_limit=work, seconds=seconds, cpu_seconds=seconds
        ),
    )
    validate(run.result, fixture)
    factor_yield = any(
        utils.valid_divisor(value, abs(fixture["n"]))
        for value in (
            *[p.value for p in run.result.factors],
            *run.result.remaining,
        )
    )
    events = run.events
    return dict(
        id=fixture["id"],
        kind=fixture["kind"],
        band=band,
        seed=seed,
        arm=arm,
        regime=regime,
        complete=run.result.complete,
        proper_factor=factor_yield,
        reason=run.reason,
        work=run.work_used,
        wall=time.perf_counter() - started,
        cpu=time.process_time() - cpu,
        cap_seconds=seconds,
        memory_cap=config.memory_bytes,
        workspace_reserve=config.workspace_reserve,
        rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        curves=sum(
            e["stage"] == "ecm" and e["outcome"] != "handoff" for e in events
        ),
        bounds=[list(tier) for tier in config.ecm_tiers],
        events=events,
        fallback_selected=any(e["stage"] == "handoff" for e in events),
        fallback_started=any(e["stage"] == "siqs" for e in events)
        or (run.checkpoint["payload"]["state"]["current"] or {}).get("stage")
        == "siqs",
        active_curve={
            key: value
            for key, value in (
                (run.checkpoint["payload"]["state"]["current"] or {}).get(
                    "job"
                )
                or {}
            ).items()
            if key
            in (
                "kind",
                "seed",
                "b1",
                "b2",
                "phase",
                "cursor",
                "start_work",
                "chain_replays",
                "chain_chunks",
            )
        },
        fallback_owned_peak=max(
            (e.get("stats", {}).get("workspace_bytes", 0) for e in events),
            default=0,
        ),
        fallback_recoveries=sum(
            e.get("stats", {}).get("recoveries", 0) for e in events
        ),
        checkpoint_bytes=len(json.dumps(run.checkpoint).encode()),
    )


def training():
    fixtures = []
    for name, kinds in (
        ("c1_fresh_corpus.json", {"balanced"}),
        (
            "b1_held_out_corpus.json",
            {"uneven_10", "power", "pm1_smooth", "close"},
        ),
    ):
        corpus = json.loads((INPUTS / "corpora" / name).read_text())
        verify_certificates(corpus["certificates"])
        chosen = set()
        for source in corpus["fixtures"]:
            kind, digits = source["kind"], source["digits"]
            key = kind, digits
            if kind not in kinds or digits not in (30, 40) or key in chosen:
                continue
            chosen.add(key)
            counts = {}
            for factor in source["factors"]:
                counts[factor] = counts.get(factor, 0) + 1
            fixtures.append(
                dict(
                    source,
                    id="training_" + source["id"],
                    factors=sorted(counts.items()),
                )
            )
    return fixtures


def measure(mode, output):
    runtime()
    protocol = verify_protocol()
    control = baseline()
    fixtures = (
        training()
        if mode in ("pilot", "train")
        else json.loads(FRESH.read_text())["fixtures"]
    )
    arms = (
        list(POLICIES)
        if mode in ("pilot", "train")
        else ["control", "selected"]
    )
    selection = json.loads(SELECTED.read_text()) if mode == "confirm" else None
    seeds = (
        protocol["training_seeds"]
        if mode != "confirm"
        else protocol["confirmation_seeds"]
    )
    if mode == "confirm":
        corpus = json.loads(FRESH.read_text())
        verify_certificates(corpus["certificates"])
        if corpus["selected_sha256"] != digest(SELECTED):
            raise ValueError("selection changed after fresh generation")
    envelope = protocol["envelopes"][mode]
    start, cpu_start = time.monotonic(), time.process_time()
    rows = []
    with machine_window():
        for fixture in fixtures:
            if mode == "pilot" and fixture["kind"] != "balanced":
                continue
            active = [
                selection["bands"][str(band_of(fixture["n"]))]
                if a == "selected"
                else a
                for a in arms
            ]
            active = list(dict.fromkeys(active))
            warm = 0 if mode == "pilot" else 3
            for arm in active:
                began = time.monotonic()
                while time.monotonic() - began < warm:
                    run_one(fixture, seeds[0], arm, control)
            sample_count = (
                1 if mode == "pilot" else (3 if mode == "train" else 9)
            )
            for target in [sample_count] if mode != "confirm" else [9, 18, 27]:
                existing = len(
                    [r for r in rows if r["id"] == fixture["id"]]
                ) // (len(active) * len(seeds))
                if target > 9:
                    unstable = False
                    for arm in active:
                        for seed in seeds:
                            values = [
                                r["wall"]
                                for r in rows
                                if r["id"] == fixture["id"]
                                and r["arm"] == arm
                                and r["seed"] == seed
                            ]
                            quartiles = statistics.quantiles(values, n=4)
                            unstable |= (
                                quartiles[2] - quartiles[0]
                            ) / statistics.median(values) > 0.15
                    if not unstable:
                        break
                    for arm in active:
                        began = time.monotonic()
                        while time.monotonic() - began < (
                            5 if target == 18 else 8
                        ):
                            run_one(fixture, seeds[0], arm, control)
                for repetition in range(existing, target):
                    for seed in seeds:
                        order = (
                            active[repetition % len(active) :]
                            + active[: repetition % len(active)]
                        )
                        if (repetition // len(active)) % 2:
                            order = list(reversed(order))
                        for arm in order:
                            if (
                                time.monotonic() - start > envelope
                                or time.process_time() - cpu_start > envelope
                            ):
                                save(
                                    output,
                                    dict(
                                        mode=mode,
                                        stopped="study_allowance",
                                        rows=rows,
                                    ),
                                )
                                return
                            check_quiet({os.getpid()})
                            row = run_one(fixture, seed, arm, control)
                            row["repetition"] = repetition
                            receipt = (
                                output.parent
                                / (output.stem + "-rows")
                                / (
                                    f"{fixture['id']}-{seed}-{repetition}-{arm}.json"
                                )
                            )
                            save(receipt, row)
                            rows.append(
                                {
                                    key: value
                                    for key, value in row.items()
                                    if key != "events"
                                }
                            )
                print(mode, fixture["id"], "samples", target, flush=True)
    save(
        output,
        dict(
            mode=mode,
            rows=rows,
            wall=time.monotonic() - start,
            cpu=time.process_time() - cpu_start,
            runtime=sys.version,
            protocol_sha256=digest(CONTROL),
        ),
    )


def select(path):
    data = json.loads(path.read_text())
    bands = {}
    for band in (30, 40):

        def rank(arm):
            rows = [
                r
                for r in data["rows"]
                if r["band"] == band and r["arm"] == arm
            ]
            complete = sum(r["complete"] for r in rows)
            cost = sum(
                r["wall"] if r["complete"] else r["cap_seconds"] for r in rows
            )
            return (-complete, cost, list(POLICIES).index(arm))

        # A tie or loss may retain the baseline.
        bands[str(band)] = min(POLICIES, key=rank)
    save(
        SELECTED,
        dict(
            bands=bands,
            training_sha256=digest(path),
            protocol_sha256=digest(CONTROL),
        ),
    )


class LimitedRandom(random.Random):
    """Bound prime-generation trials without changing candidate arithmetic."""

    remaining = 200_000

    def getrandbits(self, count):
        self.remaining -= 1
        if self.remaining < 0:
            raise RuntimeError("fresh generation allowance exhausted")
        return super().getrandbits(count)


def generate():
    protocol = verify_protocol()
    selected_hash = digest(SELECTED)
    generator = LimitedRandom(protocol["generation_seed"])
    certificates, fixtures = {}, []

    def prime(bits):
        return certified_prime(bits, generator, certificates)

    def add(kind, factors):
        n = math.prod(p**e for p, e in factors)
        fixtures.append(
            dict(
                id=f"c3_{kind}_{len(fixtures)}",
                kind=kind,
                n=n,
                factors=sorted(factors),
                digits=len(str(n)),
            )
        )

    for digits in (30, 40):
        for small_digits, count in (
            (digits // 2, 2),
            (8, 1),
            (12, 1),
            (16, 1),
        ):
            for _ in range(count):
                for attempt in range(1000):
                    p = prime(small_digits * 3322 // 1000 + 1)
                    q = prime((digits - small_digits) * 3322 // 1000 + 1)
                    if p != q and len(str(p * q)) == digits:
                        break
                else:
                    raise RuntimeError("decimal-band generation exhausted")
                add(
                    "balanced"
                    if small_digits == digits // 2
                    else f"uneven_{small_digits}",
                    [(p, 1), (q, 1)],
                )
    add("recursive", [(prime(20), 1), (prime(33), 1), (prime(48), 1)])
    add("power", [(prime(50), 2)])
    add("prime", [(prime(132), 1)])
    for p in (13001, 13003, 65521):
        certificates[str(p)] = dict(kind="trial")
    add("close", [(13001, 1), (13003, 1)])
    add("pm1_smooth", [(65521, 1), (prime(84), 1)])
    verify_certificates(certificates)
    save(
        FRESH,
        dict(
            fixtures=fixtures,
            certificates=certificates,
            selected_sha256=selected_hash,
            seed=protocol["generation_seed"],
            bias="Pocklington p=kq+1 with large q; nonuniform primes",
        ),
    )


def freeze():
    """Freeze source identity, historical inputs and all finite study rules."""
    ref = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], text=True
    ).strip()
    names = subprocess.check_output(
        ["git", "ls-tree", "-r", "--name-only", "b3b3cfb", "v2"], text=True
    ).splitlines()
    sources = {
        name: subprocess.check_output(
            ["git", "show", "b3b3cfb:" + name]
        ).decode()
        for name in names
        if name.endswith(".py")
        and not name.startswith(("v2/tests/", "v2/benchmarks/"))
    }
    if not BASELINE.exists():
        save(
            BASELINE,
            dict(
                commit="b3b3cfbea6105f08db0a8484ec09c8260d310e28",
                source=sources,
                sha256={
                    name: hashlib.sha256(value.encode()).hexdigest()
                    for name, value in sources.items()
                },
            ),
        )
    paths = [ROOT / name for name in sources]
    paths += list((ROOT / "v2/benchmarks/ecm/c3").glob("*.py"))
    paths += [
        ROOT / "v2/execution/allocation.py",
        BASELINE,
        INPUTS / "controls/c1_selected.json",
        INPUTS / "corpora/c1_fresh_corpus.json",
        INPUTS / "corpora/b1_held_out_corpus.json",
        ROOT / "v2/benchmarks/ecm/c3/protocol.md",
        ROOT / "v2/benchmarks/qs/c1/c1_implementation.py",
        ROOT / "v2/benchmarks/suites/build_phase_two_corpus.py",
        ROOT / "v2/benchmarks/ecm/b4/b4_common.py",
        ROOT / "v2/benchmarks/ecm/p52/p52_realistic.py",
        INPUTS / "controls/c6_fast_records.json",
        INPUTS / "controls/b3_cf_records.json",
    ]
    save(
        CONTROL,
        dict(
            source_commit=ref,
            policies=POLICIES,
            training_seeds=[7, 29],
            confirmation_seeds=[47, 71, 101],
            generation_seed=2026101103,
            sha256={str(p.relative_to(ROOT)): digest(p) for p in paths},
            envelopes=dict(pilot=900, train=2400, confirm=7200),
            sampling=[3, 9, 5, 18, 8, 27],
            max_relative_iqr=0.15,
        ),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "mode",
        choices=("freeze", "pilot", "train", "select", "generate", "confirm"),
    )
    parser.add_argument("--output", type=Path)
    parser.add_argument("--training", type=Path)
    args = parser.parse_args()
    runtime()
    if args.mode == "freeze":
        freeze()
    elif args.mode == "select":
        verify_protocol()
        select(args.training)
    elif args.mode == "generate":
        generate()
    else:
        if args.output is None:
            parser.error("measurement requires --output")
        measure(args.mode, args.output)


if __name__ == "__main__":
    main()
