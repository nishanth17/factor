"""Bounded per-cell v4 calibration worker; no implicit stage execution."""

import argparse
import hashlib
import importlib
import json
import resource
import signal
import time
from dataclasses import asdict

from ....common import utils
from ....execution.budget import BudgetExhaustedError
from ....qs import CountedSIQSConfig, GrowingDoubleLargeSieveConfig, SIQSJob
from ....qs.sieve_collector import SieveCollector
from . import control_v4, plan_v4

OUTPUT = plan_v4.ROOT / "v2/benchmarks/results/c10a01/calibration-v4"


def rss():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss


def identity(row, corpus):
    matches = [r for r in corpus["fixtures"] if r["id"] == row["fixture"]]
    if len(matches) != 1:
        raise ValueError("fixture identity mismatch")
    fixture = matches[0]
    assert len(str(fixture["n"])) == row["digits"]
    assert fixture["factors"][0] * fixture["factors"][1] == fixture["n"]
    return fixture


def instantiate(row, fixture, modules, *, seconds, work, observer=False):
    qs, budgets, _ = modules
    values = dict(row["configuration"])
    values["collector"] = qs.GrowingDoubleLargeSieveConfig(
        **values["collector"]
    )
    cfg = qs.CountedSIQSConfig(**values)
    cap_rss = plan_v4.RSS[row["tier"]]
    last_poll = 0.0

    def cancel():
        nonlocal last_poll
        now = time.monotonic()
        if now - last_poll < 0.1:
            return False
        last_poll = now
        return rss() > cap_rss

    budget = budgets.Budget(
        work_limit=work, seconds=seconds, cpu_seconds=seconds, cancelled=cancel
    )
    job = qs.SIQSJob(fixture["n"], seed=7, config=cfg, budget=budget)
    return job, budget, cfg


def validate_result(fixture, result, budget, common):
    n = fixture["n"]
    assert (result.divisor or 1) * result.cofactor == n
    factors, certainty = [], []
    if result.divisor is not None:
        assert 1 < result.divisor < n and n % result.divisor == 0
        factors = [result.divisor, result.cofactor]
        assert sorted(factors) == sorted(fixture["factors"])
        for factor in factors:
            try:
                budget.consume(factor.bit_length() ** 2)
                certainty.append(common.classify_prime(factor).value)
            except BudgetExhaustedError:
                certainty.append("unclassified")
    return dict(
        divisor=result.divisor,
        cofactor=result.cofactor,
        factors=factors,
        certainty=certainty,
        unresolved=[
            x
            for x, label in zip(factors, certainty)
            if label in ("composite", "unclassified")
        ]
        if factors
        else [n],
    )


def sieve(row, fixture, modules, *, seconds, work, blocks, observer=False):
    job, budget, cfg = instantiate(
        row, fixture, modules, seconds=seconds, work=work
    )
    collector = modules[0].SieveCollector
    original = collector._sieve
    records = []

    def recorded(self, lo, hi, stats):
        if len(records) >= 8192:
            raise RuntimeError("recorder_exhaustion")
        value = original(self, lo, hi, stats)
        records.append((self.polynomial.identity, lo, hi))
        return value

    started, cpu = time.monotonic(), time.process_time()
    if observer:
        collector._sieve = recorded
    try:
        result = job.run(max_blocks=blocks)
    finally:
        if observer:
            collector._sieve = original
    validated = validate_result(fixture, result, budget, modules[2])
    positions = (
        0 if job.engine is None else job.engine.stats["collected_positions"]
    )
    digest = hashlib.sha256()
    digest.update(str((job.multiplier, job.family_index, positions)).encode())
    if job.engine is not None:
        engine = job.engine
        for key in sorted(engine.collector._atoms):
            digest.update(str(key).encode())
    return dict(
        **validated,
        reason=result.reason,
        work=budget.used,
        seconds=time.monotonic() - started,
        cpu_seconds=time.process_time() - cpu,
        positions=positions,
        polynomials=job.stats["polynomials"],
        families=job.family_index,
        multiplier=job.multiplier,
        d_mod8=None if job.base is None else job.base.n_prime % 8,
        stats=result.stats,
        rss_high_water_bytes=rss(),
        semantic_digest=digest.hexdigest(),
        observer_records=len(records),
        configuration=asdict(cfg),
    )


def warm(row, fixture, modules, *, pre, work):
    start, cpu = time.monotonic(), time.process_time()
    samples, used, calls = [], 0, 0
    while pre and min(time.monotonic() - start, time.process_time() - cpu) < 3:
        if (
            calls >= 26
            or max(time.monotonic() - start, time.process_time() - cpu) >= 15
        ):
            raise RuntimeError("warmup allowance exhausted")
        out = sieve(
            row, fixture, modules, seconds=15, work=work - used, blocks=8 * 33
        )
        if out["reason"] not in ("paused", "factor_found"):
            raise RuntimeError("warmup semantic refusal: " + out["reason"])
        used += out["work"]
        calls += 1
    for _ in range(3):
        if max(time.monotonic() - start, time.process_time() - cpu) >= (
            17 if pre else 8
        ):
            raise RuntimeError("sentinel allowance exhausted")
        out = sieve(
            row, fixture, modules, seconds=8, work=work - used, blocks=8 * 33
        )
        if out["reason"] not in ("paused", "factor_found"):
            raise RuntimeError("sentinel semantic refusal: " + out["reason"])
        used += out["work"]
        samples.append(
            (
                out["seconds"],
                out["cpu_seconds"],
                out["positions"],
                out["semantic_digest"],
            )
        )
    assert calls + 3 <= 32
    return dict(
        wall=time.monotonic() - start,
        cpu=time.process_time() - cpu,
        work=used,
        calls=calls + 3,
        samples=samples,
    )


def diagnostic(row, fixture, modules):
    qs, budgets, _ = modules
    start, cpu = time.monotonic(), time.process_time()
    used = 0
    samples = []
    for _ in range(9):
        job, budget, _ = instantiate(
            row,
            fixture,
            modules,
            seconds=row["seconds"],
            work=row["work"] - used,
        )
        job._setup()
        if job.divisor is not None:
            assert 1 < job.divisor < fixture["n"]
            assert fixture["n"] % job.divisor == 0
        if row["kind"] == "roots" and job.divisor is None:
            state = hashlib.sha256()
            for index in range(8 * 64):
                step = job._next_step()
                if step is None:
                    raise RuntimeError(
                        "root family exhausted before 512 steps"
                    )
                polynomial, roots = step
                assert (
                    polynomial.a * polynomial.value(0)
                    == polynomial.b * polynomial.b - job.base.n_prime
                )
                if index % 64 == 63:
                    state.update(
                        str((polynomial.a, polynomial.b, len(roots))).encode()
                    )
            fingerprint = state.hexdigest()
        else:
            fingerprint = None
        used += budget.used
        samples.append(
            dict(k=job.multiplier, work=budget.used, fingerprint=fingerprint)
        )
    return dict(
        reason="diagnostic_complete",
        samples=samples,
        work=used,
        seconds=time.monotonic() - start,
        cpu_seconds=time.process_time() - cpu,
        rss_high_water_bytes=rss(),
    )


def runtime_modules(row):
    """Keep the owned mainline module/type graph separate from v4."""
    if row["arm"] == "control":
        qs, budgets, common = control_v4.load_control(OUTPUT)
        collector_module = importlib.import_module(
            control_v4.PACKAGE + ".qs.sieve_collector"
        )
        runtime = type(
            "Runtime",
            (),
            dict(
                SIQSJob=qs.SIQSJob,
                CountedSIQSConfig=qs.CountedSIQSConfig,
                GrowingDoubleLargeSieveConfig=qs.GrowingDoubleLargeSieveConfig,
                SieveCollector=collector_module.SieveCollector,
            ),
        )
        return runtime, budgets, common
    from ....execution import budget as budget_module

    runtime = type(
        "Runtime",
        (),
        dict(
            SIQSJob=SIQSJob,
            CountedSIQSConfig=CountedSIQSConfig,
            GrowingDoubleLargeSieveConfig=GrowingDoubleLargeSieveConfig,
            SieveCollector=SieveCollector,
        ),
    )
    return runtime, budget_module, utils


def run_cell(row, corpus):
    fixture = identity(row, corpus)
    OUTPUT.mkdir(parents=True, exist_ok=True)
    modules = runtime_modules(row)
    alarm = row["seconds"] + 30
    signal.signal(
        signal.SIGALRM,
        lambda *_: (_ for _ in ()).throw(
            TimeoutError("cell hard wall allowance exhausted")
        ),
    )
    signal.alarm(alarm)
    limits = json.loads(plan_v4.PLAN.read_text())["calibration"]["limits"]
    warm_work = limits["warm_slot_work"]
    before = warm(row, fixture, modules, pre=True, work=warm_work)
    if row["kind"] in ("roots", "setup"):
        measured = diagnostic(row, fixture, modules)
    else:
        measured = sieve(
            row,
            fixture,
            modules,
            seconds=row["seconds"],
            work=row["work"],
            blocks=row["blocks"],
            observer=row.get("observer", False),
        )
    after = warm(
        row, fixture, modules, pre=False, work=warm_work - before["work"]
    )
    signal.alarm(0)
    if (
        measured["seconds"] > row["seconds"] + 5
        or measured["cpu_seconds"] > row["seconds"] + 5
    ):
        raise RuntimeError("measured cap overshoot")
    return dict(
        id=row["id"],
        fixture=row["fixture"],
        before=before,
        result=measured,
        after=after,
        total_work=before["work"] + measured["work"] + after["work"],
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cell")
    args = parser.parse_args()
    plan = json.loads(plan_v4.PLAN.read_text())
    plan_v4.verify_plan(plan)
    rows = [r for r in plan["calibration"]["calls"] if r["id"] == args.cell]
    if len(rows) != 1:
        raise ValueError("unknown cell")
    corpus = json.loads(plan_v4.CORPUS.read_text())
    result = run_cell(rows[0], corpus)
    encoded = json.dumps(result, sort_keys=True, default=str)
    if len(encoded.encode()) > 64 * 1024:
        raise RuntimeError("recorder_exhaustion: worker stdout")
    print(encoded, flush=True)


if __name__ == "__main__":
    main()
