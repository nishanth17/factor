"""Bounded C10a-01 calibration; complete pilots require a separate freeze."""

import time

IMPORT_WALL = time.monotonic()
IMPORT_CPU = time.process_time()

import argparse  # noqa: E402
import fcntl  # noqa: E402
import hashlib  # noqa: E402
import json  # noqa: E402
import resource  # noqa: E402
import subprocess  # noqa: E402
import sys  # noqa: E402
from contextlib import contextmanager  # noqa: E402
from dataclasses import asdict  # noqa: E402
from math import prod  # noqa: E402
from pathlib import Path  # noqa: E402

from ....common import utils  # noqa: E402
from ....execution.budget import Budget, BudgetExhaustedError  # noqa: E402
from ....qs import SIQSJob  # noqa: E402
from ....qs.sieve_collector import SieveCollector  # noqa: E402
from ....qs.sign_schedule import schedule_for  # noqa: E402
from ..c10a04.accounting import Ledger, digest, output_lock  # noqa: E402
from .configuration import RSS, TIERS, configuration  # noqa: E402
from .inputs import CORPUS, validate  # noqa: E402

ROOT = Path(__file__).resolve().parents[4]
BOOTSTRAP = Path(__file__).with_name("bootstrap.json")
PLAN = Path(__file__).with_name("calibration.json")
OUTPUT = ROOT / "v2/benchmarks/results/c10a01/calibration-v1"
LAST_JOB = None
IMPORT_COST = (
    time.monotonic() - IMPORT_WALL,
    time.process_time() - IMPORT_CPU,
)


def runtime():
    if sys.implementation.name != "pypy" or sys.version_info[:2] != (3, 11):
        raise RuntimeError("PyPy implementing Python 3.11 is required")


def host_admission(digits):
    """Admit a sequential tier using actual host headroom."""
    if sys.platform != "darwin":
        raise RuntimeError("this freeze requires a separately checked host")
    ram = int(
        subprocess.check_output(["sysctl", "-n", "hw.memsize"], text=True)
    )
    pressure = subprocess.check_output(["memory_pressure", "-Q"], text=True)
    percent = int(
        pressure.split("System-wide memory free percentage: ")[1].split("%")[0]
    )
    available = ram * percent // 100
    required = RSS[TIERS[digits]] + 4 * 2**30
    if available < required:
        raise RuntimeError("host headroom cannot admit the selected tier")
    return dict(
        ram_bytes=ram,
        available_estimate_bytes=available,
        required_headroom_bytes=required,
        pressure=pressure,
    )


def rss():
    """Process lifetime high water, including warmup and previous calls."""
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(value if sys.platform == "darwin" else value * 1024)


def state_digest(job):
    """Exclude clocks; bind exact provenance, cursors and accounting."""
    engine = job.engine
    value = dict(
        work=job.budget.used,
        family=job.family_index,
        ordinal=None if job.family is None else job.family.next_index,
        next_position=None if engine is None else engine.next_position,
        atoms=[] if engine is None else sorted(engine.collector._atoms),
        rows=[]
        if engine is None
        else [
            (
                getattr(r, "relation_id", None),
                getattr(r, "atom_ids", ()),
                r.sign,
                r.exponents,
                r.u,
                getattr(r, "square_correction", 1),
            )
            for r in engine.collector.matrix_relations
        ],
    )
    return hashlib.sha256(
        json.dumps(value, sort_keys=True).encode()
    ).hexdigest()


@contextmanager
def observe(enabled):
    """A bounded block recorder; exhaustion is a harness error, not memory."""
    records = []
    original = SieveCollector._sieve

    def recorded(collector, lo, hi, stats):
        if len(records) == 8192:
            raise RuntimeError("recorder_exhaustion: block record allowance")
        result = original(collector, lo, hi, stats)
        p = collector.polynomial
        records.append(
            (p.identity, lo, hi, abs(p.value((lo + hi) // 2)).bit_length())
        )
        return result

    if enabled:
        SieveCollector._sieve = recorded
    try:
        yield records
    finally:
        SieveCollector._sieve = original


def run_one(
    n, digits, arm, seconds, work, blocks, *, multiplier=0, observer=False
):
    """Give the algorithm N only; validate every result."""
    global LAST_JOB
    LAST_JOB = None
    started, cpu = time.monotonic(), time.process_time()
    cap = RSS[TIERS[digits]]
    rss_stop = False
    last_poll = started

    def poll():
        nonlocal rss_stop, last_poll
        now = time.monotonic()
        if now - last_poll >= 0.1:
            rss_stop = rss() > cap
            last_poll = now
        return rss_stop

    budget = Budget(
        work_limit=work, seconds=seconds, cpu_seconds=seconds, cancelled=poll
    )
    cfg = configuration(digits, arm, multiplier=multiplier)
    job = SIQSJob(n, seed=7, config=cfg, budget=budget)
    with observe(observer) as records:
        result = job.run(max_blocks=blocks)
    assert (result.divisor or 1) * result.cofactor == n
    factors, labels, unresolved = [], [], [n]
    if result.divisor is not None:
        assert utils.valid_divisor(result.divisor, n)
        factors = [result.divisor, result.cofactor]
        for p in factors:
            try:
                budget.consume(p.bit_length() ** 2)
                labels.append(utils.classify_prime(p).value)
            except BudgetExhaustedError:
                labels.append("unclassified")
        unresolved = [
            p
            for p, label in zip(factors, labels)
            if label in ("composite", "unclassified")
        ]
    engine = job.engine
    stats = result.stats
    population = {}
    if engine is not None:
        collector = engine.collector
        rows = collector.matrix_relations
        population = dict(
            atoms=len(collector._atoms),
            atomic_exponent_pairs=sum(
                len(a.exponents) for a in collector._atoms.values()
            ),
            rows=len(rows),
            row_exponent_pairs=sum(len(r.exponents) for r in rows),
            provenance_references=sum(
                len(getattr(r, "atom_ids", ())) for r in rows
            ),
            collector_workspace=collector._workspace,
            atom_and_combination_reserves=sum(collector._atom_bytes.values()),
            preparation_private_bound=32768
            + sum(
                2048
                + 256 * len(r.exponents)
                + 128 * len(getattr(r, "atom_ids", ()))
                for r in rows
            ),
        )
    row = dict(
        n=n,
        digits=digits,
        arm=arm,
        forced_multiplier=multiplier,
        reason="rss_limit" if rss_stop else result.reason,
        complete=bool(factors) and not unresolved,
        divisor=result.divisor,
        cofactor=result.cofactor,
        factors=factors,
        certainty=labels,
        unresolved=unresolved,
        work=budget.used,
        work_limit=work,
        seconds_cap=seconds,
        seconds=time.monotonic() - started,
        cpu_seconds=time.process_time() - cpu,
        blocks_target=blocks,
        configuration=asdict(cfg),
        stats=stats,
        population=population,
        positions=0 if engine is None else engine.stats["collected_positions"],
        multiplier=job.multiplier,
        d_mod8=None if job.base is None else job.base.n_prime % 8,
        parity=job.family is not None and job.family.parity,
        rss_high_water_bytes=rss(),
        rss_scope="process including warmup",
        state_digest=state_digest(job),
        observer=observer,
        observation_records=records,
    )
    if row["complete"]:
        assert prod(factors) == n
    assert row["work"] <= work
    LAST_JOB = job
    return row


def describe(job):
    """Sample visited polynomials offline using direct CRT."""
    start, cpu = time.monotonic(), time.process_time()
    if job.base is None or job.family is None:
        return dict(work=0, reason="no_family")
    base, cfg = job.base, job.config
    count = job.stats["polynomials"]
    if count > 4096:
        raise RuntimeError(
            "recorder_exhaustion: polynomial metadata allowance"
        )
    parity = cfg.polynomial_policy == "parity2" and base.n_prime % 4 == 1
    schedule = schedule_for(cfg.factor_count, cfg.sign_policy)
    samples, records, charged = {}, [], 0
    old_budget = job.assignments.budget
    analysis_budget = Budget(work_limit=10**10, seconds=2.5, cpu_seconds=2.5)
    job.assignments.budget = analysis_budget
    try:
        for ordinal in range(count):
            index, gray_index = divmod(ordinal, cfg.gray_limit)
            primes = job.assignments[index]
            a = prod(primes) * (2 if parity else 1)
            terms = [
                (a // p)
                * pow(a // p, -1, p)
                * base.entries[base._columns[p] - 1].square_roots[0]
                % a
                for p in primes
            ]
            ordered = [terms[i] for i in (0,) + schedule.variable_order]
            gray = gray_index ^ (gray_index >> 1)
            raw = ordered[0] + (a // 2 if parity else 0)
            raw += sum(
                -v if gray & (1 << i) else v for i, v in enumerate(ordered[1:])
            )
            b = (raw + a // 2) % a - a // 2
            c, remainder = divmod(b * b - base.n_prime, a)
            assert remainder == 0
            hi = cfg.half_width + 1
            if ordinal == count - 1:
                hi = job.engine.next_position
                assert (a, b) == (
                    job.engine.collector.polynomial.a,
                    job.engine.collector.polynomial.b,
                )
            lo = -cfg.half_width
            if hi > lo:
                for j in range(33):
                    x = lo + (hi - lo - 1) * j // 32
                    value = (a * x + 2 * b) * x + c
                    assert (a * x + b) ** 2 - base.n_prime == a * value
                    bits = abs(value).bit_length()
                    samples[str(bits)] = samples.get(str(bits), 0) + 1
                    analysis_budget.consume(bits + 1)
                records.append(
                    dict(
                        family=index,
                        ordinal=gray_index,
                        a=a,
                        b=b,
                        lo=lo,
                        hi=hi,
                        u_first=a * lo + b,
                        u_last=a * (hi - 1) + b,
                        u_stride=a,
                    )
                )
        charged = analysis_budget.used
    finally:
        job.assignments.budget = old_budget
    assert (
        sum(r["hi"] - r["lo"] for r in records)
        == job.engine.stats["collected_positions"]
    )
    return dict(
        work=charged,
        records=records,
        sampled_norm_bit_histogram=samples,
        sampling="33 equally spaced coordinates per visited polynomial",
        wall=time.monotonic() - start,
        cpu=time.process_time() - cpu,
    )


def warm(n, digits, arm, *, multiplier=0):
    """Reuse one process until both validated execution clocks exceed 3 s."""
    start, cpu = time.monotonic(), time.process_time()
    work, calls = 0, 0
    while min(time.monotonic() - start, time.process_time() - cpu) < 3:
        if max(time.monotonic() - start, time.process_time() - cpu) >= 25:
            raise RuntimeError("warmup allowance exhausted")
        row = run_one(
            n,
            digits,
            arm,
            25 - max(time.monotonic() - start, time.process_time() - cpu),
            290000000000 - work,
            264,
            multiplier=multiplier,
        )
        if row["reason"] not in ("paused", "factor_found"):
            raise RuntimeError(
                "warmup did not complete its semantic work: " + row["reason"]
            )
        work += row["work"]
        calls += 1
    return dict(
        work=work,
        calls=calls,
        wall=time.monotonic() - start,
        cpu=time.process_time() - cpu,
        last=row,
        import_wall=IMPORT_COST[0],
        import_cpu=IMPORT_COST[1],
    )


def verify_outcome(row, fixture):
    if row["divisor"] is not None:
        assert set(row["factors"]) == set(fixture["factors"])
    assert row["n"] == fixture["n"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("warm", "prefix", "diagnostics"))
    parser.add_argument(
        "--digits", type=int, choices=(60, 62, 65, 70), required=True
    )
    parser.add_argument("--arm", choices=("control", "parity2"), required=True)
    args = parser.parse_args()
    runtime()
    corpus = json.loads(CORPUS.read_text())
    validate(corpus)
    bootstrap = json.loads(BOOTSTRAP.read_text())
    plan = json.loads(PLAN.read_text())
    fixtures = [r for r in corpus["fixtures"] if r["digits"] == args.digits]
    source = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
    ).strip()
    subprocess.run(
        ["git", "diff", "--exit-code", "HEAD", "--", "v2"],
        cwd=ROOT,
        check=True,
        stdout=subprocess.DEVNULL,
    )
    identity = dict(source=source, plan=digest(PLAN), corpus=digest(CORPUS))
    host = host_admission(args.digits)
    print(json.dumps(dict(event="host_admission", **host)), flush=True)
    with open("/private/tmp/factor-performance.lock", "a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        with output_lock(OUTPUT):
            ledger = Ledger(OUTPUT, digest(BOOTSTRAP), bootstrap["limits"])
            tag = f"{args.action}-{args.digits}-{args.arm}"
            if (OUTPUT / ("warm-" + tag + ".json")).exists():
                raise RuntimeError(
                    "process warmup cannot be replayed from a receipt"
                )
            warming = ledger.call(
                "warmup",
                "warm-" + tag,
                35,
                290000000000,
                identity,
                lambda: warm(fixtures[0]["n"], args.digits, args.arm),
            )
            print(
                json.dumps(
                    dict(
                        event="warmup",
                        tag=tag,
                        wall=warming["wall"],
                        cpu=warming["cpu"],
                    )
                ),
                flush=True,
            )
            capacity = warming["last"]["stats"].get("capacity", {})
            if not capacity.get("target_in_product_envelope"):
                raise RuntimeError(
                    "A target is outside the attainable product envelope"
                )
            if args.action == "warm":
                return
            if args.action == "prefix":
                grant = bootstrap["prefix"][str(args.digits)]
                calls = [
                    (
                        f,
                        False,
                        plan["prefix_polynomials"][str(args.digits)] * 33,
                    )
                    for f in fixtures
                ]
                phase = "prefix"
            else:
                if args.digits not in (60, 70):
                    raise ValueError("observer diagnostics only at 60 and 70")
                grant = bootstrap["diagnostics"]["observer"]
                calls = [
                    (
                        fixtures[0],
                        flag,
                        plan["observer_polynomials"][str(args.digits)] * 33,
                    )
                    for flag in (False, True)
                ]
                phase = "observer"
            previous = None
            for fixture, observer, blocks in calls:
                name = (
                    tag + "-" + fixture["id"] + ("-on" if observer else "-off")
                )
                call_identity = dict(
                    identity,
                    fixture=fixture["id"],
                    blocks=blocks,
                    observer=observer,
                )
                row = ledger.call(
                    phase,
                    name,
                    grant["seconds"] + 5,
                    grant["work"],
                    call_identity,
                    lambda: run_one(
                        fixture["n"],
                        args.digits,
                        args.arm,
                        grant["seconds"],
                        grant["work"],
                        blocks,
                        observer=observer,
                    ),
                )
                verify_outcome(row, fixture)
                print(
                    json.dumps(
                        dict(
                            event=phase,
                            name=name,
                            seconds=row["seconds"],
                            cpu=row["cpu_seconds"],
                            reason=row["reason"],
                            positions=row["positions"],
                            work=row["work"],
                        )
                    ),
                    flush=True,
                )
                if phase == "observer":
                    if previous is not None:
                        assert (
                            previous["state_digest"] == row["state_digest"]
                        ), "observer changed semantic work"
                    previous = row
                if LAST_JOB is not None:
                    ledger.call(
                        "analysis",
                        name + "-coordinates",
                        3,
                        10**10,
                        call_identity,
                        lambda: describe(LAST_JOB),
                    )


if __name__ == "__main__":
    main()
