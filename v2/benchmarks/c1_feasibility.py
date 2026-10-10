"""Bounded, diagnostic-only residual census beside the unchanged SLP job."""

import argparse
import fcntl
import hashlib
import json
import os
import platform
import random
import subprocess
import time
from collections import Counter
from contextlib import contextmanager
from math import gcd, isqrt, prod
from pathlib import Path
from unittest.mock import patch

from .. import utils
from ..budget import Budget, BudgetExhaustedError
from ..pollard_rho import RhoStats, factorize_rho
from ..qs import SIQSJob
from ..qs.linear_algebra import filter_matrix
from ..qs.sieve_collector import SieveCollector
from .b1_calibration import load_corpus, run_one, upper_config
from .p38_r1_capacity import decode_config
from .performance_audit import fingerprint
from .phase_three_reference import _rss_bytes

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).parent
CONTROL = HERE / "inputs/controls/c1_frozen.json"
LOCK = Path("/private/tmp/factor-performance.lock")
OWNER = Path("/private/tmp/factor-performance-owner.json")


def load_control():
    """Load immutable protocol and the already versioned certified corpus."""
    control = json.loads(CONTROL.read_text())
    for name, expected in control["files"].items():
        if hashlib.sha256((HERE / name).read_bytes()).hexdigest() != expected:
            raise ValueError("changed C1 input: " + name)
    corpus = load_corpus(HERE / "inputs/corpora/p38_r1_training_corpus.json")
    fixtures = [
        f
        for f in corpus["fixtures"]
        if f["kind"] == "balanced" and f["digits"] in control["bands"]
    ]
    configs = {}
    for digits, filename in ((30, "b1_selected"), (40, "b1_40d_selected")):
        selected = json.loads(
            (HERE / f"inputs/controls/{filename}.json").read_text()
        )
        configs[digits] = decode_config(
            selected["configurations"][selected["selected"]["siqs"]]
        )
    configs[60] = upper_config(60, "siqs")
    return control, fixtures, configs


@contextmanager
def machine_window():
    """Acquire the same nonblocking lease as B3/A7, including heavy QA."""
    with LOCK.open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        OWNER.write_text(
            json.dumps(dict(owner="C1", pid=os.getpid(), cwd=str(ROOT)))
        )
        try:
            yield
        finally:
            OWNER.unlink(missing_ok=True)


def save(path, value):
    data = json.dumps(value, indent=2) + "\n"
    if len(data.encode()) > 16 * 2**20:
        raise ValueError("capture byte cap exceeded")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x") as stream:
        stream.write(data)


def recover(collector, position, offset, *, full=False):
    """Recover A*F; full division independently checks root coverage."""
    polynomial = collector.polynomial
    value = polynomial.value(position)
    if not value:
        return 0, ()
    remaining, exponents = abs(value), []
    if full:
        indices = range(len(collector.factor_base.entries))
    else:
        bits = collector._hits[offset]
        for index in collector._a_support:
            bits |= 1 << index
        indices = []
        while bits:
            bit = bits & -bits
            indices.append(bit.bit_length() - 1)
            bits ^= bit
    for index in indices:
        prime = collector.factor_base.entries[index].prime
        exponent = collector._a_exponents[index]
        while remaining % prime == 0:
            remaining //= prime
            exponent += 1
        if exponent:
            exponents.append((prime, exponent))
    reconstructed = prod(p**e for p, e in exponents)
    assert reconstructed * remaining * polynomial.square_coefficient**2 == abs(
        polynomial.a * value
    )
    return remaining, tuple(exponents)


def split_residual(residual, prime_limit, budget, *, seed=7):
    """Finite diagnostic splitter; certify both endpoints independently."""
    bits = residual.bit_length()
    budget.consume(bits**2)
    if utils.classify_prime(residual) is utils.Primality.PROVEN:
        return (), "prime", 0
    root = isqrt(residual)
    if root * root == residual:
        divisor, evaluations = root, 0
    else:
        # Reserve both complete walks, including saturation recovery, first.
        budget.consume(4096 * bits**2)
        stats = RhoStats()
        divisor = factorize_rho(
            residual,
            seed=seed,
            max_attempts=2,
            max_evaluations=2048,
            batch_size=32,
            recovery_limit=64,
            stats=stats,
            _known_composite=True,
        )
        evaluations = stats.evaluations
    if divisor is None:
        return (), "split_failed", evaluations
    if not utils.valid_divisor(divisor, residual):
        raise AssertionError("invalid diagnostic split")
    endpoints = tuple(sorted((divisor, residual // divisor)))
    budget.consume(sum(p.bit_length() ** 2 for p in endpoints))
    if any(
        utils.classify_prime(p) is not utils.Primality.PROVEN
        for p in endpoints
    ):
        return (), "not_two_primes", evaluations
    if endpoints[-1] > prime_limit:
        return (), "endpoint_bound", evaluations
    assert prod(endpoints) == residual
    return endpoints, "dlp", evaluations


def kernel(rows, budget):
    """Complete row-kernel oracle without graph connectivity rules."""
    pivots, dependencies = {}, []
    for index, row in enumerate(rows):
        mask = 1 << index
        budget.consume(1 + (row.bit_length() + index) // 64)
        while row:
            pivot = row & -row
            if pivot not in pivots:
                pivots[pivot] = row, mask
                break
            other, provenance = pivots[pivot]
            budget.consume(1 + (row.bit_length() + index) // 64)
            row ^= other
            mask ^= provenance
        if not row:
            dependencies.append(mask)
    for mask in dependencies:
        assert xor_selected(rows, mask) == 0
    return tuple(dependencies)


def xor_selected(rows, mask):
    result = 0
    while mask:
        bit = mask & -mask
        result ^= rows[bit.bit_length() - 1]
        mask ^= bit
    return result


def verify_square(records, mask, n):
    """Rebuild original identities, square corrections and both GCD signs."""
    powers, sign, x, correction = Counter(), 1, 1, 1
    while mask:
        bit = mask & -mask
        record = records[bit.bit_length() - 1]
        mask ^= bit
        exponents = record["exponents"]
        assert record["u"] ** 2 - n == record["sign"] * (
            record["square"] ** 2
            * record["residual"]
            * prod(p**e for p, e in exponents)
        )
        assert prod(record["lp"]) == record["residual"]
        sign *= record["sign"]
        x = x * record["u"] % n
        correction = correction * record["square"] % n
        for prime, exponent in exponents:
            powers[prime] += exponent
        powers.update(record["lp"])
    assert sign == 1 and all(e % 2 == 0 for e in powers.values())
    y = correction
    for prime, exponent in powers.items():
        y = y * pow(prime, exponent // 2, n) % n
    assert (x * x - y * y) % n == 0
    divisors = [gcd(x - y, n), gcd(x + y, n)]
    for divisor in divisors:
        assert n % divisor == 0
    return {d for d in divisors if utils.valid_divisor(d, n)}


def incidence_report(records, n):
    """Optimistic no-eviction LP elimination followed by existing filtering."""
    budget = Budget(work_limit=10**10, seconds=20, cpu_seconds=20)
    labels = sorted({p for r in records for p in r["lp"]})
    base_labels = sorted({p for r in records for p, _ in r["exponents"]})
    # Reserve simultaneous LP, FB, pivot and provenance integers before work.
    count = len(records)
    reserve = 32768 + 512 * (len(labels) + len(base_labels))
    reserve += count * (
        4096 + 16 * ((3 * count + len(labels) + len(base_labels) + 7) // 8)
    )
    if reserve > 128 * 2**20:
        return dict(censored="matrix_memory", reserve=reserve)
    try:
        lp_columns = {p: 1 << i for i, p in enumerate(labels)}
        fb_columns = {p: 1 << (i + 1) for i, p in enumerate(base_labels)}
        lp_rows, fb_rows = [], []
        for record in records:
            lp, fb = 0, int(record["sign"] < 0)
            for prime in record["lp"]:
                lp ^= lp_columns[prime]
            for prime, exponent in record["exponents"]:
                if exponent % 2:
                    fb ^= fb_columns[prime]
            lp_rows.append(lp)
            fb_rows.append(fb)
        cycles = kernel(lp_rows, budget)
        rows = tuple(xor_selected(fb_rows, mask) for mask in cycles)
        filtered = filter_matrix(
            rows,
            weight_two=True,
            budget=budget,
            memory_bytes=128 * 2**20 - reserve,
        )
        dependencies = kernel(rows, budget)
        divisors = set()
        for dependency in dependencies:
            lifted = xor_selected(cycles, dependency)
            assert lifted and xor_selected(lp_rows, lifted) == 0
            assert xor_selected(fb_rows, lifted) == 0
            budget.consume(count + 1)
            divisors.update(verify_square(records, lifted, n))
        return dict(
            records=count,
            lp_vertices=len(labels),
            independent_lp_constraints=len(cycles),
            post_filter=filtered.stats,
            dependencies=len(dependencies),
            proper_divisors=sorted(divisors),
            work=budget.used,
            reserve=reserve,
            seconds=budget.wall_used,
            cpu_seconds=budget.cpu_used,
        )
    except (BudgetExhaustedError, MemoryError) as error:
        return dict(censored=str(error), work=budget.used, reserve=reserve)


class ResidualAudit:
    """Bounded side-channel census; never admit diagnostic rows to the job."""

    def __init__(self, config, seed):
        self.config, self.seed = config, seed
        self.prime_limit = 100 * config.base_bound
        self.product_limit = self.prime_limit**2
        self.budget = Budget(work_limit=10**11, seconds=60, cpu_seconds=60)
        self.rng = random.Random(380054 + seed)
        self.counts, self.bits, self.samples, self.sample_seen = (
            Counter(),
            Counter(),
            {},
            Counter(),
        )
        self.records = []
        self.retained_bytes = 0
        self.stop = None
        self.split_cpu = self.classify_cpu = 0.0
        self.prefixes = []

    def sample(self, reason, residual):
        """Algorithm R keeps every encountered rejection equally eligible."""
        self.sample_seen[reason] += 1
        samples = self.samples.setdefault(reason, [])
        value = dict(residual=residual, bits=residual.bit_length())
        if len(samples) < 128:
            samples.append(value)
        else:
            index = self.rng.randrange(self.sample_seen[reason])
            if index < 128:
                samples[index] = value

    def classify(self, residual):
        if residual == 1:
            return (), "full"
        if residual > self.product_limit:
            return (), "product_bound"
        started = time.process_time()
        self.budget.consume(residual.bit_length() ** 2)
        prime = utils.classify_prime(residual) is utils.Primality.PROVEN
        self.classify_cpu += time.process_time() - started
        if prime:
            if residual <= self.config.collector.residual_bound:
                return (residual,), "slp"
            return (), "prime_above_slp"
        if self.counts["split_attempts"] >= 8192:
            self.stop = "split_attempt_limit"
            return (), "split_unexamined"
        self.counts["split_attempts"] += 1
        started = time.process_time()
        endpoints, reason, evaluations = split_residual(
            residual, self.prime_limit, self.budget, seed=self.seed
        )
        self.split_cpu += time.process_time() - started
        self.counts["rho_evaluations"] += evaluations
        return endpoints, reason

    def observe(self, collector, lo, hi, threshold):
        """Census widened candidates plus uniform samples of all positions."""
        try:
            self._observe(collector, lo, hi, threshold)
        except BudgetExhaustedError:
            self.stop = "diagnostic_" + str(self.budget.reason)

    def _observe(self, collector, lo, hi, threshold):
        if self.stop:
            return
        self.budget.consume(hi - lo)
        self.counts["blocks"] += 1
        self.counts["positions"] += hi - lo
        lower, _ = collector._bounds(lo, hi)
        wide_threshold = max(
            0, lower.bit_length() - 1 - (self.product_limit - 1).bit_length()
        )
        sampled = set(self.rng.sample(range(hi - lo), min(16, hi - lo)))
        for offset in range(hi - lo):
            self.budget.consume(1)
            score = collector._scores[offset]
            slp_block = score >= threshold
            self.counts[
                "slp_block_pass" if slp_block else "slp_block_reject"
            ] += 1
            wide = score >= wide_threshold
            if not wide and offset not in sampled:
                continue
            position = lo + offset
            value = collector.polynomial.value(position)
            slp_refine = collector._candidate_passes(value, offset)
            wide_refine = score >= max(
                0,
                abs(value).bit_length()
                - 1
                - (self.product_limit - 1).bit_length(),
            )
            if wide:
                self.counts[
                    "wide_refine_pass" if wide_refine else "wide_refine_reject"
                ] += 1
            if not (wide and wide_refine) and offset not in sampled:
                continue
            residual, exponents = recover(collector, position, offset)
            self.budget.consume((len(exponents) + 1) * abs(value).bit_length())
            if offset in sampled:
                self.budget.consume(
                    len(collector.factor_base.entries)
                    * max(1, abs(value).bit_length())
                )
                assert (residual, exponents) == recover(
                    collector, position, offset, full=True
                )
                self.counts["uniform_positions"] += 1
                sample_class = (
                    "uniform_pass"
                    if slp_block and slp_refine
                    else "uniform_reject"
                )
                self.sample(sample_class, residual)
                self.bits[f"{sample_class}:{residual.bit_length()}"] += 1
                if residual and residual <= self.product_limit:
                    assert wide and wide_refine
            if not (wide and wide_refine):
                continue
            self.bits[f"census:{residual.bit_length()}"] += 1
            if not residual:
                self.counts["zero"] += 1
                continue
            endpoints, reason = self.classify(residual)
            self.counts[reason] += 1
            if reason not in ("full", "slp", "dlp"):
                self.sample(reason, residual)
                if self.stop:
                    return
                continue
            if reason in ("full", "slp"):
                assert slp_block and slp_refine
            elif not slp_block:
                self.counts["dlp_lost_at_slp_block"] += 1
            elif not slp_refine:
                self.counts["dlp_lost_at_slp_refinement"] += 1
            else:
                self.counts["dlp_lost_at_slp_residual_policy"] += 1
            reserve = 4096 + 256 * len(exponents)
            if (
                len(self.records) >= 16384
                or self.retained_bytes + reserve > 64 * 2**20
            ):
                self.stop = "diagnostic_record_cap"
                return
            polynomial = collector.polynomial
            self.records.append(
                dict(
                    kind=reason,
                    polynomial=polynomial.identity,
                    position=position,
                    u=polynomial.u_value(position),
                    sign=-1 if value < 0 else 1,
                    square=polynomial.square_coefficient,
                    exponents=exponents,
                    residual=residual,
                    lp=endpoints,
                )
            )
            self.retained_bytes += reserve
        if self.counts["blocks"] in (128, 512, 4096):
            self.prefixes.append(
                dict(
                    blocks=self.counts["blocks"],
                    counts=dict(self.counts),
                    retained_bytes=self.retained_bytes,
                )
            )

    def report(self, n):
        slp = [r for r in self.records if r["kind"] != "dlp"]
        endpoints = Counter(
            p for r in self.records if r["kind"] == "dlp" for p in r["lp"]
        )
        return dict(
            counts=dict(self.counts),
            residual_bits=dict(self.bits),
            rejection_samples=self.samples,
            sample_populations=dict(self.sample_seen),
            split_cpu=self.split_cpu,
            classify_cpu=self.classify_cpu,
            diagnostic_work=self.budget.used,
            stopped=self.stop,
            retained_bytes=self.retained_bytes,
            prefixes=self.prefixes,
            repeated_endpoints=sum(v - 1 for v in endpoints.values()),
            slp_ideal=incidence_report(slp, n),
            dlp_ideal=incidence_report(self.records, n),
        )


def decision(report, control, digits):
    slp, dlp = report["slp_ideal"], report["dlp_ideal"]
    if report["stopped"] or "censored" in slp or "censored" in dlp:
        return False
    added = (
        dlp["independent_lp_constraints"] - slp["independent_lp_constraints"]
    )
    useful = dlp["dependencies"] or dlp["post_filter"]["output_rows"]
    if digits == 60:
        useful = dlp["dependencies"]
    return bool(
        added >= max(16, slp["independent_lp_constraints"] / 10)
        and useful
        and report["split_cpu"] + report["classify_cpu"]
        <= control["cpu_seconds"] / 4
    )


def probe(fixture, seed, config, control):
    audit = ResidualAudit(config, seed)
    budget = Budget(
        work_limit=10**13,
        seconds=90,
        cpu_seconds=90,
        cancelled=lambda: audit.stop is not None,
    )
    job = SIQSJob(fixture["n"], seed=seed, config=config, budget=budget)
    original = SieveCollector._sieve

    def observed(collector, lo, hi, stats):
        threshold = original(collector, lo, hi, stats)
        audit.observe(collector, lo, hi, threshold)
        return threshold

    with patch.object(SieveCollector, "_sieve", observed):
        result = job.run(max_blocks=512)
        report = audit.report(fixture["n"])
        extend = bool(
            not audit.stop
            and "censored" not in report["slp_ideal"]
            and "censored" not in report["dlp_ideal"]
            and report["counts"].get("dlp", 0) >= 32
            and report["repeated_endpoints"] >= 16
            and report["split_cpu"] + report["classify_cpu"]
            < control["cpu_seconds"]
            and result.divisor is None
        )
        initial = report
        if extend:
            budget.seconds = budget.cpu_seconds = 180
            audit.budget.seconds = audit.budget.cpu_seconds = 120
            result = job.run(max_blocks=4096 - 512)
            report = audit.report(fixture["n"])
    assert (result.divisor or 1) * result.cofactor == fixture["n"]
    report.update(
        id=fixture["id"],
        n=fixture["n"],
        seed=seed,
        digits=fixture["digits"],
        extended=extend,
        initial=initial if extend else None,
        slp_result=dict(
            divisor=result.divisor,
            cofactor=result.cofactor,
            reason=result.reason,
            stats=result.stats,
        ),
        instrumented_seconds=budget.wall_used,
        instrumented_cpu=budget.cpu_used,
        production_work=budget.used,
        rss_bytes=_rss_bytes(),
        go=decision(report, control, fixture["digits"]),
    )
    return report, audit.records


def validate_control(row, fixture):
    """Use A10's accepted strict range without editing historical B1 evidence."""
    assert prod(row["factors"]) * prod(row["remaining"]) == fixture["n"]
    expected = Counter(fixture["factors"])
    assert not Counter(row["factors"]) - expected
    assert row["complete"] == (not row["remaining"])
    if row["complete"]:
        assert Counter(row["factors"]) == expected
    if row["divisor"] is not None:
        assert utils.valid_divisor(row["divisor"], fixture["n"])
    assert len(row["certainty"]) == len(row["factors"])
    for factor, label in zip(row["factors"], row["certainty"]):
        # Certificates in the independently verified corpus prove these factors.
        # Runtime guarantees still stop strictly before the A10 endpoint.
        wanted = (
            "proven_prime"
            if factor < 3317044064679887385961981
            else "probable_prime"
        )
        assert label == wanted
    assert row["work"] <= 10**13
    assert row["stats"].get("workspace_bytes", 0) <= 256 * 2**20


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if (
        platform.python_implementation() != "PyPy"
        or platform.python_version_tuple()[:2] != ("3", "11")
    ):
        raise RuntimeError("C1 requires PyPy implementing Python 3.11")
    control, fixtures, configs = load_control()
    args.output.mkdir(parents=True, exist_ok=True)
    with machine_window():
        started, cpu = time.monotonic(), time.process_time()
        summary = dict(
            protocol_sha256=hashlib.sha256(CONTROL.read_bytes()).hexdigest(),
            sources=fingerprint(ROOT),
            driver_sha256=hashlib.sha256(
                Path(__file__).read_bytes()
            ).hexdigest(),
            runtime=platform.python_version(),
            implementation=platform.python_implementation(),
            commit=subprocess.check_output(
                ["git", "rev-parse", "HEAD"], text=True
            ).strip(),
            cells=[],
        )
        if (args.output / "manifest.json").exists():
            original_manifest = json.loads(
                (args.output / "manifest.json").read_text()
            )
            if original_manifest["sources"] != summary["sources"]:
                raise ValueError(
                    "changed production source during continuation"
                )
            save(args.output / "continuation.json", summary)
        else:
            save(args.output / "manifest.json", summary)
        for fixture in fixtures:
            for seed in control["seeds"]:
                # Reserve a whole possible cell before launch.
                if (
                    max(time.monotonic() - started, time.process_time() - cpu)
                    + 310
                    > 2000
                ):
                    summary["stopped"] = "campaign_limit"
                    break
                name = f"{fixture['digits']}-{seed}"
                prior = args.output / (name + "-probe.json")
                if prior.exists():
                    report = json.loads(prior.read_text())
                    summary["cells"].append(
                        dict(
                            name=name,
                            go=report["go"],
                            stopped=report["stopped"],
                        )
                    )
                    continue
                with patch(
                    "v2.benchmarks.b1_calibration.validate", validate_control
                ):
                    baseline = run_one(
                        fixture,
                        seed,
                        configs[fixture["digits"]],
                        5 if fixture["digits"] == 30 else 30,
                    )
                save(args.output / (name + "-control.json"), baseline)
                report, records = probe(
                    fixture, seed, configs[fixture["digits"]], baseline
                )
                save(args.output / (name + "-probe.json"), report)
                save(args.output / (name + "-records.json"), records)
                summary["cells"].append(
                    dict(name=name, go=report["go"], stopped=report["stopped"])
                )
                print(name, json.dumps(summary["cells"][-1]), flush=True)
        summary["seconds"] = time.monotonic() - started
        summary["cpu_seconds"] = time.process_time() - cpu
        summary["go_bands"] = [
            d
            for d in control["bands"]
            if all(
                any(
                    c["name"] == f"{d}-{s}" and c["go"]
                    for c in summary["cells"]
                )
                for s in control["seeds"]
            )
        ]
        save(args.output / "decision.json", summary)


if __name__ == "__main__":
    main()
