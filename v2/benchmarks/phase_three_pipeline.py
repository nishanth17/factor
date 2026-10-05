"""Matched P3.3 diagnosis, complete QS splits, filtering and recovery costs."""

import argparse
import cProfile
import hashlib
import json
import math
import random
import resource
import statistics
import subprocess
import sys
import time
from dataclasses import asdict
from pathlib import Path

from .. import utils
from ..budget import Budget
from ..tests.test_qs import reference_positions
from .phase_one import environment
from .phase_three_collector import TRAINING, _signature
from .phase_three_reference import BASELINE, _cases, _measure, _rss_bytes
from .qs_snapshot import ROOT, load_qs_arm

CHANGES = ("qs.factor_base", "qs.sieve_collector")


def _budget():
    """Use identical declared finite allowances for all candidate arms."""
    return Budget(work_limit=200_000_000, seconds=10, cpu_seconds=10)


def _corpus(seed, count):
    """Generate known prime products; factors only validate outputs."""
    generator = random.Random(seed)
    fixtures = []
    while len(fixtures) < count:
        primes = []
        for lo, hi in ((1200, 4500), (4600, 9000)):
            candidate = generator.randrange(lo, hi) | 1
            while utils.classify_prime(candidate) != utils.Primality.PROVEN:
                candidate += 2
            primes.append(candidate)
        p, q = primes
        fixtures.append({"n": p * q, "factors": sorted(primes)})
    return fixtures


def _config(arm, **options):
    """Hold store caps/memory fixed; old sources lack only new controls."""
    values = dict(
        block_width=256,
        metadata_chunk=64,
        residual_bound=500,
        max_partials=4096,
        max_relations=4096,
        max_atoms=4096,
        memory_bytes=32 * 1024 * 1024,
    )
    values.update(options)
    config = arm.qs.sieve_collector.SieveConfig
    return config(
        **{
            key: value
            for key, value in values.items()
            if key in config.__dataclass_fields__
        }
    )


def _collection(arm, fixture, config, *, exhaustive=False):
    """Return consumed checked payloads, including common partial matching."""
    budget = _budget()
    base = arm.qs.factor_base.build_factor_base(
        fixture["n"],
        multiplier=fixture["h"],
        bound=fixture["bound"],
        budget=budget,
        memory_bytes=config.memory_bytes,
    ).factor_base
    polynomial = arm.qs.polynomial.Polynomial(
        fixture["n"],
        fixture["h"],
        fixture["a"],
        fixture["b"],
    )
    worker = arm.qs.sieve_collector.SieveCollector(
        polynomial,
        base,
        config=config,
        budget=budget,
    )
    if exhaustive:
        result = arm.qs.reference_collector.collect_block(
            polynomial,
            base,
            -128,
            129,
            residual_bound=fixture["residual"],
            budget=budget,
            memory_bytes=config.memory_bytes,
        )
        if result.reason != "complete":
            raise AssertionError("exhaustive control exhausted")
        stats = dict(
            duplicates=0,
            dropped_partials=0,
            evictions=0,
            admitted_atoms=0,
            matches=0,
        )
        for atom in result.relations:
            if worker._admit(atom, stats):
                raise AssertionError("exhaustive store refused")
        atoms, full, combined = (
            tuple(worker._atoms.values()),
            worker._full,
            worker._combined,
        )
    else:
        result = worker.collect(-128, 129)
        if result.reason != "complete":
            raise AssertionError("sieve exhausted")
        atoms, full, combined = (
            result.atoms,
            result.full_relations,
            result.combined_relations,
        )
    store = {atom.relation_id: atom for atom in atoms}
    for relation in full:
        arm.qs.relations.verify_atomic(relation, base, budget=budget)
    for relation in combined:
        arm.qs.relations.verify_combined(
            relation,
            base,
            store,
            budget=budget,
            memory_bytes=config.memory_bytes,
        )
    return _signature(atoms), len(full), len(combined)


def _complete(
    arm,
    corpus,
    config,
    *,
    weight_two=False,
    pivot="highest",
    row_excess=2,
    batch_width=256,
    details=False,
    exhaustive=False,
):
    """Complete each balanced QS fixture and prove both factors prime."""
    signatures, metadata = [], []
    for fixture in corpus:
        budget = _budget()
        n = fixture["n"]
        setup = arm.qs.factor_base.build_factor_base(
            n,
            bound=200,
            budget=budget,
            memory_bytes=config.memory_bytes,
        )
        if setup.divisor:
            raise AssertionError("fixture bypasses the relation pipeline")
        base = setup.factor_base
        job = arm.qs.pipeline.QSJob(
            arm.qs.polynomial.qs_polynomial(base),
            base,
            -512,
            2049,
            config=config,
            budget=budget,
            weight_two=weight_two,
            pivot=pivot,
            row_excess=row_excess,
            batch_width=batch_width,
            collector_class=(
                _exhaustive_collector(arm)
                if exhaustive
                else arm.qs.sieve_collector.SieveCollector
            ),
        )
        result = job.run()
        if result.reason != "factor_found":
            raise AssertionError(f"QS did not complete: {n} {result}")
        factors = sorted((result.divisor, result.cofactor))
        if factors != fixture["factors"] or factors[0] * factors[1] != n:
            raise AssertionError("complete QS reconstruction failed")
        if any(
            utils.classify_prime(value) != utils.Primality.PROVEN
            for value in factors
        ):
            raise AssertionError("QS leaves an unresolved cofactor")
        signatures.append(factors)
        metadata.append(result.stats)
    return metadata if details else signatures


def _exhaustive_collector(arm):
    """Adapt checked exhaustive atoms to the common matching/extraction store.

    This benchmark control charges the live transfer copy as additional
    workspace. Refused/failed output never qualifies as a faster completion.
    """

    class ExhaustiveCollector(arm.qs.sieve_collector.SieveCollector):
        """Collect supplied batches without score marking or root exclusion."""

        def collect(self, lo, hi):
            """Transfer one complete checked batch under the same budget."""
            raw = arm.qs.reference_collector.collect_block(
                self.polynomial,
                self.factor_base,
                lo,
                hi,
                residual_bound=self.config.residual_bound,
                budget=self.budget,
                memory_bytes=self.config.memory_bytes,
            )
            if raw.reason != "complete":
                raise AssertionError("exhaustive pipeline control refused")
            stats = dict(
                scanned=raw.scanned,
                duplicates=0,
                dropped_partials=0,
                evictions=0,
                admitted_atoms=0,
                matches=0,
            )
            for atom in raw.relations:
                if (
                    self._workspace + raw.workspace_bytes
                    > self.config.memory_bytes
                ):
                    raise AssertionError(
                        "exhaustive transfer memory exhausted"
                    )
                if self._admit(atom, stats):
                    raise AssertionError("exhaustive pipeline store refused")
            return arm.qs.sieve_collector.SieveResult(
                tuple(self._atoms.values()),
                tuple(self._full),
                tuple(self._combined),
                tuple(self._pending.values()),
                None,
                hi,
                "complete",
                stats,
                self._workspace + raw.workspace_bytes,
            )

    return ExhaustiveCollector


def _recovery(arm, config):
    """Time exact recovery on identical candidates, including resieving."""
    fixture = TRAINING[2]
    base = arm.qs.factor_base.build_factor_base(
        fixture["n"],
        multiplier=fixture["h"],
        bound=fixture["bound"],
        budget=_budget(),
    ).factor_base
    polynomial = arm.qs.polynomial.Polynomial(
        fixture["n"],
        fixture["h"],
        fixture["a"],
        fixture["b"],
    )
    worker = arm.qs.sieve_collector.SieveCollector(
        polynomial,
        base,
        config=config,
        budget=_budget(),
    )
    stats = dict(
        root_hits=0,
        slice_bytes=0,
        slice_allocations=0,
        blocks=0,
        threshold_min=2**63,
        threshold_max=0,
    )
    worker._sieve(-128, 129, stats)

    def call():
        """Consume exact payloads after residual certainty and verification."""
        worker.budget = _budget()
        stats = dict(
            candidates=0,
            zeros=0,
            division_primes=0,
            division_steps=0,
            composite_residuals=0,
        )
        if config.division == "resieve":
            if not worker._resieve(-128, 129, 0, stats):
                raise AssertionError("resieve utility memory refused")
        atoms = []
        for position in range(-128, 129):
            atom, divisor = worker._divide(position, position + 128, stats)
            if divisor:
                raise AssertionError(
                    "recovery utility found an unexpected GCD"
                )
            if atom:
                atoms.append(atom)
        return _signature(atoms)

    return call, reference_positions(polynomial, base, -128, 129, 500)


def _record(name, call, expected, args, **metadata):
    """Validate warmups and every consumed timing output; preserve samples."""
    measured = _measure(
        call, expected, 10, args.repetitions, args.warmup_seconds
    )
    print(
        name,
        round(measured["median_seconds"] * 1000, 3),
        "ms",
        "stable",
        measured["stable"],
        flush=True,
    )
    return {"name": name, "measurement": measured, **metadata}


def main():
    """Tune on training inputs, freeze, then open fresh held-out evaluation."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--warmup-seconds", type=float, default=3)
    parser.add_argument("--repetitions", type=int, default=9)
    parser.add_argument("--cold-worker", action="store_true")
    parser.add_argument(
        "--cold-arm", choices=("before", "after"), default="after"
    )
    args = parser.parse_args()
    if args.cold_worker:
        arm = load_qs_arm(
            "_qs_cold", CHANGES if args.cold_arm == "after" else ()
        )
        frozen = json.loads(
            args.output.with_suffix(".frozen.json").read_text()
        )
        values = dict(frozen["config"])
        if args.cold_arm == "before":
            values["division"] = "roots"
        config = _config(arm, **values)
        started, cpu = time.perf_counter(), time.process_time()
        factors = _complete(
            arm,
            _corpus(frozen["held_out_seed"], frozen["held_out_count"]),
            config,
            **frozen["controls"],
        )
        usage = resource.getrusage(resource.RUSAGE_SELF)
        print(
            json.dumps(
                {
                    "factors": factors,
                    "work_wall_seconds": time.perf_counter() - started,
                    "work_cpu_seconds": time.process_time() - cpu,
                    "peak_rss_bytes": _rss_bytes(),
                    "total_cpu_seconds": usage.ru_utime + usage.ru_stime,
                }
            )
        )
        return
    if (
        not math.isfinite(args.warmup_seconds)
        or args.warmup_seconds < 3
        or args.repetitions < 9
    ):
        parser.error("require three seconds warmup and nine samples")
    if (
        args.output.exists()
        or args.output.with_suffix(".frozen.json").exists()
    ):
        parser.error("use a unique capture filename")
    frozen_production = json.loads(BASELINE.read_text())
    production_unchanged = all(
        hashlib.sha256((ROOT.parent / name).read_bytes()).hexdigest() == digest
        for name, digest in frozen_production[
            "production_source_sha256"
        ].items()
    )
    if not production_unchanged:
        raise AssertionError("M17 production sources changed")
    arms = {
        name: load_qs_arm("_qs_" + name, changes)
        for name, changes in (
            ("before", ()),
            ("collector_only", ("qs.sieve_collector",)),
            ("accounting_only", ("qs.sieve_collector",)),
            ("cache_only", ("qs.factor_base",)),
            ("after", CHANGES),
        )
    }
    rows = []
    for fixture in TRAINING:
        # The independent generic oracle factors complete norms, not score
        # or root metadata. Exact full/combined yields are matched as well.
        baseline_config = _config(
            arms["before"], residual_bound=fixture["residual"]
        )
        expected = _collection(arms["before"], fixture, baseline_config)
        base = (
            arms["before"]
            .qs.factor_base.build_factor_base(
                fixture["n"], multiplier=fixture["h"], bound=fixture["bound"]
            )
            .factor_base
        )
        poly = arms["before"].qs.polynomial.Polynomial(
            fixture["n"], fixture["h"], fixture["a"], fixture["b"]
        )
        oracle = reference_positions(
            poly, base, -128, 129, fixture["residual"]
        )
        assert expected[0] == oracle
        for name, arm in arms.items():
            config = _config(arm, residual_bound=fixture["residual"])
            if name == "accounting_only":
                config = _config(
                    arm,
                    residual_bound=fixture["residual"],
                    score_policy="conservative",
                )
            rows.append(
                _record(
                    fixture["name"] + "_" + name,
                    lambda: _collection(arm, fixture, config),
                    expected,
                    args,
                    scope="matched_collection",
                    config=asdict(config),
                )
            )
        arm = arms["after"]
        config = _config(arm, residual_bound=fixture["residual"])
        rows.append(
            _record(
                fixture["name"] + "_exhaustive",
                lambda: _collection(arm, fixture, config, exhaustive=True),
                expected,
                args,
                scope="matched_collection",
            )
        )
    for division in ("full", "roots", "bucket", "resieve"):
        arm = arms["after"]
        config = _config(arm, division=division, block_width=512)
        call, expected = _recovery(arm, config)
        rows.append(
            _record(
                "recovery_" + division,
                call,
                expected,
                args,
                scope="exact_recovery_utility",
                candidate_assignment=[-128, 129],
                config=asdict(config),
            )
        )
    training = _corpus(329, 4)
    expected = [fixture["factors"] for fixture in training]
    tuning = []
    for policy, division in (
        (policy, division)
        for policy in ("adaptive", "candidate")
        for division in ("full", "roots", "bucket", "resieve")
    ):
        arm = arms["after"]
        config = _config(
            arm, division=division, residual_bound=1, score_policy=policy
        )
        row = _record(
            "training_" + policy + "_" + division,
            lambda: _complete(arm, training, config),
            expected,
            args,
            scope="complete_training",
            config=asdict(config),
            details=_complete(arm, training, config, details=True),
        )
        rows.append(row)
        if row["measurement"]["stable"]:
            tuning.append(row)
    winner = min(tuning, key=lambda row: row["measurement"]["median_seconds"])
    # Filter/pivot/stopping experiments use training only. Their kernels are
    # separately checked with the independent dense oracle in acceptance.
    arm = arms["after"]
    config = _config(arm, **winner["config"])
    filter_choices = []
    for weight_two, pivot, excess, batch in (
        (False, "highest", 2, 256),
        (False, "highest", 0, 64),
        (False, "lowest", 2, 256),
        (True, "highest", 2, 256),
        (True, "lowest", 2, 256),
    ):
        kwargs = dict(
            weight_two=weight_two,
            pivot=pivot,
            row_excess=excess,
            batch_width=batch,
        )
        expected_training = [fixture["factors"] for fixture in training]
        row = _record(
            f"filter_{weight_two}_{pivot}_{excess}_{batch}",
            lambda: _complete(arm, training, config, **kwargs),
            expected_training,
            args,
            scope="filter_stopping",
            controls=kwargs,
            details=_complete(arm, training, config, details=True, **kwargs),
        )
        rows.append(row)
        if row["measurement"]["stable"]:
            filter_choices.append(row)
    chosen = min(
        filter_choices, key=lambda row: row["measurement"]["median_seconds"]
    )
    frozen = {
        "training_seed": 329,
        "held_out_seed": 334,
        "held_out_count": 16,
        "config": winner["config"],
        "winner": winner["name"],
        "controls": chosen["controls"],
        "filter_winner": chosen["name"],
        "freeze_unix": time.time(),
    }
    args.output.with_suffix(".frozen.json").write_text(
        json.dumps(frozen, indent=2) + "\n"
    )
    held_out = _corpus(frozen["held_out_seed"], frozen["held_out_count"])
    expected = [fixture["factors"] for fixture in held_out]
    for name in ("before", "after", "after_roots", "exhaustive"):
        arm = arms["before"] if name == "before" else arms["after"]
        values = dict(frozen["config"])
        if name in ("before", "after_roots"):
            values["division"] = "roots"
        config = _config(arm, **values)
        controls = dict(frozen["controls"], exhaustive=name == "exhaustive")
        rows.append(
            _record(
                "held_out_" + name,
                lambda: _complete(arm, held_out, config, **controls),
                expected,
                args,
                scope="complete_held_out",
                config=asdict(config),
                controls=controls,
                details=_complete(
                    arm, held_out, config, details=True, **controls
                ),
            )
        )
    cold = []
    for name in ("before", "after"):
        for _ in range(9):
            started = time.perf_counter()
            result = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "v2.benchmarks.phase_three_pipeline",
                    "--output",
                    str(args.output),
                    "--cold-worker",
                    "--cold-arm",
                    name,
                ],
                capture_output=True,
                text=True,
                check=True,
                timeout=30,
            )
            sample = json.loads(result.stdout)
            assert sample["factors"] == expected
            sample["arm"] = name
            sample["lifecycle_seconds"] = time.perf_counter() - started
            cold.append(sample)
    profile = args.output.with_suffix(".prof")
    profiler = cProfile.Profile()
    profiler.runcall(
        lambda: [_complete(arm, training, config) for _ in range(100)]
    )
    profiler.dump_stats(profile)
    cases, control_corpus = _cases()
    control_name, functions, control_expected, _ = cases[-1]
    control = _record(
        control_name,
        functions["m17_runtime"],
        control_expected,
        args,
        scope="unchanged_production_control",
        seed=7,
    )
    args.output.write_text(
        json.dumps(
            {
                "milestone": "M25 / P3.3",
                "environment": environment(),
                "command": sys.orig_argv,
                "training": training,
                "held_out": held_out,
                "freeze": frozen,
                "results": rows,
                "cold_samples": cold,
                "complete_factorization_control": control,
                "control_inputs": control_corpus["complete_inputs"],
                "production_m17_hashes_unchanged": production_unchanged,
                "cold_median_seconds": statistics.median(
                    sample["lifecycle_seconds"] for sample in cold
                ),
                "profile": str(profile),
                "profile_sha256": hashlib.sha256(
                    profile.read_bytes()
                ).hexdigest(),
                "profile_policy": "diagnostic; instrumented JIT differs",
                "budgets": dict(
                    work=200_000_000,
                    wall_seconds=10,
                    cpu_seconds=10,
                    owned_bytes=32 * 1024 * 1024,
                ),
                "baseline_snapshot": str(
                    ROOT / "audit/inputs/m25_p33_before_sources.json"
                ),
                "scope": "fixed QS; no SIQS or dispatcher promotion",
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
