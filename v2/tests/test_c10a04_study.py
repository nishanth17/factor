"""Independent count arithmetic and interrupted-study contracts."""

import json
import tempfile
import unittest
from dataclasses import asdict
from math import prod
from pathlib import Path
from unittest.mock import patch

from v2.benchmarks.qs.c10a04 import configuration, selection, study
from v2.benchmarks.qs.c10a04.accounting import Ledger, digest, receipt
from v2.execution.budget import BudgetExhaustedError
from v2.qs import SIQSJob
from v2.qs.factor_base import build_factor_base
from v2.qs.multiplier import select_multiplier
from v2.tests.resource_stops import InterruptBudget as Budget


def allowance():
    return Budget(work_limit=10**9, seconds=30, cpu_seconds=30)


def limits():
    return dict(
        total=dict(wall=100, cpu=100, work=1000),
        train=dict(wall=100, cpu=100, work=1000),
        disk_bytes=2**20 * 8,
    )


class CountCalibrationTests(unittest.TestCase):
    def test_count_endpoint_matches_exhaustive_residues(self):
        # Trial division and residue enumeration share neither the production
        # prime sieve nor its Euler/Tonelli root algorithms.
        primes = [p for p in range(2, 100) if all(p % d for d in range(2, p))]
        for other in (1009, 1013, 1019, 1031):
            n = 1009 * other
            for multiplier in (1, 2, 3, 5, 7, 15):
                eligible = [
                    p
                    for p in primes
                    if any(x * x % p == n * multiplier % p for x in range(p))
                ]
                for count in (1, 3, len(eligible)):
                    with self.subTest(n=n, k=multiplier, count=count):
                        bound, divisor = configuration.endpoint_for_count(
                            n, multiplier, count, allowance(), ceiling=100
                        )
                        self.assertIsNone(divisor)
                        self.assertEqual(bound, eligible[count - 1] + 1)
                        built = build_factor_base(
                            n,
                            multiplier=multiplier,
                            bound=bound,
                            budget=allowance(),
                        )
                        self.assertEqual(
                            list(built.factor_base.primes), eligible[:count]
                        )

    def test_gcd_refusal_and_finite_budget(self):
        bound, divisor = configuration.endpoint_for_count(
            77, 7, 10, allowance(), ceiling=100
        )
        self.assertIsNone(bound)
        self.assertEqual(divisor, 7)

        with self.assertRaises(ValueError):
            configuration.endpoint_for_count(77, 77, 1, allowance())
        with self.assertRaises(ValueError):
            configuration.endpoint_for_count(
                1009 * 1013, 1, 100, allowance(), ceiling=10
            )
        with self.assertRaises(BudgetExhaustedError):
            configuration.endpoint_for_count(
                1009 * 1013, 1, 10, Budget(work_limit=0)
            )
        with self.assertRaises(ValueError):
            select_multiplier(1009 * 1013, candidates=(1, 2))

    def test_accepted_collector_and_pacing_are_fixed(self):
        for digits in configuration.BANDS:
            control = configuration.template(digits, "control")
            for arm in ("compact", "wide"):
                candidate = configuration.template(digits, arm)
                self.assertEqual(
                    asdict(control.collector), asdict(candidate.collector)
                )
                for name in (
                    "memory_bytes",
                    "checkpoint_bytes",
                    "solve_policy",
                    "row_excess",
                    "filter_row_growth",
                    "assignment_policy",
                    "family_count",
                    "max_stalled",
                    "max_trivial",
                ):
                    self.assertEqual(
                        getattr(control, name), getattr(candidate, name)
                    )
                self.assertEqual(candidate.gray_limit, 64)
                self.assertEqual(candidate.pool_size, 64)

    def test_exact_complete_and_unresolved_validation(self):
        fixture = dict(id="test", digits=55, n=10403, factors=[101, 103])
        row = study.run_one(fixture, 7, "control", 1)
        self.assertTrue(row["complete"])
        self.assertEqual(prod(row["factors"]), fixture["n"])
        self.assertEqual(row["certainty"], ["proven_prime", "proven_prime"])

        partial = dict(
            row,
            complete=False,
            divisor=None,
            factors=[],
            certainty=[],
            remaining=[fixture["n"]],
        )
        study.validate(partial, fixture)
        with self.assertRaises(ValueError):
            study.validate(dict(partial, remaining=[77]), fixture)
        with self.assertRaises(ValueError):
            study.validate(
                dict(row, certainty=["probable_prime"] * 2), fixture
            )

    def test_configuration_refusal_preserves_an_unresolved_outcome(self):
        fixture = dict(id="test", digits=55, n=10403, factors=[101, 103])
        with patch.object(
            study,
            "resolve",
            side_effect=configuration.ConfigurationRefusedError(
                "unreachable count"
            ),
        ):
            row = study.run_one(fixture, 7, "compact", 1)
        self.assertFalse(row["complete"])
        self.assertEqual(row["remaining"], [fixture["n"]])
        self.assertEqual(row["reason"], "configuration_refusal")

    def test_new_joint_config_checkpoint_preserves_charges(self):
        corpus = json.loads(study.TRAINING.read_text())
        fixture = next(f for f in corpus["fixtures"] if f["digits"] == 55)
        budget = allowance()
        config, divisor = configuration.resolve(
            fixture["n"], 55, "compact", budget
        )
        self.assertIsNone(divisor)
        job = SIQSJob(fixture["n"], seed=29, config=config, budget=budget)
        job.run(max_blocks=1)
        snapshot = json.loads(json.dumps(job.checkpoint()))
        used = budget.used
        restored_budget = allowance()
        restored = SIQSJob.from_checkpoint(
            snapshot, budget=restored_budget, config=config
        )
        result = restored.run(max_blocks=1)

        self.assertGreater(restored_budget.used, used)
        self.assertGreater(restored_budget.prior_wall, 0)
        self.assertEqual((result.divisor or 1) * result.cofactor, fixture["n"])
        self.assertEqual(restored.config, config)
        self.assertEqual(restored.base.primes, job.base.primes)


class SelectionUncertaintyTests(unittest.TestCase):
    def test_one_input_completion_gain_remains_uncertain(self):
        rows = []
        for sample in range(9):
            for arm in ("control", "compact"):
                rows.append(
                    dict(
                        id=f"t-{sample % 3}",
                        sample=sample,
                        arm=arm,
                        complete=(arm == "compact" or sample % 3 != 0),
                    )
                )
        low, high = selection.completion_interval(rows, "compact")
        self.assertEqual(low, 0)
        self.assertGreater(high, 0)
        for row in rows:
            if row["arm"] == "control":
                row["complete"] = False
        low, high = selection.completion_interval(rows, "compact")
        self.assertEqual((low, high), (1, 1))


class ReceiptTests(unittest.TestCase):
    def test_idempotency_conflicts_and_tampered_receipts(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            ledger = Ledger(output, "source", limits())
            identity = dict(arm="control")
            calls = []

            def execute():
                calls.append(1)
                return dict(work=7)

            first = ledger.call("train", "one", 5, 10, identity, execute)
            used = ledger.used()
            second = ledger.call("train", "one", 5, 10, identity, execute)
            self.assertEqual(first, second)
            self.assertEqual(len(calls), 1)
            self.assertEqual(ledger.used(), used)
            with self.assertRaises(ValueError):
                receipt(output / "one.json", {})
            with self.assertRaises(ValueError):
                ledger.call("train", "one", 5, 10, {}, execute)
            (output / "one.json").write_text("{}")
            with self.assertRaises(ValueError):
                Ledger(output, "source", limits())

    def test_crash_after_publication_settles_once(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            ledger = Ledger(output, "source", limits())
            identity = dict(arm="control")
            ledger.call("train", "one", 5, 10, identity, lambda: dict(work=7))
            row = ledger.data["calls"][0]
            row.update(status="reserved", wall=5, cpu=5, work=10)
            ledger.write()
            recovered = Ledger(output, "source", limits())
            self.assertEqual(recovered.data["calls"][0]["status"], "complete")
            self.assertEqual(
                recovered.data["calls"][0]["sha256"],
                digest(output / "one.json"),
            )
            self.assertEqual(recovered.used(field="work"), 7)

    def test_lost_measured_call_stays_censored_without_retry(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            ledger = Ledger(output, "source", limits())
            identity = dict(
                sample=0,
                digits=55,
                fixture="test",
                n=77,
                seed=7,
                arm="control",
                cap=5,
            )
            with self.assertRaises(KeyboardInterrupt):
                ledger.call(
                    "train",
                    "one",
                    6,
                    10,
                    identity,
                    lambda: (_ for _ in ()).throw(KeyboardInterrupt()),
                )
            resumed = Ledger(output, "source", limits())
            row = resumed.call(
                "train",
                "one",
                6,
                10,
                identity,
                lambda: self.fail("must not retry a lost sample"),
            )
            self.assertFalse(row["complete"])
            self.assertEqual(row["remaining"], [77])
            self.assertEqual(row["reason"], "interrupted")
            self.assertEqual(resumed.used(field="work"), 10)
            self.assertEqual(row["capped_seconds"], resumed.used())

    def test_refusal_and_failed_correctness_are_not_reset(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            ledger = Ledger(output, "source", limits())
            with self.assertRaises(RuntimeError):
                ledger.call("train", "too-large", 101, 1, {}, lambda: {})
            self.assertEqual(ledger.data["calls"], [])
            with self.assertRaises(ValueError):
                ledger.call(
                    "train",
                    "bad",
                    5,
                    10,
                    dict(sample=0),
                    lambda: (_ for _ in ()).throw(
                        ValueError("invalid reconstruction")
                    ),
                )
            self.assertEqual(ledger.used(field="work"), 10)
            with self.assertRaises(RuntimeError):
                Ledger(output, "source", limits())

    def test_repeated_yields_and_partial_group_resume(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            fixtures = [
                dict(id=f"t-{i}", n=prod(factors), factors=factors, digits=55)
                for i, factors in enumerate(([7, 11], [7, 13], [11, 13]))
            ]
            corpus = dict(
                fixtures=fixtures,
                certificates={
                    "7": dict(kind="trial"),
                    "11": dict(kind="trial"),
                    "13": dict(kind="trial"),
                },
            )
            path = output / "corpus.json"
            path.write_text(json.dumps(corpus))
            control = dict(caps={"55": 1}, arm_order_seed=1)
            calls = []

            def execute(fixture, seed, arm, _cap):
                calls.append(arm)
                if len(calls) in (2, 29):
                    raise KeyboardInterrupt()
                return dict(
                    id=fixture["id"],
                    digits=55,
                    seed=seed,
                    arm=arm,
                    complete=True,
                    divisor=fixture["factors"][0],
                    factors=fixture["factors"],
                    remaining=[],
                    certainty=["proven_prime"] * 2,
                    work=1,
                    stats={},
                    reason="factor_found",
                    seconds=0.1,
                    cpu_seconds=0.1,
                    capped_seconds=0.1,
                    cap_seconds=1,
                )

            with (
                patch.object(study, "load", return_value=control),
                patch.object(study, "CONTROL", path),
                patch.object(study, "TRAINING", path),
                patch.object(study, "warm"),
                patch.object(study, "cold", return_value=dict(work=0)),
                patch.object(study, "run_one", side_effect=execute),
            ):
                study.measure(output, "train", 55, 0)
                study.measure(output, "train", 55, 0)
                self.assertEqual(
                    len(list(output.glob("train-55-yield-*.json"))), 2
                )
                with self.assertRaises(KeyboardInterrupt):
                    study.measure(output, "train", 55, 1000)
                first = next(output.glob("train-55-00-*.json"))
                saved = first.read_bytes()
                with self.assertRaises(KeyboardInterrupt):
                    study.measure(output, "train", 55, 1000)
                marker = output / "train-55-stability-9.json"
                saved_marker = marker.read_bytes()
                study.measure(output, "train", 55, 1000)
                self.assertEqual(marker.read_bytes(), saved_marker)
                self.assertEqual(first.read_bytes(), saved)
                self.assertEqual(len(calls), 54)
                self.assertEqual(len(study.rows(output, "train", 55)), 54)
                self.assertEqual(
                    sum(
                        r["reason"] == "interrupted"
                        for r in study.rows(output, "train", 55)
                    ),
                    2,
                )


if __name__ == "__main__":
    unittest.main()
