"""Protocol identity, parameter coverage and expanded replay interchange."""

import json
import tempfile
import unittest
from dataclasses import asdict, replace
from pathlib import Path
from unittest.mock import patch

from v2 import qs
from v2.benchmarks.qs.c9a import combined, extended
from v2.benchmarks.qs.c9a.__main__ import execution_cells
from v2.benchmarks.qs.c9a.accounting import StudyLedger
from v2.benchmarks.qs.c9a.configuration import configuration, decode
from v2.benchmarks.qs.c9a.confirmation import input_prefix
from v2.benchmarks.qs.c9a.journal import atom_record, canonical
from v2.benchmarks.qs.c9a.plan import limits
from v2.benchmarks.qs.c9a.preparation import INPUTS
from v2.tests.test_c9a_extended import atom_at, collector


class ExpandedProtocolTests(unittest.TestCase):
    def test_initial_stop_continuation_and_exact_receipt_reuse(self):
        corpus = json.loads((INPUTS / "corpus.json").read_text())
        resolved = json.loads((INPUTS / "structural.json").read_text())[
            "resolved"
        ]
        plan = extended.build(corpus, resolved)
        initial, suffix = execution_cells(plan, "A", False, {}, None)
        self.assertEqual((len(initial), suffix), (4, "-initial"))
        with patch(
            "v2.benchmarks.qs.c9a.__main__.outcomes", return_value=[{}] * 3
        ):
            with self.assertRaises(RuntimeError):
                execution_cells(
                    plan, "A", True, {"feasibility_review": "checked"}, None
                )
        with patch(
            "v2.benchmarks.qs.c9a.__main__.outcomes", return_value=[{}] * 4
        ):
            with self.assertRaises(RuntimeError):
                execution_cells(plan, "A", True, {}, None)
            continued, suffix = execution_cells(
                plan, "A", True, {"feasibility_review": "checked"}, None
            )
        self.assertEqual(continued[:4], initial)
        self.assertEqual((len(continued), suffix), (32, ""))
        with tempfile.TemporaryDirectory() as directory:
            ledger = StudyLedger(Path(directory), "test-freeze", limits(plan))
            cell = initial[0]
            result = dict(
                work=0,
                measured_work=0,
                warm_work=0,
                observer_work=0,
                wall=0,
                cpu=0,
            )
            ledger.call_cell("A", "first", cell, lambda: result)
            calls = len(ledger.data["calls"])
            with patch(
                "builtins.print", side_effect=AssertionError("repeated work")
            ):
                self.assertEqual(
                    ledger.call_cell(
                        "A", "first", cell, lambda: print("rerun")
                    ),
                    result,
                )
            self.assertEqual(len(ledger.data["calls"]), calls)
            with self.assertRaises((ValueError, RuntimeError)):
                ledger.call_cell(
                    "A",
                    "first",
                    dict(cell, seed=cell["seed"] + 1),
                    lambda: result,
                )

    def test_primary_plan_keeps_four_cell_stop_and_old_identity(self):
        corpus = json.loads((INPUTS / "corpus.json").read_text())
        resolved = json.loads((INPUTS / "structural.json").read_text())[
            "resolved"
        ]
        old = combined.build(corpus, resolved)
        plan = extended.build(corpus, resolved)
        self.assertEqual(extended.validate(plan)["A"], 32)
        self.assertEqual(plan["initial_cells"], 4)
        self.assertNotEqual(plan["identity"], old["identity"])
        self.assertEqual({cell["arm"] for cell in plan["cells"]}, {"C0", "C1"})
        for cell in plan["cells"]:
            restored = decode(cell["configuration"], qs)
            self.assertEqual(asdict(restored), cell["configuration"])
            if cell["arm"] == "C0":
                self.assertEqual(restored, configuration(cell["digits"], "C0"))
            if cell["digits"] == 80:
                self.assertEqual(cell["kind"], "prefix")

    def test_joint_profiles_separate_reporting_and_slp_and_populate_three_lp(
        self,
    ):
        incumbent = configuration(60, "C0").collector
        for bound in (30000, 103814, 238724, 600000):
            choices = extended.profiles(bound, incumbent)
            self.assertEqual(len(choices), 6)
            self.assertNotEqual(
                choices["fixed-report"]["candidate_bound"],
                choices["wide-report"]["candidate_bound"],
            )
            self.assertLessEqual(
                choices["wide-report"]["residual_bound"],
                choices["wide-report"]["large_prime_bound"],
            )
            self.assertGreater(
                choices["wide-slp"]["residual_bound"],
                choices["wide-slp"]["large_prime_bound"],
            )
            self.assertGreaterEqual(
                choices["three-prime"]["large_product_bound"], 8 * bound**3
            )
            for values in choices.values():
                base = asdict(incumbent)
                base.update(values)
                qs.ExtendedYieldSieveConfig(**base)

    def test_unobserved_wider_domain_is_not_ranked_as_free_yield(self):
        choices = extended.profiles(30000, configuration(60, "C0").collector)
        observed = choices["fixed-report"]["candidate_bound"]
        rows = [
            dict(
                residual=15,
                primes=[3, 5],
                report_requirement=15,
                complete_cost=100,
                useful_credit=1,
            )
        ]
        ranked = extended.rank_profiles(
            choices, rows, observed_report_bound=observed
        )
        by_name = {row["name"]: row for row in ranked}
        self.assertEqual(
            by_name["three-prime"]["coverage"], "unobserved-width"
        )
        self.assertFalse(by_name["three-prime"]["dominated"])
        with self.assertRaises(ValueError):
            extended.rank_profiles(
                choices,
                [dict(rows[0], complete_cost=-1)],
                observed_report_bound=observed,
            )

    def test_integration_identity_binds_the_same_underlying_configuration(
        self,
    ):
        cell = dict(
            digits=60, arm="C0", fixture="sample", resolved=dict(bound=30000)
        )
        plan = dict(cells=[cell, dict(cell, arm="C1")])
        incumbent = configuration(60, "C0")
        result = extended.integrate(plan, "1" * 40, {"sample": incumbent})
        self.assertEqual(
            result["control_source"], result["candidate_base_source"]
        )
        self.assertEqual(
            result["cells"][0]["configuration"], asdict(incumbent)
        )
        self.assertNotEqual(result["identity"], extended.IDENTITY)
        self.assertNotIn("protocol", plan["cells"][0])
        self.assertEqual({c["source"] for c in result["cells"]}, {"candidate"})
        self.assertIn("60", result["integrated_incumbents"])
        self.assertNotEqual(
            input_prefix(dict(protocol=result["identity"])),
            input_prefix(dict(protocol=extended.IDENTITY)),
        )

    def test_geometry_candidates_preserve_finite_count_and_interval_domains(
        self,
    ):
        incumbent = configuration(60, "C0")
        variants = extended.geometry_candidates(incumbent)
        self.assertGreater(
            variants["wider-base"].base_count, incumbent.base_count
        )
        self.assertLess(
            variants["wider-base"].base_count, incumbent.base_bound
        )
        self.assertGreater(
            variants["wider-interval"].half_width, incumbent.half_width
        )
        near_ceiling = replace(incumbent, base_count=incumbent.base_bound - 1)
        self.assertEqual(
            extended.geometry_candidates(near_ceiling)["wider-base"],
            near_ceiling,
        )

    def test_opaque_journal_has_an_explicit_contract(self):
        worker = collector()
        record = json.loads(
            canonical(atom_record(atom_at(worker, 1, opaque=True)))
        )
        self.assertEqual(record["residual_contract"], "extended-v1")
        self.assertEqual(record["opaque"], 15)
        self.assertEqual(record["large_primes"], [])


if __name__ == "__main__":
    unittest.main()
