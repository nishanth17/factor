"""Prospective plan/configuration checks; no campaign cells are executed."""

import json
import unittest

from v2 import qs
from v2.benchmarks.infrastructure.resource_control import (
    historical_asdict as asdict,
)
from v2.benchmarks.qs.c9a import combined
from v2.benchmarks.qs.c9a.configuration import configuration, decode
from v2.benchmarks.qs.c9a.journal import atom_record, canonical
from v2.benchmarks.qs.c9a.plan import limits
from v2.benchmarks.qs.c9a.preparation import INPUTS
from v2.tests.test_c9a_combined import collector
from v2.tests.test_qs_dlp import admit


class CombinedPlanTests(unittest.TestCase):
    def test_primary_two_arm_plan_and_finite_feasibility_prefix(self):
        corpus = json.loads((INPUTS / "corpus.json").read_text())
        resolved = json.loads((INPUTS / "structural.json").read_text())[
            "resolved"
        ]
        plan = combined.build(corpus, resolved)
        self.assertEqual(combined.validate(plan)["A"], 32)
        self.assertEqual({c["arm"] for c in plan["cells"]}, {"C0", "C1"})
        self.assertEqual(
            {c["fixture"] for c in plan["cells"][:4]},
            {"exposed_60_0", "exposed_65_0"},
        )
        self.assertEqual(
            sum(
                c["seconds"] + c["margin"] + c["warm_seconds"]
                for c in plan["cells"][:4]
            ),
            410,
        )
        self.assertEqual(limits(plan)["disk_bytes"], 2**31)
        for cell in plan["cells"]:
            decoded = decode(cell["configuration"], qs)
            self.assertEqual(asdict(decoded), cell["configuration"])
            if cell["arm"] == "C0":
                self.assertEqual(decoded, configuration(cell["digits"], "C0"))
            if cell["digits"] == 80:
                self.assertEqual(cell["kind"], "prefix")

    def test_joint_domain_keeps_incumbent_and_has_distinct_boundaries(self):
        previous = configuration(60).collector
        for bound in (30000, 103814, 238724, 400000, 600000):
            values = combined.envelope(bound, previous)
            self.assertGreaterEqual(
                values["large_product_bound"], previous.large_product_bound
            )
            self.assertGreaterEqual(
                values["large_prime_bound"], previous.large_prime_bound
            )
            self.assertGreaterEqual(
                values["residual_bound"], previous.residual_bound
            )
            self.assertGreaterEqual(
                values["small_prime_bound"] ** 2
                + 2 * values["small_prime_bound"]
                + 1,
                previous.large_product_bound,
            )
            self.assertLessEqual(
                values["large_product_bound"],
                values["small_prime_bound"] * values["large_prime_bound"],
            )
            self.assertLess(values["large_product_bound"], bound**3)
        frontier = combined.envelope(600000, previous)
        self.assertEqual(frontier["large_prime_bound"], 9600000)
        self.assertEqual(frontier["small_prime_bound"], 2400000)
        self.assertEqual(frontier["large_product_bound"], 11520000000000)
        self.assertEqual(
            frontier["candidate_bound"], previous.large_product_bound
        )

    def test_compact_atoms_have_plain_replay_interchange(self):
        worker = collector()
        admit(worker, 1, (3, 5))
        atom = next(iter(worker._atoms.values()))
        record = json.loads(canonical(atom_record(atom)))
        self.assertEqual(record["id"], atom.relation_id)
        self.assertEqual(record["exponents"], [[2, 1]])


if __name__ == "__main__":
    unittest.main()
