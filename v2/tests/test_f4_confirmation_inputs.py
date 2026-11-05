"""Independent integrity checks for immutable F4 loader inputs."""

import hashlib
import json
import math
import unittest
from pathlib import Path

from v2.benchmarks.suites.build_phase_two_corpus import verify_certificates

INPUTS = Path(__file__).resolve().parents[1] / "benchmarks" / "inputs"


def read(relative):
    return json.loads((INPUTS / relative).read_text())


class ConfirmationInputTests(unittest.TestCase):
    def test_frozen_bindings_and_no_hidden_dispatch_labels(self):
        plan_path = INPUTS / "controls/b12_f4_confirmation_v1.json"
        plan = json.loads(plan_path.read_text())
        profile_path = INPUTS / "controls/b12_f4_selection_v1.json"
        profile = json.loads(profile_path.read_text())
        corpus = read("corpora/b12_f4_confirmation_v1.json")

        self.assertEqual(
            hashlib.sha256(profile_path.read_bytes()).hexdigest(),
            plan["profile_sha256"],
        )
        self.assertEqual(
            hashlib.sha256(plan_path.read_bytes()).hexdigest(),
            corpus["policy_freeze_sha256"],
        )
        self.assertEqual(len(profile["rows"]), 33)
        for row in profile["rows"]:
            self.assertEqual(
                set(row), {"digits", "options", "tiers", "cold", "warm"}
            )
        values = [r["n"] for r in corpus["fixtures"] + corpus["uniform"]]
        self.assertEqual(len(values), len(set(values)))
        self.assertFalse(set(values) & set(plan["exclusion_moduli"]))

    def test_certificates_reconstruction_and_independent_selected_moduli(self):
        corpus = read("corpora/b12_f4_confirmation_v1.json")
        verify_certificates(corpus["certificates"])
        primes = []
        for row in corpus["fixtures"]:
            self.assertEqual(
                math.prod(p**e for p, e in row["factors"]), row["n"]
            )
            primes.extend(p for p, _ in row["factors"])
        self.assertEqual(len(primes), len(set(primes)))
        selected = [
            row
            for row in corpus["fixtures"]
            if row["population"] == "selected_ecm_enriched"
        ]
        self.assertEqual(len(selected), 12)
        self.assertEqual(
            {(r["digits"], r["smaller_digits"]) for r in selected},
            {
                (digits, smaller)
                for digits in (60, 80, 100)
                for smaller in (15, 17, 18, 20)
            },
        )
        self.assertEqual(len(corpus["uniform"]), 52)
        for digits in (30, 40, 50, 60, 70, 80, 100):
            rows = [r for r in corpus["uniform"] if r["digits"] == digits]
            self.assertEqual(len(rows), 12 if digits in (60, 80, 100) else 4)
            self.assertTrue(all(len(str(r["n"])) == digits for r in rows))

    def test_portable_math_configs_bind_the_recorded_profile(self):
        reviewed = read("controls/b12_f4_selection_reviewed_v1.json")
        profile = read("controls/b12_f4_selection_v1.json")
        keys = set()
        for row in reviewed["rows"]:
            options = dict(row["config"])
            options.update(
                siqs=None, sss=None, allocation=None, memory_bytes=64 * 2**20
            )
            encoded = json.dumps(
                options,
                sort_keys=True,
                separators=(",", ":"),
                allow_nan=False,
            )
            key = row["digits"], hashlib.sha256(encoded.encode()).hexdigest()
            self.assertNotIn(key, keys)
            keys.add(key)
            found = [
                r
                for r in profile["rows"]
                if (r["digits"], r["options"]) == key
            ]
            self.assertEqual(len(found), 1)


if __name__ == "__main__":
    unittest.main()
