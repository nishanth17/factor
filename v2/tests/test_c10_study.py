"""C10 comparisons isolate policy changes and preserve unsuccessful work."""

import tempfile
import unittest
from dataclasses import asdict
from pathlib import Path

from v2.benchmarks.qs.c10 import selection, study
from v2.qs import GrowingDoubleLargeSieveConfig, SIQSConfig


def row(arm, sample, seconds, complete=True):
    return dict(
        id=f"input-{sample % 3}",
        sample=sample,
        seed=study.SEEDS[(sample // 3) % 3],
        arm=arm,
        seconds=seconds,
        capped_seconds=seconds if complete else 60,
        complete=complete,
    )


class C10StudyTests(unittest.TestCase):
    def test_arms_change_only_one_declared_policy(self):
        config = SIQSConfig(
            assignment_policy="flyer",
            collector=GrowingDoubleLargeSieveConfig(
                large_prime_bound=1000,
                large_product_bound=1000000,
            ),
        )
        arms = {a: asdict(c) for a, c in study.configurations(config).items()}
        quality = arms["quality"]
        self.assertEqual(quality.pop("assignment_policy"), "quality")
        control = dict(arms["control"])
        control.pop("assignment_policy")
        self.assertEqual(quality, control)
        plan = arms["plan"]
        self.assertEqual(plan["collector"].pop("power_plan_bytes"), 2**21)
        control = arms["control"]
        control["collector"].pop("power_plan_bytes")
        self.assertEqual(plan, control)

    def test_completions_outrank_cheaper_unsuccessful_work(self):
        rows = []
        for i in range(9):
            rows += [
                row("control", i, 20),
                row("plan", i, 1, False),
                row("quality", i, 21),
            ]
        self.assertEqual(selection.rank(rows)[0]["arm"], "control")
        report = selection.comparison(rows, "plan")
        self.assertFalse(report["passes_scoped_gate"])
        self.assertEqual(report["complete_challenger"], 0)
        self.assertEqual(report["cost_change"], 2)

    def test_small_stable_effect_has_no_arbitrary_promotion_floor(self):
        rows = [
            row(a, i, 10 if a == "control" else 9.99)
            for i in range(9)
            for a in ("control", "plan")
        ]
        result = selection.comparison(rows, "plan")
        self.assertEqual(result["independent_inputs"], 3)
        self.assertTrue(result["passes_scoped_gate"])
        rows[-1]["capped_seconds"] = 11
        self.assertTrue(selection.instability(rows))
        self.assertFalse(
            selection.comparison(rows, "plan")["passes_scoped_gate"]
        )

    def test_failed_calls_remain_charged_and_phase_cap_precedes_call(self):
        with tempfile.TemporaryDirectory() as directory:
            expense = study.expenses(Path(directory))

            def fail():
                raise ValueError("interrupted")

            with self.assertRaises(ValueError):
                study.charged(expense, "train", "failure", 20, fail)
            self.assertGreaterEqual(study.phase_used(expense, "train"), 20)
            with self.assertRaises(RuntimeError):
                study.charged(expense, "train", "too-large", 5400, fail)
            restored = study.expenses(Path(directory))
            self.assertGreaterEqual(restored.used(), 920)


if __name__ == "__main__":
    unittest.main()
