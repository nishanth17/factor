"""Independent-process sampling and frozen-input contracts for the bakeoff."""

import unittest

from v2.benchmarks import b4_bakeoff as bakeoff


class BakeoffContracts(unittest.TestCase):
    def test_counterbalanced_positions(self):
        arms = ("baseline", "one", "two", "three", "four", "five")
        positions = []
        for block in range(6):
            order = bakeoff.block_order(arms, "python-int", block, 17)
            self.assertEqual(set(order), set(arms))
            positions.append(order.index("baseline"))
        self.assertEqual(sorted(positions), list(range(6)))

    def test_paired_interval_and_sample_floor(self):
        gain, interval = bakeoff.paired_gain([100] * 9, [90] * 9)
        self.assertAlmostEqual(gain, 10)
        for endpoint in interval:
            self.assertAlmostEqual(endpoint, 10)
        with self.assertRaises(ValueError):
            bakeoff.paired_gain([100] * 8, [90] * 8)
        with self.assertRaises(ValueError):
            bakeoff.paired_gain([100] * 9, [90] * 10)

    def test_selection_rejects_regression_noise_and_loss(self):
        def row(arm, gain, stable=True, completion_ok=True):
            return dict(
                backend="python-int",
                arm=arm,
                reduction_percent=gain,
                stable=stable,
                completion_ok=completion_ok,
            )

        rows = [
            row("small_gain", 2),
            row("noise", 20, stable=False),
            row("regression", 30, completion_ok=False),
            row("loss", -5),
        ]
        self.assertEqual(bakeoff.select(rows), {"python-int": "small_gain"})
        self.assertEqual(
            bakeoff.select([row("loss", -1)]), {"python-int": None}
        )

    def test_comparison_checks_cross_process_work(self):
        captures = []
        for block in range(9):
            for arm, seconds in (("baseline", 100), ("candidate", 90)):
                captures.append(
                    dict(
                        backend="python-int",
                        arm=arm,
                        block=block,
                        seconds=seconds,
                        cpu_seconds=seconds,
                        rows=[
                            dict(
                                case="64",
                                seconds=seconds,
                                complete=True,
                                work=5,
                            )
                        ],
                    )
                )
        protocol = {"relative_iqr_limit": 0.15}
        report = bakeoff.compare(captures, "python-int", "candidate", protocol)
        self.assertTrue(report["matched_outcomes_and_work"])
        self.assertTrue(report["timing_evidence_passes"])

        captures[-1]["rows"][0]["work"] = 6
        with self.assertRaises(AssertionError):
            bakeoff.compare(captures, "python-int", "candidate", protocol)

    def test_source_inputs_disjoint_and_certified(self):
        protocol, corpus, pins = bakeoff.verify_inputs()
        self.assertEqual(protocol["seeds"], [11, 29, 53])
        self.assertEqual(len(corpus["fixtures"]), 18)
        self.assertEqual(len(pins["arms"]), 6)
        self.assertEqual(
            {fixture["split"] for fixture in corpus["fixtures"]},
            {"screen", "confirmation"},
        )


if __name__ == "__main__":
    unittest.main()
