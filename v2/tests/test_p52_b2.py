"""Frozen B2 controls and protocol gates, without running performance tests."""

import copy
import unittest
from unittest.mock import patch

from v2 import ecm, portfolio
from v2.benchmarks import p52_b2


class PairedProtocolTests(unittest.TestCase):
    def test_independent_corpus_controls_and_admitted_grid(self):
        protocol, corpus = p52_b2.inputs()
        ids = [f["id"] for f in corpus["fixtures"]]
        self.assertEqual(len(ids), len(set(ids)))
        self.assertEqual(len({f["n"] for f in corpus["fixtures"]}), len(ids))
        self.assertEqual(
            {f["split"] for f in corpus["fixtures"]}, {"training", "held_out"}
        )
        for case in protocol["cases"]:
            for arm in protocol["training_arms"]:
                config = portfolio.PortfolioConfig(
                    **p52_b2.options(protocol, case, arm)
                )
                self.assertLess(config.workspace_reserve, config.memory_bytes)
                self.assertEqual(config.backend, "python-int")
        control = p52_b2.load_control()
        self.assertNotIn(
            "ecm_pair_distance", control.PortfolioConfig.__dataclass_fields__
        )
        expected = control.advance_job.__globals__["ecm"].scalar_multiply(
            97, 3, 1, 1009, 2
        )
        with patch.object(
            ecm, "scalar_multiply", side_effect=AssertionError("active ECM")
        ):
            self.assertEqual(
                control.advance_job.__globals__["ecm"].scalar_multiply(
                    97, 3, 1, 1009, 2
                ),
                expected,
            )
        self.assertEqual(protocol["sampling"], [[3, 9], [5, 31], [8, 63]])

    def test_selection_and_uncertainty_use_fixtures_not_reruns(self):
        protocol, _ = p52_b2.inputs()
        captures = []
        for case in protocol["cases"]:
            for index, arm in enumerate(protocol["training_arms"]):
                rows = [
                    dict(
                        fixture="a",
                        seed=seed,
                        complete=True,
                        reason="complete",
                    )
                    for seed in (7, 19, 41)
                ]
                captures.append(
                    dict(
                        case=case,
                        arm=arm,
                        deterministic=True,
                        samples=[
                            dict(seconds=1 / (index + 1), rows=rows)
                            for _ in range(9)
                        ],
                    )
                )
        selected = p52_b2.select(captures, protocol)
        self.assertEqual(set(selected.values()), {"paired_2"})
        self.assertEqual(
            p52_b2.completion_interval(captures[0], captures[1]), [0, 0]
        )
        bad = copy.deepcopy(captures)
        for capture in bad:
            if capture["arm"] == "paired_2":
                capture["deterministic"] = False
        self.assertEqual(
            set(p52_b2.select(bad, protocol).values()), {"paired_1"}
        )
        bad[0]["samples"][0]["rows"][0]["reason"] = "wall_limit"
        self.assertIsNone(p52_b2.select(bad, protocol)["small"])
