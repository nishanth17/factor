"""Finite neighborhoods, compatible ablations and uncensored ranking."""

import copy
import json
import unittest
from collections import Counter
from math import gcd, prod

from v2.benchmarks.qs.a7a.cell import settings
from v2.benchmarks.qs.a7a.combined import corpus_identity, fixtures
from v2.benchmarks.qs.a7a.protocol import INPUTS
from v2.benchmarks.qs.a7a.regime_search import (
    BANDS,
    COUNTS,
    ablations,
    candidate,
    cells,
    coarse,
    refinements,
    shortlist,
)


class RegimeSearchTests(unittest.TestCase):
    def test_regime_panel_is_certified_disjoint_and_split_before_tuning(self):
        from v2.benchmarks.suites.build_phase_two_corpus import (
            verify_certificates,
        )

        data = json.loads((INPUTS / "regime_corpus.json").read_text())
        verify_certificates(data["certificates"])
        self.assertEqual(
            Counter(
                (row["purpose"], row["digits"]) for row in data["fixtures"]
            ),
            Counter(
                {
                    (purpose, digits): 3
                    for purpose in ("calibration", "confirmation")
                    for digits in BANDS
                }
            ),
        )
        previous = [r["n"] for r in fixtures("all").values() if r["n"] > 1]
        for row in data["fixtures"]:
            self.assertEqual(prod(row["factors"]), row["n"])
            self.assertEqual(len(str(row["n"])), row["digits"])
            self.assertLessEqual(
                abs(row["factor_bits"][0] - row["factor_bits"][1]), 1
            )
            self.assertTrue(
                all(gcd(row["n"], value) == 1 for value in previous)
            )
            previous.append(row["n"])
        self.assertNotEqual(corpus_identity("all"), corpus_identity("regimes"))
        self.assertEqual(len(fixtures("regimes")), len(fixtures("all")) + 54)

    def test_coarse_search_covers_counts_packets_and_serial_for_each_band(
        self,
    ):
        candidates = coarse()
        self.assertEqual(len(candidates), 1 + 2 * len(COUNTS))
        for digits in BANDS:
            fixture = dict(id=f"band_{digits}", digits=digits, factors=[2, 3])
            rows = cells(fixture, candidates)
            self.assertEqual(
                {r["arm"] for r in rows}, {"S0", *(f"C{n}" for n in COUNTS)}
            )
            altered = dict(fixture, factors=[5, 7])
            self.assertEqual(rows, cells(altered, candidates))
            self.assertTrue(all(r["memory"] == 22 * 2**30 for r in rows))

    def test_refinements_and_all_removals_remain_compatible_and_independent(
        self,
    ):
        selected = candidate(8, features=dict(coefficient_cache=True))
        alternatives = refinements(selected) + ablations(
            selected, "presieve+unroll4+parity_two+packed_queue+wide_hits"
        )
        before = copy.deepcopy(selected)
        for choice in alternatives:
            worker = settings(
                choice["arm"],
                22 * 2**30,
                False,
                choice["polynomials_per_task"],
                2048,
                choice["features"],
            )
            self.assertLessEqual(worker.lookahead, 4 * worker.workers)
        self.assertEqual(selected, before)
        self.assertEqual(
            {
                int(r["arm"][1:])
                for r in refinements(selected)
                if r["name"].startswith("w")
            },
            {7, 8, 9},
        )
        removed = next(
            r
            for r in alternatives
            if r["name"] == "without_lease_and_private_stop"
        )
        self.assertFalse(removed["features"]["private_stop"])
        self.assertEqual(removed["features"]["stop_backend"], "event")

    def test_incomplete_factors_and_partial_prefixes_never_win_by_short_time(
        self,
    ):
        def row(name, seconds, **changes):
            values = dict(
                status="ok",
                variant=name,
                seconds=seconds,
                candidate=candidate(4),
                endpoint="complete",
                complete=True,
            )
            values.update(changes)
            return values

        rows = [
            row("censored", 0.1, complete=False),
            row("winner", 2),
            row("near_tie", 2.1),
            row("slow", 8),
            row("serial", 9, candidate=candidate(0)),
        ]
        winners = shortlist(rows, target=1)
        self.assertEqual(
            {r["variant"] for r in winners}, {"winner", "near_tie", "serial"}
        )
        self.assertEqual(shortlist([row("censored", 0.1, complete=False)]), [])
        prefix = row("partial", 0.1, endpoint="prefix", coverage_valid=False)
        self.assertEqual(shortlist([prefix]), [])
        with self.assertRaises(ValueError):
            shortlist(rows + [prefix])


if __name__ == "__main__":
    unittest.main()
