"""Independent checks of bounded prospective early-stop decisions."""

import unittest

from v2.benchmarks.qs.c10a04.short_decisions import decision, schedule


def outcomes(count, *, compact_seconds=70, control_seconds=100):
    result = []
    for index in range(count):
        for arm, seconds in (
            ("control", control_seconds),
            ("compact", compact_seconds),
        ):
            result.append(
                dict(
                    sample=index,
                    id=f"input-{index % 3}",
                    arm=arm,
                    complete=True,
                    capped_seconds=seconds,
                )
            )
    return result


class ShortStudyDecisionTests(unittest.TestCase):
    def test_first_three_favorable_pairs_do_not_trigger_positive_stop(self):
        self.assertEqual(decision(outcomes(3)), "continue")

    def test_six_pairs_must_improve_every_input_without_a_regression(self):
        rows = outcomes(6)
        self.assertEqual(
            decision(rows), "stop-preliminary-consistent-compact-trend"
        )
        rows[1]["capped_seconds"] = 105
        self.assertEqual(decision(rows), "continue")

        rows = outcomes(6)
        rows[1]["capped_seconds"] = 95
        self.assertEqual(decision(rows), "continue")

    def test_partial_completion_cannot_trigger_positive_stop(self):
        rows = outcomes(6)
        rows[1]["complete"] = False
        self.assertEqual(decision(rows), "continue")

    def test_zero_compact_completion_stops_without_budget_extension(self):
        rows = outcomes(3)
        for row in rows:
            if row["arm"] == "compact":
                row["complete"] = False
        self.assertEqual(decision(rows), "stop-no-compact-completions")

    def test_nine_pairs_end_without_automatic_extension(self):
        self.assertEqual(
            decision(outcomes(9, compact_seconds=110)), "finished-nine-pairs"
        )

    def test_incomplete_or_duplicate_pairs_are_rejected(self):
        with self.assertRaises(ValueError):
            decision(outcomes(3)[:-1])
        with self.assertRaises(ValueError):
            decision(outcomes(3) + [outcomes(3)[0]])

    def test_schedule_covers_three_inputs_and_three_seeds(self):
        corpus = dict(
            fixtures=[
                dict(id=f"input-{index}", digits=60) for index in range(3)
            ]
        )
        planned = schedule(corpus)
        self.assertEqual(len(planned), 9)
        self.assertEqual(
            {(fixture["id"], seed) for _, fixture, seed, _ in planned},
            {
                (f"input-{index}", seed)
                for index in range(3)
                for seed in (7, 29, 47)
            },
        )
        for _, _, _, order in planned:
            self.assertEqual(set(order), {"control", "compact"})
