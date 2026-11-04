"""Independent receipt oracles; no ECM execution or worker startup."""

import copy
import unittest

from v2.benchmarks.ecm.b12_f4.headroom import summarise


class HeadroomTests(unittest.TestCase):
    def fixture(self):
        cells = [
            dict(id=arm, arm=arm, input="n", mode="warm", workers=workers)
            for arm, workers in (("M0", 0), ("P2", 2))
        ]
        freeze = dict(
            identity="test",
            cells=cells,
            candidate_root="candidate",
            control_root="control",
            inputs={"n": dict(n=35)},
            config=dict(ecm_tiers=[[2, 10, 1]]),
        )
        records = []
        for cell in cells:
            duration = 4 if cell["arm"] == "M0" else 2
            result = dict(
                id="curve",
                sigma=19,
                factor=5,
                reason="complete",
                phase="stage_two",
                chain_chunks=1,
            )
            child = dict(
                complete=True,
                n=35,
                reconstruct=35,
                execution_checkpoint={
                    "search": {"payload": {"completed": {"0": result}}}
                },
                warmup_seconds=3,
                warm_call_seconds=duration,
                startup_and_shutdown_stats=[dict(parent_cpu=1)],
                aggregate_cpu=99,
                pool_lifetime_stats={},
                parent_peak_rss=0,
            )
            records.append(
                dict(
                    identity="test",
                    cell=cell,
                    source_root="candidate",
                    child=child,
                    timed_out=False,
                    returncode=0,
                    failure=None,
                    wall=100,
                )
            )
        return freeze, records

    def test_warm_ratio_uses_call_duration_and_keeps_whole_cell_cpu(self):
        freeze, records = self.fixture()
        result = summarise(freeze, records)
        self.assertEqual(result["comparisons"][0]["exploratory_ratio"], 2)
        self.assertEqual(result["rows"][0]["whole_cell_observed_cpu"], 99)
        self.assertEqual(result["rows"][0]["call_parent_cpu"], 1)

    def test_bad_divisor_missing_curve_and_changed_search_are_rejected(self):
        freeze, records = self.fixture()
        for field, value in (("factor", 6), ("sigma", 23)):
            changed = copy.deepcopy(records)
            changed[1]["child"]["execution_checkpoint"]["search"]["payload"][
                "completed"
            ]["0"][field] = value
            with self.assertRaises(ValueError):
                summarise(freeze, changed)
        records[1]["child"]["execution_checkpoint"]["search"]["payload"][
            "completed"
        ] = {}
        with self.assertRaises(ValueError):
            summarise(freeze, records)

    def test_invalid_durations_cannot_create_a_speedup(self):
        freeze, records = self.fixture()
        for duration in (0, -1, float("nan"), float("inf")):
            changed = copy.deepcopy(records)
            changed[1]["child"]["warm_call_seconds"] = duration
            with self.assertRaises(ValueError):
                summarise(freeze, changed)

    def test_incomplete_and_missing_cells_stay_visible(self):
        freeze, records = self.fixture()
        records[0]["child"]["complete"] = False
        result = summarise(freeze, records[:1])
        self.assertEqual(result["failures"], ["M0"])
        self.assertEqual(result["missing"], ["P2"])
        self.assertEqual(result["comparisons"], [])
