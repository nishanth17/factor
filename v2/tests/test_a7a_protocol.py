"""Finite plan sums, endpoint analysis and corpus independence contracts."""

import unittest
from collections import defaultdict

from v2.benchmarks.qs.a7a.analysis import paired_interval, quantile_bound
from v2.benchmarks.qs.a7a.protocol import corpus, schedule


class ProtocolTests(unittest.TestCase):
    def test_combined_panel_and_memory_cover_amended_stages(self):
        from math import gcd, prod

        from v2.benchmarks.infrastructure.resource_control import (
            historical_replace as replace,
        )
        from v2.benchmarks.qs.a7a.cell import configuration, settings
        from v2.benchmarks.qs.a7a.combined import cell, fixtures

        panel = list(fixtures().values())
        self.assertEqual([r["digits"] for r in panel].count(60), 3)
        self.assertEqual([r["digits"] for r in panel].count(70), 2)
        for index, fixture in enumerate(panel):
            n = fixture["n"]
            self.assertEqual(prod(fixture["factors"]), n)
            self.assertFalse(fixture["uneven"])
            self.assertTrue(all(gcd(n, r["n"]) == 1 for r in panel[:index]))
            for arm in ("S0", "P4", "C4", "C6"):
                spec = cell(arm)
                config = replace(
                    configuration(n, "Y0"),
                    memory_bytes=spec["qs_memory_bytes"],
                    checkpoint_bytes=spec["checkpoint_bytes"],
                )
                self.assertEqual(config.backend, "python-int")
                workers = settings(
                    arm, spec["memory"], False, 64, 2048, spec["features"]
                )
                if workers is not None:
                    self.assertLessEqual(workers.workers, 10)
                    self.assertFalse(workers.arrival_order)
                if arm.startswith("P"):
                    self.assertFalse(workers.lease_budget)

    def test_combined_stage_zero_brackets_complete_local_packet_cost(self):
        plan = schedule("J0")
        self.assertEqual(plan["wall"], 900)
        self.assertEqual(
            [r["arm"] for r in plan["cells"]], ["M0", "S0", "S0", "M0"]
        )
        self.assertEqual(len(plan["cells"]), 4)
        for row in plan["cells"]:
            self.assertEqual(row["endpoint"], "complete")
            self.assertEqual(row["lifetime"], "cold")
            self.assertNotIn("cpu_guard", row)
            self.assertNotIn("cpu_profile", row)
            self.assertEqual(row["polynomials_per_task"], 64)
            self.assertEqual(row["qs_memory_bytes"], 4 * 2**30)

    def test_wide_guard_probe_brackets_full_calls_without_repricing(self):
        plan = schedule("H14")
        self.assertEqual(
            [r["cpu_guard"] for r in plan["cells"]],
            ["peer", "parent", "parent", "peer"],
        )
        self.assertEqual(
            (plan["reserve_wall"], plan["reserve_cpu"]), (300, 3000)
        )
        reference = schedule("H13")["cells"][3]
        for row in plan["cells"]:
            self.assertEqual(
                {
                    k: v
                    for k, v in row.items()
                    if k not in ("stage", "id", "cpu_guard")
                },
                {
                    k: v
                    for k, v in reference.items()
                    if k not in ("stage", "id")
                },
            )
            self.assertNotIn("cpu_profile", row)

    def test_wider_counts_share_a_funded_complete_cost_configuration(self):
        from v2.benchmarks.infrastructure.resource_control import (
            historical_replace as replace,
        )
        from v2.benchmarks.qs.a7a.cell import configuration, settings

        plan = schedule("H13")
        self.assertEqual(
            [r["arm"] for r in plan["cells"]],
            ["S0", "P4", "P6", "P10", "P8", "P4"],
        )
        self.assertEqual(
            (plan["reserve_wall"], plan["reserve_cpu"]), (300, 3000)
        )
        fixtures = {r["id"]: r for r in corpus()["fixtures"]}
        for size in (60, 70, 80):
            n = fixtures[f"exposed_{size}_0"]["n"]
            for row in plan["cells"]:
                config = replace(
                    configuration(n, row["policy"]),
                    checkpoint_bytes=row["checkpoint_bytes"],
                    memory_bytes=row["qs_memory_bytes"],
                )
                self.assertEqual(config.backend, "python-int")
                workers = settings(
                    row["arm"],
                    row["memory"],
                    row["observe"],
                    row["polynomials_per_task"],
                    row["max_batch_atoms"],
                )
                if workers is not None:
                    self.assertLessEqual(workers.workers, 10)
                self.assertEqual(row["endpoint"], "complete")
                self.assertEqual(row["samples"], 1)
                self.assertNotIn("cpu_profile", row)

    def test_counter_storage_change_compares_both_complete_serial_and_workers(
        self,
    ):
        plan = schedule("H12")
        self.assertEqual(
            (plan["reserve_wall"], plan["reserve_cpu"]), (300, 3000)
        )
        self.assertEqual(
            [(r["arm"], r["implementation"]) for r in plan["cells"]],
            [
                ("S0", "budget-control"),
                ("S0", "candidate"),
                ("P4", "candidate"),
                ("P4", "budget-control"),
            ],
        )
        for row in plan["cells"]:
            self.assertEqual(row["endpoint"], "complete")
            self.assertNotIn("cpu_profile", row)
            self.assertNotIn("cpu_guard", row)
            self.assertEqual(row["memory"], 20 * 2**30)

    def test_parent_guard_full_calls_bracket_clean_cost_before_late_profile(
        self,
    ):
        plan = schedule("H11")
        self.assertEqual(len(plan["cells"]), 5)
        self.assertEqual(
            (plan["reserve_wall"], plan["reserve_cpu"]), (300, 3000)
        )
        self.assertEqual(
            [r["cpu_guard"] for r in plan["cells"]],
            ["peer", "parent", "parent", "peer", "peer"],
        )
        for row in plan["cells"]:
            self.assertEqual(row["endpoint"], "complete")
            self.assertEqual(row["qs_memory_bytes"], 12 * 2**30)
            self.assertEqual(row["memory"], 20 * 2**30)
        self.assertFalse(any("cpu_profile" in r for r in plan["cells"][:4]))
        self.assertEqual(plan["cells"][-1]["cpu_profile"], [60, 75])

    def test_parent_guard_probe_brackets_both_sizes_with_production_controls(
        self,
    ):
        plan = schedule("H10")
        self.assertEqual(len(plan["cells"]), 8)
        self.assertEqual(
            (plan["reserve_wall"], plan["reserve_cpu"]), (120, 1200)
        )
        for fixture in ("exposed_60_0", "exposed_70_0"):
            rows = [r for r in plan["cells"] if r["fixture"] == fixture]
            self.assertEqual(
                [r["cpu_guard"] for r in rows],
                ["peer", "parent", "parent", "peer"],
            )
            self.assertEqual({r["polynomials"] for r in rows}, {768})
            self.assertEqual({r["cpu"] for r in rows}, {450})
            self.assertEqual({r["endpoint"] for r in rows}, {"prefix"})

    def test_cpu_read_isolation_preserves_controls_and_fixed_work(self):
        plan = schedule("H9")
        self.assertEqual(len(plan["cells"]), 6)
        self.assertEqual(
            (plan["reserve_wall"], plan["reserve_cpu"]), (120, 1200)
        )
        for fixture in ("exposed_60_0", "exposed_70_0"):
            rows = [r for r in plan["cells"] if r["fixture"] == fixture]
            self.assertEqual(
                {r["implementation"] for r in rows},
                {"candidate", "joint-control", "cpu-read-control"},
            )
            for row in rows:
                self.assertEqual(row["samples"], 9)
                self.assertEqual(row["lifetime"], "warm")
                self.assertEqual(row["polynomials"], 256)
                self.assertEqual(row["polynomials_per_task"], 64)

    def test_cpu_read_ablation_compares_full_cost_with_strong_serial(self):
        plan = schedule("H8")
        self.assertEqual([r["arm"] for r in plan["cells"]], ["S0", "P4", "P4"])
        self.assertEqual(
            [r["implementation"] for r in plan["cells"]],
            ["candidate", "cpu-read-control", "candidate"],
        )
        self.assertEqual(
            (plan["reserve_wall"], plan["reserve_cpu"]), (150, 1500)
        )
        for row in plan["cells"]:
            self.assertEqual(row["endpoint"], "complete")
            self.assertEqual(row["memory"], 20 * 2**30)
            self.assertEqual(row["qs_memory_bytes"], 12 * 2**30)

    def test_statistics_ablation_funds_control_and_stronger_serial(self):
        plan = schedule("H7")
        self.assertEqual(len(plan["cells"]), 3)
        self.assertEqual((plan["wall"], plan["cpu"]), (1050, 10500))
        self.assertEqual(
            (plan["reserve_wall"], plan["reserve_cpu"]), (150, 1500)
        )
        self.assertEqual([r["arm"] for r in plan["cells"]], ["S0", "P4", "S0"])
        self.assertEqual(
            [r["implementation"] for r in plan["cells"]],
            ["statistics-control", "candidate", "candidate"],
        )
        comparable = []
        for row in plan["cells"]:
            self.assertEqual(row["qs_memory_bytes"], 12 * 2**30)
            comparable.append(
                {
                    k: v
                    for k, v in row.items()
                    if k not in ("id", "arm", "implementation")
                }
            )
        self.assertEqual(comparable[0], comparable[1])
        self.assertEqual(comparable[0], comparable[2])

    def test_revised_pilot_prices_snapshot_reserve_without_retuning(self):
        from v2.benchmarks.infrastructure.resource_control import (
            historical_replace as replace,
        )
        from v2.benchmarks.qs.a7a.cell import configuration, settings

        plan = schedule("H6")
        self.assertEqual(len(plan["cells"]), 2)
        self.assertEqual(plan["machine_memory"], 20 * 2**30)
        fixture = next(
            r for r in corpus()["fixtures"] if r["id"] == "exposed_60_0"
        )
        original = configuration(fixture["n"], "Y0")
        for row in plan["cells"]:
            self.assertEqual(row["qs_memory_bytes"], 12 * 2**30)
            self.assertEqual(row["memory"], plan["machine_memory"])
            revised = replace(
                original,
                checkpoint_bytes=row["checkpoint_bytes"],
                memory_bytes=row["qs_memory_bytes"],
            )
            self.assertFalse(hasattr(revised, "memory_limit"))
            self.assertEqual(revised.collector, original.collector)
            self.assertEqual(revised.base_count, original.base_count)
            workers = settings(row["arm"], row["memory"], False, 64)
            if workers is not None:
                self.assertLessEqual(workers.workers, 10)

    def test_complete_pilot_funds_both_arms_without_prefix_quota(self):
        plan = schedule("H5")
        self.assertEqual(
            (plan["reserved_wall"], plan["reserved_cpu"]), (600, 6000)
        )
        self.assertEqual(
            (plan["reserve_wall"], plan["reserve_cpu"]), (300, 3000)
        )
        self.assertEqual({r["arm"] for r in plan["cells"]}, {"S0", "P4"})
        self.assertEqual(
            {r["fixture"] for r in plan["cells"]}, {"exposed_60_0"}
        )
        for row in plan["cells"]:
            self.assertEqual(row["endpoint"], "complete")
            self.assertEqual(row["lifetime"], "cold")
            self.assertEqual(row["polynomials_per_task"], 64)
            self.assertEqual(row["checkpoint_bytes"], 64 * 2**20)
            self.assertNotIn("polynomials", row)

    def test_coarse_search_funds_wide_packets_and_identical_coverage(self):
        plan = schedule("H4")
        self.assertEqual(len(plan["cells"]), 20)
        self.assertEqual(
            (plan["reserved_wall"], plan["reserved_cpu"]), (600, 6000)
        )
        self.assertEqual(
            (plan["reserve_wall"], plan["reserve_cpu"]), (150, 1500)
        )
        for fixture in ("exposed_60_0", "exposed_70_0"):
            rows = [r for r in plan["cells"] if r["fixture"] == fixture]
            self.assertEqual({r["polynomials"] for r in rows}, {256})
            for arm in ("P1", "P2", "P4"):
                widths = {
                    r["polynomials_per_task"] for r in rows if r["arm"] == arm
                }
                self.assertEqual(widths, {4, 16, 64})
                workers = int(arm[1:])
                self.assertGreaterEqual(256 // max(widths), workers)

    def test_polling_comparison_keeps_immutable_control_and_matched_work(self):
        plan = schedule("H3")
        self.assertEqual(len(plan["cells"]), 10)
        self.assertEqual(
            (plan["reserved_wall"], plan["reserved_cpu"]), (450, 4500)
        )
        for fixture in ("exposed_60_0", "exposed_70_0"):
            rows = [r for r in plan["cells"] if r["fixture"] == fixture]
            for arm in ("P1", "P4"):
                pair = [r for r in rows if r["arm"] == arm]
                self.assertEqual(
                    {r["implementation"] for r in pair},
                    {"candidate", "packet-control"},
                )
                self.assertEqual({r["polynomials"] for r in pair}, {128})
                self.assertEqual(
                    {r["polynomials_per_task"] for r in pair}, {4}
                )

    def test_packet_screen_pins_equal_work_instead_of_equal_task_counts(self):
        plan = schedule("H2")
        self.assertEqual(len(plan["cells"]), 14)
        self.assertEqual((plan["wall"], plan["cpu"]), (600, 6000))
        self.assertEqual(plan["reserved_wall"], 420)
        for fixture in ("exposed_60_0", "exposed_70_0"):
            rows = [r for r in plan["cells"] if r["fixture"] == fixture]
            self.assertEqual({r["polynomials"] for r in rows}, {128})
            for arm in ("M0", "P1"):
                self.assertEqual(
                    {
                        r["polynomials_per_task"]
                        for r in rows
                        if r["arm"] == arm
                    },
                    {4, 16, 64},
                )

    def test_warm_diagnostic_matches_work_and_bounds_observation(self):
        plan = schedule("H1")
        self.assertEqual(len(plan["cells"]), 18)
        self.assertEqual((plan["wall"], plan["cpu"]), (600, 6000))
        self.assertEqual(
            (plan["reserved_wall"], plan["reserved_cpu"]), (540, 5400)
        )
        self.assertEqual(plan["storage"], 32 * 2**20)
        for fixture in ("exposed_60_0", "exposed_70_0"):
            rows = [r for r in plan["cells"] if r["fixture"] == fixture]
            self.assertEqual({r["tasks"] for r in rows}, {16})
            self.assertEqual({r["seed"] for r in rows}, {104729})
            warm = [r for r in rows if r["lifetime"] == "warm"]
            self.assertEqual(
                {r["arm"] for r in warm}, {"S0", "M0", "P1", "P2", "P4"}
            )
            observed = [r for r in rows if r["observe"]]
            self.assertEqual(len(observed), 1)
            self.assertEqual(observed[0]["arm"], "P4")
            self.assertEqual(observed[0]["lifetime"], "warm")
        self.assertTrue(all(r["endpoint"] == "prefix" for r in plan["cells"]))

    def test_current_mainline_screen_is_small_clean_and_source_independent(
        self,
    ):
        plan = schedule("H0")
        self.assertEqual(len(plan["cells"]), 14)
        self.assertEqual((plan["wall"], plan["cpu"]), (600, 6000))
        self.assertEqual(plan["reserved_wall"], 360)
        self.assertEqual(plan["reserved_cpu"], 3600)
        fixtures = {row["id"] for row in corpus()["fixtures"]}
        for row in plan["cells"]:
            self.assertIn(row["fixture"], fixtures)
            self.assertFalse(row["observe"])
            self.assertEqual(row["lifetime"], "cold")
            self.assertEqual(row["cpu"], 10 * row["wall"])
            if row["fixture"] == "exposed_70_0":
                self.assertEqual(row["endpoint"], "prefix")

    def test_stage_a_has_all_named_cells_and_inclusive_budget(self):
        plan = schedule("A")
        self.assertEqual(len(plan["cells"]), 60)
        self.assertEqual(
            (plan["reserved_wall"], plan["reserved_cpu"]), (2430, 9720)
        )
        self.assertEqual(
            (plan["reserve_wall"], plan["reserve_cpu"]), (1170, 4680)
        )
        for row in plan["cells"]:
            self.assertFalse(
                row["fixture"] == "exposed_70_0"
                and row["endpoint"] == "complete"
            )

    def test_all_stage_sums_and_worker_scaling_observer_coverage(self):
        for stage in "ABCDE":
            plan = schedule(stage)
            self.assertLessEqual(plan["reserved_wall"], plan["wall"])
            self.assertLessEqual(plan["reserved_cpu"], plan["cpu"])
            self.assertEqual(
                len({r["id"] for r in plan["cells"]}), len(plan["cells"])
            )
        observed = defaultdict(set)
        for row in schedule("B")["cells"]:
            if row["observe"]:
                observed[row["fixture"]].add(row["arm"])
        self.assertTrue(observed)
        self.assertTrue(
            all(arms == {"P1", "P2", "P4"} for arms in observed.values())
        )

    def test_every_preconfirmation_cell_has_a_committed_input(self):
        fixtures = {r["id"]: r for r in corpus()["fixtures"]}
        for stage in "ABCD":
            for row in schedule(stage)["cells"]:
                self.assertIn(row["fixture"], fixtures)
        random = [
            r
            for r in fixtures.values()
            if r.get("population") == "uniform_integer"
        ]
        self.assertEqual(len(random), 36)
        self.assertEqual({r["digits"] for r in random}, {60, 80, 100})
        self.assertTrue(all(len(str(r["n"])) == r["digits"] for r in random))

    def test_quantiles_include_censored_calls_in_denominator(self):
        rows = [
            dict(seconds=2, complete=True, cap=10),
            dict(seconds=10, complete=False, cap=10),
        ]
        self.assertEqual(quantile_bound(rows, 0.5)["value"], 2)
        self.assertFalse(quantile_bound(rows, 0.9)["identified"])
        self.assertEqual(quantile_bound(rows, 0.9)["lower_bound"], 10)

    def test_repeated_seeds_do_not_inflate_independent_moduli(self):
        rows = []
        for fixture in ("a", "b", "c"):
            for seed in (1, 2, 3):
                for arm, seconds in (("S0", 8), ("P2", 6)):
                    rows.append(
                        dict(
                            fixture=fixture,
                            seed=seed,
                            arm=arm,
                            seconds=seconds,
                            complete=True,
                            cap=10,
                        )
                    )
        result = paired_interval(rows, nominee="P2", resamples=100)
        self.assertEqual(result["independent_moduli"], 3)
        self.assertEqual(result["saving_interval95"], [2, 2])
        self.assertGreater(result["zero_loss_upper95"], 0.5)

    def test_one_shot_result_formatter_preserves_certainty_and_reconstruction(
        self,
    ):
        from v2.benchmarks.infrastructure.resource_control import Budget
        from v2.benchmarks.qs.a7a.cell import complete, configuration

        budget = Budget(work_limit=1000000, seconds=2, cpu_seconds=2)
        budget.a7a_seed = 7
        result = complete(
            21,
            configuration(21, "G0"),
            None,
            budget,
            None,
            "portfolio",
            memory_bytes=20 * 2**30,
        )
        self.assertTrue(result["complete"])
        self.assertEqual(result["remaining"], [])
        self.assertEqual([factor[0] for factor in result["factors"]], [3, 7])
        self.assertIsNone(result["active_stage"])
        self.assertTrue(
            all(
                value >= 0
                for value in result["portfolio_stage_seconds"].values()
            )
        )

    def test_committed_corpus_sources_products_and_certificates(self):
        from v2.benchmarks.qs.a7a.inputs import validate_committed

        validate_committed()
