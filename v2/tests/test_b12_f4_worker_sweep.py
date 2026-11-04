"""Worker selection includes serial paths and distinguishes cold endpoints."""

import unittest

from v2.benchmarks.ecm.b12_f4.worker_sweep import select


class WorkerSweepTests(unittest.TestCase):
    def test_cold_cli_winner_can_differ_from_warm_component_winner(self):
        rows = []
        for mode in ("cold", "warm"):
            for arm, workers, cli, call in (
                ("M0", 1, 0.4, 0.2),
                ("P4", 4, 0.7, 0.1),
                ("P16", 16, 1.4, 0.05),
            ):
                rows.append(
                    dict(
                        input="C60",
                        mode=mode,
                        arm=arm,
                        workers=workers,
                        cli_wall=cli,
                        median=call,
                        observer="none",
                        complete=True,
                    )
                )
        result = select(rows)
        self.assertEqual([row["arm"] for row in result], ["M0", "P16"])
        self.assertEqual(
            [row["metric"] for row in result], ["cli_wall", "median"]
        )

    def test_instrumented_or_incomplete_arm_cannot_win_selection(self):
        base = dict(
            input="H100",
            mode="cold",
            workers=4,
            cli_wall=1,
            median=1,
            observer="none",
            complete=True,
        )
        rows = [
            dict(base, arm="P4"),
            dict(base, arm="P16", cli_wall=0.1, observer="cpu"),
            dict(base, arm="P10", cli_wall=0.01, complete=False),
        ]
        self.assertEqual(select(rows)[0]["arm"], "P4")
