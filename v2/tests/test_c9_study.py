"""Independent incidence and domain checks for the C9 study harness."""

import unittest
from dataclasses import asdict
from types import SimpleNamespace

from v2.benchmarks.qs.b1.b1_calibration import common_config
from v2.benchmarks.qs.c9.c9_study import configuration, graph_diagnostics


class C9StudyTests(unittest.TestCase):
    def test_alignment_preserves_dlp_candidate_and_product_domains(self):
        base = common_config(30000, 65536, 6, rows=8192, partials=8192)
        control = asdict(configuration(base, "dlp"))
        aligned = asdict(configuration(base, "aligned"))
        self.assertEqual(control["collector"]["residual_bound"], 900000000)
        self.assertEqual(aligned["collector"]["residual_bound"], 3000000)
        aligned["collector"]["residual_bound"] = 900000000
        self.assertEqual(control, aligned)

    def test_multigraph_core_and_disconnected_cycles(self):
        # A triangle, a parallel pair, a loop and a two-edge hanging path.
        pairs = [
            (3, 5),
            (5, 7),
            (7, 3),
            (11, 13),
            (11, 13),
            (17, 17),
            (7, 19),
            (19, 23),
        ]
        collector = SimpleNamespace(
            _atoms={
                i: SimpleNamespace(large_primes=pair, residual=1)
                for i, pair in enumerate(pairs)
            },
            config=SimpleNamespace(large_prime_bound=100),
        )
        result = graph_diagnostics(collector)
        # Enumerate every even-incidence subset independently of peeling.
        kernel = []
        for mask in range(1 << len(pairs)):
            odd = set()
            for i, pair in enumerate(pairs):
                if mask & (1 << i):
                    for vertex in pair:
                        odd.symmetric_difference_update((vertex,))
            if not odd:
                kernel.append(mask)
        self.assertEqual(len(kernel), 1 << result["retained_cycle_rank"])
        self.assertEqual(result["retained_cycle_rank"], 3)
        self.assertEqual(result["core_edges"], 6)
        self.assertEqual(result["components"], 3)
        self.assertEqual(result["core_vertices"], 6)


if __name__ == "__main__":
    unittest.main()
