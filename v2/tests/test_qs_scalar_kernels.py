"""Default specializations retain generic arithmetic and paid work."""

import ast
import inspect
import unittest
from collections import Counter
from types import MethodType

from v2.qs import scalar_kernels, sparse_queue
from v2.qs.sieve_collector import SieveCollector
from v2.tests.test_qs_hit_words import make_worker
from v2.tests.test_qs_sieve import signature


class ScalarKernelTests(unittest.TestCase):
    def test_default_loops_have_no_disabled_mechanism_guards(self):
        for function in (
            scalar_kernels.direct_scalar,
            scalar_kernels.collect_scalar,
            scalar_kernels.build_native,
            scalar_kernels.sieve_native,
        ):
            tree = ast.parse(inspect.getsource(function))
            forbidden = {
                "grouped_budget",
                "power_period",
                "candidates",
                "packed",
            }
            self.assertFalse(
                {n.id for n in ast.walk(tree) if isinstance(n, ast.Name)}
                & forbidden
            )
            self.assertFalse(
                any(
                    isinstance(n, ast.Constant)
                    and n.value in ("scan_fused", "packed_queue")
                    for n in ast.walk(tree)
                )
            )

    def test_buffers_charges_and_complete_relations_match_generic(self):
        for sparse in (0, 2**20):
            kernels = ("presieve", "unroll4", "parity_two")
            fast = make_worker(64, kernels, sparse)
            generic = make_worker(64, kernels, sparse)
            generic._sieve_direct = MethodType(
                SieveCollector._sieve_direct, generic
            )
            generic._collect = MethodType(SieveCollector._collect, generic)
            generic._sparse_sieve = sparse_queue.sieve
            for lo, hi in ((-79, -15), (3, 16), (-31, 33)):
                a, b = Counter(), Counter()
                self.assertEqual(
                    fast._sieve(lo, hi, a), generic._sieve(lo, hi, b)
                )
                self.assertEqual(fast._scores, generic._scores)
                self.assertEqual(fast._hits[:], generic._hits[:])
                self.assertEqual(a, b)
                self.assertEqual(fast.budget.used, generic.budget.used)
            a, b = fast.collect(-79, 82), generic.collect(-79, 82)
            self.assertEqual(signature(a), signature(b))
            self.assertEqual(a.next_position, b.next_position)
            self.assertEqual(fast.budget.used, generic.budget.used)
