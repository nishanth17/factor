"""A coherent private CPU view without arithmetic under the shared lock."""

import unittest
from types import SimpleNamespace

from v2.execution.budget import Budget
from v2.qs.packet_lease import ParentCPUSnapshot


class ParentSnapshotTests(unittest.TestCase):
    def test_copy_is_locked_reduction_is_private_and_cpu_is_charged_once(self):
        owner = self

        class Lock:
            active = False
            acquisitions = 0

            def __enter__(self):
                owner.assertFalse(self.active)
                self.active = True
                self.acquisitions += 1

            def __exit__(self, *args):
                self.active = False
                if self.acquisitions == 1:
                    # A publication after release must affect the next view,
                    # not the private values currently being reduced.
                    values[:] = [Sample(10), Sample(20)]

        lock = Lock()

        class Sample(float):
            def __sub__(self, prior):
                owner.assertFalse(lock.active)
                return float(self) - prior

        class Values(list):
            def __getitem__(self, index):
                owner.assertTrue(lock.active)
                return super().__getitem__(index)

        values = Values([Sample(1), Sample(2)])
        counters = SimpleNamespace(
            get_obj=lambda: values, get_lock=lambda: lock
        )
        pool = SimpleNamespace(
            state=dict(
                cpu_base=[0.0, 0.0],
                cpu=counters,
                parent_cpu=SimpleNamespace(value=0.0),
            ),
            accounted_cpu=0.0,
        )
        snapshot = ParentCPUSnapshot(pool)
        budget = Budget(work_limit=100, seconds=None, cpu_seconds=None)

        snapshot(budget)
        self.assertEqual((budget.prior_cpu, pool.accounted_cpu), (3, 3))
        snapshot(budget)
        self.assertEqual((budget.prior_cpu, pool.accounted_cpu), (30, 30))
        snapshot(budget)
        self.assertEqual((budget.prior_cpu, pool.accounted_cpu), (30, 30))
        self.assertEqual(lock.acquisitions, 3)


if __name__ == "__main__":
    unittest.main()
