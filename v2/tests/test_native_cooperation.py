"""Check bounded slices and canonical outcomes with optional unit yields."""

import copy
import unittest

from v2.ecm.native_packet import SOURCE
from v2.execution import stage_jobs
from v2.execution.portfolio_config import (
    _execution_context,
    close_execution_context,
)
from v2.tests.test_native_correctness import deterministic_budget, selected
from v2.tests.test_native_full import drive


@unittest.skipUnless(
    (SOURCE / "build/manifest.json").exists(), "scalar native build required"
)
class NativeCooperationTests(unittest.TestCase):
    def test_build_quantum_and_smaller_override_preserve_resume(self):
        case = dict(n=2**127 - 1, seed=7, b1=11000, b2=1100000)
        expected, _, _ = drive(case, "python")
        for explicit in (None, 1):
            config = selected("c-gmp-full", case["b1"], case["b2"])
            budget = deterministic_budget()
            context, programs = _execution_context({}, config, budget)
            source = programs or context
            packet = source.native_packet
            packet.actions = explicit
            job = stage_jobs.new_job(
                "ecm", case["n"], case["seed"], case["b1"], case["b2"]
            )
            previous = 0
            try:
                for _ in range(10000):
                    stage_jobs.advance_job(job, budget, source, config)
                    info = packet._read()
                    actions = int(
                        info.setup_ops
                        + info.stage_one_ops
                        + info.stage_two_ops
                    )
                    limit = (
                        packet.manifest["slice_actions"]
                        if explicit is None
                        else explicit
                    )
                    self.assertIn(packet.lib.b4_slice_actions(), (128, 256))
                    self.assertEqual(
                        packet.lib.b4_slice_actions(),
                        packet.manifest["slice_actions"],
                    )
                    self.assertLessEqual(actions - previous, limit)
                    self.assertGreaterEqual(actions - previous, 0)
                    if actions - previous == limit:
                        break
                    previous = actions
                else:
                    self.fail("selected finite action quantum was not reached")
            finally:
                close_execution_context(source)
            resumed, _, _ = drive(case, "python", job=job)
            self.assertEqual(resumed, expected)

    def test_bounded_actions_and_canonical_result_with_nonatomic_yields(self):
        case = dict(n=2**127 - 1, seed=7, b1=31, b2=997)
        expected, _, _ = drive(case, "python")
        config = selected("c-gmp-full", case["b1"], case["b2"])
        budget = deterministic_budget()
        context, programs = _execution_context({}, config, budget)
        source = programs or context
        packet = source.native_packet
        job = stage_jobs.new_job(
            "ecm", case["n"], case["seed"], case["b1"], case["b2"]
        )
        observed, previous = [], 0
        try:
            for _ in range(10000):
                stage_jobs.advance_job(job, budget, source, config)
                info = packet._read()
                actions = sum(
                    int(getattr(info, name))
                    for name in (
                        "setup_ops",
                        "stage_one_ops",
                        "stage_two_ops",
                    )
                )
                observed.append(actions - previous)
                previous = actions
                limit = (
                    packet.manifest.get("slice_actions", 128)
                    if packet.actions is None
                    else packet.actions
                )
                self.assertLessEqual(observed[-1], limit)
                self.assertGreaterEqual(observed[-1], 0)
                if job["done"]:
                    packet.sync(allow_result=True)
                    break
            else:
                self.fail("bounded complete native curve did not finish")

            actual = copy.deepcopy(job)
            self.assertEqual(actual, expected)
            self.assertTrue(any(observed))
            self.assertIn(packet.lib.b4_coalesced_unit_yields(), (0, 1))
            self.assertFalse(packet.atomic_actions)
        finally:
            close_execution_context(source)
