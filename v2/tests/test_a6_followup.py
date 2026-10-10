"""Independent p-1 paired coverage, arithmetic and bounded-resume checks."""

import copy
import importlib.util
import io
import json
import math
import sys
import unittest
from contextlib import redirect_stdout
from dataclasses import replace
from pathlib import Path

from v2 import pm1_bounded, utils
from v2.pm1_tuning import PM1TuningConfig
from v2.schedules import SieveContext
from v2.tests.test_a6_pm1 import allowance, resign


def configured(**options):
    return PM1TuningConfig(bounds=((10, 101),), segment_size=3, **options)


def prime_oracle(lo, hi):
    return [
        n
        for n in range(lo, hi + 1)
        if all(n % d for d in range(2, math.isqrt(n) + 1))
    ]


def frozen_control():
    name = "v2._a6_followup_frozen"
    if name not in sys.modules:
        source = (
            Path(__file__).parents[1]
            / "benchmarks/inputs/baselines/a6_followup_pm1_32b3c65.py"
        )
        spec = importlib.util.spec_from_file_location(name, source)
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
    return sys.modules[name]


class PairedCoverageTests(unittest.TestCase):
    def test_every_eligible_prime_and_unit_scaled_batch_product(self):
        n = 1000003
        for wheel in (30, 210):
            for b1, b2 in (
                (2, 7),
                (7, 29),
                (10, 31),
                (14, 97),
                (29, 210),
                (30, 211),
                (31, 239),
            ):
                for batch in (1, 3, 64):
                    cfg = configured(wheel=wheel, gcd_batch=batch)
                    cfg = replace(cfg, bounds=((b1, b2),))
                    state = pm1_bounded._initial(n, 2, cfg)
                    context = SieveContext(b2 + 1, segment_size=3)
                    ledger, covered = allowance(), []

                    while not state["done"]:
                        phase, before = state["phase"], len(state["terms"])
                        pm1_bounded._advance(state, ledger, context, cfg)
                        if phase != "wheel_terms":
                            continue
                        value = pow(2, math.lcm(*range(1, b1 + 1)), n)
                        self.assertEqual(state["value"], value)
                        for primes, term in zip(
                            state["wheel_term_primes"][before:],
                            state["terms"][before:],
                        ):
                            direct = (
                                math.prod(pow(value, q, n) - 1 for q in primes)
                                % n
                            )
                            if len(primes) == 2:
                                center = sum(primes) // 2
                                direct = direct * pow(value, -center, n) % n
                            self.assertEqual(term, direct)
                            covered.extend(primes)
                        self.assertEqual(
                            state["product"], math.prod(state["terms"]) % n
                        )

                    self.assertEqual(sorted(covered), prime_oracle(b1 + 1, b2))
                    self.assertEqual(len(covered), len(set(covered)))
                    self.assertEqual(state["reason"], "exhausted")

    def test_trace_identity_over_composite_rings_and_prime_powers(self):
        for n in (9, 15, 35, 49, 77):
            for value in range(2, min(n, 12)):
                if math.gcd(value, n) != 1:
                    continue
                for center, distance in ((15, 2), (30, 1), (210, 19)):
                    giant = pow(value, center, n)
                    baby = pow(value, distance, n)
                    trace = (
                        giant + pow(giant, -1, n) - baby - pow(baby, -1, n)
                    ) % n
                    direct = (
                        (pow(value, center - distance, n) - 1)
                        * (pow(value, center + distance, n) - 1)
                        % n
                    )
                    self.assertEqual(
                        trace, direct * pow(value, -center, n) % n
                    )
                    self.assertEqual(math.gcd(trace, n), math.gcd(direct, n))

    def test_independent_coverage_after_increased_b1_and_b2_only(self):
        bounds = ((7, 31), (10, 43), (10, 61))
        for wheel in (30, 210):
            cfg = replace(configured(wheel=wheel, gcd_batch=2), bounds=bounds)
            state = pm1_bounded._initial(1000003, 2, cfg)
            context, ledger = SieveContext(62, segment_size=3), allowance()
            covered = {0: [], 1: [], 2: []}
            while not state["done"]:
                phase, rung = state["phase"], state["rung"]
                before = len(state["terms"])
                pm1_bounded._advance(state, ledger, context, cfg)
                if phase == "wheel_terms":
                    value = pow(
                        2, math.lcm(*range(1, bounds[rung][0] + 1)), 1000003
                    )
                    self.assertEqual(state["value"], value)
                    for primes in state["wheel_term_primes"][before:]:
                        covered[rung].extend(primes)
            for rung, (lower, upper) in enumerate(
                ((8, 31), (11, 43), (44, 61))
            ):
                self.assertEqual(
                    sorted(covered[rung]), prime_oracle(lower, upper)
                )

    def test_singleton_excludes_partner_just_outside_bound(self):
        for wheel in (30, 210):
            cfg = configured(wheel=wheel)
            outside = pm1_bounded.factorize_pm1_bounded(
                311 * 1019,
                config=replace(cfg, bounds=((10, 30),)),
                budget=allowance(),
            )
            inside = pm1_bounded.factorize_pm1_bounded(
                311 * 1019,
                config=replace(cfg, bounds=((10, 31),)),
                budget=allowance(),
            )
            self.assertIsNone(outside.divisor)
            self.assertEqual(inside.divisor, 311)

    def test_same_pair_mixed_factors_replay_individual_primes(self):
        # Orders after M(10) are 29 and 31, paired around center 30.
        n = 59 * 311
        cfg = configured(wheel=30, gcd_batch=1)
        cfg = replace(cfg, bounds=((10, 31),))
        outcome = pm1_bounded.factorize_pm1_bounded(
            n, config=cfg, budget=allowance()
        )
        self.assertEqual(outcome.divisor, 59)
        self.assertEqual(
            outcome.checkpoint["payload"]["state"]["phase"], "term_replay"
        )
        self.assertEqual(outcome.checkpoint["payload"]["state"]["recovery"], 1)
        refused = pm1_bounded.factorize_pm1_bounded(
            n, config=replace(cfg, recovery_limit=0), budget=allowance()
        )
        self.assertEqual(refused.reason, "saturated")
        self.assertEqual(refused.result.remaining, (n,))

    def test_cross_pair_saturation_and_nonunits(self):
        for wheel in (30, 210):
            cfg = configured(wheel=wheel, gcd_batch=64)
            outcome = pm1_bounded.factorize_pm1_bounded(
                607 * 103, config=cfg, budget=allowance()
            )
            self.assertTrue(utils.valid_divisor(outcome.divisor, 607 * 103))
            nonunit = pm1_bounded.factorize_pm1_bounded(
                35, base=5, config=cfg, budget=allowance()
            )
            self.assertEqual(nonunit.divisor, 5)
            rejected = pm1_bounded.factorize_pm1_bounded(
                35, base=35, config=cfg, budget=allowance()
            )
            self.assertEqual(rejected.reason, "nonunit")


class ExecutionTests(unittest.TestCase):
    def test_bit_chunks_and_even_powers_against_small_direct_oracles(self):
        for cfg in (
            configured(chunk_size=256, chunk_bits=32),
            configured(chunk_size=64, gap_mode="recurrence", gap_entries=1),
            configured(chunk_size=256, chunk_bits=64, gap_mode="recurrence"),
        ):
            state = pm1_bounded._initial(1000003, 2, cfg)
            context = SieveContext(102, segment_size=3)
            while not state["done"]:
                pm1_bounded._advance(state, allowance(), context, cfg)
                if cfg.chunk_bits:
                    self.assertLessEqual(
                        sum(p.bit_length() for _, p in state["powers"]),
                        cfg.chunk_bits,
                    )
                if state["phase"] == "stage_two":
                    value = pow(2, math.lcm(*range(1, 11)), 1000003)
                    self.assertEqual(state["value"], value)
                    self.assertEqual(
                        state["stage_two_value"],
                        pow(value, state["previous_prime"], 1000003),
                    )
                    for i, power in enumerate(state.get("even_powers", []), 1):
                        self.assertEqual(power, pow(value, 2 * i, 1000003))

    def test_legacy_actions_state_work_and_checkpoint_remain_compatible(self):
        frozen = frozen_control()
        cfg = pm1_bounded.PM1Config(bounds=((7, 13), (10, 31)), segment_size=3)
        old_cfg = frozen.PM1Config(bounds=cfg.bounds, segment_size=3)
        state = pm1_bounded._initial(1000003, 2, cfg)
        old = frozen._initial(1000003, 2, old_cfg)
        ledger, old_ledger = allowance(), allowance()
        context = SieveContext(32, segment_size=3)
        while not old["done"]:
            frozen._advance(old, old_ledger, context, old_cfg)
            pm1_bounded._advance(state, ledger, context, cfg)
            self.assertEqual(state, old)
            self.assertEqual(ledger.used, old_ledger.used)
        checkpoint = frozen.factorize_pm1_bounded(
            1000003, config=old_cfg, budget=allowance(), max_actions=12
        ).checkpoint
        resumed = pm1_bounded.factorize_pm1_bounded(
            1000003, config=cfg, budget=allowance(), checkpoint=checkpoint
        )
        self.assertEqual(resumed.reason, "exhausted")
        self.assertEqual(
            resumed.checkpoint["payload"]["execution"], "pm1-campaign-v1"
        )

    def test_every_boundary_resume_atomic_refusal_and_cancellation(self):
        for cfg in (
            configured(chunk_bits=32, chunk_size=256, gap_mode="recurrence"),
            configured(wheel=30, gcd_batch=2),
            configured(wheel=210, gcd_batch=2),
        ):
            cfg = replace(cfg, bounds=((7, 43), (10, 101), (10, 127)))
            complete = pm1_bounded.factorize_pm1_bounded(
                1000003, config=cfg, budget=allowance()
            )
            expected = complete.checkpoint["payload"]["state"]
            state = pm1_bounded._initial(1000003, 2, cfg)
            context = SieveContext(128, segment_size=3)
            ledger = allowance()
            for actions in range(expected["steps"] + 1):
                before = copy.deepcopy(state)
                for refused in (
                    allowance(0),
                    allowance(cancelled=lambda: True),
                ):
                    with self.assertRaises(pm1_bounded.BudgetExhaustedError):
                        pm1_bounded._advance(state, refused, context, cfg)
                    self.assertEqual(state, before)
                paused = pm1_bounded.factorize_pm1_bounded(
                    1000003,
                    config=cfg,
                    budget=allowance(),
                    max_actions=actions,
                )
                checkpoint = json.loads(json.dumps(paused.checkpoint))
                resumed = pm1_bounded.factorize_pm1_bounded(
                    1000003,
                    config=cfg,
                    budget=allowance(),
                    checkpoint=checkpoint,
                )
                self.assertEqual(
                    resumed.checkpoint["payload"]["state"], expected
                )
                self.assertGreaterEqual(
                    resumed.work_used,
                    paused.work_used + resumed.verification_work,
                )
                if not state["done"]:
                    pm1_bounded._advance(state, ledger, context, cfg)

    def test_every_saturation_replay_boundary_resumes(self):
        cases = (
            (13 * 19, configured(chunk_bits=32, chunk_size=256)),
            (
                59 * 311,
                replace(configured(wheel=30, gcd_batch=1), bounds=((10, 31),)),
            ),
            (607 * 103, configured(wheel=210, gcd_batch=64)),
        )
        for n, cfg in cases:
            complete = pm1_bounded.factorize_pm1_bounded(
                n, config=cfg, budget=allowance()
            )
            expected = complete.checkpoint["payload"]["state"]
            for actions in range(expected["steps"] + 1):
                paused = pm1_bounded.factorize_pm1_bounded(
                    n, config=cfg, budget=allowance(), max_actions=actions
                )
                resumed = pm1_bounded.factorize_pm1_bounded(
                    n,
                    config=cfg,
                    budget=allowance(),
                    checkpoint=json.loads(json.dumps(paused.checkpoint)),
                )
                self.assertEqual(
                    resumed.checkpoint["payload"]["state"], expected
                )
                self.assertEqual(resumed.result.reconstruct(), n)

    def test_increased_old_power_and_table_invalidation(self):
        cfg = configured(wheel=30, gap_mode="recurrence")
        cfg = replace(cfg, bounds=((15, 15), (16, 16)))
        outcome = pm1_bounded.factorize_pm1_bounded(
            17 * 1019, base=3, config=cfg, budget=allowance()
        )
        self.assertEqual(outcome.divisor, 17)
        cfg = replace(cfg, bounds=((7, 43), (10, 101)))
        state = pm1_bounded._initial(1000003, 2, cfg)
        context = SieveContext(102, segment_size=3)
        while state["rung"] == 0:
            pm1_bounded._advance(state, allowance(), context, cfg)
        self.assertFalse(any(k.startswith("wheel_") for k in state))
        self.assertEqual(state["old_b1"], 7)

    def test_resigned_corruption_configuration_identity_and_cumulative_grant(
        self,
    ):
        cfg = configured(wheel=30)
        paused = pm1_bounded.factorize_pm1_bounded(
            1000003, config=cfg, budget=allowance(), max_actions=20
        )
        for changed in (
            replace(cfg, wheel=210),
            pm1_bounded.PM1Config(bounds=cfg.bounds, segment_size=3),
        ):
            with self.assertRaises(ValueError):
                pm1_bounded.factorize_pm1_bounded(
                    1000003,
                    config=changed,
                    budget=allowance(),
                    checkpoint=paused.checkpoint,
                )
        corrupt = copy.deepcopy(paused.checkpoint)
        corrupt["payload"]["state"]["value"] += 1
        with self.assertRaises(ValueError):
            pm1_bounded.factorize_pm1_bounded(
                1000003,
                config=cfg,
                budget=allowance(),
                checkpoint=resign(corrupt),
            )
        limited = pm1_bounded.factorize_pm1_bounded(
            1000003,
            config=cfg,
            budget=allowance(paused.work_used + 1),
            checkpoint=paused.checkpoint,
        )
        self.assertEqual(limited.reason, "work_limit")
        self.assertEqual(
            limited.checkpoint["payload"]["state"],
            paused.checkpoint["payload"]["state"],
        )
        self.assertGreater(limited.work_used, paused.work_used)

    def test_resigned_table_corruption_is_rejected(self):
        cfg = configured(wheel=30)
        for actions in range(1, 100):
            paused = pm1_bounded.factorize_pm1_bounded(
                1000003, config=cfg, budget=allowance(), max_actions=actions
            )
            if "wheel_forward" in paused.checkpoint["payload"]["state"]:
                break
        else:
            self.fail("paired setup was never reached")
        corrupt = copy.deepcopy(paused.checkpoint)
        corrupt["payload"]["state"]["wheel_forward"][0] += 1
        with self.assertRaises(ValueError):
            pm1_bounded.factorize_pm1_bounded(
                1000003,
                config=cfg,
                budget=allowance(),
                checkpoint=resign(corrupt),
            )

    def test_quiet_result_reconstruction_and_finite_workspace(self):
        output = io.StringIO()
        with redirect_stdout(output):
            result = pm1_bounded.factorize_pm1_bounded(
                607 * 1019, config=configured(wheel=30), budget=allowance()
            )
        self.assertEqual(output.getvalue(), "")
        self.assertEqual(result.result.reconstruct(), 607 * 1019)
        self.assertEqual(result.result.factors, ())
        self.assertFalse(result.result.proven)
        for options in (
            {"wheel": 60},
            {"gap_entries": 257},
            {"chunk_bits": 31},
            {"memory_bytes": 65536, "wheel": 210},
        ):
            with self.assertRaises((ValueError, MemoryError)):
                configured(**options)


if __name__ == "__main__":
    unittest.main()
