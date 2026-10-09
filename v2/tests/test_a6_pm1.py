"""Independent exponent, boundary, recovery and cumulative-resume checks."""

import copy
import hashlib
import io
import json
import math
import unittest
from contextlib import redirect_stdout
from dataclasses import replace

from v2 import pm1_bounded, schedules, utils
from v2.budget import Budget
from v2.pm1_bounded import PM1Config, factorize_pm1_bounded


def allowance(work=2_000_000, **options):
    return Budget(work_limit=work, seconds=None, cpu_seconds=None, **options)


def config(bounds=((10, 200),), **options):
    return PM1Config(bounds=bounds, segment_size=8, **options)


def resign(checkpoint):
    encoded = json.dumps(
        checkpoint["payload"], sort_keys=True, separators=(",", ":")
    )
    checkpoint["sha256"] = hashlib.sha256(encoded.encode()).hexdigest()
    return checkpoint


class RatioTests(unittest.TestCase):
    def test_exhaustive_lcm_identity_without_sieve_or_prime_power_oracle(self):
        for new in range(2, 65):
            new_lcm = math.lcm(*range(1, new + 1))
            for old in range(1, new + 1):
                old_lcm = math.lcm(*range(1, old + 1))
                actions = list(
                    schedules.prime_power_ratios(old, new, segment_size=3)
                )

                self.assertEqual(
                    old_lcm * math.prod(power for _, power in actions), new_lcm
                )
                self.assertTrue(all(power > 1 for _, power in actions))
                self.assertEqual(
                    [prime for prime, _ in actions],
                    sorted(set(prime for prime, _ in actions)),
                )

    def test_old_prime_power_increases_and_new_primes(self):
        self.assertEqual(
            list(schedules.prime_power_ratios(7, 10)), [(2, 2), (3, 3)]
        )
        self.assertEqual(
            list(schedules.prime_power_ratios(15, 17)), [(2, 2), (17, 17)]
        )
        self.assertEqual(list(schedules.prime_power_ratios(16, 16)), [])
        self.assertEqual(schedules.prime_power_ratio(2, 2**80 - 1, 2**80), 2)
        for old, new in ((0, 10), (10, 9)):
            with self.assertRaises((TypeError, ValueError)):
                list(schedules.prime_power_ratios(old, new))

    def test_direct_small_exponents_include_nonunits(self):
        for modulus in (9, 15, 35, 101):
            for base in range(modulus):
                for old, new in ((1, 7), (7, 10), (15, 17)):
                    old_lcm = math.lcm(*range(1, old + 1))
                    new_lcm = math.lcm(*range(1, new + 1))
                    value = pow(base, old_lcm, modulus)
                    for _, ratio in schedules.prime_power_ratios(old, new):
                        value = pow(value, ratio, modulus)
                    self.assertEqual(value, pow(base, new_lcm, modulus))


class CampaignTests(unittest.TestCase):
    def test_stage_one_and_every_stage_two_boundary_against_direct_powers(
        self,
    ):
        cfg = config(((7, 13), (10, 31), (10, 43)), chunk_size=3, gcd_batch=2)
        n, base = 1000003, 2
        state = pm1_bounded._initial(n, base, cfg)
        context = schedules.SieveContext(44, segment_size=8)
        budget = allowance()

        while not state["done"]:
            pm1_bounded._advance(state, budget, context, cfg)
            if state["phase"] == "stage_two":
                exponent = math.lcm(*range(1, state["b1"] + 1))
                value = pow(base, exponent, n)
                self.assertEqual(state["value"], value)
                self.assertEqual(
                    state["stage_two_value"],
                    pow(value, state["previous_prime"], n),
                )
                self.assertEqual(
                    state["product"], math.prod(state["terms"]) % n
                )
                for gap, power in state["gap_powers"].items():
                    self.assertEqual(power, pow(value, int(gap), n))
        self.assertEqual(state["reason"], "exhausted")

    def test_old_power_boundary_just_outside_and_inside(self):
        n = 17 * 1019
        self.assertEqual(pow(3, 16, 17), 1)
        self.assertNotEqual(pow(3, 8, 17), 1)
        outside = factorize_pm1_bounded(
            n, base=3, config=config(((15, 15),)), budget=allowance()
        )
        inside = factorize_pm1_bounded(
            n, base=3, config=config(((15, 15), (16, 16))), budget=allowance()
        )

        self.assertIsNone(outside.divisor)
        self.assertEqual(inside.divisor, 17)
        self.assertEqual(inside.result.reconstruct(), n)
        self.assertFalse(inside.result.proven)
        self.assertEqual(inside.result.factors, ())

    def test_stage_two_only_factor_just_inside_and_outside(self):
        n = 607 * 1019
        outside = factorize_pm1_bounded(
            n, config=config(((10, 100),)), budget=allowance()
        )
        inside = factorize_pm1_bounded(
            n, config=config(((10, 100), (10, 101))), budget=allowance()
        )

        self.assertIsNone(outside.divisor)
        self.assertEqual(inside.divisor, 607)
        self.assertEqual(inside.checkpoint["payload"]["state"]["rung"], 1)

    def test_increased_b1_rebuilds_stage_two_arithmetic(self):
        cfg = config(((7, 13), (10, 31)))
        state = pm1_bounded._initial(1000003, 2, cfg)
        context = schedules.SieveContext(32, segment_size=8)
        while state["rung"] == 0:
            pm1_bounded._advance(state, allowance(), context, cfg)

        self.assertEqual(state["old_b1"], 7)
        self.assertEqual(state["gap_powers"], {})
        self.assertEqual(state["previous_prime"], 0)
        self.assertEqual(state["stage_two_value"], 1)
        self.assertEqual(state["cursor"]["hi"], 11)

    def test_mixed_factor_saturation_and_finite_recovery(self):
        n = 13 * 19
        recovered = factorize_pm1_bounded(
            n, config=config(((10, 10),)), budget=allowance()
        )
        refused = factorize_pm1_bounded(
            n,
            config=config(((10, 10), (16, 16)), recovery_limit=0),
            budget=allowance(),
        )

        self.assertTrue(utils.valid_divisor(recovered.divisor, n))
        self.assertEqual(refused.reason, "saturated")
        self.assertIsNone(refused.divisor)
        self.assertEqual(refused.checkpoint["payload"]["state"]["rung"], 0)
        self.assertEqual(refused.result.remaining, (n,))

    def test_stage_two_mixed_batch_recovery(self):
        # 607-1=6*101 and 103-1=6*17: neither divides M(10), while
        # different stage-two relations annihilate the two factors.
        n = 607 * 103
        outcome = factorize_pm1_bounded(
            n, config=config(((10, 101),), gcd_batch=64), budget=allowance()
        )
        self.assertTrue(utils.valid_divisor(outcome.divisor, n))
        self.assertEqual(
            outcome.checkpoint["payload"]["state"]["phase"], "term_replay"
        )

    def test_stage_two_recovery_cap_retains_unresolved_input(self):
        n = 607 * 103
        outcome = factorize_pm1_bounded(
            n,
            config=config(((10, 101),), gcd_batch=64, recovery_limit=0),
            budget=allowance(),
        )
        state = outcome.checkpoint["payload"]["state"]

        self.assertEqual(outcome.reason, "saturated")
        self.assertEqual(state["phase"], "term_replay")
        self.assertEqual(state["recovery"], 0)
        self.assertIsNone(outcome.divisor)
        self.assertEqual(outcome.result.remaining, (n,))
        self.assertEqual(outcome.result.reconstruct(), n)

    def test_every_atomic_refusal_preserves_serializable_state(self):
        cfg = config(((7, 13), (10, 31)), chunk_size=2, gcd_batch=2)
        state = pm1_bounded._initial(1000003, 2, cfg)
        context = schedules.SieveContext(32, segment_size=8)
        from v2.budget import BudgetExhaustedError

        while not state["done"]:
            before = copy.deepcopy(state)
            with self.assertRaises(BudgetExhaustedError):
                pm1_bounded._advance(state, allowance(0), context, cfg)
            self.assertEqual(state, before)
            pm1_bounded._advance(state, allowance(), context, cfg)

    def test_cancellation_during_each_execution_phase_is_resumable(self):
        cfg = config(((7, 13), (10, 31)), chunk_size=2, gcd_batch=2)
        for poll in range(2, 25):
            calls = 0

            def cancel():
                nonlocal calls
                calls += 1
                return calls == poll

            paused = factorize_pm1_bounded(
                1000003, config=cfg, budget=allowance(cancelled=cancel)
            )
            resumed = factorize_pm1_bounded(
                1000003,
                config=cfg,
                budget=allowance(),
                checkpoint=paused.checkpoint,
            )
            self.assertEqual(resumed.reason, "exhausted")
            self.assertEqual(resumed.result.reconstruct(), 1000003)

    def test_nonunits_even_input_and_no_certainty_upgrade(self):
        for n, base, divisor in ((35, 5, 5), (14, 3, 2)):
            result = factorize_pm1_bounded(
                n, base=base, config=config(), budget=allowance()
            )
            self.assertEqual(result.divisor, divisor)
            self.assertEqual(result.result.reconstruct(), n)
        result = factorize_pm1_bounded(
            35, base=35, config=config(), budget=allowance()
        )
        self.assertEqual(result.reason, "nonunit")
        self.assertEqual(result.result.remaining, (35,))

    def test_library_is_quiet_and_exhaustion_is_final(self):
        cfg = config(((2, 2),))
        output = io.StringIO()
        with redirect_stdout(output):
            first = factorize_pm1_bounded(
                1019 * 1237, config=cfg, budget=allowance()
            )
            resumed = factorize_pm1_bounded(
                1019 * 1237,
                config=cfg,
                budget=allowance(),
                checkpoint=first.checkpoint,
            )
        self.assertEqual(output.getvalue(), "")
        self.assertEqual(first.reason, "exhausted")
        self.assertEqual(
            resumed.checkpoint["payload"]["state"],
            first.checkpoint["payload"]["state"],
        )

    def test_finite_input_and_storage_configuration(self):
        for kwargs in (
            {"bounds": ()},
            {"bounds": ((10, 9),)},
            {"bounds": ((10, 20), (9, 30))},
            {"bounds": ((10, 20), (10, 20))},
            {"chunk_size": 257},
        ):
            with self.assertRaises((ValueError, MemoryError)):
                PM1Config(**kwargs)
        with self.assertRaises(MemoryError):
            PM1Config(memory_bytes=65536)
        with self.assertRaises(ValueError):
            factorize_pm1_bounded(
                2**100, config=replace(config(), max_input_bits=10)
            )


class ResumeTests(unittest.TestCase):
    def test_every_boundary_reconstructs_same_state_and_charges_replay(self):
        cfg = config(((7, 13), (10, 31)), chunk_size=2, gcd_batch=2)
        n = 1000003
        full = factorize_pm1_bounded(n, config=cfg, budget=allowance())
        total = full.checkpoint["payload"]["state"]["steps"]
        for boundary in range(total + 1):
            paused = factorize_pm1_bounded(
                n, config=cfg, budget=allowance(), max_actions=boundary
            )
            saved = json.loads(json.dumps(paused.checkpoint))
            resumed = factorize_pm1_bounded(
                n, config=cfg, budget=allowance(), checkpoint=saved
            )

            self.assertEqual(
                resumed.checkpoint["payload"]["state"],
                full.checkpoint["payload"]["state"],
            )
            self.assertEqual(resumed.divisor, full.divisor)
            self.assertGreaterEqual(resumed.work_used, full.work_used)
            self.assertEqual(
                resumed.verification_work,
                paused.checkpoint["payload"]["state"]["execution_work"],
            )
            self.assertEqual(resumed.result.reconstruct(), n)

    def test_refused_actions_and_repeated_resume_cannot_reset_allowances(self):
        cfg = config()
        first = factorize_pm1_bounded(
            607 * 1019, config=cfg, budget=allowance(40)
        )
        second = factorize_pm1_bounded(
            607 * 1019,
            config=cfg,
            budget=allowance(40),
            checkpoint=first.checkpoint,
        )
        self.assertEqual(first.reason, "work_limit")
        self.assertEqual(second.reason, "work_limit")
        self.assertLessEqual(second.work_used, 40)
        self.assertEqual(
            second.checkpoint["payload"]["state"],
            first.checkpoint["payload"]["state"],
        )
        final = factorize_pm1_bounded(
            607 * 1019,
            config=cfg,
            budget=allowance(),
            checkpoint=second.checkpoint,
        )
        self.assertEqual(final.divisor, 607)
        with self.assertRaises(ValueError):
            factorize_pm1_bounded(
                607 * 1019,
                config=cfg,
                budget=allowance(first.work_used - 1),
                checkpoint=first.checkpoint,
            )

    def test_cancellation_and_cumulative_wall_cpu_limits(self):
        cfg = config()
        first = factorize_pm1_bounded(
            607 * 1019, config=cfg, budget=allowance(), max_actions=3
        )
        for budget, reason in (
            (allowance(cancelled=lambda: True), "cancelled"),
            (Budget(seconds=0, cpu_seconds=None), "wall_limit"),
            (Budget(seconds=None, cpu_seconds=0), "cpu_limit"),
        ):
            resumed = factorize_pm1_bounded(
                607 * 1019,
                config=cfg,
                budget=budget,
                checkpoint=first.checkpoint,
            )
            self.assertEqual(resumed.reason, reason)
            self.assertEqual(resumed.result.remaining, (607 * 1019,))
            self.assertGreaterEqual(resumed.wall_seconds, first.wall_seconds)
            self.assertGreaterEqual(resumed.cpu_seconds, first.cpu_seconds)

    def test_resigned_corrupt_arithmetic_is_rejected(self):
        cfg = config()
        first = factorize_pm1_bounded(
            1000003, config=cfg, budget=allowance(), max_actions=8
        )
        for field, value in (
            ("value", 1),
            ("rung", 20),
            ("base", 3),
            ("steps", -1),
            ("factor", 17),
            ("execution_work", 0),
        ):
            bad = copy.deepcopy(first.checkpoint)
            bad["payload"]["state"][field] = value
            with self.assertRaises(ValueError, msg=field):
                factorize_pm1_bounded(
                    1000003,
                    config=cfg,
                    budget=allowance(),
                    checkpoint=resign(bad),
                )
        bad = copy.deepcopy(first.checkpoint)
        bad["payload"]["state"]["cursor"]["values"] = [4, 9]
        with self.assertRaises(ValueError):
            factorize_pm1_bounded(
                1000003, config=cfg, budget=allowance(), checkpoint=resign(bad)
            )

    def test_corrupt_or_incompatible_identity_is_rejected(self):
        cfg = config()
        first = factorize_pm1_bounded(
            1000003, config=cfg, budget=allowance(), max_actions=8
        )
        for field, value in (
            ("version", 2),
            ("version", True),
            ("schedule", "future"),
            ("execution", "future"),
            ("base", 3),
            ("n", 1000033),
            ("wall_used", -1),
        ):
            bad = copy.deepcopy(first.checkpoint)
            bad["payload"][field] = value
            with self.assertRaises(ValueError, msg=field):
                factorize_pm1_bounded(
                    1000003,
                    config=cfg,
                    budget=allowance(),
                    checkpoint=resign(bad),
                )
        for options in (
            {"base": 3},
            {"config": replace(cfg, bounds=((11, 200),))},
            {"config": replace(cfg, chunk_size=1)},
        ):
            kwargs = {"config": cfg, **options}
            with self.assertRaises(ValueError):
                factorize_pm1_bounded(
                    1000003,
                    budget=allowance(),
                    checkpoint=first.checkpoint,
                    **kwargs,
                )
        bad = copy.deepcopy(first.checkpoint)
        bad["payload"]["state"]["value"] = 1
        with self.assertRaises(ValueError):
            factorize_pm1_bounded(1000003, config=cfg, checkpoint=bad)


if __name__ == "__main__":
    unittest.main()
