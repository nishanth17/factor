"""Independent CF minimum, scalar, projective and recovery verification."""

import json
import unittest
from dataclasses import replace

from v2 import ecm, prac
from v2.benchmarks import c6_cf, c6_chains, c6_fast, c6_study
from v2.benchmarks.p41_campaign import PYTHON_BACKEND, gmp_backend
from v2.benchmarks.prac_oracle import (
    affine_multiply,
    historical_points,
    matches,
    twist_point,
)


class CFTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.catalog = c6_cf.load_catalog()

    def test_independent_restricted_minima_and_metadata(self):
        data = json.loads(c6_cf.DATA.read_text())
        self.assertEqual(len(data["primes"]), 303)
        for text, values in data["primes"].items():
            if int(text) > 2:
                c6_cf.verify_minimum(int(text), len(values) - 1)
        c6_cf.verify_minimum(29, 7)
        with self.assertRaises(ValueError):
            c6_cf.verify_minimum(29, 8)
        with self.assertRaises(ValueError):
            c6_cf.verify_minimum(2001, 15)
        record = self.catalog[13]
        for changes in (dict(bits=b"\x02"), dict(exponent=2), dict(prime=11)):
            with self.assertRaises(ValueError):
                replace(record, **changes)
        for values in ([1, 3], [1, 2, 3, 6], [1, 2, 3, True]):
            with self.assertRaises(ValueError):
                c6_cf.decode(values)

    def test_every_record_and_three_point_action(self):
        controls = tuple(historical_points())
        for scalar, record in self.catalog.items():
            actions = (
                c6_fast.FastRecord(record.certified, PYTHON_BACKEND),
                c6_cf.ThreePoint(record, PYTHON_BACKEND),
            )
            self.assertLessEqual(record.certified.record.slots, 3)
            for modulus, curve_a, point in controls:
                a24 = (curve_a + 2) * pow(4, -1, modulus) % modulus
                expected = affine_multiply(scalar, point, modulus, curve_a)
                for action in actions:
                    self.assertTrue(
                        matches(
                            action((point[0], 1), modulus, a24),
                            expected,
                            modulus,
                        )
                    )

    def test_composites_prime_squares_and_saturated_product(self):
        recovered = 0
        for n in (35, 101 * 103, 1009**2):
            points = [((5, 7), 2)] if n == 35 else []
            for sigma in range(6, 12):
                setup = ecm.setup_curve(n, sigma)
                if setup.point is not None:
                    points.append((setup.point, setup.a24))
            for point, a24 in points:
                for record in list(self.catalog.values())[:60]:
                    reference = c6_chains.Executor(
                        record.certified.record, PYTHON_BACKEND
                    )
                    actions = (
                        c6_fast.FastRecord(record.certified, PYTHON_BACKEND),
                        c6_cf.ThreePoint(record, PYTHON_BACKEND),
                    )
                    try:
                        expected = reference(point, n, a24)
                    except prac.NonunitPointError as found:
                        for action in actions:
                            with self.assertRaises(
                                prac.NonunitPointError
                            ) as other:
                                action(point, n, a24)
                            self.assertEqual(
                                other.exception.factor, found.factor
                            )
                        recovered += 1
                    else:
                        for action in actions:
                            self.assertEqual(action(point, n, a24), expected)
        self.assertGreater(recovered, 0)

    def test_full_stage_false_infinity_and_finite_replay(self):
        n, sigma = 33554520197234177, 2046841451
        setup = ecm.setup_curve(n, sigma)
        point, curve_a, curve_b = twist_point(setup.point, setup.a24, n)
        expected = affine_multiply(
            ecm.stage_one_scalar(373), point, n, curve_a, curve_b
        )
        for mode in ("tuple", "three"):
            for batch in (1, 16, 64):
                program = c6_cf.build_program(
                    373, mode, PYTHON_BACKEND, batch, self.catalog
                )
                point, factor = program(setup.point, n, setup.a24)
                self.assertIsNone(factor)
                self.assertTrue(matches(point, expected, n))
                self.assertNotEqual(point[1], 0)
                extra = {}
                self.assertEqual(program((5, 7), 35, 2, extra), (None, 7))
                self.assertEqual(extra["block_replays"], 1)
                self.assertLessEqual(extra["record_replay_allowance"], batch)
        case = next(
            c
            for c in c6_study.control.load_corpus()["fixtures"]
            if c["id"] == "balanced_40d"
        )
        from v2 import constants, utils

        rng = utils.resolve_rng(41001, None)
        setup = ecm.setup_curve(
            case["n"], rng.randint(6, constants.MAX_RANDOM_ECM)
        )
        expected = c6_study.affine_controls(
            case["n"], tuple(case["factors"]), 41001, 2000
        )
        for mode in ("tuple", "three"):
            program = c6_cf.build_program(
                2000, mode, PYTHON_BACKEND, 16, self.catalog
            )
            point, factor = program(setup.point, case["n"], setup.a24)
            self.assertIsNone(factor)
            for modulus, target in expected:
                self.assertTrue(matches(point, target, modulus))

    def test_gmp_matches_int(self):
        try:
            backend = gmp_backend()
        except ImportError:
            self.skipTest("optional GMP unavailable")
        for record in list(self.catalog.values())[:60]:
            other = c6_cf.ThreePoint(record, backend)
            native = c6_cf.ThreePoint(record, PYTHON_BACKEND)
            expected = native((3, 1), 1009, 2)
            actual = other(
                (backend.integer(3), backend.integer(1)),
                backend.integer(1009),
                backend.integer(2),
            )
            self.assertEqual(tuple(map(int, actual)), expected)


if __name__ == "__main__":
    unittest.main()
