"""Independent affine, composite and resume controls for bounded B4 kernels."""

import json
import unittest
from math import gcd

from v2 import arithmetic, ecm, portfolio
from v2.benchmarks import b4_common as common
from v2.benchmarks import b4_kernels as kernels
from v2.benchmarks.prac_oracle import (
    affine_multiply,
    historical_points,
    matches,
)


def primitive(point, n):
    return gcd(gcd(*map(int, point)), int(n)) == 1


class B4KernelTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.engines = {arm: kernels.engine(arm) for arm in kernels.ARMS}
        cls.reference = kernels.ecm_module(cls.engines["baseline"])
        cls.engines["production"] = portfolio

    def test_independent_affine_fields_and_prime_powers(self):
        for prime, curve_a, point in historical_points():
            for modulus in (prime, prime * prime):
                x, y = point
                if modulus != prime:
                    polynomial = x**3 + curve_a * x * x + x
                    correction = (polynomial - y * y) // prime
                    y += prime * (correction * pow(2 * y, -1, prime) % prime)
                a24 = (curve_a + 2) * pow(4, -1, modulus) % modulus

                for scalar in (0, 1, 2, 3, 5, 7, 11, 17, 31, 63):
                    try:
                        expected = affine_multiply(
                            scalar, (x, y), modulus, curve_a
                        )
                    except ValueError:
                        # The affine oracle cannot divide a nonunit over p².
                        # Its field reductions remain independent controls.
                        continue
                    baseline = self.reference.scalar_multiply(
                        scalar, x, 1, modulus, a24
                    )
                    for arm, engine in self.engines.items():
                        candidate = kernels.ecm_module(engine).scalar_multiply(
                            scalar, x, 1, modulus, a24
                        )
                        with self.subTest(arm=arm, n=modulus, k=scalar):
                            if primitive(baseline, modulus):
                                self.assertTrue(
                                    matches(candidate, expected, modulus)
                                )
                            else:
                                # Degeneracy is an exceptional outcome,
                                # never accepted via a zero cross product.
                                self.assertFalse(primitive(candidate, modulus))

    def test_independent_crt_composite_and_unit_scaling(self):
        p, q, curve_a = 1009, 10007, 6
        points = {
            prime: point
            for prime, a, point in historical_points()
            if prime in (p, q) and a == curve_a
        }
        n = p * q

        def crt(a, b):
            return (a + p * ((b - a) * pow(p, -1, q) % q)) % n

        point = (
            crt(points[p][0], points[q][0]),
            crt(points[p][1], points[q][1]),
        )
        a24 = (curve_a + 2) * pow(4, -1, n) % n
        for scalar in (0, 1, 2, 3, 7, 19, 41, 127):
            expected = {
                prime: affine_multiply(scalar, points[prime], prime, curve_a)
                for prime in (p, q)
            }
            for arm, engine in self.engines.items():
                for scale in (1, 2, 17):
                    candidate = kernels.ecm_module(engine).scalar_multiply(
                        scalar, point[0] * scale % n, scale, n, a24
                    )
                    for prime in (p, q):
                        self.assertTrue(
                            matches(candidate, expected[prime], prime),
                            (arm, scalar, scale, prime),
                        )

    def test_convention_and_nonunit_normalization(self):
        # AA=BB+E proves the (A-2)/4 convention needs AA in its bracket.
        for x, z, a24, n in ((39, 17, 29, 1009), (107, 51, 81, 1009**2)):
            aa, bb = (x + z) ** 2, (x - z) ** 2
            delta = aa - bb
            self.assertEqual(
                delta * (bb + a24 * delta) % n,
                delta * (aa + (a24 - 1) * delta) % n,
            )
            self.assertNotEqual(
                delta * (aa + a24 * delta) % n, delta * (bb + a24 * delta) % n
            )

        for n, z, divisor in ((35, 5, 5), (49, 7, 7), (35, 0, 35)):
            with self.assertRaises(
                kernels.NormalizationNonunitError
            ) as caught:
                kernels.normalize_difference(3, z, n)
            self.assertEqual(caught.exception.divisor, divisor)
            self.assertEqual(
                caught.exception.factor, divisor if divisor < n else None
            )
        self.assertEqual(kernels.normalize_difference(9, 3, 35), (3, 1))
        self.assertFalse(matches((0, 0), None, 1009))
        self.assertFalse(matches((0, 0), (13, 15), 1009))
        self.assertFalse(matches((7, 7), (1, 1), 49))

    def test_both_representations_and_exact_formula_oracle(self):
        for name in ("python-int", "gmpy2-mpz"):
            try:
                integer = arithmetic.get_backend(name).integer
            except arithmetic.BackendUnavailableError:
                continue
            for n in (1009, 1009**2, 1009 * 1013, 2**127 - 1):
                setup = ecm.setup_curve(n, 7)
                if setup.point is None:
                    continue
                for scalar in (0, 1, 2, 3, 7, 19, 127, 2**64 + 19):
                    expected = self.reference.scalar_multiply(
                        scalar, *setup.point, n, setup.a24
                    )
                    for arm, engine in self.engines.items():
                        candidate = kernels.ecm_module(engine).scalar_multiply(
                            scalar,
                            *map(integer, setup.point),
                            integer(n),
                            integer(setup.a24),
                        )
                        self.assertEqual(
                            gcd(int(candidate[1]), n), gcd(expected[1], n)
                        )
                        if arm != "normalized":
                            self.assertEqual(
                                tuple(map(int, candidate)), expected
                            )
                        elif primitive(expected, n):
                            self.assertTrue(primitive(candidate, n))
                            self.assertEqual(
                                (
                                    candidate[0] * expected[1]
                                    - expected[0] * candidate[1]
                                )
                                % n,
                                0,
                            )
                        else:
                            self.assertFalse(primitive(candidate, n))
                        self.assertTrue(
                            all(arithmetic.is_mpz(c) for c in candidate)
                            if name == "gmpy2-mpz"
                            else all(type(c) is int for c in candidate)
                        )

    def test_validation_and_low_order_degeneracy(self):
        for engine in self.engines.values():
            multiply = kernels.ecm_module(engine).scalar_multiply
            for scalar, n, x, z in (
                (-1, 101, 2, 1),
                (1, 1, 2, 1),
                (1, 101, 0, 0),
            ):
                with self.assertRaises(ValueError):
                    multiply(scalar, x, z, n, 2)
            with self.assertRaises(TypeError):
                multiply(True, 2, 1, 101, 2)
            self.assertEqual(multiply(0, 7, 0, 101, 2), (1, 0))
            low_order = multiply(3, 0, 1, 101, 2)
            self.assertFalse(primitive(low_order, 101))

    def test_canonical_resume_and_cancellation(self):
        protocol, corpus = common.inputs()
        fixture = next(f for f in corpus["fixtures"] if f["id"] == "stage_128")
        for backend in ("python-int", "gmpy2-mpz"):
            try:
                arithmetic.get_backend(backend)
            except arithmetic.BackendUnavailableError:
                continue
            baseline = None
            for arm, engine in self.engines.items():
                config = common.config(
                    engine, protocol, fixture["case"], backend
                )
                module = kernels.ecm_module(engine)
                original = module.scalar_multiply
                progress = []

                def stop_after_scalar(*args):
                    point = original(*args)
                    progress.append(point)
                    return point

                module.scalar_multiply = stop_after_scalar
                try:
                    stopped = engine.factorize_bounded(
                        fixture["n"],
                        seed=7,
                        config=config,
                        budget=engine.Budget(
                            work_limit=protocol["work_limit"],
                            cancelled=lambda: bool(progress),
                        ),
                    )
                finally:
                    module.scalar_multiply = original
                self.assertEqual(stopped.reason, "cancelled")
                self.assertEqual(len(progress), 1)
                job = stopped.checkpoint["payload"]["state"]["current"]["job"]
                self.assertEqual(job["phase"], "stage_one")
                self.assertEqual(job["value"], list(map(int, progress[0])))
                self.assertEqual(stopped.result.reconstruct(), fixture["n"])
                checkpoint = json.loads(json.dumps(stopped.checkpoint))
                resumed = engine.factorize_bounded(
                    fixture["n"],
                    config=config,
                    checkpoint=checkpoint,
                    budget=engine.Budget(work_limit=protocol["work_limit"]),
                )
                uninterrupted = engine.factorize_bounded(
                    fixture["n"],
                    seed=7,
                    config=config,
                    budget=engine.Budget(work_limit=protocol["work_limit"]),
                )
                self.assertEqual(
                    common.validate_run(resumed, fixture),
                    common.validate_run(uninterrupted, fixture),
                )
                # Canonical X:Z and the existing a24 convention also resume in
                # the unmodified engine; no candidate identity is serialized.
                for plain in (self.engines["baseline"], portfolio):
                    cross = plain.factorize_bounded(
                        fixture["n"],
                        config=common.config(
                            plain, protocol, fixture["case"], backend
                        ),
                        checkpoint=checkpoint,
                        budget=plain.Budget(work_limit=protocol["work_limit"]),
                    )
                    self.assertEqual(
                        common.validate_run(cross, fixture),
                        common.validate_run(uninterrupted, fixture),
                    )
                cancelled = engine.factorize_bounded(
                    fixture["n"],
                    config=config,
                    budget=engine.Budget(
                        work_limit=5000, cancelled=lambda: True
                    ),
                )
                self.assertEqual(cancelled.result.reconstruct(), fixture["n"])
                self.assertFalse(cancelled.result.complete)
                signature = common.validate_run(uninterrupted, fixture)
                if baseline is None:
                    baseline = signature
                elif arm != "normalized":
                    self.assertEqual(signature, baseline)

    def test_nonunit_is_recovered_after_reserved_scalar_action(self):
        engine = self.engines["normalized"]
        package = __import__(engine.__package__, fromlist=["stage_jobs"])
        config = engine.PortfolioConfig(
            trial_bound=5,
            rho_attempts=0,
            pm1_attempts=0,
            ecm_tiers=((2, 2, 1),),
        )
        for n, x, z, factor in ((35, 3, 5, 5), (49, 3, 7, 7)):
            job = package.stage_jobs.new_job("ecm", n, 7, 2, 2)
            job.update(phase="stage_one", value=[x, z], a24=2, powers=[[2, 2]])
            job["cursor"].update(next=3)
            budget = engine.Budget(work_limit=100)
            package.stage_jobs.advance_job(job, budget, None, config)
            self.assertEqual(job["factor"], factor)
            self.assertTrue(job["done"])
            self.assertEqual(budget.used, 3)
            self.assertEqual(job["value"], [x, z])

        for engine in self.engines.values():
            package = __import__(engine.__package__, fromlist=["stage_jobs"])
            job = package.stage_jobs.new_job("ecm", 35, 7, 2, 2)
            job.update(phase="stage_two", terms=[5, 7], product=0)
            budget = engine.Budget(work_limit=100)

            # Saturation must replay and retain the proper term factor.
            package.stage_jobs._batch_check(job, budget)
            self.assertEqual(job["phase"], "term_replay")
            self.assertIsNone(job["factor"])
            package.stage_jobs._batch_check(job, budget)
            self.assertEqual(job["factor"], 5)
            self.assertTrue(job["done"])

    def test_isolation(self):
        self.assertIsNot(kernels.ecm_module(self.engines["squares"]), ecm)
        self.assertIsNot(self.engines["baseline"], portfolio)
        self.assertEqual(ecm.scalar_multiply.__module__, "v2.ecm")


if __name__ == "__main__":
    unittest.main()
