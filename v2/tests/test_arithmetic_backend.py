"""Shared exact oracles and interrupted engines on each available backend."""

import copy
import importlib.util
import io
import json
import math
import random
import subprocess
import sys
import unittest
from contextlib import redirect_stdout
from dataclasses import replace
from unittest.mock import patch

from v2 import arithmetic, ecm, pollard_pm1, pollard_rho, portfolio, utils
from v2.budget import Budget
from v2.factor import factorize
from v2.portfolio import PortfolioConfig, _canonical, factorize_bounded
from v2.qs import SIQSJob, build_factor_base
from v2.qs.checkpoint import _solver_digest
from v2.qs.families import PolynomialFamily, _checksum, family_assignments
from v2.qs.linear_algebra import DependencySolver, filter_matrix
from v2.qs.parallel import CollectionPool, ParallelSIQSJob
from v2.qs.smooth_batch import SmoothBatch
from v2.qs.sss import SSSConfig, SSSJob
from v2.schedules import SieveContext
from v2.stage_jobs import advance_job, new_job, promote_job
from v2.tests.test_ecm_programs import configuration as program_config
from v2.tests.test_phase_two import reseal
from v2.tests.test_qs_parallel import configuration as parallel_config
from v2.tests.test_siqs import configuration as siqs_config

BACKENDS = ("python-int",)
if importlib.util.find_spec("gmpy2") is not None:
    BACKENDS += ("gmpy2-mpz",)


def allowance(work=10**9):
    return Budget(work_limit=work, seconds=None, cpu_seconds=None)


def finish(kind, n, backend, seed=7):
    config = PortfolioConfig(
        backend=backend,
        trial_bound=10,
        rho_evaluations=300,
        pm1_b1=10,
        pm1_b2=200,
        ecm_tiers=((10, 200, 2),),
        chunk_size=2,
        gcd_batch=3,
    )
    job = new_job(
        kind, arithmetic.get_backend(backend).integer(n), seed, 10, 200
    )
    context, budget = SieveContext(201, segment_size=16), allowance()
    while not job["done"]:
        advance_job(job, budget, context, config)
    return arithmetic.canonical(job), budget.used


class ArithmeticTests(unittest.TestCase):
    def test_exact_contracts_and_nonunits(self):
        generator = random.Random(431)
        pairs = [(0, 0), (-15, 0), (0, -7), (-35, 21)]
        pairs += [
            (generator.getrandbits(256), -generator.getrandbits(256))
            for _ in range(32)
        ]

        for name in BACKENDS:
            integer = arithmetic.get_backend(name).integer
            for a, b in pairs:
                a, b = integer(a), integer(b)
                g, x, y = utils.extended_gcd(a, b)
                self.assertEqual(g, math.gcd(int(a), int(b)))
                self.assertEqual(a * x + b * y, g)
                self.assertEqual(arithmetic.gcd(a, b), g)

            for modulus in (15, 101, 2**127 - 1):
                for value in (-9, 0, 1, 7, 15, 23):
                    a, n = integer(value), integer(modulus)
                    divisor = math.gcd(value, modulus)
                    if divisor == 1:
                        self.assertEqual(
                            utils.modular_inverse(a, n),
                            pow(value, -1, modulus),
                        )
                        self.assertEqual(
                            arithmetic.pow(a, 31, n), pow(value, 31, modulus)
                        )
                        self.assertEqual(
                            arithmetic.pow(a, -3, n), pow(value, -3, modulus)
                        )
                    else:
                        with self.assertRaises(
                            arithmetic.NonInvertibleError
                        ) as caught:
                            utils.modular_inverse(a, n)
                        self.assertEqual(caught.exception.divisor, divisor)

            for a, b in ((0, 7), (-63, 7), (63, -7), (2**1024, 16)):
                self.assertEqual(
                    arithmetic.divexact(integer(a), integer(b)), a // b
                )
            with self.assertRaises(ValueError):
                arithmetic.divexact(integer(25), integer(7))
            with self.assertRaises(ZeroDivisionError):
                arithmetic.divexact(integer(1), integer(0))
            for value in (False, True, 1.0, "7"):
                with self.assertRaises(TypeError):
                    integer(value)
            with self.assertRaises(ValueError):
                arithmetic.pow(integer(7), -1)

    def test_root_boundaries_and_certainty(self):
        for name in BACKENDS:
            integer = arithmetic.get_backend(name).integer
            for exponent in (1, 2, 3, 7, 127, 4097):
                values = (0, 1, 2, (1 << 1024) - 1, 1 << 1024)
                for value in values:
                    root = arithmetic.integer_root(integer(value), exponent)
                    self.assertLessEqual(root**exponent, value)
                    self.assertGreater((root + 1) ** exponent, value)

            for value in (0, 1, 15, 16, 17, 2**4096 - 1):
                self.assertEqual(
                    arithmetic.isqrt(integer(value)), math.isqrt(value)
                )
            for value in (97, 3215031751, 2**64 - 59, 2**127 - 1):
                native_rng, candidate_rng = random.Random(9), random.Random(9)
                expected = utils.classify_prime(value, rng=native_rng)
                actual = utils.classify_prime(
                    integer(value), rng=candidate_rng
                )
                self.assertEqual(actual, expected)
                self.assertEqual(
                    candidate_rng.getstate(), native_rng.getstate()
                )

    def test_missing_dependency_is_explicit_and_int_import_is_lazy(self):
        process = subprocess.run(
            [
                sys.executable,
                "-c",
                "import sys; from v2.factor import factorize; "
                "assert factorize(45).reconstruct() == 45; "
                "assert 'gmpy2' not in sys.modules",
            ],
            check=False,
            capture_output=True,
            text=True,
        )
        self.assertEqual(process.returncode, 0, process.stderr)
        with (
            patch.object(arithmetic, "_gmp", None),
            patch.object(
                arithmetic.importlib, "import_module", side_effect=ImportError
            ),
        ):
            self.assertEqual(arithmetic.get_backend().integer(7), 7)
            with self.assertRaisesRegex(RuntimeError, "unavailable"):
                arithmetic.get_backend("gmpy2-mpz")


class EngineBackendTests(unittest.TestCase):
    def test_canonical_divisors_reenter_the_selected_backend(self):
        # Public splitters return native divisors. Both pending children must
        # retain the caller's backend before subsequent primality/arithmetic.
        for name in BACKENDS:
            with patch(
                "v2.factor.utils.classify_prime", wraps=utils.classify_prime
            ) as classify:
                result = factorize(626100403, level=2, seed=7, backend=name)

            self.assertTrue(result.complete)
            self.assertGreaterEqual(classify.call_count, 3)
            for call in classify.call_args_list:
                self.assertEqual(
                    arithmetic.is_mpz(call.args[0]), name == "gmpy2-mpz"
                )

            config = PortfolioConfig(
                backend=name,
                trial_bound=2,
                rho_attempts=0,
                pm1_attempts=0,
                ecm_tiers=(),
                memory_bytes=64 * 2**20,
                siqs=siqs_config(backend=name),
            )
            with patch(
                "v2.portfolio._classify_step",
                wraps=portfolio._classify_step,
            ) as classify:
                run = factorize_bounded(
                    4001 * 5003, seed=7, config=config, budget=allowance()
                )

            self.assertTrue(run.result.complete)
            self.assertGreaterEqual(classify.call_count, 3)
            for call in classify.call_args_list:
                self.assertEqual(
                    arithmetic.is_mpz(call.args[0]["n"]),
                    name == "gmpy2-mpz",
                )

    def test_factor_results_reconstruct_and_stay_canonical(self):
        inputs = (-1, 1, -72, 1009**3 * 1013**2, 626100403, 2**127 - 1)
        for name in BACKENDS:
            output = io.StringIO()
            with redirect_stdout(output):
                for n in inputs:
                    result = factorize(n, seed=7, backend=name)
                    self.assertEqual(result, factorize(n, seed=7))
                    self.assertEqual(result.reconstruct(), n)
                    self.assertIs(type(result.original), int)
                    for factor in result.factors:
                        self.assertIs(type(factor.value), int)
                partial = factorize(626100403, level=0, backend=name)
                self.assertEqual(partial.remaining, (626100403,))
                self.assertIs(type(partial.remaining[0]), int)
            self.assertEqual(output.getvalue(), "")

    def test_rho_pm1_ecm_and_complete_stage_ledgers_match(self):
        for name in BACKENDS:
            for function, options in (
                (pollard_rho.factorize_rho, dict(seed=7)),
                (pollard_pm1.factorize_pm1, dict(b1=10, b2=200)),
                (
                    ecm.factorize_ecm,
                    dict(seed=7, b1=50, b2=1000, max_curves=32),
                ),
            ):
                for n in (35, 1009 * 1013, 607 * 1019):
                    actual = function(n, backend=name, **options)
                    self.assertEqual(actual, function(n, **options))
                    if actual is not None:
                        self.assertIs(type(actual), int)
                        self.assertTrue(utils.valid_divisor(actual, n))
            for kind in ("rho", "pm1", "ecm"):
                for n in (13 * 19, 1009 * 1013, 1000003 * 1000033):
                    self.assertEqual(
                        finish(kind, n, name), finish(kind, n, "python-int")
                    )

    def test_ecm_ladder_and_mixed_batch_recovery(self):
        for name in BACKENDS:
            integer = arithmetic.get_backend(name).integer
            setup = ecm.setup_curve(integer(1009 * 1013), 6)
            for scalar in (0, 1, 2, 7, 2520):
                point = ecm.scalar_multiply(
                    scalar, *setup.point, integer(1009 * 1013), setup.a24
                )
                native = ecm.setup_curve(1009 * 1013, 6)
                expected = ecm.scalar_multiply(
                    scalar, *native.point, 1009 * 1013, native.a24
                )
                self.assertEqual(point, expected)
            self.assertEqual(
                utils.batch_factor([integer(5), integer(7)], integer(35)),
                (5, True),
            )

    def test_bounded_resume_retains_rng_work_and_backend(self):
        for name in BACKENDS:
            config = PortfolioConfig(
                backend=name,
                trial_bound=10,
                rho_attempts=0,
                pm1_b1=50,
                pm1_b2=1000,
                ecm_tiers=((50, 1000, 4),),
            )
            expected = factorize_bounded(
                1009 * 1013, seed=73, config=config, budget=allowance()
            )
            run = factorize_bounded(
                1009 * 1013, seed=73, config=config, budget=allowance(500)
            )
            self.assertEqual(run.result.reconstruct(), 1009 * 1013)
            checkpoint = json.loads(json.dumps(run.checkpoint))
            resumed = factorize_bounded(
                1009 * 1013,
                config=config,
                checkpoint=checkpoint,
                budget=allowance(),
            )
            self.assertEqual(
                (resumed.result, resumed.work_used),
                (expected.result, expected.work_used),
            )
            self.assertEqual(resumed.events, expected.events)
            damaged = copy.deepcopy(checkpoint)
            damaged["payload"]["backend"] = "wrong-backend"
            import hashlib

            damaged["sha256"] = hashlib.sha256(
                _canonical(damaged["payload"]).encode()
            ).hexdigest()
            with self.assertRaisesRegex(ValueError, "incompatible"):
                factorize_bounded(
                    1009 * 1013,
                    config=config,
                    checkpoint=damaged,
                    budget=allowance(),
                )

    def test_ecm_programs_resume_on_each_backend(self):
        n = 25013 * 25031
        reference = None
        for name in BACKENDS:
            config = program_config(backend=name)
            whole = factorize_bounded(
                n, seed=7, config=config, budget=allowance()
            )
            signature = (whole.result, whole.reason, whole.work_used)
            if reference is None:
                reference = signature
            self.assertEqual(signature, reference)

            paused = factorize_bounded(
                n, seed=7, config=config, budget=allowance(1000)
            )
            self.assertEqual(
                paused.checkpoint["payload"]["version"],
                5 if name == "python-int" else 6,
            )
            resumed = factorize_bounded(
                n,
                config=config,
                budget=allowance(),
                checkpoint=json.loads(json.dumps(paused.checkpoint)),
            )
            self.assertEqual(resumed.result, whole.result)
            self.assertEqual(resumed.reason, whole.reason)
            self.assertGreaterEqual(resumed.work_used, paused.work_used)
            self.assertEqual(
                [
                    (event["seed"], event["outcome"])
                    for event in resumed.events
                ],
                [(event["seed"], event["outcome"]) for event in whole.events],
            )

            # Pre-integration GMP version 5 had no program identity. Merely
            # relabeling a new program checkpoint must not make it legacy.
            if name == "gmpy2-mpz":
                damaged = copy.deepcopy(paused.checkpoint)
                damaged["payload"]["version"] = 5
                with self.assertRaisesRegex(ValueError, "incompatible"):
                    factorize_bounded(
                        n,
                        config=config,
                        budget=allowance(),
                        checkpoint=reseal(damaged),
                    )

    def test_preintegration_backend_checkpoints_remain_readable(self):
        n = 25013 * 25031
        for name in BACKENDS:
            config = program_config(backend=name, ecm_program_bytes=0)
            whole = factorize_bounded(
                n, seed=7, config=config, budget=allowance()
            )
            paused = factorize_bounded(
                n, seed=7, config=config, budget=allowance(1000)
            )
            legacy = copy.deepcopy(paused.checkpoint)
            legacy["payload"]["version"] = 5
            legacy["payload"]["config"]["backend"] = name
            resumed = factorize_bounded(
                n,
                config=config,
                budget=allowance(),
                checkpoint=reseal(legacy),
            )
            self.assertEqual(
                (resumed.result, resumed.reason, resumed.work_used),
                (whole.result, whole.reason, whole.work_used),
            )

    def test_verified_prac_uses_selected_arithmetic(self):
        from v2.benchmarks.prac_oracle import affine_multiply, matches

        for name in BACKENDS:
            integer = arithmetic.get_backend(name).integer
            for scalar in (13, 97):
                point = ecm.multiply_prac(
                    scalar, integer(3), integer(1), integer(1009), integer(2)
                )
                expected = affine_multiply(scalar, (3, 293), 1009, 6)
                self.assertTrue(matches(point, expected, 1009))
                self.assertTrue(
                    all(
                        arithmetic.is_mpz(value) == (name == "gmpy2-mpz")
                        for value in point
                    )
                )

    def test_backend_keyword_preserves_positional_configurations(self):
        from v2.qs.parallel import ParallelConfig
        from v2.qs.siqs import SIQSConfig

        self.assertEqual(PortfolioConfig(200).trial_bound, 200)
        self.assertEqual(SIQSConfig("mpqs").mode, "mpqs")
        self.assertEqual(SSSConfig("sssf").mode, "sssf")
        self.assertEqual(ParallelConfig(400).base_bound, 400)

    def test_retained_stage_values_are_rehydrated_without_touching_indices(
        self,
    ):
        for name in BACKENDS:
            backend = arithmetic.get_backend(name)
            job = new_job("ecm", 1022117, 7, 10, 200)
            job.update(value=[7, 9], baby=[None, [7, 9]], a24=11)
            promote_job(job, backend)
            self.assertEqual(arithmetic.is_mpz(job["n"]), name == "gmpy2-mpz")
            self.assertEqual(
                arithmetic.is_mpz(job["value"][0]), name == "gmpy2-mpz"
            )
            self.assertIs(type(job["b1"]), int)
            self.assertIs(type(job["cursor"]["index"]), int)


class QSBackendTests(unittest.TestCase):
    def test_qs_mpqs_siqs_and_sss_resume(self):
        for name in BACKENDS:
            for mode in ("qs", "mpqs", "siqs"):
                config = siqs_config(backend=name, mode=mode)
                job = SIQSJob(
                    4001 * 5003, seed=73, config=config, budget=allowance()
                )
                job.run(max_blocks=1)
                checkpoint = json.loads(json.dumps(job.checkpoint()))
                restored = SIQSJob.from_checkpoint(
                    checkpoint, budget=allowance(), config=config
                )
                expected, actual = job.run(), restored.run()
                self.assertEqual(
                    (actual.reason, actual.divisor),
                    (expected.reason, expected.divisor),
                )
                self.assertEqual(
                    actual.divisor * actual.cofactor
                    if actual.divisor
                    else actual.cofactor,
                    4001 * 5003,
                )
                self.assertEqual(
                    arithmetic.is_mpz(restored.n), name == "gmpy2-mpz"
                )
                if restored.engine is not None:
                    rows = (
                        restored.engine.prepared.rows
                        if restored.engine.prepared
                        else ()
                    )
                    for row in rows:
                        self.assertEqual(
                            arithmetic.is_mpz(row), name == "gmpy2-mpz"
                        )
                other = "gmpy2-mpz" if name == "python-int" else "python-int"
                if other in BACKENDS:
                    with self.assertRaises(ValueError):
                        SIQSJob.from_checkpoint(
                            checkpoint,
                            budget=allowance(),
                            config=replace(config, backend=other),
                        )

            for mode in ("sss", "sssf"):
                config = SSSConfig(backend=name, mode=mode, base_bound=400)
                job = SSSJob(4001 * 4003, config=config, budget=allowance())
                job.run(batch_limit=1)
                restored = SSSJob.from_checkpoint(
                    json.loads(json.dumps(job.checkpoint())),
                    budget=allowance(),
                    config=config,
                )
                actual, expected = restored.run(), job.run()
                self.assertEqual(
                    (actual.reason, actual.divisor),
                    (expected.reason, expected.divisor),
                )

    def test_external_square_and_streamed_siqs(self):
        for name in BACKENDS:
            for options in (
                dict(mode="mpqs", external_coefficients=True),
                dict(assignment_policy="flyer"),
            ):
                config = siqs_config(backend=name, factor_count=1, **options)
                job = SIQSJob(4001 * 5003, config=config, budget=allowance())
                result = job.run()
                native = SIQSJob(
                    4001 * 5003,
                    config=replace(config, backend="python-int"),
                    budget=allowance(),
                ).run()
                self.assertEqual(
                    (result.reason, result.divisor),
                    (native.reason, native.divisor),
                )

    def test_family_identity_and_backend_mismatch(self):
        for name in BACKENDS:
            base = build_factor_base(
                4001 * 5003, bound=200, backend=name
            ).factor_base
            primes = family_assignments(base, 256, family_count=1)[0]
            family = PolynomialFamily(base, primes, budget=allowance())
            family.next()
            checkpoint = json.loads(json.dumps(family.checkpoint()))
            restored = PolynomialFamily.from_checkpoint(
                base, checkpoint, budget=family.budget
            )
            self.assertEqual(restored.next(), family.next())
            damaged = copy.deepcopy(checkpoint)
            damaged["payload"]["backend"] = "wrong-backend"
            damaged["sha256"] = _checksum(damaged["payload"])
            with self.assertRaisesRegex(ValueError, "backend"):
                PolynomialFamily.from_checkpoint(
                    base, damaged, budget=allowance()
                )

    def test_smooth_trees_and_bitsets_use_the_selected_representation(self):
        generator = random.Random(43)
        rows = tuple(generator.getrandbits(80) for _ in range(100))
        rows += (0, rows[0] ^ rows[1])
        native = filter_matrix(rows, weight_two=True, budget=allowance())
        native_solver = DependencySolver(native, budget=allowance())
        expected = native_solver.run()

        for name in BACKENDS:
            integer = arithmetic.get_backend(name).integer
            values = tuple(
                integer(v)
                for v in (1, 2**257 * 3**81, 2**257 * 13, 17**32, 45)
            )
            detector = SmoothBatch(
                (2, 3, 5, 7, 11), backend=name, budget=allowance()
            )
            self.assertEqual(detector.residuals(values), (1, 1, 13, 17**32, 1))
            matrix = filter_matrix(
                tuple(integer(row) for row in rows),
                weight_two=True,
                budget=allowance(),
            )
            solver = DependencySolver(matrix, budget=allowance())
            self.assertEqual(solver.run(), expected)
            self.assertEqual(
                _solver_digest(solver, encoding="hex-v1"),
                _solver_digest(native_solver, encoding="hex-v1"),
            )
            for mask in matrix.masks:
                self.assertEqual(arithmetic.is_mpz(mask), name == "gmpy2-mpz")

    def test_parallel_serial_and_spawned_worker_contracts(self):
        for name in BACKENDS:
            expected = None
            for mode, workers in (("serial", 1), ("process", 2)):
                config = parallel_config(backend=name)
                job = ParallelSIQSJob(
                    4513 * 7717, config=config, budget=allowance()
                )
                with CollectionPool(mode, workers) as pool:
                    result = job.run(pool=pool, fixed_work=True)
                signature = (
                    result.divisor,
                    result.stats["scanned"],
                    result.stats["work_used"],
                )
                if expected is None:
                    expected = signature
                self.assertEqual(signature, expected)
                checkpoint = json.loads(json.dumps(job.checkpoint()))
                restored = ParallelSIQSJob.from_checkpoint(
                    checkpoint, budget=allowance()
                )
                self.assertEqual(
                    arithmetic.is_mpz(restored.n), name == "gmpy2-mpz"
                )


if __name__ == "__main__":
    unittest.main()
