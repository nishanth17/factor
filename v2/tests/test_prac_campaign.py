"""Full-campaign evidence contracts and isolated optional GMP arithmetic."""

import unittest
from dataclasses import asdict
from unittest.mock import patch

from v2.benchmarks.ecm.p41 import p41_campaign as campaign
from v2.benchmarks.support.prac_oracle import (
    affine_multiply,
    historical_points,
    matches,
)
from v2.common import utils
from v2.ecm import core as ecm
from v2.ecm import prac


class CampaignTests(unittest.TestCase):
    def test_certified_corpus_sizes_and_schedule_caps(self):
        data = campaign.load_corpus()
        self.assertEqual(
            {c["digits"] for c in data["fixtures"]}, {40, 50, 60, 70, 80}
        )
        self.assertEqual(len(data["seeds"]), 9)
        with self.assertRaises(ValueError):
            campaign.build_program(11001)
        program = campaign.build_program(11000)
        self.assertLessEqual(2 * len(program), campaign.MAX_PROGRAM_RECORDS)
        for prime, power, record, unit in program:
            self.assertEqual(record.scalar, power)
            self.assertEqual(unit.scalar, prime)
            self.assertTrue(prac.verify_chain(record))

    def test_recovery_propagates_factors_and_stops_on_saturation(self):
        program = campaign.build_program(8)[:1]
        extra = {"prime_power_replays": 0, "prime_units_replayed": 0}
        with patch.object(
            campaign,
            "apply_record",
            side_effect=[
                (1, 0),
                prac.NonunitPointError(7),
            ],
        ) as multiply:
            self.assertEqual(
                campaign.candidate_stage_one(
                    (2, 1),
                    77,
                    3,
                    program,
                    extra,
                ),
                (None, 7),
            )
            self.assertEqual(multiply.call_count, 2)
            self.assertEqual(extra["prime_units_replayed"], 1)
        with patch.object(campaign, "apply_record", return_value=(1, 0)):
            self.assertEqual(
                campaign.candidate_stage_one(
                    (2, 1),
                    77,
                    3,
                    program,
                    extra,
                ),
                (None, None),
            )

    def test_watchdog_retains_unresolved_cofactor(self):
        case = {"n": 77, "factors": [7, 11]}
        with patch.object(
            ecm, "factorize_ecm", side_effect=campaign.CampaignTimeoutError
        ):
            result = campaign.attempt("production_ladder", 77, 7, "current", 1)
        self.assertTrue(result["timed_out"])
        campaign.validate_result(result, case)
        for bad in (dict(result, unresolved=7), dict(result, factor=77)):
            with self.assertRaises(AssertionError):
                campaign.validate_result(bad, case)


class GmpCampaignTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            cls.backend = campaign.gmp_backend()
        except ImportError:
            raise unittest.SkipTest("optional PyPy gmpy2 is unavailable")

    def test_private_bindings_keep_native_points_and_original_modules(self):
        backend = self.backend
        self.assertIs(backend.ecm.point_add.__code__, ecm.point_add.__code__)
        self.assertIsNot(backend.ecm.utils, utils)
        self.assertIs(ecm.utils, utils)
        setup = backend.ecm.setup_curve(backend.integer(1009 * 1013), 11)
        self.assertIsNone(setup.factor)
        self.assertTrue(all(type(x) is backend.integer for x in setup.point))
        point = backend.ecm.scalar_multiply(
            31, *setup.point, backend.integer(1009 * 1013), setup.a24
        )
        self.assertTrue(all(type(x) is backend.integer for x in point))

    def test_gmp_checked_prac_against_independent_affine_oracle(self):
        backend = self.backend
        for prime, curve_a, point in historical_points():
            a24 = backend.integer((curve_a + 2) * pow(4, -1, prime) % prime)
            for scalar in range(129):
                result = backend.ecm.multiply_prac(
                    scalar,
                    backend.integer(point[0]),
                    backend.integer(1),
                    backend.integer(prime),
                    a24,
                )
                self.assertTrue(
                    matches(
                        result,
                        affine_multiply(scalar, point, prime, curve_a),
                        prime,
                    )
                )

    def test_complete_campaign_backend_equivalence(self):
        backend = self.backend
        for n in (1009 * 1013, 10007 * 10009, 1000003 * 1000033):
            for seed in (7, 41001):
                ordinary, native = ecm.EcmStats(), ecm.EcmStats()
                options = dict(
                    b1=31,
                    b2=401,
                    max_curves=4,
                    seed=seed,
                    _known_composite=True,
                )
                left = ecm.factorize_ecm(n, stats=ordinary, **options)
                right = backend.ecm.factorize_ecm(
                    backend.integer(n), stats=native, **options
                )
                self.assertEqual(left, right)
                self.assertEqual(asdict(ordinary), asdict(native))
                outputs = []
                for chosen in (campaign.PYTHON_BACKEND, backend):
                    stats = ecm.EcmStats()
                    extra = {
                        "prime_power_replays": 0,
                        "prime_units_replayed": 0,
                    }
                    factor = campaign.candidate(
                        chosen.integer(n),
                        b1=31,
                        b2=401,
                        seed=seed,
                        max_curves=4,
                        stats=stats,
                        extra=extra,
                        backend=chosen,
                    )
                    outputs.append((factor, asdict(stats), extra))
                self.assertEqual(outputs[0], outputs[1])


if __name__ == "__main__":
    unittest.main()
