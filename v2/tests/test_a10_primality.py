"""Independent proofs, strict endpoints and checked classification resume."""

import copy
import io
import random
import unittest
from contextlib import redirect_stdout
from unittest.mock import patch

from v2 import portfolio
from v2.benchmarks.primality.a10.a10_inputs import (
    REPORTED_INPUT,
    REPORTED_PRIME,
    load_corpus,
)
from v2.benchmarks.primality.a10.a10_primality import (
    load_control,
    load_protocol,
    load_v1_adapter,
)
from v2.common import arithmetic, utils
from v2.execution.budget import Budget
from v2.factor import factorize
from v2.tests.test_phase_two import reseal


def allowance(work=1000000, **kwargs):
    return Budget(work_limit=work, seconds=None, cpu_seconds=None, **kwargs)


def config(module=portfolio, **kwargs):
    return module.PortfolioConfig(
        trial_bound=5,
        rho_attempts=0,
        pm1_attempts=0,
        pm1_b1=2,
        pm1_b2=2,
        ecm_tiers=(),
        segment_size=16,
        **kwargs,
    )


def reference_sprp(n, base):
    """Independent direct congruences, avoiding production decomposition."""
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    return pow(base, d, n) == 1 or any(
        pow(base, d * 2**r, n) == n - 1 for r in range(s)
    )


class Draws:
    def __init__(self, base=2):
        self.base, self.calls = base, 0

    def randint(self, lo, hi):
        if not lo <= self.base <= hi:
            raise AssertionError("invalid test base")
        self.calls += 1
        return self.base


class A10Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = load_corpus()
        cls.control = load_control()
        load_protocol()

    def test_proof_backed_classifications_and_booleans(self):
        for fixture in self.data["fixtures"]:
            n = fixture["n"]
            with self.subTest(n=n):
                result = utils.classify_prime(n, rng=random.Random(7))
                expected = (
                    utils.Primality.PROVEN
                    if n < 3317044064679887385961981
                    else utils.Primality.PROBABLE
                )
                if not fixture["prime"]:
                    expected = utils.Primality.COMPOSITE
                self.assertEqual(result, expected)
                self.assertIs(
                    utils.is_prime(n, rng=random.Random(7)), fixture["prime"]
                )

    def test_strict_dispatch_boundaries(self):
        for n in (-1, 0, 1):
            self.assertIsNone(utils.deterministic_bases(n))
        for n in (True, 3.0, "3"):
            with self.assertRaises(TypeError):
                utils.deterministic_bases(n)
        sets = (
            (31, 73),
            (2, 7, 61),
            (2, 325, 9375, 28178, 450775, 9780504, 1795265022),
            (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37),
            (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41),
            None,
        )
        bounds = (
            9080191,
            4759123141,
            2**64,
            318665857834031151167461,
            3317044064679887385961981,
        )
        for index, bound in enumerate(bounds):
            self.assertEqual(utils.deterministic_bases(bound - 1), sets[index])
            self.assertEqual(utils.deterministic_bases(bound), sets[index + 1])
            self.assertEqual(
                utils.deterministic_bases(bound + 1), sets[index + 1]
            )

    def test_optional_gmp_proofs_and_new_policy_resume(self):
        try:
            backend = arithmetic.get_backend("gmpy2-mpz")
        except arithmetic.BackendUnavailableError:
            self.skipTest("optional gmpy2 backend unavailable")
        for fixture in self.data["fixtures"]:
            n = fixture["n"]
            self.assertEqual(
                utils.classify_prime(backend.integer(n), rng=random.Random(7)),
                utils.classify_prime(n, rng=random.Random(7)),
            )
        cfg = config(backend="gmpy2-mpz", ecm_program_bytes=65536)
        whole = portfolio.factorize_bounded(
            REPORTED_PRIME, seed=7, config=cfg, budget=allowance()
        )
        paused = portfolio.factorize_bounded(
            REPORTED_PRIME, seed=7, config=cfg, budget=allowance(200)
        )
        resumed = portfolio.factorize_bounded(
            REPORTED_PRIME,
            config=cfg,
            budget=allowance(),
            checkpoint=paused.checkpoint,
        )
        self.assertEqual(paused.checkpoint["payload"]["version"], 6)
        self.assertEqual(
            (resumed.result, resumed.work_used),
            (whole.result, whole.work_used),
        )
        self.assertTrue(resumed.result.proven)

    def test_threshold_counterexamples_with_independent_congruences(self):
        primes = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41)
        for n, count, factors in (
            (341550071728321, 7, (10670053, 32010157)),
            (318665857834031151167461, 12, (399165290221, 798330580441)),
            (3317044064679887385961981, 13, (1287836182261, 2575672364521)),
        ):
            self.assertEqual(factors[0] * factors[1], n)
            self.assertTrue(all(reference_sprp(n, a) for a in primes[:count]))
            witness = next(
                a for a in range(2, 100) if not reference_sprp(n, a)
            )
            self.assertEqual(
                utils.classify_prime(n, rng=Draws(witness)),
                utils.Primality.COMPOSITE,
            )
        # The outer endpoint can survive adversarial random bases. Such an
        # answer remains probable; finite random MR is not a certificate.
        self.assertEqual(
            utils.classify_prime(n, tolerance=3, rng=Draws()),
            utils.Primality.PROBABLE,
        )

    def test_requested_mode_and_no_rng_on_fixed_ranges(self):
        for n in (1009, REPORTED_PRIME, 2**127 - 1):
            for rounds in (1, 3, 9):
                rng = Draws()
                self.assertEqual(
                    utils.classify_prime(
                        n, use_probabilistic=True, tolerance=rounds, rng=rng
                    ),
                    utils.Primality.PROBABLE,
                )
                self.assertEqual(rng.calls, rounds)
        rng = Draws()
        self.assertEqual(
            utils.classify_prime(REPORTED_PRIME, rng=rng),
            utils.Primality.PROVEN,
        )
        self.assertEqual(rng.calls, 0)
        rng = Draws(2)
        self.assertEqual(
            utils.classify_prime(
                41**2, use_probabilistic=True, tolerance=9, rng=rng
            ),
            utils.Primality.COMPOSITE,
        )
        self.assertEqual(rng.calls, 1)
        fixed_work = []
        for rounds in (1, 3, 9):
            run = portfolio.factorize_bounded(
                REPORTED_PRIME,
                seed=7,
                config=config(primality_rounds=rounds),
                budget=allowance(),
            )
            self.assertTrue(run.result.proven)
            fixed_work.append(run.work_used)
        self.assertEqual(len(set(fixed_work)), 1)

    def test_new_resume_work_rng_and_cancellation(self):
        cfg = config()
        whole = portfolio.factorize_bounded(
            REPORTED_PRIME, seed=7, config=cfg, budget=allowance()
        )
        initial_rng = random.Random(7).getstate()
        paused = portfolio.factorize_bounded(
            REPORTED_PRIME, seed=7, config=cfg, budget=allowance(200)
        )
        payload = paused.checkpoint["payload"]
        self.assertEqual(payload["primality"], "mr13-strict-v1")
        self.assertEqual(payload["state"]["current"]["prime_job"]["index"], 1)
        self.assertEqual(portfolio._tuples(payload["rng"]), initial_rng)
        self.assertEqual(paused.reason, "work_limit")
        self.assertEqual(paused.result.remaining, (REPORTED_PRIME,))
        resumed = portfolio.factorize_bounded(
            REPORTED_PRIME,
            config=cfg,
            budget=allowance(),
            checkpoint=paused.checkpoint,
        )
        self.assertEqual(
            (resumed.result, resumed.work_used),
            (whole.result, whole.work_used),
        )
        cancelled = portfolio.factorize_bounded(
            REPORTED_PRIME,
            config=cfg,
            budget=allowance(cancelled=lambda: True),
            checkpoint=paused.checkpoint,
        )
        self.assertEqual(cancelled.reason, "cancelled")
        self.assertEqual(cancelled.work_used, paused.work_used)
        self.assertEqual(cancelled.result.reconstruct(), REPORTED_PRIME)
        self.assertEqual(
            cancelled.checkpoint["payload"]["state"],
            paused.checkpoint["payload"]["state"],
        )
        for changes in ({"seconds": 0}, {"cpu_seconds": 0}):
            budget = Budget(work_limit=1000000, **changes)
            stopped = portfolio.factorize_bounded(
                REPORTED_PRIME, config=cfg, budget=budget
            )
            self.assertEqual(stopped.result.remaining, (REPORTED_PRIME,))
            self.assertEqual(stopped.work_used, 0)

    def test_legacy_witnesses_and_terminal_labels_remain_conservative(self):
        old = self.control
        for backend, program_bytes in (
            ("python-int", 0),
            ("python-int", 65536),
            ("gmpy2-mpz", 0),
            ("gmpy2-mpz", 65536),
        ):
            if backend == "gmpy2-mpz":
                try:
                    arithmetic.get_backend(backend)
                except arithmetic.BackendUnavailableError:
                    continue
            cfg = config(
                old.portfolio,
                backend=backend,
                ecm_program_bytes=program_bytes,
            )
            new_cfg = config(backend=backend, ecm_program_bytes=program_bytes)
            kwargs = dict(seed=7, config=cfg)
            whole = old.portfolio.factorize_bounded(
                REPORTED_PRIME,
                budget=old.budget.Budget(
                    work_limit=1000000, seconds=None, cpu_seconds=None
                ),
                **kwargs,
            )
            paused = old.portfolio.factorize_bounded(
                REPORTED_PRIME,
                budget=old.budget.Budget(
                    work_limit=200, seconds=None, cpu_seconds=None
                ),
                **kwargs,
            )
            self.assertNotIn("primality", paused.checkpoint["payload"])
            resumed = portfolio.factorize_bounded(
                REPORTED_PRIME,
                config=new_cfg,
                budget=allowance(),
                checkpoint=paused.checkpoint,
            )
            self.assertEqual(
                (resumed.result.factors[0].certainty.value, resumed.work_used),
                ("probable_prime", whole.work_used),
            )
            self.assertEqual(
                resumed.checkpoint["payload"]["primality"], "mr64-strict-v1"
            )
            again = portfolio.factorize_bounded(
                REPORTED_PRIME,
                config=new_cfg,
                budget=allowance(),
                checkpoint=resumed.checkpoint,
            )
            self.assertEqual(again.result, resumed.result)

    def test_legacy_policy_covers_pending_children_and_rng(self):
        old = self.control
        options = dict(
            trial_bound=30000, rho_attempts=0, pm1_attempts=0, ecm_tiers=()
        )
        cfg = old.portfolio.PortfolioConfig(**options)
        kwargs = dict(seed=104729, config=cfg)
        whole = old.portfolio.factorize_bounded(
            REPORTED_INPUT,
            budget=old.budget.Budget(
                work_limit=1000000, seconds=None, cpu_seconds=None
            ),
            **kwargs,
        )
        paused = old.portfolio.factorize_bounded(
            REPORTED_INPUT,
            budget=old.budget.Budget(
                work_limit=200, seconds=None, cpu_seconds=None
            ),
            **kwargs,
        )
        resumed = portfolio.factorize_bounded(
            REPORTED_INPUT,
            config=portfolio.PortfolioConfig(**options),
            budget=allowance(),
            checkpoint=paused.checkpoint,
        )
        self.assertTrue(whole.result.complete)
        self.assertTrue(resumed.result.complete)
        self.assertEqual(
            [f.value for f in resumed.result.factors],
            [61, 27103, REPORTED_PRIME],
        )
        self.assertEqual(resumed.work_used, whole.work_used)
        self.assertEqual(
            resumed.checkpoint["payload"]["rng"],
            whole.checkpoint["payload"]["rng"],
        )
        self.assertEqual(
            [
                (f.value, f.exponent, f.certainty.value)
                for f in resumed.result.factors
            ],
            [
                (f.value, f.exponent, f.certainty.value)
                for f in whole.result.factors
            ],
        )
        self.assertFalse(resumed.result.proven)

    def test_primality_policies_compose_with_paired_checkpoint_schemas(self):
        for backend in ("python-int", "gmpy2-mpz"):
            if backend == "gmpy2-mpz":
                try:
                    arithmetic.get_backend(backend)
                except arithmetic.BackendUnavailableError:
                    continue
            old_cfg = config(
                self.control.portfolio,
                backend=backend,
                ecm_program_bytes=65536,
            )
            old_whole = self.control.portfolio.factorize_bounded(
                REPORTED_PRIME,
                seed=7,
                config=old_cfg,
                budget=self.control.budget.Budget(
                    work_limit=1000000, seconds=None, cpu_seconds=None
                ),
            )
            old_paused = self.control.portfolio.factorize_bounded(
                REPORTED_PRIME,
                seed=7,
                config=old_cfg,
                budget=self.control.budget.Budget(
                    work_limit=200, seconds=None, cpu_seconds=None
                ),
            )
            for version, pairing in (
                (7, {"ecm_pair_distance": 0}),
                (8, {"ecm_pair_wheel": 6}),
            ):
                with self.subTest(backend=backend, version=version):
                    cfg = config(
                        backend=backend, ecm_program_bytes=65536, **pairing
                    )
                    whole = portfolio.factorize_bounded(
                        REPORTED_PRIME,
                        seed=7,
                        config=cfg,
                        budget=allowance(),
                    )
                    paused = portfolio.factorize_bounded(
                        REPORTED_PRIME,
                        seed=7,
                        config=cfg,
                        budget=allowance(200),
                    )
                    metadata = paused.checkpoint["payload"]
                    self.assertEqual(metadata["version"], version)
                    self.assertEqual(metadata["primality"], "mr13-strict-v1")
                    resumed = portfolio.factorize_bounded(
                        REPORTED_PRIME,
                        config=cfg,
                        budget=allowance(),
                        checkpoint=paused.checkpoint,
                    )
                    self.assertTrue(resumed.result.proven)
                    self.assertEqual(
                        (resumed.result, resumed.work_used),
                        (whole.result, whole.work_used),
                    )

                    # No ECM job exists: transplant the independently frozen
                    # random-witness state into B2's configuration envelope.
                    legacy = copy.deepcopy(old_paused.checkpoint)
                    for name in ("version", "schedule", "config"):
                        legacy["payload"][name] = copy.deepcopy(metadata[name])
                    legacy = reseal(legacy)
                    restored = portfolio.factorize_bounded(
                        REPORTED_PRIME,
                        config=cfg,
                        budget=allowance(),
                        checkpoint=legacy,
                    )
                    self.assertEqual(
                        restored.result.factors[0].certainty,
                        utils.Primality.PROBABLE,
                    )
                    self.assertEqual(restored.work_used, old_whole.work_used)
                    self.assertEqual(
                        restored.checkpoint["payload"]["primality"],
                        "mr64-strict-v1",
                    )
                    self.assertEqual(
                        restored.result.reconstruct(), REPORTED_PRIME
                    )
                    bad = copy.deepcopy(legacy)
                    bad["payload"]["primality"] = "mr13-strict-v1"
                    with self.assertRaises(ValueError):
                        portfolio.factorize_bounded(
                            REPORTED_PRIME,
                            config=cfg,
                            budget=allowance(),
                            checkpoint=reseal(bad),
                        )

    def test_above_final_bound_stays_probable_across_resume(self):
        n = next(
            f["n"]
            for f in self.data["fixtures"]
            if f["prime"] and f["n"] > 3317044064679887385961981
        )
        cfg = config(primality_rounds=5)
        whole = portfolio.factorize_bounded(
            n, seed=7, config=cfg, budget=allowance()
        )
        shifts = ((n - 1) & -(n - 1)).bit_length() - 1
        initial = portfolio.factorize_bounded(
            n, seed=7, config=cfg, budget=allowance(200)
        )
        self.assertEqual(
            initial.checkpoint["payload"]["state"]["current"]["prime_job"][
                "index"
            ],
            0,
        )
        first_witness_work = initial.work_used + n.bit_length() + shifts
        paused = portfolio.factorize_bounded(
            n,
            config=cfg,
            budget=allowance(first_witness_work),
            checkpoint=initial.checkpoint,
        )
        job = paused.checkpoint["payload"]["state"]["current"]["prime_job"]
        self.assertIsNone(job["bases"])
        self.assertEqual(job["index"], 1)
        resumed = portfolio.factorize_bounded(
            n, config=cfg, budget=allowance(), checkpoint=paused.checkpoint
        )
        self.assertTrue(resumed.result.complete)
        self.assertFalse(resumed.result.proven)
        self.assertEqual(
            resumed.result.factors[0].certainty, utils.Primality.PROBABLE
        )
        self.assertEqual(
            (resumed.result, resumed.work_used),
            (whole.result, whole.work_used),
        )
        self.assertEqual(
            resumed.checkpoint["payload"]["rng"],
            whole.checkpoint["payload"]["rng"],
        )

    def test_every_fixed_witness_boundary_and_upper_range_resume(self):
        upper = next(
            f["n"]
            for f in self.data["fixtures"]
            if f["prime"]
            and 318665857834031151167461 <= f["n"] < 3317044064679887385961981
        )
        for n in (REPORTED_PRIME, upper):
            cfg = config()
            whole = portfolio.factorize_bounded(
                n, seed=7, config=cfg, budget=allowance()
            )
            paused = portfolio.factorize_bounded(
                n, seed=7, config=cfg, budget=allowance(200)
            )
            job = paused.checkpoint["payload"]["state"]["current"]["prime_job"]
            charge = n.bit_length() + job["s"]
            initial = paused.work_used - job["index"] * charge
            for index in range(len(job["bases"])):
                stopped = portfolio.factorize_bounded(
                    n,
                    seed=7,
                    config=cfg,
                    budget=allowance(initial + index * charge),
                )
                state = stopped.checkpoint["payload"]["state"]
                self.assertEqual(state["current"]["prime_job"]["index"], index)
                self.assertEqual(stopped.result.remaining, (n,))
                resumed = portfolio.factorize_bounded(
                    n,
                    config=cfg,
                    budget=allowance(),
                    checkpoint=stopped.checkpoint,
                )
                self.assertEqual(
                    (resumed.result, resumed.work_used),
                    (whole.result, whole.work_used),
                )

    def test_resealed_corrupt_policy_prefix_and_proof_are_rejected(self):
        paused = portfolio.factorize_bounded(
            REPORTED_PRIME, seed=7, config=config(), budget=allowance(200)
        )
        for change in ("policy", "base", "prefix", "index"):
            bad = copy.deepcopy(paused.checkpoint)
            payload = bad["payload"]
            job = payload["state"]["current"]["prime_job"]
            if change == "policy":
                payload["primality"] = "unknown"
            elif change == "base":
                job["bases"][-1] = 43
            elif change == "prefix":
                job["tested"][0] = 3
            else:
                job["index"] = True
            with self.assertRaises((ValueError, TypeError)):
                portfolio.factorize_bounded(
                    REPORTED_PRIME,
                    config=config(),
                    budget=allowance(),
                    checkpoint=reseal(bad),
                )

        bad = copy.deepcopy(paused.checkpoint)
        bad["payload"]["primality"] = "mr64-strict-v1"
        with self.assertRaises(ValueError):
            portfolio.factorize_bounded(
                REPORTED_PRIME,
                config=config(),
                budget=allowance(),
                checkpoint=reseal(bad),
            )

        complete = portfolio.factorize_bounded(
            REPORTED_PRIME, seed=7, config=config(), budget=allowance()
        )
        bad = copy.deepcopy(complete.checkpoint)
        bad["payload"]["primality"] = "mr64-strict-v1"
        with self.assertRaises(ValueError):
            portfolio.factorize_bounded(
                REPORTED_PRIME,
                config=config(),
                budget=allowance(),
                checkpoint=reseal(bad),
            )

    def test_composite_resume_and_prime_power_multiplicity(self):
        n = 318665857834031151167461
        whole = portfolio.factorize_bounded(
            n, seed=7, config=config(), budget=allowance()
        )
        paused = portfolio.factorize_bounded(
            n, seed=7, config=config(), budget=allowance(200)
        )
        resumed = portfolio.factorize_bounded(
            n,
            config=config(),
            budget=allowance(),
            checkpoint=paused.checkpoint,
        )
        self.assertEqual(
            (resumed.result, resumed.work_used),
            (whole.result, whole.work_used),
        )
        self.assertEqual(resumed.result.remaining, (n,))
        self.assertEqual(
            resumed.checkpoint["payload"]["state"]["classifications"][str(n)],
            "composite",
        )
        for exponent in (2, 3, 5):
            n = REPORTED_PRIME**exponent
            result = portfolio.factorize_bounded(
                n, seed=7, config=config(), budget=allowance()
            ).result
            self.assertEqual(result.reconstruct(), n)
            self.assertTrue(result.proven)
            self.assertEqual(
                [(f.value, f.exponent) for f in result.factors],
                [(REPORTED_PRIME, exponent)],
            )

    def test_reported_recursive_result_quiet_and_unresolved(self):
        output = io.StringIO()
        with redirect_stdout(output):
            direct = factorize(REPORTED_INPUT, seed=7)
            bounded = portfolio.factorize_bounded(
                REPORTED_INPUT,
                seed=7,
                config=portfolio.PortfolioConfig(
                    trial_bound=30000, ecm_tiers=()
                ),
                budget=allowance(),
            )
        self.assertEqual(output.getvalue(), "")
        for result in (direct, bounded.result):
            self.assertEqual(result.reconstruct(), REPORTED_INPUT)
            self.assertTrue(result.complete)
            self.assertEqual(
                [f.value for f in result.factors], [61, 27103, REPORTED_PRIME]
            )
            self.assertTrue(
                all(
                    f.certainty is utils.Primality.PROVEN
                    for f in result.factors
                )
            )
        n = 1000003 * REPORTED_PRIME
        limited = portfolio.factorize_bounded(
            n, seed=7, config=config(), budget=allowance()
        )
        self.assertEqual(limited.result.remaining, (n,))
        self.assertEqual(limited.result.reconstruct(), n)
        with patch("v2.common.utils.resolve_rng", side_effect=AssertionError):
            self.assertEqual(
                factorize(REPORTED_PRIME).factors[0].certainty,
                utils.Primality.PROVEN,
            )

    def test_v1_adapter_is_validated_only_in_its_supported_subset(self):
        adapter = load_v1_adapter()
        self.assertEqual(adapter.gcd(18, 24), 6)
        for fixture in self.data["fixtures"]:
            if 2 <= fixture["n"] < 3317044064679887385961981:
                self.assertEqual(
                    bool(adapter.is_prime(fixture["n"])), fixture["prime"]
                )


if __name__ == "__main__":
    unittest.main()
