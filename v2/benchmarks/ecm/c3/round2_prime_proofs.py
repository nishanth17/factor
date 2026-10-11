"""Bounded uniform candidate sampling with independently checked proofs."""

import math

from ....common import utils
from ....execution.budget import Budget
from ....portfolio import PortfolioConfig, factorize_bounded
from .round2_corpus import GenerationRandom

MAX_NODES = 10_000
MAX_BITS = 256


def verify_proofs(certificates):
    """Check trial, Pocklington and full n-1 proofs with integers."""
    if len(certificates) > MAX_NODES:
        raise ValueError("prime certificate storage allowance exceeded")
    verified = set()

    def verify(n):
        if n in verified:
            return
        if type(n) is not int or n < 2 or n.bit_length() > MAX_BITS:
            raise ValueError("invalid certificate integer")
        proof = certificates[str(n)]
        kind = proof["kind"]
        if kind == "trial":
            if n >= 2**16 or any(
                n % d == 0 for d in range(2, math.isqrt(n) + 1)
            ):
                raise ValueError("invalid trial proof")
        elif kind in ("pocklington", "lucas"):
            witness = proof["witness"]
            if type(witness) is not int or not 1 < witness < n:
                raise ValueError("invalid prime witness")
            if kind == "pocklington":
                q = proof["q"]
                factors = [(q, 1)]
                if type(q) is not int or not 2 <= q < n:
                    raise ValueError("invalid Pocklington dependency")
                if (n - 1) % q or q * q <= n:
                    raise ValueError("insufficient Pocklington factor")
            else:
                factors = proof["factors"]
                if not factors or len(factors) > MAX_BITS:
                    raise ValueError("invalid full n-1 factorization")
                seen, product = set(), 1
                for q, exponent in factors:
                    if (
                        type(q) is not int
                        or not 2 <= q < n
                        or q in seen
                        or type(exponent) is not int
                        or not 1 <= exponent <= MAX_BITS
                    ):
                        raise ValueError("invalid full n-1 dependency")
                    seen.add(q)
                    product *= q**exponent
                if product != n - 1:
                    raise ValueError("full n-1 factors do not reconstruct")
            for q, _ in factors:
                verify(q)
            if pow(witness, n - 1, n) != 1 or any(
                math.gcd(pow(witness, (n - 1) // q, n) - 1, n) != 1
                for q, _ in factors
            ):
                raise ValueError("invalid n-1 prime proof")
        else:
            raise ValueError("unknown prime certificate kind")
        verified.add(n)

    for key in certificates:
        if str(int(key)) != key:
            raise ValueError("noncanonical certificate key")
        verify(int(key))


class ProofRandom(GenerationRandom):
    """Bound candidate draws and sampling by the shared deadline."""

    def __init__(self, seed, budget):
        self.budget = budget
        super().__init__(seed)

    def getrandbits(self, count):
        self.budget.consume(0)
        return super().getrandbits(count)


class UniformPrimeSource:
    """Sample odd integers in exact decimal bands, then prove survivors.

    Candidate filtering never establishes primality. Full n-1 factorization
    uses existing bounded engines, but every terminal factor is recursively
    proved and all certificates have an independent integer-only checker.
    Failure to finish a proof aborts generation instead of skipping a hard
    candidate and silently preferring easy n-1 structure.
    """

    def __init__(self, seed, *, budget=None):
        self.budget = (
            budget
            if budget is not None
            else Budget(
                work_limit=150_000_000_000, seconds=180, cpu_seconds=180
            )
        )
        self.generator = ProofRandom(seed, self.budget)
        self.certificates = {}
        self.config = PortfolioConfig(
            ecm_tiers=(),
            rho_attempts=8,
            rho_evaluations=2_000_000,
            trace_limit=0,
        )

    def _record(self, n, proof):
        # Recursive child proofs can fill the final slot before this parent
        # returns, so recheck storage at the actual insertion boundary.
        if len(self.certificates) >= MAX_NODES:
            raise RuntimeError("prime proof node allowance exhausted")
        self.budget.consume(1)
        self.certificates[str(n)] = proof

    def _factor_predecessor(self, n):
        self.budget.consume(0)
        remaining_work = self.budget.work_limit - self.budget.used
        remaining_wall = (
            None
            if self.budget.seconds is None
            else max(0, self.budget.seconds - self.budget.wall_used)
        )
        remaining_cpu = (
            None
            if self.budget.cpu_seconds is None
            else max(0, self.budget.cpu_seconds - self.budget.cpu_used)
        )
        child = Budget(
            work_limit=remaining_work,
            seconds=remaining_wall,
            cpu_seconds=remaining_cpu,
            cancelled=self.budget.cancelled,
        )
        try:
            # Fresh portfolio calls require an unused ledger. Give this call
            # only the parent's remaining grant; the parent clocks keep
            # running across every proof and all intervening verification.
            return factorize_bounded(
                n - 1, seed=n, config=self.config, budget=child
            )
        finally:
            # This is completed work, even if a deadline was reached during
            # recovery. Record it before checking the parent deadline so a
            # failed proof cannot lose charges or reset its allowance.
            self.budget.used += child.used
            self.budget.consume(0)

    def prove(self, n):
        """Prove one candidate or stop the entire generation attempt."""
        self.budget.consume(0)
        if type(n) is not int or not 2 <= n < 10**20:
            raise ValueError("uniform proof input exceeds 20 digits")
        if str(n) in self.certificates:
            return
        if len(self.certificates) >= MAX_NODES:
            raise RuntimeError("prime proof node allowance exhausted")
        if n < 2**16:
            self.budget.consume(math.isqrt(n))
            if any(n % d == 0 for d in range(2, math.isqrt(n) + 1)):
                raise ValueError("composite proof candidate")
            self._record(n, dict(kind="trial"))
            return

        result = self._factor_predecessor(n).result
        if not result.complete or result.reconstruct() != n - 1:
            raise RuntimeError("full n-1 proof search did not complete")
        factors = [(int(f.value), f.exponent) for f in result.factors]
        if math.prod(q**e for q, e in factors) != n - 1:
            raise ValueError("proof factors do not reconstruct n-1")
        for q, _ in factors:
            self.prove(q)

        # A common witness gives order n-1 modulo n. This proves primality
        # independently of the factorer's probable/proven terminal labels.
        for witness in range(2, min(n, 258)):
            self.budget.consume((len(factors) + 1) * n.bit_length() ** 2)
            if pow(witness, n - 1, n) == 1 and all(
                math.gcd(pow(witness, (n - 1) // q, n) - 1, n) == 1
                for q, _ in factors
            ):
                self._record(
                    n, dict(kind="lucas", factors=factors, witness=witness)
                )
                return
        raise RuntimeError("finite prime witness search did not finish")

    def prime(self, digits):
        """Reject composites; never discard an unfinished proof."""
        if type(digits) is not int or not 1 <= digits <= 20:
            raise ValueError("uniform decimal band must be 1..20")
        first, upper = (10 ** (digits - 1)) | 1, 10**digits
        count = (upper - first + 1) // 2
        while True:
            candidate = first + 2 * self.generator.randrange(count)
            if not utils.is_prime(candidate, rng=self.generator):
                continue
            self.prove(candidate)
            return candidate
