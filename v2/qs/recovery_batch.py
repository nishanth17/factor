"""Dense-candidate smooth-part recovery and optional one-worker hints."""

from math import gcd

from .smooth_batch import product_tree


class BaseRecovery:
    """Reuse a base radical; recover valuations in the ordinary divider."""

    def __init__(self, base, memory_bytes, budget):
        bits = sum(p.bit_length() for p in base.primes)
        if bits > 2**21 or 65536 + 16 * bits > memory_bytes:
            raise MemoryError("factor-base product exceeds batch scratch")
        product = 1
        for prime in base.primes:
            budget.consume(product.bit_length() + prime.bit_length())
            product *= prime
        self.radical, self.memory_bytes = product, memory_bytes

    def residuals(self, values, budget):
        if len(values) > 256 or any(v < 1 for v in values):
            raise ValueError("invalid factor-base recovery batch")
        tree = product_tree(
            values,
            budget=budget,
            max_bits=2**20,
            memory_bytes=self.memory_bytes,
        )
        if not tree:
            return ()
        budget.consume(self.radical.bit_length() + tree[-1][0].bit_length())
        remainders = (self.radical % tree[-1][0],)
        for level in reversed(tree[:-1]):
            output = []
            for index, modulus in enumerate(level):
                parent = remainders[index // 2]
                budget.consume(parent.bit_length() + modulus.bit_length())
                output.append(parent % modulus)
            remainders = tuple(output)
        result = []
        for value, remainder in zip(values, remainders):
            for _ in range((value.bit_length() - 1).bit_length()):
                budget.consume(2 * value.bit_length() ** 2)
                remainder = remainder * remainder % value
            budget.consume(value.bit_length() ** 2)
            result.append(value // gcd(value, remainder))
        return tuple(result)


class PendingHints:
    """One finite background operation with deterministic ordered publication.

    The parent reserves the entire operation before submit. The worker has a
    separate bounded allowance; it cannot mutate the collector or its budget.
    Closing joins the one finite task, including on interruption. This is an
    opt-in overlap mechanism, not a native-Python parallel speedup claim.
    """

    def __init__(self, residuals, budget):
        from concurrent.futures import ThreadPoolExecutor

        from ..execution.budget import Budget
        from .cofactor_policy import batch_divisors

        count = len(residuals)
        bits = sum(r.bit_length() for r in residuals)
        reserve = 16 * (bits + 1) ** 2 * (count + 1)
        budget.consume(reserve)
        self.executor = ThreadPoolExecutor(
            max_workers=1, thread_name_prefix="c9a-hints"
        )
        allowance = Budget(work_limit=reserve, seconds=None, cpu_seconds=None)
        try:
            self.future = self.executor.submit(
                batch_divisors, residuals, allowance
            )
        except BaseException:
            self.executor.shutdown(wait=True)
            raise

    def result(self):
        try:
            return self.future.result()
        finally:
            self.executor.shutdown(wait=True)
