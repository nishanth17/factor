"""Bounded collector challengers; no dispatcher or default changes."""

import importlib
import inspect
import textwrap
from functools import lru_cache


@lru_cache(maxsize=16)
def experiment_collector(module, variant):
    """Build a collector in its runtime's own relation/base type universe."""
    smooth = importlib.import_module(module.__package__ + ".smooth_batch")

    class ExperimentCollector(module.SieveCollector):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self._batch_allowed = set()
            self._batch_values = {}
            if variant == "batch":
                # Reserve the complete simultaneous tree and leaf storage,
                # including construction, rather than subtracting base bytes.
                reserve = 8 * 2**20 + 512 * self.config.block_width
                if self._workspace + reserve > self.config.memory_bytes:
                    raise MemoryError("batch coexistence exceeds memory cap")
                self._workspace += reserve
                self._smooth = smooth.SmoothBatch(
                    tuple(entry.prime for entry in self.factor_base.entries),
                    budget=self.budget,
                    memory_bytes=8 * 2**20,
                )

        def _power_marks(self, index, maximum, lo, hi, stats):
            if variant != "tiny" or self._roots[index].prime > 3:
                return super()._power_marks(index, maximum, lo, hi, stats)
            roots, weight = self._roots[index], self._logs[index]

            def direct_marks():
                prime = roots.prime
                # The zero-weight base mark preserves bucket hit metadata.
                yield (
                    1 if roots.all_positions else prime,
                    (0,) if roots.all_positions else roots.roots,
                    0,
                )
                residues = (0,) if roots.all_positions else roots.roots
                step = 1 if roots.all_positions else prime
                modulus = hi - lo + 1

                for root in residues:
                    for position in range(lo + (root - lo) % step, hi, step):
                        self.budget.consume(
                            self.polynomial.n_prime.bit_length() + 1
                        )
                        value = abs(self.polynomial.value(position))
                        exponent = 0
                        if not value:
                            exponent = maximum.bit_length()
                        while value and value % prime == 0:
                            self.budget.consume(value.bit_length() + 1)
                            exponent += 1
                            value //= prime
                        stats["tiny_evaluations"] = (
                            stats.get("tiny_evaluations", 0) + 1
                        )
                        yield modulus, (position % modulus,), exponent * weight

            return direct_marks()

        def _sieve(self, lo, hi, stats):
            threshold = super()._sieve(lo, hi, stats)
            if variant != "batch":
                return threshold
            self._batch_allowed.clear()
            self._batch_values.clear()
            leaves, offsets = [], []

            def flush():
                residuals = self._smooth.residuals(tuple(leaves))
                stats["tree_calls"] = stats.get("tree_calls", 0) + 1
                stats["tree_leaves"] = stats.get("tree_leaves", 0) + len(
                    leaves
                )
                for offset, residual in zip(offsets, residuals):
                    self._batch_values[offset] = residual
                    if residual <= self.config.residual_bound:
                        self._batch_allowed.add(offset)
                leaves.clear()
                offsets.clear()

            for offset in range(hi - lo):
                if self._scores[offset] < threshold:
                    continue
                self.budget.consume(self.polynomial.n_prime.bit_length() + 1)
                value = self.polynomial.value(lo + offset)
                if not value:
                    self._batch_allowed.add(offset)
                    continue
                if not self._candidate_passes(value, offset):
                    continue
                leaves.append(abs(value))
                offsets.append(offset)
                if len(leaves) == 64:
                    flush()

            if leaves:
                flush()
            return threshold

        def _divide(self, position, offset, stats):
            if variant == "batch" and offset not in self._batch_allowed:
                self.budget.consume()
                stats["candidates"] += 1
                stats["batch_rejections"] = (
                    stats.get("batch_rejections", 0) + 1
                )
                return None, None

            atom, divisor = super()._divide(position, offset, stats)
            if variant == "batch" and atom is not None:
                if atom.residual != self._batch_values[offset]:
                    raise AssertionError("batch/scalar residual mismatch")
            return atom, divisor

    if variant == "tiny":
        source = textwrap.dedent(
            inspect.getsource(module.SieveCollector._sieve)
        )
        needle = "if self.config.power_plan_bytes"
        if source.count(needle) != 1:
            raise ValueError("tiny experiment requires the frozen mark loop")
        namespace = {}
        exec(
            compile(source.replace(needle, "if True"), __file__, "exec"),
            module.__dict__,
            namespace,
        )
        ExperimentCollector._sieve = namespace["_sieve"]

    if variant == "chunks":
        source = textwrap.dedent(
            inspect.getsource(module.SieveCollector._sieve)
        )
        # These literal fragments select the charged mark loop exactly;
        # preserve their spelling and spacing when formatting the runner.
        old = """for root in lifted:
                        hits = range((root - lo) % modulus, width, modulus)
                        self.budget.consume(len(hits) + 1)"""
        new = """for root_index, root in enumerate(lifted):
                        hits = range((root - lo) % modulus, width, modulus)
                        if root_index % 4 == 0:
                            amount = sum(
                                1 + len(range(
                                    (r - lo) % modulus, width, modulus))
                                for r in lifted[root_index:root_index + 4]
                            )
                            self.budget.consume(amount)"""
        if source.count(old) != 1:
            raise ValueError("chunk experiment requires the frozen mark loop")
        namespace = {}
        exec(
            compile(source.replace(old, new), __file__, "exec"),
            module.__dict__,
            namespace,
        )
        ExperimentCollector._sieve = namespace["_sieve"]

    return ExperimentCollector
