"""Process-local QS representation experiments; no production policy switch."""

import sys
from dataclasses import replace

from ....common import utils
from ....qs import extraction, factor_base, siqs
from ....qs.sieve_collector import SieveCollector


def install_native_small():
    """Keep mpz polynomials and native roots, offsets and matrix masks."""
    original_roots = factor_base.modular_square_roots
    original_inverse = utils.modular_inverse

    def modular_roots(value, prime, *, budget=None):
        return original_roots(int(value % prime), prime, budget=budget)

    def inversion(value, modulus):
        if modulus <= factor_base.MAX_FACTOR_BASE_BOUND:
            return original_inverse(int(value % modulus), int(modulus))
        return original_inverse(value, modulus)

    for imported in tuple(sys.modules.values()):
        if imported is None or not imported.__name__.startswith("v2"):
            continue
        if getattr(imported, "modular_square_roots", None) is original_roots:
            imported.modular_square_roots = modular_roots
        if getattr(imported, "modular_inverse", None) is original_inverse:
            imported.modular_inverse = inversion

    original_parity = extraction.parity_bits

    def parity(relation, base):
        bits = int(relation.sign < 0)
        for prime, exponent in relation.exponents:
            if prime not in base._columns:
                raise ValueError("exponent prime is outside the factor base")
            if exponent % 2:
                bits |= 1 << base._columns[prime]
        return bits

    for imported in tuple(sys.modules.values()):
        if imported is not None and imported.__name__.startswith("v2"):
            if getattr(imported, "parity_bits", None) is original_parity:
                imported.parity_bits = parity

    def native_roots(roots):
        return tuple(
            replace(root, roots=tuple(map(int, root.roots))) for root in roots
        )

    class NativeOffsetsCollector(SieveCollector):
        def __init__(self, *args, precomputed_roots=None, **kwargs):
            if precomputed_roots is not None:
                precomputed_roots = native_roots(precomputed_roots)
            super().__init__(
                *args, precomputed_roots=precomputed_roots, **kwargs
            )
            self._roots = native_roots(self._roots)

        def set_polynomial(self, polynomial, roots, **kwargs):
            super().set_polynomial(polynomial, native_roots(roots), **kwargs)

    siqs.SieveCollector = NativeOffsetsCollector
