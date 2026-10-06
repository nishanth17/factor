"""Optional, isolated GMP bindings for the P4.1 four-arm experiment.

Reuse the exact production function code with private globals: point values
remain mpz throughout setup, both stages and recovery. Only integer checks,
GCD and inversion change. Neither module globals nor production defaults are
patched. This adapter does not replace P4.3's production backend contract.
"""

import types

from .. import ecm, prac, utils


def _bind(module, replacements):
    namespace = vars(module).copy()
    namespace.update(replacements)
    for name, value in vars(module).items():
        if (
            name not in replacements
            and isinstance(value, types.FunctionType)
            and value.__module__ == module.__name__
        ):
            clone = types.FunctionType(
                value.__code__,
                namespace,
                name,
                value.__defaults__,
                value.__closure__,
            )
            clone.__kwdefaults__ = value.__kwdefaults__
            namespace[name] = clone
    return types.SimpleNamespace(**namespace)


def load():
    """Fail explicitly if this interpreter lacks the optional dependency."""
    import gmpy2

    def require_integer(value, name="n", minimum=None):
        if isinstance(value, bool) or not isinstance(value, (int, gmpy2.mpz)):
            raise TypeError(f"{name} must be an integer")
        if minimum is not None and value < minimum:
            raise ValueError(f"{name} must be at least {minimum}")
        return value

    def valid_divisor(divisor, n):
        return (
            isinstance(divisor, (int, gmpy2.mpz))
            and not isinstance(divisor, bool)
            and 1 < divisor < n
            and n % divisor == 0
        )

    def modular_inverse(value, modulus):
        require_integer(value, "value")
        require_integer(modulus, "modulus", 2)
        return gmpy2.invert(value, modulus)

    gmp_utils = _bind(
        utils,
        {
            "gcd": gmpy2.gcd,
            "require_integer": require_integer,
            "valid_divisor": valid_divisor,
            "modular_inverse": modular_inverse,
        },
    )
    gmp_prac = _bind(prac, {"gcd": gmpy2.gcd, "utils": gmp_utils})
    gmp_ecm = _bind(
        ecm,
        {
            "gcd": gmpy2.gcd,
            "utils": gmp_utils,
            "prac": gmp_prac,
        },
    )
    return types.SimpleNamespace(
        ecm=gmp_ecm,
        prac=gmp_prac,
        gcd=gmpy2.gcd,
        integer=gmpy2.mpz,
        identity=f"gmpy2 {gmpy2.version()} / {gmpy2.mp_version()}",
    )
