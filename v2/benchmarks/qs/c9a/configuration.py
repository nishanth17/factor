"""Independent full-arm policies with input-resolved numerical admission."""

from dataclasses import asdict, replace
from math import isqrt

from v2.benchmarks.infrastructure.resource_control import (
    adapt_resource_configuration,
)

from ....qs import (
    CountedSIQSConfig,
    GrowingDoubleLargeSieveConfig,
    SIQSConfig,
    YieldSieveConfig,
)
from .control import control_modules

BASE_COMMIT = "b0813230a05396874c41c649fa01738484a1f55e"
IDENTITY = "c9a-full-yield-v2"
SIZES = (50, 55, 60, 62, 65, 70, 75, 80)
SEEDS = (104729, 130363)
HELPER_BYTES = 32 * 2**20
COUNTS = {55: 2500, 60: 5000, 62: 6000, 65: 8500, 70: 10500}
FACTORS = {62: 9, 65: 9, 70: 10, 75: 11, 80: 12}
ROWS = {62: 16000, 65: 24000, 70: 32000, 75: 48000, 80: 60000}
ATOMS = {62: 98304, 65: 131072, 70: 196608, 75: 262144, 80: 327680}
PARTIALS = {key: value // 2 for key, value in ATOMS.items()}
PREFIX_CAPS = dict(zip(SIZES, (20, 30, 45, 60, 90, 120, 180, 300)))
COMPLETE_CAPS = dict(
    zip(SIZES[:-1], (180, 600, 1800, 2400, 3600, 7200, 14400))
)
# Historical validated 50/55/60 calls peaked at 3.001e9 charged units/s.
# Round up to 4e9 and allow eightfold headroom for a different operation mix.
# This is a finite guard, not a prediction of useful factoring throughput.
WORK_RATE_GUARD = 32 * 10**9
WORK_CAPS = {
    d: WORK_RATE_GUARD * COMPLETE_CAPS.get(d, PREFIX_CAPS[d]) for d in SIZES
}
ARMS = {
    "C0": (),
    "C1": ("P", "I", "E", "G", "D", "S"),
    "wide": ("P", "I", "E", "G", "D", "S"),
    "noE": ("P", "I", "G", "D", "S"),
    "support": ("P", "G", "D", "S"),
    "noP": ("I", "E", "G", "D", "S"),
    "noG": ("P", "I", "E", "D", "S"),
    "noD": ("P", "I", "E", "G", "S"),
    "noS": ("P", "I", "E", "G", "D"),
    "P": ("P",),
    "I": ("I",),
    "E": ("E",),
    "G": ("G",),
    "D": ("D",),
    "S": ("S",),
    "off": (),
}


def owned_bytes(digits):
    return (1 if digits <= 60 else 2 if digits <= 65 else 4) * 2**30


def rss_bytes(digits):
    return (3 if digits <= 60 else 5 if digits <= 65 else 8) * 2**30


def admission_envelope(base_bound, *, multiple=16):
    """A base-relative screening envelope, not a measured optimum.

    Preserve the old domain. The smaller cap bisects its enlarged logarithmic
    endpoint interval; the 16/64 alternatives test density versus extra intake.
    All decisions below use exact integers, including the fourth root.
    """
    if type(base_bound) is not int or not 3 <= base_bound <= 1000000:
        raise ValueError("admission needs the resolved factor-base bound")
    if multiple not in (16, 64):
        raise ValueError("unfrozen intake multiple")
    old, large = 115200000000, 3000000
    return dict(
        candidate_bound=old,
        large_product_bound=min(large**2, max(old, multiple * base_bound**2)),
        small_prime_bound=min(
            large, max(isqrt(old), isqrt(isqrt(multiple * base_bound**4)))
        ),
    )


def configuration(digits, arm="C0", *, admission_base_bound=None):
    """No routing changes; larger bands are explicit experimental policies."""
    if digits not in SIZES + (40, 51, 52, 56, 59, 61) or arm not in ARMS:
        raise ValueError("unfrozen C9a size or arm")
    owned = owned_bytes(digits)
    if digits <= 60 or digits == 61:
        # The exact 1 GiB factory input selects marking, plan and pacing.
        # Subtracting helpers before this call would silently change C0.
        # Historical v2/v3/v4 cells bind the immutable factory. A current
        # mainline comparison supplies explicit incumbents to the separate
        # integration builder; it must not silently rewrite these controls.
        from .... import qs

        control, _ = control_modules()
        config = decode(
            asdict(
                adapt_resource_configuration(
                    control.default_siqs_config,
                    10 ** (digits - 1),
                    recipe="calibrated"
                    if owned in (None, 2**30)
                    else "reference",
                )
            ),
            qs,
        )
    else:
        values = dict(
            base_bound=250000
            if digits <= 70
            else 400000
            if digits == 75
            else 600000,
            multiplier=0,
            half_width=65536,
            max_half_width=65536,
            factor_count=FACTORS[digits],
            family_count=100000,
            pool_size=64,
            assignment_policy="flyer",
            polynomials_per_family=64,
            max_stalled=100000,
            max_trivial=100000,
            row_excess=32,
            batch_width=4096,
            filter_row_growth=32,
            sign_policy="swap_high",
            solve_policy="collection",
            memory_bytes=owned,
            checkpoint_bytes=1048576,
            collector=GrowingDoubleLargeSieveConfig(
                sparse_queue_bytes=32 * 2**20,
                marking="sparse",
                block_width=4096,
                division="bucket",
                score_policy="powers",
                power_plan_bytes=2 * 2**20,
                residual_bound=3000000,
                large_prime_bound=3000000,
                large_product_bound=115200000000,
                max_partials=PARTIALS[digits],
                split_call_limit=1048576,
                preparation_cache_bytes=2 * 2**20
                if ROWS[digits] is None
                else min(2 * 2**20, 2048 * ROWS[digits]),
                target_rows=ROWS[digits],
            ),
        )
        # The existing count resolver caps its sieve at 250000. Frontier
        # bands instead declare literal bounds; no hidden ceiling extension.
        config = (
            adapt_resource_configuration(
                CountedSIQSConfig, base_count=COUNTS[digits], **values
            )
            if digits <= 70
            else adapt_resource_configuration(SIQSConfig, **values)
        )

    config = replace(config)
    if digits not in SIZES:
        return config

    flags = ARMS[arm]
    collector = config.collector
    envelope = (
        admission_envelope(
            admission_base_bound, multiple=64 if arm == "wide" else 16
        )
        if ("I" in flags or "E" in flags)
        else {}
    )
    if "I" in flags:
        collector = replace(
            collector,
            candidate_bound=envelope["candidate_bound"],
            large_product_bound=envelope["large_product_bound"],
        )
    if set(flags) - {"I"}:
        collector = adapt_resource_configuration(
            YieldSieveConfig,
            **asdict(collector),
            small_prime_bound=envelope["small_prime_bound"]
            if "E" in flags
            else None,
            preparation_policy="target-v1" if "P" in flags else "fixed",
            cycle_policy="short-v1" if "G" in flags else "forest",
            duplicate_policy="pair-v1" if "D" in flags else "identity",
            record_policy="packed-v1" if "S" in flags else "objects",
        )
    return replace(config, collector=collector)


def decode(values, qs):
    """Decode the same literal mapping in the candidate or immutable module."""
    values = dict(values)
    collector = values["collector"]
    name = (
        "BoundedYieldSieveConfig"
        if collector.get("extension_policy") == "extended-v2"
        else "ExtendedYieldSieveConfig"
        if collector.get("extension_policy") == "extended-v1"
        else "CombinedYieldSieveConfig"
        if collector.get("yield_policy") == "yield-v2"
        else "YieldSieveConfig"
        if "yield_policy" in collector
        else "GrowingDoubleLargeSieveConfig"
        if "storage_policy" in collector
        else "DoubleLargeSieveConfig"
        if "large_prime_bound" in collector
        else "SieveConfig"
    )
    values["collector"] = adapt_resource_configuration(
        getattr(qs, name), **collector
    )
    cls = qs.CountedSIQSConfig if "base_count" in values else qs.SIQSConfig
    for key in (
        "polynomial_policy",
        "multiplier_policy",
        "root_update_policy",
    ):
        if key in values and key not in cls.__dataclass_fields__:
            if values[key] != "standard":
                raise ValueError(
                    "target runtime cannot implement coupled policy " + key
                )
            values.pop(key)
    return cls(**values)
