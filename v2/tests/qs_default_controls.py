"""Explicit historical controls for tests of pre-C10a-01 preset parameters."""

from dataclasses import replace

from v2.qs import CountedSIQSConfig, default_siqs_config
from v2.qs.collector_kernels import combine, components
from v2.qs.defaults import uses_siqs_defaults


def previous_siqs_config(n, **options):
    """Remove only the documented promotion delta, preserving old presets."""
    config = default_siqs_config(n, **options)
    parts = tuple(
        p for p in components(config.collector.marking) if p != "parity_two"
    )
    return replace(
        config,
        polynomial_policy="standard",
        multiplier_policy="standard",
        root_update_policy="standard",
        assignment_policy="flyer"
        if uses_siqs_defaults(n, backend=config.backend)
        else "reference",
        multiplier=0 if isinstance(config, CountedSIQSConfig) else 1,
        collector=replace(config.collector, marking=combine(parts)),
    )
