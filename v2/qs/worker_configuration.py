"""Concrete current SIQS configuration interchange without executor imports."""

from dataclasses import asdict

from .extended_policy import BoundedYieldSieveConfig, ExtendedYieldSieveConfig
from .sieve_collector import (
    CombinedYieldSieveConfig,
    DoubleLargeSieveConfig,
    GrowingDoubleLargeSieveConfig,
    SieveConfig,
    YieldSieveConfig,
)
from .siqs import CountedSIQSConfig, SIQSConfig

COLLECTORS = {
    c.__name__: c
    for c in (
        SieveConfig,
        DoubleLargeSieveConfig,
        GrowingDoubleLargeSieveConfig,
        YieldSieveConfig,
        CombinedYieldSieveConfig,
        ExtendedYieldSieveConfig,
        BoundedYieldSieveConfig,
    )
}


def encode_config(config):
    return dict(
        kind=type(config).__name__,
        collector=type(config.collector).__name__,
        values=asdict(config),
    )


def decode_config(record):
    kinds = {c.__name__: c for c in (SIQSConfig, CountedSIQSConfig)}
    if record["kind"] not in kinds or record["collector"] not in COLLECTORS:
        raise ValueError("unknown current-worker configuration type")
    values = dict(record["values"])
    values["collector"] = COLLECTORS[record["collector"]](
        **values["collector"]
    )
    return kinds[record["kind"]](**values)
