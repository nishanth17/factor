"""Translate pre-H3 experiment recipes at the recorder boundary only."""

from ....qs.siqs import CountedSIQSConfig, SIQSConfig
from ....qs.worker_configuration import COLLECTORS
from ...infrastructure.resource_control import adapt_resource_configuration


def decode_policy(record):
    """Keep frozen numerical values while retiring recorder resource grants."""
    kinds = {c.__name__: c for c in (SIQSConfig, CountedSIQSConfig)}
    values = dict(record["values"])
    values["collector"] = adapt_resource_configuration(
        COLLECTORS[record["collector"]], **values["collector"]
    )
    return adapt_resource_configuration(kinds[record["kind"]], **values)
