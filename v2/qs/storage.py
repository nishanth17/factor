"""Optional storage ceilings for allocator-grown relation containers."""

import sys

from ..common import utils


def storage_limit(value):
    """Resolve None to addressable size without reserving that storage.

    Python containers grow from actual retained data. The integer sentinel
    lets existing exact reservation arithmetic keep reporting owned bytes
    while an absent user quota leaves allocation to the runtime.
    """
    if value is None:
        return sys.maxsize
    return utils.require_integer(value, "storage limit", 0)


class RelationPublicationError(RuntimeError):
    """Allocation interrupted publication; restart a prior checked snapshot."""


def require_published(collector):
    """Never reuse or serialize a store interrupted during publication."""
    if getattr(collector, "_publication_failed", False):
        raise RelationPublicationError(
            "relation publication was interrupted; restore a prior checkpoint"
        )
