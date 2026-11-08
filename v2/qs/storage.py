"""Optional storage ceilings for allocator-grown relation containers."""

from ..common import utils
from ..execution.resources import AllocationFailureError


def storage_limit(value):
    """Validate an optional retained policy capacity without a sentinel."""
    if value is not None:
        return utils.require_integer(value, "capacity", 0)
    return None


class RelationPublicationError(AllocationFailureError):
    """Allocation interrupted publication; restart a prior checked snapshot."""


def require_published(collector):
    """Never reuse or serialize a store interrupted during publication."""
    if getattr(collector, "_publication_failed", False):
        raise RelationPublicationError(
            "relation publication was interrupted; restore a prior checkpoint"
        )
