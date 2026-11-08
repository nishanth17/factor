"""Bounded FIFO ownership for opaque atoms outside the incidence history."""

from collections import OrderedDict

from ..common import utils
from .residual_kernel import opaque_key
from .storage import RelationPublicationError


def context_bits(atom, base):
    poly = atom.polynomial
    return (
        base.n_prime.bit_length()
        + poly.a.bit_length()
        + abs(poly.b).bit_length()
        + abs(atom.position).bit_length()
        + poly.square_coefficient.bit_length()
    )


def entry_bytes(atom, base):
    # Includes a private polynomial/payload, FIFO keys and transient index
    # pairs; unlike published records these factors have not entered slabs.
    return 4096 + 512 * len(atom.exponents) + 32 * context_bits(atom, base)


class OpaqueStaging:
    def __init__(self, config, base):
        self.config, self.base = config, base
        self.atoms = OrderedDict()
        self.used = 0

    def stage(self, atom, budget, stats, collector):
        cost = entry_bytes(atom, self.base)
        if (
            not self.config.opaque_cache_entries
            or cost > self.config.opaque_cache_bytes
        ):
            stats["opaque_staging_rejections"] = (
                stats.get("opaque_staging_rejections", 0) + 1
            )
            return
        budget.consume(len(atom.exponents) + 1)
        if atom.opaque in self.atoms:
            raise ValueError("opaque staging requires a new residual")
        used, count, retired = self.used + cost, len(self.atoms) + 1, []
        for residual, old in self.atoms.items():
            if count <= self.config.opaque_cache_entries and (
                used <= self.config.opaque_cache_bytes
            ):
                break
            budget.consume(1)
            used -= entry_bytes(old, self.base)
            count -= 1
            retired.append(residual)
        # Planning/refusal leaves the FIFO unchanged. Publication performs
        # only the insertion and actual evictions, without copying the map.
        # Unexpected allocation failure poisons the owner, as retained-atom
        # publication does; recovery then starts from a checked checkpoint.
        collector.budget.publication()
        collector._publication_failed = True
        try:
            self.atoms[atom.opaque] = atom
            for residual in retired:
                del self.atoms[residual]
            self.used = used
            stats["opaque_staged"] = stats.get("opaque_staged", 0) + 1
            stats["opaque_staging_evictions"] = stats.get(
                "opaque_staging_evictions", 0
            ) + len(retired)
        except MemoryError as error:
            raise RelationPublicationError(
                "opaque staging publication failed"
            ) from error
        else:
            collector._publication_failed = False

    def discard(self, residual):
        atom = self.atoms.pop(residual)
        self.used -= entry_bytes(atom, self.base)


def route(collector, atom, stats):
    """Return True when the incoming atom was staged or was a duplicate.

    Promotion publishes the older atom first. A later refusal preserves that
    substep and retries only the incoming atom, just as opaque conversion does.
    """
    from .extended_collector import admit_retained

    staging = collector._opaque_staging
    opaque = getattr(atom, "opaque", None)
    residual = atom.residual_product
    old = staging.atoms.get(residual)
    if old is not None:
        if old.relation_id == atom.relation_id:
            stats["duplicates"] += 1
            return True, None
        collector._checked.atom(old, collector.budget)
        if old.opaque != residual:
            raise ValueError("staged opaque identity changed")
        key = None
        if collector.config.duplicate_policy == "congruence-v2":
            key, divisor = collector._checked.congruence(old, collector.budget)
            if divisor is not None:
                raise ValueError("staged opaque unit changed")
        incoming_key = collector._incoming_key
        collector._incoming_key = key
        try:
            refusal = admit_retained(collector, old, stats)
        finally:
            collector._incoming_key = incoming_key
        if refusal:
            return True, refusal
        staging.discard(residual)
        stats["opaque_promotions"] = stats.get("opaque_promotions", 0) + 1
    if (
        opaque is not None
        and residual not in collector._graph.resolutions
        and opaque_key(residual) not in collector._graph.direct
    ):
        staging.stage(atom, collector.budget, stats, collector)
        return True, None
    return False, None


def stage_unmatched(collector, atom, stats):
    """Stage an already verified atom before building an unused exact key."""
    opaque = getattr(atom, "opaque", None)
    if (
        opaque is None
        or opaque in collector._opaque_staging.atoms
        or opaque in collector._graph.resolutions
        or opaque_key(opaque) in collector._graph.direct
    ):
        return False
    # Congruence normalization would first expose this nonunit. Preserve that
    # check even though an unmatched staged record needs no canonical key.
    if collector.config.duplicate_policy == "congruence-v2":
        divisor = collector._checked.unit(
            atom.square_correction, collector.budget
        )
        if utils.valid_divisor(divisor, collector.factor_base.n):
            raise ValueError(
                "proper divisor must be returned before admission"
            )
    collector._incoming_key = None
    collector._opaque_staging.stage(atom, collector.budget, stats, collector)
    return True
