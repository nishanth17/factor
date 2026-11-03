"""Bounded FIFO ownership for opaque atoms outside the incidence history."""

from collections import OrderedDict


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

    def stage(self, atom, budget, stats):
        cost = entry_bytes(atom, self.base)
        if (
            not self.config.opaque_cache_entries
            or cost > self.config.opaque_cache_bytes
        ):
            stats["opaque_staging_rejections"] = (
                stats.get("opaque_staging_rejections", 0) + 1
            )
            return
        budget.consume(len(atom.exponents) + len(self.atoms) + 1)
        pending = OrderedDict(self.atoms)
        used = self.used
        pending[atom.opaque] = atom
        used += cost
        evictions = 0
        while (
            len(pending) > self.config.opaque_cache_entries
            or used > self.config.opaque_cache_bytes
        ):
            _, old = pending.popitem(last=False)
            used -= entry_bytes(old, self.base)
            evictions += 1
        # The configured cache grant covers retained records. The global
        # construction scratch covers these two small coexisting FIFO maps.
        self.atoms, self.used = pending, used
        stats["opaque_staged"] = stats.get("opaque_staged", 0) + 1
        stats["opaque_staging_evictions"] = (
            stats.get("opaque_staging_evictions", 0) + evictions
        )

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
        and frozenset(((True, residual),)) not in collector._graph.direct
    ):
        staging.stage(atom, collector.budget, stats)
        return True, None
    return False, None
