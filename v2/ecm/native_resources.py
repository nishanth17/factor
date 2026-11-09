"""Exact native root observations at the existing Python resource polls."""

import errno
import os

from ..execution.resources import (
    MacResources,
    ProcessSample,
    ResourceProviderError,
)
from .native_packet import poll_budget


class NativeRootProvider:
    """Use the original tree provider whenever a childless root is unproved."""

    def __init__(self, owner):
        self.ffi, self.lib = owner.ffi, owner.lib
        self.output = self.ffi.new("b4_root_observation *")
        self.fallback = None
        self.pid = os.getpid()

    def sample(self, roots):
        if not isinstance(roots, dict) or set(roots) != {self.pid}:
            return self.fallback.sample(roots)

        status = self.lib.b4_observe_root(self.pid, self.output)
        if status == 1:
            return self.fallback.sample(roots)
        if status == -errno.ESRCH:
            return []
        if status:
            raise ResourceProviderError(
                "native root resource observation: " + str(status)
            )

        value = self.output[0]
        if roots[self.pid] is not None and roots[self.pid] != value.birth:
            return []
        return [
            ProcessSample(
                self.pid,
                int(value.birth),
                int(value.rss),
                (int(value.user) + int(value.system)) / 1_000_000_000,
            )
        ]

    def close(self):
        self.ffi.release(self.output)
        self.fallback = None


def poll_native_budget(budget, owner):
    """Preserve every policy check, replacing only the root OS transport."""
    from ..execution.ecm_workers import _base_budget

    base = _base_budget(budget)
    original = getattr(base, "provider", None)
    if type(original) is not MacResources or owner.lib is None:
        return poll_budget(budget)

    provider = getattr(owner, "_resource_provider", None)
    if provider is None:
        provider = owner._resource_provider = NativeRootProvider(owner)
    provider.fallback = original
    base.provider = provider
    try:
        return poll_budget(budget)
    finally:
        base.provider = original
        provider.fallback = None


def close_resource_checks(owner):
    provider = getattr(owner, "_resource_provider", None)
    if provider is not None:
        provider.close()
        owner._resource_provider = None
