"""Explicit whole-curve OpenMP owner for checked portfolio ECM calls."""

import hashlib
from pathlib import Path

from ..ecm import native_packet
from ..ecm.native_curve_team import OpenMPCurveTeam
from .ecm_native import NativeECMExecutor
from .worker_contract import digest

VERSION = "b4a-openmp-owner-h3-v1"


class OpenMPExecutor(NativeECMExecutor):
    """Retain private native contexts under a backend-bound checkpoint."""

    _full_team_type = OpenMPCurveTeam

    def _identity(self):
        manifest = native_packet.SOURCE / "build/openmp/manifest.json"
        try:
            build = manifest.read_bytes()
        except FileNotFoundError:
            # Cheap portfolio exits can finish before the optional backend is
            # needed. Eligible work still performs strict native admission.
            build = b"unavailable"
        backend = hashlib.sha256(
            Path(__file__).read_bytes() + build
        ).hexdigest()
        return digest([super()._identity(), VERSION, backend])

    def _team(self):
        if self.team is None:
            self.team = OpenMPCurveTeam(
                workers=self.settings.workers,
                schedule=self.schedule,
                capacity=self.capacity,
                chunk_size=self.chunk_size,
            )
        return self.team

    def _run(self, n, config, budget, lineage, checkpoint, fixed_work):
        if config.ecm_kernel != "c-gmp-full":
            raise ValueError("OpenMP owner requires c-gmp-full ECM")
        return super()._run(n, config, budget, lineage, checkpoint, fixed_work)
