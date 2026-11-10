"""Scoped whole-curve OpenMP experiment with explicit selection."""

import hashlib
from pathlib import Path

from ....ecm import native_packet
from ....ecm.native_curve_team import NativeCurveTeam
from ....execution.ecm_native import NativeECMExecutor
from ....execution.worker_contract import digest


class OpenMPCurveTeam(NativeCurveTeam):
    """Use private contexts and the admitted native OpenMP owner API."""

    variant = "openmp"
    owner_prefix = "b4_omp_"


class OpenMPExecutor(NativeECMExecutor):
    """Reuse exact assignment/history with a distinct backend identity."""

    _full_team_type = OpenMPCurveTeam

    def __init__(self, **kwargs):
        manifest = native_packet.SOURCE / "build/openmp/manifest.json"
        if not manifest.is_file():
            raise RuntimeError("explicit whole-curve OpenMP build required")
        self._backend_identity = hashlib.sha256(
            Path(__file__).read_bytes() + manifest.read_bytes()
        ).hexdigest()
        super().__init__(**kwargs)

    def _identity(self):
        return digest(
            [
                super()._identity(),
                "b4a-openmp-experiment-v1",
                self._backend_identity,
            ]
        )
