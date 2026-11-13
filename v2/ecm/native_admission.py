"""ECM-specific pre-arithmetic admission and retained-build invalidation."""

import hashlib
from dataclasses import dataclass
from pathlib import Path


class NativeUnavailableError(RuntimeError):
    """An optional native build could not be admitted before arithmetic."""

    def __init__(self, reason, detail):
        self.reason = reason
        super().__init__(detail)


class _AdmissionBudget:
    """Keep resource/provider faults outside optional-build fallback."""

    def __init__(self, budget):
        self.budget = budget
        self.failure = None

    def consume(self, amount=1):
        try:
            return self.budget.consume(amount)
        except BaseException as error:
            self.failure = error
            raise

    def poll(self):
        return self.consume(0)


def _fingerprint(path):
    value = Path(path).stat()
    return (
        value.st_dev,
        value.st_ino,
        value.st_size,
        value.st_mtime_ns,
        value.st_ctime_ns,
    )


def build_present(source, variant):
    """Check presence only; this never establishes trust or loads CFFI."""
    directory = Path(source) / "build"
    if variant != "scalar":
        directory /= variant
    return all(
        (directory / name).is_file()
        for name in (
            "manifest.json",
            "libfactor_ecm.dylib",
        )
    )


@dataclass(frozen=True)
class NativeAdmission:
    """One owner's validated immutable image; no global availability cache.

    The existing packet owns its context, buffers and library handle. This
    record owns no foreign resources. Fingerprints invalidate reuse; they
    never replace the packet loader's initial actual-byte/ABI validation.
    """

    identity: dict
    files: tuple

    def check(self):
        """Reject changed/missing files after admission without rerouting."""
        for name, expected in self.files:
            try:
                current = _fingerprint(name)
            except OSError as error:
                raise ValueError(
                    "admitted native ECM build disappeared"
                ) from error
            if current != expected:
                raise ValueError("admitted native ECM build changed")


def admit_packet(packet, n, budget, *, expected=None):
    """Admit a fresh whole-native packet before any curve action.

    Callers may handle NativeUnavailableError only for a new automatic
    decision. A restored expected identity and all post-admission faults
    refuse rather than silently starting another mathematical search.
    """
    from ..common.arithmetic import BackendUnavailableError
    from . import native_packet

    if not getattr(packet, "full_ecm", False) or packet.calls:
        raise RuntimeError("admission requires an unused whole-native packet")
    source = native_packet.SOURCE
    directory = source / "build"
    if packet.variant != "scalar":
        directory /= packet.variant
    manifest_path = directory / "manifest.json"
    library_path = directory / "libfactor_ecm.dylib"
    checked = _AdmissionBudget(budget)
    try:
        before = (_fingerprint(manifest_path), _fingerprint(library_path))
        if not packet.available(n, checked):
            raise BackendUnavailableError("native ECM build unavailable")
        if int(packet.lib.b4_whole_ecm()) != 1:
            raise ValueError("native build lacks whole-ECM execution")
        after = (_fingerprint(manifest_path), _fingerprint(library_path))
        if before != after:
            raise ValueError("native ECM build changed during admission")
    except (
        BackendUnavailableError,
        OSError,
        ValueError,
        KeyError,
        TypeError,
        AttributeError,
    ) as error:
        if checked.failure is not None:
            raise checked.failure
        if expected is not None:
            raise ValueError(
                "saved native ECM build is unavailable"
            ) from error
        reason = (
            "native_unavailable"
            if isinstance(error, (BackendUnavailableError, OSError))
            else "incompatible_native_build"
        )
        raise NativeUnavailableError(reason, str(error)) from error

    try:
        return _loaded_admission(packet, expected)
    except (OSError, ValueError, KeyError, TypeError, AttributeError) as error:
        if expected is not None:
            raise ValueError(
                "saved native ECM identity is unavailable"
            ) from error
        raise NativeUnavailableError(
            "incompatible_native_build", str(error)
        ) from error


def _loaded_admission(packet, expected):
    from . import native_packet

    source = native_packet.SOURCE
    directory = source / "build"
    if packet.variant != "scalar":
        directory /= packet.variant
    manifest_path, library_path = (
        directory / "manifest.json",
        directory / "libfactor_ecm.dylib",
    )
    manifest = packet.manifest
    paths = [manifest_path, library_path]
    binding = manifest.get("abi_binding")
    if binding is not None:
        paths.append(directory / binding["file"])
    paths.append(Path(manifest["gmp_library"]))
    paths.extend(source / name for name in manifest["source"])
    paths.extend(
        Path(native_packet.__file__).with_name(name)
        for name in (
            "native_packet.py",
            "native_curve.py",
            "native_curve_team.py",
            "native_resources.py",
            "native_metadata.py",
            "native_catalogue.py",
            "native_segment.py",
            "chains.py",
            "chain_records.py",
            "chain_data.py",
        )
    )
    paths.append(
        Path(native_packet.__file__).parents[1]
        / "execution/ecm_native_full.py"
    )
    activation_paths = [
        Path(__file__),
        Path(__file__).with_name("native_auto.py"),
        Path(__file__).with_name("execution_identity.py"),
        Path(__file__).with_name("auto.py"),
        Path(__file__).parents[1] / "execution/portfolio_config.py",
        Path(__file__).parents[1] / "portfolio.py",
    ]
    paths.extend(activation_paths)
    current_source = native_packet.source_identity.__wrapped__()
    if current_source != native_packet.source_identity():
        raise ValueError(
            "native adapter changed in this process; restart required"
        )
    activation = hashlib.sha256()
    for path in activation_paths:
        activation.update(path.read_bytes())
    identity = dict(
        source=current_source,
        variant=packet.variant,
        activation=activation.hexdigest(),
        library_sha256=manifest["library_sha256"],
        abi=manifest["abi"],
        binding_sha256=None if binding is None else binding["sha256"],
        gmp_sha256=manifest["gmp_library_sha256"],
        native_team_max_workers=manifest["native_team_max_workers"],
        slice_actions=manifest.get("slice_actions"),
    )
    if expected is not None and identity != expected:
        raise ValueError("saved native ECM build identity changed")
    return NativeAdmission(
        identity,
        tuple(
            (str(path), _fingerprint(path)) for path in dict.fromkeys(paths)
        ),
    )


def verify_saved_builds(records, budget):
    """Revalidate native-bound resumes even when the next stage is QS/prime."""
    from .native_curve import NativeCurvePacket

    checked = set()
    for record in records:
        identity = record.get("build")
        if record.get("backend") != "c-gmp-full" or identity is None:
            continue
        key = (identity.get("variant"), identity.get("library_sha256"))
        if key in checked:
            continue
        packet = NativeCurvePacket(mode="c-gmp", _variant=identity["variant"])
        try:
            admit_packet(packet, record["n"], budget, expected=identity)
            checked.add(key)
        finally:
            packet.close()


def admit_team(team, n, budget, *, expected=None):
    """Validate the first binding before any native worker is started."""
    from ..common.arithmetic import BackendUnavailableError

    checked = _AdmissionBudget(budget)
    try:
        team._grow(1, n, checked)
    except (
        BackendUnavailableError,
        OSError,
        ValueError,
        KeyError,
        TypeError,
        AttributeError,
    ) as error:
        if checked.failure is not None:
            raise checked.failure
        if expected is not None:
            raise ValueError("saved native ECM team is unavailable") from error
        raise NativeUnavailableError(
            "native_unavailable", str(error)
        ) from error
    return admit_packet(team.packets[0], n, budget, expected=expected)
