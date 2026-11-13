"""Lazy whole-native decisions, independent of the mathematical seed domain."""

from dataclasses import replace
from pathlib import Path

VERSION = "b4a-auto-full-v1"
SOURCE = Path(__file__).with_name("native")


def eligibility(n, config):
    """Reject unsupported or negligible work before native preparation."""
    if config.backend != "python-int":
        return "unsupported_automatic_backend"
    if type(n) is not int or not n % 2 or not 3 <= n < 10**100:
        return "unsupported_automatic_input"
    if n < 10**29:
        return "unmeasured_native_size"
    if (
        config.ecm_pair_distance is not None
        or config.ecm_pair_wheel is not None
        or config.ecm_chain_family != "auto"
    ):
        return "unsupported_native_formula"
    if sum(b1 * curves for b1, _, curves in config.ecm_tiers) < 16000:
        return "below_native_work_threshold"
    from .auto import _runtime_supported

    if not _runtime_supported():
        return "unsupported_automatic_runtime"
    return None


def decision(*, backend, reason, stream, build=None, workers=1):
    """Describe intent separately from execution and its actual entry count."""
    native = backend == "c-gmp-full"
    return dict(
        version=VERSION,
        requested="auto-full",
        backend=backend,
        implementation="whole-native" if native else "python",
        recipe=(
            "pthread/dynamic/chunk1/cap32"
            if workers > 1
            else "scalar/256"
            if native
            else "incumbent"
        ),
        topology="threads" if workers > 1 else "scalar",
        workers=workers,
        reason=reason,
        mathematical_stream=stream,
        build=build,
        native_entries=0,
    )


def prepare_serial(n, config, budget, context, *, saved=None):
    """Attach one admitted packet without changing parent RNG consumption.

    A saved fallback remains fallback. Once bound to C, resume requires that
    same build. Resource/provider errors and faults after admission propagate.
    Context teardown owns the packet exactly as for explicit whole-native use.
    """
    from .native_admission import (
        NativeUnavailableError,
        admit_packet,
        build_present,
    )

    if saved is not None and saved["backend"] == "python":
        old = getattr(context, "native_packet", None)
        if old is not None:
            old.close()
            context.native_packet = None
            context.automatic_native_admission = None
        return saved
    reason = eligibility(n, config)
    if reason is None and saved is None:
        reason = resource_headroom(budget)
    if reason is not None:
        if saved is not None:
            raise ValueError(
                "saved whole-native ECM is unavailable: " + reason
            )
        return prepare_serial(
            n,
            config,
            budget,
            context,
            saved=decision(
                backend="python", reason=reason, stream="parent-rng"
            ),
        )
    admission = getattr(context, "automatic_native_admission", None)
    if admission is not None:
        admission.check()
        if saved is not None and saved["build"] != admission.identity:
            raise ValueError("saved automatic native ECM build changed")
        return saved or decision(
            backend="c-gmp-full",
            reason="retained_native_scalar",
            stream="parent-rng",
            build=admission.identity,
        )
    if not build_present(SOURCE, "scalar"):
        reason = "native_unavailable"
    if reason is not None:
        if saved is not None:
            raise ValueError(
                "saved whole-native ECM is unavailable: " + reason
            )
        return decision(backend="python", reason=reason, stream="parent-rng")

    from .native_curve import NativeCurvePacket

    packet = NativeCurvePacket(mode="c-gmp")
    try:
        admission = admit_packet(
            packet,
            n,
            budget,
            expected=None if saved is None else saved["build"],
        )
    except NativeUnavailableError as error:
        packet.close()
        return decision(
            backend="python",
            reason=error.reason,
            stream="parent-rng",
        )
    except BaseException:
        packet.close()
        raise
    context.native_packet = packet
    context.automatic_native_admission = admission
    return saved or decision(
        backend="c-gmp-full",
        reason="eligible_native_scalar",
        stream="parent-rng",
        build=admission.identity,
    )


def entries(context):
    """Count actual whole-native calls, excluding admission/stage-one loads."""
    packet = getattr(context, "native_packet", None)
    return packet.calls if getattr(packet, "full_ecm", False) else 0


def _incumbent_manifest(owner, n, config, budget, lineage):
    """Determine the incumbent seed contract without running any curves.

    The B12 cost policy determines only whether the old call used Manifest
    seeds or parent RNG draws. Its process economics do not veto native scalar
    engagement or choose a native team size.
    """
    from .auto import _size_supported, cold_serial_default, selected_context

    if not _size_supported(n) or cold_serial_default(n, config):
        return False
    import multiprocessing

    if multiprocessing.parent_process() is not None:
        return False
    if owner.context is None:
        owner.context = selected_context(seed=owner.seed)
        try:
            owner.context._load(budget)
        except (OSError, ValueError):
            owner.context.close()
            owner.context = None
            owner.disabled_reason = "unavailable_automatic_policy"
            return False
    plan = owner.context._select(n, config, None, None, portfolio=True)
    if plan is None:
        return False
    settings, rule, _ = plan
    result = owner.context._portfolio_decision(
        n, config, budget, lineage, settings, rule
    )
    return result["engage"] and result["funded"]


def candidate_workers(config):
    """Provisional shortlist; freeze thresholds after exposed screening."""
    work = sum(b1 * curves for b1, _, curves in config.ecm_tiers)
    curves = sum(tier[2] for tier in config.ecm_tiers)
    if work >= 500000 and curves >= 32:
        return 8
    if work >= 50000 and curves >= 8:
        return 4
    return 1


def resource_headroom(budget):
    """Avoid native preparation when the remaining root grant is negligible."""
    base = getattr(budget, "budget", budget)
    for limit, used in (
        (base.seconds, base.wall_used),
        (base.cpu_seconds, base.cpu_used),
    ):
        if limit is not None and limit - used < 0.02:
            return "insufficient_native_headroom"
    return None


def _prepare_owner(owner, n, config, budget, workers, expected):
    from ..execution.worker_contract import Manifest, WorkerSettings
    from .native_admission import _AdmissionBudget, admit_packet, admit_team
    from .workers import ECMExecutor, arithmetic_options

    settings = WorkerSettings(
        workers=workers,
        allow_wide_workers=workers > 4,
        worker_imports="ecm",
        budget_mode="prospective",
    )
    if owner.native_context is None:
        if workers == 1:
            from .worker import CurveEngine

            selected = ECMExecutor(seed=owner.seed, settings=settings)
            selected.serial_engine = CurveEngine(settings)
        else:
            from ..execution.ecm_native import NativeECMExecutor

            selected = NativeECMExecutor(
                seed=owner.seed,
                settings=settings,
                capacity=32,
                schedule="dynamic",
                chunk_size=1,
            )
        owner.native_context = selected
    selected = owner.native_context
    if selected.settings != settings:
        raise ValueError("retained automatic native recipe changed")
    options = arithmetic_options(config, settings)
    # Only metadata is prepared here. Admission must precede curve actions.
    manifest = Manifest(n, owner.seed, owner._native_lineage, options)
    engine = selected.serial_engine if workers == 1 else selected.engine
    engine._prepare(
        dict(schedule=manifest.metadata, options=options),
        _AdmissionBudget(budget),
    )
    if owner.native_admission is None:
        if workers == 1:
            packet = (engine.programs or engine.context).native_packet
            admitted = admit_packet(packet, n, budget, expected=expected)
        else:
            from .native_curve_team import NativeCurveTeam

            selected.team = NativeCurveTeam(
                workers=workers,
                schedule="dynamic",
                capacity=32,
                chunk_size=1,
            )
            admitted = admit_team(selected.team, n, budget, expected=expected)
        owner.native_admission = admitted
    owner.native_admission.check()
    if expected is not None and owner.native_admission.identity != expected:
        raise ValueError("saved automatic native build changed")
    return selected


def run_portfolio(
    owner, n, *, config, budget, lineage="root", checkpoint=None
):
    """Select native execution while preserving the incumbent search stream."""
    from ..execution.worker_contract import digest
    from .native_admission import NativeUnavailableError, build_present

    if owner.closed:
        raise RuntimeError("automatic ECM owner is closed")
    budget.consume(0)
    incumbent = replace(config, ecm_kernel="python")
    saved = None
    if checkpoint is not None and "native_auto" in checkpoint:
        saved = checkpoint["native_auto"]
        if (
            saved.get("version") != VERSION
            or checkpoint.get("sha256") != digest(saved)
            or saved.get("inner") != checkpoint["inner"].get("sha256")
            or saved.get("workers") not in (1, 4, 8)
        ):
            raise ValueError("invalid automatic whole-native checkpoint")
    elif checkpoint is not None:
        # A previously selected Python fallback retains its exact checkpoint.
        owner.native_disabled_reason = "saved_python_fallback"

    reason = owner.native_disabled_reason or eligibility(n, config)
    if reason is None and saved is None:
        reason = resource_headroom(budget)
    if reason is None and saved is None:
        if not any(
            build_present(SOURCE, variant) for variant in ("scalar", "pthread")
        ):
            reason = "native_unavailable"
    if reason is not None:
        if saved is not None:
            raise ValueError("saved native ECM unavailable: " + reason)
        owner.native_disabled_reason = reason
        result = owner._run_incumbent(
            n,
            config=incumbent,
            budget=budget,
            lineage=lineage,
            checkpoint=checkpoint,
        )
        execution = decision(
            backend="python",
            reason=reason,
            stream="parent-rng" if result is None else "b12-manifest",
        )
        owner.last_stats = dict(owner.last_stats, ecm_execution=execution)
        return (
            None if result is None else replace(result, stats=owner.last_stats)
        )

    if saved is None and not _incumbent_manifest(
        owner, n, incumbent, budget, lineage
    ):
        owner.last_stats = dict(
            engaged=False,
            legacy_serial=True,
            engagement=dict(engage=False, reason="parent_rng_contract"),
        )
        return None

    workers = saved["workers"] if saved else candidate_workers(config)
    if saved is None and workers > 1 and not build_present(SOURCE, "pthread"):
        workers = 1
    if owner.native_context is not None:
        workers = owner.native_context.settings.workers
    concrete = replace(config, ecm_kernel="c-gmp-full")
    owner._native_lineage = lineage
    try:
        selected = _prepare_owner(
            owner,
            n,
            concrete,
            budget,
            workers,
            None if saved is None else saved["build"],
        )
    except NativeUnavailableError as error:
        # No curve action has occurred. Construction remains reachable by
        # close() even when validation failed partway through preparation.
        owner.native_context.close()
        owner.native_context = owner.native_admission = None
        owner.native_disabled_reason = error.reason
        return run_portfolio(
            owner,
            n,
            config=config,
            budget=budget,
            lineage=lineage,
        )
    before = (
        0
        if workers > 1
        else entries(
            selected.serial_engine.programs or selected.serial_engine.context
        )
    )
    result = selected.run(
        n,
        config=concrete,
        budget=budget,
        lineage=lineage,
        checkpoint=None if saved is None else checkpoint["inner"],
    )
    count = (
        result.stats.get("native_slices", 0)
        if workers > 1
        else entries(
            selected.serial_engine.programs or selected.serial_engine.context
        )
        - before
    )
    execution = decision(
        backend="c-gmp-full",
        reason="eligible_native_manifest",
        stream="b12-manifest",
        build=owner.native_admission.identity,
        workers=workers,
    )
    execution["native_entries"] = count
    execution["tiers"] = [list(tier) for tier in config.ecm_tiers]
    header = dict(
        version=VERSION,
        workers=workers,
        build=owner.native_admission.identity,
        inner=result.checkpoint["sha256"],
    )
    owner.last_stats = dict(result.stats, ecm_execution=execution)
    return replace(
        result,
        stats=owner.last_stats,
        checkpoint=dict(
            native_auto=header, sha256=digest(header), inner=result.checkpoint
        ),
    )
