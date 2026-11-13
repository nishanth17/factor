"""Actual ECM execution receipts and non-forcing benchmark assertions."""


def append_decision(state, record):
    """Retain bounded receipts and expose every omitted decision."""
    history = state.setdefault("ecm_execution", dict(jobs=[], dropped=0))
    if len(history["jobs"]) >= 4096:
        history["dropped"] += 1
        return None
    history["jobs"].append(record)
    return len(history["jobs"]) - 1


def executor_packet(executor):
    """Inspect existing ownership without preparing or importing a backend."""
    team = getattr(executor, "team", None)
    if team is not None and getattr(team, "packets", ()):
        return team.packets[0]
    engine = getattr(executor, "serial_engine", None)
    context = None if engine is None else engine.programs or engine.context
    return getattr(context, "native_packet", None)


def executor_entries(executor):
    packet = executor_packet(executor)
    if not getattr(packet, "full_ecm", False):
        return 0
    team = getattr(executor, "team", None)
    return team.slices if team is not None else packet.calls


def packet_identity(packet):
    """Report a binding already loaded by explicit native execution."""
    if packet is None or getattr(packet, "manifest", None) is None:
        return None
    from .native_packet import source_identity

    manifest = packet.manifest
    return dict(
        source=source_identity(),
        variant=packet.variant,
        library_sha256=manifest["library_sha256"],
        abi=manifest["abi"],
        gmp_sha256=manifest["gmp_library_sha256"],
        binding_sha256=(manifest.get("abi_binding") or {}).get("sha256"),
    )


def summary(state, requested):
    """Summarize actual paths while preserving every retained job decision."""
    history = state.get("ecm_execution", dict(jobs=[], dropped=0))
    jobs = history["jobs"]
    backends = sorted({record["backend"] for record in jobs})
    return dict(
        schema="ecm-execution-v1",
        requested=requested,
        backend=(
            backends[0]
            if len(backends) == 1
            else "mixed"
            if backends
            else "not-engaged"
        ),
        reason="ecm_entered" if jobs else "completed_or_stopped_before_ecm",
        native_entries=sum(record.get("native_entries", 0) for record in jobs),
        jobs=jobs,
        dropped=history["dropped"],
    )


def validate_history(state):
    """Reject malformed saved receipts before they can become resume state."""
    history = state.get("ecm_execution")
    if history is None:
        return
    if (
        not isinstance(history, dict)
        or set(history) != {"jobs", "dropped"}
        or not isinstance(history["jobs"], list)
        or len(history["jobs"]) > 4096
        or type(history["dropped"]) is not int
        or history["dropped"] < 0
    ):
        raise ValueError("invalid ECM execution history")
    for record in history["jobs"]:
        if (
            not isinstance(record, dict)
            or record.get("backend")
            not in ("python", "c-gmp", "c-gmp-full", "auto")
            or type(record.get("native_entries")) is not int
            or record["native_entries"] < 0
            or type(record.get("n")) is not int
            or record["n"] < 3
            or abs(state["original"]) % record["n"]
        ):
            raise ValueError("invalid saved ECM execution decision")
        if (
            record["backend"] == "c-gmp-full"
            and record.get("requested") == "auto-full"
        ):
            if not isinstance(record.get("build"), dict):
                raise ValueError("unbound automatic native execution")


def require_native(run, *, workers=None):
    """Assert ordinary whole-native execution after a measured call.

    This assertion never selects, imports or preloads a backend. Runners must
    save the returned observation and paid costs before invoking it, so a
    fallback or partial result remains in the intended-native denominator.
    """
    record = getattr(run, "ecm_execution", None)
    if (
        not record
        or not record["jobs"]
        or record["backend"] != "c-gmp-full"
        or not record["native_entries"]
        or record["dropped"]
    ):
        raise ValueError(
            "require-native: whole-native ECM did not exclusively execute"
        )
    for job in record["jobs"]:
        if (
            job.get("implementation") != "whole-native"
            or job.get("backend") != "c-gmp-full"
        ):
            raise ValueError("require-native: wrong native implementation")
        if workers is not None and job.get("workers") != workers:
            raise ValueError("require-native: wrong scalar/thread recipe")
    return record
