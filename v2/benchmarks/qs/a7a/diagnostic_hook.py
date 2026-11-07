"""Opt-in diagnostic initializer, absent from ordinary factoring calls."""

import atexit
import gc
import json
import os
import resource
import signal
import sys
import time
from pathlib import Path

_folder = None
_role = None
_profile = None
_cpu_file_limit = None
_events = []
_cpu_segments = 0
_closed = False


def _record(kind, **values):
    _events.append(
        dict(
            kind=kind, wall=time.monotonic(), cpu=time.process_time(), **values
        )
    )
    temporary = _folder / "events.tmp"
    temporary.write_text(json.dumps(_events, indent=2) + "\n")
    temporary.replace(_folder / "events.json")


def _frames(frame):
    """Retain names and cheap scalars, without serializing arbitrary locals."""
    rows = []
    while frame is not None:
        row = dict(
            function=frame.f_code.co_name,
            file=frame.f_code.co_filename,
            line=frame.f_lineno,
        )
        owner = frame.f_locals.get("self")
        if owner is not None:
            row["class"] = type(owner).__name__
            state = getattr(owner, "__dict__", {})
            row["scalars"] = {
                k: v
                for k, v in state.items()
                if type(v) in (int, float, bool)
                and any(
                    s in k for s in ("count", "polynomial", "block", "used")
                )
            }
            collector = state.get("collector")
            if collector is None and state.get("engine") is not None:
                collector = getattr(state["engine"], "collector", None)
            if collector is not None:
                atoms = getattr(collector, "_atoms", None)
                if atoms is not None:
                    row["retained_atoms"] = len(atoms)
        rows.append(row)
        frame = frame.f_back
    return rows


def _cpu_signal(signum, frame):
    global _profile, _cpu_file_limit, _cpu_segments
    import _vmprof

    try:
        if _profile is None:
            _cpu_file_limit = resource.getrlimit(resource.RLIMIT_FSIZE)
            previous_soft = _cpu_file_limit[0]
            soft = (
                64 * 2**20
                if previous_soft == resource.RLIM_INFINITY
                else min(64 * 2**20, previous_soft)
            )
            resource.setrlimit(
                resource.RLIMIT_FSIZE,
                (soft, _cpu_file_limit[1]),
            )
            _cpu_segments += 1
            name = (
                "cpu.vmprof"
                if _cpu_segments == 1
                else f"cpu-{_cpu_segments:03}.vmprof"
            )
            interval = (
                0.01
                if os.environ.get("FACTOR_A7A_PROFILE_ALL") == "1"
                else 0.002
            )
            _profile = (_folder / name).open("xb")
            _vmprof.enable(
                _profile.fileno(), interval, False, False, False, False
            )
            _record(
                "cpu_started",
                location=_location(frame),
                file=name,
                interval=interval,
            )
        else:
            _vmprof.disable()
            _profile.close()
            _profile = None
            resource.setrlimit(resource.RLIMIT_FSIZE, _cpu_file_limit)
            _cpu_file_limit = None
            _record("cpu_stopped", location=_location(frame))
    except Exception as error:
        try:
            if _vmprof.is_enabled():
                _vmprof.disable()
            if _profile is not None:
                _profile.close()
        finally:
            _profile = None
            if _cpu_file_limit is not None:
                resource.setrlimit(resource.RLIMIT_FSIZE, _cpu_file_limit)
                _cpu_file_limit = None
        _record("recorder_failure", operation="cpu", error=repr(error))


def _location(frame):
    # CPU stacks are already in vmprof. Reading every frame's f_locals here
    # can force PyPy JIT temporaries and distort the rest of the whole call.
    return (
        None
        if frame is None
        else dict(
            function=frame.f_code.co_name,
            file=frame.f_code.co_filename,
            line=frame.f_lineno,
        )
    )


def _heap_signal(signum, frame):
    # The supervisor stops CPU sampling first. A heap pause must never become
    # evidence about production arithmetic, memory refusal or speed.
    automatic = os.environ.get("FACTOR_A7A_PROFILE_ALL") == "1"
    resume_cpu = automatic and _profile is not None
    if _profile is not None and not automatic:
        _record("recorder_failure", operation="heap", error="CPU still active")
        return
    if resume_cpu:
        # Pauses are absent from CPU attribution. Resume into a new immutable
        # segment; never truncate or append over an earlier sample stream.
        _cpu_signal(None, frame)
    started, cpu = time.monotonic(), time.process_time()
    frames = _frames(frame)
    if os.environ.get("FACTOR_A7A_LIVE_HEAPS") == "1":
        active = any(r["function"] == "collect_packet" for r in frames)
        graph = any(r.get("retained_atoms", 0) > 0 for r in frames)
        if not (active if _role == "worker" else graph):
            _record("heap_deferred", frames=frames)
            return
    owners = _logical_owners(frame)
    before = str(gc.get_stats())
    previous = resource.getrlimit(resource.RLIMIT_FSIZE)
    limit = (2 if _role == "parent" else 1) * 2**30
    try:
        resource.setrlimit(resource.RLIMIT_FSIZE, (limit, previous[1]))
        gc.dump_rpy_heap(str(_folder / "live.heap"))
        _record(
            "heap_finished",
            seconds=time.monotonic() - started,
            cpu_seconds=time.process_time() - cpu,
            bytes=(_folder / "live.heap").stat().st_size,
            gc_before=before,
            gc_after=str(gc.get_stats()),
            rss_hwm=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            frames=frames,
            logical_owners=owners,
            phase="teardown" if frame is None else "live",
        )
    except Exception as error:
        _record("recorder_failure", operation="heap", error=repr(error))
    finally:
        resource.setrlimit(resource.RLIMIT_FSIZE, previous)
        if resume_cpu and not _closed:
            _cpu_signal(None, frame)


def _logical_owners(frame):
    """Bounded owner counts, separate from GC retained-size analysis."""
    owners, seen = [], set()
    while frame is not None and len(owners) < 64:
        values = list(frame.f_locals.items())
        for name, owner in values:
            cls = type(owner)
            if not cls.__module__.startswith("v2.qs.") or id(owner) in seen:
                continue
            candidates = [(name, owner)]
            engine = getattr(owner, "engine", None)
            collector = getattr(owner, "collector", None)
            if engine is not None:
                candidates.append((name + ".engine", engine))
                collector = getattr(engine, "collector", None)
            if collector is not None:
                candidates.append((name + ".collector", collector))
            for label, candidate in candidates:
                if id(candidate) in seen or len(owners) >= 64:
                    continue
                seen.add(id(candidate))
                state = getattr(candidate, "__dict__", {})
                containers = {
                    key: len(value)
                    for key, value in state.items()
                    if type(value)
                    in (dict, list, tuple, set, bytearray, bytes)
                }
                scalars = {
                    key: value
                    for key, value in state.items()
                    if type(value) is int
                    and any(
                        s in key for s in ("bytes", "workspace", "capacity")
                    )
                }
                owners.append(
                    dict(
                        owner=label,
                        type=type(candidate).__name__,
                        containers=containers,
                        accounted_bytes=scalars,
                    )
                )
                if type(candidate).__name__ == "SieveCollector":
                    derived = {}
                    schedule = state.get("_sparse_queue")
                    if schedule is not None:
                        queues = schedule[1]
                        derived["queue_blocks"] = len(queues)
                        derived["queue_payload_bytes"] = sum(
                            len(queue) * queue.itemsize for queue in queues
                        )
                        derived["queue_numeric_words"] = sum(
                            len(queue) for queue in queues
                        )
                    words = state.get("_word_hits")
                    derived["hit_backend"] = (
                        "list" if words is None else "WordHits"
                    )
                    inverses = state.get("_shared_lift_inverses")
                    derived["shared_inverse_entries"] = (
                        None if inverses is None else len(inverses.values)
                    )
                    transfers = state.get("_power_transfers")
                    if transfers is not None:
                        derived["root_transfer_entries"] = len(
                            transfers.values
                        )
                        derived["root_transfer_reserved"] = transfers.used
                        derived["root_transfer_builds"] = transfers.builds
                        derived["root_transfer_replays"] = transfers.replays
                    if words is not None:
                        derived["hit_word_groups"] = len(words.words)
                        derived["hit_word_cells"] = sum(map(len, words.words))
                        derived["hit_stamp_cells"] = sum(
                            map(len, words.stamps)
                        )
                    owners[-1]["derived_buffers"] = derived
        frame = frame.f_back
    return owners


def install(role):
    global _folder, _role
    _role = role
    _folder = Path(os.environ["FACTOR_A7A_DIAGNOSTIC_DIR"]) / str(os.getpid())
    _folder.mkdir()
    signal.signal(signal.SIGUSR1, _cpu_signal)
    signal.signal(signal.SIGUSR2, _heap_signal)
    _record(
        "ready",
        pid=os.getpid(),
        ppid=os.getppid(),
        role=role,
        runtime=sys.version,
        startup_gc=str(gc.get_stats()),
    )

    def close():
        global _closed
        if _closed:
            return
        _closed = True
        if _profile is not None:
            _cpu_signal(None, None)
        if os.environ.get("FACTOR_A7A_PROFILE_ALL") == "1":
            targets = _folder.parent / "heap-targets.json"
            selected = role == "parent"
            if targets.exists():
                chosen = json.loads(targets.read_text())
                selected = (
                    selected or os.getpid() in chosen or len(chosen) == 1
                )
            elif role == "worker":
                selected = True
            if selected and not any(
                e["kind"] == "heap_finished" for e in _events
            ):
                _heap_signal(None, None)

    atexit.register(close)
    if role == "worker" and os.environ.get("FACTOR_A7A_PROFILE_ALL") == "1":
        from multiprocessing.util import Finalize

        # Spawned workers can exit through multiprocessing finalizers rather
        # than ordinary interpreter atexit. The shared close is idempotent.
        Finalize(None, close, exitpriority=0)
        _cpu_signal(None, None)
    return close


def initialize_worker(state):
    from ....qs.parallel import _initialize_process

    _initialize_process(state)
    install("worker")


def main():
    from ....qs import parallel
    from . import cell

    close = install("parent")
    # Only this diagnostic process selects the wrapper. Spawn imports its
    # original production initializer, then installs identical signal hooks.
    parallel._initialize_process = initialize_worker
    if os.environ.get("FACTOR_A7A_PROFILE_ALL") == "1":
        execute = cell.execute

        def profiled_execute(request):
            # cell.main has installed its finite file limit by this point.
            # Saving/restoring that limit cannot accidentally raise a hard
            # limit that the request subsequently lowered.
            _cpu_signal(None, None)
            try:
                return execute(request)
            finally:
                close()

        cell.execute = profiled_execute
    cell.main()


if __name__ == "__main__":
    main()
