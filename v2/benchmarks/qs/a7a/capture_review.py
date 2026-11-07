"""Bounded review of every diagnostic companion before later selections."""

import json
import time
from collections import Counter
from pathlib import Path

from .diagnostic_analysis import heap_census


def review(
    capture, result, *, source, input_id, expected_pids=None, output=None
):
    """Review the declared CPU targets and their heap/buffer snapshots.

    The optional vmprof reader belongs to the analysis environment, never the
    factorer or worker environment. Heap records are streamed; no retained
    object graph or allocation-rate estimate is inferred from their census.
    Default to all recorded processes; bounded windows must pass their exact
    targets so an unsampled worker is distinct from a missing required stream.
    """
    import vmprof

    started, cpu = time.monotonic(), time.process_time()
    folders = {
        int(folder.name): folder
        for folder in Path(capture).iterdir()
        if folder.is_dir() and (folder / "events.json").exists()
    }
    expected = set(folders) if expected_pids is None else set(expected_pids)
    if not expected or not expected.issubset(folders):
        raise ValueError("missing requested capture processes")
    processes = []
    aggregate = {"parent": Counter(), "worker": Counter()}
    callers = {"parent": Counter(), "worker": Counter()}
    for pid in sorted(expected):
        folder = folders[pid]
        events = json.loads((folder / "events.json").read_text())
        ready = events[0]
        kinds = [event["kind"] for event in events]
        if (
            ready["pid"] != pid
            or ready["role"] not in aggregate
            or "recorder_failure" in kinds
            or kinds.count("cpu_started") != kinds.count("cpu_stopped")
            or "cpu_stopped" not in kinds
        ):
            raise ValueError("requested CPU capture lacks a valid flush")
        segments, leaves, samples = [], Counter(), 0
        for path in sorted(folder.glob("*.vmprof")):
            stats = vmprof.read_profile(str(path))
            segment = Counter()
            for trace, weight, *_ in stats.profiles:
                names = [
                    stats.adr_dict[a] for a in trace if a in stats.adr_dict
                ]
                segment[names[-1] if names else "unresolved"] += weight
                if len(names) >= 2:
                    edge = names[-1], names[-2]
                    callers[ready["role"]][edge] += weight
            leaves.update(segment)
            samples += stats.n_samples
            segments.append(
                dict(
                    file=path.name,
                    samples=stats.n_samples,
                    bytes=path.stat().st_size,
                )
            )
        aggregate[ready["role"]].update(leaves)
        row = dict(
            pid=ready["pid"],
            role=ready["role"],
            samples=samples,
            segments=segments,
            top_last_resolved=leaves.most_common(20),
        )
        heap_events = [
            event for event in events if event["kind"] == "heap_finished"
        ]
        if heap_events:
            row["heap"] = heap_census(folder / "live.heap")
            row["heap_events"] = heap_events
        processes.append(row)
    if not processes or any(not row["segments"] for row in processes):
        raise ValueError("missing owned CPU sample streams")
    rss = result.get("rss_samples", [])
    record = dict(
        source=source,
        input=input_id,
        processes=processes,
        expected_pids=sorted(expected),
        unreviewed_pids=sorted(set(folders) - expected),
        cpu_top_by_role={
            key: counter.most_common(30) for key, counter in aggregate.items()
        },
        owned_rss_sample_max=max((row["bytes"] for row in rss), default=None),
        callers_by_role={
            key: counter.most_common(30) for key, counter in callers.items()
        },
        rss_samples=len(rss),
        observer_failure=result.get("observer_failure"),
        wall=time.monotonic() - started,
        cpu=time.process_time() - cpu,
        limitations=[
            (
                "Last-resolved CPU samples can include unresolved JIT "
                "descendants; not exclusive or removable wall shares."
            ),
            (
                "Heap census is live GC object data, not dominators, "
                "allocation churn or unique physical memory."
            ),
            (
                "Owned RSS samples count shared pages and can miss peaks "
                "between observations."
            ),
        ],
    )
    path = (
        Path(output) if output is not None else Path(capture) / "analysis.json"
    )
    with path.open("x") as stream:
        json.dump(record, stream, indent=2)
        stream.write("\n")
    return record
