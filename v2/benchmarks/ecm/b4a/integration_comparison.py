"""Check common canonical QS work and sustained native-owner quiescence."""


def compare_prefixes(left, right, left_path, right_path):
    if left["jobs"] != right["jobs"]:
        raise ArithmeticError("QS configuration/search identities differ")
    a = {(r["job"], r["index"]): r for r in left["assignments"]}
    b = {(r["job"], r["index"]): r for r in right["assignments"]}
    shared = sorted(a.keys() & b.keys())
    common = 0
    with left_path.open("rb") as first, right_path.open("rb") as second:
        for key in shared:
            count = min(a[key]["rows"], b[key]["rows"])
            first.seek(a[key]["offset"])
            second.seek(b[key]["offset"])
            one, two = first.read(count * 32), second.read(count * 32)
            if len(one) != count * 32 or one != two:
                raise ArithmeticError("canonical published QS prefix differs")
            common += count
    return dict(
        common_assignments=len(shared),
        common_published_rows=common,
        left_published_rows=sum(r["rows"] for r in a.values()),
        right_published_rows=sum(r["rows"] for r in b.values()),
        left_extra_published_rows=sum(r["rows"] for r in a.values()) - common,
        right_extra_published_rows=sum(r["rows"] for r in b.values()) - common,
        scope="Common published atom prefixes match; unmerged peer work is "
        "separately visible in production QS counters.",
    )


def native_handoff(evidence):
    active = [
        s for s in evidence["snapshots"] if s["native_threads"] is not None
    ]
    closed = [
        s for s in evidence["snapshots"] if s["label"] == "after-ecm-close"
    ]
    if not closed or not closed[-1]["native"]["closed"]:
        raise RuntimeError("explicit ECM owner did not close")
    state = closed[-1]["native"]
    if state.get("native"):
        if (
            state["contexts"]
            or state["live_contexts"]
            or state["shared_native_blocks"]
        ):
            raise RuntimeError(
                "native contexts/blocks remain live after close"
            )
        lifecycle = state["lifecycle"]
        if lifecycle is not None and lifecycle[0] != lifecycle[1]:
            raise RuntimeError("native thread creation/join count differs")
    if not active:
        scalar = any(
            s["label"] == "qs-entry"
            and s["native"].get("topology") == "scalar"
            and s["native"].get("native")
            for s in evidence["snapshots"]
        )
        return dict(
            reached_native_QS=scalar,
            topology="scalar" if scalar else "not-engaged",
            closed=True,
            contexts_released=True,
            final_threads_joined=0,
            scope="Scalar ECM owns no background workers; its context is "
            "accounted during QS and released by final call cleanup."
            if scalar
            else "No native ECM owner observed during QS.",
        )
    before, after = active[0]["native_threads"], active[-1]["native_threads"]
    old = {r["identity"]: r for r in before["threads"]}
    new = {r["identity"]: r for r in after["threads"]}
    if old.keys() != new.keys():
        raise RuntimeError("native identities changed during QS")
    deltas = [new[key]["cpu"] - old[key]["cpu"] for key in old]
    if any(not -0.000001 <= value <= 0.001 for value in deltas):
        raise RuntimeError("native worker CPU advanced during QS")
    for sample in active[1:]:
        workers = sample["native_threads"]["activity"]["workers"]
        first = before["activity"]["workers"]
        if workers != first:
            raise RuntimeError("native ECM task counters changed during QS")
    duration = after["at"] - before["at"]
    return dict(
        reached_native_QS=True,
        elapsed_observation=duration,
        native_cpu_deltas=deltas,
        observed_run_states=[
            [r["run_state"] for r in s["native_threads"]["threads"]]
            for s in active
        ],
        sustained_interval=duration >= 3,
        steady_snapshot=any(s["label"] == "qs-steady" for s in active),
        no_native_slice_or_curve_work_during_QS=True,
        final_threads_joined=state["lifecycle"][1],
        contexts_released=True,
        scope="Mach CPU uses a 1 ms per-worker refusal bound; native task "
        "counts are unchanged through QS entry, work and completion. "
        "Mapped code/allocator RSS is not a live-context leak.",
    )
