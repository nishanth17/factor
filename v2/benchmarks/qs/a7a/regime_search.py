"""Finite QS candidate neighborhoods; rankings never route by known factors."""

import copy
import math
import statistics

from .combined import FEATURES, cell

BANDS = (40, 50, 55, 60, 62, 65, 70, 75, 80)
COUNTS = (1, 2, 4, 6, 8, 10)


def candidate(count, packet=64, *, features=None, collector=None, name=None):
    """Keep numerical packet width and compatible features fully explicit."""
    if type(count) is not int or not 0 <= count <= 14:
        raise ValueError("invalid worker count")
    if packet not in (4, 8, 16, 32, 64):
        raise ValueError("invalid packet width")
    if not count:
        return dict(
            name="serial",
            arm="S0",
            polynomials_per_task=packet,
            features={},
            collector_overrides={} if collector is None else dict(collector),
        )
    options = dict(
        FEATURES,
        private_stop=True,
        stop_backend="semaphore",
        cpu_snapshot=False,
        coefficient_cache=False,
        whole_family=packet == 64,
    )
    if features:
        options.update(features)
    # A lease removal necessarily uses the ordinary synchronized stop path.
    # Keep that coupled removal visible instead of pretending it isolates L.
    if not options["lease_budget"]:
        options.update(
            private_stop=False, stop_backend="event", cpu_snapshot=False
        )
    if not options["private_stop"]:
        options["stop_backend"] = "event"
    return dict(
        name=name or f"w{count}_p{packet}",
        arm=f"C{count}",
        polynomials_per_task=packet,
        features=options,
        collector_overrides={} if collector is None else dict(collector),
    )


def coarse(*, features=None, collector=None, gray_limit=64):
    """Two packet scales per count, plus the direct serial control."""
    if gray_limit not in (16, 32, 64):
        raise ValueError("a separately admitted family geometry is required")
    packets = (min(16, gray_limit // 2), gray_limit)
    return [candidate(0, collector=collector)] + [
        candidate(
            count,
            packet,
            features=dict(features or {}, whole_family=packet == gray_limit),
            collector=collector,
        )
        for count in COUNTS
        for packet in packets
    ]


def refinements(selected):
    """Check adjacent counts/packet widths and bounded lookahead choices."""
    if selected["arm"] == "S0":
        return []
    count = int(selected["arm"][1:])
    packet = selected["polynomials_per_task"]
    packets = (4, 8, 16, 32, 64)
    location = packets.index(packet)
    choices = []
    for neighbor in (count - 1, count + 1):
        if 1 <= neighbor <= 10:
            options = dict(selected["features"])
            options["lookahead"] = min(options["lookahead"], 4 * neighbor)
            choices.append(
                candidate(
                    neighbor,
                    packet,
                    features=options,
                    collector=selected["collector_overrides"],
                )
            )
    for neighbor in packets[max(0, location - 1) : location + 2]:
        if neighbor != packet:
            options = dict(
                selected["features"],
                whole_family=neighbor == 64
                and selected["features"]["whole_family"],
            )
            choices.append(
                candidate(
                    count,
                    neighbor,
                    features=options,
                    collector=selected["collector_overrides"],
                )
            )
    for ahead in sorted({0, 2, 4, min(8, 4 * count)}):
        if ahead != selected["features"]["lookahead"]:
            options = dict(selected["features"], lookahead=ahead)
            choices.append(
                candidate(
                    count,
                    packet,
                    features=options,
                    collector=selected["collector_overrides"],
                    name=f"lookahead_{ahead}",
                )
            )
    return choices


def ablations(selected, marking):
    """Compatible removals, hot-loop alternatives and named interactions."""
    from ....qs.collector_kernels import KERNELS

    if selected["arm"] == "S0":
        return []
    count = int(selected["arm"][1:])
    packet = selected["polynomials_per_task"]
    base = selected["features"]
    choices = []
    alternatives = [
        ("without_lease_and_private_stop", dict(lease_budget=False)),
        ("without_state_reuse", dict(reuse_state=False)),
        ("without_packed_transport", dict(packed_transport=False)),
        ("verify_parent", dict(verification="parent")),
        (
            "verify_batch64",
            dict(verification="batch", verification_batch_atoms=64),
        ),
        (
            "verify_batch256",
            dict(verification="batch", verification_batch_atoms=256),
        ),
        ("without_lookahead_family", dict(lookahead=0, whole_family=False)),
        ("private_event", dict(stop_backend="event")),
        ("shared_event", dict(private_stop=False, stop_backend="event")),
        ("cpu_snapshot", dict(cpu_snapshot=True)),
        (
            "toggle_coefficient_cache",
            dict(coefficient_cache=not base.get("coefficient_cache", False)),
        ),
        ("arrival_admission", dict(arrival_order=True)),
        (
            "interaction_no_reuse_parent_verify",
            dict(reuse_state=False, verification="parent"),
        ),
        (
            "interaction_no_transport_parent_verify",
            dict(packed_transport=False, verification="parent"),
        ),
        (
            "interaction_no_lease_no_schedule",
            dict(lease_budget=False, lookahead=0, whole_family=False),
        ),
    ]
    for name, changes in alternatives:
        choices.append(
            candidate(
                count,
                packet,
                features=dict(base, **changes),
                collector=selected["collector_overrides"],
                name=name,
            )
        )
    pieces = marking.split("+")
    collector_alternatives = [
        ("without_sparse_queue", {"sparse_queue_bytes": 0}),
        ("without_power_plan", {"power_plan_bytes": 0}),
    ]
    for piece in (
        "presieve",
        "unroll4",
        "scan_fused",
        "packed_queue",
        "wide_hits",
        "wide_inverses",
        "wide_transfers",
        "pruned_transfers",
    ):
        altered = [p for p in pieces if p != piece]
        if piece not in pieces:
            altered = [p for p in KERNELS if p in (*pieces, piece)]
        collector_alternatives.append(
            (f"toggle_{piece}", {"marking": "+".join(altered) or "sparse"})
        )
    for name, changes in collector_alternatives:
        choices.append(
            candidate(
                count,
                packet,
                features=base,
                collector=dict(selected["collector_overrides"], **changes),
                name=name,
            )
        )
    return choices


def cells(fixture, candidates, *, repetition=1, endpoint=None, samples=1):
    """Freeze equal grants and numerical work before measuring a band."""
    digits = int(fixture["digits"])
    if digits not in BANDS:
        raise ValueError(
            "band requires a separately admitted feasibility plan"
        )
    endpoint = endpoint or ("complete" if digits <= 55 else "prefix")
    if endpoint not in ("prefix", "complete"):
        raise ValueError("invalid optimization endpoint")
    quota = 640 if digits >= 75 else 1280
    rows = []
    for item in candidates:
        values = copy.deepcopy(item)
        name = values.pop("name")
        arm = values.pop("arm")
        row = cell(
            arm,
            fixture["id"],
            stage="M2-regime",
            wall=120,
            cpu=1200,
            work=96 * 10**12,
            memory=22 * 2**30,
            storage=8 * 2**20,
            endpoint=endpoint,
            lifetime="warm",
            samples=samples,
            polynomials=quota,
            max_batch_atoms=2048,
            checkpoint_bytes=64 * 2**20,
            qs_memory_bytes=4 * 2**30,
            variant=name,
            rep=repetition,
            **values,
        )
        rows.append(row)
    return rows


def shortlist(rows, *, target=3, relative_band=0.15):
    """Keep near ties and serial; exclude incomplete factor runs."""
    if type(target) is not int or not 1 <= target <= 16:
        raise ValueError("finite shortlist target required")
    if not 0 <= relative_band <= 1:
        raise ValueError("invalid exploration band")
    endpoints = {row["endpoint"] for row in rows if row.get("status") == "ok"}
    if len(endpoints) > 1:
        raise ValueError("different endpoints cannot share a ranking")
    grouped = {}
    for row in rows:
        if row.get("status") != "ok":
            continue
        value = row["seconds"]
        if (
            type(value) not in (float, int)
            or not math.isfinite(value)
            or value <= 0
        ):
            raise ValueError(
                "finite positive complete-cost measurement required"
            )
        grouped.setdefault(row["variant"], []).append(row)
    ranked = []
    for name, observations in grouped.items():
        endpoint = {r["endpoint"] for r in observations}
        if len(endpoint) != 1:
            raise ValueError("different endpoints cannot share a ranking")
        if "complete" in endpoint and not all(
            r["complete"] for r in observations
        ):
            continue
        if "prefix" in endpoint and not all(
            r["coverage_valid"] for r in observations
        ):
            continue
        ranked.append(
            dict(
                variant=name,
                seconds=statistics.median(r["seconds"] for r in observations),
                candidate=observations[0]["candidate"],
                n=len(observations),
            )
        )
    ranked.sort(key=lambda r: (r["seconds"], r["variant"]))
    if not ranked:
        return []
    result = [
        r
        for i, r in enumerate(ranked)
        if i < target
        or r["seconds"] <= ranked[0]["seconds"] * (1 + relative_band)
        or r["candidate"]["arm"] == "S0"
    ]
    return result
